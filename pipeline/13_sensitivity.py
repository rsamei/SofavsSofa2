"""Step 13: further sensitivity analyses. Only aggregates leave BigQuery; counts of 1-4 are suppressed.

1. Component AUROCs: each day-1 organ component of SOFA and SOFA-2 alone (24-hour maxima at the end of the day-1
   window), with the paired SOFA-2 minus SOFA difference per organ, primary set.
2. Sepsis-3 with a within-window baseline: in each suspicion window the baseline is the score at the window's first
   evaluated hour, and Sepsis-3 is met when a later evaluated hour of that window reaches baseline + 2. Rates,
   paired difference, kappa, discordant counts and their hospital mortality, primary set.
3. Classification (any window, baseline 0) in the set whose first suspicion fell by the end of ICU hour 23.
4. Alternative outcomes by concordance group (any window, baseline 0): death in hospital or discharge to hospice, and
   death within 30 days of ICU admission (patients.dod); crude mortality with Wilson CIs and crude OR vs Neither.
"""
import json
import math
import sys

from bq_client import BUILD, OUT, PROJECT, BQError, query
from delong import run_delong

MIN_CELL = 5
PRIM = "first_elig_stay AND NOT not_assessable_any"
GRP = ("CASE WHEN any_s1_std AND any_s2_std THEN 'BOTH' WHEN any_s2_std THEN 'SOFA2_ONLY' "
       "WHEN any_s1_std THEN 'SOFA_ONLY' ELSE 'NEITHER' END")
SUPPRESSED = []

COMP_SOURCE = f"""(SELECT c.stay_id, c.hospital_expire_flag, c.first_elig_stay, c.not_assessable_any, z.*
  EXCEPT (stay_id) FROM `{BUILD}.cohort` c JOIN `{BUILD}.hz_scores` z
  ON z.stay_id = c.stay_id AND z.hz = 'D1' AND z.hr = c.d1_last_hr) q"""
ORGANS = {"respiration": ("s1_resp24", "s2_resp24"), "coagulation": ("s1_coag24", "s2_hemo24"),
          "liver": ("s1_liver24", "s2_liver24"), "cardiovascular": ("s1_cv24", "s2_cv24"),
          "brain": ("s1_cns24", "s2_brain24"), "kidney": ("s1_renal24", "s2_kidney24")}

DELTA_SQL = f"""
WITH w AS (
  SELECT stay_id, hz_key, hr, s1_std, s2_std,
    FIRST_VALUE(s1_std) OVER (PARTITION BY stay_id, hz_key ORDER BY hr) AS b1,
    FIRST_VALUE(s2_std) OVER (PARTITION BY stay_id, hz_key ORDER BY hr) AS b2,
    ROW_NUMBER() OVER (PARTITION BY stay_id, hz_key ORDER BY hr) AS k
  FROM `{BUILD}.hz_scores` WHERE hz = 'W' AND evaluate),
s AS (SELECT stay_id, LOGICAL_OR(k > 1 AND s1_std >= b1 + 2) AS d1m, LOGICAL_OR(k > 1 AND s2_std >= b2 + 2) AS d2m
      FROM w GROUP BY stay_id),
p AS (SELECT COALESCE(s.d1m, FALSE) AS a, COALESCE(s.d2m, FALSE) AS b, c.hospital_expire_flag AS d
      FROM `{BUILD}.cohort` c LEFT JOIN s USING (stay_id) WHERE {PRIM})
SELECT COUNT(*) n, COUNTIF(a AND b) n11, COUNTIF(a AND NOT b) n10, COUNTIF(NOT a AND b) n01, COUNTIF(NOT a AND NOT b) n00,
  SUM(IF(a AND NOT b, d, 0)) d10, SUM(IF(NOT a AND b, d, 0)) d01
FROM p
"""

SUSP23_SQL = f"""
SELECT COUNT(*) n, COUNTIF(any_s1_std AND any_s2_std) n11, COUNTIF(any_s1_std AND NOT any_s2_std) n10,
  COUNTIF(NOT any_s1_std AND any_s2_std) n01, COUNTIF(NOT any_s1_std AND NOT any_s2_std) n00,
  SUM(IF(any_s1_std AND NOT any_s2_std, hospital_expire_flag, 0)) d10,
  SUM(IF(NOT any_s1_std AND any_s2_std, hospital_expire_flag, 0)) d01
FROM `{BUILD}.cohort` WHERE {PRIM} AND first_susp_hr <= 23
"""

OUTCOME_SQL = f"""
SELECT {GRP} grp, COUNT(*) n, SUM(death_or_hospice) dh, SUM(death_30d) d30, SUM(hospital_expire_flag) dhosp
FROM `{BUILD}.cohort` WHERE {PRIM} GROUP BY 1 ORDER BY 1
"""


def run(sql):
    f, rows, _ = query(PROJECT, sql)
    return [dict(zip(f, r)) for r in rows]


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [p, max(0.0, c - h), c + h]


def cnt(v, label):
    v = int(v)
    if 0 < v < MIN_CELL:
        SUPPRESSED.append(label)
        return "<5"
    return v


def paired(r, label):
    """Rates, paired difference (Wald), kappa, discordant counts and mortality from a 2x2 (a = SOFA, b = SOFA-2)."""
    n, n11, n10, n01, n00 = (int(r[k]) for k in ("n", "n11", "n10", "n01", "n00"))
    d = (n01 - n10) / n
    se = math.sqrt(((n01 + n10) / n - d * d) / n)
    po = (n11 + n00) / n
    pe = ((n11 + n10) * (n11 + n01) + (n00 + n10) * (n00 + n01)) / n / n
    out = {"n": n, "rate_sofa": wilson(n11 + n10, n), "rate_sofa2": wilson(n11 + n01, n),
           "difference_sofa2_minus_sofa": [d, d - 1.96 * se, d + 1.96 * se], "kappa": (po - pe) / (1 - pe)}
    if 0 < n10 < MIN_CELL or 0 < n01 < MIN_CELL:
        out.update({"sofa_only": "withheld", "sofa2_only": "withheld"})
        SUPPRESSED.append(f"{label} discordant pair")
        return out
    out.update({"sofa_only": n10, "sofa2_only": n01})
    for lab, k, dk in (("sofa_only", n10, int(r["d10"])), ("sofa2_only", n01, int(r["d01"]))):
        if dk >= MIN_CELL and k - dk >= MIN_CELL:
            out[f"{lab}_hosp_mortality"] = wilson(dk, k)
        else:
            SUPPRESSED.append(f"{label} {lab} deaths")
    return out


def main():
    res = {}
    try:
        sys_expr, contrasts = {}, {}
        for organ, (a, b) in ORGANS.items():
            sys_expr[f"sofa_{organ}"], sys_expr[f"sofa2_{organ}"] = a, b
            contrasts[f"{organ}_sofa2_minus_sofa"] = {f"sofa2_{organ}": 1, f"sofa_{organ}": -1}
        res["component_auroc"] = run_delong(COMP_SOURCE, {"d1_hosp": PRIM}, sys_expr, contrasts)
        res["within_window_baseline"] = paired(run(DELTA_SQL)[0], "within-window")
        res["susp_by_hr23_classification"] = paired(run(SUSP23_SQL)[0], "susp_by_hr23")
        rows = {r["grp"]: r for r in run(OUTCOME_SQL)}
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    res["alternative_outcomes"] = {}
    ne = rows["NEITHER"]
    for key, lab in (("dh", "death_or_hospice"), ("d30", "death_30d"), ("dhosp", "hospital_death")):
        out = {}
        n0, k0 = int(ne["n"]), int(ne[key])
        for g, r in rows.items():
            n, k = int(r["n"]), int(r[key])
            if 0 < k < MIN_CELL or 0 < n - k < MIN_CELL:
                out[g] = {"n": n, "events": "<5"}
                SUPPRESSED.append(f"{lab} {g}")
                continue
            row = {"n": n, "events": k, "risk": wilson(k, n)}
            if g != "NEITHER" and k and n0 - k0 and k0 and n - k:
                lor = math.log((k / (n - k)) / (k0 / (n0 - k0)))
                se = math.sqrt(1 / k + 1 / (n - k) + 1 / k0 + 1 / (n0 - k0))
                row["crude_or_vs_neither"] = [math.exp(lor), math.exp(lor - 1.96 * se), math.exp(lor + 1.96 * se)]
            out[g] = row
        res["alternative_outcomes"][lab] = out
    res["suppressed"] = SUPPRESSED
    (OUT / "sensitivity.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"within_window": {k: res["within_window_baseline"][k] for k in ("rate_sofa", "rate_sofa2")},
                      "components": {r["contrast"]: round(r["est"], 3) for r in res["component_auroc"]
                                     if r["contrast"].endswith("sofa2_minus_sofa")}}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
