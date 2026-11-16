"""Step 17 (Supplementary Methods S1): ICU-onset sepsis measured from an observed ICU baseline. Reads
BUILD_DATASET.s1_cohort (step 16). Descriptive; added after the main results had been seen.

  1. classification (baseline = lowest total over ICU hours 23-47; sensitivity: the total at hour 47): 2x2, rates
     (Wilson 95% CI), paired rate difference, Cohen's kappa and the prevalence- and bias-adjusted kappa (PABAK), both
     with Poisson-bootstrap 95% CIs as in step 5 (1,000 replicates in BigQuery), McNemar with continuity correction;
  2. concordance groups (Both, SOFA-2 only, SOFA only, Neither): hospital and ICU mortality counted from the landmark
     (Wilson CI), crude OR (Woolf CI) and Mantel-Haenszel OR (strata: age decade [top >= 80] x sex x cardiac-vascular
     ICU; Robins-Breslow-Greenland CI) against Neither;
  3. group descriptives: medians [IQR] (exact percentiles) of age, baselines, onset from suspicion, ICU and hospital
     length of stay after the landmark in hospital survivors; new organ support after the suspicion time (invasive
     ventilation, vasopressors, renal replacement therapy), with no support of the same type in the 24 hours before
     that start;
  4. drivers of discordance: the rise at onset by component for the detecting score, and for discordant stays the
     non-detecting score's rise at the detecting score's onset hour, as points short of 2;
  5. secondary: the largest rise from baseline within the window, paired DeLong AUROC for in-hospital mortality.
Output control: counts 1-4 are suppressed, and estimates resting on them are omitted.
"""
import importlib.util
import json
import math
import sys
from pathlib import Path

from bq_client import BUILD, OUT, PROJECT, BQError, query
from delong import run_delong

_spec = importlib.util.spec_from_file_location("classification", Path(__file__).with_name("05_classification.py"))
_c05 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_c05)
wilson, chi2_sf, B = _c05.wilson, _c05.chi2_sf, _c05.B

C = f"{BUILD}.s1_cohort"
MIN_CELL = 5
SET = "NOT not_assessable"
COMPONENTS = ["resp", "coag", "liver", "cv", "brain", "kidney"]
SUPPRESSED = []


def grp(p1, p2):
    return (f"CASE WHEN {p1} AND {p2} THEN 'BOTH' WHEN {p2} THEN 'SOFA2_ONLY' "
            f"WHEN {p1} THEN 'SOFA_ONLY' ELSE 'NEITHER' END")


def run(sql):
    fields, rows, _ = query(PROJECT, sql)
    return [dict(zip(fields, r)) for r in rows]


def cnt(v, label):
    v = int(v)
    if 0 < v < MIN_CELL:
        SUPPRESSED.append(label)
        return "<5"
    return v


def rate(k, n, label):
    k, n = int(k), int(n)
    if 0 < k < MIN_CELL or 0 < n - k < MIN_CELL:
        SUPPRESSED.append(label)
        return {"n": cnt(n, label + " n"), "k": "<5" if 0 < k < MIN_CELL else "withheld (complement < 5)"}
    return {"n": n, "k": k, "rate": wilson(k, n)}


def classification(p1, p2, label):
    t = run(f"""SELECT COUNTIF({p1} AND {p2}) n11, COUNTIF({p1} AND NOT {p2}) n10,
                       COUNTIF(NOT {p1} AND {p2}) n01, COUNTIF(NOT {p1} AND NOT {p2}) n00 FROM `{C}` WHERE {SET}""")[0]
    n11, n10, n01, n00 = (int(t[k]) for k in ("n11", "n10", "n01", "n00"))
    n = n11 + n10 + n01 + n00
    po = (n11 + n00) / n
    pe = ((n11 + n10) * (n11 + n01) + (n00 + n10) * (n00 + n01)) / n / n
    kappa = (po - pe) / (1 - pe)
    # Poisson bootstrap with the weights and seeds of step 5 (agreement in the main analysis)
    boot = run(f"""
      WITH base AS (SELECT stay_id, {p1} a, {p2} b FROM `{C}` WHERE {SET}),
      w AS (
        SELECT r, a, b,
          (SELECT COUNTIF(u > th) FROM UNNEST([0.367879, 0.735759, 0.919699, 0.981012, 0.996340, 0.999406, 0.999917]) th) AS wt
        FROM base, UNNEST(GENERATE_ARRAY(1, {B})) r,
          UNNEST([MOD(ABS(FARM_FINGERPRINT(CONCAT(CAST(stay_id AS STRING), '-', CAST(r AS STRING)))), 1000000000) / 1e9]) u),
      c AS (SELECT r, SUM(IF(a AND b, wt, 0)) n11, SUM(IF(a AND NOT b, wt, 0)) n10,
                   SUM(IF(NOT a AND b, wt, 0)) n01, SUM(IF(NOT a AND NOT b, wt, 0)) n00 FROM w GROUP BY r),
      k AS (SELECT r, (n11 + n00) / (n11 + n10 + n01 + n00) po,
              ((n11 + n10) * (n11 + n01) + (n00 + n10) * (n00 + n01)) / POW(n11 + n10 + n01 + n00, 2) pe FROM c)
      SELECT APPROX_QUANTILES((po - pe) / (1 - pe), 1000)[OFFSET(25)] k_lo, APPROX_QUANTILES((po - pe) / (1 - pe), 1000)[OFFSET(975)] k_hi,
             APPROX_QUANTILES(2 * po - 1, 1000)[OFFSET(25)] p_lo, APPROX_QUANTILES(2 * po - 1, 1000)[OFFSET(975)] p_hi FROM k""")[0]
    d = (n01 - n10) / n
    se_d = math.sqrt((n10 + n01) - (n01 - n10) ** 2 / n) / n
    chi2 = (abs(n10 - n01) - 1) ** 2 / (n10 + n01) if n10 + n01 else 0.0
    return {
        "n": n,
        "cells": {k: cnt(v, f"{label} {k}") for k, v in
                  (("both", n11), ("sofa_only", n10), ("sofa2_only", n01), ("neither", n00))},
        "rate_sofa": wilson(n11 + n10, n), "rate_sofa2": wilson(n11 + n01, n),
        "paired_rate_difference_sofa2_minus_sofa": [d, d - 1.96 * se_d, d + 1.96 * se_d],
        "kappa": [kappa, float(boot["k_lo"]), float(boot["k_hi"])],
        "pabak": [2 * po - 1, float(boot["p_lo"]), float(boot["p_hi"])],
        "bootstrap_reps": B,
        "mcnemar_chi2_cc": round(chi2, 2), "mcnemar_p": float(f"{chi2_sf(chi2, 1):.3g}") if chi2 else 1.0,
    }


def odds_ratios(g_expr, label):
    strata = run(f"""SELECT {g_expr} grp, CAST(FLOOR(LEAST(age, 89) / 10) AS INT64) ab, gender, cvicu,
                            COUNT(*) n, SUM(hospital_expire_flag) dh, SUM(icu_death) di
                     FROM `{C}` WHERE {SET} GROUP BY 1, 2, 3, 4""")
    tot = {}
    for s in strata:
        g = tot.setdefault(s["grp"], {"n": 0, "dh": 0, "di": 0})
        for k in ("n", "dh", "di"):
            g[k] += int(s[k])
    groups = {g: {"n": cnt(v["n"], f"{label} {g} n"),
                  "hosp": rate(v["dh"], v["n"], f"{label} {g} hospital deaths"),
                  "icu": rate(v["di"], v["n"], f"{label} {g} ICU deaths")} for g, v in tot.items()}
    ors = {}
    for g in ("BOTH", "SOFA2_ONLY", "SOFA_ONLY"):
        for o, k in (("hosp", "dh"), ("icu", "di")):
            key = f"{g}|{o}"
            if g not in tot or "NEITHER" not in tot:
                continue
            a, n1, c, n0 = tot[g][k], tot[g]["n"], tot["NEITHER"][k], tot["NEITHER"]["n"]
            if min(a, c, n1 - a, n0 - c) < MIN_CELL:
                ors[key] = "omitted (< 5 in a cell)"
                SUPPRESSED.append(f"{label} OR {key}")
                continue
            crude = a * (n0 - c) / ((n1 - a) * c)
            se = math.sqrt(1 / a + 1 / (n1 - a) + 1 / c + 1 / (n0 - c))
            R = S = Pp = Q = T = 0.0
            for st in {(s["ab"], s["gender"], s["cvicu"]) for s in strata}:
                x = next((s for s in strata if s["grp"] == g and (s["ab"], s["gender"], s["cvicu"]) == st), None)
                r = next((s for s in strata if s["grp"] == "NEITHER" and (s["ab"], s["gender"], s["cvicu"]) == st), None)
                if not x or not r:
                    continue
                aa, bb = int(x[k]), int(x["n"]) - int(x[k])
                cc, dd = int(r[k]), int(r["n"]) - int(r[k])
                nn = aa + bb + cc + dd
                R += aa * dd / nn
                S += bb * cc / nn
                Pp += (aa + dd) / nn * aa * dd / nn
                Q += (aa + dd) / nn * bb * cc / nn + (bb + cc) / nn * aa * dd / nn
                T += (bb + cc) / nn * bb * cc / nn
            row = {"crude": [crude, math.exp(math.log(crude) - 1.96 * se), math.exp(math.log(crude) + 1.96 * se)]}
            if R > 0 and S > 0:
                mh = R / S
                se_mh = math.sqrt(Pp / (2 * R * R) + Q / (2 * R * S) + T / (2 * S * S))
                row["mantel_haenszel"] = [mh, math.exp(math.log(mh) - 1.96 * se_mh), math.exp(math.log(mh) + 1.96 * se_mh)]
            row["label"] = "descriptive association, not causal"
            ors[key] = row
    return groups, ors


def descriptives(g_expr, label):
    cont = {"age": "age", "sofa_baseline": "s1_base", "sofa2_baseline": "s2_base",
            "sofa_onset_from_suspicion_hr": "s1_onset_from_susp_hr", "sofa2_onset_from_suspicion_hr": "s2_onset_from_susp_hr",
            "icu_los_after_landmark_days_survivors": "IF(hospital_expire_flag = 0, icu_los_after_landmark_days, NULL)",
            "hosp_los_after_landmark_days_survivors": "IF(hospital_expire_flag = 0, hosp_los_after_landmark_days, NULL)"}
    q = ",\n".join(
        f"PERCENTILE_CONT({e}, {p}) OVER (PARTITION BY g) AS {k}__{int(p * 100)}"
        for k, e in cont.items() for p in (0.25, 0.5, 0.75))
    cn = ",\n".join(f"COUNT({e}) OVER (PARTITION BY g) AS {k}__n" for k, e in cont.items())
    rows = run(f"""SELECT DISTINCT g, COUNT(*) OVER (PARTITION BY g) AS n,
                     COUNTIF(new_vent) OVER (PARTITION BY g) AS new_vent, COUNTIF(new_vasopressor) OVER (PARTITION BY g) AS new_vasopressor,
                     COUNTIF(new_rrt) OVER (PARTITION BY g) AS new_rrt, COUNTIF(gender = 'F') OVER (PARTITION BY g) AS female,
                     COUNTIF(cvicu) OVER (PARTITION BY g) AS cvicu,
                     {q},
                     {cn}
                   FROM (SELECT *, {g_expr} g FROM `{C}` WHERE {SET})""")
    out = {}
    for r in rows:
        g, n = r["g"], int(r["n"])
        row = {"n": cnt(n, f"{label} {g} n")}
        if n >= MIN_CELL:
            for k in cont:
                m = int(r[f"{k}__n"])
                if m >= MIN_CELL:
                    row[k] = {"n": m, "median_iqr": [float(r[f"{k}__50"]), float(r[f"{k}__25"]), float(r[f"{k}__75"])]}
            for k in ("female", "cvicu", "new_vent", "new_vasopressor", "new_rrt"):
                row[k] = rate(r[k], n, f"{label} {g} {k}")
        out[g] = row
    return out


def drivers():
    res = {}
    for g, scores in (("BOTH", ("s1", "s2")), ("SOFA2_ONLY", ("s2",)), ("SOFA_ONLY", ("s1",))):
        cond = {"BOTH": "s1_pos AND s2_pos", "SOFA2_ONLY": "s2_pos AND NOT s1_pos", "SOFA_ONLY": "s1_pos AND NOT s2_pos"}[g]
        for s in scores:
            sel = ",\n".join(
                f"AVG({s}_d_{c}) mean_{c}, COUNTIF({s}_d_{c} >= 1) up_{c}, COUNTIF({s}_d_{c} <= -1) down_{c}"
                for c in COMPONENTS)
            r = run(f"""SELECT COUNT(*) n, ANY_VALUE(med_d) med_d, AVG({s}_d) mean_d, {sel}
                        FROM (SELECT *, PERCENTILE_CONT({s}_d, 0.5) OVER () med_d
                              FROM `{C}` WHERE {SET} AND {cond})""")[0]
            n = int(r["n"])
            row = {"n": cnt(n, f"drivers {g} {s} n")}
            if n >= MIN_CELL:
                row["rise_at_onset"] = {"median": float(r["med_d"]), "mean": float(r["mean_d"])}
                row["components"] = {c: {"mean_rise": float(r[f"mean_{c}"]),
                                         "rose_ge1": rate(r[f"up_{c}"], n, f"drivers {g} {s} {c} rose"),
                                         "fell_ge1": rate(r[f"down_{c}"], n, f"drivers {g} {s} {c} fell")}
                                     for c in COMPONENTS}
            res[f"{g}|{'sofa' if s == 's1' else 'sofa2'}"] = row
    short = {}
    for g, col, cond in (("SOFA2_ONLY", "s1_d_at_s2_onset", "s2_pos AND NOT s1_pos"),
                         ("SOFA_ONLY", "s2_d_at_s1_onset", "s1_pos AND NOT s2_pos")):
        rows = run(f"""SELECT CASE WHEN {col} >= 1 THEN '1' WHEN {col} = 0 THEN '2' WHEN {col} = -1 THEN '3'
                                   ELSE '>=4' END AS short_of_2, COUNT(*) n
                       FROM `{C}` WHERE {SET} AND {cond} GROUP BY 1 ORDER BY 1""")
        short[g] = {r["short_of_2"]: cnt(r["n"], f"shortfall {g} {r['short_of_2']}") for r in rows}
    res["non_detecting_score_points_short_of_2"] = short
    return res


def main():
    res = {}
    try:
        res["classification"] = classification("s1_pos", "s2_pos", "primary")
        res["classification_hour47_baseline"] = classification("s1_pos47", "s2_pos47", "hour-47 baseline")
        res["groups"], res["odds_ratios_vs_neither"] = odds_ratios(grp("s1_pos", "s2_pos"), "primary")
        res["groups_hour47_baseline"], res["odds_ratios_vs_neither_hour47_baseline"] = odds_ratios(
            grp("s1_pos47", "s2_pos47"), "hour-47 baseline")
        res["descriptives"] = descriptives(grp("s1_pos", "s2_pos"), "descriptives")
        res["drivers"] = drivers()
        res["largest_rise_discrimination"] = run_delong(
            f"`{C}`", {"primary": SET},
            {"s1_rise": "s1_max_rise", "s2_rise": "s2_max_rise", "s1_rise47": "s1_max_rise47", "s2_rise47": "s2_max_rise47"},
            {"s2_minus_s1": {"s2_rise": 1, "s1_rise": -1}, "s2_minus_s1_hour47": {"s2_rise47": 1, "s1_rise47": -1}})
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    floor_ok = all(r["n_pos"] >= 100 and r["n_neg"] >= 100 for r in res["largest_rise_discrimination"])
    res["event_floor_ok"] = floor_ok
    res["suppressed"] = SUPPRESSED
    (OUT / "icu_onset.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    c = res["classification"]
    print(json.dumps({"n": c["n"], "cells": c["cells"], "rate_sofa": c["rate_sofa"], "rate_sofa2": c["rate_sofa2"],
                      "kappa": c["kappa"], "pabak": c["pabak"], "event_floor_ok": floor_ok, "suppressed": SUPPRESSED}))
    return 0 if floor_ok else 1


if __name__ == "__main__":
    sys.exit(main())
