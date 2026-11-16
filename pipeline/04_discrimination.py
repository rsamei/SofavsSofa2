"""Step 4: day-1 discrimination. For each analysis set, the AUROC of every score and variant and the paired
contrasts, by DeLong's method with mid-rank placement values. For a contrast sum_k c_k * AUC_k, the per-subject
combined placement value sum_k c_k * V_k gives Var = Var(V10)/m + Var(V01)/n, the DeLong covariance applied to the
contrast, so two-score differences and four-AUROC contrasts use the same machinery.

Analysis sets (primary set = first eligible stays, excluding stays with no evaluated window hour): d1_hosp (primary),
d1_icu, fw_hosp (first-window maximum, descriptive), d1_hosp_noncvicu, d1_hosp_allstays, d1_hosp_inicu,
d1_hosp_landmark, d1_hosp_complete, d1_hosp_stayid (suspicions matched by stay identifier), d1_hosp_ts_consistent
(stays with inconsistent timestamps excluded), d1_hosp_full_day1 (stays with a full 24-hour day 1).
Only aggregates leave BigQuery (every set has >= 100 deaths and >= 100 survivors, checked below)."""
import json
import math
import sys

from bq_client import BUILD, OUT, PROJECT, BQError, query

COHORT = f"{BUILD}.cohort"

SYSTEMS = ["s1_std", "s1_gcsalt", "s1_kcr", "s2_std", "s2_e", "s2_gcs15", "s2_kcr", "s2_kuo22", "s2_g", "s2_s", "s2_w", "s2_h", "s2_i"]

BASE = "first_elig_stay AND NOT not_assessable_any"
SETS = {
    "d1_hosp": ("hospital_expire_flag", "d1", BASE),
    "d1_icu": ("icu_death", "d1", BASE),
    "fw_hosp": ("hospital_expire_flag", "fw", "first_elig_stay AND NOT not_assessable_first"),
    "d1_hosp_noncvicu": ("hospital_expire_flag", "d1", BASE + " AND NOT cvicu"),
    "d1_hosp_allstays": ("hospital_expire_flag", "d1", "NOT not_assessable_any"),
    "d1_hosp_inicu": ("hospital_expire_flag", "d1", BASE + " AND first_susp_hr >= 0"),
    "d1_hosp_landmark": ("hospital_expire_flag", "d1", BASE + " AND landmark"),
    "d1_hosp_complete": ("hospital_expire_flag", "d1", BASE + " AND complete_case"),
    "d1_hosp_stayid": ("hospital_expire_flag", "d1", "first_elig_stay_old AND NOT not_assessable_old"),
    "d1_hosp_ts_consistent": ("hospital_expire_flag", "d1", BASE + " AND NOT ts_inconsistent"),
    "d1_hosp_susp_by_hr23": ("hospital_expire_flag", "d1", BASE + " AND first_susp_hr <= 23"),
    "d1_hosp_resp_measured": ("hospital_expire_flag", "d1", BASE + " AND resp_measured_both_d1"),
    "d1_death_or_hospice": ("death_or_hospice", "d1", BASE),
    "d1_death_30d": ("death_30d", "d1", BASE),
    "d1_hosp_full_day1": ("hospital_expire_flag", "d1", BASE + " AND d1_last_hr = 23"),
}

# contrasts: name -> coefficients over systems (AUC scale)
PRIMARY = {"s2_std": 1, "s1_std": -1}
FULL = {
    "primary_s2_minus_s1": PRIMARY,
    "matched_gcs_a_s2_minus_s1gcsalt": {"s2_std": 1, "s1_gcsalt": -1},
    "matched_gcs_b_s2gcs15_minus_s1": {"s2_gcs15": 1, "s1_std": -1},
    "h3a_s1gcsalt_minus_s1": {"s1_gcsalt": 1, "s1_std": -1},
    "h3b_s2_minus_s2gcs15": {"s2_std": 1, "s2_gcs15": -1},
    "matched_kidney_s2kcr_minus_s1kcr": {"s2_kcr": 1, "s1_kcr": -1},
    "uo_contribution": {"s2_std": 1, "s1_std": -1, "s2_kcr": -1, "s1_kcr": 1},
    "variant_c_s2kcr_minus_s1": {"s2_kcr": 1, "s1_std": -1},
    "variant_d_s2kuo22_minus_s1": {"s2_kuo22": 1, "s1_std": -1},
    "variant_e_s2e_minus_s1": {"s2_e": 1, "s1_std": -1},
    "direct_s2_minus_s2kcr": {"s2_std": 1, "s2_kcr": -1},
    "direct_s2_minus_s2kuo22": {"s2_std": 1, "s2_kuo22": -1},
    "direct_s2e_minus_s2": {"s2_e": 1, "s2_std": -1},
    "direct_s1kcr_minus_s1": {"s1_kcr": 1, "s1_std": -1},
    "variant_g_s2g_minus_s1": {"s2_g": 1, "s1_std": -1},
    "variant_h_s2h_minus_s1": {"s2_h": 1, "s1_std": -1},
    "variant_i_s2i_minus_s1": {"s2_i": 1, "s1_std": -1},
    "direct_s2g_minus_s2": {"s2_g": 1, "s2_std": -1},
    "variant_s_s2s_minus_s1": {"s2_s": 1, "s1_std": -1},
    "direct_s2s_minus_s2": {"s2_s": 1, "s2_std": -1},
    "direct_s2w_minus_s2": {"s2_w": 1, "s2_std": -1},
    "direct_s2h_minus_s2": {"s2_h": 1, "s2_std": -1},
    "direct_s2i_minus_s2": {"s2_i": 1, "s2_std": -1},
}


def build_sql():
    sets = ",\n      ".join(
        f"STRUCT('{an}' AS an, {y} AS y, ({where}) AS inc, '{kind}' AS kind)"
        for an, (y, kind, where) in SETS.items())
    systems = ",\n      ".join(
        f"STRUCT('{s}' AS sys, IF(st.kind = 'd1', d1_{s}, fw_{s}) AS x)" for s in SYSTEMS)
    pv_cols = ",\n    ".join(f"MAX(IF(sys = '{s}', v, NULL)) AS v_{s}" for s in SYSTEMS)
    defs = []
    for name, coef in [(f"auc_{s}", {s: 1}) for s in SYSTEMS] + list(FULL.items()):
        full_only = not (name.startswith("auc_") or name == "primary_s2_minus_s1")
        coefs = ", ".join(f"{float(coef.get(s, 0))} AS k_{s}" for s in SYSTEMS)
        defs.append(f"STRUCT('{name}' AS contrast, {str(full_only).upper()} AS full_only, {coefs})")
    defs_sql = ",\n      ".join(defs)
    combo = " + ".join(f"d.k_{s} * p.v_{s}" for s in SYSTEMS)
    return f"""
WITH long AS (
  SELECT st.an, sy.sys, stay_id, st.y, sy.x
  FROM `{COHORT}`,
    UNNEST([
      {sets}
    ]) AS st,
    UNNEST([
      {systems}
    ]) AS sy
  WHERE st.inc
),
cc AS (SELECT an, stay_id FROM long GROUP BY an, stay_id HAVING COUNTIF(x IS NULL) = 0),
r AS (
  SELECT l.an, l.sys, l.stay_id, l.y, l.x,
    RANK() OVER (PARTITION BY l.an, l.sys ORDER BY l.x) + (COUNT(*) OVER (PARTITION BY l.an, l.sys, l.x) - 1) / 2.0 AS pr,
    RANK() OVER (PARTITION BY l.an, l.sys, l.y ORDER BY l.x) + (COUNT(*) OVER (PARTITION BY l.an, l.sys, l.y, l.x) - 1) / 2.0 AS cr,
    COUNTIF(l.y = 1) OVER (PARTITION BY l.an, l.sys) AS m, COUNTIF(l.y = 0) OVER (PARTITION BY l.an, l.sys) AS n
  FROM long l JOIN cc USING (an, stay_id)
),
pl AS (SELECT an, sys, stay_id, y, IF(y = 1, (pr - cr) / n, 1 - (pr - cr) / m) AS v FROM r),
pv AS (
  SELECT an, stay_id, ANY_VALUE(y) AS y,
    {pv_cols}
  FROM pl GROUP BY an, stay_id
),
defs AS (SELECT * FROM UNNEST([
      {defs_sql}
    ])),
cv AS (
  SELECT p.an, d.contrast, p.y, {combo} AS c
  FROM pv p CROSS JOIN defs d
  WHERE NOT d.full_only OR p.an IN ('d1_hosp', 'd1_icu', 'd1_hosp_susp_by_hr23', 'd1_hosp_resp_measured')
)
SELECT an, contrast, COUNTIF(y = 1) AS n_pos, COUNTIF(y = 0) AS n_neg,
  AVG(IF(y = 1, c, NULL)) AS est,
  SQRT(VAR_SAMP(IF(y = 1, c, NULL)) / COUNTIF(y = 1) + VAR_SAMP(IF(y = 0, c, NULL)) / COUNTIF(y = 0)) AS se
FROM cv GROUP BY an, contrast ORDER BY an, contrast
"""


def main():
    try:
        fields, rows, job = query(PROJECT, build_sql())
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    out = []
    for r in rows:
        d = dict(zip(fields, r))
        est, se = float(d["est"]), float(d["se"])
        z = est / se if se > 0 else float("nan")
        rec = {"set": d["an"], "contrast": d["contrast"], "n_pos": int(d["n_pos"]), "n_neg": int(d["n_neg"]),
               "est": est, "se": se, "lo": est - 1.96 * se, "hi": est + 1.96 * se}
        if not d["contrast"].startswith("auc_"):
            rec["z"] = round(z, 2)
            rec["p"] = float(f"{math.erfc(abs(z) / math.sqrt(2)):.3g}")
        out.append(rec)
    by = {(o["set"], o["contrast"]): o for o in out}
    floor_ok = all(o["n_pos"] >= 100 and o["n_neg"] >= 100 for o in out)
    prim = by[("d1_hosp", "primary_s2_minus_s1")]
    h3a, h3b = by[("d1_hosp", "h3a_s1gcsalt_minus_s1")], by[("d1_hosp", "h3b_s2_minus_s2gcs15")]
    decisions = {
        "H2_sofa2_better": prim["lo"] > 0,
        "H3_depends_on_gcs_handling": h3a["lo"] > 0 and h3b["lo"] > 0,
        # descriptive decomposition shares (not estimands; frozen plan)
        "share_gcs_a": round(1 - by[("d1_hosp", "matched_gcs_a_s2_minus_s1gcsalt")]["est"] / prim["est"], 3),
        "share_gcs_b": round(1 - by[("d1_hosp", "matched_gcs_b_s2gcs15_minus_s1")]["est"] / prim["est"], 3),
        "share_uo": round(1 - by[("d1_hosp", "matched_kidney_s2kcr_minus_s1kcr")]["est"] / prim["est"], 3),
        "event_floor_ok_all_sets": floor_ok,
    }
    res = {"bq_job_id": job, "decisions": decisions, "rows": out}
    (OUT / "discrimination.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(decisions))
    print(f"{len(out)} rows written")
    return 0 if floor_ok else 1


if __name__ == "__main__":
    sys.exit(main())
