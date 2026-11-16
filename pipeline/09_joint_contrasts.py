"""Step 9: jointly aligned contrasts and further descriptions. Reads BUILD_DATASET.variant_stay (step 8) and writes no
table. Only aggregates leave BigQuery; counts below 5 are suppressed.
  1. Jointly aligned AUROC contrasts (paired DeLong, as in step 4), in the primary, full-day-1 and
     suspicion-by-hour-23 sets and each era:
       joint_a = SOFA-2 with PaO2/FiO2-only respiration (variant j) minus SOFA with SOFA-2 brain handling (a);
       joint_b = SOFA-2 with the mimic-code brain convention and PaO2/FiO2-only respiration (b + j) minus SOFA.
     b + j is d1_s2_gcs15 + (d1_s2_j - d1_s2_std): the two variants change different components, so the respiration
     change carries over exactly.
  2. Kidney reference counting only events after an anchor hour: the first hour Sepsis-3 was met (SOFA-2 for Both
     and SOFA-2-only, SOFA for SOFA-only) or, for Neither, the end of the first evaluated hour. Where the number of
     stays with an event at or before the anchor is 1-4, the after-anchor count of that group is withheld.
  3. Brain-rule counts on ICU day 1, primary set: GCS not assessable; sedation-block hours; SOFA-2 carried forward;
     sedated with no earlier assessable score; beyond the 72-hour cap; motor rule; delirium point raising a zero
     brain score; stays whose day-1 SOFA-2 differs between the SOFA-2 and mimic-code brain conventions.
  4. Classification in the complete-case set.
  5. Negative hospital length of stay in the primary set."""
import json
import math
import sys

from bq_client import BUILD, OUT, PROJECT, BQError, query

V4 = BUILD
T = f"{V4}.variant_stay"
MIN_CELL = 5
PRIM = "first_elig_stay AND NOT not_assessable_any"
SUPPRESSED = []

SYSTEMS = ["s1_std", "s1_gcsalt", "s2_std", "s2_gcs15", "s2_j", "s2_gcs15_j"]
SYS_EXPR = {"s2_gcs15_j": "d1_s2_gcs15 + d1_s2_j - d1_s2_std"}
SETS = {
    "primary": PRIM,
    "full_day1": PRIM + " AND d1_last_hr = 23",
    "susp_by_hr23": PRIM + " AND first_susp_hr <= 23",
    "resp_measured": PRIM + " AND resp_measured_both_d1",
    "era_2008_2013": PRIM + " AND era = '2008-2013'",
    "era_2014_2016": PRIM + " AND era = '2014-2016'",
    "era_2017_2022": PRIM + " AND era = '2017-2022'",
}
CONTRASTS = {
    "joint_a": {"s2_j": 1, "s1_gcsalt": -1},
    "joint_b": {"s2_gcs15_j": 1, "s1_std": -1},
    "s2_minus_s1": {"s2_std": 1, "s1_std": -1},
    "matched_gcs_a": {"s2_std": 1, "s1_gcsalt": -1},
    "matched_gcs_b": {"s2_gcs15": 1, "s1_std": -1},
    "j_minus_s1": {"s2_j": 1, "s1_std": -1},
}


def delong_sql():
    sets = ",\n      ".join(f"STRUCT('{k}' AS an, ({v}) AS inc)" for k, v in SETS.items())
    systems = ",\n      ".join(f"STRUCT('{s}' AS sys, {SYS_EXPR.get(s, 'd1_' + s)} AS x)" for s in SYSTEMS)
    pv_cols = ",\n    ".join(f"MAX(IF(sys = '{s}', v, NULL)) AS v_{s}" for s in SYSTEMS)
    defs = []
    for name, coef in [(f"auc_{s}", {s: 1}) for s in SYSTEMS] + list(CONTRASTS.items()):
        coefs = ", ".join(f"{float(coef.get(s, 0))} AS k_{s}" for s in SYSTEMS)
        defs.append(f"STRUCT('{name}' AS contrast, {coefs})")
    combo = " + ".join(f"d.k_{s} * p.v_{s}" for s in SYSTEMS)
    return f"""
WITH long AS (
  SELECT st.an, sy.sys, stay_id, hospital_expire_flag AS y, sy.x
  FROM `{T}`, UNNEST([
      {sets}
    ]) AS st, UNNEST([
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
pv AS (SELECT an, stay_id, ANY_VALUE(y) AS y, {pv_cols} FROM pl GROUP BY an, stay_id),
defs AS (SELECT * FROM UNNEST([
      {",".join(defs)}
    ])),
cv AS (SELECT p.an, d.contrast, p.y, {combo} AS c FROM pv p CROSS JOIN defs d)
SELECT an, contrast, COUNTIF(y = 1) AS n_pos, COUNTIF(y = 0) AS n_neg,
  AVG(IF(y = 1, c, NULL)) AS est,
  SQRT(VAR_SAMP(IF(y = 1, c, NULL)) / COUNTIF(y = 1) + VAR_SAMP(IF(y = 0, c, NULL)) / COUNTIF(y = 0)) AS se
FROM cv GROUP BY an, contrast ORDER BY an, contrast
"""


KIDNEY_SQL = f"""
WITH p AS (
  SELECT stay_id, hadm_id, icu_intime, icu_outtime, onset_s1_std, onset_s2_std,
    CASE WHEN any_s1_std AND any_s2_std THEN 'BOTH' WHEN any_s1_std THEN 'SOFA_ONLY'
         WHEN NOT any_s2_std THEN 'NEITHER'
         WHEN NOT any_s2_kcr THEN 'SOFA2_ONLY_uo_dependent' ELSE 'SOFA2_ONLY_not_uo' END AS grp
  FROM `{T}` WHERE {PRIM}),
fe AS (
  SELECT stay_id, MIN(endtime) AS first_eval_end
  FROM `{V4}.hz_scores` WHERE hz = 'W' AND evaluate AND stay_id IN (SELECT stay_id FROM p) GROUP BY stay_id),
a AS (
  SELECT p.*, CASE WHEN grp IN ('BOTH') OR grp LIKE 'SOFA2_ONLY%' THEN onset_s2_std
                   WHEN grp = 'SOFA_ONLY' THEN onset_s1_std ELSE fe.first_eval_end END AS anchor
  FROM p LEFT JOIN fe USING (stay_id)),
cr AS (
  SELECT l.hadm_id, l.charttime, l.valuenum,
    MIN(l.valuenum) OVER (PARTITION BY l.hadm_id ORDER BY UNIX_SECONDS(TIMESTAMP(l.charttime))
                          RANGE BETWEEN 172800 PRECEDING AND 1 PRECEDING) AS low48,
    MIN(l.valuenum) OVER (PARTITION BY l.hadm_id ORDER BY UNIX_SECONDS(TIMESTAMP(l.charttime))
                          RANGE BETWEEN 604800 PRECEDING AND 1 PRECEDING) AS low7d
  FROM `physionet-data.mimiciv_3_1_hosp.labevents` l
  WHERE l.itemid = 50912 AND l.valuenum > 0 AND l.valuenum < 150
    AND l.hadm_id IN (SELECT hadm_id FROM p)),
kd AS (
  SELECT a.stay_id, MIN(c.charttime) AS first_kdigo
  FROM a JOIN cr c ON c.hadm_id = a.hadm_id AND c.charttime BETWEEN a.icu_intime AND a.icu_outtime
  WHERE c.valuenum >= c.low48 + 0.3 OR c.valuenum >= 1.5 * c.low7d
  GROUP BY a.stay_id),
rr AS (
  SELECT a.stay_id, MIN(r.charttime) AS first_rrt
  FROM a JOIN `@UPSTREAM@.rrt` r
    ON r.stay_id = a.stay_id AND (r.dialysis_active = 1 OR r.dialysis_present = 1)
   AND r.charttime BETWEEN a.icu_intime AND a.icu_outtime
  GROUP BY a.stay_id),
s AS (
  SELECT a.grp, a.anchor, kd.first_kdigo, rr.first_rrt,
    LEAST(COALESCE(kd.first_kdigo, DATETIME '9999-12-31'), COALESCE(rr.first_rrt, DATETIME '9999-12-31')) AS first_any
  FROM a LEFT JOIN kd USING (stay_id) LEFT JOIN rr USING (stay_id))
SELECT grp, COUNT(*) AS n, COUNTIF(anchor IS NULL) AS no_anchor,
  COUNTIF(first_kdigo > anchor) AS kdigo_after, COUNTIF(first_rrt > anchor) AS rrt_after,
  COUNTIF(first_any < DATETIME '9999-12-31' AND first_any > anchor) AS either_after,
  COUNTIF(first_any <= anchor) AS either_at_or_before
FROM s GROUP BY grp ORDER BY grp
"""

BRAIN_SQL = f"""
WITH p AS (SELECT stay_id, d1_last_hr, d1_s2_std, d1_s2_gcs15 FROM `{T}` WHERE {PRIM}),
h AS (
  SELECT s.stay_id, s.hr, p.d1_last_hr, COALESCE(s.unable, 0) AS unable, s.motor, s.gcs_total_min, s.brain_raw,
    COALESCE(s.delirium_treated_any, 0) AS delir,
    IF(COALESCE(s.sedated_hour, 0) = 1
       OR (s.prev_sed_hr IS NOT NULL AND s.next_sed_hr IS NOT NULL AND s.next_sed_hr <= 23), 1, 0) AS sedb
  FROM `{V4}.sofa2_hourly` s JOIN p USING (stay_id)
  WHERE s.hr <= 23),
f AS (
  SELECT h.*,
    LAST_VALUE(IF(sedb = 0 AND brain_raw IS NOT NULL, hr, NULL) IGNORE NULLS)
      OVER (PARTITION BY stay_id ORDER BY hr ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) AS presed_hr,
    LAST_VALUE(IF(sedb = 0, brain_raw, NULL) IGNORE NULLS)
      OVER (PARTITION BY stay_id ORDER BY hr ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) AS presed
  FROM h),
c AS (
  SELECT stay_id, hr,
    unable = 1 AS gcs_not_assessable,
    sedb = 1 AS sedation_block,
    sedb = 1 AND presed_hr IS NOT NULL AND hr - presed_hr <= 72 AS s2_carried_forward,
    sedb = 1 AND presed_hr IS NULL AS s2_sedated_no_prior_score,
    sedb = 1 AND presed_hr IS NOT NULL AND hr - presed_hr > 72 AS s2_sedated_beyond_cap,
    sedb = 0 AND NOT (unable = 0 AND gcs_total_min IS NOT NULL) AND motor IS NOT NULL AS s2_motor_rule,
    sedb = 0 AND brain_raw IS NULL AS s2_unscored_not_sedated,
    delir = 1 AND (CASE WHEN sedb = 1 THEN IF(presed_hr IS NOT NULL AND hr - presed_hr <= 72, COALESCE(presed, 0), 0)
                        ELSE COALESCE(brain_raw, 0) END) = 0 AS s2_delirium_floor
  FROM f WHERE hr >= 0 AND hr <= d1_last_hr),
hours AS (
  SELECT COUNT(*) AS day1_hours, COUNT(DISTINCT stay_id) AS stays,
    COUNTIF(gcs_not_assessable) h_gcs_not_assessable, COUNTIF(sedation_block) h_sedation_block,
    COUNTIF(s2_carried_forward) h_s2_carried_forward, COUNTIF(s2_sedated_no_prior_score) h_s2_sedated_no_prior_score,
    COUNTIF(s2_sedated_beyond_cap) h_s2_sedated_beyond_cap, COUNTIF(s2_motor_rule) h_s2_motor_rule,
    COUNTIF(s2_unscored_not_sedated) h_s2_unscored_not_sedated, COUNTIF(s2_delirium_floor) h_s2_delirium_floor
  FROM c),
stays AS (
  SELECT COUNTIF(a) s_gcs_not_assessable, COUNTIF(b) s_sedation_block, COUNTIF(cf) s_s2_carried_forward,
    COUNTIF(np) s_s2_sedated_no_prior_score, COUNTIF(bc) s_s2_sedated_beyond_cap, COUNTIF(mr) s_s2_motor_rule,
    COUNTIF(un) s_s2_unscored_not_sedated, COUNTIF(dl) s_s2_delirium_floor
  FROM (SELECT stay_id, LOGICAL_OR(gcs_not_assessable) a, LOGICAL_OR(sedation_block) b, LOGICAL_OR(s2_carried_forward) cf,
          LOGICAL_OR(s2_sedated_no_prior_score) np, LOGICAL_OR(s2_sedated_beyond_cap) bc, LOGICAL_OR(s2_motor_rule) mr,
          LOGICAL_OR(s2_unscored_not_sedated) un, LOGICAL_OR(s2_delirium_floor) dl
        FROM c GROUP BY stay_id)),
diff AS (
  SELECT COUNTIF(d1_s2_std > d1_s2_gcs15) AS s2_higher_than_mimic_brain,
         COUNTIF(d1_s2_std < d1_s2_gcs15) AS s2_lower_than_mimic_brain,
         COUNTIF(d1_s2_std = d1_s2_gcs15) AS s2_equal_mimic_brain
  FROM p)
SELECT * FROM hours, stays, diff
"""

CC_SQL = f"""
SELECT COUNT(*) n, COUNTIF(any_s1_std) s1, COUNTIF(any_s2_std) s2,
  COUNTIF(any_s1_std AND any_s2_std) n11, COUNTIF(any_s1_std AND NOT any_s2_std) n10,
  COUNTIF(NOT any_s1_std AND any_s2_std) n01, COUNTIF(NOT any_s1_std AND NOT any_s2_std) n00
FROM `{T}` WHERE {PRIM} AND complete_case
"""

LOS_SQL = f"""
SELECT COUNTIF(a.dischtime < a.admittime) AS dischtime_before_admittime,
       COUNTIF(r.hosp_los_days < 0) AS hour_truncated_los_negative,
       COUNTIF(a.dischtime < a.admittime AND DATETIME_DIFF(a.dischtime, a.admittime, HOUR) = 0) AS negative_under_one_hour_boundary
FROM `{T}` r JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a USING (hadm_id)
WHERE {PRIM.replace('first_elig_stay', 'r.first_elig_stay').replace('not_assessable_any', 'r.not_assessable_any')}
"""

def run(sql):
    fields, rows, _ = query(PROJECT, sql)
    return [dict(zip(fields, r)) for r in rows]


def cnt(v, label):
    v = int(v)
    if 0 < v < MIN_CELL:
        SUPPRESSED.append(label)
        return "<5"
    return v


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [p, max(0.0, c - h), c + h]


def main():
    res = {"label": "exploratory"}
    try:
        rows = []
        for r in run(delong_sql()):
            est, se = float(r["est"]), float(r["se"])
            rows.append({"set": r["an"], "contrast": r["contrast"], "n_pos": int(r["n_pos"]), "n_neg": int(r["n_neg"]),
                         "est": est, "se": se, "lo": est - 1.96 * se, "hi": est + 1.96 * se})
        res["discrimination"] = rows
        res["event_floor_ok"] = all(r["n_pos"] >= 100 and r["n_neg"] >= 100 for r in rows)

        kid = {}
        for r in run(KIDNEY_SQL):
            g, n = r["grp"], int(r["n"])
            before = int(r["either_at_or_before"])
            row = {"n": n, "no_anchor": cnt(r["no_anchor"], f"kidney {g} no_anchor"),
                   "either_at_or_before_anchor": cnt(before, f"kidney {g} at_or_before")}
            withhold = 0 < before < MIN_CELL
            for k in ("kdigo_after", "rrt_after", "either_after"):
                v = int(r[k])
                if withhold:
                    row[k] = "withheld"
                    SUPPRESSED.append(f"kidney {g} {k} (at-or-before count 1-4)")
                elif 0 < v < MIN_CELL or 0 < n - v < MIN_CELL:
                    row[k] = "<5" if v < MIN_CELL else "suppressed"
                    SUPPRESSED.append(f"kidney {g} {k}")
                else:
                    row[k] = [v] + wilson(v, n)
            kid[g] = row
        res["kidney_reference_after_anchor"] = kid
        res["kidney_reference_note"] = (
            "Anchor: first hour Sepsis-3 met (SOFA-2 for Both and SOFA-2-only, SOFA for SOFA-only), or end of the "
            "first evaluated hour for Neither. Event counted when the first KDIGO-positive creatinine (rise >= 0.3 "
            "mg/dL over the lowest in the preceding 48 h or >= 1.5 x the lowest in the preceding 7 days) or the first "
            "RRT hour of the ICU stay is after the anchor; stays with an event at or before the anchor are not "
            "counted as after-anchor events; exploratory")

        b = run(BRAIN_SQL)[0]
        res["brain_rule_counts_day1"] = {k: cnt(v, f"brain {k}") for k, v in b.items()}

        c = run(CC_SQL)[0]
        n = int(c["n"])
        n11, n10, n01, n00 = (int(c[k]) for k in ("n11", "n10", "n01", "n00"))
        po = (n11 + n00) / n
        pe = ((n11 + n10) * (n11 + n01) + (n00 + n10) * (n00 + n01)) / n / n
        cc = {"n": n, "rate_sofa": wilson(int(c["s1"]), n), "rate_sofa2": wilson(int(c["s2"]), n),
              "kappa": (po - pe) / (1 - pe)}
        if 0 < n10 < MIN_CELL or 0 < n01 < MIN_CELL:
            cc.update({"sofa_only": "withheld", "sofa2_only": "withheld", "mcnemar_chi2_cc": "withheld"})
            SUPPRESSED.append("complete case discordant pair and McNemar")
        else:
            d = (n01 - n10) / n
            se = math.sqrt(((n01 + n10) / n - d * d) / n)
            cc.update({"sofa_only": n10, "sofa2_only": n01,
                       "mcnemar_chi2_cc": (abs(n10 - n01) - 1) ** 2 / (n10 + n01) if n10 + n01 else None,
                       "paired_difference_sofa2_minus_sofa": [d, d - 1.96 * se, d + 1.96 * se]})
        res["complete_case_classification"] = cc

        los = run(LOS_SQL)[0]
        res["negative_hosp_los_primary"] = {k: cnt(v, f"los {k}") for k, v in los.items()}
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    res["suppressed"] = SUPPRESSED
    (OUT / "joint_contrasts.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"event_floor_ok": res["event_floor_ok"], "suppressed": SUPPRESSED}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
