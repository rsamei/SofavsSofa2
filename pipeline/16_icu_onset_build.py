"""Step 16 (Supplementary Methods S1): build the ICU-onset analysis tables in BUILD_DATASET.

  1. s1_elig (sql/supplement_s1/01_elig.sql): landmark at the end of ICU hour 47, no suspicion of infection up to the
     landmark, first suspicion in ICU hours 48-239, first eligible stay per patient.
  2. s1_scores: the hourly SOFA and SOFA-2 of the main analysis at two horizons per eligible stay, built from
     sql/scores/05_hz_scores.sql with only the scored rows replaced, so every implementation rule is the main one:
       'B' baseline horizon, ICU hours 0-47 with look-ahead truncated at hour 47, evaluated at hours 23-47;
       'E' evaluation window, suspicion - 48 h to + 24 h, evaluated from hour 48 to ICU discharge, with look-ahead
           truncated at the window end (as the 'W' horizon of the main analysis).
  3. s1_cohort (sql/supplement_s1/02_cohort.sql): one row per analysed stay.
Checks (only counts leave BigQuery): 25 baseline hours per stay and score, component sums equal totals, hour-47
baseline present. Counts below 5 are suppressed. Reads the main tables (sofa2_hourly, sofa1_hourly, cv_hourly,
sofa2_withintime) and changes none of them.
"""
import json
import sys

from bq_client import OUT, PROJECT, SQL_DIR, BQError, dry_run, execute, query

MIN_CELL = 5
S1_DIR = SQL_DIR / "supplement_s1"

ROWS = """s1_w AS (
  SELECT stay_id, intime_hr, outtime_hr, DATETIME_SUB(susp_time, INTERVAL 48 HOUR) AS w_start,
         DATETIME_ADD(susp_time, INTERVAL 24 HOUR) AS w_end
  FROM `@BUILD@.s1_elig` WHERE first_elig_stay
),
s1_wcut AS (
  SELECT w.stay_id, MAX(s2.hr) AS t_cut
  FROM s1_w w JOIN s2 ON s2.stay_id = w.stay_id AND s2.endtime <= w.w_end
  GROUP BY 1
),
rows_s1b AS (  -- baseline horizon: ICU hours 0-47, look-ahead truncated at hour 47
  SELECT s2.stay_id, 'B' AS hz, 0 AS hz_key, 47 AS t_cut, s2.hr, CAST(NULL AS BOOL) AS old_rule,
         s2.hr BETWEEN 23 AND 47 AS evaluate, FALSE AS is_base
  FROM s2 JOIN s1_w USING (stay_id)
  WHERE s2.hr BETWEEN 0 AND 47
),
rows_s1e AS (  -- evaluation window, from hour 48 to ICU discharge, look-ahead truncated at the window end
  SELECT w.stay_id, 'E' AS hz, 0 AS hz_key, c.t_cut, s2.hr, CAST(NULL AS BOOL) AS old_rule,
         (s2.hr >= 48 AND s2.endtime BETWEEN w.w_start AND w.w_end
          AND s2.endtime BETWEEN w.intime_hr AND w.outtime_hr) AS evaluate,
         FALSE AS is_base
  FROM s1_w w JOIN s1_wcut c USING (stay_id)
  JOIN s2 ON s2.stay_id = w.stay_id AND s2.endtime <= w.w_end
         AND s2.endtime >= DATETIME_SUB(w.w_start, INTERVAL 100 HOUR)
),
r AS (SELECT * FROM rows_s1b UNION ALL SELECT * FROM rows_s1e),"""


def scores_sql():
    """sql/scores/05_hz_scores.sql with the target table and the scored rows replaced; fails if either anchor moves."""
    sql = (SQL_DIR / "scores" / "05_hz_scores.sql").read_text(encoding="utf-8")
    subs = [("CREATE OR REPLACE TABLE `@BUILD@.hz_scores` AS", "CREATE OR REPLACE TABLE `@BUILD@.s1_scores` AS"),
            ("r AS (SELECT * FROM rows_d1 UNION ALL SELECT * FROM rows_w),", ROWS)]
    for old, new in subs:
        if sql.count(old) != 1:
            raise BQError(f"05_hz_scores.sql anchor not found once: {old!r}")
        sql = sql.replace(old, new)
    return sql


CHECK_SQL = """
WITH e AS (SELECT * FROM `@BUILD@.s1_elig`), c AS (SELECT * FROM `@BUILD@.s1_cohort`),
z AS (SELECT * FROM `@BUILD@.s1_scores` WHERE evaluate)
SELECT
  (SELECT COUNT(*) FROM `@UPSTREAM@.cohort_base` WHERE age >= 18) AS base_stays,
  (SELECT COUNT(*) FROM e) AS alive_in_icu_at_landmark,
  (SELECT COUNTIF(early_susp) FROM e) AS excluded_suspicion_by_hour_47,
  (SELECT COUNTIF(NOT early_susp AND susp_time IS NULL) FROM e) AS excluded_no_later_suspicion,
  (SELECT COUNTIF(NOT early_susp AND susp_time IS NOT NULL AND NOT eligible) FROM e) AS excluded_suspicion_after_hour_239_or_discharge,
  (SELECT COUNTIF(eligible) FROM e) AS eligible_stays,
  (SELECT COUNTIF(eligible AND NOT first_elig_stay) FROM e) AS later_eligible_stays_same_patient,
  (SELECT COUNT(*) FROM c) AS first_eligible_stays,
  (SELECT COUNTIF(not_assessable) FROM c) AS not_assessable,
  (SELECT COUNTIF(NOT not_assessable) FROM c) AS analysed,
  (SELECT COUNTIF(s1_n_base_hours != 25 OR s2_n_base_hours != 25) FROM c) AS base_hours_not_25,
  (SELECT COUNTIF(s1_base47 IS NULL OR s2_base47 IS NULL) FROM c) AS missing_hour47_baseline,
  (SELECT COUNTIF(s1_std != s1_resp24 + s1_coag24 + s1_liver24 + s1_cv24 + s1_cns24 + s1_renal24) FROM z) AS s1_component_sum_mismatch,
  (SELECT COUNTIF(s2_std != s2_resp24 + s2_hemo24 + s2_liver24 + s2_cv24 + s2_brain24 + s2_kidney24) FROM z) AS s2_component_sum_mismatch,
  (SELECT COUNTIF(chk_s2_old IS DISTINCT FROM chk_s2_new) FROM `@BUILD@.s1_scores`) AS recompute_mismatch_rows
"""


def main():
    try:
        for label, sql in [("s1_elig", (S1_DIR / "01_elig.sql").read_text(encoding="utf-8")),
                           ("s1_scores", scores_sql()),
                           ("s1_cohort", (S1_DIR / "02_cohort.sql").read_text(encoding="utf-8"))]:
            dry_run(PROJECT, sql)
            execute(PROJECT, sql)
            print(f"built {label}", flush=True)
        fields, rows, _ = query(PROJECT, CHECK_SQL)
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    res = {k: (None if v is None else int(float(v))) for k, v in zip(fields, rows[0])}
    passed = all(res[k] == 0 for k in ("base_hours_not_25", "missing_hour47_baseline", "s1_component_sum_mismatch",
                                       "s2_component_sum_mismatch", "recompute_mismatch_rows"))
    suppressed = [k for k, v in res.items() if v is not None and 0 < v < MIN_CELL]
    for k in suppressed:
        res[k] = "<5"
    out = {"checks": res, "pass": passed, "suppressed": suppressed}
    (OUT / "icu_onset_build.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"pass": passed, **res}))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
