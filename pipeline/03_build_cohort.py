"""Step 3: build the horizon scores and the per-stay cohort in BUILD_DATASET (sql/scores/05-06), then report the
analysis-set sizes.

Checks (only counts leave BigQuery): every first eligible stay has both day-1 scores, and SOFA-2 rebuilt per horizon
from its components equals the total from the hourly table in every row. Counts below 5 are suppressed.
"""
import json
import sys

from bq_client import OUT, PROJECT, SQL_DIR, BQError, dry_run, execute, query

STEPS = ["05_hz_scores.sql", "06_cohort.sql"]
MIN_CELL = 5

CHECK_SQL = """
WITH n AS (SELECT * FROM `@BUILD@.cohort`),
f AS (SELECT * FROM n WHERE first_elig_stay)
SELECT
  (SELECT COUNT(*) FROM `@UPSTREAM@.cohort_base`) AS base_stays,
  (SELECT COUNT(*) FROM n) AS elig_stays,
  (SELECT COUNT(*) FROM f) AS first_elig_stays,
  (SELECT COUNTIF(chk_s2_old IS DISTINCT FROM chk_s2_new) FROM `@BUILD@.hz_scores`) AS recompute_mismatch_rows,
  (SELECT COUNTIF(d1_s1_std IS NULL OR d1_s2_std IS NULL) FROM f) AS missing_day1,
  (SELECT COUNTIF(not_assessable_any) FROM f) AS not_assessable_any,
  (SELECT COUNTIF(NOT not_assessable_any) FROM f) AS primary_n,
  (SELECT COUNTIF(landmark AND NOT not_assessable_any) FROM f) AS landmark_n,
  (SELECT COUNTIF(complete_case AND NOT not_assessable_any) FROM f) AS complete_case_n,
  (SELECT COUNTIF(ind_baseline_available AND NOT not_assessable_any) FROM f) AS ind_baseline_available_n,
  (SELECT COUNTIF(suspicion_before_icu AND NOT not_assessable_any) FROM f) AS suspicion_before_icu_n,
  (SELECT COUNTIF(era = '2008-2013' AND NOT not_assessable_any) FROM f) AS era_2008_2013_n,
  (SELECT COUNTIF(era = '2014-2016' AND NOT not_assessable_any) FROM f) AS era_2014_2016_n,
  (SELECT COUNTIF(era = '2017-2022' AND NOT not_assessable_any) FROM f) AS era_2017_2022_n,
  (SELECT COUNTIF(emar_adm AND NOT not_assessable_any) FROM f) AS emar_adm_n
"""

GROUP_SQL = """
SELECT
  CASE WHEN any_s1_std AND any_s2_std THEN 'BOTH' WHEN any_s2_std THEN 'SOFA2_ONLY'
       WHEN any_s1_std THEN 'SOFA_ONLY' ELSE 'NEITHER' END AS grp,
  COUNT(*) AS n,
  COUNTIF(one_volume_s2) AS s2_cases_resting_on_one_volume,
  COUNTIF(any_s2_std AND NOT COALESCE(any_s2_g, FALSE)) AS s2_cases_lost_under_g,
  COUNTIF(any_s2_std AND NOT COALESCE(any_s2_kcr, FALSE)) AS s2_cases_lost_without_uo
FROM `@BUILD@.cohort`
WHERE first_elig_stay AND NOT not_assessable_any
GROUP BY 1 ORDER BY 1
"""


def suppress(d, label, out):
    for k, v in d.items():
        if isinstance(v, int) and 0 < v < MIN_CELL:
            d[k] = "<5"
            out.append(f"{label}:{k}")


def main():
    try:
        for name in STEPS:
            sql = (SQL_DIR / "scores" / name).read_text(encoding="utf-8")
            dry_run(PROJECT, sql)
            execute(PROJECT, sql)
            print(f"built {name}", flush=True)
        fields, rows, _ = query(PROJECT, CHECK_SQL)
        gfields, grows, _ = query(PROJECT, GROUP_SQL)
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    res = {k: (None if v is None else int(float(v))) for k, v in zip(fields, rows[0])}
    groups = {}
    for r in grows:
        d = dict(zip(gfields, r))
        g = d.pop("grp")
        groups[g] = {k: int(v) for k, v in d.items()}
    passed = res["recompute_mismatch_rows"] == 0 and res["missing_day1"] == 0
    suppressed = []
    suppress(res, "checks", suppressed)
    for g, d in groups.items():
        suppress(d, g, suppressed)
    out = {"checks": res, "pass": passed, "by_concordance_group": groups, "suppressed": suppressed}
    (OUT / "cohort_check.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"pass": passed, **res}))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
