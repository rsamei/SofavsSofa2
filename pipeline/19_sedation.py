"""Step 19: sensitivity of the brain rule to the sedation definition. Post hoc, descriptive. Only aggregates leave
BigQuery.

Sedation definitions (sql/presubmission/19_sedation_flags.sql): std (primary, recomputed), nodex (dexmedetomidine
excluded), rass (sedative infusion and last RASS in the preceding 4 hours <= -2; hours without a RASS not
sedated) and light (primary, except hours whose last RASS is -1 or higher; hours without a RASS stay sedated). For each, the day-1 SOFA-2 and
SOFA with SOFA-2's brain handling (variant a) are rescored with sql/scores/03_sofa2_hourly.sql and
sql/scores/05_hz_scores.sql, unchanged except for the sedation flag they read, the tables they write and, in 05, the
day-1 horizon only. The std rescoring must reproduce the day-1 scores of the cohort exactly (checked).

Reports: day-1 sedated stays and hours under each definition; paired DeLong AUROCs of SOFA-2 under each definition
against SOFA, against SOFA-2 (primary), and against SOFA with the same brain handling (variant a).
"""
import json
import sys

from bq_client import BUILD, OUT, PROJECT, SQL_DIR, BQError, dry_run, execute, query
from delong import run_delong

VARIANTS = ["std", "nodex", "rass", "light"]
PRIM = "first_elig_stay AND NOT not_assessable_any"
MIN_CELL = 5


def sub(sql, pairs, name):
    for old, new, count in pairs:
        if sql.count(old) != count:
            raise BQError(f"{name}: anchor {old!r} found {sql.count(old)} times, expected {count}")
        sql = sql.replace(old, new)
    return sql


def build_sql(v):
    sw = f"""CREATE OR REPLACE VIEW `@BUILD@.s19_sw_{v}` AS
SELECT w.* REPLACE (f.sed_{v} AS sedated_hour)
FROM `@BUILD@.sofa2_withintime` w JOIN `@BUILD@.s19_sed` f USING (stay_id, hr)"""
    h = sub((SQL_DIR / "scores" / "03_sofa2_hourly.sql").read_text(encoding="utf-8"), [
        ("CREATE OR REPLACE TABLE `@BUILD@.sofa2_hourly` AS", f"CREATE OR REPLACE TABLE `@BUILD@.s19_sofa2_hourly_{v}` AS", 1),
        ("`@BUILD@.sofa2_withintime`", f"`@BUILD@.s19_sw_{v}`", 2)], "03_sofa2_hourly.sql")
    z = sub((SQL_DIR / "scores" / "05_hz_scores.sql").read_text(encoding="utf-8"), [
        ("CREATE OR REPLACE TABLE `@BUILD@.hz_scores` AS", f"CREATE OR REPLACE TABLE `@BUILD@.s19_hz_{v}` AS", 1),
        ("  FROM `@BUILD@.sofa2_hourly`\n", f"  FROM `@BUILD@.s19_sofa2_hourly_{v}`\n", 1),
        ("r AS (SELECT * FROM rows_d1 UNION ALL SELECT * FROM rows_w),", "r AS (SELECT * FROM rows_d1),", 1)],
        "05_hz_scores.sql")
    return [(f"s19_sw_{v}", sw), (f"s19_sofa2_hourly_{v}", h), (f"s19_hz_{v}", z)]


D1_SQL = "CREATE OR REPLACE TABLE `@BUILD@.s19_d1` AS\nWITH " + ",\n".join(
    f"""d_{v} AS (
  SELECT stay_id, s2_std AS d1_s2sed_{v}, s1_gcsalt AS d1_s1ased_{v}
  FROM `@BUILD@.s19_hz_{v}` WHERE hz = 'D1'
  QUALIFY ROW_NUMBER() OVER (PARTITION BY stay_id ORDER BY hr DESC) = 1)""" for v in VARIANTS) + """,
sh AS (
  SELECT c.stay_id, """ + ", ".join(f"MAX(f.sed_{v}) AS any_sed_{v}, SUM(f.sed_{v}) AS hrs_sed_{v}" for v in VARIANTS) + """
  FROM `@BUILD@.variant_stay` c JOIN `@BUILD@.s19_sed` f ON f.stay_id = c.stay_id AND f.hr BETWEEN 0 AND c.d1_last_hr
  GROUP BY 1)
SELECT c.*, """ + ", ".join(f"d_{v}.* EXCEPT (stay_id)" for v in VARIANTS) + """, sh.* EXCEPT (stay_id)
FROM `@BUILD@.variant_stay` c """ + " ".join(f"LEFT JOIN d_{v} USING (stay_id)" for v in VARIANTS) + """
LEFT JOIN sh USING (stay_id)"""

CHECK_SQL = f"""SELECT
  (SELECT COUNTIF(sed_std IS DISTINCT FROM sedated_hour) FROM `@BUILD@.s19_sed`) AS sed_std_mismatch_hours,
  (SELECT COUNTIF(f.sed_std = 1 AND f.last_rass IS NULL) FROM `@BUILD@.s19_d1` c
   JOIN `@BUILD@.s19_sed` f ON f.stay_id = c.stay_id AND f.hr BETWEEN 0 AND c.d1_last_hr
   WHERE c.first_elig_stay AND NOT c.not_assessable_any) AS hours_sed_std_without_rass,
  COUNTIF(d1_s2sed_std IS DISTINCT FROM d1_s2_std) AS s2_std_rescore_mismatch,
  COUNTIF(d1_s1ased_std IS DISTINCT FROM d1_s1_gcsalt) AS s1a_std_rescore_mismatch,
  COUNT(*) AS n,
  {", ".join(f"COUNTIF(any_sed_{v} = 1) AS stays_sed_{v}, SUM(hrs_sed_{v}) AS hours_sed_{v}" for v in VARIANTS)},
  {", ".join(f"COUNTIF(d1_s2sed_{v} != d1_s2_std) AS s2_changed_{v}, COUNTIF(d1_s2sed_{v} > d1_s2_std) AS s2_higher_{v}"
             for v in VARIANTS[1:])}
FROM `@BUILD@.s19_d1` WHERE {PRIM}"""


def main():
    try:
        steps = [("s19_sed", (SQL_DIR / "presubmission" / "19_sedation_flags.sql").read_text(encoding="utf-8"))]
        for v in VARIANTS:
            steps += build_sql(v)
        steps.append(("s19_d1", D1_SQL))
        for name, sql in steps:
            dry_run(PROJECT, sql)
            execute(PROJECT, sql)
            print(f"built {name}", flush=True)
        f, r, _ = query(PROJECT, CHECK_SQL)
        checks = {k: int(float(v)) for k, v in zip(f, r[0])}
        systems = {"s1": "d1_s1_std", "s2": "d1_s2_std", "s1a": "d1_s1_gcsalt"}
        contrasts = {}
        for v in VARIANTS[1:]:
            systems[f"s2_{v}"] = f"d1_s2sed_{v}"
            systems[f"s1a_{v}"] = f"d1_s1ased_{v}"
            contrasts[f"s2{v}_minus_s1"] = {f"s2_{v}": 1, "s1": -1}
            contrasts[f"s2{v}_minus_s2"] = {f"s2_{v}": 1, "s2": -1}
            contrasts[f"s2{v}_minus_s1a{v}"] = {f"s2_{v}": 1, f"s1a_{v}": -1}
        contrasts["s2_minus_s1"] = {"s2": 1, "s1": -1}
        contrasts["s2_minus_s1a"] = {"s2": 1, "s1a": -1}
        disc = run_delong(f"`{BUILD}.s19_d1`", {"primary": PRIM}, systems, contrasts)
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    passed = all(checks[k] == 0 for k in ("sed_std_mismatch_hours", "s2_std_rescore_mismatch", "s1a_std_rescore_mismatch"))
    suppressed = [k for k, v in checks.items() if 0 < v < MIN_CELL]
    for k in suppressed:
        checks[k] = "<5"
    res = {"pass": passed, "checks": checks, "discrimination": disc, "suppressed": suppressed,
           "definitions": {"std": "propofol, midazolam, lorazepam or dexmedetomidine infusion overlapping the hour",
                           "nodex": "as std without dexmedetomidine",
                           "rass": "as std and last RASS (item 228096) in the 4 hours up to the hour end <= -2",
                           "light": "as std except hours whose last RASS in the 4 hours up to the hour end is >= -1"}}
    (OUT / "sedation.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"pass": passed, **checks}))
    for x in disc:
        print(f"{x['contrast']:22s} {x['est']:.4f} ({x['lo']:.4f} to {x['hi']:.4f})")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
