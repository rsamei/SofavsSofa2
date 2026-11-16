"""Step 1: build the upstream tables in UPSTREAM_DATASET from MIMIC-IV v3.1.

  - 26 mimic-code concepts (sql/mimic_code, commit 4582f39), in dependency order, with physionet-data.mimiciv_hosp /
    mimiciv_icu mapped to the v3.1 datasets and mimiciv_derived mapped to UPSTREAM_DATASET;
  - cohort_base (sql/upstream/cohort_base.sql);
  - sofa_withintime, the mimic-code hourly SOFA (sql/upstream/sofa_withintime.sql).

Pass --skip-existing to keep tables that already exist (e.g. to resume after an interruption).
"""
import sys

from bq_client import PROJECT, SQL_DIR, UPSTREAM, UPSTREAM_DATASET, BQError, ensure_dataset, execute, query

CONCEPTS = [  # dependency order
    "demographics/age", "comorbidity/charlson", "demographics/icustay_times", "demographics/icustay_hourly", "demographics/weight_durations",
    "medication/antibiotic", "measurement/oxygen_delivery", "measurement/ventilator_setting", "treatment/ventilation",
    "measurement/bg", "measurement/chemistry", "measurement/complete_blood_count", "measurement/enzyme",
    "measurement/gcs", "measurement/urine_output", "measurement/urine_output_rate", "measurement/vitalsign",
    "medication/dobutamine", "medication/dopamine", "medication/epinephrine", "medication/milrinone",
    "medication/norepinephrine", "medication/phenylephrine", "medication/vasopressin", "treatment/rrt",
    "sepsis/suspicion_of_infection",
]


def map_mimic_code(sql):
    return (sql.replace("physionet-data.mimiciv_hosp.", "physionet-data.mimiciv_3_1_hosp.")
               .replace("physionet-data.mimiciv_icu.", "physionet-data.mimiciv_3_1_icu.")
               .replace("physionet-data.mimiciv_derived.", UPSTREAM + "."))


def create(table, sql):
    body = sql.strip()
    while body.endswith(";"):
        body = body[:-1].rstrip()
    return f"CREATE OR REPLACE TABLE `{UPSTREAM}.{table}` AS\n{body}\n"


def existing():
    _, rows, _ = query(PROJECT, f"SELECT table_name FROM `{UPSTREAM}.INFORMATION_SCHEMA.TABLES`")
    return {r[0] for r in rows}


def main():
    skip = "--skip-existing" in sys.argv
    try:
        ensure_dataset(PROJECT, UPSTREAM_DATASET)
        have = existing() if skip else set()
        for c in CONCEPTS:
            name = c.split("/")[1]
            if name in have:
                continue
            execute(PROJECT, create(name, map_mimic_code((SQL_DIR / "mimic_code" / f"{c}.sql").read_text(encoding="utf-8"))))
            print(f"built {name}", flush=True)
        for name in ("cohort_base", "sofa_withintime"):
            if name in have:
                continue
            execute(PROJECT, (SQL_DIR / "upstream" / f"{name}.sql").read_text(encoding="utf-8"))
            print(f"built {name}", flush=True)
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
