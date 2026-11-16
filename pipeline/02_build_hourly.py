"""Step 2: build the hourly score tables in BUILD_DATASET (sql/scores/01-04), then check them and describe eMAR
coverage.

Checks (only counts leave BigQuery):
  - the kidney column equals kidney_b after COALESCE(., 0);
  - the coverage-restricted kidney scores (variant g) never exceed the single-chart ones (variant s).
Descriptives (no outcomes), by anchor_year_group: admissions and ICU stays with any eMAR record, stays with recorded
delirium treatment or dexmedetomidine, and ICU stays with an antipsychotic or dexmedetomidine in icu.inputevents and
no eMAR record of that class. Counts below 5 are suppressed.
"""
import json
import sys

from bq_client import BUILD_DATASET, OUT, PROJECT, SQL_DIR, BQError, dry_run, ensure_dataset, execute, query

STEPS = ["01_sofa2_withintime.sql", "02_cv_hourly.sql", "03_sofa2_hourly.sql", "04_sofa1_hourly.sql"]
MIN_CELL = 5

CHECK_SQL = """
SELECT COUNT(*) AS sofa2_withintime_rows,
  COUNTIF(COALESCE(kidney, 0) != COALESCE(kidney_b, 0)) AS kidney_ne_kidney_b,
  COUNTIF(kidney_gb > COALESCE(kidney_sb, 0)) AS kidney_gb_exceeds_sb,
  COUNTIF(kidney_go > COALESCE(kidney_so, 0)) AS kidney_go_exceeds_so,
  COUNTIF(hr BETWEEN 0 AND 23 AND weight_backapplied = 1) AS d1_hours_backapplied,
  COUNTIF(hr BETWEEN 0 AND 23 AND weight_backapplied = 1 AND wt_adm_hr > 23) AS d1_hours_backapplied_after_cut,
  COUNT(DISTINCT IF(hr BETWEEN 0 AND 23 AND weight_backapplied = 1 AND wt_adm_hr > 23, stay_id, NULL)) AS d1_stays_backapplied_after_cut
FROM `@BUILD@.sofa2_withintime`
"""

EMAR_SQL = """
WITH pop AS (
  SELECT c.subject_id, c.hadm_id, c.stay_id, c.icu_intime, c.icu_outtime, p.anchor_year_group
  FROM `@UPSTREAM@.cohort_base` c
  JOIN `physionet-data.mimiciv_3_1_hosp.patients` p USING (subject_id)
  WHERE c.age >= 18
),
em_adm AS (SELECT DISTINCT hadm_id FROM `physionet-data.mimiciv_3_1_hosp.emar` WHERE hadm_id IS NOT NULL),
em_stay AS (
  SELECT DISTINCT p.stay_id
  FROM pop p JOIN `physionet-data.mimiciv_3_1_hosp.emar` e
    ON e.hadm_id = p.hadm_id AND e.charttime BETWEEN p.icu_intime AND p.icu_outtime
),
tx AS (
  SELECT stay_id, MAX(delirium_treated_any) AS del_emar, MAX(dex_treated_any) AS dex_emar
  FROM `@BUILD@.sofa2_withintime` GROUP BY 1
),
items AS (
  SELECT itemid,
    REGEXP_CONTAINS(LOWER(label), r'haloperidol|quetiapine|olanzapine|risperidone|ziprasidone|chlorpromazine') AS is_ap,
    REGEXP_CONTAINS(LOWER(label), r'dexmedetomidine|precedex') AS is_dex
  FROM `physionet-data.mimiciv_3_1_icu.d_items`
  WHERE REGEXP_CONTAINS(LOWER(label), r'haloperidol|quetiapine|olanzapine|risperidone|ziprasidone|chlorpromazine|dexmedetomidine|precedex')
),
ie AS (
  SELECT i.stay_id, LOGICAL_OR(it.is_ap) AS ap_ie, LOGICAL_OR(it.is_dex) AS dex_ie
  FROM `physionet-data.mimiciv_3_1_icu.inputevents` i JOIN items it USING (itemid)
  WHERE i.amount > 0
  GROUP BY 1
)
SELECT p.anchor_year_group,
  COUNT(DISTINCT p.hadm_id) AS admissions,
  COUNT(DISTINCT IF(ea.hadm_id IS NOT NULL, p.hadm_id, NULL)) AS admissions_with_emar,
  COUNT(*) AS icu_stays,
  COUNTIF(es.stay_id IS NOT NULL) AS stays_with_emar_in_icu,
  COUNTIF(COALESCE(tx.del_emar, 0) = 1) AS stays_delirium_tx_emar,
  COUNTIF(COALESCE(tx.dex_emar, 0) = 1) AS stays_dex_emar,
  COUNTIF(COALESCE(ie.ap_ie, FALSE)) AS stays_antipsychotic_inputevents,
  COUNTIF(COALESCE(ie.ap_ie, FALSE) AND COALESCE(tx.del_emar, 0) = 0) AS stays_antipsychotic_ie_no_emar,
  COUNTIF(COALESCE(ie.dex_ie, FALSE)) AS stays_dex_inputevents,
  COUNTIF(COALESCE(ie.dex_ie, FALSE) AND COALESCE(tx.dex_emar, 0) = 0) AS stays_dex_ie_no_emar
FROM pop p
LEFT JOIN em_adm ea ON ea.hadm_id = p.hadm_id
LEFT JOIN em_stay es ON es.stay_id = p.stay_id
LEFT JOIN tx ON tx.stay_id = p.stay_id
LEFT JOIN ie ON ie.stay_id = p.stay_id
GROUP BY 1 ORDER BY 1
"""


def main():
    try:
        ensure_dataset(PROJECT, BUILD_DATASET)
        for name in STEPS:
            sql = (SQL_DIR / "scores" / name).read_text(encoding="utf-8")
            dry_run(PROJECT, sql)
            execute(PROJECT, sql)
            print(f"built {name}", flush=True)
        fields, rows, _ = query(PROJECT, CHECK_SQL)
        efields, erows, _ = query(PROJECT, EMAR_SQL)
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    res = {k: (None if v is None else int(v)) for k, v in zip(fields, rows[0])}
    suppressed = []
    by_year = []
    for r in erows:
        d = dict(zip(efields, r))
        row = {"anchor_year_group": d.pop("anchor_year_group")}
        for k, v in d.items():
            v = int(v)
            if 0 < v < MIN_CELL:
                suppressed.append(f"{row['anchor_year_group']}:{k}")
                v = "<5"
            row[k] = v
        for num, den in [("admissions_with_emar", "admissions"), ("stays_with_emar_in_icu", "icu_stays"),
                         ("stays_delirium_tx_emar", "icu_stays"), ("stays_dex_emar", "icu_stays"),
                         ("stays_antipsychotic_inputevents", "icu_stays"), ("stays_dex_inputevents", "icu_stays")]:
            if isinstance(row[num], int) and row[den]:
                row[f"share_{num}"] = round(row[num] / row[den], 4)
        by_year.append(row)
    passed = res["kidney_ne_kidney_b"] == 0 and res["kidney_gb_exceeds_sb"] == 0 and res["kidney_go_exceeds_so"] == 0
    out = {"checks": res, "pass": passed, "emar_inputevents_by_year_group": by_year, "suppressed": suppressed}
    (OUT / "hourly_check.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"pass": passed, **res}))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
