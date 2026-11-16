-- Base cohort: first ICU stay of each hospital admission, age at least 18 years, ICU length of stay at least 4 hours.
-- Age is the age at hospital admission (mimic-code age concept).
CREATE OR REPLACE TABLE `@UPSTREAM@.cohort_base` AS
WITH icu AS (
  SELECT i.subject_id, i.hadm_id, i.stay_id, i.first_careunit, i.last_careunit,
         i.intime AS icu_intime, i.outtime AS icu_outtime, i.los AS icu_los_days,
         ROW_NUMBER() OVER (PARTITION BY i.hadm_id ORDER BY i.intime) AS stay_rank
  FROM `physionet-data.mimiciv_3_1_icu.icustays` i
),
dx AS (
  -- primary (seq_num = 1) discharge diagnosis
  SELECT hadm_id, icd_code AS primary_icd, icd_version AS primary_icd_version
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd` WHERE seq_num = 1
)
SELECT icu.subject_id, icu.hadm_id, icu.stay_id, icu.first_careunit, icu.last_careunit,
       icu.icu_intime, icu.icu_outtime, icu.icu_los_days, LN(icu.icu_los_days) AS log_icu_los,
       a.age, p.gender,
       ad.insurance AS insurance_raw,
       CASE WHEN ad.insurance = 'Medicare' THEN 'Medicare' WHEN ad.insurance = 'Medicaid' THEN 'Medicaid'
            WHEN ad.insurance = 'Private' THEN 'Private' ELSE 'Other' END AS payer_group,
       ad.admission_type AS admission_type_raw,
       IF(ad.admission_type = 'URGENT', 'ED', 'Direct') AS admission_source,
       ad.discharge_location, ad.hospital_expire_flag, ad.race, ad.admittime, ad.dischtime,
       dx.primary_icd, dx.primary_icd_version,
       CAST(ad.hospital_expire_flag = 0 AS INT64) AS in_main_analysis   -- hospital survivors (75,472)
FROM icu
JOIN `physionet-data.mimiciv_3_1_hosp.admissions` ad USING (hadm_id)
JOIN `physionet-data.mimiciv_3_1_hosp.patients` p ON p.subject_id = icu.subject_id
JOIN `@UPSTREAM@.age` a ON a.hadm_id = icu.hadm_id
LEFT JOIN dx ON dx.hadm_id = icu.hadm_id
WHERE icu.stay_rank = 1 AND a.age >= 18 AND icu.icu_los_days >= 4/24;
