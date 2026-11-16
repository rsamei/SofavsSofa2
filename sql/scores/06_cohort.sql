-- One row per eligible ICU stay: eligibility and analysis-set flags, outcomes, day-1 scores and any-window Sepsis-3
-- flags for SOFA, SOFA-2 and the variants.
CREATE OR REPLACE TABLE `@BUILD@.cohort` AS
WITH pop AS (
  SELECT c.subject_id, c.hadm_id, c.stay_id, t.intime_hr, t.outtime_hr, c.icu_intime, c.icu_outtime, pt.anchor_year_group,
         c.admittime, c.dischtime, c.hospital_expire_flag, c.icu_los_days, c.age, c.gender, c.first_careunit
  FROM `@UPSTREAM@.cohort_base` c JOIN `@UPSTREAM@.icustay_times` t USING (stay_id)
  JOIN `physionet-data.mimiciv_3_1_hosp.patients` pt ON pt.subject_id = c.subject_id
  WHERE c.age >= 18
),
susp AS (
  SELECT hadm_id, suspected_infection_time, ARRAY_AGG(DISTINCT stay_id IGNORE NULLS) AS s_stays
  FROM `@UPSTREAM@.suspicion_of_infection`
  WHERE suspected_infection = 1 AND hadm_id IS NOT NULL
  GROUP BY 1, 2
),
win AS (
  SELECT p.stay_id, s.suspected_infection_time, p.stay_id IN UNNEST(s.s_stays) AS old_rule
  FROM pop p JOIN susp s USING (hadm_id)
  WHERE DATETIME_ADD(s.suspected_infection_time, INTERVAL 24 HOUR) >= p.intime_hr
    AND DATETIME_SUB(s.suspected_infection_time, INTERVAL 48 HOUR) <= p.outtime_hr
),
elig AS (
  SELECT p.*, w.first_susp_time, w.n_windows, w.elig_old,
         DATETIME_DIFF(w.first_susp_time, p.intime_hr, HOUR) AS first_susp_hr,
         ROW_NUMBER() OVER (PARTITION BY p.subject_id ORDER BY p.icu_intime, p.stay_id) = 1 AS first_elig_stay,
         w.elig_old AND ROW_NUMBER() OVER (PARTITION BY p.subject_id, w.elig_old ORDER BY p.icu_intime, p.stay_id) = 1
           AS first_elig_stay_old
  FROM pop p
  JOIN (SELECT stay_id, MIN(suspected_infection_time) AS first_susp_time, COUNT(*) AS n_windows,
               LOGICAL_OR(old_rule) AS elig_old
        FROM win GROUP BY 1) w
    USING (stay_id)
),
hzs AS (SELECT * FROM `@BUILD@.hz_scores`),
d1 AS (
  SELECT * EXCEPT (rn) FROM (
    SELECT stay_id, hr AS d1_last_hr,
      s1_std AS d1_s1_std, s1_gcsalt AS d1_s1_gcsalt, s1_kcr AS d1_s1_kcr,
      s2_std AS d1_s2_std, s2_e AS d1_s2_e, s2_gcs15 AS d1_s2_gcs15, s2_kcr AS d1_s2_kcr, s2_kuo22 AS d1_s2_kuo22,
      s2_g AS d1_s2_g, s2_s AS d1_s2_s, s2_w AS d1_s2_w, s1_pipe AS d1_s1_pipe, s1_orig AS d1_s1_orig, s1_orig_n AS d1_s1_orig_n, s2_gcs15_h AS d1_s2_gcs15_h, s2_h AS d1_s2_h, s2_i AS d1_s2_i, s2_kidney24 AS d1_s2_kidney, s2_kidney_cr24 AS d1_s2_kidney_cr,
      s1_cns24 AS d1_s1_cns, s1_cns_alt24 AS d1_s1_cns_alt, s1_cns_orig24 AS d1_s1_cns_orig, s2_brain24 AS d1_s2_brain,
      ROW_NUMBER() OVER (PARTITION BY stay_id ORDER BY hr DESC) AS rn
    FROM hzs WHERE hz = 'D1')
  WHERE rn = 1
),
ev AS (SELECT * FROM hzs WHERE hz = 'W' AND evaluate),
anyw AS (
  SELECT stay_id,
    COUNT(*) AS n_eval_hours_any,
    LOGICAL_OR(s1_std >= 2) AS any_s1_std, LOGICAL_OR(s1_gcsalt >= 2) AS any_s1_gcsalt, LOGICAL_OR(s1_kcr >= 2) AS any_s1_kcr,
    LOGICAL_OR(s1_pipe >= 2) AS any_s1_pipe, LOGICAL_OR(s1_orig >= 2) AS any_s1_orig, LOGICAL_OR(s1_orig_n >= 2) AS any_s1_orig_n,
    LOGICAL_OR(s2_std >= 2) AS any_s2_std, LOGICAL_OR(s2_e >= 2) AS any_s2_e, LOGICAL_OR(s2_gcs15 >= 2) AS any_s2_gcs15,
    LOGICAL_OR(s2_kcr >= 2) AS any_s2_kcr, LOGICAL_OR(s2_kuo22 >= 2) AS any_s2_kuo22,
    LOGICAL_OR(s2_g >= 2) AS any_s2_g, LOGICAL_OR(s2_s >= 2) AS any_s2_s, LOGICAL_OR(s2_w >= 2) AS any_s2_w, LOGICAL_OR(s2_h >= 2) AS any_s2_h, LOGICAL_OR(s2_i >= 2) AS any_s2_i,
    LOGICAL_OR(s2_1v >= 2) AS any_s2_1v,
    MIN(IF(s1_std >= 2, endtime, NULL)) AS onset_s1_std, MIN(IF(s2_std >= 2, endtime, NULL)) AS onset_s2_std,
    -- stay_id eligibility rule (sensitivity row)
    COUNTIF(old_rule) AS n_eval_hours_old,
    LOGICAL_OR(old_rule AND s1_std >= 2) AS old_any_s1_std, LOGICAL_OR(old_rule AND s2_std >= 2) AS old_any_s2_std,
    -- SOFA-2 brain component >= 1 at any evaluated hour (era description, )
    LOGICAL_OR(s2_brain24 >= 1) AS any_s2_brain_pos
  FROM ev GROUP BY stay_id
),
fw AS (
  SELECT stay_id,
    COUNT(*) AS n_eval_hours_first,
    MAX(s1_std) AS fw_s1_std, MAX(s1_gcsalt) AS fw_s1_gcsalt, MAX(s1_kcr) AS fw_s1_kcr,
    MAX(s2_std) AS fw_s2_std, MAX(s2_e) AS fw_s2_e, MAX(s2_gcs15) AS fw_s2_gcs15, MAX(s2_kcr) AS fw_s2_kcr, MAX(s2_kuo22) AS fw_s2_kuo22,
    MAX(s2_g) AS fw_s2_g, MAX(s2_s) AS fw_s2_s, MAX(s2_w) AS fw_s2_w, MAX(s2_h) AS fw_s2_h, MAX(s2_i) AS fw_s2_i,
    MAX(s1_resp24) AS fw_s1_resp, MAX(s1_coag24) AS fw_s1_coag, MAX(s1_liver24) AS fw_s1_liver, MAX(s1_cv24) AS fw_s1_cv,
    MAX(s1_cns24) AS fw_s1_cns, MAX(s1_renal24) AS fw_s1_renal, MAX(s1_renal_cr24) AS fw_s1_renal_cr,
    MAX(s2_resp24) AS fw_s2_resp, MAX(s2_hemo24) AS fw_s2_hemo, MAX(s2_liver24) AS fw_s2_liver, MAX(s2_cv24) AS fw_s2_cv,
    MAX(s2_brain24) AS fw_s2_brain, MAX(s2_kidney24) AS fw_s2_kidney, MAX(s2_kidney_cr24) AS fw_s2_kidney_cr
  FROM ev WHERE hz_key = 1 GROUP BY stay_id
),
base AS (
  SELECT stay_id, hz_key, ANY_VALUE(s1_std) AS base_s1, ANY_VALUE(s2_std) AS base_s2
  FROM hzs WHERE hz = 'W' AND is_base GROUP BY 1, 2
),
indw AS (
  SELECT ev.stay_id,
    LOGICAL_OR(ev.s1_std >= COALESCE(b.base_s1, 0) + 2) AS ind_s1_std,
    LOGICAL_OR(ev.s2_std >= COALESCE(b.base_s2, 0) + 2) AS ind_s2_std
  FROM ev LEFT JOIN base b USING (stay_id, hz_key)
  GROUP BY ev.stay_id
),
base_first AS (SELECT stay_id, base_s1 AS ind_base_s1, base_s2 AS ind_base_s2 FROM base WHERE hz_key = 1),
meas AS (
  SELECT h1.stay_id,
    MAX(h1.m1_resp) = 1 AND MAX(h1.m1_coag) = 1 AND MAX(h1.m1_liver) = 1 AND MAX(h1.m1_cv) = 1
      AND MAX(h1.m1_cns) = 1 AND MAX(h1.m1_renal) = 1 AS s1_complete_d1,
    MAX(h2.m2_resp) = 1 AND MAX(h2.m2_hemo) = 1 AND MAX(h2.m2_liver) = 1 AND MAX(h2.m2_cv) = 1
      AND MAX(h2.m2_brain) = 1 AND MAX(h2.m2_kidney) = 1 AS s2_complete_d1,
    MAX(h1.m1_resp) = 1 AND MAX(h2.m2_resp) = 1 AS resp_measured_both_d1
  FROM `@BUILD@.sofa1_hourly` h1
  JOIN `@BUILD@.sofa2_hourly` h2 USING (stay_id, hr)
  WHERE hr BETWEEN 0 AND 23
  GROUP BY h1.stay_id
)
SELECT e.subject_id, e.hadm_id, e.stay_id, e.first_elig_stay, e.elig_old, e.first_elig_stay_old,
  COALESCE(a.n_eval_hours_old, 0) = 0 AS not_assessable_old,
  COALESCE(a.any_s2_std AND NOT a.any_s2_1v, FALSE) AS one_volume_s2,
  e.anchor_year_group,
  CASE WHEN e.anchor_year_group IN ('2008 - 2010', '2011 - 2013') THEN '2008-2013'
       WHEN e.anchor_year_group = '2014 - 2016' THEN '2014-2016' ELSE '2017-2022' END AS era,
  em.hadm_id IS NOT NULL AS emar_adm,
  e.intime_hr, e.outtime_hr, e.icu_intime, e.icu_outtime,
  e.first_susp_time, e.first_susp_hr, e.n_windows,
  e.first_susp_hr < 0 AS suspicion_before_icu, e.first_susp_hr >= 48 AS late_suspicion,
  COALESCE(a.n_eval_hours_any, 0) = 0 AS not_assessable_any,
  COALESCE(f.n_eval_hours_first, 0) = 0 AS not_assessable_first,
  -- landmark set: alive and in the ICU at the end of hour 23, first suspicion by hour 23
  (e.icu_outtime >= DATETIME_ADD(e.intime_hr, INTERVAL 24 HOUR)
   AND (adm.deathtime IS NULL OR adm.deathtime >= DATETIME_ADD(e.intime_hr, INTERVAL 24 HOUR))
   AND e.first_susp_hr <= 23) AS landmark,
  -- complete-case set: all six components measured on day 1 under both scores
  COALESCE(m.s1_complete_d1 AND m.s2_complete_d1, FALSE) AS complete_case,
  d.* EXCEPT (stay_id),
  a.* EXCEPT (stay_id),
  f.* EXCEPT (stay_id),
  i.ind_s1_std, i.ind_s2_std, bf.ind_base_s1, bf.ind_base_s2, bf.ind_base_s1 IS NOT NULL AS ind_baseline_available,
  e.hospital_expire_flag,
  CAST(adm.deathtime IS NOT NULL AND adm.deathtime BETWEEN e.icu_intime AND DATETIME_ADD(e.icu_outtime, INTERVAL 6 HOUR) AS INT64) AS icu_death,
  e.icu_los_days, DATETIME_DIFF(e.dischtime, e.admittime, HOUR) / 24.0 AS hosp_los_days,
  e.age, e.gender, e.first_careunit, e.first_careunit LIKE '%Cardiac Vascular%' AS cvicu,
  (e.icu_intime < adm.admittime OR COALESCE(adm.deathtime > adm.dischtime, FALSE) OR adm.dischtime < adm.admittime) AS ts_inconsistent,
  COALESCE(m.resp_measured_both_d1, FALSE) AS resp_measured_both_d1,
  CAST(e.hospital_expire_flag = 1 OR COALESCE(adm.discharge_location = 'HOSPICE', FALSE) AS INT64) AS death_or_hospice,
  CAST(pt2.dod IS NOT NULL AND pt2.dod <= DATE_ADD(DATE(e.icu_intime), INTERVAL 30 DAY) AS INT64) AS death_30d
FROM elig e
LEFT JOIN d1 d USING (stay_id)
LEFT JOIN anyw a USING (stay_id)
LEFT JOIN fw f USING (stay_id)
LEFT JOIN indw i USING (stay_id)
LEFT JOIN base_first bf USING (stay_id)
LEFT JOIN meas m USING (stay_id)
LEFT JOIN `physionet-data.mimiciv_3_1_hosp.admissions` adm ON adm.hadm_id = e.hadm_id
LEFT JOIN `physionet-data.mimiciv_3_1_hosp.patients` pt2 ON pt2.subject_id = e.subject_id
LEFT JOIN (SELECT DISTINCT hadm_id FROM `physionet-data.mimiciv_3_1_hosp.emar` WHERE hadm_id IS NOT NULL) em ON em.hadm_id = e.hadm_id;
