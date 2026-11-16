-- Supplementary Methods S1: ICU-onset sepsis measured from an observed ICU baseline. Eligible stays, one row per stay
-- of the base cohort that is alive and in the ICU at the end of ICU hour 47 (the landmark), with the eligibility flags.
-- Hours follow the main analysis: hour h = DATETIME_DIFF(time, intime_hr, HOUR), and the end of hour h is
-- intime_hr + (h + 1) hours (as the landmark of the main analysis, cohort.sql).
--   early_susp:      any suspicion of infection of the admission up to the end of ICU hour 47 (excluded);
--   susp_time:       otherwise the first suspicion, eligible when in hours 48 to 239 and before ICU discharge;
--   first_elig_stay: first eligible stay of the patient (the analysed stay).
CREATE OR REPLACE TABLE `@BUILD@.s1_elig` AS
WITH pop AS (
  SELECT c.subject_id, c.hadm_id, c.stay_id, c.icu_intime, c.icu_outtime, t.intime_hr, t.outtime_hr,
         DATETIME_ADD(t.intime_hr, INTERVAL 48 HOUR) AS landmark_time
  FROM `@UPSTREAM@.cohort_base` c JOIN `@UPSTREAM@.icustay_times` t USING (stay_id)
  WHERE c.age >= 18
),
lm AS (  -- alive and in the ICU at the end of hour 47
  SELECT p.*
  FROM pop p JOIN `physionet-data.mimiciv_3_1_hosp.admissions` adm ON adm.hadm_id = p.hadm_id
  WHERE p.icu_outtime >= p.landmark_time AND (adm.deathtime IS NULL OR adm.deathtime >= p.landmark_time)
),
susp AS (
  SELECT hadm_id, suspected_infection_time
  FROM `@UPSTREAM@.suspicion_of_infection`
  WHERE suspected_infection = 1 AND hadm_id IS NOT NULL
),
s AS (
  SELECT l.stay_id,
    LOGICAL_OR(s.suspected_infection_time < l.landmark_time) AS early_susp,
    MIN(IF(s.suspected_infection_time >= l.landmark_time, s.suspected_infection_time, NULL)) AS susp_time
  FROM lm l JOIN susp s USING (hadm_id)
  GROUP BY 1
),
e AS (
  SELECT l.*, COALESCE(s.early_susp, FALSE) AS early_susp, s.susp_time,
         DATETIME_DIFF(s.susp_time, l.intime_hr, HOUR) AS susp_hr,
         (NOT COALESCE(s.early_susp, FALSE) AND s.susp_time IS NOT NULL
          AND DATETIME_DIFF(s.susp_time, l.intime_hr, HOUR) BETWEEN 48 AND 239
          AND s.susp_time <= l.outtime_hr) AS eligible
  FROM lm l LEFT JOIN s USING (stay_id)
)
SELECT e.*,
  e.eligible AND ROW_NUMBER() OVER (PARTITION BY e.subject_id, e.eligible ORDER BY e.icu_intime, e.stay_id) = 1
    AS first_elig_stay
FROM e;
