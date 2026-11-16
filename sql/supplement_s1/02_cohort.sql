-- Supplementary Methods S1: one row per analysed stay (first eligible stay of each patient, s1_elig) with, for SOFA
-- (s1) and SOFA-2 (s2) separately: the baseline (lowest 24-hour total over ICU hours 23-47, earliest hour on ties, with
-- its components) and the hour-47 baseline (sensitivity); Sepsis-3 (total >= baseline + 2 at an evaluated window hour)
-- and onset (first such hour); the rise at onset by component; the other score's rise at this score's onset hour; the
-- largest rise within the window; outcomes counted from the landmark; organ support (invasive ventilation,
-- vasopressors, renal replacement therapy) started after the suspicion time, with no support of the same type in the
-- 24 hours before that start.
CREATE OR REPLACE TABLE `@BUILD@.s1_cohort` AS
WITH el AS (SELECT * FROM `@BUILD@.s1_elig` WHERE first_elig_stay),
long AS (
  SELECT stay_id, hz, hr, 's1' AS sc, s1_std AS tot, s1_resp24 AS resp, s1_coag24 AS coag, s1_liver24 AS liver,
         s1_cv24 AS cv, s1_cns24 AS brain, s1_renal24 AS kidney
  FROM `@BUILD@.s1_scores` WHERE evaluate
  UNION ALL
  SELECT stay_id, hz, hr, 's2', s2_std, s2_resp24, s2_hemo24, s2_liver24, s2_cv24, s2_brain24, s2_kidney24
  FROM `@BUILD@.s1_scores` WHERE evaluate
),
base AS (
  SELECT stay_id, sc, COUNT(*) AS n_base_hours,
    ARRAY_AGG(STRUCT(hr, tot, resp, coag, liver, cv, brain, kidney) ORDER BY tot, hr LIMIT 1)[OFFSET(0)] AS bmin,
    ARRAY_AGG(IF(hr = 47, STRUCT(hr, tot, resp, coag, liver, cv, brain, kidney), NULL) IGNORE NULLS LIMIT 1)[SAFE_OFFSET(0)] AS b47
  FROM long WHERE hz = 'B' GROUP BY 1, 2
),
win AS (SELECT l.*, b.bmin, b.b47 FROM long l JOIN base b USING (stay_id, sc) WHERE l.hz = 'E'),
ws AS (
  SELECT stay_id, sc, COUNT(*) AS n_eval,
    MIN(IF(tot >= bmin.tot + 2, hr, NULL)) AS onset_hr,
    MIN(IF(tot >= b47.tot + 2, hr, NULL)) AS onset_hr47,
    MAX(tot - bmin.tot) AS max_rise, MAX(tot - b47.tot) AS max_rise47
  FROM win GROUP BY 1, 2
),
on_c AS (  -- rise at onset, total and by component (component at onset minus component at the baseline hour)
  SELECT w.stay_id, w.sc, w.tot - w.bmin.tot AS d_tot, w.resp - w.bmin.resp AS d_resp, w.coag - w.bmin.coag AS d_coag,
         w.liver - w.bmin.liver AS d_liver, w.cv - w.bmin.cv AS d_cv, w.brain - w.bmin.brain AS d_brain,
         w.kidney - w.bmin.kidney AS d_kidney
  FROM win w JOIN ws USING (stay_id, sc) WHERE w.hr = ws.onset_hr
),
cross_d AS (  -- the other score's rise, against its own baseline, at this score's onset hour
  SELECT ws.stay_id, ws.sc, o.tot - o.bmin.tot AS other_d
  FROM ws JOIN win o ON o.stay_id = ws.stay_id AND o.sc != ws.sc AND o.hr = ws.onset_hr
),
per AS (
  SELECT b.stay_id, b.sc, b.n_base_hours, b.bmin, b.b47, ws.n_eval, ws.onset_hr, ws.onset_hr47, ws.max_rise,
         ws.max_rise47, oc.d_tot, oc.d_resp, oc.d_coag, oc.d_liver, oc.d_cv, oc.d_brain, oc.d_kidney, cd.other_d
  FROM base b LEFT JOIN ws USING (stay_id, sc) LEFT JOIN on_c oc USING (stay_id, sc) LEFT JOIN cross_d cd USING (stay_id, sc)
),
sup AS (
  SELECT 'vent' AS k, stay_id, starttime, endtime FROM `@UPSTREAM@.ventilation`
  WHERE ventilation_status = 'InvasiveVent'
  UNION ALL SELECT 'vaso', stay_id, starttime, endtime FROM `@UPSTREAM@.norepinephrine`
  UNION ALL SELECT 'vaso', stay_id, starttime, endtime FROM `@UPSTREAM@.epinephrine`
  UNION ALL SELECT 'vaso', stay_id, starttime, endtime FROM `@UPSTREAM@.dopamine`
  UNION ALL SELECT 'vaso', stay_id, starttime, endtime FROM `@UPSTREAM@.phenylephrine`
  UNION ALL SELECT 'vaso', stay_id, starttime, endtime FROM `@UPSTREAM@.vasopressin`
  UNION ALL SELECT 'rrt', stay_id, charttime, charttime FROM `@UPSTREAM@.rrt` WHERE dialysis_active = 1
),
cand AS (  -- support starts after the suspicion time, in the ICU
  SELECT el.stay_id, a.k, a.starttime
  FROM el JOIN sup a ON a.stay_id = el.stay_id AND a.starttime > el.susp_time AND a.starttime <= el.icu_outtime
),
newsup AS (  -- a start counts when there was no support of the same type in the 24 hours before that start
  SELECT stay_id, k FROM (
    SELECT c.stay_id, c.k, c.starttime, LOGICAL_OR(b.stay_id IS NOT NULL) AS had_prior
    FROM cand c LEFT JOIN sup b
      ON b.stay_id = c.stay_id AND b.k = c.k AND b.starttime < c.starttime
     AND b.endtime >= DATETIME_SUB(c.starttime, INTERVAL 24 HOUR)
    GROUP BY 1, 2, 3)
  WHERE NOT had_prior
  GROUP BY 1, 2
),
ns AS (
  SELECT stay_id, LOGICAL_OR(k = 'vent') AS new_vent, LOGICAL_OR(k = 'vaso') AS new_vasopressor, LOGICAL_OR(k = 'rrt') AS new_rrt
  FROM newsup GROUP BY 1
)
SELECT el.subject_id, el.hadm_id, el.stay_id, el.susp_time, el.susp_hr, el.landmark_time,
  COALESCE(p1.n_eval, 0) = 0 OR COALESCE(p2.n_eval, 0) = 0 AS not_assessable,
  p1.n_base_hours AS s1_n_base_hours, p2.n_base_hours AS s2_n_base_hours, p1.n_eval AS n_eval_hours,
  -- SOFA
  p1.bmin.tot AS s1_base, p1.bmin.hr AS s1_base_hr, p1.b47.tot AS s1_base47,
  p1.onset_hr IS NOT NULL AS s1_pos, p1.onset_hr AS s1_onset_hr, p1.onset_hr - el.susp_hr AS s1_onset_from_susp_hr,
  p1.onset_hr47 IS NOT NULL AS s1_pos47, p1.max_rise AS s1_max_rise, p1.max_rise47 AS s1_max_rise47,
  p1.d_tot AS s1_d, p1.d_resp AS s1_d_resp, p1.d_coag AS s1_d_coag, p1.d_liver AS s1_d_liver, p1.d_cv AS s1_d_cv,
  p1.d_brain AS s1_d_brain, p1.d_kidney AS s1_d_kidney, p1.other_d AS s2_d_at_s1_onset,
  -- SOFA-2
  p2.bmin.tot AS s2_base, p2.bmin.hr AS s2_base_hr, p2.b47.tot AS s2_base47,
  p2.onset_hr IS NOT NULL AS s2_pos, p2.onset_hr AS s2_onset_hr, p2.onset_hr - el.susp_hr AS s2_onset_from_susp_hr,
  p2.onset_hr47 IS NOT NULL AS s2_pos47, p2.max_rise AS s2_max_rise, p2.max_rise47 AS s2_max_rise47,
  p2.d_tot AS s2_d, p2.d_resp AS s2_d_resp, p2.d_coag AS s2_d_coag, p2.d_liver AS s2_d_liver, p2.d_cv AS s2_d_cv,
  p2.d_brain AS s2_d_brain, p2.d_kidney AS s2_d_kidney, p2.other_d AS s1_d_at_s2_onset,
  -- outcomes, counted from the landmark
  cb.hospital_expire_flag,
  CAST(adm.deathtime IS NOT NULL AND adm.deathtime BETWEEN cb.icu_intime AND DATETIME_ADD(cb.icu_outtime, INTERVAL 6 HOUR) AS INT64) AS icu_death,
  DATETIME_DIFF(cb.icu_outtime, el.landmark_time, MINUTE) / 1440.0 AS icu_los_after_landmark_days,
  DATETIME_DIFF(cb.dischtime, el.landmark_time, MINUTE) / 1440.0 AS hosp_los_after_landmark_days,
  cb.age, cb.gender, cb.first_careunit, cb.first_careunit LIKE '%Cardiac Vascular%' AS cvicu,
  COALESCE(ns.new_vent, FALSE) AS new_vent, COALESCE(ns.new_vasopressor, FALSE) AS new_vasopressor,
  COALESCE(ns.new_rrt, FALSE) AS new_rrt
FROM el
JOIN `@UPSTREAM@.cohort_base` cb USING (stay_id)
LEFT JOIN `physionet-data.mimiciv_3_1_hosp.admissions` adm ON adm.hadm_id = el.hadm_id
LEFT JOIN per p1 ON p1.stay_id = el.stay_id AND p1.sc = 's1'
LEFT JOIN per p2 ON p2.stay_id = el.stay_id AND p2.sc = 's2'
LEFT JOIN ns ON ns.stay_id = el.stay_id;
