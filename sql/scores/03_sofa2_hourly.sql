-- Hourly SOFA-2 with the brain component (GCS bands, motor rule when verbal is not assessable, pre-sedation score
-- carried forward through sedated hours with a 72-hour cap, delirium point), the cardiovascular component of
-- 02_cv_hourly.sql, 24-hour component maxima, and the per-hour inputs needed for the horizon scores and variants.
CREATE OR REPLACE TABLE `@BUILD@.sofa2_hourly` AS
WITH g AS (
  SELECT s.stay_id, s.hr, MIN(g.gcs) AS gcs_total_min, MAX(g.gcs_unable) AS gcs_unable_any2, MAX(g.gcs_motor) AS gcs_motor_best
  FROM `@BUILD@.sofa2_withintime` s
  JOIN `@UPSTREAM@.gcs` g ON g.stay_id = s.stay_id AND g.charttime > s.starttime AND g.charttime <= s.endtime
  GROUP BY 1,2
),
h AS (
  SELECT s.*, g.gcs_total_min,
    COALESCE(g.gcs_unable_any2, s.gcs_unable_any, 0) AS unable,
    g.gcs_motor_best AS motor,
    CASE
      WHEN COALESCE(g.gcs_unable_any2, s.gcs_unable_any, 0) = 0 AND g.gcs_total_min IS NOT NULL THEN
        CASE WHEN g.gcs_total_min >= 15 THEN 0 WHEN g.gcs_total_min >= 13 THEN 1 WHEN g.gcs_total_min >= 9 THEN 2 WHEN g.gcs_total_min >= 6 THEN 3 ELSE 4 END
      WHEN g.gcs_motor_best IS NOT NULL THEN
        CASE WHEN g.gcs_motor_best >= 6 THEN 0 WHEN g.gcs_motor_best = 5 THEN 1 WHEN g.gcs_motor_best = 4 THEN 2 WHEN g.gcs_motor_best = 3 THEN 3 ELSE 4 END
      ELSE NULL END AS brain_raw,
    -- MIMIC SOFA-1 convention: intubated -> GCS 15 -> 0
    CASE WHEN COALESCE(g.gcs_unable_any2, s.gcs_unable_any, 0) = 1 THEN 0
         WHEN g.gcs_total_min IS NOT NULL THEN
           CASE WHEN g.gcs_total_min >= 15 THEN 0 WHEN g.gcs_total_min >= 13 THEN 1 WHEN g.gcs_total_min >= 9 THEN 2 WHEN g.gcs_total_min >= 6 THEN 3 ELSE 4 END
         ELSE NULL END AS brain_raw_gcs15,
    CASE WHEN COALESCE(s.on_rrt_any,0)=1 OR COALESCE(s.rrt_eligible_flag,0)=1 THEN 4
         WHEN s.creatinine_max > 3.5 THEN 3 WHEN s.creatinine_max > 2.0 THEN 2 WHEN s.creatinine_max > 1.2 THEN 1 ELSE 0 END AS kidney_cr,
    CASE WHEN COALESCE(s.on_rrt_any,0)=1 OR COALESCE(s.rrt_eligible_flag,0)=1 THEN 4
         WHEN s.hr >= 22 AND COALESCE(s.uo_gt24_or_anuria_flag,0)=1 THEN 3
         WHEN s.hr >= 22 AND COALESCE(s.uo_gt12_flag,0)=1 THEN 2
         WHEN s.hr >= 22 AND COALESCE(s.uo_6_12_flag,0)=1 THEN 1 ELSE 0 END AS kidney_uo_late
  FROM `@BUILD@.sofa2_withintime` s
  LEFT JOIN g USING (stay_id, hr)
),
sed1 AS (
  SELECT h.*,
    MAX(IF(COALESCE(sedated_hour,0)=1, hr, NULL)) OVER (PARTITION BY stay_id ORDER BY hr ROWS BETWEEN 6 PRECEDING AND CURRENT ROW) AS prev_sed_hr,
    MIN(IF(COALESCE(sedated_hour,0)=1, hr, NULL)) OVER (PARTITION BY stay_id ORDER BY hr ROWS BETWEEN CURRENT ROW AND 6 FOLLOWING) AS next_sed_hr
  FROM h
),
sed2 AS (
  SELECT *,
    IF(COALESCE(sedated_hour,0)=1 OR (prev_sed_hr IS NOT NULL AND next_sed_hr IS NOT NULL), 1, 0) AS sed_block
  FROM sed1
),
ff AS (
  SELECT sed2.*,
    LAST_VALUE(IF(sed_block=0, brain_raw, NULL) IGNORE NULLS) OVER (PARTITION BY stay_id ORDER BY hr ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) AS brain_presed,
    LAST_VALUE(IF(sed_block=0 AND brain_raw IS NOT NULL, hr, NULL) IGNORE NULLS) OVER (PARTITION BY stay_id ORDER BY hr ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) AS brain_presed_hr
  FROM sed2
),
b AS (
  SELECT ff.*,
    GREATEST(
      CASE WHEN sed_block=1 THEN IF(hr - brain_presed_hr <= 72, COALESCE(brain_presed, 0), 0) ELSE COALESCE(brain_raw, 0) END,
      IF(COALESCE(delirium_treated_any,0)=1, 1, 0)) AS brain_v2,
    GREATEST(CASE WHEN sed_block=1 THEN 0 ELSE COALESCE(brain_raw, 0) END, IF(COALESCE(delirium_treated_any,0)=1, 1, 0)) AS brain_v2_nocarry,
    GREATEST(CASE WHEN sed_block=1 THEN 0 ELSE COALESCE(brain_raw_gcs15, 0) END, IF(COALESCE(delirium_treated_any,0)=1, 1, 0)) AS brain_gcs15,
    GREATEST(kidney_cr, kidney_uo_late) AS kidney_uo22
  FROM ff
),
w AS (
  SELECT b.*, c.cardiovascular_v2, c.cardiovascular_v2_cut, c.meanbp_min AS meanbp_min_v2, c.ne_eq, c.other_agent,
    MAX(brain_v2) OVER (PARTITION BY b.stay_id ORDER BY b.hr ROWS BETWEEN 23 PRECEDING AND CURRENT ROW) AS brain_v2_24h,
    MAX(brain_v2_nocarry) OVER (PARTITION BY b.stay_id ORDER BY b.hr ROWS BETWEEN 23 PRECEDING AND CURRENT ROW) AS brain_nc_24h,
    MAX(brain_gcs15) OVER (PARTITION BY b.stay_id ORDER BY b.hr ROWS BETWEEN 23 PRECEDING AND CURRENT ROW) AS brain_gcs15_24h,
    MAX(c.cardiovascular_v2) OVER (PARTITION BY b.stay_id ORDER BY b.hr ROWS BETWEEN 23 PRECEDING AND CURRENT ROW) AS cv_v2_24h,
    MAX(kidney_cr) OVER (PARTITION BY b.stay_id ORDER BY b.hr ROWS BETWEEN 23 PRECEDING AND CURRENT ROW) AS kidney_cr_24h,
    MAX(kidney_uo22) OVER (PARTITION BY b.stay_id ORDER BY b.hr ROWS BETWEEN 23 PRECEDING AND CURRENT ROW) AS kidney_uo22_24h
  FROM b JOIN `@BUILD@.cv_hourly` c USING (stay_id, hr)
)
SELECT stay_id, hr, starttime, endtime, intime,
       respiration, hemostasis, liver, brain AS brain_orig, brain_v2, brain_v2_nocarry, brain_gcs15,
       cardiovascular AS cardiovascular_orig, cardiovascular_v2 AS cardiovascular, kidney, kidney_cr, kidney_uo22,
       respiration_24hours, hemostasis_24hours, liver_24hours,
       brain_24hours AS brain_orig_24hours, brain_v2_24h AS brain_24hours, brain_nc_24h AS brain_nocarry_24hours, brain_gcs15_24h AS brain_gcs15_24hours,
       cardiovascular_24hours AS cardiovascular_orig_24hours, COALESCE(cv_v2_24h, 0) AS cardiovascular_24hours,
       kidney_24hours, kidney_cr_24h AS kidney_cr_24hours, kidney_uo22_24h AS kidney_uo22_24hours,
       sofa2_24hours AS sofa2_orig_24hours,
       COALESCE(respiration_24hours,0) + COALESCE(hemostasis_24hours,0) + COALESCE(liver_24hours,0) + brain_v2_24h + COALESCE(cv_v2_24h,0) + COALESCE(kidney_24hours,0) AS sofa2_24hours,
       COALESCE(respiration_24hours,0) + COALESCE(hemostasis_24hours,0) + COALESCE(liver_24hours,0) + brain_nc_24h + COALESCE(cv_v2_24h,0) + COALESCE(kidney_24hours,0) AS sofa2_nocarry_24hours,
       COALESCE(respiration_24hours,0) + COALESCE(hemostasis_24hours,0) + COALESCE(liver_24hours,0) + brain_gcs15_24h + COALESCE(cv_v2_24h,0) + COALESCE(kidney_24hours,0) AS sofa2_gcs15_24hours,
       COALESCE(respiration_24hours,0) + COALESCE(hemostasis_24hours,0) + COALESCE(liver_24hours,0) + brain_v2_24h + COALESCE(cv_v2_24h,0) + kidney_cr_24h AS sofa2_kcr_24hours,
       COALESCE(respiration_24hours,0) + COALESCE(hemostasis_24hours,0) + COALESCE(liver_24hours,0) + brain_v2_24h + COALESCE(cv_v2_24h,0) + kidney_uo22_24h AS sofa2_kuo22_24hours,
       gcs_total_min, unable, motor, sedated_hour, sed_block, delirium_treated_any, meanbp_min_v2, ne_eq, other_agent,
       -- inputs for the truncated brain scores (05_hz_scores.sql) and variant e
       brain_raw, brain_raw_gcs15, prev_sed_hr, next_sed_hr, dex_treated_any,
       -- 24-h totals without the brain component, per kidney rule
       COALESCE(respiration_24hours,0) + COALESCE(hemostasis_24hours,0) + COALESCE(liver_24hours,0) + COALESCE(cv_v2_24h,0) + COALESCE(kidney_24hours,0) AS s2_nobrain_24hours,
       COALESCE(respiration_24hours,0) + COALESCE(hemostasis_24hours,0) + COALESCE(liver_24hours,0) + COALESCE(cv_v2_24h,0) + kidney_cr_24h AS s2_nobrain_kcr_24hours,
       COALESCE(respiration_24hours,0) + COALESCE(hemostasis_24hours,0) + COALESCE(liver_24hours,0) + COALESCE(cv_v2_24h,0) + kidney_uo22_24h AS s2_nobrain_kuo22_24hours,
       -- component measured in this hour (complete-case set)
       CAST(pao2fio2ratio_min IS NOT NULL OR spo2fio2ratio_min IS NOT NULL OR COALESCE(ecmo_any,0) = 1 AS INT64) AS m2_resp,
       CAST(platelet_min IS NOT NULL AS INT64) AS m2_hemo,
       CAST(bilirubin_max IS NOT NULL AS INT64) AS m2_liver,
       CAST(meanbp_min_v2 IS NOT NULL OR COALESCE(ne_eq,0) > 0 OR COALESCE(other_agent,FALSE) OR COALESCE(mechanical_support_flag,0) = 1 AS INT64) AS m2_cv,
       CAST((unable = 0 AND gcs_total_min IS NOT NULL) OR motor IS NOT NULL AS INT64) AS m2_brain,
       CAST(creatinine_max IS NOT NULL OR uo_6_12_flag IS NOT NULL OR uo_gt12_flag IS NOT NULL OR uo_gt24_or_anuria_flag IS NOT NULL OR COALESCE(on_rrt_any,0) = 1 AS INT64) AS m2_kidney,
       -- inputs for the horizon scores (05_hz_scores.sql)
       COALESCE(respiration_24hours,0) + COALESCE(hemostasis_24hours,0) + COALESCE(liver_24hours,0) AS s2_rhl_24hours,
       cardiovascular_v2_cut AS cardiovascular_cut,
       kidney_o, kidney_b, kidney_cr_o, kidney_cr_b, kidney_uo22_o, kidney_uo22_b, kidney_go, kidney_gb, kidney_so, kidney_sb, kidney_pre24,
       weight_own_kg, weight_ba_kg, weight_backapplied, wt_adm_hr, uo_n_6h, uo_n_12h, uo_n_24h, cov_6h, cov_12h, cov_24h,
       uo_6_12_flag, uo_gt12_flag, uo_gt24_or_anuria_flag
FROM w;
