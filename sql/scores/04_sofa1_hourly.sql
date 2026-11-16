-- Hourly inputs for SOFA (mimic-code sofa_withintime, unchanged scoring), plus the SOFA brain band applied to the
-- SOFA-2 GCS handling (variant a) and the SOFA renal component by creatinine only (variant f).
CREATE OR REPLACE TABLE `@BUILD@.sofa1_hourly` AS
WITH a AS (SELECT * FROM `@UPSTREAM@.sofa_withintime`),
s2 AS (SELECT stay_id, hr, gcs_total_min, unable, motor FROM `@BUILD@.sofa2_hourly`),
h AS (
  SELECT a.stay_id, a.hr, a.endtime, a.sofa_24hours, a.cns_24hours, a.renal_24hours,
    CASE WHEN s2.unable = 0 AND s2.gcs_total_min IS NOT NULL THEN
           CASE WHEN s2.gcs_total_min >= 15 THEN 0 WHEN s2.gcs_total_min >= 13 THEN 1 WHEN s2.gcs_total_min >= 10 THEN 2 WHEN s2.gcs_total_min >= 6 THEN 3 ELSE 4 END
         WHEN s2.motor IS NOT NULL THEN
           CASE WHEN s2.motor >= 6 THEN 0 WHEN s2.motor = 5 THEN 1 WHEN s2.motor = 4 THEN 2 WHEN s2.motor = 3 THEN 3 ELSE 4 END
         ELSE NULL END AS cns_raw,
    CASE WHEN a.creatinine_max >= 5.0 THEN 4 WHEN a.creatinine_max >= 3.5 THEN 3 WHEN a.creatinine_max >= 2.0 THEN 2
         WHEN a.creatinine_max >= 1.2 THEN 1 WHEN a.creatinine_max IS NULL THEN NULL ELSE 0 END AS renal_cr,
    CAST(a.pao2fio2ratio_vent IS NOT NULL OR a.pao2fio2ratio_novent IS NOT NULL AS INT64) AS m1_resp,
    CAST(a.platelet_min IS NOT NULL AS INT64) AS m1_coag,
    CAST(a.bilirubin_max IS NOT NULL AS INT64) AS m1_liver,
    CAST(COALESCE(a.meanbp_min, a.rate_dopamine, a.rate_dobutamine, a.rate_epinephrine, a.rate_norepinephrine) IS NOT NULL AS INT64) AS m1_cv,
    CAST(a.gcs_min IS NOT NULL AS INT64) AS m1_cns,
    CAST(a.creatinine_max IS NOT NULL OR a.uo_24hr IS NOT NULL AS INT64) AS m1_renal
  FROM a LEFT JOIN s2 USING (stay_id, hr)
)
SELECT h.*,
  COALESCE(MAX(renal_cr) OVER (PARTITION BY stay_id ORDER BY hr ROWS BETWEEN 23 PRECEDING AND CURRENT ROW), 0) AS renal_cr_24hours,
  sofa_24hours - COALESCE(cns_24hours,0) AS s1_nocns_24hours,
  sofa_24hours - COALESCE(renal_24hours,0)
    + COALESCE(MAX(renal_cr) OVER (PARTITION BY stay_id ORDER BY hr ROWS BETWEEN 23 PRECEDING AND CURRENT ROW), 0) AS sofa_kcr_24hours
FROM h;
