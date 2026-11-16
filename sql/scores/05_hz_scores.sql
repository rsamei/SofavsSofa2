-- Scores at each evaluation horizon (hz): 'D' = ICU day 1 (end of hour 23 or last scored hour), 'W' = every hour of a
-- suspicion-of-infection window, with look-ahead truncated at the horizon. One row per stay, horizon and hour, with
-- SOFA, SOFA-2 and every implementation variant.
CREATE OR REPLACE TABLE `@BUILD@.hz_scores` AS
WITH pop AS (
  SELECT c.stay_id, c.hadm_id, t.intime_hr, t.outtime_hr
  FROM `@UPSTREAM@.cohort_base` c JOIN `@UPSTREAM@.icustay_times` t USING (stay_id)
  WHERE c.age >= 18
),
susp AS (  -- by admission; old_rule = the stay_id match also holds this suspicion time for the stay
  SELECT hadm_id, suspected_infection_time, ARRAY_AGG(DISTINCT stay_id IGNORE NULLS) AS s_stays
  FROM `@UPSTREAM@.suspicion_of_infection`
  WHERE suspected_infection = 1 AND hadm_id IS NOT NULL
  GROUP BY 1, 2
),
winx AS (
  SELECT p.stay_id, p.intime_hr, p.outtime_hr, s.suspected_infection_time AS susp_time,
         DATETIME_SUB(s.suspected_infection_time, INTERVAL 48 HOUR) AS w_start,
         DATETIME_ADD(s.suspected_infection_time, INTERVAL 24 HOUR) AS w_end,
         ROW_NUMBER() OVER (PARTITION BY p.stay_id ORDER BY s.suspected_infection_time) AS hz_key,
         p.stay_id IN UNNEST(s.s_stays) AS old_rule
  FROM pop p JOIN susp s USING (hadm_id)
  WHERE DATETIME_ADD(s.suspected_infection_time, INTERVAL 24 HOUR) >= p.intime_hr
    AND DATETIME_SUB(s.suspected_infection_time, INTERVAL 48 HOUR) <= p.outtime_hr
),
s2 AS (
  SELECT stay_id, hr, starttime, endtime, sedated_hour, prev_sed_hr, next_sed_hr, brain_raw, brain_raw_gcs15,
         gcs_total_min AS gcs_hr_min, unable AS gcs_unable_hr,
         delirium_treated_any, dex_treated_any,
         s2_nobrain_24hours, s2_nobrain_kcr_24hours, s2_nobrain_kuo22_24hours,
         respiration_24hours AS s2_resp24, hemostasis_24hours AS s2_hemo24, liver_24hours AS s2_liver24,
         s2_rhl_24hours, cardiovascular AS cv_full, cardiovascular_cut AS cv_cut,
         kidney_o, kidney_b, kidney_cr_o, kidney_cr_b, kidney_uo22_o, kidney_uo22_b, kidney_go, kidney_gb, kidney_so, kidney_sb, kidney_pre24,
         weight_own_kg, weight_ba_kg, weight_backapplied, wt_adm_hr, uo_n_6h, uo_n_12h, uo_n_24h, uo_6_12_flag, uo_gt12_flag, uo_gt24_or_anuria_flag
  FROM `@BUILD@.sofa2_hourly`
),
ap_items AS (
  SELECT itemid FROM `physionet-data.mimiciv_3_1_icu.d_items`
  WHERE REGEXP_CONTAINS(LOWER(label), r'haloperidol|quetiapine|olanzapine|risperidone|ziprasidone|chlorpromazine')
),
ie_hr AS (  -- variant i: primary-list antipsychotic in inputevents overlapping the hour
  SELECT DISTINCT s2.stay_id, s2.hr
  FROM s2 JOIN `physionet-data.mimiciv_3_1_icu.inputevents` i
    ON i.stay_id = s2.stay_id AND i.amount > 0
   AND i.starttime < s2.endtime AND COALESCE(i.endtime, i.starttime) >= s2.starttime
  WHERE i.itemid IN (SELECT itemid FROM ap_items)
),
s1 AS (
  SELECT h.stay_id, h.hr, h.cns_raw, h.sofa_24hours, h.s1_nocns_24hours, h.sofa_kcr_24hours,
         h.cns_24hours AS s1_cns24, h.renal_24hours AS s1_renal24, h.renal_cr_24hours AS s1_renal_cr24,
         a.respiration_24hours AS s1_resp24, a.coagulation_24hours AS s1_coag24,
         a.liver_24hours AS s1_liver24, a.cardiovascular_24hours AS s1_cv24
  FROM `@BUILD@.sofa1_hourly` h
  JOIN `@UPSTREAM@.sofa_withintime` a USING (stay_id, hr)
),
s1p AS (  -- original SOFA definitions (Vincent 1996) scored on the SOFA-2 pipeline's hourly inputs
  SELECT w.stay_id, w.hr,
    CASE WHEN w.pao2fio2ratio_min < 100 AND w.advanced_vent_any = 1 THEN 4
         WHEN w.pao2fio2ratio_min < 200 AND w.advanced_vent_any = 1 THEN 3
         WHEN w.pao2fio2ratio_min < 300 THEN 2 WHEN w.pao2fio2ratio_min < 400 THEN 1 ELSE 0 END AS p_resp,
    CASE WHEN w.platelet_min < 20 THEN 4 WHEN w.platelet_min < 50 THEN 3 WHEN w.platelet_min < 100 THEN 2
         WHEN w.platelet_min < 150 THEN 1 ELSE 0 END AS p_coag,
    CASE WHEN w.bilirubin_max >= 12 THEN 4 WHEN w.bilirubin_max >= 6 THEN 3 WHEN w.bilirubin_max >= 2 THEN 2
         WHEN w.bilirubin_max >= 1.2 THEN 1 ELSE 0 END AS p_liver,
    GREATEST(CASE WHEN w.creatinine_max >= 5 THEN 4 WHEN w.creatinine_max >= 3.5 THEN 3 WHEN w.creatinine_max >= 2 THEN 2
                  WHEN w.creatinine_max >= 1.2 THEN 1 ELSE 0 END,
             CASE WHEN w.hr >= 23 AND w.iv_24h AND w.uo_24h_ml < 200 THEN 4
                  WHEN w.hr >= 23 AND w.iv_24h AND w.uo_24h_ml < 500 THEN 3 ELSE 0 END) AS p_renal,
    CASE WHEN COALESCE(c.dopa, 0) > 15 OR COALESCE(c.epi, 0) > 0.1 OR COALESCE(c.ne, 0) > 0.1 THEN 4
         WHEN COALESCE(c.dopa, 0) > 5 OR COALESCE(c.epi, 0) > 0 OR COALESCE(c.ne, 0) > 0 THEN 3
         WHEN COALESCE(c.dopa, 0) > 0 OR COALESCE(c.dobu, 0) > 0 THEN 2
         WHEN c.meanbp_min < 70 THEN 1 ELSE 0 END AS p_cv_full,
    CASE WHEN COALESCE(c.dopa_cut, 0) > 15 OR COALESCE(c.epi_cut, 0) > 0.1 OR COALESCE(c.ne_cut, 0) > 0.1 THEN 4
         WHEN COALESCE(c.dopa_cut, 0) > 5 OR COALESCE(c.epi_cut, 0) > 0 OR COALESCE(c.ne_cut, 0) > 0 THEN 3
         WHEN COALESCE(c.dopa_cut, 0) > 0 OR COALESCE(c.dobu_cut, 0) > 0 THEN 2
         WHEN c.meanbp_min < 70 THEN 1 ELSE 0 END AS p_cv_cut
  FROM `@BUILD@.sofa2_withintime` w JOIN `@BUILD@.cv_hourly` c USING (stay_id, hr)
),
rows_d1 AS (
  SELECT s2.stay_id, 'D1' AS hz, 0 AS hz_key, 23 AS t_cut, s2.hr, CAST(NULL AS BOOL) AS old_rule, TRUE AS evaluate, FALSE AS is_base
  FROM s2 JOIN (SELECT DISTINCT stay_id FROM pop) USING (stay_id)
  WHERE s2.hr BETWEEN 0 AND 23
),
wcut AS (
  SELECT w.stay_id, w.hz_key, MAX(s2.hr) AS t_cut
  FROM winx w JOIN s2 ON s2.stay_id = w.stay_id AND s2.endtime <= w.w_end
  GROUP BY 1, 2
),
rows_w AS (
  SELECT w.stay_id, 'W' AS hz, w.hz_key, c.t_cut, s2.hr, w.old_rule,
         (s2.endtime BETWEEN w.w_start AND w.w_end AND s2.endtime BETWEEN w.intime_hr AND w.outtime_hr) AS evaluate,
         (s2.endtime >= DATETIME_SUB(w.w_start, INTERVAL 1 HOUR) AND s2.endtime < w.w_start AND s2.hr >= 23) AS is_base
  FROM winx w JOIN wcut c USING (stay_id, hz_key)
  JOIN s2 ON s2.stay_id = w.stay_id AND s2.endtime <= w.w_end
         AND s2.endtime >= DATETIME_SUB(w.w_start, INTERVAL 100 HOUR)
),
r AS (SELECT * FROM rows_d1 UNION ALL SELECT * FROM rows_w),
x AS (
  SELECT r.*, s2.* EXCEPT (stay_id, hr), s1.* EXCEPT (stay_id, hr), s1p.* EXCEPT (stay_id, hr),
    IF(COALESCE(s2.sedated_hour, 0) = 1
       OR (s2.prev_sed_hr IS NOT NULL AND s2.next_sed_hr IS NOT NULL AND s2.next_sed_hr <= r.t_cut), 1, 0) AS sedb,
    (COALESCE(s2.weight_backapplied, 0) = 0 OR s2.wt_adm_hr <= r.t_cut) AS wt_ok,
    ie_hr.stay_id IS NOT NULL AS ap_ie
  FROM r JOIN s2 USING (stay_id, hr) JOIN s1 USING (stay_id, hr) JOIN s1p USING (stay_id, hr)
  LEFT JOIN ie_hr USING (stay_id, hr)
),
ff AS (
  SELECT x.*,
    LAST_VALUE(IF(sedb = 0, brain_raw, NULL) IGNORE NULLS) OVER wp AS b_presed,
    LAST_VALUE(IF(sedb = 0 AND brain_raw IS NOT NULL, hr, NULL) IGNORE NULLS) OVER wp AS b_presed_hr,
    LAST_VALUE(IF(sedb = 0, cns_raw, NULL) IGNORE NULLS) OVER wp AS c_presed,
    LAST_VALUE(IF(sedb = 0 AND cns_raw IS NOT NULL, hr, NULL) IGNORE NULLS) OVER wp AS c_presed_hr,
    -- original SOFA brain rule: only a complete (assessable) GCS is carried forward, no best-motor substitute
    LAST_VALUE(IF(sedb = 0, cns_full, NULL) IGNORE NULLS) OVER wp AS o_prev,
    LAST_VALUE(IF(sedb = 0 AND cns_full IS NOT NULL, hr, NULL) IGNORE NULLS) OVER wp AS o_prev_hr
  FROM (
    SELECT x.*,
      CASE WHEN COALESCE(gcs_unable_hr, 0) = 0 AND gcs_hr_min IS NOT NULL THEN
             CASE WHEN gcs_hr_min >= 15 THEN 0 WHEN gcs_hr_min >= 13 THEN 1 WHEN gcs_hr_min >= 10 THEN 2
                  WHEN gcs_hr_min >= 6 THEN 3 ELSE 4 END
           ELSE NULL END AS cns_full
    FROM x) x
  WINDOW wp AS (PARTITION BY stay_id, hz, hz_key ORDER BY hr ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING)
),
b AS (
  SELECT ff.*,
    CASE WHEN sedb = 1 THEN IF(hr - b_presed_hr <= 72, COALESCE(b_presed, 0), 0) ELSE COALESCE(brain_raw, 0) END AS b_core,
    CASE WHEN sedb = 1 THEN IF(hr - c_presed_hr <= 72, COALESCE(c_presed, 0), 0) ELSE COALESCE(cns_raw, 0) END AS cns_alt,
    CASE WHEN sedb = 1 THEN 0 ELSE COALESCE(brain_raw_gcs15, 0) END AS b_gcs15_core,
    -- cns_orig: sedated or not-assessable hours take the last complete GCS (72-hour cap, none -> 0);
    -- cns_orig_n: as cns_orig, but a not-assessable hour outside sedation is scored as normal (0)
    CASE WHEN sedb = 1 OR COALESCE(gcs_unable_hr, 0) = 1 THEN IF(hr - o_prev_hr <= 72, COALESCE(o_prev, 0), 0)
         ELSE COALESCE(cns_full, 0) END AS cns_orig,
    CASE WHEN sedb = 1 THEN IF(hr - o_prev_hr <= 72, COALESCE(o_prev, 0), 0)
         WHEN COALESCE(gcs_unable_hr, 0) = 1 THEN 0
         ELSE COALESCE(cns_full, 0) END AS cns_orig_n
  FROM ff
),
b2 AS (
  SELECT b.*,
    GREATEST(b_core, IF(COALESCE(delirium_treated_any, 0) = 1, 1, 0)) AS brain_std,
    GREATEST(b_core, IF(COALESCE(delirium_treated_any, 0) = 1 OR COALESCE(dex_treated_any, 0) = 1, 1, 0)) AS brain_e,
    GREATEST(b_gcs15_core, IF(COALESCE(delirium_treated_any, 0) = 1, 1, 0)) AS brain_gcs15,
    b_core AS brain_h,
    b_gcs15_core AS brain_gcs15_h,
    GREATEST(b_core, IF(COALESCE(delirium_treated_any, 0) = 1 OR ap_ie, 1, 0)) AS brain_i,
    -- per-horizon choices
    IF(wt_ok, kidney_b, kidney_o) AS k_std,
    IF(wt_ok, kidney_cr_b, kidney_cr_o) AS k_cr,
    IF(wt_ok, kidney_uo22_b, kidney_uo22_o) AS k_uo22,
    IF(wt_ok, kidney_gb, kidney_go) AS k_g,
    IF(wt_ok, kidney_sb, kidney_so) AS k_s,
    -- variant w: weights outside 30-300 kg are not used, so those hours lose the urine criteria
    IF(wt_ok, IF(weight_ba_kg BETWEEN 30 AND 300, kidney_b, kidney_cr_b),
       IF(weight_own_kg BETWEEN 30 AND 300, kidney_o, kidney_cr_o)) AS k_w,
    IF(wt_ok, GREATEST(kidney_cr_b,
         CASE WHEN uo_gt24_or_anuria_flag = 1 AND uo_n_24h >= 2 THEN 3
              WHEN uo_gt12_flag = 1 AND uo_n_12h >= 2 THEN 2
              WHEN uo_6_12_flag = 1 AND uo_n_6h >= 2 THEN 1 ELSE 0 END),
       kidney_cr_o) AS k_1v,
    IF(hr = t_cut, cv_cut, cv_full) AS cv_hz,
    IF(hr = t_cut, p_cv_cut, p_cv_full) AS p_cv_hz
  FROM b
),
w24 AS (
  SELECT b2.*,
    MAX(brain_std) OVER w AS brain_std_24h, MAX(brain_e) OVER w AS brain_e_24h,
    MAX(brain_gcs15) OVER w AS brain_gcs15_24h, MAX(cns_alt) OVER w AS cns_alt_24h,
    MAX(cns_orig) OVER w AS cns_orig_24h, MAX(cns_orig_n) OVER w AS cns_orig_n_24h,
    MAX(brain_h) OVER w AS brain_h_24h, MAX(brain_i) OVER w AS brain_i_24h, MAX(brain_gcs15_h) OVER w AS brain_gcs15_h_24h,
    COALESCE(MAX(cv_hz) OVER w, 0) AS cv24,
    COALESCE(MAX(cv_full) OVER w, 0) AS cv_full24,
    GREATEST(COALESCE(kidney_pre24, 0), COALESCE(MAX(k_std) OVER w, 0)) AS k24,
    GREATEST(COALESCE(kidney_pre24, 0), COALESCE(MAX(k_cr) OVER w, 0)) AS kcr24,
    GREATEST(COALESCE(kidney_pre24, 0), COALESCE(MAX(k_uo22) OVER w, 0)) AS kuo22_24,
    GREATEST(COALESCE(kidney_pre24, 0), COALESCE(MAX(k_g) OVER w, 0)) AS kg24,
    GREATEST(COALESCE(kidney_pre24, 0), COALESCE(MAX(k_s) OVER w, 0)) AS ks24,
    GREATEST(COALESCE(kidney_pre24, 0), COALESCE(MAX(k_w) OVER w, 0)) AS kw24,
    GREATEST(COALESCE(kidney_pre24, 0), COALESCE(MAX(k_1v) OVER w, 0)) AS k1v24,
    GREATEST(COALESCE(kidney_pre24, 0), COALESCE(MAX(kidney_b) OVER w, 0)) AS kb24,
    MAX(p_resp) OVER w + MAX(p_coag) OVER w + MAX(p_liver) OVER w + MAX(p_renal) OVER w + MAX(p_cv_hz) OVER w AS s1p_nocns24
  FROM b2
  WINDOW w AS (PARTITION BY stay_id, hz, hz_key ORDER BY hr ROWS BETWEEN 23 PRECEDING AND CURRENT ROW)
)
SELECT stay_id, hz, hz_key, hr, endtime, t_cut, old_rule, evaluate, is_base,
  s2_rhl_24hours + cv24 + k24 + brain_std_24h AS s2_std,
  s2_rhl_24hours + cv24 + k24 + brain_e_24h AS s2_e,
  s2_rhl_24hours + cv24 + k24 + brain_gcs15_24h AS s2_gcs15,
  s2_rhl_24hours + cv24 + kcr24 + brain_std_24h AS s2_kcr,
  s2_rhl_24hours + cv24 + kuo22_24 + brain_std_24h AS s2_kuo22,
  s2_rhl_24hours + cv24 + kg24 + brain_std_24h AS s2_g,
  s2_rhl_24hours + cv24 + ks24 + brain_std_24h AS s2_s,
  s2_rhl_24hours + cv24 + kw24 + brain_std_24h AS s2_w,
  s2_rhl_24hours + cv24 + k24 + brain_h_24h AS s2_h,
  s2_rhl_24hours + cv24 + k24 + brain_gcs15_h_24h AS s2_gcs15_h,
  s2_rhl_24hours + cv24 + k24 + brain_i_24h AS s2_i,
  s2_rhl_24hours + cv24 + k1v24 + brain_std_24h AS s2_1v,
  s2_nobrain_24hours + brain_std_24h AS chk_s2_old,
  s2_rhl_24hours + cv_full24 + kb24 + brain_std_24h AS chk_s2_new,
  sofa_24hours AS s1_std,
  s1p_nocns24 + cns_alt_24h AS s1_pipe,
  s1p_nocns24 + cns_orig_24h AS s1_orig,
  s1p_nocns24 + cns_orig_n_24h AS s1_orig_n,
  s1_nocns_24hours + cns_alt_24h AS s1_gcsalt,
  sofa_kcr_24hours AS s1_kcr,
  s2_resp24, s2_hemo24, s2_liver24, cv24 AS s2_cv24, brain_std_24h AS s2_brain24, k24 AS s2_kidney24, kcr24 AS s2_kidney_cr24,
  s1_resp24, s1_coag24, s1_liver24, s1_cv24, s1_cns24, cns_alt_24h AS s1_cns_alt24, cns_orig_24h AS s1_cns_orig24, s1_renal24, s1_renal_cr24
FROM w24
WHERE evaluate OR is_base;
