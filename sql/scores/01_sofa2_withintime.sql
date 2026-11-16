-- Hourly SOFA-2 inputs and component scores for every ICU hour (one row per stay_id, hr).
-- Built from the mimic-code concept tables in @UPSTREAM@. Look-ahead is limited to what is known by each hour:
--   weight: the weight at hour h is the last weight charted by the end of hour h; before the first charted weight the
--     first admission weight is back-applied, and its charting hour (wt_adm_hr) is kept so the day-1 score can use it
--     only when it was charted by hour 23;
--   SpO2/FiO2: FiO2 is matched to SpO2 only if charted at or before the SpO2 and at most 60 minutes earlier;
--   urine output: besides the primary rule (volume charted in the trailing window / window length / weight), the
--     coverage rule of variant g requires at least two charted volumes and charting intervals covering at least
--     5 of 6, 10 of 12 or 22 of 24 hours;
--   ties in the weight and FiO2 matches are broken explicitly so the table is deterministic.
CREATE OR REPLACE TABLE `@BUILD@.sofa2_withintime` AS
WITH co AS (
    SELECT ih.stay_id, ie.hadm_id,ie.intime
        , hr
        , TIMESTAMP_SUB(ih.endtime, INTERVAL '1' HOUR) AS starttime
        , ih.endtime
    FROM `@UPSTREAM@.icustay_hourly` ih
    INNER JOIN `physionet-data.mimiciv_3_1_icu.icustays` ie
        ON ih.stay_id = ie.stay_id
)
-- ========= SpO2/FiO2 fallback
, spo2_obs AS (
  SELECT
    stay_id,
    charttime,
    spo2
  FROM `@UPSTREAM@.vitalsign`
  WHERE charttime IS NOT NULL
    AND spo2 IS NOT NULL
    AND spo2 < 98   -- SOFA-2 rule: only use S/F when SpO2 < 98%
)

, fio2_obs AS (
  SELECT
    stay_id,
    charttime,
    -- normalize FiO2 to a fraction (0.21–1.0)
    CASE
      WHEN fio2 IS NULL THEN NULL
      WHEN fio2 > 1.0 THEN fio2/100.0
      ELSE fio2
    END AS fio2_frac
  FROM `@UPSTREAM@.ventilator_setting`
  WHERE charttime IS NOT NULL
    AND fio2 IS NOT NULL
)

, sf_inst AS (
  -- pair SpO2 with an FiO2 charted at or before it, at most 60 minutes earlier
  SELECT
    s.stay_id,
    s.charttime AS spo2_time,
    f.charttime AS fio2_time,
    s.spo2,
    f.fio2_frac,
    SAFE_DIVIDE(s.spo2, f.fio2_frac) AS spo2fio2ratio
  FROM spo2_obs s
  JOIN fio2_obs f
    ON s.stay_id = f.stay_id
   AND f.charttime <= s.charttime
   AND f.charttime >= DATETIME_SUB(s.charttime, INTERVAL 60 MINUTE)
  WHERE f.fio2_frac BETWEEN 0.21 AND 1.00 --feasible values
)

, sf AS (
  -- hourly min S/F ratio (lower is worse)
  SELECT
    co.stay_id,
    co.hr,
    MIN(sf_inst.spo2fio2ratio) AS spo2fio2ratio_min
  FROM co
  LEFT JOIN sf_inst
    ON co.stay_id = sf_inst.stay_id
   AND co.starttime <= sf_inst.spo2_time
   AND co.endtime > sf_inst.spo2_time
  GROUP BY co.stay_id, co.hr
)
-- ========= Respiratory support flags
, vent_support AS (
    SELECT
        vd.stay_id
        , vd.starttime
        , vd.endtime
        ,CASE
            WHEN vd.ventilation_status IN ('InvasiveVent', 'NonInvasiveVent', 'HFNC', 'Tracheostomy') THEN 1
            ELSE 0
          END AS advanced_vent -- THERE IS ALSO SupplementalOxygen IN THE VALUES WHICH SEEMS NOT MECH VENT
    FROM `@UPSTREAM@.ventilation` vd
)
-- advanced respiratory support in the hour (any overlapping episode), for the SpO2/FiO2 branch
, adv_hourly AS (
  SELECT co.stay_id, co.hr, MAX(vs.advanced_vent) AS advanced_vent_hour
  FROM co
  JOIN vent_support vs
    ON co.stay_id = vs.stay_id AND vs.starttime < co.endtime AND vs.endtime > co.starttime
  GROUP BY co.stay_id, co.hr
)
-- =============== ECMO support

, vv_ecmo_events AS (
  SELECT
    ce.stay_id,
    ce.charttime
  FROM `physionet-data.mimiciv_3_1_icu.chartevents` ce
  WHERE ce.charttime IS NOT NULL
    AND ce.itemid = 229268              -- Circuit Configuration (ECMO)
    AND LOWER(COALESCE(ce.value,'')) IN ('vv','va','vav')
)

, ecmo_hourly AS (
  SELECT
    co.stay_id,
    co.hr,
    1 AS ecmo_any
  FROM co
  JOIN vv_ecmo_events e
    ON co.stay_id = e.stay_id
   AND e.charttime >= co.starttime
   AND e.charttime <  co.endtime
  GROUP BY co.stay_id, co.hr
)
-- ========= PaO2/FiO2 extraction
, pafi AS (
    SELECT ie.stay_id
        , bg.charttime
        , bg.pao2fio2ratio
        , MAX(CASE WHEN vs.advanced_vent = 1 THEN 1 ELSE 0 END) AS advanced_vent_now
    FROM `physionet-data.mimiciv_3_1_icu.icustays` ie
    INNER JOIN `@UPSTREAM@.bg` bg --blood gas table
        ON ie.hadm_id = bg.hadm_id
    LEFT JOIN vent_support vs
        ON ie.stay_id = vs.stay_id
        AND bg.charttime >= vs.starttime
        AND bg.charttime < vs.endtime
    WHERE bg.specimen = 'ART.' -- arterial blood sample
    GROUP BY ie.stay_id, bg.charttime, bg.pao2fio2ratio
)

, pf AS ( -- defining flags
    SELECT co.stay_id, co.hr
        , MIN(pafi.pao2fio2ratio) AS pao2fio2ratio_min
        , MAX(pafi.advanced_vent_now) AS advanced_vent_any
    FROM co
    LEFT JOIN pafi
        ON co.stay_id = pafi.stay_id
        AND pafi.charttime >= co.starttime
        AND pafi.charttime <  co.endtime
    GROUP BY co.stay_id, co.hr
)
-- ========= Vitals /labs
, vs AS (
    SELECT co.stay_id, co.hr
        , MIN(vs.mbp) AS meanbp_min
    FROM co
    LEFT JOIN `@UPSTREAM@.vitalsign` vs
        ON co.stay_id = vs.stay_id
        AND co.starttime <= vs.charttime
        AND co.endtime > vs.charttime
    GROUP BY co.stay_id, co.hr
)
-- ========= Sedation timing for SOFA-2 brain (footnote c)
-- list of meds to be checked
, sedative_itemids AS (
  SELECT itemid
  FROM `physionet-data.mimiciv_3_1_icu.d_items`
  WHERE REGEXP_CONTAINS(LOWER(label), r'\b(propofol|midazolam|lorazepam|dexmedetomidine)\b') --exact words / list should be checked with clinicians
)

, sedation_intervals AS ( -- taking all sedatives in the stay
  SELECT
    ie.stay_id,
    ie.starttime,
    ie.endtime
  FROM `physionet-data.mimiciv_3_1_icu.inputevents` ie
  WHERE ie.stay_id IS NOT NULL
    AND ie.starttime IS NOT NULL
    AND ie.endtime IS NOT NULL
    AND ie.itemid IN (SELECT itemid FROM sedative_itemids)
)

, sed_start AS ( -- start time of sedation
  SELECT
    stay_id,
    MIN(starttime) AS sed_starttime
  FROM sedation_intervals
  GROUP BY stay_id
)

, sedation_hourly AS ( --binary flag for sedated hour
  SELECT
    co.stay_id,
    co.hr,
    CASE WHEN COUNTIF(s.starttime < co.endtime AND s.endtime > co.starttime) > 0 THEN 1 ELSE 0 END AS sedated_hour
  FROM co
  LEFT JOIN sedation_intervals s
    ON co.stay_id = s.stay_id
  GROUP BY co.stay_id, co.hr
)
, gcs AS (
  -- Raw component events
  WITH gcs_components AS (
    SELECT
      ce.stay_id,
      ce.charttime,
      CASE WHEN ce.itemid = 223901 THEN SAFE_CAST(ce.valuenum AS INT64) END AS motor, --motor
      CASE WHEN ce.itemid = 223900 THEN --verbal
        CASE
          WHEN LOWER(ce.value) = 'no response-ett' THEN NULL  -- cannot evaluate verbal
          ELSE SAFE_CAST(ce.valuenum AS INT64)
        END
      END AS verbal,
      CASE WHEN ce.itemid = 220739 THEN SAFE_CAST(ce.valuenum AS INT64) END AS eyes, --eye opening
      CASE WHEN ce.itemid = 223900 AND LOWER(ce.value) = 'no response-ett' THEN 1 ELSE 0 END AS ett_flag --intubated
    FROM `physionet-data.mimiciv_3_1_icu.chartevents` ce
    WHERE ce.stay_id IS NOT NULL
      AND ce.charttime IS NOT NULL
      AND ce.itemid IN (223901, 223900, 220739)
  ),

  -- Hourly observed total + best motor in the hour
  gcs_hourly_obs AS (
    SELECT
      co.stay_id,
      co.hr,
      MIN(
        CASE
          WHEN gc.motor IS NOT NULL AND gc.verbal IS NOT NULL AND gc.eyes IS NOT NULL
          THEN (gc.motor + gc.verbal + gc.eyes)
        END
      ) AS total_gcs_obs,
      MAX(gc.motor) AS motor_best, -- fotnote d
      MAX(gc.ett_flag) AS gcs_unable_any --footnote d
    FROM co
    LEFT JOIN gcs_components gc
      ON co.stay_id = gc.stay_id
     AND gc.charttime >= co.starttime
     AND gc.charttime <  co.endtime
    GROUP BY co.stay_id, co.hr
  ),

  -- Last recorded TOTAL GCS before sedation starttime (per stay)
  pre_sed AS (
    SELECT
      ss.stay_id,
      ARRAY_AGG(
        (gc.motor + gc.verbal + gc.eyes)
        ORDER BY gc.charttime DESC
        LIMIT 1
      )[SAFE_OFFSET(0)] AS pre_sed_total_gcs --incase the array is empty
    FROM sed_start ss
    JOIN gcs_components gc
      ON ss.stay_id = gc.stay_id
     AND gc.charttime < ss.sed_starttime
    WHERE gc.motor IS NOT NULL AND gc.verbal IS NOT NULL AND gc.eyes IS NOT NULL
    GROUP BY ss.stay_id
  )

  SELECT
    o.stay_id,
    o.hr,
    COALESCE(sh.sedated_hour, 0) AS sedated_hour,


    CASE
      WHEN COALESCE(sh.sedated_hour,0) = 1 THEN ps.pre_sed_total_gcs
      ELSE o.total_gcs_obs
    END AS gcs_min,

    -- motor fallback (footnote d): best achieved motor in the hour
    o.motor_best AS gcs_motor_min,

    o.gcs_unable_any
  FROM gcs_hourly_obs o
  LEFT JOIN sedation_hourly sh
    ON o.stay_id = sh.stay_id AND o.hr = sh.hr
  LEFT JOIN pre_sed ps
    ON o.stay_id = ps.stay_id
)

, bili AS (
    SELECT co.stay_id, co.hr
        , MAX(enz.bilirubin_total) AS bilirubin_max
    FROM co
    LEFT JOIN `@UPSTREAM@.enzyme` enz
        ON co.hadm_id = enz.hadm_id
        AND co.starttime <= enz.charttime
        AND co.endtime > enz.charttime
    GROUP BY co.stay_id, co.hr

)

, cr AS (
    SELECT co.stay_id, co.hr
        , MAX(chem.creatinine) AS creatinine_max
    FROM co
    LEFT JOIN `@UPSTREAM@.chemistry` chem
        ON co.hadm_id = chem.hadm_id
        AND co.starttime <= chem.charttime
        AND co.endtime > chem.charttime
    GROUP BY co.stay_id, co.hr
)
-- ========= RRT eligibility labs

, k_hco3 AS (
  SELECT
    co.stay_id,
    co.hr,
    MAX(chem.potassium)   AS potassium_max,
    MIN(chem.bicarbonate) AS bicarbonate_min
  FROM co
  LEFT JOIN `@UPSTREAM@.chemistry` chem
    ON co.hadm_id = chem.hadm_id
   AND co.starttime <= chem.charttime
   AND co.endtime > chem.charttime
  GROUP BY co.stay_id, co.hr
)

, ph_lab AS (
  SELECT
    co.stay_id,
    co.hr,
    MIN(bg.ph) AS ph_min
  FROM co
  LEFT JOIN `@UPSTREAM@.bg` bg
    ON co.hadm_id = bg.hadm_id
   AND co.starttime <= bg.charttime
   AND co.endtime > bg.charttime
  WHERE bg.ph IS NOT NULL
  GROUP BY co.stay_id, co.hr
)
, plt AS (
    SELECT co.stay_id, co.hr
        , MIN(cbc.platelet) AS platelet_min
    FROM co
    LEFT JOIN `@UPSTREAM@.complete_blood_count` cbc
        ON co.hadm_id = cbc.hadm_id
        AND cbc.charttime >= co.starttime
        AND cbc.charttime <  co.endtime
    GROUP BY co.stay_id, co.hr
)

-- ========= Urine output

, wt_charted AS (
  SELECT ce.stay_id, ce.itemid, ce.charttime, ce.valuenum AS weight
  FROM `physionet-data.mimiciv_3_1_icu.chartevents` ce
  WHERE ce.itemid IN (226512, 224639) -- admission weight, daily weight (kg)
    AND ce.valuenum > 0
)

, wt_adm AS (  -- first charted admission weight (back-applied only when charted by the horizon cut)
  SELECT stay_id, ARRAY_AGG(STRUCT(charttime, weight) ORDER BY charttime, weight LIMIT 1)[OFFSET(0)] AS a
  FROM wt_charted WHERE itemid = 226512
  GROUP BY stay_id
)

, wt_hour AS (
  SELECT
    co.stay_id,
    co.hr,
    co.endtime,
    (ARRAY_AGG(w.weight IGNORE NULLS ORDER BY w.charttime DESC, w.weight DESC LIMIT 1))[SAFE_OFFSET(0)] AS weight_hour_kg
  FROM co
  LEFT JOIN wt_charted w
    ON co.stay_id = w.stay_id
   AND w.charttime <= co.endtime
  GROUP BY co.stay_id, co.hr, co.endtime
)

, wt AS (
  SELECT
    wh.stay_id,
    wh.hr,
    wh.weight_hour_kg AS weight_own_kg,
    COALESCE(wh.weight_hour_kg, wa.a.weight) AS weight_ba_kg,
    IF(wh.weight_hour_kg IS NULL AND wa.a.weight IS NOT NULL, 1, 0) AS weight_backapplied,
    -- grid hour whose (starttime, endtime] holds the admission-weight chart time
    wh.hr + CAST(CEIL(DATETIME_DIFF(wa.a.charttime, wh.endtime, SECOND) / 3600.0) AS INT64) AS wt_adm_hr
  FROM wt_hour wh
  LEFT JOIN wt_adm wa
    ON wh.stay_id = wa.stay_id
)

, uo_tm AS (  -- minutes since the previous urine chart (first chart: since intime_hr), as mimic-code urine_output_rate
  SELECT u.stay_id, u.charttime, u.urineoutput,
    DATETIME_DIFF(u.charttime,
      COALESCE(LAG(u.charttime) OVER (PARTITION BY u.stay_id ORDER BY u.charttime), t.intime_hr), MINUTE) AS tm_min
  FROM `@UPSTREAM@.urine_output` u
  JOIN `@UPSTREAM@.icustay_times` t USING (stay_id)
)

,uo AS (
  SELECT
    co.stay_id,
    co.hr,
    wt.weight_own_kg,
    wt.weight_ba_kg,
    wt.weight_backapplied,
    wt.wt_adm_hr,

    SUM(CASE WHEN u.charttime >= DATETIME_SUB(co.endtime, INTERVAL 6 HOUR)
             AND u.charttime < co.endtime
             THEN u.urineoutput ELSE 0 END) AS uo_6h_ml,

    SUM(CASE WHEN u.charttime >= DATETIME_SUB(co.endtime, INTERVAL 12 HOUR)
             AND u.charttime < co.endtime
             THEN u.urineoutput ELSE 0 END) AS uo_12h_ml,

    SUM(CASE WHEN u.charttime >= DATETIME_SUB(co.endtime, INTERVAL 24 HOUR)
             AND u.charttime < co.endtime
             THEN u.urineoutput ELSE 0 END) AS uo_24h_ml,

    MAX(CASE WHEN u.charttime >= DATETIME_SUB(co.endtime, INTERVAL 6 HOUR)
             AND u.charttime < co.endtime THEN 1 ELSE 0 END) AS has_uo_meas_6h,

    MAX(CASE WHEN u.charttime >= DATETIME_SUB(co.endtime, INTERVAL 12 HOUR)
             AND u.charttime < co.endtime THEN 1 ELSE 0 END) AS has_uo_meas_12h,

    MAX(CASE WHEN u.charttime >= DATETIME_SUB(co.endtime, INTERVAL 24 HOUR)
             AND u.charttime < co.endtime THEN 1 ELSE 0 END) AS has_uo_meas_24h,

    -- variant g: charts and charting-interval coverage (minutes) inside each window
    COUNTIF(u.charttime >= DATETIME_SUB(co.endtime, INTERVAL 6 HOUR) AND u.charttime < co.endtime) AS uo_n_6h,
    COUNTIF(u.charttime >= DATETIME_SUB(co.endtime, INTERVAL 12 HOUR) AND u.charttime < co.endtime) AS uo_n_12h,
    COUNTIF(u.charttime >= DATETIME_SUB(co.endtime, INTERVAL 24 HOUR) AND u.charttime < co.endtime) AS uo_n_24h,
    SUM(CASE WHEN u.charttime >= DATETIME_SUB(co.endtime, INTERVAL 6 HOUR)
             AND u.charttime < co.endtime THEN u.tm_min ELSE 0 END) AS uo_tm_6h,
    SUM(CASE WHEN u.charttime >= DATETIME_SUB(co.endtime, INTERVAL 12 HOUR)
             AND u.charttime < co.endtime THEN u.tm_min ELSE 0 END) AS uo_tm_12h,
    SUM(CASE WHEN u.charttime >= DATETIME_SUB(co.endtime, INTERVAL 24 HOUR)
             AND u.charttime < co.endtime THEN u.tm_min ELSE 0 END) AS uo_tm_24h

  FROM co
  LEFT JOIN uo_tm u
    ON co.stay_id = u.stay_id
   AND u.charttime >= DATETIME_SUB(co.endtime, INTERVAL 24 HOUR)
   AND u.charttime < co.endtime
  LEFT JOIN wt
    ON co.stay_id = wt.stay_id
   AND co.hr = wt.hr
  GROUP BY co.stay_id, co.hr, wt.weight_own_kg, wt.weight_ba_kg, wt.weight_backapplied, wt.wt_adm_hr
),

uo_cov AS (
  SELECT uo.*,
    (uo_n_6h >= 2 AND uo_tm_6h >= 5 * 60) AS cov_6h,
    (uo_n_12h >= 2 AND uo_tm_12h >= 10 * 60) AS cov_12h,
    (uo_n_24h >= 2 AND uo_tm_24h >= 22 * 60) AS cov_24h,
    -- primary rule: the charting intervals of the charts in the window cover the whole window (mimic-code urine_output_rate)
    (uo_tm_6h >= 6 * 60) AS iv_6h,
    (uo_tm_12h >= 12 * 60) AS iv_12h,
    (uo_tm_24h >= 24 * 60) AS iv_24h
  FROM uo
),

uo_flags AS (
  SELECT
    stay_id,
    hr,
    uo_6h_ml,
    uo_24h_ml,
    iv_24h,
    weight_own_kg,
    weight_ba_kg,
    weight_backapplied,
    wt_adm_hr,
    uo_n_6h, uo_n_12h, uo_n_24h, cov_6h, cov_12h, cov_24h,
    -- primary rule: rate = volume charted in the window / summed charting intervals of those charts / weight, scored
    -- only when the intervals cover the window; anuria needs at least two charts
    CASE WHEN weight_own_kg IS NULL OR weight_own_kg <= 0 THEN NULL
      WHEN NOT (iv_6h) THEN NULL
      WHEN hr < 5 THEN NULL
      WHEN SAFE_DIVIDE(uo_6h_ml, weight_own_kg * uo_tm_6h / 60.0) < 0.5
       AND (hr < 11 OR NOT (iv_12h) OR SAFE_DIVIDE(uo_12h_ml, weight_own_kg * uo_tm_12h / 60.0) >= 0.5) THEN 1
      ELSE 0 END AS uo_6_12_o,
    CASE WHEN weight_own_kg IS NULL OR weight_own_kg <= 0 THEN NULL
      WHEN NOT (iv_12h) THEN NULL
      WHEN hr < 11 THEN NULL
      WHEN SAFE_DIVIDE(uo_12h_ml, weight_own_kg * uo_tm_12h / 60.0) < 0.5 THEN 1
      ELSE 0 END AS uo_gt12_o,
    CASE WHEN weight_own_kg IS NULL OR weight_own_kg <= 0 THEN NULL
      WHEN hr >= 11 AND (iv_12h) AND uo_n_12h >= 2 AND uo_12h_ml = 0 THEN 1
      WHEN hr >= 23 AND (iv_24h) AND SAFE_DIVIDE(uo_24h_ml, weight_own_kg * uo_tm_24h / 60.0) < 0.3 THEN 1
      ELSE 0 END AS uo_gt24_o,
    (weight_own_kg IS NOT NULL AND weight_own_kg > 0 AND hr >= 5 AND (iv_6h)
       AND SAFE_DIVIDE(uo_6h_ml, weight_own_kg * uo_tm_6h / 60.0) < 0.3) AS olig6_o,
    CASE WHEN weight_ba_kg IS NULL OR weight_ba_kg <= 0 THEN NULL
      WHEN NOT (iv_6h) THEN NULL
      WHEN hr < 5 THEN NULL
      WHEN SAFE_DIVIDE(uo_6h_ml, weight_ba_kg * uo_tm_6h / 60.0) < 0.5
       AND (hr < 11 OR NOT (iv_12h) OR SAFE_DIVIDE(uo_12h_ml, weight_ba_kg * uo_tm_12h / 60.0) >= 0.5) THEN 1
      ELSE 0 END AS uo_6_12_b,
    CASE WHEN weight_ba_kg IS NULL OR weight_ba_kg <= 0 THEN NULL
      WHEN NOT (iv_12h) THEN NULL
      WHEN hr < 11 THEN NULL
      WHEN SAFE_DIVIDE(uo_12h_ml, weight_ba_kg * uo_tm_12h / 60.0) < 0.5 THEN 1
      ELSE 0 END AS uo_gt12_b,
    CASE WHEN weight_ba_kg IS NULL OR weight_ba_kg <= 0 THEN NULL
      WHEN hr >= 11 AND (iv_12h) AND uo_n_12h >= 2 AND uo_12h_ml = 0 THEN 1
      WHEN hr >= 23 AND (iv_24h) AND SAFE_DIVIDE(uo_24h_ml, weight_ba_kg * uo_tm_24h / 60.0) < 0.3 THEN 1
      ELSE 0 END AS uo_gt24_b,
    (weight_ba_kg IS NOT NULL AND weight_ba_kg > 0 AND hr >= 5 AND (iv_6h)
       AND SAFE_DIVIDE(uo_6h_ml, weight_ba_kg * uo_tm_6h / 60.0) < 0.3) AS olig6_b,
    -- variant s: single-chart rule (volume charted in the trailing window / window length / weight)
    CASE WHEN weight_own_kg IS NULL OR weight_own_kg <= 0 THEN NULL
      WHEN NOT (has_uo_meas_6h = 1) THEN NULL
      WHEN hr < 5 THEN NULL
      WHEN SAFE_DIVIDE(uo_6h_ml, weight_own_kg * 6.0) < 0.5
       AND (hr < 11 OR NOT (has_uo_meas_12h = 1) OR SAFE_DIVIDE(uo_12h_ml, weight_own_kg * 12.0) >= 0.5) THEN 1
      ELSE 0 END AS uo_6_12_so,
    CASE WHEN weight_own_kg IS NULL OR weight_own_kg <= 0 THEN NULL
      WHEN NOT (has_uo_meas_12h = 1) THEN NULL
      WHEN hr < 11 THEN NULL
      WHEN SAFE_DIVIDE(uo_12h_ml, weight_own_kg * 12.0) < 0.5 THEN 1
      ELSE 0 END AS uo_gt12_so,
    CASE WHEN weight_own_kg IS NULL OR weight_own_kg <= 0 THEN NULL
      WHEN hr >= 11 AND (has_uo_meas_12h = 1) AND uo_12h_ml = 0 THEN 1
      WHEN hr >= 23 AND (has_uo_meas_24h = 1) AND SAFE_DIVIDE(uo_24h_ml, weight_own_kg * 24.0) < 0.3 THEN 1
      ELSE 0 END AS uo_gt24_so,
    (weight_own_kg IS NOT NULL AND weight_own_kg > 0 AND hr >= 5 
       AND SAFE_DIVIDE(uo_6h_ml, weight_own_kg * 6.0) < 0.3) AS olig6_so,
    CASE WHEN weight_ba_kg IS NULL OR weight_ba_kg <= 0 THEN NULL
      WHEN NOT (has_uo_meas_6h = 1) THEN NULL
      WHEN hr < 5 THEN NULL
      WHEN SAFE_DIVIDE(uo_6h_ml, weight_ba_kg * 6.0) < 0.5
       AND (hr < 11 OR NOT (has_uo_meas_12h = 1) OR SAFE_DIVIDE(uo_12h_ml, weight_ba_kg * 12.0) >= 0.5) THEN 1
      ELSE 0 END AS uo_6_12_sb,
    CASE WHEN weight_ba_kg IS NULL OR weight_ba_kg <= 0 THEN NULL
      WHEN NOT (has_uo_meas_12h = 1) THEN NULL
      WHEN hr < 11 THEN NULL
      WHEN SAFE_DIVIDE(uo_12h_ml, weight_ba_kg * 12.0) < 0.5 THEN 1
      ELSE 0 END AS uo_gt12_sb,
    CASE WHEN weight_ba_kg IS NULL OR weight_ba_kg <= 0 THEN NULL
      WHEN hr >= 11 AND (has_uo_meas_12h = 1) AND uo_12h_ml = 0 THEN 1
      WHEN hr >= 23 AND (has_uo_meas_24h = 1) AND SAFE_DIVIDE(uo_24h_ml, weight_ba_kg * 24.0) < 0.3 THEN 1
      ELSE 0 END AS uo_gt24_sb,
    (weight_ba_kg IS NOT NULL AND weight_ba_kg > 0 AND hr >= 5 
       AND SAFE_DIVIDE(uo_6h_ml, weight_ba_kg * 6.0) < 0.3) AS olig6_sb,
    CASE WHEN weight_own_kg IS NULL OR weight_own_kg <= 0 THEN NULL
      WHEN NOT (cov_6h) THEN NULL
      WHEN hr < 5 THEN NULL
      WHEN SAFE_DIVIDE(uo_6h_ml, weight_own_kg * 6.0) < 0.5
       AND (hr < 11 OR NOT (cov_12h) OR SAFE_DIVIDE(uo_12h_ml, weight_own_kg * 12.0) >= 0.5) THEN 1
      ELSE 0 END AS uo_6_12_go,
    CASE WHEN weight_own_kg IS NULL OR weight_own_kg <= 0 THEN NULL
      WHEN NOT (cov_12h) THEN NULL
      WHEN hr < 11 THEN NULL
      WHEN SAFE_DIVIDE(uo_12h_ml, weight_own_kg * 12.0) < 0.5 THEN 1
      ELSE 0 END AS uo_gt12_go,
    CASE WHEN weight_own_kg IS NULL OR weight_own_kg <= 0 THEN NULL
      WHEN hr >= 11 AND (cov_12h) AND uo_12h_ml = 0 THEN 1
      WHEN hr >= 23 AND (cov_24h) AND SAFE_DIVIDE(uo_24h_ml, weight_own_kg * 24.0) < 0.3 THEN 1
      ELSE 0 END AS uo_gt24_go,
    (weight_own_kg IS NOT NULL AND weight_own_kg > 0 AND hr >= 5 AND (cov_6h)
       AND SAFE_DIVIDE(uo_6h_ml, weight_own_kg * 6.0) < 0.3) AS olig6_go,
    CASE WHEN weight_ba_kg IS NULL OR weight_ba_kg <= 0 THEN NULL
      WHEN NOT (cov_6h) THEN NULL
      WHEN hr < 5 THEN NULL
      WHEN SAFE_DIVIDE(uo_6h_ml, weight_ba_kg * 6.0) < 0.5
       AND (hr < 11 OR NOT (cov_12h) OR SAFE_DIVIDE(uo_12h_ml, weight_ba_kg * 12.0) >= 0.5) THEN 1
      ELSE 0 END AS uo_6_12_gb,
    CASE WHEN weight_ba_kg IS NULL OR weight_ba_kg <= 0 THEN NULL
      WHEN NOT (cov_12h) THEN NULL
      WHEN hr < 11 THEN NULL
      WHEN SAFE_DIVIDE(uo_12h_ml, weight_ba_kg * 12.0) < 0.5 THEN 1
      ELSE 0 END AS uo_gt12_gb,
    CASE WHEN weight_ba_kg IS NULL OR weight_ba_kg <= 0 THEN NULL
      WHEN hr >= 11 AND (cov_12h) AND uo_12h_ml = 0 THEN 1
      WHEN hr >= 23 AND (cov_24h) AND SAFE_DIVIDE(uo_24h_ml, weight_ba_kg * 24.0) < 0.3 THEN 1
      ELSE 0 END AS uo_gt24_gb,
    (weight_ba_kg IS NOT NULL AND weight_ba_kg > 0 AND hr >= 5 AND (cov_6h)
       AND SAFE_DIVIDE(uo_6h_ml, weight_ba_kg * 6.0) < 0.3) AS olig6_gb
  FROM uo_cov
)

-- ========= RRT flag
, rrt AS (
  SELECT
      stay_id
    , charttime
    , CASE
        WHEN dialysis_active = 1 THEN 1
        WHEN dialysis_present = 1 THEN 1
      END AS on_rrt
    , dialysis_type
  FROM `@UPSTREAM@.rrt`
  WHERE charttime IS NOT NULL
)

-- =============== Delirium treatment from EMAR
-- the list of meds should be consulted with clinicians
, delirium_meds AS (
  SELECT 'haloperidol' AS token UNION ALL
  SELECT 'quetiapine' UNION ALL
  SELECT 'olanzapine' UNION ALL
  SELECT 'risperidone' UNION ALL
  SELECT 'ziprasidone' UNION ALL
  SELECT 'chlorpromazine'
)
-- variant e only: dexmedetomidine counted as delirium treatment
, delirium_meds_dex AS (
  SELECT 'dexmedetomidine' AS token UNION ALL
  SELECT 'precedex'
)

, emar_filtered AS ( -- can be improved  with better filtering
  SELECT
    i.stay_id,
    e.charttime,
    LOWER(COALESCE(e.medication, '')) AS med_text,
    LOWER(COALESCE(e.event_txt, '')) AS event_txt
  FROM `physionet-data.mimiciv_3_1_icu.icustays` i
  INNER JOIN `physionet-data.mimiciv_3_1_hosp.emar` e
    ON i.hadm_id = e.hadm_id
  WHERE e.charttime BETWEEN i.intime AND i.outtime
    AND (
      LOWER(COALESCE(e.event_txt,'')) LIKE '%admin%'
      OR LOWER(COALESCE(e.event_txt,'')) LIKE '%started%'
      OR LOWER(COALESCE(e.event_txt,'')) LIKE '%new bag%'
      OR LOWER(COALESCE(e.event_txt,'')) LIKE '%confirmed%'
    )
    AND LOWER(COALESCE(e.event_txt,'')) NOT LIKE '%not given%'
    AND LOWER(COALESCE(e.event_txt,'')) NOT LIKE '%not started%'
    AND LOWER(COALESCE(e.event_txt,'')) NOT LIKE '%not confirmed%'
    AND LOWER(COALESCE(e.event_txt,'')) NOT LIKE '%held%'
    AND LOWER(COALESCE(e.event_txt,'')) NOT LIKE '%hold%'
    AND LOWER(COALESCE(e.event_txt,'')) NOT LIKE '%refused%'
)

, delirium_tx AS (
  SELECT
    stay_id,
    charttime,
    IF(EXISTS (SELECT 1 FROM delirium_meds d
               WHERE LOWER(m.med_text) LIKE CONCAT('%', LOWER(d.token), '%')), 1, 0) AS delirium_treated,
    IF(EXISTS (SELECT 1 FROM delirium_meds_dex d
               WHERE LOWER(m.med_text) LIKE CONCAT('%', LOWER(d.token), '%')), 1, 0) AS dex_treated
  FROM emar_filtered m
  WHERE EXISTS (SELECT 1 FROM delirium_meds d WHERE LOWER(m.med_text) LIKE CONCAT('%', LOWER(d.token), '%'))
     OR EXISTS (SELECT 1 FROM delirium_meds_dex d WHERE LOWER(m.med_text) LIKE CONCAT('%', LOWER(d.token), '%'))
)

-- ========== Other vaso hourly flag
, other_vaso AS (
  SELECT
    stay_id,
    hr,
    MAX(flag) AS other_vaso_flag
  FROM (

    -- Vasopressin
    SELECT
      co.stay_id,
      co.hr,
      MAX(CASE WHEN vp.vaso_rate IS NOT NULL AND vp.vaso_rate > 0 THEN 1 ELSE 0 END) AS flag
    FROM co
    LEFT JOIN `@UPSTREAM@.vasopressin` vp
      ON co.stay_id = vp.stay_id
      AND TIMESTAMP_DIFF(
       LEAST(vp.endtime, co.endtime),
       GREATEST(vp.starttime, co.starttime),
       MINUTE
     ) >= 60
    GROUP BY co.stay_id, co.hr

    UNION ALL

    -- Phenylephrine
    SELECT
      co.stay_id,
      co.hr,
      MAX(CASE WHEN phe.vaso_rate IS NOT NULL AND phe.vaso_rate > 0 THEN 1 ELSE 0 END) AS flag
    FROM co
    LEFT JOIN `@UPSTREAM@.phenylephrine` phe
      ON co.stay_id = phe.stay_id
      AND TIMESTAMP_DIFF(
       LEAST(phe.endtime, co.endtime),
       GREATEST(phe.starttime, co.starttime),
       MINUTE
     ) >= 60
    GROUP BY co.stay_id, co.hr

    UNION ALL

    SELECT
      co.stay_id,
      co.hr,
      MAX(CASE WHEN mil.vaso_rate IS NOT NULL AND mil.vaso_rate > 0 THEN 1 ELSE 0 END) AS flag
    FROM co
    LEFT JOIN `@UPSTREAM@.milrinone` mil
      ON co.stay_id = mil.stay_id
      AND TIMESTAMP_DIFF(
       LEAST(mil.endtime, co.endtime),
       GREATEST(mil.starttime, co.starttime),
       MINUTE
     ) >= 60
    GROUP BY co.stay_id, co.hr

  )
  GROUP BY stay_id, hr
)

-- ========== VA-ECMO detection for cardiovascular mechanical support
, ecmo_circuit_item AS (
  SELECT itemid
  FROM `physionet-data.mimiciv_3_1_icu.d_items`
  WHERE LOWER(label) = 'circuit configuration (ecmo)'
  LIMIT 1
)

, va_ecmo_events AS (
  SELECT
    ce.stay_id,
    ce.charttime
  FROM `physionet-data.mimiciv_3_1_icu.chartevents` ce
  WHERE ce.charttime IS NOT NULL
    AND ce.itemid IN (SELECT itemid FROM ecmo_circuit_item)
    AND LOWER(COALESCE(ce.value,'')) IN ('va', 'vav')
)

, va_ecmo_hourly AS (
  SELECT
    co.stay_id,
    co.hr,
    1 AS va_ecmo_flag
  FROM co
  JOIN va_ecmo_events e
    ON co.stay_id = e.stay_id
   AND co.starttime <= e.charttime
   AND co.endtime   > e.charttime
  GROUP BY co.stay_id, co.hr
)
-- ========== Mechanical circulatory support

, mcs_itemids AS ( --very doubtful about this one as the distinct values
--are alot and I'm not sure about many of them
  SELECT itemid
  FROM `physionet-data.mimiciv_3_1_icu.d_items`
  WHERE
       LOWER(label) LIKE '%iabp%'
    OR LOWER(label) LIKE '%impella%'
    OR LOWER(label) LIKE '%lvad%'
    OR LOWER(label) LIKE '% vad%'
)

, mcs_device_events AS (
  SELECT
    ce.stay_id,
    ce.charttime
  FROM `physionet-data.mimiciv_3_1_icu.chartevents` ce
  WHERE ce.charttime IS NOT NULL
    AND ce.itemid IN (SELECT itemid FROM mcs_itemids)
)

, mcs_hourly AS (
  SELECT
    co.stay_id,
    co.hr,
    1 AS mechanical_support_flag
  FROM co
  LEFT JOIN mcs_device_events d
    ON co.stay_id = d.stay_id
   AND d.charttime >= co.starttime AND d.charttime < co.endtime
  LEFT JOIN va_ecmo_hourly va
    ON co.stay_id = va.stay_id
   AND co.hr = va.hr
  WHERE d.charttime IS NOT NULL OR va.va_ecmo_flag = 1
  GROUP BY co.stay_id, co.hr
)
, vaso AS (
    SELECT
        co.stay_id
        , co.hr
        , MAX(epi.vaso_rate) AS rate_epinephrine_ugkgmin
        , MAX(nor.vaso_rate) AS rate_norepinephrine_ugkgmin
        , MAX(dop.vaso_rate) AS rate_dopamine_mcgkgmin
        , MAX(dob.vaso_rate) AS rate_dobutamine
        , COALESCE(MAX(ov.other_vaso_flag), 0) AS other_vaso_flag
        , COALESCE(MAX(mcs.mechanical_support_flag), 0) AS mechanical_support_flag
    FROM co
    LEFT JOIN `@UPSTREAM@.epinephrine` epi
        ON co.stay_id = epi.stay_id AND TIMESTAMP_DIFF(
       LEAST(epi.endtime, co.endtime),
       GREATEST(epi.starttime, co.starttime),
       MINUTE
     ) >= 60
    LEFT JOIN `@UPSTREAM@.norepinephrine` nor
        ON co.stay_id = nor.stay_id AND TIMESTAMP_DIFF(
       LEAST(nor.endtime, co.endtime),
       GREATEST(nor.starttime, co.starttime),
       MINUTE
     ) >= 60
    LEFT JOIN `@UPSTREAM@.dopamine` dop
        ON co.stay_id = dop.stay_id AND TIMESTAMP_DIFF(
       LEAST(dop.endtime, co.endtime),
       GREATEST(dop.starttime, co.starttime),
       MINUTE
     ) >= 60
    LEFT JOIN `@UPSTREAM@.dobutamine` dob
        ON co.stay_id = dob.stay_id AND TIMESTAMP_DIFF(
       LEAST(dob.endtime, co.endtime),
       GREATEST(dob.starttime, co.starttime),
       MINUTE
     ) >= 60
    LEFT JOIN other_vaso ov ON co.stay_id = ov.stay_id AND co.hr = ov.hr
    LEFT JOIN mcs_hourly mcs ON co.stay_id = mcs.stay_id AND co.hr = mcs.hr
    GROUP BY co.stay_id, co.hr
)

, scorecomp AS (
    SELECT
        co.stay_id, co.hr, co.starttime, co.endtime, co.intime
        , pf.pao2fio2ratio_min
        , pf.advanced_vent_any
        , MAX(adv.advanced_vent_hour) AS advanced_vent_hour
        -- , pf.ecmo_any
        , vs.meanbp_min
        , gcs.gcs_min
        , gcs.gcs_motor_min
        , gcs.gcs_unable_any
        , bili.bilirubin_max
        , cr.creatinine_max
        , plt.platelet_min
        , vaso.rate_epinephrine_ugkgmin
        , vaso.rate_norepinephrine_ugkgmin
        , vaso.rate_dopamine_mcgkgmin
        , vaso.rate_dobutamine
        , vaso.other_vaso_flag
        , vaso.mechanical_support_flag
        , MAX(rrt.on_rrt) AS on_rrt_any
        , COALESCE(MAX(del.delirium_treated), 0) AS delirium_treated_any
        , COALESCE(MAX(del.dex_treated), 0) AS dex_treated_any
        , MAX(uo.uo_6_12_b) AS uo_6_12_flag
        , MAX(uo.uo_gt12_b) AS uo_gt12_flag
        , MAX(uo.uo_gt24_b) AS uo_gt24_or_anuria_flag
        , MAX(uo.uo_6_12_b) AS uo_6_12_b, MAX(uo.uo_gt12_b) AS uo_gt12_b, MAX(uo.uo_gt24_b) AS uo_gt24_b
        , MAX(uo.uo_6_12_o) AS uo_6_12_o, MAX(uo.uo_gt12_o) AS uo_gt12_o, MAX(uo.uo_gt24_o) AS uo_gt24_o
        , MAX(uo.uo_6_12_go) AS uo_6_12_go, MAX(uo.uo_gt12_go) AS uo_gt12_go, MAX(uo.uo_gt24_go) AS uo_gt24_go
        , MAX(uo.uo_6_12_gb) AS uo_6_12_gb, MAX(uo.uo_gt12_gb) AS uo_gt12_gb, MAX(uo.uo_gt24_gb) AS uo_gt24_gb
        , LOGICAL_OR(uo.olig6_o) AS olig6_o, LOGICAL_OR(uo.olig6_b) AS olig6_b
        , LOGICAL_OR(uo.olig6_go) AS olig6_go, LOGICAL_OR(uo.olig6_gb) AS olig6_gb
        , MAX(uo.uo_6_12_so) AS uo_6_12_so, MAX(uo.uo_gt12_so) AS uo_gt12_so, MAX(uo.uo_gt24_so) AS uo_gt24_so
        , MAX(uo.uo_6_12_sb) AS uo_6_12_sb, MAX(uo.uo_gt12_sb) AS uo_gt12_sb, MAX(uo.uo_gt24_sb) AS uo_gt24_sb
        , LOGICAL_OR(uo.olig6_so) AS olig6_so, LOGICAL_OR(uo.olig6_sb) AS olig6_sb
        , MAX(uo.weight_own_kg) AS weight_own_kg, MAX(uo.weight_ba_kg) AS weight_ba_kg, MAX(uo.wt_adm_hr) AS wt_adm_hr
        , MAX(uo.uo_n_6h) AS uo_n_6h, MAX(uo.uo_n_12h) AS uo_n_12h, MAX(uo.uo_n_24h) AS uo_n_24h
        , LOGICAL_OR(uo.cov_6h) AS cov_6h, LOGICAL_OR(uo.cov_12h) AS cov_12h, LOGICAL_OR(uo.cov_24h) AS cov_24h
        , MAX(uo.uo_24h_ml) AS uo_24h_ml, LOGICAL_OR(uo.iv_24h) AS iv_24h
        , (MAX(cr.creatinine_max) > 1.2) AS cr_trig
        , (MAX(kh.potassium_max) >= 6.0 OR (MIN(ph.ph_min) <= 7.20 AND MIN(kh.bicarbonate_min) <= 12.0)) AS rrt_sev
        , COALESCE(MAX(eh.ecmo_any), 0) AS ecmo_any
        , MIN(sf.spo2fio2ratio_min) AS spo2fio2ratio_min
        , MAX(kh.potassium_max) AS f0_ -- unnamed in the upstream query; saved table named it f0_
        , MIN(kh.bicarbonate_min) AS f1_
        , MIN(ph.ph_min) AS f2_
        , MAX(uo.uo_6h_ml) AS uo_6h_ml
        , MAX(uo.weight_backapplied) AS weight_backapplied
        , gcs.sedated_hour

        , CASE

            WHEN (MAX(cr.creatinine_max) IS NULL AND MAX(uo.uo_6h_ml) IS NULL) THEN 0

            WHEN (
                -- renal dysfunction trigger (SOFA-2 footnote p): Cr > 1.2 OR oliguria <0.3 mL/kg/h for >6h
                (MAX(cr.creatinine_max) > 1.2)
                OR LOGICAL_OR(uo.olig6_b)
                )
            AND (
                -- severity trigger: K >= 6 OR (pH <= 7.20 AND HCO3 <= 12)
                MAX(kh.potassium_max) >= 6.0
                OR (MIN(ph.ph_min) <= 7.20 AND MIN(kh.bicarbonate_min) <= 12.0)
                )
            THEN 1
            ELSE 0
        END AS rrt_eligible_flag
    FROM co
    LEFT JOIN pf  ON co.stay_id = pf.stay_id  AND co.hr = pf.hr
    LEFT JOIN ecmo_hourly eh ON co.stay_id = eh.stay_id AND co.hr = eh.hr
    LEFT JOIN vs  ON co.stay_id = vs.stay_id  AND co.hr = vs.hr
    LEFT JOIN gcs ON co.stay_id = gcs.stay_id AND co.hr = gcs.hr
    LEFT JOIN bili ON co.stay_id = bili.stay_id AND co.hr = bili.hr
    LEFT JOIN cr   ON co.stay_id = cr.stay_id   AND co.hr = cr.hr
    LEFT JOIN plt  ON co.stay_id = plt.stay_id  AND co.hr = plt.hr
    LEFT JOIN vaso ON co.stay_id = vaso.stay_id AND co.hr = vaso.hr
    LEFT JOIN rrt ON co.stay_id = rrt.stay_id AND rrt.charttime >= co.starttime
                                            AND rrt.charttime < co.endtime
    LEFT JOIN delirium_tx del
        ON co.stay_id = del.stay_id AND del.charttime >= co.starttime
                                    AND del.charttime < co.endtime
    LEFT JOIN uo_flags uo ON co.stay_id = uo.stay_id AND co.hr = uo.hr
    LEFT JOIN sf ON co.stay_id = sf.stay_id AND co.hr = sf.hr
    LEFT JOIN adv_hourly adv ON co.stay_id = adv.stay_id AND co.hr = adv.hr
    LEFT JOIN k_hco3 kh ON co.stay_id = kh.stay_id AND co.hr = kh.hr
    LEFT JOIN ph_lab ph ON co.stay_id = ph.stay_id AND co.hr = ph.hr
    GROUP BY
        co.stay_id, co.hr, co.starttime, co.endtime,co.intime,
        pf.pao2fio2ratio_min, pf.advanced_vent_any, --pf.ecmo_any,
        vs.meanbp_min, gcs.gcs_min, gcs.gcs_motor_min, gcs.gcs_unable_any, gcs.sedated_hour,
        bili.bilirubin_max, cr.creatinine_max, plt.platelet_min,
        vaso.rate_epinephrine_ugkgmin, vaso.rate_norepinephrine_ugkgmin,
        vaso.rate_dopamine_mcgkgmin, vaso.rate_dobutamine,
        vaso.other_vaso_flag, vaso.mechanical_support_flag
)

, scorecalc AS (
    SELECT s.*

    -- ===================== Respiratory
    , CASE

    -- ECMO always respiratory = 4
    WHEN ecmo_any = 1 THEN 4

    -- Prefer PaO2/FiO2 when available
    WHEN pao2fio2ratio_min IS NOT NULL THEN
      CASE
        WHEN advanced_vent_any = 1 AND pao2fio2ratio_min <= 75 THEN 4
        WHEN advanced_vent_any = 1 AND pao2fio2ratio_min <= 150 THEN 3
        WHEN pao2fio2ratio_min <= 225 THEN 2
        WHEN pao2fio2ratio_min <= 300 THEN 1
        ELSE 0
      END

    -- Fallback SpO2/FiO2
    WHEN spo2fio2ratio_min IS NOT NULL THEN
      CASE
        WHEN advanced_vent_hour = 1 AND spo2fio2ratio_min <= 120 THEN 4
        WHEN advanced_vent_hour = 1 AND spo2fio2ratio_min <= 200 THEN 3
        WHEN spo2fio2ratio_min <= 250 THEN 2
        WHEN spo2fio2ratio_min <= 300 THEN 1
        ELSE 0
      END

    ELSE NULL
  END AS respiration
    -- ===================== Hemostasis
    , CASE
        WHEN platelet_min <= 50 THEN 4
        WHEN platelet_min <= 80 THEN 3
        WHEN platelet_min <= 100 THEN 2
        WHEN platelet_min <= 150 THEN 1
        WHEN platelet_min IS NULL THEN NULL
        ELSE 0
      END AS hemostasis

    -- ===================== Liver
    , CASE
        WHEN bilirubin_max > 12.0 THEN 4
        WHEN bilirubin_max > 6.0 THEN 3
        WHEN bilirubin_max > 3.0 THEN 2
        WHEN bilirubin_max > 1.2 THEN 1
        WHEN bilirubin_max IS NULL THEN NULL
        ELSE 0
      END AS liver

    -- ===================== Brain
    , CASE
        -- Footnote c: if sedated and previous GCS unknown -> score 0
        WHEN sedated_hour = 1 AND gcs_min IS NULL THEN 0

        -- Footnote e: delirium treatment -> score at least 1 (even if GCS 15 / motor 6)
        WHEN delirium_treated_any = 1 THEN
        CASE
            WHEN gcs_min BETWEEN 3 AND 5 THEN 4
            WHEN gcs_min BETWEEN 6 AND 8 THEN 3
            WHEN gcs_min BETWEEN 9 AND 12 THEN 2
            WHEN gcs_min BETWEEN 13 AND 14 THEN 1
            WHEN gcs_min = 15 THEN 1

            -- Footnote d motor fallback if total unavailable
            WHEN gcs_motor_min <= 2 THEN 4
            WHEN gcs_motor_min = 3 THEN 3
            WHEN gcs_motor_min = 4 THEN 2
            WHEN gcs_motor_min = 5 THEN 1
            WHEN gcs_motor_min = 6 THEN 1
            ELSE 1
        END

        -- Total GCS available
        WHEN gcs_min IS NOT NULL THEN
        CASE
            WHEN gcs_min BETWEEN 3 AND 5 THEN 4
            WHEN gcs_min BETWEEN 6 AND 8 THEN 3
            WHEN gcs_min BETWEEN 9 AND 12 THEN 2
            WHEN gcs_min BETWEEN 13 AND 14 THEN 1
            WHEN gcs_min = 15 THEN 0
            ELSE NULL
        END

        -- Footnote d: cannot evaluate all 3 domains -> best achieved motor
        WHEN gcs_motor_min IS NOT NULL THEN
        CASE
            WHEN gcs_motor_min <= 2 THEN 4
            WHEN gcs_motor_min = 3 THEN 3
            WHEN gcs_motor_min = 4 THEN 2
            WHEN gcs_motor_min = 5 THEN 1
            WHEN gcs_motor_min = 6 THEN 0
            ELSE NULL
        END

        ELSE NULL
    END AS brain
    -- ===================== Cardiovascular
    , CASE
        WHEN mechanical_support_flag = 1 THEN 4
        WHEN (
            -- dopamine-only path (must be sole vaso)
            rate_dopamine_mcgkgmin > 0
            AND COALESCE(rate_norepinephrine_ugkgmin,0) = 0
            AND COALESCE(rate_epinephrine_ugkgmin,0) = 0
            AND COALESCE(rate_dobutamine,0) = 0
            AND COALESCE(other_vaso_flag,0) = 0
        ) THEN
            CASE
              WHEN rate_dopamine_mcgkgmin > 40 THEN 4
              WHEN rate_dopamine_mcgkgmin > 20 THEN 3
              ELSE 2
            END
        ELSE
          -- norepi+epi sum path
          CASE
            WHEN (COALESCE(rate_norepinephrine_ugkgmin,0) + COALESCE(rate_epinephrine_ugkgmin,0)) > 0.4
              THEN 4
            WHEN (
              (COALESCE(rate_norepinephrine_ugkgmin,0) + COALESCE(rate_epinephrine_ugkgmin,0)) > 0.2
              AND (COALESCE(rate_norepinephrine_ugkgmin,0) + COALESCE(rate_epinephrine_ugkgmin,0)) <= 0.4
              AND (COALESCE(other_vaso_flag,0) = 1 OR COALESCE(rate_dobutamine,0) > 0)
            ) THEN 4
            WHEN (
              (COALESCE(rate_norepinephrine_ugkgmin,0) + COALESCE(rate_epinephrine_ugkgmin,0)) > 0.2
              AND (COALESCE(rate_norepinephrine_ugkgmin,0) + COALESCE(rate_epinephrine_ugkgmin,0)) <= 0.4
            ) THEN 3
            WHEN (
              (COALESCE(rate_norepinephrine_ugkgmin,0) + COALESCE(rate_epinephrine_ugkgmin,0)) <= 0.2
              AND (COALESCE(rate_norepinephrine_ugkgmin,0) + COALESCE(rate_epinephrine_ugkgmin,0)) > 0
              AND (COALESCE(other_vaso_flag,0) = 1 OR COALESCE(rate_dobutamine,0) > 0)
            ) THEN 3
            WHEN (
              (COALESCE(rate_norepinephrine_ugkgmin,0) + COALESCE(rate_epinephrine_ugkgmin,0)) <= 0.2
              AND (COALESCE(rate_norepinephrine_ugkgmin,0) + COALESCE(rate_epinephrine_ugkgmin,0)) > 0
            ) THEN 2
            WHEN (COALESCE(other_vaso_flag,0) = 1 OR COALESCE(rate_dobutamine,0) > 0) THEN 2
            WHEN meanbp_min < 70
                AND rate_norepinephrine_ugkgmin IS NULL
                AND rate_epinephrine_ugkgmin IS NULL
                AND rate_dopamine_mcgkgmin IS NULL
                AND rate_dobutamine IS NULL
                AND COALESCE(other_vaso_flag,0) = 0
                THEN 1
            WHEN meanbp_min < 70 THEN 1
            WHEN meanbp_min IS NULL THEN NULL
            ELSE 0
          END
      END AS cardiovascular

    -- ===================== Kidney
    , CASE
        WHEN (on_rrt_any = 1 OR rrt_eligible_flag = 1) THEN 4
        WHEN uo_gt24_or_anuria_flag = 1 THEN 3
        WHEN creatinine_max > 3.5 THEN 3
        WHEN uo_gt12_flag = 1 THEN 2
        WHEN creatinine_max > 2.0 AND creatinine_max <= 3.5 THEN 2
        WHEN uo_6_12_flag = 1 THEN 1
        WHEN creatinine_max > 1.2 AND creatinine_max <= 2.0 THEN 1
        WHEN creatinine_max IS NULL AND on_rrt_any IS NULL
            AND uo_6_12_flag IS NULL AND uo_gt12_flag IS NULL AND uo_gt24_or_anuria_flag IS NULL
            AND rrt_eligible_flag IS NULL
            THEN NULL
        ELSE 0
        END AS kidney

    -- kidney scores: rrt_eligible and urine-output flags per weight rule (o, b) and coverage (g)
    , CASE
        WHEN (COALESCE(on_rrt_any,0) = 1 OR ((COALESCE(cr_trig,FALSE) OR COALESCE(olig6_o,FALSE)) AND COALESCE(rrt_sev,FALSE))) THEN 4
        WHEN uo_gt24_o = 1 THEN 3
        WHEN creatinine_max > 3.5 THEN 3
        WHEN uo_gt12_o = 1 THEN 2
        WHEN creatinine_max > 2.0 THEN 2
        WHEN uo_6_12_o = 1 THEN 1
        WHEN creatinine_max > 1.2 THEN 1
        WHEN creatinine_max IS NULL AND on_rrt_any IS NULL
            AND uo_6_12_o IS NULL AND uo_gt12_o IS NULL AND uo_gt24_o IS NULL
            THEN NULL
        ELSE 0 END AS kidney_o
    , CASE WHEN (COALESCE(on_rrt_any,0) = 1 OR ((COALESCE(cr_trig,FALSE) OR COALESCE(olig6_o,FALSE)) AND COALESCE(rrt_sev,FALSE))) THEN 4
        WHEN creatinine_max > 3.5 THEN 3 WHEN creatinine_max > 2.0 THEN 2 WHEN creatinine_max > 1.2 THEN 1 ELSE 0 END AS kidney_cr_o
    , GREATEST(
        CASE WHEN (COALESCE(on_rrt_any,0) = 1 OR ((COALESCE(cr_trig,FALSE) OR COALESCE(olig6_o,FALSE)) AND COALESCE(rrt_sev,FALSE))) THEN 4
          WHEN creatinine_max > 3.5 THEN 3 WHEN creatinine_max > 2.0 THEN 2 WHEN creatinine_max > 1.2 THEN 1 ELSE 0 END,
        CASE WHEN hr >= 22 AND COALESCE(uo_gt24_o,0) = 1 THEN 3
             WHEN hr >= 22 AND COALESCE(uo_gt12_o,0) = 1 THEN 2
             WHEN hr >= 22 AND COALESCE(uo_6_12_o,0) = 1 THEN 1 ELSE 0 END) AS kidney_uo22_o
    , CASE
        WHEN (COALESCE(on_rrt_any,0) = 1 OR ((COALESCE(cr_trig,FALSE) OR COALESCE(olig6_b,FALSE)) AND COALESCE(rrt_sev,FALSE))) THEN 4
        WHEN uo_gt24_b = 1 THEN 3
        WHEN creatinine_max > 3.5 THEN 3
        WHEN uo_gt12_b = 1 THEN 2
        WHEN creatinine_max > 2.0 THEN 2
        WHEN uo_6_12_b = 1 THEN 1
        WHEN creatinine_max > 1.2 THEN 1
        WHEN creatinine_max IS NULL AND on_rrt_any IS NULL
            AND uo_6_12_b IS NULL AND uo_gt12_b IS NULL AND uo_gt24_b IS NULL
            THEN NULL
        ELSE 0 END AS kidney_b
    , CASE WHEN (COALESCE(on_rrt_any,0) = 1 OR ((COALESCE(cr_trig,FALSE) OR COALESCE(olig6_b,FALSE)) AND COALESCE(rrt_sev,FALSE))) THEN 4
        WHEN creatinine_max > 3.5 THEN 3 WHEN creatinine_max > 2.0 THEN 2 WHEN creatinine_max > 1.2 THEN 1 ELSE 0 END AS kidney_cr_b
    , GREATEST(
        CASE WHEN (COALESCE(on_rrt_any,0) = 1 OR ((COALESCE(cr_trig,FALSE) OR COALESCE(olig6_b,FALSE)) AND COALESCE(rrt_sev,FALSE))) THEN 4
          WHEN creatinine_max > 3.5 THEN 3 WHEN creatinine_max > 2.0 THEN 2 WHEN creatinine_max > 1.2 THEN 1 ELSE 0 END,
        CASE WHEN hr >= 22 AND COALESCE(uo_gt24_b,0) = 1 THEN 3
             WHEN hr >= 22 AND COALESCE(uo_gt12_b,0) = 1 THEN 2
             WHEN hr >= 22 AND COALESCE(uo_6_12_b,0) = 1 THEN 1 ELSE 0 END) AS kidney_uo22_b
    , CASE
        WHEN (COALESCE(on_rrt_any,0) = 1 OR ((COALESCE(cr_trig,FALSE) OR COALESCE(olig6_go,FALSE)) AND COALESCE(rrt_sev,FALSE))) THEN 4
        WHEN uo_gt24_go = 1 THEN 3
        WHEN creatinine_max > 3.5 THEN 3
        WHEN uo_gt12_go = 1 THEN 2
        WHEN creatinine_max > 2.0 THEN 2
        WHEN uo_6_12_go = 1 THEN 1
        WHEN creatinine_max > 1.2 THEN 1
        WHEN creatinine_max IS NULL AND on_rrt_any IS NULL
            AND uo_6_12_go IS NULL AND uo_gt12_go IS NULL AND uo_gt24_go IS NULL
            THEN NULL
        ELSE 0 END AS kidney_go
    , CASE
        WHEN (COALESCE(on_rrt_any,0) = 1 OR ((COALESCE(cr_trig,FALSE) OR COALESCE(olig6_gb,FALSE)) AND COALESCE(rrt_sev,FALSE))) THEN 4
        WHEN uo_gt24_gb = 1 THEN 3
        WHEN creatinine_max > 3.5 THEN 3
        WHEN uo_gt12_gb = 1 THEN 2
        WHEN creatinine_max > 2.0 THEN 2
        WHEN uo_6_12_gb = 1 THEN 1
        WHEN creatinine_max > 1.2 THEN 1
        WHEN creatinine_max IS NULL AND on_rrt_any IS NULL
            AND uo_6_12_gb IS NULL AND uo_gt12_gb IS NULL AND uo_gt24_gb IS NULL
            THEN NULL
        ELSE 0 END AS kidney_gb
    , CASE
        WHEN (COALESCE(on_rrt_any,0) = 1 OR ((COALESCE(cr_trig,FALSE) OR COALESCE(olig6_so,FALSE)) AND COALESCE(rrt_sev,FALSE))) THEN 4
        WHEN uo_gt24_so = 1 THEN 3
        WHEN creatinine_max > 3.5 THEN 3
        WHEN uo_gt12_so = 1 THEN 2
        WHEN creatinine_max > 2.0 THEN 2
        WHEN uo_6_12_so = 1 THEN 1
        WHEN creatinine_max > 1.2 THEN 1
        WHEN creatinine_max IS NULL AND on_rrt_any IS NULL
            AND uo_6_12_so IS NULL AND uo_gt12_so IS NULL AND uo_gt24_so IS NULL
            THEN NULL
        ELSE 0 END AS kidney_so
    , CASE
        WHEN (COALESCE(on_rrt_any,0) = 1 OR ((COALESCE(cr_trig,FALSE) OR COALESCE(olig6_sb,FALSE)) AND COALESCE(rrt_sev,FALSE))) THEN 4
        WHEN uo_gt24_sb = 1 THEN 3
        WHEN creatinine_max > 3.5 THEN 3
        WHEN uo_gt12_sb = 1 THEN 2
        WHEN creatinine_max > 2.0 THEN 2
        WHEN uo_6_12_sb = 1 THEN 1
        WHEN creatinine_max > 1.2 THEN 1
        WHEN creatinine_max IS NULL AND on_rrt_any IS NULL
            AND uo_6_12_sb IS NULL AND uo_gt12_sb IS NULL AND uo_gt24_sb IS NULL
            THEN NULL
        ELSE 0 END AS kidney_sb

    FROM scorecomp s
)

, score_final AS (
    SELECT s.*
        , COALESCE(MAX(respiration)     OVER w, 0) AS respiration_24hours
        , COALESCE(MAX(hemostasis)      OVER w, 0) AS hemostasis_24hours
        , COALESCE(MAX(liver)           OVER w, 0) AS liver_24hours
        , COALESCE(MAX(cardiovascular)  OVER w, 0) AS cardiovascular_24hours
        , COALESCE(MAX(brain)           OVER w, 0) AS brain_24hours
        , COALESCE(MAX(kidney)          OVER w, 0) AS kidney_24hours
        -- pre-ICU part of the kidney 24-h maximum (urine-output criteria cannot fire before hour 5)
        , MAX(IF(hr < 0, kidney_b, NULL)) OVER w AS kidney_pre24

        , COALESCE(MAX(respiration)     OVER w, 0)
        + COALESCE(MAX(hemostasis)      OVER w, 0)
        + COALESCE(MAX(liver)           OVER w, 0)
        + COALESCE(MAX(cardiovascular)  OVER w, 0)
        + COALESCE(MAX(brain)           OVER w, 0)
        + COALESCE(MAX(kidney)          OVER w, 0)
        AS sofa2_24hours
    FROM scorecalc s
    WINDOW w AS (
        PARTITION BY stay_id
        ORDER BY hr
        ROWS BETWEEN 23 PRECEDING AND 0 FOLLOWING -- 23 hours before including the current row
    )
)
SELECT *
FROM score_final
WHERE hr >= 0;
