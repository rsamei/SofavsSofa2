-- Hourly SOFA-2 cardiovascular component (Ranzani et al., JAMA 2025, Table 2, footnotes i-n).
-- Vasoactive rate per hour is the maximum over all infusion rows overlapping the hour. cardiovascular_v2_cut counts an
-- infusion in the hour only once it has lasted 60 minutes by the end of that hour (used at the prediction horizon).
CREATE OR REPLACE TABLE `@BUILD@.cv_hourly` AS
WITH hrs AS (SELECT stay_id, hr, starttime, endtime FROM `@BUILD@.sofa2_withintime`),
agents AS (
  SELECT 'ne' ag, stay_id, starttime, endtime, vaso_rate FROM `@UPSTREAM@.norepinephrine` WHERE vaso_rate > 0
  UNION ALL SELECT 'epi', stay_id, starttime, endtime, vaso_rate FROM `@UPSTREAM@.epinephrine` WHERE vaso_rate > 0
  UNION ALL SELECT 'dopa', stay_id, starttime, endtime, vaso_rate FROM `@UPSTREAM@.dopamine` WHERE vaso_rate > 0
  UNION ALL SELECT 'dobu', stay_id, starttime, endtime, vaso_rate FROM `@UPSTREAM@.dobutamine` WHERE vaso_rate > 0
  UNION ALL SELECT 'phenyl', stay_id, starttime, endtime, vaso_rate FROM `@UPSTREAM@.phenylephrine` WHERE vaso_rate > 0
  UNION ALL SELECT 'vaso', stay_id, starttime, endtime, vaso_rate FROM `@UPSTREAM@.vasopressin` WHERE vaso_rate > 0
  UNION ALL SELECT 'milr', stay_id, starttime, endtime, vaso_rate FROM `@UPSTREAM@.milrinone` WHERE vaso_rate > 0
),
-- continuous infusion episodes per agent (rows merged when the next row starts within 1 min of the previous end)
ep0 AS (
  SELECT *, IF(starttime > DATETIME_ADD(MAX(endtime) OVER (PARTITION BY ag, stay_id ORDER BY starttime, endtime ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING), INTERVAL 1 MINUTE), 1, 0) AS new_ep
  FROM agents),
ep AS (SELECT *, SUM(new_ep) OVER (PARTITION BY ag, stay_id ORDER BY starttime, endtime ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW) AS ep_id FROM ep0),
epdur AS (SELECT ag, stay_id, ep_id, MIN(starttime) ep_start, DATETIME_DIFF(MAX(endtime), MIN(starttime), MINUTE) ep_min FROM ep GROUP BY 1,2,3),
hourly AS (
  SELECT h.stay_id, h.hr,
    MAX(IF(e.ag='ne', e.vaso_rate, NULL)) ne, MAX(IF(e.ag='epi', e.vaso_rate, NULL)) epi, MAX(IF(e.ag='dopa', e.vaso_rate, NULL)) dopa,
    MAX(IF(e.ag='dobu', e.vaso_rate, NULL)) dobu, MAX(IF(e.ag='phenyl', e.vaso_rate, NULL)) phenyl, MAX(IF(e.ag='vaso', e.vaso_rate, NULL)) vaso, MAX(IF(e.ag='milr', e.vaso_rate, NULL)) milr
  FROM hrs h
  JOIN ep e ON e.stay_id = h.stay_id AND e.starttime < h.endtime AND e.endtime > h.starttime
  JOIN epdur d ON d.ag = e.ag AND d.stay_id = e.stay_id AND d.ep_id = e.ep_id AND d.ep_min >= 60
  GROUP BY 1,2),
hourly_cut AS (
  SELECT h.stay_id, h.hr,
    MAX(IF(e.ag='ne', e.vaso_rate, NULL)) ne, MAX(IF(e.ag='epi', e.vaso_rate, NULL)) epi, MAX(IF(e.ag='dopa', e.vaso_rate, NULL)) dopa,
    MAX(IF(e.ag='dobu', e.vaso_rate, NULL)) dobu, MAX(IF(e.ag='phenyl', e.vaso_rate, NULL)) phenyl, MAX(IF(e.ag='vaso', e.vaso_rate, NULL)) vaso, MAX(IF(e.ag='milr', e.vaso_rate, NULL)) milr
  FROM hrs h
  JOIN ep e ON e.stay_id = h.stay_id AND e.starttime < h.endtime AND e.endtime > h.starttime
  JOIN epdur d ON d.ag = e.ag AND d.stay_id = e.stay_id AND d.ep_id = e.ep_id AND d.ep_min >= 60
   AND DATETIME_ADD(d.ep_start, INTERVAL 60 MINUTE) <= h.endtime
  GROUP BY 1,2),
mapv AS (
  SELECT h.stay_id, h.hr, MIN(v.mbp) meanbp_min
  FROM hrs h JOIN `@UPSTREAM@.vitalsign` v ON v.stay_id = h.stay_id AND v.charttime > h.starttime AND v.charttime <= h.endtime
  WHERE v.mbp IS NOT NULL GROUP BY 1,2)
SELECT h.stay_id, h.hr, h.starttime, h.endtime, m.meanbp_min, v.ne, v.epi, v.dopa, v.dobu, v.phenyl, v.vaso, v.milr,
  vc.ne AS ne_cut, vc.epi AS epi_cut, vc.dopa AS dopa_cut, vc.dobu AS dobu_cut,
  s.mechanical_support_flag, s.ecmo_any,
  COALESCE(v.ne,0) + COALESCE(v.epi,0) AS ne_eq,
  (COALESCE(v.phenyl,0) > 0 OR COALESCE(v.vaso,0) > 0 OR COALESCE(v.dobu,0) > 0 OR COALESCE(v.milr,0) > 0
   OR (COALESCE(v.dopa,0) > 0 AND COALESCE(v.ne,0) + COALESCE(v.epi,0) > 0)) AS other_agent,
  CASE
    WHEN COALESCE(s.mechanical_support_flag,0) = 1 THEN 4
    WHEN COALESCE(v.ne,0) + COALESCE(v.epi,0) > 0.4 THEN 4
    WHEN COALESCE(v.ne,0) + COALESCE(v.epi,0) > 0.2 THEN IF(COALESCE(v.phenyl,0) > 0 OR COALESCE(v.vaso,0) > 0 OR COALESCE(v.dobu,0) > 0 OR COALESCE(v.milr,0) > 0 OR COALESCE(v.dopa,0) > 0, 4, 3)
    WHEN COALESCE(v.ne,0) + COALESCE(v.epi,0) > 0 THEN IF(COALESCE(v.phenyl,0) > 0 OR COALESCE(v.vaso,0) > 0 OR COALESCE(v.dobu,0) > 0 OR COALESCE(v.milr,0) > 0 OR COALESCE(v.dopa,0) > 0, 3, 2)
    WHEN COALESCE(v.dopa,0) > 40 THEN 4
    WHEN COALESCE(v.dopa,0) > 20 THEN 3
    WHEN COALESCE(v.dopa,0) > 0 THEN 2
    WHEN COALESCE(v.phenyl,0) > 0 OR COALESCE(v.vaso,0) > 0 OR COALESCE(v.dobu,0) > 0 OR COALESCE(v.milr,0) > 0 THEN 2
    WHEN m.meanbp_min < 70 THEN 1
    WHEN m.meanbp_min IS NULL THEN NULL
    ELSE 0 END AS cardiovascular_v2,
  CASE
    WHEN COALESCE(s.mechanical_support_flag,0) = 1 THEN 4
    WHEN COALESCE(vc.ne,0) + COALESCE(vc.epi,0) > 0.4 THEN 4
    WHEN COALESCE(vc.ne,0) + COALESCE(vc.epi,0) > 0.2 THEN IF(COALESCE(vc.phenyl,0) > 0 OR COALESCE(vc.vaso,0) > 0 OR COALESCE(vc.dobu,0) > 0 OR COALESCE(vc.milr,0) > 0 OR COALESCE(vc.dopa,0) > 0, 4, 3)
    WHEN COALESCE(vc.ne,0) + COALESCE(vc.epi,0) > 0 THEN IF(COALESCE(vc.phenyl,0) > 0 OR COALESCE(vc.vaso,0) > 0 OR COALESCE(vc.dobu,0) > 0 OR COALESCE(vc.milr,0) > 0 OR COALESCE(vc.dopa,0) > 0, 3, 2)
    WHEN COALESCE(vc.dopa,0) > 40 THEN 4
    WHEN COALESCE(vc.dopa,0) > 20 THEN 3
    WHEN COALESCE(vc.dopa,0) > 0 THEN 2
    WHEN COALESCE(vc.phenyl,0) > 0 OR COALESCE(vc.vaso,0) > 0 OR COALESCE(vc.dobu,0) > 0 OR COALESCE(vc.milr,0) > 0 THEN 2
    WHEN m.meanbp_min < 70 THEN 1
    WHEN m.meanbp_min IS NULL THEN NULL
    ELSE 0 END AS cardiovascular_v2_cut
FROM hrs h
JOIN `@BUILD@.sofa2_withintime` s USING (stay_id, hr)
LEFT JOIN hourly v USING (stay_id, hr)
LEFT JOIN hourly_cut vc USING (stay_id, hr)
LEFT JOIN mapv m USING (stay_id, hr);
