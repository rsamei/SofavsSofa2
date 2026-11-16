"""Step 10: implementation audit. A second, separately written implementation re-derives the hourly SOFA-2 inputs and
component scores from the mimic-code concept tables for a 1% random sample of ICU hours 0-47 of primary-set stays
(sample fixed by a hash of stay and hour), and compares them with the pipeline's tables. Only agreement counts leave
BigQuery; disagreement patterns are shown only when they hold at least 5 hours.

Audited per hour: platelet, bilirubin and creatinine values and the hemostasis and liver scores; respiration (ECMO,
PaO2/FiO2 with support at the blood gas, SpO2/FiO2 with support in the hour); the brain score before sedation and
delirium rules (GCS bands or best motor score); the cardiovascular score; the charted weight; and the urine-output
flags of the primary rule (given the pipeline's weight).
"""
import json
import sys

from bq_client import BUILD, OUT, PROJECT, UPSTREAM, BQError, query

MIN_CELL = 5
PRIM = "c.first_elig_stay AND NOT c.not_assessable_any"
HOSP, ICU = "physionet-data.mimiciv_3_1_hosp", "physionet-data.mimiciv_3_1_icu"

SQL = f"""
WITH smp AS (
  SELECT w.stay_id, c.hadm_id, w.hr, w.starttime AS t0, w.endtime AS t1
  FROM `{BUILD}.sofa2_withintime` w JOIN `{BUILD}.cohort` c USING (stay_id)
  WHERE {PRIM} AND w.hr BETWEEN 0 AND 47
    AND MOD(ABS(FARM_FINGERPRINT(CONCAT('audit-', CAST(w.stay_id AS STRING), '-', CAST(w.hr AS STRING)))), 100) = 0
),
lab AS (
  SELECT s.stay_id, s.hr, MIN(x.platelet) AS plt
  FROM smp s LEFT JOIN `{UPSTREAM}.complete_blood_count` x
    ON x.hadm_id = s.hadm_id AND x.charttime >= s.t0 AND x.charttime < s.t1
  GROUP BY 1, 2
),
lab2 AS (
  SELECT s.stay_id, s.hr, MAX(x.bilirubin_total) AS bili
  FROM smp s LEFT JOIN `{UPSTREAM}.enzyme` x ON x.hadm_id = s.hadm_id AND x.charttime >= s.t0 AND x.charttime < s.t1
  GROUP BY 1, 2
),
lab3 AS (
  SELECT s.stay_id, s.hr, MAX(x.creatinine) AS creat
  FROM smp s LEFT JOIN `{UPSTREAM}.chemistry` x ON x.hadm_id = s.hadm_id AND x.charttime >= s.t0 AND x.charttime < s.t1
  GROUP BY 1, 2
),
adv_ep AS (
  SELECT stay_id, starttime, endtime FROM `{UPSTREAM}.ventilation`
  WHERE ventilation_status IN ('InvasiveVent', 'NonInvasiveVent', 'HFNC', 'Tracheostomy')
),
gas AS (
  SELECT s.stay_id, s.hr, g.charttime, g.pao2fio2ratio, LOGICAL_OR(a.stay_id IS NOT NULL) AS sup
  FROM smp s JOIN `{UPSTREAM}.bg` g
    ON g.hadm_id = s.hadm_id AND g.specimen = 'ART.' AND g.charttime >= s.t0 AND g.charttime < s.t1
  LEFT JOIN adv_ep a ON a.stay_id = s.stay_id AND g.charttime >= a.starttime AND g.charttime < a.endtime
  GROUP BY 1, 2, 3, 4
),
pf AS (SELECT stay_id, hr, MIN(pao2fio2ratio) AS pf, LOGICAL_OR(sup) AS pf_sup FROM gas GROUP BY 1, 2),
fio2 AS (
  SELECT stay_id, charttime, IF(fio2 > 1, fio2 / 100, fio2) AS f
  FROM `{UPSTREAM}.ventilator_setting` WHERE fio2 IS NOT NULL
),
sfr AS (
  SELECT s.stay_id, s.hr, MIN(v.spo2 / f.f) AS sf
  FROM smp s
  JOIN `{UPSTREAM}.vitalsign` v ON v.stay_id = s.stay_id AND v.charttime >= s.t0 AND v.charttime < s.t1 AND v.spo2 < 98
  JOIN fio2 f ON f.stay_id = s.stay_id AND f.charttime <= v.charttime
             AND f.charttime >= DATETIME_SUB(v.charttime, INTERVAL 60 MINUTE) AND f.f BETWEEN 0.21 AND 1.0
  GROUP BY 1, 2
),
ecmo_h AS (
  SELECT DISTINCT s.stay_id, s.hr FROM smp s JOIN `{ICU}.chartevents` ce
    ON ce.stay_id = s.stay_id AND ce.itemid = 229268 AND LOWER(ce.value) IN ('vv', 'va', 'vav')
   AND ce.charttime >= s.t0 AND ce.charttime < s.t1
),
sup_h AS (
  SELECT DISTINCT s.stay_id, s.hr FROM smp s JOIN adv_ep a
    ON a.stay_id = s.stay_id AND a.starttime < s.t1 AND a.endtime > s.t0
),
resp AS (
  SELECT s.stay_id, s.hr, e.stay_id IS NOT NULL AS ecmo, u.stay_id IS NOT NULL AS sup_hr, pf.pf, pf.pf_sup, sfr.sf
  FROM smp s LEFT JOIN pf USING (stay_id, hr) LEFT JOIN sfr USING (stay_id, hr)
  LEFT JOIN ecmo_h e USING (stay_id, hr) LEFT JOIN sup_h u USING (stay_id, hr)
),
resp_score AS (
  SELECT stay_id, hr,
    CASE WHEN ecmo THEN 4
         WHEN pf IS NOT NULL THEN CASE WHEN pf_sup AND pf <= 75 THEN 4 WHEN pf_sup AND pf <= 150 THEN 3
                                       WHEN pf <= 225 THEN 2 WHEN pf <= 300 THEN 1 ELSE 0 END
         WHEN sf IS NOT NULL THEN CASE WHEN sup_hr AND sf <= 120 THEN 4 WHEN sup_hr AND sf <= 200 THEN 3
                                       WHEN sf <= 250 THEN 2 WHEN sf <= 300 THEN 1 ELSE 0 END
    END AS resp
  FROM resp
),
gcs_h AS (
  SELECT s.stay_id, s.hr, MIN(g.gcs) AS gcs, MAX(g.gcs_unable) AS unable, MAX(g.gcs_motor) AS motor
  FROM smp s JOIN `{UPSTREAM}.gcs` g ON g.stay_id = s.stay_id AND g.charttime > s.t0 AND g.charttime <= s.t1
  GROUP BY 1, 2
),
brain AS (
  SELECT stay_id, hr,
    CASE WHEN COALESCE(unable, 0) = 0 AND gcs IS NOT NULL THEN
           CASE WHEN gcs = 15 THEN 0 WHEN gcs >= 13 THEN 1 WHEN gcs >= 9 THEN 2 WHEN gcs >= 6 THEN 3 ELSE 4 END
         WHEN motor IS NOT NULL THEN
           CASE WHEN motor = 6 THEN 0 WHEN motor = 5 THEN 1 WHEN motor = 4 THEN 2 WHEN motor = 3 THEN 3 ELSE 4 END
    END AS brain
  FROM gcs_h
),
inf AS (
  SELECT 'ne' AS ag, stay_id, starttime, endtime, vaso_rate FROM `{UPSTREAM}.norepinephrine`
  UNION ALL SELECT 'epi', stay_id, starttime, endtime, vaso_rate FROM `{UPSTREAM}.epinephrine`
  UNION ALL SELECT 'dopa', stay_id, starttime, endtime, vaso_rate FROM `{UPSTREAM}.dopamine`
  UNION ALL SELECT 'dobu', stay_id, starttime, endtime, vaso_rate FROM `{UPSTREAM}.dobutamine`
  UNION ALL SELECT 'phenyl', stay_id, starttime, endtime, vaso_rate FROM `{UPSTREAM}.phenylephrine`
  UNION ALL SELECT 'vaso', stay_id, starttime, endtime, vaso_rate FROM `{UPSTREAM}.vasopressin`
  UNION ALL SELECT 'milr', stay_id, starttime, endtime, vaso_rate FROM `{UPSTREAM}.milrinone`
),
inf_run AS (  -- continuous runs per agent: a new run starts when a row begins > 1 min after every earlier row ended
  SELECT *, SUM(brk) OVER (PARTITION BY ag, stay_id ORDER BY starttime, endtime) AS run_id
  FROM (
    SELECT i.*,
      IF(starttime > DATETIME_ADD(MAX(endtime) OVER (PARTITION BY ag, stay_id ORDER BY starttime, endtime
                                                    ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING), INTERVAL 1 MINUTE)
         OR LAG(endtime) OVER (PARTITION BY ag, stay_id ORDER BY starttime, endtime) IS NULL, 1, 0) AS brk
    FROM inf i WHERE vaso_rate > 0
  )
),
long_runs AS (
  SELECT ag, stay_id, run_id FROM inf_run GROUP BY 1, 2, 3
  HAVING DATETIME_DIFF(MAX(endtime), MIN(starttime), MINUTE) >= 60
),
dose AS (
  SELECT s.stay_id, s.hr,
    MAX(IF(r.ag = 'ne', r.vaso_rate, NULL)) AS ne, MAX(IF(r.ag = 'epi', r.vaso_rate, NULL)) AS epi,
    MAX(IF(r.ag = 'dopa', r.vaso_rate, NULL)) AS dopa, LOGICAL_OR(r.ag IN ('dobu', 'phenyl', 'vaso', 'milr')) AS other
  FROM smp s
  JOIN inf_run r ON r.stay_id = s.stay_id AND r.starttime < s.t1 AND r.endtime > s.t0
  JOIN long_runs l ON l.ag = r.ag AND l.stay_id = r.stay_id AND l.run_id = r.run_id
  GROUP BY 1, 2
),
mcs_items AS (
  SELECT itemid FROM `{ICU}.d_items`
  WHERE LOWER(label) LIKE '%iabp%' OR LOWER(label) LIKE '%impella%' OR LOWER(label) LIKE '%lvad%' OR LOWER(label) LIKE '% vad%'
),
map_h AS (
  SELECT s.stay_id, s.hr, MIN(v.mbp) AS map
  FROM smp s JOIN `{UPSTREAM}.vitalsign` v ON v.stay_id = s.stay_id AND v.charttime > s.t0 AND v.charttime <= s.t1
  GROUP BY 1, 2
),
mcs_h AS (
  SELECT DISTINCT s.stay_id, s.hr FROM smp s JOIN `{ICU}.chartevents` ce
    ON ce.stay_id = s.stay_id AND ce.charttime >= s.t0 AND ce.charttime < s.t1
  WHERE ce.itemid IN (SELECT itemid FROM mcs_items) OR (ce.itemid = 229268 AND LOWER(ce.value) IN ('va', 'vav'))
),
cv AS (
  SELECT s.stay_id, s.hr, d.ne, d.epi, d.dopa, COALESCE(d.other, FALSE) AS other, m.map, mc.stay_id IS NOT NULL AS mcs
  FROM smp s LEFT JOIN dose d USING (stay_id, hr) LEFT JOIN map_h m USING (stay_id, hr) LEFT JOIN mcs_h mc USING (stay_id, hr)
),
cv_score AS (
  SELECT stay_id, hr,
    CASE WHEN mcs THEN 4
         WHEN COALESCE(ne, 0) + COALESCE(epi, 0) > 0.4 THEN 4
         WHEN COALESCE(ne, 0) + COALESCE(epi, 0) > 0.2 THEN IF(other OR COALESCE(dopa, 0) > 0, 4, 3)
         WHEN COALESCE(ne, 0) + COALESCE(epi, 0) > 0 THEN IF(other OR COALESCE(dopa, 0) > 0, 3, 2)
         WHEN COALESCE(dopa, 0) > 40 THEN 4 WHEN COALESCE(dopa, 0) > 20 THEN 3 WHEN COALESCE(dopa, 0) > 0 THEN 2
         WHEN other THEN 2
         WHEN map < 70 THEN 1 WHEN map IS NULL THEN NULL ELSE 0 END AS cv
  FROM cv
),
wt AS (
  SELECT s.stay_id, s.hr, ARRAY_AGG(ce.valuenum IGNORE NULLS ORDER BY ce.charttime DESC, ce.valuenum DESC LIMIT 1)[SAFE_OFFSET(0)] AS weight
  FROM smp s LEFT JOIN `{ICU}.chartevents` ce
    ON ce.stay_id = s.stay_id AND ce.itemid IN (226512, 224639) AND ce.valuenum > 0 AND ce.charttime <= s.t1
  GROUP BY 1, 2
),
uo_iv AS (
  SELECT u.stay_id, u.charttime, u.urineoutput,
    DATETIME_DIFF(u.charttime, IFNULL(LAG(u.charttime) OVER (PARTITION BY u.stay_id ORDER BY u.charttime), t.intime_hr), MINUTE) AS iv
  FROM `{UPSTREAM}.urine_output` u JOIN `{UPSTREAM}.icustay_times` t USING (stay_id)
  WHERE u.stay_id IN (SELECT stay_id FROM smp)
),
uo_win AS (
  SELECT s.stay_id, s.hr, h.weight_ba_kg AS wkg,
    SUM(IF(u.charttime >= DATETIME_SUB(s.t1, INTERVAL 6 HOUR), u.urineoutput, 0)) AS v6,
    SUM(IF(u.charttime >= DATETIME_SUB(s.t1, INTERVAL 6 HOUR), u.iv, 0)) AS m6,
    SUM(IF(u.charttime >= DATETIME_SUB(s.t1, INTERVAL 12 HOUR), u.urineoutput, 0)) AS v12,
    SUM(IF(u.charttime >= DATETIME_SUB(s.t1, INTERVAL 12 HOUR), u.iv, 0)) AS m12,
    COUNTIF(u.charttime >= DATETIME_SUB(s.t1, INTERVAL 12 HOUR)) AS n12,
    SUM(u.urineoutput) AS v24, SUM(u.iv) AS m24
  FROM smp s
  JOIN `{BUILD}.sofa2_withintime` h ON h.stay_id = s.stay_id AND h.hr = s.hr
  LEFT JOIN uo_iv u ON u.stay_id = s.stay_id AND u.charttime >= DATETIME_SUB(s.t1, INTERVAL 24 HOUR) AND u.charttime < s.t1
  GROUP BY 1, 2, 3
),
uo AS (
  SELECT stay_id, hr,
    CASE WHEN wkg IS NULL OR wkg <= 0 OR m6 < 360 OR hr < 5 THEN NULL
         WHEN v6 / (m6 / 60) / wkg < 0.5 AND (hr < 11 OR m12 < 720 OR v12 / (m12 / 60) / wkg >= 0.5) THEN 1 ELSE 0 END AS f612,
    CASE WHEN wkg IS NULL OR wkg <= 0 OR m12 < 720 OR hr < 11 THEN NULL
         WHEN v12 / (m12 / 60) / wkg < 0.5 THEN 1 ELSE 0 END AS f12,
    CASE WHEN wkg IS NULL OR wkg <= 0 THEN NULL
         WHEN hr >= 11 AND m12 >= 720 AND n12 >= 2 AND v12 = 0 THEN 1
         WHEN hr >= 23 AND m24 >= 1440 AND v24 / (m24 / 60) / wkg < 0.3 THEN 1 ELSE 0 END AS f24
  FROM uo_win
),
cmp AS (
  SELECT s.stay_id, s.hr,
    STRUCT(l.plt AS a, w.platelet_min AS p) AS platelet, STRUCT(l2.bili AS a, w.bilirubin_max AS p) AS bilirubin,
    STRUCT(l3.creat AS a, w.creatinine_max AS p) AS creatinine,
    STRUCT(CASE WHEN l.plt <= 50 THEN 4 WHEN l.plt <= 80 THEN 3 WHEN l.plt <= 100 THEN 2 WHEN l.plt <= 150 THEN 1
                WHEN l.plt IS NOT NULL THEN 0 END AS a, w.hemostasis AS p) AS hemostasis,
    STRUCT(CASE WHEN l2.bili > 12 THEN 4 WHEN l2.bili > 6 THEN 3 WHEN l2.bili > 3 THEN 2 WHEN l2.bili > 1.2 THEN 1
                WHEN l2.bili IS NOT NULL THEN 0 END AS a, w.liver AS p) AS liver,
    STRUCT(rs.resp AS a, w.respiration AS p) AS respiration,
    STRUCT(b.brain AS a, h.brain_raw AS p) AS brain,
    STRUCT(cs.cv AS a, c.cardiovascular_v2 AS p) AS cardiovascular,
    STRUCT(wt.weight AS a, w.weight_own_kg AS p) AS weight,
    STRUCT(uo.f612 AS a, w.uo_6_12_b AS p) AS uo_6_12,
    STRUCT(uo.f12 AS a, w.uo_gt12_b AS p) AS uo_gt12,
    STRUCT(uo.f24 AS a, w.uo_gt24_b AS p) AS uo_gt24
  FROM smp s
  JOIN `{BUILD}.sofa2_withintime` w USING (stay_id, hr)
  JOIN `{BUILD}.sofa2_hourly` h USING (stay_id, hr)
  JOIN `{BUILD}.cv_hourly` c USING (stay_id, hr)
  LEFT JOIN lab l USING (stay_id, hr) LEFT JOIN lab2 l2 USING (stay_id, hr) LEFT JOIN lab3 l3 USING (stay_id, hr)
  LEFT JOIN resp_score rs USING (stay_id, hr) LEFT JOIN brain b USING (stay_id, hr)
  LEFT JOIN cv_score cs USING (stay_id, hr) LEFT JOIN wt USING (stay_id, hr) LEFT JOIN uo USING (stay_id, hr)
)
SELECT item, COUNT(*) AS hours, COUNTIF(a IS NOT NULL OR p IS NOT NULL) AS hours_with_value,
  COUNTIF(ROUND(CAST(a AS FLOAT64), 4) IS NOT DISTINCT FROM ROUND(CAST(p AS FLOAT64), 4)) AS agree
FROM cmp, UNNEST([
  STRUCT('platelet' AS item, platelet.a AS a, platelet.p AS p), ('bilirubin', bilirubin.a, bilirubin.p),
  ('creatinine', creatinine.a, creatinine.p), ('hemostasis', hemostasis.a, hemostasis.p), ('liver', liver.a, liver.p),
  ('respiration', respiration.a, respiration.p), ('brain', brain.a, brain.p), ('cardiovascular', cardiovascular.a, cardiovascular.p),
  ('weight', weight.a, weight.p), ('uo_6_12', uo_6_12.a, uo_6_12.p), ('uo_gt12', uo_gt12.a, uo_gt12.p), ('uo_gt24', uo_gt24.a, uo_gt24.p)])
GROUP BY item ORDER BY item
"""


def main():
    try:
        fields, rows, _ = query(PROJECT, SQL)
        dfields, drows, _ = query(PROJECT, DIS_SQL)
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    res = {"sample": "1% of ICU hours 0-47 of primary-set stays", "items": {}, "disagreement_patterns": {}}
    for r in rows:
        d = dict(zip(fields, r))
        res["items"][d["item"]] = {"hours": int(d["hours"]), "hours_with_value": int(d["hours_with_value"]),
                                   "agree": int(d["agree"]), "disagree": int(d["hours"]) - int(d["agree"])}
    for r in drows:
        d = dict(zip(dfields, r))
        if int(d["n"]) >= MIN_CELL:
            res["disagreement_patterns"].setdefault(d["item"], []).append(
                {"audit": d["audit"], "pipeline": d["pipeline"], "hours": int(d["n"])})
    (OUT / "audit.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    for k, v in res["items"].items():
        print(f"{k:15s} hours {v['hours']:6d}  agree {v['agree']:6d}  disagree {v['disagree']:5d}")
    return 0


# disagreement patterns for categorical items (scores and flags), counts only
DIS_SQL = SQL.split("SELECT item, COUNT(*)")[0] + """
SELECT item, audit, pipeline, COUNT(*) AS n
FROM (
  SELECT x.item, CAST(x.a AS STRING) AS audit, CAST(x.p AS STRING) AS pipeline
  FROM cmp, UNNEST([
    STRUCT('hemostasis' AS item, hemostasis.a AS a, hemostasis.p AS p), ('liver', liver.a, liver.p),
    ('respiration', respiration.a, respiration.p), ('brain', brain.a, brain.p),
    ('cardiovascular', cardiovascular.a, cardiovascular.p), ('uo_6_12', uo_6_12.a, uo_6_12.p),
    ('uo_gt12', uo_gt12.a, uo_gt12.p), ('uo_gt24', uo_gt24.a, uo_gt24.p)]) x
  WHERE x.a IS DISTINCT FROM x.p)
GROUP BY 1, 2, 3 ORDER BY 1, 4 DESC
"""


if __name__ == "__main__":
    sys.exit(main())
