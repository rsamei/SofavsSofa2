-- Step 19: alternative sedation definitions for SOFA-2's brain rule (footnote c), one row per hour of sofa2_withintime.
--   sed_std    the primary definition (sofa2_withintime.sedated_hour): a propofol, midazolam, lorazepam or
--              dexmedetomidine infusion in icu.inputevents overlapping the hour; recomputed here as a check;
--   sed_nodex  the same without dexmedetomidine;
--   sed_rass   sed_std and the last Richmond Agitation-Sedation Scale (item 228096) charted in the 4 hours up to the
--              end of the hour is -2 or lower (moderate or deeper sedation); an hour without such a RASS is not sedated;
--   sed_light  sed_std except hours whose last RASS in the 4 hours up to the end of the hour is -1 or higher (light
--              sedation, GCS presumed assessable); hours without a RASS stay sedated.
-- No value after the end of the hour is used.
CREATE OR REPLACE TABLE `@BUILD@.s19_sed` AS
WITH co AS (SELECT stay_id, hr, starttime, endtime, sedated_hour FROM `@BUILD@.sofa2_withintime`),
items AS (
  SELECT itemid, REGEXP_CONTAINS(LOWER(label), r'\bdexmedetomidine\b') AS is_dex
  FROM `physionet-data.mimiciv_3_1_icu.d_items`
  WHERE REGEXP_CONTAINS(LOWER(label), r'\b(propofol|midazolam|lorazepam|dexmedetomidine)\b')
),
iv AS (
  SELECT ie.stay_id, ie.starttime, ie.endtime, it.is_dex
  FROM `physionet-data.mimiciv_3_1_icu.inputevents` ie JOIN items it USING (itemid)
  WHERE ie.stay_id IS NOT NULL AND ie.starttime IS NOT NULL AND ie.endtime IS NOT NULL
),
inf AS (
  SELECT co.stay_id, co.hr, COUNT(iv.stay_id) > 0 AS any_sed, COUNTIF(NOT iv.is_dex) > 0 AS nodex
  FROM co LEFT JOIN iv ON iv.stay_id = co.stay_id AND iv.starttime < co.endtime AND iv.endtime > co.starttime
  GROUP BY 1, 2
),
rass AS (
  SELECT co.stay_id, co.hr, ARRAY_AGG(c.valuenum ORDER BY c.charttime DESC LIMIT 1)[OFFSET(0)] AS last_rass
  FROM co JOIN `physionet-data.mimiciv_3_1_icu.chartevents` c
    ON c.stay_id = co.stay_id AND c.itemid = 228096 AND c.valuenum IS NOT NULL
   AND c.charttime > DATETIME_SUB(co.endtime, INTERVAL 4 HOUR) AND c.charttime <= co.endtime
  GROUP BY 1, 2
)
SELECT co.stay_id, co.hr, co.sedated_hour,
  CAST(COALESCE(inf.any_sed, FALSE) AS INT64) AS sed_std,
  CAST(COALESCE(inf.nodex, FALSE) AS INT64) AS sed_nodex,
  CAST(COALESCE(inf.any_sed, FALSE) AND COALESCE(rass.last_rass <= -2, FALSE) AS INT64) AS sed_rass,
  CAST(COALESCE(inf.any_sed, FALSE) AND NOT COALESCE(rass.last_rass > -2, FALSE) AS INT64) AS sed_light,
  rass.last_rass
FROM co LEFT JOIN inf USING (stay_id, hr) LEFT JOIN rass USING (stay_id, hr);
