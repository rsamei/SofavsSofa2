"""Step 6: descriptive tables for the primary set.
  - Table 1 by concordance group: n, age, female, suspicion before ICU, day-1 SOFA / SOFA-2, ICU and hospital
    length of stay (median and IQR), stays with a negative hospital length of stay;
  - in-hospital mortality by day-1 score point for both scores (Wilson CI; points >= 16 merged);
  - day-1 measurement of each component by score and concordance group (output-controlled counts; SOFA-only merged
    with Neither);
  - urine-output documentation for stays whose first-window SOFA-2 kidney point comes from urine output alone;
  - day-1 weight: back-applied admission weight hours used and excluded;
  - eMAR coverage by era and anchor_year_group.
Output control: counts < 5 suppressed; distribution summaries only for groups >= 5."""
import json
import math
import sys

from bq_client import BUILD, OUT, PROJECT, UPSTREAM, BQError, query

C = f"{BUILD}.cohort"
H = BUILD
MIN_CELL = 5
PRIM = "first_elig_stay AND NOT not_assessable_any"
GRP = ("CASE WHEN any_s1_std AND any_s2_std THEN 'BOTH' WHEN any_s1_std THEN 'SOFA_ONLY' "
       "WHEN any_s2_std THEN 'SOFA2_ONLY' ELSE 'NEITHER' END")


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(p, 6), round(c - h, 6), round(c + h, 6)]


def run(sql):
    fields, rows, _ = query(PROJECT, sql)
    return [dict(zip(fields, r)) for r in rows]


def f(v, nd=4):
    return None if v is None else round(float(v), nd)


def main():
    res = {"suppressed": []}
    try:
        t1 = run(f"""
          WITH g AS (SELECT {GRP} grp, * FROM `{C}` WHERE {PRIM}),
          a AS (SELECT 'ALL' grp, * EXCEPT (grp) FROM g UNION ALL SELECT * FROM g),
          q AS (
            SELECT grp,
              PERCENTILE_CONT(d1_s1_std, 0.25) OVER w AS s1_q1, PERCENTILE_CONT(d1_s1_std, 0.5) OVER w AS s1_q2,
              PERCENTILE_CONT(d1_s1_std, 0.75) OVER w AS s1_q3,
              PERCENTILE_CONT(d1_s2_std, 0.25) OVER w AS s2_q1, PERCENTILE_CONT(d1_s2_std, 0.5) OVER w AS s2_q2,
              PERCENTILE_CONT(d1_s2_std, 0.75) OVER w AS s2_q3,
              PERCENTILE_CONT(icu_los_days, 0.25) OVER w AS icu_q1, PERCENTILE_CONT(icu_los_days, 0.5) OVER w AS icu_q2,
              PERCENTILE_CONT(icu_los_days, 0.75) OVER w AS icu_q3,
              PERCENTILE_CONT(hosp_los_days, 0.25) OVER w AS hosp_q1, PERCENTILE_CONT(hosp_los_days, 0.5) OVER w AS hosp_q2,
              PERCENTILE_CONT(hosp_los_days, 0.75) OVER w AS hosp_q3
            FROM a WINDOW w AS (PARTITION BY grp)),
          qq AS (SELECT grp, ANY_VALUE(s1_q1) AS s1_q1, ANY_VALUE(s1_q2) AS s1_q2, ANY_VALUE(s1_q3) AS s1_q3, ANY_VALUE(s2_q1) AS s2_q1, ANY_VALUE(s2_q2) AS s2_q2, ANY_VALUE(s2_q3) AS s2_q3, ANY_VALUE(icu_q1) AS icu_q1, ANY_VALUE(icu_q2) AS icu_q2, ANY_VALUE(icu_q3) AS icu_q3, ANY_VALUE(hosp_q1) AS hosp_q1, ANY_VALUE(hosp_q2) AS hosp_q2, ANY_VALUE(hosp_q3) AS hosp_q3 FROM q GROUP BY grp)
          SELECT grp, COUNT(*) n, AVG(age) age_mean, STDDEV(age) age_sd, AVG(IF(gender = 'F', 1, 0)) female,
            AVG(IF(suspicion_before_icu, 1, 0)) susp_before_icu,
            AVG(d1_s1_std) d1_sofa_mean, STDDEV(d1_s1_std) d1_sofa_sd,
            AVG(d1_s2_std) d1_sofa2_mean, STDDEV(d1_s2_std) d1_sofa2_sd,
            SUM(hospital_expire_flag) hosp_deaths,
            COUNTIF(hosp_los_days < 0) neg_hosp_los, ANY_VALUE(qq.s1_q1) s1_q1, ANY_VALUE(qq.s1_q2) s1_q2, ANY_VALUE(qq.s1_q3) s1_q3, ANY_VALUE(qq.s2_q1) s2_q1, ANY_VALUE(qq.s2_q2) s2_q2, ANY_VALUE(qq.s2_q3) s2_q3, ANY_VALUE(qq.icu_q1) icu_q1, ANY_VALUE(qq.icu_q2) icu_q2, ANY_VALUE(qq.icu_q3) icu_q3, ANY_VALUE(qq.hosp_q1) hosp_q1, ANY_VALUE(qq.hosp_q2) hosp_q2, ANY_VALUE(qq.hosp_q3) hosp_q3
          FROM a JOIN qq USING (grp) GROUP BY grp""")
        res["table1"] = {}
        for r in t1:
            q = {k: r[k] for k in r if k[:3] in ("s1_", "s2_", "icu", "hos") and k[-3:-1] == "_q"}
            n = int(r["n"])
            res["table1"][r["grp"]] = {
                "n": n, "age": [f(r["age_mean"], 1), f(r["age_sd"], 1)], "female": f(r["female"], 3),
                "suspicion_before_icu": f(r["susp_before_icu"], 3),
                "hospital_mortality": f(int(r["hosp_deaths"]) / n, 3),
                "day1_sofa": {"mean": f(r["d1_sofa_mean"]), "sd": f(r["d1_sofa_sd"]),
                              "median_iqr": [f(q.get("s1_q2")), f(q.get("s1_q1")), f(q.get("s1_q3"))]},
                "day1_sofa2": {"mean": f(r["d1_sofa2_mean"]), "sd": f(r["d1_sofa2_sd"]),
                               "median_iqr": [f(q.get("s2_q2")), f(q.get("s2_q1")), f(q.get("s2_q3"))]},
                "icu_los_days_median_iqr": [f(q.get("icu_q2")), f(q.get("icu_q1")), f(q.get("icu_q3"))],
                "hosp_los_days_median_iqr": [f(q.get("hosp_q2")), f(q.get("hosp_q1")), f(q.get("hosp_q3"))],
                "negative_hosp_los_stays": (lambda k: k if k == 0 or k >= MIN_CELL else "<5")(int(r["neg_hosp_los"]))}
        # ---- case mix: shares, output-controlled (share withheld when the count or its complement is 1-4) -----
        cm = run(f"""
          WITH g AS (SELECT {GRP} grp, c.* FROM `{C}` c WHERE {PRIM}),
          a AS (SELECT 'ALL' grp, * EXCEPT (grp) FROM g UNION ALL SELECT * FROM g),
          vent AS (
            SELECT DISTINCT v.stay_id FROM `{UPSTREAM}.ventilation` v JOIN `{C}` c USING (stay_id)
            WHERE v.ventilation_status = 'InvasiveVent' AND v.starttime < DATETIME_ADD(c.intime_hr, INTERVAL 24 HOUR)
              AND v.endtime > c.intime_hr),
          vaso AS (
            SELECT DISTINCT stay_id FROM `{BUILD}.cv_hourly`
            WHERE hr BETWEEN 0 AND 23 AND (COALESCE(ne, 0) > 0 OR COALESCE(epi, 0) > 0 OR COALESCE(dopa, 0) > 0
                  OR COALESCE(phenyl, 0) > 0 OR COALESCE(vaso, 0) > 0)),
          x AS (
            SELECT a.grp, a.stay_id, b.admission_type_raw, a.first_careunit, ch.charlson_comorbidity_index AS cci,
              vent.stay_id IS NOT NULL AS imv_d1, vaso.stay_id IS NOT NULL AS vaso_d1,
              CASE WHEN a.first_careunit LIKE '%CVICU%' THEN 'cardiac_vascular' WHEN a.first_careunit LIKE '%(CCU)%' THEN 'coronary'
                   WHEN a.first_careunit LIKE '%Neuro%' THEN 'neuro' WHEN a.first_careunit LIKE '%MICU/SICU%' THEN 'medical_surgical'
                   WHEN a.first_careunit LIKE '%MICU%' THEN 'medical' WHEN a.first_careunit LIKE '%SICU%' THEN 'surgical_trauma'
                   ELSE 'medical_surgical' END AS icu_type  -- the few other units are grouped with medical/surgical
            FROM a JOIN `{UPSTREAM}.cohort_base` b USING (stay_id)
            LEFT JOIN `{UPSTREAM}.charlson` ch ON ch.hadm_id = a.hadm_id
            LEFT JOIN vent USING (stay_id) LEFT JOIN vaso USING (stay_id)),
          qc AS (SELECT grp, ANY_VALUE(q1) q1, ANY_VALUE(q2) q2, ANY_VALUE(q3) q3 FROM (
            SELECT grp, PERCENTILE_CONT(cci, 0.25) OVER w q1, PERCENTILE_CONT(cci, 0.5) OVER w q2,
                   PERCENTILE_CONT(cci, 0.75) OVER w q3 FROM x WINDOW w AS (PARTITION BY grp)) GROUP BY grp)
          SELECT grp, COUNT(*) n, COUNTIF(imv_d1) imv_d1, COUNTIF(vaso_d1) vaso_d1,
            COUNTIF(admission_type_raw IN ('EW EMER.', 'DIRECT EMER.', 'URGENT')) emergency_admission,
            COUNTIF(icu_type = 'medical') icu_medical, COUNTIF(icu_type = 'surgical_trauma') icu_surgical_trauma,
            COUNTIF(icu_type = 'medical_surgical') icu_medical_surgical, COUNTIF(icu_type = 'cardiac_vascular') icu_cardiac_vascular,
            COUNTIF(icu_type = 'coronary') icu_coronary, COUNTIF(icu_type = 'neuro') icu_neuro,
            ANY_VALUE(qc.q1) cci_q1, ANY_VALUE(qc.q2) cci_q2, ANY_VALUE(qc.q3) cci_q3
          FROM x JOIN qc USING (grp) GROUP BY grp""")
        for r in cm:
            n = int(r["n"])
            row = {}
            for k in ("imv_d1", "vaso_d1", "emergency_admission", "icu_medical", "icu_surgical_trauma",
                      "icu_medical_surgical", "icu_cardiac_vascular", "icu_coronary", "icu_neuro"):
                k_ = int(r[k])
                if 0 < k_ < MIN_CELL or 0 < n - k_ < MIN_CELL:
                    row[k] = "withheld"
                    res["suppressed"].append(f"case mix {r['grp']} {k}")
                else:
                    row[k] = f(k_ / n, 3)
            icu = [k for k in row if k.startswith("icu_")]
            if any(row[k] == "withheld" for k in icu):  # shares sum to 1: one withheld cell is recoverable
                for k in icu:
                    row[k] = "withheld"
            row["charlson_median_iqr"] = [f(r["cci_q2"]), f(r["cci_q1"]), f(r["cci_q3"])]
            res["table1"][r["grp"]]["case_mix"] = row
        # ---- mortality by day-1 score point ------------------------------------------------
        sp = run(f"""
          SELECT sys, LEAST(x, 16) pt, COUNT(*) n, SUM(hospital_expire_flag) d FROM (
            SELECT 'SOFA' sys, d1_s1_std x, hospital_expire_flag FROM `{C}` WHERE {PRIM}
            UNION ALL SELECT 'SOFA2', d1_s2_std, hospital_expire_flag FROM `{C}` WHERE {PRIM})
          GROUP BY 1, 2 ORDER BY 1, 2""")
        res["mortality_by_day1_point"] = {}
        for r in sp:
            n, d = int(r["n"]), int(r["d"])
            key = f"{r['sys']}|{'>=16' if int(r['pt']) == 16 else int(r['pt'])}"
            if n < MIN_CELL or d < MIN_CELL or n - d < MIN_CELL:
                res["mortality_by_day1_point"][key] = {"n": n if n >= MIN_CELL else "<5", "mortality": "suppressed"}
                res["suppressed"].append(f"score point {key}")
            else:
                res["mortality_by_day1_point"][key] = {"n": n, "deaths": d, "mortality": wilson(d, n)}
        # ---- day-1 component missingness by group -------------------------------------------
        ms = run(f"""
          WITH m AS (
            SELECT h1.stay_id,
              MAX(m1_resp) m1_resp, MAX(m1_coag) m1_coag, MAX(m1_liver) m1_liver, MAX(m1_cv) m1_cv, MAX(m1_cns) m1_cns, MAX(m1_renal) m1_renal,
              MAX(m2_resp) m2_resp, MAX(m2_hemo) m2_hemo, MAX(m2_liver) m2_liver, MAX(m2_cv) m2_cv, MAX(m2_brain) m2_brain, MAX(m2_kidney) m2_kidney,
              COUNTIF(h2.hr BETWEEN 0 AND 23) day1_hours
            FROM `{H}.sofa1_hourly` h1 JOIN `{H}.sofa2_hourly` h2 USING (stay_id, hr)
            WHERE hr BETWEEN 0 AND 23 GROUP BY h1.stay_id),
          g AS (SELECT {GRP} grp, stay_id FROM `{C}` WHERE {PRIM}),
          a AS (SELECT 'ALL' grp, stay_id FROM g UNION ALL SELECT grp, stay_id FROM g)
          SELECT a.grp, COUNT(*) n,
            SUM(m1_resp) m1_resp, SUM(m1_coag) m1_coag, SUM(m1_liver) m1_liver, SUM(m1_cv) m1_cv, SUM(m1_cns) m1_cns, SUM(m1_renal) m1_renal,
            SUM(m2_resp) m2_resp, SUM(m2_hemo) m2_hemo, SUM(m2_liver) m2_liver, SUM(m2_cv) m2_cv, SUM(m2_brain) m2_brain, SUM(m2_kidney) m2_kidney
          FROM a JOIN m USING (stay_id) GROUP BY a.grp""")
        # shares are computed from counts here and output-controlled; the SOFA-only group is merged with Neither.
        cnts = {r["grp"]: {k: int(v) for k, v in r.items() if k.startswith("m")} | {"n": int(r["n"])} for r in ms}
        res["day1_measured_share"] = {g: {k: f(c[k] / c["n"], 3) for k in ("m1_resp", "m2_resp")} | {"n": c["n"]}
                                      for g, c in cnts.items()}
        merged = dict(cnts)
        so, ne = merged.pop("SOFA_ONLY"), merged.pop("NEITHER")
        merged["SOFA2_NEGATIVE"] = {k: so[k] + ne[k] for k in so}
        pairs = [("resp", "m1_resp", "m2_resp"), ("coag", "m1_coag", "m2_hemo"), ("liver", "m1_liver", "m2_liver"),
                 ("cv", "m1_cv", "m2_cv"), ("brain", "m1_cns", "m2_brain"), ("kidney", "m1_renal", "m2_kidney")]
        table = {}
        for g, c in merged.items():
            n = c["n"]
            row = {"n": n}
            for comp, a, b in pairs:
                cell = {}
                for lab, k in (("sofa", a), ("sofa2", b)):
                    v = c[k]
                    if 0 < v < MIN_CELL:
                        cell[lab] = "<5 stays"
                    elif 0 < n - v < MIN_CELL:
                        cell[lab] = "all but <5"
                    else:
                        cell[lab] = f(v / n, 3)
                    if isinstance(cell[lab], str):
                        res["suppressed"].append(f"day1_measured {g} {comp} {lab}")
                d = abs(c[a] - c[b])
                if 0 < d < MIN_CELL and not any(isinstance(x, str) for x in cell.values()):
                    cell = {"sofa": "pair", "sofa2": "pair", "both": f(c[b] / n, 3), "note": "scores differ by 1-4 stays"}
                    res["suppressed"].append(f"day1_measured {g} {comp} between-score difference")
                row[comp] = cell
            table[g] = row
        res["day1_measured_table"] = table
        # ---- urine-output documentation -------------------------------------------------------
        uo = run(f"""
          WITH g AS (SELECT {GRP} grp, stay_id FROM `{C}` WHERE {PRIM} AND fw_s2_kidney >= 1 AND fw_s2_kidney_cr = 0),
          fire AS (
            SELECT z.stay_id, ARRAY_AGG(STRUCT(s.endtime, IF(s.uo_gt24_or_anuria_flag = 1, 24, IF(s.uo_gt12_flag = 1, 12, 6)) AS len)
                   ORDER BY z.hr LIMIT 1)[OFFSET(0)] q
            FROM `{H}.hz_scores` z
            JOIN `{H}.sofa2_withintime` s USING (stay_id, hr)
            JOIN g USING (stay_id)
            WHERE z.hz = 'W' AND z.hz_key = 1 AND z.evaluate
              AND (s.uo_6_12_flag = 1 OR s.uo_gt12_flag = 1 OR s.uo_gt24_or_anuria_flag = 1)
            GROUP BY z.stay_id),
          vols AS (
            SELECT f.stay_id, f.q.len, u.charttime,
              DATETIME_SUB(f.q.endtime, INTERVAL f.q.len HOUR) ws, f.q.endtime we
            FROM fire f LEFT JOIN `@UPSTREAM@.urine_output` u
              ON u.stay_id = f.stay_id AND u.charttime >= DATETIME_SUB(f.q.endtime, INTERVAL f.q.len HOUR) AND u.charttime < f.q.endtime),
          gaps AS (
            SELECT stay_id, ANY_VALUE(len) len, COUNT(charttime) n_vol,
              MAX(DATETIME_DIFF(COALESCE(nxt, we), prv, MINUTE)) / 60.0 AS max_gap_h
            FROM (SELECT stay_id, len, we, charttime,
                    COALESCE(LAG(charttime) OVER (PARTITION BY stay_id ORDER BY charttime), ws) prv,
                    LEAD(charttime) OVER (PARTITION BY stay_id ORDER BY charttime) nxt
                  FROM vols) GROUP BY stay_id),
          x AS (SELECT g.grp, gaps.* FROM gaps JOIN g USING (stay_id)),
          a AS (SELECT 'ALL_UO_ONLY' grp, * EXCEPT (grp) FROM x UNION ALL SELECT grp, * EXCEPT (grp) FROM x WHERE grp = 'SOFA2_ONLY')
          SELECT grp, COUNT(*) n, APPROX_QUANTILES(n_vol, 4) nvol_q, AVG(IF(n_vol <= 1, 1, 0)) single_volume_share,
                 APPROX_QUANTILES(max_gap_h, 4) max_gap_q, AVG(IF(len = 6, 1, 0)) share_6h_criterion
          FROM a GROUP BY grp""")
        res["uo_documentation"] = {r["grp"]: {
            "n": int(r["n"]), "charted_volumes_median_iqr": [f(x["v"] if isinstance(x, dict) else x, 1) for x in r["nvol_q"]][1:4],
            "single_volume_share": f(r["single_volume_share"], 3),
            "longest_gap_h_median_iqr": [f(x["v"] if isinstance(x, dict) else x, 1) for x in r["max_gap_q"]][1:4],
            "share_first_hit_6_12h_criterion": f(r["share_6h_criterion"], 3)} for r in uo}
        # ---- weight back-applied on day 1 -------------------------------------------------------
        w = run(f"""SELECT COUNT(*) hours, COUNTIF(s.weight_backapplied = 1) backapplied, COUNT(DISTINCT IF(s.weight_backapplied = 1, s.stay_id, NULL)) stays,
                           COUNTIF(s.weight_backapplied = 1 AND s.wt_adm_hr <= 23) backapplied_used_day1,
                           COUNTIF(s.weight_backapplied = 1 AND s.wt_adm_hr > 23) backapplied_excluded_day1
                    FROM `{H}.sofa2_withintime` s JOIN `{C}` c USING (stay_id) WHERE c.{PRIM.replace(' AND NOT ', ' AND NOT c.')} AND s.hr BETWEEN 0 AND 23""")[0]
        res["weight_backapplied_day1"] = {k: int(v) for k, v in w.items()}
        # ---- eMAR coverage by anchor_year_group ---------------------------------------------
        em = run(f"""
          WITH d AS (SELECT stay_id, MAX(delirium_treated_any) del, MAX(dex_treated_any) dex FROM `{H}.sofa2_withintime` GROUP BY 1)
          SELECT yg, COUNT(*) n, SUM(del) n_del, SUM(dex) n_dex FROM (
            SELECT c.anchor_year_group yg, d.del, d.dex FROM `{C}` c JOIN d USING (stay_id)
            WHERE c.first_elig_stay AND NOT c.not_assessable_any
            UNION ALL SELECT CONCAT('era ', c.era), d.del, d.dex FROM `{C}` c JOIN d USING (stay_id)
            WHERE c.first_elig_stay AND NOT c.not_assessable_any)
          GROUP BY 1 ORDER BY 1""")
        res["emar_by_year_group"] = {}
        for r in em:
            n, a, b = int(r["n"]), int(r["n_del"]), int(r["n_dex"])
            res["emar_by_year_group"][r["yg"]] = {
                "n": n, "primary_delirium_treatment": wilson(a, n) if a >= MIN_CELL else "suppressed",
                "dexmedetomidine": wilson(b, n) if b >= MIN_CELL else "suppressed"}
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    (OUT / "descriptives.json").write_text(json.dumps(res, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps({"groups": list(res["table1"]), "uo_documentation": res["uo_documentation"], "suppressed": res["suppressed"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
