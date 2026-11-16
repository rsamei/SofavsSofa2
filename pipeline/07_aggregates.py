"""Step 7: aggregates for the flow table, operating characteristics, weight and data-quality checks.
  1. flow: cohort_base -> adults -> eligible -> first eligible -> not assessable (with deaths) -> primary; first
     suspicion after ICU hours 23 and 48;
  2. operating characteristics of the Sepsis-3 threshold for in-hospital death: sensitivity, specificity (paired
     differences, Wald CI), PPV and NPV (Wilson CI);
  3. weight: day-1 weight distribution, counts < 30 kg and > 300 kg, stays without a usable day-1 weight, and the
     size of the admission-weight item restriction;
  4. SOFA-2 kidney 'measured' flag on day 1 by concordance group;
  5. stays whose last day-1 hour is before hour 23: n, deaths, and day-1 kidney scores raised by pre-ICU hours;
  6. timestamp inconsistencies in the primary set;
  7. SOFA-2-only cases: which component scores higher under SOFA-2 than SOFA at the first qualifying hour.
Counts below 5 are suppressed (0 is reported); distribution summaries are exact median and IQR."""
import json
import math
import sys

from bq_client import BUILD, OUT, PROJECT, BQError, query

V4 = BUILD
C = f"{V4}.cohort"
MIN_CELL = 5
PRIM = "first_elig_stay AND NOT not_assessable_any"
GRP = ("CASE WHEN any_s1_std AND any_s2_std THEN 'BOTH' WHEN any_s1_std THEN 'SOFA_ONLY' "
       "WHEN any_s2_std THEN 'SOFA2_ONLY' ELSE 'NEITHER' END")
SUPPRESSED = []


def run(sql):
    fields, rows, _ = query(PROJECT, sql)
    return [dict(zip(fields, r)) for r in rows]


def cnt(v, label):
    v = int(v)
    if 0 < v < MIN_CELL:
        SUPPRESSED.append(label)
        return "<5"
    return v


def wilson(k, n, z=1.96):
    if n == 0:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(p, 4), round(c - h, 4), round(c + h, 4)]


def paired(b, c, n):
    d = (c - b) / n
    se = math.sqrt((b + c) - (c - b) ** 2 / n) / n
    return [round(d, 4), round(d - 1.96 * se, 4), round(d + 1.96 * se, 4)]


def main():
    res = {}
    try:
        # 1. flow ---------------------------------------------------------------------------
        f = run(f"""
          SELECT
            (SELECT COUNT(*) FROM `@UPSTREAM@.cohort_base`) cohort_base,
            (SELECT COUNT(*) FROM `@UPSTREAM@.cohort_base` c JOIN `@UPSTREAM@.icustay_times` t
               USING (stay_id) WHERE c.age >= 18) adults,
            (SELECT COUNT(*) FROM `{C}`) eligible,
            (SELECT COUNTIF(first_elig_stay) FROM `{C}`) first_eligible,
            (SELECT COUNTIF(first_elig_stay AND not_assessable_any) FROM `{C}`) not_assessable,
            (SELECT COUNTIF(first_elig_stay AND not_assessable_any AND hospital_expire_flag = 1) FROM `{C}`) not_assessable_deaths,
            (SELECT COUNTIF({PRIM}) FROM `{C}`) primary_n,
            (SELECT COUNTIF({PRIM} AND first_susp_hr > 23) FROM `{C}`) first_susp_after_hr23,
            (SELECT COUNTIF({PRIM} AND first_susp_hr >= 48) FROM `{C}`) first_susp_hr48_or_later,
            (SELECT COUNTIF(not_assessable_any) FROM `{C}`) all_stays_not_assessable""")[0]
        res["flow"] = {k: cnt(v, f"flow {k}") for k, v in f.items()}
        if isinstance(res["flow"]["adults"], int) and isinstance(res["flow"]["eligible"], int):
            res["flow"]["adults_without_overlapping_suspicion"] = res["flow"]["adults"] - res["flow"]["eligible"]
        # 2. operating characteristics -------------------------------------------------------
        o = run(f"""
          SELECT hospital_expire_flag y,
            COUNTIF(any_s1_std AND any_s2_std) n11, COUNTIF(any_s1_std AND NOT any_s2_std) n10,
            COUNTIF(NOT any_s1_std AND any_s2_std) n01, COUNTIF(NOT any_s1_std AND NOT any_s2_std) n00
          FROM `{C}` WHERE {PRIM} GROUP BY 1""")
        t = {int(r["y"]): {k: int(r[k]) for k in ("n11", "n10", "n01", "n00")} for r in o}
        dth, srv = t[1], t[0]
        D, S = sum(dth.values()), sum(srv.values())
        pos1_d, pos2_d = dth["n11"] + dth["n10"], dth["n11"] + dth["n01"]
        neg1_s, neg2_s = srv["n00"] + srv["n01"], srv["n00"] + srv["n10"]
        pos1, pos2 = pos1_d + srv["n11"] + srv["n10"], pos2_d + srv["n11"] + srv["n01"]
        neg1, neg2 = D + S - pos1, D + S - pos2
        res["operating_characteristics_hosp_death"] = {
            "deaths": D, "survivors": S,
            "sensitivity_sofa": wilson(pos1_d, D), "sensitivity_sofa2": wilson(pos2_d, D),
            "sensitivity_diff_sofa2_minus_sofa": paired(dth["n10"], dth["n01"], D),
            "specificity_sofa": wilson(neg1_s, S), "specificity_sofa2": wilson(neg2_s, S),
            "specificity_diff_sofa2_minus_sofa": paired(srv["n01"], srv["n10"], S),
            "ppv_sofa": wilson(pos1_d, pos1), "ppv_sofa2": wilson(pos2_d, pos2),
            "npv_sofa": wilson(neg1_s, neg1), "npv_sofa2": wilson(neg2_s, neg2),
            "note": "threshold Sepsis-3 (score rise >= 2 from baseline 0, any window); net benefit not assessed; post hoc"}
        # 3. weight ------------------------------------------------------------------------
        w = run(f"""
          WITH p AS (SELECT stay_id, {GRP} grp FROM `{C}` WHERE {PRIM}),
          h AS (
            SELECT s.stay_id, s.hr, s.weight_own_kg, s.weight_ba_kg, s.weight_backapplied, s.wt_adm_hr
            FROM `{V4}.sofa2_withintime` s JOIN p USING (stay_id) WHERE s.hr BETWEEN 0 AND 23),
          anyw AS (
            SELECT DISTINCT ce.stay_id FROM `physionet-data.mimiciv_3_1_icu.chartevents` ce JOIN p USING (stay_id)
            WHERE ce.itemid IN (226512, 224639) AND ce.valuenum > 0),
          per AS (
            SELECT h.stay_id,
              ARRAY_AGG(IF(h.weight_own_kg IS NOT NULL, h.weight_own_kg,
                           IF(h.weight_backapplied = 1 AND h.wt_adm_hr <= 23, h.weight_ba_kg, NULL)) IGNORE NULLS
                        ORDER BY h.hr DESC LIMIT 1)[SAFE_OFFSET(0)] d1_weight,
              COUNTIF(h.weight_own_kg IS NULL AND a.stay_id IS NOT NULL) hours_v3_rule_backapplied,
              COUNTIF(h.weight_backapplied = 1) hours_v4_backapplied,
              COUNTIF(h.weight_backapplied = 1 AND h.wt_adm_hr <= 23) hours_v4_backapplied_used
            FROM h LEFT JOIN anyw a USING (stay_id) GROUP BY h.stay_id)
          SELECT p.grp, COUNT(*) n,
            COUNTIF(per.d1_weight IS NULL) no_usable_d1_weight,
            COUNTIF(anyw.stay_id IS NULL) no_weight_charted,
            COUNTIF(per.d1_weight < 30) under_30kg, COUNTIF(per.d1_weight > 300) over_300kg,
            SUM(per.hours_v3_rule_backapplied) hours_v3_rule, COUNTIF(per.hours_v3_rule_backapplied > 0) stays_v3_rule,
            SUM(per.hours_v4_backapplied) hours_v4, COUNTIF(per.hours_v4_backapplied > 0) stays_v4,
            SUM(per.hours_v4_backapplied_used) hours_v4_used
          FROM p JOIN per USING (stay_id) LEFT JOIN anyw USING (stay_id) GROUP BY 1""")
        tot = {}
        res["weight"] = {"by_group": {}}
        for r in w:
            # < 30 kg and > 300 kg only for all stays: by group they would make a suppressed cell recoverable
            row = {k: cnt(v, f"weight {r['grp']} {k}") for k, v in r.items() if k not in ("grp", "under_30kg", "over_300kg")}
            res["weight"]["by_group"][r["grp"]] = row
            for k, v in r.items():
                if k != "grp":
                    tot[k] = tot.get(k, 0) + int(v)
        res["weight"]["all"] = {k: cnt(v, f"weight all {k}") for k, v in tot.items()}
        q = run(f"""
          WITH p AS (SELECT stay_id FROM `{C}` WHERE {PRIM}),
          per AS (
            SELECT s.stay_id, ARRAY_AGG(IF(s.weight_own_kg IS NOT NULL, s.weight_own_kg,
                IF(s.weight_backapplied = 1 AND s.wt_adm_hr <= 23, s.weight_ba_kg, NULL)) IGNORE NULLS
                ORDER BY s.hr DESC LIMIT 1)[SAFE_OFFSET(0)] wt
            FROM `{V4}.sofa2_withintime` s JOIN p USING (stay_id) WHERE s.hr BETWEEN 0 AND 23 GROUP BY 1)
          SELECT PERCENTILE_CONT(wt, 0.25) OVER () q1, PERCENTILE_CONT(wt, 0.5) OVER () q2, PERCENTILE_CONT(wt, 0.75) OVER () q3
          FROM per WHERE wt IS NOT NULL LIMIT 1""")[0]
        res["weight"]["d1_weight_kg_median_iqr"] = [round(float(q["q2"]), 1), round(float(q["q1"]), 1), round(float(q["q3"]), 1)]
        # 4. corrected kidney measured flag and complete case --------------------------------
        k = run(f"""
          WITH p AS (SELECT stay_id, {GRP} grp, hospital_expire_flag y, complete_case FROM `{C}` WHERE {PRIM}),
          m AS (
            SELECT s.stay_id,
              MAX(CAST(s.creatinine_max IS NOT NULL OR COALESCE(s.on_rrt_any, 0) = 1
                OR (s.uo_n_6h >= 1 AND (COALESCE(s.weight_own_kg, 0) > 0
                    OR (s.weight_backapplied = 1 AND s.wt_adm_hr <= 23 AND s.weight_ba_kg > 0))) AS INT64)) m2k_fixed,
              MAX(h2.m2_kidney) m2k_old,
              MAX(h1.m1_resp) = 1 AND MAX(h1.m1_coag) = 1 AND MAX(h1.m1_liver) = 1 AND MAX(h1.m1_cv) = 1
                AND MAX(h1.m1_cns) = 1 AND MAX(h1.m1_renal) = 1 AS s1_complete,
              MAX(h2.m2_resp) = 1 AND MAX(h2.m2_hemo) = 1 AND MAX(h2.m2_liver) = 1 AND MAX(h2.m2_cv) = 1
                AND MAX(h2.m2_brain) = 1 AS s2_complete_nokidney
            FROM `{V4}.sofa2_withintime` s
            JOIN `{V4}.sofa2_hourly` h2 USING (stay_id, hr) JOIN `{V4}.sofa1_hourly` h1 USING (stay_id, hr)
            JOIN p USING (stay_id)
            WHERE s.hr BETWEEN 0 AND 23 GROUP BY s.stay_id)
          SELECT p.grp, COUNT(*) n, AVG(m.m2k_old) m2_kidney_old, AVG(m.m2k_fixed) m2_kidney_fixed,
            COUNTIF(p.complete_case) cc_old, COUNTIF(p.complete_case AND p.y = 1) cc_old_deaths,
            COUNTIF(m.s1_complete AND m.s2_complete_nokidney AND m.m2k_fixed = 1) cc_fixed,
            COUNTIF(m.s1_complete AND m.s2_complete_nokidney AND m.m2k_fixed = 1 AND p.y = 1) cc_fixed_deaths
          FROM p JOIN m USING (stay_id) GROUP BY 1""")
        res["kidney_measured_day1"] = {}
        cc = {"cc_old": 0, "cc_old_deaths": 0, "cc_fixed": 0, "cc_fixed_deaths": 0}
        for r in k:
            res["kidney_measured_day1"][r["grp"]] = {"n": int(r["n"]), "m2_kidney_old": round(float(r["m2_kidney_old"]), 3),
                                                     "m2_kidney_fixed": round(float(r["m2_kidney_fixed"]), 3)}
            for key in cc:
                cc[key] += int(r[key])
        res["complete_case_set"] = {k2: cnt(v, f"complete case {k2}") for k2, v in cc.items()}
        # 5. last day-1 hour before 23 -------------------------------------------------------
        e = run(f"""
          WITH p AS (SELECT stay_id, d1_last_hr, hospital_expire_flag y FROM `{C}` WHERE {PRIM} AND d1_last_hr < 23),
          k AS (
            SELECT p.stay_id, MAX(IF(s.hr = p.d1_last_hr, s.kidney_pre24, NULL)) pre24,
                   MAX(IF(s.hr BETWEEN 0 AND p.d1_last_hr, s.kidney_b, NULL)) icu_max
            FROM p JOIN `{V4}.sofa2_hourly` s USING (stay_id) GROUP BY 1)
          SELECT COUNT(*) n, SUM(p.y) deaths, COUNTIF(COALESCE(k.pre24, 0) > COALESCE(k.icu_max, 0)) kidney_raised_by_pre_icu
          FROM p JOIN k USING (stay_id)""")[0]
        res["d1_last_hr_before_23"] = {k2: cnt(v, f"short day 1 {k2}") for k2, v in e.items()}
        # 6. timestamps ------------------------------------------------------------------------
        ts = run(f"""
          SELECT COUNTIF(c.icu_intime < a.admittime) intime_before_admit, COUNTIF(c.icu_intime > a.dischtime) intime_after_disch,
            COUNTIF(a.deathtime < c.icu_intime) death_before_icu, COUNTIF(a.deathtime > a.dischtime) death_after_disch,
            COUNTIF(a.dischtime < a.admittime) negative_hosp_los
          FROM `{C}` c JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a USING (hadm_id) WHERE c.{PRIM.replace(' AND NOT ', ' AND NOT c.')}""")[0]
        res["timestamps_primary"] = {k2: cnt(v, f"timestamps {k2}") for k2, v in ts.items()}
        # 7. component carrying the qualifying point in SOFA-2-only cases ---------------------
        comp = run(f"""
          WITH p AS (
            SELECT stay_id,
              CASE WHEN NOT any_s2_kcr THEN 'uo_dependent' ELSE 'not_uo_dependent' END cls
            FROM `{C}` WHERE {PRIM} AND any_s2_std AND NOT any_s1_std),
          o AS (
            SELECT z.* FROM `{V4}.hz_scores` z JOIN p USING (stay_id)
            WHERE z.hz = 'W' AND z.evaluate AND z.s2_std >= 2
            QUALIFY ROW_NUMBER() OVER (PARTITION BY z.stay_id ORDER BY z.endtime, z.hz_key) = 1),
          x AS (
            SELECT p.cls,
              o.s2_resp24 > o.s1_resp24 resp, o.s2_hemo24 > o.s1_coag24 hemo, o.s2_liver24 > o.s1_liver24 liver,
              o.s2_cv24 > o.s1_cv24 cv, o.s2_brain24 > o.s1_cns24 brain, o.s2_kidney24 > o.s1_renal24 kidney
            FROM p JOIN o USING (stay_id))
          SELECT cls, COUNT(*) n,
            COUNTIF(resp) resp, COUNTIF(hemo) hemo, COUNTIF(liver) liver, COUNTIF(cv) cv, COUNTIF(brain) brain, COUNTIF(kidney) kidney,
            COUNTIF(resp AND NOT (hemo OR liver OR cv OR brain OR kidney)) only_resp,
            COUNTIF(brain AND NOT (resp OR hemo OR liver OR cv OR kidney)) only_brain,
            COUNTIF(cv AND NOT (resp OR hemo OR liver OR brain OR kidney)) only_cv,
            COUNTIF(kidney AND NOT (resp OR hemo OR liver OR cv OR brain)) only_kidney,
            COUNTIF(NOT (resp OR hemo OR liver OR cv OR brain OR kidney)) none_higher
          FROM x GROUP BY 1""")
        res["sofa2_only_components_at_first_qualifying_hour"] = {}
        for r in comp:
            n = int(r["n"])
            row = {"n": cnt(n, f"components {r['cls']} n")}
            for key in ("resp", "hemo", "liver", "cv", "brain", "kidney", "only_resp", "only_brain", "only_cv",
                        "only_kidney", "none_higher"):
                v = cnt(r[key], f"components {r['cls']} {key}")
                row[key] = [v, round(int(r[key]) / n, 3)] if isinstance(v, int) and v > 0 else v
            res["sofa2_only_components_at_first_qualifying_hour"][r["cls"]] = row
        res["sofa2_only_components_note"] = ("component scores higher under SOFA-2 than under SOFA (24-h maxima) at the "
                                             "first evaluated hour with SOFA-2 >= 2; several can be higher; post hoc")
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    res["suppressed"] = SUPPRESSED
    (OUT / "aggregates.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"flow": res["flow"], "suppressed": SUPPRESSED}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
