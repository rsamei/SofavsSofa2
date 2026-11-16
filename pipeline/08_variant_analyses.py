"""Step 8: variant j and further sensitivity analyses. Builds BUILD_DATASET.variant_stay (one row per eligible stay):
  - variant j: SOFA-2 respiration from PaO2/FiO2 only (ECMO scores 4), no SpO2/FiO2 fallback, and j-ref, the same
    construction keeping the SpO2/FiO2-based respiration; both use a 24-hour respiratory maximum over ICU hours
    (hr >= 0), so j - j-ref isolates the input;
  - any-window Sepsis-3 flags for j, j-ref and SOFA-2 at thresholds 3, 4 and 5.
Then reports (only aggregates leave BigQuery; counts below 5 are suppressed):
  1. paired DeLong AUROC contrasts (as in step 4) for j and j-ref, and SOFA-2 vs SOFA in the full-24-hour-day-1 and
     suspicion-by-hour-23 sets and in each era (with the two matched-GCS contrasts);
  2. Sepsis-3 rates for j and j-ref; matched positivity: the SOFA-2 threshold whose rate is closest to SOFA's at 2;
  3. paired baseline-0 vs individual-baseline classification on the stays with an individual baseline;
  4. kidney reference (KDIGO creatinine criterion or RRT during the ICU stay) for SOFA-2-only stays by urine
     dependence and for Neither;
  5. overlap of the admission-matched and stay-matched primary sets."""
import json
import math
import sys

from bq_client import BUILD, OUT, PROJECT, BQError, execute, query

V4 = BUILD
MIN_CELL = 5
PRIM = "first_elig_stay AND NOT not_assessable_any"
SUPPRESSED = []

BUILD_SQL = f"""
CREATE OR REPLACE TABLE `{V4}.variant_stay` AS
WITH r AS (
  SELECT stay_id, hr,
    CASE WHEN COALESCE(ecmo_any, 0) = 1 THEN 4
         WHEN pao2fio2ratio_min IS NOT NULL THEN
           CASE WHEN advanced_vent_any = 1 AND pao2fio2ratio_min <= 75 THEN 4
                WHEN advanced_vent_any = 1 AND pao2fio2ratio_min <= 150 THEN 3
                WHEN pao2fio2ratio_min <= 225 THEN 2
                WHEN pao2fio2ratio_min <= 300 THEN 1 ELSE 0 END
         ELSE NULL END AS resp_pf,
    respiration AS resp_ref
  FROM `{V4}.sofa2_withintime`
),
r24 AS (
  SELECT stay_id, hr,
    COALESCE(MAX(resp_pf) OVER w, 0) AS rpf24, COALESCE(MAX(resp_ref) OVER w, 0) AS rref24
  FROM r WINDOW w AS (PARTITION BY stay_id ORDER BY hr ROWS BETWEEN 23 PRECEDING AND CURRENT ROW)
),
z AS (
  SELECT z.stay_id, z.hz, z.hr, z.evaluate, z.s2_std,
    z.s2_std - z.s2_resp24 + r24.rpf24 AS s2_j,
    z.s2_std - z.s2_resp24 + r24.rref24 AS s2_jref
  FROM `{V4}.hz_scores` z JOIN r24 USING (stay_id, hr)
),
d1 AS (
  SELECT z.stay_id, z.s2_j AS d1_s2_j, z.s2_jref AS d1_s2_jref
  FROM z JOIN `{V4}.cohort` c ON c.stay_id = z.stay_id AND z.hz = 'D1' AND z.hr = c.d1_last_hr
),
w AS (
  SELECT stay_id,
    LOGICAL_OR(s2_j >= 2) AS any_s2_j, LOGICAL_OR(s2_jref >= 2) AS any_s2_jref,
    LOGICAL_OR(s2_std >= 3) AS any_s2_t3, LOGICAL_OR(s2_std >= 4) AS any_s2_t4, LOGICAL_OR(s2_std >= 5) AS any_s2_t5
  FROM z WHERE hz = 'W' AND evaluate GROUP BY stay_id
)
SELECT c.*, d1.d1_s2_j, d1.d1_s2_jref,
  COALESCE(w.any_s2_j, FALSE) AS any_s2_j, COALESCE(w.any_s2_jref, FALSE) AS any_s2_jref,
  COALESCE(w.any_s2_t3, FALSE) AS any_s2_t3, COALESCE(w.any_s2_t4, FALSE) AS any_s2_t4,
  COALESCE(w.any_s2_t5, FALSE) AS any_s2_t5
FROM `{V4}.cohort` c LEFT JOIN d1 USING (stay_id) LEFT JOIN w USING (stay_id)
"""

T = f"{V4}.variant_stay"
SYSTEMS = ["s1_std", "s1_gcsalt", "s2_std", "s2_gcs15", "s2_j", "s2_jref"]
SETS = {
    "primary": PRIM,
    "full_day1": PRIM + " AND d1_last_hr = 23",
    "susp_by_hr23": PRIM + " AND first_susp_hr <= 23",
    "era_2008_2013": PRIM + " AND era = '2008-2013'",
    "era_2014_2016": PRIM + " AND era = '2014-2016'",
    "era_2017_2022": PRIM + " AND era = '2017-2022'",
}
CONTRASTS = {
    "s2_minus_s1": {"s2_std": 1, "s1_std": -1},
    "matched_gcs_a": {"s2_std": 1, "s1_gcsalt": -1},
    "matched_gcs_b": {"s2_gcs15": 1, "s1_std": -1},
    "j_minus_s1": {"s2_j": 1, "s1_std": -1},
    "jref_minus_s1": {"s2_jref": 1, "s1_std": -1},
    "j_minus_jref": {"s2_j": 1, "s2_jref": -1},
    "s2_minus_jref": {"s2_std": 1, "s2_jref": -1},
}


def delong_sql():
    sets = ",\n      ".join(f"STRUCT('{k}' AS an, ({v}) AS inc)" for k, v in SETS.items())
    systems = ",\n      ".join(f"STRUCT('{s}' AS sys, d1_{s} AS x)" for s in SYSTEMS)
    pv_cols = ",\n    ".join(f"MAX(IF(sys = '{s}', v, NULL)) AS v_{s}" for s in SYSTEMS)
    defs = []
    for name, coef in [(f"auc_{s}", {s: 1}) for s in SYSTEMS] + list(CONTRASTS.items()):
        coefs = ", ".join(f"{float(coef.get(s, 0))} AS k_{s}" for s in SYSTEMS)
        defs.append(f"STRUCT('{name}' AS contrast, {coefs})")
    combo = " + ".join(f"d.k_{s} * p.v_{s}" for s in SYSTEMS)
    return f"""
WITH long AS (
  SELECT st.an, sy.sys, stay_id, hospital_expire_flag AS y, sy.x
  FROM `{T}`, UNNEST([
      {sets}
    ]) AS st, UNNEST([
      {systems}
    ]) AS sy
  WHERE st.inc
),
cc AS (SELECT an, stay_id FROM long GROUP BY an, stay_id HAVING COUNTIF(x IS NULL) = 0),
r AS (
  SELECT l.an, l.sys, l.stay_id, l.y, l.x,
    RANK() OVER (PARTITION BY l.an, l.sys ORDER BY l.x) + (COUNT(*) OVER (PARTITION BY l.an, l.sys, l.x) - 1) / 2.0 AS pr,
    RANK() OVER (PARTITION BY l.an, l.sys, l.y ORDER BY l.x) + (COUNT(*) OVER (PARTITION BY l.an, l.sys, l.y, l.x) - 1) / 2.0 AS cr,
    COUNTIF(l.y = 1) OVER (PARTITION BY l.an, l.sys) AS m, COUNTIF(l.y = 0) OVER (PARTITION BY l.an, l.sys) AS n
  FROM long l JOIN cc USING (an, stay_id)
),
pl AS (SELECT an, sys, stay_id, y, IF(y = 1, (pr - cr) / n, 1 - (pr - cr) / m) AS v FROM r),
pv AS (SELECT an, stay_id, ANY_VALUE(y) AS y, {pv_cols} FROM pl GROUP BY an, stay_id),
defs AS (SELECT * FROM UNNEST([
      {",".join(defs)}
    ])),
cv AS (SELECT p.an, d.contrast, p.y, {combo} AS c FROM pv p CROSS JOIN defs d)
SELECT an, contrast, COUNTIF(y = 1) AS n_pos, COUNTIF(y = 0) AS n_neg,
  AVG(IF(y = 1, c, NULL)) AS est,
  SQRT(VAR_SAMP(IF(y = 1, c, NULL)) / COUNTIF(y = 1) + VAR_SAMP(IF(y = 0, c, NULL)) / COUNTIF(y = 0)) AS se
FROM cv GROUP BY an, contrast ORDER BY an, contrast
"""


def run(sql):
    fields, rows, _ = query(PROJECT, sql)
    return [dict(zip(fields, r)) for r in rows]


def cnt(v, label):
    v = int(v)
    if 0 < v < MIN_CELL:
        SUPPRESSED.append(label)
        return "<5"
    return v


def wilson(k, n, z=1.96, nd=6):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(p, nd), max(0.0, round(c - h, nd)), round(c + h, nd)]


def kappa(n11, n10, n01, n00):
    n = n11 + n10 + n01 + n00
    po = (n11 + n00) / n
    pe = ((n11 + n10) * (n11 + n01) + (n00 + n10) * (n00 + n01)) / n / n
    return round((po - pe) / (1 - pe), 3)


def mcnemar(b, c):
    return round((abs(b - c) - 1) ** 2 / (b + c), 1) if b + c else None


def discordant(b, c, label):
    """Discordant counts and McNemar. With both marginal rates reported, b - c is known, so when either
    count is 1-4 both counts and the statistic are withheld."""
    b, c = int(b), int(c)
    if 0 < b < MIN_CELL or 0 < c < MIN_CELL:
        SUPPRESSED.append(f"{label} discordant pair and McNemar")
        return {"sofa_only": "withheld", "sofa2_only": "withheld", "mcnemar_chi2_cc": "withheld",
                "one_discordant_count_below_5": True}
    return {"sofa_only": b, "sofa2_only": c, "mcnemar_chi2_cc": mcnemar(b, c)}


def main():
    res = {}
    try:
        execute(PROJECT, BUILD_SQL)
        # 1. DeLong ----------------------------------------------------------------------------
        rows = []
        for r in run(delong_sql()):
            est, se = float(r["est"]), float(r["se"])
            rows.append({"set": r["an"], "contrast": r["contrast"], "n_pos": int(r["n_pos"]), "n_neg": int(r["n_neg"]),
                         "est": round(est, 6), "lo": round(est - 1.96 * se, 6), "hi": round(est + 1.96 * se, 6)})
        res["discrimination"] = rows
        res["event_floor_ok"] = all(r["n_pos"] >= 100 and r["n_neg"] >= 100 for r in rows)
        # 2. classification: j, j-ref, thresholds ------------------------------------------------
        c = run(f"""SELECT COUNT(*) n, COUNTIF(any_s1_std) s1, COUNTIF(any_s2_std) s2, COUNTIF(any_s2_j) j, COUNTIF(any_s2_jref) jref,
                           COUNTIF(any_s2_t3) t3, COUNTIF(any_s2_t4) t4, COUNTIF(any_s2_t5) t5,
                           COUNTIF(any_s2_j AND NOT any_s1_std) j_only, COUNTIF(any_s2_jref AND NOT any_s1_std) jref_only
                    FROM `{T}` WHERE {PRIM}""")[0]
        n = int(c["n"])
        res["classification"] = {"n": n}
        for k in ("s1", "s2", "j", "jref", "t3", "t4", "t5"):
            res["classification"]["rate_" + k] = wilson(int(c[k]), n)
        res["classification"]["sofa2_only_count_j"] = int(c["j_only"])
        res["classification"]["sofa2_only_count_jref"] = int(c["jref_only"])
        s1 = int(c["s1"])
        best = min((3, 4, 5), key=lambda t: abs(int(c[f"t{t}"]) - s1))
        g = run(f"""SELECT COUNTIF(any_s1_std AND any_s2_t{best}) n11, COUNTIF(any_s1_std AND NOT any_s2_t{best}) n10,
                           COUNTIF(NOT any_s1_std AND any_s2_t{best}) n01, COUNTIF(NOT any_s1_std AND NOT any_s2_t{best}) n00,
                           SUM(IF(any_s1_std AND NOT any_s2_t{best}, hospital_expire_flag, 0)) d10,
                           SUM(IF(NOT any_s1_std AND any_s2_t{best}, hospital_expire_flag, 0)) d01
                    FROM `{T}` WHERE {PRIM}""")[0]
        n11, n10, n01, n00 = (int(g[k]) for k in ("n11", "n10", "n01", "n00"))
        mp = {"threshold": best, "rate_sofa2": wilson(n11 + n01, n), "rate_sofa": wilson(n11 + n10, n),
              "kappa": kappa(n11, n10, n01, n00), **discordant(n10, n01, "mp")}
        for lab, d, nn in (("sofa_only", int(g["d10"]), n10), ("sofa2_only", int(g["d01"]), n01)):
            if d >= MIN_CELL and nn - d >= MIN_CELL:
                mp[f"{lab}_hosp_mortality"] = wilson(d, nn)
            else:
                SUPPRESSED.append(f"mp {lab} deaths")
        res["matched_positivity"] = mp
        # 3. individual baseline, paired on the same stays ----------------------------------------
        b = run(f"""SELECT COUNT(*) n,
            COUNTIF(any_s1_std) s1_b0, COUNTIF(ind_s1_std) s1_ind, COUNTIF(any_s2_std) s2_b0, COUNTIF(ind_s2_std) s2_ind,
            COUNTIF(any_s1_std AND NOT any_s2_std) b0_sofa_only, COUNTIF(any_s2_std AND NOT any_s1_std) b0_sofa2_only,
            COUNTIF(ind_s1_std AND NOT ind_s2_std) ind_sofa_only, COUNTIF(ind_s2_std AND NOT ind_s1_std) ind_sofa2_only,
            COUNTIF(any_s1_std AND NOT ind_s1_std) s1_lost, COUNTIF(any_s2_std AND NOT ind_s2_std) s2_lost,
            AVG(ind_base_s1) base_s1, AVG(ind_base_s2) base_s2
          FROM `{T}` WHERE {PRIM} AND ind_baseline_available""")[0]
        nb = int(b["n"])
        res["individual_baseline_paired"] = {
            "n": nb,
            "baseline0": {"rate_sofa": wilson(int(b["s1_b0"]), nb), "rate_sofa2": wilson(int(b["s2_b0"]), nb),
                          **discordant(b["b0_sofa_only"], b["b0_sofa2_only"], "ib b0")},
            "individual": {"rate_sofa": wilson(int(b["s1_ind"]), nb), "rate_sofa2": wilson(int(b["s2_ind"]), nb),
                           **discordant(b["ind_sofa_only"], b["ind_sofa2_only"], "ib ind")},
            "lost_with_individual_baseline": {"sofa": cnt(b["s1_lost"], "ib s1 lost"), "sofa2": cnt(b["s2_lost"], "ib s2 lost")},
            "mean_baseline": {"sofa": round(float(b["base_s1"]), 2), "sofa2": round(float(b["base_s2"]), 2)},
            "note": "selected subset: at least 24 h of ICU data before the window"}
        # 4. kidney reference ---------------------------------------------------------------------
        k = run(f"""
          WITH p AS (
            SELECT stay_id, hadm_id, icu_intime, icu_outtime,
              CASE WHEN any_s1_std AND any_s2_std THEN 'BOTH' WHEN any_s1_std THEN 'SOFA_ONLY'
                   WHEN NOT any_s2_std THEN 'NEITHER'
                   WHEN NOT any_s2_kcr THEN 'SOFA2_ONLY_uo_dependent' ELSE 'SOFA2_ONLY_not_uo' END AS grp
            FROM `{T}` WHERE {PRIM}),
          cr AS (
            SELECT l.hadm_id, l.charttime, l.valuenum,
              MIN(l.valuenum) OVER (PARTITION BY l.hadm_id ORDER BY UNIX_SECONDS(TIMESTAMP(l.charttime))
                                    RANGE BETWEEN 172800 PRECEDING AND 1 PRECEDING) AS low48,
              MIN(l.valuenum) OVER (PARTITION BY l.hadm_id ORDER BY UNIX_SECONDS(TIMESTAMP(l.charttime))
                                    RANGE BETWEEN 604800 PRECEDING AND 1 PRECEDING) AS low7d
            FROM `physionet-data.mimiciv_3_1_hosp.labevents` l
            WHERE l.itemid = 50912 AND l.valuenum > 0 AND l.valuenum < 150
              AND l.hadm_id IN (SELECT hadm_id FROM p)),
          kd AS (
            SELECT p.stay_id, LOGICAL_OR(c.valuenum >= c.low48 + 0.3 OR c.valuenum >= 1.5 * c.low7d) AS kdigo_cr
            FROM p JOIN cr c ON c.hadm_id = p.hadm_id AND c.charttime BETWEEN p.icu_intime AND p.icu_outtime
            GROUP BY p.stay_id),
          rr AS (
            SELECT DISTINCT p.stay_id FROM p JOIN `@UPSTREAM@.rrt` r
              ON r.stay_id = p.stay_id AND (r.dialysis_active = 1 OR r.dialysis_present = 1)
             AND r.charttime BETWEEN p.icu_intime AND p.icu_outtime)
          SELECT p.grp, COUNT(*) n, COUNTIF(COALESCE(kd.kdigo_cr, FALSE)) kdigo, COUNTIF(rr.stay_id IS NOT NULL) rrt,
                 COUNTIF(COALESCE(kd.kdigo_cr, FALSE) OR rr.stay_id IS NOT NULL) either
          FROM p LEFT JOIN kd USING (stay_id) LEFT JOIN rr USING (stay_id) GROUP BY 1""")
        res["kidney_reference"] = {}
        for r in k:
            nn = int(r["n"])
            row = {"n": nn}
            for key in ("kdigo", "rrt", "either"):
                v = int(r[key])
                if 0 < v < MIN_CELL or 0 < nn - v < MIN_CELL:
                    row[key] = "<5" if v < MIN_CELL else "suppressed"
                    SUPPRESSED.append(f"kidney {r['grp']} {key}")
                else:
                    row[key] = [v] + wilson(v, nn)
            res["kidney_reference"][r["grp"]] = row
        res["kidney_reference_note"] = ("KDIGO creatinine: rise >= 0.3 mg/dL over the lowest value in the preceding 48 h or "
                                        ">= 1.5 x the lowest in the preceding 7 days (same admission), during the ICU stay; RRT from "
                                        "2525.rrt during the ICU stay; exploratory")
        # 5. overlap and unrounded urine shares -------------------------------------------------
        o = run(f"""SELECT COUNTIF({PRIM} AND first_elig_stay_old AND NOT not_assessable_old) both_primary,
                           COUNTIF({PRIM} AND NOT (first_elig_stay_old AND NOT not_assessable_old)) new_only,
                           COUNTIF(first_elig_stay_old AND NOT not_assessable_old AND NOT ({PRIM})) old_only
                    FROM `{T}`""")[0]
        res["primary_set_overlap"] = {k2: int(v) for k2, v in o.items()}
        u = run(f"""SELECT COUNT(*) n, COUNTIF(NOT any_s2_kcr) c, COUNTIF(NOT any_s2_g) g, COUNTIF(one_volume_s2) v
                    FROM `{T}` WHERE {PRIM} AND any_s2_std AND NOT any_s1_std""")[0]
        un = int(u["n"])
        res["uo_share_unrounded"] = {k2: wilson(int(u[k2]), un, nd=6) for k2 in ("c", "g", "v")}
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    res["suppressed"] = SUPPRESSED
    (OUT / "variant_analyses.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"event_floor_ok": res["event_floor_ok"], "matched_threshold": res["matched_positivity"]["threshold"],
                      "suppressed": SUPPRESSED}))
    return 0 if res["event_floor_ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
