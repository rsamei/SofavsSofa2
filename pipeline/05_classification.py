"""Step 5: Sepsis-3 classification (any suspicion window, baseline 0) in the primary set.
  - agreement: 2x2, rates (Wilson 95% CI), paired rate difference, observed/positive/negative agreement, Cohen kappa,
    PABAK (Poisson-bootstrap 95% CI, 1,000 replicates computed in BigQuery), McNemar with continuity correction;
  - concordance groups: n, hospital and ICU deaths, mortality (Wilson CI), crude OR (Woolf CI) and Mantel-Haenszel OR
    vs Neither (strata: age decade [top >= 80] x sex x cardiac-vascular ICU; Robins-Breslow-Greenland CI) with
    Breslow-Day homogeneity; ICU deaths are reported for SOFA-only and Neither merged as SOFA-2 negative, which is
    also the reference of the ICU odds ratios;
  - variants in the first suspicion window; urine-output dependence of SOFA-2-only stays (variants c, d, g) and
    stays resting on a single charted urine volume, with mortality and crude risk differences vs Neither (Newcombe);
  - era description (rates and SOFA-2 brain-component prevalence); individual baseline; subgroup rates.
Output control: counts < 5 and any complementary cell that would reveal them are suppressed; estimates resting on
< 5 events are omitted."""
import json
import math
import sys

from bq_client import BUILD, OUT, PROJECT, BQError, query

C = f"{BUILD}.cohort"
B = 1000  # bootstrap replicates
MIN_CELL = 5


def wilson(k, n, z=1.96):
    if n == 0:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [p, c - h, c + h]


def chi2_sf(x, df):
    """Survival function of chi-square (regularised upper incomplete gamma, series/continued fraction)."""
    if x <= 0:
        return 1.0
    a, xx = df / 2.0, x / 2.0
    if xx < a + 1:
        s, t, k = 1.0 / a, 1.0 / a, a
        for _ in range(500):
            k += 1
            t *= xx / k
            s += t
        return max(0.0, 1 - s * math.exp(-xx + a * math.log(xx) - math.lgamma(a)))
    b, c, d = xx + 1 - a, 1e300, 1 / (xx + 1 - a)
    h = d
    for i in range(1, 500):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = 1 / d if d != 0 else 1e300
        c = b + an / c
        h *= d * c
    return h * math.exp(-xx + a * math.log(xx) - math.lgamma(a))


def run(sql):
    fields, rows, _ = query(PROJECT, sql)
    return [dict(zip(fields, r)) for r in rows]


def num(v):
    return None if v is None else float(v)


PRIM = "first_elig_stay AND NOT not_assessable_any"
GRP = ("CASE WHEN any_s1_std AND any_s2_std THEN 'BOTH' WHEN any_s1_std THEN 'SOFA_ONLY' "
       "WHEN any_s2_std THEN 'SOFA2_ONLY' ELSE 'NEITHER' END")


def main():
    res = {"suppressed": []}
    try:
        # ---- 2x2 and agreement (primary) -------------------------------------------------
        t = run(f"""SELECT COUNTIF(any_s1_std AND any_s2_std) n11, COUNTIF(any_s1_std AND NOT any_s2_std) n10,
                           COUNTIF(NOT any_s1_std AND any_s2_std) n01, COUNTIF(NOT any_s1_std AND NOT any_s2_std) n00,
                           (SELECT COUNTIF(first_elig_stay AND not_assessable_any) FROM `{C}`) not_assessable
                    FROM `{C}` WHERE {PRIM}""")[0]
        n11, n10, n01, n00 = (int(t[k]) for k in ("n11", "n10", "n01", "n00"))
        n = n11 + n10 + n01 + n00
        po = (n11 + n00) / n
        pe = ((n11 + n10) * (n11 + n01) + (n00 + n10) * (n00 + n01)) / n / n
        boot = run(f"""
          WITH base AS (SELECT stay_id, any_s1_std a, any_s2_std b FROM `{C}` WHERE {PRIM}),
          w AS (
            SELECT r, a, b,
              (SELECT COUNTIF(u > th) FROM UNNEST([0.367879, 0.735759, 0.919699, 0.981012, 0.996340, 0.999406, 0.999917]) th) AS wt
            FROM base, UNNEST(GENERATE_ARRAY(1, {B})) r,
              UNNEST([MOD(ABS(FARM_FINGERPRINT(CONCAT(CAST(stay_id AS STRING), '-', CAST(r AS STRING)))), 1000000000) / 1e9]) u),
          c AS (SELECT r, SUM(IF(a AND b, wt, 0)) n11, SUM(IF(a AND NOT b, wt, 0)) n10,
                       SUM(IF(NOT a AND b, wt, 0)) n01, SUM(IF(NOT a AND NOT b, wt, 0)) n00 FROM w GROUP BY r),
          k AS (SELECT r, (n11 + n00) / (n11 + n10 + n01 + n00) po,
                  ((n11 + n10) * (n11 + n01) + (n00 + n10) * (n00 + n01)) / POW(n11 + n10 + n01 + n00, 2) pe FROM c)
          SELECT APPROX_QUANTILES((po - pe) / (1 - pe), 1000)[OFFSET(25)] k_lo, APPROX_QUANTILES((po - pe) / (1 - pe), 1000)[OFFSET(975)] k_hi,
                 APPROX_QUANTILES(2 * po - 1, 1000)[OFFSET(25)] p_lo, APPROX_QUANTILES(2 * po - 1, 1000)[OFFSET(975)] p_hi FROM k""")[0]
        res["agreement"] = {
            "n": n, "not_assessable_excluded": int(t["not_assessable"]),
            "cells": {"both": n11, "sofa_only": n10, "sofa2_only": n01, "neither": n00},
            "rate_sofa": wilson(n11 + n10, n), "rate_sofa2": wilson(n11 + n01, n),
            "observed_agreement": round(po, 6),
            "positive_agreement": round(2 * n11 / (2 * n11 + n10 + n01), 6),
            "negative_agreement": round(2 * n00 / (2 * n00 + n10 + n01), 6),
            "kappa": [round((po - pe) / (1 - pe), 6), round(float(boot["k_lo"]), 6), round(float(boot["k_hi"]), 6)],
            "pabak": [round(2 * po - 1, 6), round(float(boot["p_lo"]), 6), round(float(boot["p_hi"]), 6)],
            "mcnemar_chi2_cc": round((abs(n10 - n01) - 1) ** 2 / (n10 + n01), 1),
            "mcnemar_p": (lambda p: float(f"{p:.3g}") if p > 0 else "< 1e-300")(chi2_sf((abs(n10 - n01) - 1) ** 2 / (n10 + n01), 1)),
            "paired_rate_difference_sofa2_minus_sofa": (lambda d, se: [round(d, 6), round(d - 1.96 * se, 6), round(d + 1.96 * se, 6)])(
                (n01 - n10) / n, math.sqrt((n10 + n01) - (n01 - n10) ** 2 / n) / n),
            "bootstrap_reps": B,
        }
        # ---- concordance groups, crude and MH ORs, Breslow-Day ----------------------------
        strata = run(f"""
          SELECT {GRP} grp, CAST(FLOOR(LEAST(age, 89) / 10) AS INT64) ab, gender, cvicu,
                 COUNT(*) n, SUM(hospital_expire_flag) dh, SUM(icu_death) di
          FROM `{C}` WHERE {PRIM} GROUP BY 1, 2, 3, 4""")
        groups = {}
        for s in strata:
            g = groups.setdefault(s["grp"], {"n": 0, "dh": 0, "di": 0})
            for k in ("n", "dh", "di"):
                g[k] += int(s[k])
        # merged SOFA-2-negative reference for ICU death (per stratum as well)
        groups["SOFA2_NEGATIVE"] = {k: groups["NEITHER"][k] + groups["SOFA_ONLY"][k] for k in ("n", "dh", "di")}
        merged = {}
        for st in strata:
            if st["grp"] in ("NEITHER", "SOFA_ONLY"):
                key = (st["ab"], st["gender"], st["cvicu"])
                m = merged.setdefault(key, {"grp": "SOFA2_NEGATIVE", "ab": st["ab"], "gender": st["gender"],
                                            "cvicu": st["cvicu"], "n": 0, "dh": 0, "di": 0})
                for k in ("n", "dh", "di"):
                    m[k] += int(st[k])
        strata = list(strata) + list(merged.values())
        res["groups"] = {}
        for g, v in groups.items():
            row = {"n": v["n"]}
            outcomes = {"NEITHER": (("hosp", "dh"),), "SOFA_ONLY": (("hosp", "dh"),)}.get(g, (("hosp", "dh"), ("icu", "di")))
            for o, k in outcomes:
                if 0 < v[k] < MIN_CELL:
                    row[f"{o}_deaths"] = "<5"
                    res["suppressed"].append(f"{g} {o} deaths")
                else:
                    row[f"{o}_deaths"] = v[k]
                    row[f"{o}_mortality"] = wilson(v[k], v["n"])
            res["groups"][g] = row
        res["odds_ratios"] = {}
        comparisons = [(g, "hosp", "dh", "NEITHER", f"{g}|hosp") for g in ("BOTH", "SOFA_ONLY", "SOFA2_ONLY")]
        comparisons += [(g, "icu", "di", "SOFA2_NEGATIVE", f"{g}|icu|vs_sofa2_negative") for g in ("BOTH", "SOFA2_ONLY")]
        comparisons += [("SOFA2_ONLY", "hosp", "dh", "SOFA2_NEGATIVE", "SOFA2_ONLY|hosp|vs_sofa2_negative")]
        res["odds_ratio_reference_note"] = (
            "Hospital-death ORs use Neither as reference (pre-specified). ICU-death ORs use SOFA-2 negative "
            "(SOFA-only + Neither) because the SOFA-only ICU-death count (1-4) is recoverable from the published "
            "ICU event total; SOFA2_ONLY|hosp|vs_sofa2_negative shows the effect of that reference change.")
        for g, o, k, ref, key in comparisons:
                a, n1 = groups[g][k], groups[g]["n"]
                c, n0 = groups[ref][k], groups[ref]["n"]
                if min(a, c, n1 - a, n0 - c) < MIN_CELL:
                    res["odds_ratios"][key] = "omitted (< 5 events)"
                    res["suppressed"].append(f"OR {key}")
                    continue
                crude = (a * (n0 - c)) / ((n1 - a) * c)
                se = math.sqrt(1 / a + 1 / (n1 - a) + 1 / c + 1 / (n0 - c))
                # Mantel-Haenszel with RBG variance, and Breslow-Day
                R = S = Pp = Q = T = 0.0
                tabs = []
                for st in {(s["ab"], s["gender"], s["cvicu"]) for s in strata}:
                    x = next((s for s in strata if s["grp"] == g and (s["ab"], s["gender"], s["cvicu"]) == st), None)
                    r = next((s for s in strata if s["grp"] == ref and (s["ab"], s["gender"], s["cvicu"]) == st), None)
                    if not x or not r:
                        continue
                    aa, bb = int(x[k]), int(x["n"]) - int(x[k])
                    cc, dd = int(r[k]), int(r["n"]) - int(r[k])
                    nn = aa + bb + cc + dd
                    R += aa * dd / nn
                    S += bb * cc / nn
                    Pp += (aa + dd) / nn * aa * dd / nn
                    Q += (aa + dd) / nn * bb * cc / nn + (bb + cc) / nn * aa * dd / nn
                    T += (bb + cc) / nn * bb * cc / nn
                    tabs.append((aa, bb, cc, dd))
                mh = R / S
                se_mh = math.sqrt(Pp / (2 * R * R) + Q / (2 * R * S) + T / (2 * S * S))
                bd, df = 0.0, 0
                for aa, bb, cc, dd in tabs:
                    n1s, n0s, m1 = aa + bb, cc + dd, aa + cc
                    A, Bq, Cq = 1 - mh, n0s - m1 + mh * (n1s + m1), -mh * n1s * m1
                    lo_b, hi_b = max(0, m1 - n0s), min(n1s, m1)
                    if abs(A) < 1e-12:
                        roots = [-Cq / Bq]
                    else:
                        disc = math.sqrt(max(Bq * Bq - 4 * A * Cq, 0))
                        roots = [(-Bq + disc) / (2 * A), (-Bq - disc) / (2 * A)]
                    ah = next((x for x in roots if lo_b - 1e-9 <= x <= hi_b + 1e-9), None)
                    if ah is None or ah <= 0 or n1s - ah <= 0 or m1 - ah <= 0 or n0s - m1 + ah <= 0:
                        continue
                    var = 1 / (1 / ah + 1 / (n1s - ah) + 1 / (m1 - ah) + 1 / (n0s - m1 + ah))
                    bd += (aa - ah) ** 2 / var
                    df += 1
                df = max(df - 1, 1)
                res["odds_ratios"][key] = {
                    "crude": [round(crude, 2), round(math.exp(math.log(crude) - 1.96 * se), 2), round(math.exp(math.log(crude) + 1.96 * se), 2)],
                    "mantel_haenszel": [round(mh, 2), round(math.exp(math.log(mh) - 1.96 * se_mh), 2), round(math.exp(math.log(mh) + 1.96 * se_mh), 2)],
                    "breslow_day": {"chi2": round(bd, 1), "df": df, "p": float(f"{chi2_sf(bd, df):.3g}")},
                    "label": "descriptive/prognostic, not causal"}
        # ---- variants in the first suspicion window (Table 3) and matched pairs ---------------
        pairs = [("s2_std", "s1_std"), ("s2_gcs15", "s1_std"), ("s2_kcr", "s1_std"), ("s2_kuo22", "s1_std"),
                 ("s2_e", "s1_std"), ("s2_std", "s1_gcsalt"), ("s2_kcr", "s1_kcr"),
                 ("s2_g", "s1_std"), ("s2_s", "s1_std"), ("s2_h", "s1_std"), ("s2_i", "s1_std")]
        sel = ",\n".join(
            f"COUNTIF(fw_{x} >= 2) {x}__{y}__pos_x, COUNTIF(fw_{y} >= 2) {x}__{y}__pos_y, "
            f"COUNTIF(fw_{x} >= 2 AND fw_{y} >= 2) {x}__{y}__both, COUNTIF(fw_{x} >= 2 AND fw_{y} < 2) {x}__{y}__only_x, "
            f"COUNTIF(fw_{x} < 2 AND fw_{y} >= 2) {x}__{y}__only_y, "
            f"SUM(IF(fw_{x} >= 2 AND fw_{y} < 2, hospital_expire_flag, 0)) {x}__{y}__only_x_deaths"
            for x, y in pairs)
        v = run(f"SELECT COUNT(*) n, {sel} FROM `{C}` WHERE first_elig_stay AND NOT not_assessable_first")[0]
        res["first_window_variants"] = {"n": int(v["n"])}
        for x, y in pairs:
            g = {m: int(v[f"{x}__{y}__{m}"]) for m in ("pos_x", "pos_y", "both", "only_x", "only_y", "only_x_deaths")}
            if g["only_x_deaths"] >= MIN_CELL and g["only_x"] - g["only_x_deaths"] >= MIN_CELL:
                g["only_x_mortality"] = wilson(g["only_x_deaths"], g["only_x"])
            else:
                g["only_x_deaths"] = "<5"
                res["suppressed"].append(f"first-window {x} vs {y} only-x deaths")
            res["first_window_variants"][f"{x}_vs_{y}"] = g
        # ---- single-chart urine rule (variant s) against SOFA, any window ------------------------
        s = run(f"""SELECT COUNT(*) n, COUNTIF(any_s2_s) pos, COUNTIF(any_s2_s AND NOT any_s1_std) s_only,
                           COUNTIF(any_s1_std AND NOT any_s2_s) sofa_only,
                           SUM(IF(any_s2_s AND NOT any_s1_std, hospital_expire_flag, 0)) s_only_deaths,
                           COUNTIF(any_s2_s AND NOT any_s1_std AND NOT any_s2_std) s_only_not_primary,
                           SUM(IF(any_s2_s AND NOT any_s1_std AND NOT any_s2_std, hospital_expire_flag, 0)) s_only_not_primary_deaths
                    FROM `{C}` WHERE {PRIM}""")[0]
        sn = int(s["n"])
        row = {"rate_s": wilson(int(s["pos"]), sn)}
        for k in ("s_only", "sofa_only", "s_only_not_primary"):
            v = int(s[k])
            row[k] = v if not 0 < v < MIN_CELL else "<5"
        for k, dk in (("s_only", "s_only_deaths"), ("s_only_not_primary", "s_only_not_primary_deaths")):
            n_, d_ = int(s[k]), int(s[dk])
            if d_ >= MIN_CELL and n_ - d_ >= MIN_CELL:
                row[f"{k}_hosp_mortality"] = wilson(d_, n_)
            else:
                res["suppressed"].append(f"variant s {k} deaths")
        d = (int(s["s_only"]) - int(s["sofa_only"])) / sn
        se = math.sqrt(((int(s["s_only"]) + int(s["sofa_only"])) / sn - d * d) / sn)
        row["paired_rate_difference_s_minus_sofa"] = [d, d - 1.96 * se, d + 1.96 * se]
        res["single_chart_rule_s"] = row
        # ---- transition: groups under the single-chart rule (s) against groups under the primary rule --------
        gs = GRP.replace("any_s2_std", "any_s2_s")
        tr = run(f"""SELECT {gs} g_s, {GRP} g_p, COUNT(*) n, SUM(hospital_expire_flag) d
                     FROM `{C}` WHERE {PRIM} GROUP BY 1, 2""")
        cells = {(r["g_s"], r["g_p"]): [int(r["n"]), int(r["d"])] for r in tr}
        tt = {}
        for (a_, b_), (n_, d_) in sorted(cells.items()):
            tt[f"{a_}__{b_}"] = n_ if not 0 < n_ < MIN_CELL else "<5"
        res["transition_s_to_primary"] = tt
        # ---- bounded weight (variant w): weights outside 30-300 kg lose the urine criteria ----------------
        w = run(f"""SELECT COUNT(*) n, COUNTIF(any_s2_w) pos, COUNTIF(any_s2_w AND NOT any_s1_std) w_only,
                           COUNTIF(any_s1_std AND NOT any_s2_w) sofa_only,
                           COUNTIF(any_s2_std AND NOT any_s2_w) lost, COUNTIF(any_s2_w AND NOT any_s2_std) gained
                    FROM `{C}` WHERE {PRIM}""")[0]
        wn = int(w["n"])
        roww = {"rate_w": wilson(int(w["pos"]), wn)}
        for k in ("w_only", "sofa_only", "lost", "gained"):
            v = int(w[k])
            roww[k] = v if not 0 < v < MIN_CELL else "<5"
        d = (int(w["w_only"]) - int(w["sofa_only"])) / wn
        se = math.sqrt(((int(w["w_only"]) + int(w["sofa_only"])) / wn - d * d) / wn)
        roww["paired_rate_difference_w_minus_sofa"] = [d, d - 1.96 * se, d + 1.96 * se]
        res["bounded_weight_w"] = roww
        # ---- urine-output share (primary any-window definition) ------------------------------
        u = run(f"""SELECT COUNT(*) n, COUNTIF(NOT any_s2_kcr) lose_c, COUNTIF(NOT any_s2_kuo22) lose_d,
                           COUNTIF(NOT any_s2_g) lose_g, COUNTIF(one_volume_s2) one_vol,
                           COUNTIF(any_s2_g) pos_g, COUNTIF(any_s2_g AND NOT any_s2_kcr) uo_dep_g
                    FROM `{C}` WHERE {PRIM} AND any_s2_std AND NOT any_s1_std""")[0]
        un = int(u["n"])
        res["uo_share_sofa2_only"] = {"n_sofa2_only": un,
                                      "lose_without_uo_c": [int(u["lose_c"])] + wilson(int(u["lose_c"]), un)[0:3],
                                      "lose_uo_from_hr22_d": [int(u["lose_d"])] + wilson(int(u["lose_d"]), un)[0:3],
                                      "lose_under_coverage_g": [int(u["lose_g"])] + wilson(int(u["lose_g"]), un)[0:3],
                                      "rest_on_one_volume": [int(u["one_vol"])] + wilson(int(u["one_vol"]), un)[0:3],
                                      # counted directly (variant g can score the kidney lower
                                      # than variant c through the RRT-criteria branch, so arithmetic on lose_c
                                      # and lose_g gives only a bound)
                                      "positive_under_g": int(u["pos_g"])}
        # Complementary control: the stays not urine-dependent under c but lost under g
        # equal (n - lose_c) - (pos_g - uo_dep_g); if that is 1-4, the count and CI would reveal it, so only
        # the share, rounded to print precision, is kept.
        pos_g, dep_g = int(u["pos_g"]), int(u["uo_dep_g"])
        comp = (un - int(u["lose_c"])) - (pos_g - dep_g)
        if 0 < comp < MIN_CELL:
            res["uo_share_sofa2_only"]["urine_dependent_under_g"] = {"share": round(dep_g / pos_g, 3),
                                                                     "count_and_ci": "withheld (complementary cell 1-4)"}
            res["suppressed"].append("urine_dependent_under_g count and CI (complementary cell)")
        else:
            res["uo_share_sofa2_only"]["urine_dependent_under_g"] = [dep_g] + wilson(dep_g, pos_g)[0:3]
        # ---- cases by urine dependence, per concordance group, with mortality ------------
        ov = run(f"""
          SELECT {GRP} grp, cls, COUNT(*) n, SUM(hospital_expire_flag) d FROM (
            SELECT *, CASE WHEN one_volume_s2 THEN 'one_volume'
                           WHEN any_s2_std AND NOT any_s2_kcr THEN 'uo_dependent_multi_volume'
                           WHEN any_s2_std THEN 'not_uo_dependent' ELSE 'sofa2_negative' END cls
            FROM `{C}` WHERE {PRIM})
          GROUP BY 1, 2""")
        res["uo_dependence_by_group"] = {}
        for r in ov:
            nn, d = int(r["n"]), int(r["d"])
            row = {"n": nn if nn >= MIN_CELL or nn == 0 else "<5"}
            if d >= MIN_CELL and nn - d >= MIN_CELL:
                row["hosp_deaths"], row["hosp_mortality"] = d, wilson(d, nn)
            elif d == 0:
                row["hosp_deaths"] = 0
            else:
                row["hosp_deaths"] = "<5"
                res["suppressed"].append(f"uo dependence {r['grp']} {r['cls']} deaths")
            res["uo_dependence_by_group"].setdefault(r["grp"], {})[r["cls"]] = row
        # One-volume stays by n only, with deaths for all urine-dependent stays combined, whenever the one-volume deaths
        # are withheld (1-4): with the group total and the other subgroups printed, they would be recoverable.
        for grp in ("BOTH", "SOFA2_ONLY"):
            bo = {r["cls"]: (int(r["n"]), int(r["d"])) for r in ov if r["grp"] == grp}
            gr = res["uo_dependence_by_group"].get(grp, {})
            if "one_volume" not in bo or isinstance(gr.get("one_volume", {}).get("hosp_deaths"), int):
                continue
            n1v = bo["one_volume"][0]
            nud = n1v + bo.get("uo_dependent_multi_volume", (0, 0))[0]
            dud = bo["one_volume"][1] + bo.get("uo_dependent_multi_volume", (0, 0))[1]
            gr.pop("uo_dependent_multi_volume", None)
            gr["one_volume"] = {"n": n1v if n1v >= MIN_CELL else "<5"}
            gr["uo_dependent_all"] = ({"n": nud, "hosp_deaths": dud, "hosp_mortality": wilson(dud, nud)}
                                      if dud >= MIN_CELL and nud - dud >= MIN_CELL else {"n": nud, "hosp_deaths": "<5"})
            res["suppressed"] = [x for x in res["suppressed"] if x != f"uo dependence {grp} one_volume deaths"]
            res["suppressed"].append(f"uo dependence {grp}: one-volume deaths merged into all urine-dependent stays")
        # ---- crude risk differences vs Neither by urine dependence (arithmetic on the counts
        # above; Newcombe hybrid-score CI from the Wilson intervals; post hoc, crude) -----------------------------
        ne = res["groups"]["NEITHER"]
        n2, d2 = int(ne["n"]), int(ne["hosp_deaths"])
        p2, l2, u2 = wilson(d2, n2)[0:3]
        res["uo_dependence_risk_difference_vs_neither"] = {}
        for cls, row in res["uo_dependence_by_group"].get("SOFA2_ONLY", {}).items():
            if not isinstance(row.get("hosp_deaths"), int) or row["hosp_deaths"] < MIN_CELL:
                continue
            p1, l1, u1 = wilson(row["hosp_deaths"], row["n"])[0:3]
            dd = p1 - p2
            res["uo_dependence_risk_difference_vs_neither"][cls] = [
                dd, dd - math.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2), dd + math.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)]
        # ---- eligibility by stay identifier (sensitivity row) --------------------------------
        o = run(f"""SELECT COUNT(*) n, COUNTIF(old_any_s1_std) p1, COUNTIF(old_any_s2_std) p2,
                           COUNTIF(old_any_s1_std AND old_any_s2_std) n11, COUNTIF(old_any_s1_std AND NOT old_any_s2_std) n10,
                           COUNTIF(NOT old_any_s1_std AND old_any_s2_std) n01, SUM(hospital_expire_flag) d
                    FROM `{C}` WHERE first_elig_stay_old AND NOT not_assessable_old""")[0]
        on = int(o["n"])
        res["stayid_rule"] = {"n": on, "rate_sofa": wilson(int(o["p1"]), on), "rate_sofa2": wilson(int(o["p2"]), on),
                              "cells": {"both": int(o["n11"]), "sofa_only": int(o["n10"]), "sofa2_only": int(o["n01"]),
                                        "neither": on - int(o["n11"]) - int(o["n10"]) - int(o["n01"])},
                              "hosp_mortality": wilson(int(o["d"]), on)}
        # ---- era description: rates and brain-component prevalence, no test ----------
        er = run(f"""SELECT era, COUNT(*) n, COUNTIF(any_s1_std) s1, COUNTIF(any_s2_std) s2, COUNTIF(any_s2_h) s2h,
                            COUNTIF(any_s2_i) s2i, COUNTIF(any_s2_brain_pos) brain_any, COUNTIF(d1_s2_brain >= 1) brain_d1,
                            COUNTIF(emar_adm) emar
                     FROM `{C}` WHERE {PRIM} GROUP BY era ORDER BY era""")
        res["era"] = {r["era"]: {"n": int(r["n"]), "rate_sofa": wilson(int(r["s1"]), int(r["n"])),
                                 "rate_sofa2": wilson(int(r["s2"]), int(r["n"])),
                                 "rate_sofa2_h_no_delirium": wilson(int(r["s2h"]), int(r["n"])),
                                 "rate_sofa2_i_emar_or_inputevents": wilson(int(r["s2i"]), int(r["n"])),
                                 "brain_component_ge1_any_window": wilson(int(r["brain_any"]), int(r["n"])),
                                 "brain_component_ge1_day1": wilson(int(r["brain_d1"]), int(r["n"])),
                                 "admission_with_emar": wilson(int(r["emar"]), int(r["n"])),
                                 "label": "descriptive, not a test"} for r in er}
        # ---- individual baseline (sensitivity 8) ---------------------------------------------
        ib = run(f"""SELECT {GRP.replace('any_s1_std', 'ind_s1_std').replace('any_s2_std', 'ind_s2_std')} grp,
                            COUNT(*) n, SUM(hospital_expire_flag) d, AVG(ind_base_s1) base_s1, AVG(ind_base_s2) base_s2
                     FROM `{C}` WHERE first_elig_stay AND ind_baseline_available GROUP BY 1""")
        res["individual_baseline"] = {}
        for r in ib:
            d, nn = int(r["d"]), int(r["n"])
            row = {"n": nn if nn >= MIN_CELL else "<5"}
            if d >= MIN_CELL and nn - d >= MIN_CELL:
                row["hosp_mortality"] = wilson(d, nn)
            else:
                res["suppressed"].append(f"individual baseline {r['grp']} deaths")
            res["individual_baseline"][r["grp"]] = row
        res["individual_baseline"]["mean_baseline"] = {
            "sofa": round(sum(num(r["base_s1"]) * int(r["n"]) for r in ib) / sum(int(r["n"]) for r in ib), 2),
            "sofa2": round(sum(num(r["base_s2"]) * int(r["n"]) for r in ib) / sum(int(r["n"]) for r in ib), 2)}
        # ---- subgroup classification rates (sensitivity 5, 6, 7; CVICU) ----------------------
        sub = run(f"""
          SELECT s, COUNT(*) n, COUNTIF(any_s1_std) pos_s1, COUNTIF(any_s2_std) pos_s2, SUM(hospital_expire_flag) d FROM (
            SELECT 'non_cvicu' s, * FROM `{C}` WHERE {PRIM} AND NOT cvicu
            UNION ALL SELECT 'cvicu', * FROM `{C}` WHERE {PRIM} AND cvicu
            UNION ALL SELECT 'all_stays', * FROM `{C}` WHERE NOT not_assessable_any
            UNION ALL SELECT 'in_icu_suspicion', * FROM `{C}` WHERE {PRIM} AND first_susp_hr >= 0
            UNION ALL SELECT 'landmark', * FROM `{C}` WHERE {PRIM} AND landmark)
          GROUP BY s""")
        res["subgroups"] = {r["s"]: {"n": int(r["n"]), "rate_sofa": wilson(int(r["pos_s1"]), int(r["n"])),
                                     "rate_sofa2": wilson(int(r["pos_s2"]), int(r["n"])),
                                     "hosp_mortality": wilson(int(r["d"]), int(r["n"]))} for r in sub}
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    (OUT / "classification.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"agreement": res["agreement"], "suppressed": res["suppressed"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
