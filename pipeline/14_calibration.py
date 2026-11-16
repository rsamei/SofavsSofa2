"""Step 14: recalibration, Brier scores, likelihood-ratio tests and external calibration of the day-1 scores, primary set.
Models are fitted to counts per score point (or per pair of points), so only aggregates leave BigQuery.

1. Logistic recalibration of each score for in-hospital death (logit p = a + b * score): Brier score, scaled Brier
   score (1 - Brier / Brier of the prevalence-only model), log-likelihood and AIC. Likelihood-ratio tests of adding
   SOFA-2 to a model with SOFA, and SOFA to a model with SOFA-2 (1 df each).
2. External calibration of the published day-1 score-to-ICU-mortality mappings (JAMA 2025 SOFA-2 paper, Figure 3C,
   reference/jama2025_fig3_counts.json) against ICU death: observed/expected, calibration-in-the-large (offset model),
   calibration slope, and bands for a calibration plot (points merged until each band has >= 5 deaths and survivors).
   Published points >= 20 are pooled because the top SOFA point has no survivors.
"""
import json
import math
import sys

from bq_client import BUILD, OUT, PROJECT, ROOT, BQError, query

MIN_CELL = 5
POOL_FROM = 20
PRIM = "first_elig_stay AND NOT not_assessable_any"


def solve(a, b):
    """Solve a x = b (small dense system, Gaussian elimination with partial pivoting)."""
    n = len(b)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(m[r][c]))
        m[c], m[p] = m[p], m[c]
        for r in range(n):
            if r != c:
                f = m[r][c] / m[c][c]
                for k in range(c, n + 1):
                    m[r][k] -= f * m[c][k]
    return [m[i][n] / m[i][i] for i in range(n)]


def inverse(a):
    n = len(a)
    cols = [solve(a, [1.0 if i == j else 0.0 for i in range(n)]) for j in range(n)]
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def expit(x):
    return 1 / (1 + math.exp(-x))


def logistic(groups, offset=None):
    """Grouped logistic regression. groups: list of (n, deaths, [covariates incl. 1 for the intercept]).
    offset: optional list of offsets (same order). Returns coefficients, covariance and log-likelihood."""
    k = len(groups[0][2])
    beta = [0.0] * k
    off = offset or [0.0] * len(groups)
    for _ in range(100):
        g = [0.0] * k
        h = [[0.0] * k for _ in range(k)]
        for (n, d, x), o in zip(groups, off):
            p = expit(o + sum(b * xi for b, xi in zip(beta, x)))
            w = n * p * (1 - p)
            for i in range(k):
                g[i] += (d - n * p) * x[i]
                for j in range(k):
                    h[i][j] += w * x[i] * x[j]
        step = solve(h, g)
        beta = [b + s for b, s in zip(beta, step)]
        if max(abs(s) for s in step) < 1e-10:
            break
    else:
        raise RuntimeError("logistic model did not converge")
    ll = 0.0
    for (n, d, x), o in zip(groups, off):
        p = expit(o + sum(b * xi for b, xi in zip(beta, x)))
        ll += d * math.log(p) + (n - d) * math.log(1 - p)
    return beta, inverse(h), ll


def chi2_sf_1df(x):
    return math.erfc(math.sqrt(x / 2))


def recalibration(pairs):
    """pairs: {(s1, s2): (n, deaths)}."""
    by = {"sofa": {}, "sofa2": {}}
    for (a, b), (n, d) in pairs.items():
        for name, s in (("sofa", a), ("sofa2", b)):
            m = by[name].setdefault(s, [0, 0])
            m[0] += n
            m[1] += d
    N = sum(n for n, d in pairs.values())
    D = sum(d for n, d in pairs.values())
    ybar = D / N
    brier0 = ybar * (1 - ybar)
    out = {"n": N, "deaths": D, "prevalence_brier": brier0}
    lls = {}
    for name in ("sofa", "sofa2"):
        groups = [(n, d, [1.0, float(s)]) for s, (n, d) in sorted(by[name].items())]
        beta, cov, ll = logistic(groups)
        brier = sum(d * (1 - expit(beta[0] + beta[1] * x[1])) ** 2 + (n - d) * expit(beta[0] + beta[1] * x[1]) ** 2
                    for n, d, x in groups) / N
        lls[name] = ll
        out[name] = {"intercept": beta[0], "slope_per_point": beta[1], "slope_se": math.sqrt(cov[1][1]),
                     "brier": brier, "scaled_brier": 1 - brier / brier0, "log_likelihood": ll, "aic": 4 - 2 * ll}
    groups = [(n, d, [1.0, float(a), float(b)]) for (a, b), (n, d) in sorted(pairs.items())]
    beta, cov, ll_both = logistic(groups)
    out["both"] = {"coef_sofa": beta[1], "coef_sofa2": beta[2], "log_likelihood": ll_both}
    for added, base in (("sofa2", "sofa"), ("sofa", "sofa2")):
        lr = 2 * (ll_both - lls[base])
        out[f"lr_add_{added}_to_{base}"] = {"chi2": lr, "df": 1, "p": chi2_sf_1df(lr)}
    out["aic_difference_sofa2_minus_sofa"] = out["sofa2"]["aic"] - out["sofa"]["aic"]
    return out


def published_risk(ref):
    e, d = ref["encounters"], ref["icu_deaths"]
    risk = {k: d[k] / e[k] for k in range(POOL_FROM)}
    top = sum(d[POOL_FROM:]) / sum(e[POOL_FROM:])
    for k in range(POOL_FROM, 25):
        risk[k] = top
    return risk


def logit(p):
    return math.log(p / (1 - p))


def wilson(k, n, z=1.96):
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [p, max(0.0, c - h), min(1.0, c + h)]


def bands(rows, risk):
    out, cur = [], None
    for k in sorted(rows, reverse=True):
        n, d = rows[k]
        if cur is None:
            cur = {"hi": k, "lo": k, "n": 0, "d": 0, "e": 0.0}
        cur["lo"], cur["n"], cur["d"], cur["e"] = k, cur["n"] + n, cur["d"] + d, cur["e"] + n * risk[k]
        if cur["d"] >= MIN_CELL and cur["n"] - cur["d"] >= MIN_CELL:
            out.append(cur)
            cur = None
    if cur is not None:
        last = out.pop()
        out.append({"hi": last["hi"], "lo": cur["lo"], "n": last["n"] + cur["n"], "d": last["d"] + cur["d"],
                    "e": last["e"] + cur["e"]})
    return [{"points": str(b["lo"]) if b["lo"] == b["hi"] else f"{b['lo']}-{b['hi']}", "n": b["n"], "deaths": b["d"],
             "observed": wilson(b["d"], b["n"]), "expected": b["e"] / b["n"]} for b in reversed(out)]


def external(rows, risk):
    o = sum(d for n, d in rows.values())
    n_all = sum(n for n, d in rows.values())
    e = sum(n * risk[k] for k, (n, d) in rows.items())
    s = 1.96 / math.sqrt(o)
    groups = [(n, d, [1.0]) for k, (n, d) in rows.items()]
    offs = [logit(risk[k]) for k in rows]
    b0, c0, _ = logistic(groups, offs)
    b1, c1, _ = logistic([(n, d, [1.0, logit(risk[k])]) for k, (n, d) in rows.items()])
    return {"n": n_all, "observed_deaths": o, "expected_deaths": e,
            "o_to_e": [o / e, o / e * math.exp(-s), o / e * math.exp(s)],
            "calibration_in_the_large": [b0[0], b0[0] - 1.96 * math.sqrt(c0[0][0]), b0[0] + 1.96 * math.sqrt(c0[0][0])],
            "slope": [b1[1], b1[1] - 1.96 * math.sqrt(c1[1][1]), b1[1] + 1.96 * math.sqrt(c1[1][1])],
            "bands": bands(rows, risk)}


def main():
    try:
        f, r, _ = query(PROJECT, f"""SELECT d1_s1_std a, d1_s2_std b, COUNT(*) n, SUM(hospital_expire_flag) d
                                     FROM `{BUILD}.cohort` WHERE {PRIM} GROUP BY 1, 2""")
        pairs = {(int(a), int(b)): (int(n), int(d)) for a, b, n, d in r}
        ext = {}
        for name, col in (("sofa1", "d1_s1_std"), ("sofa2", "d1_s2_std")):
            f, r, _ = query(PROJECT, f"""SELECT {col} x, COUNT(*) n, SUM(icu_death) d FROM `{BUILD}.cohort`
                                         WHERE {PRIM} GROUP BY 1""")
            ext[name] = {int(x): (int(n), int(d)) for x, n, d in r}
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    ref = json.loads((ROOT / "reference" / "jama2025_fig3_counts.json").read_text(encoding="utf-8"))
    res = {"recalibration_hospital_death": recalibration(pairs),
           "external_calibration_icu_death": {k: external(v, published_risk(ref[k])) for k, v in ext.items()},
           "published_source": ref["source"]}
    for k in res["external_calibration_icu_death"].values():
        for b in k["bands"]:
            assert b["deaths"] >= MIN_CELL and b["n"] - b["deaths"] >= MIN_CELL
    (OUT / "calibration.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    rc = res["recalibration_hospital_death"]
    print(json.dumps({"brier": {k: round(rc[k]["brier"], 5) for k in ("sofa", "sofa2")},
                      "scaled_brier": {k: round(rc[k]["scaled_brier"], 4) for k in ("sofa", "sofa2")},
                      "lr": {k: round(v["chi2"], 1) for k, v in rc.items() if k.startswith("lr_")},
                      "o_to_e": {k: [round(x, 2) for x in v["o_to_e"]] for k, v in res["external_calibration_icu_death"].items()}}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
