"""Step 18: does SOFA-2's day-1 advantage reflect respiratory support rather than organ dysfunction? Reads
BUILD_DATASET.variant_stay (step 8), sofa2_hourly (step 2) and the mimic-code ventilation concept. Post hoc,
descriptive. Only aggregates leave BigQuery.

Day-1 indicators (ICU hours 0 to the day-1 scoring hour, as the day-1 score):
  imv_d1  invasive ventilation (ventilation_status 'InvasiveVent') overlapping day 1;
  adv_d1  advanced respiratory support as in the SOFA-2 respiratory rule (invasive or non-invasive ventilation,
          high-flow nasal oxygen, tracheostomy) overlapping day 1;
  ett_d1  a day-1 hour with the GCS verbal response not assessable (the trigger of SOFA-2's best-motor rule).

1. Likelihood-ratio tests (grouped logistic regression, linear score terms, as step 14) of adding SOFA-2 to SOFA plus
   an indicator, and SOFA to SOFA-2 plus an indicator, for each indicator; the unadjusted tests for reference.
2. Paired DeLong AUROCs of the fitted linear predictors: SOFA + indicator against SOFA-2 + indicator.
3. Day-1 AUROC of SOFA and SOFA-2 and the main rule contrasts within strata of each indicator: primary difference,
   brain rule aligned (variants a and b), respiratory input aligned (variant j), both aligned (joint contrasts).
"""
import importlib.util
import json
import math
import sys
from pathlib import Path

from bq_client import BUILD, OUT, PROJECT, UPSTREAM, BQError, query
from delong import run_delong

_spec = importlib.util.spec_from_file_location("calibration", Path(__file__).with_name("14_calibration.py"))
_c14 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_c14)
logistic, chi2_sf_1df = _c14.logistic, _c14.chi2_sf_1df

PRIM = "first_elig_stay AND NOT not_assessable_any"
INDICATORS = ["imv_d1", "adv_d1", "ett_d1"]

SOURCE = f"""(
  WITH c AS (SELECT * FROM `{BUILD}.variant_stay` WHERE {PRIM}),
  win AS (
    SELECT c.stay_id, MIN(h.starttime) AS t0, MAX(h.endtime) AS t1, MAX(COALESCE(h.unable, 0)) AS ett
    FROM c JOIN `{BUILD}.sofa2_hourly` h ON h.stay_id = c.stay_id AND h.hr BETWEEN 0 AND c.d1_last_hr
    GROUP BY 1
  ),
  v AS (
    SELECT w.stay_id,
      LOGICAL_OR(e.ventilation_status = 'InvasiveVent') AS imv,
      LOGICAL_OR(e.ventilation_status IN ('InvasiveVent', 'NonInvasiveVent', 'HFNC', 'Tracheostomy')) AS adv
    FROM win w JOIN `{UPSTREAM}.ventilation` e
      ON e.stay_id = w.stay_id AND e.starttime < w.t1 AND e.endtime > w.t0
    GROUP BY 1
  )
  SELECT c.*, CAST(COALESCE(v.imv, FALSE) AS INT64) AS imv_d1, CAST(COALESCE(v.adv, FALSE) AS INT64) AS adv_d1,
         CAST(COALESCE(w.ett, 0) = 1 AS INT64) AS ett_d1
  FROM c LEFT JOIN win w USING (stay_id) LEFT JOIN v USING (stay_id)
) z"""


def lr_tests(ind):
    f, r, _ = query(PROJECT, f"""SELECT d1_s1_std a, d1_s2_std b, {ind} v, COUNT(*) n, SUM(hospital_expire_flag) d
                                 FROM {SOURCE} GROUP BY 1, 2, 3""")
    cells = [(int(n), int(d), float(a), float(b), float(v)) for a, b, v, n, d in r]

    def fit(cols):
        groups = [(n, d, [1.0] + [{"a": a, "b": b, "v": v}[k] for k in cols]) for n, d, a, b, v in cells]
        beta, cov, ll = logistic(groups)
        return {"coef": dict(zip(["intercept"] + cols, beta)),
                "se": dict(zip(["intercept"] + cols, (math.sqrt(cov[i][i]) for i in range(len(beta))))), "ll": ll}

    m = {k: fit(cols) for k, cols in {"sofa": ["a"], "sofa2": ["b"], "sofa_ind": ["a", "v"], "sofa2_ind": ["b", "v"],
                                       "both_ind": ["a", "b", "v"], "both": ["a", "b"], "ind": ["v"]}.items()}
    lr = lambda big, small: {"chi2": 2 * (m[big]["ll"] - m[small]["ll"]), "df": 1,
                             "p": chi2_sf_1df(2 * (m[big]["ll"] - m[small]["ll"]))}
    return {"n": sum(c[0] for c in cells), "deaths": sum(c[1] for c in cells),
            "with_indicator": {"prevalence": sum(c[0] for c in cells if c[4] == 1) / sum(c[0] for c in cells),
                               "deaths": sum(c[1] for c in cells if c[4] == 1)},
            "lr_add_sofa2_to_sofa": lr("both", "sofa"), "lr_add_sofa_to_sofa2": lr("both", "sofa2"),
            "lr_add_sofa2_to_sofa_plus_indicator": lr("both_ind", "sofa_ind"),
            "lr_add_sofa_to_sofa2_plus_indicator": lr("both_ind", "sofa2_ind"),
            "lr_add_indicator_to_sofa": lr("sofa_ind", "sofa"), "lr_add_indicator_to_sofa2": lr("sofa2_ind", "sofa2"),
            "models": {k: {"coef": v["coef"], "se": v["se"], "log_likelihood": v["ll"]} for k, v in m.items()}}


def main():
    res = {"indicators": {"imv_d1": "invasive ventilation on day 1",
                          "adv_d1": "advanced respiratory support on day 1 (SOFA-2 respiratory rule)",
                          "ett_d1": "GCS verbal response not assessable on day 1 (best-motor rule trigger)"},
           "likelihood_ratio": {}, "model_auroc": {}, "stratified": None}
    try:
        for ind in INDICATORS:
            t = lr_tests(ind)
            res["likelihood_ratio"][ind] = t
            ca, cb = t["models"]["sofa_ind"]["coef"], t["models"]["sofa2_ind"]["coef"]
            res["model_auroc"][ind] = run_delong(
                SOURCE, {"primary": "TRUE"},
                # linear predictors on an integer scale (DeLong ranks by value; order and ties are unchanged)
                {"sofa_ind": f"CAST(ROUND(1e9 * ({ca['a']} * d1_s1_std + {ca['v']} * {ind})) AS INT64)",
                 "sofa2_ind": f"CAST(ROUND(1e9 * ({cb['b']} * d1_s2_std + {cb['v']} * {ind})) AS INT64)",
                 "ind": f"{ind}", "s1": "d1_s1_std", "s2": "d1_s2_std"},
                {"sofa2ind_minus_sofaind": {"sofa2_ind": 1, "sofa_ind": -1}})
        sets = {"all": "TRUE"}
        for ind in INDICATORS:
            sets[f"{ind}_yes"] = f"{ind} = 1"
            sets[f"{ind}_no"] = f"{ind} = 0"
        res["stratified"] = run_delong(
            SOURCE, sets,
            {"s1": "d1_s1_std", "s2": "d1_s2_std", "s1a": "d1_s1_gcsalt", "s2b": "d1_s2_gcs15", "s2j": "d1_s2_j",
             "s2bj": "d1_s2_gcs15 + d1_s2_j - d1_s2_std"},
            {"primary_s2_minus_s1": {"s2": 1, "s1": -1},
             "brain_a_s2_minus_s1a": {"s2": 1, "s1a": -1},
             "brain_b_s2b_minus_s1": {"s2b": 1, "s1": -1},
             "resp_j_s2j_minus_s1": {"s2j": 1, "s1": -1},
             "joint_a_s2j_minus_s1a": {"s2j": 1, "s1a": -1},
             "joint_b_s2bj_minus_s1": {"s2bj": 1, "s1": -1}})
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    floor_ok = all(r["n_pos"] >= 100 and r["n_neg"] >= 100
                   for rows in [res["stratified"]] + list(res["model_auroc"].values()) for r in rows)
    res["event_floor_ok"] = floor_ok
    (OUT / "ventilation.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    for ind, t in res["likelihood_ratio"].items():
        print(ind, f"prev {t['with_indicator']['prevalence']:.3f}",
              {k: round(v["chi2"], 1) for k, v in t.items() if k.startswith("lr_")})
        for r in res["model_auroc"][ind]:
            print("  ", r["contrast"], f"{r['est']:.4f} ({r['lo']:.4f} to {r['hi']:.4f})")
    for r in res["stratified"]:
        if not r["contrast"].startswith("auc_s1") or r["contrast"] in ("auc_s1", "auc_s2"):
            print(f"{r['set']:12s} {r['contrast']:24s} {r['est']:.4f} ({r['lo']:.4f} to {r['hi']:.4f}) {r['n_pos']}/{r['n_neg']}")
    return 0 if floor_ok else 1


if __name__ == "__main__":
    sys.exit(main())
