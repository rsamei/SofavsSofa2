"""Step 15: SOFA-2 against the original SOFA, both scored as published on the same hourly inputs.

The comparator (d1_s1_orig) is the original SOFA definitions (Vincent 1996) on the SOFA-2 pipeline's hourly inputs
(as d1_s1_pipe, step 11) with the original brain rule instead of SOFA-2's: a GCS that cannot be assessed (sedation,
or verbal not assessable) takes the last complete GCS before it, carried forward for up to 72 hours (none -> 0 points);
there is no best-motor substitute, which is a SOFA-2 rule. Sensitivity (d1_s1_orig_n): a not-assessable GCS outside
sedation is scored as normal, as in mimic-code. d1_s1_origbrain keeps the mimic-code SOFA and changes only its brain
rule to the original one. Because the day-1 score is a 24-hour worst value, a carried-forward GCS cannot exceed the
window maximum; the original rule can only remove points charted during sedation. Reads BUILD_DATASET.cohort (step 3).
Only aggregates leave BigQuery; counts below 5 are suppressed.

Reports: day-1 AUROC and paired DeLong contrasts; any-window Sepsis-3 rates (SOFA >= 2); day-1 brain points by rule.
"""
import json
import sys

from bq_client import BUILD, OUT, PROJECT, BQError, query
from delong import run_delong

T = f"{BUILD}.cohort"
PRIM = "first_elig_stay AND NOT not_assessable_any"
MIN_CELL = 5
SETS = {"d1_hosp": PRIM, "d1_hosp_susp_by_hr23": PRIM + " AND first_susp_hr <= 23",
        "d1_hosp_resp_measured": PRIM + " AND resp_measured_both_d1"}
SYS_EXPR = {"s2_std": "d1_s2_std", "s1_std": "d1_s1_std", "s1_gcsalt": "d1_s1_gcsalt", "s1_pipe": "d1_s1_pipe",
            "s1_orig": "d1_s1_orig", "s1_orig_n": "d1_s1_orig_n",
            "s1_origbrain": "d1_s1_std - d1_s1_cns + d1_s1_cns_orig"}
CONTRASTS = {
    "s2_minus_s1orig": {"s2_std": 1, "s1_orig": -1},
    "s2_minus_s1orign": {"s2_std": 1, "s1_orig_n": -1},
    "s2_minus_s1pipe": {"s2_std": 1, "s1_pipe": -1},
    "s2_minus_s1": {"s2_std": 1, "s1_std": -1},
    "s1orig_minus_s1": {"s1_orig": 1, "s1_std": -1},
    "s1orig_minus_s1pipe": {"s1_orig": 1, "s1_pipe": -1},
    "s1orign_minus_s1orig": {"s1_orig_n": 1, "s1_orig": -1},
    "s2_minus_s1origbrain": {"s2_std": 1, "s1_origbrain": -1},
    "s1origbrain_minus_s1": {"s1_origbrain": 1, "s1_std": -1},
    "s1gcsalt_minus_s1origbrain": {"s1_gcsalt": 1, "s1_origbrain": -1},
}

CLASS_SQL = f"""
SELECT COUNT(*) AS n,
  COUNTIF(any_s1_std) AS s1_std, COUNTIF(any_s1_pipe) AS s1_pipe, COUNTIF(any_s1_orig) AS s1_orig,
  COUNTIF(any_s1_orig_n) AS s1_orig_n, COUNTIF(any_s2_std) AS s2_std,
  COUNTIF(any_s2_std AND any_s1_orig) AS both_s2_s1orig, COUNTIF(any_s2_std AND NOT any_s1_orig) AS s2_only_s1orig,
  COUNTIF(NOT any_s2_std AND any_s1_orig) AS s1orig_only
FROM `{T}` WHERE {PRIM}
"""

BRAIN_SQL = f"""
SELECT COUNT(*) AS n,
  COUNTIF(d1_s1_cns_orig <> d1_s1_cns_alt) AS orig_differs_from_sofa2_rule,
  COUNTIF(d1_s1_cns_orig > d1_s1_cns_alt) AS orig_higher, COUNTIF(d1_s1_cns_orig < d1_s1_cns_alt) AS orig_lower,
  COUNTIF(d1_s1_cns_orig <> d1_s1_cns) AS orig_differs_from_mimic_code,
  COUNTIF(d1_s1_cns_orig > d1_s1_cns) AS orig_higher_than_mimic_code,
  COUNTIF(d1_s1_cns_orig < d1_s1_cns) AS orig_lower_than_mimic_code,
  AVG(d1_s1_cns) AS mean_cns_mimic_code, AVG(d1_s1_cns_alt) AS mean_cns_sofa2_rule, AVG(d1_s1_cns_orig) AS mean_cns_orig
FROM `{T}` WHERE {PRIM}
"""


def one_row(sql):
    f, rows, _ = query(PROJECT, sql)
    d = dict(zip(f, rows[0]))
    out = {}
    for k, v in d.items():
        if v is None:
            out[k] = None
        elif "." in str(v):
            out[k] = float(v)
        else:
            n = int(v)
            out[k] = n if n == 0 or n >= MIN_CELL else "<5"
    return out


def main():
    try:
        disc = run_delong(f"`{T}`", SETS, SYS_EXPR, CONTRASTS)
        cls = one_row(CLASS_SQL)
        brain = one_row(BRAIN_SQL)
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    floor_ok = all(r["n_pos"] >= 100 and r["n_neg"] >= 100 for r in disc)
    res = {"event_floor_ok": floor_ok, "discrimination": disc, "sepsis3_any_window": cls, "day1_brain_points": brain}
    (OUT / "original_sofa.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    for r in disc:
        print(f"{r['set']:24s} {r['contrast']:22s} {r['est']:.4f} ({r['lo']:.4f} to {r['hi']:.4f})")
    print("sepsis3", cls)
    print("brain", brain)
    return 0 if floor_ok else 1


if __name__ == "__main__":
    sys.exit(main())
