"""Step 11: order-independent decomposition of the day-1 AUROC difference, and SOFA scored on the SOFA-2 pipeline.
Reads BUILD_DATASET.variant_stay (step 8). Only aggregates leave BigQuery.

1. Shapley decomposition. Four rules are aligned to SOFA's handling, in every combination (16 coalitions):
     B  brain: the mimic-code convention for a non-assessable or sedated GCS (variant b);
     D  delirium: no delirium point (variant h);
     R  respiration: PaO2/FiO2 only, no SpO2/FiO2 fallback (variant j);
     U  urine output: kidney by creatinine or RRT only on SOFA-2 (variant c) and on SOFA (variant f).
   For a coalition S, delta(S) = AUROC(SOFA-2 with the rules in S aligned) - AUROC(SOFA, with its urine term removed when
   U is in S). The rules change separate components (brain and delirium both change the brain component and are
   combined explicitly), so each coalition's day-1 score is the primary score plus the component differences. The
   Shapley value of a rule is its average marginal change in delta over all orderings; the four values sum to
   delta(all) - delta(none). Each is a linear contrast of AUROCs, so its CI comes from the paired DeLong variance.
2. SOFA on the SOFA-2 pipeline (d1_s1_pipe): the original SOFA definitions scored from the SOFA-2 pipeline's hourly
   inputs, with SOFA-2's handling of a non-assessable GCS. SOFA-2 minus it isolates the definitions on one pipeline;
   it minus SOFA with the same brain handling (variant a) isolates the pipeline.
"""
import itertools
import json
import math
import sys

from bq_client import BUILD, OUT, BQError
from delong import run_delong

T = f"{BUILD}.variant_stay"
PRIM = "first_elig_stay AND NOT not_assessable_any"
SETS = {"d1_hosp": PRIM, "susp_by_hr23": PRIM + " AND first_susp_hr <= 23",
        "resp_measured": PRIM + " AND resp_measured_both_d1"}
PLAYERS = ["B", "D", "R", "U"]


def s2_expr(S):
    e = ["d1_s2_std"]
    if "B" in S and "D" in S:
        e.append("(d1_s2_gcs15_h - d1_s2_std)")
    elif "B" in S:
        e.append("(d1_s2_gcs15 - d1_s2_std)")
    elif "D" in S:
        e.append("(d1_s2_h - d1_s2_std)")
    if "R" in S:
        e.append("(d1_s2_j - d1_s2_std)")
    if "U" in S:
        e.append("(d1_s2_kcr - d1_s2_std)")
    return " + ".join(e)


def name(S):
    return "s2_" + ("".join(sorted(S)).lower() or "none")


COALITIONS = [frozenset(c) for k in range(5) for c in itertools.combinations(PLAYERS, k)]
SYS_EXPR = {name(S): s2_expr(S) for S in COALITIONS}
SYS_EXPR.update({"s1_std": "d1_s1_std", "s1_kcr": "d1_s1_kcr", "s1_gcsalt": "d1_s1_gcsalt", "s1_pipe": "d1_s1_pipe"})


def delta(S):
    """Coefficients of delta(S) = AUROC(S2(S)) - AUROC(S1(S))."""
    return {name(S): 1.0, ("s1_kcr" if "U" in S else "s1_std"): -1.0}


def add(a, b, w=1.0):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0.0) + w * v
    return out


CONTRASTS = {f"delta_{name(S)[3:]}": delta(S) for S in COALITIONS}
n = len(PLAYERS)
for p in PLAYERS:
    phi = {}
    for S in COALITIONS:
        if p in S:
            continue
        w = math.factorial(len(S)) * math.factorial(n - len(S) - 1) / math.factorial(n)
        phi = add(phi, delta(S | {p}), w)
        phi = add(phi, delta(S), -w)
    CONTRASTS[f"shapley_{p}"] = phi
CONTRASTS["s2_minus_s1pipe"] = {"s2_none": 1.0, "s1_pipe": -1.0}
CONTRASTS["s1pipe_minus_s1gcsalt"] = {"s1_pipe": 1.0, "s1_gcsalt": -1.0}
CONTRASTS["s1pipe_minus_s1"] = {"s1_pipe": 1.0, "s1_std": -1.0}


def main():
    try:
        out = run_delong(f"`{T}`", SETS, SYS_EXPR, CONTRASTS)
    except BQError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    floor_ok = all(r["n_pos"] >= 100 and r["n_neg"] >= 100 for r in out)
    by = {(r["set"], r["contrast"]): r for r in out}
    checks = {}
    for s in SETS:
        total = sum(by[(s, f"shapley_{p}")]["est"] for p in PLAYERS)
        target = by[(s, "delta_bdru")]["est"] - by[(s, "delta_none")]["est"]
        checks[s] = {"shapley_sum": total, "delta_all_minus_delta_none": target, "ok": abs(total - target) < 1e-9}
    res = {"players": {"B": "brain convention", "D": "delirium point", "R": "respiratory input", "U": "urine output"},
           "event_floor_ok": floor_ok, "checks": checks, "discrimination": out}
    (OUT / "decomposition.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    for s in SETS:
        print(s, {p: round(by[(s, f"shapley_{p}")]["est"], 4) for p in PLAYERS},
              "none", round(by[(s, "delta_none")]["est"], 4), "all", round(by[(s, "delta_bdru")]["est"], 4),
              "s2-s1pipe", round(by[(s, "s2_minus_s1pipe")]["est"], 4), "s1pipe-s1gcsalt", round(by[(s, "s1pipe_minus_s1gcsalt")]["est"], 4))
    return 0 if floor_ok and all(c["ok"] for c in checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
