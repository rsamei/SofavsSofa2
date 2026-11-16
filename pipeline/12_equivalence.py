"""Step 12: equivalence tests (two one-sided tests, TOST) for the AUROC contrasts, from the estimates written by steps 4,
9 and 11 (no BigQuery).

Margin: +/-0.01 AUROC, chosen after the results were seen and reported as such. For each contrast: estimate, 90% CI,
the two one-sided p values and whether the 90% CI lies inside the margin (equivalence at alpha = 0.05), plus whether
the 95% CI excludes 0.
"""
import json
import math

from bq_client import OUT

MARGIN = 0.01
Z90 = 1.6448536269514722
CONTRASTS = [
    ("discrimination.json", "rows", "d1_hosp", "primary_s2_minus_s1"),
    ("discrimination.json", "rows", "d1_hosp", "matched_gcs_a_s2_minus_s1gcsalt"),
    ("discrimination.json", "rows", "d1_hosp", "matched_gcs_b_s2gcs15_minus_s1"),
    ("discrimination.json", "rows", "d1_hosp", "variant_s_s2s_minus_s1"),
    ("joint_contrasts.json", "discrimination", "primary", "joint_a"),
    ("joint_contrasts.json", "discrimination", "primary", "joint_b"),
    ("joint_contrasts.json", "discrimination", "full_day1", "joint_a"),
    ("joint_contrasts.json", "discrimination", "full_day1", "joint_b"),
    ("joint_contrasts.json", "discrimination", "susp_by_hr23", "joint_a"),
    ("joint_contrasts.json", "discrimination", "susp_by_hr23", "joint_b"),
    ("joint_contrasts.json", "discrimination", "resp_measured", "joint_a"),
    ("joint_contrasts.json", "discrimination", "resp_measured", "joint_b"),
    ("discrimination.json", "rows", "d1_hosp_susp_by_hr23", "primary_s2_minus_s1"),
    ("decomposition.json", "discrimination", "d1_hosp", "s2_minus_s1pipe"),
    ("decomposition.json", "discrimination", "d1_hosp", "s1pipe_minus_s1"),
    ("decomposition.json", "discrimination", "d1_hosp", "delta_bdru"),
]


def phi(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def main():
    res = {"margin": MARGIN, "note": "margin chosen after results were seen", "contrasts": []}
    for fname, key, set_, contrast in CONTRASTS:
        path = OUT / fname
        if not path.exists():
            continue
        rows = [r for r in json.loads(path.read_text(encoding="utf-8"))[key]
                if r["set"] == set_ and r["contrast"] == contrast]
        if not rows:
            continue
        est, se = float(rows[0]["est"]), float(rows[0]["se"])
        p_low = 1 - phi((est + MARGIN) / se)   # H0: diff <= -margin
        p_high = phi((est - MARGIN) / se)      # H0: diff >= +margin
        res["contrasts"].append({
            "set": set_, "contrast": contrast, "est": est, "se": se,
            "ci90": [est - Z90 * se, est + Z90 * se], "ci95": [est - 1.96 * se, est + 1.96 * se],
            "p_tost": max(p_low, p_high), "equivalent": max(p_low, p_high) < 0.05,
            "ci95_excludes_0": est - 1.96 * se > 0 or est + 1.96 * se < 0})
    (OUT / "equivalence.json").write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    for c in res["contrasts"]:
        print(f"{c['set']:22s} {c['contrast']:34s} {c['est']:+.4f}  90% CI {c['ci90'][0]:+.4f} to {c['ci90'][1]:+.4f}"
              f"  equivalent={c['equivalent']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
