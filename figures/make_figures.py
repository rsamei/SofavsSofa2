"""Figures for the manuscript, drawn from the pipeline's JSON outputs (OUTPUT_DIR, default output/) into
OUTPUT_DIR/figures. Needs matplotlib.

  fig1_flow         participant flow
  fig2_roc          day-1 ROC curves for in-hospital death, SOFA and SOFA-2 (from counts per score point)
  fig3_forest       AUROC contrasts for the implementation variants, with the +/-0.01 equivalence margin
  fig4_concordance  Sepsis-3 concordance groups: size and in-hospital mortality
  fig5_calibration  external calibration of the published day-1 score-to-ICU-mortality mappings
"""
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402
import matplotlib.ticker  # noqa: E402,F401

OUT = Path(os.environ.get("OUTPUT_DIR", "output"))
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
SOFA, SOFA2 = "#2a78d6", "#eb6834"          # categorical slots 1-2 (validated pair)
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "white", "axes.facecolor": "white", "savefig.dpi": 300})


def load(name):
    return json.loads((OUT / name).read_text(encoding="utf-8"))


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)


def fig1_flow():
    f = load("aggregates.json")["flow"]
    c = load("cohort_check.json")["checks"]
    steps = [
        (f"First ICU stay of each admission, adults,\nICU stay ≥ 4 h: {f['cohort_base']:,}", None),
        (f"Eligible stays: {f['eligible']:,}", f"No suspicion-of-infection window\noverlapping the stay: {f['adults_without_overlapping_suspicion']:,}"),
        (f"First eligible stay per patient: {f['first_eligible']:,}", f"Later eligible stays of the\nsame patient: {f['eligible'] - f['first_eligible']:,}"),
        (f"Primary set: {f['primary_n']:,}", f"Not assessable (no scored hour\nin any window): {f['not_assessable']:,}"),
    ]
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ys = [9.0, 6.6, 4.2, 1.8]
    for (main, side), y in zip(steps, ys):
        ax.add_patch(FancyBboxPatch((0.4, y - 0.5), 4.8, 1.0, boxstyle="round,pad=0.04,rounding_size=0.08",
                                    fc=SURF, ec=INK2, lw=0.8))
        ax.text(2.8, y, main, ha="center", va="center", color=INK, fontsize=8.4)
        if side:
            ym = y + 1.2
            ax.add_patch(FancyBboxPatch((5.9, ym - 0.45), 3.9, 0.9, boxstyle="round,pad=0.04,rounding_size=0.08",
                                        fc="white", ec=INK2, lw=0.6))
            ax.text(7.85, ym, side, ha="center", va="center", color=INK2, fontsize=7.6)
            ax.plot([2.8, 5.9], [ym, ym], color=INK2, lw=0.8)
        if y < 9.0:
            ax.annotate("", xy=(2.8, y + 0.52), xytext=(2.8, y + 1.88), arrowprops=dict(arrowstyle="->", color=INK2, lw=0.8))
    ax.text(0.4, 0.55, f"Subsets of the primary set: suspicion by ICU hour 23 {f['primary_n'] - f['first_susp_after_hr23']:,}; "
                       f"landmark {c['landmark_n']:,}; complete case {c['complete_case_n']:,}\n"
                       f"Eras 2008–2013 / 2014–2016 / 2017–2022: {c['era_2008_2013_n']:,} / {c['era_2014_2016_n']:,} / "
                       f"{c['era_2017_2022_n']:,}", ha="left", va="center", color=INK2, fontsize=7.4)
    save(fig, "fig1_flow")


def roc_points(counts):
    """counts: {score: (n, deaths)} -> ROC coordinates for thresholds score >= k."""
    pos = sum(d for n, d in counts.values())
    neg = sum(n - d for n, d in counts.values())
    xs, ys = [0.0], [0.0]
    tp = fp = 0
    for k in sorted(counts, reverse=True):
        n, d = counts[k]
        tp += d
        fp += n - d
        xs.append(fp / neg)
        ys.append(tp / pos)
    return xs, ys


def fig2_roc():
    m = load("descriptives.json")["mortality_by_day1_point"]
    rows = {r["set"] + "|" + r["contrast"]: r for r in load("discrimination.json")["rows"]}
    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    ax.plot([0, 1], [0, 1], color=GRID, lw=1)
    for sys, col, key, lab in (("SOFA", SOFA, "auc_s1_std", "SOFA"), ("SOFA2", SOFA2, "auc_s2_std", "SOFA-2")):
        counts = {int(k.split("|")[1].lstrip(">=")): (v["n"], v["deaths"]) for k, v in m.items() if k.startswith(sys + "|")}
        xs, ys = roc_points(counts)
        auc = rows[f"d1_hosp|{key}"]
        ax.plot(xs, ys, color=col, lw=2, marker="o", ms=3,
                label=f"{lab}: AUROC {auc['est']:.3f} ({auc['lo']:.3f}–{auc['hi']:.3f})")
    ax.set_xlabel("1 − specificity")
    ax.set_ylabel("Sensitivity")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.grid(color=GRID, lw=0.5)
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    ax.set_title("Day-1 score, in-hospital death (primary set)", fontsize=9, color=INK, loc="left")
    save(fig, "fig2_roc")


def fig3_forest():
    rows = {r["set"] + "|" + r["contrast"]: r for r in load("discrimination.json")["rows"]}
    joint = {r["set"] + "|" + r["contrast"]: r for r in load("joint_contrasts.json")["discrimination"]}
    dec = {r["set"] + "|" + r["contrast"]: r for r in load("decomposition.json")["discrimination"]}
    items = [
        ("SOFA-2 − SOFA (primary)", rows["d1_hosp|primary_s2_minus_s1"], "Unaligned"),
        ("SOFA-2 − SOFA, suspicion by hour 23", rows["d1_hosp_susp_by_hr23|primary_s2_minus_s1"], "Unaligned"),
        ("SOFA-2 − SOFA with SOFA-2 brain handling (a)", rows["d1_hosp|matched_gcs_a_s2_minus_s1gcsalt"], "Brain rule aligned"),
        ("SOFA-2 with mimic-code brain convention (b) − SOFA", rows["d1_hosp|matched_gcs_b_s2gcs15_minus_s1"], "Brain rule aligned"),
        ("j − a (PaO2/FiO2-only respiration, brain aligned)", joint["primary|joint_a"], "Brain and respiratory input aligned"),
        ("(b + j) − SOFA", joint["primary|joint_b"], "Brain and respiratory input aligned"),
        ("SOFA-2 − SOFA on the SOFA-2 pipeline", dec["d1_hosp|s2_minus_s1pipe"], "Same pipeline"),
        ("Without urine criteria (c) − SOFA-2", None, "Single rules within SOFA-2"),
        ("Single-chart urine rule (s) − SOFA-2", rows["d1_hosp|direct_s2s_minus_s2"], "Single rules within SOFA-2"),
        ("No delirium point (h) − SOFA-2", rows["d1_hosp|direct_s2h_minus_s2"], "Single rules within SOFA-2"),
    ]
    # the pipeline reports SOFA-2 minus c; flip it to c minus SOFA-2
    r = rows["d1_hosp|direct_s2_minus_s2kcr"]
    items[7] = ("Without urine criteria (c) − SOFA-2", {"est": -r["est"], "lo": -r["hi"], "hi": -r["lo"]}, items[7][2])
    def fmt(x):
        return f"{x:+.4f}".replace("-", "−") if abs(x) >= 0.00005 else "0.0000"

    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.axvspan(-0.01, 0.01, color="#f0efec", lw=0)
    ax.axvline(0, color=INK2, lw=0.8)
    ticks, labels, y, last = [], [], 0, None
    for lab, rr, grp in items:
        if grp != last:
            y -= 0.4
            ticks.append(y)
            labels.append(grp)
            y -= 1
            last = grp
        ax.plot([rr["lo"], rr["hi"]], [y, y], color=INK, lw=1.2)
        ax.plot(rr["est"], y, "s", color=INK, ms=4.5)
        ax.text(1.02, y, f"{fmt(rr['est'])} ({fmt(rr['lo'])} to {fmt(rr['hi'])})", transform=ax.get_yaxis_transform(),
                fontsize=7.4, color=INK, va="center")
        ticks.append(y)
        labels.append("   " + lab)
        y -= 1
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=7.6)
    for tl, lab in zip(ax.get_yticklabels(), labels):
        if not lab.startswith("   "):
            tl.set_fontweight("bold")
            tl.set_color(INK2)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(y + 0.3, 0)
    ax.set_xlim(-0.02, 0.04)
    ax.set_xticks([-0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Difference in day-1 AUROC for in-hospital death (95% CI); shaded band ±0.01")
    save(fig, "fig3_forest")


def fig4_concordance():
    g = load("classification.json")["groups"]
    order = [("BOTH", "Both"), ("SOFA2_ONLY", "SOFA-2 only"), ("SOFA_ONLY", "SOFA only"), ("NEITHER", "Neither")]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.6), gridspec_kw={"width_ratios": [1, 1.2]}, sharey=True)
    ys = list(range(len(order)))[::-1]
    for y, (k, lab) in zip(ys, order):
        n = g[k]["n"]
        a1.barh(y, n, height=0.55, color=INK2)
        a1.text(n, y, f" {n:,}", va="center", fontsize=8, color=INK)
        m = g[k].get("hosp_mortality")
        if isinstance(m, list):
            a2.plot([m[1], m[2]], [y, y], color=INK, lw=1.2)
            a2.plot(m[0], y, "o", color=INK, ms=5)
            a2.text(m[2] + 0.004, y, f"{m[0]:.1%}", va="center", fontsize=8, color=INK)
    a1.set_yticks(ys)
    a1.set_yticklabels([lab for _, lab in order])
    a1.set_xscale("log")
    a1.set_xlim(100, 200000)
    a1.set_xlabel("Stays (log scale)")
    a1.tick_params(axis="y", length=0)
    a2.set_xlabel("In-hospital mortality (95% CI)")
    a2.set_xlim(0, 0.2)
    a2.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
    for a in (a1, a2):
        a.grid(axis="x", color=GRID, lw=0.5)
    a1.set_title("Sepsis-3 concordance groups (any window, baseline 0)", fontsize=9, color=INK, loc="left")
    save(fig, "fig4_concordance")


def fig5_calibration():
    ext = load("calibration.json")["external_calibration_icu_death"]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.5), sharey=True)
    for ax, (k, lab, colr) in zip(axes, (("sofa1", "SOFA", SOFA), ("sofa2", "SOFA-2", SOFA2))):
        ax.plot([0, 1], [0, 1], color=GRID, lw=1, ls="--")
        for b in ext[k]["bands"]:
            o = b["observed"]
            ax.plot([b["expected"], b["expected"]], [o[1], o[2]], color=colr, lw=1)
            ax.plot(b["expected"], o[0], "o", color=colr, ms=4)
        e = ext[k]
        ax.text(0.03, 0.95, f"{lab}\nO/E {e['o_to_e'][0]:.2f} ({e['o_to_e'][1]:.2f}–{e['o_to_e'][2]:.2f})\n"
                            f"slope {e['slope'][0]:.2f} ({e['slope'][1]:.2f}–{e['slope'][2]:.2f})",
                va="top", fontsize=8, color=INK)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal")
        ax.grid(color=GRID, lw=0.5)
        ax.set_xlabel("Published ICU mortality")
    axes[0].set_ylabel("Observed ICU mortality")
    save(fig, "fig5_calibration")


if __name__ == "__main__":
    fig1_flow()
    fig2_roc()
    fig3_forest()
    fig4_concordance()
    fig5_calibration()
    print(f"figures written to {FIG}")
