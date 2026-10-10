"""Supplementary Figure 10: disease-axis and distance sensitivity computed from exported statistics.

Reads the CSVs written by 2_Disease_Axis/06_Sensitivity_From_Exports.py and draws
a, lambda by gene set; b, raw and noise-corrected pairwise distances; c, observed t-profile correlations
against the correlation expected from the shared no-COPD reference alone.
Usage: python SFig10_Axis_Distance_Sensitivity.py <sensitivity_dir> <out_stem>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as ns

NEW, EXC = "#0072B2", "#D55E00"


def build(sens_dir, out_stem):
    ns.apply()
    sd = Path(sens_dir)
    lam = pd.read_csv(sd / "lambda_gene_set_sensitivity.csv")
    dist = pd.read_csv(sd / "distance_noise_corrected.csv")
    cor = pd.read_csv(sd / "tprofile_correlation_vs_shared_reference.csv")
    fig, axs = plt.subplots(1, 3, figsize=(ns.DOUBLE, 62 * ns.MM), gridspec_kw=dict(width_ratios=[1.25, 1.1, 1.0]))
    fig.subplots_adjust(left=0.17, right=0.99, top=0.80, bottom=0.17, wspace=0.62)

    ax = axs[0]
    labels = {"Axis transcripts (FDR < 0.05, primary)": "Axis, FDR < 0.05\n(primary)", "All transcripts": "All transcripts",
              "Axis, FDR < 1e-6": "Axis, FDR < 10⁻⁶", "Axis, 1e-6 <= FDR < 0.05": "Axis, 10⁻⁶ ≤ FDR < 0.05"}
    y = np.arange(len(lam))[::-1]
    ax.scatter(lam.lambda_new, y, color=NEW, s=14, zorder=3, label="Newly diagnosed")
    ax.scatter(lam.lambda_excluded, y, color=EXC, s=14, zorder=3, label="Potentially excluded")
    for yi, a, b in zip(y, lam.lambda_new, lam.lambda_excluded):
        ax.plot([b, a], [yi, yi], color="#BDBDBD", lw=0.6, zorder=1)
    ax.set_yticks(y); ax.set_yticklabels([labels[g] for g in lam.gene_set])
    for yi, n in zip(y, lam.n):
        ax.text(1.02, yi, f"{n:,}", va="center", ha="left", fontsize=5, color="#555", transform=ax.get_yaxis_transform(), clip_on=False)
    ax.set_xlim(0, 0.7); ax.set_xlabel("Disease-axis coefficient λ"); ax.set_ylim(-0.6, len(lam) - 0.4)
    ax.text(1.02, len(lam) - 0.45, "Genes", ha="left", va="bottom", fontsize=5, color="#555", transform=ax.get_yaxis_transform(), clip_on=False)
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, handletextpad=0.1, columnspacing=0.8)
    ns.panel_label(ax, "a", x=-0.62, y=1.13)

    ax = axs[1]
    u = dist[dist.gene_set.str.startswith("Union")].set_index("contrast")
    comps = [("2vs0", "Excluded–\nno COPD", EXC), ("2vs3", "Excluded–\nestablished", EXC),
             ("1vs0", "Newly–\nno COPD", NEW), ("1vs3", "Newly–\nestablished", NEW)]
    x = np.arange(len(comps)); w = 0.36
    for i, (c, lab, col) in enumerate(comps):
        ax.bar(i - w / 2, u.loc[c, "rms_raw"], w, color="white", edgecolor=col, lw=0.6, hatch="////")
        ax.bar(i + w / 2, u.loc[c, "rms_noise_corrected"], w, color=col, edgecolor=col, lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels([c[1] for c in comps], fontsize=5)
    ax.set_ylabel("RMS log₂ fold change"); ax.set_ylim(0, 0.1)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor="white", edgecolor="#555", hatch="////", lw=0.6, label="Uncorrected"),
                       Patch(facecolor="#777", edgecolor="#777", label="Noise floor removed")],
              frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, handletextpad=0.3, columnspacing=0.8)
    ns.panel_label(ax, "b", x=-0.30, y=1.13)

    ax = axs[2]
    pairs = [("1vs0 with 3vs0", "Newly vs\nestablished"), ("2vs0 with 3vs0", "Excluded vs\nestablished"),
             ("1vs0 with 2vs0", "Newly vs\nexcluded")]
    for i, (p, lab) in enumerate(pairs):
        r = cor[cor.pair == p].set_index("gene_set")
        ax.plot([i - 0.3, i + 0.3], [r.expected_from_shared_reference.iloc[0]] * 2, color="#595959", ls="--", lw=0.7)
        ax.scatter(i - 0.1, r.loc["Union", "spearman_observed"], marker="o", s=14, color="#333", zorder=3)
        ax.scatter(i + 0.1, r.loc["All transcripts", "spearman_observed"], marker="s", s=12, facecolor="white", edgecolor="#333", lw=0.6, zorder=3)
    ax.axhline(0, color="#BDBDBD", lw=0.4)
    ax.set_xticks(range(3)); ax.set_xticklabels([p[1] for p in pairs], fontsize=5)
    ax.set_ylabel("Spearman ρ of t-profiles"); ax.set_ylim(-0.2, 0.8); ax.set_xlim(-0.5, 2.5)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color="#333", ms=3, label="Union genes"),
                       Line2D([], [], marker="s", ls="", mfc="white", mec="#333", ms=3, mew=0.6, label="All genes"),
                       Line2D([], [], ls="--", color="#595959", lw=0.7, label="Expected from\nshared reference")],
              frameon=False, loc="upper right", ncol=1, handletextpad=0.3, borderaxespad=0.2)
    ns.panel_label(ax, "c", x=-0.32, y=1.13)
    Path(out_stem).parent.mkdir(parents=True, exist_ok=True)
    ns.save(fig, str(out_stem)); plt.close(fig)


if __name__ == "__main__":
    build(sys.argv[1], sys.argv[2])
