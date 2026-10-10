"""Figure 3: volcano plots for each disease group versus no COPD (Python port of Fig3_Volcano_Plots.R).

Same selection as the R version: common limits (|log2 FC| <= 1, -log10 P <= 25), capped points as triangles,
dashed line at the largest unadjusted P among BH FDR < 0.05 transcripts, and the three most significant
transcripts per direction labelled. Drawn at Nature Communications specifications (183 mm, 5-7 pt type).
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as ns

#### Edit paths to match your file system ####
SOURCE = Path("path/to/Schema_omics_csv_bhatt/copdschema/extracted/gene_de_Model_1_DisCOPDAI_allresults_2026.04.25.csv")
OUT_STEM = Path(__file__).resolve().parents[1] / "figures" / "Fig3_volcano"

CONTRASTS = {"DisCOPDAI_1vs0": "Newly diagnosed vs no COPD", "DisCOPDAI_2vs0": "Potentially excluded vs no COPD",
             "DisCOPDAI_3vs0": "Established COPD vs no COPD"}
XL, YL = 1.0, 25.0
COL = {"up": "#B2182B", "down": "#2166AC", "ns": "#BDBDBD"}


def place_labels(ax, pts, labels):
    """Greedy leader-line placement: push each label away from the origin until it clears earlier labels."""
    placed = []
    for (x, y), lab in sorted(zip(pts, labels), key=lambda t: -t[0][1]):
        dx = 0.18 if x >= 0 else -0.18; tx, ty = np.clip(x + dx, -0.92, 0.92), min(y + 1.2, 23.5)
        for _ in range(40):
            if all(abs(tx - px) > 0.24 or abs(ty - py) > 1.3 for px, py in placed): break
            ty -= 1.4
        placed.append((tx, ty))
        ax.annotate(lab, (x, y), (tx, ty), fontsize=5, style="italic", ha="left" if dx > 0 else "right", va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.3, color="#555", shrinkA=0.5, shrinkB=0.5), zorder=5)


def load_de(source=SOURCE):
    de = pd.read_csv(source)
    de = de[de.Contrast.isin(CONTRASTS)].copy()
    de["mlp"] = -np.log10(de["P.Value"].clip(lower=np.finfo(float).tiny))
    de["status"] = np.where(de["adj.P.Val"] < 0.05, np.where(de.logFC > 0, "up", "down"), "ns")
    de["capped"] = (de.logFC.abs() > XL) | (de.mlp > YL)
    de["dx"], de["dy"] = de.logFC.clip(-XL, XL), de.mlp.clip(upper=YL)
    return de


def draw_volcano(ax, d, title):
    """One volcano panel; returns the labelled transcripts."""
    nsig = int((d["adj.P.Val"] < 0.05).sum())
    for st in ("ns", "down", "up"):
        s = d[(d.status == st) & ~d.capped]
        ax.scatter(s.dx, s.dy, s=0.6, color=COL[st], alpha=0.6, lw=0, rasterized=True, zorder=2 if st == "ns" else 3)
        c = d[(d.status == st) & d.capped]
        ax.scatter(c.dx, c.dy, s=4, marker="^", facecolors="white", edgecolors=COL[st], lw=0.4, zorder=4, clip_on=False)
    if nsig:
        ax.axhline(-np.log10(d.loc[d["adj.P.Val"] < 0.05, "P.Value"].max()), color="#595959", ls="--", lw=0.5, zorder=1)
    ax.axvline(0, color="#B3B3B3", lw=0.4, zorder=1)
    lab = (d[(d["adj.P.Val"] < 0.05) & d.gene_label.notna() & (d.gene_label != "")]
           .assign(a=lambda x: -x.logFC.abs()).sort_values(["status", "adj.P.Val", "a"]).groupby("status").head(3))
    if len(lab):
        place_labels(ax, list(zip(lab.dx, lab.dy)), list(lab.gene_label))
    ax.set_xlim(-XL, XL); ax.set_ylim(0, YL); ax.set_xticks([-0.5, 0, 0.5]); ax.set_yticks(range(0, 26, 5))
    ax.set_title(f"{title}\n{nsig:,} transcripts at FDR < 0.05", loc="center", fontsize=6)
    return lab.assign(contrast=title)


def legend_handles():
    return [Line2D([], [], marker="o", ls="", color=COL["up"], ms=2.5, label="Higher expression (FDR < 0.05)"),
            Line2D([], [], marker="o", ls="", color=COL["down"], ms=2.5, label="Lower expression (FDR < 0.05)"),
            Line2D([], [], marker="o", ls="", color=COL["ns"], ms=2.5, label="Not significant"),
            Line2D([], [], marker="^", ls="", mfc="white", mec="#555", ms=3, mew=0.4, label="Beyond axis limits")]


def build(source=SOURCE, out_stem=OUT_STEM):
    ns.apply()
    de = load_de(source)
    fig, axs = plt.subplots(1, 3, figsize=(ns.DOUBLE, 64 * ns.MM), sharey=True)
    fig.subplots_adjust(left=0.065, right=0.99, top=0.84, bottom=0.24, wspace=0.08)
    manifest = []
    for j, (k, title) in enumerate(CONTRASTS.items()):
        manifest.append(draw_volcano(axs[j], de[de.Contrast == k], title))
        ns.panel_label(axs[j], "abc"[j], x=-0.04 if j else -0.16, y=1.13)
    axs[1].set_xlabel("log₂ fold change"); axs[0].set_ylabel(r"$\mathsf{−log₁₀}\ P$")
    fig.legend(handles=legend_handles(), frameon=False, loc="lower center", ncol=4, bbox_to_anchor=(0.53, 0.0),
               handletextpad=0.2, columnspacing=1.2)
    out_stem.parent.mkdir(parents=True, exist_ok=True); ns.save(fig, str(out_stem)); plt.close(fig)
    return pd.concat(manifest)[["contrast", "gene_label", "logFC", "P.Value", "adj.P.Val", "status"]]


if __name__ == "__main__":
    print(build())
