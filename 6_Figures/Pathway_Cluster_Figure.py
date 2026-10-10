#!/usr/bin/env python3
"""Render the objective Model 1 pathway-profile cluster figure."""

from __future__ import annotations

from pathlib import Path
import re
import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"

CONTRASTS = ["Newly diagnosed", "Potentially excluded", "Established COPD"]
COLORS = {
    "Newly diagnosed": "#0072B2",
    "Potentially excluded": "#D55E00",
    "Established COPD": "#4D4D4D",
}


def clean_pathway(value: str) -> str:
    value = re.sub(r"^(HALLMARK|REACTOME)_", "", value)
    value = value.replace("_", " ").title()
    replacements = {
        "Dna": "DNA", "Rna": "RNA", "Mrna": "mRNA", "Rrna": "rRNA",
        "Nfkb": "NF-kB", "Tnfa": "TNF-alpha", "Mtorc1": "mTORC1",
        "Ecm": "ECM", "P53": "p53", "Myc": "MYC", "Iii": "III",
        "Ca2 Pathway": "Ca2+ Pathway",
        "Asparagine N Linked": "Asparagine N-Linked",
        "Mapk1 Mapk3": "MAPK1/MAPK3", "Snrnp": "snRNP",
        "Trna": "tRNA", "Hox": "HOX", "Sirt1": "SIRT1",
        "Interferon Alpha Beta": "Interferon-alpha/beta",
    }
    for old, new in replacements.items():
        value = value.replace(old, new)
    return value


DISPLAY_PATHWAY_LABELS = {
    "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION":
        "Neurotransmitter receptor/postsynaptic signaling",
    "REACTOME_MAPK1_MAPK3_SIGNALING": "MAPK1/MAPK3 signaling",
    "REACTOME_CA2_PATHWAY": "Ca2+ pathway",
    "REACTOME_SNRNP_ASSEMBLY": "snRNP assembly",
    "HALLMARK_MYC_TARGETS_V2": "MYC targets V2",
    "REACTOME_TRNA_MODIFICATION_IN_THE_NUCLEUS_AND_CYTOSOL":
        "tRNA modification (nucleus/cytosol)",
    "REACTOME_ACTIVATION_OF_ANTERIOR_HOX_GENES_IN_HINDBRAIN_DEVELOPMENT_DURING_EARLY_EMBRYOGENESIS":
        "Anterior HOX activation in hindbrain",
    "REACTOME_SIRT1_NEGATIVELY_REGULATES_RRNA_EXPRESSION":
        "SIRT1 repression of rRNA expression",
    "HALLMARK_ADIPOGENESIS": "Adipogenesis",
    "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING": "Interferon-alpha/beta signaling",
    "REACTOME_ASPARAGINE_N_LINKED_GLYCOSYLATION":
        "Asparagine N-linked glycosylation",
    "REACTOME_DISEASES_OF_SIGNAL_TRANSDUCTION_BY_GROWTH_FACTOR_RECEPTORS_AND_SECOND_MESSENGERS":
        "Growth-factor/second-messenger signaling",
    "REACTOME_TRANSCRIPTIONAL_REGULATION_OF_BROWN_AND_BEIGE_ADIPOCYTE_DIFFERENTIATION":
        "Brown/beige adipocyte differentiation",
    "REACTOME_HCMV_EARLY_EVENTS": "HCMV early events",
    "REACTOME_INTERFERON_SIGNALING": "Interferon signaling",
}


def display_pathway(value: str) -> str:
    label = DISPLAY_PATHWAY_LABELS.get(value, clean_pathway(value))
    return "\n".join(textwrap.wrap(label, width=43, break_long_words=False))


def main() -> None:
    summary = pd.read_csv(TABLES / "pathway_profile_cluster_summary_model1.csv")
    membership = pd.read_csv(TABLES / "pathway_profile_cluster_membership_model1.csv")
    reps = pd.read_csv(TABLES / "pathway_profile_cluster_representatives_model1.csv")
    cluster_drivers = pd.read_csv(TABLES / "pathway_profile_cluster_drivers_model1.csv")
    shared_genes = pd.read_csv(TABLES / "pathway_profile_shared_significant_genes_model1.csv")
    clusters = sorted(summary["cluster_number"].unique())
    median_matrix = (
        summary.pivot(index="cluster_number", columns="contrast_label", values="median_NES")
        .loc[clusters, CONTRASTS]
    )
    significant_matrix = (
        summary.pivot(index="cluster_number", columns="contrast_label", values="fdr_significant_pathways")
        .loc[clusters, CONTRASTS]
    )
    sizes = membership.groupby("cluster_number").size().reindex(clusters)

    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })

    fig = plt.figure(figsize=(11.0, 8.4), constrained_layout=False, facecolor="white")
    outer = fig.add_gridspec(
        2, 1, height_ratios=[1.38, 1.08],
        left=0.055, right=0.975, top=0.945, bottom=0.105, hspace=0.27,
    )
    ax_table = fig.add_subplot(outer[0])
    ax_table.set_xlim(0, 1)
    ax_table.set_ylim(0, 1)
    ax_table.axis("off")
    limit = float(np.ceil(np.nanmax(np.abs(median_matrix.to_numpy())) * 2) / 2)
    nes_cmap = LinearSegmentedColormap.from_list(
        "nes_diverging", ["#2166AC", "#F7F7F7", "#B2182B"]
    )
    nes_norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    profile_x = np.array([0.115, 0.225, 0.335])
    profile_left, profile_right = 0.085, 0.385
    tile_width, tile_height = 0.094, 0.092
    representative_role_x, representative_text_x = 0.405, 0.475
    driver_text_x = 0.805
    row_top, row_height = 0.79, 0.135

    ax_table.text(0.0, 1.045, "A", fontsize=12, fontweight="bold", va="top", clip_on=False)
    ax_table.text(
        0.045, 1.045, "Pathway-cluster profiles and representative programs",
        fontsize=10.2, fontweight="bold", va="top", clip_on=False,
    )
    ax_table.plot([0.0, 1.0], [0.985, 0.985], color="#1B1B1B", lw=1.05)
    ax_table.text((profile_left + profile_right) / 2, 0.947, "Median FGSEA NES vs no COPD", ha="center", va="center", fontsize=7.5, fontweight="bold")
    ax_table.text(0.69, 0.947, "Cluster annotation", ha="center", va="center", fontsize=7.5, fontweight="bold")
    ax_table.plot([profile_left, profile_right], [0.918, 0.918], color="#777777", lw=0.55)
    ax_table.plot([0.405, 0.975], [0.918, 0.918], color="#777777", lw=0.55)
    ax_table.text(0.035, 0.855, "Cluster", ha="center", va="center", fontsize=7.1, fontweight="bold")
    ax_table.text(0.405, 0.855, "Selection\nrole", ha="left", va="center", fontsize=6.6, fontweight="bold", linespacing=0.95)
    ax_table.text(0.475, 0.855, "Representative pathways", ha="left", va="center", fontsize=7.1, fontweight="bold")
    ax_table.text(driver_text_x, 0.855, "Leading-edge drivers", ha="left", va="center", fontsize=7.1, fontweight="bold")
    for x_pos, contrast in zip(profile_x, CONTRASTS):
        ax_table.text(
            x_pos, 0.855, contrast.replace(" ", "\n", 1),
            ha="center", va="center", fontsize=7.1, fontweight="bold",
            color="#222222", linespacing=1.0,
        )
    ax_table.text(
        profile_x[1], 0.803, "0 genes at gene-level FDR < .05",
        ha="center", va="center", fontsize=5.9, color=COLORS["Potentially excluded"],
    )
    ax_table.plot([0.0, 1.0], [0.79, 0.79], color="#1B1B1B", lw=0.75)
    for row, cluster in enumerate(clusters):
        row_center = row_top - row * row_height - row_height / 2
        row_bottom = row_top - (row + 1) * row_height
        cluster_reps = reps.loc[reps["cluster_number"] == cluster].sort_values("representative_rank")
        profile_rep = cluster_reps.loc[cluster_reps["representative_role"] == "Profile exemplar"].iloc[0]
        evidence_rep = cluster_reps.loc[cluster_reps["representative_role"] == "Central evidence"].iloc[0]
        genes = (
            cluster_drivers.loc[cluster_drivers["cluster_number"] == cluster]
            .sort_values("driver_rank")["gene"].tolist()
        )

        ax_table.text(
            0.035, row_center, f"C{cluster}\n(n={sizes.loc[cluster]})",
            ha="center", va="center", fontsize=8.7, fontweight="bold", linespacing=1.0,
        )
        values = median_matrix.loc[cluster, CONTRASTS].to_numpy(dtype=float)
        significant_counts = significant_matrix.loc[cluster, CONTRASTS].to_numpy(dtype=int)
        for x_pos, value, count in zip(profile_x, values, significant_counts):
            fill = nes_cmap(nes_norm(value))
            text_color = "white" if abs(value) >= 1.05 else "#222222"
            ax_table.add_patch(Rectangle(
                (x_pos - tile_width / 2, row_center - tile_height / 2),
                tile_width, tile_height, facecolor=fill, edgecolor="white", linewidth=0.8,
            ))
            ax_table.text(
                x_pos, row_center + 0.014, f"{value:+.2f}",
                ha="center", va="center", fontsize=7.5, fontweight="bold", color=text_color,
            )
            ax_table.text(
                x_pos, row_center - 0.022, f"{count}/{sizes.loc[cluster]} FDR",
                ha="center", va="center", fontsize=5.5, color=text_color,
            )

        source_profile = "H" if profile_rep["source"] == "Hallmark" else "R"
        source_evidence = "H" if evidence_rep["source"] == "Hallmark" else "R"
        representative_rows = [
            (f"PROFILE [{source_profile}]", display_pathway(profile_rep["pathway"]), row_center + 0.032),
            (f"EVIDENCE [{source_evidence}]", display_pathway(evidence_rep["pathway"]), row_center - 0.031),
        ]
        for role, value, y_pos in representative_rows:
            ax_table.text(
                representative_role_x, y_pos, role,
                ha="left", va="center", fontsize=6.6,
                fontstyle="italic", color="#111111",
            )
            ax_table.text(
                representative_text_x, y_pos, value,
                ha="left", va="center", fontsize=7.3, color="#222222",
                linespacing=0.92,
            )
        ax_table.text(
            representative_text_x, row_center + 0.005,
            f"standardized-profile distance {profile_rep['profile_distance_to_centroid']:.2f}",
            ha="left", va="center", fontsize=5.9, color="#666666",
        )
        ax_table.text(
            representative_text_x, row_center - 0.054,
            f"aggregate FDR evidence {evidence_rep['aggregate_evidence']:.1f}",
            ha="left", va="center", fontsize=5.7, color="#666666",
        )

        for gene_row, gene_slice in enumerate((genes[:3], genes[3:6])):
            ranked = " · ".join(gene_slice)
            ax_table.text(
                driver_text_x, row_center + 0.023 - gene_row * 0.046, ranked,
                ha="left", va="center", fontsize=7.6, fontstyle="italic", color="#222222",
            )
        if row < len(clusters) - 1:
            ax_table.plot([0.0, 1.0], [row_bottom, row_bottom], color="#DDDDDD", lw=0.45)

    table_bottom = row_top - len(clusters) * row_height
    ax_table.plot([0.0, 1.0], [table_bottom, table_bottom], color="#1B1B1B", lw=0.9)
    cax = ax_table.inset_axes([0.085, 0.036, 0.20, 0.026])
    colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=nes_norm, cmap=nes_cmap), cax=cax, orientation="horizontal")
    colorbar.set_ticks([-limit, 0, limit])
    colorbar.ax.tick_params(labelsize=5.5, length=2, pad=1)
    colorbar.outline.set_linewidth(0.5)
    ax_table.text(
        0.0, 0.049, "Median FGSEA NES\n<0: lower expression\n>0: higher expression",
        ha="left", va="center", fontsize=5.2, color="#333333", linespacing=0.95,
    )
    ax_table.text(
        0.405, 0.075,
        "PROFILE = pathway nearest the cluster centroid; EVIDENCE = highest FDR score among the nearest 25%.",
        ha="left", va="center", fontsize=5.1, color="#111111", fontstyle="italic",
    )
    ax_table.text(
        0.405, 0.052,
        "FDR score = Σ min[−log10(FDR), 10] across 3 contrasts; PROFILE is excluded from EVIDENCE selection.",
        ha="left", va="center", fontsize=5.1, color="#111111", fontstyle="italic",
    )
    ax_table.text(
        0.405, 0.029,
        "Cells show median original FGSEA NES; lower text gives pathways at FDR < .05.",
        ha="left", va="center", fontsize=5.1, color="#111111", fontstyle="italic",
    )
    ax_table.text(
        0.405, 0.006,
        "[H], Hallmark; [R], Reactome. Drivers are cluster-specific leading-edge genes.",
        ha="left", va="center", fontsize=5.1, color="#111111", fontstyle="italic",
    )

    bottom = outer[1].subgridspec(2, 3, height_ratios=[0.24, 0.76], hspace=0.03, wspace=0.46)
    ax_bottom_title = fig.add_subplot(bottom[0, :])
    ax_bottom_title.axis("off")
    ax_bottom_title.text(0.0, 0.96, "B", fontsize=12, fontweight="bold", va="top")
    ax_bottom_title.text(
        0.045, 0.96,
        "Selected shared leading-edge genes across all three contrasts",
        fontsize=9, fontweight="bold", va="top",
    )
    axes = [fig.add_subplot(bottom[1, i]) for i in range(3)]
    t_limit = float(np.ceil(shared_genes["moderated_t"].abs().max()))
    for ax, contrast in zip(axes, CONTRASTS):
        data = shared_genes.loc[shared_genes["contrast_label"] == contrast].copy()
        data = data.sort_values("display_order", ascending=False).reset_index(drop=True)
        y = np.arange(len(data))
        color = COLORS[contrast]
        ax.axvline(0, color="#737373", lw=0.7, zorder=0)
        ax.axhline(5.5, color="#B8B8B8", lw=0.55, ls=(0, (3, 2)), zorder=0)
        ax.hlines(y, 0, data["moderated_t"], color=color, lw=1.0, alpha=0.72, zorder=1)
        significant = data["gene_fdr"] < 0.05
        ax.scatter(
            data.loc[~significant, "moderated_t"], y[~significant],
            s=21, facecolor="white", edgecolor=color, linewidth=0.9, zorder=3,
        )
        ax.scatter(
            data.loc[significant, "moderated_t"], y[significant],
            s=21, facecolor=color, edgecolor=color, linewidth=0.9, zorder=3,
        )
        ax.set_yticks(y, data["gene"])
        ax.set_xlim(-t_limit, t_limit)
        ax.grid(axis="x", color="#E6E6E6", lw=0.45)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(
            f"{contrast}\n{int(significant.sum())} of {len(data)} genes at FDR < .05",
            color=color, fontweight="bold", fontsize=7.5, pad=6, linespacing=1.15,
        )
        ax.set_xlabel("Moderated t statistic vs no COPD")

    legend = [
        Line2D([0], [0], marker="o", color="#222222", markerfacecolor="#222222", lw=0, markersize=5, label="Gene FDR < .05"),
        Line2D([0], [0], marker="o", color="#222222", markerfacecolor="white", lw=0, markersize=5, label="Gene FDR ≥ .05"),
    ]
    fig.legend(handles=legend, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.025), fontsize=7)

    FIGURES.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf", "svg"):
        path = FIGURES / f"Fig5_model1_pathway_clusters.{suffix}"
        fig.savefig(path, dpi=600 if suffix == "png" else None, facecolor="white")
    plt.close(fig)
    print(f"Wrote Figure 5 as PNG, PDF, and SVG in {FIGURES}")


if __name__ == "__main__":
    main()

