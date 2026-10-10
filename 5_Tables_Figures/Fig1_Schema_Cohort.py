#!/usr/bin/env python3
"""Build Figure 1 from the observed AI-CT schema partition and cohort flow."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = ROOT / "figures"
PARTITION_SOURCE = Path(
    "path/to/"
    "Schema_omics_csv_bhatt/reproduce_2026-09-01/"
    "figures_from_review/partition_groups_2026.05.16.csv"
)

INK = "#202124"
MUTED = "#5F6368"
GRID = "#B8BDC3"
NO_COPD = "#8D969D"
NEW = "#0072B2"
EXCLUDED = "#D55E00"
ESTABLISHED = "#3B7C44"


def add_cell(
    ax: plt.Axes,
    *,
    xy: tuple[float, float],
    width: float,
    height: float,
    color: str,
    group: str,
    count: int,
    reclassified: bool = False,
) -> None:
    ax.add_patch(
        Rectangle(
            xy,
            width,
            height,
            facecolor=color,
            edgecolor=color,
            alpha=0.12,
            linewidth=0,
        )
    )
    ax.add_patch(
        Rectangle(
            xy,
            width,
            height,
            facecolor="none",
            edgecolor=color if reclassified else GRID,
            linewidth=1.6 if reclassified else 0.8,
        )
    )
    x, y = xy
    cx = x + width / 2
    cy = y + height / 2
    ax.text(
        cx,
        cy + 0.040,
        group,
        ha="center",
        va="center",
        fontsize=6.5,
        fontweight="bold",
        color=INK,
        linespacing=1.05,
    )
    ax.text(
        cx,
        cy - 0.050,
        f"n = {count:,}",
        ha="center",
        va="center",
        fontsize=6,
        color=INK,
    )
    if reclassified:
        ax.text(
            cx,
            y + 0.035,
            "RECLASSIFIED",
            ha="center",
            va="bottom",
            fontsize=5.5,
            fontweight="bold",
            color=color,
        )


def add_flow_box(
    ax: plt.Axes,
    *,
    xy: tuple[float, float],
    width: float,
    height: float,
    count: str,
    lines: tuple[str, ...],
    fill: str = "white",
    edge: str = INK,
    count_size: float = 7,
    text_size: float = 6,
) -> None:
    x, y = xy
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.008,rounding_size=0.012",
            facecolor=fill,
            edgecolor=edge,
            linewidth=0.9,
        )
    )
    ax.text(
        x + width / 2,
        y + height * 0.68,
        count,
        ha="center",
        va="center",
        fontsize=count_size,
        fontweight="bold",
        color=INK,
    )
    ax.text(
        x + width / 2,
        y + height * 0.36,
        "\n".join(lines),
        ha="center",
        va="center",
        fontsize=text_size,
        fontweight="bold" if fill == "white" else "normal",
        color=INK,
        linespacing=1.05,
    )


def build_figure() -> None:
    parts = pd.read_csv(PARTITION_SOURCE, index_col=0)
    row = parts.loc["DisCOPDAI (AI-CT)"]
    total = int(row["Total"])
    counts = {
        "No COPD": int(row["Grp_0 No COPD"]),
        "Newly diagnosed": int(row["Grp_1 New COPD"]),
        "Potentially excluded": int(row["Grp_2 Excluded"]),
        "Established COPD": int(row["Grp_3 Existing COPD"]),
    }

    with plt.rc_context(
        {
            "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
            "font.size": 6,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    ):
        fig = plt.figure(figsize=(7.2, 5.0), facecolor="white")
        grid = fig.add_gridspec(
            2,
            height_ratios=(1.08, 0.72),
            left=0.065,
            right=0.985,
            top=0.965,
            bottom=0.09,
            hspace=0.16,
        )
        ax_a = fig.add_subplot(grid[0])
        ax_b = fig.add_subplot(grid[1])

        for ax in (ax_a, ax_b):
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis("off")

        # Panel A: a direct 2 x 2 definition of bidirectional reclassification.
        ax_a.text(-0.01, 1.01, "a", fontsize=8, fontweight="bold", color=INK, clip_on=False)
        x0, xm, x1 = 0.30, 0.50, 0.70
        y0, ym, y1 = 0.14, 0.48, 0.82
        ax_a.text(
            (x0 + x1) / 2,
            0.955,
            "Multidimensional COPD diagnosis",
            ha="center",
            fontsize=6.5,
            fontweight="bold",
            color=INK,
        )
        ax_a.text((x0 + xm) / 2, 0.875, "No", ha="center", fontsize=6, color=MUTED)
        ax_a.text((xm + x1) / 2, 0.875, "Yes", ha="center", fontsize=6, color=MUTED)
        ax_a.text(
            0.205,
            (y0 + y1) / 2,
            "Spirometric airflow\nobstruction",
            ha="center",
            va="center",
            rotation=90,
            fontsize=6,
            fontweight="bold",
            color=INK,
            linespacing=1.05,
        )
        ax_a.text(0.27, (ym + y1) / 2, "No", ha="right", va="center", fontsize=6, color=MUTED)
        ax_a.text(0.27, (y0 + ym) / 2, "Yes", ha="right", va="center", fontsize=6, color=MUTED)

        add_cell(
            ax_a,
            xy=(x0, ym),
            width=xm - x0,
            height=y1 - ym,
            color=NO_COPD,
            group="No COPD",
            count=counts["No COPD"],
        )
        add_cell(
            ax_a,
            xy=(xm, ym),
            width=x1 - xm,
            height=y1 - ym,
            color=NEW,
            group="Newly\ndiagnosed COPD",
            count=counts["Newly diagnosed"],
            reclassified=True,
        )
        add_cell(
            ax_a,
            xy=(x0, y0),
            width=xm - x0,
            height=ym - y0,
            color=EXCLUDED,
            group="Potentially\nexcluded",
            count=counts["Potentially excluded"],
            reclassified=True,
        )
        add_cell(
            ax_a,
            xy=(xm, y0),
            width=x1 - xm,
            height=ym - y0,
            color=ESTABLISHED,
            group="Established COPD",
            count=counts["Established COPD"],
        )
        # Panel B: observed analytic-cohort flow reconstructed from the source pipeline.
        ax_b.text(-0.005, 0.97, "b", fontsize=8, fontweight="bold", color=INK, clip_on=False)
        main_y, main_h, main_w = 0.57, 0.27, 0.19
        main_x = (0.015, 0.275, 0.535, 0.795)
        main = (
            ("5,113", ("RNA-seq masterfile", "records")),
            ("4,146", ("Final-analyzed", "RNA-seq samples")),
            ("4,142", ("Linked visit-2 phenotype", "and schema data")),
            ("3,309", ("Complete-case", "analytic cohort")),
        )
        for x, (count, lines) in zip(main_x, main):
            add_flow_box(
                ax_b,
                xy=(x, main_y),
                width=main_w,
                height=main_h,
                count=count,
                lines=lines,
            )

        for left, right in zip(main_x[:-1], main_x[1:]):
            ax_b.add_patch(
                FancyArrowPatch(
                    (left + main_w, main_y + main_h / 2),
                    (right, main_y + main_h / 2),
                    arrowstyle="-|>",
                    mutation_scale=8,
                    linewidth=0.8,
                    color=MUTED,
                )
            )

        exclusion_x = (0.145, 0.405, 0.665)
        exclusions = (
            ("−967", ("Not marked", "final analyzed")),
            ("−4", ("No linked visit-2", "phenotype")),
            ("−833", ("Incomplete classification", "or model covariates")),
        )
        for i, (x, (count, lines)) in enumerate(zip(exclusion_x, exclusions)):
            branch_x = (main_x[i] + main_w + main_x[i + 1]) / 2
            ax_b.plot([branch_x, branch_x], [main_y + main_h / 2, 0.43], color=MUTED, lw=0.7)
            ax_b.add_patch(
                FancyArrowPatch(
                    (branch_x, 0.43),
                    (branch_x, 0.39),
                    arrowstyle="-|>",
                    mutation_scale=7,
                    linewidth=0.7,
                    color=MUTED,
                )
            )
            add_flow_box(
                ax_b,
                xy=(x, 0.12),
                width=0.19,
                height=0.24,
                count=count,
                lines=lines,
                fill="#F1F2F3",
                edge=GRID,
                count_size=7,
                text_size=5.5,
            )

        fig.text(
            0.065,
            0.025,
            "Primary model adjusted for age, race, sex and smoking status. Secondary model additionally adjusted for leukocyte percentages (Supplementary Information).",
            fontsize=6,
            color="black",
            fontstyle="italic",
        )

        FIGURE_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            FIGURE_DIR / "Fig1_analysis_lineage.pdf",
            facecolor="white",
            metadata={"Title": "Multidimensional COPD Schema and RNA-seq Analytic Cohort in COPDGene"},
        )
        fig.savefig(
            FIGURE_DIR / "Fig1_analysis_lineage.png",
            dpi=600,
            facecolor="white",
        )
        plt.close(fig)


if __name__ == "__main__":
    build_figure()

