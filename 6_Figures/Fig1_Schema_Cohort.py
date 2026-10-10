#!/usr/bin/env python3
"""Build Figure 1 from the observed AI-CT schema partition and cohort flow."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import sys
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import nature_style as ns  # noqa: E402
FIGURE_DIR = ROOT / "figures"
# Aggregate Table 1 (group medians / percentages) written by build_table1_ai.py; no participant-level data.
CHARACTERISTICS_SOURCE = Path("path/to/natcomms_2026-09-11/tables/Table1_AIvisual_Phase2.csv")

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


# Characteristics shown in panel c: (section, label, Table 1 row, value kind).
# kind "pct" takes the percentage in parentheses; "med" takes the median.
PANEL_C_ROWS = (
    ("Demographics", "Age, years", "Age, years", "med"),
    ("Demographics", "Female, %", "Female, n (%)", "pct"),
    ("Demographics", "African American, %", "African American, n (%)", "pct"),
    ("Demographics", "Current smoker, %", "Current smoker, n (%)", "pct"),
    ("Demographics", "Pack-years", "Pack-years", "med"),
    ("Symptoms and function", "SGRQ total score", "SGRQ total score", "med"),
    ("Symptoms and function", "mMRC dyspnea ≥2, %", "mMRC dyspnea >=2, n (%)", "pct"),
    ("Symptoms and function", "Chronic bronchitis, %", "Chronic bronchitis, n (%)", "pct"),
    ("Symptoms and function", "6MWD, m", "6MWD, m", "med"),
    ("Spirometry", "FEV₁, % predicted", "FEV1 % predicted", "med"),
    ("Spirometry", "FEV₁/FVC", "FEV1/FVC", "med"),
    ("Spirometry", "GOLD 0, %", "GOLD 0", "pct"),
    ("Spirometry", "GOLD 1, %", "GOLD 1", "pct"),
    ("Spirometry", "GOLD 2, %", "GOLD 2", "pct"),
    ("Spirometry", "GOLD 3–4, %", "GOLD 3-4", "pct"),
    ("CT", "Emphysema, % LAA−950", "Emphysema % (LAA<-950 HU)", "med"),
    ("CT", "Pi10 airway wall, mm", "Pi10 airway wall, mm", "med"),
    ("Blood", "White blood cells, ×10⁹ per l", "White blood cells, x10^9/L", "med"),
)
GROUP_ORDER = (("No COPD", NO_COPD), ("Newly\ndiagnosed", NEW), ("Potentially\nexcluded", EXCLUDED),
               ("Established\nCOPD", ESTABLISHED))


def _cell_value(text: str, kind: str) -> float:
    text = str(text).replace(",", "").strip()
    if kind == "pct":
        return float(text.split("(")[1].split(")")[0])
    return float(text.split()[0])


def load_characteristics(source: Path) -> tuple[list[tuple[str, str]], np.ndarray]:
    table = pd.read_csv(source, dtype=str).set_index("Characteristic")
    group_cols = list(table.columns[:4])
    rows, values = [], []
    for section, label, key, kind in PANEL_C_ROWS:
        if key not in table.index:
            continue
        rows.append((section, label))
        values.append([_cell_value(table.loc[key, c], kind) for c in group_cols])
    return rows, np.array(values), group_cols


def draw_characteristics(ax: plt.Axes, source: Path) -> None:
    """Row-scaled heat map of Phase 2 characteristics; each cell prints the median or percentage."""
    rows, values, group_cols = load_characteristics(source)
    lo, hi = values.min(axis=1, keepdims=True), values.max(axis=1, keepdims=True)
    scaled = np.divide(values - lo, hi - lo, out=np.zeros_like(values), where=(hi - lo) > 0)
    cmap = plt.get_cmap("PuBu")
    n = len(rows)
    sections = [r[0] for r in rows]
    # Insert a small gap before each new section.
    ypos, y, prev = [], 0.0, None
    for sec in sections:
        if prev is not None and sec != prev:
            y += 0.45
        ypos.append(y); y += 1.0; prev = sec
    total = y
    for i, yy in enumerate(ypos):
        for j in range(4):
            v = values[i, j]
            ax.add_patch(Rectangle((j, yy), 1, 1, facecolor=cmap(0.08 + 0.72 * scaled[i, j]), edgecolor="white", lw=0.8))
            row_max = np.abs(values[i]).max()
            txt = f"{v:.2f}" if row_max < 3 else (f"{v:.1f}" if row_max < 100 else f"{v:.0f}")
            ax.text(j + 0.5, yy + 0.5, txt, ha="center", va="center", fontsize=5.5,
                    color="white" if scaled[i, j] > 0.62 else INK)
        ax.text(-0.08, yy + 0.5, rows[i][1], ha="right", va="center", fontsize=5.5, color=INK)
    prev = None
    for i, sec in enumerate(sections):
        if sec != prev:
            ax.text(-2.55, ypos[i] - 0.12, sec, ha="left", va="bottom", fontsize=5.5, fontweight="bold", color=MUTED)
            prev = sec
    for j, (name, color) in enumerate(GROUP_ORDER):
        n_txt = group_cols[j].split("(n=")[-1].rstrip(")") if "(n=" in group_cols[j] else ""
        ax.text(j + 0.5, -0.35, name, ha="center", va="bottom", fontsize=5.5, fontweight="bold", color=color,
                linespacing=1.0)
        ax.text(j + 0.5, -1.3, f"n = {n_txt}", ha="center", va="bottom", fontsize=5, color=MUTED)
    ax.set_xlim(-2.6, 4.05)
    ax.set_ylim(total + 0.1, -2.15)
    ax.axis("off")
    ax.text(4.0, total + 0.55, "Shading scales each row from lowest (light) to highest (dark) group value",
            ha="right", va="top", fontsize=5, color=MUTED, fontstyle="italic")


def build_figure() -> None:
    parts = pd.read_csv(PARTITION_SOURCE, index_col=0)
    row = parts.loc["DisCOPDAI (AI-CT)"]
    counts = {
        "No COPD": int(row["Grp_0 No COPD"]),
        "Newly diagnosed": int(row["Grp_1 New COPD"]),
        "Potentially excluded": int(row["Grp_2 Excluded"]),
        "Established COPD": int(row["Grp_3 Existing COPD"]),
    }
    MM = 1 / 25.4
    with plt.rc_context(
        {"font.family": ["Arial", "Arial Unicode MS", "Liberation Sans", "DejaVu Sans"],
         "font.size": 6, "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42}
    ):
        fig = plt.figure(figsize=(183 * MM, 140 * MM), facecolor="white")
        grid = fig.add_gridspec(2, 2, height_ratios=(1.55, 0.62), width_ratios=(0.43, 0.57),
                                left=0.01, right=0.99, top=0.975, bottom=0.045, hspace=0.10, wspace=0.04)
        ax_a = fig.add_subplot(grid[0, 0])
        ax_c = fig.add_subplot(grid[0, 1])
        ax_b = fig.add_subplot(grid[1, :])
        for ax in (ax_a, ax_b):
            ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

        # Panel a: 2 x 2 definition of bidirectional reclassification.
        ax_a.text(0.0, 0.985, "a", fontsize=8, fontweight="bold", color=INK, va="top")
        x0, xm, x1 = 0.17, 0.585, 1.0
        y0, ym, y1 = 0.14, 0.50, 0.86
        ax_a.text((x0 + x1) / 2, 0.955, "Multidimensional COPD diagnosis", ha="center", fontsize=6.5,
                  fontweight="bold", color=INK)
        ax_a.text((x0 + xm) / 2, 0.895, "No", ha="center", fontsize=6, color=MUTED)
        ax_a.text((xm + x1) / 2, 0.895, "Yes", ha="center", fontsize=6, color=MUTED)
        ax_a.text(0.055, (y0 + y1) / 2, "Spirometric airflow\nobstruction", ha="center", va="center", rotation=90,
                  fontsize=6, fontweight="bold", color=INK, linespacing=1.05)
        ax_a.text(0.155, (ym + y1) / 2, "No", ha="right", va="center", fontsize=6, color=MUTED)
        ax_a.text(0.155, (y0 + ym) / 2, "Yes", ha="right", va="center", fontsize=6, color=MUTED)
        add_cell(ax_a, xy=(x0, ym), width=xm - x0, height=y1 - ym, color=NO_COPD, group="No COPD",
                 count=counts["No COPD"])
        add_cell(ax_a, xy=(xm, ym), width=x1 - xm, height=y1 - ym, color=NEW, group="Newly\ndiagnosed COPD",
                 count=counts["Newly diagnosed"], reclassified=True)
        add_cell(ax_a, xy=(x0, y0), width=xm - x0, height=ym - y0, color=EXCLUDED, group="Potentially\nexcluded",
                 count=counts["Potentially excluded"], reclassified=True)
        add_cell(ax_a, xy=(xm, y0), width=x1 - xm, height=ym - y0, color=ESTABLISHED, group="Established COPD",
                 count=counts["Established COPD"])
        ax_a.text((x0 + x1) / 2, 0.055, "Groups assigned at Phase 1 (AI visual CT);\ncharacteristics at the Phase 2 blood draw",
                  ha="center", va="center", fontsize=5.5, color=MUTED, linespacing=1.1)

        # Panel c: Phase 2 participant characteristics by group.
        ax_c.text(-0.005, 0.985, "c", transform=ax_c.transAxes, fontsize=8, fontweight="bold", color=INK, va="top")
        draw_characteristics(ax_c, CHARACTERISTICS_SOURCE)

        # Panel b: analytic-cohort flow.
        ax_b.text(0.0, 1.0, "b", fontsize=8, fontweight="bold", color=INK, va="top")
        main_y, main_h, main_w = 0.56, 0.36, 0.19
        main_x = (0.03, 0.285, 0.54, 0.795)
        main = (("5,113", ("RNA-seq masterfile", "records")), ("4,146", ("Final-analyzed", "RNA-seq samples")),
                ("4,142", ("Linked visit-2 phenotype", "and schema data")), ("3,309", ("Complete-case", "analytic cohort")))
        for x, (count, lines) in zip(main_x, main):
            add_flow_box(ax_b, xy=(x, main_y), width=main_w, height=main_h, count=count, lines=lines)
        for left, right in zip(main_x[:-1], main_x[1:]):
            ax_b.add_patch(FancyArrowPatch((left + main_w, main_y + main_h / 2), (right, main_y + main_h / 2),
                                           arrowstyle="-|>", mutation_scale=7, linewidth=0.7, color=MUTED))
        exclusions = (("−967", ("Not marked", "final analyzed")), ("−4", ("No linked visit-2", "phenotype")),
                      ("−833", ("Incomplete classification", "or model covariates")))
        for i, (count, lines) in enumerate(exclusions):
            branch_x = (main_x[i] + main_w + main_x[i + 1]) / 2
            ax_b.plot([branch_x, branch_x], [main_y + main_h / 2, 0.44], color=MUTED, lw=0.6)
            ax_b.add_patch(FancyArrowPatch((branch_x, 0.44), (branch_x, 0.40), arrowstyle="-|>", mutation_scale=6,
                                           linewidth=0.6, color=MUTED))
            add_flow_box(ax_b, xy=(branch_x - 0.095, 0.08), width=0.19, height=0.31, count=count, lines=lines,
                         fill="#F1F2F3", edge=GRID, count_size=6.5, text_size=5.5)
        fig.text(0.015, 0.012, "Primary model adjusted for age, race, sex and smoking status; secondary model "
                 "additionally adjusted for leukocyte percentages.", fontsize=5.5, color=MUTED, fontstyle="italic")

        FIGURE_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIGURE_DIR / "Fig1_schema_cohort.pdf", facecolor="white",
                    metadata={"Title": "Multidimensional COPD schema, analytic cohort and participant characteristics"})
        fig.savefig(FIGURE_DIR / "Fig1_schema_cohort.png", dpi=450, facecolor="white")
        ns.export_panels(fig, str(FIGURE_DIR / "Fig1_schema_cohort"))
        plt.close(fig)


if __name__ == "__main__":
    build_figure()
