#!/usr/bin/env python3
"""Build the primary Model 1 effect-profile geometry figure with Matplotlib.

This script visualizes existing exported limma and scaling estimates. It does
not refit the expression model, λ estimator, or bootstrap.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
TABLES = ROOT / "tables"

DE_SOURCE = Path(
    "path/to/"
    "Schema_omics_csv_bhatt/copdschema/extracted/"
    "gene_de_Model_1_DisCOPDAI_allresults_2026.04.25.csv"
)
BOOTSTRAP_SOURCE = Path(
    "path/to/"
    "Schema_omics_csv_bhatt/reproduce_2026-09-01/outputs/tables/"
    "bootstrap_lambda.csv"
)

OUTPUT_STEM = FIGURES / "Fig2_model1_axis"
PROVENANCE_FILE = TABLES / "Fig2_model1_axis_provenance.json"

ESTABLISHED = "DisCOPDAI_3vs0"
NEW = "DisCOPDAI_1vs0"
EXCLUDED = "DisCOPDAI_2vs0"
GRID_SIZE = 55

GROUP_COLORS = {
    "Newly diagnosed": "#0072B2",
    "Potentially excluded": "#D55E00",
}


def load_effect_profiles() -> pd.DataFrame:
    de = pd.read_csv(DE_SOURCE)
    required = {"id", "gene_label", "logFC", "adj.P.Val", "Contrast"}
    missing = required.difference(de.columns)
    if missing:
        raise ValueError(f"Missing required DE columns: {sorted(missing)}")

    selected = de.loc[
        de["Contrast"].isin([NEW, EXCLUDED, ESTABLISHED]),
        ["id", "gene_label", "logFC", "adj.P.Val", "Contrast"],
    ]
    effects = selected.pivot(index="id", columns="Contrast", values="logFC")
    fdr = selected.pivot(index="id", columns="Contrast", values="adj.P.Val")
    labels = selected.drop_duplicates("id").set_index("id")["gene_label"]

    axis_genes = fdr[ESTABLISHED] < 0.05
    profiles = effects.loc[axis_genes, [ESTABLISHED, NEW, EXCLUDED]].dropna().copy()
    profiles["gene_label"] = labels.reindex(profiles.index)
    return profiles


def panel_label(ax: mpl.axes.Axes, label: str, y: float = 1.08) -> None:
    ax.text(
        -0.18,
        y,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
        ha="left",
    )


def style_axis(ax: mpl.axes.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3, width=0.6)
    ax.grid(False)


def draw_density_panel(
    ax: mpl.axes.Axes,
    profiles: pd.DataFrame,
    outcome: str,
    group: str,
    lam: float,
    ci_lo: float,
    ci_hi: float,
    limit: float,
) -> mpl.collections.PolyCollection:
    x = profiles[ESTABLISHED].to_numpy()
    y = profiles[outcome].to_numpy()
    hb = ax.hexbin(
        x,
        y,
        gridsize=GRID_SIZE,
        extent=(-limit, limit, -limit, limit),
        mincnt=1,
        cmap="viridis",
        norm=LogNorm(vmin=1),
        linewidths=0,
        rasterized=True,
    )

    identity_x = np.array([-limit, limit])
    observed_x = np.array([np.min(x), np.max(x)])
    ax.plot(identity_x, identity_x, color="0.55", lw=0.8, ls=(0, (4, 3)), label="Equal effect")
    ax.plot(
        observed_x,
        lam * observed_x,
        color=GROUP_COLORS[group],
        lw=1.7,
        label="Corrected disease-axis slope",
    )
    ax.axhline(0, color="0.82", lw=0.5, zorder=0)
    ax.axvline(0, color="0.82", lw=0.5, zorder=0)
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Established COPD vs no COPD\ngene-expression $\\log_{2}$ fold change")
    ax.set_ylabel(f"{group} vs no COPD\ngene-expression $\\log_{{2}}$ fold change")
    ax.set_title(
        f"{group} effect profile",
        loc="left",
        y=1.08,
        pad=0,
        fontsize=9.5,
        fontweight="bold",
    )
    ax.text(
        0.04,
        0.96,
        f"λ = {lam:.2f} ({ci_lo:.2f}-{ci_hi:.2f})\n"
        f"n = {len(profiles):,} axis genes",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=7.3,
        color="0.12",
    )
    style_axis(ax)
    return hb


def build_figure() -> dict[str, object]:
    FIGURES.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)

    profiles = load_effect_profiles()
    bootstrap = pd.read_csv(BOOTSTRAP_SOURCE).set_index("group")

    required_groups = {"Newly diagnosed", "Potentially excluded"}
    if not required_groups.issubset(bootstrap.index):
        raise ValueError("Bootstrap table is missing a primary reclassified group")

    finite_values = profiles[[ESTABLISHED, NEW, EXCLUDED]].to_numpy().ravel()
    limit = np.ceil(np.nanmax(np.abs(finite_values)) * 10) / 10

    style = {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    }

    with mpl.rc_context(style):
        fig = plt.figure(figsize=(7.25, 3.55))
        upper_grid = fig.add_gridspec(
            1,
            3,
            width_ratios=[1, 1, 0.05],
            left=0.13,
            right=0.96,
            bottom=0.20,
            top=0.85,
            wspace=0.58,
        )
        axes = {
            "A": fig.add_subplot(upper_grid[0, 0]),
            "B": fig.add_subplot(upper_grid[0, 1]),
        }
        colorbar_axis = fig.add_subplot(upper_grid[0, 2])

        hb_a = draw_density_panel(
            axes["A"],
            profiles,
            NEW,
            "Newly diagnosed",
            float(bootstrap.loc["Newly diagnosed", "lambda"]),
            float(bootstrap.loc["Newly diagnosed", "ci_lo"]),
            float(bootstrap.loc["Newly diagnosed", "ci_hi"]),
            limit,
        )
        hb_b = draw_density_panel(
            axes["B"],
            profiles,
            EXCLUDED,
            "Potentially excluded",
            float(bootstrap.loc["Potentially excluded", "lambda"]),
            float(bootstrap.loc["Potentially excluded", "ci_lo"]),
            float(bootstrap.loc["Potentially excluded", "ci_hi"]),
            limit,
        )
        density_max = max(float(np.max(hb_a.get_array())), float(np.max(hb_b.get_array())))
        hb_a.set_clim(1, density_max)
        hb_b.set_clim(1, density_max)
        colorbar = fig.colorbar(hb_b, cax=colorbar_axis)
        colorbar.ax.yaxis.set_label_position("left")
        colorbar.set_label("Genes per hexagon (log scale)", fontsize=7.5, labelpad=9)
        colorbar.ax.tick_params(labelsize=6.5, length=2)

        axes["A"].legend(
            loc="lower right",
            frameon=False,
            handlelength=2.6,
            borderaxespad=0.4,
        )

        for label, key in zip("AB", ["A", "B"], strict=True):
            panel_label(axes[key], label, y=1.13)

        fig.align_xlabels([axes["A"], axes["B"]])

        for extension, kwargs in {
            "pdf": {},
            "svg": {},
            "png": {"dpi": 600},
        }.items():
            fig.savefig(OUTPUT_STEM.with_suffix(f".{extension}"), **kwargs)
        plt.close(fig)

    provenance = {
        "figure": "Fig2_model1_axis",
        "purpose": "Primary Model 1 effect-profile geometry and disease-axis position",
        "source_data": [str(DE_SOURCE), str(BOOTSTRAP_SOURCE)],
        "analysis_boundary": "Visualization of existing exported estimates; no model refitting",
        "gene_selection": "Genes with BH FDR < 0.05 in established COPD versus no COPD under Model 1",
        "n_axis_genes": int(len(profiles)),
        "density_encoding": {
            "type": "hexbin",
            "gridsize": GRID_SIZE,
            "common_extent": [-float(limit), float(limit), -float(limit), float(limit)],
            "normalization": "logarithmic gene count",
            "min_count": 1,
        },
        "line_encoding": {
            "gray_dashed": "identity line y=x",
            "colored_solid": "corrected λ over the observed established-COPD effect range",
        },
        "uncertainty": "95% paired gene-bootstrap confidence intervals from bootstrap_lambda.csv",
        "software": {
            "matplotlib": mpl.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "outputs": [str(OUTPUT_STEM.with_suffix(f".{ext}")) for ext in ("pdf", "png", "svg")],
    }
    PROVENANCE_FILE.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    return provenance


if __name__ == "__main__":
    print(json.dumps(build_figure(), indent=2))

