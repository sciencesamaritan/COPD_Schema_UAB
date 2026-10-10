"""Figure 2: disease-axis geometry (a, b) and genome-wide volcano plots (c-e) in one 183 mm figure.

Combines the panels of 2_Disease_Axis/Fig2_Disease_Axis_Geometry.py and Fig3_Volcano_Plots.py without refitting
anything. Point clouds and hexbins are rasterized at 450 dpi inside a vector PDF; text and lines stay vector.
"""
import importlib.util, sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nature_style as ns
import Fig3_Volcano_Plots as vol

_spec = importlib.util.spec_from_file_location("axisfig", HERE.parent / "2_Disease_Axis" / "Fig2_Disease_Axis_Geometry.py")
axisfig = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(axisfig)

#### Edit paths to match your file system ####
DE_SOURCE = axisfig.DE_SOURCE
BOOTSTRAP_SOURCE = axisfig.BOOTSTRAP_SOURCE
OUT_STEM = HERE.parent / "figures" / "Fig2_axis_volcano"


def build(de_source=DE_SOURCE, bootstrap_source=BOOTSTRAP_SOURCE, out_stem=OUT_STEM):
    ns.apply()
    axisfig.DE_SOURCE = Path(de_source)
    profiles = axisfig.load_effect_profiles()
    boot = pd.read_csv(bootstrap_source).set_index("group")
    vals = profiles[[axisfig.ESTABLISHED, axisfig.NEW, axisfig.EXCLUDED]].to_numpy().ravel()
    limit = np.ceil(np.nanmax(np.abs(vals)) * 10) / 10
    de = vol.load_de(de_source)

    fig = plt.figure(figsize=(ns.DOUBLE, 150 * ns.MM))
    top = fig.add_gridspec(1, 2, left=0.09, right=0.84, top=0.94, bottom=0.535, wspace=0.30)
    bot = fig.add_gridspec(1, 3, left=0.065, right=0.99, top=0.405, bottom=0.10, wspace=0.08)

    ax_a, ax_b = (fig.add_subplot(top[0, i]) for i in range(2))
    hbs = []
    for ax, key, grp in ((ax_a, axisfig.NEW, "Newly diagnosed"), (ax_b, axisfig.EXCLUDED, "Potentially excluded")):
        r = boot.loc[grp]
        hbs.append(axisfig.draw_density_panel(ax, profiles, key, grp, float(r["lambda"]), float(r["ci_lo"]),
                                              float(r["ci_hi"]), limit))
        ax.title.set_fontsize(6.5)
        ax.tick_params(length=2, width=0.5)
        ax.set_xticks(np.arange(-1, 1.01, 0.5)); ax.set_yticks(np.arange(-1, 1.01, 0.5))
    vmax = max(float(np.max(h.get_array())) for h in hbs)
    for h in hbs:
        h.set_clim(1, vmax)
    fig.canvas.draw()
    pb = ax_b.get_position()
    cax = fig.add_axes([pb.x1 + 0.015, pb.y0, 0.011, pb.height])
    cb = fig.colorbar(hbs[1], cax=cax)
    cb.set_label("Genes per hexagon (log scale)", fontsize=6, labelpad=3)
    cb.ax.tick_params(labelsize=5.5, length=2)
    cb.ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    cb.ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax_a.legend(loc="lower right", frameon=False, handlelength=2.2, borderaxespad=0.3, fontsize=5.5)
    for ax, letter in ((ax_a, "a"), (ax_b, "b")):
        ns.panel_label(ax, letter, x=-0.20, y=1.045)

    vaxs = [fig.add_subplot(bot[0, 0])]
    vaxs += [fig.add_subplot(bot[0, i], sharey=vaxs[0]) for i in (1, 2)]
    manifest = []
    for j, (k, title) in enumerate(vol.CONTRASTS.items()):
        manifest.append(vol.draw_volcano(vaxs[j], de[de.Contrast == k], title))
        if j:
            plt.setp(vaxs[j].get_yticklabels(), visible=False)
        ns.panel_label(vaxs[j], "cde"[j], x=-0.04 if j else -0.16, y=1.13)
    vaxs[1].set_xlabel("log₂ fold change"); vaxs[0].set_ylabel(r"$\mathsf{−log₁₀}\ P$")
    fig.legend(handles=vol.legend_handles(), frameon=False, loc="lower center", ncol=4, bbox_to_anchor=(0.53, 0.0),
               handletextpad=0.2, columnspacing=1.2)

    out_stem = Path(out_stem); out_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out_stem}.pdf", dpi=450)
    fig.savefig(f"{out_stem}.png", dpi=450)
    ns.export_panels(fig, str(out_stem))
    plt.close(fig)
    return pd.concat(manifest)


if __name__ == "__main__":
    build()
