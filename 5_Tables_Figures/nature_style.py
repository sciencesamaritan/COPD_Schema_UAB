"""Shared figure style for Nature Communications submission.

Nature portfolio figure specifications: sans-serif type (Arial or Helvetica), 5-7 pt at final size,
panel labels 8 pt bold upright lowercase, vector PDF with embedded TrueType fonts, RGB.
Widths: 89 mm (single column), 120-136 mm (1.5 column), 183 mm (double column); maximum height 170 mm.

Arial is used when installed; Liberation Sans (metric-identical to Arial) is the fallback.
"""
import matplotlib as mpl

MM = 1 / 25.4
SINGLE, ONE_HALF, DOUBLE, MAX_H = 89 * MM, 120 * MM, 183 * MM, 170 * MM
SANS = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]

GROUP = {"New": "Newly diagnosed", "Excluded": "Potentially excluded", "Existing": "Established COPD",
         "New vs Existing": "Newly diagnosed vs established"}
GROUP_2L = {"New": "Newly\ndiagnosed", "Excluded": "Potentially\nexcluded", "Existing": "Established\nCOPD",
            "New vs Existing": "Newly diagnosed\nvs established"}

RC = {
    "font.family": "sans-serif", "font.sans-serif": SANS, "font.size": 6,
    "axes.titlesize": 6.5, "axes.titleweight": "normal", "axes.titlelocation": "left", "axes.titlepad": 3,
    "axes.labelsize": 6, "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5,
    "legend.title_fontsize": 6, "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4, "xtick.major.size": 2, "ytick.major.size": 2,
    "xtick.minor.size": 1.2, "ytick.minor.size": 1.2, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
    "axes.labelpad": 2, "lines.linewidth": 0.8, "patch.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "mathtext.fontset": "custom", "mathtext.rm": "sans", "mathtext.it": "sans:italic", "mathtext.bf": "sans:bold",
    "mathtext.sf": "sans", "mathtext.default": "it",
    "savefig.dpi": 450, "figure.dpi": 150, "image.interpolation": "bilinear",
}


def apply():
    mpl.rcParams.update(RC)


def panel_label(ax, letter, x=-0.02, y=1.0, fig=None):
    """8 pt bold upright lowercase panel label placed in axes-fraction coordinates."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8, fontweight="bold", ha="right", va="bottom")


def save(fig, stem):
    fig.savefig(f"{stem}.pdf")
    fig.savefig(f"{stem}.png", dpi=300)
