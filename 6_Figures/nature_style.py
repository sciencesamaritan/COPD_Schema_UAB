"""Shared figure style for Nature Communications submission.

Nature portfolio figure specifications: sans-serif type (Arial or Helvetica), 5-7 pt at final size,
panel labels 8 pt bold upright lowercase, vector PDF with embedded TrueType fonts, RGB.
Widths: 89 mm (single column), 120-136 mm (1.5 column), 183 mm (double column); maximum height 170 mm.

Arial matches the manuscript body text. Arial Unicode MS supplies glyphs Arial lacks (subscript digits);
Liberation Sans is a last resort on systems without Arial.
"""
import matplotlib as mpl

MM = 1 / 25.4
SINGLE, ONE_HALF, DOUBLE, MAX_H = 89 * MM, 120 * MM, 183 * MM, 170 * MM
SANS = ["Arial", "Arial Unicode MS", "Helvetica", "Liberation Sans", "DejaVu Sans"]

GROUP = {"New": "Newly diagnosed", "Excluded": "Potentially excluded", "Existing": "Established COPD",
         "New vs Existing": "Newly diagnosed vs established"}
GROUP_2L = {"New": "Newly\ndiagnosed", "Excluded": "Potentially\nexcluded", "Existing": "Established\nCOPD",
            "New vs Existing": "Newly diagnosed\nvs established"}

RC = {
    "font.family": SANS, "font.sans-serif": SANS, "font.size": 6,
    "axes.titlesize": 6.5, "axes.titleweight": "normal", "axes.titlelocation": "left", "axes.titlepad": 3,
    "axes.labelsize": 6, "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5,
    "legend.title_fontsize": 6, "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4, "xtick.major.size": 2, "ytick.major.size": 2,
    "xtick.minor.size": 1.2, "ytick.minor.size": 1.2, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
    "axes.labelpad": 2, "lines.linewidth": 0.8, "patch.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "mathtext.sf": "Arial Unicode MS",  # \mathsf{} carries Unicode sub/superscript digits in Arial design
    "mathtext.default": "it",
    "savefig.dpi": 450, "figure.dpi": 150, "image.interpolation": "bilinear",
}


def apply():
    mpl.rcParams.update(RC)


def panel_label(ax, letter, x=-0.02, y=1.0, fig=None):
    """8 pt bold upright lowercase panel label placed in axes-fraction coordinates."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8, fontweight="bold", ha="right", va="bottom")


def _panel_labels(fig):
    """Panel labels are 8 pt bold single lowercase letters, drawn on the figure or on any axes."""
    out = []
    for t in list(fig.texts) + [t for ax in fig.axes for t in ax.texts]:
        s = t.get_text().strip()
        if len(s) == 1 and s.islower() and t.get_fontsize() == 8 and t.get_fontweight() in ("bold", 700):
            out.append((s, t))
    return out


def export_panels(fig, stem, subdir="panels"):
    """Write each lettered panel as its own vector PDF: <figure dir>/<subdir>/<stem>_<letter>.pdf.

    Axes, figure legends and figure-level notes are assigned to the nearest panel label above and to their left
    (an axes that carries its own label keeps it). Everything outside the panel is hidden before saving, so no
    neighbouring content leaks into the crop. A colour bar that spans more than one panel is written separately as
    <stem>_colorbar.pdf. Rasterized layers keep savefig.dpi.
    """
    from pathlib import Path
    from matplotlib.transforms import Bbox
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    labels = _panel_labels(fig)
    if not labels:
        return []
    lab_bb = {k: t.get_window_extent(r) for k, t in labels}
    label_texts = {t for _, t in labels}

    def nearest(bb):
        cx, cy = bb.x0 + bb.width / 2, bb.y0 + bb.height / 2
        cand = [(k, lb) for k, lb in lab_bb.items() if lb.x0 <= cx and lb.y0 >= cy - 2] or list(lab_bb.items())
        return min(cand, key=lambda kb: (kb[1].x0 - bb.x0) ** 2 + (kb[1].y1 - bb.y1) ** 2)[0]

    items = {k: [] for k in lab_bb}
    cbars = []
    for ax in fig.axes:
        if not ax.get_visible():
            continue
        own = [k for k, t in labels if t in ax.texts]
        if hasattr(ax, "_colorbar") and not own:
            cbars.append(ax); continue
        items[own[0] if own else nearest(ax.get_window_extent(r))].append(ax)
    for art in list(fig.legends) + [t for t in fig.texts if t not in label_texts and t.get_visible()]:
        items[nearest(art.get_window_extent(r))].append(art)
    loose = []
    for cb in cbars:  # attach a colour bar to the panel whose vertical span contains it
        cbb = cb.get_window_extent(r); home = None
        for k, arts in items.items():
            axs = [a for a in arts if hasattr(a, "get_tightbbox") and a in fig.axes]
            if not axs:
                continue
            u = Bbox.union([a.get_window_extent(r) for a in axs])
            if cbb.y0 >= u.y0 - 0.1 * u.height and cbb.y1 <= u.y1 + 0.1 * u.height and abs(cbb.x0 - u.x1) < 0.25 * fig.bbox.width:
                home = k; break
        (items[home].append(cb) if home else loose.append(cb))

    every = list(fig.axes) + list(fig.legends) + list(fig.texts)
    state = {a: a.get_visible() for a in every}
    out_dir = Path(stem).parent / subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    groups = dict(items)
    if loose:
        groups["colorbar"] = loose
    for key, arts in sorted(groups.items()):
        keep = set(arts)
        if key in lab_bb:
            keep.add(dict(labels)[key])
        for a in every:
            a.set_visible(state[a] and (a in keep or (a in fig.axes and False)))
        for k, t in labels:  # labels drawn inside axes follow their axes; hide other panels' figure-level labels
            if t in fig.texts:
                t.set_visible(k == key)
        fig.canvas.draw()
        boxes = [a.get_tightbbox(r) if hasattr(a, "get_tightbbox") else a.get_window_extent(r) for a in keep]
        boxes = [b for b in boxes if b is not None and b.width > 0]
        bb = Bbox.union(boxes).padded(3).transformed(fig.dpi_scale_trans.inverted())
        f = out_dir / f"{Path(stem).name}_{key}.pdf"
        fig.savefig(f, bbox_inches=bb)
        written.append(f)
    for a, v in state.items():
        a.set_visible(v)
    fig.canvas.draw()
    return written


def save(fig, stem, panels=True):
    fig.savefig(f"{stem}.pdf")
    fig.savefig(f"{stem}.png", dpi=450)
    if panels:
        export_panels(fig, stem)
