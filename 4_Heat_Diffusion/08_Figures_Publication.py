"""Publication versions of the propagation figures at Nature Communications specifications.

Main: Fig. 4 (pathway landscape), Fig. 5 (calibration and agreement), Fig. 6 (gene modules).
Supplementary: Supplementary Figs. 6-9 (persistence, gene landscape, gene-level calibration, gene-pathway bridge).
Reads the outputs of scripts 01-05 (same inputs as 06 and 07); only the drawing differs.
"""
import sys, os, json, re
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "5_Tables_Figures"))
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter
from scipy.stats import spearmanr
import heat_kernel as core, calibrated_null as core2
from rds_reader import read_rds, names
import nature_style as ns
from paths import MSIGDB_CACHE, OUT, FIG, CLS, MODEL

ns.apply()
CM = LinearSegmentedColormap.from_list("d", ["#2a78d6", "#f0efec", "#e34948"]); INK = "#222222"; GR = "#9A9890"
BLUE, MAG = "#1f5fa8", "#B5446E"
PANELS = ["New", "Excluded", "Existing", "New vs Existing"]
G = ns.GROUP

# ============================ pathway-level state ============================
S = np.load(f"{OUT}/state_{CLS}_{MODEL}.npz", allow_pickle=True); V = np.load(f"{OUT}/state_v2_{CLS}_{MODEL}.npz")
M1 = json.load(open(f"{OUT}/diffusion_meta_{CLS}_{MODEL}.json")); M2 = json.load(open(f"{OUT}/diffusion_meta_v2_{CLS}_{MODEL}.json"))
pw = list(S["pw"]); P = S["P"]; lam, U = S["lam"], S["U"]; prog_of = V["prog_of"]; PROG = M2["programs"]
C = V["C"].astype(float); F = core2.null_factor(C)
SUMM = pd.read_csv(f"{OUT}/diffusion_summary_v2_{CLS}_{MODEL}.csv").set_index(["method", "panel"])
gx, gy = np.mgrid[0:1:220j, 0:1:220j]
d2 = (gx.ravel()[:, None] - P[:, 0]) ** 2 + (gy.ravel()[:, None] - P[:, 1]) ** 2; K = np.exp(-d2 / (2 * 0.022 ** 2)); Ks = K.sum(1); dens = Ks.reshape(gx.shape)
def pfield(v): return (K @ v / (Ks + 1e-12)).reshape(gx.shape)
H1 = core.heat(lam, U, 1.0); s1 = np.sqrt(((H1 @ F) ** 2).sum(1))
ABBR = {"Growth-factor RTK / PI3K": "RTK/PI3K", "MAPK / Ras / ERK": "MAPK/Ras/ERK", "Proteasome / ubiquitin": "Proteasome", "Innate / neutrophil": "Innate",
        "Cytokine / NF-κB": "Cytokine/NF-κB", "Chromatin / methylation": "Chromatin", "DNA repair / TP53": "DNA repair", "Translation / ribosome": "Translation",
        "GPCR / Ca2+ / NFAT": "GPCR/calcium", "Trafficking / vesicles": "Trafficking", "Mitochondria / OXPHOS": "OXPHOS", "Cell cycle / mitosis": "Cell cycle", "RNA processing": "RNA processing"}
def short_lab(l): return " + ".join(ABBR.get(x, x) for x in l.split(" + "))
top_prog = [p["program"] for p in sorted(PROG, key=lambda p: -p["size"])[:7]]


def p_landscape(ax, m, g, calibrated=True):
    z = V[f"z|{m}|{g}"]; Tf = np.clip((H1 @ z) / s1, -6, 6); Fv = pfield(Tf); fade = np.clip(dens / np.median(dens), 0, 1)
    rgba = CM((np.clip(Fv, -6, 6) + 6) / 12); rgba[..., 3] = 0.15 + 0.8 * fade
    ax.imshow(np.transpose(rgba, (1, 0, 2)), origin="lower", extent=[0, 1, 0, 1], interpolation="bilinear", rasterized=True)
    call = V[f"call|{m}|{g}"].astype(bool); sign = V[f"sign|{m}|{g}"]
    for sgn, col in [(1, "#7a1010"), (-1, "#0b3a75")]:
        ind = pfield((call & (sign == sgn)).astype(float)) * (fade > 0.45)
        if ind.max() > 0.42: ax.contour(gx, gy, ind, levels=[0.42], colors=col if calibrated else GR, linewidths=0.6, linestyles="-" if calibrated else "--")
    fdr = S[f"fdr|{m}|{g}"]; sz = 0.4 + 2.2 * np.minimum(-np.log10(np.clip(fdr, 1e-10, 1)), 10) / 10
    o = np.argsort(np.abs(z))
    ax.scatter(P[o, 0], P[o, 1], s=sz[o], c=np.clip(z[o], -6, 6), cmap=CM, vmin=-6, vmax=6, ec=np.where(fdr[o] < .05, INK, "none"), lw=0.2, zorder=3, rasterized=True)
    s_ = SUMM.loc[(m, g)]
    note = (f"Programs {int(s_.calls_overlap_fdr)} · raw {int(s_.raw_fdr05)} · inflation {s_.inflation_lambda:.1f}" if calibrated
            else f"Raw {int(s_.raw_fdr05)} · inflation {s_.inflation_lambda:.1f}")
    ax.text(0.0, -0.02, note, transform=ax.transAxes, fontsize=5, color=INK, va="top")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([]); ax.set_yticks([]); [s.set_visible(False) for s in ax.spines.values()]


def p_labels(ax):
    placed = []
    for k in top_prog:
        idx = np.where(prog_of == k)[0]
        if len(idx) < 5: continue
        c = np.median(P[idx], 0); d = c - 0.5; d = d / (np.linalg.norm(d) + 1e-9); tgt = np.clip(c + d * 0.2, [0.22, 0.05], [0.78, 0.96])
        for _ in range(25):
            if all(np.hypot(*(tgt - q)) > 0.1 for q in placed): break
            tgt = tgt + np.array([0, 0.04])
        lab = short_lab(PROG[k]["label"]); hw = 0.0105 * len(lab) + 0.01; tgt[0] = np.clip(tgt[0], hw, 1 - hw)
        placed.append(tgt)
        ax.annotate(lab, c, tgt, fontsize=5, ha="center", va="center", color=INK, zorder=6,
                    arrowprops=dict(arrowstyle="-", lw=0.35, color=INK), bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.85))


# ---------------- Fig. 4: pathway landscape ----------------
fig = plt.figure(figsize=(ns.DOUBLE, 112 * ns.MM))
gs = fig.add_gridspec(2, 5, width_ratios=[1, 1, 1, 1, 0.045], left=0.075, right=0.925, top=0.93, bottom=0.06, wspace=0.06, hspace=0.2)
for r, (m, cal) in enumerate([("singscore", True), ("FGSEA", False)]):
    for j, g in enumerate(PANELS):
        ax = fig.add_subplot(gs[r, j]); p_landscape(ax, m, g, cal)
        if r == 0: ax.set_title(G[g], fontsize=6.5, loc="center", pad=3)
        if cal and g in ("New", "Existing"): p_labels(ax)
        if j == 0:
            ax.text(-0.09, 0.5, "Participant-level singscore\n(calibrated)" if cal else "Preranked FGSEA\n(descriptive)", transform=ax.transAxes,
                    rotation=90, ha="center", va="center", fontsize=6)
            ns.panel_label(ax, "ab"[r], x=-0.13, y=1.02)
cax = fig.add_subplot(gs[:, 4]); sm = plt.cm.ScalarMappable(cmap=CM, norm=plt.Normalize(-6, 6))
cb = fig.colorbar(sm, cax=cax); cb.set_label("Standardized score", fontsize=5.5); cb.ax.tick_params(labelsize=5, width=0.4, length=1.5); cb.outline.set_linewidth(0.4)
cax.set_position([0.94, 0.3, 0.012, 0.4])
ns.save(fig, f"{FIG}/Fig4_pathway_landscape"); plt.close(fig)


# ---------------- Fig. 5: calibration and agreement ----------------
cal = M2["calibration_overlap"]; mor = M1["moran"]
GS_ = {"New": "newly diagnosed", "Excluded": "potentially excluded", "Existing": "established", "New vs Existing": "newly diagnosed vs established"}
fig = plt.figure(figsize=(ns.DOUBLE, 120 * ns.MM))
top = fig.add_gridspec(1, 2, width_ratios=[1, 1.9], wspace=0.06, left=0.215, right=0.985, top=0.95, bottom=0.6)
gsp = fig.add_gridspec(1, 3, wspace=0.42, left=0.075, right=0.985, top=0.42, bottom=0.12)
keys = [(m, g) for m in ("singscore", "FGSEA") for g in PANELS]
ax = fig.add_subplot(top[0, 0])
ax.axvspan(cal["moran_I_overlap_null_mean"] - 0.03, cal["moran_I_overlap_null_q95"], color="#E6E4DD", lw=0); ax.axvline(cal["moran_I_overlap_null_mean"], color=GR, lw=0.6)
for i, (m, g) in enumerate(keys): ax.plot(mor[f"{m}|{g}"][0], i, "o", color=BLUE if m == "singscore" else "#999", ms=3)
ax.set_yticks(range(len(keys))); ax.set_yticklabels([f"{m}, {GS_[g]}" for m, g in keys]); ax.set_ylim(len(keys) - 0.5, -1.4)
ax.set_xlabel("Moran's $I$ of raw pathway scores")
ax.text(cal["moran_I_overlap_null_mean"] + 0.005, -0.95, "Overlap noise", ha="center", va="center", fontsize=5, color="#555")
fig.text(0.01, 0.955, "a", fontsize=8, fontweight="bold", va="bottom")
ax = fig.add_subplot(top[0, 1], sharey=ax)
for i, (m, g) in enumerate(keys):
    r = SUMM.loc[(m, g)]
    ax.plot([r.calls_overlap_fdr, r.calls_label_null_fdr], [i, i], color="#DDD", lw=1.6, zorder=1, solid_capstyle="butt")
    ax.plot(r.calls_label_null_fdr, i, "o", color=MAG, ms=3, zorder=2); ax.plot(r.calls_overlap_fdr, i, "o", color=BLUE, ms=3, zorder=3)
    ax.plot(r.raw_fdr05, i, "|", color=INK, ms=5, mew=0.6, zorder=4)
ax.tick_params(axis="y", labelleft=False); ax.set_xlabel("Pathways called (FDR < 0.05)")
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=MAG, ms=3, label=f"Label-permutation null ({cal['label_null_fwer_rate']:.0%} FWER on overlap noise)"),
                   Line2D([], [], marker="o", ls="", color=BLUE, ms=3, label=f"Overlap-preserving null ({cal['overlap_null_fwer_rate']:.0%} FWER)"),
                   Line2D([], [], marker="|", ls="", color=INK, ms=5, mew=0.6, label="Raw FDR < 0.05")],
          frameon=False, loc="upper center", bbox_to_anchor=(0.4, -0.17), ncol=3, handletextpad=0.3, columnspacing=0.8)
ns.panel_label(ax, "b", x=-0.03, y=1.0)
SW = pd.read_csv(f"{OUT}/diffusion_tsweep_{CLS}_{MODEL}.csv")
for j, (m, a, b) in enumerate([("singscore", "New", "Existing"), ("singscore", "Excluded", "Existing"), ("FGSEA", "Excluded", "Existing")]):
    ax = fig.add_subplot(gsp[0, j]); d = SW[(SW.method == m) & (SW.a == a) & (SW.b == b)]
    ax.axvspan(M1["cap"], 4, color="#EEE", lw=0); ax.fill_between(d.t, d.msr_lo, d.msr_hi, color="#C3DBF2", lw=0, label="Spectral null, 95%")
    ax.plot(d.t, d.lab_lo, ":", color=GR, lw=0.7, label="Label null, 95%"); ax.plot(d.t, d.lab_hi, ":", color=GR, lw=0.7)
    ax.plot(d.t, d.r, color=MAG, lw=1.1, label="Observed"); ax.axvline(1, color=INK, ls="--", lw=0.5); ax.axhline(0, color="#CCC", lw=0.4)
    r1 = d[d.t == 1.0].iloc[0]
    ax.text(0.04, 0.05, f"$r$ = {r1.r:.2f}, spectral $P$ = {r1.p_msr:.3f}", transform=ax.transAxes, fontsize=5)
    ax.set_ylim(-1, 1.02); ax.set_xlim(0, 4); ax.set_xlabel("Diffusion time $t$"); ax.set_ylabel("Correlation of diffused maps" if j == 0 else "")
    ax.set_title(f"{m}: {GS_[a]} vs {GS_[b]}", fontsize=6); ax.text(3.12, 0.88, "Beyond cap", fontsize=5, color="#777", ha="center")
    if j == 0: fig.legend(*ax.get_legend_handles_labels(), frameon=False, loc="lower center", bbox_to_anchor=(0.53, 0.0), ncol=3, handlelength=1.6, handletextpad=0.3, columnspacing=1.5)
    ns.panel_label(ax, "cde"[j], x=-0.22 if j == 0 else -0.14, y=1.03)
ns.save(fig, f"{FIG}/Fig5_calibration_agreement"); plt.close(fig)

# ============================ gene-level state ============================
GSV = np.load(f"{OUT}/gene_state_{CLS}_{MODEL}.npz"); META = json.load(open(f"{OUT}/gene_meta_{CLS}_{MODEL}.json"))
MOD = pd.read_csv(f"{OUT}/gene_modules_{CLS}_{MODEL}.csv"); SUM = pd.read_csv(f"{OUT}/gene_diffusion_summary_{CLS}_{MODEL}.csv")
POW = json.load(open(f"{OUT}/gene_power_coexpr.json"))
Y = GSV["Y"]; genes = list(GSV["genes"]); n = len(genes); mod_of = GSV["mod_of"]; PC = GSV["P_coexpr"]

ACR = {"dna": "DNA", "rna": "RNA", "rrna": "rRNA", "trna": "tRNA", "tcf": "TCF", "wnt": "WNT", "bmal1": "BMAL1", "clock": "CLOCK", "arntl": "ARNTL",
       "er": "ER", "golgi": "Golgi", "slc": "SLC", "rac1": "RAC1", "fcgr": "FCGR", "gtpase": "GTPase", "atp": "ATP", "rtk": "RTK", "ras": "RAS", "erk": "ERK",
       "gf": "GF", "escrt": "ESCRT", "nf": "NF", "kb": "κB"}
MANUAL = {"HALLMARK_INTERFERON_GAMMA_RESPONSE": "Interferon-γ response", "REACTOME_TOLL_LIKE_RECEPTOR_CASCADES": "Toll-like receptor cascades",
          "REACTOME_CLASS_A_1_RHODOPSIN_LIKE_RECEPTORS": "Class A/1 (rhodopsin-like) receptors",
          "REACTOME_PHOSPHORYLATED_BMAL1_CLOCK_ARNTL_CLOCK_ACTIVATES_EXPRESSION_OF_CORE_CLOCK_GENES": "BMAL1:CLOCK activates core clock genes",
          "REACTOME_RRNA_MODIFICATION_IN_THE_NUCLEUS_AND_CYTOSOL": "rRNA modification in nucleus and cytosol",
          "REACTOME_TCF_DEPENDENT_SIGNALING_IN_RESPONSE_TO_WNT": "TCF-dependent signaling in response to WNT",
          "REACTOME_ER_TO_GOLGI_ANTEROGRADE_TRANSPORT": "ER-to-Golgi anterograde transport",
          "REACTOME_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT": "SLC-mediated transmembrane transport",
          "REACTOME_RAC1_GTPASE_CYCLE": "RAC1 GTPase cycle", "REACTOME_FCGR_ACTIVATION": "FCGR activation"}
def set_label(s):
    if s in MANUAL: return MANUAL[s]
    w = re.sub(r"^(REACTOME|HALLMARK|KEGG|BIOCARTA)_", "", s).lower().split("_")
    w = [ACR.get(x, x) for x in w]; t = " ".join(w); return t[0].upper() + t[1:]
MOD["disp"] = MOD.best_set.map(set_label)
MOD["sel"] = (MOD["frac_called_New"] >= 0.25) | (MOD["frac_called_New vs Existing"] >= 0.25)

# ---------------- Fig. 6: module matrix ----------------
show = pd.concat([MOD[MOD.sel], MOD[~MOD.sel].sort_values("size", ascending=False).head(14)]).drop_duplicates("module").sort_values("mean_T_New", ascending=False)
nr = len(show)
fig = plt.figure(figsize=(ns.DOUBLE, (4.0 * nr + 30) * ns.MM))
ax = fig.add_axes([0.29, 30 / (4.0 * nr + 30) * 0.62, 0.34, 1 - 30 / (4.0 * nr + 30)])
for i, (_, r) in enumerate(show.iterrows()):
    for j, g in enumerate(PANELS):
        mt, fc = r[f"mean_T_{g}"], r[f"frac_called_{g}"]
        ax.scatter(j, i, s=2 + 52 * fc, c=[np.clip(mt, -6, 6)], cmap=CM, vmin=-6, vmax=6, ec=INK if fc >= 0.5 else "none", lw=0.4)
    ax.text(-0.55, i, r.disp, ha="right", va="center", fontsize=5.5, weight="bold" if r.sel else "normal")
    ax.text(4.0, i, f"{r['size']}", ha="right", va="center", fontsize=5.5, color="#444", transform=ax.transData)
    ax.text(4.15, i, ", ".join(eval(r.hubs)[:3]), ha="left", va="center", fontsize=5.5, color="#444", style="italic")
ax.set_xticks(range(4)); ax.set_xticklabels([ns.GROUP_2L[g] for g in PANELS], fontsize=5.5); ax.xaxis.tick_top(); ax.tick_params(axis="x", length=0, pad=2)
ax.set_yticks([]); ax.set_ylim(nr - 0.5, -0.7); ax.set_xlim(-0.45, 3.4); [s.set_visible(False) for s in ax.spines.values()]
ax.text(4.0, -0.95, "Genes", fontsize=5.5, ha="right", va="bottom", weight="bold", color="#444"); ax.text(4.15, -0.95, "Hub genes", fontsize=5.5, ha="left", va="bottom", weight="bold", color="#444")
lax = fig.add_axes([0.33, 0.0, 0.36, 22 / (4.0 * nr + 30)]); lax.axis("off"); lax.set_xlim(0, 1); lax.set_ylim(0, 1)
for k, f in enumerate([0.1, 0.5, 1.0]):
    lax.scatter(0.05 + k * 0.3, 0.55, s=2 + 52 * f, c="#bbb", ec=INK if f >= 0.5 else "none", lw=0.4, clip_on=False); lax.text(0.1 + k * 0.3, 0.55, f"{f:.0%} called", va="center", fontsize=5.5)
cax = fig.add_axes([0.76, 0.55 * 22 / (4.0 * nr + 30) - 0.004, 0.18, 0.008]); sm = plt.cm.ScalarMappable(cmap=CM, norm=plt.Normalize(-6, 6))
cb = fig.colorbar(sm, cax=cax, orientation="horizontal"); cb.set_label("Mean diffused score", fontsize=5.5, labelpad=1); cb.ax.tick_params(labelsize=5, width=0.4, length=1.5, pad=1); cb.outline.set_linewidth(0.4)
ns.save(fig, f"{FIG}/Fig6_gene_modules"); plt.close(fig)

# ---------------- Supplementary Fig. 7: gene landscape ----------------
NG = 300; SIG = 4.0
ix = np.clip((Y * NG).astype(int), 0, NG - 1)
def gsmooth(v):
    a = np.zeros((NG, NG)); np.add.at(a, (ix[:, 0], ix[:, 1]), v); return gaussian_filter(a, SIG)
DEN = gsmooth(np.ones(n)); gfade = np.clip(DEN / np.quantile(DEN[DEN > 1e-3], 0.5), 0, 1)
def gfield(v): return gsmooth(v) / (DEN + 1e-9)
LABELS = {"New": list(MOD[MOD["frac_called_New"] >= 0.25].module), "New vs Existing": list(MOD[MOD["frac_called_New vs Existing"] >= 0.25].module),
          "Existing": list(MOD.sort_values("size", ascending=False).module[:4]), "Excluded": []}
def short(s): return s if len(s) <= 26 else s[:24] + "…"
def g_landscape(ax, g):
    T = np.clip(GSV[f"T1|{g}"], -6, 6); Fv = np.clip(gfield(T), -6, 6)
    rgba = CM((Fv + 6) / 12); rgba[..., 3] = 0.1 + 0.6 * gfade
    ax.imshow(np.transpose(rgba, (1, 0, 2)), origin="lower", extent=[0, 1, 0, 1], interpolation="bilinear", rasterized=True)
    call = GSV[f"call|{g}"].astype(bool); raw = GSV[f"fdr|{g}"] < .05; xs = (np.arange(NG) + 0.5) / NG
    ax.scatter(Y[~call, 0], Y[~call, 1], s=0.08, color="#bbb", lw=0, alpha=0.5, zorder=2, rasterized=True)
    o = np.where(call)[0]; o = o[np.argsort(np.abs(T[o]))]
    ax.scatter(Y[o, 0], Y[o, 1], s=0.6, c=T[o], cmap=CM, vmin=-6, vmax=6, lw=0, zorder=3, rasterized=True)
    ax.scatter(Y[raw, 0], Y[raw, 1], s=1.2, facecolors="none", edgecolors=INK, lw=0.15, zorder=4, rasterized=True)
    items = []
    for m in LABELS[g]:
        idx = np.where(mod_of == m)[0]
        if len(idx) < 8: continue
        ind = gfield((mod_of == m).astype(float)) * (gfade > 0.2)
        if ind.max() > 0.5: ax.contour(xs, xs, ind.T, levels=[0.5], colors="#333", linewidths=0.4, linestyles="--", zorder=5)
        items.append((m, Y[idx][np.argmax(ind[ix[idx, 0], ix[idx, 1]])]))
    for side in (0, 1):
        its = sorted([it for it in items if (it[1][0] >= 0.5) == bool(side)], key=lambda it: -it[1][1]); ys = []
        for _, c in its:
            y = min(c[1], ys[-1] - 0.11) if ys else min(c[1], 0.93); ys.append(y)
        if ys and ys[-1] < 0.08: ys = list(np.array(ys) + (0.08 - ys[-1]))
        for (m, c), y in zip(its, ys):
            row = MOD[MOD.module == m].iloc[0]; x = 0.99 if side else 0.01
            ax.annotate(short(row.disp), c, (x, y), fontsize=5, ha="right" if side else "left", va="center", color=INK, zorder=6,
                        arrowprops=dict(arrowstyle="-", lw=0.3, color=INK), bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="#ccc", lw=0.25, alpha=0.92))
    r = SUM[(SUM.null == "coexpr") & (SUM.panel == g)].iloc[0]; ri = SUM[(SUM.null == "iid") & (SUM.panel == g)].iloc[0]
    ax.text(0.0, -0.02, f"Calls {r.calls_fdr} (rescued {r.rescued}) · raw {r.raw_fdr05}\ninflation {r.lambda_gc:.1f} · independent-gene null {ri.calls_fdr}",
            transform=ax.transAxes, fontsize=5, va="top")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([]); ax.set_yticks([]); [s_.set_visible(False) for s_ in ax.spines.values()]
fig = plt.figure(figsize=(ns.DOUBLE, 60 * ns.MM))
gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 1, 1, 0.045], left=0.01, right=0.955, top=0.9, bottom=0.17, wspace=0.05)
for j, g in enumerate(PANELS):
    ax = fig.add_subplot(gs[0, j]); g_landscape(ax, g); ax.set_title(G[g], loc="center", fontsize=6.5); ns.panel_label(ax, "abcd"[j], x=0.04, y=1.0)
cax = fig.add_subplot(gs[0, 4]); cb = fig.colorbar(plt.cm.ScalarMappable(cmap=CM, norm=plt.Normalize(-6, 6)), cax=cax)
cb.set_label("Standardized score", fontsize=5.5); cb.ax.tick_params(labelsize=5, width=0.4, length=1.5); cb.outline.set_linewidth(0.4); cax.set_position([0.965, 0.3, 0.012, 0.45])
ns.save(fig, f"{FIG}/SFig7_gene_landscape"); plt.close(fig)

# ---------------- Supplementary Fig. 8: gene-level calibration ----------------
fig = plt.figure(figsize=(ns.DOUBLE, 62 * ns.MM)); gs = fig.add_gridspec(1, 3, wspace=0.75, left=0.07, right=0.985, top=0.9, bottom=0.3, width_ratios=[1, 1, 1.1])
ax = fig.add_subplot(gs[0]); k = np.arange(1, n + 1); zh2 = (GSV["z|Excluded"] @ np.load(f"{OUT}/gene_graph_{CLS}_{MODEL}.npz")["U"]) ** 2
ax.scatter(k, zh2, s=0.15, color="#ccc", lw=0, rasterized=True); ax.plot(k, PC, color=BLUE, lw=0.9, label="Co-expression null")
ax.axhline(1, color=MAG, ls="--", lw=0.7, label="Independent-gene null"); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-3, 60)
from matplotlib.ticker import FuncFormatter as _FF; _plain=_FF(lambda v,_: f"{v:g}"); ax.xaxis.set_major_formatter(_plain); ax.yaxis.set_major_formatter(_plain)
ax.set_xlabel("Graph eigenvector (smooth to rough)"); ax.set_ylabel("Power, potentially excluded $t$ map"); ax.legend(frameon=False, loc="upper right", handlelength=1.4)
ns.panel_label(ax, "a", x=-0.3, y=1.02)
ax = fig.add_subplot(gs[1])
for i, g in enumerate(PANELS):
    a = SUM[(SUM.null == "iid") & (SUM.panel == g)].iloc[0]; b = SUM[(SUM.null == "coexpr") & (SUM.panel == g)].iloc[0]
    ax.plot([b.calls_fdr, a.calls_fdr], [i, i], color="#DDD", lw=1.6, solid_capstyle="butt"); ax.plot(a.calls_fdr, i, "o", color=MAG, ms=3); ax.plot(b.calls_fdr, i, "o", color=BLUE, ms=3)
    ax.plot(a.raw_fdr05, i, "|", color=INK, ms=5, mew=0.6)
ax.set_yticks(range(4)); ax.set_yticklabels([ns.GROUP_2L[g] for g in PANELS]); ax.invert_yaxis(); ax.set_xscale("symlog", linthresh=10); ax.set_xlim(-1, 1.2e4)
ax.xaxis.set_major_formatter(_FF(lambda v,_: f"{v:,.0f}"))
ax.set_xlabel("Genes called (FDR < 0.05)")
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=MAG, ms=3, label="Independent-gene null"), Line2D([], [], marker="o", ls="", color=BLUE, ms=3, label="Co-expression null"),
                   Line2D([], [], marker="|", ls="", color=INK, ms=5, mew=0.6, label="Raw FDR < 0.05")], frameon=False, loc="upper center", bbox_to_anchor=(0.4, -0.25), ncol=1)
ns.panel_label(ax, "b", x=-0.62, y=1.02)
ax = fig.add_subplot(gs[2]); keys = list(POW)
for i, kk in enumerate(keys):
    ax.barh(i - 0.18, POW[kk]["raw_power"], 0.34, color="#bbb"); ax.barh(i + 0.18, POW[kk]["power"], 0.34, color=BLUE)
    ax.text(POW[kk]["power"] + 0.015, i + 0.18, f"{POW[kk]['power']:.0%}", va="center", fontsize=5)
def pl(kk):
    m_ = re.match(r"(\d+)g_shift([\d.]+)", kk); return f"{m_.group(1)} genes, {m_.group(2)} s.d." if m_ else kk
ax.set_yticks(range(len(keys))); ax.set_yticklabels([pl(kk) for kk in keys]); ax.invert_yaxis(); ax.set_xlim(0, 0.8)
ax.set_xlabel("Planted genes recovered")
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color="#bbb", label="Gene-wise FDR"), plt.Rectangle((0, 0), 1, 1, color=BLUE, label="Propagation, co-expression null")],
          frameon=False, loc="upper center", bbox_to_anchor=(0.4, -0.25), ncol=1)
ns.panel_label(ax, "c", x=-0.5, y=1.02)
ns.save(fig, f"{FIG}/SFig8_gene_calibration"); plt.close(fig)

# ---------------- Supplementary Fig. 9: gene-pathway bridge ----------------
cache = read_rds(MSIGDB_CACHE); so = cache.value[names(cache).index("sets")]
GSETS = {k: v.value for k, v in zip(names(so), so.value)}; gi = {g: i for i, g in enumerate(genes)}
memb = [np.array([gi[g] for g in set(GSETS[p]) if g in gi]) for p in pw]
fig, axs = plt.subplots(1, 4, figsize=(ns.DOUBLE, 52 * ns.MM)); fig.subplots_adjust(left=0.065, right=0.99, top=0.88, bottom=0.2, wspace=0.3)
for j, g in enumerate(PANELS):
    gm = np.array([GSV[f"T1|{g}"][m].mean() if len(m) else np.nan for m in memb]); pz = V[f"z|singscore|{g}"]; ok = np.isfinite(gm)
    rho = spearmanr(gm[ok], pz[ok])[0]; pc = V[f"call|singscore|{g}"].astype(bool)
    ax = axs[j]; ax.scatter(gm[~pc], pz[~pc], s=0.6, color="#ccc", lw=0, rasterized=True)
    ax.scatter(gm[pc], pz[pc], s=1, c=np.where(V[f"sign|singscore|{g}"][pc] > 0, "#e34948", "#2a78d6"), lw=0, rasterized=True)
    ax.axhline(0, color="#ddd", lw=0.4); ax.axvline(0, color="#ddd", lw=0.4)
    ax.set_title(G[g], fontsize=6); ax.set_xlabel("Mean diffused score of member genes"); ax.set_ylabel("Pathway singscore $z$" if j == 0 else "")
    ax.text(0.04, 0.95, f"$\\rho$ = {rho:.2f}", transform=ax.transAxes, va="top", fontsize=5.5)
    ns.panel_label(ax, "abcd"[j], x=-0.2 if j == 0 else -0.12, y=1.04)
ns.save(fig, f"{FIG}/SFig9_gene_pathway_bridge"); plt.close(fig)
# ---------------- Supplementary Fig. 6: persistence ----------------
TG = S["TG"]
cols = [("singscore", g) for g in PANELS] + [("FGSEA", g) for g in PANELS]
rows = []
for k in top_prog[:6]:
    idx = np.where(prog_of == k)[0]; sc = np.zeros(len(idx))
    for (m, g) in cols[:4]: sc += V[f"call|{m}|{g}"][idx]
    best = idx[np.argsort(-(sc + 0.01 * np.abs(V["z|singscore|Existing"][idx])))[:5]]
    rows.append(("hdr", k)); rows += [("pw", i) for i in best]
ACR2 = {"tcr": "TCR", "plcg": "PLCG", "itpr": "ITPR", "pips": "PIPs", "edg1": "EDG1", "par1": "PAR1", "fcer1": "FCER1", "ets": "ETS", "nfat": "NFAT", "tnf": "TNF", "tp53": "TP53",
        "cd28": "CD28", "pi3k": "PI3K", "akt": "AKT", "hdr": "HDR", "hrr": "HRR", "mapk": "MAPK", "kit": "KIT", "g1": "G1", "s": "S", "clrs": "CLRs", "o": "O", "n": "N", "map": "MAP"}
def pw_label(s):
    s = re.sub(r"^(KEGG_MEDICUS_REFERENCE|KEGG_MEDICUS|REACTOME|BIOCARTA|KEGG)_", "", s)
    w = [ACR2.get(x, ACR.get(x, x)) for x in s.lower().split("_")]; s = " ".join(w); s = s[0].upper() + s[1:]
    for a, b in [("G1 S ", "G1/S "), ("C type lectin receptors CLRs", "C-type lectin receptors"), ("N glycan", "N-glycan"), ("O linked", "O-linked"),
                 ("PIPs at the early", "PIPs at the early"), ("Hypusinylation hydroxylation", "hypusinylation, hydroxylation"), ("HDR through homologous recombination HRR", "HDR through homologous recombination"),
                 ("Gamma carboxylation", "γ-Carboxylation"), ("Interleukin 1 family", "Interleukin-1 family"), ("Interleukin 1 signaling", "Interleukin-1 signaling"), ("Interleukin 10", "Interleukin-10"),
                 ("CD28 dependent", "CD28-dependent"), ("TCR PLCG ITPR signaling pathway", "TCR, PLCG and ITPR signaling")]:
        s = s.replace(a, b)
    return s if len(s) <= 52 else s[:50] + "…"
hgt = (3.0 * len(rows) + 22) * ns.MM
fig = plt.figure(figsize=(ns.DOUBLE, hgt)); ax = fig.add_axes([0.36, 4 * ns.MM / hgt, 0.62, 1 - 18 * ns.MM / hgt])
y = 0
for kind, v in rows:
    if kind == "hdr":
        y += 0.4; ax.text(-0.01, y, PROG[v]["label"].replace("Ca2+", "calcium"), transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6, weight="bold"); y += 1; continue
    ax.text(-0.01, y, pw_label(pw[v]), transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=5, color="#333")
    for c, (m, g) in enumerate(cols):
        per = V[f"persist|{m}|{g}"][v]; call = V[f"call|{m}|{g}"][v]; sgn = V[f"sign|{m}|{g}"][v]; col = "#c0392b" if sgn > 0 else BLUE
        if call and per.any():
            ts = TG[per]; x0 = c + 0.06 + ts.min() / TG.max() * 0.84; x1 = c + 0.06 + ts.max() / TG.max() * 0.84 + 0.02
            ax.add_patch(plt.Rectangle((x0, y - 0.32), x1 - x0, 0.64, color=col, alpha=0.9 if m == "singscore" else 0.45, lw=0))
            xb = c + 0.06 + V[f"best_t|{m}|{g}"][v] / TG.max() * 0.84; ax.plot([xb, xb], [y - 0.4, y + 0.4], color=INK, lw=0.5)
        if S[f"fdr|{m}|{g}"][v] < .05: ax.plot(c + 0.95, y, "o", ms=1.5, color=INK)
    y += 1
for c in range(len(cols) + 1): ax.axvline(c, color="#DDD", lw=0.4)
ax.axvline(4, color=INK, lw=0.6)
for c, (m, g) in enumerate(cols): ax.text(c + 0.5, -0.5, ns.GROUP_2L[g], ha="center", va="bottom", fontsize=5, linespacing=1.0)
ax.text(2, -2.6, "Participant-level singscore (calibrated)", ha="center", fontsize=6, weight="bold"); ax.text(6, -2.6, "Preranked FGSEA (descriptive)", ha="center", fontsize=6, weight="bold", color="#777")
ax.set_xlim(0, len(cols)); ax.set_ylim(y, -3.2); ax.set_xticks([]); ax.set_yticks([]); [s.set_visible(False) for s in ax.spines.values()]
ns.save(fig, f"{FIG}/SFig6_persistence"); plt.close(fig)
print("publication figures done")
