"""Gene-level diffusion figures: G1 landscape, G2 module matrix, G3 null calibration, G4 gene-pathway bridge."""
import sys, os, json, re
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.figure as _mf
_mf.Figure.suptitle=lambda self,*a,**k: None
_mf.Figure.text=lambda self,*a,**k: None

from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter
from scipy.stats import spearmanr, norm
from rds_reader import read_rds, names

from paths import DE_FILE, SINGSCORE_FILE, FGSEA_FILE, MSIGDB_CACHE, OUT, FIG, CLS, MODEL
S = np.load(f"{OUT}/gene_state_{CLS}_{MODEL}.npz"); META = json.load(open(f"{OUT}/gene_meta_{CLS}_{MODEL}.json"))
MOD = pd.read_csv(f"{OUT}/gene_modules_{CLS}_{MODEL}.csv"); SUM = pd.read_csv(f"{OUT}/gene_diffusion_summary_{CLS}_{MODEL}.csv")
POW = json.load(open(f"{OUT}/gene_power_coexpr.json"))
Y = S["Y"]; genes = list(S["genes"]); n = len(genes); mod_of = S["mod_of"]; lam = S["lam"]; P = S["P_coexpr"]
PANELS = ["New", "Excluded", "Existing", "New vs Existing"]
TITLE = {"New": "New vs No COPD", "Excluded": "Excluded vs No COPD", "Existing": "Existing vs No COPD", "New vs Existing": "New vs Existing (direct)"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.titlesize": 8.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
                     "axes.spines.top": False, "axes.spines.right": False})
CM = LinearSegmentedColormap.from_list("d", ["#2a78d6", "#f0efec", "#e34948"]); INK = "#222"; GR = "#9A9890"
def short(s): s = str(s); return s if len(s) <= 30 else s[:28] + "…"

# ---------- smoothing on a grid ----------
NG = 300; SIG = 4.0
ix = np.clip((Y * NG).astype(int), 0, NG - 1)
def gsmooth(v):
    a = np.zeros((NG, NG)); np.add.at(a, (ix[:, 0], ix[:, 1]), v); return gaussian_filter(a, SIG)
DEN = gsmooth(np.ones(n)); fade = np.clip(DEN / np.quantile(DEN[DEN > 1e-3], 0.5), 0, 1)
def field(v): return gsmooth(v) / (DEN + 1e-9)

# modules to label: discriminating for New or the direct contrast, plus the largest Existing modules
MOD["sel"] = (MOD["frac_called_New"] >= 0.25) | (MOD["frac_called_New vs Existing"] >= 0.25)
lab_mods = list(MOD[MOD.sel].module) + [m for m in MOD.sort_values("size", ascending=False).module if m not in set(MOD[MOD.sel].module)][:3]

LABELS = {"New": list(MOD[MOD["frac_called_New"] >= 0.25].module), "New vs Existing": list(MOD[MOD["frac_called_New vs Existing"] >= 0.25].module),
          "Existing": list(MOD.sort_values("size", ascending=False).module[:4]), "Excluded": []}
def landscape(ax, g):
    T = np.clip(S[f"T1|{g}"], -6, 6); Fv = np.clip(field(T), -6, 6)
    rgba = CM((Fv + 6) / 12); rgba[..., 3] = 0.1 + 0.6 * fade
    ax.imshow(np.transpose(rgba, (1, 0, 2)), origin="lower", extent=[0, 1, 0, 1], interpolation="bilinear")
    call = S[f"call|{g}"].astype(bool); raw = S[f"fdr|{g}"] < .05; xs = (np.arange(NG) + 0.5) / NG
    ax.scatter(Y[~call, 0], Y[~call, 1], s=0.5, color="#bbb", lw=0, alpha=0.5, zorder=2)
    o = np.where(call)[0]; o = o[np.argsort(np.abs(T[o]))]
    ax.scatter(Y[o, 0], Y[o, 1], s=3.2, c=T[o], cmap=CM, vmin=-6, vmax=6, lw=0, zorder=3)
    ax.scatter(Y[raw, 0], Y[raw, 1], s=6, facecolors="none", edgecolors=INK, lw=0.3, zorder=4)
    items = []
    for m in LABELS[g]:
        idx = np.where(mod_of == m)[0]
        if len(idx) < 8: continue
        ind = field((mod_of == m).astype(float)) * (fade > 0.2)
        if ind.max() > 0.5: ax.contour(xs, xs, ind.T, levels=[0.5], colors="#333", linewidths=0.6, linestyles="--", zorder=5)
        items.append((m, Y[idx][np.argmax(ind[ix[idx, 0], ix[idx, 1]])]))
    for side in (0, 1):   # stacked callouts in the left and right margins, ordered by height
        its = sorted([it for it in items if (it[1][0] >= 0.5) == bool(side)], key=lambda it: -it[1][1])
        ys = []
        for _, c in its:
            y = min(c[1], ys[-1] - 0.085) if ys else min(c[1], 0.93); ys.append(y)
        if ys and ys[-1] < 0.07: ys = list(np.array(ys) + (0.07 - ys[-1]))
        for (m, c), y in zip(its, ys):
            row = MOD[MOD.module == m].iloc[0]; x = 0.985 if side else 0.015
            ax.annotate(f"{short(row.label)}\n{', '.join(eval(row.hubs)[:2])}", c, (x, y), fontsize=5.4, ha="right" if side else "left", va="center", color=INK, zorder=6,
                        arrowprops=dict(arrowstyle="-", lw=0.4, color=INK), bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="#ccc", lw=0.3, alpha=0.92))
    r = SUM[(SUM.null == "coexpr") & (SUM.panel == g)].iloc[0]; ri = SUM[(SUM.null == "iid") & (SUM.panel == g)].iloc[0]
    ax.text(0.01, 0.01, f"calls {r.calls_fdr} (rescued {r.rescued})   raw {r.raw_fdr05}   λ {r.lambda_gc:.1f}   [iid null: {ri.calls_fdr}]",
            transform=ax.transAxes, fontsize=6.2, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.2))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([]); ax.set_yticks([]); [s_.set_visible(False) for s_ in ax.spines.values()]

fig, axs = plt.subplots(1, 4, figsize=(16, 5.3)); fig.subplots_adjust(left=0.01, right=0.95, top=0.78, bottom=0.03, wspace=0.03)
for j, g in enumerate(PANELS): landscape(axs[j], g); axs[j].set_title(TITLE[g])
fig.suptitle("Gene-level diffusion: New COPD carries interferon, heme and histone modules shared with Existing, plus a New-specific ribosomal module",
             x=0.01, ha="left", fontsize=12, weight="bold")
fig.text(0.01, 0.86, f"Network propagation of limma moderated t over {n:,} genes linked by shared curated-set membership (Hallmark, Reactome, BioCarta, KEGG; Jaccard kNN k = 10). Layout: t-SNE of diffusion-map coordinates.\n"
         f"Field = diffused score at t = 1 scaled by its null SD; coloured points = multiscale calls (FDR 0.05, t in [0, {META['cap']}]) against a co-expression null estimated from the Excluded contrast (0 DE genes);\n"
         "black rings = raw limma FDR < 0.05; dashed outlines = labelled gene modules (hub genes listed).", fontsize=7, color="#444")
sm = plt.cm.ScalarMappable(cmap=CM, norm=plt.Normalize(-6, 6)); cb = fig.colorbar(sm, ax=axs, fraction=0.01, pad=0.01); cb.set_label("standardized score", fontsize=7)
fig.savefig(f"{FIG}/G1_gene_heat_landscape.png", dpi=300); fig.savefig(f"{FIG}/G1_gene_heat_landscape.pdf"); plt.close(fig)

# ---------- G2: module matrix ----------
M = MOD.copy(); M["spec"] = M["mean_T_New"] - M["mean_T_Existing"]
show = pd.concat([M[M.sel], M[~M.sel].sort_values("size", ascending=False).head(14)]).drop_duplicates("module")
show = show.sort_values("mean_T_New", ascending=False)
fig, ax = plt.subplots(figsize=(11.5, 0.27 * len(show) + 2.2)); fig.subplots_adjust(left=0.3, right=0.74, top=0.86, bottom=0.03)
for i, (_, r) in enumerate(show.iterrows()):
    for j, g in enumerate(PANELS):
        mt, fc = r[f"mean_T_{g}"], r[f"frac_called_{g}"]
        ax.scatter(j, i, s=12 + 260 * fc, c=[np.clip(mt, -6, 6)], cmap=CM, vmin=-6, vmax=6, ec=INK if fc >= 0.5 else "none", lw=0.5)
    ax.text(-0.6, i, f"{r.label}", ha="right", va="center", fontsize=7, weight="bold" if r.sel else "normal")
    ax.text(3.6, i, f"{r['size']}  {', '.join(eval(r.hubs)[:3])}", ha="left", va="center", fontsize=6.2, color="#555")
ax.set_xticks(range(4)); ax.set_xticklabels([TITLE[g].replace(" vs No COPD", "").replace(" (direct)", "") for g in PANELS], fontsize=7.5)
ax.xaxis.tick_top(); ax.set_yticks([]); ax.set_ylim(len(show) + 1.3, -0.8); ax.set_xlim(-0.5, 3.5); [s.set_visible(False) for s in ax.spines.values()]
ax.text(3.6, -0.8, "genes  hubs", fontsize=6.4, color="#555", weight="bold")
fig.suptitle("Gene modules: which are shared, which are New-specific, and which are Existing-only", x=0.02, ha="left", fontsize=11.5, weight="bold")
fig.text(0.02, 0.905, "Colour: mean diffused score of module genes at the scale cap; size: fraction of module genes called (FDR 0.05, co-expression null); ring: at least half called. "
         "Bold: module discriminating New (≥ 25% called in New or New vs Existing).", fontsize=6.8, color="#444", wrap=True)
for k, f in enumerate([0.1, 0.5, 1.0]): ax.scatter(0.6 + k * 0.9, len(show) + 0.7, s=12 + 260 * f, c="#bbb", clip_on=False, ec=INK if f >= 0.5 else "none", lw=0.5); ax.text(0.78 + k * 0.9, len(show) + 0.7, f"{f:.0%} called", va="center", fontsize=6.4)
fig.savefig(f"{FIG}/G2_gene_modules.png", dpi=300); fig.savefig(f"{FIG}/G2_gene_modules.pdf"); plt.close(fig)

# ---------- G3: calibration ----------
fig = plt.figure(figsize=(13.5, 4.4)); gs = fig.add_gridspec(1, 3, wspace=0.42, left=0.06, right=0.985, top=0.78, bottom=0.2)
ax = fig.add_subplot(gs[0]); k = np.arange(1, n + 1); zh2 = (S["z|Excluded"] @ np.load(f"{OUT}/gene_graph_{CLS}_{MODEL}.npz")["U"]) ** 2
ax.scatter(k, zh2, s=0.6, color="#ccc", lw=0, rasterized=True); ax.plot(k, P, color="#1f5fa8", lw=1.6, label="co-expression null spectrum (smoothed)")
ax.axhline(1, color="#B5446E", ls="--", lw=1.2, label="iid null (flat)"); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-3, 60)
ax.set_xlabel("graph eigenvector (smooth → rough)"); ax.set_ylabel("power of Excluded t-map"); ax.legend(frameon=False, fontsize=6.4, loc="upper right")
ax.set_title("a  Without any DE, t-scores are smooth on the graph")
ax = fig.add_subplot(gs[1])
for i, g in enumerate(PANELS):
    a = SUM[(SUM.null == "iid") & (SUM.panel == g)].iloc[0]; b = SUM[(SUM.null == "coexpr") & (SUM.panel == g)].iloc[0]
    ax.plot([b.calls_fdr, a.calls_fdr], [i, i], color="#DDD", lw=3); ax.plot(a.calls_fdr, i, "o", color="#B5446E", ms=6); ax.plot(b.calls_fdr, i, "o", color="#1f5fa8", ms=6)
    ax.plot(a.raw_fdr05, i, "|", color=INK, ms=9)
ax.set_yticks(range(4)); ax.set_yticklabels([TITLE[g] for g in PANELS], fontsize=7); ax.invert_yaxis(); ax.set_xscale("symlog", linthresh=10); ax.set_xlim(-1, 1.2e4); ax.set_xlabel("genes called (FDR 0.05, symlog)")
ax.legend(handles=[Line2D([], [], marker="o", ls="", color="#B5446E", label="iid null"), Line2D([], [], marker="o", ls="", color="#1f5fa8", label="co-expression null"),
                   Line2D([], [], marker="|", ls="", color=INK, ms=9, label="raw FDR < 0.05")], frameon=False, fontsize=6.4, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3)
ax.set_title("b  The iid null calls 254 genes where none are DE")
ax = fig.add_subplot(gs[2]); keys = list(POW)
for i, kk in enumerate(keys):
    ax.barh(i - 0.18, POW[kk]["raw_power"], 0.34, color="#bbb"); ax.barh(i + 0.18, POW[kk]["power"], 0.34, color="#1f5fa8")
    ax.text(POW[kk]["power"] + 0.01, i + 0.18, f"{POW[kk]['power']:.0%} vs raw {POW[kk]['raw_power']:.0%}  (false {POW[kk]['false_call_prop']:.0%})", va="center", fontsize=6.4)
ax.set_yticks(range(len(keys))); ax.set_yticklabels([kk.replace("g_shift", " genes, shift ") for kk in keys], fontsize=7); ax.invert_yaxis(); ax.set_xlim(0, 0.95)
ax.set_xlabel("fraction of planted genes recovered"); ax.set_title("c  Power on a planted module")
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color="#bbb", label="raw limma FDR"), plt.Rectangle((0, 0), 1, 1, color="#1f5fa8", label="diffusion, co-expression null")],
          frameon=False, fontsize=6.4, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
fig.suptitle("Gene-level calibration: genes are co-expressed, so the null must be too", x=0.06, ha="left", fontsize=12, weight="bold")
fig.text(0.06, 0.865, "The co-expression null uses the graph-Fourier power of the Excluded vs No COPD t-map (0 DE genes, λ 1.13). Under H0 the gene-gene correlation of moderated t equals residual "
         "expression correlation, which is shared by all contrasts.\nExcluded = 0 calls is therefore by construction; the evidence is panel b's iid row and panel c's power at ~4% false calls.", fontsize=6.9, color="#444")
fig.savefig(f"{FIG}/G3_gene_null_calibration.png", dpi=300); fig.savefig(f"{FIG}/G3_gene_null_calibration.pdf"); plt.close(fig)

# ---------- G4: gene-to-pathway bridge ----------
pst = np.load(f"{OUT}/state_{CLS}_{MODEL}.npz", allow_pickle=True); V = np.load(f"{OUT}/state_v2_{CLS}_{MODEL}.npz"); pw = list(pst["pw"])
cache = read_rds(MSIGDB_CACHE); so = cache.value[names(cache).index("sets")]
GS = {k: v.value for k, v in zip(names(so), so.value)}; gi = {g: i for i, g in enumerate(genes)}
memb = [np.array([gi[g] for g in set(GS[p]) if g in gi]) for p in pw]
fig, axs = plt.subplots(1, 4, figsize=(15, 4.2)); fig.subplots_adjust(left=0.05, right=0.99, top=0.76, bottom=0.14, wspace=0.28)
BR = {}
for j, g in enumerate(PANELS):
    gm = np.array([S[f"T1|{g}"][m].mean() if len(m) else np.nan for m in memb]); pz = V[f"z|singscore|{g}"]; ok = np.isfinite(gm)
    rho = spearmanr(gm[ok], pz[ok])[0]; pc = V[f"call|singscore|{g}"].astype(bool)
    agree = np.mean(np.sign(gm[pc & ok]) == V[f"sign|singscore|{g}"][pc & ok]) if (pc & ok).any() else np.nan; BR[g] = dict(rho=float(rho), sign_agree=float(agree), n_called=int(pc.sum()))
    ax = axs[j]; ax.scatter(gm[~pc], pz[~pc], s=3, color="#ccc", lw=0); ax.scatter(gm[pc], pz[pc], s=5, c=np.where(V[f"sign|singscore|{g}"][pc] > 0, "#e34948", "#2a78d6"), lw=0)
    ax.axhline(0, color="#ddd", lw=0.6); ax.axvline(0, color="#ddd", lw=0.6)
    ax.set_title(TITLE[g]); ax.set_xlabel("mean diffused gene score of members"); ax.set_ylabel("pathway singscore z" if j == 0 else "")
    ax.text(0.03, 0.95, f"Spearman ρ = {rho:.2f}" + (f"\nsign agreement in program calls {agree:.0%}" if np.isfinite(agree) else ""), transform=ax.transAxes, va="top", fontsize=6.6, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5), zorder=5)
fig.suptitle("Gene and pathway diffusion agree: pathway programs are carried by their member genes", x=0.05, ha="left", fontsize=12, weight="bold")
fig.text(0.05, 0.85, f"{len(pw):,} pathways (Reactome, BioCarta, KEGG). Coloured: pathway program calls (singscore, overlap-preserving null; red higher, blue lower). Grey: not called.", fontsize=7, color="#444")
fig.savefig(f"{FIG}/G4_gene_pathway_bridge.png", dpi=300); fig.savefig(f"{FIG}/G4_gene_pathway_bridge.pdf"); plt.close(fig)
json.dump(BR, open(f"{OUT}/gene_pathway_bridge_{CLS}_{MODEL}.json", "w"), indent=1)
print(BR); print("gene figures done")
