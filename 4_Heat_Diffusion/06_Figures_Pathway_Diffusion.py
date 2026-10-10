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
import heat_kernel as core, calibrated_null as core2

from paths import DE_FILE, SINGSCORE_FILE, FGSEA_FILE, MSIGDB_CACHE, OUT, FIG, CLS, MODEL
S = np.load(f"{OUT}/state_{CLS}_{MODEL}.npz", allow_pickle=True); V = np.load(f"{OUT}/state_v2_{CLS}_{MODEL}.npz")
M1 = json.load(open(f"{OUT}/diffusion_meta_{CLS}_{MODEL}.json")); M2 = json.load(open(f"{OUT}/diffusion_meta_v2_{CLS}_{MODEL}.json"))
pw = list(S["pw"]); n = len(pw); P = S["P"]; lam, U, TG, W = S["lam"], S["U"], S["TG"], S["W"]; prog_of = V["prog_of"]; PROG = M2["programs"]
C = V["C"].astype(float); F = core2.null_factor(C)
PANELS = ["New", "Excluded", "Existing", "New vs Existing"]; TITLE = {"New": "New vs No COPD", "Excluded": "Excluded vs No COPD", "Existing": "Existing vs No COPD", "New vs Existing": "New vs Existing (direct)"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.titlesize": 8.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
                     "axes.spines.top": False, "axes.spines.right": False})
CM = LinearSegmentedColormap.from_list("d", ["#2a78d6", "#f0efec", "#e34948"]); INK = "#222"; GR = "#9A9890"
gx, gy = np.mgrid[0:1:220j, 0:1:220j]
d2 = (gx.ravel()[:, None] - P[:, 0]) ** 2 + (gy.ravel()[:, None] - P[:, 1]) ** 2; K = np.exp(-d2 / (2 * 0.022 ** 2)); Ks = K.sum(1); dens = Ks.reshape(gx.shape)
def field(v): return (K @ v / (Ks + 1e-12)).reshape(gx.shape)
H1 = core.heat(lam, U, 1.0); s1 = np.sqrt(((H1 @ F) ** 2).sum(1))
SUMM = pd.read_csv(f"{OUT}/diffusion_summary_v2_{CLS}_{MODEL}.csv").set_index(["method", "panel"])
ABBR = {"Growth-factor RTK / PI3K": "RTK/PI3K", "MAPK / Ras / ERK": "MAPK/Ras/ERK", "Proteasome / ubiquitin": "Proteasome", "Innate / neutrophil": "Innate",
        "Cytokine / NF-κB": "Cytokine/NF-κB", "Chromatin / methylation": "Chromatin", "DNA repair / TP53": "DNA repair", "Translation / ribosome": "Translation",
        "GPCR / Ca2+ / NFAT": "GPCR/Ca²⁺", "Trafficking / vesicles": "Trafficking", "Mitochondria / OXPHOS": "OXPHOS", "Cell cycle / mitosis": "Cell cycle", "RNA processing": "RNA proc."}
def short_lab(l): return " + ".join(ABBR.get(x, x) for x in l.split(" + "))
top_prog = [p["program"] for p in sorted(PROG, key=lambda p: -p["size"])[:7]]

def landscape(ax, m, g, calibrated=True):
    z = V[f"z|{m}|{g}"]; Tf = np.clip((H1 @ z) / s1, -6, 6); Fv = field(Tf); fade = np.clip(dens / np.median(dens), 0, 1)
    rgba = CM((np.clip(Fv, -6, 6) + 6) / 12); rgba[..., 3] = 0.15 + 0.8 * fade
    ax.imshow(np.transpose(rgba, (1, 0, 2)), origin="lower", extent=[0, 1, 0, 1], interpolation="bilinear")
    call = V[f"call|{m}|{g}"].astype(bool); sign = V[f"sign|{m}|{g}"]
    for sgn, col in [(1, "#7a1010"), (-1, "#0b3a75")]:
        ind = field((call & (sign == sgn)).astype(float)) * (fade > 0.45)
        if ind.max() > 0.42: ax.contour(gx, gy, ind, levels=[0.42], colors=col if calibrated else GR, linewidths=0.9, linestyles="-" if calibrated else "--")
    fdr = S[f"fdr|{m}|{g}"]; sz = 1.5 + 5 * np.minimum(-np.log10(np.clip(fdr, 1e-10, 1)), 10) / 10 * 2
    o = np.argsort(np.abs(z))
    ax.scatter(P[o, 0], P[o, 1], s=sz[o], c=np.clip(z[o], -6, 6), cmap=CM, vmin=-6, vmax=6, ec=np.where(fdr[o] < .05, INK, "none"), lw=0.35, zorder=3)
    s_ = SUMM.loc[(m, g)]
    note = f"program calls {int(s_.calls_overlap_fdr)}   raw {int(s_.raw_fdr05)}   λ {s_.inflation_lambda:.1f}" if calibrated else f"raw {int(s_.raw_fdr05)}   λ {s_.inflation_lambda:.1f}   (outlines not calibrated)"
    ax.text(0.01, 0.01, note, transform=ax.transAxes, fontsize=6.3, color=INK, bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([]); ax.set_yticks([]); [s.set_visible(False) for s in ax.spines.values()]

def labels(ax):
    placed = []
    for k in top_prog:
        idx = np.where(prog_of == k)[0]
        if len(idx) < 5: continue
        c = np.median(P[idx], 0); d = c - 0.5; d = d / (np.linalg.norm(d) + 1e-9); tgt = np.clip(c + d * 0.17, [0.1, 0.06], [0.9, 0.95])
        for _ in range(20):
            if all(np.hypot(*(tgt - q)) > 0.09 for q in placed): break
            tgt = tgt + np.array([0, 0.035])
        placed.append(tgt); lab = short_lab(PROG[k]["label"])
        ax.annotate(lab, c, tgt, fontsize=5.6, ha="center", color=INK, zorder=6, arrowprops=dict(arrowstyle="-", lw=0.45, color=INK),
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))

# ---------------- Figure D1: landscape ----------------
fig, axs = plt.subplots(2, 4, figsize=(15, 8.2)); fig.subplots_adjust(left=0.05, right=0.95, top=0.845, bottom=0.07, wspace=0.05, hspace=0.16)
for j, g in enumerate(PANELS):
    landscape(axs[0, j], "singscore", g, True); axs[0, j].set_title(TITLE[g])
    landscape(axs[1, j], "FGSEA", g, False)
    if j == 0: labels(axs[0, j])
    if g == "Existing": labels(axs[0, j])
axs[0, 0].text(-0.08, 0.5, "singscore\n(calibrated)", transform=axs[0, 0].transAxes, rotation=90, ha="center", va="center", fontsize=8.5, weight="bold")
axs[1, 0].text(-0.08, 0.5, "FGSEA\n(descriptive)", transform=axs[1, 0].transAxes, rotation=90, ha="center", va="center", fontsize=8.5, weight="bold")
fig.suptitle("New and Existing COPD share programs (innate/NF-κB and proteasome up; MAPK, DNA repair, cell cycle down); Excluded has none", x=0.05, ha="left", fontsize=12.5, weight="bold")
fig.text(0.05, 0.9, f"Pathway heat-diffusion landscape, {CLS} {MODEL}. n = {n} sets (Reactome, BioCarta, KEGG) scored by both methods; Jaccard kNN graph (k = 10); scale cap t = {M1['cap']}.\n"
         "Field = diffused score at t = 1; points = raw score (ring: raw FDR < 0.05); outlines = multiscale program calls (FDR 0.05) against an overlap-preserving null.", fontsize=7.2, color="#444")
sm = plt.cm.ScalarMappable(cmap=CM, norm=plt.Normalize(-6, 6)); cb = fig.colorbar(sm, ax=axs, fraction=0.012, pad=0.01); cb.set_label("standardized score (− lower, + higher)", fontsize=7)
fig.text(0.05, 0.02, "FGSEA row: preranked p-values are inflated (λ 4 to 7.5, including Excluded, where singscore λ = 0.71 and no gene is DE), so its outlines are shown dashed and are not inference. "
         "Program labels: themes enriched among member sets relative to all 1,232 tested sets.", fontsize=6.8, color="#555")
fig.savefig(f"{FIG}/D1_pathway_heat_landscape.png", dpi=300); fig.savefig(f"{FIG}/D1_pathway_heat_landscape.pdf"); plt.close(fig)

# ---------------- Figure D2: persistence barcode ----------------
cols = [("singscore", g) for g in PANELS] + [("FGSEA", g) for g in PANELS]
rows = []
for k in top_prog[:6]:
    idx = np.where(prog_of == k)[0]
    sc = np.zeros(len(idx))
    for (m, g) in cols[:4]: sc += V[f"call|{m}|{g}"][idx]
    best = idx[np.argsort(-(sc + 0.01 * np.abs(V["z|singscore|Existing"][idx])))[:5]]
    rows.append(("hdr", k)); rows += [("pw", i) for i in best]
fig, ax = plt.subplots(figsize=(12.5, 0.2 * len(rows) + 1.6)); fig.subplots_adjust(left=0.33, right=0.985, top=0.86, bottom=0.08)
wcol = 1.0; y = 0
for kind, v in rows:
    if kind == "hdr":
        y += 0.5; ax.text(-0.02, y, PROG[v]["label"], transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7.4, weight="bold"); y += 1; continue
    nm = re.sub(r"^(REACTOME|BIOCARTA|KEGG_MEDICUS_REFERENCE|KEGG_MEDICUS)_", "", pw[v]).replace("_", " ").lower()
    ax.text(-0.02, y, nm[:58], transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.4, color="#333")
    for c, (m, g) in enumerate(cols):
        per = V[f"persist|{m}|{g}"][v]; call = V[f"call|{m}|{g}"][v]; sgn = V[f"sign|{m}|{g}"][v]
        col = "#c0392b" if sgn > 0 else "#1f5fa8"
        if call and per.any():
            ts = TG[per]; x0 = c + 0.06 + ts.min() / TG.max() * 0.84; x1 = c + 0.06 + ts.max() / TG.max() * 0.84 + 0.02
            ax.add_patch(plt.Rectangle((x0, y - 0.32), x1 - x0, 0.64, color=col, alpha=0.9 if m == "singscore" else 0.45, lw=0))
            bt = V[f"best_t|{m}|{g}"][v]; xb = c + 0.06 + bt / TG.max() * 0.84; ax.plot([xb, xb], [y - 0.4, y + 0.4], color=INK, lw=0.8)
        if S[f"fdr|{m}|{g}"][v] < .05: ax.plot(c + 0.95, y, "o", ms=2.2, color=INK)
    y += 1
for c in range(len(cols) + 1): ax.axvline(c, color="#DDD", lw=0.6)
ax.axvline(4, color=INK, lw=1)
for c, (m, g) in enumerate(cols): ax.text(c + 0.5, -0.6, TITLE[g].replace(" (direct)", "").replace(" vs No COPD", ""), ha="center", fontsize=6.6, rotation=0)
ax.text(2, -1.6, "singscore (calibrated)", ha="center", fontsize=8, weight="bold"); ax.text(6, -1.6, "FGSEA (descriptive)", ha="center", fontsize=8, weight="bold", color="#777")
ax.set_xlim(0, len(cols)); ax.set_ylim(y, -2.2); ax.set_xticks([]); ax.set_yticks([]); [s.set_visible(False) for s in ax.spines.values()]
fig.suptitle("Programs persist across diffusion scales in New and Existing, never in Excluded", x=0.02, ha="left", fontsize=11.5, weight="bold")
fig.text(0.02, 0.9, f"Bars: scales t in [0, {M1['cap']}] where the program call holds (red higher, blue lower than reference); tick: most significant scale; dot: raw FDR < 0.05. "
         "A bar from t = 0 is strong alone; a bar starting later is weak alone but part of a coherent program.", fontsize=7, color="#444")
fig.savefig(f"{FIG}/D2_persistence_barcode.png", dpi=300); fig.savefig(f"{FIG}/D2_persistence_barcode.pdf"); plt.close(fig)

# ---------------- Figure D3: null comparison + t-sweep ----------------
cal = M2["calibration_overlap"]; mor = M1["moran"]
fig = plt.figure(figsize=(13.5, 7.6)); gsp = fig.add_gridspec(2, 3, hspace=0.75, wspace=0.55, left=0.1, right=0.985, top=0.86, bottom=0.08)
ax = fig.add_subplot(gsp[0, 0]); keys = [f"{m}|{g}" for m in ("singscore", "FGSEA") for g in PANELS]
ax.axvspan(cal["moran_I_overlap_null_mean"] - 0.03, cal["moran_I_overlap_null_q95"], color="#E6E4DD"); ax.axvline(cal["moran_I_overlap_null_mean"], color=GR, lw=1)
for i, k in enumerate(keys): ax.plot(mor[k][0], i, "o", color="#1f5fa8" if k.startswith("sing") else "#999", ms=6)
ax.set_yticks(range(len(keys))); ax.set_yticklabels([k.replace("|", " · ").replace(" vs Existing", " v Exist.") for k in keys], fontsize=6.6); ax.invert_yaxis()
ax.set_xlabel("Moran's I of raw scores on the Jaccard graph"); ax.set_title("a  Smoothness is mostly gene sharing")
ax.set_ylim(len(keys) - 0.5, -1.4); ax.text(cal["moran_I_overlap_null_mean"], -0.75, "pure overlap noise (mean to q95)", ha="center", va="center", fontsize=6.2, color="#666")
ax = fig.add_subplot(gsp[0, 1:]); Ssum = SUMM.reset_index()
for i, r in Ssum.iterrows():
    ax.plot([r.calls_overlap_fdr, r.calls_label_null_fdr], [i, i], color="#DDD", lw=3, zorder=1)
    ax.plot(r.calls_label_null_fdr, i, "o", color="#B5446E", ms=6, zorder=2); ax.plot(r.calls_overlap_fdr, i, "o", color="#1f5fa8", ms=6, zorder=3); ax.plot(r.raw_fdr05, i, "|", color=INK, ms=9, zorder=4)
ax.set_yticks(range(len(Ssum))); ax.set_yticklabels([f"{r.method} · {r.panel}" for r in Ssum.itertuples()], fontsize=6.8); ax.invert_yaxis(); ax.set_xlabel("pathways called (FDR 0.05)")
ax.legend(handles=[Line2D([], [], marker="o", ls="", color="#B5446E", label=f"label-permutation null (recipe as written): {cal['label_null_fwer_rate']:.0%} FWER, {cal['label_null_mean_fdr_calls']:.0f} calls on pure overlap noise"),
                   Line2D([], [], marker="o", ls="", color="#1f5fa8", label=f"overlap-preserving null: {cal['overlap_null_fwer_rate']:.0%} FWER, {cal['overlap_null_mean_fdr_calls']:.2f} calls on pure overlap noise"),
                   Line2D([], [], marker="|", ls="", color=INK, ms=9, label="raw FDR < 0.05")], fontsize=6.6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
ax.set_title("b  The label-permutation null calls programs in noise; the overlap null does not")
SW = pd.read_csv(f"{OUT}/diffusion_tsweep_{CLS}_{MODEL}.csv")
for j, (m, a, b) in enumerate([("singscore", "New", "Existing"), ("singscore", "Excluded", "Existing"), ("FGSEA", "Excluded", "Existing")]):
    ax = fig.add_subplot(gsp[1, j]); d = SW[(SW.method == m) & (SW.a == a) & (SW.b == b)]
    ax.axvspan(M1["cap"], 4, color="#EEE", lw=0); ax.fill_between(d.t, d.msr_lo, d.msr_hi, color="#C3DBF2", lw=0, label="spectral null 95%")
    ax.plot(d.t, d.lab_lo, ":", color=GR, lw=1, label="label null 95%"); ax.plot(d.t, d.lab_hi, ":", color=GR, lw=1)
    ax.plot(d.t, d.r, color="#B5446E", lw=2, label="observed r(t)"); ax.axvline(1, color=INK, ls="--", lw=0.8); ax.axhline(0, color="#DDD", lw=0.6)
    r1 = d[d.t == 1.0].iloc[0]; ax.text(0.03, 0.06, f"r(1) = {r1.r:.2f}, spectral p = {r1.p_msr:.3f}, raw r = {r1.r_raw:.2f}", transform=ax.transAxes, fontsize=6.6)
    ax.set_ylim(-1, 1.02); ax.set_xlim(0, 4); ax.set_xlabel("diffusion time t"); ax.set_ylabel("correlation of diffused maps" if j == 0 else "")
    ax.set_title(f"{'cde'[j]}  {m}: {a} vs {b}"); ax.text(3.1, 0.9, "beyond cap", fontsize=6.2, color=GR, ha="center")
    if j == 0: ax.legend(fontsize=6.2, frameon=False, loc="lower right", bbox_to_anchor=(1, 0.12))
fig.suptitle("Which null is right: pathway scores are already smooth on the gene-set graph before any biology", x=0.06, ha="left", fontsize=12, weight="bold")
fig.text(0.06, 0.905, f"Calibration on {cal['reps']} simulated null maps z ~ N(0, C), C = cosine gene overlap. Spectral agreement null on the same noise: Type I error {cal['agreement_typeI_spectral_on_overlap_noise']:.1%}. "
         f"Planted program: {cal['overlap_null_planted_false_call_prop']:.1%} of calls outside the program and its neighbours.", fontsize=7.2, color="#444")
fig.savefig(f"{FIG}/D3_null_calibration_tsweep.png", dpi=300); fig.savefig(f"{FIG}/D3_null_calibration_tsweep.pdf"); plt.close(fig)
print("figures done")
