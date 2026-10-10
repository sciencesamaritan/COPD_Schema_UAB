"""Gene-level heat diffusion (network propagation) of limma moderated t on a gene graph built from
shared curated-set membership (Hallmark, Reactome, BioCarta, KEGG_MEDICUS; sets of 10 to 200 measured genes).
Null: z_i ~ N(0, 1) independent (moderated t, large df). Excluded vs No COPD (0 DE genes) is the negative control
that tests this independence assumption; co-expression would show up there as spurious module calls."""
import sys, os, json, time, re
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd, networkx as nx
from scipy.stats import norm
import heat_kernel as core, calibrated_null as core2
from rds_reader import read_rds, names

from paths import DE_FILE, SINGSCORE_FILE, FGSEA_FILE, MSIGDB_CACHE, OUT, FIG, CLS, MODEL
B = 1000
t0 = time.time(); log = lambda *a: print(f"[{time.time()-t0:6.1f}s]", *a, flush=True)
CT = {"New": f"{CLS}_1vs0", "Excluded": f"{CLS}_2vs0", "Existing": f"{CLS}_3vs0", "New vs Existing": f"{CLS}_1vs3"}

de = pd.read_csv(DE_FILE, usecols=["gene_label", "Contrast", "t", "adj.P.Val", "logFC"])
de = de[de.Contrast.isin(CT.values()) & de.gene_label.notna()]
de = de.loc[de.groupby(["Contrast", "gene_label"])["t"].apply(lambda s: s.abs().idxmax()).values]
cache = read_rds(MSIGDB_CACHE); so = cache.value[0]; src = cache.value[1].value
GS = {k: v.value for k, v in zip(names(so), so.value)}; measured = set(de.gene_label)
sets = [k for k, s in zip(names(so), src) if s in ("Hallmark", "Reactome", "BioCarta", "KEGG_MEDICUS") and 10 <= len(set(GS[k]) & measured) <= 200]
mem = {}
for k in sets:
    for g in set(GS[k]) & measured: mem.setdefault(g, []).append(k)
genes = sorted(mem); n = len(genes); log("genes", n, "sets", len(sets))
GC = f"{OUT}/gene_graph_{CLS}_{MODEL}.npz"
if os.path.exists(GC):
    q = np.load(GC); W, lam, U, cap = q["W"], q["lam"], q["U"], float(q["cap"]); log("graph cached")
else:
    W, L, lam, U = core.graph(genes, mem, k=10); cap = core.scale_cap(lam, U, target=20)
    np.savez(GC, W=W, lam=lam, U=U, cap=cap); log("graph built")
G = nx.from_numpy_array(W); ncomp = nx.number_connected_components(G); log("components", ncomp, "cap", cap)
TG = np.round(np.arange(0, cap + 1e-9, 0.25), 2)

Z, FDR, LFC = {}, {}, {}
for g, ct in CT.items():
    d = de[de.Contrast == ct].set_index("gene_label").loc[genes]; Z[g] = d.t.values; FDR[g] = d["adj.P.Val"].values; LFC[g] = d.logFC.values
lamGC = {g: core2.inflation(z) for g, z in Z.items()}; log("gene lambda", lamGC)
MOR = {g: core.moran(z, W, B=499, rng=np.random.default_rng(3)) for g, z in Z.items()}; log("moran", MOR)
ZH = {g: U.T @ z for g, z in Z.items()}

# Null spectra. iid: P = 1 (independent genes). coexpr: graph-Fourier power of the negative-control contrast
# (Excluded vs No COPD, 0 DE genes), smoothed over eigenvalue order. Under H0 the correlation of moderated t
# between genes equals residual expression correlation, which does not depend on the contrast, so one
# no-signal contrast estimates the null spectrum for all. Conservative: any real Excluded signal inflates P.
def smooth_power(zh, w=201):
    p = zh ** 2; k = np.ones(w) / w; pad = np.pad(p, (w // 2, w // 2), mode="reflect"); return np.convolve(pad, k, mode="valid")
NULLS = {"iid": np.ones(n), "coexpr": smooth_power(ZH["Excluded"])}
log("coexpr spectrum: low-freq (first 200) mean power %.2f, high-freq (last 2000) %.2f" % (NULLS["coexpr"][:200].mean(), NULLS["coexpr"][-2000:].mean()))

def null_node(P, B, rng):
    Nh = np.sqrt(P)[:, None] * rng.normal(size=(n, B)); out = np.ones((n, B)); S_ = {}
    for t in TG:
        f = np.exp(-t * lam); s = np.sqrt((U ** 2) @ (f ** 2 * P)); S_[t] = s
        out = np.minimum(out, 2 * norm.sf(np.abs((U @ (f[:, None] * Nh)) / s[:, None])))
    return out, S_

def multiscale_spec(zh, nul, S_):
    nT = len(TG); T = np.empty((n, nT)); obs = np.empty((n, nT))
    for j, t in enumerate(TG):
        f = np.exp(-t * lam); To = (U @ (f * zh)) / S_[t]; T[:, j] = To; obs[:, j] = 2 * norm.sf(np.abs(To))
    minp = obs.min(1); bi = obs.argmin(1); pool = np.sort(nul.ravel()); M = pool.size
    p = (np.searchsorted(pool, minp, side="right") + 1) / (M + 1); q = core2._bh(p); na = np.sort(nul.min(0))
    pf = (np.searchsorted(na, minp, side="right") + 1) / (len(na) + 1); call = q < .05; cut = minp[call].max() if call.any() else 0
    return dict(T=T, minp=minp, q=q, call=call, call_fwer=pf < .05, sign=np.sign(T[np.arange(n), bi]), best_t=TG[bi], persist=obs <= cut)

ALL = {}
for nm, P in NULLS.items():
    nul, S_ = null_node(P, B, np.random.default_rng(7)); log("null draws", nm)
    for g in CT:
        r = multiscale_spec(ZH[g], nul, S_); ALL[(nm, g)] = r
        log(f"[{nm:6s}] {g:16s} FDR {int(r['call'].sum()):5d} FWER {int(r['call_fwer'].sum()):5d} raw {int((FDR[g]<.05).sum()):5d} rescued {int((r['call'] & (FDR[g]>=.05)).sum()):5d} up {int((r['call']&(r['sign']>0)).sum())} down {int((r['call']&(r['sign']<0)).sum())}")
    del nul
RES = {g: ALL[("coexpr", g)] for g in CT}

# planted-module power under the co-expression null: shift 40 or 100 graph-adjacent genes of the
# negative-control (Excluded) map and record recovery, recovery by gene-wise FDR, and the proportion of calls
# outside the module (strict) and outside the module and its direct graph neighbours (as for pathway programs)
nul, S_ = null_node(NULLS["coexpr"], B, np.random.default_rng(7)); rp = np.random.default_rng(13); POW = {}
for m_size in (40, 100):
    for shift in (1.0, 1.5, 2.0):
        pw_, raw_, fcp, fcn = [], [], [], []
        for rep in range(10):
            seed = int(rp.integers(n)); mod = list(nx.single_source_shortest_path_length(G, seed))[:m_size]
            z = Z["Excluded"].copy(); z[mod] += shift
            r = multiscale_spec(U.T @ z, nul, S_); q_raw = core2._bh(2 * norm.sf(np.abs(z)))
            pw_.append(r["call"][mod].mean()); raw_.append((q_raw[mod] < .05).mean())
            c = np.where(r["call"])[0]; nb = set(mod).union(*[set(G.neighbors(i)) for i in mod])
            fcp.append(0.0 if len(c) == 0 else float(np.mean(~np.isin(c, mod))))
            fcn.append(0.0 if len(c) == 0 else float(np.mean(~np.isin(c, list(nb)))))
        POW[f"{m_size}g_shift{shift}"] = dict(power=float(np.mean(pw_)), raw_power=float(np.mean(raw_)),
                                              false_call_prop=float(np.mean(fcp)), false_call_prop_beyond_neighbours=float(np.mean(fcn)))
        log("power", m_size, shift, POW[f"{m_size}g_shift{shift}"])
json.dump(POW, open(f"{OUT}/gene_power_coexpr.json", "w"), indent=1)
del nul

# modules: communities among genes called in any vs-0 panel or the direct contrast
called = np.zeros(n, bool)
for r in RES.values(): called |= r["call"]
idx = np.where(called)[0]; Gs = G.subgraph(idx).copy()
comms = [sorted(c) for c in nx.community.greedy_modularity_communities(Gs, weight="weight", resolution=1.0) if len(c) >= 8]
comms.sort(key=len, reverse=True)
from scipy.stats import hypergeom
setgenes = {k: set(GS[k]) & set(genes) for k in sets}
def annotate(c):
    cg = {genes[i] for i in c}; best = None
    for k, sg in setgenes.items():
        ov = len(cg & sg)
        if ov < 3: continue
        p = hypergeom.sf(ov - 1, n, len(sg), len(cg))
        if best is None or p < best[1]: best = (k, p, ov)
    return best
MOD = []; mod_of = -np.ones(n, int)
for k, c in enumerate(comms[:40]):
    mod_of[c] = k; sub = Gs.subgraph(c); dg = dict(sub.degree(weight="weight")); hubs = sorted(c, key=lambda i: -dg[i])[:4]
    a = annotate(c); lab = re.sub(r"^(REACTOME|HALLMARK|BIOCARTA|KEGG_MEDICUS_REFERENCE|KEGG_MEDICUS)_", "", a[0]).replace("_", " ").lower() if a else ""
    MOD.append(dict(module=k, size=len(c), hubs=[genes[i] for i in hubs], best_set=a[0] if a else "", best_p=float(a[1]) if a else 1.0, overlap=a[2] if a else 0,
                    label=(lab[:1].upper() + lab[1:])[:40],
                    **{f"frac_called_{g}": float(RES[g]["call"][c].mean()) for g in CT}, **{f"mean_T_{g}": float(np.mean(RES[g]["T"][c, -1])) for g in CT}))
log("modules", [(m["label"], m["size"], m["hubs"][:2]) for m in MOD[:15]])

# embedding: diffusion-map coordinates (t = cap) -> t-SNE
from sklearn.manifold import TSNE
k = 40; Dm = U[:, 1:k + 1] * np.exp(-cap * lam[1:k + 1])
Y = TSNE(n_components=2, perplexity=40, init="pca", random_state=1, learning_rate="auto").fit_transform(Dm / np.linalg.norm(Dm, axis=1, keepdims=True).clip(1e-12))
Y = (Y - Y.min(0)) / (Y.max(0) - Y.min(0)) * 0.9 + 0.05; log("embedding done")

rows = []
for g, r in RES.items():
    for i in np.where(r["call"] | (FDR[g] < .05))[0]:
        rows.append(dict(panel=g, gene=genes[i], module=int(mod_of[i]), module_label=MOD[mod_of[i]]["label"] if mod_of[i] >= 0 else "", t=round(float(Z[g][i]), 3), logFC=float(LFC[g][i]),
                         raw_padj=float(FDR[g][i]), min_p=float(r["minp"][i]), q_multiscale=float(r["q"][i]), call=bool(r["call"][i]), call_fwer=bool(r["call_fwer"][i]),
                         rescued=bool(r["call"][i] and FDR[g][i] >= .05), sign=int(r["sign"][i]), best_t=float(r["best_t"][i])))
pd.DataFrame(rows).to_csv(f"{OUT}/gene_diffusion_calls_{CLS}_{MODEL}.csv", index=False)
pd.DataFrame(MOD).to_csv(f"{OUT}/gene_modules_{CLS}_{MODEL}.csv", index=False)
summ = pd.DataFrame([dict(null=nm, panel=g, raw_fdr05=int((FDR[g] < .05).sum()), calls_fdr=int(r["call"].sum()), calls_fwer=int(r["call_fwer"].sum()), up=int((r["call"] & (r["sign"] > 0)).sum()),
                          down=int((r["call"] & (r["sign"] < 0)).sum()), rescued=int((r["call"] & (FDR[g] >= .05)).sum()), lambda_gc=round(lamGC[g], 3), moran_I=round(MOR[g][0], 3)) for (nm, g), r in ALL.items()])
summ.to_csv(f"{OUT}/gene_diffusion_summary_{CLS}_{MODEL}.csv", index=False); print(summ.to_string(index=False))
json.dump(dict(n=n, sets=len(sets), cap=cap, TG=list(map(float, TG)), null_spectrum_coexpr=NULLS['coexpr'][::50].tolist(), components=ncomp, lambda_gc=lamGC, moran={g: list(v) for g, v in MOR.items()}, modules=MOD), open(f"{OUT}/gene_meta_{CLS}_{MODEL}.json", "w"), indent=1, default=float)
np.savez_compressed(f"{OUT}/gene_state_{CLS}_{MODEL}.npz", Y=Y, genes=np.array(genes), mod_of=mod_of, TG=TG, W_edges=np.array(np.nonzero(np.triu(W > 0.05))).T,
                    **{f"z|{g}": Z[g] for g in CT}, **{f"fdr|{g}": FDR[g] for g in CT}, **{f"{k}|{g}": RES[g][k] for g in CT for k in ("call", "sign", "persist", "q", "best_t")},
                    **{f"Tcap|{g}": RES[g]["T"][:, -1] for g in CT}, **{f"T1|{g}": RES[g]["T"][:, list(TG).index(1.0)] for g in CT},
                    **{f"iidcall|{g}": ALL[("iid", g)]["call"] for g in CT}, lam=lam, P_coexpr=NULLS["coexpr"])
log("done")
