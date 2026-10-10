"""Pathway diffusion, v2: overlap-preserving null (primary) vs label-permutation null (recipe as written).
Reuses graph/layout/label-null results from run_pathway.py (state_*.npz)."""
import sys, os, json, time, re
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd, networkx as nx
from scipy.stats import norm
import heat_kernel as core, calibrated_null as core2
from rds_reader import read_rds, names
from pathway_themes import label_program_enriched

from paths import DE_FILE, SINGSCORE_FILE, FGSEA_FILE, MSIGDB_CACHE, OUT, FIG, CLS, MODEL
B = 2000
t0 = time.time(); log = lambda *a: print(f"[{time.time()-t0:6.1f}s]", *a, flush=True)
S = np.load(f"{OUT}/state_{CLS}_{MODEL}.npz", allow_pickle=True)
pw = list(S["pw"]); n = len(pw); lam, U, W, TG, P = S["lam"], S["U"], S["W"], S["TG"], S["P"]
CT = {"New": f"{CLS}_1vs0", "Excluded": f"{CLS}_2vs0", "Existing": f"{CLS}_3vs0", "New vs Existing": f"{CLS}_1vs3"}
VS0 = ["New", "Excluded", "Existing"]

fg = pd.read_csv(FGSEA_FILE, usecols=["pathway", "pval", "padj", "NES", "Contrast"])
ss = pd.read_csv(SINGSCORE_FILE)
de = pd.read_csv(DE_FILE, usecols=["gene_label"])
cache = read_rds(MSIGDB_CACHE); so = cache.value[names(cache).index("sets")]
GS = {k: v.value for k, v in zip(names(so), so.value)}; measured = set(de.gene_label.dropna().astype(str))
gs = {p: sorted(set(GS[p]) & measured) for p in pw}

def piv(df, val): return df[df.Contrast.isin(CT.values())].pivot_table(index="pathway", columns="Contrast", values=val, aggfunc="first").loc[pw]
fgz = fg.assign(z=np.sign(fg.NES) * norm.isf(fg.pval.clip(1e-300) / 2)); ssz = ss.assign(z=np.sign(ss.logFC) * norm.isf(ss["P.Value"].clip(1e-300) / 2))
Zf, Zs = piv(fgz, "z"), piv(ssz, "z")
Z = {("FGSEA", g): Zf[CT[g]].values for g in CT} | {("singscore", g): Zs[CT[g]].values for g in CT}
FDR = {(m, g): S[f"fdr|{m}|{g}"] for (m, g) in Z}

C = core2.overlap_corr(pw, gs); F = core2.null_factor(C); log("overlap corr built; median offdiag", float(np.median(C[np.triu_indices(n, 1)])))
LAMBDA = {f"{m}|{g}": core2.inflation(z) for (m, g), z in Z.items()}; log("inflation", {k: round(v, 2) for k, v in LAMBDA.items()})

# realistic null: how smooth is pure overlap noise?
rr = np.random.default_rng(5); mor_null = [core.moran(core2.draws(F, 1, rr)[:, 0], W, B=199, rng=rr)[0] for _ in range(40)]
log("Moran I under overlap null: mean %.3f, 95%% %.3f" % (np.mean(mor_null), np.quantile(mor_null, .95)))

# ---------------- calibration: both nulls on overlap-correlated noise ----------------
cal = {"moran_I_overlap_null_mean": float(np.mean(mor_null)), "moran_I_overlap_null_q95": float(np.quantile(mor_null, .95))}
lab_any, ov_any, lab_n, ov_n = [], [], [], []
for k in range(40):
    z0 = core2.draws(F, 1, rr)[:, 0]
    a = core.multiscale(z0, lam, U, TG, B=400, rng=rr); b = core2.multiscale_overlap(z0, lam, U, F, TG, B=400, rng=rr)
    lab_any.append(a["call_fwer"].any()); ov_any.append(b["call_fwer"].any()); lab_n.append(int(a["call"].sum())); ov_n.append(int(b["call"].sum()))
cal.update(label_null_fwer_rate=float(np.mean(lab_any)), label_null_mean_fdr_calls=float(np.mean(lab_n)),
           overlap_null_fwer_rate=float(np.mean(ov_any)), overlap_null_mean_fdr_calls=float(np.mean(ov_n)), reps=40)
log("calibration on overlap noise", cal)
G = nx.from_numpy_array(W); fdp = []
for k in range(25):
    seed = rr.integers(n); prog = list(nx.single_source_shortest_path_length(G, seed, cutoff=2))[:20]
    z0 = core2.draws(F, 1, rr)[:, 0]; z0[prog] += 1.5
    c = core2.multiscale_overlap(z0, lam, U, F, TG, B=400, rng=rr)["call"]; called = np.where(c)[0]
    # a call counts as false only if it is neither planted nor a direct neighbour of a planted set
    near = set(prog) | {j for i in prog for j in G.neighbors(i)}
    fdp.append(0 if len(called) == 0 else np.mean([i not in near for i in called]))
cal["overlap_null_planted_false_call_prop"] = float(np.mean(fdp)); cal["planted_reps"] = 25
rej = 0
for _ in range(200):
    a, b = core2.draws(F, 2, rr).T; q = core.agreement(a, b, lam, U, 1.0, B=200, rng=rr); rej += q["p_msr"] < .05
cal["agreement_typeI_spectral_on_overlap_noise"] = rej / 200
log("calibration done", cal)

# ---------------- real calls ----------------
RES = {}
for (m, g), z in Z.items():
    r = core2.multiscale_overlap(z, lam, U, F, TG, B=B, rng=np.random.default_rng(7)); RES[(m, g)] = r
    log(f"overlap-null calls {m:9s} {g:16s} FDR {int(r['call'].sum()):4d} FWER {int(r['call_fwer'].sum()):4d} | label-null {int(S[f'call|{m}|{g}'].sum()):4d} | raw {int((FDR[(m,g)]<.05).sum()):4d} | lambda {LAMBDA[f'{m}|{g}']:.2f}")

# ---------------- programs: communities among overlap-null calls, labelled by themes enriched among member sets ----------------
called_any = np.zeros(n, bool)
for r in RES.values(): called_any |= r["call"]
idx = np.where(called_any)[0]; Gs = G.subgraph(idx).copy()
comms = sorted(nx.community.greedy_modularity_communities(Gs, weight="weight", resolution=1.2), key=len, reverse=True)
def short(p):
    p = re.sub(r"^(REACTOME|BIOCARTA|KEGG_MEDICUS_(REFERENCE|VARIANT|PATHOGEN|ENV_FACTOR)?)_", "", p); p = p.replace("_PATHWAY", "").replace("_", " ").lower()
    return p[:1].upper() + p[1:]
prog_of = -np.ones(n, int); PROG = []
for k, c in enumerate(comms):
    c = sorted(c); prog_of[c] = k; sub = Gs.subgraph(c); dg = dict(sub.degree(weight="weight")); hub = max(c, key=lambda i: dg[i])
    PROG.append(dict(program=k, label=label_program_enriched([pw[i] for i in c], pw), hub=pw[hub], size=len(c), central_set_label=short(pw[hub])[:42]))
log("programs", [(p["label"], p["size"]) for p in PROG[:14]])

# ---------------- agreement at t = 1 under realistic overlap noise already calibrated; reuse v1 sweep ----------------
rows = []
for (m, g), r in RES.items():
    for i in np.where(r["call"] | (FDR[(m, g)] < .05) | S[f"call|{m}|{g}"])[0]:
        rows.append(dict(method=m, panel=g, pathway=pw[i], program=int(prog_of[i]), program_label=PROG[prog_of[i]]["label"] if prog_of[i] >= 0 else "",
                         z=round(float(Z[(m, g)][i]), 3), raw_padj=float(FDR[(m, g)][i]), min_p=float(r["minp"][i]), q_overlap=float(r["q"][i]),
                         p_fwer_overlap=float(r["p_fwer"][i]), call_overlap=bool(r["call"][i]), call_label_null=bool(S[f"call|{m}|{g}"][i]),
                         sign=int(r["sign"][i]), best_t=float(r["best_t"][i]), persistence_t=";".join(str(t) for t, f in zip(TG, r["persist"][i]) if f)))
pd.DataFrame(rows).to_csv(f"{OUT}/diffusion_calls_v2_{CLS}_{MODEL}.csv", index=False)
summ = pd.DataFrame([dict(method=m, panel=g, raw_fdr05=int((FDR[(m, g)] < .05).sum()), calls_overlap_fdr=int(r["call"].sum()), calls_overlap_fwer=int(r["call_fwer"].sum()),
                          up=int((r["call"] & (r["sign"] > 0)).sum()), down=int((r["call"] & (r["sign"] < 0)).sum()), calls_label_null_fdr=int(S[f"call|{m}|{g}"].sum()),
                          inflation_lambda=round(LAMBDA[f"{m}|{g}"], 3)) for (m, g), r in RES.items()])
summ.to_csv(f"{OUT}/diffusion_summary_v2_{CLS}_{MODEL}.csv", index=False); print(summ.to_string(index=False))
json.dump(dict(calibration_overlap=cal, inflation=LAMBDA, programs=PROG), open(f"{OUT}/diffusion_meta_v2_{CLS}_{MODEL}.json", "w"), indent=1, default=float)
np.savez_compressed(f"{OUT}/state_v2_{CLS}_{MODEL}.npz", prog_of=prog_of, C=C.astype(np.float32),
                    **{f"z|{m}|{g}": z for (m, g), z in Z.items()},
                    **{f"{k}|{m}|{g}": RES[(m, g)][k] for (m, g) in RES for k in ("call", "call_fwer", "sign", "best_t", "persist", "minp", "q")})
log("done")
