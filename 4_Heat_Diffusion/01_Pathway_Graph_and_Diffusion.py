"""Recipe 1 (plot-atlas): pathway heat-diffusion on the COPD-schema enrichment results.
Primary: DisCOPDAI, Model_1. Universe: pathways tested by FGSEA and singscore in every panel
(BioCarta, KEGG_MEDICUS, Reactome; Hallmark absent from singscore; ImmuneSigDB excluded by default)."""
import sys, os, json, time, re
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd, networkx as nx
from scipy.stats import norm
import heat_kernel as core
from rds_reader import read_rds, names

from paths import DE_FILE, SINGSCORE_FILE, FGSEA_FILE, MSIGDB_CACHE, OUT, FIG, CLS, MODEL
INCLUDE_IMMUNESIG = os.environ.get("IMMUNESIG", "0") == "1"
B = int(os.environ.get("B", "2000"))
rng = np.random.default_rng(20261008)
t0 = time.time(); log = lambda *a: print(f"[{time.time()-t0:6.1f}s]", *a, flush=True)

# ---------------- inputs ----------------
fg = pd.read_csv(FGSEA_FILE, usecols=["pathway", "source", "padj", "NES", "size", "Contrast", "leadingEdge"])
ss = pd.read_csv(SINGSCORE_FILE)
de = pd.read_csv(DE_FILE, usecols=["gene_label", "Contrast", "t"])
cache = read_rds(MSIGDB_CACHE)
sets_obj = cache.value[names(cache).index("sets")]
GS = {n: v.value for n, v in zip(names(sets_obj), sets_obj.value)}
measured = set(de["gene_label"].dropna().astype(str))
log("sets", len(GS), "measured genes", len(measured))

CT = {"New": f"{CLS}_1vs0", "Excluded": f"{CLS}_2vs0", "Existing": f"{CLS}_3vs0", "New vs Existing": f"{CLS}_1vs3"}
VS0 = ["New", "Excluded", "Existing"]
srcs = {"BioCarta", "KEGG_MEDICUS", "Reactome"} | ({"ImmuneSigDB"} if INCLUDE_IMMUNESIG else set())
fgu = fg[fg.source.isin(srcs) & fg.Contrast.isin(CT.values())]
ssu = ss[ss.Contrast.isin(CT.values())]
common = set(GS)
for ct in CT.values():
    common &= set(fgu[fgu.Contrast == ct].pathway) & set(ssu[ssu.Contrast == ct].pathway)
pw = sorted(p for p in common if len(set(GS[p]) & measured) >= 10)
n = len(pw); log("universe", n)
src = fg.drop_duplicates("pathway").set_index("pathway").loc[pw, "source"].values
gs = {p: sorted(set(GS[p]) & measured) for p in pw}

# ---------------- signals ----------------
def piv(df, val):
    return df.pivot_table(index="pathway", columns="Contrast", values=val, aggfunc="first").loc[pw]
nes, fdr_fg = piv(fgu, "NES"), piv(fgu, "padj")
ssz = ssu.assign(z=np.sign(ssu.logFC) * norm.isf(ssu["P.Value"].clip(1e-300) / 2))
zss, fdr_ss = piv(ssz, "z"), piv(ssu, "adj.P.Val")
sd_fg = np.nanstd(nes[[CT[g] for g in VS0]].values, ddof=1); sd_ss = np.nanstd(zss[[CT[g] for g in VS0]].values, ddof=1)
Z = {("FGSEA", g): nes[CT[g]].values / sd_fg for g in CT} | {("singscore", g): zss[CT[g]].values / sd_ss for g in CT}
FDR = {("FGSEA", g): fdr_fg[CT[g]].values for g in CT} | {("singscore", g): fdr_ss[CT[g]].values for g in CT}

# ---------------- graph + checks ----------------
W, L, lam, U = core.graph(pw, gs, k=10)
ncomp = int(np.sum(lam < 1e-8)); isolated = int(np.sum((W > 1e-3).sum(1) == 0))
G = nx.from_numpy_array(W); log("graph components", nx.number_connected_components(G), "zero eig", ncomp, "isolated", isolated)
cap = core.scale_cap(lam, U); TG = np.round(np.arange(0, cap + 1e-9, 0.25), 2); log("scale cap", cap, "grid", TG)
from scipy.stats import spearmanr
A3 = {g: float(spearmanr(Z[("FGSEA", g)], Z[("singscore", g)]).statistic) for g in CT}
MOR = {f"{m}|{g}": core.moran(Z[(m, g)], W, B=999, rng=np.random.default_rng(3)) for m in ("FGSEA", "singscore") for g in CT}
log("moran", {k: (round(v[0], 3), round(v[1], 4)) for k, v in MOR.items()})

# ---------------- multiscale calls ----------------
RES = {}
for (m, g), z in Z.items():
    r = core.multiscale(z, lam, U, TG, B=B, q=0.05, rng=np.random.default_rng(7))
    RES[(m, g)] = r; log(f"calls {m:9s} {g:16s} FDR {int(r['call'].sum()):4d}  FWER {int(r['call_fwer'].sum()):4d}  raw {int((FDR[(m,g)]<.05).sum()):4d}")

# ---------------- agreement + t-sweep ----------------
PAIRS = [("FGSEA", "New", "Existing"), ("FGSEA", "Excluded", "Existing"), ("FGSEA", "New", "Excluded"),
         ("singscore", "New", "Existing"), ("singscore", "Excluded", "Existing"), ("singscore", "New", "Excluded")]
METHOD_PAIRS = [("FGSEA vs singscore", g) for g in VS0]
AGR = []; SWEEP = []
tsw = np.round(np.arange(0, 4.0001, 0.25), 2)
for m, a, b in PAIRS:
    for t in tsw:
        r = core.agreement(Z[(m, a)], Z[(m, b)], lam, U, t, B=1000, rng=np.random.default_rng(11))
        SWEEP.append(dict(method=m, a=a, b=b, t=t, r=r["r_t"], msr_lo=r["band_msr"][0], msr_hi=r["band_msr"][1], p_msr=r["p_msr"],
                          lab_lo=r["band_label"][0], lab_hi=r["band_label"][1], p_label=r["p_label"], r_raw=r["r_raw"]))
        if t == 1.0: AGR.append(dict(kind="between groups", method=m, a=a, b=b, **{k: (v if not isinstance(v, np.ndarray) else list(np.round(v, 4))) for k, v in r.items()}))
for lab, g in METHOD_PAIRS:
    r = core.agreement(Z[("FGSEA", g)], Z[("singscore", g)], lam, U, 1.0, B=1000, rng=np.random.default_rng(12))
    AGR.append(dict(kind="between methods", method=lab, a=g, b=g, **{k: (v if not isinstance(v, np.ndarray) else list(np.round(v, 4))) for k, v in r.items()}))
log("agreement done")

# ---------------- calibration on the real graph ----------------
cal = {}
H1 = core.heat(lam, U, 1.0); rr = np.random.default_rng(99); rej_m = rej_l = 0; NP = 300
for _ in range(NP):
    a, b = H1 @ rr.normal(size=n), H1 @ rr.normal(size=n)
    q = core.agreement(a, b, lam, U, 1.0, B=200, rng=rr); rej_m += q["p_msr"] < .05; rej_l += q["p_label"] < .05
cal["agreement_typeI_spectral"] = rej_m / NP; cal["agreement_typeI_label"] = rej_l / NP
fw = []
for k in range(60):
    z0 = rr.normal(size=n) if k % 2 == 0 else rr.standard_t(3, size=n)
    fw.append(core.multiscale(z0, lam, U, TG, B=400, rng=rr)["call_fwer"].any())
cal["fwer_any_call_rate_noise"] = float(np.mean(fw)); cal["fwer_reps"] = len(fw)
# FDR: plant one program = a connected neighbourhood of ~20 sets
fdp = []
for k in range(30):
    seed = rr.integers(n); prog = list(nx.single_source_shortest_path_length(G, seed, cutoff=2))[:20]
    z0 = rr.normal(size=n); z0[prog] += 1.0
    c = core.multiscale(z0, lam, U, TG, B=400, rng=rr)["call"]; called = np.where(c)[0]
    fdp.append(0 if len(called) == 0 else np.mean(~np.isin(called, prog)))
cal["fdr_false_call_proportion"] = float(np.mean(fdp)); cal["fdr_reps"] = len(fdp)
log("calibration", cal)

# ---------------- programs (communities among called sets) ----------------
called_any = np.zeros(n, bool)
for r in RES.values(): called_any |= r["call"]
idx = np.where(called_any)[0]
Gs = G.subgraph(idx).copy()
comms = sorted(nx.community.greedy_modularity_communities(Gs, weight="weight", resolution=1.0), key=len, reverse=True) if len(idx) else []
STOP = set("REACTOME BIOCARTA KEGG MEDICUS PATHWAY PATHWAYS OF THE AND IN BY TO SIGNALING SIGNALLING VIA REGULATION ACTIVATION FROM WITH ON FOR A RESPONSE".split())
def label_comm(c):
    from collections import Counter
    cnt = Counter(w for i in c for w in re.split(r"[_\s]+", pw[i]) if w and w not in STOP and not w.isdigit() and len(w) > 2)
    return " / ".join(w.title() for w, _ in cnt.most_common(2))
prog_of = -np.ones(n, int); PROG = []
for k, c in enumerate(comms):
    c = sorted(c); prog_of[c] = k; PROG.append(dict(program=k, label=label_comm(c), size=len(c)))
log("programs", [(p["label"], p["size"]) for p in PROG[:12]])

# ---------------- layout ----------------
pos = nx.kamada_kawai_layout(G, weight=None) if n <= 2500 else nx.spectral_layout(G)
pos = nx.spring_layout(G, pos=pos, weight="weight", seed=1, iterations=50)
P = np.array([pos[i] for i in range(n)]); P = (P - P.min(0)) / (P.max(0) - P.min(0)) * 0.88 + 0.06
log("layout done")

# ---------------- tables ----------------
rows = []
for (m, g), r in RES.items():
    for i in np.where(r["call"] | (FDR[(m, g)] < .05))[0]:
        rows.append(dict(method=m, panel=g, pathway=pw[i], source=src[i], program=int(prog_of[i]),
                         program_label=PROG[prog_of[i]]["label"] if prog_of[i] >= 0 else "", z=round(float(Z[(m, g)][i]), 4),
                         raw_padj=float(FDR[(m, g)][i]), min_p=float(r["minp"][i]), q_multiscale=float(r["q"][i]), p_fwer=float(r["p_fwer"][i]),
                         call=bool(r["call"][i]), call_fwer=bool(r["call_fwer"][i]), sign=int(r["sign"][i]), best_t=float(r["best_t"][i]),
                         persistence_t=";".join(str(t) for t, f in zip(TG, r["persist"][i]) if f)))
calls = pd.DataFrame(rows); calls.to_csv(f"{OUT}/diffusion_calls_{CLS}_{MODEL}.csv", index=False)
pd.DataFrame(AGR).to_csv(f"{OUT}/diffusion_agreement_t1_{CLS}_{MODEL}.csv", index=False)
pd.DataFrame(SWEEP).to_csv(f"{OUT}/diffusion_tsweep_{CLS}_{MODEL}.csv", index=False)
summ = pd.DataFrame([dict(method=m, panel=g, contrast=CT[g], raw_fdr05=int((FDR[(m, g)] < .05).sum()), program_calls_fdr=int(r["call"].sum()),
                          program_calls_fwer=int(r["call_fwer"].sum()), up=int((r["call"] & (r["sign"] > 0)).sum()), down=int((r["call"] & (r["sign"] < 0)).sum()),
                          moran_I=round(MOR[f"{m}|{g}"][0], 4), moran_p=MOR[f"{m}|{g}"][1]) for (m, g), r in RES.items()])
summ.to_csv(f"{OUT}/diffusion_summary_{CLS}_{MODEL}.csv", index=False)
meta = dict(classification=CLS, model=MODEL, n=n, k=10, cap=cap, TG=list(map(float, TG)), B=B, components=ncomp, isolated=isolated,
            sources=sorted(srcs), sd_fgsea_NES=sd_fg, sd_singscore_z=sd_ss, A3_spearman=A3, moran={k: list(v) for k, v in MOR.items()}, calibration=cal,
            programs=PROG)
json.dump(meta, open(f"{OUT}/diffusion_meta_{CLS}_{MODEL}.json", "w"), indent=1, default=float)
np.savez_compressed(f"{OUT}/state_{CLS}_{MODEL}.npz", P=P, W=W, lam=lam, U=U, TG=TG, prog_of=prog_of, src=src, pw=np.array(pw),
                    **{f"z|{m}|{g}": z for (m, g), z in Z.items()}, **{f"fdr|{m}|{g}": f for (m, g), f in FDR.items()},
                    **{f"{k}|{m}|{g}": RES[(m, g)][k] for (m, g) in RES for k in ("call", "call_fwer", "sign", "best_t", "persist", "minp", "q")})
print(summ.to_string(index=False)); log("done")
