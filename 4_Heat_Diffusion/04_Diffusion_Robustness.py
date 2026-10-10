"""Robustness of the key diffusion claims: graph choices, shared-reference null, singscore null smoothness."""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import heat_kernel as core, calibrated_null as core2
from scipy.stats import norm
from rds_reader import read_rds, names
from paths import DE_FILE, SINGSCORE_FILE, FGSEA_FILE, MSIGDB_CACHE, OUT, FIG, CLS, MODEL
S=np.load(f"{OUT}/state_{CLS}_{MODEL}.npz",allow_pickle=True); V=np.load(f"{OUT}/state_v2_{CLS}_{MODEL}.npz")
pw=list(S["pw"]); n=len(pw)
de=pd.read_csv(DE_FILE,usecols=["gene_label"])
c=read_rds(MSIGDB_CACHE); so=c.value[names(c).index("sets")]
GS={k:v.value for k,v in zip(names(so),so.value)}; meas=set(de.gene_label.dropna().astype(str)); gs={p:sorted(set(GS[p])&meas) for p in pw}
Z={g:V[f"z|singscore|{g}"] for g in ["New","Excluded","Existing","New vs Existing"]}
rng=np.random.default_rng(11); res={}
def r_t(a,b,lam,U,t=1.0):
    f=np.exp(-t*lam); ha=U@(f*(U.T@a)); hb=U@(f*(U.T@b)); return float(np.corrcoef(ha,hb)[0,1])
C=V["C"].astype(float); F=core2.null_factor(C)
for k,sim in [(5,"jaccard"),(10,"jaccard"),(20,"jaccard"),(10,"overlap")]:
    W,L,lam,U=core.graph(pw,gs,k=k,sim=sim); row={}
    for a,b in [("New","Existing"),("Excluded","Existing"),("New","Excluded")]:
        obs=r_t(Z[a],Z[b],lam,U)
        # spectral (MSR) null
        cb=U.T@Z[b]; f=np.exp(-lam); ha=U@(f*(U.T@Z[a]))
        nm=[np.corrcoef(ha,U@(f*cb*rng.choice([-1,1],len(cb))))[0,1] for _ in range(500)]
        row[f"{a}-{b}"]=dict(r=round(obs,3),p_msr=round((1+np.sum(np.abs(nm)>=abs(obs)))/501,4))
    res[f"k{k}_{sim}"]=row; print(k,sim,row,flush=True)
# shared-reference null: pairs from N(0,C) with cross-correlation rho_ref (both contrasts share No COPD)
W,L,lam,U=core.graph(pw,gs,k=10); f=np.exp(-lam)
for rho in [0.2,0.3,0.4]:
    rs=[]
    for _ in range(300):
        e0=core2.draws(F,1,rng)[:,0]; e1=core2.draws(F,1,rng)[:,0]; e2=core2.draws(F,1,rng)[:,0]
        a=np.sqrt(rho)*e0+np.sqrt(1-rho)*e1; b=np.sqrt(rho)*e0+np.sqrt(1-rho)*e2
        rs.append(np.corrcoef(U@(f*(U.T@a)),U@(f*(U.T@b)))[0,1])
    res[f"shared_ref_null_rho{rho}"]=dict(mean=round(float(np.mean(rs)),3),q975=round(float(np.quantile(rs,.975)),3)); print(rho,res[f"shared_ref_null_rho{rho}"],flush=True)
# empirical shared-reference correlation: raw r of Excluded vs Existing and New vs Excluded
res["raw_r"]={k:round(float(np.corrcoef(Z[a],Z[b])[0,1]),3) for k,(a,b) in {"New-Existing":("New","Existing"),"Excluded-Existing":("Excluded","Existing"),"New-Excluded":("New","Excluded")}.items()}
# singscore null smoothness: Excluded Moran's I vs overlap-noise band
mI=core.moran(Z["Excluded"],W,B=199)[0]; band=core2.__dict__.get("x")
nb=[core.moran(core2.draws(F,1,rng)[:,0],W,B=1)[0] for _ in range(40)]
res["moran_singscore_Excluded"]=dict(I=round(mI,3),overlap_noise_mean=round(float(np.mean(nb)),3),overlap_noise_q05=round(float(np.quantile(nb,.05)),3))
# program-call stability across k: Jaccard of New call sets
calls={}
for k in [5,10,20]:
    W,L,lam,U=core.graph(pw,gs,k=k); cap=core.scale_cap(lam,U); TG=np.round(np.arange(0,cap+1e-9,0.25),2)
    calls[k]={g:core2.multiscale_overlap(Z[g],lam,U,F,TG,B=500,rng=np.random.default_rng(7))["call"] for g in ["New","Excluded","Existing"]}
    print("k",k,"cap",cap,{g:int(v.sum()) for g,v in calls[k].items()},flush=True)
jac=lambda a,b:float((a&b).sum()/max(1,(a|b).sum()))
res["calls_by_k"]={str(k):{g:int(v.sum()) for g,v in calls[k].items()} for k in calls}
res["call_jaccard_vs_k10"]={f"k{k}":{g:round(jac(calls[k][g],calls[10][g]),3) for g in ["New","Existing"]} for k in [5,20]}
json.dump(res,open(f"{OUT}/robustness.json","w"),indent=1); print(json.dumps(res,indent=1))
