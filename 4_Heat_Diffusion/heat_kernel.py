import numpy as np, scipy.sparse as sp

def graph(pw, gs, k=10, sim="jaccard"):
    """pw: list of set names; gs: dict name -> list of genes. Returns W, L, lam, U."""
    genes=sorted({g for p in pw for g in gs[p]}); gi={g:i for i,g in enumerate(genes)}; n=len(pw)
    rows=[i for i,p in enumerate(pw) for _ in set(gs[p])]; cols=[gi[g] for p in pw for g in set(gs[p])]
    A=sp.csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(n,len(genes)))
    I=(A@A.T).toarray(); s=np.diag(I).copy()
    S=I/(s[:,None]+s[None,:]-I) if sim=="jaccard" else I/np.minimum(s[:,None],s[None,:])
    np.fill_diagonal(S,0)
    nb=np.argsort(-S,1)[:,:k]; ii=np.repeat(np.arange(n),k); jj=nb.ravel(); ww=S[ii,jj].copy(); ww[ww<=0]=1e-3
    W=sp.csr_matrix((ww,(ii,jj)),shape=(n,n)); W=W.maximum(W.T).toarray()
    Dh=1/np.sqrt(W.sum(1)); L=np.eye(n)-Dh[:,None]*W*Dh[None,:]
    lam,U=np.linalg.eigh(L); return W, L, np.clip(lam,0,None), U

def heat(lam,U,t): return (U*np.exp(-t*lam))@U.T
def n_eff(lam,U,t): H=heat(lam,U,t); return (H.sum(1)**2)/(H*H).sum(1)
def scale_cap(lam,U,target=20,step=0.25,tmax=6):
    for t in np.arange(step,tmax+1e-9,step):
        if np.median(n_eff(lam,U,t))>=target: return float(t)
    return float(tmax)

def cor_cols(A,B):
    A=A-A.mean(0); B=B-B.mean(0); return (A*B).sum(0)/np.sqrt((A*A).sum(0)*(B*B).sum(0))

def msr(b,U,B,rng):
    c=U.T@b; S=rng.choice([-1.,1.],size=(len(c),B)); return U@(c[:,None]*S)

def agreement(a,b,lam,U,t,B=1000,rng=None):
    rng=rng or np.random.default_rng(1); H=heat(lam,U,t); ha=H@a; hb=H@b
    r=cor_cols(ha[:,None],hb[:,None])[0]; HA=np.repeat(ha[:,None],B,1)
    nm=cor_cols(HA,H@msr(b,U,B,rng)); nl=cor_cols(HA,H@np.column_stack([rng.permutation(b) for _ in range(B)]))
    p=lambda nul:(1+np.sum(np.abs(nul)>=abs(r)))/(B+1)
    return dict(r_raw=float(np.corrcoef(a,b)[0,1]),r_t=float(r),p_msr=p(nm),p_label=p(nl),
                band_msr=np.quantile(nm,[.025,.975]),band_label=np.quantile(nl,[.025,.975]))

def _scales(z,lam,U,tgrid,B,rng):
    zc=z-z.mean(); sdz=zc.std(ddof=1); n=len(z); P=np.column_stack([rng.permutation(zc) for _ in range(B)])
    obs=np.empty((n,len(tgrid))); T=np.empty((n,len(tgrid))); nul_node=np.ones((n,B)); nul_any=np.ones(B)
    for j,t in enumerate(tgrid):
        H=heat(lam,U,t); nr=np.sqrt((H*H).sum(1)); To=(H@zc)/(sdz*nr); Tn=np.abs((H@P)/(sdz*nr[:,None]))
        pool=np.sort(Tn.ravel()); M=pool.size; pv=lambda x:(M-np.searchsorted(pool,x,side="left")+1)/(M+1)
        T[:,j]=To; obs[:,j]=pv(np.abs(To)); pn=pv(Tn); nul_node=np.minimum(nul_node,pn); nul_any=np.minimum(nul_any,pn.min(0))
    return T,obs,nul_node,nul_any

def multiscale(z,lam,U,tgrid,B=2000,q=0.05,rng=None):
    """Scale-balanced min-p program calls with FDR (primary) and FWER."""
    rng=rng or np.random.default_rng(7); n=len(z); T,obs,nul_node,nul_any=_scales(z,lam,U,tgrid,B,rng)
    minp=obs.min(1); bi=obs.argmin(1)
    pool=np.sort(nul_node.ravel()); M=pool.size; p=(np.searchsorted(pool,minp,side="right")+1)/(M+1)
    o=np.argsort(p); bh=np.empty(n); bh[o]=np.minimum.accumulate((p[o]*n/np.arange(1,n+1))[::-1])[::-1]; bh=np.minimum(bh,1)
    p_fwer=np.array([(1+np.sum(nul_any<=m))/(B+1) for m in minp]); call=bh<q
    cut=minp[call].max() if call.any() else 0
    return dict(T=T,p_scale=obs,minp=minp,q=bh,call=call,p_fwer=p_fwer,call_fwer=p_fwer<q,
                best_t=np.asarray(tgrid)[bi],sign=np.sign(T[np.arange(n),bi]),persist=obs<=cut)

def moran(z,W,B=999,rng=None):
    rng=rng or np.random.default_rng(3); Wr=W/W.sum(1,keepdims=True); zc=z-z.mean()
    I=lambda v:(v@Wr@v)/(v@v); obs=I(zc); nul=np.array([I(rng.permutation(zc)) for _ in range(B)])
    return float(obs),float((1+np.sum(nul>=obs))/(B+1))
