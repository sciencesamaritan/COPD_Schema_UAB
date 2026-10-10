"""Overlap-preserving null for diffusion program calls (extends plot-atlas Recipe 1 core.py).

Problem: pathway statistics are correlated under the null because overlapping sets share genes,
so the score map is already smooth on the Jaccard graph with no biology. Permuting pathway labels
destroys that smoothness and makes diffused neighbourhood scores look extreme (anti-conservative).

Fix: model null scores as z ~ N(0, C), C_ij = |Gi n Gj| / sqrt(|Gi| |Gj|) (cosine overlap; exact for
mean-type statistics of independent genes). Then (H_t z)_i ~ N(0, (H_t C H_t)_ii) exactly, so each scale
has an analytic per-pathway p, and simulated draws from N(0, C) calibrate the min-p over scales,
BH-FDR across pathways, and FWER. Inputs z must be on a null-unit scale (p-derived z per method)."""
import numpy as np, scipy.sparse as sp
from scipy.stats import norm
from heat_kernel import heat

def overlap_corr(pw, gs):
    genes = sorted({g for p in pw for g in gs[p]}); gi = {g: i for i, g in enumerate(genes)}; n = len(pw)
    rows = [i for i, p in enumerate(pw) for _ in set(gs[p])]; cols = [gi[g] for p in pw for g in set(gs[p])]
    A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, len(genes)))
    I = (A @ A.T).toarray(); s = np.diag(I).copy(); return I / np.sqrt(s[:, None] * s[None, :])

def null_factor(C, eps=1e-9):
    w, V = np.linalg.eigh(C); w = np.clip(w, 0, None); return V * np.sqrt(w)      # C = F F'

def draws(F, B, rng): return F @ rng.normal(size=(F.shape[1], B))

def _bh(p):
    n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]; return np.minimum(q, 1)

def multiscale_overlap(z, lam, U, F, tgrid, B=2000, q=0.05, rng=None, Zn=None):
    """Program calls against N(0, C). Returns the same keys as core.multiscale."""
    rng = rng or np.random.default_rng(7); n = len(z)
    Zn = draws(F, B, rng) if Zn is None else Zn
    T = np.empty((n, len(tgrid))); obs = np.empty((n, len(tgrid))); nul_node = np.ones((n, Zn.shape[1]))
    for j, t in enumerate(tgrid):
        H = heat(lam, U, t); HF = H @ F; s = np.sqrt((HF * HF).sum(1))
        To = (H @ z) / s; T[:, j] = To; obs[:, j] = 2 * norm.sf(np.abs(To))
        nul_node = np.minimum(nul_node, 2 * norm.sf(np.abs((H @ Zn) / s[:, None])))
    minp = obs.min(1); bi = obs.argmin(1)
    pool = np.sort(nul_node.ravel()); M = pool.size; p = (np.searchsorted(pool, minp, side="right") + 1) / (M + 1)
    qv = _bh(p); nul_any = nul_node.min(0); p_fwer = np.array([(1 + np.sum(nul_any <= m)) / (len(nul_any) + 1) for m in minp])
    call = qv < q; cut = minp[call].max() if call.any() else 0
    return dict(T=T, p_scale=obs, minp=minp, p=p, q=qv, call=call, p_fwer=p_fwer, call_fwer=p_fwer < q,
                best_t=np.asarray(tgrid)[bi], sign=np.sign(T[np.arange(n), bi]), persist=obs <= cut)

def inflation(z): return float(np.median(z ** 2) / 0.4549364)   # genomic-control style lambda
