"""Calibration diagnostic for preranked pathway enrichment: break-even inter-gene correlation.

Genomic-control inflation (lambda_GC) is computed on the diffusion universes in
02_Pathway_Calibrated_Programs.py (1,232 pathways) and 03_Gene_Network_Propagation.py (8,046 genes).

Break-even inter-gene correlation for each pathway: the competitive set-mean z of limma moderated t under
   gene independence, z0 = (mean t_set - mean t_all) / (sd t_all / sqrt(m)), and the average within-set
   correlation rho* at which |z0| falls to 1.96 given the variance inflation 1 + (m - 1) rho (Wu & Smyth 2012):
   rho* = ((z0 / 1.96)^2 - 1) / (m - 1), floored at 0.

Output (OUT): break_even_intergene_correlation.csv (all MSigDB sets tested by FGSEA with >= 15 measured genes)
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd
from scipy.stats import norm
from rds_reader import read_rds, names
from paths import DE_FILE, FGSEA_FILE, MSIGDB_CACHE, OUT, CLS

CT = {"New": f"{CLS}_1vs0", "Excluded": f"{CLS}_2vs0", "Existing": f"{CLS}_3vs0", "New vs Existing": f"{CLS}_1vs3"}
Z975 = norm.isf(0.025)

de = pd.read_csv(DE_FILE, usecols=["gene_label", "Contrast", "t"])
fg = pd.read_csv(FGSEA_FILE, usecols=["pathway"])

# ---------------- break-even inter-gene correlation per pathway ----------------
cache = read_rds(MSIGDB_CACHE); so = cache.value[names(cache).index("sets")]
GS = {k: set(v.value) for k, v in zip(names(so), so.value)}
de = de[de.gene_label.notna()]
rows = []
for g, ct in CT.items():
    d = de[de.Contrast == ct]
    d = d.loc[d.groupby("gene_label")["t"].apply(lambda x: x.abs().idxmax()).values].set_index("gene_label")["t"]  # collapse duplicate symbols by max |t|
    mu, sd = d.mean(), d.std(ddof=1)
    for p in fg.pathway.unique():
        if p not in GS: continue
        tt = d.reindex(sorted(GS[p] & set(d.index))).dropna(); m = len(tt)
        if m < 15: continue
        z0 = (tt.mean() - mu) / (sd / np.sqrt(m))
        rho = max(0.0, ((z0 / Z975) ** 2 - 1) / (m - 1))
        rows.append(dict(pathway=p, contrast=g, genes=m, mean_t=round(tt.mean(), 3), z_independent=round(z0, 2), break_even_rho=round(rho, 4)))
be = pd.DataFrame(rows); be.to_csv(f"{OUT}/break_even_intergene_correlation.csv", index=False)
print(f"break-even rho written for {be.pathway.nunique()} pathways x {be.contrast.nunique()} contrasts")
