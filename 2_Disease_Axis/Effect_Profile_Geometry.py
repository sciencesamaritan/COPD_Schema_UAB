"""Axis-independent geometry of the four schema-group effect profiles.

Places no COPD (0), newly diagnosed (1), potentially excluded (2) and established COPD (3) relative to one
another without privileging the established-COPD axis:
  * pairwise distance between groups a and b = root-mean-square log2 fold change of the a-vs-b contrast
    (no COPD is the reference, so d(g, 0) uses the g-vs-0 contrast);
  * distance ratio = d(g, established) / d(g, no COPD);
  * classical multidimensional scaling of the 4 x 4 distance matrix;
  * Spearman correlation among the three moderated-t profiles versus no COPD.
Computed on the union of transcripts at FDR < 0.05 in any disease-versus-no-COPD contrast (primary), on the
established-COPD signature transcripts, and on all transcripts.

Input: gene-level limma results (all six contrasts) from
       1_Differential_Expression/COPDGene_RNAseq_DE_FGSEA_Workflow.R
Output: effect_profile_distances.csv, effect_profile_mds.csv, effect_profile_spearman.csv
"""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr

#### Edit paths to match your file system ####
DE_FILE = "path/to/gene_de_Model_1_DisCOPDAI_allresults.csv"
OUT = "path/to/disease_axis_outputs/tables"
CLS = "DisCOPDAI"
FDR = 0.05

os.makedirs(OUT, exist_ok=True)
d = pd.read_csv(DE_FILE, usecols=["id", "Contrast", "logFC", "t", "adj.P.Val"])
d["c"] = d.Contrast.str.replace(f"{CLS}_", "", regex=False)
L = d.pivot_table(index="id", columns="c", values="logFC")
T = d.pivot_table(index="id", columns="c", values="t")
Q = d.pivot_table(index="id", columns="c", values="adj.P.Val")

GROUPS = {"0": "No COPD", "1": "Newly diagnosed", "2": "Potentially excluded", "3": "Established COPD"}
def contrast(a, b):
    """Column holding the a-vs-b effect (sign is irrelevant for distances)."""
    return f"{a}vs{b}" if f"{a}vs{b}" in L.columns else f"{b}vs{a}"

masks = {
    "union_any_vs_noCOPD": Q[["1vs0", "2vs0", "3vs0"]].lt(FDR).any(axis=1),
    "established_signature": Q["3vs0"].lt(FDR),
    "all_transcripts": pd.Series(True, index=L.index),
}

dist_rows, mds_rows, sp_rows = [], [], []
for name, m in masks.items():
    rms = lambda c: float(np.sqrt(np.nanmean(L.loc[m, c] ** 2)))
    keys = list(GROUPS)
    D = np.zeros((4, 4))
    for i, a in enumerate(keys):
        for j, b in enumerate(keys):
            if i < j: D[i, j] = D[j, i] = rms(contrast(a, b))
    for g in ["1", "2"]:
        dist_rows.append(dict(transcripts=name, n=int(m.sum()), group=GROUPS[g],
                              d_to_established=round(D[int(g), 3], 4), d_to_noCOPD=round(D[int(g), 0], 4),
                              distance_ratio=round(D[int(g), 3] / D[int(g), 0], 2)))
    # classical MDS
    J = np.eye(4) - np.ones((4, 4)) / 4; Bm = -0.5 * J @ (D ** 2) @ J
    w, v = np.linalg.eigh(Bm); o = np.argsort(w)[::-1][:2]
    X = v[:, o] * np.sqrt(np.clip(w[o], 0, None))
    mds_rows += [dict(transcripts=name, group=GROUPS[k], dim1=round(X[i, 0], 4), dim2=round(X[i, 1], 4)) for i, k in enumerate(keys)]
    for a, b in [("1vs0", "3vs0"), ("2vs0", "3vs0"), ("1vs0", "2vs0")]:
        sp_rows.append(dict(transcripts=name, profile_a=a, profile_b=b, spearman_rho=round(spearmanr(T.loc[m, a], T.loc[m, b])[0], 3)))

pd.DataFrame(dist_rows).to_csv(f"{OUT}/effect_profile_distances.csv", index=False)
pd.DataFrame(mds_rows).to_csv(f"{OUT}/effect_profile_mds.csv", index=False)
pd.DataFrame(sp_rows).to_csv(f"{OUT}/effect_profile_spearman.csv", index=False)
print(pd.DataFrame(dist_rows).to_string(index=False)); print(pd.DataFrame(sp_rows).to_string(index=False))
