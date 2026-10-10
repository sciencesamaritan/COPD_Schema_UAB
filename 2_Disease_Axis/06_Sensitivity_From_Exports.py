"""Sensitivity analyses that need only the exported Model 1 limma table (no participant-level data).

1. Disease-axis lambda on alternative gene sets (axis transcripts, all transcripts, strong/moderate axis strata),
   with the same dilution and shared-reference corrections as 02_scaling.R.
2. Pairwise effect-profile distances with the measurement-noise floor removed:
   d_corrected = sqrt(mean(L^2) - mean(se^2)) over the chosen transcripts.
3. Expected Spearman/Pearson correlation of two contrasts' t-profiles from the shared no-COPD reference alone.
4. Per-transcript equivalence (two one-sided tests) of the potentially excluded effect within +/-25% of the
   established-COPD effect on the axis transcripts, including the shared-reference covariance; BH across transcripts.

Usage: python 06_Sensitivity_From_Exports.py <gene_de_Model_1_DisCOPDAI_allresults.csv> <out_dir>
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm

N = {"0": 1736, "1": 181, "2": 150, "3": 1242}
src, out = Path(sys.argv[1]), Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
d = pd.read_csv(src)
d["se"] = np.where(np.abs(d.t) > 1e-6, d.logFC / d.t, (d["CI.R"] - d["CI.L"]) / (2 * 1.959964))
L = d.pivot_table(index="id", columns="Contrast", values="logFC")
S = d.pivot_table(index="id", columns="Contrast", values="se")
Q = d.pivot_table(index="id", columns="Contrast", values="adj.P.Val")
T = d.pivot_table(index="id", columns="Contrast", values="t")
C = lambda g, h="0": f"DisCOPDAI_{g}vs{h}"


def rho(g, h, ref="0"):
    """Error correlation of contrasts g-vs-ref and h-vs-ref from the shared reference group."""
    return (1 / N[ref]) / np.sqrt((1 / N[g] + 1 / N[ref]) * (1 / N[h] + 1 / N[ref]))


def lam(g, mask):
    x, y, sx, sy = L.loc[mask, C(3)], L.loc[mask, C(g)], S.loc[mask, C(3)], S.loc[mask, C(g)]
    num = np.mean(x * y) - x.mean() * y.mean() - rho(g, "3") * np.mean(sx * sy)
    den = np.mean(x ** 2) - x.mean() ** 2 - np.mean(sx ** 2)
    return num / den


axis = Q[C(3)] < 0.05
strong = Q[C(3)] < 1e-6
sets = {"Axis transcripts (FDR < 0.05, primary)": axis, "All transcripts": pd.Series(True, index=Q.index),
        "Axis, FDR < 1e-6": strong, "Axis, 1e-6 <= FDR < 0.05": axis & ~strong}
lam_rows = [dict(gene_set=k, n=int(m.sum()), lambda_new=lam("1", m), lambda_excluded=lam("2", m),
                 difference=lam("1", m) - lam("2", m)) for k, m in sets.items()]
lam_tab = pd.DataFrame(lam_rows)

union = (Q[[C(1), C(2), C(3)]] < 0.05).any(axis=1)
dsets = {"Union of disease-vs-no-COPD transcripts": union, "Established-COPD axis transcripts": axis,
         "All transcripts": pd.Series(True, index=Q.index)}
dist_rows = []
for k, m in dsets.items():
    for c in [C(1), C(2), C(1, 3), C(2, 3), C(1, 2)]:
        raw = np.sqrt(np.mean(L.loc[m, c] ** 2))
        cor = np.sqrt(max(np.mean(L.loc[m, c] ** 2) - np.mean(S.loc[m, c] ** 2), 0))
        dist_rows.append(dict(gene_set=k, contrast=c.replace("DisCOPDAI_", ""), n=int(m.sum()), rms_raw=raw, rms_noise_corrected=cor))
dist = pd.DataFrame(dist_rows)
ratio = (dist.pivot_table(index="gene_set", columns="contrast", values=["rms_raw", "rms_noise_corrected"]))
ratios = pd.DataFrame({
    "excluded_raw": ratio[("rms_raw", "2vs3")] / ratio[("rms_raw", "2vs0")],
    "excluded_corrected": ratio[("rms_noise_corrected", "2vs3")] / ratio[("rms_noise_corrected", "2vs0")],
    "new_raw": ratio[("rms_raw", "1vs3")] / ratio[("rms_raw", "1vs0")],
    "new_corrected": ratio[("rms_noise_corrected", "1vs3")] / ratio[("rms_noise_corrected", "1vs0")]})

cor_rows = []
for a, b in [("1", "3"), ("2", "3"), ("1", "2")]:
    for k, m in (("Union", union), ("All transcripts", pd.Series(True, index=Q.index))):
        cor_rows.append(dict(pair=f"{a}vs0 with {b}vs0", gene_set=k,
                             spearman_observed=T.loc[m, C(a)].corr(T.loc[m, C(b)], method="spearman"),
                             expected_from_shared_reference=rho(a, b)))
cors = pd.DataFrame(cor_rows)

# Equivalence: |b2| < 0.25 |b3| on axis transcripts, TOST with shared-reference covariance
m = axis
b2, b3, s2, s3 = L.loc[m, C(2)], L.loc[m, C(3)], S.loc[m, C(2)], S.loc[m, C(3)]
r = rho("2", "3"); sg = np.sign(b3); a = 0.25
up = b2 * sg - a * np.abs(b3)                    # want < 0
lo = b2 * sg + a * np.abs(b3)                    # want > 0
v_up = s2 ** 2 + a ** 2 * s3 ** 2 - 2 * a * r * s2 * s3
v_lo = s2 ** 2 + a ** 2 * s3 ** 2 + 2 * a * r * s2 * s3
p_up = norm.cdf(up / np.sqrt(v_up)); p_lo = norm.sf(lo / np.sqrt(v_lo))
p_tost = pd.Series(np.maximum(p_up, p_lo), index=b2.index)
order = np.argsort(p_tost.values); n = len(p_tost)
bh = np.empty(n); bh[order] = np.minimum.accumulate((p_tost.values[order] * n / np.arange(1, n + 1))[::-1])[::-1]
eq = dict(margin_fraction=a, n_axis=int(n), n_equivalent_p05=int((p_tost < 0.05).sum()),
          n_equivalent_fdr05=int((np.minimum(bh, 1) < 0.05).sum()), median_p_tost=float(np.median(p_tost)),
          median_margin_over_se=float(np.median(a * np.abs(b3) / s2)))

lam_tab.to_csv(out / "lambda_gene_set_sensitivity.csv", index=False)
dist.to_csv(out / "distance_noise_corrected.csv", index=False)
ratios.to_csv(out / "distance_ratios.csv")
cors.to_csv(out / "tprofile_correlation_vs_shared_reference.csv", index=False)
json.dump(eq, open(out / "equivalence_tost.json", "w"), indent=1)
print(lam_tab.round(3).to_string(index=False)); print(ratios.round(2)); print(cors.round(3).to_string(index=False)); print(eq)
