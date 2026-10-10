"""Edit these paths to match your file system before running the heat-diffusion scripts.

Inputs are summary statistics produced by the R workflows in this repository:
  DE_FILE         1_Differential_Expression/COPDGene_RNAseq_DE_FGSEA_Workflow.R  (gene-level limma results, all contrasts)
  SINGSCORE_FILE  1_Differential_Expression/COPDGene_singscore_Workflow.R      (participant-level pathway-score regressions)
  FGSEA_FILE      2_Disease_Axis/run_reproduction.R                            (fgsea_recomputed.csv)
  MSIGDB_CACHE    2_Disease_Axis/03_enrichment.R :: build_pathways()          (msigdb_cache.rds, msigdbr 25.1.1)
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))

CLS = "DisCOPDAI"      # case definition: AI-derived visual CT (primary); "DisCOPDQuant" for the sensitivity analysis
MODEL = "Model_1"      # Model_1: age, race, sex, smoking status; Model_2 adds leukocyte percentages

DE_FILE = f"path/to/gene_de_{MODEL}_{CLS}_allresults.csv"
SINGSCORE_FILE = f"path/to/singscore_{MODEL}_{CLS}_allresults.csv"
FGSEA_FILE = f"path/to/fgsea_recomputed_{MODEL}_{CLS}.csv"
MSIGDB_CACHE = "path/to/msigdb_cache.rds"

OUT = os.path.join(HERE, "out")          # intermediate states and tables
FIG = os.path.join(HERE, "figures")      # rendered figures
os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)
