# Whole-Blood Transcriptomic Signatures Distinguish Multidimensional COPD Schema Reclassified Groups

### Whole-blood RNA-seq analysis and calibrated pathway heat diffusion in COPDGene

---

## 1. Summary

This repository holds the **R scripts** (v.4.5) and **Python scripts** (v.3.11+) behind the manuscript *"Whole-Blood Transcriptomic Signatures Distinguish Multidimensional COPD Schema Reclassified Groups"*.

The multidimensional COPD diagnostic schema (Bhatt et al., JAMA 2025) moves people in two directions. It brings in participants without airflow obstruction who carry imaging and symptom evidence of disease, and it moves out participants with obstruction who carry none. We asked whether the blood transcriptome moves with them. COPDGene participants fall into four groups:

- **No COPD**
- **Newly diagnosed COPD:** no airflow obstruction, three or more minor criteria
- **Potentially excluded:** airflow obstruction, no minor criterion
- **Established COPD**

The code covers four layers of evidence:

- **Differential expression and enrichment.** Gene-level differential expression (**limma-voom**), preranked pathway enrichment (**fgsea**) and participant-level pathway scores (**singscore**).
- **Disease axis.** Each reclassified group is placed on the established-COPD axis by a scaling coefficient **λ**, corrected for regression dilution and for the error correlation that the shared no-COPD reference induces.
- **Calibration.** Genomic-control inflation (**λ<sub>GC</sub>**) and the **break-even inter-gene correlation** separate methods that hold their size under the null from methods that do not.
- **Heat diffusion on pathway and gene graphs.** Programs of pathways and modules of genes are tested against nulls that keep gene overlap and co-expression, and maps are compared by **Moran spectral randomization**.

No participant-level data are committed. COPDGene data are available to qualified investigators through the COPDGene Study (www.copdgene.org) and dbGaP.

---

## 2. Pipeline Outline

```mermaid
flowchart TD
    D0["COPDGene Freeze 5<br/>whole-blood RNA-seq, Phase 2<br/>n = 3,309"] --> S1
    P0["Phenotypes and<br/>schema classification"] --> S1
    P0 --> S5

    subgraph S1["1_Differential_Expression"]
        direction LR
        A1["filterByExpr, TMM, voom<br/>limma Model 1 and Model 2"] --> A2["fgsea on moderated t"]
        A1 --> A3["singscore per participant,<br/>limma regression"]
    end

    S1 --> S2
    S1 --> S3

    subgraph S2["2_Disease_Axis"]
        direction LR
        B1["λ: through-origin scaling<br/>dilution and shared-reference corrections"] --> B2["Bootstrap CIs,<br/>variance partition"]
        B1 --> B3["Axis-independent geometry:<br/>RMS distances, MDS, Spearman"]
    end

    subgraph S3["3_Pathway_Enrichment"]
        direction LR
        C1["Pathway evidence engine:<br/>redundancy collapse, Leiden communities"] --> C2["Profile clusters,<br/>Model 2 sensitivity"]
    end

    S1 --> S4
    S3 --> S4

    subgraph S4["4_Heat_Diffusion"]
        direction LR
        E1["Pathway graph, heat kernel,<br/>overlap-preserving null"] --> E2["Calibrated program calls,<br/>Moran spectral randomization"]
        E1 --> E3["Gene graph,<br/>co-expression null, modules"]
        E2 --> E4["λGC, break-even ρ*,<br/>robustness"]
    end

    S5["5_Tables<br/>Table 1, Phase 1 to Phase 2 migration"]
    S2 --> S6
    S3 --> S6
    S4 --> S6
    S6["6_Figures<br/>Figures 1 to 6, supplementary figures,<br/>Cytoscape network"]
```

1. **Cohort and classification**
   - COPDGene Freeze 5 whole-blood RNA-seq, Phase 2 (visit 2); de-identified (randomized) IDs.
   - Complete case on the three schema classifications, smoking status and blood-cell differentials: **n = 3,309**.
   - Schema: AI-derived visual CT classification (**primary**); quantitative CT classification (emphysema ≥5%, Pi10 ≥2.24 mm; **sensitivity**).

2. **Differential expression** (`limma` v.3.66.0, `edgeR`)
   - **filterByExpr** (min.count 50, min.prop 0.9): 16,554 transcripts. **TMM** normalization.
   - **voom** + **lmFit** + **eBayes** per model and classification:
     - Model 1: library batch, age, race, sex, smoking status.
     - Model 2: Model 1 + eosinophil, neutrophil, lymphocyte and monocyte %.
   - Contrasts: each disease group vs no COPD (primary); pairwise disease-group contrasts (secondary).
   - Benjamini–Hochberg FDR.

3. **Pathway enrichment and participant-level scores**
   - **fgseaMultilevel** (v.1.36.2) on moderated t against MSigDB collections (`msigdbr` v.25.1.1; Hallmark, Reactome, BioCarta, KEGG_MEDICUS, ImmuneSigDB; set size 15–500).
   - **singscore** (undirected) on batch-corrected log-CPM, regressed on the same contrasts.
   - Pathway evidence engine: redundancy collapse, Leiden communities, profile clustering (Ward.D2, silhouette).
   - HGNC symbol validation (`biobouncer` v.0.1.4).

4. **Disease-axis scaling**
   - λ is the through-origin regression of a group's log2 fold changes on the established-COPD log2 fold changes over the 7,515 established-COPD transcripts (FDR < 0.05).
   - Standard errors recovered from limma confidence intervals.
   - Corrections for regression dilution and for the shared-reference covariance.
   - 2,000 transcript bootstraps; paired bootstrap of the difference; ρ sensitivity; on- and off-axis variance partition; per-gene Wald departures.
   - Axis-independent geometry: pairwise RMS distances, classical MDS, Spearman correlation of profiles.

5. **Calibrated heat diffusion** (Python: `numpy`, `scipy`, `pandas`, `networkx`, `scikit-learn`, `matplotlib`)
   - **Pathway graph.** 1,232 Reactome, BioCarta and KEGG pathways scored by both methods; Jaccard k-nearest-neighbor graph (k = 10); normalized-Laplacian heat kernel; scale cap at an effective neighborhood of 20 sets.
   - **Overlap-preserving null.** N(0, C), with C the cosine gene overlap. Multiscale minimum-P calibrated with 2,000 null draws; BH-FDR and FWER.
   - **Agreement between contrasts.** Correlation of diffused maps, tested by Moran spectral randomization; shared-reference null.
   - **Gene propagation.** 8,046 genes on a shared-membership graph; co-expression-preserving graph-Fourier null estimated from the potentially excluded contrast; planted-module power; module annotation.
   - **Calibration.** λ<sub>GC</sub> = median(z²)/0.455; break-even inter-gene correlation ρ* = ((z/1.96)² − 1)/(m − 1).
   - **Robustness.** k = 5/10/20, Jaccard vs overlap coefficient.

6. **Tables**
   - Table 1 (Phase 2 characteristics) and the Phase 1 → Phase 2 migration table.

7. **Figures**
   - Figures 1–6 and supplementary figures.
   - Cytoscape pathway-overlap network.

---

## 3. Script Descriptions

### `1_Differential_Expression/`
- **`COPDGene_RNAseq_DE_FGSEA_Workflow.R`**: Main pipeline. Steps:
  - loading, cohort assembly, filtering, TMM and QC (library size, PCA, voom mean–variance);
  - limma-voom fits for both models and the three classifications, all contrasts;
  - fgsea.

  Writes `gene_de_<Model>_<Classification>_allresults_<date>.csv` and `gsea_..._allresults_<date>.csv`.
- **`COPDGene_singscore_Workflow.R`**: Participant-level singscore pathway scores and their limma regression. Writes `singscore_..._allresults_<date>.csv`.

### `2_Disease_Axis/`
- **`config.yml`**: Every path and parameter for this folder (group sizes, FDR, bootstrap replicates, seed).
- **`00_config.R`**: Path resolution, run log, check recorder.
- **`01_load.R`**: Loads gene-level results; recovers standard errors from confidence intervals and verifies them against the moderated t.
- **`02_scaling.R`**: The λ estimator:
  - design-implied shared-reference correlation;
  - regression-dilution correction;
  - on- and off-axis variance partition;
  - paired bootstrap, ρ sensitivity, per-gene departures and their detection power.
- **`03_enrichment.R`**: Builds the MSigDB cache (`msigdb_cache.rds`) and recomputes fgsea from the moderated t.
- **`05_reproduce.R`** / **`run_reproduction.R`**: Runs the scaling analysis and the fgsea recompute and writes tables and an analysis log. Checks against reference tables when they are present.
- **`Effect_Profile_Geometry.py`**: Pairwise RMS distances, distance ratios, classical MDS and Spearman profile correlations.
- **`Fig2_Disease_Axis_Geometry.py`**: Figure 2 (hexbin effect plots with λ and 95% CI).

### `3_Pathway_Enrichment/`
- **`Pathway_Evidence_Engine.R`**: Integrates fgsea, singscore and gene-level evidence; collapses redundant pathways; builds Leiden communities of overlap; produces the competitive-test concordance audit.
- **`Pathway_Profile_Clusters.R`**: Clusters the 263 non-redundant pathways by their standardized three-contrast NES profiles (Ward.D2, average silhouette); extracts leading-edge genes.
- **`Pathway_Enrichment_singscore_Figure.R`**: FGSEA representatives and the singscore matrix figure.
- **`Pathway_Robustness_Audit.R`**: Graph-threshold stability and rank-test concordance figure.
- **`Model2_Sensitivity.R`**: Model 1 vs Model 2 effect concordance and the λ difference.
- **`Hallmark_singscore_Extension.R`**: Optional Hallmark singscore extension, run from de-identified counts.
- **`GOATEA_Input_Preparation.R`**: Exports ranked gene lists for interactive enrichment review.
- **`HGNC_Nomenclature.R`**, **`Nomenclature_Audit.R`**: HGNC symbol validation helpers and the audit.

### `4_Heat_Diffusion/`
- **`paths.py`**: Input and output paths for every diffusion script.
- **`heat_kernel.py`**: Shared functions:
  - kNN gene-set graph, heat kernel, scale cap;
  - multiscale calls, Moran's I;
  - Moran spectral randomization agreement.
- **`calibrated_null.py`**: Shared functions:
  - overlap correlation C and its null factor;
  - overlap-preserving multiscale minimum-P calls with FDR/FWER;
  - λ<sub>GC</sub>.
- **`rds_reader.py`**: Minimal reader for the R `msigdb_cache.rds`.
- **`pathway_themes.py`**: Theme lexicon and enrichment-based labels for programs.
- **`01_Pathway_Graph_and_Diffusion.py`**: Builds the pathway graph and layout; runs a label-permutation reference analysis; produces scale sweeps and map agreement.
- **`02_Pathway_Calibrated_Programs.py`**: Overlap-preserving null; calibration on simulated overlap noise; calibrated program calls; program communities and labels.
- **`03_Gene_Network_Propagation.py`**: Gene graph; co-expression-preserving null; gene calls; planted-module power; modules; gene–pathway bridge.
- **`04_Diffusion_Robustness.py`**: Graph-construction sensitivity, shared-reference null and null smoothness.
- **`05_Calibration_BreakEven_Correlation.py`**: Break-even inter-gene correlation for every tested pathway.
- **`06_Figures_Pathway_Diffusion.py`**, **`07_Figures_Gene_Propagation.py`**, **`08_Figures_Publication.py`**: Diffusion figures (Figures 4–6, Supplementary Figures 6–9) and their publication layout.

### `5_Tables/`
- **`copdgene_phenotypes.py`**, **`multidimensional_schema.py`**, **`Table1_Phase2_Characteristics.py`**: The multidimensional schema computed from phenotype variables. The scripts:
  - validate against the published Phase 1 counts;
  - classify participants at Phase 2;
  - tabulate Phase 1 → Phase 2 migration;
  - build Table 1.

  The analysis scripts write the other manuscript tables: the disease-axis tables come from `2_Disease_Axis/run_reproduction.R`, and the calibration, program, break-even, robustness and gene-module tables from `4_Heat_Diffusion/02_` to `05_`.

### `6_Figures/`
- **`nature_style.py`**: Shared figure style (fonts, sizes, palette), imported by the Python figure scripts here and in `4_Heat_Diffusion/`.
- **`Fig1_Schema_Cohort.py`**: Figure 1 (schema cross-classification and cohort flow).
- **`Fig2_Axis_Volcano_Combined.py`**: Combined axis-geometry and volcano figure, assembled from `2_Disease_Axis/Fig2_Disease_Axis_Geometry.py` and `Fig3_Volcano_Plots.py` without refitting.
- **`Fig3_Volcano_Plots.R`**, **`Fig3_Volcano_Plots.py`**: Figure 3 (three-contrast volcano plots).
- **`Pathway_Cluster_Figure.py`**: Pathway-profile cluster figure.
- **`Cytoscape_Network_Export.py`**, **`Cytoscape_Render.py`**: Pathway-overlap network export and CyREST rendering.

---

## 4. Usage

1. **Obtain data access.** Request COPDGene phenotype and Freeze 5 RNA-seq data through the COPDGene Study and dbGaP, along with the de-identified schema classification files.
2. **Install packages.**
   - **R:** `data.table`, `edgeR`, `limma`, `fgsea`, `msigdbr`, `singscore`, `GSEABase`, `AnnotationHub`, `ensembldb`, `igraph`, `ComplexHeatmap`, `circlize`, `ggplot2`, `ggrepel`, `patchwork`, `hexbin`, `matrixStats`, `org.Hs.eg.db`, `yaml`, `jsonlite`, `dplyr`, `DT`, `biobouncer`.
   - **Python:** `numpy`, `scipy`, `pandas`, `networkx`, `scikit-learn`, `matplotlib`.
3. **Edit paths.** Each script opens with a path block: `DATA_DIR`/`RES_DIR` in the R workflows, `2_Disease_Axis/config.yml`, `4_Heat_Diffusion/paths.py`, and the `path/to/` entries elsewhere. Point them at your file system.
   The path blocks keep the file names of the original analysis. They correspond to these outputs:

   | File named in path blocks | Produced by |
   |---|---|
   | `gene_de_<Model>_<Classification>_allresults_*.csv` | `1_Differential_Expression/COPDGene_RNAseq_DE_FGSEA_Workflow.R` |
   | `gsea_<Model>_<Classification>_allresults_*.csv`; `fgsea_differential_all.csv` (all models, classifications and contrasts stacked, with `Model` and `Classification` columns) | `1_Differential_Expression/COPDGene_RNAseq_DE_FGSEA_Workflow.R` |
   | `singscore_<Model>_<Classification>_allresults_*.csv` | `1_Differential_Expression/COPDGene_singscore_Workflow.R` |
   | `fgsea_recomputed_<Model>_<Classification>.csv`, `msigdb_cache.rds`, `bootstrap_lambda.csv` | `2_Disease_Axis/run_reproduction.R` |
   | `pathway_integrated_evidence_model1.csv`, `igraph_pathway_clusters_model1.csv` | `3_Pathway_Enrichment/Pathway_Evidence_Engine.R` |
   | `partition_groups_*.csv` (counts of airflow obstruction by schema group) | cross-tabulation of the de-identified classification file with post-bronchodilator spirometry |

4. **Run in order:**
   ```
   Rscript 1_Differential_Expression/COPDGene_RNAseq_DE_FGSEA_Workflow.R
   Rscript 1_Differential_Expression/COPDGene_singscore_Workflow.R
   Rscript 2_Disease_Axis/run_reproduction.R
   python  2_Disease_Axis/Effect_Profile_Geometry.py
   Rscript 3_Pathway_Enrichment/Pathway_Evidence_Engine.R
   Rscript 3_Pathway_Enrichment/Pathway_Profile_Clusters.R
   Rscript 3_Pathway_Enrichment/Model2_Sensitivity.R
   python  4_Heat_Diffusion/01_Pathway_Graph_and_Diffusion.py      # ~5 min
   python  4_Heat_Diffusion/02_Pathway_Calibrated_Programs.py
   python  4_Heat_Diffusion/03_Gene_Network_Propagation.py
   python  4_Heat_Diffusion/04_Diffusion_Robustness.py
   python  4_Heat_Diffusion/05_Calibration_BreakEven_Correlation.py
   python  4_Heat_Diffusion/06_Figures_Pathway_Diffusion.py
   python  4_Heat_Diffusion/07_Figures_Gene_Propagation.py
   python  5_Tables/Table1_Phase2_Characteristics.py
   ```
   The figure scripts in `2_`, `3_`, `4_` and `6_Figures/` run after the analyses they plot. Seeds are fixed in each script.

---

## 5. References & Tools

- Bhatt SP, et al. Multidimensional COPD diagnostic schema. *JAMA* 2025. doi:10.1001/jama.2025.7358
- Regan EA, et al. Genetic epidemiology of COPD (COPDGene) study design. *COPD* 2010.
- Ritchie ME, et al. limma. *Nucleic Acids Res* 2015. Law CW, et al. voom. *Genome Biol* 2014.
- Korotkevich G, et al. Fast gene set enrichment analysis (fgsea). *bioRxiv*.
- Foroutan M, et al. singscore. *BMC Bioinformatics* 2018.
- Liberzon A, et al. MSigDB Hallmark gene sets. *Cell Syst* 2015.
- Wu D, Smyth GK. camera: competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res* 2012.
- Devlin B, Roeder K. Genomic control for association studies. *Biometrics* 1999.
- Kondor RI, Lafferty J. Diffusion kernels on graphs. *ICML* 2002. Cowen L, et al. Network propagation. *Nat Rev Genet* 2017.
- Wagner HH, Dray S. Moran spectral randomization. *Methods Ecol Evol* 2015.

---

## 6. Authors & Contact

- **Lead author:** Stéphanie Marie Aguilera (Center for Lung Analytics and Imaging Research, University of Alabama at Birmingham)
- **Principal Investigators:**
  - Surya P. Bhatt, University of Alabama at Birmingham: sbhatt@uabmc.edu
  - Peter J. Castaldi, Channing Division of Network Medicine, Brigham and Women's Hospital: peter.castaldi@channing.harvard.edu
- **Co-authors:** Jose S. Cabrera, Adel Boueiz, Craig P. Hersh, Arie Nakhmani, Sandeep Bodduluri, Richard Casaburi, James D. Crapo, Edwin K. Silverman, for the COPDGene Investigators.
