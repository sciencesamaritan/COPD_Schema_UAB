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

### Gene expression data processing

Whole-blood RNA sequencing data from COPDGene Study participants were obtained from Freeze 5, with raw gene-level counts quantified using Salmon against the GENCODE v37 annotation on the GRCh38.p13 reference genome. Samples were restricted to Visit 2 participants who passed quality control (`final.analyzed = 1`) and had complete data for all three COPD classification variables, smoking status, and peripheral blood cell differential percentages (n = 3,309). Gene-level counts were assembled into a `DGEList` object using edgeR. Lowly expressed genes were filtered once across all analyses using `filterByExpr` with a minimum count threshold of 50 and a minimum sample proportion of 0.9, applied under a covariate-only design so that every downstream model shares the same 16,554 transcripts. Library sizes were normalized using the trimmed mean of M-values (TMM) method via `calcNormFactors`. Log-transformed counts per million (log-CPM) were computed with a prior count of 2, and technical batch effects attributable to library construction batch (`Lc.Batch`) were removed from the log-CPM matrix using `removeBatchEffect` from limma. Gene expression values were subsequently mean-centered per gene across all samples. `Lc.Batch` was excluded from all downstream design matrices to avoid double adjustment.

### COPD classification variables

Three COPD classification schemes were compared, each encoding disease status as a four-level factor (0, 1, 2, 3) with level 0 as the reference: the original multidimensional schema (`DisCOPD`), the AI-derived visual CT classification (`DisCOPDAI`; primary), and the quantitative CT classification (`DisCOPDQuant`; emphysema ≥5%, Pi10 ≥2.24 mm; sensitivity). All three variables were merged onto the phenotype data, and only participants with non-missing values for all three classifications were retained so that every comparison runs on an identical participant set. Groups were assigned at Phase 1, the only visit with visual CT reads; blood was drawn at Phase 2.

### Differential gene expression analysis

Two linear models were specified for each classification variable. Model 1 adjusted for age, race, sex, and smoking status (`SmokCigNow`). Model 2 additionally included peripheral blood cell percentages (eosinophils, neutrophils, lymphocytes, and monocytes). For each model, a design matrix was constructed and precision weights were estimated using the `voom` transformation applied to the uncorrected TMM-normalized `DGEList` with a batch-inclusive design, to model the mean-variance relationship accurately. These weights were then applied during model fitting via `lmFit` on the batch-corrected, mean-centered log-CPM matrix. Six contrasts were tested per classification variable: each disease group versus the control group (1 vs 0, 2 vs 0, 3 vs 0) and all pairwise comparisons among disease groups (1 vs 3, 2 vs 3, 1 vs 2). Contrasts were estimated using `contrasts.fit`, and moderated t-statistics and p-values were obtained via the empirical Bayes framework implemented in `eBayes`. Gene-level results were extracted using `topTable` with Benjamini-Hochberg false discovery rate (FDR) correction.

### Gene set enrichment analysis

For each contrast, genes were ranked by their moderated t-statistics and mapped from ENSEMBL versioned identifiers to Entrez gene IDs. When multiple ENSEMBL identifiers mapped to the same Entrez ID, the gene with the largest absolute t-statistic was retained. Pre-ranked gene set enrichment analysis was performed using `fgseaMultilevel` from the fgsea package with Hallmark and Reactome gene sets obtained from MSigDB via msigdbr. Gene sets were filtered to include between 15 and 500 genes, and exact p-values were computed (`eps = 0`). Pathway-level results were adjusted for multiple testing using the Benjamini-Hochberg method. Gene symbols were validated against HGNC with biobouncer. Preranked enrichment is reported as a descriptive layer; the calibration below gives the reason.

### Single-sample pathway scoring

Single-sample pathway activity scores were computed using the singscore package. Four pathway collections were obtained from MSigDB: Reactome, KEGG (KEGG_MEDICUS), ImmuneSigDB, and BioCarta. The batch-corrected, mean-centered log-CPM expression matrix was mapped from ENSEMBL versioned identifiers to gene symbols using a pre-built gene annotation map; when multiple ENSEMBL identifiers mapped to the same symbol, the gene with the highest average expression was retained. Genes were ranked within each sample using `rankGenes` with ties resolved by the minimum method. Pathway gene sets were filtered to retain those with 15 to 500 genes overlapping the measured transcriptome, then converted to `GeneSet` objects. All pathways were scored simultaneously across all samples using `multiScore` with `knownDirection = FALSE` (MSigDB pathway gene sets do not specify regulatory direction) and `centerScore = TRUE`.

The resulting pathway score matrix (pathways × samples) was treated as a feature matrix analogous to gene expression, with each pathway score representing a single-sample measure of pathway activity. This matrix was then regressed against the same COPD classification contrasts using limma: for each model (Model 1 and Model 2) and each classification variable, design matrices were constructed without `Lc.Batch`, and the standard limma pipeline (`lmFit` → `makeContrasts` → `contrasts.fit` → `eBayes` → `topTable`) was applied. No voom weights were used for pathway score regression, as the scores are continuous and approximately normally distributed. P-values were adjusted using the Benjamini-Hochberg FDR method.

### Pathway evidence engine

Pathway evidence from fgsea, singscore, and gene-level results was integrated for the 263 candidate pathways. Redundant pathways were collapsed, Leiden communities of gene overlap were identified, and pathway profiles were clustered using Ward.D2 on standardized three-contrast NES profiles, with the number of clusters selected by silhouette analysis. Leading-edge genes were extracted for each profile cluster, and graph-threshold stability of the partition was assessed by the adjusted Rand index.

### Cross-classification comparisons

To assess concordance across the three classification systems, pairwise overlaps of FDR-significant genes and pathways (adj. P < 0.05) were computed for each model. Effect-size concordance was evaluated by computing Pearson correlation coefficients of log-fold-change estimates between all pairs of classification variables for each matched contrast under Model 1. Model 1 and Model 2 effects were compared the same way.

### Disease-axis scaling

For the 7,515 transcripts significant in the established-COPD contrast (FDR < 0.05), the scaling coefficient λ was estimated as the through-origin regression of a group's log2 fold changes on the established-COPD log2 fold changes, with λ = 0 the no-COPD reference and λ = 1 the established-COPD profile. Standard errors were recovered from limma confidence intervals and verified against the moderated t-statistics. The denominator was corrected for regression dilution by subtracting the mean error variance, and the numerator for the shared-reference error correlation ρ = (1/n0)/√((1/ng + 1/n0)(1/nh + 1/n0)), which every contrast against the same no-COPD group induces. Confidence intervals were obtained from 2,000 transcript bootstraps, and the two coefficients were compared by a paired bootstrap of 10,000 replicates. Systematic variance was partitioned into on-axis and orthogonal components, per-transcript departures from the scaled expectation were tested by Wald tests with Benjamini-Hochberg control, and a sign-permutation null checked the estimator. Axis-independent geometry used pairwise root-mean-square log2 fold-change distances over the 7,781 transcripts significant in any disease contrast, classical multidimensional scaling, and Spearman correlations among the three moderated-t profiles.

### Calibration

P-values were converted to signed z-scores, z = sign(effect) · Φ⁻¹(1 − P/2), and the genomic-control inflation factor λGC = median(z²)/0.4549 was computed for fgsea, singscore, and limma in every contrast. The potentially excluded contrast, which has no differentially expressed transcripts, served as the null reference. For sets of m measured genes, the competitive set-mean z of limma t under independence was z0 = (t̄set − t̄all)/(st/√m), and the break-even inter-gene correlation ρ* = ((z0/1.96)² − 1)/(m − 1) followed from the variance inflation 1 + (m − 1)ρ.

### Calibrated heat diffusion on the pathway graph

The 1,232 Reactome, BioCarta, and KEGG pathways scored by both methods in all four contrasts, with at least 10 measured genes, were joined to their 10 nearest neighbors by Jaccard similarity of gene membership (symmetrized). With L = I − D^(−1/2) W D^(−1/2) = U diag(λ) Uᵀ, the heat kernel Ht = U diag(e^(−tλ)) Uᵀ was applied to each contrast's signed z-scores, and the scale cap was the first t (step 0.25) at which the median effective neighborhood (Σj Hij)²/Σj Hij² reached 20 sets (t = 2.25). Under no group difference, set statistics of independent genes have correlation Cij = |Gi ∩ Gj|/√(|Gi||Gj|); for z ~ N(0, C) the program statistic Ti(t) = (Ht z)i/√((Ht C Ht)ii) is standard normal at every scale. Each pathway's statistic was the minimum two-sided P over t ∈ {0, 0.25, …, 2.25}, referred to the pooled null distribution of the same minimum from 2,000 draws Fε (C = FFᵀ) through identical kernels, with Benjamini-Hochberg FDR across pathways and family-wise error from the null minimum over pathways and scales. Calibration used 40 simulated null maps and 25 maps with a planted 20-set program. Programs were greedy-modularity communities among called pathways, labeled by themes enriched among member sets relative to the universe.

### Map agreement

Agreement between contrasts was the Pearson correlation of H1 za and H1 zb, tested by Moran spectral randomization: random signs on the graph-Fourier coefficients of zb, rebuilt and diffused (1,000 draws), which preserves each map's power spectrum. The shared-reference null simulated pairs √ρr e0 + √(1 − ρr) e1,2 with e ~ N(0, C) and ρr ∈ {0.2, 0.3, 0.4}. Graph sensitivity used k = 5, 10, 20 and the overlap coefficient in place of Jaccard.

### Gene-level network propagation

Genes in Hallmark, Reactome, BioCarta, or KEGG_MEDICUS sets of 10 to 200 measured genes (8,046 genes; 1,756 sets) were joined by Jaccard similarity of their set-membership vectors (k = 10; scale cap 2.5), with z the moderated t. Under no differential expression z ~ N(0, R), with R the residual expression correlation common to all contrasts; its graph-Fourier power Pk = ukᵀ R uk was estimated from the potentially excluded contrast as (ukᵀ z)² smoothed over 201 adjacent eigenvectors. Null standard deviations were si(t) = (Σk Uik² e^(−2tλk) Pk)^(1/2), and 1,000 null draws U(√P ∘ ε) calibrated multiscale calls as above. An independent-gene null (Pk = 1) is reported for comparison. Power used shifts of 1.0 to 2.0 added to 40 or 100 graph-adjacent genes of the potentially excluded map (10 replicates per setting). Modules were greedy-modularity communities of at least 8 called genes, annotated by hypergeometric overlap with curated sets; layouts used t-SNE of diffusion-map coordinates and serve for display only.

### Tables and figures

Table 1 (Phase 2 characteristics) and the Phase 1 → Phase 2 migration table are built from the phenotype variables by the scripts in `5_Tables/`. Figures 1 to 6, the supplementary figures, and the Cytoscape pathway-overlap network are built by the scripts in `6_Figures/` and `4_Heat_Diffusion/`.

### Software

All R analyses were performed in R (v4.5.1) using the Bioconductor packages edgeR, limma, fgsea, singscore, GSEABase, and msigdbr, with igraph for community detection and biobouncer for HGNC symbol validation. Volcano plots were generated with EnhancedVolcano and ggplot2, and interactive tables were rendered with DT. The disease-axis geometry, heat diffusion, and gene propagation were implemented in Python (v3.11+) with numpy, scipy, pandas, networkx, scikit-learn, and matplotlib. Seeds are fixed in each script.

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

**Cohort and diagnostic schema**
- Bhatt SP, Abadi E, Anzueto A, et al. A multidimensional diagnostic approach for chronic obstructive pulmonary disease. *JAMA*. 2025;333(24):2164–2175. [doi:10.1001/jama.2025.7358](https://doi.org/10.1001/jama.2025.7358)
- Regan EA, Hokanson JE, Murphy JR, et al. Genetic epidemiology of COPD (COPDGene) study design. *COPD*. 2010;7(1):32–43. [doi:10.3109/15412550903499522](https://doi.org/10.3109/15412550903499522)

**Differential expression**
- Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*. 2010;26(1):139–140. [doi:10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616)
- Robinson MD, Oshlack A. A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biol*. 2010;11(3):R25. [doi:10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25)
- Law CW, Chen Y, Shi W, Smyth GK. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biol*. 2014;15(2):R29. [doi:10.1186/gb-2014-15-2-r29](https://doi.org/10.1186/gb-2014-15-2-r29)
- Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res*. 2015;43(7):e47. [doi:10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007)

**Pathway enrichment and participant-level scores**
- Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc Natl Acad Sci USA*. 2005;102(43):15545–15550. [doi:10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)
- Korotkevich G, Sukhov V, Budin N, et al. Fast gene set enrichment analysis. *bioRxiv*. 2021:060012. [doi:10.1101/060012](https://doi.org/10.1101/060012)
- Foroutan M, Bhuva DD, Lyu R, et al. Single sample scoring of molecular phenotypes. *BMC Bioinformatics*. 2018;19(1):404. [doi:10.1186/s12859-018-2435-4](https://doi.org/10.1186/s12859-018-2435-4)
- Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Syst*. 2015;1(6):417–425. [doi:10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004)
- Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. *Sci Rep*. 2019;9(1):5233. [doi:10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z)

**Calibration**
- Devlin B, Roeder K. Genomic control for association studies. *Biometrics*. 1999;55(4):997–1004. [doi:10.1111/j.0006-341x.1999.00997.x](https://doi.org/10.1111/j.0006-341x.1999.00997.x)
- Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res*. 2012;40(17):e133. [doi:10.1093/nar/gks461](https://doi.org/10.1093/nar/gks461)
- Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B*. 1995;57(1):289–300.
- Westfall PH, Young SS. *Resampling-Based Multiple Testing: Examples and Methods for p-Value Adjustment*. Wiley; 1993.

**Network propagation and heat diffusion**
- Kondor RI, Lafferty J. Diffusion kernels on graphs and other discrete input spaces. In: *Proceedings of the 19th International Conference on Machine Learning (ICML)*. 2002:315–322.
- Coifman RR, Lafon S, Lee AB, et al. Geometric diffusions as a tool for harmonic analysis and structure definition of data: diffusion maps. *Proc Natl Acad Sci USA*. 2005;102(21):7426–7431. [doi:10.1073/pnas.0500334102](https://doi.org/10.1073/pnas.0500334102)
- Cowen L, Ideker T, Raphael BJ, Sharan R. Network propagation: a universal amplifier of genetic associations. *Nat Rev Genet*. 2017;18(9):551–562. [doi:10.1038/nrg.2017.38](https://doi.org/10.1038/nrg.2017.38)
- Wagner HH, Dray S. Generating spatially constrained null models for irregularly spaced data using Moran spectral randomization methods. *Methods Ecol Evol*. 2015;6(10):1169–1178.
- Vos de Wael R, Benkarim O, Paquola C, et al. BrainSpace: a toolbox for the analysis of macroscale gradients in neuroimaging and connectomics datasets. *Commun Biol*. 2020;3(1):103. [doi:10.1038/s42003-020-0794-7](https://doi.org/10.1038/s42003-020-0794-7)
- Clauset A, Newman MEJ, Moore C. Finding community structure in very large networks. *Phys Rev E*. 2004;70(6):066111.

**Network visualization**
- Shannon P, Markiel A, Ozier O, et al. Cytoscape: a software environment for integrated models of biomolecular interaction networks. *Genome Res*. 2003;13(11):2498–2504. [doi:10.1101/gr.1239303](https://doi.org/10.1101/gr.1239303)

**Software**
- R (v.4.5): [limma](https://bioconductor.org/packages/limma), [edgeR](https://bioconductor.org/packages/edgeR), [fgsea](https://bioconductor.org/packages/fgsea), [singscore](https://bioconductor.org/packages/singscore), [msigdbr](https://cran.r-project.org/package=msigdbr), [igraph](https://r.igraph.org), [ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap), [ggplot2](https://ggplot2.tidyverse.org), biobouncer (HGNC symbol validation).
- Python (v.3.11+): [NumPy](https://numpy.org), [SciPy](https://scipy.org), [pandas](https://pandas.pydata.org), [NetworkX](https://networkx.org), [scikit-learn](https://scikit-learn.org), [Matplotlib](https://matplotlib.org).
- [Cytoscape](https://cytoscape.org) (CyREST) for the pathway-overlap network.
- Data: [COPDGene](https://www.copdgene.org) and [dbGaP](https://www.ncbi.nlm.nih.gov/gap/); gene sets from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb), [Reactome](https://reactome.org) and [KEGG](https://www.genome.jp/kegg/); symbols from [HGNC](https://www.genenames.org).

---

## 6. Authors & Contact

- **Lead author:** Stéphanie Marie Aguilera (Center for Lung Analytics and Imaging Research, University of Alabama at Birmingham)
- **Principal Investigators:**
  - Surya P. Bhatt, University of Alabama at Birmingham: sbhatt@uabmc.edu
  - Peter J. Castaldi, Channing Division of Network Medicine, Brigham and Women's Hospital: peter.castaldi@channing.harvard.edu
- **Co-authors:** Jose S. Cabrera, Adel Boueiz, Craig P. Hersh, Arie Nakhmani, Sandeep Bodduluri, Richard Casaburi, James D. Crapo, Edwin K. Silverman, for the COPDGene Investigators.
