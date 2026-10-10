# ============================================================================
# COPDGene_RNAseq_DE_FGSEA_Workflow.R
#
# Whole-blood RNA-seq (COPDGene Freeze 5, Phase 2) differential expression and
# preranked enrichment for the multidimensional COPD diagnostic schema.
#   1. Load de-identified counts, RNA master file, Phase 2 phenotypes and the
#      three schema classifications (original visual, AI-derived visual CT,
#      quantitative CT); complete-case cohort (n = 3,309)
#   2. filterByExpr (min.count 50, min.prop 0.9), TMM normalization, QC
#   3. limma-voom per Model x Classification:
#        Model_1: batch + age + race + sex + smoking status
#        Model_2: Model_1 + eosinophil, neutrophil, lymphocyte, monocyte %
#      contrasts: each disease group vs no COPD, and pairwise disease groups
#   4. fgsea on moderated t (Hallmark, Reactome, KEGG_MEDICUS, ImmuneSigDB,
#      BioCarta; set size 15-500)
#
# Outputs (RES_DIR): gene_de/*_allresults_<date>.csv, gsea/*_allresults_<date>.csv,
#   gene_de_summary_<date>.csv, gsea_summary_<date>.csv
#
# Code recovered from the executed notebook bhatt_repro_genede_gsea_2026.
# R 4.5; limma 3.66.0; edgeR; fgsea 1.36.2; msigdbr 25.1.1.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma); library(fgsea)
  library(msigdbr); library(dplyr); library(tibble); library(ggplot2); library(DT)
})

#### Edit paths to match your file system ####
DATA_DIR  <- "path/to/COPDGene/deidentified"   # *_randid.csv and *_deid.csv inputs
CACHE_DIR <- "path/to/cache"
RES_DIR   <- "path/to/results"
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RES_DIR,   showWarnings = FALSE, recursive = TRUE)
DATE_STAMP <- format(Sys.Date(), "%Y.%m.%d")

#### Schema classifications (0 = no COPD, 1 = newly diagnosed,
#### 2 = potentially excluded, 3 = established COPD) ####
cls_vars   <- c("DisCOPD", "DisCOPDAI", "DisCOPDQuant")
cls_labels <- c(DisCOPD      = "Original visual classification",
                DisCOPDAI    = "AI-derived visual CT classification (primary)",
                DisCOPDQuant = "Quantitative CT classification")

# ============================================================================
# Source:     bhatt_repro_genede_gsea_2026.Rmd / chunk: preprocess
# Inputs:     counts_raw_freeze5_randid.csv,
#             masterfile_rna_freeze5_randid.csv,
#             COPDGene_P1P2P3_SM_NS_ILDBR_Long_Sep24_randid.csv,
#             {MultidimensionalSchema, NewClassQualitativeAI, NewClassQuantitative}_deid.csv
# Population: visit-2 RNA samples, complete-case on all 3 classifications +
#             CBC + smoking → n = 3309 expected
# Output:     pheno_df, dge (DGEList), gmap (ENSEMBL→{symbol, entrez})
# Methods:    METHODS.md §2-3
# ============================================================================
counts_rds <- file.path(CACHE_DIR, "counts_raw.rds")
if (!file.exists(counts_rds)) {
  cnt <- fread(file.path(DATA_DIR, "counts_raw_freeze5_randid.csv"), header = TRUE)
  gids <- cnt$GENEID; cnt[, GENEID := NULL]
  mat <- as.matrix(cnt); rownames(mat) <- gids
  storage.mode(mat) <- "double"
  saveRDS(mat, counts_rds); rm(cnt); gc()
} else {
  mat <- readRDS(counts_rds)
}
n_counts_cols <- ncol(mat); n_genes_raw <- nrow(mat)

mf <- fread(file.path(DATA_DIR, "masterfile_rna_freeze5_randid.csv"))
n_master_total <- nrow(mf)
mf <- mf[final.analyzed == 1, .(rand_id, Lc.Batch)]
mf[, rand_id := as.character(rand_id)]
n_master_keep <- nrow(mf)

ph <- fread(file.path(DATA_DIR, "COPDGene_P1P2P3_SM_NS_ILDBR_Long_Sep24_randid.csv"),
            select = c("rand_id","visitnum","age_visit","race","gender","SmokCigNow",
                       "eosinphl_pct","neutrophl_pct","lymphcyt_pct","monocyt_pct"))
n_pheno_total <- nrow(ph)
ph <- ph[visitnum == 2]; ph[, rand_id := as.character(rand_id)]
n_pheno_v2 <- nrow(ph)

mds  <- fread(file.path(DATA_DIR, "MultidimensionalSchema_deid.csv"))[, .(rand_id, DisCOPD)]
cqai <- fread(file.path(DATA_DIR, "NewClassQualitativeAI_deid.csv"))[, .(rand_id, DisCOPDAI)]
cqt  <- fread(file.path(DATA_DIR, "NewClassQuantitative_deid.csv"))[, .(rand_id, DisCOPDQuant)]
for (x in list(mds, cqai, cqt)) x[, rand_id := as.character(rand_id)]

pheno <- ph |>
  merge(mf,   by = "rand_id", all.x = FALSE, all.y = FALSE) |>
  merge(mds,  by = "rand_id", all.x = TRUE) |>
  merge(cqai, by = "rand_id", all.x = TRUE) |>
  merge(cqt,  by = "rand_id", all.x = TRUE)
setDT(pheno)
n_after_master_join <- nrow(pheno)
cc_cols <- c("age_visit","race","gender",
             "DisCOPD","DisCOPDAI","DisCOPDQuant",
             "SmokCigNow","eosinphl_pct","neutrophl_pct","lymphcyt_pct","monocyt_pct")
pheno <- pheno[complete.cases(pheno[, ..cc_cols])]
n_after_cc <- nrow(pheno)
pheno <- pheno[rand_id %in% colnames(mat)]
n_analytic <- nrow(pheno)
mat_ana <- mat[, pheno$rand_id, drop = FALSE]
stopifnot(identical(colnames(mat_ana), pheno$rand_id))

pheno[, `:=`(
  age_visit    = as.numeric(age_visit),
  race         = factor(race),
  gender       = factor(gender),
  SmokCigNow   = factor(SmokCigNow),
  Lc.Batch     = factor(Lc.Batch),
  DisCOPD      = factor(DisCOPD,      levels = c("0","1","2","3")),
  DisCOPDAI    = factor(DisCOPDAI,    levels = c("0","1","2","3")),
  DisCOPDQuant = factor(DisCOPDQuant, levels = c("0","1","2","3"))
)]
pheno_df <- as.data.frame(pheno)
rownames(pheno_df) <- pheno_df$rand_id
pheno_df$sid <- pheno_df$rand_id

dge <- DGEList(counts = mat_ana, samples = pheno_df,
               genes  = data.frame(id = rownames(mat_ana), stringsAsFactors = FALSE))
n_lc_batch_levels <- length(unique(pheno_df$Lc.Batch))

# ============================================================================
# Source:     chunk: libsize
# Inputs:     dge
# Population: n_analytic subjects
# Methods:    METHODS.md §3 — pre-normalization QC
# ============================================================================
libsizes <- data.frame(rand_id = colnames(dge),
                       lib_size = colSums(dge$counts),
                       Lc.Batch = dge$samples$Lc.Batch)
ggplot(libsizes, aes(x = lib_size / 1e6)) +
  geom_histogram(bins = 60, fill = "steelblue", color = "white") +
  scale_x_continuous(labels = scales::comma) +
  labs(title = sprintf("Library size distribution (n = %d samples)", nrow(libsizes)),
       subtitle = sprintf("median = %.1f M reads, range %.1f – %.1f M",
                          median(libsizes$lib_size) / 1e6,
                          min(libsizes$lib_size) / 1e6,
                          max(libsizes$lib_size) / 1e6),
       x = "Library size (millions of reads)", y = "Sample count") +
  theme_minimal(base_size = 11)

# ============================================================================
# Source:     chunk: filter-norm
# Inputs:     dge (analytic), pheno_df
# Output:     dge (filtered + TMM-normalized)
# Methods:    METHODS.md §3
# ============================================================================
design_filt <- model.matrix(~ age_visit + race + gender + SmokCigNow + Lc.Batch,
                            data = dge$samples)
keep_g <- filterByExpr(dge, design = design_filt,
                       min.count = 50, min.prop = 0.9)
dge <- dge[keep_g, , keep.lib.sizes = TRUE]
n_genes_kept <- nrow(dge)
dge <- calcNormFactors(dge, method = "TMM")

nf <- data.frame(norm_factor = dge$samples$norm.factors)
ggplot(nf, aes(x = norm_factor)) +
  geom_histogram(bins = 50, fill = "darkorange", color = "white") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  labs(title = sprintf("TMM normalization factors (n = %d samples)", nrow(nf)),
       subtitle = sprintf("median = %.3f, IQR [%.3f, %.3f]",
                          median(nf$norm_factor),
                          quantile(nf$norm_factor, 0.25),
                          quantile(nf$norm_factor, 0.75)),
       x = "Norm factor", y = "Sample count") +
  theme_minimal(base_size = 11)

# ============================================================================
# Source:     chunk: gmap
# Inputs:     rownames(dge); cache/gmap_v103.rds (built from
#             AnnotationHub AH89426 = Ensembl v103 EnsDb, matching the
#             original analysis's gmap_ENSGv103_2025.6.09.rda).
# Output:     gmap (ENSEMBL → SYMBOL, ENTREZ); id_map (ENSEMBL → ENTREZ for ranks)
# Methods:    METHODS.md §3 — gene annotation: Ensembl v103 (matches the
#             GENCODE v37 reference upstream).
# ============================================================================
ens_full <- rownames(dge)
ens_bare <- sub("\\.[0-9]+$", "", ens_full)
v103_rds <- file.path(CACHE_DIR, "gmap_v103.rds")
if (file.exists(v103_rds)) {
  gmap <- readRDS(v103_rds)$gmap
  if (!identical(gmap$ensembl_gene_id_version, ens_full)) {
    file.remove(v103_rds); gmap <- NULL
  }
} else gmap <- NULL

if (is.null(gmap)) {
  suppressPackageStartupMessages({
    library(AnnotationHub); library(ensembldb)
  })
  ah <- AnnotationHub()
  edb <- ah[["AH89426"]]
  genes_df <- as.data.frame(genes(edb,
    columns = c("gene_id","gene_name","entrezid","gene_biotype")))
  genes_df$entrez_first <- vapply(genes_df$entrezid, function(x) {
    x <- as.character(x); x <- x[!is.na(x) & x != ""]
    if (length(x)) x[1] else NA_character_
  }, character(1))
  hits <- genes_df[match(ens_bare, genes_df$gene_id), ]
  gmap <- data.frame(
    ensembl_gene_id_version = ens_full,
    ensembl_gene_id         = ens_bare,
    external_gene_name      = hits$gene_name,
    entrezgene_id           = hits$entrez_first,
    stringsAsFactors        = FALSE)
  saveRDS(list(source = "AnnotationHub AH89426 (Ensembl v103 EnsDb)",
               gmap = gmap), v103_rds)
}
n_with_symbol <- sum(!is.na(gmap$external_gene_name) & gmap$external_gene_name != "")
n_with_entrez <- sum(!is.na(gmap$entrezgene_id) & gmap$entrezgene_id != "")

id_map <- as.data.table(gmap)[!is.na(entrezgene_id) & entrezgene_id != "",
                              .(ensembl = ensembl_gene_id_version,
                                entrez  = as.character(entrezgene_id))]
setkey(id_map, ensembl)

# ============================================================================
# Source:     chunk: pca
# Inputs:     dge (TMM-normalized)
# Population: n_analytic subjects
# Methods:    METHODS.md §4.2 — diagnostic only, not used in DE
# Note:       PCA on the top 2000 most variable AUTOSOMAL genes — X/Y/MT
#             excluded so the diagnostic reflects batch structure rather
#             than gender / mitochondrial content.
# ============================================================================
logcpm_diag <- cpm(dge, log = TRUE, prior.count = 2)
ens_bare_pca <- sub("\\.[0-9]+$", "", rownames(logcpm_diag))
# Use the same Ensembl v103 EnsDb that built `gmap` for chromosome lookup —
# consistent annotation source across the whole report (METHODS.md §3).
chr_rds <- file.path(CACHE_DIR, "chr_v103.rds")
chr_src <- "AnnotationHub AH89426 (Ensembl v103 EnsDb)"
chr_obj <- if (file.exists(chr_rds)) readRDS(chr_rds) else NULL
# Validate cache: must be a v103-tagged list with a sufficient gene universe.
if (!is.list(chr_obj) || !identical(chr_obj$source, chr_src) ||
    !all(ens_bare_pca %in% names(chr_obj$chr))) chr_obj <- NULL

if (is.null(chr_obj)) {
  suppressPackageStartupMessages({ library(AnnotationHub); library(ensembldb) })
  edb_pca <- AnnotationHub()[["AH89426"]]
  genes_pca <- as.data.frame(genes(edb_pca, columns = c("gene_id")))
  chr_lookup <- setNames(as.character(genes_pca$seqnames), genes_pca$gene_id)
  chr_obj <- list(source = chr_src, chr = chr_lookup)
  saveRDS(chr_obj, chr_rds)
}
chr_v <- chr_obj$chr[ens_bare_pca]
autosomal_mask <- !is.na(chr_v) & chr_v %in% as.character(1:22)
n_autosomal <- sum(autosomal_mask)

vars_all <- matrixStats::rowVars(logcpm_diag)
vars_all[!autosomal_mask] <- -Inf
top_var_idx <- order(vars_all, decreasing = TRUE)[1:2000]
pca_diag <- prcomp(t(logcpm_diag[top_var_idx, ]), center = TRUE, scale. = FALSE)
pca_df <- data.frame(PC1 = pca_diag$x[, 1], PC2 = pca_diag$x[, 2],
                     Lc.Batch = pheno_df$Lc.Batch,
                     DisCOPD  = pheno_df$DisCOPD)

p_batch <- ggplot(pca_df, aes(PC1, PC2, color = Lc.Batch)) +
  geom_point(alpha = 0.5, size = 0.7) +
  guides(color = "none") +
  labs(title = sprintf("Lc.Batch (n = %d, %d batches)", n_analytic, n_lc_batch_levels)) +
  theme_minimal(base_size = 10)
p_class <- ggplot(pca_df, aes(PC1, PC2, color = DisCOPD)) +
  geom_point(alpha = 0.5, size = 0.7) +
  scale_color_brewer(palette = "Set1", na.value = "grey80") +
  labs(title = sprintf("DisCOPD class (n = %d)", n_analytic)) +
  theme_minimal(base_size = 10)
(p_batch | p_class) +
  plot_annotation(title = "PCA of log-CPM colored by Lc.Batch and by DisCOPD class",
                  subtitle = sprintf("Top 2000 most variable AUTOSOMAL genes (X/Y/MT excluded; %d autosomal genes available)",
                                     n_autosomal))

# ============================================================================
# Source:     chunk: voom-demo
# Inputs:     dge (TMM-normalized), pheno_df
# Population: n_analytic subjects
# Output:     voom_demo (voom output, feeds the pre-fit mean-variance plot);
#             demo_fit2 (lmFit + eBayes on the same demo design, feeds the
#             post-fit residual-SD diagnostic plot below).
# Methods:    METHODS.md §4.2, §6 — voom precision weights
# Note:       Representative voom + lmFit + eBayes call (Model 1, DisCOPD
#             design) used by both the pre-fit and post-fit diagnostic
#             plots only. The de-gsea chunk refits voom per
#             (Model × Classification) cell for the actual DE results.
# ============================================================================
demo_design <- model.matrix(~ Lc.Batch + age_visit + race + gender + SmokCigNow + DisCOPD,
                            data = pheno_df)
voom_demo   <- voom(dge, demo_design, plot = FALSE, save.plot = TRUE)

# Also fit lmFit + eBayes on the same demo design so the post-fit mean-variance
# diagnostic (plotSA) below uses the actual fitted object.
demo_fit  <- lmFit(voom_demo, demo_design)
demo_fit2 <- eBayes(demo_fit)

# Reproduce voom's mean-variance plot using the saved data
mv <- voom_demo$voom.xy
ml <- voom_demo$voom.line
ggplot(data.frame(x = mv$x, y = mv$y), aes(x, y)) +
  geom_point(alpha = 0.10, size = 0.4, color = "darkblue") +
  geom_line(data = data.frame(x = ml$x, y = ml$y), aes(x, y),
            color = "red", linewidth = 0.7) +
  labs(title = sprintf("voom mean-variance trend (n = %d genes, %d samples; Model 1 / DisCOPD design)",
                       nrow(dge), n_analytic),
       subtitle = "After voom precision weights, residual heteroscedasticity is captured by the trend (red curve)",
       x = expression(log[2]~"(count size + 0.5)"),
       y = expression(sqrt(sigma))) +
  theme_minimal(base_size = 11)

# ============================================================================
# Source:     chunk: voom-postfit-plot
# Inputs:     demo_fit2 (lmFit fit object; we read $Amean and $sigma which
#             are the pre-eBayes residual statistics — eBayes runs after
#             but doesn't change either quantity).
# Output:     post-fit sqrt(residual SD) vs AveExpr plot (limma::plotSA convention)
# Methods:    METHODS.md §6 — voom + lmFit; this plot confirms the
#             mean-variance trend is corrected after voom-weighted lmFit
#             (independent of any subsequent eBayes shrinkage).
# ============================================================================
sa_df <- data.frame(AveExpr = demo_fit2$Amean,
                    sqrt_sigma = sqrt(demo_fit2$sigma))
ggplot(sa_df, aes(x = AveExpr, y = sqrt_sigma)) +
  geom_point(alpha = 0.10, size = 0.4, color = "darkgreen") +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 0.7) +
  labs(title = sprintf("Post-fit sqrt(residual SD) vs average expression (n = %d genes; Model 1 / DisCOPD design)",
                       nrow(demo_fit2)),
       subtitle = "After voom-weighted lmFit (pre-eBayes); flat trend confirms heteroscedasticity is corrected — same y-axis convention as limma::plotSA",
       x = expression("AveExpr (mean " * log[2] ~ "CPM)"),
       y = expression(sqrt(sigma) * " = " * sqrt("residual SD"))) +
  theme_minimal(base_size = 11)

# ============================================================================
# Source:     chunk: pathways
# Inputs:     msigdbr (Hallmark + Reactome + KEGG_MEDICUS + ImmuneSigDB + BioCarta)
# Output:     all_pathways (named list of Entrez vectors), pathway_source
# Methods:    METHODS.md §7
# ============================================================================
hall  <- msigdbr(species = "Homo sapiens", collection = "H")
react <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
kegg  <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS")
imm   <- msigdbr(species = "Homo sapiens", collection = "C7", subcollection = "IMMUNESIGDB")
bio   <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:BIOCARTA")

build_list <- function(df) split(as.character(df$ncbi_gene), df$gs_name)
hall_list  <- build_list(hall)
react_list <- build_list(react)
kegg_list  <- build_list(kegg)
imm_list   <- build_list(imm)
bio_list   <- build_list(bio)

# Combined pathway list, with source tracked. Dedup by gs_name (Reactome
# appears in original DE GSEA and in singscore's collections — keep first).
pw_src_pairs <- c(
  setNames(rep("Hallmark",     length(hall_list)),  names(hall_list)),
  setNames(rep("Reactome",     length(react_list)), names(react_list)),
  setNames(rep("KEGG_MEDICUS", length(kegg_list)),  names(kegg_list)),
  setNames(rep("ImmuneSigDB",  length(imm_list)),   names(imm_list)),
  setNames(rep("BioCarta",     length(bio_list)),   names(bio_list))
)
pw_src_pairs <- pw_src_pairs[!duplicated(names(pw_src_pairs))]

all_pathways_full <- c(hall_list, react_list, kegg_list, imm_list, bio_list)
all_pathways      <- all_pathways_full[!duplicated(names(all_pathways_full))]
pathway_source    <- pw_src_pairs[names(all_pathways)]
n_pathways_loaded <- length(all_pathways)

# ============================================================================
# Source:     chunk: de-gsea
# Inputs:     dge (TMM-normalized counts), pheno_df, all_pathways, id_map
# Population: n_analytic subjects per cell
# Output:     all_de, all_gsea, de_summary, gsea_summary
# Methods:    METHODS.md §6 (DE), §7 (GSEA)
# ============================================================================

# ---- 1. Covariate strings for the two models ----
# Right-hand side of the design formula EXCLUDING the classification
# variable (which is appended per cell in the inner loop). Lc.Batch is
# included as a covariate here (NOT pre-removed from expression) — voom
# needs to see batch in the design so the precision weights are
# calibrated against post-batch residuals, not the gross sample-level
# variance.
cov_m1 <- "Lc.Batch + age_visit + race + gender + SmokCigNow"
cov_m2 <- paste(cov_m1, "+ eosinphl_pct + neutrophl_pct + lymphcyt_pct + monocyt_pct")
model_configs <- list(Model_1 = cov_m1, Model_2 = cov_m2)

# ---- 2. Helper: per-contrast topTable extraction ----
# Builds the per-cell DE results frame: topTable + symbol annotation
# + display rounding + sort by FDR. Note `t` is intentionally NOT rounded
# (the de-concordance plot reads it at full precision; rounding to 4 sig
# figs degrades Pearson r at the 5th decimal).
extract_table <- function(fit_obj, coef_name, gmap_df, n = Inf) {
  tt <- topTable(fit_obj, coef = coef_name, number = n,
                 confint = TRUE, sort.by = "none")
  tt$id <- rownames(tt)
  tt$gene_label <- gmap_df$external_gene_name[match(tt$id, gmap_df$ensembl_gene_id_version)]
  tt <- as_tibble(tt) |>
    mutate(across(c(logFC, P.Value, adj.P.Val, AveExpr, CI.L, CI.R),
                  ~ signif(as.numeric(.x), 4))) |>
    arrange(adj.P.Val) |>
    dplyr::select(gene_label, id, logFC, AveExpr, t, P.Value, adj.P.Val, CI.L, CI.R)
  tt
}

# ---- 3. Helper: build the GSEA rank vector ----
# fgsea needs a named numeric vector of statistics, one entry per Entrez
# ID, sorted decreasing. We start from per-gene moderated t-statistics
# (named by ENSEMBL versioned ID), join to id_map (ENSEMBL → Entrez),
# and resolve Entrez collisions (multi-ENSEMBL → same Entrez) by
# keeping the gene with the largest |t|. Stable order ties pick the
# first ENSEMBL in id_map's order, matching the original analysis.
build_ranks <- function(t_named, id_map) {
  dt <- data.table(ensembl = names(t_named), stat = t_named)
  setkey(dt, ensembl)
  m <- id_map[dt, nomatch = 0]
  m_unique <- m[order(-abs(stat)), .SD[1], by = entrez]
  ranks <- m_unique$stat
  names(ranks) <- m_unique$entrez
  sort(ranks, decreasing = TRUE)
}

# ---- 4. Per-cell loop ----
# Two outer loops (model × classification) build the fit; an inner loop
# extracts each of the 6 contrasts. Per-cell results accumulate in 4
# named lists keyed by "Model|Classification|Contrast".
all_de <- list(); all_gsea <- list()
all_de_sum <- list(); all_gsea_sum <- list()

for (model_label in names(model_configs)) {
  cov <- model_configs[[model_label]]

  for (cls_var in cls_vars) {

    # ---- 4a. Subset to subjects with non-NA classification ----
    # Complete-case upstream guarantees this is a no-op (all 3
    # classifications non-NA in the analytic cohort). Kept for safety.
    keep_s <- !is.na(pheno_df[[cls_var]])
    ph_s   <- droplevels(pheno_df[keep_s, , drop = FALSE])
    dge_s  <- dge[, keep_s]

    # ---- 4b. Design matrix ----
    # ~ <model covariates> + <classification factor>. R's treatment-
    # coding produces dummy columns named `<cls_var>1`, `<cls_var>2`,
    # `<cls_var>3` (level 0 is the reference). `make.names` ensures
    # column names are syntactically valid R identifiers.
    design_mat <- model.matrix(as.formula(paste("~", cov, "+", cls_var)),
                               data = ph_s)
    colnames(design_mat) <- make.names(colnames(design_mat))

    # ---- 4c. voom refit + lmFit ----
    # voom() estimates the count mean-variance trend and produces
    # precision weights per (gene, sample). lmFit then does
    # weighted-least-squares per gene. We refit voom per cell so the
    # weights reflect THIS design's residuals — refitting matters when
    # the design changes between cells (e.g., Model 2 adds CBC).
    v   <- voom(dge_s, design_mat, plot = FALSE)
    fit <- lmFit(v, design_mat)

    # ---- 4d. Contrast matrix ----
    # Same 6 contrasts as the singscore report — three "level vs ref 0"
    # taken directly from coefficients, three pairwise differences via
    # subtraction. Each is added only if the involved levels are
    # present in this cohort (defensive against partial cohorts).
    lvls    <- levels(ph_s[[cls_var]])
    present <- intersect(c("1","2","3"), lvls)
    cs <- list()
    for (k in present) cs[[paste0("Grp", k, "vs0")]] <- paste0(cls_var, k)
    if (all(c("1","3") %in% present)) cs[["Grp1vs3"]] <- paste0(cls_var,"1 - ", cls_var,"3")
    if (all(c("2","3") %in% present)) cs[["Grp2vs3"]] <- paste0(cls_var,"2 - ", cls_var,"3")
    if (all(c("1","2") %in% present)) cs[["Grp1vs2"]] <- paste0(cls_var,"1 - ", cls_var,"2")
    contrast_mtx <- do.call(makeContrasts, c(cs, list(levels = design_mat)))
    # Rename "Grp1vs0" → "DisCOPD_1vs0" etc., so the per-cell key
    # uniquely identifies model + classification + contrast downstream.
    colnames(contrast_mtx) <- paste0(cls_var, "_", sub("^Grp","", names(cs)))

    # ---- 4e. Apply contrasts + empirical Bayes ----
    # contrasts.fit rotates SEs into the contrast basis without
    # refitting; eBayes shrinks per-gene residual variances toward a
    # common trend, which is what makes a per-gene t-test on a small
    # number of df powerful enough to be useful.
    fit2 <- eBayes(contrasts.fit(fit, contrast_mtx))

    # ---- 4f. Per-contrast extraction (DE + GSEA) ----
    for (cn in colnames(contrast_mtx)) {
      key <- paste(model_label, cls_var, cn, sep = "|")

      # ---- DE topTable (display + per-cell summary) ----
      tt <- extract_table(fit2, cn, gmap)
      all_de[[key]] <- tt
      all_de_sum[[key]] <- data.frame(
        Model = model_label, Classification = cls_var, Contrast = cn,
        N_samples    = nrow(ph_s),
        Genes_Tested = nrow(tt),
        Genes_Nominal_05 = sum(tt$P.Value   < 0.05, na.rm = TRUE),
        Genes_FDR_05     = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
        stringsAsFactors = FALSE)

      # ---- GSEA on full-precision t-statistic ranks ----
      # Read from fit2$t (NOT the rounded tt$t) so the rank order is
      # at full precision. set.seed(42) inside the loop gives per-cell
      # reproducibility (each cell uses the same RNG sequence).
      t_named <- setNames(as.numeric(fit2$t[, cn]), rownames(fit2$t))
      ranks   <- build_ranks(t_named, id_map)
      set.seed(42)
      gsea <- fgsea(pathways = all_pathways, stats = ranks,
                    minSize = 15, maxSize = 500, eps = 0)
      gsea[, source := pathway_source[pathway]]
      # leadingEdge is a list-column — flatten to ";" delimited string
      # so the result fits in a CSV row.
      gsea[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
      setorder(gsea, padj)
      all_gsea[[key]] <- gsea
      all_gsea_sum[[key]] <- data.frame(
        Model = model_label, Classification = cls_var, Contrast = cn,
        N_pathways_tested  = nrow(gsea),
        Pathways_Nominal_05 = sum(gsea$pval < 0.05, na.rm = TRUE),
        Pathways_FDR_05     = sum(gsea$padj < 0.05, na.rm = TRUE),
        stringsAsFactors = FALSE)
    }
  }
}

# ---- 5. Collapse per-cell summaries ----
# `de_summary` and `gsea_summary` each have one row per cell; the
# `short_contrast` helper column strips the `<cls_var>_` prefix so
# downstream bar plots can use a clean x-axis.
de_summary   <- bind_rows(all_de_sum)   |>
  mutate(short_contrast = sub("^[^_]+_", "", Contrast))
gsea_summary <- bind_rows(all_gsea_sum) |>
  mutate(short_contrast = sub("^[^_]+_", "", Contrast))

DT::datatable(de_summary[, c("Model","Classification","Contrast",
                             "N_samples","Genes_Tested","Genes_FDR_05")],
              caption = sprintf("Gene-level DE significance — %d cells (n = %d subjects, %d genes tested per cell)",
                                nrow(de_summary), n_analytic, n_genes_kept),
              options = list(pageLength = 36, scrollX = TRUE),
              rownames = FALSE, filter = "top")

ggplot(de_summary, aes(x = short_contrast, y = Genes_FDR_05, fill = Classification)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.75) +
  facet_wrap(~ Model, ncol = 1) +
  labs(title = sprintf("FDR-significant DE genes per cell (n = %d subjects, %d genes tested)",
                       n_analytic, n_genes_kept),
       subtitle = "adj.P.Val < 0.05",
       x = "Contrast", y = "Number of genes") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")

phist_df <- bind_rows(lapply(names(all_de), function(k) {
  parts <- strsplit(k, "[|]")[[1]]
  data.frame(P.Value = all_de[[k]]$P.Value,
             Model   = parts[1], Classification = parts[2], Contrast = parts[3],
             stringsAsFactors = FALSE)
}))
ggplot(phist_df, aes(x = P.Value)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  facet_grid(Contrast ~ Model + Classification, scales = "free_y") +
  labs(title = sprintf("Per-cell p-value distributions — DE (n = %d genes)", n_genes_kept),
       subtitle = "Uniform left tail with spike at 0 = well-specified",
       x = "P.Value", y = "Genes") +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(size = 7))

DT::datatable(gsea_summary[, c("Model","Classification","Contrast",
                               "N_pathways_tested","Pathways_FDR_05")],
              caption = sprintf("GSEA significance — %d cells (n = %d subjects; per cell, fgsea tests pathways with overlap-size in [15, 500] from a library of %d loaded — column N_pathways_tested has the per-cell value)",
                                nrow(gsea_summary), n_analytic, n_pathways_loaded),
              options = list(pageLength = 36, scrollX = TRUE),
              rownames = FALSE, filter = "top")

ggplot(gsea_summary, aes(x = short_contrast, y = Pathways_FDR_05, fill = Classification)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.75) +
  facet_wrap(~ Model, ncol = 1) +
  labs(title = sprintf("FDR-significant GSEA pathways per cell (n = %d subjects, %d pathways)",
                       n_analytic, n_pathways_loaded),
       subtitle = "padj < 0.05",
       x = "Contrast", y = "Number of pathways") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")

gs_breakdown <- bind_rows(lapply(names(all_gsea), function(k) {
  parts <- strsplit(k, "[|]")[[1]]
  sig <- as.data.frame(all_gsea[[k]]) |> dplyr::filter(padj < 0.05)
  if (!nrow(sig)) return(NULL)
  sig |> dplyr::count(source, name = "n_FDR_sig") |>
    mutate(Model = parts[1], Classification = parts[2], Contrast = parts[3])
}))
if (nrow(gs_breakdown)) {
  DT::datatable(gs_breakdown[, c("Model","Classification","Contrast",
                                 "source","n_FDR_sig")],
                caption = "FDR-significant GSEA pathways by collection",
                options = list(pageLength = 25, scrollX = TRUE),
                rownames = FALSE, filter = "top")
} else {
  cat("(Data) No FDR-significant GSEA pathways across any cell.\n")
}

de_sets <- lapply(cls_vars, function(v) {
  pat <- paste0("^Model_1[|]", v, "[|]")
  ks  <- grep(pat, names(all_de), value = TRUE)
  unique(unlist(lapply(ks, function(k) {
    x <- all_de[[k]]; x$id[x$adj.P.Val < 0.05]
  })))
})
names(de_sets) <- cls_vars
mat_de <- outer(seq_along(cls_vars), seq_along(cls_vars),
                Vectorize(function(i, j) length(intersect(de_sets[[i]], de_sets[[j]]))))
dimnames(mat_de) <- list(cls_vars, cls_vars)
kable(mat_de, caption = "Pairwise overlap of FDR-significant DE genes (Model 1, union over the 6 contrasts per classification). Diagonal = total set size for that classification (any contrast); off-diagonals = gene intersections, counted across any contrast pair (not within-matched-contrast).")

gs_sets <- lapply(cls_vars, function(v) {
  pat <- paste0("^Model_1[|]", v, "[|]")
  ks  <- grep(pat, names(all_gsea), value = TRUE)
  unique(unlist(lapply(ks, function(k) {
    x <- as.data.frame(all_gsea[[k]]); x$pathway[x$padj < 0.05]
  })))
})
names(gs_sets) <- cls_vars
mat_gs <- outer(seq_along(cls_vars), seq_along(cls_vars),
                Vectorize(function(i, j) length(intersect(gs_sets[[i]], gs_sets[[j]]))))
dimnames(mat_gs) <- list(cls_vars, cls_vars)
kable(mat_gs, caption = "Pairwise overlap of FDR-significant GSEA pathways (Model 1, union over the 6 contrasts per classification). Diagonal = total set size for that classification (any contrast); off-diagonals = pathway intersections, counted across any contrast pair (not within-matched-contrast).")

t_lookup <- lapply(cls_vars, function(v) {
  k <- paste0("Model_1|", v, "|", v, "_1vs0")
  x <- all_de[[k]]
  setNames(x$t, x$id)
})
names(t_lookup) <- cls_vars
pairs_l <- combn(cls_vars, 2, simplify = FALSE)
plots <- lapply(pairs_l, function(pr) {
  ids <- intersect(names(t_lookup[[pr[1]]]), names(t_lookup[[pr[2]]]))
  d <- data.frame(x = t_lookup[[pr[1]]][ids],
                  y = t_lookup[[pr[2]]][ids])
  d <- d[is.finite(d$x) & is.finite(d$y), ]
  r <- cor(d$x, d$y, use = "complete.obs")
  ggplot(d, aes(x, y)) +
    geom_point(alpha = 0.10, size = 0.4) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = sprintf("%s vs %s (n = %d genes)", pr[1], pr[2], nrow(d)),
         subtitle = sprintf("Pearson r = %.3f", r),
         x = paste("t-stat", pr[1]), y = paste("t-stat", pr[2])) +
    theme_minimal(base_size = 10)
})
wrap_plots(plots, ncol = 3) +
  plot_annotation(title = sprintf("Gene t-stat concordance — Model 1, contrast 1vs0 (n = %d subjects)",
                                  n_analytic))

for (model_label in unique(de_summary$Model)) {
  cat("\n\n## ", model_label, " {.tabset}\n\n", sep = "")
  for (v in cls_vars) {
    cat("\n\n### ", cls_labels[v], " (`", v, "`) {.tabset}\n\n", sep = "")
    pat  <- paste0("^", model_label, "[|]", v, "[|]")
    keys <- grep(pat, names(all_de), value = TRUE)
    for (k in keys) {
      cn <- strsplit(k, "[|]")[[1]][3]
      cat("\n\n#### ", cn, "\n\n", sep = "")

      tt <- all_de[[k]]
      cat("\n**Top 20 differentially expressed genes (n = ", n_analytic, ")**\n\n", sep = "")
      print(kable(head(tt, 20),
                  digits = 3,
                  caption = sprintf("DE — %s | %s | %s", model_label, v, cn)))

      gg2 <- as.data.frame(all_gsea[[k]]) |> arrange(padj)
      cat("\n\n**Top 20 GSEA pathways (n = ", n_analytic, ")**\n\n", sep = "")
      print(kable(head(gg2[, c("pathway","source","ES","NES","pval","padj","size")], 20),
                  digits = 3,
                  caption = sprintf("GSEA — %s | %s | %s", model_label, v, cn)))
      cat("\n")
    }
  }
}

# ============================================================================
# Source:     chunk: exports
# Inputs:     all_de, all_gsea, de_summary, gsea_summary, all_pathways
# Output:     CSVs in results/ and results/{gene_de,gsea}/;
#             standalone analyzed-pathways RDS in results/;
#             full results bundle RDS in cache/
# ============================================================================
fwrite(de_summary[, c("Model","Classification","Contrast",
                      "N_samples","Genes_Tested",
                      "Genes_Nominal_05","Genes_FDR_05")],
       file.path(RES_DIR, paste0("gene_de_summary_", DATE_STAMP, ".csv")))
fwrite(gsea_summary[, c("Model","Classification","Contrast",
                        "N_pathways_tested",
                        "Pathways_Nominal_05","Pathways_FDR_05")],
       file.path(RES_DIR, paste0("gsea_summary_", DATE_STAMP, ".csv")))

for (model_label in names(model_configs)) {
  for (cls_var in cls_vars) {
    pat  <- paste0("^", model_label, "[|]", cls_var, "[|]")
    keys <- grep(pat, names(all_de), value = TRUE)
    if (!length(keys)) next
    combined_de <- bind_rows(lapply(keys, function(k) {
      parts <- strsplit(k, "[|]")[[1]]
      all_de[[k]] |> mutate(Model = parts[1], Contrast = parts[3])
    }))
    fwrite(combined_de, file.path(RES_DIR, "gene_de",
            sprintf("gene_de_%s_%s_allresults_%s.csv",
                    model_label, cls_var, DATE_STAMP)))
    combined_gs <- rbindlist(lapply(keys, function(k) {
      parts <- strsplit(k, "[|]")[[1]]
      cbind(all_gsea[[k]], Model = parts[1], Contrast = parts[3])
    }))
    fwrite(combined_gs, file.path(RES_DIR, "gsea",
            sprintf("gsea_%s_%s_allresults_%s.csv",
                    model_label, cls_var, DATE_STAMP)))
  }
}

# Standalone RDS of the analyzed pathway library — convenient for
# downstream consumers who want the Entrez membership lists without
# loading the full results bundle.
gsea_analyzed_pathways <- list(
  pathway_list    = all_pathways,
  pathway_source  = pathway_source,
  n_pathways      = length(all_pathways),
  collections     = c("Hallmark (H)",
                      "Reactome (C2:CP:REACTOME)",
                      "KEGG Medicus (C2:CP:KEGG_MEDICUS)",
                      "ImmuneSigDB (C7:IMMUNESIGDB)",
                      "BioCarta (C2:CP:BIOCARTA)"),
  msigdbr_version = as.character(packageVersion("msigdbr")),
  filter_rule     = paste("fgsea minSize=15 maxSize=500;",
                          "library deduplicated by gs_name across collections"),
  date_stamp      = DATE_STAMP
)
saveRDS(gsea_analyzed_pathways,
        file.path(RES_DIR, paste0("gsea_analyzed_pathways_",
                                  DATE_STAMP, ".rds")))

saveRDS(list(all_de = all_de, all_gsea = all_gsea,
             de_summary = de_summary, gsea_summary = gsea_summary,
             pathway_source = pathway_source, pheno = pheno_df,
             cls_vars = cls_vars, cls_labels = cls_labels,
             date_stamp = DATE_STAMP),
        file.path(CACHE_DIR, "gene_de_gsea_results.rds"))