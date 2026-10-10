# ============================================================================
# COPDGene_singscore_Workflow.R
#
# Participant-level pathway scores (singscore) for the multidimensional COPD
# diagnostic schema, regressed on the same models and contrasts.
#   1. Load de-identified counts, RNA master file, Phase 2 phenotypes and the
#      three schema classifications (original visual, AI-derived visual CT,
#      quantitative CT); complete-case cohort (n = 3,309)
#   2. filterByExpr (min.count 50, min.prop 0.9), TMM normalization, QC
#   3. batch-corrected log-CPM; singscore (undirected; sets of 15-500
#      measured genes; knownDirection = FALSE); scores regressed with limma per Model x Classification:
#        Model_1: batch + age + race + sex + smoking status
#        Model_2: Model_1 + eosinophil, neutrophil, lymphocyte, monocyte %
#      contrasts: each disease group vs no COPD, and pairwise disease groups
#   4. Benjamini-Hochberg control across 5,995 sets
#
# Outputs (RES_DIR): singscore/*_allresults_<date>.csv, singscore_summary_<date>.csv,
#   singscore_score_matrix_<date>.csv
#
# Code recovered from the executed notebook bhatt_repro_singscore_2026.
# R 4.5; limma 3.66.0; edgeR; fgsea 1.36.2; msigdbr 25.1.1.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma); library(fgsea)
  library(msigdbr); library(singscore); library(GSEABase); library(dplyr); library(tibble); library(ggplot2); library(DT)
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
# Source:     bhatt_repro_singscore_2026.Rmd / chunk: preprocess
# Inputs:     counts_raw_freeze5_randid.csv,
#             masterfile_rna_freeze5_randid.csv,
#             COPDGene_P1P2P3_SM_NS_ILDBR_Long_Sep24_randid.csv,
#             {MultidimensionalSchema, NewClassQualitativeAI, NewClassQuantitative}_deid.csv
# Population: visit-2 RNA samples, complete-case on all 3 classifications +
#             CBC + smoking → n = 3309 expected
# Output:     pheno_df, dge (DGEList), gmap (ENSEMBL→symbol)
# Methods:    METHODS.md §2-3
# ============================================================================

# 1. Counts (cached as RDS so we don't re-read the 770 MB CSV)
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
n_counts_cols    <- ncol(mat)
n_genes_raw      <- nrow(mat)

# 2. Masterfile filter
mf <- fread(file.path(DATA_DIR, "masterfile_rna_freeze5_randid.csv"))
n_master_total   <- nrow(mf)
mf <- mf[final.analyzed == 1, .(rand_id, Lc.Batch)]
mf[, rand_id := as.character(rand_id)]
n_master_keep    <- nrow(mf)

# 3. Phenotype, visit 2
ph <- fread(file.path(DATA_DIR, "COPDGene_P1P2P3_SM_NS_ILDBR_Long_Sep24_randid.csv"),
            select = c("rand_id","visitnum","age_visit","race","gender","SmokCigNow",
                       "eosinphl_pct","neutrophl_pct","lymphcyt_pct","monocyt_pct"))
n_pheno_total    <- nrow(ph)
ph <- ph[visitnum == 2]; ph[, rand_id := as.character(rand_id)]
n_pheno_v2       <- nrow(ph)

# 4. Classification files
mds  <- fread(file.path(DATA_DIR, "MultidimensionalSchema_deid.csv"))[, .(rand_id, DisCOPD)]
cqai <- fread(file.path(DATA_DIR, "NewClassQualitativeAI_deid.csv"))[, .(rand_id, DisCOPDAI)]
cqt  <- fread(file.path(DATA_DIR, "NewClassQuantitative_deid.csv"))[, .(rand_id, DisCOPDQuant)]
for (x in list(mds, cqai, cqt)) x[, rand_id := as.character(rand_id)]

# 5. Merge + complete-case on classifications + CBC + smoking
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
n_after_cc       <- nrow(pheno)

# 6. Restrict to samples present in counts
pheno <- pheno[rand_id %in% colnames(mat)]
n_analytic       <- nrow(pheno)
mat_ana <- mat[, pheno$rand_id, drop = FALSE]
stopifnot(identical(colnames(mat_ana), pheno$rand_id))

# 7. Coerce variable types
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

# 8. DGEList (no filter / norm yet — those are the next chunks)
dge <- DGEList(counts = mat_ana, samples = pheno_df,
               genes  = data.frame(id = rownames(mat_ana), stringsAsFactors = FALSE))
n_lc_batch_levels <- length(unique(pheno_df$Lc.Batch))

# ============================================================================
# Source:     chunk: libsize
# Inputs:     dge (analytic DGEList)
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
# Source:     chunk: filter
# Inputs:     dge
# Population: n_analytic subjects, n_genes_raw genes
# Output:     dge (filtered), n_genes_kept
# Methods:    METHODS.md §3 — filterByExpr min.count=50, min.prop=0.9
# ============================================================================
design_filt <- model.matrix(~ age_visit + race + gender + SmokCigNow + Lc.Batch,
                            data = dge$samples)
keep_g <- filterByExpr(dge, design = design_filt,
                       min.count = 50, min.prop = 0.9)
dge <- dge[keep_g, , keep.lib.sizes = TRUE]
n_genes_kept <- nrow(dge)

# ============================================================================
# Source:     chunk: normalize
# Inputs:     dge (filtered)
# Output:     dge (TMM-normalized), logcpm, logcpm_corrected
# Methods:    METHODS.md §3 (TMM, log-CPM) + §4.1 (singscore batch handling)
# ============================================================================
dge <- calcNormFactors(dge, method = "TMM")
logcpm           <- cpm(dge, log = TRUE, prior.count = 2)
logcpm_corrected <- removeBatchEffect(logcpm, batch = dge$samples$Lc.Batch)

# Mean-centering note: the post-4.08 shared preprocessing in the original
# analysis added a per-gene mean-centering step here (subtract rowMeans).
# We deliberately do NOT apply mean-centering in this *main* analysis
# because — combined with `knownDirection = FALSE` scoring — it suppresses
# detectable pathway signal in extreme-vs-control contrasts (e.g.,
# DisCOPD 3vs0) where gene-level DE shows thousands of FDR-significant
# genes but mean-centered singscore returns only a handful. A faithful
# original-matching reproduction with mean-centering is in the parallel
# report `bhatt_repro_singscore_meancentered_2026.Rmd`. See METHODS.md
# §4.1 and `old_vs_new_comparison.md` for the full rationale.

nf <- data.frame(norm_factor = dge$samples$norm.factors,
                 Lc.Batch    = dge$samples$Lc.Batch)
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

set.seed(1)
samp_idx <- sample(ncol(logcpm), 50)
plot_df <- data.frame(
  rand_id = rep(colnames(logcpm)[samp_idx], each = nrow(logcpm)),
  logCPM  = as.vector(logcpm[, samp_idx])
)
ggplot(plot_df, aes(x = logCPM, group = rand_id)) +
  geom_density(alpha = 0.05, color = "steelblue") +
  labs(title = sprintf("log-CPM density (random subsample of 50 / %d samples)", ncol(logcpm)),
       subtitle = "After TMM normalization, prior.count = 2",
       x = "log2(CPM + prior)", y = "Density") +
  theme_minimal(base_size = 11)

# ============================================================================
# Source:     chunk: pca
# Inputs:     logcpm (pre-batch-correction), logcpm_corrected
# Population: n_analytic subjects
# Methods:    METHODS.md §4.1 — singscore batch correction
# Note:       PCA on the top 2000 most variable AUTOSOMAL genes — X/Y/MT
#             excluded so the diagnostic reflects batch structure rather
#             than gender / mitochondrial content.
# ============================================================================
ens_bare_pca <- sub("\\.[0-9]+$", "", rownames(logcpm))
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

vars_all <- matrixStats::rowVars(logcpm)
vars_all[!autosomal_mask] <- -Inf  # exclude non-autosomal from ranking
top_var_idx <- order(vars_all, decreasing = TRUE)[1:2000]

pca_pre  <- prcomp(t(logcpm[top_var_idx, ]), center = TRUE, scale. = FALSE)
pca_post <- prcomp(t(logcpm_corrected[top_var_idx, ]), center = TRUE, scale. = FALSE)

pca_df <- function(p, label) {
  data.frame(PC1 = p$x[, 1], PC2 = p$x[, 2],
             Lc.Batch = pheno_df$Lc.Batch,
             DisCOPD  = pheno_df$DisCOPD,
             stage    = label,
             stringsAsFactors = FALSE)
}
pcd <- bind_rows(pca_df(pca_pre,  "Before batch correction"),
                 pca_df(pca_post, "After removeBatchEffect"))

p_batch <- ggplot(pcd, aes(PC1, PC2, color = Lc.Batch)) +
  geom_point(alpha = 0.5, size = 0.7) +
  facet_wrap(~ stage) +
  guides(color = "none") +
  labs(title = sprintf("PCA of log-CPM colored by Lc.Batch (n = %d, top 2000 autosomal variable genes)",
                       n_analytic),
       subtitle = sprintf("%d distinct batches; X/Y/MT excluded from gene ranking",
                          n_lc_batch_levels)) +
  theme_minimal(base_size = 11)
print(p_batch)

p_class <- ggplot(pcd, aes(PC1, PC2, color = DisCOPD)) +
  geom_point(alpha = 0.5, size = 0.7) +
  facet_wrap(~ stage) +
  scale_color_brewer(palette = "Set1", na.value = "grey80") +
  labs(title = sprintf("PCA of log-CPM colored by DisCOPD class (n = %d, autosomal only)", n_analytic),
       subtitle = "X/Y/MT excluded from gene ranking") +
  theme_minimal(base_size = 11)
print(p_class)

# ============================================================================
# Source:     chunk: gmap
# Inputs:     rownames(dge); cache/gmap_v103.rds (built from
#             AnnotationHub AH89426 = Ensembl v103 EnsDb, matching the
#             original analysis's gmap_ENSGv103_2025.6.09.rda).
# Output:     gmap (4-col data.frame: ENS versioned, ENS bare, SYMBOL, ENTREZ)
# Methods:    METHODS.md §3 — gene annotation: Ensembl v103 via
#             AnnotationHub (release matched the GENCODE v37 reference used
#             upstream in the COPDGene salmon counts pipeline).
# ============================================================================
ens_full <- rownames(dge)
ens_bare <- sub("\\.[0-9]+$", "", ens_full)
v103_rds <- file.path(CACHE_DIR, "gmap_v103.rds")
if (file.exists(v103_rds)) {
  gmap <- readRDS(v103_rds)$gmap
  # The cached v103 gmap was built against the same filtered ENSEMBL set;
  # if rownames were to differ, regenerate.
  if (!identical(gmap$ensembl_gene_id_version, ens_full)) {
    file.remove(v103_rds)  # force rebuild
    gmap <- NULL
  }
} else gmap <- NULL

if (is.null(gmap)) {
  suppressPackageStartupMessages({
    library(AnnotationHub); library(ensembldb)
  })
  ah <- AnnotationHub()
  edb <- ah[["AH89426"]]   # Ensembl v103 EnsDb for Homo sapiens
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

# ============================================================================
# Source:     chunk: pathways
# Inputs:     msigdbr (collections C2:CP:REACTOME, C2:CP:KEGG_MEDICUS,
#                       C7:IMMUNESIGDB, C2:CP:BIOCARTA)
# Output:     pathway_list (named list of gene-symbol vectors), pathway_source
# Methods:    METHODS.md §5
# ============================================================================
pulls <- list(
  Reactome    = list(coll = "C2", sub = "CP:REACTOME"),
  KEGG        = list(coll = "C2", sub = "CP:KEGG_MEDICUS"),
  ImmuneSigDB = list(coll = "C7", sub = "IMMUNESIGDB"),
  BioCarta    = list(coll = "C2", sub = "CP:BIOCARTA")
)
all_df <- bind_rows(lapply(names(pulls), function(src) {
  cfg <- pulls[[src]]
  df  <- msigdbr(species = "Homo sapiens",
                 collection = cfg$coll, subcollection = cfg$sub)
  data.frame(gs_name = df$gs_name, gene_symbol = df$gene_symbol,
             source  = src, stringsAsFactors = FALSE)
}))
all_df <- all_df[!duplicated(all_df[, c("gs_name","gene_symbol")]), ]
pathway_list   <- split(all_df$gene_symbol, all_df$gs_name)
pathway_source <- setNames(all_df$source[!duplicated(all_df$gs_name)],
                           all_df$gs_name[!duplicated(all_df$gs_name)])
n_pathways_loaded <- length(pathway_list)

# ============================================================================
# Source:     chunk: score
# Inputs:     logcpm_corrected (batch-corrected log-CPM, gene SYMBOL keyed),
#             pathway_list
# Population: n_analytic samples × n_genes_kept genes
# Output:     score_mat (pathways × samples)
# Methods:    METHODS.md §5 — singscore knownDirection=FALSE, vectorized
# ============================================================================

# ENSEMBL → SYMBOL on logcpm_corrected; dedup by max mean expression
sym_of <- gmap$external_gene_name[match(rownames(logcpm_corrected),
                                        gmap$ensembl_gene_id_version)]
mask   <- !is.na(sym_of) & sym_of != ""
expr   <- logcpm_corrected[mask, , drop = FALSE]
sym_of <- sym_of[mask]
ord <- order(sym_of, -rowMeans(expr))
expr <- expr[ord, , drop = FALSE]; sym_of <- sym_of[ord]
expr_sym <- expr[!duplicated(sym_of), , drop = FALSE]
rownames(expr_sym) <- sym_of[!duplicated(sym_of)]
n_genes_with_symbol_dedup <- nrow(expr_sym)

rankData <- rankGenes(expr_sym)
avail <- rownames(rankData)
# Two-part pathway filter (METHODS.md §5):
#   (a) 15-500 overlap with available genes (standard MSigDB practice)
#   (b) ≥50% of the pathway's annotated genes must be present in the
#       expression data, ensuring the score is a faithful representation
#       of pathway activity rather than a fragment.
sizes_overlap <- vapply(pathway_list, function(g) length(intersect(g, avail)),
                        integer(1))
sizes_total   <- lengths(pathway_list)
frac_present  <- sizes_overlap / sizes_total
pass_size     <- sizes_overlap >= 15 & sizes_overlap <= 500
pass_coverage <- frac_present >= 0.5
keep_pw       <- names(sizes_overlap)[pass_size & pass_coverage]
pathway_list_filt <- lapply(pathway_list[keep_pw], function(g) intersect(g, avail))
n_pathways_kept       <- length(pathway_list_filt)
n_pathways_drop_size  <- sum(!pass_size)
n_pathways_drop_cover <- sum(pass_size & !pass_coverage)

# Per-source × per-filter-step cascade for the diagnostic table below.
# Important: `pathway_source` is named in msigdbr-load order (Reactome,
# KEGG, ImmuneSigDB, BioCarta), but `pass_size` / `pass_coverage` are
# named over `pathway_list` which `split()` returns alphabetically. We
# must align by pathway name (not positional binding) before tabulating
# per-source counts; otherwise sources mix with the wrong filter results.
stopifnot(setequal(names(pathway_source), names(pass_size)))
pw_cascade <- data.frame(
  Source         = pathway_source[names(pass_size)],
  pass_size      = pass_size,
  pass_size_cov  = pass_size & pass_coverage,
  stringsAsFactors = FALSE
)
pw_cascade_tab <- as.data.frame(table(Source = pw_cascade$Source))
colnames(pw_cascade_tab) <- c("Source", "Loaded")
pw_cascade_tab$Pass_size_15_500 <- as.integer(table(pw_cascade$Source[pw_cascade$pass_size]))[
  match(pw_cascade_tab$Source, names(table(pw_cascade$Source[pw_cascade$pass_size])))]
pw_cascade_tab$Pass_size_15_500[is.na(pw_cascade_tab$Pass_size_15_500)] <- 0L
pw_cascade_tab$Pass_coverage_50pct <- as.integer(table(pw_cascade$Source[pw_cascade$pass_size_cov]))[
  match(pw_cascade_tab$Source, names(table(pw_cascade$Source[pw_cascade$pass_size_cov])))]
pw_cascade_tab$Pass_coverage_50pct[is.na(pw_cascade_tab$Pass_coverage_50pct)] <- 0L
pw_cascade_tab <- rbind(pw_cascade_tab,
  data.frame(Source = "TOTAL",
             Loaded              = sum(pw_cascade_tab$Loaded),
             Pass_size_15_500    = sum(pw_cascade_tab$Pass_size_15_500),
             Pass_coverage_50pct = sum(pw_cascade_tab$Pass_coverage_50pct)))

# Vectorized re-implementation of singscore::multiScore for
# knownDirection=FALSE (METHODS.md §5)
multi_score_undirected <- function(rankData, pathway_list_filt) {
  G <- nrow(rankData); N <- ncol(rankData)
  col_med  <- ceiling(matrixStats::colMedians(rankData))
  folded   <- abs(rankData - rep(col_med, each = G))
  gene_idx <- setNames(seq_len(G), rownames(rankData))
  pw_names <- names(pathway_list_filt); P <- length(pw_names)
  triplet_i <- integer(sum(lengths(pathway_list_filt)))
  triplet_j <- integer(length(triplet_i))
  pos <- 1L; sizes_v <- integer(P)
  for (pi in seq_along(pathway_list_filt)) {
    g_ids <- gene_idx[pathway_list_filt[[pi]]]
    k <- length(g_ids); sizes_v[pi] <- k
    triplet_i[pos:(pos + k - 1L)] <- pi
    triplet_j[pos:(pos + k - 1L)] <- g_ids
    pos <- pos + k
  }
  membership <- sparseMatrix(i = triplet_i, j = triplet_j,
                             x = 1, dims = c(P, G))
  rownames(membership) <- pw_names
  pw_sum  <- as.matrix(membership %*% folded)
  pw_mean <- pw_sum / sizes_v
  k_half <- floor(sizes_v / 2); bg_half <- ceiling(G / 2)
  low <- (k_half + 1) / 2
  up  <- (2 * bg_half - k_half + 1) / 2
  scores <- (pw_mean - low) / (up - low)
  rownames(scores) <- pw_names
  colnames(scores) <- colnames(rankData)
  scores
}
score_mat <- multi_score_undirected(rankData, pathway_list_filt)
stopifnot(identical(colnames(score_mat), pheno_df$rand_id))

# In-Rmd validation: score 3 random pathways with singscore::simpleScore and
# verify our vectorized output matches bit-exact. Breaks the knit if the
# vectorized formula ever drifts from the canonical singscore implementation.
# Note: `centerScore = TRUE` is silently a no-op for `knownDirection = FALSE`
#       in singscore (see singscore source: `multiSingscore`, the
#       `if (!knownDirection)` branch sets `center_offset <- 0` with a
#       warning). Our vectorized impl also subtracts no offset, so the
#       comparison is bit-exact under matching parameters.
set.seed(20260425)
.probe_pw <- sample(names(pathway_list_filt), 3)
.simple <- suppressWarnings(vapply(.probe_pw, function(nm) {
  gs <- GeneSet(pathway_list_filt[[nm]], setName = nm)
  simpleScore(rankData, upSet = gs, knownDirection = FALSE,
              centerScore = TRUE)$TotalScore
}, numeric(ncol(rankData))))
.max_diff <- max(abs(t(.simple) - score_mat[.probe_pw, , drop = FALSE]))
stopifnot(.max_diff < 1e-10)

# Diagnostic: pathway-score distribution per sample (subset)
set.seed(1)
samp_idx <- sample(ncol(score_mat), 50)
sd_df <- data.frame(
  sample = rep(colnames(score_mat)[samp_idx], each = nrow(score_mat)),
  score  = as.vector(score_mat[, samp_idx])
)
ggplot(sd_df, aes(x = score, group = sample)) +
  geom_density(alpha = 0.05, color = "purple") +
  labs(title = sprintf("Pathway score density per sample (subsample of 50 / %d)", ncol(score_mat)),
       subtitle = sprintf("%d pathways scored", nrow(score_mat)),
       x = "Singscore (centered, normalized)", y = "Density") +
  theme_minimal(base_size = 11)

# ============================================================================
# Source:     chunk: regress
# Inputs:     score_mat (pathways × samples), pheno_df
# Population: n_analytic subjects
# Output:     all_pw (per-cell topTables), pw_summary
# Methods:    METHODS.md §4.1 (Lc.Batch in design for score regression) +
#             §8 (statistical testing — eBayes, BH FDR per cell)
# ============================================================================

# ---- 1. Covariate strings for the two models ----
# Each string is the right-hand side of the design formula, EXCLUDING the
# classification variable (which is appended per cell in the inner loop).
# Lc.Batch is included as a fixed-effect covariate even though batch was
# also removed from the expression matrix upstream — defense in depth so
# any residual batch-correlated score variance is absorbed here too.
cov_m1 <- "Lc.Batch + age_visit + race + gender + SmokCigNow"
cov_m2 <- paste(cov_m1, "+ eosinphl_pct + neutrophl_pct + lymphcyt_pct + monocyt_pct")
model_configs <- list(Model_1 = cov_m1, Model_2 = cov_m2)

# ---- 2. Per-cell loop ----
# Two outer loops (model × classification) build the fit; an inner loop
# extracts each of the 6 contrasts from that fit. Per-cell results are
# accumulated in two named lists keyed by "Model|Classification|Contrast".
all_pw <- list(); all_pw_sig <- list()

for (model_label in names(model_configs)) {
  cov <- model_configs[[model_label]]

  for (cls_var in cls_vars) {

    # ---- 2a. Subset to subjects with non-NA classification ----
    # Complete-case upstream guarantees this is a no-op (all 3 classifications
    # are non-NA in the analytic cohort), but the guard is kept for safety
    # in case a future cohort has subjects missing one classification but
    # not the others.
    keep_s <- !is.na(pheno_df[[cls_var]])
    ph_s   <- droplevels(pheno_df[keep_s, , drop = FALSE])
    scr_s  <- score_mat[, keep_s, drop = FALSE]

    # ---- 2b. Build the design matrix ----
    # Formula = ~ <model covariates> + <classification factor>. R's default
    # treatment-coding for the factor produces dummy columns named
    # `<cls_var>1`, `<cls_var>2`, `<cls_var>3` (level 0 is the reference).
    # `make.names` ensures column names are syntactically valid R identifiers
    # (e.g., factor level "2" becomes column "DisCOPD2", not "DisCOPD 2").
    design_mat <- model.matrix(as.formula(paste("~", cov, "+", cls_var)),
                               data = ph_s)
    colnames(design_mat) <- make.names(colnames(design_mat))

    # ---- 2c. Build the contrast matrix ----
    # We want six contrasts: the three "level vs reference 0" coefficients
    # taken directly, plus three pairwise differences (1 vs 3, 2 vs 3,
    # 1 vs 2) constructed by subtraction. Each contrast is added only if
    # the involved levels are actually present in this cohort.
    lvls    <- levels(ph_s[[cls_var]])
    present <- intersect(c("1","2","3"), lvls)
    cs <- list()
    for (k in present) {
      # "Grp1vs0" -> the dummy column for level 1 vs reference (level 0)
      cs[[paste0("Grp", k, "vs0")]] <- paste0(cls_var, k)
    }
    if (all(c("1","3") %in% present)) cs[["Grp1vs3"]] <- paste0(cls_var,"1 - ", cls_var,"3")
    if (all(c("2","3") %in% present)) cs[["Grp2vs3"]] <- paste0(cls_var,"2 - ", cls_var,"3")
    if (all(c("1","2") %in% present)) cs[["Grp1vs2"]] <- paste0(cls_var,"1 - ", cls_var,"2")
    contrast_mtx <- do.call(makeContrasts, c(cs, list(levels = design_mat)))
    # Rename columns from "Grp1vs0" → "DisCOPD_1vs0" etc., so the per-cell
    # downstream key uniquely identifies model + classification + contrast.
    colnames(contrast_mtx) <- paste0(cls_var, "_", sub("^Grp","", names(cs)))

    # ---- 2d. Fit + apply contrasts + empirical Bayes ----
    # `lmFit(scores, design)` runs ordinary least squares per pathway
    # (each row of scr_s = one pathway across samples). `contrasts.fit`
    # rotates the per-coefficient SEs to the per-contrast SEs without
    # refitting. `eBayes` then borrows variance information across
    # pathways to shrink per-pathway residual variances toward a common
    # trend — this is what makes a small per-pathway sample more powerful
    # than a naive t-test would be.
    fit  <- lmFit(scr_s, design_mat)
    fit2 <- eBayes(contrasts.fit(fit, contrast_mtx))

    # ---- 2e. Per-contrast extraction ----
    # `topTable(coef = cn, sort.by = "none")` returns one row per pathway
    # with logFC, t, P.Value, adj.P.Val (BH-adjusted within this contrast),
    # AveExpr, B, plus 95% CI bounds (because confint = TRUE). We attach
    # the pathway name and source collection for downstream reporting,
    # round display columns to 3 sig figs, and re-sort by FDR for table
    # readability. P-value counts are then accumulated for the summary.
    for (cn in colnames(contrast_mtx)) {
      key <- paste(model_label, cls_var, cn, sep = "|")
      tt  <- topTable(fit2, coef = cn, number = Inf, confint = TRUE,
                      sort.by = "none")
      tt$pathway <- rownames(tt)
      tt$source  <- pathway_source[tt$pathway]
      tt <- as_tibble(tt) |>
        mutate(across(c(logFC, P.Value, adj.P.Val, AveExpr, CI.L, CI.R),
                      ~ signif(as.numeric(.x), 3))) |>
        arrange(adj.P.Val) |>
        dplyr::select(pathway, source, logFC, P.Value, adj.P.Val, AveExpr, CI.L, CI.R)
      all_pw[[key]] <- tt
      all_pw_sig[[key]] <- data.frame(
        Model = model_label, Classification = cls_var, Contrast = cn,
        N_samples = nrow(ph_s),
        PW_Nominal_05 = sum(tt$P.Value   < 0.05, na.rm = TRUE),
        PW_FDR_05     = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
        stringsAsFactors = FALSE)
    }
  }
}

# ---- 3. Collapse the per-cell summaries into a single data.frame ----
# `pw_summary` has one row per cell; the helper column `short_contrast`
# strips the "<cls_var>_" prefix so the bar plots can use a clean x-axis.
pw_summary <- bind_rows(all_pw_sig) |>
  mutate(short_contrast = sub("^[^_]+_", "", Contrast))

DT::datatable(pw_summary[, c("Model","Classification","Contrast",
                             "N_samples","PW_FDR_05")],
              caption = sprintf("Singscore pathway significance — %d cells (Model × Classification × Contrast)",
                                nrow(pw_summary)),
              options = list(pageLength = 36, scrollX = TRUE),
              rownames = FALSE, filter = "top")

ggplot(pw_summary, aes(x = short_contrast, y = PW_FDR_05, fill = Classification)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.75) +
  facet_wrap(~ Model, ncol = 1) +
  labs(title = sprintf("FDR-significant singscore pathways per cell (n = %d subjects, %d pathways tested)",
                       n_analytic, n_pathways_kept),
       subtitle = "adj.P.Val < 0.05",
       x = "Contrast", y = "Number of pathways") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")

src_breakdown <- bind_rows(lapply(names(all_pw), function(key) {
  parts <- strsplit(key, "[|]")[[1]]
  sig <- all_pw[[key]] |> dplyr::filter(adj.P.Val < 0.05)
  if (!nrow(sig)) return(NULL)
  sig |> dplyr::count(source, name = "n_FDR_sig") |>
    mutate(Model = parts[1], Classification = parts[2], Contrast = parts[3])
}))
if (nrow(src_breakdown)) {
  DT::datatable(src_breakdown[, c("Model","Classification","Contrast",
                                  "source","n_FDR_sig")],
                caption = "FDR-significant pathways by collection",
                options = list(pageLength = 25, scrollX = TRUE),
                rownames = FALSE, filter = "top")
} else {
  cat("(Data) No FDR-significant pathways across any cell.\n")
}

sig_sets <- lapply(cls_vars, function(v) {
  pat <- paste0("^Model_1[|]", v, "[|]")
  ks  <- grep(pat, names(all_pw), value = TRUE)
  unique(unlist(lapply(ks, function(k) {
    x <- all_pw[[k]]; x$pathway[x$adj.P.Val < 0.05]
  })))
})
names(sig_sets) <- cls_vars
mat_overlap <- outer(seq_along(cls_vars), seq_along(cls_vars),
             Vectorize(function(i, j) length(intersect(sig_sets[[i]], sig_sets[[j]]))))
dimnames(mat_overlap) <- list(cls_vars, cls_vars)
kable(mat_overlap, caption = "Pairwise overlap of FDR-significant pathways (Model 1, union over the 6 contrasts per classification). Diagonal = total set size for that classification (any contrast); off-diagonals = pathway intersections, counted across any contrast pair (not within-matched-contrast).")

lfc_lookup <- lapply(cls_vars, function(v) {
  k <- paste0("Model_1|", v, "|", v, "_1vs0")
  setNames(all_pw[[k]]$logFC, all_pw[[k]]$pathway)
})
names(lfc_lookup) <- cls_vars
pairs_l <- combn(cls_vars, 2, simplify = FALSE)
plots <- lapply(pairs_l, function(pr) {
  paths <- intersect(names(lfc_lookup[[pr[1]]]), names(lfc_lookup[[pr[2]]]))
  d <- data.frame(x = lfc_lookup[[pr[1]]][paths],
                  y = lfc_lookup[[pr[2]]][paths])
  r <- cor(d$x, d$y, use = "complete.obs")
  ggplot(d, aes(x, y)) +
    geom_point(alpha = 0.18, size = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = sprintf("%s vs %s (n = %d pathways)", pr[1], pr[2], nrow(d)),
         subtitle = sprintf("Pearson r = %.3f", r),
         x = paste("logFC", pr[1]), y = paste("logFC", pr[2])) +
    theme_minimal(base_size = 10)
})
wrap_plots(plots, ncol = 3) +
  plot_annotation(title = sprintf("Singscore logFC concordance — Model 1, contrast 1vs0 (n = %d subjects)",
                                  n_analytic))

for (model_label in unique(pw_summary$Model)) {
  cat("\n\n## ", model_label, " {.tabset}\n\n", sep = "")
  for (v in cls_vars) {
    cat("\n\n### ", cls_labels[v], " (`", v, "`) {.tabset}\n\n", sep = "")
    pat <- paste0("^", model_label, "[|]", v, "[|]")
    keys <- grep(pat, names(all_pw), value = TRUE)
    for (k in keys) {
      cn <- strsplit(k, "[|]")[[1]][3]
      cat("\n\n#### ", cn, "\n\n", sep = "")
      tt <- all_pw[[k]]
      print(kable(head(tt, 20), digits = 4,
                  caption = sprintf("Top 20 pathways — %s | %s | %s (n = %d)",
                                    model_label, v, cn, n_analytic)))
      cat("\n")
    }
  }
}

# ============================================================================
# Source:     chunk: exports
# Inputs:     all_pw, pw_summary, score_mat
# Output:     CSVs in results/, RDS in cache/
# ============================================================================
fwrite(pw_summary[, c("Model","Classification","Contrast",
                      "N_samples","PW_Nominal_05","PW_FDR_05")],
       file.path(RES_DIR, paste0("singscore_summary_", DATE_STAMP, ".csv")))
for (model_label in names(model_configs)) {
  for (cls_var in cls_vars) {
    pat  <- paste0("^", model_label, "[|]", cls_var, "[|]")
    keys <- grep(pat, names(all_pw), value = TRUE)
    if (!length(keys)) next
    combined <- bind_rows(lapply(keys, function(k) {
      parts <- strsplit(k, "[|]")[[1]]
      all_pw[[k]] |> mutate(Model = parts[1], Contrast = parts[3])
    }))
    fwrite(combined, file.path(RES_DIR,
            sprintf("singscore_%s_%s_allresults_%s.csv",
                    model_label, cls_var, DATE_STAMP)))
  }
}
# Score matrix export (samples × pathways) — supports downstream
# per-sample analyses. First column is `rand_id` (the deidentified
# sample key used throughout this project); subsequent columns are
# pathway names. Original analysis's score matrix used `sid` as the
# column name with non-deidentified IDs; we use `rand_id` to match the
# deidentified-data convention.
score_export <- as.data.frame(t(score_mat))
score_export$rand_id <- rownames(score_export)
score_export <- score_export[, c("rand_id",
                                 setdiff(colnames(score_export), "rand_id"))]
fwrite(score_export,
       file.path(RES_DIR, paste0("singscore_score_matrix_", DATE_STAMP, ".csv")))

saveRDS(list(all_pw = all_pw, pw_summary = pw_summary,
             pathway_source = pathway_source, pheno = pheno_df,
             cls_vars = cls_vars, cls_labels = cls_labels,
             score_mat = score_mat, pathway_list_filt = pathway_list_filt,
             date_stamp = DATE_STAMP),
        file.path(CACHE_DIR, "singscore_results.rds"))

# Standalone RDS of the analyzed pathway library — convenient for
# downstream consumers who want the gene-symbol membership lists without
# loading the full results bundle.
analyzed_pathways <- list(
  pathway_list    = pathway_list_filt,
  pathway_source  = pathway_source[names(pathway_list_filt)],
  n_pathways      = length(pathway_list_filt),
  collections     = c("Reactome (C2:CP:REACTOME)",
                      "KEGG Medicus (C2:CP:KEGG_MEDICUS)",
                      "ImmuneSigDB (C7:IMMUNESIGDB)",
                      "BioCarta (C2:CP:BIOCARTA)"),
  msigdbr_version = as.character(packageVersion("msigdbr")),
  filter_rule     = paste("15 <= overlap <= 500 AND >= 50% of pathway",
                          "annotated genes present in expression matrix"),
  date_stamp      = DATE_STAMP
)
saveRDS(analyzed_pathways,
        file.path(RES_DIR, paste0("singscore_analyzed_pathways_",
                                  DATE_STAMP, ".rds")))