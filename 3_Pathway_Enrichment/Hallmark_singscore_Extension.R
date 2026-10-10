#!/usr/bin/env Rscript

# Compute the missing Hallmark participant-level singscores and refit the
# complete 6,045-set Model 1 DisCOPDAI analysis. This extends the same analysis;
# it does not substitute GSEA statistics for single-sample scores.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(Matrix)
  library(matrixStats)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]))
root <- dirname(dirname(script_file))
table_dir <- file.path(root, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = "") {
  hit <- match(flag, args)
  if (is.na(hit)) return(default)
  if (hit == length(args)) stop(flag, " requires a value")
  args[[hit + 1L]]
}

data_dir <- arg_value("--data-dir", Sys.getenv("COPDGENE_DEID_DIR", ""))
if (!nzchar(data_dir)) {
  data_dir <- "path/to/COPDGene/deidentified"
}

required_files <- c(
  counts = "counts_raw_freeze5_randid.csv",
  master = "masterfile_rna_freeze5_randid.csv",
  pheno = "COPDGene_P1P2P3_SM_NS_ILDBR_Long_Sep24_randid.csv",
  original_class = "MultidimensionalSchema_deid.csv",
  ai_class = "NewClassQualitativeAI_deid.csv",
  quantitative_class = "NewClassQuantitative_deid.csv"
)
input_paths <- file.path(data_dir, required_files)
missing <- input_paths[!file.exists(input_paths)]
if (length(missing)) {
  stop(
    "Hallmark singscore requires the original deidentified participant-level inputs.\n",
    "Set --data-dir or COPDGENE_DEID_DIR to a directory containing:\n  ",
    paste(required_files, collapse = "\n  "),
    "\nMissing under ", data_dir, ":\n  ", paste(basename(missing), collapse = "\n  ")
  )
}

schema_root <- paste0(
  "path/to/",
  "Schema_omics_csv_bhatt"
)
de_file <- file.path(
  schema_root, "copdschema/extracted/",
  "gene_de_Model_1_DisCOPDAI_allresults_2026.04.25.csv"
)
msigdb_file <- file.path(
  schema_root, "reproduce_2026-09-01/outputs/tables/msigdb_cache.rds"
)
singscore_zip <- file.path(schema_root, "copdschema/singscore_archive.zip")
singscore_member <- "singscore_score_matrix_2026.04.26.csv"
singscore_pathways_member <- "singscore_analyzed_pathways_2026.04.26.rds"
existing_results_file <- file.path(
  schema_root, "reanalysis_2026/tables/singscore_differential_all.csv"
)

stopifnot(
  file.exists(de_file), file.exists(msigdb_file),
  file.exists(singscore_zip), file.exists(existing_results_file)
)

cat("Reading deidentified count matrix...\n")
counts_dt <- fread(input_paths[["counts"]], header = TRUE)
gene_id_column <- if ("GENEID" %chin% names(counts_dt)) "GENEID" else names(counts_dt)[[1]]
gene_ids <- counts_dt[[gene_id_column]]
counts_dt[, (gene_id_column) := NULL]
counts <- as.matrix(counts_dt)
rm(counts_dt)
storage.mode(counts) <- "double"
rownames(counts) <- gene_ids

cat("Reconstructing the original 3,309-participant analytic cohort...\n")
master <- fread(input_paths[["master"]])[final.analyzed == 1, .(rand_id, Lc.Batch)]
pheno <- fread(
  input_paths[["pheno"]],
  select = c(
    "rand_id", "visitnum", "age_visit", "race", "gender", "SmokCigNow",
    "eosinphl_pct", "neutrophl_pct", "lymphcyt_pct", "monocyt_pct"
  )
)[visitnum == 2]
original_class <- fread(input_paths[["original_class"]])[, .(rand_id, DisCOPD)]
ai_class <- fread(input_paths[["ai_class"]])[, .(rand_id, DisCOPDAI)]
quantitative_class <- fread(input_paths[["quantitative_class"]])[, .(rand_id, DisCOPDQuant)]

for (object in list(master, pheno, original_class, ai_class, quantitative_class)) {
  object[, rand_id := as.character(rand_id)]
}
pheno <- Reduce(
  function(x, y) merge(x, y, by = "rand_id", all = FALSE),
  list(pheno, master, original_class, ai_class, quantitative_class)
)
complete_columns <- c(
  "DisCOPD", "DisCOPDAI", "DisCOPDQuant",
  "SmokCigNow", "eosinphl_pct", "neutrophl_pct", "lymphcyt_pct", "monocyt_pct"
)
pheno <- pheno[complete.cases(pheno[, ..complete_columns])]
pheno <- pheno[rand_id %chin% colnames(counts)]
pheno[, sample_order := match(rand_id, colnames(counts))]
setorder(pheno, sample_order)
pheno[, sample_order := NULL]
counts <- counts[, pheno$rand_id, drop = FALSE]
stopifnot(identical(colnames(counts), pheno$rand_id), nrow(pheno) == 3309L)

pheno[, `:=`(
  age_visit = as.numeric(age_visit),
  race = factor(race),
  gender = factor(gender),
  SmokCigNow = factor(SmokCigNow),
  Lc.Batch = factor(Lc.Batch),
  DisCOPD = factor(DisCOPD, levels = c("0", "1", "2", "3")),
  DisCOPDAI = factor(DisCOPDAI, levels = c("0", "1", "2", "3")),
  DisCOPDQuant = factor(DisCOPDQuant, levels = c("0", "1", "2", "3"))
)]
pheno_df <- as.data.frame(pheno)
rownames(pheno_df) <- pheno_df$rand_id

cat("Applying the original filter, TMM normalization, and batch correction...\n")
dge <- DGEList(
  counts = counts,
  samples = pheno_df,
  genes = data.frame(id = rownames(counts), stringsAsFactors = FALSE)
)
rm(counts)
filter_design <- model.matrix(
  ~ age_visit + race + gender + SmokCigNow + Lc.Batch,
  data = dge$samples
)
keep <- filterByExpr(dge, design = filter_design, min.count = 50, min.prop = 0.9)
dge <- dge[keep, , keep.lib.sizes = TRUE]
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 2)
logcpm_corrected <- removeBatchEffect(logcpm, batch = dge$samples$Lc.Batch)
stopifnot(nrow(logcpm_corrected) == 16554L)
filtered_gene_n <- nrow(logcpm_corrected)

cat("Mapping Ensembl identifiers to the exported Model 1 gene symbols...\n")
de_map <- unique(
  fread(de_file, select = c("id", "gene_label"))[, .(id, gene_label)]
)
symbol <- de_map$gene_label[match(rownames(logcpm_corrected), de_map$id)]
if (mean(!is.na(symbol) & nzchar(symbol)) < 0.9) {
  symbol <- de_map$gene_label[
    match(sub("\\.[0-9]+$", "", rownames(logcpm_corrected)), sub("\\.[0-9]+$", "", de_map$id))
  ]
}
mapped <- !is.na(symbol) & nzchar(symbol)
if (mean(mapped) < 0.9) stop("Fewer than 90% of filtered genes mapped to symbols")
expression <- logcpm_corrected[mapped, , drop = FALSE]
symbol <- symbol[mapped]
ordering <- order(symbol, -rowMeans(expression))
expression <- expression[ordering, , drop = FALSE]
symbol <- symbol[ordering]
expression <- expression[!duplicated(symbol), , drop = FALSE]
rownames(expression) <- symbol[!duplicated(symbol)]
rm(logcpm, logcpm_corrected, dge)

cat("Loading the 50 Hallmark sets used by the current pathway workflow...\n")
msigdb <- readRDS(msigdb_file)
hallmark_names <- names(msigdb$sets)[msigdb$source[names(msigdb$sets)] == "Hallmark"]
hallmark <- msigdb$sets[hallmark_names]
available <- rownames(expression)
overlap_n <- vapply(hallmark, function(x) length(intersect(x, available)), integer(1))
coverage <- overlap_n / lengths(hallmark)
keep_hallmark <- overlap_n >= 15 & overlap_n <= 500 & coverage >= 0.5
if (!all(keep_hallmark)) {
  stop("One or more Hallmark sets failed the prespecified overlap/coverage criteria")
}
hallmark <- lapply(hallmark, intersect, y = available)
stopifnot(length(hallmark) == 50L)

cat("Ranking genes within each participant and computing undirected singscores...\n")
rank_data <- apply(expression, 2, rank, ties.method = "min")
rownames(rank_data) <- rownames(expression)
colnames(rank_data) <- colnames(expression)
rm(expression)

multi_score_undirected <- function(rankData, pathways) {
  genes_n <- nrow(rankData)
  sample_median <- ceiling(matrixStats::colMedians(rankData))
  folded <- abs(rankData - rep(sample_median, each = genes_n))
  gene_index <- setNames(seq_len(genes_n), rownames(rankData))
  pathway_names <- names(pathways)
  pathway_sizes <- lengths(pathways)
  row_index <- rep(seq_along(pathways), pathway_sizes)
  column_index <- unname(unlist(lapply(pathways, function(x) gene_index[x]), use.names = FALSE))
  membership <- sparseMatrix(
    i = row_index, j = column_index, x = 1,
    dims = c(length(pathways), genes_n)
  )
  rownames(membership) <- pathway_names
  pathway_mean <- as.matrix(membership %*% folded) / pathway_sizes
  pathway_half <- floor(pathway_sizes / 2)
  background_half <- ceiling(genes_n / 2)
  lower <- (pathway_half + 1) / 2
  upper <- (2 * background_half - pathway_half + 1) / 2
  scores <- (pathway_mean - lower) / (upper - lower)
  rownames(scores) <- pathway_names
  colnames(scores) <- colnames(rankData)
  scores
}
hallmark_scores <- multi_score_undirected(rank_data, hallmark)
stopifnot(nrow(hallmark_scores) == 50L, ncol(hallmark_scores) == 3309L)

cat("Loading archived participant-level scores and refitting all 6,045 sets together...\n")
temporary_directory <- tempfile("singscore_archive_")
dir.create(temporary_directory)
on.exit(unlink(temporary_directory, recursive = TRUE), add = TRUE)
unzip(
  singscore_zip,
  files = c(singscore_member, singscore_pathways_member),
  exdir = temporary_directory
)
existing_scores_dt <- fread(file.path(temporary_directory, singscore_member))
sample_column <- if ("rand_id" %chin% names(existing_scores_dt)) "rand_id" else names(existing_scores_dt)[[1]]
existing_ids <- as.character(existing_scores_dt[[sample_column]])
existing_scores_dt[, (sample_column) := NULL]
existing_scores <- t(as.matrix(existing_scores_dt))
rm(existing_scores_dt)
storage.mode(existing_scores) <- "double"
colnames(existing_scores) <- existing_ids
sample_match <- match(pheno$rand_id, colnames(existing_scores))
stopifnot(!anyNA(sample_match))
existing_scores <- existing_scores[, sample_match, drop = FALSE]
stopifnot(identical(colnames(existing_scores), colnames(hallmark_scores)))

archived_pathways <- readRDS(
  file.path(temporary_directory, singscore_pathways_member)
)
validation_names <- unlist(lapply(
  unique(archived_pathways$pathway_source),
  function(collection) {
    candidates <- names(archived_pathways$pathway_source)[
      archived_pathways$pathway_source == collection
    ]
    intersect(candidates, rownames(existing_scores))[[1]]
  }
), use.names = FALSE)
validation_sets <- lapply(
  archived_pathways$pathway_list[validation_names], intersect, y = available
)
validation_scores <- multi_score_undirected(rank_data, validation_sets)
validation_error <- max(abs(
  validation_scores - existing_scores[validation_names, , drop = FALSE]
))
if (!is.finite(validation_error) || validation_error > 1e-8) {
  stop(
    "Source-equivalence check failed against archived singscores; maximum ",
    "absolute difference = ", signif(validation_error, 5)
  )
}
cat(
  "Archived singscore equivalence check passed; maximum absolute difference = ",
  signif(validation_error, 4), "\n", sep = ""
)
rm(rank_data, archived_pathways, validation_scores)

score_matrix <- rbind(existing_scores, hallmark_scores)
stopifnot(
  !anyDuplicated(rownames(score_matrix)),
  nrow(existing_scores) == 5995L,
  nrow(score_matrix) == 6045L
)

# Lc.Batch was removed from log-CPM before scoring and is not entered twice.
design <- model.matrix(
  ~ age_visit + race + gender + SmokCigNow + DisCOPDAI,
  data = pheno_df
)
colnames(design) <- make.names(colnames(design))
contrast_matrix <- makeContrasts(
  DisCOPDAI_1vs0 = DisCOPDAI1,
  DisCOPDAI_2vs0 = DisCOPDAI2,
  DisCOPDAI_3vs0 = DisCOPDAI3,
  levels = design
)
fit <- eBayes(contrasts.fit(lmFit(score_matrix, design), contrast_matrix))

source_map <- unique(
  fread(existing_results_file, select = c("pathway", "source"))[, .(pathway, source)]
)
source_lookup <- setNames(source_map$source, source_map$pathway)
source_lookup[hallmark_names] <- "Hallmark"

results <- rbindlist(lapply(colnames(contrast_matrix), function(contrast) {
  table <- as.data.table(
    topTable(fit, coef = contrast, number = Inf, confint = TRUE, sort.by = "none"),
    keep.rownames = "pathway"
  )
  table[, `:=`(
    source = unname(source_lookup[pathway]),
    Model = "Model_1",
    Classification = "DisCOPDAI",
    Contrast = contrast,
    contrast_label = c(
      DisCOPDAI_1vs0 = "Newly diagnosed vs no COPD",
      DisCOPDAI_2vs0 = "Potentially excluded vs no COPD",
      DisCOPDAI_3vs0 = "Established COPD vs no COPD"
    )[[contrast]],
    score_diff = logFC,
    pval = P.Value,
    FDR = adj.P.Val,
    direction = fifelse(logFC >= 0, "up", "down"),
    sig = adj.P.Val < 0.05
  )]
  table[]
}))
setcolorder(
  results,
  c(
    "pathway", "source", "Model", "Classification", "Contrast",
    "contrast_label", "score_diff", "pval", "FDR", "direction", "sig",
    "AveExpr", "t", "B", "CI.L", "CI.R"
  )
)
setorder(results, Contrast, FDR, pathway)

fwrite(results, file.path(table_dir, "singscore_differential_hallmark_extended_model1.csv"))
fwrite(
  results[source == "Hallmark"],
  file.path(table_dir, "singscore_hallmark_model1_all.csv")
)
saveRDS(
  list(
    scores = hallmark_scores,
    pathway_membership = hallmark,
    pathway_overlap = overlap_n,
    pathway_coverage = coverage,
    sample_ids = colnames(hallmark_scores)
  ),
  file.path(table_dir, "singscore_hallmark_score_matrix_model1.rds")
)

provenance <- data.table(
  field = c(
    "analysis", "participants", "filtered_genes", "mapped_unique_symbols",
    "existing_pathways", "hallmark_pathways", "combined_test_family",
    "reference", "model", "rank_ties", "known_direction",
    "archived_score_max_abs_error", "gene_set_cache", "msigdbr_version"
  ),
  value = c(
    "Extension of archived singscore analysis", nrow(pheno), filtered_gene_n,
    length(available), nrow(existing_scores), nrow(hallmark_scores), nrow(score_matrix),
    "No COPD (DisCOPDAI level 0)",
    "age_visit + race + gender + SmokCigNow + DisCOPDAI; Lc.Batch removed upstream",
    "minimum rank", "FALSE (undirected)", signif(validation_error, 8),
    msigdb_file, as.character(msigdb$msigdbr_version)
  )
)
fwrite(provenance, file.path(table_dir, "singscore_hallmark_provenance.csv"))

summary <- results[, .(
  pathways = .N,
  fdr_05 = sum(FDR < 0.05),
  hallmark_fdr_05 = sum(source == "Hallmark" & FDR < 0.05)
), by = Contrast]
print(summary)
cat("Hallmark singscore extension complete.\n")

