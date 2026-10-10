#!/usr/bin/env Rscript

# Cluster collapsed Model 1 pathways by their three no-COPD contrast profiles.
# The statistical inputs are the existing FGSEA and limma exports; no model is
# refit here. Cluster number is chosen by average silhouette width for k=2..8.

suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]))
root <- dirname(dirname(script_file))
table_dir <- file.path(root, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(dirname(script_file), "HGNC_Nomenclature.R"))

analysis_root <- paste0(
  "path/to/",
  "Schema_omics_csv_bhatt"
)
fgsea_file <- file.path(
  analysis_root, "reanalysis_2026", "tables", "fgsea_differential_all.csv"
)
de_file <- file.path(
  analysis_root, "copdschema", "extracted",
  "gene_de_Model_1_DisCOPDAI_allresults_2026.04.25.csv"
)
candidate_file <- file.path(table_dir, "igraph_pathway_clusters_model1.csv")

contrasts <- c("DisCOPDAI_1vs0", "DisCOPDAI_2vs0", "DisCOPDAI_3vs0")
contrast_labels <- c(
  DisCOPDAI_1vs0 = "Newly diagnosed",
  DisCOPDAI_2vs0 = "Potentially excluded",
  DisCOPDAI_3vs0 = "Established COPD"
)

stopifnot(file.exists(fgsea_file), file.exists(de_file), file.exists(candidate_file))

fg <- fread(fgsea_file)[
  Model == "Model_1" & Classification == "DisCOPDAI" &
    Contrast %chin% contrasts & source %chin% c("Hallmark", "Reactome")
]
de <- fread(de_file)[Contrast %chin% contrasts & is.finite(t)]
candidates <- fread(candidate_file)$pathway
fg <- fg[pathway %chin% candidates]

wide <- dcast(
  fg[, .(pathway, Contrast, NES)],
  pathway ~ Contrast,
  value.var = "NES"
)
wide <- wide[complete.cases(wide[, ..contrasts])]
stopifnot(nrow(wide) == length(candidates))

# Standardize each contrast, then normalize each pathway vector. This clusters
# relative direction patterns without allowing the highest-variance contrast
# or the largest-magnitude pathway to determine the partition by itself.
profile_matrix <- scale(as.matrix(wide[, ..contrasts]))
profile_norm <- sqrt(rowSums(profile_matrix^2))
profile_matrix <- profile_matrix / profile_norm
rownames(profile_matrix) <- wide$pathway

distance_matrix <- as.matrix(dist(profile_matrix, method = "euclidean"))
tree <- hclust(as.dist(distance_matrix), method = "ward.D2")

mean_silhouette <- function(cluster) {
  widths <- vapply(seq_along(cluster), function(i) {
    same <- which(cluster == cluster[[i]] & seq_along(cluster) != i)
    if (!length(same)) return(0)
    within_distance <- mean(distance_matrix[i, same])
    other_clusters <- setdiff(unique(cluster), cluster[[i]])
    nearest_distance <- min(vapply(other_clusters, function(group) {
      mean(distance_matrix[i, cluster == group])
    }, numeric(1)))
    (nearest_distance - within_distance) / max(nearest_distance, within_distance)
  }, numeric(1))
  mean(widths)
}

k_selection <- rbindlist(lapply(2:8, function(k) {
  groups <- cutree(tree, k = k)
  sizes <- sort(table(groups), decreasing = TRUE)
  data.table(
    k = k,
    average_silhouette = mean_silhouette(groups),
    smallest_cluster = min(sizes),
    largest_cluster = max(sizes),
    cluster_sizes = paste(as.integer(sizes), collapse = ";")
  )
}))
chosen_k <- k_selection[which.max(average_silhouette), k]
fwrite(k_selection, file.path(table_dir, "pathway_profile_cluster_k_selection_model1.csv"))

raw_cluster <- cutree(tree, k = chosen_k)
dendrogram_order <- unique(raw_cluster[tree$order])
cluster_map <- setNames(seq_along(dendrogram_order), dendrogram_order)
cluster_number <- unname(cluster_map[as.character(raw_cluster)])
cluster_id <- paste0("C", cluster_number)
names(cluster_id) <- names(raw_cluster)

membership <- copy(wide)
membership[, `:=`(
  cluster_number = cluster_number,
  cluster_id = cluster_id,
  dendrogram_position = match(seq_len(.N), tree$order)
)]

evidence <- fg[, .(
  source = first(source),
  significant_contrasts = sum(padj < 0.05),
  minimum_fdr = min(padj),
  aggregate_evidence = sum(pmin(-log10(pmax(padj, .Machine$double.xmin)), 10))
), by = pathway]
membership <- merge(membership, evidence, by = "pathway", all.x = TRUE)

centroid_distance <- numeric(nrow(membership))
names(centroid_distance) <- membership$pathway
for (cid in unique(cluster_id)) {
  members <- names(cluster_id)[cluster_id == cid]
  centroid <- colMeans(profile_matrix[members, , drop = FALSE])
  centroid_distance[members] <- sqrt(rowSums(
    (profile_matrix[members, , drop = FALSE] -
      matrix(centroid, nrow = length(members), ncol = length(centroid), byrow = TRUE))^2
  ))
}
membership[, profile_distance_to_centroid := centroid_distance[pathway]]
setorder(membership, cluster_number, profile_distance_to_centroid, -aggregate_evidence)
fwrite(membership, file.path(table_dir, "pathway_profile_cluster_membership_model1.csv"))

# Give the two displayed pathways distinct roles. The profile exemplar is the
# centroid-nearest member. The central evidence member has the strongest
# aggregate FGSEA evidence among the closest quartile, excluding the exemplar.
profile_representatives <- membership[
  order(cluster_number, profile_distance_to_centroid, -aggregate_evidence),
  .SD[1],
  by = cluster_number
]
profile_representatives[, `:=`(
  representative_rank = 1L,
  representative_role = "Profile exemplar",
  selection_basis = paste0(
    "Minimum Euclidean distance to centroid in standardized, ",
    "unit-length 3-contrast space"
  )
)]

centrality_cutoffs <- membership[, .(
  centrality_cutoff = quantile(profile_distance_to_centroid, 0.25)
), by = cluster_number]
central_candidates <- merge(membership, centrality_cutoffs, by = "cluster_number")
central_candidates <- merge(
  central_candidates,
  profile_representatives[, .(cluster_number, profile_pathway = pathway)],
  by = "cluster_number"
)[
  profile_distance_to_centroid <= centrality_cutoff & pathway != profile_pathway
]
evidence_representatives <- central_candidates[
  order(cluster_number, -aggregate_evidence, profile_distance_to_centroid),
  .SD[1],
  by = cluster_number
]
evidence_representatives[, `:=`(
  representative_rank = 2L,
  representative_role = "Central evidence",
  selection_basis = paste0(
    "Maximum sum(min(-log10(FDR),10)) across contrasts among centroid-nearest ",
    "quartile, excluding PROFILE"
  )
)]

representatives <- rbindlist(
  list(profile_representatives, evidence_representatives),
  use.names = TRUE,
  fill = TRUE
)
setorder(representatives, cluster_number, representative_rank)
fwrite(representatives, file.path(table_dir, "pathway_profile_cluster_representatives_model1.csv"))

cluster_long <- melt(
  membership,
  id.vars = c("pathway", "cluster_number", "cluster_id"),
  measure.vars = contrasts,
  variable.name = "Contrast",
  value.name = "NES"
)
cluster_long <- merge(
  cluster_long,
  fg[, .(pathway, Contrast, padj)],
  by = c("pathway", "Contrast"),
  all.x = TRUE
)
cluster_summary <- cluster_long[, .(
  mean_NES = mean(NES),
  sd_NES = sd(NES),
  median_NES = median(NES),
  q1_NES = quantile(NES, 0.25),
  q3_NES = quantile(NES, 0.75),
  fdr_significant_pathways = sum(padj < 0.05),
  pathway_count = .N
), by = .(cluster_number, cluster_id, Contrast)]
cluster_summary[, contrast_label := unname(contrast_labels[Contrast])]

clean_pathway <- function(x) {
  x <- sub("^(HALLMARK|REACTOME)_", "", x)
  x <- tools::toTitleCase(tolower(gsub("_", " ", x, fixed = TRUE)))
  replacements <- c(
    "Dna" = "DNA", "Rna" = "RNA", "Mrna" = "mRNA", "Rrna" = "rRNA",
    "Nfkb" = "NF-kB", "Tnfa" = "TNF-alpha", "Mtorc1" = "mTORC1",
    "Ecm" = "ECM", "Il " = "IL ", "P53" = "p53", "Myc" = "MYC",
    "Iii" = "III", "Ca2 Pathway" = "Ca2+ Pathway",
    "Asparagine N Linked" = "Asparagine N-Linked",
    "Mapk1 Mapk3" = "MAPK1/MAPK3", "Snrnp" = "snRNP",
    "Trna" = "tRNA", "Hox" = "HOX", "Sirt1" = "SIRT1",
    "Interferon Alpha Beta" = "Interferon-alpha/beta"
  )
  for (from in names(replacements)) x <- gsub(from, replacements[[from]], x, fixed = TRUE)
  x
}

program_labels <- representatives[order(cluster_number, representative_rank), .(
  program_label = paste(
    paste0(representative_role, ": ", vapply(pathway, clean_pathway, character(1))),
    collapse = " | "
  )
), by = cluster_number]
cluster_summary <- merge(cluster_summary, program_labels, by = "cluster_number", all.x = TRUE)
fwrite(cluster_summary, file.path(table_dir, "pathway_profile_cluster_summary_model1.csv"))

# Leading-edge genes are counted only in pathways reaching FGSEA FDR < .05.
leading <- copy(fg[padj < 0.05 & nzchar(leadingEdge)])
leading <- merge(
  leading,
  membership[, .(pathway, cluster_number, cluster_id)],
  by = "pathway"
)
leading[, gene := strsplit(leadingEdge, ";", fixed = TRUE)]
leading <- leading[, .(gene = unique(trimws(unlist(gene)))), by = .(
  pathway, Contrast, cluster_number, cluster_id
)]
leading <- leading[nzchar(gene)]
leading[, gene_archived := gene]
leading_edge_audit <- audit_hgnc_symbols(
  leading$gene_archived,
  "FDR-significant candidate-pathway leading-edge pool",
  "Figure 5"
)
if (any(leading_edge_audit$status == "unrecognised")) {
  stop("Figure 5 leading-edge pool contains an unrecognised HGNC label")
}
leading[, gene := normalize_hgnc_symbols(gene_archived, leading_edge_audit)]
leading <- unique(leading[, .(pathway, Contrast, cluster_number, cluster_id, gene)])
fwrite(
  leading_edge_audit,
  file.path(table_dir, "figure5_leading_edge_hgnc_validation_biobouncer.csv")
)

cluster_hits <- unique(leading[, .(pathway, Contrast, cluster_number, cluster_id, gene)])[
  , .(cluster_hits = .N), by = .(cluster_number, cluster_id, gene)
]
all_hits <- unique(leading[, .(pathway, Contrast, cluster_number, gene)])[
  , .(all_hits = .N), by = gene
]
gene_strength <- de[, .(maximum_absolute_t = max(abs(t), na.rm = TRUE)), by = .(gene = gene_label)]
cluster_drivers <- merge(cluster_hits, all_hits, by = "gene")
cluster_drivers <- merge(cluster_drivers, gene_strength, by = "gene", all.x = TRUE)
cluster_drivers[, specificity := cluster_hits / all_hits]
cluster_drivers[, driver_score := cluster_hits * specificity * maximum_absolute_t]
cluster_drivers <- cluster_drivers[
  order(cluster_number, -driver_score, -cluster_hits, gene),
  head(.SD, 6),
  by = cluster_number
]
cluster_drivers[, driver_rank := seq_len(.N), by = cluster_number]
fwrite(cluster_drivers, file.path(table_dir, "pathway_profile_cluster_drivers_model1.csv"))

contrast_hits <- unique(leading[, .(pathway, Contrast, gene)])[
  , .(pathway_hits = .N), by = .(Contrast, gene)
]
contrast_de <- de[, .(
  Contrast, gene = gene_label, moderated_t = t, log2_fold_change = logFC,
  gene_fdr = adj.P.Val
)]
contrast_drivers <- merge(contrast_hits, contrast_de, by = c("Contrast", "gene"))
contrast_drivers[, `:=`(
  direction = fifelse(moderated_t >= 0, "Higher", "Lower"),
  gene_fdr_status = fifelse(gene_fdr < 0.05, "FDR < 0.05", "FDR >= 0.05"),
  driver_score = pathway_hits * abs(moderated_t)
)]
contrast_drivers <- contrast_drivers[
  order(Contrast, direction, -driver_score, gene),
  head(.SD, 6),
  by = .(Contrast, direction)
]
contrast_drivers[, contrast_label := unname(contrast_labels[Contrast])]
contrast_drivers[, display_rank := frank(-driver_score, ties.method = "first"), by = Contrast]
fwrite(contrast_drivers, file.path(table_dir, "pathway_profile_contrast_drivers_model1.csv"))

# Use one shared gene set for direct visual comparison across all 3 contrasts.
# Eligibility requires gene-level FDR < .05 and membership in an
# FDR-significant pathway leading edge in both newly diagnosed and established
# COPD. Rank by the weaker absolute t statistic across those 2 contrasts, with
# 6 concordantly higher and 6 concordantly lower genes retained.
shared_reference_contrasts <- c("DisCOPDAI_1vs0", "DisCOPDAI_3vs0")
de_unique <- de[
  order(Contrast, gene_label, adj.P.Val, -abs(t)),
  .SD[1],
  by = .(Contrast, gene_label)
]
shared_gene_fdr <- Reduce(intersect, lapply(shared_reference_contrasts, function(contrast) {
  de_unique[Contrast == contrast & adj.P.Val < 0.05, gene_label]
}))
shared_leading_edge <- Reduce(intersect, lapply(shared_reference_contrasts, function(contrast) {
  unique(leading[Contrast == contrast, gene])
}))
shared_candidates <- intersect(shared_gene_fdr, shared_leading_edge)

shared_wide <- dcast(
  de_unique[
    gene_label %chin% shared_candidates & Contrast %chin% shared_reference_contrasts,
    .(gene_label, Contrast, t)
  ],
  gene_label ~ Contrast,
  value.var = "t"
)
shared_wide <- shared_wide[
  DisCOPDAI_1vs0 * DisCOPDAI_3vs0 > 0
]
shared_wide[, `:=`(
  shared_direction = fifelse(DisCOPDAI_1vs0 > 0, "Higher", "Lower"),
  joint_min_abs_t = pmin(abs(DisCOPDAI_1vs0), abs(DisCOPDAI_3vs0))
)]
selected_shared <- shared_wide[
  order(-joint_min_abs_t, gene_label),
  head(.SD, 6),
  by = shared_direction
]
selected_shared[, direction_order := fifelse(shared_direction == "Higher", 1L, 2L)]
setorder(selected_shared, direction_order, -joint_min_abs_t, gene_label)
selected_shared[, display_order := seq_len(.N)]

shared_plot <- merge(
  de_unique[
    gene_label %chin% selected_shared$gene_label & Contrast %chin% contrasts,
    .(
      gene = gene_label, Contrast, moderated_t = t,
      log2_fold_change = logFC, gene_fdr = adj.P.Val
    )
  ],
  selected_shared[, .(
    gene = gene_label, shared_direction, joint_min_abs_t, display_order
  )],
  by = "gene"
)
shared_plot[, `:=`(
  contrast_label = unname(contrast_labels[Contrast]),
  gene_fdr_status = fifelse(gene_fdr < 0.05, "FDR < 0.05", "FDR >= 0.05")
)]
setorder(shared_plot, Contrast, display_order)
stopifnot(
  nrow(shared_plot) == 36L,
  all(shared_plot[, .N, by = Contrast]$N == 12L),
  shared_plot[Contrast == "DisCOPDAI_2vs0", sum(gene_fdr < 0.05)] == 0L
)

# Validate every symbol printed in Figure 5 against the same pinned HGNC
# snapshot used by the source repository. Fail rather than silently relabeling
# a figure if a future driver selection introduces an obsolete or unknown name.
displayed_genes <- sort(unique(c(cluster_drivers$gene, shared_plot$gene)))
symbol_validation <- audit_hgnc_symbols(
  displayed_genes,
  "Printed cluster and shared-effect labels",
  "Figure 5"
)
fwrite(
  symbol_validation,
  file.path(table_dir, "figure5_hgnc_symbol_validation_biobouncer.csv")
)
if (!all(symbol_validation$valid %in% TRUE)) {
  failed <- symbol_validation[!(valid %in% TRUE), paste(input, status, sep = ": ")]
  stop("Figure 5 contains non-current HGNC symbols: ", paste(failed, collapse = "; "))
}

fwrite(
  shared_plot,
  file.path(table_dir, "pathway_profile_shared_significant_genes_model1.csv")
)

cat(sprintf(
  paste0(
    "Pathway profile clustering complete: %d pathways, k=%d, ",
    paste0(
      "mean silhouette=%.3f; %d shared significant leading-edge genes displayed; ",
      "%d/%d displayed symbols current by biobouncer %s (HGNC %s).\n"
    )
  ),
  nrow(membership), chosen_k,
  k_selection[k == chosen_k, average_silhouette],
  nrow(selected_shared),
  sum(symbol_validation$valid %in% TRUE), nrow(symbol_validation),
  unique(symbol_validation$biobouncer_version), unique(symbol_validation$version)
))

