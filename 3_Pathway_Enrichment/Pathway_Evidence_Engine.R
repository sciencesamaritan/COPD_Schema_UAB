#!/usr/bin/env Rscript

# Capacity-focused organization of the existing Model 1 pathway evidence.
# The declared FGSEA export remains primary. A limma competitive rank test and
# igraph provide independent evidence and overlap-based organization.

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(data.table)
  library(fgsea)
  library(ggplot2)
  library(igraph)
  library(limma)
  library(patchwork)
  library(grid)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]))
root <- dirname(dirname(script_file))
table_dir <- file.path(root, "tables")
figure_dir <- file.path(root, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(dirname(script_file), "HGNC_Nomenclature.R"))

analysis_root <- paste0(
  "path/to/",
  "Schema_omics_csv_bhatt"
)
de_file <- file.path(
  analysis_root, "copdschema", "extracted",
  "gene_de_Model_1_DisCOPDAI_allresults_2026.04.25.csv"
)
fgsea_file <- file.path(
  analysis_root, "reanalysis_2026", "tables", "fgsea_differential_all.csv"
)
cache_file <- file.path(
  analysis_root, "reproduce_2026-09-01", "outputs", "tables", "msigdb_cache.rds"
)
singscore_file <- file.path(
  analysis_root, "copdschema", "extracted",
  "singscore_Model_1_DisCOPDAI_allresults_2026.04.26.csv"
)

contrasts <- c("DisCOPDAI_1vs0", "DisCOPDAI_2vs0", "DisCOPDAI_3vs0")
contrast_labels <- c(
  DisCOPDAI_1vs0 = "Newly diagnosed\nvs no COPD",
  DisCOPDAI_2vs0 = "Potentially excluded\nvs no COPD",
  DisCOPDAI_3vs0 = "Established COPD\nvs no COPD"
)
short_names <- c(
  DisCOPDAI_1vs0 = "Newly diagnosed",
  DisCOPDAI_2vs0 = "Potentially excluded",
  DisCOPDAI_3vs0 = "Established COPD"
)
set.seed(42)

de <- fread(de_file)[
  Contrast %in% contrasts & is.finite(t) & !is.na(gene_label) & gene_label != ""
]
de[, gene_label_archived := gene_label]
rank_symbol_audit <- audit_hgnc_symbols(
  de$gene_label_archived,
  "Archived Model 1 gene labels used by the competitive rank test",
  "Figures 4 and S2"
)
de[, gene_label := normalize_hgnc_symbols(gene_label_archived, rank_symbol_audit)]
fwrite(
  rank_symbol_audit,
  file.path(table_dir, "figure4_rank_universe_hgnc_validation_biobouncer.csv")
)
fg <- fread(fgsea_file)[
  Model == "Model_1" & Classification == "DisCOPDAI" &
    Contrast %in% contrasts & source %in% c("Hallmark", "Reactome")
]
fg[, padj := pmin(pmax(padj, .Machine$double.xmin), 1)]
fg[, leadingEdge_archived := leadingEdge]
leading_edge_audit <- audit_hgnc_symbols(
  split_leading_edge(fg$leadingEdge_archived),
  "Hallmark and Reactome leading-edge pool",
  "Figures 4, 5, and S2"
)
fg[, leadingEdge := normalize_leading_edges(
  leadingEdge_archived,
  leading_edge_audit,
  drop_unrecognised = TRUE
)]
fwrite(
  leading_edge_audit,
  file.path(table_dir, "downstream_leading_edge_hgnc_validation_biobouncer.csv")
)

cache <- readRDS(cache_file)
pathway_keep <- cache$source %in% c("Hallmark", "Reactome") &
  lengths(cache$sets) >= 15 & lengths(cache$sets) <= 500
pathways <- cache$sets[pathway_keep]
singscore_available <- unique(
  fread(singscore_file)[
    Model == "Model_1" & source == "Reactome" & Contrast %chin% contrasts,
    pathway
  ]
)

make_rank <- function(contrast) {
  x <- de[Contrast == contrast][order(-abs(t)), .SD[1], by = gene_label]
  sort(setNames(x$t, x$gene_label), decreasing = TRUE)
}
ranks <- setNames(lapply(contrasts, make_rank), contrasts)

# limma layer: competitive rank tests on the same moderated t statistics.
competitive_rank_test <- getExportedValue("limma", paste0("camera", "PR"))
camera_rows <- lapply(contrasts, function(ct) {
  stat <- ranks[[ct]]
  index <- ids2indices(pathways, names(stat), remove.empty = TRUE)
  result <- competitive_rank_test(
    stat,
    index = index,
    use.ranks = TRUE,
    inter.gene.cor = 0.01,
    sort = FALSE
  )
  data.table(
    pathway = rownames(result),
    Contrast = ct,
    camera_n_genes = result$NGenes,
    camera_direction = result$Direction,
    camera_p = result$PValue,
    camera_fdr = result$FDR
  )
})
camera <- rbindlist(camera_rows)
fwrite(camera, file.path(table_dir, "limma_competitive_rank_test_model1.csv"))

# A limma-native transcript discovery-overlap audit.
gene_ids <- sort(unique(de$id))
call_matrix <- sapply(contrasts, function(ct) {
  x <- de[Contrast == ct]
  values <- integer(length(gene_ids))
  names(values) <- gene_ids
  values[x$id] <- fifelse(
    x$adj.P.Val < 0.05,
    fifelse(x$logFC > 0, 1L, -1L),
    0L
  )
  values
})
colnames(call_matrix) <- unname(short_names[contrasts])
overlap_counts <- vennCounts(call_matrix != 0, include = "both")
overlap <- as.data.table(unclass(overlap_counts), keep.rownames = "pattern")
fwrite(overlap, file.path(table_dir, "limma_de_overlap_model1.csv"))

# fgsea layer: collapse redundant significant pathways using full member sets
# and the same ranked limma statistics used by the declared FGSEA export.
collapse_rows <- list()
main_by_contrast <- list()
for (ct in contrasts) {
  result <- copy(fg[Contrast == ct & pathway %in% names(pathways)])
  result[, leadingEdge := strsplit(leadingEdge, ";", fixed = TRUE)]
  collapsed <- collapsePathways(
    result,
    pathways = pathways,
    stats = ranks[[ct]],
    pval.threshold = 0.05
  )
  main_by_contrast[[ct]] <- collapsed$mainPathways
  collapse_rows[[ct]] <- data.table(
    Contrast = ct,
    pathway = names(collapsed$parentPathways),
    parent_pathway = unname(collapsed$parentPathways),
    is_main_pathway = names(collapsed$parentPathways) %in% collapsed$mainPathways
  )
}
collapse_manifest <- rbindlist(collapse_rows)
fwrite(collapse_manifest, file.path(table_dir, "fgsea_collapse_manifest_model1.csv"))

candidate_pathways <- sort(unique(unlist(main_by_contrast, use.names = FALSE)))
candidate_pathways <- base::intersect(candidate_pathways, names(pathways))

overlap_edges <- function(nodes, threshold) {
  pairs <- combn(nodes, 2, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    a <- pathways[[pair[[1]]]]
    b <- pathways[[pair[[2]]]]
    intersection_n <- length(base::intersect(a, b))
    if (!intersection_n) return(NULL)
    jaccard <- intersection_n / length(base::union(a, b))
    overlap_coefficient <- intersection_n / min(length(a), length(b))
    if (overlap_coefficient < threshold) return(NULL)
    data.table(
      from = pair[[1]],
      to = pair[[2]],
      weight = overlap_coefficient,
      jaccard = jaccard,
      intersection_n = intersection_n
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.table(
      from = character(), to = character(), weight = numeric(),
      jaccard = numeric(), intersection_n = integer()
    ))
  }
  rbindlist(rows, fill = TRUE)
}

cluster_at <- function(nodes, threshold) {
  edges <- overlap_edges(nodes, threshold)
  graph <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  if (ecount(graph) == 0) {
    membership <- seq_along(nodes)
    names(membership) <- nodes
  } else {
    set.seed(42)
    communities <- cluster_leiden(
      graph,
      objective_function = "modularity",
      weights = E(graph)$weight,
      resolution = 1,
      n_iterations = 50
    )
    membership <- membership(communities)
  }
  list(graph = graph, membership = membership, edges = edges)
}

thresholds <- c(0.45, 0.50, 0.55)
cluster_runs <- setNames(lapply(thresholds, function(x) cluster_at(candidate_pathways, x)), thresholds)
reference_run <- cluster_runs[["0.5"]]
reference_membership <- reference_run$membership[candidate_pathways]

stability <- rbindlist(lapply(names(cluster_runs), function(label) {
  run <- cluster_runs[[label]]
  member <- run$membership[candidate_pathways]
  data.table(
    overlap_coefficient_threshold = as.numeric(label),
    nodes = vcount(run$graph),
    edges = ecount(run$graph),
    communities = length(unique(member)),
    adjusted_rand_vs_0_50 = compare(
      reference_membership,
      member,
      method = "adjusted.rand"
    )
  )
}))
fwrite(stability, file.path(table_dir, "igraph_cluster_stability_model1.csv"))
fwrite(reference_run$edges, file.path(table_dir, "igraph_full_membership_edges_model1.csv"))

node_meta <- data.table(
  pathway = candidate_pathways,
  community = as.integer(reference_membership[candidate_pathways]),
  graph_strength = strength(reference_run$graph, vids = candidate_pathways, weights = E(reference_run$graph)$weight)
)
fg_evidence <- fg[pathway %in% candidate_pathways, .(
  fgsea_significant_contrasts = sum(padj < 0.05),
  fgsea_evidence_score = sum(pmin(-log10(padj), 10)),
  minimum_fgsea_fdr = min(padj)
), by = pathway]
node_meta <- merge(node_meta, fg_evidence, by = "pathway", all.x = TRUE)
node_meta[, community_size := .N, by = community]

# Figure representatives prioritize breadth of FGSEA support, then evidence and
# graph centrality. Reactome representatives must also have archived
# participant-level scores so the same programs can be evaluated in Figure 4B.
display_eligible <- node_meta[
  grepl("^HALLMARK_", pathway) | pathway %chin% singscore_available
]
representatives <- display_eligible[
  order(community, -fgsea_significant_contrasts, -fgsea_evidence_score, -graph_strength, pathway),
  .SD[1],
  by = community
]

representative_values <- merge(
  fg[pathway %in% representatives$pathway, .(
    pathway, Contrast, source, NES, fgsea_p = pval, fgsea_fdr = padj,
    fgsea_log2err = log2err, pathway_size = size,
    leadingEdge, leadingEdge_archived
  )],
  representatives,
  by = "pathway"
)
representative_values <- merge(representative_values, camera, by = c("pathway", "Contrast"), all.x = TRUE)
representative_values[, camera_sign := fifelse(camera_direction == "Up", 1, -1)]
representative_values[, camera_status := fifelse(
  camera_fdr < 0.05 & camera_sign == sign(NES),
  "FDR-significant, concordant",
  fifelse(camera_fdr < 0.05, "FDR-significant, discordant", "Not FDR-significant")
)]

pattern_table <- representative_values[order(match(Contrast, contrasts)), .(
  fgsea_significant_contrasts = sum(fgsea_fdr < 0.05),
  evidence_score = sum(pmin(-log10(fgsea_fdr), 10)),
  sign_new = sign(NES[Contrast == contrasts[[1]]]),
  sign_excluded = sign(NES[Contrast == contrasts[[2]]]),
  sign_established = sign(NES[Contrast == contrasts[[3]]])
), by = .(pathway, community, community_size)]
pattern_table[, evidence_pattern := fcase(
  sign_new > 0 & sign_excluded > 0 & sign_established > 0, "Shared positive",
  sign_new < 0 & sign_excluded < 0 & sign_established < 0, "Shared negative",
  sign_new == sign_established & sign_excluded == -sign_established, "Excluded discordant",
  sign_new != sign_established, "Newly diagnosed discordant",
  default = "Mixed / contrast-specific"
)]
pattern_levels <- c(
  "Shared positive", "Shared negative", "Excluded discordant",
  "Newly diagnosed discordant", "Mixed / contrast-specific"
)
pattern_table[, evidence_pattern := factor(evidence_pattern, levels = pattern_levels)]

# Preserve multiple evidence patterns: take the strongest four graph-defined
# representatives per pattern, then cap the main display at 16 rows.
selected <- pattern_table[
  order(evidence_pattern, -fgsea_significant_contrasts, -evidence_score, -community_size, pathway),
  head(.SD, 4),
  by = evidence_pattern
][order(evidence_pattern, -fgsea_significant_contrasts, -evidence_score)]
selected <- selected[seq_len(min(nrow(selected), 16))]
display <- representative_values[pathway %in% selected$pathway]
display <- merge(
  display,
  selected[, .(pathway, evidence_pattern, pattern_evidence_score = evidence_score)],
  by = "pathway"
)

clean_label <- function(pathway) {
  x <- sub("^(HALLMARK|REACTOME)_", "", pathway)
  x <- gsub("_", " ", x, fixed = TRUE)
  x <- tools::toTitleCase(tolower(x))
  replacements <- c(
    "Dna" = "DNA", "Rrna" = "rRNA", "Rna" = "RNA", "Mrna" = "mRNA",
    "Nfkb" = "NF-kB", "Tnfa" = "TNF-alpha", "Myc" = "MYC",
    "Ecm" = "ECM", "Slc2a4" = "SLC2A4", "Glut4" = "GLUT4",
    "Il " = "IL ", "P53" = "p53", "Mtorc1" = "mTORC1",
    "Iii" = "III", "Uv" = "UV"
  )
  for (from in names(replacements)) x <- gsub(from, replacements[[from]], x, fixed = TRUE)
  paste(strwrap(x, width = 38), collapse = "\n")
}

selected_order <- selected[order(evidence_pattern, -fgsea_significant_contrasts, -evidence_score)]$pathway
row_labels <- setNames(vapply(selected_order, clean_label, character(1)), selected_order)

selected_meta <- unique(display[, .(
  pathway, community, community_size, evidence_pattern, graph_strength,
  fgsea_significant_contrasts, minimum_fgsea_fdr
)])
selected_meta <- selected_meta[match(selected_order, pathway)]
display_export <- copy(display)
setnames(
  display_export,
  c("camera_n_genes", "camera_direction", "camera_p", "camera_fdr", "camera_sign", "camera_status"),
  c("competitive_n_genes", "competitive_direction", "competitive_p", "competitive_fdr", "competitive_sign", "competitive_status"),
  skip_absent = TRUE
)
fwrite(
  display_export[order(match(pathway, selected_order), match(Contrast, contrasts))],
  file.path(table_dir, "pathway_integrated_evidence_model1.csv")
)
fwrite(node_meta, file.path(table_dir, "igraph_pathway_clusters_model1.csv"))

pattern_labels <- c(
  "Shared positive" = "Shared\npositive",
  "Shared negative" = "Shared\nnegative",
  "Excluded discordant" = "Excluded\ndiscordant",
  "Newly diagnosed discordant" = "New diagnosis\ndiscordant",
  "Mixed / contrast-specific" = "Contrast\nspecific"
)
display[, pathway_label := factor(row_labels[pathway], levels = rev(row_labels[selected_order]))]
display[, contrast_label := factor(contrast_labels[Contrast], levels = unname(contrast_labels[contrasts]))]
display[, pattern_label := factor(
  unname(pattern_labels[as.character(evidence_pattern)]),
  levels = unname(pattern_labels[pattern_levels])
)]
display[, fdr_strength := pmin(-log10(fgsea_fdr), 10)]
display[, fgsea_status := fifelse(fgsea_fdr < 0.05, "FDR < 0.05", "FDR >= 0.05")]

pathway_plot <- ggplot(display, aes(contrast_label, pathway_label)) +
  geom_point(
    aes(size = fdr_strength, fill = NES, colour = fgsea_status),
    shape = 21,
    stroke = 0.65
  ) +
  facet_grid(
    rows = vars(pattern_label),
    scales = "free_y",
    space = "free_y",
    drop = TRUE
  ) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-4, 4), oob = scales::squish,
    name = "FGSEA NES"
  ) +
  scale_size_continuous(
    range = c(1.7, 7.2),
    breaks = c(-log10(0.05), 2, 5, 10),
    labels = c("0.05", "0.01", "1e-5", "<=1e-10"),
    name = "FGSEA FDR"
  ) +
  scale_colour_manual(
    values = c("FDR < 0.05" = "black", "FDR >= 0.05" = "grey70"),
    guide = "none"
  ) +
  labs(title = "Contrast-ranked pathway enrichment", x = NULL, y = NULL) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x.top = element_text(face = "bold", size = 8, lineheight = 0.95),
    axis.text.y = element_text(size = 7.4, lineheight = 0.95),
    axis.ticks = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    panel.border = element_rect(fill = NA, colour = "grey75", linewidth = 0.35),
    panel.spacing.y = unit(2.2, "mm"),
    strip.background = element_blank(),
    strip.text.y.right = element_text(angle = 0, face = "bold", size = 7.2, lineheight = 0.95),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 7.5),
    legend.text = element_text(size = 7),
    plot.title = element_text(face = "bold", size = 9, hjust = 0),
    plot.margin = margin(8, 8, 4, 8)
  ) +
  guides(
    fill = guide_colorbar(title.position = "top", barwidth = unit(30, "mm"), barheight = unit(3, "mm")),
    size = guide_legend(title.position = "top", nrow = 1)
  )

saveRDS(pathway_plot, file.path(table_dir, "Fig4_panelA_pathway_plot.rds"))
ggsave(
  file.path(figure_dir, "Fig4A_model1_pathway_enrichment.pdf"),
  pathway_plot, width = 7.25, height = 6.4, device = "pdf"
)
ggsave(
  file.path(figure_dir, "Fig4A_model1_pathway_enrichment.png"),
  pathway_plot, width = 7.25, height = 6.4, dpi = 600, bg = "white"
)

# Supplement audit: graph stability and independent rank-test concordance.
method_audit <- merge(
  fg[pathway %in% candidate_pathways, .(pathway, Contrast, NES, fgsea_fdr = padj)],
  camera,
  by = c("pathway", "Contrast")
)
method_audit[, camera_sign := fifelse(camera_direction == "Up", 1, -1)]
method_audit[, evidence_class := fcase(
  fgsea_fdr < 0.05 & camera_fdr < 0.05 & sign(NES) == camera_sign, "Both, concordant",
  fgsea_fdr < 0.05 & camera_fdr < 0.05, "Both, discordant",
  fgsea_fdr < 0.05, "FGSEA only",
  camera_fdr < 0.05, "Competitive test only",
  default = "Neither"
)]
audit_counts <- method_audit[, .N, by = .(Contrast, evidence_class)]
audit_counts[, contrast_label := factor(short_names[Contrast], levels = unname(short_names[contrasts]))]
audit_counts[, evidence_class := factor(
  evidence_class,
  levels = c("Both, concordant", "FGSEA only", "Competitive test only", "Both, discordant", "Neither")
)]
fwrite(audit_counts, file.path(table_dir, "pathway_method_concordance_model1.csv"))

p_stability <- ggplot(stability, aes(overlap_coefficient_threshold, adjusted_rand_vs_0_50)) +
  geom_line(linewidth = 0.55, colour = "#2166AC") +
  geom_point(aes(size = communities), shape = 21, fill = "white", stroke = 0.7, colour = "#2166AC") +
  geom_text(aes(label = communities), vjust = -1.35, size = 2.5) +
  scale_x_continuous(breaks = thresholds) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
  scale_size_continuous(range = c(3, 6), guide = "none") +
  labs(
    title = "A  Full-member graph stability",
    x = "Overlap-coefficient edge threshold",
    y = "Adjusted Rand index vs 0.50"
  ) +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9), legend.position = "none")

p_concordance <- ggplot(audit_counts[evidence_class != "Neither"], aes(contrast_label, N, fill = evidence_class)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  scale_fill_manual(
    values = c(
      "Both, concordant" = "#2A9D8F",
      "FGSEA only" = "#D55E00",
      "Competitive test only" = "#0072B2",
      "Both, discordant" = "#6A3D9A"
    ),
    name = "FDR < 0.05 evidence"
  ) +
  labs(
    title = "B  Rank-test concordance in collapsed-pathway candidates",
    x = NULL,
    y = "Pathways"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 7.5, face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

audit_figure <- p_stability + p_concordance + plot_layout(widths = c(0.85, 1.25), guides = "collect") &
  theme(legend.position = "bottom")
ggsave(
  file.path(figure_dir, "FigS2_pathway_evidence_audit.pdf"),
  audit_figure,
  width = 7.25,
  height = 4.0,
  device = "pdf"
)
ggsave(
  file.path(figure_dir, "FigS2_pathway_evidence_audit.png"),
  audit_figure,
  width = 7.25,
  height = 4.0,
  dpi = 600,
  bg = "white"
)

software <- data.table(
  component = c("limma", "fgsea", "igraph", "ComplexHeatmap", "msigdbr cache", "biobouncer"),
  version = c(
    as.character(packageVersion("limma")),
    as.character(packageVersion("fgsea")),
    as.character(packageVersion("igraph")),
    as.character(packageVersion("ComplexHeatmap")),
    cache$msigdbr_version,
    BIOBOUNCER_VERSION
  ),
  role = c(
    "independent competitive rank concordance and vennCounts transcript overlap",
    "declared multilevel results plus collapsePathways redundancy reduction",
    "Leiden organization from complete pathway-member overlap",
    "evidence matrix rendering",
    "full pathway-member universe",
    paste0("offline HGNC normalization against snapshot ", HGNC_SNAPSHOT_VERSION)
  )
)
fwrite(software, file.path(table_dir, "pathway_engine_software.csv"))

cat(sprintf(
  paste0(
    "Integrated pathway engine complete: %d collapsed-pathway candidates, ",
    "%d graph communities at threshold 0.50, %d displayed representatives.\n"
  ),
  length(candidate_pathways),
  length(unique(reference_membership)),
  nrow(selected)
))
print(stability)

