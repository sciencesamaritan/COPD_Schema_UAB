#!/usr/bin/env Rscript

# Consolidated nomenclature audit for gene symbols displayed by, or used as
# gene-set keys in, the downstream main and supplementary figures.

suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]))
root <- dirname(dirname(script_file))
table_dir <- file.path(root, "tables")
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

de <- fread(de_file)
fg <- fread(fgsea_file)[
  Model == "Model_1" & Classification == "DisCOPDAI" &
    Contrast %chin% c("DisCOPDAI_1vs0", "DisCOPDAI_2vs0", "DisCOPDAI_3vs0")
]
representatives <- fread(file.path(table_dir, "pathway_integrated_evidence_model1.csv"))
membership <- fread(file.path(table_dir, "pathway_profile_cluster_membership_model1.csv"))
figure3 <- fread(file.path(table_dir, "enhanced_volcano_label_manifest.csv"))
figure5_drivers <- fread(file.path(table_dir, "pathway_profile_cluster_drivers_model1.csv"))
figure5_shared <- fread(file.path(table_dir, "pathway_profile_shared_significant_genes_model1.csv"))
figure5_pool <- fread(file.path(
  table_dir, "figure5_leading_edge_hgnc_validation_biobouncer.csv"
))

scopes <- list(
  list(
    figure = "Figures 2, 3, 4, S1, and S2",
    scope = "Archived Model 1 DE gene-label universe",
    genes = unique(de$gene_label)
  ),
  list(
    figure = "Figure 3",
    scope = "Printed volcano labels",
    genes = unique(figure3$gene_label)
  ),
  list(
    figure = "Figure 4",
    scope = "Archived representative-pathway leading-edge labels",
    genes = split_leading_edge(representatives$leadingEdge_archived)
  ),
  list(
    figure = "Figure 4",
    scope = "Normalized representative-pathway leading-edge genes",
    genes = split_leading_edge(representatives$leadingEdge)
  ),
  list(
    figure = "Figure 5",
    scope = "Archived candidate-pathway leading-edge labels",
    genes = figure5_pool$input
  ),
  list(
    figure = "Figure 5",
    scope = "Normalized candidate-pathway leading-edge genes used for driver selection",
    genes = figure5_pool$approved_symbol
  ),
  list(
    figure = "Figure 5",
    scope = "Printed cluster and shared-effect labels",
    genes = unique(c(figure5_drivers$gene, figure5_shared$gene))
  )
)

audits <- rbindlist(lapply(scopes, function(item) {
  audit_hgnc_symbols(item$genes, item$scope, item$figure)
}), use.names = TRUE, fill = TRUE)

fwrite(
  audits,
  file.path(table_dir, "downstream_figure_hgnc_validation_biobouncer.csv")
)
summary <- summarize_hgnc_audit(audits)
setorder(summary, figure, scope)
fwrite(
  summary,
  file.path(table_dir, "downstream_figure_nomenclature_summary.csv")
)
print(summary)

