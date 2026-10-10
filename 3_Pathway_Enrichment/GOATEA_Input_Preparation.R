#!/usr/bin/env Rscript

# Prepare the three primary Model 1 contrasts for local GOATEA review.
# This script changes file format and identifier namespace only; it does not
# refit differential-expression or pathway models.

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]))
root <- dirname(dirname(script_file))

source_file <- paste0(
  "path/to/",
  "Schema_omics_csv_bhatt/copdschema/extracted/",
  "gene_de_Model_1_DisCOPDAI_allresults_2026.04.25.csv"
)
out_dir <- file.path(root, "R", "goatea_inputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

contrast_names <- c(
  DisCOPDAI_1vs0 = "newly_diagnosed_vs_no_copd",
  DisCOPDAI_2vs0 = "potentially_excluded_vs_no_copd",
  DisCOPDAI_3vs0 = "established_copd_vs_no_copd"
)

de <- read.csv(source_file, check.names = FALSE)
required <- c("gene_label", "id", "logFC", "P.Value", "adj.P.Val", "Contrast")
stopifnot(all(required %in% names(de)))

ensembl <- sub("\\..*$", "", de$id)
entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(ensembl),
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)
de$gene <- unname(entrez[ensembl])

manifest <- list()
for (contrast in names(contrast_names)) {
  x <- de[de$Contrast == contrast & !is.na(de$gene), ]
  x <- x[order(x$adj.P.Val, -abs(x$logFC)), ]
  x <- x[!duplicated(x$gene), ]

  goatea_input <- data.frame(
    gene = as.integer(x$gene),
    symbol = x$gene_label,
    effectsize = x$logFC,
    pvalue = x$P.Value,
    signif = x$adj.P.Val < 0.05,
    fdr = x$adj.P.Val,
    moderated_t = x$t,
    stringsAsFactors = FALSE
  )

  output_file <- file.path(out_dir, paste0(contrast_names[[contrast]], ".tsv"))
  write.table(
    goatea_input,
    output_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )

  manifest[[length(manifest) + 1]] <- data.frame(
    contrast = contrast,
    label = contrast_names[[contrast]],
    source_rows = sum(de$Contrast == contrast),
    mapped_unique_entrez = nrow(goatea_input),
    fdr_significant = sum(goatea_input$signif),
    output = output_file,
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, manifest)
write.csv(manifest, file.path(out_dir, "goatea_input_manifest.csv"), row.names = FALSE)
print(manifest[, c("label", "source_rows", "mapped_unique_entrez", "fdr_significant")])


