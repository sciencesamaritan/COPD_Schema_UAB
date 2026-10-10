# =============================================================================
# 01_load.R -- differential expression, with standard errors recovered from CI
# =============================================================================

#' Load the gene-level DE table for the configured model/case definition.
#' Recovers per-gene SE from the limma 95% CI half-width and *verifies* the
#' recovery by reproducing the reported moderated t. The check is relative,
#' not absolute: the source rounds CI bounds to 4 decimals, which floors the
#' achievable precision.
load_de <- function(verbose = TRUE) {
  f <- cfg_path("gene_de")
  if (verbose) log_msg("reading ", basename(f))
  de <- fread(f)
  de[, short := sub(".*_(\\d+vs\\d+)$", "\\1", Contrast)]
  de[, se := (CI.R - CI.L) / (2 * CFG$params$z975)]

  chk <- de[is.finite(t) & is.finite(se) & se > 0]
  rel <- max(abs(chk$logFC / chk$se - chk$t) / abs(chk$t), na.rm = TRUE)
  abs_worst <- max(abs(chk$logFC / chk$se - chk$t), na.rm = TRUE)
  if (verbose) log_msg(sprintf(
    "SE recovery: max relative error %.2e (worst absolute t discrepancy %.4f)",
    rel, abs_worst))
  if (rel >= 5e-3) stop(sprintf("SE recovery failed: max relative error %.3e", rel))

  wide <- function(v) {
    m <- dcast(de, id ~ short, value.var = v)
    ids <- m$id; m[, id := NULL]
    as.matrix(m) -> mm; rownames(mm) <- ids; mm
  }
  sym <- de[!duplicated(id), .(id, gene_label)]

  list(de = de, lfc = wide("logFC"), padj = wide("adj.P.Val"),
       se = wide("se"), t = wide("t"),
       sym = setNames(sym$gene_label, sym$id),
       se_rel_error = rel, se_abs_worst = abs_worst,
       contrasts = sort(unique(de$short)), n_genes = uniqueN(de$id))
}

#' Per-contrast significant-gene counts.
de_counts <- function(D, fdr = CFG$params$fdr) {
  ct <- names(CFG$params$contrasts)
  rbindlist(lapply(ct, function(g) {
    p <- D$padj[, g]
    data.table(contrast = g, group = CFG$params$contrasts[[g]],
               transcripts_FDR05 = sum(p < fdr, na.rm = TRUE),
               unique_symbols = uniqueN(D$sym[names(p)[which(p < fdr)]]),
               min_FDR = round(min(p, na.rm = TRUE), 4))
  }))
}
