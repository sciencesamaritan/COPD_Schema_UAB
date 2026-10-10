# =============================================================================
# 03_enrichment.R -- FGSEA recompute, fidelity against the distributed tables,
# and the symbol-keying diagnostic that separates the two candidate causes of
# the FDR-call disagreement.
# =============================================================================
suppressPackageStartupMessages({ library(fgsea); library(msigdbr) })

MSIG_SPEC <- list(
  Hallmark     = list("H",  NULL),
  Reactome     = list("C2", "CP:REACTOME"),
  BioCarta     = list("C2", "CP:BIOCARTA"),
  KEGG_MEDICUS = list("C2", "CP:KEGG_MEDICUS"),
  ImmuneSigDB  = list("C7", "IMMUNESIGDB")
)

#' Build the five MSigDB collections keyed by HGNC symbol. Cached to disk:
#' msigdbr is slow and the content is version-pinned by the installed package.
build_pathways <- function(cache = file.path(cfg_out("tables"), "msigdb_cache.rds")) {
  ver <- as.character(packageVersion("msigdbr"))
  if (file.exists(cache)) {
    cached <- readRDS(cache)
    # A cache built under a different msigdbr is a different gene-set universe.
    if (identical(cached$msigdbr_version, ver)) {
      log_msg("pathway cache hit: ", basename(cache), " (msigdbr ", ver, ")")
      return(cached)
    }
    log_msg("pathway cache was built under msigdbr ",
            cached$msigdbr_version %||% "unrecorded", " but ", ver,
            " is installed; rebuilding", level = "WARN")
  }
  log_msg("building MSigDB collections (msigdbr ", ver, ")")
  sets <- list(); src <- character(0)
  for (nm in names(MSIG_SPEC)) {
    s <- MSIG_SPEC[[nm]]
    df <- as.data.table(msigdbr(species = "Homo sapiens", collection = s[[1]],
                                subcollection = s[[2]]))
    l <- lapply(split(as.character(df$gene_symbol), df$gs_name), unique)
    log_msg(sprintf("  %-13s %5d sets", nm, length(l)))
    sets <- c(sets, l); src <- c(src, setNames(rep(nm, length(l)), names(l)))
  }
  out <- list(sets = sets, source = src, msigdbr_version = ver,
              built = as.character(Sys.Date()))
  saveRDS(out, cache)
  out
}

#' Pre-ranked FGSEA for one contrast from the provided moderated t-statistics.
#' `key` selects the gene identifier: "symbol" reproduces the 2026 recompute,
#' "symbol_validated" drops labels that are not current HGNC symbols.
run_fgsea_contrast <- function(D, contrast, P, key = "symbol", drop_ids = NULL,
                               minSize = 15, maxSize = 500, seed = CFG$params$seed) {
  set.seed(seed)
  sub <- D$de[short == contrast & is.finite(t) & !is.na(gene_label) & gene_label != ""]
  if (!is.null(drop_ids)) sub <- sub[!gene_label %in% drop_ids]
  sub <- sub[order(-abs(t))][, .SD[1], by = gene_label]     # collapse dups, keep max|t|
  ranks <- setNames(sub$t, sub$gene_label)
  ranks <- sort(ranks, decreasing = TRUE)
  res <- as.data.table(fgseaMultilevel(pathways = P$sets, stats = ranks,
                                       minSize = minSize, maxSize = maxSize, eps = 0))
  if (nrow(res)) {
    res[, leadingEdge := vapply(leadingEdge, paste, "", collapse = ";")]
    res[, source := P$source[pathway]]
  }
  # set unconditionally so an empty result still carries its identity columns
  res[, `:=`(Model = CFG$meta$model, Contrast = contrast, n_ranked = length(ranks))]
  res[]
}

#' Concordance of a recompute against the distributed GSEA table for one contrast.
fidelity <- function(recomputed, shipped, contrast, fdr = CFG$params$fdr) {
  a <- recomputed[Contrast == contrast, .(pathway, NES, padj)]
  b <- shipped[short == contrast, .(pathway, NES_ship = NES, padj_ship = padj)]
  m <- merge(a, b, by = "pathway")
  if (!nrow(m)) return(data.table(contrast = contrast, shared = 0L))
  sa <- m$padj < fdr; sb <- m$padj_ship < fdr
  data.table(contrast = contrast, shared = nrow(m),
             NES_rho = round(cor(m$NES, m$NES_ship, method = "spearman"), 4),
             fdr_agree = round(mean(sa == sb), 4),
             fdr_jaccard = round(sum(sa & sb) / max(sum(sa | sb), 1), 4),
             n_sig_recomputed = sum(sa), n_sig_shipped = sum(sb))
}

#' Load the distributed GSEA table, normalising the contrast label.
load_shipped_gsea <- function() {
  f <- cfg_path("gsea")
  log_msg("reading ", basename(f))
  g <- fread(f)
  g[, short := sub(".*_(\\d+vs\\d+)$", "\\1", Contrast)]
  g[]
}

#' Pathway-level comparison between groups: which of the reference signature's
#' pathways are shared, absent, or reversed in each reclassified group.
pathway_partition <- function(gs, ref = CFG$params$reference_contrast,
                              groups = c("1vs0", "2vs0"), fdr = CFG$params$fdr) {
  R <- gs[short == ref & padj < fdr, .(pathway, NES_ref = NES)]
  rbindlist(lapply(groups, function(g) {
    G <- gs[short == g, .(pathway, NES_g = NES, padj_g = padj)]
    m <- merge(R, G, by = "pathway", all.x = TRUE)
    m[, status := fifelse(is.na(padj_g) | padj_g >= fdr, "not significant",
                  fifelse(sign(NES_g) == sign(NES_ref), "shared, same direction",
                          "significant, reversed"))]
    m[, .(contrast = g, group = CFG$params$contrasts[[g]], n_reference = nrow(R),
          n = .N, pct = round(100 * .N / nrow(R), 1)), by = status]
  }))
}
