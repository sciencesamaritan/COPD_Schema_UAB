# =============================================================================
# 05_reproduce.R -- disease-axis scaling (lambda), bootstrap, rho sensitivity, per-gene
# departures and FGSEA recompute. Writes outputs/tables/*.csv and a check table.
# =============================================================================

#' Disease-axis scaling from gene-level summary statistics. Records a pass/fail check
#' per reported number when reference tables are supplied.
reproduce_all <- function(with_fgsea = TRUE, fgsea_contrasts = c("1vs0", "2vs0", "3vs0")) {
  check_reset()
  T <- tryCatch(cfg_path("shipped_tables"), error = function(e) NULL)
  OUT <- cfg_out("tables")

  log_step("1. Load differential expression and recover standard errors")
  D <- load_de()
  cnt <- de_counts(D)
  print(cnt); fwrite(cnt, file.path(OUT, "de_counts.csv"))
  record_check("SE recovery: max relative error", D$se_rel_error, 0, 5e-3)

  log_step("2. Scaling along the established-disease axis")
  n_axis <- length(axis_mask(D))
  log_msg("reference signature (established COPD, FDR<0.05): ", n_axis, " transcripts")
  dec <- rbindlist(lapply(c("1vs0", "2vs0"), function(g) decompose(D, g)))
  print(dec); fwrite(dec, file.path(OUT, "scaling_coefficients.csv"))

  # Optional: check against previously reported values when the reference tables are supplied
  has_ref <- !is.null(T) && file.exists(file.path(T, "Table_E1_scaling_coefficients.csv"))
  if (has_ref) {
  E1 <- fread(file.path(T, "Table_E1_scaling_coefficients.csv"))
  E1p <- E1[Schema == "Qualitative AI"]
  record_check("n_transcripts on axis", n_axis, E1p$n_transcripts[1], 0)
  }
  if (has_ref) for (i in seq_len(nrow(dec))) {
    g <- dec$group[i]; s <- E1p[group == g]
    record_check(sprintf("lambda [%s]", g), dec$lambda[i], s$lambda, 5e-3)
    record_check(sprintf("rho [%s]", g), dec$rho[i], s$rho, 1e-3)
    record_check(sprintf("on-axis SD [%s]", g), dec$on_axis_sd[i], s$on_axis_sd, 1e-3)
    record_check(sprintf("%% variance on axis [%s]", g),
                 dec$pct_variance_on_axis[i], s$pct_variance_on_axis, 2.0, "pct pts")
  }

  log_step("3. Paired bootstrap and rho sensitivity")
  bt <- bootstrap_lambda(D)
  dd <- attr(bt, "diff")
  print(bt)
  log_msg(sprintf("observed difference %.3f | bootstrap SE %.3f | 95%% CI %.3f to %.3f",
                  dd$observed, dd$se, dd$ci[1], dd$ci[2]))
  log_msg(sprintf("replicates with difference <= 0: %d of %d | bootstrap r between coefficients %+.3f",
                  dd$n_le_zero, dd$B, dd$boot_cor))
  log_msg(sprintf("separation: %.1f-fold", dd$fold))
  fwrite(bt, file.path(OUT, "bootstrap_lambda.csv"))
  if (has_ref) for (i in seq_len(nrow(bt))) {
    s <- E1p[group == bt$group[i]]
    record_check(sprintf("bootstrap CI low [%s]", bt$group[i]), bt$ci_lo[i], s$ci_low, 0.02)
    record_check(sprintf("bootstrap CI high [%s]", bt$group[i]), bt$ci_hi[i], s$ci_high, 0.02)
  }

  sw <- rho_sweep(D); fwrite(sw, file.path(OUT, "rho_sensitivity.csv"))
  log_msg("rho sweep written; separation stays ",
          sprintf("%.1f-%.1f fold across rho in [0, 0.5]",
                  min(sw$separation_fold), max(sw$separation_fold)))

  log_step("4. Per-gene departures from proportional scaling")
  off <- departures(D, "1vs0")
  pw <- departure_power(D, "1vs0")
  log_msg(sprintf("departing at FDR<0.05: %d (propagated SE), %d (response SE alone), 31 reported",
                  pw$n_departing_propagated, pw$n_departing_response))
  log_msg(sprintf("smallest detectable departure %.3f log2; typical expectation tested %.3f log2",
                  pw$smallest_detectable, pw$typical_expectation))
  log_msg(sprintf("-> detection floor is %.1fx the typical effect: this count does NOT bound the off-axis component",
                  pw$floor_ratio))
  fwrite(off[q < CFG$params$fdr], file.path(OUT, "off_axis_departing_genes.csv"))

  log_step("5. Group gene signatures")
  sig <- rbindlist(lapply(names(CFG$params$contrasts), function(g) {
    d <- data.table(id = rownames(D$lfc), gene = unname(D$sym[rownames(D$lfc)]),
                    logFC = D$lfc[, g], padj = D$padj[, g])
    d <- d[is.finite(logFC)][order(-logFC)]
    rbind(head(d, 8)[, dir := "positive"], tail(d, 8)[, dir := "negative"])[
      , `:=`(contrast = g, group = CFG$params$contrasts[[g]],
             fdr_sig = padj < CFG$params$fdr)][]
  }))
  fwrite(sig, file.path(OUT, "group_gene_signatures.csv"))
  log_msg("wrote group_gene_signatures.csv (", nrow(sig), " rows)")

  res <- list(D = D, counts = cnt, decomposition = dec, bootstrap = bt,
              boot_diff = dd, rho_sweep = sw, departures = off, power = pw,
              signatures = sig, n_axis = n_axis)

  if (with_fgsea) {
    log_step("6. FGSEA recompute and fidelity against the distributed tables")
    P <- build_pathways()
    log_msg("total gene sets: ", length(P$sets))
    gs <- tryCatch(load_shipped_gsea(), error = function(e) NULL)
    rec <- rbindlist(lapply(fgsea_contrasts, function(g) {
      log_msg("fgseaMultilevel on ", g, " ...")
      r <- run_fgsea_contrast(D, g, P)
      log_msg(sprintf("  %s: %d ranked genes, %d pathways tested, %d at FDR<0.05",
                      g, r$n_ranked[1], nrow(r), sum(r$padj < CFG$params$fdr, na.rm = TRUE)))
      r
    }))
    fwrite(rec, file.path(OUT, "fgsea_recomputed.csv"))
    res$fgsea <- rec; res$pathways <- P
    if (!is.null(gs)) {   # optional fidelity against the originally distributed FGSEA tables
      fid <- rbindlist(lapply(fgsea_contrasts, function(g) fidelity(rec, gs, g)))
      print(fid); fwrite(fid, file.path(OUT, "fgsea_fidelity.csv"))
      log_msg(sprintf("NES rho: median %.3f, min %.3f | FDR Jaccard: median %.3f, min %.3f",
                      median(fid$NES_rho), min(fid$NES_rho),
                      median(fid$fdr_jaccard), min(fid$fdr_jaccard)))
      part <- pathway_partition(gs); fwrite(part, file.path(OUT, "pathway_partition.csv"))
      print(part)
      res$fidelity <- fid; res$partition <- part
    }
  }

  log_step("Reproduction checks")
  ck <- checks_table()
  print(ck)
  fwrite(ck, file.path(OUT, "reproduction_checks.csv"))
  log_msg(sprintf("%d of %d checks passed", sum(ck$pass), nrow(ck)),
          level = if (all(ck$pass)) "OK" else "WARN")
  res$checks <- ck
  res
}

#' Write the analysis log to disk with provenance and the data-tier table.
write_analysis_log <- function(res, path = file.path(cfg_out("logs"),
                                                     sprintf("ANALYSIS_LOG_%s.md", Sys.Date()))) {
  ck <- res$checks
  tiers <- data_tiers()
  fmt_tbl <- function(dt) {
    if (!nrow(dt)) return("_(empty)_")
    hdr <- paste("|", paste(names(dt), collapse = " | "), "|")
    sep <- paste("|", paste(rep("---", ncol(dt)), collapse = " | "), "|")
    rows <- apply(dt, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))
    paste(c(hdr, sep, rows), collapse = "\n")
  }
  txt <- c(
    sprintf("# Analysis log -- %s", Sys.Date()),
    "",
    sprintf("Run started %s, elapsed %.1f s.",
            format(.LOG$t0, "%Y-%m-%d %H:%M:%S"), log_elapsed()),
    sprintf("Case definition `%s`, model `%s` (%s).",
            CFG$meta$case_definition, CFG$meta$model,
            paste(CFG$params$covariates_used, collapse = ", ")),
    "",
    "## Data availability", "",
    fmt_tbl(tiers[, .(tier, item, exists, size_mb)]),
    "",
    "All results below are computed from the gene-level summary statistics.",
    "",
    "## Reproduction checks", "",
    fmt_tbl(ck[, .(check, recomputed = round(recomputed, 4),
                   shipped = round(shipped, 4), abs_diff = signif(abs_diff, 3),
                   tol, pass)]),
    "",
    sprintf("**%d of %d checks passed.**", sum(ck$pass), nrow(ck)),
    "",
    "## Run log", "", "```", .LOG$lines, "```", "",
    "## Session", "", "```", session_block(), "```")
  writeLines(txt, path)
  log_msg("wrote ", path)
  path
}
