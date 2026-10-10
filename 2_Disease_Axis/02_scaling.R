# =============================================================================
# 02_scaling.R -- effect-size scaling along the established-disease axis
#
# Ports the estimator from the submission notebook (06_notebook/build_notebook.py)
# to R, with both bias corrections intact:
#   1. regression dilution  -- subtract mean error variance from the denominator
#   2. shared-reference     -- subtract rho * E[sx*sy] from the numerator, with
#                              rho determined by the design (group sizes)
# =============================================================================

#' Split "1vs0" into its two group labels.
.contrast_groups <- function(g) {
  p <- strsplit(g, "vs", fixed = TRUE)[[1]]
  if (length(p) != 2L) stop("cannot parse contrast: ", g)
  p
}

#' Design-implied correlation between two contrast estimators.
#'
#' For contrasts (a - b) and (c - d) under a common-variance linear model,
#'   Cov = sigma^2 * [ [a==c]/n_a - [a==d]/n_a - [b==c]/n_b + [b==d]/n_b ]
#'   Var = sigma^2 * ( 1/n_a + 1/n_b )
#' so sigma^2 cancels in the correlation. The earlier implementation read only
#' the first group of each contrast and assumed the second was always the shared
#' reference, which returned 0.094 rather than 1 for a contrast against itself
#' and a non-zero value for contrasts sharing no group at all.
rho_design <- function(g, ref = CFG$params$reference_contrast, n = CFG$params$group_n) {
  n  <- lapply(n, as.numeric)
  ab <- .contrast_groups(g); cd <- .contrast_groups(ref)
  a <- ab[1]; b <- ab[2]; c_ <- cd[1]; d <- cd[2]
  for (k in c(a, b, c_, d))
    if (is.null(n[[k]])) stop("no group size for level '", k, "' in config params$group_n")
  cov <- (a == c_) / n[[a]] - (a == d) / n[[a]] - (b == c_) / n[[b]] + (b == d) / n[[b]]
  v1  <- 1 / n[[a]] + 1 / n[[b]]
  v2  <- 1 / n[[c_]] + 1 / n[[d]]
  cov / sqrt(v1 * v2)
}

#' Restrict to the reference (established-disease) signature: genes at FDR<0.05.
axis_mask <- function(D, ref = CFG$params$reference_contrast, fdr = CFG$params$fdr) {
  which(D$padj[, ref] < fdr)
}

#' Scaling coefficient of group `g`'s effect profile onto the `ref` axis.
lambda_hat <- function(D, g, ref = CFG$params$reference_contrast,
                       rho = NULL, idx = NULL) {
  m <- axis_mask(D, ref)
  y <- D$lfc[m, g]; x <- D$lfc[m, ref]
  sy <- D$se[m, g]; sx <- D$se[m, ref]
  if (!is.null(idx)) { x <- x[idx]; y <- y[idx]; sx <- sx[idx]; sy <- sy[idx] }
  r <- if (is.null(rho)) rho_design(g, ref) else rho
  # population moments (denominator matches numpy var(), ddof=0)
  num <- mean(x * y) - mean(x) * mean(y) - r * mean(sx * sy)
  den <- mean(x^2) - mean(x)^2 - mean(sx^2)
  # The denominator is the predictor's variance net of its measurement error.
  # If error exceeds observed spread it goes non-positive and lambda is
  # meaningless rather than merely imprecise, so refuse instead of returning it.
  if (!is.finite(den) || den <= 0)
    stop(sprintf("lambda_hat: predictor variance net of measurement error is %.3g; ",
                 den), "the axis carries less signal than noise, so the ratio is undefined")
  num / den
}

#' Split systematic profile variance into on-axis and off-axis components.
#' Off-axis is (total systematic) - (on-axis), NOT the raw residual variance:
#' the raw residual still carries the predictor's measurement error scaled by
#' lambda, which inflates the off-axis term.
decompose <- function(D, g, ref = CFG$params$reference_contrast) {
  m <- axis_mask(D, ref)
  y <- D$lfc[m, g]; x <- D$lfc[m, ref]
  sy <- D$se[m, g]; sx <- D$se[m, ref]
  l <- lambda_hat(D, g, ref)
  total_sys <- (mean(y^2) - mean(y)^2) - mean(sy^2)
  pred_sys  <- (mean(x^2) - mean(x)^2) - mean(sx^2)
  on_var  <- max(l^2 * pred_sys, 0)
  off_var <- max(total_sys - on_var, 0)
  tot <- on_var + off_var
  data.table(contrast = g, group = CFG$params$contrasts[[g]],
             lambda = l, rho = rho_design(g, ref),
             on_axis_sd = sqrt(on_var), off_axis_sd = sqrt(off_var),
             # both components vanish only if the group has no systematic
             # variance at all; report NA rather than 0/0 = NaN
             pct_variance_on_axis = if (tot > 0) 100 * on_var / tot else NA_real_)
}

#' Paired bootstrap over the genes of the reference signature. Both coefficients
#' are recomputed on each resample, so the shared reference and shared axis enter
#' the null exactly as they do in the data.
bootstrap_lambda <- function(D, groups = c("1vs0", "2vs0"),
                             B = CFG$params$bootstrap_B, seed = CFG$params$seed) {
  # Seeding is deliberate, but the caller's stream must not be disturbed: any
  # random draw after this call would otherwise silently depend on it.
  if (exists(".Random.seed", envir = globalenv())) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  n_g <- length(axis_mask(D))
  out <- matrix(NA_real_, B, length(groups), dimnames = list(NULL, groups))
  for (b in seq_len(B)) {
    idx <- sample.int(n_g, n_g, replace = TRUE)
    for (g in groups) out[b, g] <- lambda_hat(D, g, idx = idx)
  }
  if (anyNA(out)) stop("bootstrap produced ", sum(is.na(out)), " non-finite replicates")
  ci <- t(apply(out, 2, quantile, c(.025, .975), names = FALSE))
  res <- data.table(contrast = groups,
                    group = unlist(CFG$params$contrasts[groups]),
                    lambda = vapply(groups, function(g) lambda_hat(D, g), 0),
                    ci_lo = ci[, 1], ci_hi = ci[, 2],
                    boot_sd = apply(out, 2, sd))
  d <- out[, 1] - out[, 2]
  attr(res, "diff") <- list(
    observed = res$lambda[1] - res$lambda[2], se = sd(d),
    ci = quantile(d, c(.025, .975), names = FALSE),
    n_le_zero = sum(d <= 0), B = B,
    boot_cor = cor(out[, 1], out[, 2]),
    fold = res$lambda[1] / res$lambda[2])
  attr(res, "draws") <- out
  res
}

#' Sensitivity of lambda to the assumed shared-reference correlation.
rho_sweep <- function(D, groups = c("1vs0", "2vs0"), rhos = seq(0, 0.5, length.out = 11)) {
  rbindlist(lapply(rhos, function(r) {
    l <- vapply(groups, function(g) lambda_hat(D, g, rho = r), 0)
    data.table(rho = r, contrast = groups, group = unlist(CFG$params$contrasts[groups]),
               lambda = l)
  }))[, separation_fold := lambda[contrast == "1vs0"] / lambda[contrast == "2vs0"], by = rho][]
}

#' Per-gene test of the observed effect against its lambda-scaled expectation.
#' propagate=TRUE uses the error of the difference (response error plus the
#' lambda-scaled predictor error); FALSE uses the response error alone.
departures <- function(D, g = "1vs0", ref = CFG$params$reference_contrast,
                       propagate = TRUE) {
  m <- axis_mask(D, ref)
  ids <- rownames(D$lfc)[m]
  y <- D$lfc[m, g]; x <- D$lfc[m, ref]
  sy <- D$se[m, g]; sx <- D$se[m, ref]
  l <- lambda_hat(D, g, ref)
  dev <- y - l * x
  sd_dev <- if (propagate) sqrt(sy^2 + (l * sx)^2) else sy
  z <- dev / sd_dev
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  out <- data.table(id = ids, gene = unname(D$sym[ids]), logFC = y,
                    expected = l * x, deviation = dev, z = z, p = p,
                    q = p.adjust(p, "BH"))
  setattr(out, "lambda", l)
  setattr(out, "sd_dev", sd_dev)
  out[order(-abs(deviation))]
}

#' The detection floor: what departure size this test could actually see,
#' relative to the effect it is testing against.
departure_power <- function(D, g = "1vs0", ref = CFG$params$reference_contrast) {
  off <- departures(D, g, ref, propagate = TRUE)
  sd_dev <- attr(off, "sd_dev"); l <- attr(off, "lambda")
  m <- axis_mask(D, ref)
  smallest <- 1.96 * median(sd_dev)
  typical <- median(abs(l * D$lfc[m, ref]))
  list(n_departing_propagated = sum(off$q < CFG$params$fdr),
       n_departing_response   = sum(departures(D, g, ref, propagate = FALSE)$q < CFG$params$fdr),
       n_axis_genes = length(m),
       smallest_detectable = smallest, typical_expectation = typical,
       floor_ratio = smallest / typical)
}
