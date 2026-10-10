# =============================================================================
# 00_config.R -- config resolution and the run log
# Every path used anywhere in this project resolves through cfg_path().
# No absolute path appears outside config.yml.
# =============================================================================
suppressPackageStartupMessages({
  library(yaml); library(data.table); library(jsonlite)
})

APP_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NULL)
if (is.null(APP_DIR) || !nzchar(APP_DIR)) APP_DIR <- getwd()
if (basename(APP_DIR) == "R") APP_DIR <- dirname(APP_DIR)

CFG <- yaml::read_yaml(file.path(APP_DIR, "config.yml"))

#' Resolve a {placeholder} path from config.
#' Placeholders draw from roots, previously-resolved input keys, meta, and outputs.
cfg_path <- function(key, section = "inputs") {
  raw <- CFG[[section]][[key]]
  if (is.null(raw)) stop("no config entry: ", section, "/", key)
  vars <- c(
    CFG$roots,
    list(
      model = CFG$meta$model,
      cd    = CFG$meta$case_definition,
      extracted_dir = NULL,
      outputs_dir   = NULL
    )
  )
  # two dependent keys, resolved first
  vars$extracted_dir <- gsub("\\{schema\\}", CFG$roots$schema, CFG$inputs$extracted_dir)
  vars$outputs_dir   <- gsub("\\{schema\\}", CFG$roots$schema, CFG$outputs$dir)
  out <- raw
  for (nm in names(vars)) {
    if (!is.null(vars[[nm]]))
      out <- gsub(sprintf("\\{%s\\}", nm), vars[[nm]], out, fixed = FALSE)
  }
  out
}

cfg_out <- function(key) cfg_path(key, "outputs")

#' Input availability for the analysis log.
data_tiers <- function() {
  rbindlist(lapply(c("gene_de", "gsea", "singscore"), function(k) {
    p <- cfg_path(k)
    data.table(tier = "summary statistics", item = k, path = p, exists = file.exists(p),
               size_mb = if (file.exists(p)) round(file.size(p) / 1e6, 1) else NA_real_)
  }))
}

# ------------------------------------------------------------------ run log --
.LOG <- new.env(parent = emptyenv())
.LOG$lines <- character(0)
.LOG$t0 <- Sys.time()

log_reset <- function() { .LOG$lines <- character(0); .LOG$t0 <- Sys.time(); invisible(NULL) }

log_msg <- function(..., level = "INFO") {
  msg <- sprintf("%-19s  %-5s  %s",
                 format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, paste0(...))
  .LOG$lines <- c(.LOG$lines, msg)
  cat(msg, "\n")
  invisible(msg)
}

log_step <- function(title) {
  bar <- strrep("-", 76)
  .LOG$lines <- c(.LOG$lines, "", bar, paste("##", title), bar)
  cat("\n", bar, "\n## ", title, "\n", bar, "\n", sep = "")
  invisible(NULL)
}

#' Record a check of a recomputed value against a previously reported value.
.CHECKS <- new.env(parent = emptyenv())
.CHECKS$rows <- list()
check_reset <- function() .CHECKS$rows <- list()
record_check <- function(what, recomputed, shipped, tol, units = "") {
  ok <- is.finite(recomputed) && is.finite(shipped) && abs(recomputed - shipped) <= tol
  .CHECKS$rows[[length(.CHECKS$rows) + 1]] <- data.table(
    check = what, recomputed = recomputed, shipped = shipped,
    abs_diff = abs(recomputed - shipped), tol = tol, units = units, pass = ok)
  log_msg(sprintf("%-46s recomputed %9.4f | shipped %9.4f | tol %.4f | %s",
                  what, recomputed, shipped, tol, if (ok) "PASS" else "FAIL"),
          level = if (ok) "CHECK" else "FAIL")
  ok
}
checks_table <- function() if (length(.CHECKS$rows)) rbindlist(.CHECKS$rows) else data.table()

log_text <- function() paste(.LOG$lines, collapse = "\n")
log_elapsed <- function() as.numeric(difftime(Sys.time(), .LOG$t0, units = "secs"))

session_block <- function() {
  si <- sessionInfo()
  pk <- sort(vapply(si$otherPkgs, function(p) paste0(p$Package, " ", p$Version), ""))
  c(paste("R version:", si$R.version$version.string),
    paste("Platform:", si$platform),
    paste("Running under:", si$running),
    "Attached packages:", paste0("  ", pk))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
