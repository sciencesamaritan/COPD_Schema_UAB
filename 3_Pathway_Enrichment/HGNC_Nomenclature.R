# Shared BioBouncer helpers for figure-facing HGNC nomenclature.

BIOBOUNCER_VERSION <- "0.1.4"
HGNC_SNAPSHOT_VERSION <- "2026-07-07"

require_biobouncer <- function() {
  if (!requireNamespace("biobouncer", quietly = TRUE)) {
    stop("biobouncer is required for downstream HGNC nomenclature checks")
  }
  installed <- as.character(utils::packageVersion("biobouncer"))
  if (!identical(installed, BIOBOUNCER_VERSION)) {
    stop("Expected biobouncer ", BIOBOUNCER_VERSION, "; found ", installed)
  }
}

audit_hgnc_symbols <- function(symbols, scope, figure) {
  require_biobouncer()
  symbols <- sort(unique(trimws(as.character(symbols))))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  audit <- data.table::as.data.table(biobouncer::check_id(
    symbols,
    source_db = "hgnc",
    how = "cache",
    version = HGNC_SNAPSHOT_VERSION
  ))
  audit[, `:=`(
    status = data.table::fifelse(
      valid %in% TRUE,
      "current",
      data.table::fifelse(!is.na(suggestion), "repairable", "unrecognised")
    ),
    approved_symbol = data.table::fifelse(
      valid %in% TRUE,
      normalized,
      data.table::fifelse(!is.na(suggestion), suggestion, input)
    ),
    scope = scope,
    figure = figure,
    biobouncer_version = BIOBOUNCER_VERSION
  )]
  audit
}

normalize_hgnc_symbols <- function(symbols, audit) {
  lookup <- setNames(audit$approved_symbol, audit$input)
  normalized <- unname(lookup[as.character(symbols)])
  normalized[is.na(normalized)] <- as.character(symbols)[is.na(normalized)]
  normalized
}

split_leading_edge <- function(values) {
  values <- values[!is.na(values) & nzchar(values)]
  sort(unique(trimws(unlist(strsplit(values, ";", fixed = TRUE)))))
}

normalize_leading_edges <- function(values, audit, drop_unrecognised = FALSE) {
  if (drop_unrecognised) {
    keep <- audit$status != "unrecognised"
    audit <- audit[keep]
  }
  vapply(values, function(value) {
    if (is.na(value) || !nzchar(value)) return("")
    genes <- trimws(strsplit(value, ";", fixed = TRUE)[[1]])
    normalized <- normalize_hgnc_symbols(genes, audit)
    if (drop_unrecognised) normalized <- normalized[genes %in% audit$input]
    paste(unique(normalized[nzchar(normalized)]), collapse = ";")
  }, character(1))
}

summarize_hgnc_audit <- function(audit) {
  audit[, .(
    total_unique_labels = .N,
    current = sum(status == "current"),
    repairable = sum(status == "repairable"),
    unrecognised = sum(status == "unrecognised"),
    all_resolved_to_approved_hgnc = all(status != "unrecognised"),
    hgnc_snapshot = unique(version),
    biobouncer_version = unique(biobouncer_version)
  ), by = .(figure, scope)]
}

