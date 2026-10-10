#!/usr/bin/env Rscript
# Disease-axis scaling and FGSEA recompute from gene-level summary statistics.
# Usage: Rscript run_reproduction.R [--no-fgsea]
# Edit config.yml first. Writes <outputs>/tables/*.csv and <outputs>/logs/ANALYSIS_LOG_<date>.md
setwd(if (nzchar(d <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))) d else ".")
for (f in sort(list.files(".", "^0[0-9]_.*\\.R$", full.names = TRUE))) source(f)
for (k in c("tables", "figures", "logs")) dir.create(cfg_out(k), showWarnings = FALSE, recursive = TRUE)
args <- commandArgs(TRUE)
res <- reproduce_all(with_fgsea = !("--no-fgsea" %in% args))
write_analysis_log(res)
saveRDS(res[setdiff(names(res), c("pathways"))], file.path(cfg_out("tables"), "reproduction.rds"))
quit(status = if (all(res$checks$pass)) 0 else 1)
