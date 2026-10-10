#!/usr/bin/env Rscript

# Direct Model 1 versus Model 2 sensitivity display from exported gene effects
# and paired-bootstrap λ-difference summaries.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(hexbin)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]))
root <- dirname(dirname(script_file))
figure_dir <- file.path(root, "figures")
table_dir <- file.path(root, "tables")

analysis_root <- paste0(
  "path/to/",
  "Schema_omics_csv_bhatt"
)
de_dir <- file.path(analysis_root, "copdschema", "extracted")
m1_file <- file.path(de_dir, "gene_de_Model_1_DisCOPDAI_allresults_2026.04.25.csv")
m2_file <- file.path(de_dir, "gene_de_Model_2_DisCOPDAI_allresults_2026.04.25.csv")
lambda_file <- file.path(
  analysis_root, "reproduce_2026-09-01", "outputs", "tables",
  "model1_vs_model2_all_classifications.csv"
)

contrasts <- c("DisCOPDAI_1vs0", "DisCOPDAI_2vs0", "DisCOPDAI_3vs0")
contrast_titles <- c(
  DisCOPDAI_1vs0 = "A  Newly diagnosed vs no COPD",
  DisCOPDAI_2vs0 = "B  Potentially excluded vs no COPD",
  DisCOPDAI_3vs0 = "C  Established COPD vs no COPD"
)

m1 <- fread(m1_file)[Contrast %chin% contrasts, .(
  id, gene = gene_label, Contrast,
  logFC_model1 = logFC,
  FDR_model1 = adj.P.Val
)]
m2 <- fread(m2_file)[Contrast %chin% contrasts, .(
  id, Contrast,
  logFC_model2 = logFC,
  FDR_model2 = adj.P.Val
)]
effects <- merge(m1, m2, by = c("id", "Contrast"))
stopifnot(nrow(effects) == 16554L * 3L)

effect_summary <- effects[, .(
  genes = .N,
  spearman_rho = cor(logFC_model1, logFC_model2, method = "spearman"),
  pearson_r = cor(logFC_model1, logFC_model2),
  model1_fdr_05 = sum(FDR_model1 < 0.05),
  model2_fdr_05 = sum(FDR_model2 < 0.05),
  significant_both = sum(FDR_model1 < 0.05 & FDR_model2 < 0.05)
), by = Contrast]
fwrite(effect_summary, file.path(table_dir, "model1_model2_gene_effect_concordance.csv"))

common_limit <- ceiling(max(abs(c(effects$logFC_model1, effects$logFC_model2))) * 10) / 10

effect_panel <- function(contrast) {
  x <- effects[Contrast == contrast]
  s <- effect_summary[Contrast == contrast]
  ggplot(x, aes(logFC_model1, logFC_model2)) +
    geom_hex(bins = 48) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#4D4D4D", linewidth = 0.45) +
    annotate(
      "text", x = -1.32, y = 1.28,
      label = sprintf("'Spearman'~rho == %.2f", s$spearman_rho),
      parse = TRUE,
      hjust = 0, vjust = 1, size = 2.6
    ) +
    annotate(
      "text", x = -1.32, y = 1.12,
      label = sprintf(
        "FDR < 0.05: M1 %s; M2 %s",
        format(s$model1_fdr_05, big.mark = ","),
        format(s$model2_fdr_05, big.mark = ",")
      ),
      hjust = 0, vjust = 1, size = 2.6
    ) +
    scale_fill_gradient(
      low = "#F2F2F2", high = "#2166AC",
      trans = "log10", name = "Genes per\nhexagon"
    ) +
    coord_equal(xlim = c(-common_limit, common_limit), ylim = c(-common_limit, common_limit), expand = FALSE) +
    labs(
      title = contrast_titles[[contrast]],
      x = expression("Model 1 gene-expression " * log[2] * " fold change"),
      y = expression("Model 2 gene-expression " * log[2] * " fold change")
    ) +
    theme_classic(base_size = 8.5) +
    theme(
      plot.title = element_text(face = "bold", size = 9),
      axis.title = element_text(size = 8),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 7),
      aspect.ratio = 1,
      plot.margin = margin(6, 6, 6, 6)
    )
}

lambda <- fread(lambda_file)[case_def == "DisCOPDAI"]
lambda <- lambda[, .(
  diff = first(diff),
  diff_lo = first(diff_lo),
  diff_hi = first(diff_hi),
  axis_n = first(axis_n),
  lambda_new = lambda[group == "Newly diagnosed"][[1]],
  lambda_excluded = lambda[group == "Potentially excluded"][[1]]
), by = model]
lambda[, model_label := factor(model, levels = c("Model_2", "Model_1"), labels = c("Model 2", "Model 1"))]
fwrite(lambda, file.path(table_dir, "model1_model2_lambda_difference.csv"))

p_lambda <- ggplot(lambda, aes(diff, model_label)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.45) +
  geom_errorbar(
    aes(xmin = diff_lo, xmax = diff_hi),
    orientation = "y", width = 0.16, linewidth = 0.65
  ) +
  geom_point(shape = 21, size = 3.2, stroke = 0.65, fill = "white") +
  geom_text(
    aes(
      x = diff,
      label = sprintf("%.3f (%.3f to %.3f)", diff, diff_lo, diff_hi)
    ),
    vjust = -1.35, size = 2.55
  ) +
  annotate(
    "text", x = 0.01, y = 1.48,
    label = sprintf(
      "'Model 1'~lambda*':'~'%.3f, %.3f (axis genes = %s)'",
      lambda[model == "Model_1"]$lambda_new,
      lambda[model == "Model_1"]$lambda_excluded,
      format(lambda[model == "Model_1"]$axis_n, big.mark = ",")
    ),
    parse = TRUE, hjust = 0, size = 2.45
  ) +
  annotate(
    "text", x = 0.01, y = 0.72,
    label = sprintf(
      "'Model 2'~lambda*':'~'%.3f, %.3f (axis genes = %s)'",
      lambda[model == "Model_2"]$lambda_new,
      lambda[model == "Model_2"]$lambda_excluded,
      format(lambda[model == "Model_2"]$axis_n, big.mark = ",")
    ),
    parse = TRUE, hjust = 0, size = 2.45
  ) +
  scale_x_continuous(limits = c(0, 0.57), breaks = seq(0, 0.5, 0.1)) +
  labs(
    title = "D  Difference in established-COPD-axis position",
    subtitle = "Newly diagnosed minus potentially excluded\nPaired gene bootstrap",
    x = expression("Difference in " * lambda * " (95% CI)"),
    y = NULL
  ) +
  theme_classic(base_size = 8.5) +
  theme(
    plot.title = element_text(face = "bold", size = 9),
    plot.subtitle = element_text(size = 7.5),
    axis.text.y = element_text(face = "bold"),
    aspect.ratio = 1,
    plot.margin = margin(6, 6, 6, 6)
  )

figure <- wrap_plots(
  effect_panel(contrasts[[1]]), effect_panel(contrasts[[2]]),
  effect_panel(contrasts[[3]]), p_lambda,
  ncol = 2,
  widths = c(1, 1),
  heights = c(1, 1),
  guides = "collect"
) &
  theme(legend.position = "bottom")

ggsave(
  file.path(figure_dir, "FigS1_model2_sensitivity.pdf"),
  figure, width = 7.25, height = 7.0, device = "pdf"
)
ggsave(
  file.path(figure_dir, "FigS1_model2_sensitivity.png"),
  figure, width = 7.25, height = 7.0, dpi = 600, bg = "white"
)

cat(sprintf(
  "Model 2 sensitivity figure complete: %s matched gene effects across 3 reference contrasts.\n",
  format(nrow(effects) / 3, big.mark = ",")
))

