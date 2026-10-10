#!/usr/bin/env Rscript

# Render the pathway-organization and rank-test robustness audit from the
# tables produced by 04_pathway_evidence_engine.R. No enrichment is recomputed.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]))
root <- dirname(dirname(script_file))
table_dir <- file.path(root, "tables")
figure_dir <- file.path(root, "figures")

stability_file <- file.path(table_dir, "igraph_cluster_stability_model1.csv")
concordance_file <- file.path(table_dir, "pathway_method_concordance_model1.csv")
stopifnot(file.exists(stability_file), file.exists(concordance_file))

stability <- fread(stability_file)
counts <- fread(concordance_file)[evidence_class != "Neither"]
stability[, community_label_y := fcase(
  overlap_coefficient_threshold == 0.45, 0.912,
  overlap_coefficient_threshold == 0.50, 1.010,
  overlap_coefficient_threshold == 0.55, 0.978
)]

contrast_order <- c("Newly diagnosed", "Potentially excluded", "Established COPD")
evidence_order <- c("Both, concordant", "FGSEA only", "Competitive test only")
evidence_labels <- c(
  "Both, concordant" = "Both tests, same direction",
  "FGSEA only" = "FGSEA only",
  "Competitive test only" = "Rank test only"
)
evidence_colours <- c(
  "Both tests, same direction" = "#2A9D8F",
  "FGSEA only" = "#D55E00",
  "Rank test only" = "#0072B2"
)

counts[, contrast_label := factor(contrast_label, levels = contrast_order)]
counts[, evidence_label := factor(
  unname(evidence_labels[evidence_class]),
  levels = unname(evidence_labels[evidence_order])
)]
totals <- counts[, .(total = sum(N)), by = contrast_label]

p_stability <- ggplot(
  stability,
  aes(overlap_coefficient_threshold, adjusted_rand_vs_0_50)
) +
  geom_hline(yintercept = 1, colour = "grey70", linewidth = 0.4, linetype = 3) +
  geom_vline(xintercept = 0.50, colour = "grey70", linewidth = 0.4, linetype = 3) +
  geom_line(linewidth = 0.7, colour = "#2166AC") +
  geom_point(shape = 21, size = 4.2, fill = "white", stroke = 0.9, colour = "#2166AC") +
  geom_text(
    aes(
      y = community_label_y,
      label = paste0(communities, " communities")
    ),
    hjust = 0.5,
    vjust = 0,
    size = 3.1
  ) +
  scale_x_continuous(
    limits = c(0.44, 0.56),
    breaks = c(0.45, 0.50, 0.55),
    labels = c("0.45\nLower threshold", "0.50\nPrimary threshold", "0.55\nHigher threshold"),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(0.875, 1.025),
    breaks = c(0.85, 0.90, 0.95, 1.00)
  ) +
  labs(
    title = "A  Stability of pathway clusters",
    subtitle = "Adjusted Rand index near 1 indicates similar pathway membership",
    x = "Full-membership overlap threshold",
    y = "Agreement with primary 0.50 partition"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, colour = "grey30", margin = margin(b = 8)),
    axis.title = element_text(size = 9.2),
    axis.text = element_text(size = 8.5),
    plot.margin = margin(10, 16, 26, 12)
  )

p_concordance <- ggplot(counts, aes(contrast_label, N, fill = evidence_label)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.3) +
  geom_text(
    aes(label = ifelse(N >= 8, N, "")),
    position = position_stack(vjust = 0.5),
    colour = "white",
    fontface = "bold",
    size = 3.2
  ) +
  geom_text(
    data = totals,
    aes(x = contrast_label, y = total, label = paste0("Total ", total)),
    inherit.aes = FALSE,
    vjust = -0.45,
    size = 3.1,
    fontface = "bold"
  ) +
  scale_fill_manual(values = evidence_colours, name = "FDR < 0.05 evidence") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "B  Agreement between pathway tests",
    subtitle = "Candidates detected by at least one test; pathways detected by neither are omitted",
    x = NULL,
    y = "Candidate pathways"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, colour = "grey30", margin = margin(b = 8)),
    axis.title.y = element_text(size = 9.2),
    axis.text.x = element_text(size = 8.7, face = "bold"),
    axis.text.y = element_text(size = 8.5),
    legend.position = "bottom",
    legend.title = element_text(size = 8.8, face = "bold"),
    legend.text = element_text(size = 8.2),
    legend.key.width = grid::unit(9, "mm"),
    legend.spacing.x = grid::unit(3, "mm"),
    plot.margin = margin(14, 16, 8, 12)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "top"))

draw_audit <- function() {
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 2,
    ncol = 1,
    heights = grid::unit(c(1.02, 1), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  print(
    p_stability,
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1),
    newpage = FALSE
  )
  print(
    p_concordance,
    vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1),
    newpage = FALSE
  )
  grid::popViewport()
}

grDevices::pdf(
  file.path(figure_dir, "FigS2_pathway_evidence_audit.pdf"),
  width = 7.25,
  height = 7.5,
  useDingbats = FALSE
)
draw_audit()
grDevices::dev.off()

ragg::agg_png(
  file.path(figure_dir, "FigS2_pathway_evidence_audit.png"),
  width = 7.25,
  height = 7.5,
  units = "in",
  res = 600,
  background = "white"
)
draw_audit()
grDevices::dev.off()

cat("Figure S2 rebuilt from pathway stability and test-concordance audit tables.\n")

