#!/usr/bin/env Rscript

# Build Figure 4B from the archived participant-level Reactome singscore
# contrasts, then combine it with the ranked enrichment panel.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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
singscore_file <- file.path(
  analysis_root, "copdschema", "extracted",
  "singscore_Model_1_DisCOPDAI_allresults_2026.04.26.csv"
)
integrated_file <- file.path(table_dir, "pathway_integrated_evidence_model1.csv")
panel_a_file <- file.path(table_dir, "Fig4_panelA_pathway_plot.rds")

stopifnot(
  file.exists(singscore_file),
  file.exists(integrated_file),
  file.exists(panel_a_file)
)

contrasts <- c("DisCOPDAI_1vs0", "DisCOPDAI_2vs0", "DisCOPDAI_3vs0")
contrast_labels <- c(
  DisCOPDAI_1vs0 = "Newly\ndiagnosed\nvs no COPD",
  DisCOPDAI_2vs0 = "Potentially\nexcluded\nvs no COPD",
  DisCOPDAI_3vs0 = "Established\nCOPD\nvs no COPD"
)

integrated <- fread(integrated_file)
selected_order <- unique(integrated$pathway)
selected_meta <- unique(integrated[, .(pathway, evidence_pattern)])

scores_all <- fread(singscore_file)[
  Model == "Model_1" & Contrast %chin% contrasts & source == "Reactome"
]
scores <- merge(
  scores_all[pathway %chin% selected_order],
  selected_meta,
  by = "pathway",
  all.x = TRUE
)
selected_reactome <- selected_order[grepl("^REACTOME_", selected_order)]
stopifnot(
  uniqueN(scores$pathway) == length(selected_reactome),
  nrow(scores) == length(selected_reactome) * length(contrasts)
)

clean_label <- function(pathway) {
  x <- sub("^(HALLMARK|REACTOME)_", "", pathway)
  x <- tools::toTitleCase(tolower(gsub("_", " ", x, fixed = TRUE)))
  replacements <- c(
    "Dna" = "DNA", "Rrna" = "rRNA", "Rna" = "RNA", "Mrna" = "mRNA",
    "Ecm" = "ECM", "Nfkb" = "NF-kB", "Tnfa" = "TNF-alpha",
    "Slc2a4" = "SLC2A4", "Glut4" = "GLUT4", "Mtorc1" = "mTORC1",
    "Iii" = "III", "Uv" = "UV"
  )
  for (from in names(replacements)) x <- gsub(from, replacements[[from]], x, fixed = TRUE)
  x
}

score_order <- selected_reactome
scores[, `:=`(
  pathway_label = factor(clean_label(pathway), levels = rev(clean_label(score_order))),
  contrast_label = factor(contrast_labels[Contrast], levels = unname(contrast_labels[contrasts])),
  score_diff = logFC,
  FDR = adj.P.Val,
  ci_low = CI.L,
  ci_high = CI.R,
  direction = fifelse(logFC > 0, "Higher", fifelse(logFC < 0, "Lower", "No difference")),
  sig = fifelse(adj.P.Val < 0.05, "FDR < 0.05", "FDR >= 0.05")
)]
fwrite(
  scores[order(match(pathway, score_order), match(Contrast, contrasts))],
  file.path(table_dir, "singscore_selected_model1.csv")
)

limit <- max(abs(c(scores$ci_low, scores$ci_high)), na.rm = TRUE)
limit <- ceiling(limit * 1000) / 1000
contrast_colours <- setNames(
  c("#0072B2", "#D55E00", "#4D4D4D"),
  unname(contrast_labels[contrasts])
)
singscore_plot <- ggplot(
  scores,
  aes(score_diff, pathway_label, colour = contrast_label)
) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = "grey45") +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    width = 0,
    linewidth = 0.55
  ) +
  geom_point(
    aes(fill = sig),
    shape = 21,
    size = 2.5,
    stroke = 0.75
  ) +
  facet_grid(cols = vars(contrast_label)) +
  scale_colour_manual(values = contrast_colours, guide = "none") +
  scale_fill_manual(
    values = c("FDR < 0.05" = "black", "FDR >= 0.05" = "white"),
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(-limit, limit),
    breaks = scales::breaks_pretty(n = 3),
    expand = expansion(mult = c(0.06, 0.06))
  ) +
  scale_y_discrete(labels = function(x) {
    vapply(x, function(y) paste(strwrap(y, 38), collapse = "\n"), character(1))
  }) +
  labs(
    title = "Participant-level pathway-score differences",
    x = "Adjusted undirected pathway-score difference",
    y = NULL,
    caption = paste0(
      "Differences from No COPD Group; bars are 95% CIs. Filled points: ",
      "Benjamini-Hochberg FDR < 0.05."
    )
  ) +
  theme_classic(base_size = 8.5) +
  theme(
    plot.caption = element_text(size = 6.6, hjust = 0, colour = "grey25"),
    plot.caption.position = "plot",
    axis.title.x = element_text(size = 7.5, margin = margin(t = 5)),
    axis.text.x = element_text(size = 6.8),
    axis.text.y = element_text(size = 6.4, lineheight = 0.92),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    panel.spacing.x = grid::unit(3.5, "mm"),
    strip.background = element_rect(fill = "grey96", colour = "grey75", linewidth = 0.35),
    strip.text = element_text(face = "bold", size = 7.2, lineheight = 0.95),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    plot.title = element_text(size = 9, face = "bold", hjust = 0),
    plot.margin = margin(5, 8, 6, 8)
  ) +
  guides(fill = guide_legend(override.aes = list(colour = "black"), nrow = 1))

ggsave(
  file.path(figure_dir, "Fig4B_model1_singscore.pdf"),
  singscore_plot, width = 7.25, height = 4.8, device = "pdf"
)
ggsave(
  file.path(figure_dir, "Fig4B_model1_singscore.png"),
  singscore_plot, width = 7.25, height = 4.8, dpi = 600, bg = "white"
)

pathway_plot <- readRDS(panel_a_file)
pathway_plot <- pathway_plot +
  scale_x_discrete(
    position = "top",
    labels = c(
      "Newly diagnosed\nvs no COPD" = "Newly\ndiagnosed\nvs no COPD",
      "Potentially excluded\nvs no COPD" = "Potentially\nexcluded\nvs no COPD",
      "Established COPD\nvs no COPD" = "Established\nCOPD\nvs no COPD"
    )
  ) +
  labs(tag = "A") +
  theme(
    axis.text.x.top = element_text(face = "bold", size = 8.2, lineheight = 0.90),
    axis.text.y = element_text(size = 8.2, lineheight = 0.95),
    strip.text.y.right = element_text(
      angle = 0, face = "bold", size = 8, lineheight = 0.95
    ),
    legend.title = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 7.8),
    plot.title = element_text(face = "bold", size = 10.5, hjust = 0),
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = "topleft",
    plot.margin = margin(8, 7, 6, 8)
  )
singscore_plot <- singscore_plot +
  labs(tag = "B") +
  theme(
    plot.caption = element_text(size = 7.2, hjust = 0, colour = "grey20"),
    axis.title.x = element_text(size = 8.2, margin = margin(t = 5)),
    axis.text.x = element_text(size = 7.4),
    axis.text.y = element_text(size = 7.5, lineheight = 0.92),
    strip.text = element_text(face = "bold", size = 7.5, lineheight = 0.90),
    legend.text = element_text(size = 7.8),
    plot.title = element_text(size = 10.5, face = "bold", hjust = 0),
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = "topleft",
    plot.margin = margin(8, 8, 6, 7)
  )

draw_horizontal_figure <- function() {
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 1,
    ncol = 2,
    widths = grid::unit(c(1, 1), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  print(
    pathway_plot,
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1),
    newpage = FALSE
  )
  print(
    singscore_plot,
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2),
    newpage = FALSE
  )
  grid::popViewport()
}

grDevices::pdf(
  file.path(figure_dir, "Fig4_model1_pathway_programs.pdf"),
  width = 13.0,
  height = 7.2,
  useDingbats = FALSE
)
draw_horizontal_figure()
grDevices::dev.off()

grDevices::png(
  file.path(figure_dir, "Fig4_model1_pathway_programs.png"),
  width = 13.0,
  height = 7.2,
  units = "in",
  res = 600,
  bg = "white",
  type = "cairo"
)
draw_horizontal_figure()
grDevices::dev.off()

cat(
  "Figure 4 complete: panel A ranked enrichment; panel B archived Reactome ",
  "participant-level singscore contrasts saved separately and combined.\n",
  sep = ""
)

