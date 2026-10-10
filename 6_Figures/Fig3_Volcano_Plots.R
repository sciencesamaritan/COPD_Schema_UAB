#!/usr/bin/env Rscript

# Traditional three-panel limma volcano plot for the primary Model 1 contrasts.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]))
root <- dirname(dirname(script_file))
source(file.path(dirname(script_file), "..", "3_Pathway_Enrichment", "HGNC_Nomenclature.R"))

source_file <- paste0(
  "path/to/",
  "Schema_omics_csv_bhatt/copdschema/extracted/",
  "gene_de_Model_1_DisCOPDAI_allresults_2026.04.25.csv"
)
figure_stem <- file.path(root, "figures", "Fig3_model1_enhanced_volcanoes")
table_dir <- file.path(root, "tables")
dir.create(dirname(figure_stem), recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

contrast_labels <- c(
  DisCOPDAI_1vs0 = "Newly diagnosed vs no COPD",
  DisCOPDAI_2vs0 = "Potentially excluded vs no COPD",
  DisCOPDAI_3vs0 = "Established COPD vs no COPD"
)

# Fixed display limits keep the three reference contrasts directly comparable.
# Values beyond a limit are retained at the boundary and marked with triangles.
x_limit <- 1
y_limit <- 25

de <- read.csv(source_file, check.names = FALSE) |>
  filter(Contrast %in% names(contrast_labels)) |>
  mutate(
    contrast = factor(contrast_labels[Contrast], levels = unname(contrast_labels)),
    fdr = pmin(pmax(adj.P.Val, .Machine$double.xmin), 1),
    nominal_p = pmin(pmax(P.Value, .Machine$double.xmin), 1),
    minus_log10_p = -log10(nominal_p),
    display_logFC = pmax(-x_limit, pmin(x_limit, logFC)),
    display_y = pmin(y_limit, minus_log10_p),
    capped = abs(logFC) > x_limit | minus_log10_p > y_limit,
    status = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Higher expression",
      adj.P.Val < 0.05 & logFC < 0 ~ "Lower expression",
      TRUE ~ "Not FDR-significant"
    )
  )

summary_table <- de |>
  group_by(Contrast, contrast) |>
  summarise(
    tested = n(),
    higher = sum(status == "Higher expression"),
    lower = sum(status == "Lower expression"),
    total_fdr_05 = sum(adj.P.Val < 0.05),
    smallest_fdr = min(adj.P.Val, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(summary_table, file.path(table_dir, "enhanced_volcano_summary.csv"), row.names = FALSE)

labels <- de |>
  filter(adj.P.Val < 0.05, !is.na(gene_label), gene_label != "") |>
  group_by(contrast, status) |>
  arrange(adj.P.Val, desc(abs(logFC)), .by_group = TRUE) |>
  slice_head(n = 3) |>
  ungroup()
label_audit <- audit_hgnc_symbols(
  labels$gene_label,
  "Printed volcano labels",
  "Figure 3"
)
if (any(label_audit$status == "unrecognised")) {
  stop("Figure 3 selected an unrecognised HGNC label")
}
labels <- labels |>
  mutate(
    gene_label_archived = gene_label,
    gene_label = normalize_hgnc_symbols(gene_label, label_audit)
  )
data.table::fwrite(
  label_audit,
  file.path(table_dir, "figure3_hgnc_symbol_validation_biobouncer.csv")
)
write.csv(
  labels |>
    select(
      Contrast, contrast, gene_label, gene_label_archived,
      logFC, P.Value, adj.P.Val, status
    ),
  file.path(table_dir, "enhanced_volcano_label_manifest.csv"),
  row.names = FALSE
)

strip_labels <- summary_table |>
  mutate(strip = paste0(as.character(contrast), "\n", total_fdr_05, " genes at FDR < 0.05")) |>
  select(contrast, strip)
strip_lookup <- setNames(strip_labels$strip, strip_labels$contrast)

fdr_boundaries <- de |>
  filter(adj.P.Val < 0.05) |>
  group_by(contrast) |>
  summarise(
    boundary_y = -log10(max(nominal_p, na.rm = TRUE)),
    .groups = "drop"
  )

p <- ggplot(de, aes(display_logFC, display_y)) +
  geom_hline(
    data = fdr_boundaries,
    aes(yintercept = boundary_y),
    inherit.aes = FALSE,
    linewidth = 0.35,
    linetype = "dashed",
    colour = "grey35"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_point(
    data = filter(de, !capped),
    aes(colour = status),
    size = 0.65,
    alpha = 0.58,
    stroke = 0
  ) +
  geom_point(
    data = filter(de, capped),
    aes(colour = status),
    shape = 24,
    size = 1.15,
    fill = "white",
    stroke = 0.4,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = labels,
    aes(label = gene_label),
    size = 2.3,
    box.padding = 0.22,
    point.padding = 0.12,
    min.segment.length = 0,
    segment.size = 0.25,
    segment.color = "grey45",
    max.overlaps = Inf,
    seed = 42,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~contrast,
    nrow = 1,
    labeller = as_labeller(strip_lookup)
  ) +
  coord_cartesian(xlim = c(-x_limit, x_limit), ylim = c(0, y_limit), expand = FALSE) +
  scale_colour_manual(
    values = c(
      "Lower expression" = "#2166AC",
      "Not FDR-significant" = "#BDBDBD",
      "Higher expression" = "#B2182B"
    ),
    breaks = c("Higher expression", "Lower expression", "Not FDR-significant")
  ) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 25, by = 5)) +
  labs(
    x = expression(log[2]~fold~change),
    y = expression(-log[10](italic(P)*"-value")),
    colour = NULL
  ) +
  theme_classic(base_size = 9, base_family = "sans") +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.45),
    strip.text = element_text(size = 8.3, face = "bold", margin = margin(5, 3, 5, 3)),
    panel.spacing = grid::unit(7, "pt"),
    axis.title = element_text(size = 9.5),
    axis.text = element_text(size = 8, colour = "black"),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.35),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.key.width = grid::unit(11, "pt"),
    plot.margin = margin(6, 7, 4, 4)
  )

ggsave(paste0(figure_stem, ".pdf"), p, width = 7.25, height = 3.45, device = "pdf")
ggsave(paste0(figure_stem, ".png"), p, width = 7.25, height = 3.45, dpi = 600, bg = "white")
ggsave(paste0(figure_stem, ".svg"), p, width = 7.25, height = 3.45, bg = "white")

print(summary_table)

