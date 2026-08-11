# 15_make_main_figure_panel_A.R
#
# Main Figure Panel A
#
# Cross-species analysis workflow

library(ggplot2)

# -------------------------------------------------------------------------
# 0. Read the shared measurable universe size from the data
# -------------------------------------------------------------------------
# Read the shared measurable universe size from the actual universe table so
# the schematic always matches the run, rather than being hardcoded.

shared_universe_file <- file.path(
  "data_processed",
  "human_mouse_shared_measurable_universe.csv"
)

if (!file.exists(shared_universe_file)) {
  stop("Shared measurable universe table not found: ", shared_universe_file)
}

shared_universe_n <- nrow(
  read.csv(shared_universe_file, check.names = FALSE)
)

# -------------------------------------------------------------------------
# 1. Define output paths
# -------------------------------------------------------------------------

figure_dir <- file.path(
  "results",
  "figures",
  "main_figure"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

panel_a_pdf <- file.path(
  figure_dir,
  "panel_A_cross_species_workflow.pdf"
)

panel_a_png <- file.path(
  figure_dir,
  "panel_A_cross_species_workflow.png"
)

# -------------------------------------------------------------------------
# 2. Helper function for workflow boxes
# -------------------------------------------------------------------------

add_box <- function(
    plot,
    xmin,
    xmax,
    ymin,
    ymax,
    label,
    font_size = 4.0,
    font_face = "bold",
    line_width = 0.55
) {
  
  plot +
    ggplot2::annotate(
      "rect",
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = "white",
      color = "black",
      linewidth = line_width
    ) +
    ggplot2::annotate(
      "text",
      x = (xmin + xmax) / 2,
      y = (ymin + ymax) / 2,
      label = label,
      size = font_size,
      fontface = font_face,
      lineheight = 1.05
    )
}

# -------------------------------------------------------------------------
# 3. Create plotting canvas
# -------------------------------------------------------------------------

panel_a <- ggplot2::ggplot() +
  ggplot2::coord_cartesian(
    xlim = c(0, 10),
    ylim = c(0, 10),
    clip = "off"
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )

# -------------------------------------------------------------------------
# 4. Add workflow boxes
# -------------------------------------------------------------------------

panel_a <- add_box(
  panel_a,
  xmin = 0.4,
  xmax = 4.6,
  ymin = 8.0,
  ymax = 9.5,
  label = paste0(
    "Human schizophrenia\n",
    "synaptic proteome\n",
    "Ranked by proteomic t statistic"
  ),
  font_size = 4.0
)

panel_a <- add_box(
  panel_a,
  xmin = 5.4,
  xmax = 9.6,
  ymin = 8.0,
  ymax = 9.5,
  label = paste0(
    "Mouse neuronal C4-OE\n",
    "transcriptome\n",
    "Ranked by DESeq2 Wald statistic"
  ),
  font_size = 4.0
)

panel_a <- add_box(
  panel_a,
  xmin = 2.2,
  xmax = 7.8,
  ymin = 5.7,
  ymax = 7.2,
  label = paste0(
    "One-to-one human-mouse ortholog mapping\n",
    paste0("Shared measurable universe: ",
           format(shared_universe_n, big.mark = ","), " genes")
  ),
  font_size = 3.7
)

panel_a <- add_box(
  panel_a,
  xmin = 2.0,
  xmax = 8.0,
  ymin = 3.4,
  ymax = 5.0,
  label = paste0(
    "Preranked gene-set enrichment analysis\n",
    "GO Biological Process and Reactome pathways"
  ),
  font_size = 3.7
)

panel_a <- add_box(
  panel_a,
  xmin = 0.3,
  xmax = 3.4,
  ymin = 0.7,
  ymax = 2.5,
  label = paste0(
    "Global pathway landscape\n",
    "GO NES correspondence\n",
    "Panel B"
  ),
  font_size = 3.9
)

panel_a <- add_box(
  panel_a,
  xmin = 3.45,
  xmax = 6.55,
  ymin = 0.7,
  ymax = 2.5,
  label = paste0(
    "RNA-localization convergence\n",
    "Shared leading-edge genes\n",
    "Panel C"
  ),
  font_size = 3.9
)

panel_a <- add_box(
  panel_a,
  xmin = 6.6,
  xmax = 9.7,
  ymin = 0.7,
  ymax = 2.5,
  label = paste0(
    "Shared Reactome pathways\n",
    "Mechanistic modules\n",
    "Panel D"
  ),
  font_size = 3.9
)

# -------------------------------------------------------------------------
# 5. Define arrow style
# -------------------------------------------------------------------------

arrow_style <- grid::arrow(
  length = grid::unit(
    0.14,
    "inches"
  ),
  type = "closed"
)

# -------------------------------------------------------------------------
# 6. Add arrows from input datasets to shared ortholog universe
# -------------------------------------------------------------------------

panel_a <- panel_a +
  ggplot2::annotate(
    "segment",
    x = 2.5,
    xend = 4.3,
    y = 8.0,
    yend = 7.35,
    linewidth = 0.65,
    arrow = arrow_style
  )

panel_a <- panel_a +
  ggplot2::annotate(
    "segment",
    x = 7.5,
    xend = 5.7,
    y = 8.0,
    yend = 7.35,
    linewidth = 0.65,
    arrow = arrow_style
  )

# -------------------------------------------------------------------------
# 7. Add arrow from ortholog universe to pathway analysis
# -------------------------------------------------------------------------

panel_a <- panel_a +
  ggplot2::annotate(
    "segment",
    x = 5.0,
    xend = 5.0,
    y = 5.7,
    yend = 5.15,
    linewidth = 0.65,
    arrow = arrow_style
  )

# -------------------------------------------------------------------------
# 8. Add branched connector from pathway analysis to Panels B-D
# -------------------------------------------------------------------------

# Vertical stem
panel_a <- panel_a +
  ggplot2::annotate(
    "segment",
    x = 5.0,
    xend = 5.0,
    y = 3.4,
    yend = 3.0,
    linewidth = 0.6
  )

# Horizontal branch
panel_a <- panel_a +
  ggplot2::annotate(
    "segment",
    x = 1.85,
    xend = 8.15,
    y = 3.0,
    yend = 3.0,
    linewidth = 0.6
  )

# Arrow to Panel B
panel_a <- panel_a +
  ggplot2::annotate(
    "segment",
    x = 1.85,
    xend = 1.85,
    y = 3.0,
    yend = 2.65,
    linewidth = 0.6,
    arrow = arrow_style
  )

# Arrow to Panel C
panel_a <- panel_a +
  ggplot2::annotate(
    "segment",
    x = 5.0,
    xend = 5.0,
    y = 3.0,
    yend = 2.65,
    linewidth = 0.6,
    arrow = arrow_style
  )

# Arrow to Panel D
panel_a <- panel_a +
  ggplot2::annotate(
    "segment",
    x = 8.15,
    xend = 8.15,
    y = 3.0,
    yend = 2.65,
    linewidth = 0.6,
    arrow = arrow_style
  )

# -------------------------------------------------------------------------
# 9. Save Panel A
# -------------------------------------------------------------------------

ggplot2::ggsave(
  filename = panel_a_pdf,
  plot = panel_a,
  width = 8.4,
  height = 7.2,
  units = "in"
)

ggplot2::ggsave(
  filename = panel_a_png,
  plot = panel_a,
  width = 8.4,
  height = 7.2,
  units = "in",
  dpi = 600
)

# -------------------------------------------------------------------------
# 10. Completion messages
# -------------------------------------------------------------------------

message(
  "Revised Main Figure Panel A completed."
)

message(
  "PDF: ",
  panel_a_pdf
)

message(
  "PNG: ",
  panel_a_png
)