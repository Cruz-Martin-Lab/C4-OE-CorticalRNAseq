# 12_make_main_figure_panel_B.R
#
# Main Figure Panel B
# Global GO Biological Process comparison between:
#   Human schizophrenia synaptic proteome
#   Mouse neuronal C4-OE transcriptome

library(readr)
library(dplyr)
library(ggplot2)

# -------------------------------------------------------------------------
# 1. Define paths
# -------------------------------------------------------------------------

input_file <- file.path(
  "results",
  "tables",
  "human_mouse_exact_go_bp_pathway_comparison.csv"
)

figure_dir <- file.path(
  "results",
  "figures",
  "main_figure"
)

table_dir <- file.path(
  "results",
  "tables"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ",
    input_file
  )
}

# -------------------------------------------------------------------------
# 2. Read and validate pathway data
# -------------------------------------------------------------------------

pathways <- readr::read_csv(
  input_file,
  show_col_types = FALSE
)

required_columns <- c(
  "pathway_name",
  "human_NES",
  "human_fdr",
  "mouse_NES",
  "mouse_fdr"
)

missing_columns <- setdiff(
  required_columns,
  names(pathways)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

# -------------------------------------------------------------------------
# 3. Prepare plot data
# -------------------------------------------------------------------------

plot_data <- pathways %>%
  dplyr::filter(
    !is.na(human_NES),
    !is.na(mouse_NES),
    is.finite(human_NES),
    is.finite(mouse_NES)
  ) %>%
  dplyr::mutate(
    significant_both =
      !is.na(human_fdr) &
      !is.na(mouse_fdr) &
      human_fdr < 0.05 &
      mouse_fdr < 0.05,
    
    is_rna_localization =
      pathway_name ==
      "GOBP_RNA_LOCALIZATION",
    
    point_group =
      dplyr::case_when(
        is_rna_localization ~
          "RNA localization",
        
        significant_both ~
          "Significant in both",
        
        TRUE ~
          "All other pathways"
      )
  )

# -------------------------------------------------------------------------
# 4. Calculate global pathway correlation
# -------------------------------------------------------------------------

spearman_result <- stats::cor.test(
  plot_data$human_NES,
  plot_data$mouse_NES,
  method = "spearman",
  exact = FALSE
)

rho <- unname(
  spearman_result$estimate
)

p_value <- spearman_result$p.value

# -------------------------------------------------------------------------
# 5. Extract RNA-localization pathway
# -------------------------------------------------------------------------

rna_point <- plot_data %>%
  dplyr::filter(
    is_rna_localization
  )

if (nrow(rna_point) != 1) {
  stop(
    paste0(
      "Expected exactly one GOBP_RNA_LOCALIZATION row, but found ",
      nrow(rna_point),
      "."
    )
  )
}

rna_mouse_nes <- rna_point$mouse_NES[[1]]
rna_human_nes <- rna_point$human_NES[[1]]

# -------------------------------------------------------------------------
# 6. Create Main Figure Panel B
# -------------------------------------------------------------------------

panel_b <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = mouse_NES,
    y = human_NES
  )
) +
  
  # Zero-reference lines
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    color = "grey70"
  ) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    color = "grey70"
  ) +
  
  # All background pathways
  ggplot2::geom_point(
    data = plot_data %>%
      dplyr::filter(
        point_group ==
          "All other pathways"
      ),
    size = 1.15,
    alpha = 0.30,
    color = "grey55"
  ) +
  
  # Other pathways significant in both datasets
  ggplot2::geom_point(
    data = plot_data %>%
      dplyr::filter(
        point_group ==
          "Significant in both"
      ),
    size = 2.1,
    alpha = 0.85,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.7
  ) +
  
  # RNA-localization pathway
  ggplot2::geom_point(
    data = rna_point,
    size = 4,
    shape = 21,
    fill = "#C43C39",
    color = "black",
    stroke = 0.8
  ) +
  
  # Connector line placed to the upper-left
  ggplot2::annotate(
    geom = "segment",
    x = rna_mouse_nes,
    y = rna_human_nes,
    xend = rna_mouse_nes - 0.28,
    yend = rna_human_nes + 0.24,
    linewidth = 0.5,
    color = "black"
  ) +
  
  # RNA-localization label
  ggplot2::annotate(
    geom = "text",
    x = rna_mouse_nes - 0.32,
    y = rna_human_nes + 0.28,
    label = "RNA localization",
    hjust = 1,
    vjust = 0,
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +
  
  # Correlation annotation
  ggplot2::annotate(
    geom = "text",
    x = -Inf,
    y = Inf,
    label = paste0(
      "Spearman rho = ",
      sprintf(
        "%.3f",
        rho
      ),
      "\nP = ",
      format(
        p_value,
        scientific = TRUE,
        digits = 2
      )
    ),
    hjust = -0.08,
    vjust = 1.15,
    size = 3.3,
    color = "black"
  ) +
  
  # Axis labels
  ggplot2::labs(
    x =
      "Mouse C4-OE normalized enrichment score",
    
    y =
      paste0(
        "Human schizophrenia synaptic-proteome\n",
        "normalized enrichment score"
      )
  ) +
  
  ggplot2::coord_cartesian(
    clip = "off"
  ) +
  
  ggplot2::theme_classic(
    base_size = 10.5
  ) +
  
  ggplot2::theme(
    axis.title = ggplot2::element_text(
      face = "bold",
      size = 10.5
    ),
    
    axis.text = ggplot2::element_text(
      color = "black",
      size = 9.5
    ),
    
    axis.line = ggplot2::element_line(
      color = "black",
      linewidth = 0.55
    ),
    
    axis.ticks = ggplot2::element_line(
      color = "black",
      linewidth = 0.45
    ),
    
    plot.margin = ggplot2::margin(
      t = 12,
      r = 12,
      b = 8,
      l = 8
    )
  )

# -------------------------------------------------------------------------
# 7. Save the revised panel
# -------------------------------------------------------------------------

pdf_file <- file.path(
  figure_dir,
  "panel_B_global_GO_pathway_landscape.pdf"
)

png_file <- file.path(
  figure_dir,
  "panel_B_global_GO_pathway_landscape.png"
)

ggplot2::ggsave(
  filename = pdf_file,
  plot = panel_b,
  width = 5.4,
  height = 4.8,
  units = "in"
)

ggplot2::ggsave(
  filename = png_file,
  plot = panel_b,
  width = 5.4,
  height = 4.8,
  units = "in",
  dpi = 600
)

# -------------------------------------------------------------------------
# 8. Save the values highlighted in Panel B
# -------------------------------------------------------------------------

readr::write_csv(
  rna_point %>%
    dplyr::select(
      pathway_name,
      human_NES,
      human_fdr,
      mouse_NES,
      mouse_fdr
    ),
  file.path(
    table_dir,
    "main_figure_panel_B_RNA_localization_values.csv"
  )
)

# -------------------------------------------------------------------------
# 9. Completion messages
# -------------------------------------------------------------------------

message(
  "Revised Main Figure Panel B completed."
)

message(
  "PDF: ",
  pdf_file
)

message(
  "PNG: ",
  png_file
)