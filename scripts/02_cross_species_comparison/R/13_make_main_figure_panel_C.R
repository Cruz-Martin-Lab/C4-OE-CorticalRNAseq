# 13_make_main_figure_panel_C.R
#
# Main Figure Panel C
#
# RNA-localization convergence between:
#   Human schizophrenia synaptic proteome
#   Mouse neuronal C4-OE transcriptome

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(msigdbr)
library(fgsea)

# Package installation belongs in scripts/00_setup/, not in an analysis script:
# a side-effecting install.packages() here would change the environment midway
# through a run and, non-interactively, block waiting for a CRAN mirror.
library(patchwork)

# -------------------------------------------------------------------------
# 1. Define paths
# -------------------------------------------------------------------------

shared_universe_file <- file.path(
  "data_processed",
  "human_mouse_shared_measurable_universe.csv"
)

rna_gene_file <- file.path(
  "results",
  "tables",
  "rna_localization_shared_leading_edge_genes.csv"
)

rna_summary_file <- file.path(
  "results",
  "logs",
  "rna_localization_global_rank_summary.csv"
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

# -------------------------------------------------------------------------
# 2. Confirm that required files exist
# -------------------------------------------------------------------------

required_files <- c(
  shared_universe_file,
  rna_gene_file,
  rna_summary_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "The following required files are missing:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

# -------------------------------------------------------------------------
# 3. Read shared ortholog universe
# -------------------------------------------------------------------------

shared <- readr::read_csv(
  shared_universe_file,
  show_col_types = FALSE
)

required_shared_columns <- c(
  "human_gene_symbol",
  "scz_t_stat",
  "mouse_stat"
)

missing_shared_columns <- setdiff(
  required_shared_columns,
  names(shared)
)

if (length(missing_shared_columns) > 0) {
  stop(
    "Missing columns in the shared universe: ",
    paste(
      missing_shared_columns,
      collapse = ", "
    )
  )
}

# -------------------------------------------------------------------------
# 4. Build native-statistic ranked vectors
# -------------------------------------------------------------------------

human_rank_table <- shared %>%
  dplyr::filter(
    !is.na(human_gene_symbol),
    !is.na(scz_t_stat),
    is.finite(scz_t_stat)
  ) %>%
  dplyr::arrange(
    dplyr::desc(scz_t_stat),
    human_gene_symbol
  )

mouse_rank_table <- shared %>%
  dplyr::filter(
    !is.na(human_gene_symbol),
    !is.na(mouse_stat),
    is.finite(mouse_stat)
  ) %>%
  dplyr::arrange(
    dplyr::desc(mouse_stat),
    human_gene_symbol
  )

human_ranks <- human_rank_table$scz_t_stat
names(human_ranks) <- human_rank_table$human_gene_symbol

mouse_ranks <- mouse_rank_table$mouse_stat
names(mouse_ranks) <- mouse_rank_table$human_gene_symbol

human_ranks <- sort(
  human_ranks,
  decreasing = TRUE
)

mouse_ranks <- sort(
  mouse_ranks,
  decreasing = TRUE
)

shared_symbols <- intersect(
  names(human_ranks),
  names(mouse_ranks)
)

# -------------------------------------------------------------------------
# 5. Retrieve GO Biological Process RNA-localization genes
# -------------------------------------------------------------------------

go_bp_raw <- msigdbr::msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "C5",
  subcollection = "GO:BP"
)

required_msigdb_columns <- c(
  "gs_name",
  "gene_symbol"
)

missing_msigdb_columns <- setdiff(
  required_msigdb_columns,
  names(go_bp_raw)
)

if (length(missing_msigdb_columns) > 0) {
  stop(
    "Missing msigdbr columns: ",
    paste(
      missing_msigdb_columns,
      collapse = ", "
    )
  )
}

rna_pathway_genes <- go_bp_raw %>%
  dplyr::filter(
    gs_name == "GOBP_RNA_LOCALIZATION"
  ) %>%
  dplyr::pull(
    gene_symbol
  ) %>%
  unique()

rna_pathway_genes <- intersect(
  rna_pathway_genes,
  shared_symbols
)

if (length(rna_pathway_genes) < 10) {
  stop(
    "Too few measurable RNA-localization genes were found."
  )
}

# -------------------------------------------------------------------------
# 6. Read RNA-localization pathway statistics
# -------------------------------------------------------------------------

rna_summary <- readr::read_csv(
  rna_summary_file,
  show_col_types = FALSE
)

required_summary_columns <- c(
  "pathway_name",
  "human_NES",
  "human_fdr",
  "mouse_NES",
  "mouse_fdr",
  "human_positive_rank",
  "mouse_positive_rank"
)

missing_summary_columns <- setdiff(
  required_summary_columns,
  names(rna_summary)
)

if (length(missing_summary_columns) > 0) {
  stop(
    "Missing RNA summary columns: ",
    paste(
      missing_summary_columns,
      collapse = ", "
    )
  )
}

rna_summary <- rna_summary %>%
  dplyr::filter(
    pathway_name == "GOBP_RNA_LOCALIZATION"
  )

if (nrow(rna_summary) != 1) {
  stop(
    "Expected exactly one RNA-localization summary row."
  )
}

human_nes <- rna_summary$human_NES[[1]]
human_fdr <- rna_summary$human_fdr[[1]]
mouse_nes <- rna_summary$mouse_NES[[1]]
mouse_fdr <- rna_summary$mouse_fdr[[1]]

human_rank <- rna_summary$human_positive_rank[[1]]
mouse_rank <- rna_summary$mouse_positive_rank[[1]]

# -------------------------------------------------------------------------
# 7. Create enrichment plots
# -------------------------------------------------------------------------

human_enrichment_plot <- fgsea::plotEnrichment(
  pathway = rna_pathway_genes,
  stats = human_ranks
) +
  ggplot2::labs(
    title = "Human schizophrenia synaptic proteome",
    
    subtitle = paste0(
      "NES = ",
      sprintf("%.2f", human_nes),
      "; FDR = ",
      format(
        human_fdr,
        scientific = FALSE,
        digits = 2
      )
    ),
    
    x = "Genes ranked by schizophrenia t statistic",
    
    y = "Enrichment score"
  ) +
  ggplot2::theme_classic(
    base_size = 10
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 10.5
    ),
    
    plot.subtitle = ggplot2::element_text(
      size = 9.5
    ),
    
    axis.title = ggplot2::element_text(
      face = "bold",
      size = 9.5
    ),
    
    axis.text = ggplot2::element_text(
      color = "black",
      size = 8.5
    ),
    
    plot.margin = ggplot2::margin(
      t = 4,
      r = 8,
      b = 4,
      l = 4
    )
  )

mouse_enrichment_plot <- fgsea::plotEnrichment(
  pathway = rna_pathway_genes,
  stats = mouse_ranks
) +
  ggplot2::labs(
    title = "Mouse neuronal C4-OE transcriptome",
    
    subtitle = paste0(
      "NES = ",
      sprintf("%.2f", mouse_nes),
      "; FDR = ",
      format(
        mouse_fdr,
        scientific = FALSE,
        digits = 2
      )
    ),
    
    x = "Genes ranked by DESeq2 Wald statistic",
    
    y = "Enrichment score"
  ) +
  ggplot2::theme_classic(
    base_size = 10
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 10.5
    ),
    
    plot.subtitle = ggplot2::element_text(
      size = 9.5
    ),
    
    axis.title = ggplot2::element_text(
      face = "bold",
      size = 9.5
    ),
    
    axis.text = ggplot2::element_text(
      color = "black",
      size = 8.5
    ),
    
    plot.margin = ggplot2::margin(
      t = 4,
      r = 4,
      b = 4,
      l = 8
    )
  )

# -------------------------------------------------------------------------
# 8. Read shared leading-edge gene table
# -------------------------------------------------------------------------

rna_genes <- readr::read_csv(
  rna_gene_file,
  show_col_types = FALSE
)

required_gene_columns <- c(
  "human_gene_symbol",
  "human_t_stat",
  "human_fdr",
  "mouse_wald_stat",
  "mouse_padj",
  "direction_relation",
  "shared_leading_edge"
)

missing_gene_columns <- setdiff(
  required_gene_columns,
  names(rna_genes)
)

if (length(missing_gene_columns) > 0) {
  stop(
    "Missing RNA gene-table columns: ",
    paste(
      missing_gene_columns,
      collapse = ", "
    )
  )
}

rna_genes <- rna_genes %>%
  dplyr::filter(
    shared_leading_edge
  ) %>%
  dplyr::mutate(
    mean_rank_statistic = (
      rank(
        human_t_stat,
        ties.method = "average"
      ) +
        rank(
          mouse_wald_stat,
          ties.method = "average"
        )
    ) / 2
  ) %>%
  dplyr::arrange(
    mean_rank_statistic
  ) %>%
  dplyr::mutate(
    human_gene_symbol = factor(
      human_gene_symbol,
      levels = human_gene_symbol
    )
  )

shared_gene_count <- nrow(rna_genes)

same_direction_count <- sum(
  rna_genes$direction_relation ==
    "same direction",
  na.rm = TRUE
)

# -------------------------------------------------------------------------
# 9. Create separate human and mouse gene plots
# -------------------------------------------------------------------------

human_gene_plot <- ggplot2::ggplot(
  rna_genes,
  ggplot2::aes(
    x = human_t_stat,
    y = human_gene_symbol
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    color = "grey70"
  ) +
  ggplot2::geom_point(
    size = 2.3,
    shape = 16,
    color = "grey20"
  ) +
  ggplot2::labs(
    title = "Human",
    
    x = "Schizophrenia t statistic",
    
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 9.5
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 10,
      hjust = 0.5
    ),
    
    axis.title.x = ggplot2::element_text(
      face = "bold",
      size = 9.5
    ),
    
    axis.text.x = ggplot2::element_text(
      color = "black",
      size = 8.5
    ),
    
    axis.text.y = ggplot2::element_text(
      color = "black",
      size = 8.8
    ),
    
    plot.margin = ggplot2::margin(
      t = 4,
      r = 8,
      b = 4,
      l = 4
    )
  )

mouse_gene_plot <- ggplot2::ggplot(
  rna_genes,
  ggplot2::aes(
    x = mouse_wald_stat,
    y = human_gene_symbol
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    color = "grey70"
  ) +
  ggplot2::geom_point(
    size = 2.4,
    shape = 17,
    color = "grey20"
  ) +
  ggplot2::labs(
    title = "Mouse",
    
    x = "DESeq2 Wald statistic",
    
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 9.5
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 10,
      hjust = 0.5
    ),
    
    axis.title.x = ggplot2::element_text(
      face = "bold",
      size = 9.5
    ),
    
    axis.text.x = ggplot2::element_text(
      color = "black",
      size = 8.5
    ),
    
    axis.text.y = ggplot2::element_blank(),
    
    axis.ticks.y = ggplot2::element_blank(),
    
    plot.margin = ggplot2::margin(
      t = 4,
      r = 4,
      b = 4,
      l = 8
    )
  )

# -------------------------------------------------------------------------
# 10. Create shared-gene header
# -------------------------------------------------------------------------

gene_header <- patchwork::wrap_elements(
  full = grid::textGrob(
    paste0(
      "Shared RNA-localization leading-edge genes\n",
      shared_gene_count,
      " shared genes; ",
      same_direction_count,
      "/",
      shared_gene_count,
      " with concordant direction"
    ),
    x = 0,
    hjust = 0,
    gp = grid::gpar(
      fontsize = 11,
      fontface = "bold"
    )
  )
)

# -------------------------------------------------------------------------
# 11. Assemble the leading-edge section
# -------------------------------------------------------------------------

gene_plots <- human_gene_plot +
  mouse_gene_plot +
  patchwork::plot_layout(
    widths = c(
      1,
      1
    )
  )

bottom_section <- gene_header /
  gene_plots +
  patchwork::plot_layout(
    heights = c(
      0.12,
      1
    )
  )

# -------------------------------------------------------------------------
# 12. Assemble revised Panel C
# -------------------------------------------------------------------------

top_row <- human_enrichment_plot +
  mouse_enrichment_plot +
  patchwork::plot_layout(
    ncol = 2
  )

panel_c <- top_row /
  bottom_section +
  patchwork::plot_layout(
    heights = c(
      0.95,
      1.15
    )
  )

# -------------------------------------------------------------------------
# 13. Save individual components
# -------------------------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "panel_C1_human_RNA_localization_enrichment.pdf"
  ),
  plot = human_enrichment_plot,
  width = 5,
  height = 4,
  units = "in"
)

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "panel_C2_mouse_RNA_localization_enrichment.pdf"
  ),
  plot = mouse_enrichment_plot,
  width = 5,
  height = 4,
  units = "in"
)

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "panel_C3_human_shared_RNA_localization_genes.pdf"
  ),
  plot = human_gene_plot,
  width = 4,
  height = 6,
  units = "in"
)

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "panel_C4_mouse_shared_RNA_localization_genes.pdf"
  ),
  plot = mouse_gene_plot,
  width = 4,
  height = 6,
  units = "in"
)

# -------------------------------------------------------------------------
# 14. Save revised composite Panel C
# -------------------------------------------------------------------------

panel_c_pdf <- file.path(
  figure_dir,
  "panel_C_RNA_localization_convergence.pdf"
)

panel_c_png <- file.path(
  figure_dir,
  "panel_C_RNA_localization_convergence.png"
)

ggplot2::ggsave(
  filename = panel_c_pdf,
  plot = panel_c,
  width = 10.5,
  height = 9.2,
  units = "in"
)

ggplot2::ggsave(
  filename = panel_c_png,
  plot = panel_c,
  width = 10.5,
  height = 9.2,
  units = "in",
  dpi = 600
)

# -------------------------------------------------------------------------
# 15. Save figure-source table
# -------------------------------------------------------------------------

readr::write_csv(
  rna_genes,
  file.path(
    table_dir,
    "main_figure_panel_C_shared_RNA_localization_genes.csv"
  )
)

# -------------------------------------------------------------------------
# 16. Completion messages
# -------------------------------------------------------------------------

message(
  "Revised Main Figure Panel C completed."
)

message(
  "PDF: ",
  panel_c_pdf
)

message(
  "PNG: ",
  panel_c_png
)