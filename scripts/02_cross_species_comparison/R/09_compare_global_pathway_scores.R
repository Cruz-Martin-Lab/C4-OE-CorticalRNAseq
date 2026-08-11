# 09_compare_global_pathway_scores.R
#
# Purpose:
# Compare normalized enrichment scores across all GO Biological Process
# pathways tested in both datasets.
#
# This analysis does not require pathways to pass an FDR threshold.
# It asks whether the overall pathway-ranking landscapes are related.

library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)

# -------------------------------------------------------------------------
# 1. Define paths
# -------------------------------------------------------------------------

comparison_file <- file.path(
  "results",
  "tables",
  "human_mouse_exact_go_bp_pathway_comparison.csv"
)

table_dir <- file.path("results", "tables")
figure_dir <- file.path("results", "figures")
log_dir <- file.path("results", "logs")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(comparison_file)) {
  stop("Pathway comparison file not found: ", comparison_file)
}

# -------------------------------------------------------------------------
# 2. Read and validate data
# -------------------------------------------------------------------------

pathways <- readr::read_csv(
  comparison_file,
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
    "Required columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

analysis_data <- pathways %>%
  dplyr::filter(
    !is.na(human_NES),
    !is.na(mouse_NES),
    is.finite(human_NES),
    is.finite(mouse_NES)
  ) %>%
  dplyr::mutate(
    pathway_label = pathway_name %>%
      stringr::str_remove("^GOBP_") %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_to_sentence(),
    
    same_nes_direction =
      sign(human_NES) == sign(mouse_NES),
    
    absolute_nes_sum =
      abs(human_NES) + abs(mouse_NES),
    
    mean_absolute_nes =
      (abs(human_NES) + abs(mouse_NES)) / 2
  )

# -------------------------------------------------------------------------
# 3. Global pathway correlation
# -------------------------------------------------------------------------

spearman_test <- cor.test(
  analysis_data$human_NES,
  analysis_data$mouse_NES,
  method = "spearman",
  exact = FALSE
)

pearson_test <- cor.test(
  analysis_data$human_NES,
  analysis_data$mouse_NES,
  method = "pearson"
)

global_summary <- tibble::tibble(
  pathways_compared = nrow(analysis_data),
  
  spearman_rho =
    unname(spearman_test$estimate),
  
  spearman_p_value =
    spearman_test$p.value,
  
  pearson_r =
    unname(pearson_test$estimate),
  
  pearson_p_value =
    pearson_test$p.value,
  
  same_direction_pathways =
    sum(
      analysis_data$same_nes_direction,
      na.rm = TRUE
    ),
  
  opposite_direction_pathways =
    sum(
      !analysis_data$same_nes_direction,
      na.rm = TRUE
    ),
  
  same_direction_proportion =
    mean(
      analysis_data$same_nes_direction,
      na.rm = TRUE
    )
)

readr::write_csv(
  global_summary,
  file.path(
    log_dir,
    "global_go_bp_nes_correlation_summary.csv"
  )
)

# -------------------------------------------------------------------------
# 4. Rank and percentile each pathway
# -------------------------------------------------------------------------

ranked_pathways <- analysis_data %>%
  dplyr::mutate(
    human_positive_rank =
      rank(
        -human_NES,
        ties.method = "average"
      ),
    
    mouse_positive_rank =
      rank(
        -mouse_NES,
        ties.method = "average"
      ),
    
    human_positive_percentile =
      1 -
      (
        human_positive_rank - 1
      ) /
      (
        n() - 1
      ),
    
    mouse_positive_percentile =
      1 -
      (
        mouse_positive_rank - 1
      ) /
      (
        n() - 1
      ),
    
    mean_positive_percentile =
      (
        human_positive_percentile +
          mouse_positive_percentile
      ) / 2,
    
    cross_species_rank_score =
      sqrt(
        human_positive_percentile *
          mouse_positive_percentile
      )
  ) %>%
  dplyr::arrange(
    dplyr::desc(cross_species_rank_score),
    dplyr::desc(mean_absolute_nes)
  )

readr::write_csv(
  ranked_pathways,
  file.path(
    table_dir,
    "all_go_bp_cross_species_nes_ranks.csv"
  )
)

# -------------------------------------------------------------------------
# 5. Extract RNA-localization family
# -------------------------------------------------------------------------

rna_family <- ranked_pathways %>%
  dplyr::filter(
    stringr::str_detect(
      pathway_name,
      "RNA_LOCALIZATION|MRNA_TRANSPORT|RNA_EXPORT|NUCLEAR_TRANSPORT"
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(cross_species_rank_score)
  )

readr::write_csv(
  rna_family,
  file.path(
    table_dir,
    "rna_localization_family_cross_species_ranks.csv"
  )
)

# -------------------------------------------------------------------------
# 6. Summarize prespecified biological modules
# -------------------------------------------------------------------------

module_pathways <- ranked_pathways %>%
  dplyr::mutate(
    biological_module = dplyr::case_when(
      stringr::str_detect(
        pathway_name,
        "RNA_LOCALIZATION|MRNA_TRANSPORT|RNA_EXPORT|NUCLEAR_TRANSPORT|NUCLEOBASE_CONTAINING_COMPOUND_TRANSPORT"
      ) ~
        "RNA localization and transport",
      
      stringr::str_detect(
        pathway_name,
        "CHOLESTEROL|STEROL|STEROID|ISOPRENOID|GLYCEROLIPID|PHOSPHOLIPID|LIPID"
      ) ~
        "Lipid and sterol metabolism",
      
      stringr::str_detect(
        pathway_name,
        "PROTEIN_FOLDING|PROTEIN_MATURATION|ENDOPLASMIC_RETICULUM_STRESS|PROTEASOM|AUTOPHAG"
      ) ~
        "Proteostasis and protein quality control",
      
      TRUE ~
        NA_character_
    )
  ) %>%
  dplyr::filter(
    !is.na(biological_module)
  )

calculate_module_summary <- function(module_data) {
  
  if (nrow(module_data) < 3) {
    return(
      tibble::tibble(
        pathways = nrow(module_data),
        spearman_rho = NA_real_,
        spearman_p_value = NA_real_,
        same_direction_proportion =
          mean(
            module_data$same_nes_direction,
            na.rm = TRUE
          ),
        median_human_NES =
          median(
            module_data$human_NES,
            na.rm = TRUE
          ),
        median_mouse_NES =
          median(
            module_data$mouse_NES,
            na.rm = TRUE
          )
      )
    )
  }
  
  correlation <- cor.test(
    module_data$human_NES,
    module_data$mouse_NES,
    method = "spearman",
    exact = FALSE
  )
  
  tibble::tibble(
    pathways = nrow(module_data),
    
    spearman_rho =
      unname(correlation$estimate),
    
    spearman_p_value =
      correlation$p.value,
    
    same_direction_proportion =
      mean(
        module_data$same_nes_direction,
        na.rm = TRUE
      ),
    
    median_human_NES =
      median(
        module_data$human_NES,
        na.rm = TRUE
      ),
    
    median_mouse_NES =
      median(
        module_data$mouse_NES,
        na.rm = TRUE
      )
  )
}

module_summary <- module_pathways %>%
  dplyr::group_by(
    biological_module
  ) %>%
  dplyr::group_modify(
    ~ calculate_module_summary(.x)
  ) %>%
  dplyr::ungroup()

readr::write_csv(
  module_summary,
  file.path(
    table_dir,
    "prespecified_module_nes_summary.csv"
  )
)

# -------------------------------------------------------------------------
# 7. Create pathway scatter plot
# -------------------------------------------------------------------------

plot_data <- ranked_pathways %>%
  dplyr::mutate(
    highlight = dplyr::case_when(
      pathway_name ==
        "GOBP_RNA_LOCALIZATION" ~
        "RNA localization",
      
      stringr::str_detect(
        pathway_name,
        "CHOLESTEROL_BIOSYNTHETIC_PROCESS_VIA_DESMOSTEROL|STEROL_BIOSYNTHETIC_PROCESS|PROTEIN_FOLDING"
      ) ~
        "Secondary pathway",
      
      TRUE ~
        "Other pathways"
    )
  )

# The RNA-localization point was being buried under the ~2,700 overlapping
# background points because a single geom_point layer draws points in data
# order. Split the data and draw the highlighted point in its OWN layer, LAST,
# as a red, black-outlined circle with a bold label -- so it is always rendered
# on top and clearly visible.
background_data <- plot_data %>%
  dplyr::filter(highlight != "RNA localization")

rna_localization_data <- plot_data %>%
  dplyr::filter(highlight == "RNA localization")

nes_scatter <- ggplot2::ggplot() +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.4
  ) +
  ggplot2::geom_point(
    data = background_data,
    ggplot2::aes(
      x = mouse_NES,
      y = human_NES,
      shape = highlight
    ),
    colour = "grey55",
    alpha = 0.45,
    size = 1.6
  ) +
  ggplot2::geom_point(
    data = rna_localization_data,
    ggplot2::aes(
      x = mouse_NES,
      y = human_NES,
      fill = "RNA localization"
    ),
    shape = 21,
    colour = "black",
    size = 4,
    stroke = 1.1
  ) +
  ggplot2::geom_text(
    data = rna_localization_data,
    ggplot2::aes(
      x = mouse_NES,
      y = human_NES,
      label = "RNA localization"
    ),
    nudge_x = 0.16,
    nudge_y = 0.18,
    fontface = "bold",
    colour = "red",
    size = 3.6
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Other pathways" = 16,
      "Secondary pathway" = 17
    )
  ) +
  ggplot2::scale_fill_manual(
    values = c("RNA localization" = "red")
  ) +
  ggplot2::labs(
    title =
      "Global GO Biological Process comparison",
    subtitle =
      paste0(
        "Spearman rho = ",
        round(
          unname(spearman_test$estimate),
          3
        ),
        "; P = ",
        signif(
          spearman_test$p.value,
          3
        )
      ),
    x =
      "Mouse C4-OE NES",
    y =
      "Human schizophrenia synaptic-proteome NES",
    shape =
      NULL,
    fill =
      NULL
  ) +
  ggplot2::theme_classic()

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "global_human_mouse_go_bp_nes_scatter.pdf"
  ),
  plot = nes_scatter,
  width = 7,
  height = 6
)

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "global_human_mouse_go_bp_nes_scatter.png"
  ),
  plot = nes_scatter,
  width = 7,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------------
# 8. Record RNA-localization rank specifically
# -------------------------------------------------------------------------

rna_rank_summary <- ranked_pathways %>%
  dplyr::filter(
    pathway_name ==
      "GOBP_RNA_LOCALIZATION"
  ) %>%
  dplyr::select(
    pathway_name,
    human_NES,
    human_fdr,
    mouse_NES,
    mouse_fdr,
    human_positive_rank,
    mouse_positive_rank,
    human_positive_percentile,
    mouse_positive_percentile,
    mean_positive_percentile,
    cross_species_rank_score,
    same_nes_direction
  )

readr::write_csv(
  rna_rank_summary,
  file.path(
    log_dir,
    "rna_localization_global_rank_summary.csv"
  )
)

# -------------------------------------------------------------------------
# 9. Record software environment
# -------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    log_dir,
    "sessionInfo_09_compare_global_pathway_scores.txt"
  )
)

message(
  "Global pathway-score comparison completed."
)

message(
  "Global summary: ",
  file.path(
    log_dir,
    "global_go_bp_nes_correlation_summary.csv"
  )
)

message(
  "RNA rank summary: ",
  file.path(
    log_dir,
    "rna_localization_global_rank_summary.csv"
  )
)

message(
  "Module summary: ",
  file.path(
    table_dir,
    "prespecified_module_nes_summary.csv"
  )
)