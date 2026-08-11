# 11_summarize_shared_reactome_modules.R
#
# Purpose:
# Collapse the Reactome pathways significant in both datasets into
# nonredundant biological modules and identify recurrent shared
# leading-edge genes.
#
# This is intended as the final pathway-level synthesis.

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)

# -------------------------------------------------------------------------
# 1. Paths
# -------------------------------------------------------------------------

input_file <- file.path(
  "results",
  "tables",
  "human_mouse_reactome_pathway_comparison.csv"
)

table_dir <- file.path("results", "tables")
figure_dir <- file.path("results", "figures")
log_dir <- file.path("results", "logs")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Reactome comparison file not found: ", input_file)
}

# -------------------------------------------------------------------------
# 2. Read and validate
# -------------------------------------------------------------------------

reactome <- readr::read_csv(
  input_file,
  show_col_types = FALSE
)

required_columns <- c(
  "pathway_name",
  "human_NES",
  "human_fdr",
  "human_leading_edge_genes",
  "mouse_NES",
  "mouse_fdr",
  "mouse_leading_edge_genes",
  "significant_both_fdr_0_05",
  "nes_direction_relation"
)

missing_columns <- setdiff(
  required_columns,
  names(reactome)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

shared_sig <- reactome %>%
  dplyr::filter(
    significant_both_fdr_0_05
  )

if (nrow(shared_sig) == 0) {
  stop("No Reactome pathways are significant in both datasets.")
}

# -------------------------------------------------------------------------
# 3. Assign nonredundant biological modules
# -------------------------------------------------------------------------

shared_sig <- shared_sig %>%
  dplyr::mutate(
    biological_module = dplyr::case_when(
      
      stringr::str_detect(
        pathway_name,
        "CARGO_TRAFFICKING_TO_THE_PERICILIARY_MEMBRANE|FORMATION_OF_TUBULIN_FOLDING_INTERMEDIATES|COOPERATION_OF_PREFOLDIN_AND_TRIC_CCT|TRANSPORT_TO_THE_GOLGI|TRANSLOCATION_OF_SLC2A4"
      ) ~
        "Cytoskeletal folding and intracellular trafficking",
      
      stringr::str_detect(
        pathway_name,
        "ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY|CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING"
      ) ~
        "ER folding, ubiquitination, and antigen processing",
      
      pathway_name ==
        "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION" ~
        "Translation elongation",
      
      pathway_name ==
        "REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS" ~
        "IGF and extracellular secretory remodeling",
      
      stringr::str_detect(
        pathway_name,
        "ADAPTIVE_IMMUNE_SYSTEM|SIGNALING_BY_INTERLEUKINS"
      ) ~
        "Stress and immune-associated signaling",
      
      TRUE ~
        "Other shared Reactome pathway"
    ),
    
    pathway_label = pathway_name %>%
      stringr::str_remove("^REACTOME_") %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_to_sentence()
  )

readr::write_csv(
  shared_sig,
  file.path(
    table_dir,
    "shared_reactome_pathways_with_modules.csv"
  )
)

# -------------------------------------------------------------------------
# 4. Parse leading-edge genes
# -------------------------------------------------------------------------

split_gene_string <- function(x) {
  
  if (
    is.na(x) ||
    x == ""
  ) {
    return(character(0))
  }
  
  unique(
    stringr::str_split(
      x,
      pattern = ";",
      simplify = FALSE
    )[[1]]
  )
}

gene_overlap_by_pathway <- shared_sig %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    human_gene_list = list(
      split_gene_string(
        human_leading_edge_genes
      )
    ),
    
    mouse_gene_list = list(
      split_gene_string(
        mouse_leading_edge_genes
      )
    ),
    
    shared_gene_list = list(
      intersect(
        human_gene_list,
        mouse_gene_list
      )
    ),
    
    human_leading_edge_count =
      length(human_gene_list),
    
    mouse_leading_edge_count =
      length(mouse_gene_list),
    
    shared_leading_edge_count =
      length(shared_gene_list),
    
    leading_edge_union_count =
      length(
        union(
          human_gene_list,
          mouse_gene_list
        )
      ),
    
    leading_edge_jaccard =
      ifelse(
        leading_edge_union_count > 0,
        shared_leading_edge_count /
          leading_edge_union_count,
        NA_real_
      ),
    
    shared_leading_edge_genes =
      paste(
        shared_gene_list,
        collapse = ";"
      )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    biological_module,
    pathway_name,
    pathway_label,
    human_NES,
    human_fdr,
    mouse_NES,
    mouse_fdr,
    nes_direction_relation,
    human_leading_edge_count,
    mouse_leading_edge_count,
    shared_leading_edge_count,
    leading_edge_union_count,
    leading_edge_jaccard,
    shared_leading_edge_genes
  ) %>%
  dplyr::arrange(
    biological_module,
    dplyr::desc(shared_leading_edge_count),
    human_fdr,
    mouse_fdr
  )

readr::write_csv(
  gene_overlap_by_pathway,
  file.path(
    table_dir,
    "shared_reactome_leading_edge_overlap.csv"
  )
)

# -------------------------------------------------------------------------
# 5. Build recurrent-gene table
# -------------------------------------------------------------------------

recurrent_genes <- shared_sig %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    shared_gene_list = list(
      intersect(
        split_gene_string(
          human_leading_edge_genes
        ),
        split_gene_string(
          mouse_leading_edge_genes
        )
      )
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    biological_module,
    pathway_name,
    shared_gene_list
  ) %>%
  tidyr::unnest_longer(
    shared_gene_list,
    values_to = "gene_symbol"
  ) %>%
  dplyr::filter(
    !is.na(gene_symbol),
    gene_symbol != ""
  ) %>%
  dplyr::group_by(
    gene_symbol
  ) %>%
  dplyr::summarise(
    shared_pathway_count =
      dplyr::n_distinct(pathway_name),
    
    module_count =
      dplyr::n_distinct(biological_module),
    
    modules =
      paste(
        sort(
          unique(biological_module)
        ),
        collapse = "; "
      ),
    
    pathways =
      paste(
        sort(
          unique(pathway_name)
        ),
        collapse = ";"
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(shared_pathway_count),
    dplyr::desc(module_count),
    gene_symbol
  )

readr::write_csv(
  recurrent_genes,
  file.path(
    table_dir,
    "recurrent_shared_reactome_leading_edge_genes.csv"
  )
)

# -------------------------------------------------------------------------
# 6. Summarize modules
# -------------------------------------------------------------------------

module_summary <- gene_overlap_by_pathway %>%
  dplyr::group_by(
    biological_module
  ) %>%
  dplyr::summarise(
    pathways =
      dplyr::n(),
    
    same_direction_pathways =
      sum(
        nes_direction_relation ==
          "same NES direction"
      ),
    
    same_direction_proportion =
      mean(
        nes_direction_relation ==
          "same NES direction"
      ),
    
    median_human_NES =
      median(
        human_NES,
        na.rm = TRUE
      ),
    
    median_mouse_NES =
      median(
        mouse_NES,
        na.rm = TRUE
      ),
    
    total_shared_leading_edge_occurrences =
      sum(
        shared_leading_edge_count,
        na.rm = TRUE
      ),
    
    median_shared_leading_edge_count =
      median(
        shared_leading_edge_count,
        na.rm = TRUE
      ),
    
    median_leading_edge_jaccard =
      median(
        leading_edge_jaccard,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(pathways),
    dplyr::desc(
      total_shared_leading_edge_occurrences
    )
  )

readr::write_csv(
  module_summary,
  file.path(
    table_dir,
    "shared_reactome_module_summary.csv"
  )
)

# -------------------------------------------------------------------------
# 7. Create recurrent-gene plot
# -------------------------------------------------------------------------

top_genes <- recurrent_genes %>%
  dplyr::filter(
    shared_pathway_count >= 2
  ) %>%
  dplyr::slice_head(
    n = 30
  )

if (nrow(top_genes) > 0) {
  
  recurrent_gene_plot <- ggplot2::ggplot(
    top_genes,
    ggplot2::aes(
      x = reorder(
        gene_symbol,
        shared_pathway_count
      ),
      y = shared_pathway_count
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title =
        "Recurrent cross-species Reactome leading-edge genes",
      
      subtitle =
        "Genes appearing in the shared leading edge of at least two pathways",
      
      x =
        NULL,
      
      y =
        "Number of shared significant Reactome pathways"
    ) +
    ggplot2::theme_classic()
  
  ggplot2::ggsave(
    filename = file.path(
      figure_dir,
      "recurrent_shared_reactome_leading_edge_genes.pdf"
    ),
    plot = recurrent_gene_plot,
    width = 7,
    height = 7
  )
  
  ggplot2::ggsave(
    filename = file.path(
      figure_dir,
      "recurrent_shared_reactome_leading_edge_genes.png"
    ),
    plot = recurrent_gene_plot,
    width = 7,
    height = 7,
    dpi = 300
  )
}

# -------------------------------------------------------------------------
# 8. Create pathway NES plot
# -------------------------------------------------------------------------

pathway_plot_data <- gene_overlap_by_pathway %>%
  dplyr::select(
    biological_module,
    pathway_label,
    human_NES,
    mouse_NES
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      human_NES,
      mouse_NES
    ),
    names_to = "dataset",
    values_to = "NES"
  ) %>%
  dplyr::mutate(
    dataset = dplyr::recode(
      dataset,
      human_NES =
        "Human SCZ synaptic proteome",
      mouse_NES =
        "Mouse C4-OE transcriptome"
    )
  )

pathway_nes_plot <- ggplot2::ggplot(
  pathway_plot_data,
  ggplot2::aes(
    x = NES,
    y = reorder(
      pathway_label,
      NES
    ),
    shape = dataset
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.4
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_dodge(
      width = 0.5
    ),
    size = 2.5
  ) +
  ggplot2::facet_grid(
    biological_module ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::labs(
    title =
      "Shared significant Reactome pathways",
    
    x =
      "Normalized enrichment score",
    
    y =
      NULL,
    
    shape =
      NULL
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    strip.text.y =
      ggplot2::element_text(
        angle = 0
      )
  )

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "shared_reactome_pathway_nes_by_module.pdf"
  ),
  plot = pathway_nes_plot,
  width = 9,
  height = 10
)

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "shared_reactome_pathway_nes_by_module.png"
  ),
  plot = pathway_nes_plot,
  width = 9,
  height = 10,
  dpi = 300
)

# -------------------------------------------------------------------------
# 9. Final analysis summary
# -------------------------------------------------------------------------

final_summary <- tibble::tibble(
  metric = c(
    "Shared Reactome pathways at FDR < 0.05",
    "Nonredundant modules represented",
    "Shared pathways with same NES direction",
    "Unique genes in at least one shared leading edge",
    "Genes recurring in at least two shared pathways",
    "Maximum shared pathway count for one gene"
  ),
  
  value = c(
    nrow(shared_sig),
    
    dplyr::n_distinct(
      shared_sig$biological_module
    ),
    
    sum(
      shared_sig$nes_direction_relation ==
        "same NES direction"
    ),
    
    nrow(recurrent_genes),
    
    sum(
      recurrent_genes$
        shared_pathway_count >= 2
    ),
    
    ifelse(
      nrow(recurrent_genes) > 0,
      max(
        recurrent_genes$
          shared_pathway_count
      ),
      0
    )
  )
)

readr::write_csv(
  final_summary,
  file.path(
    log_dir,
    "shared_reactome_module_analysis_summary.csv"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    log_dir,
    "sessionInfo_11_summarize_shared_reactome_modules.txt"
  )
)

message(
  "Shared Reactome module analysis completed."
)

message(
  "Final summary: ",
  file.path(
    log_dir,
    "shared_reactome_module_analysis_summary.csv"
  )
)

message(
  "Module summary: ",
  file.path(
    table_dir,
    "shared_reactome_module_summary.csv"
  )
)

message(
  "Recurrent genes: ",
  file.path(
    table_dir,
    "recurrent_shared_reactome_leading_edge_genes.csv"
  )
)