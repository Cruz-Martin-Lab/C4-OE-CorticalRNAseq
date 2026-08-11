# 14_make_main_figure_panel_D.R
#
# Main Figure Panel D
#
# Shared Reactome pathways grouped into nonredundant biological modules.
#
# This version uses a compact paired-dot layout with:
#   - pathway labels on the y-axis
#   - human and mouse NES shown side by side
#   - connecting lines between datasets
#   - module headers embedded in the pathway labels
#   - legend placed below the plot

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# -------------------------------------------------------------------------
# 1. Define paths
# -------------------------------------------------------------------------

input_file <- file.path(
  "results",
  "tables",
  "shared_reactome_pathways_with_modules.csv"
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
# 2. Read and validate data
# -------------------------------------------------------------------------

reactome <- readr::read_csv(
  input_file,
  show_col_types = FALSE
)

required_columns <- c(
  "biological_module",
  "pathway_name",
  "human_NES",
  "human_fdr",
  "mouse_NES",
  "mouse_fdr",
  "nes_direction_relation"
)

missing_columns <- setdiff(
  required_columns,
  names(reactome)
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
# 3. Create concise pathway labels
# -------------------------------------------------------------------------

plot_data <- reactome %>%
  dplyr::mutate(
    pathway_label = dplyr::case_when(
      
      pathway_name ==
        "REACTOME_CARGO_TRAFFICKING_TO_THE_PERICILIARY_MEMBRANE" ~
        "Periciliary cargo trafficking",
      
      pathway_name ==
        "REACTOME_FORMATION_OF_TUBULIN_FOLDING_INTERMEDIATES_BY_CCT_TRIC" ~
        "Tubulin folding by CCT/TRiC",
      
      pathway_name ==
        "REACTOME_COOPERATION_OF_PREFOLDIN_AND_TRIC_CCT_IN_ACTIN_AND_TUBULIN_FOLDING" ~
        "Prefoldin-CCT/TRiC folding",
      
      pathway_name ==
        "REACTOME_TRANSPORT_TO_THE_GOLGI_AND_SUBSEQUENT_MODIFICATION" ~
        "Golgi transport and modification",
      
      pathway_name ==
        "REACTOME_TRANSLOCATION_OF_SLC2A4_GLUT4_TO_THE_PLASMA_MEMBRANE" ~
        "Plasma-membrane cargo translocation",
      
      pathway_name ==
        "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC" ~
        "ER folding and peptide loading",
      
      pathway_name ==
        "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION" ~
        "Ubiquitination and protein processing",
      
      pathway_name ==
        "REACTOME_ADAPTIVE_IMMUNE_SYSTEM" ~
        "Stress/immune-associated machinery",
      
      pathway_name ==
        "REACTOME_SIGNALING_BY_INTERLEUKINS" ~
        "Interleukin-associated signaling",
      
      pathway_name ==
        "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION" ~
        "Translation elongation",
      
      stringr::str_detect(
        pathway_name,
        "REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR"
      ) ~
        "IGF transport and extracellular remodeling",
      
      TRUE ~
        pathway_name %>%
        stringr::str_remove("^REACTOME_") %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_to_sentence()
    )
  )

# -------------------------------------------------------------------------
# 4. Define module order
# -------------------------------------------------------------------------

# Curated module order for the biologically-grouped modules. Any shared
# pathway assigned a module NOT in this list (e.g. "Other shared Reactome
# pathway") would otherwise be coerced to NA by factor() and silently dropped
# from the panel. A full run surfaced 9 such pathways
# (synaptic transmission, neuronal system, potassium channels, rRNA
# processing, SRP-dependent targeting, cytokine signalling, ...). Append every
# module actually present, so the panel shows all shared pathways rather than
# hiding those that fall outside the curated groups. The biological grouping /
# relabelling of the appended pathways is a presentation decision left open.
module_order <- c(
  "Cytoskeletal folding and intracellular trafficking",
  "ER folding, ubiquitination, and antigen processing",
  "Stress and immune-associated signaling",
  "Translation elongation",
  "IGF and extracellular secretory remodeling"
)

extra_modules <- setdiff(
  unique(as.character(plot_data$biological_module)),
  module_order
)
if (length(extra_modules) > 0) {
  message(
    "Panel D: ", length(extra_modules),
    " module(s) not in the curated order will be appended so no shared ",
    "pathway is dropped: ", paste(extra_modules, collapse = "; ")
  )
  module_order <- c(module_order, sort(extra_modules))
}

plot_data <- plot_data %>%
  dplyr::mutate(
    biological_module = factor(
      biological_module,
      levels = module_order
    )
  )

# -------------------------------------------------------------------------
# 5. Define pathway order within modules
# -------------------------------------------------------------------------

plot_data <- plot_data %>%
  dplyr::group_by(
    biological_module
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      (human_NES + mouse_NES) / 2
    ),
    .by_group = TRUE
  ) %>%
  dplyr::ungroup()

# -------------------------------------------------------------------------
# 6. Add module-header rows
# -------------------------------------------------------------------------

module_headers <- plot_data %>%
  dplyr::distinct(
    biological_module
  ) %>%
  dplyr::mutate(
    pathway_name = NA_character_,
    pathway_label = paste0(
      "MODULE: ",
      as.character(biological_module)
    ),
    human_NES = NA_real_,
    human_fdr = NA_real_,
    mouse_NES = NA_real_,
    mouse_fdr = NA_real_,
    nes_direction_relation = NA_character_,
    is_module_header = TRUE
  )

plot_data <- plot_data %>%
  dplyr::mutate(
    is_module_header = FALSE
  )

combined_rows <- dplyr::bind_rows(
  module_headers,
  plot_data
)

# -------------------------------------------------------------------------
# 7. Build explicit row order
# -------------------------------------------------------------------------

ordered_labels <- character(0)

for (module_name in module_order) {
  
  module_header_label <- paste0(
    "MODULE: ",
    module_name
  )
  
  module_pathways <- plot_data %>%
    dplyr::filter(
      as.character(biological_module) ==
        module_name
    ) %>%
    dplyr::pull(
      pathway_label
    )
  
  ordered_labels <- c(
    ordered_labels,
    module_header_label,
    module_pathways
  )
}

combined_rows <- combined_rows %>%
  dplyr::mutate(
    display_label = factor(
      pathway_label,
      levels = rev(ordered_labels)
    )
  )

# -------------------------------------------------------------------------
# 8. Create long-format point data
# -------------------------------------------------------------------------

point_data <- plot_data %>%
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
      human_NES = "Human schizophrenia",
      mouse_NES = "Mouse C4-OE"
    ),
    
    display_label = factor(
      pathway_label,
      levels = rev(ordered_labels)
    )
  )

# -------------------------------------------------------------------------
# 9. Create segment data
# -------------------------------------------------------------------------

segment_data <- plot_data %>%
  dplyr::mutate(
    display_label = factor(
      pathway_label,
      levels = rev(ordered_labels)
    )
  )

# -------------------------------------------------------------------------
# 10. Create module-divider positions
# -------------------------------------------------------------------------

label_index <- tibble::tibble(
  label = levels(combined_rows$display_label),
  y_position = seq_along(
    levels(combined_rows$display_label)
  )
)

module_header_positions <- label_index %>%
  dplyr::filter(
    stringr::str_starts(
      label,
      "MODULE:"
    )
  )

# -------------------------------------------------------------------------
# 11. Create custom y-axis labels
# -------------------------------------------------------------------------

label_lookup <- setNames(
  levels(combined_rows$display_label),
  levels(combined_rows$display_label)
)

label_function <- function(x) {
  
  labels <- label_lookup[x]
  
  labels <- ifelse(
    stringr::str_starts(
      labels,
      "MODULE:"
    ),
    paste0(
      "Module: ",
      stringr::str_remove(
        labels,
        "^MODULE: "
      )
    ),
    paste0(
      "   ",
      labels
    )
  )
  
  labels
}

# -------------------------------------------------------------------------
# 12. Create Panel D
# -------------------------------------------------------------------------

panel_d <- ggplot2::ggplot() +
  
  # Zero reference line
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.45,
    color = "grey65"
  ) +
  
  # Connecting lines
  ggplot2::geom_segment(
    data = segment_data,
    ggplot2::aes(
      x = human_NES,
      xend = mouse_NES,
      y = display_label,
      yend = display_label
    ),
    linewidth = 0.8,
    color = "grey55"
  ) +
  
  # Human and mouse points
  ggplot2::geom_point(
    data = point_data,
    ggplot2::aes(
      x = NES,
      y = display_label,
      shape = dataset
    ),
    size = 3.0,
    stroke = 0.9,
    position = ggplot2::position_dodge(
      width = 0.4
    )
  ) +
  
  # Horizontal separators above module headers
  ggplot2::geom_hline(
    data = module_header_positions,
    ggplot2::aes(
      yintercept = y_position - 0.5
    ),
    linewidth = 0.35,
    color = "grey85"
  ) +
  
  ggplot2::scale_shape_manual(
    values = c(
      "Human schizophrenia" = 16,
      "Mouse C4-OE" = 17
    )
  ) +
  
  ggplot2::scale_y_discrete(
    labels = label_function,
    drop = FALSE
  ) +
  
  ggplot2::labs(
    x = "Normalized enrichment score",
    y = NULL,
    shape = NULL
  ) +
  
  ggplot2::theme_classic(
    base_size = 10
  ) +
  
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(
      face = "bold",
      size = 10.5
    ),
    
    axis.text.x = ggplot2::element_text(
      color = "black",
      size = 9
    ),
    
    axis.text.y = ggplot2::element_text(
      color = "black",
      size = 9,
      hjust = 1
    ),
    
    axis.ticks.y = ggplot2::element_blank(),
    
    legend.position = "bottom",
    
    legend.direction = "horizontal",
    
    legend.text = ggplot2::element_text(
      size = 9
    ),
    
    legend.key.width = grid::unit(
      0.9,
      "lines"
    ),
    
    plot.margin = ggplot2::margin(
      t = 8,
      r = 10,
      b = 8,
      l = 8
    )
  )

# -------------------------------------------------------------------------
# 13. Bold module-header labels
# -------------------------------------------------------------------------

panel_d <- panel_d +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(
      color = "black",
      size = 9
    )
  )

# -------------------------------------------------------------------------
# 14. Save Panel D
# -------------------------------------------------------------------------

panel_d_pdf <- file.path(
  figure_dir,
  "panel_D_shared_Reactome_modules.pdf"
)

panel_d_png <- file.path(
  figure_dir,
  "panel_D_shared_Reactome_modules.png"
)

ggplot2::ggsave(
  filename = panel_d_pdf,
  plot = panel_d,
  width = 8.2,
  height = 7.6,
  units = "in"
)

ggplot2::ggsave(
  filename = panel_d_png,
  plot = panel_d,
  width = 8.2,
  height = 7.6,
  units = "in",
  dpi = 600
)

# -------------------------------------------------------------------------
# 15. Save figure-source table
# -------------------------------------------------------------------------

readr::write_csv(
  plot_data,
  file.path(
    table_dir,
    "main_figure_panel_D_shared_Reactome_modules.csv"
  )
)

# -------------------------------------------------------------------------
# 16. Completion messages
# -------------------------------------------------------------------------

message(
  "Revised Main Figure Panel D completed."
)

message(
  "PDF: ",
  panel_d_pdf
)

message(
  "PNG: ",
  panel_d_png
)