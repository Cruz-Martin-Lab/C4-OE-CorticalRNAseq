# 08_build_cross_species_summary.R
#
# Purpose:
# Build presentation-ready summaries of cross-species convergence between:
#   1. Human schizophrenia synaptic proteomics
#   2. Mouse neuronal C4-OE transcriptomics
#
# Primary result:
#   RNA localization is significantly enriched in both datasets and is
#   supported by a shared leading-edge gene core.
#
# Secondary modules:
#   - lipid and sterol metabolism
#   - proteostasis and protein quality control
#
# Direction is reported descriptively but is not used as a requirement
# for defining biological convergence.

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(msigdbr)
library(fgsea)
library(ggplot2)

# -------------------------------------------------------------------------
# 1. Define paths
# -------------------------------------------------------------------------

shared_file <- file.path(
  "data_processed",
  "human_mouse_shared_measurable_universe.csv"
)

human_gsea_file <- file.path(
  "results",
  "tables",
  "human_scz_go_bp_preranked_gsea.csv"
)

mouse_gsea_file <- file.path(
  "results",
  "tables",
  "mouse_c4oe_go_bp_preranked_gsea.csv"
)

semantic_file <- file.path(
  "results",
  "tables",
  "best_related_human_go_pathway_per_mouse_pathway.csv"
)

table_dir <- file.path("results", "tables")
figure_dir <- file.path("results", "figures")
log_dir <- file.path("results", "logs")

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

required_files <- c(
  shared_file,
  human_gsea_file,
  mouse_gsea_file,
  semantic_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Required files are missing: ",
    paste(missing_files, collapse = ", ")
  )
}

# -------------------------------------------------------------------------
# 2. Read inputs
# -------------------------------------------------------------------------

shared <- readr::read_csv(
  shared_file,
  show_col_types = FALSE
)

human_gsea <- readr::read_csv(
  human_gsea_file,
  show_col_types = FALSE
)

mouse_gsea <- readr::read_csv(
  mouse_gsea_file,
  show_col_types = FALSE
)

semantic_matches <- readr::read_csv(
  semantic_file,
  show_col_types = FALSE
)

required_shared_columns <- c(
  "human_gene_symbol",
  "mouse_gene_symbol",
  "scz_log2fc",
  "scz_t_stat",
  "scz_p_value",
  "scz_fdr",
  "mouse_log2fc",
  "mouse_stat",
  "mouse_p_value",
  "mouse_padj"
)

missing_shared_columns <- setdiff(
  required_shared_columns,
  names(shared)
)

if (length(missing_shared_columns) > 0) {
  stop(
    "Missing shared-universe columns: ",
    paste(missing_shared_columns, collapse = ", ")
  )
}

# -------------------------------------------------------------------------
# 3. Extract RNA-localization pathway results
# -------------------------------------------------------------------------

rna_pathway_name <- "GOBP_RNA_LOCALIZATION"

human_rna <- human_gsea %>%
  dplyr::filter(
    pathway_name == rna_pathway_name
  )

mouse_rna <- mouse_gsea %>%
  dplyr::filter(
    pathway_name == rna_pathway_name
  )

if (nrow(human_rna) != 1) {
  stop(
    "Expected exactly one human RNA-localization pathway row."
  )
}

if (nrow(mouse_rna) != 1) {
  stop(
    "Expected exactly one mouse RNA-localization pathway row."
  )
}

split_leading_edge <- function(x) {
  
  if (
    length(x) == 0 ||
    is.na(x) ||
    x == ""
  ) {
    return(character(0))
  }
  
  stringr::str_split(
    x,
    pattern = ";",
    simplify = FALSE
  )[[1]] %>%
    unique()
}

human_rna_leading <- split_leading_edge(
  human_rna$leading_edge_genes[[1]]
)

mouse_rna_leading <- split_leading_edge(
  mouse_rna$leading_edge_genes[[1]]
)

shared_rna_leading <- intersect(
  human_rna_leading,
  mouse_rna_leading
) %>%
  sort()

# -------------------------------------------------------------------------
# 4. Build shared RNA-localization gene table
# -------------------------------------------------------------------------

rna_gene_table <- shared %>%
  dplyr::filter(
    human_gene_symbol %in%
      shared_rna_leading
  ) %>%
  dplyr::transmute(
    human_gene_symbol,
    mouse_gene_symbol,
    
    human_log2fc =
      as.numeric(scz_log2fc),
    
    human_t_stat =
      as.numeric(scz_t_stat),
    
    human_p_value =
      as.numeric(scz_p_value),
    
    human_fdr =
      as.numeric(scz_fdr),
    
    human_significant_fdr_0_05 =
      !is.na(human_fdr) &
      human_fdr < 0.05,
    
    human_significant_fdr_0_10 =
      !is.na(human_fdr) &
      human_fdr < 0.10,
    
    mouse_log2fc =
      as.numeric(mouse_log2fc),
    
    mouse_wald_stat =
      as.numeric(mouse_stat),
    
    mouse_p_value =
      as.numeric(mouse_p_value),
    
    mouse_padj =
      as.numeric(mouse_padj),
    
    mouse_significant_padj_0_05 =
      !is.na(mouse_padj) &
      mouse_padj < 0.05,
    
    mouse_significant_padj_0_10 =
      !is.na(mouse_padj) &
      mouse_padj < 0.10,
    
    direction_relation = dplyr::case_when(
      sign(human_log2fc) ==
        sign(mouse_log2fc) ~
        "same direction",
      
      sign(human_log2fc) !=
        sign(mouse_log2fc) ~
        "opposite direction",
      
      TRUE ~
        "unclassified"
    ),
    
    individually_significant_both_0_05 =
      human_significant_fdr_0_05 &
      mouse_significant_padj_0_05,
    
    individually_significant_either_0_05 =
      human_significant_fdr_0_05 |
      mouse_significant_padj_0_05,
    
    pathway =
      "RNA localization",
    
    shared_leading_edge =
      TRUE
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      abs(human_t_stat) +
        abs(mouse_wald_stat)
    ),
    human_gene_symbol
  )

if (nrow(rna_gene_table) !=
    length(shared_rna_leading)) {
  
  warning(
    "Not every shared leading-edge symbol was found ",
    "in the shared ortholog universe."
  )
}

readr::write_csv(
  rna_gene_table,
  file.path(
    table_dir,
    "rna_localization_shared_leading_edge_genes.csv"
  )
)

# -------------------------------------------------------------------------
# 5. Build RNA-localization pathway summary
# -------------------------------------------------------------------------

rna_pathway_summary <- tibble::tibble(
  pathway = "RNA localization",
  go_term = "GO:0006403",
  
  human_NES =
    human_rna$NES[[1]],
  
  human_p_value =
    human_rna$pval[[1]],
  
  human_fdr =
    human_rna$padj[[1]],
  
  mouse_NES =
    mouse_rna$NES[[1]],
  
  mouse_p_value =
    mouse_rna$pval[[1]],
  
  mouse_fdr =
    mouse_rna$padj[[1]],
  
  human_leading_edge_n =
    length(human_rna_leading),
  
  mouse_leading_edge_n =
    length(mouse_rna_leading),
  
  shared_leading_edge_n =
    length(shared_rna_leading),
  
  leading_edge_jaccard =
    length(shared_rna_leading) /
    length(
      union(
        human_rna_leading,
        mouse_rna_leading
      )
    ),
  
  nes_direction_relation =
    dplyr::case_when(
      sign(human_rna$NES[[1]]) ==
        sign(mouse_rna$NES[[1]]) ~
        "same NES direction",
      
      TRUE ~
        "opposite NES direction"
    ),
  
  shared_leading_edge_genes =
    paste(
      shared_rna_leading,
      collapse = ";"
    )
)

readr::write_csv(
  rna_pathway_summary,
  file.path(
    table_dir,
    "rna_localization_cross_species_summary.csv"
  )
)

# -------------------------------------------------------------------------
# 6. Build selected module-level summary
# -------------------------------------------------------------------------
#
# Each row represents a biologically interpretable comparison.
# Direction is included descriptively and is not used as a criterion
# for biological convergence.

selected_module_matches <- semantic_matches %>%
  dplyr::filter(
    best_match_type ==
      "best among significant human pathways"
  ) %>%
  dplyr::filter(
    mouse_pathway_name %in% c(
      "GOBP_RNA_LOCALIZATION",
      "GOBP_ESTABLISHMENT_OF_RNA_LOCALIZATION",
      "GOBP_CHOLESTEROL_BIOSYNTHETIC_PROCESS_VIA_DESMOSTEROL",
      "GOBP_STEROL_BIOSYNTHETIC_PROCESS",
      "GOBP_ISOPRENOID_BIOSYNTHETIC_PROCESS",
      "GOBP_PROTEIN_FOLDING"
    )
  ) %>%
  dplyr::mutate(
    biological_module =
      dplyr::case_when(
        stringr::str_detect(
          mouse_pathway_name,
          "RNA_LOCALIZATION"
        ) ~
          "RNA localization and transport",
        
        stringr::str_detect(
          mouse_pathway_name,
          "CHOLESTEROL|STEROL|ISOPRENOID"
        ) ~
          "Lipid and sterol metabolism",
        
        mouse_pathway_name ==
          "GOBP_PROTEIN_FOLDING" ~
          "Proteostasis and protein quality control",
        
        TRUE ~
          "Other"
      ),
    
    mouse_pathway_label =
      mouse_pathway_name %>%
      stringr::str_remove("^GOBP_") %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_to_sentence(),
    
    human_pathway_label =
      human_pathway_name %>%
      stringr::str_remove("^GOBP_") %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_to_sentence(),
    
    direction_is_descriptive_only =
      TRUE
  ) %>%
  dplyr::select(
    biological_module,
    
    mouse_pathway_name,
    mouse_pathway_label,
    mouse_go_id,
    mouse_NES,
    mouse_fdr,
    
    human_pathway_name,
    human_pathway_label,
    human_go_id,
    human_NES,
    human_fdr,
    
    semantic_similarity,
    leading_edge_shared_n,
    leading_edge_jaccard,
    shared_leading_edge_genes,
    nes_direction_relation,
    direction_is_descriptive_only
  ) %>%
  dplyr::arrange(
    factor(
      biological_module,
      levels = c(
        "RNA localization and transport",
        "Lipid and sterol metabolism",
        "Proteostasis and protein quality control"
      )
    ),
    mouse_fdr
  )

readr::write_csv(
  selected_module_matches,
  file.path(
    table_dir,
    "cross_species_selected_pathway_modules.csv"
  )
)

# -------------------------------------------------------------------------
# 7. Reconstruct native ranked vectors
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
names(human_ranks) <-
  human_rank_table$human_gene_symbol

mouse_ranks <- mouse_rank_table$mouse_stat
names(mouse_ranks) <-
  mouse_rank_table$human_gene_symbol

human_ranks <- sort(
  human_ranks,
  decreasing = TRUE
)

mouse_ranks <- sort(
  mouse_ranks,
  decreasing = TRUE
)

# -------------------------------------------------------------------------
# 8. Retrieve the RNA-localization gene set
# -------------------------------------------------------------------------

go_bp <- msigdbr::msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "C5",
  subcollection = "GO:BP"
)

rna_gene_set <- go_bp %>%
  dplyr::filter(
    gs_name == rna_pathway_name
  ) %>%
  dplyr::pull(gene_symbol) %>%
  unique()

rna_gene_set_human <- intersect(
  rna_gene_set,
  names(human_ranks)
)

rna_gene_set_mouse <- intersect(
  rna_gene_set,
  names(mouse_ranks)
)

# -------------------------------------------------------------------------
# 9. Create enrichment plots
# -------------------------------------------------------------------------

human_enrichment_plot <-
  fgsea::plotEnrichment(
    pathway = rna_gene_set_human,
    stats = human_ranks
  ) +
  ggplot2::labs(
    title =
      "Human schizophrenia synaptic proteome",
    subtitle =
      paste0(
        "RNA localization: NES = ",
        round(human_rna$NES[[1]], 2),
        "; FDR = ",
        signif(human_rna$padj[[1]], 2)
      ),
    x =
      "Genes ranked by proteomic t statistic",
    y =
      "Enrichment score"
  ) +
  ggplot2::theme_classic()

mouse_enrichment_plot <-
  fgsea::plotEnrichment(
    pathway = rna_gene_set_mouse,
    stats = mouse_ranks
  ) +
  ggplot2::labs(
    title =
      "Mouse neuronal C4-OE transcriptome",
    subtitle =
      paste0(
        "RNA localization: NES = ",
        round(mouse_rna$NES[[1]], 2),
        "; FDR = ",
        signif(mouse_rna$padj[[1]], 2)
      ),
    x =
      "Genes ranked by DESeq2 Wald statistic",
    y =
      "Enrichment score"
  ) +
  ggplot2::theme_classic()

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "human_rna_localization_enrichment.pdf"
  ),
  plot = human_enrichment_plot,
  width = 7,
  height = 5
)

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "human_rna_localization_enrichment.png"
  ),
  plot = human_enrichment_plot,
  width = 7,
  height = 5,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "mouse_rna_localization_enrichment.pdf"
  ),
  plot = mouse_enrichment_plot,
  width = 7,
  height = 5
)

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "mouse_rna_localization_enrichment.png"
  ),
  plot = mouse_enrichment_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# -------------------------------------------------------------------------
# 10. Create paired shared-gene statistic plot
# -------------------------------------------------------------------------

rna_plot_long <- rna_gene_table %>%
  dplyr::select(
    human_gene_symbol,
    human_t_stat,
    mouse_wald_stat
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      human_t_stat,
      mouse_wald_stat
    ),
    names_to = "dataset",
    values_to = "native_statistic"
  ) %>%
  dplyr::mutate(
    dataset = dplyr::recode(
      dataset,
      human_t_stat =
        "Human proteomic t statistic",
      mouse_wald_stat =
        "Mouse DESeq2 Wald statistic"
    )
  )

gene_order <- rna_gene_table %>%
  dplyr::mutate(
    maximum_absolute_statistic =
      pmax(
        abs(human_t_stat),
        abs(mouse_wald_stat),
        na.rm = TRUE
      )
  ) %>%
  dplyr::arrange(
    maximum_absolute_statistic
  ) %>%
  dplyr::pull(human_gene_symbol)

rna_plot_long <- rna_plot_long %>%
  dplyr::mutate(
    human_gene_symbol = factor(
      human_gene_symbol,
      levels = gene_order
    )
  )

rna_gene_plot <- ggplot2::ggplot(
  rna_plot_long,
  ggplot2::aes(
    x = native_statistic,
    y = human_gene_symbol,
    shape = dataset
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.4
  ) +
  ggplot2::geom_point(
    size = 2.5,
    alpha = 0.85,
    position = ggplot2::position_dodge(
      width = 0.5
    )
  ) +
  ggplot2::labs(
    title =
      "Shared RNA-localization leading-edge genes",
    subtitle =
      paste0(
        length(shared_rna_leading),
        " genes shared between human and mouse leading edges"
      ),
    x =
      "Native test statistic",
    y =
      NULL,
    shape =
      "Dataset"
  ) +
  ggplot2::theme_classic()

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "rna_localization_shared_gene_statistics.pdf"
  ),
  plot = rna_gene_plot,
  width = 8,
  height = 8
)

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "rna_localization_shared_gene_statistics.png"
  ),
  plot = rna_gene_plot,
  width = 8,
  height = 8,
  dpi = 300
)

# -------------------------------------------------------------------------
# 11. Create concise analysis summary
# -------------------------------------------------------------------------

analysis_summary <- tibble::tibble(
  metric = c(
    "Human RNA-localization NES",
    "Human RNA-localization FDR",
    "Mouse RNA-localization NES",
    "Mouse RNA-localization FDR",
    "Human RNA-localization leading-edge genes",
    "Mouse RNA-localization leading-edge genes",
    "Shared RNA-localization leading-edge genes",
    "RNA-localization leading-edge Jaccard index",
    "Shared genes individually significant in both at 0.05",
    "Shared genes individually significant in either dataset at 0.05",
    "Selected module comparisons retained"
  ),
  
  value = c(
    human_rna$NES[[1]],
    human_rna$padj[[1]],
    mouse_rna$NES[[1]],
    mouse_rna$padj[[1]],
    length(human_rna_leading),
    length(mouse_rna_leading),
    length(shared_rna_leading),
    rna_pathway_summary$leading_edge_jaccard[[1]],
    
    sum(
      rna_gene_table$
        individually_significant_both_0_05,
      na.rm = TRUE
    ),
    
    sum(
      rna_gene_table$
        individually_significant_either_0_05,
      na.rm = TRUE
    ),
    
    nrow(selected_module_matches)
  )
)

readr::write_csv(
  analysis_summary,
  file.path(
    log_dir,
    "cross_species_summary_analysis.csv"
  )
)

# -------------------------------------------------------------------------
# 12. Record software environment
# -------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    log_dir,
    "sessionInfo_08_build_cross_species_summary.txt"
  )
)

message(
  "Cross-species summary completed."
)

message(
  "RNA gene table: ",
  file.path(
    table_dir,
    "rna_localization_shared_leading_edge_genes.csv"
  )
)

message(
  "Module summary: ",
  file.path(
    table_dir,
    "cross_species_selected_pathway_modules.csv"
  )
)

message(
  "Review summary: ",
  file.path(
    log_dir,
    "cross_species_summary_analysis.csv"
  )
)