# 10_compare_hallmark_reactome.R
#
# Purpose:
# Validate the cross-species pathway results using two less GO-dependent
# MSigDB collections:
#
#   1. Hallmark pathways
#   2. Reactome pathways
#
# Ranking statistics:
#   Human: schizophrenia-versus-control proteomic t statistic
#   Mouse: DESeq2 Wald statistic
#
# Direction is reported descriptively and is not required for defining
# biological convergence.

library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
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

table_dir <- file.path("results", "tables")
figure_dir <- file.path("results", "figures")
log_dir <- file.path("results", "logs")
processed_dir <- "data_processed"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(shared_file)) {
  stop("Shared measurable universe not found: ", shared_file)
}

# -------------------------------------------------------------------------
# 2. Read and validate the shared ortholog universe
# -------------------------------------------------------------------------

shared <- readr::read_csv(
  shared_file,
  show_col_types = FALSE
)

required_columns <- c(
  "human_gene_symbol",
  "mouse_gene_symbol",
  "scz_t_stat",
  "mouse_stat"
)

missing_columns <- setdiff(
  required_columns,
  names(shared)
)

if (length(missing_columns) > 0) {
  stop(
    "Required columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (anyDuplicated(shared$human_gene_symbol) > 0) {
  stop("Duplicated human gene symbols remain in the shared universe.")
}

# -------------------------------------------------------------------------
# 3. Build native-statistic ranked vectors
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
# 4. Retrieve Hallmark and Reactome collections
# -------------------------------------------------------------------------

hallmark_raw <- msigdbr::msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "H"
)

reactome_raw <- msigdbr::msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:REACTOME"
)

standardize_collection <- function(
    raw_data,
    collection_name
) {
  
  required_gene_set_columns <- c(
    "gs_name",
    "gene_symbol"
  )
  
  missing_gene_set_columns <- setdiff(
    required_gene_set_columns,
    names(raw_data)
  )
  
  if (length(missing_gene_set_columns) > 0) {
    stop(
      "Missing gene-set columns for ",
      collection_name,
      ": ",
      paste(
        missing_gene_set_columns,
        collapse = ", "
      )
    )
  }
  
  raw_data %>%
    dplyr::transmute(
      collection = collection_name,
      pathway_name = as.character(gs_name),
      gene_symbol = as.character(gene_symbol)
    ) %>%
    dplyr::filter(
      !is.na(pathway_name),
      !is.na(gene_symbol),
      pathway_name != "",
      gene_symbol != "",
      gene_symbol %in% shared_symbols
    ) %>%
    dplyr::distinct()
}

hallmark_standardized <- standardize_collection(
  hallmark_raw,
  "Hallmark"
)

reactome_standardized <- standardize_collection(
  reactome_raw,
  "Reactome"
)

readr::write_csv(
  hallmark_standardized,
  file.path(
    processed_dir,
    "msigdb_hallmark_shared_universe.csv"
  )
)

readr::write_csv(
  reactome_standardized,
  file.path(
    processed_dir,
    "msigdb_reactome_shared_universe.csv"
  )
)

# -------------------------------------------------------------------------
# 5. Convert gene sets into fgsea pathway lists
# -------------------------------------------------------------------------

make_pathway_list <- function(collection_data) {
  
  pathway_list <- split(
    collection_data$gene_symbol,
    collection_data$pathway_name
  )
  
  lapply(
    pathway_list,
    unique
  )
}

hallmark_pathways <- make_pathway_list(
  hallmark_standardized
)

reactome_pathways <- make_pathway_list(
  reactome_standardized
)

# -------------------------------------------------------------------------
# 6. Run fgsea
# -------------------------------------------------------------------------

run_collection_fgsea <- function(
    pathways,
    ranks,
    minimum_size,
    maximum_size
) {
  
  set.seed(20260723)
  
  fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = minimum_size,
    maxSize = maximum_size,
    eps = 0
  )
}

human_hallmark_raw <- run_collection_fgsea(
  hallmark_pathways,
  human_ranks,
  minimum_size = 10,
  maximum_size = 500
)

mouse_hallmark_raw <- run_collection_fgsea(
  hallmark_pathways,
  mouse_ranks,
  minimum_size = 10,
  maximum_size = 500
)

human_reactome_raw <- run_collection_fgsea(
  reactome_pathways,
  human_ranks,
  minimum_size = 15,
  maximum_size = 500
)

mouse_reactome_raw <- run_collection_fgsea(
  reactome_pathways,
  mouse_ranks,
  minimum_size = 15,
  maximum_size = 500
)

# -------------------------------------------------------------------------
# 7. Standardize fgsea outputs
# -------------------------------------------------------------------------

standardize_fgsea <- function(
    fgsea_result,
    dataset_name,
    collection_name
) {
  
  as.data.frame(fgsea_result) %>%
    tibble::as_tibble() %>%
    dplyr::rename(
      pathway_name = pathway
    ) %>%
    dplyr::mutate(
      dataset = dataset_name,
      collection = collection_name,
      
      enrichment_direction =
        dplyr::case_when(
          NES > 0 ~ "positive",
          NES < 0 ~ "negative",
          TRUE ~ "zero"
        ),
      
      significant_fdr_0_05 =
        !is.na(padj) &
        padj < 0.05,
      
      significant_fdr_0_10 =
        !is.na(padj) &
        padj < 0.10,
      
      leading_edge_genes =
        vapply(
          leadingEdge,
          function(x) {
            paste(
              x,
              collapse = ";"
            )
          },
          character(1)
        )
    ) %>%
    dplyr::select(
      dataset,
      collection,
      pathway_name,
      size,
      ES,
      NES,
      pval,
      padj,
      log2err,
      enrichment_direction,
      significant_fdr_0_05,
      significant_fdr_0_10,
      leading_edge_genes
    ) %>%
    dplyr::arrange(
      padj,
      dplyr::desc(abs(NES)),
      pathway_name
    )
}

human_hallmark <- standardize_fgsea(
  human_hallmark_raw,
  "Human SCZ synaptic proteome",
  "Hallmark"
)

mouse_hallmark <- standardize_fgsea(
  mouse_hallmark_raw,
  "Mouse C4-OE transcriptome",
  "Hallmark"
)

human_reactome <- standardize_fgsea(
  human_reactome_raw,
  "Human SCZ synaptic proteome",
  "Reactome"
)

mouse_reactome <- standardize_fgsea(
  mouse_reactome_raw,
  "Mouse C4-OE transcriptome",
  "Reactome"
)

readr::write_csv(
  human_hallmark,
  file.path(
    table_dir,
    "human_scz_hallmark_preranked_gsea.csv"
  )
)

readr::write_csv(
  mouse_hallmark,
  file.path(
    table_dir,
    "mouse_c4oe_hallmark_preranked_gsea.csv"
  )
)

readr::write_csv(
  human_reactome,
  file.path(
    table_dir,
    "human_scz_reactome_preranked_gsea.csv"
  )
)

readr::write_csv(
  mouse_reactome,
  file.path(
    table_dir,
    "mouse_c4oe_reactome_preranked_gsea.csv"
  )
)

# -------------------------------------------------------------------------
# 8. Match exact pathways between datasets
# -------------------------------------------------------------------------

compare_collection <- function(
    human_result,
    mouse_result,
    collection_name
) {
  
  human_result %>%
    dplyr::select(
      pathway_name,
      human_size = size,
      human_NES = NES,
      human_p_value = pval,
      human_fdr = padj,
      human_significant_fdr_0_05 =
        significant_fdr_0_05,
      human_significant_fdr_0_10 =
        significant_fdr_0_10,
      human_leading_edge_genes =
        leading_edge_genes
    ) %>%
    dplyr::inner_join(
      mouse_result %>%
        dplyr::select(
          pathway_name,
          mouse_size = size,
          mouse_NES = NES,
          mouse_p_value = pval,
          mouse_fdr = padj,
          mouse_significant_fdr_0_05 =
            significant_fdr_0_05,
          mouse_significant_fdr_0_10 =
            significant_fdr_0_10,
          mouse_leading_edge_genes =
            leading_edge_genes
        ),
      by = "pathway_name"
    ) %>%
    dplyr::mutate(
      collection = collection_name,
      
      significant_both_fdr_0_05 =
        human_significant_fdr_0_05 &
        mouse_significant_fdr_0_05,
      
      significant_both_fdr_0_10 =
        human_significant_fdr_0_10 &
        mouse_significant_fdr_0_10,
      
      nes_direction_relation =
        dplyr::case_when(
          sign(human_NES) ==
            sign(mouse_NES) ~
            "same NES direction",
          
          TRUE ~
            "opposite NES direction"
        ),
      
      absolute_nes_sum =
        abs(human_NES) +
        abs(mouse_NES)
    ) %>%
    dplyr::arrange(
      dplyr::desc(
        significant_both_fdr_0_05
      ),
      dplyr::desc(
        significant_both_fdr_0_10
      ),
      human_fdr,
      mouse_fdr,
      dplyr::desc(
        absolute_nes_sum
      )
    )
}

hallmark_comparison <- compare_collection(
  human_hallmark,
  mouse_hallmark,
  "Hallmark"
)

reactome_comparison <- compare_collection(
  human_reactome,
  mouse_reactome,
  "Reactome"
)

readr::write_csv(
  hallmark_comparison,
  file.path(
    table_dir,
    "human_mouse_hallmark_pathway_comparison.csv"
  )
)

readr::write_csv(
  reactome_comparison,
  file.path(
    table_dir,
    "human_mouse_reactome_pathway_comparison.csv"
  )
)

# -------------------------------------------------------------------------
# 9. Test overlap of significant pathways
# -------------------------------------------------------------------------

run_overlap_test <- function(
    comparison,
    human_flag,
    mouse_flag,
    collection_name,
    threshold_name
) {
  
  human_sig <- comparison[[human_flag]]
  mouse_sig <- comparison[[mouse_flag]]
  
  both <- sum(
    human_sig &
      mouse_sig,
    na.rm = TRUE
  )
  
  human_only <- sum(
    human_sig &
      !mouse_sig,
    na.rm = TRUE
  )
  
  mouse_only <- sum(
    !human_sig &
      mouse_sig,
    na.rm = TRUE
  )
  
  neither <- sum(
    !human_sig &
      !mouse_sig,
    na.rm = TRUE
  )
  
  contingency <- matrix(
    c(
      both,
      human_only,
      mouse_only,
      neither
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  fisher_result <- fisher.test(
    contingency,
    alternative = "greater"
  )
  
  pathways_tested <- nrow(comparison)
  human_total <- sum(human_sig, na.rm = TRUE)
  mouse_total <- sum(mouse_sig, na.rm = TRUE)
  
  expected_overlap <-
    human_total *
    mouse_total /
    pathways_tested
  
  tibble::tibble(
    collection = collection_name,
    threshold = threshold_name,
    pathways_tested = pathways_tested,
    human_significant_pathways =
      human_total,
    mouse_significant_pathways =
      mouse_total,
    observed_shared_pathways =
      both,
    expected_shared_pathways =
      expected_overlap,
    enrichment_ratio =
      ifelse(
        expected_overlap > 0,
        both / expected_overlap,
        NA_real_
      ),
    fisher_odds_ratio =
      unname(fisher_result$estimate),
    fisher_p_value =
      fisher_result$p.value
  )
}

overlap_summary <- dplyr::bind_rows(
  run_overlap_test(
    hallmark_comparison,
    "human_significant_fdr_0_05",
    "mouse_significant_fdr_0_05",
    "Hallmark",
    "FDR < 0.05"
  ),
  
  run_overlap_test(
    hallmark_comparison,
    "human_significant_fdr_0_10",
    "mouse_significant_fdr_0_10",
    "Hallmark",
    "FDR < 0.10"
  ),
  
  run_overlap_test(
    reactome_comparison,
    "human_significant_fdr_0_05",
    "mouse_significant_fdr_0_05",
    "Reactome",
    "FDR < 0.05"
  ),
  
  run_overlap_test(
    reactome_comparison,
    "human_significant_fdr_0_10",
    "mouse_significant_fdr_0_10",
    "Reactome",
    "FDR < 0.10"
  )
)

readr::write_csv(
  overlap_summary,
  file.path(
    table_dir,
    "hallmark_reactome_pathway_overlap_tests.csv"
  )
)

# -------------------------------------------------------------------------
# 10. Calculate global NES correlations
# -------------------------------------------------------------------------

calculate_global_correlation <- function(
    comparison,
    collection_name
) {
  
  analysis_data <- comparison %>%
    dplyr::filter(
      !is.na(human_NES),
      !is.na(mouse_NES),
      is.finite(human_NES),
      is.finite(mouse_NES)
    )
  
  spearman_result <- cor.test(
    analysis_data$human_NES,
    analysis_data$mouse_NES,
    method = "spearman",
    exact = FALSE
  )
  
  tibble::tibble(
    collection = collection_name,
    pathways_compared =
      nrow(analysis_data),
    
    spearman_rho =
      unname(
        spearman_result$estimate
      ),
    
    spearman_p_value =
      spearman_result$p.value,
    
    same_direction_pathways =
      sum(
        sign(analysis_data$human_NES) ==
          sign(analysis_data$mouse_NES),
        na.rm = TRUE
      ),
    
    same_direction_proportion =
      mean(
        sign(analysis_data$human_NES) ==
          sign(analysis_data$mouse_NES),
        na.rm = TRUE
      )
  )
}

global_correlation_summary <- dplyr::bind_rows(
  calculate_global_correlation(
    hallmark_comparison,
    "Hallmark"
  ),
  
  calculate_global_correlation(
    reactome_comparison,
    "Reactome"
  )
)

readr::write_csv(
  global_correlation_summary,
  file.path(
    log_dir,
    "hallmark_reactome_global_nes_correlations.csv"
  )
)

# -------------------------------------------------------------------------
# 11. Extract biologically relevant pathway families
# -------------------------------------------------------------------------

relevant_pattern <- paste(
  c(
    "RNA",
    "MRNA",
    "RIBONUCLEOPROTEIN",
    "NUCLEAR_EXPORT",
    "NUCLEAR_IMPORT",
    "TRANSPORT_OF_RNA",
    "CHOLESTEROL",
    "STEROL",
    "LIPID",
    "FATTY_ACID",
    "UNFOLDED_PROTEIN",
    "PROTEIN_FOLDING",
    "CHAPERON",
    "AUTOPHAG",
    "PROTEASOM"
  ),
  collapse = "|"
)

focused_pathways <- dplyr::bind_rows(
  hallmark_comparison,
  reactome_comparison
) %>%
  dplyr::filter(
    stringr::str_detect(
      pathway_name,
      relevant_pattern
    )
  ) %>%
  dplyr::arrange(
    collection,
    human_fdr,
    mouse_fdr,
    dplyr::desc(
      absolute_nes_sum
    )
  )

readr::write_csv(
  focused_pathways,
  file.path(
    table_dir,
    "hallmark_reactome_focused_cross_species_pathways.csv"
  )
)

# -------------------------------------------------------------------------
# 12. Create concise summary
# -------------------------------------------------------------------------

analysis_summary <- tibble::tibble(
  metric = c(
    "Hallmark pathways tested",
    "Human Hallmark pathways FDR < 0.05",
    "Mouse Hallmark pathways FDR < 0.05",
    "Shared Hallmark pathways FDR < 0.05",
    "Reactome pathways tested",
    "Human Reactome pathways FDR < 0.05",
    "Mouse Reactome pathways FDR < 0.05",
    "Shared Reactome pathways FDR < 0.05",
    "Focused RNA/lipid/proteostasis comparisons retained"
  ),
  
  value = c(
    nrow(hallmark_comparison),
    
    sum(
      hallmark_comparison$
        human_significant_fdr_0_05,
      na.rm = TRUE
    ),
    
    sum(
      hallmark_comparison$
        mouse_significant_fdr_0_05,
      na.rm = TRUE
    ),
    
    sum(
      hallmark_comparison$
        significant_both_fdr_0_05,
      na.rm = TRUE
    ),
    
    nrow(reactome_comparison),
    
    sum(
      reactome_comparison$
        human_significant_fdr_0_05,
      na.rm = TRUE
    ),
    
    sum(
      reactome_comparison$
        mouse_significant_fdr_0_05,
      na.rm = TRUE
    ),
    
    sum(
      reactome_comparison$
        significant_both_fdr_0_05,
      na.rm = TRUE
    ),
    
    nrow(focused_pathways)
  )
)

readr::write_csv(
  analysis_summary,
  file.path(
    log_dir,
    "hallmark_reactome_analysis_summary.csv"
  )
)

# -------------------------------------------------------------------------
# 13. Create scatter plots
# -------------------------------------------------------------------------

make_scatter <- function(
    comparison,
    collection_name
) {
  
  correlation <- cor.test(
    comparison$human_NES,
    comparison$mouse_NES,
    method = "spearman",
    exact = FALSE
  )
  
  ggplot2::ggplot(
    comparison,
    ggplot2::aes(
      x = mouse_NES,
      y = human_NES
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.4
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      alpha = 0.7,
      size = 2
    ) +
    ggplot2::labs(
      title =
        paste0(
          collection_name,
          " pathway comparison"
        ),
      
      subtitle =
        paste0(
          "Spearman rho = ",
          round(
            unname(
              correlation$estimate
            ),
            3
          ),
          "; P = ",
          signif(
            correlation$p.value,
            3
          )
        ),
      
      x =
        "Mouse C4-OE NES",
      
      y =
        "Human schizophrenia synaptic-proteome NES"
    ) +
    ggplot2::theme_classic()
}

hallmark_plot <- make_scatter(
  hallmark_comparison,
  "Hallmark"
)

reactome_plot <- make_scatter(
  reactome_comparison,
  "Reactome"
)

ggplot2::ggsave(
  file.path(
    figure_dir,
    "human_mouse_hallmark_nes_scatter.pdf"
  ),
  hallmark_plot,
  width = 7,
  height = 6
)

ggplot2::ggsave(
  file.path(
    figure_dir,
    "human_mouse_hallmark_nes_scatter.png"
  ),
  hallmark_plot,
  width = 7,
  height = 6,
  dpi = 300
)

ggplot2::ggsave(
  file.path(
    figure_dir,
    "human_mouse_reactome_nes_scatter.pdf"
  ),
  reactome_plot,
  width = 7,
  height = 6
)

ggplot2::ggsave(
  file.path(
    figure_dir,
    "human_mouse_reactome_nes_scatter.png"
  ),
  reactome_plot,
  width = 7,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------------
# 14. Record software environment
# -------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    log_dir,
    "sessionInfo_10_compare_hallmark_reactome.txt"
  )
)

message(
  "Hallmark and Reactome comparison completed."
)

message(
  "Analysis summary: ",
  file.path(
    log_dir,
    "hallmark_reactome_analysis_summary.csv"
  )
)

message(
  "Global correlations: ",
  file.path(
    log_dir,
    "hallmark_reactome_global_nes_correlations.csv"
  )
)

message(
  "Focused pathways: ",
  file.path(
    table_dir,
    "hallmark_reactome_focused_cross_species_pathways.csv"
  )
)