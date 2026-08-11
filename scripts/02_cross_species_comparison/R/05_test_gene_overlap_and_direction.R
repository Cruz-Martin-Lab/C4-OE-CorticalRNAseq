# 05_test_gene_overlap_and_direction.R
#
# Purpose:
# Test cross-species overlap between:
#   1. Human schizophrenia synaptic proteomic changes
#   2. Mouse C4-OE cortical transcriptomic changes
#
# Primary analysis:
#   Human SCZ DEP: FDR < 0.10
#   Mouse C4-OE DEG: padj < 0.05
#
# Sensitivity analysis:
#   Human SCZ DEP: FDR < 0.05
#   Mouse C4-OE DEG: padj < 0.05
#
# The analysis uses the high-confidence one-to-one shared
# measurable ortholog universe created in Script 04.

library(readr)
library(dplyr)
library(tibble)
library(ggplot2)

# -------------------------------------------------------------------------
# 1. Define paths
# -------------------------------------------------------------------------

input_file <- file.path(
  "data_processed",
  "human_mouse_shared_measurable_universe.csv"
)

table_dir <- file.path("results", "tables")
figure_dir <- file.path("results", "figures")
log_dir <- file.path("results", "logs")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Shared measurable universe not found: ", input_file)
}

# -------------------------------------------------------------------------
# 2. Read the shared measurable universe
# -------------------------------------------------------------------------

shared <- readr::read_csv(
  input_file,
  show_col_types = FALSE
)

required_columns <- c(
  "human_gene_symbol",
  "mouse_gene_symbol",
  "scz_log2fc",
  "scz_fdr",
  "scz_dep_fdr_0_10",
  "scz_dep_fdr_0_05",
  "scz_direction",
  "mouse_log2fc",
  "mouse_padj",
  "mouse_deg_padj_0_05",
  "mouse_direction"
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
  stop("Shared universe contains duplicated human gene symbols.")
}

if (anyDuplicated(shared$mouse_gene_symbol) > 0) {
  stop("Shared universe contains duplicated mouse gene symbols.")
}

# -------------------------------------------------------------------------
# 3. Add primary and sensitivity significance classifications
# -------------------------------------------------------------------------

shared_classified <- shared %>%
  mutate(
    human_sig_primary =
      !is.na(scz_dep_fdr_0_10) &
      scz_dep_fdr_0_10,
    
    human_sig_sensitivity =
      !is.na(scz_dep_fdr_0_05) &
      scz_dep_fdr_0_05,
    
    mouse_sig =
      !is.na(mouse_deg_padj_0_05) &
      mouse_deg_padj_0_05,
    
    overlap_primary =
      human_sig_primary & mouse_sig,
    
    overlap_sensitivity =
      human_sig_sensitivity & mouse_sig,
    
    direction_relation = case_when(
      is.na(scz_log2fc) |
        is.na(mouse_log2fc) ~
        NA_character_,
      
      sign(scz_log2fc) ==
        sign(mouse_log2fc) ~
        "concordant",
      
      sign(scz_log2fc) !=
        sign(mouse_log2fc) ~
        "discordant",
      
      TRUE ~
        "unclassified"
    ),
    
    directional_category = case_when(
      scz_log2fc > 0 & mouse_log2fc > 0 ~
        "human up / mouse up",
      
      scz_log2fc < 0 & mouse_log2fc < 0 ~
        "human down / mouse down",
      
      scz_log2fc > 0 & mouse_log2fc < 0 ~
        "human up / mouse down",
      
      scz_log2fc < 0 & mouse_log2fc > 0 ~
        "human down / mouse up",
      
      TRUE ~
        "unclassified"
    )
  )

readr::write_csv(
  shared_classified,
  file.path(
    table_dir,
    "shared_universe_with_overlap_classification.csv"
  )
)

# -------------------------------------------------------------------------
# 4. Function for Fisher exact overlap testing
# -------------------------------------------------------------------------

run_overlap_test <- function(
    data,
    human_sig_column,
    analysis_name
) {
  
  human_sig <- data[[human_sig_column]]
  mouse_sig <- data$mouse_sig
  
  both <- sum(
    human_sig & mouse_sig,
    na.rm = TRUE
  )
  
  human_only <- sum(
    human_sig & !mouse_sig,
    na.rm = TRUE
  )
  
  mouse_only <- sum(
    !human_sig & mouse_sig,
    na.rm = TRUE
  )
  
  neither <- sum(
    !human_sig & !mouse_sig,
    na.rm = TRUE
  )
  
  contingency_matrix <- matrix(
    c(
      both,
      human_only,
      mouse_only,
      neither
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Human = c(
        "Significant",
        "Not_significant"
      ),
      Mouse = c(
        "Significant",
        "Not_significant"
      )
    )
  )
  
  fisher_result <- fisher.test(
    contingency_matrix,
    alternative = "greater"
  )
  
  human_total <- sum(
    human_sig,
    na.rm = TRUE
  )
  
  mouse_total <- sum(
    mouse_sig,
    na.rm = TRUE
  )
  
  universe_size <- nrow(data)
  
  expected_overlap <- (
    human_total * mouse_total
  ) / universe_size
  
  enrichment_ratio <- ifelse(
    expected_overlap > 0,
    both / expected_overlap,
    NA_real_
  )
  
  list(
    contingency = as.data.frame.matrix(
      contingency_matrix
    ) %>%
      rownames_to_column("human_status"),
    
    summary = tibble(
      analysis = analysis_name,
      universe_size = universe_size,
      human_significant = human_total,
      mouse_significant = mouse_total,
      observed_overlap = both,
      expected_overlap = expected_overlap,
      enrichment_ratio = enrichment_ratio,
      fisher_odds_ratio =
        unname(fisher_result$estimate),
      fisher_p_value =
        fisher_result$p.value,
      fisher_confidence_lower =
        fisher_result$conf.int[[1]],
      fisher_confidence_upper =
        fisher_result$conf.int[[2]]
    )
  )
}

# -------------------------------------------------------------------------
# 5. Run primary overlap test
# -------------------------------------------------------------------------

primary_test <- run_overlap_test(
  data = shared_classified,
  human_sig_column = "human_sig_primary",
  analysis_name =
    "Human FDR < 0.10; mouse padj < 0.05"
)

# -------------------------------------------------------------------------
# 6. Run sensitivity overlap test
# -------------------------------------------------------------------------

sensitivity_test <- run_overlap_test(
  data = shared_classified,
  human_sig_column = "human_sig_sensitivity",
  analysis_name =
    "Human FDR < 0.05; mouse padj < 0.05"
)

overlap_test_summary <- bind_rows(
  primary_test$summary,
  sensitivity_test$summary
)

readr::write_csv(
  overlap_test_summary,
  file.path(
    table_dir,
    "gene_overlap_fisher_test_summary.csv"
  )
)

readr::write_csv(
  primary_test$contingency,
  file.path(
    log_dir,
    "gene_overlap_primary_contingency_table.csv"
  )
)

readr::write_csv(
  sensitivity_test$contingency,
  file.path(
    log_dir,
    "gene_overlap_sensitivity_contingency_table.csv"
  )
)

# -------------------------------------------------------------------------
# 7. Export complete primary overlap table
# -------------------------------------------------------------------------

primary_overlap <- shared_classified %>%
  filter(overlap_primary) %>%
  dplyr::select(
    human_gene_symbol,
    mouse_gene_symbol,
    human_ensembl_gene_id,
    mouse_ensembl_gene_id,
    scz_log2fc,
    scz_fdr,
    scz_direction,
    mouse_log2fc,
    mouse_padj,
    mouse_direction,
    direction_relation,
    directional_category,
    base_mean,
    num_spectra,
    num_unique_peptides,
    percent_coverage,
    accession_number,
    entry_name
  ) %>%
  arrange(
    direction_relation,
    scz_fdr,
    mouse_padj,
    human_gene_symbol
  )

readr::write_csv(
  primary_overlap,
  file.path(
    table_dir,
    "primary_cross_species_gene_overlap.csv"
  )
)

# -------------------------------------------------------------------------
# 8. Export complete sensitivity overlap table
# -------------------------------------------------------------------------

sensitivity_overlap <- shared_classified %>%
  filter(overlap_sensitivity) %>%
  dplyr::select(
    human_gene_symbol,
    mouse_gene_symbol,
    human_ensembl_gene_id,
    mouse_ensembl_gene_id,
    scz_log2fc,
    scz_fdr,
    scz_direction,
    mouse_log2fc,
    mouse_padj,
    mouse_direction,
    direction_relation,
    directional_category,
    base_mean,
    num_spectra,
    num_unique_peptides,
    percent_coverage,
    accession_number,
    entry_name
  ) %>%
  arrange(
    direction_relation,
    scz_fdr,
    mouse_padj,
    human_gene_symbol
  )

readr::write_csv(
  sensitivity_overlap,
  file.path(
    table_dir,
    "sensitivity_cross_species_gene_overlap.csv"
  )
)

# -------------------------------------------------------------------------
# 9. Summarize directional concordance
# -------------------------------------------------------------------------

summarize_direction <- function(
    overlap_table,
    analysis_name
) {
  
  overlap_table %>%
    count(
      direction_relation,
      directional_category,
      name = "n_genes"
    ) %>%
    mutate(
      analysis = analysis_name,
      total_overlap = nrow(overlap_table),
      proportion_of_overlap =
        n_genes / total_overlap
    ) %>%
    dplyr::select(
      analysis,
      direction_relation,
      directional_category,
      n_genes,
      total_overlap,
      proportion_of_overlap
    )
}

direction_summary <- bind_rows(
  summarize_direction(
    primary_overlap,
    "Human FDR < 0.10; mouse padj < 0.05"
  ),
  
  summarize_direction(
    sensitivity_overlap,
    "Human FDR < 0.05; mouse padj < 0.05"
  )
)

readr::write_csv(
  direction_summary,
  file.path(
    table_dir,
    "cross_species_direction_summary.csv"
  )
)

# -------------------------------------------------------------------------
# 10. Test whether concordance exceeds 50%
# -------------------------------------------------------------------------
#
# This is a descriptive directional test among overlapping genes.
# It asks whether concordant direction occurs more often than a
# simple 50:50 concordant-versus-discordant expectation.

run_direction_binomial_test <- function(
    overlap_table,
    analysis_name
) {
  
  n_concordant <- sum(
    overlap_table$direction_relation ==
      "concordant",
    na.rm = TRUE
  )
  
  n_discordant <- sum(
    overlap_table$direction_relation ==
      "discordant",
    na.rm = TRUE
  )
  
  n_directional <- n_concordant + n_discordant
  
  if (n_directional == 0) {
    return(
      tibble(
        analysis = analysis_name,
        concordant = 0,
        discordant = 0,
        proportion_concordant = NA_real_,
        binomial_p_value = NA_real_,
        confidence_lower = NA_real_,
        confidence_upper = NA_real_
      )
    )
  }
  
  binom_result <- binom.test(
    x = n_concordant,
    n = n_directional,
    p = 0.5,
    alternative = "two.sided"
  )
  
  tibble(
    analysis = analysis_name,
    concordant = n_concordant,
    discordant = n_discordant,
    proportion_concordant =
      n_concordant / n_directional,
    binomial_p_value =
      binom_result$p.value,
    confidence_lower =
      binom_result$conf.int[[1]],
    confidence_upper =
      binom_result$conf.int[[2]]
  )
}

direction_test_summary <- bind_rows(
  run_direction_binomial_test(
    primary_overlap,
    "Human FDR < 0.10; mouse padj < 0.05"
  ),
  
  run_direction_binomial_test(
    sensitivity_overlap,
    "Human FDR < 0.05; mouse padj < 0.05"
  )
)

readr::write_csv(
  direction_test_summary,
  file.path(
    table_dir,
    "cross_species_direction_binomial_test.csv"
  )
)

# -------------------------------------------------------------------------
# 11. Create scatter plot of overlapping genes
# -------------------------------------------------------------------------

if (nrow(primary_overlap) > 0) {
  
  overlap_plot <- ggplot(
    primary_overlap,
    aes(
      x = mouse_log2fc,
      y = scz_log2fc,
      shape = direction_relation
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.4
    ) +
    geom_point(
      size = 2.5,
      alpha = 0.8
    ) +
    labs(
      title =
        "Cross-species overlap of C4-OE DEGs and SCZ synaptic DEPs",
      subtitle =
        "Human FDR < 0.10; mouse padj < 0.05",
      x =
        "Mouse C4-OE log2 fold change",
      y =
        "Human SCZ synaptic-proteome log2 fold change",
      shape =
        "Direction"
    ) +
    theme_classic()
  
  ggsave(
    filename = file.path(
      figure_dir,
      "primary_cross_species_overlap_scatter.pdf"
    ),
    plot = overlap_plot,
    width = 7,
    height = 6
  )
  
  ggsave(
    filename = file.path(
      figure_dir,
      "primary_cross_species_overlap_scatter.png"
    ),
    plot = overlap_plot,
    width = 7,
    height = 6,
    dpi = 300
  )
}

# -------------------------------------------------------------------------
# 12. Create concise analysis summary
# -------------------------------------------------------------------------

primary_direction_counts <- primary_overlap %>%
  count(
    direction_relation,
    name = "n"
  )

sensitivity_direction_counts <- sensitivity_overlap %>%
  count(
    direction_relation,
    name = "n"
  )

analysis_summary <- tibble(
  metric = c(
    "Shared measurable ortholog universe",
    "Human SCZ DEPs at FDR < 0.10",
    "Mouse C4-OE DEGs at padj < 0.05",
    "Primary observed overlap",
    "Primary expected overlap",
    "Primary overlap enrichment ratio",
    "Primary Fisher exact P value",
    "Primary concordant overlap genes",
    "Primary discordant overlap genes",
    "Human SCZ DEPs at FDR < 0.05",
    "Sensitivity observed overlap",
    "Sensitivity expected overlap",
    "Sensitivity overlap enrichment ratio",
    "Sensitivity Fisher exact P value",
    "Sensitivity concordant overlap genes",
    "Sensitivity discordant overlap genes"
  ),
  
  value = c(
    nrow(shared_classified),
    
    sum(
      shared_classified$human_sig_primary,
      na.rm = TRUE
    ),
    
    sum(
      shared_classified$mouse_sig,
      na.rm = TRUE
    ),
    
    primary_test$summary$observed_overlap,
    
    primary_test$summary$expected_overlap,
    
    primary_test$summary$enrichment_ratio,
    
    primary_test$summary$fisher_p_value,
    
    sum(
      primary_overlap$direction_relation ==
        "concordant",
      na.rm = TRUE
    ),
    
    sum(
      primary_overlap$direction_relation ==
        "discordant",
      na.rm = TRUE
    ),
    
    sum(
      shared_classified$human_sig_sensitivity,
      na.rm = TRUE
    ),
    
    sensitivity_test$summary$observed_overlap,
    
    sensitivity_test$summary$expected_overlap,
    
    sensitivity_test$summary$enrichment_ratio,
    
    sensitivity_test$summary$fisher_p_value,
    
    sum(
      sensitivity_overlap$direction_relation ==
        "concordant",
      na.rm = TRUE
    ),
    
    sum(
      sensitivity_overlap$direction_relation ==
        "discordant",
      na.rm = TRUE
    )
  )
)

readr::write_csv(
  analysis_summary,
  file.path(
    log_dir,
    "gene_overlap_analysis_summary.csv"
  )
)

# -------------------------------------------------------------------------
# 13. Record software environment
# -------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    log_dir,
    "sessionInfo_05_test_gene_overlap_and_direction.txt"
  )
)

message(
  "Cross-species gene-overlap analysis completed."
)

message(
  "Primary overlap table: ",
  file.path(
    table_dir,
    "primary_cross_species_gene_overlap.csv"
  )
)

message(
  "Statistical summary: ",
  file.path(
    table_dir,
    "gene_overlap_fisher_test_summary.csv"
  )
)

message(
  "Review summary: ",
  file.path(
    log_dir,
    "gene_overlap_analysis_summary.csv"
  )
)