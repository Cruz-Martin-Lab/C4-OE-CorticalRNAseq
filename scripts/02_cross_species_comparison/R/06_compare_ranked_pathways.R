# 06_compare_ranked_pathways.R
#
# Purpose:
# Compare pathway-level enrichment between:
#   1. Human schizophrenia synaptic proteomics
#   2. Mouse C4-OE cortical transcriptomics
#
# Both datasets are analyzed using:
#   - the same high-confidence one-to-one ortholog universe;
#   - the same human GO Biological Process gene sets;
#   - the same signed significance ranking formula.
#
# Ranking score:
#   sign(log2 fold change) * -log10(nominal P value)
#
# Positive NES:
#   pathway genes are concentrated toward the positively changed end.
#
# Negative NES:
#   pathway genes are concentrated toward the negatively changed end.

library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(msigdbr)
library(fgsea)
library(ggplot2)

# -------------------------------------------------------------------------
# 1. Define paths
# -------------------------------------------------------------------------

input_file <- file.path(
  "data_processed",
  "human_mouse_shared_measurable_universe.csv"
)

processed_dir <- "data_processed"
table_dir <- file.path("results", "tables")
figure_dir <- file.path("results", "figures")
log_dir <- file.path("results", "logs")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Shared measurable universe not found: ", input_file)
}

# -------------------------------------------------------------------------
# 2. Read and validate the shared universe
# -------------------------------------------------------------------------

shared <- readr::read_csv(
  input_file,
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

if (anyDuplicated(shared$mouse_gene_symbol) > 0) {
  stop("Duplicated mouse gene symbols remain in the shared universe.")
}

# -------------------------------------------------------------------------
# 3. Create matched human and mouse rank vectors
# -------------------------------------------------------------------------


rank_table <- shared %>%
  dplyr::transmute(
    human_gene_symbol,
    mouse_gene_symbol,
    
    human_rank_score =
      as.numeric(scz_t_stat),
    
    mouse_rank_score =
      as.numeric(mouse_stat)
  )

readr::write_csv(
  rank_table,
  file.path(
    processed_dir,
    "human_mouse_matched_pathway_rank_scores.csv"
  )
)

human_rank_table <- rank_table %>%
  filter(
    !is.na(human_gene_symbol),
    !is.na(human_rank_score),
    is.finite(human_rank_score)
  ) %>%
  arrange(
    desc(human_rank_score),
    human_gene_symbol
  )

mouse_rank_table <- rank_table %>%
  filter(
    !is.na(human_gene_symbol),
    !is.na(mouse_rank_score),
    is.finite(mouse_rank_score)
  ) %>%
  arrange(
    desc(mouse_rank_score),
    human_gene_symbol
  )

human_ranks <- human_rank_table$human_rank_score
names(human_ranks) <- human_rank_table$human_gene_symbol

mouse_ranks <- mouse_rank_table$mouse_rank_score
names(mouse_ranks) <- mouse_rank_table$human_gene_symbol

human_ranks <- sort(human_ranks, decreasing = TRUE)
mouse_ranks <- sort(mouse_ranks, decreasing = TRUE)

if (anyDuplicated(names(human_ranks)) > 0) {
  stop("Human ranked vector contains duplicated gene symbols.")
}

if (anyDuplicated(names(mouse_ranks)) > 0) {
  stop("Mouse ranked vector contains duplicated human ortholog symbols.")
}

# Audit tied ranking values.
human_tie_fraction <- mean(duplicated(human_ranks))
mouse_tie_fraction <- mean(duplicated(mouse_ranks))

rank_summary <- tibble(
  dataset = c(
    "Human SCZ synaptic proteome",
    "Mouse C4-OE transcriptome"
  ),
  
  genes_ranked = c(
    length(human_ranks),
    length(mouse_ranks)
  ),
  
  minimum_score = c(
    min(human_ranks),
    min(mouse_ranks)
  ),
  
  maximum_score = c(
    max(human_ranks),
    max(mouse_ranks)
  ),
  
  fraction_of_tied_scores = c(
    human_tie_fraction,
    mouse_tie_fraction
  )
)

readr::write_csv(
  rank_summary,
  file.path(
    log_dir,
    "ranked_pathway_input_summary.csv"
  )
)

# -------------------------------------------------------------------------
# 4. Retrieve GO Biological Process gene sets
# -------------------------------------------------------------------------

collection_metadata <- msigdbr::msigdbr_collections(
  db_species = "HS"
)

readr::write_csv(
  collection_metadata,
  file.path(
    log_dir,
    "msigdbr_collection_metadata.csv"
  )
)

go_bp_raw <- msigdbr::msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "C5",
  subcollection = "GO:BP"
)

required_geneset_columns <- c(
  "gs_name",
  "gene_symbol"
)

missing_geneset_columns <- setdiff(
  required_geneset_columns,
  names(go_bp_raw)
)

if (length(missing_geneset_columns) > 0) {
  stop(
    "Required msigdbr columns are missing: ",
    paste(missing_geneset_columns, collapse = ", ")
  )
}

# Preserve GO identifiers where available.
go_id_column <- intersect(
  c("gs_id", "gs_exact_source"),
  names(go_bp_raw)
)

if (length(go_id_column) == 0) {
  go_bp_standardized <- go_bp_raw %>%
    dplyr::transmute(
      pathway_name = as.character(gs_name),
      go_id = NA_character_,
      gene_symbol = as.character(gene_symbol)
    )
} else {
  selected_go_column <- go_id_column[[1]]
  
  go_bp_standardized <- go_bp_raw %>%
    dplyr::transmute(
      pathway_name = as.character(gs_name),
      go_id = as.character(.data[[selected_go_column]]),
      gene_symbol = as.character(gene_symbol)
    )
}

go_bp_standardized <- go_bp_standardized %>%
  filter(
    !is.na(pathway_name),
    !is.na(gene_symbol),
    pathway_name != "",
    gene_symbol != ""
  ) %>%
  distinct(
    pathway_name,
    go_id,
    gene_symbol
  )

readr::write_csv(
  go_bp_standardized,
  file.path(
    processed_dir,
    "msigdb_human_go_biological_process_gene_sets.csv"
  )
)

# -------------------------------------------------------------------------
# 5. Restrict gene sets to the shared measurable universe
# -------------------------------------------------------------------------

shared_human_symbols <- intersect(
  names(human_ranks),
  names(mouse_ranks)
)

go_bp_measurable <- go_bp_standardized %>%
  filter(gene_symbol %in% shared_human_symbols) %>%
  distinct(
    pathway_name,
    go_id,
    gene_symbol
  )

pathway_size_audit <- go_bp_measurable %>%
  count(
    pathway_name,
    go_id,
    name = "measured_gene_count"
  ) %>%
  arrange(
    desc(measured_gene_count),
    pathway_name
  )

readr::write_csv(
  pathway_size_audit,
  file.path(
    log_dir,
    "go_bp_shared_universe_pathway_sizes.csv"
  )
)

# Use pathway names as fgsea list identifiers.
# GO identifiers remain attached through pathway_metadata.
pathways <- split(
  go_bp_measurable$gene_symbol,
  go_bp_measurable$pathway_name
)

pathways <- lapply(
  pathways,
  unique
)

pathway_metadata <- go_bp_measurable %>%
  dplyr::select(
    pathway_name,
    go_id
  ) %>%
  distinct()

# -------------------------------------------------------------------------
# 6. Run preranked GSEA
# -------------------------------------------------------------------------

set.seed(20260723)

human_fgsea_raw <- fgsea::fgseaMultilevel(
  pathways = pathways,
  stats = human_ranks,
  minSize = 15,
  maxSize = 500,
  eps = 0
)

set.seed(20260723)

mouse_fgsea_raw <- fgsea::fgseaMultilevel(
  pathways = pathways,
  stats = mouse_ranks,
  minSize = 15,
  maxSize = 500,
  eps = 0
)

# -------------------------------------------------------------------------
# 7. Standardize and save complete GSEA results
# -------------------------------------------------------------------------

standardize_fgsea <- function(
    fgsea_result,
    dataset_name
) {
  
  as.data.frame(fgsea_result) %>%
    as_tibble() %>%
    dplyr::rename(
      pathway_name = pathway
    ) %>%
    left_join(
      pathway_metadata,
      by = "pathway_name"
    ) %>%
    mutate(
      dataset = dataset_name,
      
      enrichment_direction = case_when(
        NES > 0 ~ "positive",
        NES < 0 ~ "negative",
        TRUE ~ "zero"
      ),
      
      significant_fdr_0_05 =
        !is.na(padj) & padj < 0.05,
      
      significant_fdr_0_10 =
        !is.na(padj) & padj < 0.10,
      
      leading_edge_genes = vapply(
        leadingEdge,
        function(x) {
          paste(x, collapse = ";")
        },
        character(1)
      )
    ) %>%
    dplyr::select(
      dataset,
      pathway_name,
      go_id,
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
    arrange(
      padj,
      desc(abs(NES)),
      pathway_name
    )
}

human_fgsea <- standardize_fgsea(
  human_fgsea_raw,
  "Human SCZ synaptic proteome"
)

mouse_fgsea <- standardize_fgsea(
  mouse_fgsea_raw,
  "Mouse C4-OE transcriptome"
)

readr::write_csv(
  human_fgsea,
  file.path(
    table_dir,
    "human_scz_go_bp_preranked_gsea.csv"
  )
)

readr::write_csv(
  mouse_fgsea,
  file.path(
    table_dir,
    "mouse_c4oe_go_bp_preranked_gsea.csv"
  )
)

# -------------------------------------------------------------------------
# 8. Match exact pathways across datasets
# -------------------------------------------------------------------------

pathway_comparison <- human_fgsea %>%
  dplyr::select(
    pathway_name,
    go_id,
    human_size = size,
    human_ES = ES,
    human_NES = NES,
    human_p_value = pval,
    human_fdr = padj,
    human_significant_fdr_0_05 = significant_fdr_0_05,
    human_significant_fdr_0_10 = significant_fdr_0_10,
    human_leading_edge_genes = leading_edge_genes
  ) %>%
  inner_join(
    mouse_fgsea %>%
      dplyr::select(
        pathway_name,
        mouse_size = size,
        mouse_ES = ES,
        mouse_NES = NES,
        mouse_p_value = pval,
        mouse_fdr = padj,
        mouse_significant_fdr_0_05 = significant_fdr_0_05,
        mouse_significant_fdr_0_10 = significant_fdr_0_10,
        mouse_leading_edge_genes = leading_edge_genes
      ),
    by = "pathway_name"
  ) %>%
  mutate(
    significant_both_fdr_0_05 =
      human_significant_fdr_0_05 &
      mouse_significant_fdr_0_05,
    
    significant_both_fdr_0_10 =
      human_significant_fdr_0_10 &
      mouse_significant_fdr_0_10,
    
    nes_direction_relation = case_when(
      sign(human_NES) == sign(mouse_NES) ~
        "same NES direction",
      
      sign(human_NES) != sign(mouse_NES) ~
        "opposite NES direction",
      
      TRUE ~
        "unclassified"
    ),
    
    absolute_nes_sum =
      abs(human_NES) + abs(mouse_NES)
  ) %>%
  arrange(
    desc(significant_both_fdr_0_05),
    desc(significant_both_fdr_0_10),
    human_fdr,
    mouse_fdr,
    desc(absolute_nes_sum)
  )

readr::write_csv(
  pathway_comparison,
  file.path(
    table_dir,
    "human_mouse_exact_go_bp_pathway_comparison.csv"
  )
)

shared_significant_pathways <- pathway_comparison %>%
  filter(significant_both_fdr_0_05)

readr::write_csv(
  shared_significant_pathways,
  file.path(
    table_dir,
    "human_mouse_shared_significant_go_bp_pathways_fdr005.csv"
  )
)

# -------------------------------------------------------------------------
# 9. Test overlap of significant pathways
# -------------------------------------------------------------------------

run_pathway_overlap_test <- function(
    comparison,
    human_flag,
    mouse_flag,
    analysis_name
) {
  
  human_sig <- comparison[[human_flag]]
  mouse_sig <- comparison[[mouse_flag]]
  
  both <- sum(human_sig & mouse_sig, na.rm = TRUE)
  human_only <- sum(human_sig & !mouse_sig, na.rm = TRUE)
  mouse_only <- sum(!human_sig & mouse_sig, na.rm = TRUE)
  neither <- sum(!human_sig & !mouse_sig, na.rm = TRUE)
  
  contingency <- matrix(
    c(
      both,
      human_only,
      mouse_only,
      neither
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Human = c("Significant", "Not_significant"),
      Mouse = c("Significant", "Not_significant")
    )
  )
  
  fisher_result <- fisher.test(
    contingency,
    alternative = "greater"
  )
  
  universe_size <- nrow(comparison)
  human_total <- sum(human_sig, na.rm = TRUE)
  mouse_total <- sum(mouse_sig, na.rm = TRUE)
  
  expected_overlap <-
    human_total * mouse_total / universe_size
  
  tibble(
    analysis = analysis_name,
    pathways_tested_in_both = universe_size,
    human_significant_pathways = human_total,
    mouse_significant_pathways = mouse_total,
    observed_shared_pathways = both,
    expected_shared_pathways = expected_overlap,
    enrichment_ratio = ifelse(
      expected_overlap > 0,
      both / expected_overlap,
      NA_real_
    ),
    fisher_odds_ratio =
      unname(fisher_result$estimate),
    fisher_p_value =
      fisher_result$p.value,
    fisher_confidence_lower =
      fisher_result$conf.int[[1]],
    fisher_confidence_upper =
      fisher_result$conf.int[[2]]
  )
}

pathway_overlap_summary <- bind_rows(
  run_pathway_overlap_test(
    pathway_comparison,
    "human_significant_fdr_0_05",
    "mouse_significant_fdr_0_05",
    "Exact GO:BP overlap at FDR < 0.05"
  ),
  
  run_pathway_overlap_test(
    pathway_comparison,
    "human_significant_fdr_0_10",
    "mouse_significant_fdr_0_10",
    "Exact GO:BP overlap at FDR < 0.10"
  )
)

readr::write_csv(
  pathway_overlap_summary,
  file.path(
    table_dir,
    "go_bp_pathway_overlap_fisher_test.csv"
  )
)

# -------------------------------------------------------------------------
# 10. Summarize NES direction among shared pathways
# -------------------------------------------------------------------------

shared_direction_summary <- pathway_comparison %>%
  filter(significant_both_fdr_0_05) %>%
  count(
    nes_direction_relation,
    name = "n_pathways"
  ) %>%
  mutate(
    total_shared_pathways =
      sum(n_pathways),
    
    proportion =
      n_pathways / total_shared_pathways
  )

readr::write_csv(
  shared_direction_summary,
  file.path(
    table_dir,
    "shared_pathway_nes_direction_summary.csv"
  )
)

# -------------------------------------------------------------------------
# 11. Plot human versus mouse pathway NES
# -------------------------------------------------------------------------

plot_data <- pathway_comparison %>%
  mutate(
    pathway_class = case_when(
      significant_both_fdr_0_05 ~
        "significant in both",
      
      human_significant_fdr_0_05 ~
        "human only",
      
      mouse_significant_fdr_0_05 ~
        "mouse only",
      
      TRUE ~
        "not significant"
    )
  )

pathway_scatter <- ggplot(
  plot_data,
  aes(
    x = mouse_NES,
    y = human_NES,
    shape = pathway_class
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
    alpha = 0.65,
    size = 1.8
  ) +
  labs(
    title =
      "GO Biological Process enrichment across species and modalities",
    subtitle =
      "Same one-to-one ortholog universe and pathway definitions",
    x =
      "Mouse C4-OE normalized enrichment score",
    y =
      "Human SCZ synaptic-proteome normalized enrichment score",
    shape =
      "FDR < 0.05"
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    figure_dir,
    "human_mouse_go_bp_nes_scatter.pdf"
  ),
  plot = pathway_scatter,
  width = 7,
  height = 6
)

ggsave(
  filename = file.path(
    figure_dir,
    "human_mouse_go_bp_nes_scatter.png"
  ),
  plot = pathway_scatter,
  width = 7,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------------
# 12. Create concise processing summary
# -------------------------------------------------------------------------

pathway_summary <- tibble(
  metric = c(
    "Shared ortholog genes available for pathway analysis",
    "Human genes in ranked vector",
    "Mouse genes in ranked vector",
    "GO Biological Process gene sets retrieved",
    "GO Biological Process gene sets with shared-universe genes",
    "GO pathways tested in human GSEA",
    "GO pathways tested in mouse GSEA",
    "Human significant pathways at FDR < 0.05",
    "Mouse significant pathways at FDR < 0.05",
    "Shared significant pathways at FDR < 0.05",
    "Shared pathways with same NES direction at FDR < 0.05",
    "Shared pathways with opposite NES direction at FDR < 0.05"
  ),
  
  value = c(
    nrow(shared),
    length(human_ranks),
    length(mouse_ranks),
    n_distinct(go_bp_standardized$pathway_name),
    length(pathways),
    nrow(human_fgsea),
    nrow(mouse_fgsea),
    sum(
      human_fgsea$significant_fdr_0_05,
      na.rm = TRUE
    ),
    sum(
      mouse_fgsea$significant_fdr_0_05,
      na.rm = TRUE
    ),
    sum(
      pathway_comparison$significant_both_fdr_0_05,
      na.rm = TRUE
    ),
    sum(
      pathway_comparison$significant_both_fdr_0_05 &
        pathway_comparison$nes_direction_relation ==
        "same NES direction",
      na.rm = TRUE
    ),
    sum(
      pathway_comparison$significant_both_fdr_0_05 &
        pathway_comparison$nes_direction_relation ==
        "opposite NES direction",
      na.rm = TRUE
    )
  )
)

readr::write_csv(
  pathway_summary,
  file.path(
    log_dir,
    "ranked_pathway_comparison_summary.csv"
  )
)

# -------------------------------------------------------------------------
# 13. Record software environment
# -------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    log_dir,
    "sessionInfo_06_compare_ranked_pathways.txt"
  )
)

message("Ranked pathway comparison completed.")

message(
  "Human GSEA results: ",
  file.path(
    table_dir,
    "human_scz_go_bp_preranked_gsea.csv"
  )
)

message(
  "Mouse GSEA results: ",
  file.path(
    table_dir,
    "mouse_c4oe_go_bp_preranked_gsea.csv"
  )
)

message(
  "Exact pathway comparison: ",
  file.path(
    table_dir,
    "human_mouse_exact_go_bp_pathway_comparison.csv"
  )
)

message(
  "Review summary: ",
  file.path(
    log_dir,
    "ranked_pathway_comparison_summary.csv"
  )
)