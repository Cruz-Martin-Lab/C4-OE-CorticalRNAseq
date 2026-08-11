# 07_compare_related_go_pathways.R
#
# Purpose:
# Compare biologically related GO Biological Process terms between:
#   1. Human schizophrenia synaptic proteomics
#   2. Mouse C4-OE cortical transcriptomics
#
# This is an exploratory ontology-aware analysis.
# It does not replace the exact GO-term overlap test from Script 06.
#
# Primary exploratory question:
# For each mouse pathway significant at FDR < 0.10, what are the most
# semantically related human pathways?
#
# Similarity measures:
#   - GO semantic similarity using the Wang method
#   - Jaccard overlap of leading-edge genes
#
# Human pathways are reported in two ways:
#   1. Best match among all tested human pathways
#   2. Best match among human pathways significant at FDR < 0.05

library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(msigdbr)
library(GOSemSim)
library(GO.db)
library(org.Hs.eg.db)
library(ggplot2)

# -------------------------------------------------------------------------
# 1. Define paths
# -------------------------------------------------------------------------

human_file <- file.path(
  "results",
  "tables",
  "human_scz_go_bp_preranked_gsea.csv"
)

mouse_file <- file.path(
  "results",
  "tables",
  "mouse_c4oe_go_bp_preranked_gsea.csv"
)

processed_dir <- "data_processed"
table_dir <- file.path("results", "tables")
figure_dir <- file.path("results", "figures")
log_dir <- file.path("results", "logs")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(human_file)) {
  stop("Human GSEA file not found: ", human_file)
}

if (!file.exists(mouse_file)) {
  stop("Mouse GSEA file not found: ", mouse_file)
}

# -------------------------------------------------------------------------
# 2. Read GSEA results
# -------------------------------------------------------------------------

human_gsea <- readr::read_csv(
  human_file,
  show_col_types = FALSE
)

mouse_gsea <- readr::read_csv(
  mouse_file,
  show_col_types = FALSE
)

required_columns <- c(
  "pathway_name",
  "NES",
  "pval",
  "padj",
  "leading_edge_genes"
)

missing_human <- setdiff(
  required_columns,
  names(human_gsea)
)

missing_mouse <- setdiff(
  required_columns,
  names(mouse_gsea)
)

if (length(missing_human) > 0) {
  stop(
    "Missing human GSEA columns: ",
    paste(missing_human, collapse = ", ")
  )
}

if (length(missing_mouse) > 0) {
  stop(
    "Missing mouse GSEA columns: ",
    paste(missing_mouse, collapse = ", ")
  )
}

# -------------------------------------------------------------------------
# 3. Retrieve genuine GO identifiers
# -------------------------------------------------------------------------
#
# In msigdbr, gs_id may be an MSigDB identifier such as M50425.
# gs_exact_source contains the original source identifier, which for
# GO collections should be a GO identifier such as GO:0006695.

go_bp_raw <- msigdbr::msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "C5",
  subcollection = "GO:BP"
)

required_msig_columns <- c(
  "gs_name",
  "gs_exact_source",
  "gene_symbol"
)

missing_msig_columns <- setdiff(
  required_msig_columns,
  names(go_bp_raw)
)

if (length(missing_msig_columns) > 0) {
  stop(
    "Required msigdbr columns are missing: ",
    paste(missing_msig_columns, collapse = ", ")
  )
}

go_metadata <- go_bp_raw %>%
  dplyr::transmute(
    pathway_name = as.character(gs_name),
    go_id = as.character(gs_exact_source),
    gene_symbol = as.character(gene_symbol)
  ) %>%
  filter(
    !is.na(pathway_name),
    !is.na(go_id),
    !is.na(gene_symbol),
    pathway_name != "",
    gene_symbol != "",
    str_detect(go_id, "^GO:[0-9]{7}$")
  ) %>%
  distinct()

readr::write_csv(
  go_metadata,
  file.path(
    processed_dir,
    "msigdb_go_bp_with_true_go_ids.csv"
  )
)

pathway_to_go <- go_metadata %>%
  dplyr::select(
    pathway_name,
    go_id
  ) %>%
  distinct()

# Audit pathways associated with more than one GO identifier.
pathway_go_counts <- pathway_to_go %>%
  count(
    pathway_name,
    name = "n_go_ids"
  )

ambiguous_pathway_ids <- pathway_to_go %>%
  semi_join(
    pathway_go_counts %>%
      filter(n_go_ids > 1),
    by = "pathway_name"
  ) %>%
  arrange(
    pathway_name,
    go_id
  )

readr::write_csv(
  ambiguous_pathway_ids,
  file.path(
    log_dir,
    "go_pathways_with_multiple_go_ids.csv"
  )
)

# Retain pathways with exactly one true GO identifier.
pathway_to_go_unique <- pathway_to_go %>%
  inner_join(
    pathway_go_counts %>%
      filter(n_go_ids == 1),
    by = "pathway_name"
  ) %>%
  dplyr::select(
    pathway_name,
    go_id
  )

# -------------------------------------------------------------------------
# 4. Add true GO identifiers to GSEA results
# -------------------------------------------------------------------------

human_go <- human_gsea %>%
  dplyr::select(
    pathway_name,
    human_NES = NES,
    human_p_value = pval,
    human_fdr = padj,
    human_leading_edge = leading_edge_genes
  ) %>%
  inner_join(
    pathway_to_go_unique,
    by = "pathway_name"
  ) %>%
  distinct(
    pathway_name,
    go_id,
    .keep_all = TRUE
  )

mouse_go <- mouse_gsea %>%
  dplyr::select(
    pathway_name,
    mouse_NES = NES,
    mouse_p_value = pval,
    mouse_fdr = padj,
    mouse_leading_edge = leading_edge_genes
  ) %>%
  inner_join(
    pathway_to_go_unique,
    by = "pathway_name"
  ) %>%
  distinct(
    pathway_name,
    go_id,
    .keep_all = TRUE
  )

mouse_focus <- mouse_go %>%
  filter(
    !is.na(mouse_fdr),
    mouse_fdr < 0.10
  ) %>%
  arrange(
    mouse_fdr,
    desc(abs(mouse_NES)),
    pathway_name
  )

human_significant <- human_go %>%
  filter(
    !is.na(human_fdr),
    human_fdr < 0.05
  )

if (nrow(mouse_focus) == 0) {
  stop("No mouse GO pathways were significant at FDR < 0.10.")
}

if (nrow(human_significant) == 0) {
  warning("No human GO pathways were significant at FDR < 0.05.")
}

readr::write_csv(
  mouse_focus,
  file.path(
    table_dir,
    "mouse_go_pathways_fdr010_with_true_go_ids.csv"
  )
)

# -------------------------------------------------------------------------
# 5. Build GO Biological Process semantic data
# -------------------------------------------------------------------------

go_semantic_data <- GOSemSim::godata(
  annoDb = "org.Hs.eg.db",
  ont = "BP",
  computeIC = FALSE
)

# -------------------------------------------------------------------------
# 6. Calculate GO semantic similarity
# -------------------------------------------------------------------------
#
# Wang similarity is topology-based and does not require information
# content. Scores generally range from 0 to 1, with larger values
# indicating closer ontology relationships.

mouse_go_ids <- unique(mouse_focus$go_id)
human_go_ids <- unique(human_go$go_id)

semantic_matrix <- GOSemSim::mgoSim(
  mouse_go_ids,
  human_go_ids,
  semData = go_semantic_data,
  measure = "Wang",
  combine = NULL
)

if (is.null(dim(semantic_matrix))) {
  semantic_matrix <- matrix(
    semantic_matrix,
    nrow = length(mouse_go_ids),
    ncol = length(human_go_ids),
    dimnames = list(
      mouse_go_ids,
      human_go_ids
    )
  )
}

semantic_long <- as.data.frame(
  as.table(semantic_matrix),
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  dplyr::rename(
    mouse_go_id = Var1,
    human_go_id = Var2,
    semantic_similarity = Freq
  ) %>%
  mutate(
    mouse_go_id = as.character(mouse_go_id),
    human_go_id = as.character(human_go_id),
    semantic_similarity =
      as.numeric(semantic_similarity)
  ) %>%
  filter(
    !is.na(semantic_similarity)
  )

# -------------------------------------------------------------------------
# 7. Join pathway statistics
# -------------------------------------------------------------------------

comparison_all <- semantic_long %>%
  inner_join(
    mouse_focus %>%
      dplyr::rename(
        mouse_pathway_name = pathway_name,
        mouse_go_id = go_id
      ),
    by = "mouse_go_id"
  ) %>%
  inner_join(
    human_go %>%
      dplyr::rename(
        human_pathway_name = pathway_name,
        human_go_id = go_id
      ),
    by = "human_go_id"
  ) %>%
  mutate(
    nes_direction_relation = case_when(
      sign(mouse_NES) == sign(human_NES) ~
        "same NES direction",
      
      sign(mouse_NES) != sign(human_NES) ~
        "opposite NES direction",
      
      TRUE ~
        "unclassified"
    ),
    
    human_significant_fdr_0_05 =
      !is.na(human_fdr) &
      human_fdr < 0.05
  )

# -------------------------------------------------------------------------
# 8. Calculate leading-edge gene overlap
# -------------------------------------------------------------------------

split_gene_string <- function(x) {
  
  if (
    is.na(x) ||
    x == ""
  ) {
    return(character(0))
  }
  
  unique(
    str_split(
      x,
      pattern = ";",
      simplify = FALSE
    )[[1]]
  )
}

calculate_gene_overlap <- function(
    mouse_string,
    human_string
) {
  
  mouse_genes <- split_gene_string(mouse_string)
  human_genes <- split_gene_string(human_string)
  
  shared_genes <- intersect(
    mouse_genes,
    human_genes
  )
  
  union_genes <- union(
    mouse_genes,
    human_genes
  )
  
  tibble(
    leading_edge_shared_n =
      length(shared_genes),
    
    leading_edge_jaccard =
      ifelse(
        length(union_genes) > 0,
        length(shared_genes) /
          length(union_genes),
        NA_real_
      ),
    
    shared_leading_edge_genes =
      paste(
        sort(shared_genes),
        collapse = ";"
      )
  )
}

gene_overlap_results <- lapply(
  seq_len(nrow(comparison_all)),
  function(i) {
    
    calculate_gene_overlap(
      comparison_all$mouse_leading_edge[[i]],
      comparison_all$human_leading_edge[[i]]
    )
  }
) %>%
  bind_rows()

comparison_all <- bind_cols(
  comparison_all,
  gene_overlap_results
) %>%
  mutate(
    combined_similarity_score =
      semantic_similarity +
      leading_edge_jaccard
  )

readr::write_csv(
  comparison_all,
  file.path(
    table_dir,
    "all_mouse_human_go_semantic_comparisons.csv"
  )
)

# -------------------------------------------------------------------------
# 9. Select top human matches for each mouse pathway
# -------------------------------------------------------------------------

top_matches_all_human <- comparison_all %>%
  group_by(
    mouse_pathway_name,
    mouse_go_id
  ) %>%
  arrange(
    desc(semantic_similarity),
    desc(leading_edge_jaccard),
    human_fdr,
    human_pathway_name,
    .by_group = TRUE
  ) %>%
  slice_head(n = 10) %>%
  mutate(
    match_rank = row_number(),
    match_scope = "all tested human pathways"
  ) %>%
  ungroup()

top_matches_significant_human <- comparison_all %>%
  filter(
    human_significant_fdr_0_05
  ) %>%
  group_by(
    mouse_pathway_name,
    mouse_go_id
  ) %>%
  arrange(
    desc(semantic_similarity),
    desc(leading_edge_jaccard),
    human_fdr,
    human_pathway_name,
    .by_group = TRUE
  ) %>%
  slice_head(n = 10) %>%
  mutate(
    match_rank = row_number(),
    match_scope =
      "human pathways significant at FDR < 0.05"
  ) %>%
  ungroup()

top_matches_combined <- bind_rows(
  top_matches_all_human,
  top_matches_significant_human
) %>%
  dplyr::select(
    match_scope,
    match_rank,
    
    mouse_pathway_name,
    mouse_go_id,
    mouse_NES,
    mouse_p_value,
    mouse_fdr,
    
    human_pathway_name,
    human_go_id,
    human_NES,
    human_p_value,
    human_fdr,
    
    semantic_similarity,
    leading_edge_shared_n,
    leading_edge_jaccard,
    shared_leading_edge_genes,
    nes_direction_relation,
    combined_similarity_score
  ) %>%
  arrange(
    mouse_fdr,
    mouse_pathway_name,
    match_scope,
    match_rank
  )

readr::write_csv(
  top_matches_combined,
  file.path(
    table_dir,
    "top_related_human_go_pathways_for_mouse_fdr010.csv"
  )
)

# -------------------------------------------------------------------------
# 10. Create one best-match table
# -------------------------------------------------------------------------

best_all_human <- top_matches_all_human %>%
  filter(match_rank == 1) %>%
  mutate(
    best_match_type =
      "best among all human pathways"
  )

best_significant_human <- top_matches_significant_human %>%
  filter(match_rank == 1) %>%
  mutate(
    best_match_type =
      "best among significant human pathways"
  )

best_match_summary <- bind_rows(
  best_all_human,
  best_significant_human
) %>%
  dplyr::select(
    best_match_type,
    
    mouse_pathway_name,
    mouse_go_id,
    mouse_NES,
    mouse_fdr,
    
    human_pathway_name,
    human_go_id,
    human_NES,
    human_fdr,
    
    semantic_similarity,
    leading_edge_shared_n,
    leading_edge_jaccard,
    shared_leading_edge_genes,
    nes_direction_relation
  ) %>%
  arrange(
    mouse_fdr,
    mouse_pathway_name,
    best_match_type
  )

readr::write_csv(
  best_match_summary,
  file.path(
    table_dir,
    "best_related_human_go_pathway_per_mouse_pathway.csv"
  )
)

# -------------------------------------------------------------------------
# 11. Assign broad mouse biological themes
# -------------------------------------------------------------------------
#
# These labels summarize the significant mouse terms only.
# They do not restrict which human pathways are considered.

best_match_themed <- best_match_summary %>%
  mutate(
    mouse_biological_theme = case_when(
      str_detect(
        mouse_pathway_name,
        "CHOLESTEROL|STEROL|STEROID|ISOPRENOID|ALCOHOL"
      ) ~
        "sterol and lipid biosynthesis",
      
      str_detect(
        mouse_pathway_name,
        "PROTEIN_FOLDING|PROTEIN_MATURATION|ENDOPLASMIC_RETICULUM|PROTEASOM"
      ) ~
        "proteostasis and ER stress",
      
      str_detect(
        mouse_pathway_name,
        "RNA_LOCALIZATION|RNA_EXPORT|MRNA_TRANSPORT|NUCLEAR_TRANSPORT"
      ) ~
        "RNA localization and transport",
      
      str_detect(
        mouse_pathway_name,
        "AXON|NEURON|SYNAP|PRESYNAPTIC"
      ) ~
        "neuronal and synaptic remodeling",
      
      TRUE ~
        "other"
    )
  )

readr::write_csv(
  best_match_themed,
  file.path(
    table_dir,
    "best_related_go_pathways_with_mouse_themes.csv"
  )
)

# -------------------------------------------------------------------------
# 12. Plot top significant-human matches
# -------------------------------------------------------------------------

plot_data <- top_matches_significant_human %>%
  filter(match_rank <= 5) %>%
  mutate(
    comparison_label = paste0(
      str_replace_all(
        mouse_pathway_name,
        "^GOBP_",
        ""
      ),
      "\n-> ",
      str_replace_all(
        human_pathway_name,
        "^GOBP_",
        ""
      )
    )
  ) %>%
  arrange(
    semantic_similarity
  ) %>%
  mutate(
    comparison_label = factor(
      comparison_label,
      levels = unique(comparison_label)
    )
  )

if (nrow(plot_data) > 0) {
  
  similarity_plot <- ggplot(
    plot_data,
    aes(
      x = semantic_similarity,
      y = comparison_label,
      shape = nes_direction_relation
    )
  ) +
    geom_point(
      size = 2.5,
      alpha = 0.8
    ) +
    labs(
      title =
        "Ontology-related human pathways for mouse C4-OE programs",
      subtitle =
        "Mouse FDR < 0.10; human pathways FDR < 0.05",
      x =
        "GO Biological Process semantic similarity (Wang)",
      y =
        NULL,
      shape =
        "NES relationship"
    ) +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 7)
    )
  
  ggsave(
    filename = file.path(
      figure_dir,
      "mouse_human_go_semantic_similarity.pdf"
    ),
    plot = similarity_plot,
    width = 10,
    height = 10
  )
  
  ggsave(
    filename = file.path(
      figure_dir,
      "mouse_human_go_semantic_similarity.png"
    ),
    plot = similarity_plot,
    width = 10,
    height = 10,
    dpi = 300
  )
}

# -------------------------------------------------------------------------
# 13. Create processing summary
# -------------------------------------------------------------------------

analysis_summary <- tibble(
  metric = c(
    "Mouse pathways significant at FDR < 0.10",
    "Human pathways tested with unique true GO IDs",
    "Human pathways significant at FDR < 0.05 with true GO IDs",
    "Total mouse-human semantic comparisons",
    "Mouse pathways with a significant-human match",
    "Best significant-human matches with semantic similarity >= 0.70",
    "Best significant-human matches sharing leading-edge genes"
  ),
  
  value = c(
    nrow(mouse_focus),
    nrow(human_go),
    nrow(human_significant),
    nrow(comparison_all),
    
    n_distinct(
      top_matches_significant_human$
        mouse_pathway_name
    ),
    
    sum(
      best_significant_human$
        semantic_similarity >= 0.70,
      na.rm = TRUE
    ),
    
    sum(
      best_significant_human$
        leading_edge_shared_n > 0,
      na.rm = TRUE
    )
  )
)

readr::write_csv(
  analysis_summary,
  file.path(
    log_dir,
    "go_semantic_comparison_summary.csv"
  )
)

# -------------------------------------------------------------------------
# 14. Record software environment
# -------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    log_dir,
    "sessionInfo_07_compare_related_go_pathways.txt"
  )
)

message(
  "Exploratory GO semantic comparison completed."
)

message(
  "Best-match table: ",
  file.path(
    table_dir,
    "best_related_human_go_pathway_per_mouse_pathway.csv"
  )
)

message(
  "Top-match table: ",
  file.path(
    table_dir,
    "top_related_human_go_pathways_for_mouse_fdr010.csv"
  )
)

message(
  "Review summary: ",
  file.path(
    log_dir,
    "go_semantic_comparison_summary.csv"
  )
)