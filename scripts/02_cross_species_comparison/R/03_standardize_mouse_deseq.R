# 03_standardize_mouse_deseq.R
#
# Purpose:
# Standardize the complete mouse C4-OE DESeq2 results table.
#
# This script preserves the full set of tested genes so it can later
# define the mouse analysis universe for cross-species overlap testing.

library(readr)
library(dplyr)
library(stringr)
library(tibble)

# -------------------------------------------------------------------------
# 1. Define paths
# -------------------------------------------------------------------------

# The mouse input is the DESeq2 table produced by strand 01
# (scripts/01_bulk_rnaseq_DE), design ~ sex + genotype. run_full_analysis.R
# copies it into data_raw/ before this script runs, so the cross-species
# comparison always uses the DEG table produced by this pipeline.
input_file <- file.path(
  "data_raw",
  "S1_Table_DESeq.csv"
)

processed_dir <- "data_processed"
log_dir <- file.path("results", "logs")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# -------------------------------------------------------------------------
# 2. Read the complete DESeq2 table
# -------------------------------------------------------------------------

mouse_raw <- readr::read_csv(
  input_file,
  show_col_types = FALSE
)

# Column names in the mouse table. The gene-symbol column is `gene`, and there
# are no per-sample count columns -- nothing downstream in R/04-15 ever read
# those, so they are simply absent rather than carried along.
required_columns <- c(
  "gene",
  "baseMean",
  "logFC",
  "lfcSE",
  "stat",
  "PValue",
  "padj",
  # Listed here so that supplying a DEG table without these columns produces the
  # clear error below rather than dying later inside transmute() with
  # "object 'sex_confound_flag' not found".
  "sex_confound_flag",
  "sex_linked"
)

missing_columns <- setdiff(
  required_columns,
  names(mouse_raw)
)

if (length(missing_columns) > 0) {
  stop(
    "Required columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

# -------------------------------------------------------------------------
# 3. Create a standardized mouse DESeq2 table
# -------------------------------------------------------------------------

mouse_standardized <- mouse_raw %>%
  transmute(
    mouse_gene_symbol = str_trim(
      as.character(gene)
    ),

    base_mean = suppressWarnings(
      as.numeric(baseMean)
    ),

    mouse_log2fc = suppressWarnings(
      as.numeric(logFC)
    ),

    mouse_lfc_se = suppressWarnings(
      as.numeric(lfcSE)
    ),

    mouse_stat = suppressWarnings(
      as.numeric(stat)
    ),

    mouse_p_value = suppressWarnings(
      as.numeric(PValue)
    ),

    mouse_padj = suppressWarnings(
      as.numeric(padj)
    ),

    # Carried through deliberately. transmute() drops anything not named, and
    # these two carry the sex diagnostics through: sex_confound_flag marks genes
    # that are also sex-differential within Control, sex_linked marks X/Y genes.
    # Downstream scripts can then exclude or annotate them rather than silently
    # treating them as clean genotype effects.
    mouse_sex_confound_flag = as.logical(sex_confound_flag),
    mouse_sex_linked        = as.logical(sex_linked)
  ) %>%
  mutate(
    mouse_gene_symbol = na_if(
      mouse_gene_symbol,
      ""
    ),
    
    mouse_gene_symbol = if_else(
      str_to_lower(mouse_gene_symbol) %in%
        c("na", "n/a", "nan"),
      NA_character_,
      mouse_gene_symbol
    ),
    
    mouse_direction = case_when(
      is.na(mouse_log2fc) ~ NA_character_,
      mouse_log2fc > 0 ~ "up",
      mouse_log2fc < 0 ~ "down",
      TRUE ~ "unchanged"
    ),
    
    mouse_deg_padj_0_05 =
      !is.na(mouse_padj) &
      mouse_padj < 0.05,
    
    mouse_deg_padj_0_10 =
      !is.na(mouse_padj) &
      mouse_padj < 0.10
  )

# Save all tested rows, including rows without usable symbols.
readr::write_csv(
  mouse_standardized,
  file.path(
    processed_dir,
    "mouse_c4oe_deseq_all_tested_rows.csv"
  )
)

# -------------------------------------------------------------------------
# 4. Separate rows with and without usable gene symbols
# -------------------------------------------------------------------------

mouse_with_gene <- mouse_standardized %>%
  filter(!is.na(mouse_gene_symbol))

mouse_without_gene <- mouse_standardized %>%
  filter(is.na(mouse_gene_symbol))

readr::write_csv(
  mouse_without_gene,
  file.path(
    log_dir,
    "mouse_deseq_rows_without_gene_symbol.csv"
  )
)

# -------------------------------------------------------------------------
# 5. Audit duplicated mouse gene symbols
# -------------------------------------------------------------------------

mouse_duplicate_gene_summary <- mouse_with_gene %>%
  count(
    mouse_gene_symbol,
    name = "n_rows"
  ) %>%
  filter(n_rows > 1) %>%
  arrange(
    desc(n_rows),
    mouse_gene_symbol
  )

mouse_duplicate_gene_rows <- mouse_with_gene %>%
  semi_join(
    mouse_duplicate_gene_summary,
    by = "mouse_gene_symbol"
  ) %>%
  left_join(
    mouse_duplicate_gene_summary,
    by = "mouse_gene_symbol"
  ) %>%
  arrange(
    mouse_gene_symbol,
    mouse_padj,
    desc(base_mean)
  )

readr::write_csv(
  mouse_duplicate_gene_summary,
  file.path(
    log_dir,
    "mouse_deseq_duplicate_gene_summary.csv"
  )
)

readr::write_csv(
  mouse_duplicate_gene_rows,
  file.path(
    log_dir,
    "mouse_deseq_duplicate_gene_rows.csv"
  )
)

# -------------------------------------------------------------------------
# 6. Create one row per mouse gene
# -------------------------------------------------------------------------
#
# Duplicate rows, if present, are resolved without selecting on statistical
# significance. The row with the highest baseMean is retained.
# The gene symbol is used as the deterministic final tie-breaker.

mouse_gene_level <- mouse_with_gene %>%
  group_by(mouse_gene_symbol) %>%
  mutate(
    n_rows_per_gene = n()
  ) %>%
  arrange(
    desc(base_mean),
    mouse_gene_symbol,
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(
    duplicate_resolution = if_else(
      n_rows_per_gene == 1,
      "single DESeq2 row",
      paste0(
        "selected highest baseMean from ",
        n_rows_per_gene,
        " rows"
      )
    )
  ) %>%
  arrange(mouse_gene_symbol)

if (
  anyDuplicated(
    mouse_gene_level$mouse_gene_symbol
  ) > 0
) {
  stop(
    "Mouse gene-level output still contains duplicated gene symbols."
  )
}

# -------------------------------------------------------------------------
# 7. Save final mouse gene-level table
# -------------------------------------------------------------------------

readr::write_csv(
  mouse_gene_level,
  file.path(
    processed_dir,
    "mouse_c4oe_deseq_gene_level.csv"
  )
)

# -------------------------------------------------------------------------
# 8. Create processing summary
# -------------------------------------------------------------------------

mouse_summary <- tibble(
  metric = c(
    "Raw DESeq2 rows",
    "Rows with a usable mouse gene symbol",
    "Rows without a usable mouse gene symbol",
    "Unique mouse genes before duplicate collapse",
    "Mouse genes with multiple rows",
    "Final mouse gene-level rows",
    "Genes with non-missing adjusted P value",
    "Mouse DEGs at padj < 0.05",
    "Mouse DEGs at padj < 0.10",
    "Mouse DEGs up at padj < 0.05",
    "Mouse DEGs down at padj < 0.05"
  ),
  
  value = c(
    nrow(mouse_raw),
    nrow(mouse_with_gene),
    nrow(mouse_without_gene),
    n_distinct(
      mouse_with_gene$mouse_gene_symbol
    ),
    nrow(mouse_duplicate_gene_summary),
    nrow(mouse_gene_level),
    sum(
      !is.na(mouse_gene_level$mouse_padj)
    ),
    sum(
      mouse_gene_level$mouse_deg_padj_0_05,
      na.rm = TRUE
    ),
    sum(
      mouse_gene_level$mouse_deg_padj_0_10,
      na.rm = TRUE
    ),
    sum(
      mouse_gene_level$mouse_deg_padj_0_05 &
        mouse_gene_level$mouse_direction == "up",
      na.rm = TRUE
    ),
    sum(
      mouse_gene_level$mouse_deg_padj_0_05 &
        mouse_gene_level$mouse_direction == "down",
      na.rm = TRUE
    )
  )
)

readr::write_csv(
  mouse_summary,
  file.path(
    log_dir,
    "mouse_deseq_standardization_summary.csv"
  )
)

# -------------------------------------------------------------------------
# 9. Save provisional QC genes for inspection only
# -------------------------------------------------------------------------

mouse_qc_genes <- mouse_gene_level %>%
  filter(
    mouse_gene_symbol %in%
      c("Hcn1", "Kcnab2")
  )

readr::write_csv(
  mouse_qc_genes,
  file.path(
    log_dir,
    "mouse_deseq_QC_Hcn1_Kcnab2.csv"
  )
)

# -------------------------------------------------------------------------
# 10. Record software environment
# -------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    log_dir,
    "sessionInfo_03_standardize_mouse_deseq.txt"
  )
)

message(
  "Mouse DESeq2 standardization completed."
)

message(
  "Final mouse gene-level table: ",
  file.path(
    processed_dir,
    "mouse_c4oe_deseq_gene_level.csv"
  )
)

message(
  "Review summary: ",
  file.path(
    log_dir,
    "mouse_deseq_standardization_summary.csv"
  )
)