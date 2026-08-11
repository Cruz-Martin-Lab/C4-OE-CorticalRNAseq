# 02_standardize_aryal_s2.R
#
# Purpose:
# Standardize Aryal et al. Table S2 and create one human SCZ
# proteomics record per gene.
#
# Duplicate protein entries are resolved using measurement support,
# not statistical significance:
#   1. Highest numSpectraProteinObserved
#   2. Highest numPepsUnique
#   3. Highest percentCoverage
#   4. Alphabetical accession_number as a final deterministic tie-breaker

library(readxl)
library(readr)
library(dplyr)
library(stringr)
library(tibble)

# -------------------------------------------------------------------------
# 1. Define paths
# -------------------------------------------------------------------------

find_one_file <- function(patterns, files, label) {
  for (pattern in patterns) {
    matches <- files[stringr::str_detect(stringr::str_to_lower(basename(files)), stringr::str_to_lower(pattern))]
    if (length(matches) == 1) {
      return(matches)
    }
  }

  stop(
    "Could not uniquely find ",
    label,
    ". Tried patterns: ",
    paste(patterns, collapse = ", ")
  )
}

processed_dir <- "data_processed"
log_dir <- file.path("results", "logs")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists("data_raw")) {
  stop("Missing project directory: data_raw")
}

raw_files <- list.files("data_raw", full.names = TRUE)
input_file <- find_one_file(
  patterns = c(
    "table s2",
    "ms_ms analysis",
    "ms:ms analysis",
    "dlpfc synapse proteomes"
  ),
  files = raw_files,
  label = "Aryal Table S2 workbook"
)

message("Using Aryal Table S2 workbook: ", basename(input_file))

# -------------------------------------------------------------------------
# 2. Read the full Results sheet
# -------------------------------------------------------------------------

aryal_raw <- read_excel(
  input_file,
  sheet = "Results",
  skip = 2,
  .name_repair = "unique"
)

# The Results sheet contains eight sample-metadata rows beneath
# the column-header row. These are not protein observations.
aryal_metadata_rows <- aryal_raw %>%
  slice_head(n = 8)

aryal_protein_raw <- aryal_raw %>%
  dplyr::slice(-(1:8))

readr::write_csv(
  aryal_metadata_rows,
  file.path(
    log_dir,
    "aryal_s2_spreadsheet_metadata_rows.csv"
  )
)

if (nrow(aryal_protein_raw) != 8996) {
  warning(
    "Expected 8,996 protein entries after removing metadata rows, but found ",
    nrow(aryal_protein_raw)
  )
}

required_columns <- c(
  "accession_number",
  "geneSymbol",
  "numSpectraProteinObserved",
  "numPepsUnique",
  "percentCoverage",
  "entry_name",
  "Schizophrenia_vs_Normal.logFC",
  "Schizophrenia_vs_Normal.P.Value",
  "Schizophrenia_vs_Normal.adj.P.Val",
  "Schizophrenia_vs_Normal.t"
)

missing_columns <- setdiff(required_columns, names(aryal_raw))

if (length(missing_columns) > 0) {
  stop(
    "Required columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

# -------------------------------------------------------------------------
# 3. Create a clean protein-level table
# -------------------------------------------------------------------------

aryal_protein_level <- aryal_protein_raw %>%
  transmute(
    accession_number = as.character(accession_number),
    gene_symbol = str_trim(as.character(geneSymbol)),
    entry_name = as.character(entry_name),
    
    num_spectra = suppressWarnings(
      as.numeric(numSpectraProteinObserved)
    ),
    
    num_unique_peptides = suppressWarnings(
      as.numeric(numPepsUnique)
    ),
    
    percent_coverage = suppressWarnings(
      as.numeric(percentCoverage)
    ),
    
    scz_log2fc = suppressWarnings(
      as.numeric(Schizophrenia_vs_Normal.logFC)
    ),
    
    scz_t_stat = suppressWarnings(
      as.numeric(Schizophrenia_vs_Normal.t)
    ),
    
    scz_p_value = suppressWarnings(
      as.numeric(Schizophrenia_vs_Normal.P.Value)
    ),
    
    scz_fdr = suppressWarnings(
      as.numeric(Schizophrenia_vs_Normal.adj.P.Val)
    )
  ) %>%
  mutate(
    gene_symbol = na_if(gene_symbol, ""),
    gene_symbol = if_else(
      str_to_lower(gene_symbol) %in% c("na", "n/a", "nan"),
      NA_character_,
      gene_symbol
    ),
    
    scz_direction = case_when(
      is.na(scz_log2fc) ~ NA_character_,
      scz_log2fc > 0 ~ "up",
      scz_log2fc < 0 ~ "down",
      TRUE ~ "unchanged"
    ),
    
    scz_dep_fdr_0_10 = !is.na(scz_fdr) & scz_fdr < 0.10,
    scz_dep_fdr_0_05 = !is.na(scz_fdr) & scz_fdr < 0.05
  )

# Save all protein entries, including entries lacking gene symbols,
# so no information is silently discarded.
readr::write_csv(
  aryal_protein_level,
  file.path(
    processed_dir,
    "aryal_scz_protein_level_all_entries.csv"
  )
)

# -------------------------------------------------------------------------
# 4. Separate rows that can be used for gene-level analysis
# -------------------------------------------------------------------------

aryal_with_gene <- aryal_protein_level %>%
  filter(!is.na(gene_symbol))

aryal_without_gene <- aryal_protein_level %>%
  filter(is.na(gene_symbol))

readr::write_csv(
  aryal_without_gene,
  file.path(
    log_dir,
    "aryal_s2_entries_without_gene_symbol.csv"
  )
)

# -------------------------------------------------------------------------
# 5. Identify genes represented by multiple protein entries
# -------------------------------------------------------------------------

duplicate_gene_summary <- aryal_with_gene %>%
  count(gene_symbol, name = "n_protein_entries") %>%
  filter(n_protein_entries > 1) %>%
  arrange(desc(n_protein_entries), gene_symbol)

duplicate_gene_entries <- aryal_with_gene %>%
  semi_join(
    duplicate_gene_summary,
    by = "gene_symbol"
  ) %>%
  left_join(
    duplicate_gene_summary,
    by = "gene_symbol"
  ) %>%
  arrange(
    gene_symbol,
    desc(num_spectra),
    desc(num_unique_peptides),
    desc(percent_coverage),
    accession_number
  )

readr::write_csv(
  duplicate_gene_summary,
  file.path(
    log_dir,
    "aryal_s2_duplicate_gene_summary.csv"
  )
)

readr::write_csv(
  duplicate_gene_entries,
  file.path(
    log_dir,
    "aryal_s2_duplicate_gene_entries.csv"
  )
)

# -------------------------------------------------------------------------
# 6. Select one representative protein entry per gene
# -------------------------------------------------------------------------

aryal_gene_level <- aryal_with_gene %>%
  group_by(gene_symbol) %>%
  mutate(
    n_protein_entries = n()
  ) %>%
  arrange(
    desc(num_spectra),
    desc(num_unique_peptides),
    desc(percent_coverage),
    accession_number,
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(
    duplicate_resolution = if_else(
      n_protein_entries == 1,
      "single protein entry",
      paste0(
        "selected highest spectral support from ",
        n_protein_entries,
        " protein entries"
      )
    )
  ) %>%
  arrange(gene_symbol)

if (anyDuplicated(aryal_gene_level$gene_symbol) > 0) {
  stop("Gene-level output still contains duplicated gene symbols.")
}

# -------------------------------------------------------------------------
# 7. Save the final standardized gene-level table
# -------------------------------------------------------------------------

readr::write_csv(
  aryal_gene_level,
  file.path(
    processed_dir,
    "aryal_scz_gene_level.csv"
  )
)

# Audit only: which row was selected for genes with duplicate entries?
selected_duplicate_entries <- aryal_gene_level %>%
  filter(n_protein_entries > 1)

readr::write_csv(
  selected_duplicate_entries,
  file.path(
    log_dir,
    "aryal_s2_selected_duplicate_entries.csv"
  )
)

# -------------------------------------------------------------------------
# 8. Create a processing summary
# -------------------------------------------------------------------------

processing_summary <- tibble(
  metric = c(
    "Raw rows read beneath the header",
    "Spreadsheet metadata rows removed",
    "Protein entries after metadata removal",
    "Protein entries with a usable gene symbol",
    "Protein entries without a usable gene symbol",
    "Unique genes before duplicate collapse",
    "Genes with multiple protein entries",
    "Final gene-level rows",
    "Gene-level SCZ DEPs at FDR < 0.10",
    "Gene-level SCZ DEPs at FDR < 0.05",
    "Gene-level SCZ DEPs up at FDR < 0.10",
    "Gene-level SCZ DEPs down at FDR < 0.10"
  ),
  
  value = c(
    nrow(aryal_raw),
    nrow(aryal_metadata_rows),
    nrow(aryal_protein_raw),
    nrow(aryal_with_gene),
    nrow(aryal_without_gene),
    n_distinct(aryal_with_gene$gene_symbol),
    nrow(duplicate_gene_summary),
    nrow(aryal_gene_level),
    sum(aryal_gene_level$scz_dep_fdr_0_10, na.rm = TRUE),
    sum(aryal_gene_level$scz_dep_fdr_0_05, na.rm = TRUE),
    sum(
      aryal_gene_level$scz_dep_fdr_0_10 &
        aryal_gene_level$scz_direction == "up",
      na.rm = TRUE
    ),
    sum(
      aryal_gene_level$scz_dep_fdr_0_10 &
        aryal_gene_level$scz_direction == "down",
      na.rm = TRUE
    )
  )
)

readr::write_csv(
  processing_summary,
  file.path(
    log_dir,
    "aryal_s2_standardization_summary.csv"
  )
)

# -------------------------------------------------------------------------
# 9. Save selected QC genes for inspection only
# -------------------------------------------------------------------------

qc_gene_rows <- aryal_gene_level %>%
  filter(gene_symbol %in% c("HCN1", "KCNAB2"))

readr::write_csv(
  qc_gene_rows,
  file.path(
    log_dir,
    "aryal_s2_gene_level_QC_HCN1_KCNAB2.csv"
  )
)

# -------------------------------------------------------------------------
# 10. Record software environment
# -------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    log_dir,
    "sessionInfo_02_standardize_aryal_s2.txt"
  )
)

message("Aryal Table S2 standardization completed.")
message(
  "Final gene-level table: ",
  file.path(processed_dir, "aryal_scz_gene_level.csv")
)
message(
  "Review summary: ",
  file.path(log_dir, "aryal_s2_standardization_summary.csv")
)