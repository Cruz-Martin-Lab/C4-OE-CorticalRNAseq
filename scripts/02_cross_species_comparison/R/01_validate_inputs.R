# 01_validate_inputs.R
# Audit raw mouse RNA-seq and Aryal human proteomics inputs.
# This script reads but never modifies files in data_raw/.

required_packages <- c(
  "readr",
  "readxl",
  "dplyr",
  "purrr",
  "stringr",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

library(readr)
library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)

# -------------------------------------------------------------------------
# 1. Confirm project folders
# -------------------------------------------------------------------------

required_directories <- c(
  "data_raw",
  "data_processed",
  "results",
  "results/logs"
)

# Create analysis/output directories if they do not yet exist.
# The raw-input folder must exist, because this script reads from it.
if (!dir.exists("data_raw")) {
  stop("Missing project directory: data_raw")
}

for (relative_dir in setdiff(required_directories, "data_raw")) {
  dir.create(relative_dir, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------------------------------------------------------
# 2. Locate raw input files
# -------------------------------------------------------------------------

raw_files <- list.files(
  path = "data_raw",
  full.names = TRUE,
  recursive = FALSE
)

if (length(raw_files) == 0) {
  stop("No files were found in data_raw/.")
}

file_manifest <- tibble(
  file_path = raw_files,
  file_name = basename(raw_files),
  extension = tools::file_ext(raw_files),
  size_bytes = file.info(raw_files)$size
) |>
  arrange(file_name)

readr::write_csv(
  file_manifest,
  "results/logs/input_file_manifest.csv"
)

message("Files detected in data_raw/:")
print(file_manifest)

# -------------------------------------------------------------------------
# 3. Identify expected files by filename pattern
# -------------------------------------------------------------------------

find_one_file <- function(pattern, files, label) {
  matches <- files[
    str_detect(
      str_to_lower(basename(files)),
      str_to_lower(pattern)
    )
  ]
  
  if (length(matches) == 0) {
    stop(
      "Could not find ",
      label,
      " using pattern: ",
      pattern
    )
  }
  
  if (length(matches) > 1) {
    stop(
      "More than one file matched ",
      label,
      ": ",
      paste(basename(matches), collapse = "; ")
    )
  }
  
  matches
}

# Anchored and exact. Two reasons this pattern is fussy:
#
#   1. A pattern like "s1_table_deseq2" would match `S1_Table_DESeq2.csv` but
#      NOT the expected `S1_Table_DESeq.csv` -- so this script would stop dead.
#   2. A loose pattern like "s1_table_deseq" would match BOTH names if such a
#      file is still lying around in data_raw/, and find_one_file() then errors
#      with "More than one file matched".
#
# run_full_analysis.R puts exactly one mouse table here: strand 01's output.
# If you are running this strand by hand, make sure no stray copies remain in
# data_raw/.
mouse_file <- find_one_file(
  pattern = "^s1_table_deseq\\.csv$",
  files = raw_files,
  label = "mouse DESeq2 table (expected: S1_Table_DESeq.csv from strand 01)"
)

aryal_s2_file <- find_one_file(
  pattern = "table s2",
  files = raw_files,
  label = "Aryal Table S2"
)

aryal_s3_file <- find_one_file(
  pattern = "table s3",
  files = raw_files,
  label = "Aryal Table S3"
)

aryal_s4_file <- find_one_file(
  pattern = "table s4",
  files = raw_files,
  label = "Aryal Table S4"
)

# -------------------------------------------------------------------------
# 4. Audit mouse DESeq2 CSV
# -------------------------------------------------------------------------

mouse <- read_csv(
  mouse_file,
  show_col_types = FALSE,
  name_repair = "unique"
)

mouse_audit <- tibble(
  dataset = "Mouse C4-OE DESeq2",
  rows = nrow(mouse),
  columns = ncol(mouse),
  duplicate_column_names = sum(duplicated(names(mouse))),
  completely_empty_rows = sum(
    apply(mouse, 1, function(x) all(is.na(x)))
  )
)

mouse_columns <- tibble(
  dataset = "Mouse C4-OE DESeq2",
  column_number = seq_along(names(mouse)),
  column_name = names(mouse),
  class = map_chr(
    mouse,
    ~ paste(class(.x), collapse = "/")
  ),
  missing_values = map_int(
    mouse,
    ~ sum(is.na(.x))
  ),
  unique_values = map_int(
    mouse,
    ~ n_distinct(.x, na.rm = TRUE)
  )
)

readr::write_csv(
  mouse_audit,
  "results/logs/mouse_dataset_summary.csv"
)

readr::write_csv(
  mouse_columns,
  "results/logs/mouse_column_audit.csv"
)

# Find likely mouse gene, effect-size, and adjusted-P columns.

mouse_gene_candidates <- names(mouse)[
  str_detect(
    str_to_lower(names(mouse)),
    "^(gene|symbol|gene_symbol|genesymbol)$"
  )
]

mouse_logfc_candidates <- names(mouse)[
  str_detect(
    str_to_lower(names(mouse)),
    "log2foldchange|log2fc|logfc"
  )
]

mouse_padj_candidates <- names(mouse)[
  str_detect(
    str_to_lower(names(mouse)),
    "padj|adj.*p|fdr"
  )
]

mouse_key_columns <- tibble(
  field = c(
    "gene",
    "effect_size",
    "adjusted_p"
  ),
  candidate_columns = c(
    paste(mouse_gene_candidates, collapse = "; "),
    paste(mouse_logfc_candidates, collapse = "; "),
    paste(mouse_padj_candidates, collapse = "; ")
  )
)

readr::write_csv(
  mouse_key_columns,
  "results/logs/mouse_key_column_candidates.csv"
)

# Extract Hcn1 and Kcnab2 rows.

if (length(mouse_gene_candidates) == 1) {
  mouse_gene_column <- mouse_gene_candidates[[1]]
  
  mouse_target_rows <- mouse |>
    filter(
      str_to_upper(
        as.character(.data[[mouse_gene_column]])
      ) %in% c("HCN1", "KCNAB2")
    )
  
  readr::write_csv(
    mouse_target_rows,
    "results/logs/mouse_Hcn1_Kcnab2_rows.csv"
  )
} else {
  warning(
    "Could not uniquely identify the mouse gene-symbol column. ",
    "Candidates: ",
    paste(mouse_gene_candidates, collapse = ", ")
  )
}

# -------------------------------------------------------------------------
# 5. Audit Aryal Excel workbooks
# -------------------------------------------------------------------------

audit_workbook <- function(file_path, dataset_label) {
  sheets <- excel_sheets(file_path)
  
  map_dfr(
    sheets,
    function(sheet_name) {
      dat <- read_excel(
        file_path,
        sheet = sheet_name,
        n_max = 25,
        .name_repair = "unique"
      )
      
      tibble(
        dataset = dataset_label,
        file_name = basename(file_path),
        sheet = sheet_name,
        preview_rows_read = nrow(dat),
        columns = ncol(dat),
        column_names = paste(
          names(dat),
          collapse = " | "
        )
      )
    }
  )
}

aryal_s2_sheets <- audit_workbook(
  aryal_s2_file,
  "Aryal Table S2"
)

aryal_s3_sheets <- audit_workbook(
  aryal_s3_file,
  "Aryal Table S3"
)

aryal_s4_sheets <- audit_workbook(
  aryal_s4_file,
  "Aryal Table S4"
)

aryal_workbook_audit <- bind_rows(
  aryal_s2_sheets,
  aryal_s3_sheets,
  aryal_s4_sheets
)

readr::write_csv(
  aryal_workbook_audit,
  "results/logs/aryal_workbook_sheet_audit.csv"
)

# -------------------------------------------------------------------------
# 6. Read Aryal Table S2 Results sheet
# -------------------------------------------------------------------------

s2_sheet_names <- excel_sheets(aryal_s2_file)

results_sheet <- s2_sheet_names[
  str_to_lower(
    str_trim(s2_sheet_names)
  ) == "results"
]

if (length(results_sheet) != 1) {
  stop(
    "Could not uniquely identify the 'Results' sheet in Aryal Table S2. ",
    "Available sheets: ",
    paste(s2_sheet_names, collapse = ", ")
  )
}

aryal_s2 <- read_excel(
  aryal_s2_file,
  sheet = results_sheet[[1]],
  skip = 2,
  .name_repair = "unique"
)

aryal_s2_summary <- tibble(
  dataset = "Aryal Table S2 Results",
  rows = nrow(aryal_s2),
  columns = ncol(aryal_s2),
  duplicate_column_names = sum(
    duplicated(names(aryal_s2))
  )
)

aryal_s2_columns <- tibble(
  dataset = "Aryal Table S2 Results",
  column_number = seq_along(names(aryal_s2)),
  column_name = names(aryal_s2),
  class = map_chr(
    aryal_s2,
    ~ paste(class(.x), collapse = "/")
  ),
  missing_values = map_int(
    aryal_s2,
    ~ sum(is.na(.x))
  ),
  unique_values = map_int(
    aryal_s2,
    ~ n_distinct(.x, na.rm = TRUE)
  )
)

readr::write_csv(
  aryal_s2_summary,
  "results/logs/aryal_s2_dataset_summary.csv"
)

readr::write_csv(
  aryal_s2_columns,
  "results/logs/aryal_s2_column_audit.csv"
)

# Identify likely human gene-symbol and SCZ-comparison columns.

human_gene_candidates <- names(aryal_s2)[
  str_detect(
    str_to_lower(names(aryal_s2)),
    "gene.*symbol|genesymbol"
  )
]

scz_logfc_candidates <- names(aryal_s2)[
  str_detect(
    str_to_lower(names(aryal_s2)),
    "schizophrenia.*normal.*logfc|scz.*ctrl.*logfc"
  )
]

scz_padj_candidates <- names(aryal_s2)[
  str_detect(
    str_to_lower(names(aryal_s2)),
    paste0(
      "schizophrenia.*normal.*adj|",
      "scz.*ctrl.*adj|",
      "schizophrenia.*normal.*fdr"
    )
  )
]

aryal_key_columns <- tibble(
  field = c(
    "human_gene",
    "SCZ_effect_size",
    "SCZ_adjusted_p"
  ),
  candidate_columns = c(
    paste(human_gene_candidates, collapse = "; "),
    paste(scz_logfc_candidates, collapse = "; "),
    paste(scz_padj_candidates, collapse = "; ")
  )
)

readr::write_csv(
  aryal_key_columns,
  "results/logs/aryal_s2_key_column_candidates.csv"
)

# Extract all HCN1 and KCNAB2 protein/isoform rows.

if (length(human_gene_candidates) == 1) {
  human_gene_column <- human_gene_candidates[[1]]
  
  aryal_target_rows <- aryal_s2 |>
    filter(
      str_to_upper(
        as.character(.data[[human_gene_column]])
      ) %in% c("HCN1", "KCNAB2")
    )
  
  readr::write_csv(
    aryal_target_rows,
    "results/logs/aryal_HCN1_KCNAB2_rows.csv"
  )
} else {
  warning(
    "Could not uniquely identify the Aryal gene-symbol column. ",
    "Candidates: ",
    paste(human_gene_candidates, collapse = ", ")
  )
}

# -------------------------------------------------------------------------
# 7. Save session information
# -------------------------------------------------------------------------

session_text <- capture.output(
  sessionInfo()
)

writeLines(
  session_text,
  "results/logs/sessionInfo_01_validate_inputs.txt"
)

message("")
message("Input audit completed.")
message("Review the files written to results/logs/.")