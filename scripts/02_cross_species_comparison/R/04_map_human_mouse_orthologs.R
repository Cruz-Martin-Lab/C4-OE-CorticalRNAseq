# 04_map_human_mouse_orthologs.R
#
# Purpose:
# Retrieve and audit Ensembl human-mouse ortholog mappings for the
# standardized Aryal human proteomics and mouse C4-OE RNA-seq datasets.
#
# Primary cross-species analyses will use Ensembl high-confidence
# one-to-one orthologs. Ambiguous and missing mappings are retained
# in audit outputs rather than silently discarded.

library(biomaRt)
library(readr)
library(dplyr)
library(stringr)
library(tibble)

# -------------------------------------------------------------------------
# 1. Define paths
# -------------------------------------------------------------------------

human_file <- file.path(
  "data_processed",
  "aryal_scz_gene_level.csv"
)

mouse_file <- file.path(
  "data_processed",
  "mouse_c4oe_deseq_gene_level.csv"
)

processed_dir <- "data_processed"
log_dir <- file.path("results", "logs")

dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(human_file)) {
  stop("Human standardized file not found: ", human_file)
}

if (!file.exists(mouse_file)) {
  stop("Mouse standardized file not found: ", mouse_file)
}

# -------------------------------------------------------------------------
# 2. Read standardized human and mouse datasets
# -------------------------------------------------------------------------

human_data <- readr::read_csv(
  human_file,
  show_col_types = FALSE
)

mouse_data <- readr::read_csv(
  mouse_file,
  show_col_types = FALSE
)

required_human_columns <- c(
  "gene_symbol",
  "scz_log2fc",
  "scz_fdr",
  "scz_dep_fdr_0_10",
  "scz_dep_fdr_0_05"
)

required_mouse_columns <- c(
  "mouse_gene_symbol",
  "mouse_log2fc",
  "mouse_padj",
  "mouse_deg_padj_0_05",
  "mouse_deg_padj_0_10"
)

missing_human_columns <- setdiff(
  required_human_columns,
  names(human_data)
)

missing_mouse_columns <- setdiff(
  required_mouse_columns,
  names(mouse_data)
)

if (length(missing_human_columns) > 0) {
  stop(
    "Missing required human columns: ",
    paste(missing_human_columns, collapse = ", ")
  )
}

if (length(missing_mouse_columns) > 0) {
  stop(
    "Missing required mouse columns: ",
    paste(missing_mouse_columns, collapse = ", ")
  )
}

human_symbols <- human_data %>%
  filter(!is.na(gene_symbol)) %>%
  distinct(gene_symbol) %>%
  pull(gene_symbol)

mouse_symbols <- mouse_data %>%
  filter(!is.na(mouse_gene_symbol)) %>%
  distinct(mouse_gene_symbol) %>%
  pull(mouse_gene_symbol)

# -------------------------------------------------------------------------
# 3. Obtain the human-mouse ortholog table
# -------------------------------------------------------------------------
#
# REPRODUCIBILITY. This is the one step in the whole project that depended on
# a live external service. `useEnsembl()` without a `version` argument returns
# whatever release Ensembl is currently serving, so the ortholog set -- and
# therefore the shared measurable universe, every Fisher denominator, and every
# GSEA gene-set size -- silently changed over time. Another group running this
# next year would not have reproduced the published numbers.
#
# Two changes fix that:
#
#   1. The Ensembl release is PINNED via ENSEMBL_RELEASE below. Archived
#      releases stay available at Ensembl indefinitely, so this is stable.
#   2. The query result is FROZEN to a CSV that is committed to the repo.
#      When that file is present it is used directly and Ensembl is never
#      contacted, so the pipeline is fully offline and deterministic.
#
# To deliberately refresh the frozen table, delete it and re-run; the new file
# is written automatically and should then be committed with a note saying why.

# Ensembl release the analysis is pinned to. Pinning fixes the human-mouse
# ortholog set (and therefore the shared gene universe) so it does not depend on
# whatever release Ensembl serves on a given day. Release 116 (Jun 2026) is
# served by BioMart at www.ensembl.org; the query result is frozen to CSV and
# committed, so every later run is offline and deterministic regardless of which
# release is reachable. The methods text quotes release 116.
ENSEMBL_RELEASE <- 116

frozen_ortholog_file <- file.path(
  "reference",
  "ensembl_human_mouse_orthologs_frozen.csv"
)

connect_to_ensembl <- function() {

  mirrors <- c(
    "useast",
    "www",
    "asia"
  )

  last_error <- NULL

  for (selected_mirror in mirrors) {

    result <- tryCatch(
      {
        biomaRt::useEnsembl(
          biomart = "genes",
          dataset = "hsapiens_gene_ensembl",
          version = ENSEMBL_RELEASE,
          mirror  = selected_mirror
        )
      },
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(result)) {
      message(
        "Connected to Ensembl release ", ENSEMBL_RELEASE,
        " via mirror: ", selected_mirror
      )
      return(result)
    }
  }

  stop(
    "Could not connect to an Ensembl BioMart mirror for release ",
    ENSEMBL_RELEASE, ". Last error: ", last_error,
    "\nIf you are offline, place the frozen ortholog table at: ",
    frozen_ortholog_file
  )
}

USING_FROZEN_ORTHOLOGS <- file.exists(frozen_ortholog_file)

if (USING_FROZEN_ORTHOLOGS) {
  message(
    "Using the FROZEN Ensembl ortholog table (Ensembl is not contacted): ",
    frozen_ortholog_file
  )
  human_mart <- NULL
} else {
  message(
    "No frozen ortholog table found at ", frozen_ortholog_file,
    " -- querying Ensembl release ", ENSEMBL_RELEASE,
    " and freezing the result for future runs."
  )
  human_mart <- connect_to_ensembl()
}

# -------------------------------------------------------------------------
# 4. Verify the required BioMart attributes
# -------------------------------------------------------------------------

requested_attributes <- c(
  "ensembl_gene_id",
  "external_gene_name",
  "mmusculus_homolog_ensembl_gene",
  "mmusculus_homolog_associated_gene_name",
  "mmusculus_homolog_orthology_type",
  "mmusculus_homolog_orthology_confidence"
)

if (!USING_FROZEN_ORTHOLOGS) {

  available_attributes <- biomaRt::listAttributes(
    human_mart
  )

  missing_attributes <- setdiff(
    requested_attributes,
    available_attributes$name
  )

  if (length(missing_attributes) > 0) {

    readr::write_csv(
      available_attributes,
      file.path(
        log_dir,
        "ensembl_available_attributes.csv"
      )
    )

    stop(
      "Required Ensembl attributes are unavailable in Ensembl release ",
      ENSEMBL_RELEASE, ": ",
      paste(missing_attributes, collapse = ", "),
      ". Available attributes were saved to ",
      file.path(
        log_dir,
        "ensembl_available_attributes.csv"
      )
    )
  }
}

# -------------------------------------------------------------------------
# 5. Query Ensembl in batches
# -------------------------------------------------------------------------
#
# Batching avoids sending all 8,444 human symbols in one request.

query_ortholog_batch <- function(symbol_batch) {

  biomaRt::getBM(
    attributes = requested_attributes,
    filters = "external_gene_name",
    values = symbol_batch,
    mart = human_mart,
    uniqueRows = TRUE
  )
}

batch_size <- 500

if (USING_FROZEN_ORTHOLOGS) {

  # Offline path: read the committed table and skip Ensembl entirely.
  ortholog_raw <- readr::read_csv(
    frozen_ortholog_file,
    show_col_types = FALSE
  )

  missing_frozen_columns <- setdiff(
    requested_attributes,
    names(ortholog_raw)
  )

  if (length(missing_frozen_columns) > 0) {
    stop(
      "The frozen ortholog table is missing required columns: ",
      paste(missing_frozen_columns, collapse = ", "),
      ". Delete ", frozen_ortholog_file,
      " and re-run to regenerate it from Ensembl release ", ENSEMBL_RELEASE, "."
    )
  }

  message(
    "Loaded ", nrow(ortholog_raw),
    " frozen ortholog rows; Ensembl was not contacted."
  )

  # The frozen table was built for a particular set of human symbols. If the
  # upstream human data changes, symbols absent from the frozen file would
  # silently fall through as "not returned by Ensembl" and quietly shrink the
  # shared universe. Warn loudly instead.
  unqueried <- setdiff(human_symbols, ortholog_raw$external_gene_name)

  if (length(unqueried) > 0) {
    warning(
      length(unqueried), " of ", length(human_symbols),
      " human symbols are absent from the frozen ortholog table, so they ",
      "cannot be mapped. This means the frozen table predates a change in ",
      "the human input. Delete ", frozen_ortholog_file,
      " and re-run to rebuild it from Ensembl release ", ENSEMBL_RELEASE, ".",
      call. = FALSE
    )

    readr::write_csv(
      tibble(human_gene_symbol = unqueried),
      file.path(log_dir, "symbols_absent_from_frozen_ortholog_table.csv")
    )
  }

} else {

human_symbol_batches <- split(
  human_symbols,
  ceiling(
    seq_along(human_symbols) / batch_size
  )
)

ortholog_batch_results <- vector(
  mode = "list",
  length = length(human_symbol_batches)
)

for (batch_index in seq_along(human_symbol_batches)) {
  
  message(
    "Querying Ensembl batch ",
    batch_index,
    " of ",
    length(human_symbol_batches)
  )
  
  ortholog_batch_results[[batch_index]] <- tryCatch(
    {
      query_ortholog_batch(
        human_symbol_batches[[batch_index]]
      )
    },
    error = function(e) {
      stop(
        "Ensembl query failed for batch ",
        batch_index,
        ": ",
        conditionMessage(e)
      )
    }
  )
}

  ortholog_raw <- bind_rows(
    ortholog_batch_results
  )

  # Freeze the query result so every later run is offline and identical.
  dir.create(dirname(frozen_ortholog_file), recursive = TRUE,
             showWarnings = FALSE)

  readr::write_csv(
    ortholog_raw,
    frozen_ortholog_file
  )

  message(
    "Froze ", nrow(ortholog_raw), " ortholog rows to ", frozen_ortholog_file,
    " -- COMMIT THIS FILE so others reproduce the published universe."
  )
}

# -------------------------------------------------------------------------
# 6. Standardize the Ensembl mapping table
# -------------------------------------------------------------------------

ortholog_standardized <- ortholog_raw %>%
  transmute(
    human_ensembl_gene_id =
      as.character(ensembl_gene_id),
    
    human_gene_symbol =
      str_trim(as.character(external_gene_name)),
    
    mouse_ensembl_gene_id =
      as.character(
        mmusculus_homolog_ensembl_gene
      ),
    
    mouse_gene_symbol =
      str_trim(
        as.character(
          mmusculus_homolog_associated_gene_name
        )
      ),
    
    ensembl_orthology_type =
      as.character(
        mmusculus_homolog_orthology_type
      ),
    
    ensembl_orthology_confidence =
      suppressWarnings(
        as.numeric(
          mmusculus_homolog_orthology_confidence
        )
      ),
    
  ) %>%
  mutate(
    human_gene_symbol =
      na_if(human_gene_symbol, ""),
    
    mouse_gene_symbol =
      na_if(mouse_gene_symbol, ""),
    
    human_ensembl_gene_id =
      na_if(human_ensembl_gene_id, ""),
    
    mouse_ensembl_gene_id =
      na_if(mouse_ensembl_gene_id, ""),
    
    mapping_has_mouse_ortholog =
      !is.na(mouse_gene_symbol) &
      !is.na(mouse_ensembl_gene_id),
    
    mapping_is_one_to_one =
      ensembl_orthology_type ==
      "ortholog_one2one",
    
    mapping_is_high_confidence =
      ensembl_orthology_confidence == 1,
    
    mapping_is_primary =
      mapping_has_mouse_ortholog &
      mapping_is_one_to_one &
      mapping_is_high_confidence
  ) %>%
  distinct()

readr::write_csv(
  ortholog_standardized,
  file.path(
    processed_dir,
    "ensembl_human_mouse_orthologs_raw.csv"
  )
)

# -------------------------------------------------------------------------
# 7. Add human genes missing from the Ensembl return
# -------------------------------------------------------------------------

all_human_mapping_audit <- tibble(
  human_gene_symbol = human_symbols
) %>%
  left_join(
    ortholog_standardized,
    by = "human_gene_symbol"
  ) %>%
  mutate(
    mapping_returned_by_ensembl =
      !is.na(human_ensembl_gene_id),
    
    mapping_status = case_when(
      !mapping_returned_by_ensembl ~
        "not returned by Ensembl",
      
      is.na(mouse_gene_symbol) ~
        "human gene returned but no mouse ortholog",
      
      mapping_is_primary ~
        "high-confidence one-to-one",
      
      mapping_is_one_to_one &
        !mapping_is_high_confidence ~
        "one-to-one, not high confidence",
      
      ensembl_orthology_type ==
        "ortholog_one2many" ~
        "one-to-many",
      
      ensembl_orthology_type ==
        "ortholog_many2many" ~
        "many-to-many",
      
      TRUE ~
        "other or unclassified"
    )
  )

readr::write_csv(
  all_human_mapping_audit,
  file.path(
    log_dir,
    "human_mouse_ortholog_mapping_audit.csv"
  )
)

# -------------------------------------------------------------------------
# 8. Identify duplicate or ambiguous mappings
# -------------------------------------------------------------------------

human_mapping_counts <- ortholog_standardized %>%
  filter(mapping_has_mouse_ortholog) %>%
  count(
    human_gene_symbol,
    name = "n_mouse_mappings"
  )

mouse_mapping_counts <- ortholog_standardized %>%
  filter(mapping_has_mouse_ortholog) %>%
  count(
    mouse_gene_symbol,
    name = "n_human_mappings"
  )

ambiguous_mapping_rows <- ortholog_standardized %>%
  left_join(
    human_mapping_counts,
    by = "human_gene_symbol"
  ) %>%
  left_join(
    mouse_mapping_counts,
    by = "mouse_gene_symbol"
  ) %>%
  filter(
    n_mouse_mappings > 1 |
      n_human_mappings > 1 |
      !mapping_is_primary
  ) %>%
  arrange(
    human_gene_symbol,
    mouse_gene_symbol
  )

readr::write_csv(
  human_mapping_counts,
  file.path(
    log_dir,
    "ortholog_human_mapping_counts.csv"
  )
)

readr::write_csv(
  mouse_mapping_counts,
  file.path(
    log_dir,
    "ortholog_mouse_mapping_counts.csv"
  )
)

readr::write_csv(
  ambiguous_mapping_rows,
  file.path(
    log_dir,
    "ortholog_ambiguous_mapping_rows.csv"
  )
)

# -------------------------------------------------------------------------
# 9. Create the primary high-confidence one-to-one mapping
# -------------------------------------------------------------------------

primary_one_to_one <- ortholog_standardized %>%
  filter(mapping_is_primary) %>%
  add_count(
    human_gene_symbol,
    name = "n_mouse_per_human"
  ) %>%
  add_count(
    mouse_gene_symbol,
    name = "n_human_per_mouse"
  ) %>%
  filter(
    n_mouse_per_human == 1,
    n_human_per_mouse == 1
  ) %>%
  dplyr::select(
    human_ensembl_gene_id,
    human_gene_symbol,
    mouse_ensembl_gene_id,
    mouse_gene_symbol,
    ensembl_orthology_type,
    ensembl_orthology_confidence
  ) %>%
  distinct() %>%
  arrange(human_gene_symbol)

if (
  anyDuplicated(
    primary_one_to_one$human_gene_symbol
  ) > 0
) {
  stop(
    "Primary mapping still contains duplicated human symbols."
  )
}

if (
  anyDuplicated(
    primary_one_to_one$mouse_gene_symbol
  ) > 0
) {
  stop(
    "Primary mapping still contains duplicated mouse symbols."
  )
}

readr::write_csv(
  primary_one_to_one,
  file.path(
    processed_dir,
    "ensembl_human_mouse_one_to_one.csv"
  )
)

# -------------------------------------------------------------------------
# 10. Create the shared measurable universe
# -------------------------------------------------------------------------
#
# A pair enters the shared universe only when:
#   1. it is a high-confidence one-to-one Ensembl ortholog;
#   2. the human gene was measured in Aryal et al.;
#   3. the mouse gene was tested in the C4-OE DESeq2 analysis.

shared_universe <- primary_one_to_one %>%
  inner_join(
    human_data,
    by = c(
      "human_gene_symbol" = "gene_symbol"
    )
  ) %>%
  inner_join(
    mouse_data,
    by = "mouse_gene_symbol"
  ) %>%
  arrange(human_gene_symbol)

readr::write_csv(
  shared_universe,
  file.path(
    processed_dir,
    "human_mouse_shared_measurable_universe.csv"
  )
)

# -------------------------------------------------------------------------
# 11. Create summary
# -------------------------------------------------------------------------

mapping_summary <- tibble(
  metric = c(
    "Human genes submitted to Ensembl",
    "Human genes returned by Ensembl",
    "Human genes not returned by Ensembl",
    "Human genes with any mouse ortholog",
    "High-confidence one-to-one mapping rows",
    "Primary unique one-to-one ortholog pairs",
    "Mouse genes in the DESeq2 gene-level table",
    "Ortholog pairs in the shared measurable universe",
    "Human SCZ DEPs in the shared universe at FDR < 0.10",
    "Mouse C4-OE DEGs in the shared universe at padj < 0.05"
  ),
  
  value = c(
    length(human_symbols),
    
    n_distinct(
      ortholog_standardized$human_gene_symbol,
      na.rm = TRUE
    ),
    
    sum(
      !all_human_mapping_audit$
        mapping_returned_by_ensembl
    ),
    
    n_distinct(
      ortholog_standardized$
        human_gene_symbol[
          ortholog_standardized$
            mapping_has_mouse_ortholog
        ],
      na.rm = TRUE
    ),
    
    sum(
      ortholog_standardized$
        mapping_is_primary,
      na.rm = TRUE
    ),
    
    nrow(primary_one_to_one),
    
    length(mouse_symbols),
    
    nrow(shared_universe),
    
    sum(
      shared_universe$
        scz_dep_fdr_0_10,
      na.rm = TRUE
    ),
    
    sum(
      shared_universe$
        mouse_deg_padj_0_05,
      na.rm = TRUE
    )
  )
)

readr::write_csv(
  mapping_summary,
  file.path(
    log_dir,
    "human_mouse_ortholog_mapping_summary.csv"
  )
)

# -------------------------------------------------------------------------
# 12. Record software environment
# -------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    log_dir,
    "sessionInfo_04_map_human_mouse_orthologs.txt"
  )
)

message(
  "Human-mouse ortholog mapping completed."
)

message(
  "Primary one-to-one mapping: ",
  file.path(
    processed_dir,
    "ensembl_human_mouse_one_to_one.csv"
  )
)

message(
  "Shared measurable universe: ",
  file.path(
    processed_dir,
    "human_mouse_shared_measurable_universe.csv"
  )
)

message(
  "Review summary: ",
  file.path(
    log_dir,
    "human_mouse_ortholog_mapping_summary.csv"
  )
)