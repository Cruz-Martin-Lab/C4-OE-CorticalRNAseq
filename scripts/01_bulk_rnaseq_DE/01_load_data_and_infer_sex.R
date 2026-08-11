# =============================================================================
# 01_load_data_and_infer_sex.R
# -----------------------------------------------------------------------------
# - Loads raw_counts.csv
# - Maps Ensembl gene IDs to MGI symbols (org.Mm.eg.db, offline)
# - Filters zero-variance / all-but-one-zero genes
# - Infers sex per sample from Xist vs. Ddx3y and cross-checks it against
#   SAMPLE_METADATA in 00_config.R (which was itself derived from this
#   inference -- this script reproduces and documents that derivation)
# - Saves cleaned, symbol-mapped counts + metadata + QC plots for downstream
#   scripts
# =============================================================================

# ---- Locate this script's own folder, independent of the working directory
# (works via Rscript, RStudio's "Source" button, or source() from anywhere) --
.get_script_dir <- function() {
  # Check the source() frame FIRST. commandArgs() reports the Rscript ENTRY
  # POINT, which is the wrong answer whenever this script is sourced by a
  # master runner (run_all.R, run_full_analysis.R) living in another folder.
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  # gsub("~+~"): R's front end encodes spaces in --file= as "~+~". This
  # repository's path contains spaces, so without decoding, every derived path
  # is wrong and dir.create(recursive=TRUE) silently builds a bogus tree.
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1) return(dirname(normalizePath(
    gsub("~+~", " ", sub("^--file=", "", file_arg), fixed = TRUE))))
  # Fallback: this happens when code is run line-by-line (e.g. RStudio's
  # "Run" rather than "Source") instead of via source()/Rscript, so there's
  # no frame to introspect. Fall back to checking known relative locations
  # from the current working directory instead.
  candidates <- c(
    getwd(),                                          # wd is already this folder
    file.path(getwd(), "scripts", "01_bulk_rnaseq_DE") # wd is the repo root
  )
  for (cand in candidates) {
    if (file.exists(file.path(cand, "00_config.R"))) return(normalizePath(cand))
  }
  stop(
    "Could not determine this script's folder automatically (checked Rscript, ",
    "source(), and the working directory). This usually happens when script ",
    "code is run line-by-line instead of sourced as a whole file.\n",
    "Fix: source() the whole script (or Rscript it) rather than running lines ",
    "individually, e.g.:\n",
    "  source(\"scripts/01_bulk_rnaseq_DE/01_load_data_and_infer_sex.R\")\n",
    "Current working directory: ", getwd()
  )
}
SCRIPT_DIR <- .get_script_dir()
source(file.path(SCRIPT_DIR, "00_config.R"))
source(file.path(SCRIPT_DIR, "R", "functions.R"))

DIR_TABLES  <- step_dir("01", "tables")
DIR_FIGURES <- step_dir("01", "figures")

# ---- 1. Load raw counts ------------------------------------------------------
raw_counts <- read_counts(FILE_RAW_COUNTS)
message(sprintf("Loaded raw_counts.csv: %d genes x %d samples",
                nrow(raw_counts), ncol(raw_counts) - 1))

# ---- 2. Sex inference (BEFORE any filtering, so marker genes are guaranteed
#         to still be present regardless of downstream filters) --------------
sex_markers <- c("Xist" = "ENSMUSG00000086503", "Ddx3y" = "ENSMUSG00000069045")
ens_base <- gsub("\\.\\d+$", "", raw_counts$Gene)

marker_counts <- purrr::map_dfr(names(sex_markers), function(g) {
  row <- raw_counts[ens_base == sex_markers[[g]], -1, drop = FALSE]
  if (nrow(row) == 0) stop("Sex marker gene not found in raw counts: ", g)
  tibble(gene = g, sample = colnames(row), count = as.numeric(row[1, ]))
})

marker_wide <- marker_counts %>% tidyr::pivot_wider(names_from = gene, values_from = count)
marker_wide <- marker_wide %>%
  mutate(sex_inferred = dplyr::case_when(
    Xist > 1000 & Ddx3y < 100 ~ "F",
    Ddx3y > 1000 & Xist < 100 ~ "M",
    TRUE ~ "AMBIGUOUS"
  ))

if (any(marker_wide$sex_inferred == "AMBIGUOUS")) {
  warning("Ambiguous sex call for sample(s): ",
          paste(marker_wide$sample[marker_wide$sex_inferred == "AMBIGUOUS"], collapse = ", "),
          " -- inspect manually before trusting SAMPLE_METADATA.")
}

metadata_check <- SAMPLE_METADATA %>%
  dplyr::left_join(marker_wide, by = "sample")

mismatches <- metadata_check %>% dplyr::filter(sex != sex_inferred)
if (nrow(mismatches) > 0) {
  warning("Sex mismatch between SAMPLE_METADATA (00_config.R) and Xist/Ddx3y-",
          "inferred sex for: ", paste(mismatches$sample, collapse = ", "),
          " -- update 00_config.R before proceeding.")
} else {
  message("Sex inference matches SAMPLE_METADATA for all samples. Good.")
}

# Written with display labels (Control / C4-OE, Control_10 / C4-OE_4) so the
# supplementary table matches the figures.
metadata_check_out <- metadata_check
metadata_check_out$genotype <- display_genotype(metadata_check_out$genotype)
metadata_check_out$sample   <- display_sample(metadata_check_out$sample)
write.csv(metadata_check_out, file.path(DIR_TABLES, "sex_inference_check.csv"), row.names = FALSE)

sex_marker_plot <- plot_sex_markers(metadata_check)
ggplot2::ggsave(file.path(DIR_FIGURES, "sex_marker_diagnostic.pdf"),
                 sex_marker_plot, width = 6, height = 5)

# ---- 3. Build metadata for DESeq2 -------------------------------------------
metadata <- build_metadata(raw_counts, SAMPLE_METADATA)
metadata_out <- metadata
metadata_out$genotype <- display_genotype(metadata_out$genotype)
metadata_out$sample   <- display_sample(metadata_out$sample)
write.csv(metadata_out[, c("sample", "genotype", "sex")],
          file.path(DIR_TABLES, "sample_metadata.csv"), row.names = FALSE)

# ---- 4. Filter zero-variance genes, then map to gene symbols -----------------
filtered_counts <- zero_var_genes(raw_counts)
message(sprintf("After zero-variance filter: %d genes retained", nrow(filtered_counts)))

symbol_counts <- map_ensembl_to_symbol(filtered_counts)
message(sprintf("After Ensembl -> symbol mapping (org.Mm.eg.db): %d genes retained",
                nrow(symbol_counts)))
write.csv(symbol_counts, file.path(DIR_TABLES, "counts_symbol_mapped_filtered.csv"),
          row.names = FALSE)

# ---- 5. Exploratory QC plots -------------------------------------------------
lib_data <- get_library_size(symbol_counts)
ggplot2::ggsave(file.path(DIR_FIGURES, "library_size_barplot.pdf"),
                 lib_size_bar(lib_data), width = 6, height = 4)

var_mean_plot <- plot_variance_vs_mean(symbol_counts)
ggplot2::ggsave(file.path(DIR_FIGURES, "variance_vs_mean.pdf"),
                 var_mean_plot, width = 6, height = 4)

# ---- 6. Cache for downstream scripts -----------------------------------------
saveRDS(list(symbol_counts = symbol_counts, metadata = metadata),
        file.path(DIR_RDATA, "01_loaded_data.rds"))

message("01_load_data_and_infer_sex.R complete.")
