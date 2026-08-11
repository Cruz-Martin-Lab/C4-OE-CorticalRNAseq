# =============================================================================
# scripts/00_setup/install_dependencies.R
# -----------------------------------------------------------------------------
# Installs every R package the three analysis strands need.
#
# PREFERRED ROUTE: if renv.lock exists at the repo root, use that instead --
# it pins exact versions, which is what makes results reproducible:
#
#   renv::restore()
#
# Use this script only to bootstrap an environment from scratch (for example
# when first creating the lockfile), or if you are deliberately not using renv.
#
#   Rscript scripts/00_setup/install_dependencies.R
# =============================================================================

# Under `Rscript`, R's default repos is "@CRAN@" and install.packages() aborts
# with "trying to use CRAN without setting a mirror". Set one explicitly.
if (is.null(getOption("repos")) || any(getOption("repos") == "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

CRAN_PACKAGES <- c(
  # --- tidyverse and general tooling ---
  "dplyr", "readr", "tibble", "tidyr", "stringr", "purrr", "ggplot2",
  "readxl", "data.table", "rstudioapi",
  # --- plotting ---
  "patchwork", "cowplot", "ggrepel", "gridExtra", "scales",
  "pheatmap", "RColorBrewer", "gplots",
  # --- analysis ---
  "WGCNA", "dynamicTreeCut", "msigdbr",
  # --- single cell (strand 03) ---
  "Seurat", "SeuratObject"
)

BIOC_PACKAGES <- c(
  "DESeq2",
  "clusterProfiler",
  "org.Mm.eg.db",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "GO.db",
  "GOSemSim",
  "fgsea",
  "EnhancedVolcano",
  "biomaRt",
  # WGCNA's Bioconductor dependencies -- installing them explicitly avoids the
  # confusing "package 'impute' is not available" failure on a fresh machine.
  "impute",
  "preprocessCore"
)

message("R version: ", R.version.string)

# -----------------------------------------------------------------------------
# CRAN
# -----------------------------------------------------------------------------

missing_cran <- CRAN_PACKAGES[
  !vapply(CRAN_PACKAGES, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  message("\nInstalling ", length(missing_cran), " CRAN package(s):")
  message("  ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, dependencies = TRUE)
} else {
  message("\nAll CRAN packages already installed.")
}

# -----------------------------------------------------------------------------
# Bioconductor
# -----------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

missing_bioc <- BIOC_PACKAGES[
  !vapply(BIOC_PACKAGES, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc) > 0) {
  message("\nInstalling ", length(missing_bioc), " Bioconductor package(s):")
  message("  ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
} else {
  message("All Bioconductor packages already installed.")
}

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------

still_missing <- c(
  CRAN_PACKAGES[!vapply(CRAN_PACKAGES, requireNamespace, logical(1), quietly = TRUE)],
  BIOC_PACKAGES[!vapply(BIOC_PACKAGES, requireNamespace, logical(1), quietly = TRUE)]
)

if (length(still_missing) > 0) {
  stop("\nThese packages could not be installed:\n  ",
       paste(still_missing, collapse = ", "),
       "\nInstall them manually before running run_full_analysis.R.")
}

# msigdbr changed its API at version 10.0: scripts in strand 02 call
# msigdbr(db_species=, collection=, subcollection=), which errors with
# "unused arguments" on 7.x. Fail here with a clear message rather than
# forty minutes into a run.
if (utils::packageVersion("msigdbr") < "10.0.0") {
  stop("msigdbr ", utils::packageVersion("msigdbr"), " is too old. ",
       "Strand 02 needs >= 10.0.0 for the collection=/subcollection= API.\n",
       "Update with: install.packages('msigdbr')")
}

message("\nAll dependencies satisfied.")
message("Next: Rscript scripts/00_setup/check_environment.R")
