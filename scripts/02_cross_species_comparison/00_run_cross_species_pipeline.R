# -----------------------------------------------------------------------------
# Master script: run the full cross-species comparison pipeline
#
# This script sources the numbered analysis scripts in order.
# Run it from R or with:
#   Rscript 00_run_cross_species_pipeline.R
# -----------------------------------------------------------------------------
##############################################################################
# Install missing packages
##############################################################################

cran_packages <- c(
  "readr",
  "readxl",
  "dplyr",
  "purrr",
  "stringr",
  "tibble",
  "tidyr",
  "ggplot2",
  "patchwork",
  "msigdbr",
  "fgsea"
)

bioc_packages <- c(
  "biomaRt",
  "GOSemSim",
  "GO.db",
  "org.Hs.eg.db"
)

# --------------------------------------------------------------------------
# Install CRAN packages
# --------------------------------------------------------------------------

missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_cran) > 0) {
  
  message("Installing CRAN packages:")
  
  print(missing_cran)
  
  install.packages(
    missing_cran,
    dependencies = TRUE
  )
}

# --------------------------------------------------------------------------
# Install Bioconductor manager if needed
# --------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  
  install.packages("BiocManager")
  
}

# --------------------------------------------------------------------------
# Install Bioconductor packages
# --------------------------------------------------------------------------

missing_bioc <- bioc_packages[
  !vapply(
    bioc_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_bioc) > 0) {
  
  message("Installing Bioconductor packages:")
  
  print(missing_bioc)
  
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
  
}

message("All required packages are installed.")


get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[[1]]))))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    path <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(path)) {
      return(normalizePath(dirname(path)))
    }
  }

  normalizePath(getwd())
}

find_project_root <- function(start_dir) {
  candidates <- c(start_dir, dirname(start_dir), dirname(dirname(start_dir)))
  candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))

  for (candidate in candidates) {
    has_data_raw <- dir.exists(file.path(candidate, "data_raw"))
    has_results <- dir.exists(file.path(candidate, "results"))
    has_rproj <- length(list.files(candidate, pattern = "\\.Rproj$")) > 0

    if ((has_data_raw && has_results) || has_rproj) {
      return(candidate)
    }
  }

  start_dir
}

find_scripts_dir <- function(project_root, scripts) {
  candidates <- c(project_root, file.path(project_root, "R"))
  candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))

  for (candidate in candidates) {
    if (all(file.exists(file.path(candidate, scripts)))) {
      return(candidate)
    }
  }

  stop(
    "Could not find all numbered scripts in either the project root or the R/ subfolder."
  )
}

script_dir <- get_script_dir()
project_root <- find_project_root(script_dir)

scripts <- c(
  "01_validate_inputs.R",
  "02_standardize_aryal_s2.R",
  "03_standardize_mouse_deseq.R",
  "04_map_human_mouse_orthologs.R",
  "05_test_gene_overlap_and_direction.R",
  "06_compare_ranked_pathways.R",
  "07_compare_related_go_pathways.R",
  "08_build_cross_species_summary.R",
  "09_compare_global_pathway_scores.R",
  "10_compare_hallmark_reactome.R",
  "11_summarize_shared_reactome_modules.R",
  "12_make_main_figure_panel_B.R",
  "13_make_main_figure_panel_C.R",
  "14_make_main_figure_panel_D.R",
  "15_make_main_figure_panel_A.R"
)

scripts_dir <- find_scripts_dir(project_root, scripts)
script_paths <- file.path(scripts_dir, scripts)

setwd(project_root)

message("Project root: ", project_root)
message("Scripts directory: ", scripts_dir)
message("Starting pipeline at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

run_one_script <- function(script_path) {
  message("\n=== Running: ", basename(script_path), " ===")
  start_time <- Sys.time()
  source(script_path, local = new.env(parent = globalenv()))
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  message("=== Finished: ", basename(script_path), " (", round(as.numeric(elapsed), 2), " min) ===")
}

for (script_path in script_paths) {
  run_one_script(script_path)
}

message("\nPipeline completed successfully at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
