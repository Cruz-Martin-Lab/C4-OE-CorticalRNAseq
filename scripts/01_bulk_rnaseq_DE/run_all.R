# =============================================================================
# run_all.R
# -----------------------------------------------------------------------------
# Runs the full C4-OE cortical bulk RNA-seq pipeline top-to-bottom.
#
# Works regardless of your R session's working directory -- run it via
# Rscript, RStudio's "Source" button, or source() from anywhere:
#
#   Rscript scripts/01_bulk_rnaseq_DE/run_all.R
#
# Requires (see also README.md):
#   raw_counts.csv and the six MSigDB .gmt files present in "Input files/"
#   R packages: tidyverse, DESeq2, clusterProfiler, org.Mm.eg.db,
#   AnnotationDbi, GO.db, fgsea, EnhancedVolcano, pheatmap, RColorBrewer, WGCNA,
#   ggrepel, gridExtra, scales
# =============================================================================

.get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
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

steps <- c(
  "01_load_data_and_infer_sex.R",
  "02_differential_expression_deseq2.R",
  "03_go_enrichment_ORA_GSEA.R",
  "04_coexpression_wgcna.R"
)

for (step in steps) {
  message("\n=============================================================")
  message("Running: ", step)
  message("=============================================================")
  source(file.path(SCRIPT_DIR, step))
}

message("\nAll steps complete. Each step wrote to its own folder:")
message("  Outputs/01_qc/                  {tables,figures}")
message("  Outputs/02_differential_expression/ {tables,figures}")
message("  Outputs/03_enrichment/          {tables,figures}")
message("  Outputs/04_wgcna/               {tables,figures}")
message("  Outputs/rdata/                  cached objects, shared between steps")
