# =============================================================================
# scripts/00_setup/check_environment.R
# -----------------------------------------------------------------------------
# Reports the exact versions of everything that affects the numbers, so two
# runs can be compared. It records rather than validates: the authoritative
# version check is renv::restore() against renv.lock.
#
# Run this before an analysis run, and paste the output into any issue report:
#
#   Rscript scripts/00_setup/check_environment.R
# =============================================================================

# Versions used for the published analysis. Update these only when the
# published results are regenerated, never to match a new machine.
REFERENCE_ENVIRONMENT <- list(
  ensembl_release = 116,        # pinned in scripts/02_cross_species_comparison/R/04
  msigdb_gmt_version = "2023.2" # the .gmt files in "Input files/"
)

# Packages whose version genuinely changes results, as opposed to merely
# changing plot cosmetics or messages.
RESULT_CRITICAL <- c(
  "DESeq2",          # dispersion estimation, independent filtering
  "fgsea",           # GSEA p values
  "clusterProfiler", # ORA
  "org.Mm.eg.db",    # gene symbol -> ID mapping; changes the tested universe
  "org.Hs.eg.db",
  "GO.db",           # GO structure; changes term membership
  "msigdbr",         # gene set contents
  "WGCNA",           # module detection
  "biomaRt",
  "GOSemSim"
)

rule <- function() cat(strrep("-", 72), "\n")

cat("\n")
rule()
cat("ENVIRONMENT REPORT\n")
rule()
cat("Date:      ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version: ", R.version.string, "\n")
cat("Platform:  ", R.version$platform, "\n")
cat("Locale:    ", Sys.getlocale("LC_COLLATE"), "\n")

# Collation order affects sort() and therefore any alphabetical tie-breaking.
if (!grepl("^C$|UTF-8", Sys.getlocale("LC_COLLATE"))) {
  cat("  WARNING: an unusual collation locale can change alphabetical\n")
  cat("           tie-breaking in sorted output.\n")
}

cat("\n")
rule()
cat("RESULT-CRITICAL PACKAGES\n")
rule()

for (p in RESULT_CRITICAL) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("  %-18s %s\n", p, as.character(utils::packageVersion(p))))
  } else {
    cat(sprintf("  %-18s NOT INSTALLED\n", p))
  }
}

cat("\n")
rule()
cat("PINNED EXTERNAL DATA\n")
rule()
cat(sprintf("  %-18s %s\n", "Ensembl release", REFERENCE_ENVIRONMENT$ensembl_release))
cat(sprintf("  %-18s %s\n", "MSigDB .gmt", REFERENCE_ENVIRONMENT$msigdb_gmt_version))

# The frozen ortholog table is the single most important reproducibility
# artefact: with it, the shared gene universe is fixed; without it, the
# cross-species numbers depend on what Ensembl serves today.
# NOTE: this path is relative, so run this script from the repository root.
frozen <- file.path("scripts", "02_cross_species_comparison", "reference",
                    "ensembl_human_mouse_orthologs_frozen.csv")

cat(sprintf("  %-18s %s\n", "Frozen orthologs",
            if (file.exists(frozen)) "present (offline, deterministic)"
            else "ABSENT -- R/04 will query Ensembl live"))

if (!file.exists(frozen)) {
  cat("\n  The cross-species gene universe will be whatever Ensembl serves\n")
  cat("  today. After the first run, commit the generated file at:\n    ",
      frozen, "\n")
}

cat("\n")
rule()
cat("renv\n")
rule()

if (file.exists("renv.lock")) {
  cat("  renv.lock present -- restore the pinned library with renv::restore()\n")
  if (requireNamespace("renv", quietly = TRUE)) {
    status <- tryCatch(capture.output(renv::status()), error = function(e) NULL)
    if (!is.null(status)) cat(paste0("  ", status, collapse = "\n"), "\n")
  }
} else {
  cat("  No renv.lock. Package versions are NOT pinned, so another machine\n")
  cat("  may produce different numbers. Create one with:\n")
  cat("    Rscript scripts/00_setup/snapshot_renv.R\n")
}

cat("\n")
rule()
cat("Done. Include this output in any report of a discrepancy.\n")
rule()
cat("\n")
