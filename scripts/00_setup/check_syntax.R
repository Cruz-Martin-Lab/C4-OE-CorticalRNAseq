# =============================================================================
# scripts/00_setup/check_syntax.R
# -----------------------------------------------------------------------------
# Parses every .R file in the project WITHOUT running any of it.
#
# This catches syntax errors -- unbalanced braces, stray commas, broken string
# literals -- in a couple of seconds, instead of thirty minutes into an
# analysis run. Worth doing after any edit, and especially after an edit made
# by someone who could not execute R at the time.
#
# Run from the repository root:
#
#   Rscript scripts/00_setup/check_syntax.R
#
# Exit status is 1 if anything failed to parse, so it also works in CI.
# =============================================================================

if (!dir.exists("scripts") || !file.exists("run_full_analysis.R")) {
  stop("Run this from the repository root (the folder containing ",
       "run_full_analysis.R). Current directory: ", getwd())
}

r_files <- c(
  "run_full_analysis.R",
  list.files("scripts", pattern = "\\.[Rr]$", recursive = TRUE,
             full.names = TRUE)
)
r_files <- unique(r_files[file.exists(r_files)])

cat("\n", strrep("=", 72), "\n", sep = "")
cat("SYNTAX CHECK -- ", length(r_files), " R files\n", sep = "")
cat(strrep("=", 72), "\n\n", sep = "")

failures <- list()

for (f in r_files) {
  result <- tryCatch({
    parse(f)
    NULL
  }, error = function(e) conditionMessage(e))

  if (is.null(result)) {
    cat(sprintf("  [ok]     %s\n", f))
  } else {
    cat(sprintf("  [FAILED] %s\n", f))
    failures[[f]] <- result
  }
}

cat("\n", strrep("=", 72), "\n", sep = "")

if (length(failures) == 0) {
  cat("All ", length(r_files), " files parse cleanly.\n", sep = "")
  cat(strrep("=", 72), "\n\n", sep = "")
} else {
  cat(length(failures), " file(s) FAILED to parse:\n\n", sep = "")
  for (f in names(failures)) {
    cat("--- ", f, " ---\n", sep = "")
    cat(failures[[f]], "\n\n")
  }
  cat(strrep("=", 72), "\n\n", sep = "")
  if (!interactive()) quit(status = 1)
}
