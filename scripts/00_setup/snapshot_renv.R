# =============================================================================
# scripts/00_setup/snapshot_renv.R
# -----------------------------------------------------------------------------
# Creates renv.lock: the record of exactly which package versions produced the
# published results. This is the single most useful thing in the repository for
# anyone trying to reproduce them.
#
# Run ONCE, from the repository root, on the machine that produced the results,
# AFTER a successful full run:
#
#   Rscript scripts/00_setup/snapshot_renv.R
#
# Then commit renv.lock (and the small renv/ scaffolding files renv creates).
#
# Anyone else then reproduces your library with:
#
#   renv::restore()
#
# Re-run this only when you deliberately change a package version and
# regenerate the results to match.
# =============================================================================

if (is.null(getOption("repos")) || any(getOption("repos") == "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

if (!requireNamespace("renv", quietly = TRUE)) {
  message("Installing renv...")
  install.packages("renv")
}

if (!dir.exists("scripts") || !file.exists("run_full_analysis.R")) {
  stop("Run this from the repository root (the folder containing ",
       "run_full_analysis.R).")
}

# Capture the libraries where the project's packages are ACTUALLY installed,
# BEFORE renv::init() runs. init(bare = TRUE) removes the system/user library
# from .libPaths() for the rest of the session, so a plain snapshot() afterwards
# sees an (almost) empty library and records next to nothing -- a lockfile with
# ~2 packages is the symptom. `.Library` / `.Library.site` are base constants
# that keep pointing at the real system libraries regardless of what renv does
# to .libPaths(), so include them explicitly; we hydrate from these below.
source_libs <- unique(c(.libPaths(), .Library, .Library.site))
source_libs <- source_libs[nzchar(source_libs) & dir.exists(source_libs)]

# Initialise renv if this project is not already using it. bare = TRUE avoids
# renv trying to install everything itself -- the library is already built, we
# only want to record it.
if (!file.exists("renv.lock") && !dir.exists("renv")) {
  message("Initialising renv for this project...")
  renv::init(bare = TRUE, restart = FALSE)
}

# Discover dependencies by scanning the code rather than snapshotting the whole
# user library, so the lockfile describes this project and nothing else.
message("Scanning scripts for package dependencies...")
deps <- unique(renv::dependencies(path = c("scripts", "run_full_analysis.R"),
                                  progress = FALSE)$Package)
message("Found ", length(deps), " packages referenced in code.")

# Bioconductor packages need their repository recorded, or restore() on another
# machine will fail to find them.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
options(renv.config.bioconductor.version = BiocManager::version())
message("Bioconductor version recorded: ", BiocManager::version())

# Populate the project library from the real installs BEFORE snapshotting.
# This is the step that makes the difference: without it, renv::snapshot() below
# records almost nothing because init(bare = TRUE) has already removed the system
# library from the search path. hydrate() copies (hard-links) the discovered
# packages and their dependencies from `source_libs` into the project library.
message("Hydrating project library from installed packages...")
renv::hydrate(packages = deps, sources = source_libs, prompt = FALSE)

# No `type =` argument: type = "explicit" records only packages listed in a
# DESCRIPTION file, and this project has none -- that combination writes a
# near-empty lockfile and silently defeats the point of the script. Passing
# `packages =` alone snapshots exactly the set discovered above.
renv::snapshot(packages = deps, prompt = FALSE)

# Guardrail: a healthy lockfile has far more than `deps` entries (it also records
# every recursive dependency). If it is suspiciously small, the hydrate step
# failed to find the installed packages -- fail loudly rather than commit junk.
n_locked <- length(renv::lockfile_read("renv.lock")$Packages)
message("\nrenv.lock written with ", n_locked, " packages.")
if (n_locked < length(deps)) {
  stop("renv.lock has only ", n_locked, " packages but ", length(deps),
       " were referenced in code. The hydrate step did not find the installed\n",
       "packages -- check that they are installed and on .libPaths(), then re-run.",
       call. = FALSE)
}
message("COMMIT IT, together with .Rprofile and renv/activate.R.")
message("Others then run: renv::restore()")
