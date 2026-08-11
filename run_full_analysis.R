# =============================================================================
# run_full_analysis.R
# -----------------------------------------------------------------------------
# Master runner for the C4-OE cortical RNA-seq analysis.
#
# Creates the full output tree, checks every required input is present, then
# runs all three analysis strands in dependency order:
#
#   01_bulk_rnaseq_DE                QC -> DESeq2 -> GO/GSEA -> WGCNA
#   02_cross_species_comparison      mouse C4-OE vs human SCZ synapse proteome
#   03_scRNAseq_reference_projection DEG projection onto an Allen sc reference
#
# Strand 02 consumes strand 01's DEG table, so the order matters.
#
# USAGE
# -----
#   Rscript run_full_analysis.R
# or from an R session anywhere:
#   source("run_full_analysis.R")
#
# To run only some strands, set this before sourcing:
#   RUN_STRANDS <- c("01")
#
# Each strand runs in its OWN R process, so one strand's attached packages can
# never mask another's functions. Set ISOLATE_STRANDS <- FALSE to run in the
# current session instead (handy when debugging, but see the warning below).
#
# The working directory does not matter; this script locates the repo itself
# and restores your original working directory when it finishes.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Locate the repository root
# -----------------------------------------------------------------------------

# R's front end encodes spaces in the --file= argument as "~+~". If that is not
# undone, every path derived from it is wrong -- and because dir.create() with
# recursive = TRUE will happily create the bogus tree, the failure stays silent
# until something tries to READ a file that was supposed to already exist.
# This repository's own path contains spaces, so this is not hypothetical.
.decode_r_path <- function(p) gsub("~+~", " ", p, fixed = TRUE)

.get_script_dir <- function() {
  # Passed explicitly by run_strands_isolated(); not subject to the ~+~ quirk.
  root_arg <- grep("^--root=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(root_arg) == 1) {
    return(normalizePath(sub("^--root=", "", root_arg[[1]])))
  }

  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(
      .decode_r_path(sub("^--file=", "", file_arg))
    )))
  }
  for (f in rev(sys.frames())) {
    if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (!is.null(p) && nzchar(p)) return(dirname(normalizePath(p)))
  }
  # Last resort: assume the working directory is the repo root.
  if (dir.exists(file.path(getwd(), "scripts", "01_bulk_rnaseq_DE"))) {
    return(normalizePath(getwd()))
  }
  stop(
    "Could not locate the repository root. Run this file with Rscript, or ",
    "source() it as a whole file rather than running lines individually.\n",
    "Current working directory: ", getwd()
  )
}

MASTER_ROOT <- .get_script_dir()

# Sanity-check the root BEFORE anything creates directories. Without this, a
# mis-resolved path (see .decode_r_path) produces a complete phantom output
# tree somewhere unexpected and the run only fails much later, confusingly.
.expected_markers <- c("run_full_analysis.R",
                       file.path("scripts", "01_bulk_rnaseq_DE", "00_config.R"))
.missing_markers <- .expected_markers[
  !file.exists(file.path(MASTER_ROOT, .expected_markers))
]
if (length(.missing_markers) > 0) {
  stop(
    "Resolved the repository root to:\n  ", MASTER_ROOT,
    "\nbut it does not contain: ", paste(.missing_markers, collapse = ", "),
    "\n\nThis is not the repository. If the path above shows '~+~' where there ",
    "should be\nspaces, R has mangled the --file= argument; report it, as the ",
    "decoder in\n.decode_r_path() should have handled that.",
    call. = FALSE
  )
}

# Remembered so the working directory can be restored. NOTE: on.exit() cannot
# be used here -- at top level there is no function context, so R would discard
# it silently. The per-strand runners register their own on.exit(), which do
# work because those are real function contexts, and main() below restores it.
.original_wd <- getwd()

MASTER_SCRIPTS <- file.path(MASTER_ROOT, "scripts")
MASTER_INPUT   <- file.path(MASTER_ROOT, "Input files")
MASTER_OUTPUT  <- file.path(MASTER_ROOT, "Outputs")

if (!exists("RUN_STRANDS")) RUN_STRANDS <- c("01", "02", "03")

# ---- Session isolation -------------------------------------------------------
# Each strand runs in its OWN R process by default.
#
# This is not fussiness. Strand 01 attaches DESeq2, AnnotationDbi, clusterProfiler
# and WGCNA; those export functions named `count`, `select`, `filter`, `rename`,
# `first`, `slice` and more. Strand 02 is written in dplyr and calls 107 of those
# verbs unqualified. Because `library(dplyr)` is a no-op when dplyr is already
# attached, dplyr never gets promoted back to the front of the search path, and
# strand 02 silently calls the wrong function -- e.g.
#   Error in count(., gene_symbol, name = "n_protein_entries") :
#     Argument 'x' is not a vector: list
# A fresh process per strand removes the whole class of problem, and makes each
# strand independently reproducible rather than dependent on what ran before it.
#
# Set ISOLATE_STRANDS <- FALSE to run everything in the current session (faster
# to debug, but you inherit whatever is already attached).
if (!exists("ISOLATE_STRANDS")) ISOLATE_STRANDS <- TRUE

# When this file re-invokes itself for a single strand it passes --strand=NN.
.strand_arg <- grep("^--strand=", commandArgs(trailingOnly = FALSE), value = TRUE)
IS_CHILD_PROCESS <- length(.strand_arg) == 1L
if (IS_CHILD_PROCESS) {
  RUN_STRANDS <- sub("^--strand=", "", .strand_arg[[1]])
}

RUN_STARTED_AT <- Sys.time()

.rule <- function(char = "=") message(strrep(char, 78))
.header <- function(txt) {
  message("")
  .rule()
  message(txt)
  .rule()
}


# -----------------------------------------------------------------------------
# 1. The full directory tree
# -----------------------------------------------------------------------------
# Outputs/ mirrors scripts/ one-for-one, so every output folder has an obvious
# matching code folder. Declared in one place so the structure is documentation
# as well as instruction.

OUTPUT_TREE <- c(
  # -- strand 01: bulk RNA-seq differential expression ------------------------
  "01_bulk_rnaseq_DE/01_qc/tables",
  "01_bulk_rnaseq_DE/01_qc/figures",
  "01_bulk_rnaseq_DE/02_differential_expression/tables",
  "01_bulk_rnaseq_DE/02_differential_expression/figures",
  "01_bulk_rnaseq_DE/03_enrichment/tables",
  "01_bulk_rnaseq_DE/03_enrichment/figures",
  "01_bulk_rnaseq_DE/04_wgcna/tables",
  "01_bulk_rnaseq_DE/04_wgcna/figures",
  "01_bulk_rnaseq_DE/rdata",

  # -- strand 02: cross-species comparison ------------------------------------
  # These names are fixed by the 15 scripts in that strand, which use paths
  # relative to their working directory. See runStrand02() below.
  "02_cross_species_comparison/data_raw",
  "02_cross_species_comparison/data_processed",
  "02_cross_species_comparison/reference",
  "02_cross_species_comparison/results/tables",
  "02_cross_species_comparison/results/figures",
  "02_cross_species_comparison/results/logs",

  # -- strand 03: scRNA-seq reference projection ------------------------------
  "03_scRNAseq_reference_projection/tables",
  "03_scRNAseq_reference_projection/figures",

  # -- misc -------------------------------------------------------------------
  "archive",      # DIR_ARCHIVE in 00_config.R: files no script produces
  "_provenance"
)

create_output_tree <- function() {
  .header("1. Creating the output tree")
  created <- 0L
  for (rel in OUTPUT_TREE) {
    d <- file.path(MASTER_OUTPUT, rel)
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
      created <- created + 1L
    }
  }
  message(sprintf("  %d folders declared, %d newly created.",
                  length(OUTPUT_TREE), created))
  message("  Root: ", MASTER_OUTPUT)
}


# -----------------------------------------------------------------------------
# 2. Required inputs
# -----------------------------------------------------------------------------
# Checked up front so a missing file fails in seconds rather than forty minutes
# into a run. Large inputs are not in the repository -- see REPRODUCIBILITY.md
# for where to obtain each one.

REQUIRED_INPUTS <- list(
  list(strand = "01", path = "raw_counts.csv",
       what = "Raw gene x sample count matrix (Ensembl IDs)"),
  list(strand = "01", path = "mh.all.v2023.2.Mm.symbols.gmt",
       what = "MSigDB mouse Hallmark"),
  list(strand = "01", path = "m2.all.v2023.2.Mm.symbols.gmt",
       what = "MSigDB mouse curated, all"),
  list(strand = "01", path = "m2.cp.v2023.2.Mm.symbols.gmt",
       what = "MSigDB mouse canonical pathways"),
  list(strand = "01", path = "m3.all.v2023.2.Mm.symbols.gmt",
       what = "MSigDB mouse regulatory targets"),
  list(strand = "01", path = "m5.go.bp.v2023.2.Mm.symbols.gmt",
       what = "MSigDB mouse GO Biological Process"),
  list(strand = "01", path = "m8.all.v2023.2.Mm.symbols.gmt",
       what = "MSigDB mouse cell type signatures"),

  list(strand = "02",
       path = "cross_species/Table S2. Results of MS_MS analysis of DLPFC synapse proteomes.xlsx",
       what = "Aryal et al. Table S2 (SCZ/BP synapse proteome)"),
  list(strand = "02",
       path = "cross_species/Table S3. Results of GSEA of changes in SCZ and BP synapse proteomes.xlsx",
       what = "Aryal et al. Table S3"),
  list(strand = "02",
       path = "cross_species/Table S4. Results of GO analyses of module proteins.xlsx",
       what = "Aryal et al. Table S4"),

  # Optional because this file is very large (~5 GB) and not shipped in the repo.
  # When absent, run_strand_03 skips strand 03 gracefully; when present, strand 03
  # runs. A missing copy therefore warns rather than aborting.
  list(strand = "03", path = "sc_data_C4OE_PCA.rds", optional = TRUE,
       what = "Allen single-cell reference Seurat object")
)

check_inputs <- function(strands) {
  .header("2. Checking required inputs")

  needed <- Filter(function(x) x$strand %in% strands, REQUIRED_INPUTS)
  missing <- character(0)

  for (item in needed) {
    full <- file.path(MASTER_INPUT, item$path)
    optional <- isTRUE(item$optional)
    if (file.exists(full)) {
      size_mb <- round(file.info(full)$size / 1024^2, 1)
      message(sprintf("  [ok]      %-58s %6.1f MB", item$path, size_mb))
    } else if (optional) {
      message(sprintf("  [absent]  %-58s  %s (optional)", item$path, item$what))
    } else {
      message(sprintf("  [MISSING] %-58s  %s", item$path, item$what))
      missing <- c(missing, item$path)
    }
  }

  if (length(missing) > 0) {
    stop(
      "\n", length(missing), " required input file(s) are missing from '",
      MASTER_INPUT, "'.\n",
      "These are too large to ship in the git repository. See ",
      "REPRODUCIBILITY.md\nfor where to download each one.\n\nMissing:\n  - ",
      paste(missing, collapse = "\n  - "),
      call. = FALSE
    )
  }

  message("  All required inputs present.")
}


# -----------------------------------------------------------------------------
# 3. Strand runners
# -----------------------------------------------------------------------------

# `isolate` controls the evaluation environment, and the correct value differs
# between the two strands:
#
#   Strand 01 scripts assign SCRIPT_DIR and then source 00_config.R, which
#   checks exists("SCRIPT_DIR"). That inner source() has no local= argument, so
#   it evaluates in globalenv() -- meaning SCRIPT_DIR must be in globalenv() for
#   the handshake to work. Sourcing them into a private environment silently
#   breaks it: the config's fallback fires and MASTER_ROOT resolves two levels
#   above the repo. So: isolate = FALSE (this is what run_all.R has always done).
#
#   Strand 02 scripts share nothing and rely only on relative paths, so they get
#   a private environment, as their own runner does.
source_step <- function(path, label, isolate = TRUE) {
  message("")
  message("  --> ", label)
  t0 <- Sys.time()
  if (isolate) {
    source(path, local = new.env(parent = globalenv()), chdir = FALSE)
  } else {
    source(path, local = FALSE, chdir = FALSE)
  }
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  message("  <-- ", label, " done (", mins, " min)")
  invisible(mins)
}


run_strand_01 <- function() {
  .header("Strand 01: bulk RNA-seq differential expression")

  d <- file.path(MASTER_SCRIPTS, "01_bulk_rnaseq_DE")
  steps <- c(
    "01_load_data_and_infer_sex.R",
    "02_differential_expression_deseq2.R",
    "03_go_enrichment_ORA_GSEA.R",
    "04_coexpression_wgcna.R"
  )
  # These scripts resolve their own paths from 00_config.R, so they write into
  # Outputs/01_bulk_rnaseq_DE/<step>/ regardless of the working directory.
  # isolate = FALSE: see the note on source_step().
  for (s in steps) source_step(file.path(d, s), s, isolate = FALSE)
}


run_strand_02 <- function() {
  .header("Strand 02: cross-species comparison")

  code_dir <- file.path(MASTER_SCRIPTS, "02_cross_species_comparison")
  work_dir <- file.path(MASTER_OUTPUT, "02_cross_species_comparison")

  # Every path inside the 15 scripts of this strand is relative to the working
  # directory ("data_raw/...", "results/tables/..."). Rather than edit ~2,000
  # lines, we point the working directory at the strand's output folder, so
  # they write where we want without modification.
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)

  # -- assemble data_raw ------------------------------------------------------
  raw_dest <- file.path(work_dir, "data_raw")

  # Remove any stray mouse tables left in data_raw/ by earlier runs or copied
  # in by hand. Leaving one in place makes 01_validate_inputs.R fail with
  # "More than one file matched"; the mouse table always comes from strand 01.
  stale_tables <- c("final_deseq_df.csv", "S1_Table_DESeq2.csv")
  for (f in stale_tables) {
    stale <- file.path(raw_dest, f)
    if (file.exists(stale)) {
      file.remove(stale)
      message("  Removed stray mouse table from data_raw/: ", f)
    }
  }

  # (a) The mouse DEG table comes from strand 01, NOT from Input files (design
  #     ~ sex + genotype). Copying it here is what makes the cross-species
  #     comparison use the DEG table produced by this pipeline.
  mouse_table <- file.path(
    MASTER_OUTPUT, "01_bulk_rnaseq_DE", "02_differential_expression",
    "tables", "S1_Table_DESeq.csv"
  )
  if (!file.exists(mouse_table)) {
    stop("Strand 02 needs strand 01's DEG table, which does not exist:\n  ",
         mouse_table, "\nRun strand 01 first.", call. = FALSE)
  }
  file.copy(mouse_table, file.path(raw_dest, "S1_Table_DESeq.csv"),
            overwrite = TRUE)
  message("  Copied DEG table into data_raw/")

  # (b) The Aryal supplementary workbooks.
  aryal_src <- file.path(MASTER_INPUT, "cross_species")
  aryal <- list.files(aryal_src, pattern = "\\.xlsx$", full.names = TRUE)
  file.copy(aryal, raw_dest, overwrite = TRUE)
  message("  Copied ", length(aryal), " Aryal workbook(s) into data_raw/")

  # (c) The frozen Ensembl ortholog table, if one has been committed. Its
  #     presence makes R/04 skip the live BioMart query entirely, which is what
  #     makes the shared gene universe reproducible. See REPRODUCIBILITY.md.
  frozen_src <- file.path(code_dir, "reference",
                          "ensembl_human_mouse_orthologs_frozen.csv")
  if (file.exists(frozen_src)) {
    file.copy(frozen_src, file.path(work_dir, "reference"), overwrite = TRUE)
    message("  Using frozen Ensembl ortholog table (offline, deterministic)")
  } else {
    message("  NOTE: no frozen ortholog table found. R/04 will query Ensembl ",
            "live and\n        write one. Commit it afterwards so others ",
            "reproduce this universe.")
  }

  # -- run --------------------------------------------------------------------
  setwd(work_dir)

  steps <- sprintf("%02d", 1:15)
  script_files <- list.files(file.path(code_dir, "R"), pattern = "\\.R$")
  ordered <- vapply(steps, function(n) {
    hit <- grep(paste0("^", n, "_"), script_files, value = TRUE)
    if (length(hit) != 1) {
      stop("Expected exactly one script starting '", n, "_' in ",
           file.path(code_dir, "R"), ", found ", length(hit), call. = FALSE)
    }
    hit
  }, character(1))

  for (s in ordered) source_step(file.path(code_dir, "R", s), s)

  # -- copy any newly frozen reference table back into the repo ---------------
  frozen_new <- file.path(work_dir, "reference",
                          "ensembl_human_mouse_orthologs_frozen.csv")
  if (file.exists(frozen_new) && !file.exists(frozen_src)) {
    dir.create(dirname(frozen_src), recursive = TRUE, showWarnings = FALSE)
    file.copy(frozen_new, frozen_src, overwrite = TRUE)
    message("\n  A frozen ortholog table was written to:\n    ", frozen_src,
            "\n  COMMIT IT -- it pins the shared gene universe for everyone else.")
  }
}


run_strand_03 <- function() {
  .header("Strand 03: scRNA-seq reference projection")

  script <- file.path(MASTER_SCRIPTS, "03_scRNAseq_reference_projection",
                      "Projection.R")

  # The single-cell reference (~5 GB) is not shipped in the repository. If it is
  # absent, skip strand 03 gracefully so strands 01-02 still complete on a
  # machine without it, rather than aborting the whole run.
  sc_ref <- file.path(MASTER_INPUT, "sc_data_C4OE_PCA.rds")
  if (!file.exists(sc_ref)) {
    message("  SKIPPED. The single-cell reference is not present:")
    message("    ", sc_ref)
    message("  Place sc_data_C4OE_PCA.rds in 'Input files/' to enable strand 03.")
    return(invisible(FALSE))
  }

  # Strand 03 projects strand 01's DEG table onto the reference, so that table
  # must already exist.
  mouse_table <- file.path(
    MASTER_OUTPUT, "01_bulk_rnaseq_DE", "02_differential_expression",
    "tables", "S1_Table_DESeq.csv"
  )
  if (!file.exists(mouse_table)) {
    stop("Strand 03 needs strand 01's DEG table, which does not exist:\n  ",
         mouse_table, "\nRun strand 01 first.", call. = FALSE)
  }

  # Projection.R resolves its input/output paths from MASTER_INPUT / MASTER_OUTPUT
  # (both visible here) and writes to Outputs/03_.../{figures,tables}. The setwd
  # is kept for parity with the other strands; the script uses absolute paths.
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(file.path(MASTER_OUTPUT, "03_scRNAseq_reference_projection"))
  source_step(script, "Projection.R")
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# 4. Provenance record
# -----------------------------------------------------------------------------

write_provenance <- function(strands) {
  .header("Writing provenance record")

  prov_dir <- file.path(MASTER_OUTPUT, "_provenance")
  stamp <- format(RUN_STARTED_AT, "%Y%m%d_%H%M%S")

  writeLines(
    capture.output(sessionInfo()),
    file.path(prov_dir, paste0("sessionInfo_", stamp, ".txt"))
  )

  # Every file produced, with size and modification time, so two runs can be
  # diffed against each other directly.
  produced <- list.files(MASTER_OUTPUT, recursive = TRUE, full.names = TRUE)
  produced <- produced[!grepl("/_provenance/", produced, fixed = TRUE)]

  if (length(produced) > 0) {
    info <- file.info(produced)
    manifest <- data.frame(
      # substring(), not sub(): MASTER_OUTPUT is a path, and this project's own
      # path contains "(Sonia)". Interpolated into a regex those parentheses
      # become a capture group, the pattern fails to match, and the manifest
      # silently records absolute paths instead of relative ones.
      file = substring(produced, nchar(MASTER_OUTPUT) + 2L),
      size_bytes = info$size,
      modified = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
      stringsAsFactors = FALSE
    )
    manifest <- manifest[order(manifest$file), ]
    write.csv(manifest,
              file.path(prov_dir, paste0("output_manifest_", stamp, ".csv")),
              row.names = FALSE)
    message("  ", nrow(manifest), " output files recorded.")
  }

  writeLines(
    c(
      paste("Run started:  ", format(RUN_STARTED_AT, "%Y-%m-%d %H:%M:%S")),
      paste("Run finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste("Strands run:  ", paste(strands, collapse = ", ")),
      paste("R version:    ", R.version.string),
      paste("Platform:     ", R.version$platform),
      paste("Project root: ", MASTER_ROOT)
    ),
    file.path(prov_dir, paste0("run_summary_", stamp, ".txt"))
  )

  message("  Provenance written to Outputs/_provenance/")
}


# -----------------------------------------------------------------------------
# 5. Go
# -----------------------------------------------------------------------------
# Wrapped in a function so on.exit() has a real context: the working directory
# is then restored even if a strand throws.

# Re-invoke this same file once per strand, each in a clean R process.
run_strands_isolated <- function(strands) {
  rscript <- file.path(
    R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  )
  if (!file.exists(rscript)) {
    stop("Could not find Rscript at ", rscript,
         ". Set ISOLATE_STRANDS <- FALSE to run in the current session ",
         "instead (see the note at the top of this file).", call. = FALSE)
  }

  self <- file.path(MASTER_ROOT, "run_full_analysis.R")

  for (s in strands) {
    .header(paste0("Launching strand ", s, " in a fresh R session"))
    # --root is passed explicitly so the child never has to recover the
    # repository path from --file=, which mangles spaces (see .decode_r_path).
    status <- system2(rscript, args = c(
      shQuote(self),
      paste0("--strand=", s),
      shQuote(paste0("--root=", MASTER_ROOT))
    ))
    if (!identical(as.integer(status), 0L)) {
      stop("Strand ", s, " failed (exit status ", status,
           "). See the output above for the error.", call. = FALSE)
    }
  }

  .header("All requested strands completed")
  message("Each strand ran in its own R session and wrote its own provenance")
  message("record to Outputs/_provenance/.")
  invisible(TRUE)
}


main <- function() {

  on.exit(setwd(.original_wd), add = TRUE)

  # Parent process: hand each strand to a clean child and stop here.
  if (ISOLATE_STRANDS && !IS_CHILD_PROCESS) {
    return(run_strands_isolated(RUN_STRANDS))
  }

  .header(paste0("C4-OE cortical RNA-seq | full analysis | strands: ",
                 paste(RUN_STRANDS, collapse = ", ")))
  message("Repository: ", MASTER_ROOT)
  message("Started:    ", format(RUN_STARTED_AT, "%Y-%m-%d %H:%M:%S"))
  if (IS_CHILD_PROCESS) {
    message("Session:    isolated child process for strand ", RUN_STRANDS)
  }

  create_output_tree()
  check_inputs(RUN_STRANDS)

  if ("01" %in% RUN_STRANDS) run_strand_01()
  if ("02" %in% RUN_STRANDS) run_strand_02()
  if ("03" %in% RUN_STRANDS) run_strand_03()

  setwd(.original_wd)
  write_provenance(RUN_STRANDS)

  total <- round(
    as.numeric(difftime(Sys.time(), RUN_STARTED_AT, units = "mins")), 1
  )
  .header(paste0("Analysis complete in ", total, " min"))
  message("Outputs: ", MASTER_OUTPUT)
  message("Provenance and the full output manifest: Outputs/_provenance/")

  invisible(TRUE)
}

main()
