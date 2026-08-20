# =============================================================================
# 00_config.R
# -----------------------------------------------------------------------------
# Central configuration for the C4-OE cortical bulk RNA-seq analysis.
# Every other script in this pipeline sources this file first, so paths,
# thresholds, and sample metadata only need to be edited in one place.
# =============================================================================

# ---- Project root -----------------------------------------------------------
# Resolved from SCRIPT_DIR (set by the bootstrap at the top of whichever
# 01-04 / run_all.R script sourced this file), NOT from getwd(). This makes
# the pipeline work regardless of your R session's working directory --
# whether you run scripts via Rscript, RStudio's "Source" button (which sets
# wd to the script's own folder), or source() from the repo root.
if (!exists("SCRIPT_DIR")) {
  # Fallback: 00_config.R was sourced directly without going through a
  # 01-04 script's bootstrap. Assume it's being run from its own folder.
  warning("SCRIPT_DIR not set; assuming 00_config.R is being sourced from ",
          "its own directory (scripts/01_bulk_rnaseq_DE/). If paths below ",
          "look wrong, source one of the numbered scripts (01-04) instead, ",
          "or run_all.R, rather than sourcing 00_config.R directly.")
  SCRIPT_DIR <- getwd()
}
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))

DIR_INPUT   <- file.path(PROJECT_ROOT, "Input files")
DIR_OUTPUT  <- file.path(PROJECT_ROOT, "Outputs")

# Outputs/ mirrors scripts/ one-for-one, so every output folder has an obvious
# matching code folder:
#
#   scripts/01_bulk_rnaseq_DE/   ->  Outputs/01_bulk_rnaseq_DE/
#   scripts/02_cross_species_comparison/ -> Outputs/02_cross_species_comparison/
#   scripts/03_scRNAseq_reference_projection/ -> Outputs/03_scRNAseq_reference_projection/
#
# This config governs the first strand only. The other two are wired up by
# run_full_analysis.R at the repo root.
DIR_OUTPUT_BULK <- file.path(DIR_OUTPUT, "01_bulk_rnaseq_DE")
DIR_RDATA   <- file.path(DIR_OUTPUT_BULK, "rdata")  # cached intermediate R objects
DIR_ARCHIVE <- file.path(DIR_OUTPUT, "archive")     # ad hoc files, NOT pipeline-produced

if (!dir.exists(DIR_RDATA)) dir.create(DIR_RDATA, recursive = TRUE)

# ---- Per-step output folders -------------------------------------------------
# Each numbered script writes into its own Outputs/<step>/{tables,figures}
# folder, so it is always obvious which script produced a given file and a
# re-run of one step cannot silently overwrite another step's output.
#
# Scripts do NOT get DIR_TABLES / DIR_FIGURES from this file. Each numbered
# script sets them itself, immediately after sourcing this config:
#
#   DIR_TABLES  <- step_dir("02", "tables")
#   DIR_FIGURES <- step_dir("02", "figures")
#
# That is deliberate: if a script forgets, it fails with "object 'DIR_TABLES'
# not found" instead of quietly writing into a shared folder.
#
# `rdata/` stays shared across steps, because that is exactly what it is for --
# step 02 reads step 01's cache, and steps 03 and 04 both read step 02's.
OUTPUT_STEP_DIRS <- c(
  "01" = "01_qc",
  "02" = "02_differential_expression",
  "03" = "03_enrichment",
  "04" = "04_wgcna"
)

step_dir <- function(step, kind = c("tables", "figures")) {
  kind <- match.arg(kind)
  step <- as.character(step)
  if (!step %in% names(OUTPUT_STEP_DIRS)) {
    stop("Unknown pipeline step: '", step, "'. Known steps: ",
         paste(names(OUTPUT_STEP_DIRS), collapse = ", "),
         ". If this is a new step, add it to OUTPUT_STEP_DIRS in 00_config.R.")
  }
  d <- file.path(DIR_OUTPUT_BULK, OUTPUT_STEP_DIRS[[step]], kind)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

# ---- Input files --------------------------------------------------------
FILE_RAW_COUNTS <- file.path(DIR_INPUT, "raw_counts.csv")

GMT_FILES <- list(
  hallmark        = file.path(DIR_INPUT, "mh.all.v2023.2.Mm.symbols.gmt"),
  canonical_paths = file.path(DIR_INPUT, "m2.all.v2023.2.Mm.symbols.gmt"),
  cp_only         = file.path(DIR_INPUT, "m2.cp.v2023.2.Mm.symbols.gmt"),
  regulatory      = file.path(DIR_INPUT, "m3.all.v2023.2.Mm.symbols.gmt"),
  go_bp           = file.path(DIR_INPUT, "m5.go.bp.v2023.2.Mm.symbols.gmt"),
  cell_type       = file.path(DIR_INPUT, "m8.all.v2023.2.Mm.symbols.gmt")
)

# ---- Sample metadata ---------------------------------------------------------
# Genotype is read from the sample name prefix (WT_ / mC4_ in the raw counts
# file), then relabelled to Control / C4-OE for all reader-facing output.
#
# Sex was NOT recorded. It is inferred here
# from Xist (female-specific, Xist+ / Ddx3y-) vs. Ddx3y (male-specific,
# Y-linked) raw counts, which are unambiguous in this dataset (each gene is
# essentially on/off per sample -- see 01_load_data_and_infer_sex.R for the
# full inference + a diagnostic plot). CONFIRM these against your animal
# records before treating them as ground truth for the manuscript.
#
# IMPORTANT CONFOUND: every C4-OE sample here is male. Control is split
# 2F/2M. Sex and genotype are therefore confounded within the C4-OE group -- adding
# sex as a covariate is the statistically appropriate way to handle this
# (it uses the within-Control sex contrast to estimate and remove the sex
# effect), but it cannot fully disentangle "genotype effect" from "effect
# only visible in males" for genes that are also strongly sex-differential.
# This caveat is written into the DE script's output and should be stated
# in the manuscript's limitations.
SAMPLE_METADATA <- data.frame(
  sample   = c("WT_10", "WT_12", "WT_7", "WT_8", "mC4_4", "mC4_5", "mC4_6", "mC4_7"),
  genotype = c("Control", "Control", "Control", "Control", "C4_OE", "C4_OE", "C4_OE", "C4_OE"),
  sex      = c("F",     "M",     "M",    "F",    "M",     "M",     "M",     "M"),
  stringsAsFactors = FALSE
)
# NOTE on naming: `sample` MUST keep the original WT_/mC4_ strings because
# those are the literal column headers in raw_counts.csv. Everything the
# reader sees is relabelled for display:
#   genotype level "Control" -> "Control", "C4_OE" -> "C4-OE"
#   sample "WT_10" -> "Control_10", "mC4_4" -> "C4-OE_4"
# The internal genotype level is "C4_OE" (underscore, not hyphen) because a
# hyphen is not a syntactically valid R name: DESeq2 silently mangles such
# factor levels when building the model matrix, which would break
# resultsNames()/contrast lookups. Display conversion happens at plot/table
# time via display_genotype() and display_sample() in functions.R.
GENOTYPE_LEVELS  <- c("Control", "C4_OE")             # reference level first
GENOTYPE_DISPLAY <- c(Control = "Control", C4_OE = "C4-OE")
DESEQ_RESULT_NAME <- "genotype_C4_OE_vs_Control"
CONTRAST_LABEL    <- "C4-OE vs Control"

# Group fill/point colours, keyed by DISPLAY label. Every genotype-coloured plot
# in the pipeline (and Figures_paper.R) reads from here, so changing these two
# values recolours the whole analysis consistently.
GENOTYPE_COLOURS <- c("Control" = "#BBD8B4", "C4-OE" = "#FEB751")

# ---- Sex-linked genes -----------------------------------------------------
# Kept as an explicit, documented QC/annotation step (flagging, not silent
# removal) -- see rationale in 02_differential_expression_deseq2.R.
SEX_LINKED_GENES <- c(
  "Xist", "Usp9y", "Ube1y1", "Kdm5d", "Eif2s3y", "Uty", "Ddx3y",
  "Usp9x", "Ube1x", "Kdm5c", "Eif2s3x", "Kdm6a", "Ddx3x"
)
# Note: several of these are legacy gene symbols
# (Smcy -> Kdm5d, Dby -> Ddx3y, Smcx -> Kdm5c, Utx -> Kdm6a, Ube1y -> Ube1y1,
# Ube1x -> retained). Both legacy and current MGI symbols are included so the
# filter is robust to which gene-symbol build produced your annotation.

# ---- Synaptic GO Biological Process definition -------------------------------
# Root GO BP terms used to define "synaptic" for the focused synaptic barplot
# in 03_go_enrichment_ORA_GSEA.R. The full synaptic term set is these roots
# PLUS all of their GO descendants (resolved via GO.db::GOBPOFFSPRING in
# get_synaptic_go_ids()), so specific child terms are captured even when their
# names don't contain the word "synapse" (e.g. "neurotransmitter secretion",
# "long-term synaptic potentiation"). This is more complete and more
# reproducible than grepping term descriptions for keywords.
SYNAPTIC_GO_BP_ROOTS <- c(
  "GO:0099536",  # synaptic signaling
  "GO:0050808",  # synapse organization
  "GO:0050803",  # regulation of synapse structure or activity
  "GO:0099003",  # vesicle-mediated transport in synapse
  "GO:0007269"   # neurotransmitter secretion
)

# ---- WGCNA soft-thresholding power ------------------------------------------
# Set WGCNA_POWER to a number to use it directly, or to NULL to auto-select
# (first power whose SIGNED R^2 reaches WGCNA_RSQ_CUTOFF, falling back to the
# WGCNA authors' sample-size recommendation if none does).
#
# Why 18 for this dataset: the soft-threshold scan on these 8 samples never
# reaches the conventional R^2 = 0.9 anywhere in powers 1-20 -- it plateaus at
# roughly 0.85-0.89 from power ~16 onward. Langfelder & Horvath's guidance for
# exactly this situation (no power attains the cutoff) is to use a default
# based on sample number and network type; for a SIGNED network with fewer
# than 20 samples that default is 18. At power 18 the signed R^2 is ~0.87 and
# the curve has flattened, so it sits on the plateau rather than on the steep
# part of the curve.
#
# For comparison, power 12 gives a signed R^2 of only ~0.78 -- below even the
# more permissive 0.8 threshold, and still on the rising part of the curve.
# Re-running with WGCNA_POWER <- 12 is worthwhile as a sensitivity check;
# module assignments should be compared rather than assumed identical.
WGCNA_POWER       <- 18
WGCNA_RSQ_CUTOFF  <- 0.90

# Number of largest modules to annotate with GO ORA (grey, the unassigned
# bin, is always excluded and does not count toward this total).
WGCNA_N_TOP_MODULES <- 8

# Minimum mean normalized count for a gene to enter the WGCNA network. This is
# a plain low-expression cutoff computed from the normalized counts alone: it
# never looks at genotype, at p-values, or at any DESeq2 model output, so the
# co-expression network is built on an input set chosen independently of the
# differential expression result. (Filtering WGCNA input by differential
# expression is strongly discouraged -- it collapses the correlation structure
# into a single trait-driven module. Filtering by low expression is standard
# and recommended.) At 20 the retained set is a superset of the genes DESeq2
# was able to test, so nothing that was previously analysed is lost.
WGCNA_MIN_MEAN_COUNT <- 20

# ---- WGCNA module resolution -----------------------------------------------
# These two set how finely the co-expression tree is cut into modules. A large
# minimum module size is appropriate for the sample size of this experiment
# (n = 8): with 8 samples every gene-gene correlation rests on 8 points, so
# fine module boundaries are not supported by the data. WGCNA's own default
# (minModuleSize 30) splits these data into 73 modules, most of them small and
# interleaved rather than forming clean contiguous blocks.
#
# Both criteria are computed from the expression correlation structure ALONE
# and never see the genotype labels, so changing them cannot bias the
# module-trait comparison -- it only changes the resolution at which modules
# are defined. Report the values used in the methods.
WGCNA_MERGE_CUT_HEIGHT <- 0.25   # merge modules whose eigengenes correlate > 0.75
WGCNA_MIN_MODULE_SIZE  <- 200    # smallest module retained, in genes

# Number of modules drawn in the manuscript eigengene heatmap (largest first,
# grey excluded). The full set is written alongside it as a separate file, so
# nothing is hidden -- this only keeps the manuscript panel readable.
WGCNA_HEATMAP_N_MODULES <- 10

# ---- Analysis thresholds ---------------------------------------------------
PADJ_THRESHOLD  <- 0.05
LOGFC_THRESHOLD <- 1.0     # used for GO/GSEA gene-list definitions and plots
GO_MIN_SET_SIZE <- 15
GO_MAX_SET_SIZE <- 500

# ---- Reproducibility ---------------------------------------------------------
set.seed(42)
