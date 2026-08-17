# =============================================================================
# Figures_paper.R
# -----------------------------------------------------------------------------
# Stand-alone, modular generator for hand-picked manuscript figures.
#
# Each figure is a self-contained function registered in PAPER_FIGURES (bottom
# of this file). Sourcing the script regenerates the requested figures WITHOUT
# touching the analysis pipeline: it only READS strand 01's cached results
# (Outputs/01_bulk_rnaseq_DE/rdata/), so it never re-runs anything.
#
#   # make every registered figure
#   source("Figures_paper.R")
#
#   # make only one (or a few) -- set this BEFORE sourcing
#   FIGURES_TO_MAKE <- "c4_expression"
#   source("Figures_paper.R")
#
# To ADD a figure: write a fig_<name>() function that returns a ggplot, then
# add an entry to PAPER_FIGURES. Nothing else needs changing.
#
# Figures are written to Outputs/paper_figures/ as PDF (vector, for the
# manuscript) and PNG (for quick viewing).
# =============================================================================

# ---- Locate this script, independent of the working directory ---------------
# (works via Rscript, RStudio's "Source" button, or source() from anywhere;
# decodes the "~+~" that R's front end substitutes for spaces in --file=).
.get_figures_dir <- function() {
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1) return(dirname(normalizePath(
    gsub("~+~", " ", sub("^--file=", "", file_arg), fixed = TRUE))))
  if (file.exists("run_full_analysis.R")) return(normalizePath(getwd()))
  stop("Could not locate Figures_paper.R automatically; setwd() to the repo root first.")
}

REPO_ROOT  <- .get_figures_dir()
SCRIPT_DIR <- file.path(REPO_ROOT, "scripts", "01_bulk_rnaseq_DE")
source(file.path(SCRIPT_DIR, "00_config.R"))   # provides DIR_RDATA, DIR_OUTPUT, ...

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# ---- Output folder + shared style -------------------------------------------
DIR_PAPER_FIGS <- file.path(DIR_OUTPUT, "paper_figures")
if (!dir.exists(DIR_PAPER_FIGS)) dir.create(DIR_PAPER_FIGS, recursive = TRUE)

# Genotype fill colours come from 00_config.R (GENOTYPE_COLOURS), so the paper
# figures and the pipeline figures share one palette. To recolour everything,
# edit GENOTYPE_COLOURS in scripts/01_bulk_rnaseq_DE/00_config.R.
stopifnot(exists("GENOTYPE_COLOURS"))

save_figure <- function(plot, name, width = 6.5, height = 4.5) {
  pdf_path <- file.path(DIR_PAPER_FIGS, paste0(name, ".pdf"))
  png_path <- file.path(DIR_PAPER_FIGS, paste0(name, ".png"))
  ggsave(pdf_path, plot, width = width, height = height)
  ggsave(png_path, plot, width = width, height = height, dpi = 300)
  message("  wrote ", basename(pdf_path), " + ", basename(png_path))
  invisible(plot)
}

# ---- Cached strand-01 results (read on demand) ------------------------------
.load_deseq_cache <- function() {
  f <- file.path(DIR_RDATA, "02_deseq2_results.rds")
  if (!file.exists(f)) {
    stop("Strand 01's DESeq2 cache not found:\n  ", f,
         "\nRun strand 01 first: source(\"run_full_analysis.R\").", call. = FALSE)
  }
  readRDS(f)   # list: results_table, vst_matrix, norm_counts, metadata
}

# Long-format expression for one gene: one row per sample, with display labels.
# `value` is "normalized" (DESeq2 median-of-ratios counts) or "vst".
.gene_expression_long <- function(gene = "C4b", value = c("normalized", "vst")) {
  value <- match.arg(value)
  d <- .load_deseq_cache()

  expr_tbl <- if (value == "normalized") {
    d$norm_counts                                            # tibble: gene + samples
  } else {
    as.data.frame(d$vst_matrix) %>% rownames_to_column("gene")
  }

  if (!gene %in% expr_tbl$gene) {
    stop("Gene '", gene, "' not found in the ", value,
         " expression table.", call. = FALSE)
  }

  sample_cols <- setdiff(names(expr_tbl), "gene")
  expr_tbl %>%
    filter(gene == !!gene) %>%
    pivot_longer(all_of(sample_cols), names_to = "sample", values_to = "expr") %>%
    left_join(d$metadata, by = "sample") %>%
    mutate(
      genotype_label = factor(genotype_label, levels = c("Control", "C4-OE")),
      # keep the metadata's sample order (Control group first, then C4-OE)
      sample_label   = factor(sample_label, levels = d$metadata$sample_label)
    )
}

# =============================================================================
# FIGURES
# =============================================================================

# --- C4 expression across samples --------------------------------------------
# Bar graph of per-sample expression for one gene (default C4b, the mouse
# complement C4 gene driven up in the C4-OE model), one bar per sample,
# coloured by genotype. `value` picks DESeq2 median-of-ratios normalized counts
# (the natural "expression level"; default) or the variance-stabilized (VST)
# values used for the PCA/heatmap.
fig_c4_expression <- function(gene = "C4b", value = c("normalized", "vst")) {
  value <- match.arg(value)
  df    <- .gene_expression_long(gene, value)
  y_lab <- if (value == "normalized") "Normalized counts" else "Expression (VST)"

  p <- ggplot(df, aes(x = sample_label, y = expr, fill = genotype_label)) +
    geom_col(width = 0.72, colour = "black", linewidth = 0.3) +
    scale_fill_manual(values = GENOTYPE_COLOURS, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = paste0(gene, " expression across samples"),
         x = NULL, y = y_lab) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      plot.title = element_text(face = "bold")
    )

  save_figure(p, paste0("fig_", tolower(gene), "_expression_", value))
}

# --- C4 expression by genotype: box-and-whiskers + one dot per animal ---------
# Box-and-whiskers of a gene's expression per genotype group, with each animal
# (sample) overlaid as a dot. Same data as fig_c4_expression, summarised by
# group. `value` picks normalized counts (default) or VST.
fig_c4_expression_box <- function(gene = "C4b", value = c("normalized", "vst")) {
  value <- match.arg(value)
  df    <- .gene_expression_long(gene, value)
  y_lab <- if (value == "normalized") "Normalized counts" else "Expression (VST)"

  p <- ggplot(df, aes(x = genotype_label, y = expr)) +
    geom_boxplot(aes(fill = genotype_label), width = 0.55, alpha = 0.45,
                 colour = "black", outlier.shape = NA) +
    # one dot per animal; fixed seed keeps the horizontal jitter reproducible
    geom_point(aes(fill = genotype_label),
               position = position_jitter(width = 0.10, height = 0, seed = 1234),
               shape = 21, size = 3, colour = "black", stroke = 0.4) +
    scale_fill_manual(values = GENOTYPE_COLOURS, guide = "none") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(title = paste0(gene, " expression by genotype"), x = NULL, y = y_lab) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  save_figure(p, paste0("fig_", tolower(gene), "_expression_box_", value),
              width = 4.2, height = 4.5)
}

# =============================================================================
# REGISTRY + DISPATCH
# =============================================================================
# name -> a zero-argument thunk that builds and saves the figure. Add new
# figures here; set FIGURES_TO_MAKE (above, before sourcing) to a subset of
# these names, or leave it as "all".
PAPER_FIGURES <- list(
  c4_expression     = function() fig_c4_expression(gene = "C4b", value = "normalized"),
  c4_expression_box = function() fig_c4_expression_box(gene = "C4b", value = "normalized")
)

if (!exists("FIGURES_TO_MAKE")) FIGURES_TO_MAKE <- "all"
.to_make <- if (identical(FIGURES_TO_MAKE, "all")) names(PAPER_FIGURES) else FIGURES_TO_MAKE

.unknown <- setdiff(.to_make, names(PAPER_FIGURES))
if (length(.unknown)) {
  stop("Unknown figure(s): ", paste(.unknown, collapse = ", "),
       ".\nAvailable: ", paste(names(PAPER_FIGURES), collapse = ", "), call. = FALSE)
}

message("Generating paper figure(s): ", paste(.to_make, collapse = ", "))
for (.nm in .to_make) {
  message("- ", .nm)
  PAPER_FIGURES[[.nm]]()
}
message("Done. Figures written to: ", DIR_PAPER_FIGS)
