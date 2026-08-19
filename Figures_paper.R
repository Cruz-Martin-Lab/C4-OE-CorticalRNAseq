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

# GO ORA / GSEA result tables produced by strand 01 step 03.
DIR_ENRICH <- file.path(DIR_OUTPUT, "01_bulk_rnaseq_DE", "03_enrichment", "tables")

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
# SHARED PLOT BUILDERS
# =============================================================================
# Each returns a ggplot (no file side effects) so panels can be saved on their
# own or composed later. Styling mirrors the pipeline: GO ORA barplots follow
# make_go_barplot() (blue->red p.adjust, theme_minimal); GSEA NES barplots
# follow plot_gsea_top() (orangered->pink p.adjust, theme_bw).

# Tidy an MSigDB-style pathway id for display: drop the collection prefix and
# turn underscores into spaces.
.clean_pathway <- function(x) {
  x <- sub("^(HALLMARK|GOBP|GOCC|GOMF|REACTOME|WP|KEGG|BIOCARTA|PID)_", "", x)
  tolower(gsub("_", " ", x))
}

# The collection prefix of an MSigDB-style id ("" if it has none).
.pathway_collection <- function(x) {
  re <- "^(HALLMARK|GOBP|GOCC|GOMF|REACTOME|WP|KEGG|BIOCARTA|PID)_.*$"
  ifelse(grepl(re, x), sub(re, "\\1", x), "")
}

# Different collections curate the same pathway under the same name (e.g.
# WP_CHOLESTEROL_BIOSYNTHESIS and REACTOME_CHOLESTEROL_BIOSYNTHESIS both clean
# to "cholesterol biosynthesis"). Duplicate labels would make ggplot STACK the
# two bars into one impossible NES, so disambiguate them with the collection.
.unique_labels <- function(pathway) {
  lab  <- .clean_pathway(pathway)
  dup  <- lab %in% lab[duplicated(lab)]
  coll <- .pathway_collection(pathway)
  lab[dup] <- paste0(lab[dup], " (", coll[dup], ")")
  lab
}

# GO ORA barplot in make_go_barplot() style. Reads one or more GO ORA CSVs
# (basenames under DIR_ENRICH), keeps significant terms (optionally matching
# `pattern` in the Description), and plots the top `num_paths` by gene Count.
build_go_ora_barplot <- function(files, title, pattern = NULL, include_terms = NULL,
                                 exclude = NULL, num_paths = 20,
                                 padj_threshold = PADJ_THRESHOLD) {
  d <- do.call(rbind, lapply(files, function(f) {
    p <- file.path(DIR_ENRICH, f)
    if (!file.exists(p)) {
      stop("GO ORA table not found:\n  ", p,
           "\nRun strand 01 first: source(\"run_full_analysis.R\").", call. = FALSE)
    }
    read.csv(p, check.names = FALSE)[, c("Description", "Count", "p.adjust")]
  }))
  d <- d %>% filter(!is.na(p.adjust), p.adjust < padj_threshold) %>%
    distinct(Description, .keep_all = TRUE)
  # `pattern` keeps terms matching a regex; `include_terms` keeps an exact set;
  # `exclude` drops terms matching a regex.
  if (!is.null(pattern))       d <- d %>% filter(grepl(pattern, Description, ignore.case = TRUE))
  if (!is.null(include_terms)) d <- d %>% filter(Description %in% include_terms)
  if (!is.null(exclude))       d <- d %>% filter(!grepl(exclude, Description, ignore.case = TRUE))
  d <- d %>% arrange(desc(Count)) %>% head(num_paths)
  if (nrow(d) == 0) stop("No significant GO terms for '", title, "'.", call. = FALSE)

  ggplot(d, aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "blue", high = "red", name = "p.adjust") +
    coord_flip() +
    labs(title = title, x = NULL, y = "Gene count") +
    theme_minimal()
}

# GSEA NES barplot in plot_gsea_top() style. Reads one or more GSEA CSVs
# (basenames under DIR_ENRICH), keeps significant sets (optionally matching
# `pattern`), optionally taking the `top_each` most extreme sets per direction,
# and plots NES coloured by padj. `include_pattern` adds back every significant
# set matching a regex AFTER the top-N trim, so a rank-based selection can never
# drop a gene set that the manuscript text cites. Long set names are wrapped at
# `wrap` characters so the bars keep their width.
build_gsea_nes_barplot <- function(files, title, top_each = NULL, pattern = NULL,
                                   include_pattern = NULL,
                                   direction = c("both", "up", "down"),
                                   padj_threshold = PADJ_THRESHOLD, wrap = 55) {
  direction <- match.arg(direction)
  all_sig <- do.call(rbind, lapply(files, function(f) {
    p <- file.path(DIR_ENRICH, f)
    if (!file.exists(p)) {
      stop("GSEA table not found:\n  ", p,
           "\nRun strand 01 first: source(\"run_full_analysis.R\").", call. = FALSE)
    }
    read.csv(p, check.names = FALSE)[, c("pathway", "NES", "padj")]
  })) %>%
    filter(!is.na(padj), padj < padj_threshold) %>%
    distinct(pathway, .keep_all = TRUE)   # theme tables repeat a set across collections

  d <- all_sig
  if (!is.null(pattern)) d <- d %>% filter(grepl(pattern, pathway, ignore.case = TRUE))
  if (direction == "up")   d <- d %>% filter(NES > 0) %>% arrange(desc(NES))
  if (direction == "down") d <- d %>% filter(NES < 0) %>% arrange(NES)
  if (direction == "both" && !is.null(top_each)) {
    d <- bind_rows(
      d %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(top_each),
      d %>% filter(NES < 0) %>% arrange(NES)       %>% head(top_each))
  }
  if (direction != "both" && !is.null(top_each)) d <- head(d, top_each)
  if (!is.null(include_pattern)) {
    d <- bind_rows(d, all_sig %>%
                     filter(grepl(include_pattern, pathway, ignore.case = TRUE))) %>%
      distinct(pathway, .keep_all = TRUE)
  }
  if (nrow(d) == 0) stop("No significant gene sets for '", title, "'.", call. = FALSE)
  d <- d %>% mutate(label = .unique_labels(pathway))
  if (!is.null(wrap)) {
    d$label <- vapply(d$label,
                      function(s) paste(strwrap(s, width = wrap), collapse = "\n"), "")
  }

  ggplot(d, aes(x = reorder(label, NES), y = NES, fill = padj)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_gradient(low = "orangered", high = "pink", name = "p.adjust") +
    coord_flip() +
    labs(title = title, x = NULL, y = "Normalized Enrichment Score (NES)") +
    theme_bw(base_size = 9)
}

# DESeq2 Wald-statistic ranking used by GSEA (NA-padj genes dropped), as a named
# decreasing numeric vector -- the `stats` argument to fgsea::plotEnrichment().
.ranked_stats <- function() {
  rt <- .load_deseq_cache()$results_table %>%
    filter(!is.na(padj), !is.na(stat), !is.na(gene)) %>%
    distinct(gene, .keep_all = TRUE)
  stats <- rt$stat
  names(stats) <- rt$gene
  sort(stats, decreasing = TRUE)
}

# Genes of one pathway, read from the committed GMT for `collection`
# (a name in GMT_FILES: hallmark, canonical_paths, cp_only, regulatory,
# go_bp, cell_type).
.gmt_pathway <- function(collection, pathway) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required for enrichment-score curves.", call. = FALSE)
  }
  gp <- fgsea::gmtPathways(GMT_FILES[[collection]])
  if (!pathway %in% names(gp)) {
    stop("Pathway '", pathway, "' not found in collection '", collection, "'.", call. = FALSE)
  }
  gp[[pathway]]
}

# GSEA running-enrichment-score curve (fgsea::plotEnrichment) for one pathway.
build_enrichment_curve <- function(collection, pathway, title = NULL) {
  if (is.null(title)) title <- .clean_pathway(pathway)
  fgsea::plotEnrichment(.gmt_pathway(collection, pathway), .ranked_stats()) +
    labs(title = title, x = "Rank in ranked gene list", y = "Enrichment score (ES)") +
    theme(plot.title = element_text(face = "bold", size = 10))
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

# --- FIGURE 3: cholesterol biosynthesis --------------------------------------

# 3A. GSEA running enrichment-score curve for WP_CHOLESTEROL_BIOSYNTHESIS.
fig3a_cholesterol_gsea_curve <- function() {
  p <- build_enrichment_curve("cp_only", "WP_CHOLESTEROL_BIOSYNTHESIS",
                              "WP cholesterol biosynthesis (GSEA)")
  save_figure(p, "fig3a_cholesterol_gsea_curve", width = 6, height = 4.2)
}

# 3B. Cholesterol / sterol gene sets across collections; NES coloured by padj.
fig3b_cholesterol_gsea_nes <- function() {
  p <- build_gsea_nes_barplot("GSEA_theme_cholesterol_lipid.csv",
                              "Cholesterol / sterol gene sets (GSEA)")
  save_figure(p, "fig3b_cholesterol_gsea_nes", width = 8.5, height = 6)
}

# 3C. Cholesterol GO:BP terms among up-regulated genes (make_go_barplot style).
#     Every significant cholesterol / sterol / steroid term in GO_ORA_BP_up.csv
#     (10 terms), which includes the 6 cited in the manuscript text. The two
#     significant glycosphingolipid terms are lipid but not sterol, so they are
#     left out of the cholesterol panel.
fig_cholesterol_go_up <- function() {
  p <- build_go_ora_barplot(
    "GO_ORA_BP_up.csv",
    "Cholesterol-related biological processes (up-regulated genes)",
    include_terms = c(
      # cited in the text
      "cholesterol biosynthetic process",
      "cholesterol biosynthetic process via desmosterol",
      "cholesterol biosynthetic process via lathosterol",
      "zymosterol biosynthetic process",
      "zymosterol metabolic process",
      "steroid biosynthetic process",
      # remaining significant sterol / steroid terms
      "steroid metabolic process",
      "cholesterol metabolic process",
      "sterol metabolic process",
      "sterol biosynthetic process"))
  save_figure(p, "fig_cholesterol_go_up", width = 8.5, height = 5)
}

# 3D. WP cholesterol biosynthesis genes as a lollipop by log2 fold change.
#     Reproduces the pipeline's plot_pathway_lollipop() exactly: coloured by
#     significance (red/blue), point size by |logFC|, genes ordered by logFC.
fig3d_wp_cholesterol_lollipop <- function() {
  d     <- .load_deseq_cache()
  genes <- .gmt_pathway("cp_only", "WP_CHOLESTEROL_BIOSYNTHESIS")
  df <- d$results_table %>%
    filter(gene %in% genes, !is.na(logFC), !is.na(padj)) %>%
    mutate(Significance = factor(
      ifelse(padj < PADJ_THRESHOLD, "Significant", "Not Significant"),
      levels = c("Not Significant", "Significant"))) %>%
    arrange(logFC)
  df$gene_f <- factor(df$gene, levels = df$gene)   # ascending -> largest at top

  p <- ggplot(df, aes(x = logFC, y = gene_f)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = logFC, y = gene_f, yend = gene_f,
                     color = Significance), linewidth = 0.5) +
    geom_point(aes(color = Significance, size = abs(logFC))) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "blue"),
                       name = "Significance") +
    scale_size_continuous(name = "|LogFC|", range = c(1.5, 5)) +
    labs(title = "WP Cholesterol Biosynthesis Genes",
         x = expression(Log[2] ~ "Fold Change (LogFC)"), y = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 9))
  save_figure(p, "fig3d_wp_cholesterol_lollipop", width = 6, height = 5)
}

# --- FIGURE 4: synaptic directional split ------------------------------------

# 4A. Synaptic + axonal GO:BP terms among up-regulated genes. All 6 significant
#     terms of GO_ORA_BP_synaptic_up.csv plus axonogenesis / axon extension,
#     which are cited in the text and live in the broader BP-up table.
#     Note: "regulation of synapse structure or activity" is the parent of
#     "regulation of synapse organization" and carries the same 47 genes; it is
#     part of the synaptic set, so both appear. Drop it here to de-duplicate.
fig4a_go_bp_synaptic_up <- function() {
  p <- build_go_ora_barplot(
    c("GO_ORA_BP_synaptic_up.csv", "GO_ORA_BP_up.csv"),
    "Synaptic and axonal biological processes (up-regulated)",
    include_terms = c(
      # the 6 significant terms of the synaptic GO:BP set
      "regulation of synapse organization",
      "regulation of synapse structure or activity",
      "synapse assembly",
      "neuromuscular junction development",
      "postsynaptic modulation of chemical synaptic transmission",
      "inhibitory synapse assembly",
      # axonal terms cited in the text
      "axonogenesis",
      "axon extension"))
  save_figure(p, "fig4a_go_bp_synaptic_up", width = 8.5, height = 4)
}

# 4B. GO:CC compartments among up-regulated genes -- all 15 significant terms
#     (the postsynaptic / dendritic compartments cited in the text plus the
#     matrix and basolateral-membrane compartments that come with them).
fig4b_go_cc_up <- function() {
  p <- build_go_ora_barplot(
    "GO_ORA_CC_up.csv", "Cellular compartments (up-regulated)")
  save_figure(p, "fig4b_go_cc_up", width = 8.5, height = 5.5)
}

# 4C. Down-regulated GO:BP + GO:CC -- every significant term in both tables
#     (5 BP + cytosolic ribosome): the translation-at-synapse terms and the
#     cytosolic ribosome cited in the text, plus the two amine-catabolism terms.
fig4c_go_down <- function() {
  p <- build_go_ora_barplot(
    c("GO_ORA_BP_down.csv", "GO_ORA_CC_down.csv"),
    "Down-regulated processes and compartments")
  save_figure(p, "fig4c_go_down", width = 8, height = 3.2)
}

# 4D. Paired GSEA curves: an up-regulated synaptic set beside the down-regulated
#     synaptic-translation set. Both GO:BP; swap via the arguments.
fig4d_synaptic_gsea_curves <- function(
    up_pathway   = "GOBP_RECEPTOR_LOCALIZATION_TO_SYNAPSE",
    down_pathway = "GOBP_TRANSLATION_AT_SYNAPSE") {
  up   <- build_enrichment_curve("go_bp", up_pathway)
  down <- build_enrichment_curve("go_bp", down_pathway)
  combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(up, down, ncol = 2)
  } else if (requireNamespace("cowplot", quietly = TRUE)) {
    cowplot::plot_grid(up, down, ncol = 2)
  } else {
    stop("Package 'patchwork' or 'cowplot' is required to pair the curves.", call. = FALSE)
  }
  save_figure(combined, "fig4d_synaptic_gsea_curves", width = 10, height = 4.2)
}

# 4E. Down-regulated canonical pathways (GSEA) -- the ribosomal / translation
#     collapse that the text cites alongside 4D (cytoplasmic ribosomal proteins,
#     NMD, eukaryotic translation initiation among the most negative sets).
fig4e_translation_gsea_nes <- function() {
  p <- build_gsea_nes_barplot(
    "GSEA_cp_only.csv", "Down-regulated canonical pathways (GSEA)",
    direction = "down", top_each = 10)
  save_figure(p, "fig4e_translation_gsea_nes", width = 9, height = 5)
}

# --- FIGURE 5: vascular, adhesion, and rank-based findings --------------------

# 5A. Vascular / endothelial GO:BP terms among up-regulated genes. Curated to
#     the terms cited in the text.
fig5a_go_bp_vascular_up <- function() {
  p <- build_go_ora_barplot(
    "GO_ORA_BP_up.csv", "Vascular / endothelial processes (up-regulated)",
    include_terms = c(
      # cited in the text
      "endothelial cell migration",
      "blood vessel endothelial cell migration",
      "regulation of vasculature development",
      "regulation of angiogenesis",
      "sprouting angiogenesis",
      "endothelium development",
      "establishment of endothelial barrier",
      # remaining significant vascular / endothelial terms
      "circulatory system process",
      "vascular process in circulatory system",
      "regulation of endothelial cell migration",
      "blood vessel diameter maintenance",
      "endothelial cell proliferation",
      "regulation of endothelial cell proliferation",
      "endothelial cell differentiation",
      "regulation of blood vessel endothelial cell migration",
      "negative regulation of angiogenesis",
      "negative regulation of blood vessel morphogenesis",
      "endothelial cell development",
      "cell migration involved in sprouting angiogenesis"))
  save_figure(p, "fig5a_go_bp_vascular_up", width = 9, height = 6.5)
}

# 5B. Molecular functions among up-regulated genes: all 12 significant GO:MF
#     terms (adhesion, integrin and ECM among them) plus the 3 GO:BP adhesion
#     terms cited in the text, which live in the BP-up table.
fig5b_go_mf_up <- function() {
  p <- build_go_ora_barplot(
    c("GO_ORA_MF_up.csv", "GO_ORA_BP_up.csv"),
    "Molecular functions and adhesion (up-regulated)",
    include_terms = c(
      # all 12 significant GO:MF terms
      "cell adhesion molecule binding",
      "signaling receptor regulator activity",
      "peptide binding",
      "integrin binding",
      "glycosaminoglycan binding",
      "extracellular matrix structural constituent",
      "protein folding chaperone",
      "cargo receptor activity",
      "ATP-dependent protein folding chaperone",
      "lipoprotein particle binding",
      "protein-lipid complex binding",
      "opsonin binding",
      # GO:BP adhesion terms cited in the text
      "cell-substrate adhesion",
      "cell-cell adhesion via plasma-membrane adhesion molecules",
      "regulation of cell junction assembly"))
  save_figure(p, "fig5b_go_mf_up", width = 9, height = 5.5)
}

# 5C. Canonical pathways, diverging NES: the 10 most positive and 10 most
#     negative significant sets, plus every significant cell-junction /
#     adhesion / collagen set (`include_pattern`) so the adhesion pathways cited
#     in the text stay on the panel even when they fall outside the top 10.
fig5c_canonical_gsea_nes <- function() {
  p <- build_gsea_nes_barplot(
    "GSEA_cp_only.csv", "Canonical pathways (GSEA)",
    direction = "both", top_each = 10,
    include_pattern = "JUNCTION|ADHERENS|ADHESION|COLLAGEN")
  save_figure(p, "fig5c_canonical_gsea_nes", width = 9, height = 8)
}

# 5D. Hallmark gene sets, diverging NES: all 24 significant Hallmark sets, plus
#     the 3 non-Hallmark interferon / cytokine sets cited in the text (response
#     to interferon beta, interferon signaling, cytokine signaling), which come
#     from the immune theme table. Overlapping sets are de-duplicated.
fig5d_hallmark_gsea_nes <- function() {
  p <- build_gsea_nes_barplot(
    c("GSEA_hallmark.csv", "GSEA_theme_immune.csv"),
    "Hallmark and interferon gene sets (GSEA)")
  save_figure(p, "fig5d_hallmark_gsea_nes", width = 8.5, height = 8)
}

# =============================================================================
# REGISTRY + DISPATCH
# =============================================================================
# name -> a zero-argument thunk that builds and saves the figure. Add new
# figures here; set FIGURES_TO_MAKE (above, before sourcing) to a subset of
# these names, or leave it as "all".
PAPER_FIGURES <- list(
  c4_expression      = function() fig_c4_expression(gene = "C4b", value = "normalized"),
  c4_expression_box  = function() fig_c4_expression_box(gene = "C4b", value = "normalized"),

  # Figure 3 -- cholesterol biosynthesis (3C = cholesterol_go_up)
  fig3a_cholesterol_gsea_curve  = fig3a_cholesterol_gsea_curve,
  fig3b_cholesterol_gsea_nes    = fig3b_cholesterol_gsea_nes,
  cholesterol_go_up             = function() fig_cholesterol_go_up(),
  fig3d_wp_cholesterol_lollipop = fig3d_wp_cholesterol_lollipop,

  # Figure 4 -- synaptic directional split
  fig4a_go_bp_synaptic_up    = fig4a_go_bp_synaptic_up,
  fig4b_go_cc_up             = fig4b_go_cc_up,
  fig4c_go_down              = fig4c_go_down,
  fig4d_synaptic_gsea_curves = fig4d_synaptic_gsea_curves,
  fig4e_translation_gsea_nes = fig4e_translation_gsea_nes,

  # Figure 5 -- vascular, adhesion, rank-based
  fig5a_go_bp_vascular_up  = fig5a_go_bp_vascular_up,
  fig5b_go_mf_up           = fig5b_go_mf_up,
  fig5c_canonical_gsea_nes = fig5c_canonical_gsea_nes,
  fig5d_hallmark_gsea_nes  = fig5d_hallmark_gsea_nes
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
