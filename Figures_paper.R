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

# Base-graphics figures cannot go through ggsave(): they draw straight to a
# device. `draw` is a zero-argument function that issues the plotting calls.
save_base_figure <- function(draw, name, width = 9, height = 5, res = 300) {
  pdf_path <- file.path(DIR_PAPER_FIGS, paste0(name, ".pdf"))
  png_path <- file.path(DIR_PAPER_FIGS, paste0(name, ".png"))
  pdf(pdf_path, width = width, height = height); draw(); dev.off()
  png(png_path, width = width * res, height = height * res, res = res); draw(); dev.off()
  message("  wrote ", basename(pdf_path), " + ", basename(png_path))
  invisible(NULL)
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

# Display labels for a set of pathway ids, tagged with their collection where
# that matters:
#   * Different collections curate the same pathway under the same name (e.g.
#     WP_CHOLESTEROL_BIOSYNTHESIS and REACTOME_CHOLESTEROL_BIOSYNTHESIS both
#     clean to "cholesterol biosynthesis"). Duplicate labels would make ggplot
#     STACK the two bars into one impossible NES, so those are always tagged.
#   * `mark_collections` additionally tags every set that is NOT from the
#     panel's dominant collection, so a mostly-Hallmark panel does not silently
#     present a GO:BP or Reactome set as if it were a Hallmark one.
.set_labels <- function(pathway, mark_collections = FALSE) {
  lab  <- .clean_pathway(pathway)
  coll <- .pathway_collection(pathway)
  tag  <- lab %in% lab[duplicated(lab)]
  if (isTRUE(mark_collections) && any(nzchar(coll))) {
    dominant <- names(sort(table(coll[nzchar(coll)]), decreasing = TRUE))[1]
    tag <- tag | (nzchar(coll) & coll != dominant)
  }
  lab[tag] <- paste0(lab[tag], " (", coll[tag], ")")
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
                                   include_pattern = NULL, exclude_sets = NULL,
                                   mark_collections = FALSE,
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

  # Drop sets shown on another figure BEFORE any top-N cut, so the next-ranked
  # sets move up to fill their places rather than leaving the panel short.
  # Matched on the pathway id AND on the cleaned display name, so the same
  # pathway curated by two collections (WP_X / REACTOME_X) is caught either way.
  if (!is.null(exclude_sets) && length(exclude_sets)) {
    all_sig <- all_sig %>%
      filter(!pathway %in% exclude_sets,
             !.clean_pathway(pathway) %in% .clean_pathway(exclude_sets))
  }

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
    # Added-back sets obey `direction` too: an up-only panel must never gain a
    # negative set just because its NAME matches (e.g. "JUNCTION" also occurs in
    # REACTOME_..._EXON_JUNCTION_COMPLEX_EJC, a strongly negative NMD set).
    extra <- all_sig %>% filter(grepl(include_pattern, pathway, ignore.case = TRUE))
    if (direction == "up")   extra <- extra %>% filter(NES > 0)
    if (direction == "down") extra <- extra %>% filter(NES < 0)
    d <- bind_rows(d, extra) %>% distinct(pathway, .keep_all = TRUE)
  }
  if (nrow(d) == 0) stop("No significant gene sets for '", title, "'.", call. = FALSE)
  d <- d %>% mutate(label = .set_labels(pathway, mark_collections))
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
fig_c4_expression <- function(gene = "C4b", value = c("normalized", "vst"),
                                save = TRUE) {
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

  if (save) save_figure(p, paste0("fig_", tolower(gene), "_expression_", value))
  p
}

# --- C4 expression by genotype: box-and-whiskers + one dot per animal ---------
# Box-and-whiskers of a gene's expression per genotype group, with each animal
# (sample) overlaid as a dot. Same data as fig_c4_expression, summarised by
# group. `value` picks normalized counts (default) or VST.
fig_c4_expression_box <- function(gene = "C4b", value = c("normalized", "vst"),
                                    save = TRUE) {
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

  if (save) save_figure(p, paste0("fig_", tolower(gene), "_expression_box_", value),
                        width = 4.2, height = 4.5)
  p
}

# --- FIGURE 3: cholesterol biosynthesis --------------------------------------

# 3A. GSEA running enrichment-score curve for WP_CHOLESTEROL_BIOSYNTHESIS.
fig3a_cholesterol_gsea_curve <- function(save = TRUE) {
  p <- build_enrichment_curve("cp_only", "WP_CHOLESTEROL_BIOSYNTHESIS",
                              "WP cholesterol biosynthesis (GSEA)")
  if (save) save_figure(p, "fig3a_cholesterol_gsea_curve", width = 6, height = 4.2)
  p
}

# 3B. Cholesterol / sterol gene sets across collections; NES coloured by padj.
fig3b_cholesterol_gsea_nes <- function(save = TRUE) {
  p <- build_gsea_nes_barplot("GSEA_theme_cholesterol_lipid.csv",
                              "Cholesterol / sterol gene sets (GSEA)")
  if (save) save_figure(p, "fig3b_cholesterol_gsea_nes", width = 8.5, height = 6)
  p
}

# 3C. Cholesterol GO:BP terms among up-regulated genes (make_go_barplot style).
#     Every significant cholesterol / sterol / steroid term in GO_ORA_BP_up.csv
#     (10 terms), which includes the 6 cited in the manuscript text. The two
#     significant glycosphingolipid terms are lipid but not sterol, so they are
#     left out of the cholesterol panel.
fig_cholesterol_go_up <- function(save = TRUE) {
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
  if (save) save_figure(p, "fig_cholesterol_go_up", width = 8.5, height = 5)
  p
}

# 3D. WP cholesterol biosynthesis genes as a lollipop by log2 fold change.
#     Reproduces the pipeline's plot_pathway_lollipop() exactly: coloured by
#     significance (red/blue), point size by |logFC|, genes ordered by logFC.
fig3d_wp_cholesterol_lollipop <- function(save = TRUE) {
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
  if (save) save_figure(p, "fig3d_wp_cholesterol_lollipop", width = 6, height = 5)
  p
}

# --- FIGURE 4: synaptic directional split ------------------------------------

# 4A. Synaptic + axonal GO:BP terms among up-regulated genes. All 6 significant
#     terms of GO_ORA_BP_synaptic_up.csv plus axonogenesis / axon extension,
#     which are cited in the text and live in the broader BP-up table.
#     Note: "regulation of synapse structure or activity" is the parent of
#     "regulation of synapse organization" and carries the same 47 genes; it is
#     part of the synaptic set, so both appear. Drop it here to de-duplicate.
fig4a_go_bp_synaptic_up <- function(save = TRUE) {
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
  if (save) save_figure(p, "fig4a_go_bp_synaptic_up", width = 8.5, height = 4)
  p
}

# 4B. GO:CC compartments among up-regulated genes -- the significant terms
#     (the postsynaptic / dendritic compartments cited in the text plus the
#     matrix and basolateral-membrane compartments that come with them).
#     "Schaffer collateral - CA1 synapse" is excluded by request.
fig4b_go_cc_up <- function(save = TRUE) {
  p <- build_go_ora_barplot(
    "GO_ORA_CC_up.csv", "Cellular compartments (up-regulated)",
    exclude = "^Schaffer collateral")
  if (save) save_figure(p, "fig4b_go_cc_up", width = 8.5, height = 5.2)
  p
}

# 4C. Down-regulated GO:BP + GO:CC -- the significant terms of both tables:
#     translation at synapse and the cytosolic ribosome cited in the text, plus
#     the two amine-catabolism terms.
#     "translation at presynapse" and "translation at postsynapse" are excluded:
#     they are the same 11 genes at the same p.adjust as "translation at
#     synapse", so they add three identical bars instead of one.
fig4c_go_down <- function(save = TRUE) {
  p <- build_go_ora_barplot(
    c("GO_ORA_BP_down.csv", "GO_ORA_CC_down.csv"),
    "Down-regulated processes and compartments",
    exclude = "translation at (pre|post)synapse")
  if (save) save_figure(p, "fig4c_go_down", width = 8, height = 2.6)
  p
}

# 4D. Paired GSEA curves: an up-regulated synaptic set beside the down-regulated
#     synaptic-translation set. Both GO:BP; swap via the arguments.
fig4d_synaptic_gsea_curves <- function(
    up_pathway   = "GOBP_RECEPTOR_LOCALIZATION_TO_SYNAPSE",
    down_pathway = "GOBP_TRANSLATION_AT_SYNAPSE", ncol = 2, save = TRUE) {
  up   <- build_enrichment_curve("go_bp", up_pathway)
  down <- build_enrichment_curve("go_bp", down_pathway)
  combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(up, down, ncol = ncol)
  } else if (requireNamespace("cowplot", quietly = TRUE)) {
    cowplot::plot_grid(up, down, ncol = ncol)
  } else {
    stop("Package 'patchwork' or 'cowplot' is required to pair the curves.", call. = FALSE)
  }
  if (save) save_figure(combined, "fig4d_synaptic_gsea_curves", width = 10, height = 4.2)
  combined
}

# 4E. Down-regulated canonical pathways (GSEA) -- the ribosomal / translation
#     collapse that the text cites alongside 4D (cytoplasmic ribosomal proteins,
#     NMD, eukaryotic translation initiation among the most negative sets).
fig4e_translation_gsea_nes <- function(save = TRUE) {
  p <- build_gsea_nes_barplot(
    "GSEA_cp_only.csv", "Down-regulated canonical pathways (GSEA)",
    direction = "down", top_each = 10)
  if (save) save_figure(p, "fig4e_translation_gsea_nes", width = 9, height = 5)
  p
}

# --- FIGURE 5: vascular, adhesion, and rank-based findings --------------------

# Gene sets that already appear on another figure. 5C is the broad "everything
# else" canonical-pathway panel, so it should not repeat what Figure 3 and panel
# 5D already show; `build_gsea_nes_barplot(exclude_sets = )` drops these before
# the top-N cut, so the next-ranked sets take their place.
#
# Read from the tables the other panels use, so this follows them automatically
# if those panels change.
#
# NOTE: fig4e_translation_gsea_nes is deliberately NOT counted here. It is a
# stand-alone panel that is not on any assembled figure, and its sets are 5C's
# entire negative half -- counting it would empty that half of the panel.
.sets_shown_elsewhere <- function(padj_threshold = PADJ_THRESHOLD) {
  tables <- c("GSEA_theme_cholesterol_lipid.csv",                 # panel 3B
              "GSEA_hallmark.csv", "GSEA_theme_immune.csv")       # panel 5D
  from_tables <- unlist(lapply(tables, function(f) {
    path <- file.path(DIR_ENRICH, f)
    if (!file.exists(path)) return(character(0))
    d <- read.csv(path, check.names = FALSE)
    d$pathway[!is.na(d$padj) & d$padj < padj_threshold]
  }))
  unique(c(from_tables,
           "WP_CHOLESTEROL_BIOSYNTHESIS",             # panel 3A curve
           "GOBP_RECEPTOR_LOCALIZATION_TO_SYNAPSE",   # panel 4D curves
           "GOBP_TRANSLATION_AT_SYNAPSE"))
}

# 5A. Vascular / endothelial GO:BP terms among up-regulated genes. Curated to
#     the terms cited in the text.
fig5a_go_bp_vascular_up <- function(save = TRUE) {
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
  if (save) save_figure(p, "fig5a_go_bp_vascular_up", width = 9, height = 6.5)
  p
}

# 5B. Molecular functions among up-regulated genes: all 12 significant GO:MF
#     terms (adhesion, integrin and ECM among them) plus the 3 GO:BP adhesion
#     terms cited in the text, which live in the BP-up table.
fig5b_go_mf_up <- function(save = TRUE) {
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
  if (save) save_figure(p, "fig5b_go_mf_up", width = 9, height = 5.5)
  p
}

# 5C. Canonical pathways, diverging NES: the 10 most positive and 10 most
#     negative significant sets, plus every significant cell-junction /
#     adhesion / collagen set (`include_pattern`) so the adhesion pathways cited
#     in the text stay on the panel even when they fall outside the top 10.
fig5c_canonical_gsea_nes <- function(save = TRUE) {
  p <- build_gsea_nes_barplot(
    "GSEA_cp_only.csv", "Canonical pathways (GSEA)",
    direction = "both", top_each = 10,
    include_pattern = "JUNCTION|ADHERENS|ADHESION|COLLAGEN",
    exclude_sets = .sets_shown_elsewhere())
  if (save) save_figure(p, "fig5c_canonical_gsea_nes", width = 9, height = 8)
  p
}

# 5D. Hallmark gene sets, diverging NES: all 24 significant Hallmark sets, plus
#     the 3 non-Hallmark interferon / cytokine sets cited in the text (response
#     to interferon beta, interferon signaling, cytokine signaling), which come
#     from the immune theme table. Overlapping sets are de-duplicated.
fig5d_hallmark_gsea_nes <- function(save = TRUE) {
  p <- build_gsea_nes_barplot(
    c("GSEA_hallmark.csv", "GSEA_theme_immune.csv"),
    "Hallmark and interferon gene sets (GSEA)",
    mark_collections = TRUE)
  if (save) save_figure(p, "fig5d_hallmark_gsea_nes", width = 8.5, height = 8)
  p
}

# --- FIGURE 6: WGCNA co-expression modules -----------------------------------
# Reads strand 01 step 04's cached network and the module tables it writes;
# nothing here re-runs WGCNA.

DIR_WGCNA <- file.path(DIR_OUTPUT, "01_bulk_rnaseq_DE", "04_wgcna")

.load_wgcna_cache <- function() {
  f <- file.path(DIR_RDATA, "04_wgcna_results.rds")
  if (!file.exists(f)) {
    stop("WGCNA cache not found:\n  ", f,
         "\nRun strand 01 step 04 first:\n",
         "  source(\"scripts/01_bulk_rnaseq_DE/04_coexpression_wgcna.R\")", call. = FALSE)
  }
  readRDS(f)   # list: net, module_df, MEs, picked_power
}

.wgcna_table <- function(file) {
  p <- file.path(DIR_WGCNA, "tables", file)
  if (!file.exists(p)) {
    stop("WGCNA table not found:\n  ", p,
         "\nRun strand 01 step 04 (and 04b for the stability column).", call. = FALSE)
  }
  read.csv(p, check.names = FALSE)
}

# Module sizes, largest first, grey excluded.
.wgcna_modules_by_size <- function(module_df) {
  tb <- sort(table(module_df$module), decreasing = TRUE)
  setdiff(names(tb), "grey")
}

# 6A. Gene clustering dendrogram with the module colour bar underneath.
#
#     NOTE: WGCNA::plotDendroAndColors() is base graphics and splits the device
#     with layout(). That call does not survive being replayed into a grid
#     viewport, so this panel CANNOT be composed with patchwork -- doing so
#     collapses the tree and draws the colour bar over it. It is therefore a
#     stand-alone figure, written straight to its own device, and figure6()
#     leaves it out. Place it manually alongside the composite, or use it as a
#     supplementary panel.
build_wgcna_dendrogram <- function() {
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' is required for the dendrogram panel.", call. = FALSE)
  }
  w    <- .load_wgcna_cache()
  net  <- w$net
  cols <- WGCNA::labels2colors(net$colors)
  function() {
    WGCNA::plotDendroAndColors(
      net$dendrograms[[1]], cols[net$blockGenes[[1]]], "Module colors",
      dendroLabels = FALSE, hang = 0.05, addGuide = TRUE, guideHang = 0.05,
      main = "Gene clustering and module assignment")
  }
}

# 6B. Module eigengenes across samples, largest `n_modules` modules only,
#     annotated by genotype and sex. pheatmap draws with grid, so the gtable is
#     wrapped for patchwork the same way.
build_wgcna_eigengene_heatmap <- function(n_modules = 10) {
  for (pkg in c("pheatmap", "RColorBrewer", "WGCNA")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for the eigengene panel.", call. = FALSE)
    }
  }
  .require_patchwork()
  w    <- .load_wgcna_cache()
  meta <- .load_deseq_cache()$metadata
  MEs  <- WGCNA::orderMEs(w$MEs)
  sid  <- rownames(MEs)

  mat <- t(scale(MEs))
  colnames(mat) <- meta$sample_label[match(sid, meta$sample)]
  ann <- data.frame(
    Genotype = factor(meta$genotype_label[match(sid, meta$sample)],
                       levels = names(GENOTYPE_COLOURS)),
    Sex      = meta$sex[match(sid, meta$sample)])
  rownames(ann) <- colnames(mat)

  mods <- .wgcna_modules_by_size(w$module_df)
  keep <- paste0("ME", head(mods, n_modules))
  keep <- keep[keep %in% rownames(mat)]

  ph <- pheatmap::pheatmap(
    mat[keep, , drop = FALSE], annotation_col = ann, silent = TRUE,
    annotation_colors = list(Genotype = GENOTYPE_COLOURS),
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100),
    main = sprintf("Module eigengenes (z-scored)\n%d largest of %d modules (grey excluded)",
                    length(keep), length(mods)))
  patchwork::wrap_elements(ph$gtable)
}

# 6C. Module-theme dot plot: each displayed module's strongest GO BP terms, with
#     a term drawn in every displayed module where it is significant -- so the
#     panel shows how module-specific each theme is.
#
#     `select` picks which modules appear. "stability" uses the power-sensitivity
#     score (04b) and is the defensible choice for the manuscript: it shows the
#     modules that survive the parameter choice. "terms" falls back to ranking by
#     number of significant terms when that column is absent. Modules surviving
#     BH correction for genotype are always kept, so a small but trait-associated
#     module is never dropped.
build_wgcna_module_dotplot <- function(n_modules = 8, n_terms = 3,
                                        select = c("stability", "terms"),
                                        wrap = 42) {
  select  <- match.arg(select)
  go      <- .wgcna_table("GO_ORA_modules_combined.csv")
  summ    <- .wgcna_table("wgcna_module_summary.csv")
  if (select == "stability" && !"stability" %in% names(summ)) {
    message("  (no stability column yet -- run 04b; ranking modules by term count instead)")
    select <- "terms"
  }

  n_sig <- table(go$module)
  annotated <- names(n_sig)
  rank_by <- if (select == "stability") {
    st <- summ$stability[match(annotated, summ$module)]
    annotated[order(-st)]
  } else {
    annotated[order(-as.integer(n_sig[annotated]))]
  }
  always <- summ$module[!is.na(summ$padj_genotype) & summ$padj_genotype < PADJ_THRESHOLD]
  show   <- unique(c(head(rank_by, n_modules), intersect(always, annotated)))
  show   <- rank_by[rank_by %in% show]                  # keep the ranking order

  terms <- unique(unlist(lapply(show, function(m) {
    d <- go[go$module == m, ]
    head(d$Description[order(d$p.adjust)], n_terms)
  })))
  dp <- go[go$module %in% show & go$Description %in% terms, ]
  if (!nrow(dp)) stop("No significant GO terms for the module dot plot.", call. = FALSE)

  wrap_term <- function(x) vapply(x, function(z) paste(strwrap(z, width = wrap), collapse = "\n"), "")
  terms_w <- wrap_term(terms)
  dp$module      <- factor(dp$module, levels = show)
  dp$Description <- factor(wrap_term(as.character(dp$Description)), levels = rev(terms_w))

  ggplot(dp, aes(x = module, y = Description, size = Count, colour = -log10(p.adjust))) +
    geom_point() +
    scale_colour_gradient(low = "steelblue", high = "firebrick",
                           name = expression(-log[10] ~ "p.adjust")) +
    scale_size_continuous(name = "Genes", range = c(1.5, 6)) +
    labs(title = "Module GO biological process enrichment", x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_line(colour = "grey92"))
}

# Base graphics: returns the draw function, and saving goes through
# save_base_figure() rather than ggsave().
fig6a_wgcna_dendrogram <- function(save = TRUE) {
  draw <- build_wgcna_dendrogram()
  if (save) save_base_figure(draw, "fig6a_wgcna_dendrogram", width = 9, height = 5)
  invisible(draw)
}

fig6b_wgcna_eigengene_heatmap <- function(save = TRUE) {
  p <- build_wgcna_eigengene_heatmap(n_modules = 10)
  if (save) save_figure(p, "fig6b_wgcna_eigengene_heatmap", width = 7.5, height = 5.5)
  p
}

fig6c_wgcna_module_dotplot <- function(save = TRUE) {
  p <- build_wgcna_module_dotplot(n_modules = 8, n_terms = 3)
  if (save) save_figure(p, "fig6c_wgcna_module_dotplot", width = 8, height = 9)
  p
}

# =============================================================================
# ASSEMBLED FIGURES
# =============================================================================
# One builder per manuscript figure: it lays the panels out on a single page and
# labels them A, B, C, ... in the style of Figure 2 -- a large, plain (not bold)
# letter at the top-left of each panel.
#
# The panels come from the SAME functions that write the stand-alone files,
# called with save = FALSE so they return the ggplot without touching the disk.
# A panel and its place in the composite therefore can never drift apart: edit
# the panel function once and both the single file and the figure update.
#
# To re-arrange a figure, edit its `design` string. Each letter is an area of
# the grid and refers to a panel BY POSITION in `panels` (A = 1st, B = 2nd, ...),
# which is also the order the tags are assigned in -- so keep `panels` in
# reading order and the letters on the page will match the panel names.

.require_patchwork <- function() {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required to assemble the figures.\n",
         "Install it, or build the panels individually instead.", call. = FALSE)
  }
}

# Lay panels out on `design` and add the A/B/C/... tags.
.assemble <- function(panels, design, heights = NULL, tag_size = 20) {
  .require_patchwork()
  (patchwork::wrap_plots(panels, design = design, heights = heights) +
     patchwork::plot_annotation(tag_levels = "A")) &
    theme(plot.tag = element_text(size = tag_size, face = "plain", hjust = 0))
}

# --- Figure 3: cholesterol biosynthesis --------------------------------------
# A the GSEA curve and B the gene sets on top; C the GO terms and D the
# per-gene lollipop below.
figure3 <- function() {
  panels <- list(
    fig3a_cholesterol_gsea_curve(save = FALSE),
    fig3b_cholesterol_gsea_nes(save = FALSE),
    fig_cholesterol_go_up(save = FALSE),
    fig3d_wp_cholesterol_lollipop(save = FALSE))
  p <- .assemble(panels, design = "AB\nCD", heights = c(1, 1))
  save_figure(p, "figure3_cholesterol", width = 15, height = 10)
}

# --- Figure 4: the synaptic directional split --------------------------------
# The four panels of the figure proposal: up-regulated on top (A synaptic and
# axonal processes, B compartments), down-regulated below (C processes and
# compartments, D the paired enrichment curves). Laying A/B above C/D makes the
# directional asymmetry the visual argument of the figure.
#
# D holds two curves; here they are STACKED (ncol = 1) so each gets the full
# width of its half-page cell -- side by side in that space their titles and
# axis labels collide. The stand-alone fig4d file keeps them side by side.
#
# fig4e_translation_gsea_nes (the canonical ribosomal collapse) is deliberately
# NOT on this page: those sets are the negative half of panel 5C. It stays
# registered as a stand-alone panel -- add it back by appending it to `panels`
# and using the design "AABB\nCCDD\nEEEE".
figure4 <- function() {
  .require_patchwork()
  panels <- list(
    fig4a_go_bp_synaptic_up(save = FALSE),   # -> A
    fig4b_go_cc_up(save = FALSE),            # -> B
    fig4c_go_down(save = FALSE),             # -> C
    # wrap_elements makes patchwork treat the curve pair as ONE panel, so it
    # gets a single tag instead of one tag per curve.
    patchwork::wrap_elements(
      fig4d_synaptic_gsea_curves(ncol = 1, save = FALSE)))   # -> D
  p <- .assemble(panels, design = "AB\nCD", heights = c(1, 1.15))
  save_figure(p, "figure4_synaptic_split", width = 15, height = 11)
}

# --- Figure 6: WGCNA co-expression modules -----------------------------------
# A the module eigengenes, B the module-theme dot plot. The dot plot carries the
# argument -- each theme occupies its own module -- so it gets the larger share.
# The dendrogram (fig6a) is base graphics and cannot be composed here; it stays
# a stand-alone file to place manually or use as supplementary.
figure6 <- function() {
  panels <- list(
    fig6b_wgcna_eigengene_heatmap(save = FALSE), # -> A
    fig6c_wgcna_module_dotplot(save = FALSE))    # -> B
  p <- .assemble(panels, design = "A\nB", heights = c(0.75, 1.6))
  save_figure(p, "figure6_wgcna_modules", width = 9, height = 14)
}

# --- Figure 5: vascular, adhesion, and rank-based findings --------------------
# A vascular processes and B molecular functions on top; C the up-regulated
# canonical pathways and D the Hallmark sets below.
figure5 <- function() {
  panels <- list(
    fig5a_go_bp_vascular_up(save = FALSE),
    fig5b_go_mf_up(save = FALSE),
    fig5c_canonical_gsea_nes(save = FALSE),
    fig5d_hallmark_gsea_nes(save = FALSE))
  p <- .assemble(panels, design = "AB\nCD", heights = c(1, 1.35))
  save_figure(p, "figure5_vascular_adhesion", width = 15, height = 15)
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
  fig5d_hallmark_gsea_nes  = fig5d_hallmark_gsea_nes,

  # Figure 6 -- WGCNA co-expression modules
  fig6a_wgcna_dendrogram        = fig6a_wgcna_dendrogram,
  fig6b_wgcna_eigengene_heatmap = fig6b_wgcna_eigengene_heatmap,
  fig6c_wgcna_module_dotplot    = fig6c_wgcna_module_dotplot,

  # Assembled, panel-labelled figures (the manuscript-ready pages)
  figure3 = figure3,
  figure4 = figure4,
  figure5 = figure5,
  figure6 = figure6
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
