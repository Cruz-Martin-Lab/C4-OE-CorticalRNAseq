# =============================================================================
# 04_coexpression_wgcna.R
# -----------------------------------------------------------------------------
# Gene co-expression network analysis (WGCNA) on the VST-transformed matrix,
# reproducing the module-detection + per-module GO/GSEA figures from the
# original manuscript pipeline.
#
# Produces (all under Outputs/04_wgcna/):
#   figures/wgcna_soft_threshold.pdf
#   figures/wgcna_dendrogram.pdf
#   figures/wgcna_module_sizes.pdf
#   figures/wgcna_module_eigengene_heatmap.pdf            largest N modules
#   figures/wgcna_module_eigengene_heatmap_all_modules.pdf every module
#   figures/GO_ORA_module_<color>_barplot.pdf  per top module WITH significant terms
#   tables/wgcna_soft_threshold_scan.csv
#   tables/wgcna_module_assignments.csv
#   tables/wgcna_module_trait_correlation.csv
#   tables/GO_ORA_module_<color>.csv           for each of the top modules
#
# NOTE: none of these are currently on disk -- only the rdata cache survived an
# earlier folder copy. Re-run this script alone to regenerate them; it reads
# Outputs/rdata/02_deseq2_results.rds and does not need steps 01-03 re-run.
# =============================================================================

.get_script_dir <- function() {
  # Check the source() frame FIRST. commandArgs() reports the Rscript ENTRY
  # POINT, which is the wrong answer whenever this script is sourced by a
  # master runner (run_all.R, run_full_analysis.R) living in another folder.
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  # gsub("~+~"): R's front end encodes spaces in --file= as "~+~". This
  # repository's path contains spaces, so without decoding, every derived path
  # is wrong and dir.create(recursive=TRUE) silently builds a bogus tree.
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1) return(dirname(normalizePath(
    gsub("~+~", " ", sub("^--file=", "", file_arg), fixed = TRUE))))
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
source(file.path(SCRIPT_DIR, "00_config.R"))
source(file.path(SCRIPT_DIR, "R", "functions.R"))

DIR_TABLES  <- step_dir("04", "tables")
DIR_FIGURES <- step_dir("04", "figures")

deseq_cached <- readRDS(file.path(DIR_RDATA, "02_deseq2_results.rds"))
results_table <- deseq_cached$results_table
vst_matrix    <- deseq_cached$vst_matrix
norm_counts   <- deseq_cached$norm_counts
metadata      <- deseq_cached$metadata

# ---- 1. Prepare input matrix (samples in rows, genes in columns) -----------
WGCNA::allowWGCNAThreads(4)

# Restrict to adequately expressed genes. Roughly a third of the matrix is
# low-count material whose sample-to-sample variation is mostly sampling noise
# -- and WGCNA builds modules purely from correlation structure, so feeding
# those genes in creates modules of co-fluctuating noise that compete with real
# biology and distort the soft-threshold scan.
#
# The cutoff is the MEAN NORMALIZED COUNT (WGCNA_MIN_MEAN_COUNT in 00_config.R),
# computed from the normalized count matrix alone. It is deliberately NOT the
# "genes DESeq2 could test" set: filtering WGCNA input by anything derived from
# the differential expression model is discouraged, because the network should
# be built without reference to the trait. This criterion never sees genotype,
# p-values, or any model output, so the modules below are genuinely unsupervised.
count_mat <- as.matrix(norm_counts[, setdiff(names(norm_counts), "gene")])
rownames(count_mat) <- norm_counts$gene
mean_norm_count <- rowMeans(count_mat)
expressed_genes <- names(mean_norm_count)[mean_norm_count >= WGCNA_MIN_MEAN_COUNT]

vst_filtered <- vst_matrix[rownames(vst_matrix) %in% expressed_genes, , drop = FALSE]
message(sprintf("WGCNA input: %d of %d genes retained (mean normalized count >= %g).",
                 nrow(vst_filtered), nrow(vst_matrix), WGCNA_MIN_MEAN_COUNT))

wgcna_input <- prepare_wgcna_matrix(vst_filtered)

good_samples_genes <- WGCNA::goodSamplesGenes(wgcna_input, verbose = 0)
if (!good_samples_genes$allOK) {
  wgcna_input <- wgcna_input[good_samples_genes$goodSamples, good_samples_genes$goodGenes]
  message("Removed genes/samples flagged by WGCNA::goodSamplesGenes().")
}

# ---- 2. Soft-threshold power selection --------------------------------------
sft <- run_soft_threshold_scan(wgcna_input, powers = 1:20)

# Soft-thresholding power. Uses WGCNA_POWER from 00_config.R when set;
# otherwise auto-selects the first power whose SIGNED R^2 (-sign(slope) *
# SFT.R.sq, the same quantity plotted in wgcna_soft_threshold.pdf) reaches
# WGCNA_RSQ_CUTOFF. Using the signed value matters: at low powers the raw
# SFT.R.sq can be high while the slope is POSITIVE, which is the opposite of
# scale-free topology -- selecting on raw R^2 alone can therefore pick power 1
# or 2 on some datasets.
#
# If no power reaches the cutoff, the fallback follows Langfelder & Horvath's
# sample-size recommendation rather than an arbitrary constant.
sft_fit <- sft$fitIndices
sft_fit$signed_rsq <- -sign(sft_fit$slope) * sft_fit$SFT.R.sq

wgcna_default_power <- function(n_samples, signed = TRUE) {
  if (signed) {
    if (n_samples < 20) 18 else if (n_samples < 30) 16 else if (n_samples < 40) 14 else 12
  } else {
    if (n_samples < 20) 9  else if (n_samples < 30) 8  else if (n_samples < 40) 7  else 6
  }
}

picked_power <- if (!is.null(WGCNA_POWER)) {
  message(sprintf("WGCNA power: using WGCNA_POWER = %d from 00_config.R.", WGCNA_POWER))
  WGCNA_POWER
} else {
  candidates <- sft_fit$Power[sft_fit$signed_rsq >= WGCNA_RSQ_CUTOFF]
  if (length(candidates) > 0) {
    p <- min(candidates)
    message(sprintf("WGCNA power: %d (first power with signed R^2 >= %.2f).", p, WGCNA_RSQ_CUTOFF))
    p
  } else {
    p <- wgcna_default_power(nrow(wgcna_input), signed = TRUE)
    message(sprintf(
      "WGCNA power: no power reached signed R^2 >= %.2f (max observed %.3f at power %d). Falling back to the WGCNA sample-size recommendation for a signed network with %d samples: %d.",
      WGCNA_RSQ_CUTOFF, max(sft_fit$signed_rsq, na.rm = TRUE),
      sft_fit$Power[which.max(sft_fit$signed_rsq)], nrow(wgcna_input), p))
    p
  }
}

# Report where the chosen power actually sits on the curve, for the record.
row_at <- sft_fit[sft_fit$Power == picked_power, ]
if (nrow(row_at) == 1) {
  message(sprintf("  -> at power %d: signed R^2 = %.3f, mean connectivity = %.1f",
                   picked_power, row_at$signed_rsq, row_at$mean.k.))
}
write.csv(sft_fit, file.path(DIR_TABLES, "wgcna_soft_threshold_scan.csv"), row.names = FALSE)

# Plot the scan with the chosen power marked (drawn after selection so the
# figure documents which power was actually used).
sft_plots <- plot_soft_threshold(sft, chosen_power = picked_power)
pdf(file.path(DIR_FIGURES, "wgcna_soft_threshold.pdf"), width = 10, height = 5)
gridExtra::grid.arrange(sft_plots$scale_independence, sft_plots$mean_connectivity, ncol = 2)
dev.off()

# ---- 3. Module detection -----------------------------------------------------
# Module resolution comes from 00_config.R (WGCNA_MIN_MODULE_SIZE,
# WGCNA_MERGE_CUT_HEIGHT).
net <- run_wgcna_modules(wgcna_input, power = picked_power,
                          min_module_size = WGCNA_MIN_MODULE_SIZE,
                          merge_cut_height = WGCNA_MERGE_CUT_HEIGHT)
message(sprintf("  -> minModuleSize = %d, mergeCutHeight = %.2f",
                 WGCNA_MIN_MODULE_SIZE, WGCNA_MERGE_CUT_HEIGHT))
module_colors <- WGCNA::labels2colors(net$colors)

pdf(file.path(DIR_FIGURES, "wgcna_dendrogram.pdf"), width = 9, height = 5)
WGCNA::plotDendroAndColors(
  net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
  "Module colors", dendroLabels = FALSE, hang = 0.05,
  addGuide = TRUE, guideHang = 0.05
)
dev.off()

module_df <- data.frame(gene_id = names(net$colors), module = module_colors,
                          stringsAsFactors = FALSE)
write.csv(module_df, file.path(DIR_TABLES, "wgcna_module_assignments.csv"), row.names = FALSE)

module_sizes <- as.data.frame(table(module_df$module))
names(module_sizes) <- c("Module", "Gene_Count")
module_sizes <- module_sizes[order(-module_sizes$Gene_Count), ]

ggplot2::ggsave(
  file.path(DIR_FIGURES, "wgcna_module_sizes.pdf"),
  ggplot2::ggplot(head(module_sizes, 20),
                   ggplot2::aes(x = reorder(Module, -Gene_Count), y = Gene_Count, fill = Module)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(title = "Module sizes", x = "Module", y = "Gene count") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)),
  width = 8, height = 5
)

# ---- 4. Module eigengenes + heatmap ------------------------------------------
me_list <- WGCNA::moduleEigengenes(wgcna_input, colors = module_colors)
MEs <- WGCNA::orderMEs(me_list$eigengenes)

# Apply the Control / C4-OE display labels (the eigengene matrix inherits the
# raw raw_counts.csv sample names via wgcna_input).
sample_ids <- rownames(MEs)          # raw sample IDs, as in raw_counts.csv
me_mat <- t(scale(MEs))
colnames(me_mat) <- display_sample(sample_ids)

annotation_col <- data.frame(
  Genotype = display_genotype(metadata[sample_ids, "genotype"]),
  Sex      = metadata[sample_ids, "sex"]
)
rownames(annotation_col) <- display_sample(sample_ids)

me_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100)

# The complete heatmap, written for the record so no module is hidden.
pdf(file.path(DIR_FIGURES, "wgcna_module_eigengene_heatmap_all_modules.pdf"),
    width = 8, height = max(6, 0.16 * nrow(me_mat) + 3))
pheatmap::pheatmap(me_mat,
                    annotation_col = annotation_col,
                    color = me_palette,
                    main = sprintf("Module eigengenes (z-scored)\nall %d modules",
                                    nrow(me_mat)))
dev.off()

# The manuscript panel: the largest WGCNA_HEATMAP_N_MODULES modules only. At
# this sample size the small modules are the least robust, and plotting every
# one makes the panel unreadable. Grey is the unassigned-gene bin, not a
# module, so it is excluded. The title states what is shown out of what.
modules_ranked <- setdiff(as.character(module_sizes$Module), "grey")
me_top         <- paste0("ME", head(modules_ranked, WGCNA_HEATMAP_N_MODULES))
me_top         <- me_top[me_top %in% rownames(me_mat)]
me_mat_top     <- me_mat[me_top, , drop = FALSE]

pdf(file.path(DIR_FIGURES, "wgcna_module_eigengene_heatmap.pdf"),
    width = 8, height = max(4.5, 0.42 * nrow(me_mat_top) + 2.5))
pheatmap::pheatmap(me_mat_top,
                    annotation_col = annotation_col,
                    color = me_palette,
                    main = sprintf("Module eigengenes (z-scored)\n%d largest of %d modules (grey excluded)",
                                    nrow(me_mat_top), length(modules_ranked)))
dev.off()
message(sprintf("  -> eigengene heatmap: %d of %d modules shown (grey excluded)",
                 nrow(me_mat_top), length(modules_ranked)))

# ---- 4b. Module-trait relationships -----------------------------------------
# The central WGCNA output for a designed experiment: correlate each module
# eigengene with genotype, and (given the confound in this dataset) with sex.
#
# IMPORTANT POWER CAVEAT: with 8 samples, each correlation is estimated from
# 8 points. |r| must exceed ~0.71 for a nominal p < 0.05, and with dozens of
# modules tested the FDR correction is severe. We treat these as descriptive /
# hypothesis-generating, not as confirmatory statistics.
trait_df <- data.frame(
  genotype = as.numeric(metadata[sample_ids, "genotype"] == GENOTYPE_LEVELS[2]),
  sex      = as.numeric(metadata[sample_ids, "sex"] == "M")
)

mt_cor  <- WGCNA::cor(MEs, trait_df, use = "p")
mt_p    <- WGCNA::corPvalueStudent(mt_cor, nrow(MEs))
module_trait <- data.frame(
  module        = sub("^ME", "", rownames(mt_cor)),
  size          = as.integer(table(module_colors)[sub("^ME", "", rownames(mt_cor))]),
  cor_genotype  = round(mt_cor[, "genotype"], 3),
  p_genotype    = signif(mt_p[, "genotype"], 3),
  padj_genotype = signif(p.adjust(mt_p[, "genotype"], method = "BH"), 3),
  cor_sex       = round(mt_cor[, "sex"], 3),
  p_sex         = signif(mt_p[, "sex"], 3),
  stringsAsFactors = FALSE
)
module_trait <- module_trait[order(module_trait$p_genotype), ]
write.csv(module_trait, file.path(DIR_TABLES, "wgcna_module_trait_correlation.csv"),
          row.names = FALSE)

n_sig_raw  <- sum(module_trait$p_genotype < 0.05, na.rm = TRUE)
n_sig_adj  <- sum(module_trait$padj_genotype < 0.05, na.rm = TRUE)
message(sprintf(
  "Module-trait: %d of %d modules correlate with genotype at nominal p < 0.05; %d survive BH correction.",
  n_sig_raw, nrow(module_trait), n_sig_adj))

# ---- 5. Per-module GO ORA for the largest few modules (excluding "grey" =
#         unassigned genes), using the SAME universe fix as script 03 -------
all_tested_genes <- results_table %>% dplyr::filter(!is.na(padj)) %>% dplyr::pull(gene)
# Profile the N largest real modules. "grey" is excluded first rather than
# subtracted afterwards, so excluding it never silently reduces the count.
real_modules <- setdiff(as.character(module_sizes$Module), "grey")
n_top <- if (!exists("WGCNA_N_TOP_MODULES") || is.null(WGCNA_N_TOP_MODULES)) {
  length(real_modules)          # NULL => annotate every module
} else {
  WGCNA_N_TOP_MODULES
}
top_modules <- head(real_modules, n_top)
message(sprintf("Annotating the %d largest modules: %s",
                 length(top_modules), paste(top_modules, collapse = ", ")))

for (mod in top_modules) {
  mod_genes <- module_df$gene_id[module_df$module == mod]
  go_res <- run_go_ora(mod_genes, universe = all_tested_genes, ont = "BP")
  go_df <- as.data.frame(go_res)
  write.csv(go_df, file.path(DIR_TABLES, sprintf("GO_ORA_module_%s.csv", mod)), row.names = FALSE)

  if (nrow(go_df) > 0) {
    mod_plot <- make_go_barplot(go_res, num_paths = 20,
                                 title = sprintf("Module '%s' GO BP enrichment", mod))
    if (!is.null(mod_plot)) {
      ggplot2::ggsave(
        file.path(DIR_FIGURES, sprintf("GO_ORA_module_%s_barplot.pdf", mod)),
        mod_plot, width = 8, height = 7
      )
    } else {
      message(sprintf("Module '%s': no significant GO BP terms; barplot skipped.", mod))
    }
  }
}

# ---- 6. Cache -----------------------------------------------------------------
saveRDS(list(net = net, module_df = module_df, MEs = MEs, picked_power = picked_power),
        file.path(DIR_RDATA, "04_wgcna_results.rds"))


# ---- 6. Supplementary tables + module-theme dot plot -------------------------
# One row per significant GO term per module, so the per-module CSVs do not have
# to be read separately.
go_all <- do.call(rbind, lapply(top_modules, function(mod) {
  f <- file.path(DIR_TABLES, sprintf("GO_ORA_module_%s.csv", mod))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f)
  d <- d[!is.na(d$p.adjust) & d$p.adjust < PADJ_THRESHOLD, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  data.frame(module = mod, d[, c("ID", "Description", "GeneRatio", "Count",
                                  "pvalue", "p.adjust", "geneID")],
             stringsAsFactors = FALSE)
}))
go_all <- go_all[order(go_all$module, go_all$p.adjust), ]
write.csv(go_all, file.path(DIR_TABLES, "GO_ORA_modules_combined.csv"), row.names = FALSE)
message(sprintf("  -> combined GO table: %d significant terms across %d modules",
                 nrow(go_all), length(unique(go_all$module))))

# One row per module: size, genotype correlation, and GO summary. 04b adds a
# `stability` column to this file when it runs.
n_sig <- vapply(module_sizes$Module, function(m) sum(go_all$module == m), 0L)
top_terms <- vapply(module_sizes$Module, function(m) {
  d <- go_all[go_all$module == m, ]
  if (!nrow(d)) "" else paste(utils::head(d$Description[order(d$p.adjust)], 3), collapse = "; ")
}, "")
mod_summary <- data.frame(
  module = as.character(module_sizes$Module),
  genes  = as.integer(module_sizes$Gene_Count),
  n_sig_GO = as.integer(n_sig), top_GO_terms = top_terms, stringsAsFactors = FALSE)
mt_join <- module_trait[, intersect(names(module_trait),
             c("module","cor_genotype","p_genotype","padj_genotype","cor_sex","p_sex"))]
mod_summary <- merge(mod_summary, mt_join, by = "module", all.x = TRUE)
mod_summary <- mod_summary[order(-mod_summary$genes), ]
write.csv(mod_summary, file.path(DIR_TABLES, "wgcna_module_summary.csv"), row.names = FALSE)

# Dot plot: the modules worth naming, each with its strongest GO terms. A term
# is drawn in every displayed module where it is significant, so the panel also
# shows how module-specific each theme is.
if (nrow(go_all)) {
  by_count <- names(sort(tapply(go_all$module, go_all$module, length), decreasing = TRUE))
  always   <- mod_summary$module[!is.na(mod_summary$padj_genotype) &
                                   mod_summary$padj_genotype < PADJ_THRESHOLD]
  show_mods <- unique(c(utils::head(by_count, WGCNA_DOTPLOT_N_MODULES),
                        intersect(always, unique(go_all$module))))
  show_mods <- by_count[by_count %in% show_mods]          # keep the ranking order
  terms <- unique(unlist(lapply(show_mods, function(m) {
    d <- go_all[go_all$module == m, ]
    utils::head(d$Description[order(d$p.adjust)], WGCNA_DOTPLOT_N_TERMS)
  })))
  dp <- go_all[go_all$module %in% show_mods & go_all$Description %in% terms, ]
  # Wrap long term names so one verbose GO description cannot squeeze the panel.
  wrap_term <- function(x) vapply(x, function(z)
    paste(strwrap(z, width = 42), collapse = "\n"), "")
  terms_w <- wrap_term(terms)
  dp$module      <- factor(dp$module, levels = show_mods)
  dp$Description <- factor(wrap_term(as.character(dp$Description)), levels = rev(terms_w))
  p_dot <- ggplot2::ggplot(dp, ggplot2::aes(x = module, y = Description,
                                             size = Count, colour = -log10(p.adjust))) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient(low = "steelblue", high = "firebrick",
                                    name = expression(-log[10] ~ "p.adjust")) +
    ggplot2::scale_size_continuous(name = "Genes", range = c(1.5, 6)) +
    ggplot2::labs(title = "Module GO biological process enrichment",
                   x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                    panel.grid.major = ggplot2::element_line(colour = "grey92"))
  ggplot2::ggsave(file.path(DIR_FIGURES, "wgcna_module_theme_dotplot.pdf"), p_dot,
                   width = 8.5, height = max(5, 0.3 * length(terms) + 2.2), limitsize = FALSE)
  message(sprintf("  -> dot plot: %d modules x %d terms", length(show_mods), length(terms)))
}

message("04_coexpression_wgcna.R complete.")
