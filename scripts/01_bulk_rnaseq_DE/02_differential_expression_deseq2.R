# =============================================================================
# 02_differential_expression_deseq2.R
# -----------------------------------------------------------------------------
# DESeq2 differential expression: design = ~ sex + genotype.
#
# Produces (all under Outputs/02_differential_expression/):
#   tables/S1_Table_DESeq.csv        the DEG table
#   figures/volcano_deseq2.pdf
#   figures/volcano_deseq2_by_category.pdf   (manuscript volcano)
#   figures/pca_deseq2.pdf
#   figures/dispersion_deseq2.pdf
#   figures/logfc_histogram_deseq2.pdf
#   figures/heatmap_top_degs.pdf
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

DIR_TABLES  <- step_dir("02", "tables")
DIR_FIGURES <- step_dir("02", "figures")

cached <- readRDS(file.path(DIR_RDATA, "01_loaded_data.rds"))
symbol_counts <- cached$symbol_counts
metadata      <- cached$metadata

# ---- 1. Run DESeq2 with sex covariate ---------------------------------------
deseq_out <- run_deseq2(symbol_counts, metadata)

# ---- 2. Flag (don't drop) sex-confounded genes and sex-linked genes ---------
results_table <- deseq_out$results_table
results_table <- flag_sex_confounded_genes(symbol_counts, metadata, results_table)
results_table <- flag_sex_linked_genes(results_table, gene_col = "gene")

n_sexlinked_sig <- results_table %>%
  dplyr::filter(sex_linked, !is.na(padj), padj < PADJ_THRESHOLD) %>% nrow()
message(sprintf(
  "%d sex-linked genes are flagged in the DE table; %d of those are significant at padj<%.2f.\n",
  sum(results_table$sex_linked), n_sexlinked_sig, PADJ_THRESHOLD),
  "They remain in the full table (flagged via `sex_linked` column) for transparency, ",
  "and are excluded from headline gene-set/GO input lists in 03_go_enrichment (see that script)."
)

# ---- 3. Save the DEG table ---------------------------------------------------
results_table <- results_table %>% dplyr::arrange(padj)
write.csv(results_table, file.path(DIR_TABLES, "S1_Table_DESeq.csv"),
          row.names = FALSE)

# ---- 4. Manuscript figures ----------------------------------------------------
ggplot2::ggsave(file.path(DIR_FIGURES, "logfc_histogram_deseq2.pdf"),
                 plot_logfc_hist(results_table, PADJ_THRESHOLD, "log2 fold change distribution"),
                 width = 6, height = 4)

pdf(file.path(DIR_FIGURES, "dispersion_deseq2.pdf"), width = 6, height = 5)
deseq_out$dispersion_plot()
dev.off()

ggplot2::ggsave(file.path(DIR_FIGURES, "pca_deseq2.pdf"),
                 plot_pca(deseq_out$vst_matrix, metadata,
                          color_by = "genotype_label", shape_by = "sex",
                          title = "PCA (VST-transformed counts)"),
                 width = 6, height = 5)

ggplot2::ggsave(file.path(DIR_FIGURES, "volcano_deseq2.pdf"),
                 plot_volcano(results_table, PADJ_THRESHOLD, LOGFC_THRESHOLD,
                              CONTRAST_LABEL, highlight_sex_confound = FALSE),
                 width = 7, height = 6)

# Category volcano: significant DEGs colored red (up) / blue (down) by padj
# alone, with labels chosen from the curated GO categories (lipid metabolism /
# synapse / immune) plus the genes named in force_label_genes. This is the
# manuscript volcano figure.
#
# NOTE: the older EnhancedVolcano-based figure was removed from the pipeline.
# EnhancedVolcano gates both point color AND label selection on |logFC| >
# FCcutoff as well as padj, which excludes the coordinated modest-effect
# biology in this dataset (e.g. every cholesterol-biosynthesis gene sits
# below |logFC| = 1). plot_enhanced_volcano() is still defined in
# functions.R if it is ever needed, but it is no longer called.
ggplot2::ggsave(
  file.path(DIR_FIGURES, "volcano_deseq2_by_category.pdf"),
  plot_category_volcano(results_table, PADJ_THRESHOLD, LOGFC_THRESHOLD,
                         title = paste0(CONTRAST_LABEL, ": DEGs by category")),
  width = 9, height = 7
)

heatmap_plot <- plot_heatmap_top_deg(results_table, deseq_out$norm_counts, metadata,
                                      top_n = 50, title = "Top 50 DEGs")
pdf(file.path(DIR_FIGURES, "heatmap_top_degs.pdf"), width = 7, height = 9)
grid::grid.newpage()
grid::grid.draw(heatmap_plot$gtable)
dev.off()

# ---- 5. Cache for downstream scripts -----------------------------------------
saveRDS(list(results_table = results_table, vst_matrix = deseq_out$vst_matrix,
             norm_counts = deseq_out$norm_counts, metadata = metadata),
        file.path(DIR_RDATA, "02_deseq2_results.rds"))

message("02_differential_expression_deseq2.R complete. DEG table written to ",
        file.path(DIR_TABLES, "S1_Table_DESeq.csv"))
