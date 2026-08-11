# =============================================================================
# 03_go_enrichment_ORA_GSEA.R
# -----------------------------------------------------------------------------
# Gene Ontology / pathway enrichment on the differential-expression table.
#
# Key methodological choices:
#   - enrichGO() receives `universe = all_tested_genes` (every gene that
#     survived filtering and was actually tested by DESeq2), instead of
#     defaulting to every annotated gene in org.Mm.eg.db. Using the whole
#     genome as background inflates enrichment significance for any gene set
#     correlated with "detectable in bulk cortical RNA-seq" (e.g. broadly
#     expressed housekeeping-adjacent categories).
#   - Sex-linked genes are excluded from the significant-gene input list (but
#     NOT from the universe) so a handful of sex-chromosome genes can't drive
#     GO terms in the headline enrichment; see SEX_LINKED_GENES in 00_config.R.
#
# Produces (all under Outputs/03_enrichment/), per ontology (BP/MF/CC)
# x direction (combined/up/down):
#   tables/GO_ORA_<ont>_<direction>.csv
#   figures/GO_ORA_<ont>_<direction>_barplot.pdf   only where terms are significant
#   tables/GO_ORA_BP_synaptic_<direction>.csv
#   figures/GO_ORA_BP_synaptic_top10_<direction>.pdf
#
# ORA is run three ways because up- and down-regulated genes frequently
# enrich for different biology; combining them can dilute both signals.
# The universe is the same for all three (every tested gene).
#
# And per MSigDB collection:
#   tables/GSEA_<collection>.csv
#   figures/GSEA_<collection>_barplot.pdf
#   tables/GSEA_theme_<theme>.csv           filtered views, no figure by design
#   figures/lollipop_WP_cholesterol_biosynthesis.pdf
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

DIR_TABLES  <- step_dir("03", "tables")
DIR_FIGURES <- step_dir("03", "figures")

cached <- readRDS(file.path(DIR_RDATA, "02_deseq2_results.rds"))
results_table <- cached$results_table

# ---- 1. Define the gene universe and the significant-gene input list --------
all_tested_genes <- results_table %>% dplyr::filter(!is.na(padj)) %>% dplyr::pull(gene)

# NOTE on the fold-change cutoff. The DEG list fed to ORA is defined by
# adjusted p-value ONLY (padj < PADJ_THRESHOLD), deliberately WITHOUT an
# additional |logFC| filter.
#
# Requiring |logFC| > 1 as well shrinks the input from 1,573 significant
# genes to 189, and it does so in a biologically biased way: it removes
# precisely the coordinated, modest-magnitude pathway responses that
# enrichment analysis exists to detect. In this dataset all 14 significant
# cholesterol-biosynthesis genes (Hmgcs1 0.64, Hmgcr 0.66, Sqle 0.36,
# Cyp51 0.58, Dhcr7 0.34, Msmo1 0.44, ...) are significant at padj < 0.05 but
# NONE reach |logFC| > 1, so the fold-change filter silently made the GO ORA
# blind to the manuscript's headline finding. This is the same failure mode
# that the fgsea ranking metric had (see make_ranked_list() in functions.R).
#
# LOGFC_THRESHOLD is still used for volcano-plot reference lines; it is not
# used to define enrichment input. Set ORA_USE_LOGFC_FILTER <- TRUE in
# 00_config.R only if you specifically want the old, stricter behaviour.
use_lfc_filter <- if (exists("ORA_USE_LOGFC_FILTER")) ORA_USE_LOGFC_FILTER else FALSE

deg_table <- results_table %>%
  dplyr::filter(!is.na(padj), padj < PADJ_THRESHOLD,
                if (use_lfc_filter) abs(logFC) > LOGFC_THRESHOLD else TRUE,
                !sex_linked, !grepl("^Gm[0-9]", gene), !grepl("Rik$", gene))

# ORA is run three ways. Up- and down-regulated genes frequently enrich for
# entirely different biology, and combining them can dilute both signals to
# non-significance (a term enriched only among down-regulated genes is
# competing against all the up-regulated genes as well). Running the
# directions separately and the combined list alongside them makes it
# explicit which terms are direction-specific and which hold overall.
#
# The universe is IDENTICAL for all three runs -- it is every gene that was
# tested, not just the tested genes of that direction. The background is the
# set of genes that could have been called, which does not change depending
# on which subset is being interrogated.
gene_sets_by_direction <- list(
  combined = deg_table$gene,
  up       = deg_table %>% dplyr::filter(logFC > 0) %>% dplyr::pull(gene),
  down     = deg_table %>% dplyr::filter(logFC < 0) %>% dplyr::pull(gene)
)

message(sprintf(
  "GO universe: %d tested genes. DEG input (padj < %.2f%s; sex-linked and Gm*/Rik excluded): combined %d, up %d, down %d.",
  length(all_tested_genes), PADJ_THRESHOLD,
  if (use_lfc_filter) sprintf(" AND |logFC| > %.1f", LOGFC_THRESHOLD) else "",
  length(gene_sets_by_direction$combined),
  length(gene_sets_by_direction$up),
  length(gene_sets_by_direction$down)))

for (dir_name in names(gene_sets_by_direction)) {
  if (length(gene_sets_by_direction[[dir_name]]) < 5) {
    warning(sprintf("Only %d gene(s) in the '%s' set -- GO ORA will be unstable or empty.",
                     length(gene_sets_by_direction[[dir_name]]), dir_name))
  }
}

# ---- 2. GO ORA: 3 ontologies x 3 directions ---------------------------------
# Results are stored as go_results[[ontology]][[direction]] and written to
# Outputs/tables/GO_ORA_<ont>_<direction>.csv with a matching barplot.
go_results <- list()
for (ont in c("BP", "MF", "CC")) {
  go_results[[ont]] <- list()
  for (dir_name in names(gene_sets_by_direction)) {
    genes_in <- gene_sets_by_direction[[dir_name]]
    if (length(genes_in) == 0) {
      message(sprintf("GO %s (%s): no input genes; skipped.", ont, dir_name))
      next
    }

    go_res <- run_go_ora(genes_in, universe = all_tested_genes, ont = ont)
    go_results[[ont]][[dir_name]] <- go_res

    go_df <- as.data.frame(go_res)
    write.csv(go_df,
              file.path(DIR_TABLES, sprintf("GO_ORA_%s_%s.csv", ont, dir_name)),
              row.names = FALSE)

    n_sig <- sum(go_df$p.adjust < PADJ_THRESHOLD, na.rm = TRUE)
    message(sprintf("GO %s (%s): %d terms returned, %d significant at padj < %.2f.",
                     ont, dir_name, nrow(go_df), n_sig, PADJ_THRESHOLD))

    if (nrow(go_df) > 0) {
      dir_label <- c(combined = "all DEGs", up = "upregulated", down = "downregulated")[[dir_name]]
      go_plot <- make_go_barplot(go_res, num_paths = 20,
                                  title = sprintf("GO %s Enrichment (%s)", ont, dir_label))
      if (!is.null(go_plot)) {
        ggplot2::ggsave(
          file.path(DIR_FIGURES, sprintf("GO_ORA_%s_%s_barplot.pdf", ont, dir_name)),
          go_plot, width = 8, height = 7
        )
      } else {
        message(sprintf("GO %s (%s): no significant terms; barplot skipped.", ont, dir_name))
      }
    }
  }
}

# ---- 2b. Focused synaptic GO BP barplots -------------------------------------
# Top 10 most significant synaptic Biological Process terms. "Synaptic" is
# defined by GO hierarchy (SYNAPTIC_GO_BP_ROOTS in 00_config.R plus all their
# GO descendants) rather than by keyword-matching term names, so child terms
# whose names don't contain "synapse" are still captured.
# Produced for all three directions, since synaptic terms enriched only among
# downregulated genes are a different claim from ones enriched overall.
synaptic_ids <- get_synaptic_go_ids()
message(sprintf("Synaptic GO BP definition resolved to %d GO IDs.", length(synaptic_ids)))

for (dir_name in names(gene_sets_by_direction)) {
  go_bp <- go_results[["BP"]][[dir_name]]
  if (is.null(go_bp) || nrow(as.data.frame(go_bp)) == 0) next

  synaptic_df <- as.data.frame(go_bp) %>%
    dplyr::filter(ID %in% synaptic_ids) %>%
    dplyr::arrange(p.adjust)
  write.csv(synaptic_df,
            file.path(DIR_TABLES, sprintf("GO_ORA_BP_synaptic_%s.csv", dir_name)),
            row.names = FALSE)

  dir_label <- c(combined = "all DEGs", up = "upregulated", down = "downregulated")[[dir_name]]
  synaptic_plot <- plot_go_subset_barplot(
    go_bp, keep_ids = synaptic_ids, num_paths = 10,
    title = sprintf("Top 10 synaptic GO BP terms (%s)", dir_label)
  )
  if (!is.null(synaptic_plot)) {
    ggplot2::ggsave(file.path(DIR_FIGURES, sprintf("GO_ORA_BP_synaptic_top10_%s.pdf", dir_name)),
                     synaptic_plot, width = 8, height = 5)
  }
}

# ---- 3. GSEA across all MSigDB collections -----------------------------------
ranked_list <- make_ranked_list(results_table %>% dplyr::filter(!sex_linked))

gsea_results <- list()
for (collection_name in names(GMT_FILES)) {
  gmt_path <- GMT_FILES[[collection_name]]
  if (!file.exists(gmt_path)) {
    warning("GMT file not found, skipping: ", gmt_path)
    next
  }
  gsea_res <- run_fgsea(gmt_path, ranked_list)
  gsea_results[[collection_name]] <- gsea_res

  # fgsea's `leadingEdge` column is a list-column; flatten for CSV export.
  gsea_export <- gsea_res
  gsea_export$leadingEdge <- vapply(gsea_export$leadingEdge, paste, collapse = ",", FUN.VALUE = character(1))
  write.csv(gsea_export, file.path(DIR_TABLES, sprintf("GSEA_%s.csv", collection_name)),
            row.names = FALSE)

  gsea_plot <- plot_gsea_top(gsea_res, num_paths = 15,
                              title = sprintf("GSEA: %s", collection_name))
  if (!is.null(gsea_plot)) {
    ggplot2::ggsave(
      file.path(DIR_FIGURES, sprintf("GSEA_%s_barplot.pdf", collection_name)),
      gsea_plot, width = 9, height = 7
    )
  } else {
    message(sprintf("GSEA (%s): no gene sets significant at padj<=%.2f.", collection_name, PADJ_THRESHOLD))
  }
}

# ---- 3b. Themed GSEA summaries across all collections ------------------------
# Pulls every gene set matching a biological theme into one table, so a
# specific pathway question can be answered without searching six CSVs.
# Reports all matching sets with their NES and padj, significant or not --
# a clearly non-significant result for a pathway of interest is itself an
# answer worth recording.
gsea_themes <- list(
  cholesterol_lipid = "CHOLESTEROL|STEROL|LIPID|FATTY_ACID|MEVALONATE|SREBP|LIPOPROTEIN",
  synaptic          = "SYNAP|NEUROTRANSMITTER|DENDRIT|AXON",
  immune            = "IMMUN|INFLAMMAT|INTERFERON|COMPLEMENT|CYTOKINE"
)

for (theme_name in names(gsea_themes)) {
  theme_tbl <- collect_gsea_by_keyword(gsea_results, gsea_themes[[theme_name]])
  if (is.null(theme_tbl) || nrow(theme_tbl) == 0) {
    message(sprintf("GSEA theme '%s': no matching gene sets found.", theme_name))
    next
  }
  theme_tbl$leadingEdge <- vapply(theme_tbl$leadingEdge, paste, collapse = ",",
                                   FUN.VALUE = character(1))
  write.csv(theme_tbl, file.path(DIR_TABLES, sprintf("GSEA_theme_%s.csv", theme_name)),
            row.names = FALSE)
  n_sig <- sum(theme_tbl$padj <= PADJ_THRESHOLD, na.rm = TRUE)
  message(sprintf("GSEA theme '%s': %d matching gene sets, %d significant at padj<=%.2f.",
                   theme_name, nrow(theme_tbl), n_sig, PADJ_THRESHOLD))
}

# ---- 3c. Per-gene lollipop plot for the WP cholesterol biosynthesis set ------
# Shows the individual fold changes behind the pathway-level result. The gene
# set is read from the MSigDB canonical-pathways GMT by name, so it stays tied
# to the same MSigDB version as the GSEA above rather than being hard-coded.
# Significance is padj-only (no fold-change cutoff) -- these genes are exactly
# the coordinated, modest-magnitude changes that a |logFC| filter would hide.
chol_set <- get_gmt_pathway_genes(GMT_FILES$cp_only, "WP_CHOLESTEROL_BIOSYNTHESIS")
if (is.null(chol_set)) {
  chol_set <- get_gmt_pathway_genes(GMT_FILES$canonical_paths, "WP_CHOLESTEROL_BIOSYNTHESIS")
}

if (!is.null(chol_set)) {
  message("Lollipop gene set: ", chol_set$name, " (", length(chol_set$genes), " genes)")
  chol_lollipop <- plot_pathway_lollipop(
    results_table, chol_set$genes,
    title = "WP Cholesterol Biosynthesis Genes"
  )
  if (!is.null(chol_lollipop)) {
    ggplot2::ggsave(file.path(DIR_FIGURES, "lollipop_WP_cholesterol_biosynthesis.pdf"),
                     chol_lollipop, width = 6, height = 5)
    write.csv(
      results_table %>%
        dplyr::filter(gene %in% chol_set$genes) %>%
        dplyr::arrange(dplyr::desc(logFC)),
      file.path(DIR_TABLES, "WP_cholesterol_biosynthesis_genes.csv"), row.names = FALSE
    )
  }
} else {
  message("WP cholesterol biosynthesis set not found in the available GMT files; ",
          "lollipop plot skipped.")
}

# ---- 4. Cache for downstream scripts -----------------------------------------
saveRDS(list(go_results = go_results, gsea_results = gsea_results,
             gene_sets_by_direction = gene_sets_by_direction,
             all_tested_genes = all_tested_genes),
        file.path(DIR_RDATA, "03_go_gsea_results.rds"))

message("03_go_enrichment_ORA_GSEA.R complete.")
