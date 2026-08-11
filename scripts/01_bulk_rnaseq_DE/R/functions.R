# =============================================================================
# functions.R
# -----------------------------------------------------------------------------
# Shared helper functions for the C4-OE cortical bulk RNA-seq analysis.
# Key methodological choices these functions implement:
#   - DESeq2 design formula `~ sex + genotype`
#   - Ensembl -> gene symbol mapping via the local org.Mm.eg.db annotation
#     package (offline, reproducible, and immune to Ensembl mirror downtime)
#   - clusterProfiler::enrichGO() takes an explicit `universe` argument
#     (the set of genes actually tested for DE)
#   - sex-linked gene handling is explicit and auditable rather than a silent
#     row filter with no confirmation step
#
# Every function is pure (input -> output, no reliance on global state other
# than the constants defined in 00_config.R) so the pipeline can be run
# top-to-bottom or function-by-function during review.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(fgsea)
  library(EnhancedVolcano)
  library(pheatmap)
  library(RColorBrewer)
  library(WGCNA)
})

# -----------------------------------------------------------------------------
# 1. Data import & preprocessing
# -----------------------------------------------------------------------------

#' Read a raw counts CSV (handles CR-only line endings from older macOS exports)
read_counts <- function(file) {
  raw_lines <- readLines(file, warn = FALSE)
  # readLines already splits on \r, \n, or \r\n, so this is robust regardless
  # of which line-ending convention the exporting machine used.
  con <- textConnection(raw_lines)
  counts <- read.csv(con, header = TRUE, check.names = FALSE)
  close(con)
  names(counts)[1] <- "Gene"
  counts
}

#' Drop genes with zero variance across all samples, and genes with fewer
#' than 2 non-zero samples (uninformative for DE testing).
zero_var_genes <- function(counts) {
  counts %>%
    dplyr::filter(!apply(dplyr::select(., -Gene), 1, function(row) var(row, na.rm = TRUE) == 0)) %>%
    dplyr::filter(rowSums(dplyr::select(., -Gene) != 0) > 1)
}

#' Convert internal genotype levels ("Control", "C4_OE") to display labels
#' ("Control", "C4-OE"). Used everywhere output is written or plotted, so the
#' hyphenated label never has to survive R's model-matrix name mangling.
display_genotype <- function(x) {
  out <- unname(GENOTYPE_DISPLAY[as.character(x)])
  factor(out, levels = unname(GENOTYPE_DISPLAY[GENOTYPE_LEVELS]))
}

#' Convert raw sample IDs (which must match raw_counts.csv column headers)
#' to display labels: WT_10 -> Control_10, mC4_4 -> C4-OE_4.
display_sample <- function(x) {
  out <- sub("^WT_", "Control_", as.character(x))
  sub("^mC4_", "C4-OE_", out)
}

#' Add display-label columns (genotype_label, sample_label) to a metadata
#' table without disturbing the internal columns DESeq2 relies on.
add_display_labels <- function(metadata) {
  metadata$genotype_label <- display_genotype(metadata$genotype)
  metadata$sample_label   <- display_sample(metadata$sample)
  metadata
}

#' Attach genotype + (inferred) sex metadata to a counts matrix's sample columns.
#' `metadata_table` should be SAMPLE_METADATA from 00_config.R (or equivalent).
build_metadata <- function(counts, metadata_table) {
  samples <- colnames(counts)[colnames(counts) != "Gene"]
  meta <- metadata_table[match(samples, metadata_table$sample), ]
  if (anyNA(meta$sample)) {
    stop("Some sample columns in the counts matrix are missing from SAMPLE_METADATA: ",
         paste(samples[is.na(meta$sample)], collapse = ", "))
  }
  # Control is the reference level, so the DESeq2 coefficient is C4-OE vs Control.
  meta$genotype <- factor(meta$genotype, levels = GENOTYPE_LEVELS)
  meta$sex      <- factor(meta$sex, levels = c("F", "M"))
  meta <- add_display_labels(meta)
  rownames(meta) <- meta$sample
  meta
}

get_library_size <- function(count_data) {
  sample_names <- colnames(count_data)[colnames(count_data) != "Gene"]
  reads <- colSums(count_data[, sample_names])
  tibble(sample = sample_names, total_reads = reads)
}

#' Map Ensembl gene IDs (with or without version suffix) to MGI gene symbols
#' using the local org.Mm.eg.db package. No internet connection required, and
#' fully reproducible across machines/dates, so mappings do not depend on
#' Ensembl mirror uptime or drift between Ensembl releases (as live
#' biomaRt::useEnsembl() calls can, silently returning different mappings on
#' different days).
map_ensembl_to_symbol <- function(counts) {
  ens_ids <- gsub("\\.\\d+$", "", counts$Gene)
  symbols <- suppressMessages(
    AnnotationDbi::mapIds(
      org.Mm.eg.db,
      keys      = ens_ids,
      column    = "SYMBOL",
      keytype   = "ENSEMBL",
      multiVals = "first"
    )
  )
  out <- counts
  out$Gene <- unname(symbols)
  out <- out[!is.na(out$Gene), ]
  # Collapse duplicate symbols (a small number of Ensembl IDs map to the same
  # symbol) by summing counts, which is the standard approach for count data.
  sample_cols <- colnames(out)[colnames(out) != "Gene"]
  out <- out %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(sample_cols), sum), .groups = "drop")
  as.data.frame(out)
}

#' Flag (do not silently drop) sex-linked genes. Returns the same data frame
#' with an added `sex_linked` logical column so downstream users can choose
#' to exclude them from headline figures while keeping them auditable in the
#' full supplementary table.
flag_sex_linked_genes <- function(df, gene_col = "gene", sex_genes = SEX_LINKED_GENES) {
  df[["sex_linked"]] <- df[[gene_col]] %in% sex_genes
  df
}

# -----------------------------------------------------------------------------
# 2. Exploratory plots
# -----------------------------------------------------------------------------

lib_size_bar <- function(library_data) {
  library_data$sample <- display_sample(library_data$sample)
  library_data$sample <- factor(library_data$sample, levels = unique(library_data$sample))
  ggplot(library_data, aes(x = sample, y = total_reads, fill = sample)) +
    geom_bar(stat = "identity") +
    labs(title = "Library Size per Sample", x = "Sample", y = "Number of Counts") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_variance_vs_mean <- function(filtered_counts) {
  gene_stats <- filtered_counts %>%
    pivot_longer(cols = -Gene, names_to = "sample", values_to = "value") %>%
    group_by(Gene) %>%
    summarize(mean_count = mean(value), variance = var(value), .groups = "drop") %>%
    arrange(mean_count) %>%
    mutate(mean_rank = rank(mean_count, ties.method = "first"))

  ggplot(gene_stats, aes(x = mean_rank, y = variance)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_smooth() +
    scale_y_continuous(trans = "log10") +
    labs(title = "Variance vs. Mean", x = "Rank(Mean)", y = "Variance") +
    theme_bw()
}

#' Sex-marker diagnostic plot: Xist vs. a Y-linked gene, colored by inferred
#' sex and shaped by genotype. Serves as the visual evidence for the sex
#' calls in SAMPLE_METADATA and belongs in supplementary material.
plot_sex_markers <- function(marker_df) {
  marker_df$genotype <- display_genotype(marker_df$genotype)
  marker_df$sample   <- display_sample(marker_df$sample)
  ggplot(marker_df, aes(x = Xist, y = Ddx3y, color = sex, shape = genotype)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(aes(label = sample), show.legend = FALSE, size = 3) +
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    labs(title = "Sex inference from Xist / Ddx3y raw counts",
         x = "Xist (raw counts)", y = "Ddx3y (raw counts)",
         color = "Inferred sex", shape = "Genotype") +
    theme_bw()
}

# -----------------------------------------------------------------------------
# 3. Differential expression (DESeq2, design = ~ sex + genotype)
# -----------------------------------------------------------------------------

#' Run DESeq2 with sex as an additive covariate.
#'
#' design = ~ sex + genotype  =>  the genotype coefficient (C4-OE vs Control) is
#' estimated after adjusting for the sex effect learned from all 8 samples
#' (using the sex contrast present within Control, since C4-OE is all-male --
#' see the confound note in 00_config.R). Control is the reference level.
run_deseq2 <- function(counts_symbols, metadata) {
  counts_symbols$Gene <- make.names(counts_symbols$Gene, unique = TRUE)
  count_mat <- as.matrix(counts_symbols[, metadata$sample])
  rownames(count_mat) <- counts_symbols$Gene
  storage.mode(count_mat) <- "integer"

  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = metadata,
    design    = ~ sex + genotype
  )
  dds <- DESeq(dds)

  vst_data   <- vst(dds, blind = FALSE)
  vst_matrix <- assay(vst_data)

  res <- results(dds, name = DESEQ_RESULT_NAME)

  results_table <- as_tibble(as.data.frame(res)) %>%
    mutate(gene = rownames(count_mat)) %>%
    relocate(gene) %>%
    dplyr::rename(logFC = log2FoldChange, PValue = pvalue) %>%
    dplyr::select(gene, baseMean, logFC, lfcSE, stat, PValue, padj)

  norm_counts <- as_tibble(as.data.frame(counts(dds, normalized = TRUE))) %>%
    mutate(gene = rownames(count_mat)) %>%
    relocate(gene)

  list(
    dds             = dds,
    results_table   = results_table,
    norm_counts     = norm_counts,
    vst_matrix      = vst_matrix,
    dispersion_plot = function() plotDispEsts(dds)
  )
}

#' Secondary model used ONLY to flag which genotype-DE genes are also
#' strongly sex-differential within Control (n=2 vs n=2, so treat as a screening
#' flag, not a confirmatory test). This directly documents the confound
#' rather than leaving it implicit.
flag_sex_confounded_genes <- function(counts_symbols, metadata, deseq_results,
                                       padj_thresh = PADJ_THRESHOLD) {
  wt_meta <- metadata[metadata$genotype == "Control", ]
  if (length(unique(wt_meta$sex)) < 2 || min(table(wt_meta$sex)) < 2) {
    warning("Not enough Control samples per sex to screen for sex-confounded genes; skipping.")
    deseq_results$sex_confound_flag <- NA
    return(deseq_results)
  }
  counts_symbols$Gene <- make.names(counts_symbols$Gene, unique = TRUE)
  count_mat <- as.matrix(counts_symbols[, wt_meta$sample])
  rownames(count_mat) <- counts_symbols$Gene
  storage.mode(count_mat) <- "integer"

  dds_sex <- DESeqDataSetFromMatrix(count_mat, colData = wt_meta, design = ~ sex)
  dds_sex <- DESeq(dds_sex)
  res_sex <- as.data.frame(results(dds_sex, name = "sex_M_vs_F"))
  res_sex$gene <- rownames(count_mat)

  sex_hits <- res_sex$gene[!is.na(res_sex$padj) & res_sex$padj < padj_thresh]
  deseq_results$sex_confound_flag <- deseq_results$gene %in% sex_hits
  deseq_results
}

# -----------------------------------------------------------------------------
# 4. Post-DEG visualization (consistent color mapping and labeled-gene
#    convention across the differential-expression figures)
# -----------------------------------------------------------------------------

plot_logfc_hist <- function(res_table, padj_threshold, title) {
  filtered <- res_table %>% dplyr::filter(padj < padj_threshold)
  ggplot(filtered, aes(x = logFC)) +
    geom_histogram(binwidth = 0.2, fill = "royalblue1", color = "black") +
    labs(title = title, x = "Log2FoldChange", y = "Count") +
    theme_bw()
}

plot_pca <- function(vst_matrix, metadata, color_by = "genotype_label", shape_by = "sex",
                      title = "PCA") {
  pca_results <- prcomp(t(vst_matrix))
  variance <- pca_results$sdev^2
  var_explained <- variance / sum(variance) * 100
  pca_df <- as.data.frame(pca_results$x)
  pca_df <- cbind(pca_df, metadata[rownames(pca_df), , drop = FALSE])

  ggplot(pca_df, aes(x = PC1, y = PC2,
                      color = .data[[color_by]], shape = .data[[shape_by]])) +
    geom_point(size = 3) +
    stat_ellipse(aes(group = .data[[color_by]]), linetype = "dashed") +
    labs(title = title, color = "Genotype", shape = "Sex",
         x = paste0("PC1: ", round(var_explained[1], 0), "% Variance"),
         y = paste0("PC2: ", round(var_explained[2], 0), "% Variance")) +
    theme_bw()
}

plot_dispersion_alt <- function(vst_matrix, gene_col_name = "gene") {
  # Simple mean-variance sanity plot on the VST matrix (complements
  # DESeq2::plotDispEsts, which needs the dds object directly).
  m <- rowMeans(vst_matrix)
  v <- apply(vst_matrix, 1, var)
  ggplot(data.frame(mean = m, variance = v), aes(mean, variance)) +
    geom_point(alpha = 0.2, size = 0.5) +
    geom_smooth() +
    labs(title = "VST mean-variance relationship", x = "Mean (VST)", y = "Variance (VST)") +
    theme_bw()
}

plot_volcano <- function(res_table, padj_threshold, logfc_threshold, title,
                          labeled_genes = c("Spink10", "Prss56", "C4b", "Tmem41b.ps",
                                             "Gdpd3", "Plin4"),
                          highlight_sex_confound = TRUE) {
  # Color reflects significance (padj) and direction only -- a gene with a
  # small fold change but a significant padj is still colored UP/DOWN, not
  # grey. logfc_threshold is kept only as a visual reference (see vline
  # below if added by the caller); it no longer gates which points get color.
  df <- res_table %>%
    mutate(status = case_when(
      is.na(padj) | padj >= padj_threshold ~ "NS",
      logFC > 0                            ~ "UP",
      logFC < 0                            ~ "DOWN",
      TRUE                                 ~ "NS"
    ))

  p <- ggplot(df, aes(x = logFC, y = -log10(padj), color = status)) +
    geom_point(size = 1, alpha = 0.7) +
    scale_color_manual(values = c(UP = "#EB3F20", DOWN = "#6C7DFF", NS = "grey80"),
                        name = "Status") +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold),
               linetype = "dashed", color = "grey40") +
    theme_bw() +
    theme(legend.position = "bottom") +
    ggtitle(title)

  if (highlight_sex_confound && "sex_confound_flag" %in% names(df)) {
    p <- p + geom_point(
      data = dplyr::filter(df, sex_confound_flag %in% TRUE, status != "NS"),
      shape = 1, size = 2.4, color = "black", stroke = 0.8
    )
  }

  labeled_data <- df[df$gene %in% labeled_genes, ]
  p + ggrepel::geom_text_repel(data = labeled_data, aes(label = gene), color = "black", size = 3)
}

plot_enhanced_volcano <- function(df, select_lab = NULL, title = NULL,
                                   fc_cutoff = LOGFC_THRESHOLD, p_cutoff = PADJ_THRESHOLD,
                                   n_labels = 25) {
  # Default: label only the top N most significant genes clearing both
  # cutoffs, not every gene in the dataset. EnhancedVolcano's `lab` argument
  # is the pool of POSSIBLE labels; `selectLab` is what actually gets drawn.
  # Leaving selectLab = every gene would produce an unreadable wall of text
  # on a genome-wide volcano plot.
  if (is.null(select_lab)) {
    sig <- df %>%
      dplyr::filter(!is.na(padj), padj < p_cutoff, abs(logFC) > fc_cutoff) %>%
      dplyr::arrange(padj)
    select_lab <- head(sig$gene, n_labels)
    if (length(select_lab) == 0) {
      warning("No genes pass the padj/logFC cutoffs for labeling; plot will have no gene labels.")
    }
  }
  EnhancedVolcano(df,
    lab = df$gene, x = "logFC", y = "padj",
    pCutoff = p_cutoff, FCcutoff = fc_cutoff,
    selectLab = select_lab, labSize = 3, labCol = "black", labFace = "bold",
    boxedLabels = TRUE, colAlpha = 0.8, legendPosition = "right",
    drawConnectors = TRUE, widthConnectors = 1.0, colConnectors = "black",
    maxoverlapsConnectors = Inf, title = title, subtitle = NULL,
    caption = paste0(nrow(df), " genes tested; top ", length(select_lab), " labeled"),
    gridlines.major = FALSE, gridlines.minor = FALSE
  )
}

#' Classify genes into curated biological categories using propagated GO
#' annotations from org.Mm.eg.db (GOALL: a gene annotated to a specific
#' child term, e.g. "cholesterol biosynthetic process", is also captured
#' under its broader parent, e.g. "lipid metabolic process"). This is
#' reproducible and auditable (anyone can re-derive the same gene lists from
#' the GO IDs below) rather than a hand-curated, undocumented gene list.
#'
#'   lipid_metabolism : GO:0006629 "lipid metabolic process" (BP)
#'   synapse          : GO:0045202 "synapse" (CC)
#'   immune           : GO:0002376 "immune system process" (BP)
annotate_gene_categories <- function(genes) {
  go_ids <- c(
    lipid_metabolism = "GO:0006629",
    synapse           = "GO:0045202",
    immune            = "GO:0002376"
  )
  ann <- suppressMessages(
    AnnotationDbi::select(
      org.Mm.eg.db, keys = unique(genes), keytype = "SYMBOL", columns = "GOALL"
    )
  )
  out <- data.frame(gene = unique(genes), stringsAsFactors = FALSE)
  for (cat_name in names(go_ids)) {
    hit_genes <- unique(ann$SYMBOL[ann$GOALL == go_ids[[cat_name]]])
    out[[cat_name]] <- out$gene %in% hit_genes
  }
  out
}

#' Volcano plot colored by significance/direction (UP = red, DOWN = blue, NS
#' = grey), where the
#' text labels are still chosen from the curated GO categories (lipid
#' metabolism / synapse / immune), capped at `max_labels_per_category` genes
#' per category (by padj), plus any genes explicitly requested via
#' `force_label_genes` regardless of category or significance. Category
#' membership drives WHICH points get labeled; it no longer drives point
#' color -- that stays a straightforward significance/direction plot.
plot_category_volcano <- function(res_table, padj_threshold = PADJ_THRESHOLD,
                                   logfc_threshold = LOGFC_THRESHOLD,
                                   title = "DEGs by biological category",
                                   max_labels_per_category = 8,
                                   force_label_genes = c("Hmgcr", "Plin4","C4b")) {
  cats <- annotate_gene_categories(res_table$gene)
  df <- res_table %>% dplyr::left_join(cats, by = "gene")
  cat_cols <- c("lipid_metabolism", "synapse", "immune")
  df[cat_cols] <- lapply(df[cat_cols], function(x) ifelse(is.na(x), FALSE, x))

  df <- df %>%
    dplyr::mutate(
      # Significance is padj-only; logfc_threshold is retained purely as a
      # visual reference line, not a color/label cutoff -- a small-fold-change
      # but padj-significant gene is still colored and eligible for labeling.
      significant = !is.na(padj) & padj < padj_threshold,
      status = dplyr::case_when(
        !significant ~ "NS",
        logFC > 0     ~ "UP",
        logFC < 0     ~ "DOWN",
        TRUE          ~ "NS"
      ),
      n_categories = lipid_metabolism + synapse + immune,
      category = dplyr::case_when(
        n_categories >= 2  ~ "Multiple categories",
        lipid_metabolism   ~ "Lipid metabolism",
        synapse            ~ "Synapse",
        immune             ~ "Immune / inflammatory",
        TRUE               ~ NA_character_
      )
    )

  status_colors <- c(UP = "#EB3F20", DOWN = "#6C7DFF", NS = "grey80")

  category_labels <- df %>%
    dplyr::filter(significant, !is.na(category)) %>%
    dplyr::group_by(category) %>%
    dplyr::arrange(padj, .by_group = TRUE) %>%
    dplyr::slice_head(n = max_labels_per_category) %>%
    dplyr::ungroup()

  forced_labels <- df %>% dplyr::filter(gene %in% force_label_genes)
  missing_forced <- setdiff(force_label_genes, df$gene)
  if (length(missing_forced) > 0) {
    warning("force_label_genes not found in results table (skipped): ",
            paste(missing_forced, collapse = ", "))
  }

  label_df <- dplyr::bind_rows(category_labels, forced_labels) %>%
    dplyr::distinct(gene, .keep_all = TRUE)

  message(sprintf(
    "plot_category_volcano: labeling %d category genes (up to %d/category) + %d forced gene(s) = %d labels total.",
    nrow(category_labels), max_labels_per_category, nrow(forced_labels), nrow(label_df)
  ))

  ggplot(df, aes(x = logFC, y = -log10(padj), color = status)) +
    geom_point(size = 1, alpha = 0.7) +
    scale_color_manual(values = status_colors, name = "Status") +
    ggrepel::geom_text_repel(
      data = label_df, aes(label = gene), size = 3, color = "black",
      max.overlaps = Inf, min.segment.length = 0, segment.size = 0.3, show.legend = FALSE
    ) +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "grey40") +
    labs(title = title, x = paste0("log2 fold change (", CONTRAST_LABEL, ")"),
         y = expression(-log[10] ~ "(adjusted p-value)")) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
}

plot_heatmap_top_deg <- function(res_table, norm_counts, metadata, top_n = 50,
                                  title = "Top DEGs (VST-normalized)") {
  sig <- res_table %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = top_n)

  # Built manually (rather than via tibble::column_to_rownames()) because
  # norm_counts/metadata may already carry a non-default rownames attribute
  # from upstream matrix conversions, which column_to_rownames() refuses to
  # overwrite. Plain base-R assignment sidesteps that entirely.
  sub_df <- as.data.frame(norm_counts) %>% dplyr::filter(gene %in% sig$gene)
  mat <- as.matrix(sub_df[, metadata$sample, drop = FALSE])
  rownames(mat) <- sub_df$gene
  mat <- log2(mat + 1)
  mat <- t(scale(t(mat)))  # z-score by gene
  # Show display labels (Control_10 / C4-OE_4) on the heatmap columns rather
  # than the raw raw_counts.csv headers (WT_10 / mC4_4).
  colnames(mat) <- display_sample(colnames(mat))

  annotation_col <- data.frame(
    Genotype = display_genotype(metadata$genotype),
    Sex      = metadata$sex
  )
  rownames(annotation_col) <- display_sample(metadata$sample)

  pheatmap(mat,
    annotation_col = annotation_col,
    show_rownames = TRUE, fontsize_row = 6,
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    main = title, silent = TRUE
  )
}

# -----------------------------------------------------------------------------
# 5. GO over-representation analysis (ORA) -- WITH gene universe
# -----------------------------------------------------------------------------

#' clusterProfiler::enrichGO with an explicit universe (the fix for the
#' original bug). `universe` should be every gene that was actually tested
#' for DE (i.e. res_table$gene after filtering), not left as the
#' clusterProfiler default (every annotated gene in org.Mm.eg.db).
run_go_ora <- function(sig_genes, universe, ont = "BP") {
  clusterProfiler::enrichGO(
    gene          = sig_genes,
    universe      = universe,
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )
}

#' Barplot of enriched GO terms.
#'
#' Terms are FILTERED to p.adjust < padj_threshold before plotting, then the
#' top `num_paths` are taken by gene Count. Because run_go_ora() deliberately
#' sets enrichGO's pvalueCutoff/qvalueCutoff to 1 so the exported CSVs contain
#' every tested term, the significance filtering has to happen here instead --
#' without it the barplots show non-significant terms, which a reader would
#' reasonably assume were enriched.
#'
#' Returns NULL (with a warning) when no term passes the threshold, so callers
#' can skip writing an empty figure.
make_go_barplot <- function(go_res, num_paths = 20, title = "GO Enrichment",
                             order_by = c("Count", "p.adjust"),
                             padj_threshold = PADJ_THRESHOLD) {
  order_by <- match.arg(order_by)
  df <- as.data.frame(go_res)

  n_total <- nrow(df)
  df <- df %>% dplyr::filter(!is.na(p.adjust), p.adjust < padj_threshold)
  if (nrow(df) == 0) {
    warning(sprintf("make_go_barplot: no terms pass p.adjust < %.2f (of %d tested) for '%s'; no plot produced.",
                     padj_threshold, n_total, title))
    return(NULL)
  }

  df <- if (order_by == "Count") {
    df %>% dplyr::arrange(desc(Count))
  } else {
    df %>% dplyr::arrange(p.adjust)
  }
  n_sig <- nrow(df)
  df <- head(df, num_paths)
  message(sprintf("make_go_barplot: %d/%d terms significant at p.adjust < %.2f; plotting top %d by %s.",
                   n_sig, n_total, padj_threshold, nrow(df), order_by))

  ggplot(df, aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "blue", high = "red", name = "p.adjust") +
    coord_flip() +
    labs(title = title, x = NULL, y = "Gene count") +
    theme_minimal()
}

#' Resolve the full set of synaptic GO BP term IDs: the roots defined in
#' SYNAPTIC_GO_BP_ROOTS (00_config.R) plus every GO descendant of those roots.
#' Using the GO hierarchy rather than keyword-matching term names means child
#' terms like "neurotransmitter secretion" or "long-term synaptic
#' potentiation" are included even though their names lack the word "synapse".
get_synaptic_go_ids <- function(roots = SYNAPTIC_GO_BP_ROOTS) {
  if (!requireNamespace("GO.db", quietly = TRUE)) {
    stop("GO.db is required to resolve synaptic GO terms. Install with: ",
         "BiocManager::install('GO.db')")
  }
  offspring <- AnnotationDbi::as.list(GO.db::GOBPOFFSPRING[roots])
  ids <- unique(c(roots, unlist(offspring, use.names = FALSE)))
  ids[!is.na(ids)]
}

#' Barplot of the top N enriched GO terms restricted to a given set of GO IDs
#' (e.g. the synaptic set from get_synaptic_go_ids()). Ordered by adjusted
#' p-value so the most significant terms are shown, and filtered to terms
#' passing `padj_threshold`.
plot_go_subset_barplot <- function(go_res, keep_ids, num_paths = 10,
                                    padj_threshold = PADJ_THRESHOLD,
                                    title = "GO enrichment (subset)") {
  df <- as.data.frame(go_res) %>%
    dplyr::filter(ID %in% keep_ids, p.adjust < padj_threshold) %>%
    dplyr::arrange(p.adjust)

  if (nrow(df) == 0) {
    warning("No terms in the requested GO subset pass padj < ", padj_threshold,
            "; no plot produced.")
    return(NULL)
  }
  df <- head(df, num_paths)
  message(sprintf("plot_go_subset_barplot: %d matching term(s) significant; plotting top %d.",
                   nrow(as.data.frame(go_res) %>%
                          dplyr::filter(ID %in% keep_ids, p.adjust < padj_threshold)),
                   nrow(df)))

  ggplot(df, aes(x = reorder(Description, -p.adjust), y = Count, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "blue", high = "red", name = "p.adjust") +
    coord_flip() +
    labs(title = title, x = NULL, y = "Gene count") +
    theme_minimal()
}

# -----------------------------------------------------------------------------
# 6. GSEA (fgsea against MSigDB collections)
# -----------------------------------------------------------------------------

#' Build the ranked gene list that fgsea consumes.
#'
#' IMPORTANT -- ranking metric. The default is the DESeq2 Wald statistic
#' (`stat` = logFC / lfcSE), NOT raw logFC. This matters enormously here:
#'
#' GSEA is driven entirely by which genes sit at the two extremes of the
#' ranked list. Ranking by raw logFC puts low-count, high-variance genes
#' there -- in this dataset, of the 400 most extreme genes by logFC, only 23
#' had a testable adjusted p-value and their median baseMean was 1.7, i.e.
#' the poles of the list were essentially noise. Meanwhile the genuinely
#' coordinated biology (cholesterol biosynthesis: Hmgcs1, Hmgcr, Sqle,
#' Msmo1, Cyp51, Dhcr7 ...) has modest fold changes (~0.3-0.66) despite very
#' small p-values, so logFC-ranking buried those genes around rank 3000-5000
#' of 22,761 -- exactly where GSEA cannot see them.
#'
#' Ranking by `stat` divides by the standard error, so uncertain low-count
#' genes collapse toward zero and reliable coordinated changes rise to the
#' top (the same cholesterol genes move to ranks 5, 17, 37, 67, 94, 102).
#' This is the standard, recommended DESeq2 -> fgsea handoff.
#'
#' `metric = "logFC"` ranks by log2 fold change instead of the Wald statistic.
#'
#' `drop_untested` removes genes with any NA field -- in practice the ~7,065
#' low-count genes DESeq2's independent filtering left with NA padj. This
#' filter is essential: without it, those noise genes occupy the extremes of
#' the ranked list (median
#' baseMean 1.7 for the 400 most extreme genes, vs 72.0 once filtered) and
#' GSEA returns essentially nothing.
make_ranked_list <- function(res_table, metric = c("stat", "logFC"),
                              drop_untested = TRUE) {
  metric <- match.arg(metric)
  if (!metric %in% names(res_table)) {
    stop("Requested ranking metric '", metric, "' is not a column of the results table.")
  }
  ranked <- res_table
  if (drop_untested) {
    ranked <- ranked %>% dplyr::filter(!is.na(padj))
  }
  ranked <- ranked %>%
    dplyr::filter(!is.na(.data[[metric]]), !is.na(gene)) %>%
    dplyr::distinct(gene, .keep_all = TRUE)
  rnk <- ranked[[metric]]
  names(rnk) <- ranked$gene
  message(sprintf("make_ranked_list: ranking %d genes by '%s'%s.",
                   length(rnk), metric,
                   if (drop_untested) " (untested/NA-padj genes dropped)" else ""))
  sort(rnk, decreasing = TRUE)
}

run_fgsea <- function(gmt_file_path, ranked_list, min_size = GO_MIN_SET_SIZE,
                       max_size = GO_MAX_SET_SIZE) {
  pathways <- fgsea::gmtPathways(gmt_file_path)
  res <- fgsea(pathways, ranked_list, minSize = min_size, maxSize = max_size)
  as_tibble(res) %>% dplyr::arrange(padj)
}

#' Look up a gene set from a .gmt file by matching its name against a regex.
#' Returns list(name = <matched set name>, genes = <character vector>), or
#' NULL if nothing matches. Reading the set from the GMT (rather than
#' hard-coding a gene list) keeps the figure reproducible and tied to the
#' same MSigDB version used everywhere else in the pipeline.
get_gmt_pathway_genes <- function(gmt_file, pattern, ignore_case = TRUE) {
  if (!file.exists(gmt_file)) {
    warning("GMT file not found: ", gmt_file)
    return(NULL)
  }
  pw <- fgsea::gmtPathways(gmt_file)
  hits <- grep(pattern, names(pw), value = TRUE, ignore.case = ignore_case)
  if (length(hits) == 0) {
    warning("No gene set in ", basename(gmt_file), " matched pattern: ", pattern)
    return(NULL)
  }
  if (length(hits) > 1) {
    message(sprintf("get_gmt_pathway_genes: %d sets matched '%s'; using '%s'. Others: %s",
                     length(hits), pattern, hits[1], paste(hits[-1], collapse = ", ")))
  }
  list(name = hits[1], genes = pw[[hits[1]]])
}

#' Lollipop plot of per-gene fold changes for a defined gene set.
#'
#' Each gene gets a horizontal stem from 0 to its log2 fold change, with a
#' point at the end coloured by significance (padj only -- no fold-change
#' cutoff) and sized by the magnitude of the change. Genes are ordered by
#' log2FC so the coordinated direction of a pathway is immediately readable.
#'
#' Note on point size: sizing is by |log2FC|, because ggplot2's size scale
#' cannot represent negative values. The legend is labelled accordingly.
plot_pathway_lollipop <- function(res_table, genes, title = "Pathway genes",
                                   padj_threshold = PADJ_THRESHOLD,
                                   drop_untested = TRUE) {
  df <- res_table %>%
    dplyr::filter(gene %in% genes, !is.na(logFC))
  if (drop_untested) df <- df %>% dplyr::filter(!is.na(padj))

  missing <- setdiff(genes, df$gene)
  if (length(missing) > 0) {
    message(sprintf("plot_pathway_lollipop: %d of %d set genes not plotted (absent or untested): %s",
                     length(missing), length(genes),
                     paste(utils::head(missing, 10), collapse = ", ")))
  }
  if (nrow(df) == 0) {
    warning("No genes from this set are present in the results table; no plot produced.")
    return(NULL)
  }

  df <- df %>%
    dplyr::mutate(
      Significance = ifelse(!is.na(padj) & padj < padj_threshold,
                             "Significant", "Not Significant"),
      Significance = factor(Significance, levels = c("Not Significant", "Significant"))
    ) %>%
    dplyr::arrange(logFC)
  df$gene_f <- factor(df$gene, levels = df$gene)  # ascending -> largest at top

  message(sprintf("plot_pathway_lollipop: %d genes plotted, %d significant at padj < %.2f.",
                   nrow(df), sum(df$Significance == "Significant"), padj_threshold))

  ggplot(df, aes(x = logFC, y = gene_f)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = logFC, y = gene_f, yend = gene_f,
                      color = Significance), linewidth = 0.5) +
    geom_point(aes(color = Significance, size = abs(logFC))) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "blue"),
                        name = "Significance") +
    scale_size_continuous(name = "|LogFC|", range = c(1.5, 5)) +
    labs(title = title, x = expression(Log[2] ~ "Fold Change (LogFC)"), y = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 9))
}

#' Extract GSEA results for pathways matching a keyword pattern, across a
#' named list of fgsea result tables (one per MSigDB collection). Returns a
#' single tidy table with a `collection` column, sorted by padj. Useful for
#' interrogating a specific biological theme (e.g. cholesterol/sterol) without
#' hunting through six separate CSVs.
collect_gsea_by_keyword <- function(gsea_results, pattern, ignore_case = TRUE) {
  purrr::imap_dfr(gsea_results, function(res, collection) {
    if (is.null(res) || nrow(res) == 0) return(NULL)
    hits <- res %>% dplyr::filter(grepl(pattern, pathway, ignore.case = ignore_case))
    if (nrow(hits) == 0) return(NULL)
    hits$collection <- collection
    hits
  }) %>%
    dplyr::arrange(padj) %>%
    dplyr::relocate(collection)
}

plot_gsea_top <- function(fgsea_res, num_paths = 15, padj_threshold = PADJ_THRESHOLD,
                           title = "GSEA results") {
  pos <- fgsea_res %>% dplyr::filter(NES > 0, padj <= padj_threshold) %>%
    dplyr::arrange(desc(NES)) %>% head(num_paths)
  neg <- fgsea_res %>% dplyr::filter(NES < 0, padj <= padj_threshold) %>%
    dplyr::arrange(NES) %>% head(num_paths)
  top <- dplyr::bind_rows(pos, neg)

  if (nrow(top) == 0) {
    warning("No gene sets pass padj <= ", padj_threshold, "; no plot produced.")
    return(NULL)
  }
  if (nrow(top) < 3) {
    message(sprintf(
      "plot_gsea_top: only %d gene set(s) significant -- a barplot with so few bars is usually not worth showing as a figure.",
      nrow(top)))
  }

  ggplot(top, aes(x = fct_reorder(pathway, NES), y = NES, fill = padj)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_gradient(low = "orangered", high = "pink", name = "p.adjust") +
    coord_flip() +
    labs(title = title, x = NULL, y = "Normalized Enrichment Score (NES)") +
    theme_bw(base_size = 9)
}

# -----------------------------------------------------------------------------
# 7. WGCNA co-expression network
# -----------------------------------------------------------------------------

prepare_wgcna_matrix <- function(vst_matrix) {
  t(vst_matrix)  # WGCNA expects samples-in-rows, genes-in-columns
}

run_soft_threshold_scan <- function(wgcna_input, powers = 1:20) {
  WGCNA::pickSoftThreshold(data = wgcna_input, powerVector = powers, verbose = 0)
}

plot_soft_threshold <- function(sft, chosen_power = NULL,
                                 rsq_cutoff = if (exists("WGCNA_RSQ_CUTOFF")) WGCNA_RSQ_CUTOFF else 0.9) {
  fit <- sft$fitIndices
  p1 <- ggplot(fit, aes(Power, -sign(slope) * SFT.R.sq)) +
    geom_point() + geom_text(aes(label = Power), vjust = -0.6, size = 3) +
    geom_hline(yintercept = rsq_cutoff, color = "red", linetype = "dashed") +
    labs(title = "Scale independence", y = "Signed R^2") + theme_bw()
  p2 <- ggplot(fit, aes(Power, mean.k.)) +
    geom_point() + geom_text(aes(label = Power), vjust = -0.6, size = 3) +
    labs(title = "Mean connectivity", y = "Mean connectivity") + theme_bw()
  # Mark the power actually used, so the figure documents the choice.
  if (!is.null(chosen_power)) {
    p1 <- p1 + geom_vline(xintercept = chosen_power, color = "blue", linetype = "dotted")
    p2 <- p2 + geom_vline(xintercept = chosen_power, color = "blue", linetype = "dotted")
  }
  list(scale_independence = p1, mean_connectivity = p2)
}

run_wgcna_modules <- function(wgcna_input, power, min_module_size = 30,
                               merge_cut_height = 0.25) {
  WGCNA::blockwiseModules(
    datExpr = wgcna_input, power = power, networkType = "signed",
    deepSplit = 2, pamRespectsDendro = FALSE, minModuleSize = min_module_size,
    maxBlockSize = ncol(wgcna_input) + 1, reassignThreshold = 0,
    mergeCutHeight = merge_cut_height, numericLabels = TRUE, verbose = 0
  )
}
