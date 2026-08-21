# =============================================================================
# 04b_wgcna_power_sensitivity.R
# -----------------------------------------------------------------------------
# How stable are the WGCNA modules across the soft-thresholding power?
#
# The soft-threshold scan on this dataset never reaches R^2 = 0.9, so several
# powers are defensible. This script rebuilds the network at each power in
# WGCNA_SENSITIVITY_POWERS, holding everything else fixed (same genes, same
# minModuleSize, same mergeCutHeight), and measures how much the resulting
# partitions agree.
#
# Produces (under Outputs/.../04_wgcna/):
#   tables/wgcna_power_sensitivity_summary.csv      per power: R^2, modules, grey
#   tables/wgcna_power_sensitivity_ari.csv          pairwise adjusted Rand index
#   tables/wgcna_power_sensitivity_jaccard.csv      pairwise mean best-match Jaccard
#   tables/wgcna_power_sensitivity_modules.csv      per-module stability at WGCNA_POWER
#   tables/wgcna_power_sensitivity_assignments.csv  gene x power module assignment
#   figures/wgcna_power_sensitivity.pdf             agreement matrix + module stability
#
# Reads Outputs/rdata/02_deseq2_results.rds; does not need steps 01-04 re-run.
# Runtime is dominated by one TOM per power (a few minutes each).
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

suppressPackageStartupMessages(library(WGCNA))
WGCNA::allowWGCNAThreads(4)

DIR_TABLES  <- step_dir("04", "tables")
DIR_FIGURES <- step_dir("04", "figures")

powers <- WGCNA_SENSITIVITY_POWERS
ref    <- WGCNA_POWER
if (!ref %in% powers) {
  stop("WGCNA_POWER (", ref, ") must be one of WGCNA_SENSITIVITY_POWERS.", call. = FALSE)
}

# ---- 1. Same input matrix as 04 ---------------------------------------------
deseq_cached <- readRDS(file.path(DIR_RDATA, "02_deseq2_results.rds"))
nc  <- deseq_cached$norm_counts
cnt <- as.matrix(nc[, setdiff(names(nc), "gene")]); rownames(cnt) <- nc$gene
expressed <- names(which(rowMeans(cnt) >= WGCNA_MIN_MEAN_COUNT))
vst_filtered <- deseq_cached$vst_matrix[rownames(deseq_cached$vst_matrix) %in% expressed, , drop = FALSE]
datExpr <- prepare_wgcna_matrix(vst_filtered)
genes   <- colnames(datExpr)
message(sprintf("Sensitivity input: %d genes x %d samples; powers %s (reference %d).",
                 length(genes), nrow(datExpr), paste(powers, collapse = ", "), ref))

# ---- 2. One partition per power ---------------------------------------------
# Each partition is built with run_wgcna_modules() -- the SAME function step 04
# uses -- so module names, sizes and the grey bin are directly comparable with
# the reported results, and the only thing varying between runs is the power.
# This recomputes the TOM inside every call, which is the cost of that
# comparability.
partitions <- list(); summary_rows <- list()
for (p in powers) {
  message(sprintf("  power %d: detecting modules ...", p))
  net  <- run_wgcna_modules(datExpr, power = p,
                            min_module_size  = WGCNA_MIN_MODULE_SIZE,
                            merge_cut_height = WGCNA_MERGE_CUT_HEIGHT)
  cols <- labels2colors(net$colors)
  names(cols) <- genes
  partitions[[as.character(p)]] <- cols
  summary_rows[[as.character(p)]] <- data.frame(
    power = p, n_modules = length(setdiff(unique(cols), "grey")),
    grey_genes = sum(cols == "grey"),
    largest_module = max(table(cols[cols != "grey"])),
    median_module  = as.integer(median(table(cols[cols != "grey"]))))
}
sft <- read.csv(file.path(DIR_TABLES, "wgcna_soft_threshold_scan.csv"))
summ <- do.call(rbind, summary_rows)
summ$signed_rsq <- round(sft$signed_rsq[match(summ$power, sft$Power)], 3)
summ$mean_k     <- round(sft$mean.k.[match(summ$power, sft$Power)], 1)
summ <- summ[, c("power","signed_rsq","mean_k","n_modules","grey_genes","largest_module","median_module")]
write.csv(summ, file.path(DIR_TABLES, "wgcna_power_sensitivity_summary.csv"), row.names = FALSE)

# ---- 3. Pairwise agreement ---------------------------------------------------
# Adjusted Rand index: agreement between two partitions of the same genes,
# corrected for the agreement expected by chance. 1 = identical, 0 = chance.
adj_rand <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  ch2 <- function(x) sum(x * (x - 1) / 2)
  idx <- ch2(as.vector(tab)); ea <- ch2(rowSums(tab)); eb <- ch2(colSums(tab))
  exp_idx <- ea * eb / ch2(n); max_idx <- (ea + eb) / 2
  (idx - exp_idx) / (max_idx - exp_idx)
}
# Mean best-match Jaccard: for each module of `a`, the best Jaccard overlap
# with any module of `b`, averaged over modules (size-weighted).
best_jaccard <- function(a, b) {
  ma <- setdiff(unique(a), "grey"); mb <- setdiff(unique(b), "grey")
  sets_b <- lapply(mb, function(k) names(b)[b == k])
  j <- vapply(ma, function(m) {
    g <- names(a)[a == m]
    max(vapply(sets_b, function(h) length(intersect(g, h)) / length(union(g, h)), 0))
  }, 0)
  sizes <- vapply(ma, function(m) sum(a == m), 0)
  sum(j * sizes) / sum(sizes)
}
np <- length(powers)
ari <- jac <- matrix(NA_real_, np, np, dimnames = list(powers, powers))
for (i in seq_len(np)) for (k in seq_len(np)) {
  a <- partitions[[i]]; b <- partitions[[k]]
  ari[i, k] <- if (i == k) 1 else adj_rand(a, b)
  jac[i, k] <- if (i == k) 1 else best_jaccard(a, b)
}
write.csv(round(ari, 3), file.path(DIR_TABLES, "wgcna_power_sensitivity_ari.csv"))
write.csv(round(jac, 3), file.path(DIR_TABLES, "wgcna_power_sensitivity_jaccard.csv"))

# ---- 4. Per-module stability at the reference power --------------------------
# For every gene, compare the set of genes sharing its module at the reference
# power with the set sharing its module at each other power (Jaccard). A gene
# scoring 1 keeps exactly the same neighbours; a gene scoring 0 keeps none.
# Averaging within a module gives that module's stability.
ref_part <- partitions[[as.character(ref)]]
others   <- setdiff(as.character(powers), as.character(ref))
gene_stab <- vapply(genes, function(g) {
  set_ref <- names(ref_part)[ref_part == ref_part[[g]]]
  mean(vapply(others, function(p) {
    q <- partitions[[p]]; set_p <- names(q)[q == q[[g]]]
    length(intersect(set_ref, set_p)) / length(union(set_ref, set_p))
  }, 0))
}, 0)
mod_stab <- data.frame(
  module    = names(tapply(gene_stab, ref_part, mean)),
  size      = as.integer(table(ref_part)[names(tapply(gene_stab, ref_part, mean))]),
  stability = round(as.numeric(tapply(gene_stab, ref_part, mean)), 3))
mod_stab <- mod_stab[order(-mod_stab$stability), ]
write.csv(mod_stab, file.path(DIR_TABLES, "wgcna_power_sensitivity_modules.csv"), row.names = FALSE)

assign_tbl <- data.frame(gene = genes, stability = round(gene_stab, 3),
                         do.call(cbind, lapply(partitions, function(x) unname(x[genes]))),
                         check.names = FALSE)
names(assign_tbl)[-(1:2)] <- paste0("power_", powers)
write.csv(assign_tbl, file.path(DIR_TABLES, "wgcna_power_sensitivity_assignments.csv"), row.names = FALSE)

# ---- 5. Figure ---------------------------------------------------------------
pdf(file.path(DIR_FIGURES, "wgcna_power_sensitivity.pdf"), width = 11, height = 5)
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 2))
image(seq_len(np), seq_len(np), ari, axes = FALSE, xlab = "power", ylab = "power",
      main = "Partition agreement\n(adjusted Rand index)", zlim = c(0, 1),
      col = colorRampPalette(c("white", "steelblue4"))(64))
axis(1, seq_len(np), powers); axis(2, seq_len(np), powers, las = 1)
for (i in seq_len(np)) for (k in seq_len(np))
  text(i, k, sprintf("%.2f", ari[i, k]), cex = 0.8,
       col = ifelse(ari[i, k] > 0.6, "white", "black"))
bp <- barplot(rev(mod_stab$stability), horiz = TRUE, xlim = c(0, 1), las = 1,
              names.arg = rev(mod_stab$module), cex.names = 0.6,
              col = rev(mod_stab$module), border = "grey30",
              xlab = "mean co-membership Jaccard across powers",
              main = sprintf("Module stability at power %d", ref))
abline(v = mean(gene_stab), lty = 2, col = "red")
par(op); dev.off()

# Fold the stability score into the module summary written by 04, so a single
# supplementary table carries size, function, genotype correlation and stability.
summary_file <- file.path(DIR_TABLES, "wgcna_module_summary.csv")
if (file.exists(summary_file)) {
  ms <- read.csv(summary_file, stringsAsFactors = FALSE)
  ms$stability <- mod_stab$stability[match(ms$module, mod_stab$module)]
  ms$stability_ref_power <- ref
  write.csv(ms, summary_file, row.names = FALSE)
  message("  -> stability column added to wgcna_module_summary.csv")
} else {
  message("  -> wgcna_module_summary.csv not found; run 04 first to get the merged table")
}

message(sprintf("Median pairwise ARI: %.3f | mean gene-level stability: %.3f",
                 median(ari[upper.tri(ari)]), mean(gene_stab)))
message("04b_wgcna_power_sensitivity.R complete.")
