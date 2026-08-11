# =============================================================================
# Combined Deconvolution / DEG Projection Script
# Reference Prep + Sham/Control Projection + Treatment/C4 Projection
# No RDS saving
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Libraries
# -----------------------------------------------------------------------------
library(data.table)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(gplots)
library(cowplot)
library(patchwork)
library(pheatmap)
library(dynamicTreeCut)
library(dplyr)
# -----------------------------------------------------------------------------
# 0.5 Set seed for reproducibility
# -----------------------------------------------------------------------------
set.seed(1234)


# -----------------------------------------------------------------------------
# 0.6 Resolve inputs and outputs (portable -- no absolute machine paths)
# -----------------------------------------------------------------------------
# When run via run_full_analysis.R the runner defines MASTER_INPUT / MASTER_OUTPUT.
# When run on its own, locate this script on disk and derive the repository root.
# The repository path contains spaces, which R's front end encodes as "~+~" in
# the --file= argument, so decode that (see CLAUDE.md).
.get_projection_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    path <- gsub("~+~", " ", sub("^--file=", "", file_arg[[1]]), fixed = TRUE)
    return(dirname(normalizePath(path, mustWork = FALSE)))
  }
  for (fr in rev(sys.frames())) {
    if (!is.null(fr$ofile)) {
      return(dirname(normalizePath(fr$ofile, mustWork = FALSE)))
    }
  }
  getwd()
}

if (exists("MASTER_INPUT") && exists("MASTER_OUTPUT")) {
  input_dir    <- MASTER_INPUT
  outputs_root <- MASTER_OUTPUT
} else {
  repo_root    <- dirname(dirname(.get_projection_dir()))  # scripts/03_.. -> root
  input_dir    <- file.path(repo_root, "Input files")
  outputs_root <- file.path(repo_root, "Outputs")
}

strand_out  <- file.path(outputs_root, "03_scRNAseq_reference_projection")
dir_figures <- file.path(strand_out, "figures")
dir_tables  <- file.path(strand_out, "tables")
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)

sc_reference_path <- file.path(input_dir, "sc_data_C4OE_PCA.rds")
mouse_deg_path <- file.path(
  outputs_root, "01_bulk_rnaseq_DE", "02_differential_expression",
  "tables", "S1_Table_DESeq.csv"
)

if (!file.exists(sc_reference_path)) {
  stop("Single-cell reference not found:\n  ", sc_reference_path,
       "\nPlace sc_data_C4OE_PCA.rds in 'Input files/'.", call. = FALSE)
}
if (!file.exists(mouse_deg_path)) {
  stop("Mouse DEG table not found:\n  ", mouse_deg_path,
       "\nRun strand 01 first so it writes S1_Table_DESeq.csv.", call. = FALSE)
}


# =============================================================================
# PART 1: LOAD AND PREPARE SINGLE-CELL REFERENCE
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Load data directly
# -----------------------------------------------------------------------------
#sc_data <- ss.seurat
sc_data <- readRDS(sc_reference_path)


print(sc_data)


# -----------------------------------------------------------------------------
# 2. Subset to neocortical / isocortical regions
# -----------------------------------------------------------------------------
regions <- c(
  "VISp", "ALM", "SSp", "MOp", "ACA", "PL-ILA",
  "ORB", "SSs", "RSP", "AI", "TEa-PERI-ECT",
  "RSPv", "PTLp", "AUD", "GU", "VIS"
)

sc_data <- subset(
  x = sc_data,
  subset = region_label %in% regions
)

subclasses_to_remove <- c(
  "", "CA1-ProS", "CA2", "CA3", "DG",
  "SUB", "SUB-ProS", "NP SUB", "CLA", "NP"
)

sc_data <- subset(
  x = sc_data,
  subset = !(subclass_label %in% subclasses_to_remove)
)

sc_data$region_label <- droplevels(factor(sc_data$region_label))
sc_data$subclass_label <- droplevels(factor(sc_data$subclass_label))

region_counts <- table(sc_data$region_label)
subclass_counts <- table(sc_data$subclass_label)

print(region_counts)
print(subclass_counts)



# -----------------------------------------------------------------------------
# 3. Rename subclass labels into grouped cell types
# -----------------------------------------------------------------------------
Idents(sc_data) <- "subclass_label"

sc_data <- RenameIdents(
  sc_data,
  
  # Non-neuronal
  "Astro"     = "Astrocytes",
  "Endo"      = "Endothelial",
  "Micro-PVM" = "Microglia",
  "Oligo"     = "Oligodendrocytes",
  "VLMC"      = "Mural",
  "SMC-Peri"  = "Mural",
  
  # Inhibitory interneurons
  "Lamp5"     = "Lamp5",
  "Pvalb"     = "Pvalb",
  "Sncg"      = "Sncg",
  "Sst"       = "Sst",
  "Vip"       = "Vip",
  "Meis2"     = "Other_GABA",
  "Sst Chodl" = "Sst_Chodl",
  
  # Excitatory / pyramidal subclasses
  "CR"              = "Cajal-Retzius",
  
  # Upper-layer IT
  "L2 IT ENTl"      = "L2/3_IT",
  "L2/3 IT CTX"     = "L2/3_IT",
  "L2/3 IT ENTl"    = "L2/3_IT",
  "L2/3 IT PPP"     = "L2/3_IT",
  "L3 IT ENT"       = "L2/3_IT",
  
  # Middle/deep IT
  "L4 RSP-ACA"      = "L4/5_IT",
  "L4/5 IT CTX"     = "L4/5_IT",
  "L5 IT CTX"       = "L4/5_IT",
  "L5/6 IT TPE-ENT" = "L4/5_IT",
  
  # Pyramidal tract
  "L5 PT CTX"       = "L5_PT",
  
  # Corticothalamic
  "L6 CT CTX"       = "L6_CT",
  "L6b CTX"         = "L6_CT",
  "L6b/CT ENT"      = "L6_CT",
  
  # Layer 6 IT / near-projecting
  "L6 IT CTX"       = "L6_IT",
  "L6 IT ENTl"      = "L6_IT",
  "Car3"            = "L6_IT",
  "L5/6 NP CTX"     = "L5/6_NP"
)

sc_data$celltype_grouped <- as.character(Idents(sc_data))
sc_data$celltype_grouped <- factor(sc_data$celltype_grouped)

Idents(sc_data) <- "celltype_grouped"

sc_data$celltype_grouped <- factor(
  sc_data$celltype_grouped,
  levels = sort(levels(sc_data$celltype_grouped))
)

Idents(sc_data) <- "celltype_grouped"

grouped_counts <- table(sc_data$celltype_grouped)

print(grouped_counts)


# -----------------------------------------------------------------------------
# 4. Normalize, variable features, scale, PCA
# -----------------------------------------------------------------------------
sc_data <- NormalizeData(
  sc_data,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

sc_data <- FindVariableFeatures(
  sc_data,
  selection.method = "vst",
  nfeatures = 2000
)

sc_data <- ScaleData(
  sc_data,
  features = VariableFeatures(sc_data)
)

sc_data <- RunPCA(
  sc_data,
  features = VariableFeatures(sc_data),
  seed.use = 1234
)

elbow_plot <- ElbowPlot(sc_data, ndims = 50)

ggsave(
  filename = file.path(dir_figures, "elbow_plot.png"),
  plot = elbow_plot,
  width = 7,
  height = 5,
  dpi = 300
)


# -----------------------------------------------------------------------------
# 5. Choose PCs and run t-SNE / UMAP
# -----------------------------------------------------------------------------
pcs <- 14

sc_data <- RunTSNE(
  sc_data,
  dims = 1:pcs,
  check_duplicates = FALSE,
  seed.use = 1234
)

sc_data <- RunUMAP(
  sc_data,
  dims = 1:pcs,
  seed.use = 1234
)


# -----------------------------------------------------------------------------
# 6. Cell-type colors
# -----------------------------------------------------------------------------
cols <- c(
  # Non-neuronal
  "Astrocytes"       = "limegreen",
  "Endothelial"      = "steelblue",
  "Microglia"        = "firebrick2",
  "Mural"            = "magenta",
  "Oligodendrocytes" = "gray52",
  
  # Inhibitory
  "Lamp5"            = "plum1",
  "Pvalb"            = "mediumorchid4",
  "Sncg"             = "orchid1",
  "Sst"              = "darkorchid",
  "Vip"              = "purple3",
  "Other_GABA"       = "purple4",
  "Sst_Chodl"        = "mediumpurple3",
  
  # Excitatory
  "Cajal-Retzius"    = "gold",
  "L2/3_IT"          = "tan1",
  "L4/5_IT"          = "tan2",
  "L5_PT"            = "chocolate2",
  "L6_CT"            = "sienna",
  "L6_IT"            = "tan3",
  "L5/6_NP"          = "wheat3"
)

cols_use <- cols[names(cols) %in% levels(sc_data$celltype_grouped)]

missing_cols <- setdiff(levels(sc_data$celltype_grouped), names(cols))

if (length(missing_cols) > 0) {
  warning("Missing colors for: ", paste(missing_cols, collapse = ", "))
}


# -----------------------------------------------------------------------------
# 7. Reference plots
# -----------------------------------------------------------------------------
pca_reference_plot <- DimPlot(
  sc_data,
  reduction = "pca",
  group.by = "celltype_grouped",
  pt.size = 0.1,
  label = TRUE,
  cols = cols_use
) +
  ggtitle("Reference PCA")

tsne_reference_plot <- DimPlot(
  sc_data,
  reduction = "tsne",
  group.by = "celltype_grouped",
  pt.size = 0.1,
  label = TRUE,
  cols = cols_use
) +
  ggtitle("Reference t-SNE")

umap_reference_plot <- DimPlot(
  sc_data,
  reduction = "umap",
  group.by = "celltype_grouped",
  pt.size = 0.1,
  label = TRUE,
  cols = cols_use
) +
  ggtitle("Reference UMAP")

ggsave(
  filename = file.path(dir_figures, "reference_pca.png"),
  plot = pca_reference_plot,
  width = 8,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(dir_figures, "reference_tsne.png"),
  plot = tsne_reference_plot,
  width = 8,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(dir_figures, "reference_umap.png"),
  plot = umap_reference_plot,
  width = 8,
  height = 7,
  dpi = 300
)


# -----------------------------------------------------------------------------
# 8. Marker validation plots
# -----------------------------------------------------------------------------
broad_marker_genes <- c(
  "Slc17a7", "Satb2",
  "Gad1", "Gad2",
  "Gfap", "Aqp4",
  "Mbp", "Plp1",
  "Cx3cr1", "P2ry12",
  "Cldn5", "Pecam1",
  "Acta2", "Pdgfrb"
)

broad_marker_genes_use <- broad_marker_genes[
  broad_marker_genes %in% rownames(sc_data)
]

broad_marker_dotplot <- DotPlot(
  sc_data,
  features = broad_marker_genes_use,
  group.by = "celltype_grouped",
  cols = c("lightgrey", "red")
) +
  RotatedAxis() +
  ggtitle("Broad Cell-Type Marker Validation") +
  xlab("Marker gene") +
  ylab("Grouped cell type")

ggsave(
  filename = file.path(dir_figures, "broad_marker_dotplot.png"),
  plot = broad_marker_dotplot,
  width = 10,
  height = 7,
  dpi = 300
)

layer_marker_genes <- c(
  "Cux2",
  "Rorb",
  "Bcl11b",
  "Fezf2",
  "Foxp2",
  "Syt6"
)

layer_marker_genes_use <- layer_marker_genes[
  layer_marker_genes %in% rownames(sc_data)
]

layer_marker_dotplot <- DotPlot(
  sc_data,
  features = layer_marker_genes_use,
  group.by = "celltype_grouped",
  cols = c("lightgrey", "red")
) +
  RotatedAxis() +
  ggtitle("Layer-Specific Marker Validation") +
  xlab("Marker gene") +
  ylab("Grouped cell type")

ggsave(
  filename = file.path(dir_figures, "layer_marker_dotplot.png"),
  plot = layer_marker_dotplot,
  width = 8,
  height = 7,
  dpi = 300
)


# =============================================================================
# PART 2: LOAD DESEQ RESULTS
# =============================================================================


# -----------------------------------------------------------------------------
# 9. Load DESeq table
# -----------------------------------------------------------------------------
all_deseq <- read.csv(
  file = mouse_deg_path,
  row.names = 1,
  check.names = FALSE
)

if (!"logFC" %in% colnames(all_deseq)) {
  stop("Column 'logFC' was not found in S1_Table_DESeq.csv.")
}

gene_names <- rownames(sc_data)

print(head(all_deseq))
print(colnames(all_deseq))

all_deseq <- all_deseq[all_deseq$padj < 0.05, ]
# =============================================================================
# PART 3: Downregulated PROJECTION
# =============================================================================


# -----------------------------------------------------------------------------
# 9. Select negative DEGs
# -----------------------------------------------------------------------------
Downregulated_deg_table <- all_deseq[all_deseq$logFC < 0, ]

Downregulated_deg_table <- Downregulated_deg_table[
  order(Downregulated_deg_table$logFC),
]

keep <- rownames(Downregulated_deg_table) %in% rownames(sc_data)
sum(keep)
Downregulated_deg_table <- Downregulated_deg_table[keep, , drop = FALSE]

Downregulated_deg_table <- Downregulated_deg_table[1:250, ]

Downregulated_genes <- rownames(Downregulated_deg_table)

write.csv(
  Downregulated_deg_table,
  file.path(dir_tables, "Downregulated_control_DEGs.csv")
)


# -----------------------------------------------------------------------------
# 10. PCA using Downregulated genes
# -----------------------------------------------------------------------------
sc_data_Downregulated <- ScaleData(
  sc_data,
  features = Downregulated_genes
)

sc_data_Downregulated <- RunPCA(
  sc_data_Downregulated,
  features = Downregulated_genes,
  reduction.name = "pca_Downregulated_deg",
  reduction.key = "DownregulatedDEGPC_",
  seed.use = 1234
)


# -----------------------------------------------------------------------------
# 11. Plot Downregulated PCA
# -----------------------------------------------------------------------------
Downregulated_pca_plot <- DimPlot(
  sc_data_Downregulated,
  reduction = "pca_Downregulated_deg",
  group.by = "celltype_grouped",
  pt.size = 0.1,
  cols = cols_use
) +
  ggtitle("Downregulated DEGs")

ggsave(
  file.path(dir_figures, "Downregulated_pca_plot.png"),
  Downregulated_pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)

make_pc_coordinate_bins <- function(
    seurat_obj,
    reduction_name,
    condition_label,
    celltype_col = "major_class",
    pc_axis = 1,
    n_bins = 4
) {
  
  # 1. Extract PCA coordinates
  if (!reduction_name %in% names(seurat_obj@reductions)) {
    stop(paste("ERROR: Reduction", reduction_name, "not found."))
  }
  emb <- Embeddings(seurat_obj, reduction = reduction_name)
  
  # 2. Safely grab the correct PC column
  pc_col <- paste0(reduction_name, "_", pc_axis)
  if (!pc_col %in% colnames(emb)) {
    pc_col <- colnames(emb)[pc_axis]
  }
  
  # 3. Build the base dataframe
  df <- data.frame(
    cell = rownames(emb),
    PC_value = emb[, pc_col],
    cell_type = seurat_obj@meta.data[rownames(emb), celltype_col]
  )
  
  # Remove missing cell types
  df <- df %>% filter(!is.na(cell_type))
  
  # 4. Bin by COORDINATE range (Equal Width)
  # cut() breaks the range from min(PC_value) to max(PC_value) into equal chunks
  df <- df %>%
    mutate(
      PC_bin = cut(PC_value, breaks = n_bins, labels = FALSE),
      
      # Optional: Capture the exact coordinate ranges as text labels for clarity
      PC_range_label = cut(PC_value, breaks = n_bins)
    )
  
  # 5. Assign human-readable names to the coordinate zones
  if (n_bins == 4) {
    df <- df %>%
      mutate(
        PC_region = case_when(
          PC_bin == 1 ~ paste0("Lowest Range PC", pc_axis),
          PC_bin == 2 ~ paste0("Low-mid Range PC", pc_axis),
          PC_bin == 3 ~ paste0("High-mid Range PC", pc_axis),
          PC_bin == 4 ~ paste0("Highest Range PC", pc_axis)
        )
      )
  } else {
    df <- df %>%
      mutate(
        PC_region = paste0("PC", pc_axis, " Coordinate Bin ", PC_bin)
      )
  }
  
  # 6. Calculate cell-type percentages inside each spatial bin
  percent_df <- df %>%
    mutate(condition = condition_label) %>%
    dplyr::count(condition, PC_region, PC_range_label, cell_type) %>%
    group_by(condition, PC_region) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    ungroup() %>%
    mutate(PC_axis = paste0("PC", pc_axis))

  return(percent_df)
}

pc1_percent <- make_pc_coordinate_bins(
  seurat_obj = sc_data_Downregulated,
  reduction_name = "pca_Downregulated_deg",
  condition_label = "Downregulated_DEGs",
  celltype_col = "celltype_grouped",
  pc_axis = 1,
  n_bins = 4
)
pc1_percent$PC_region <- factor(
  pc1_percent$PC_region,
  levels = c(
    "Lowest Range PC1",
    "Low-mid Range PC1",
    "High-mid Range PC1",
    "Highest Range PC1"
  )
)
PC1 <- ggplot(
  pc1_percent,
  aes(x = PC_region, y = percent, fill = cell_type)
) +
  scale_fill_manual(values = cols, drop = FALSE) +
  geom_bar(stat = "identity") +
  facet_wrap(~PC_axis, nrow = 1) +
  theme_bw() +
  labs(
    title = "Cell-type composition across PCA coordinate bins Down reg genes",
    x = "Coordinate bin",
    y = "Percent of cells"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

PC1
ggsave(
  filename = file.path(dir_figures, "Down_PC1_celltype_composition.png"),
  plot = PC1,
  width = 8,
  height = 6,
  dpi = 300
)

pc2_percent <- make_pc_coordinate_bins(
  seurat_obj = sc_data_Downregulated,
  reduction_name = "pca_Downregulated_deg",
  condition_label = "Downregulated_DEGs",
  celltype_col = "celltype_grouped",
  pc_axis = 2,
  n_bins = 4
)
pc2_percent$PC_region <- factor(
  pc2_percent$PC_region,
  levels = c(
    "Lowest Range PC2",
    "Low-mid Range PC2",
    "High-mid Range PC2",
    "Highest Range PC2"
  )
)

PC2 <- ggplot(
  pc2_percent,
  aes(x = PC_region, y = percent, fill = cell_type)
) +
  scale_fill_manual(values = cols, drop = FALSE) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(
    title = "Cell-type composition across PC2 coordinate bins Down reg genes",
    x = "Coordinate bin",
    y = "Percent of cells"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

PC2

ggsave(
  file.path(dir_figures, "Down_PC2_celltype_composition.png"),
  PC2,
  width = 10,
  height = 6,
  dpi = 300
)
# -----------------------------------------------------------------------------
# 12. Downregulated correlation heatmap
# -----------------------------------------------------------------------------
Downregulated_matrix <- as.matrix(
  GetAssayData(sc_data_Downregulated, layer = "data")[Downregulated_genes, ]
)

Downregulated_variance <- apply(Downregulated_matrix, 1, var)

Downregulated_valid_genes <- names(Downregulated_variance[Downregulated_variance > 0])

Downregulated_matrix_valid <- as.matrix(
  GetAssayData(sc_data_Downregulated, layer = "data")[Downregulated_valid_genes, ]
)

Downregulated_correlations <- cor(
  (t(Downregulated_matrix_valid) + 1),
  method = "pearson"
)

write.csv(
  Downregulated_correlations,
  file.path(dir_tables, "Downregulated_correlations.csv")
)

# -----------------------------------------------------------------------------
# Hierarchical clustering
# -----------------------------------------------------------------------------

gene_distance <- as.dist(1 - Downregulated_correlations)

gene_clustering <- hclust(
  gene_distance,
  method = "average"
)

# -----------------------------------------------------------------------------
# Dynamic tree cutting
# -----------------------------------------------------------------------------

dynamic_modules <- cutreeDynamic(
  dendro = gene_clustering,
  distM = as.matrix(gene_distance),
  
  deepSplit = 2,
  minClusterSize = 10
)

# -----------------------------------------------------------------------------
# Inspect module sizes
# -----------------------------------------------------------------------------

table(dynamic_modules)

# -----------------------------------------------------------------------------
# Create gene-module table
# -----------------------------------------------------------------------------

gene_module_table <- data.frame(
  Gene = rownames(Downregulated_correlations),
  Module = dynamic_modules
)

# -----------------------------------------------------------------------------
# Save module assignments
# -----------------------------------------------------------------------------

write.csv(gene_module_table, file.path(dir_tables, "Downregulated_gene_modules.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# Create module annotation for heatmap in the oreder of the modules + order correlation matrix
# -----------------------------------------------------------------------------

annotation_df <- data.frame(
  Module = factor(dynamic_modules)
)

rownames(annotation_df) <- rownames(Downregulated_correlations)


# -----------------------------------------------------------------------------
# Plot heatmap
# -----------------------------------------------------------------------------
# 1. Define the scale limits and center
min_val <- min(
  Downregulated_correlations,
  na.rm = TRUE
)

max_val <- max(
  Downregulated_correlations,
  na.rm = TRUE
)

center_val <- 0

# Inspect values
min_val
max_val

# 2. Create asymmetric breaks
# Ensure the number of colors on each side of the center is proportional to the range
palette_length <- 100
my_breaks <- c(seq(min_val, center_val, length.out = ceiling(palette_length/2) + 1), 
               seq(max_val/palette_length, max_val, length.out = floor(palette_length/2)))

# 3. Define colors
my_colors <- colorRampPalette(c("#0571b0", "#f7f7f7", "#ca0020"))(palette_length)


p_down<-pheatmap(
  Downregulated_correlations,
  
  cluster_rows = gene_clustering,
  cluster_cols = gene_clustering,
  
  annotation_row = annotation_df,
  annotation_col = annotation_df,
  
  clustering_method = "average",
  
  color=my_colors,
  breaks = my_breaks,
  
  border_color = NA,

  fontsize_row = 4,
  fontsize_col = 4,

  # pheatmap writes the file itself (ggsave() does not work on a pheatmap object).
  filename = file.path(dir_figures, "Down_heatmap_gene_gene_correlations.png"),
  width = 20,
  height = 20
)
# =============================================================================
# PART 4: Upregulated PROJECTION
# =============================================================================


# -----------------------------------------------------------------------------
# 13. Select negative DEGs
# -----------------------------------------------------------------------------
Upregulated_deg_table <- all_deseq[all_deseq$logFC > 0, ]

keep <- rownames(Upregulated_deg_table) %in% rownames(sc_data)
Upregulated_deg_table <- Upregulated_deg_table[keep, , drop = FALSE]

Upregulated_deg_table <- Upregulated_deg_table[
  order(Upregulated_deg_table$logFC, decreasing = TRUE),
  ,
  drop = FALSE
]

Upregulated_deg_table <- head(Upregulated_deg_table, 250)

Upregulated_genes <- rownames(Upregulated_deg_table)


write.csv(
  Upregulated_deg_table,
  file.path(dir_tables, "Upregulated_control_DEGs.csv")
)


# -----------------------------------------------------------------------------
# 14. PCA using Upregulated genes
# -----------------------------------------------------------------------------
sc_data_Upregulated <- ScaleData(
  sc_data,
  features = Upregulated_genes
)

sc_data_Upregulated <- RunPCA(
  sc_data_Upregulated,
  features = Upregulated_genes,
  reduction.name = "pca_Upregulated_deg",
  reduction.key = "UpregulatedDEGPC_",
  seed.use = 1234
)


# -----------------------------------------------------------------------------
# 15. Plot Upregulated PCA
# -----------------------------------------------------------------------------
Upregulated_pca_plot <- DimPlot(
  sc_data_Upregulated,
  reduction = "pca_Upregulated_deg",
  group.by = "celltype_grouped",
  pt.size = 0.1,
  cols = cols_use
) +
  ggtitle("Upregulated DEGs")

ggsave(
  file.path(dir_figures, "Upregulated_pca_plot.png"),
  Upregulated_pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)

make_pc_coordinate_bins <- function(
    seurat_obj,
    reduction_name,
    condition_label,
    celltype_col = "major_class",
    pc_axis = 1,
    n_bins = 4
) {
  
  # 1. Extract PCA coordinates
  if (!reduction_name %in% names(seurat_obj@reductions)) {
    stop(paste("ERROR: Reduction", reduction_name, "not found."))
  }
  emb <- Embeddings(seurat_obj, reduction = reduction_name)
  
  # 2. Safely grab the correct PC column
  pc_col <- paste0(reduction_name, "_", pc_axis)
  if (!pc_col %in% colnames(emb)) {
    pc_col <- colnames(emb)[pc_axis]
  }
  
  # 3. Build the base dataframe
  df <- data.frame(
    cell = rownames(emb),
    PC_value = emb[, pc_col],
    cell_type = seurat_obj@meta.data[rownames(emb), celltype_col]
  )
  
  # Remove missing cell types
  df <- df %>% filter(!is.na(cell_type))
  
  # 4. Bin by COORDINATE range (Equal Width)
  # cut() breaks the range from min(PC_value) to max(PC_value) into equal chunks
  df <- df %>%
    mutate(
      PC_bin = cut(PC_value, breaks = n_bins, labels = FALSE),
      
      # Optional: Capture the exact coordinate ranges as text labels for clarity
      PC_range_label = cut(PC_value, breaks = n_bins)
    )
  
  # 5. Assign human-readable names to the coordinate zones
  if (n_bins == 4) {
    df <- df %>%
      mutate(
        PC_region = case_when(
          PC_bin == 1 ~ paste0("Lowest Range PC", pc_axis),
          PC_bin == 2 ~ paste0("Low-mid Range PC", pc_axis),
          PC_bin == 3 ~ paste0("High-mid Range PC", pc_axis),
          PC_bin == 4 ~ paste0("Highest Range PC", pc_axis)
        )
      )
  } else {
    df <- df %>%
      mutate(
        PC_region = paste0("PC", pc_axis, " Coordinate Bin ", PC_bin)
      )
  }
  
  # 6. Calculate cell-type percentages inside each spatial bin
  percent_df <- df %>%
    mutate(condition = condition_label) %>%
    dplyr::count(condition, PC_region, PC_range_label, cell_type) %>%
    group_by(condition, PC_region) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    ungroup() %>%
    mutate(PC_axis = paste0("PC", pc_axis))

  return(percent_df)
}

pc1_percent <- make_pc_coordinate_bins(
  seurat_obj = sc_data_Upregulated,
  reduction_name = "pca_Upregulated_deg",
  condition_label = "Upregulated_DEGs",
  celltype_col = "celltype_grouped",
  pc_axis = 1,
  n_bins = 4
)
pc1_percent$PC_region <- factor(
  pc1_percent$PC_region,
  levels = c(
    "Lowest Range PC1",
    "Low-mid Range PC1",
    "High-mid Range PC1",
    "Highest Range PC1"
  )
)
PC1 <- ggplot(
  pc1_percent,
  aes(x = PC_region, y = percent, fill = cell_type)
) +
  scale_fill_manual(values = cols, drop = FALSE) +
  geom_bar(stat = "identity") +
  facet_wrap(~PC_axis, nrow = 1) +
  theme_bw() +
  labs(
    title = "Cell-type composition across PCA coordinate bins Up reg genes",
    x = "Coordinate bin",
    y = "Percent of cells"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

PC1
ggsave(
  filename = file.path(dir_figures, "Up_PC1_celltype_composition.png"),
  plot = PC1,
  width = 8,
  height = 6,
  dpi = 300
)

pc2_percent <- make_pc_coordinate_bins(
  seurat_obj = sc_data_Upregulated,
  reduction_name = "pca_Upregulated_deg",
  condition_label = "Upregulated_DEGs",
  celltype_col = "celltype_grouped",
  pc_axis = 2,
  n_bins = 4
)
pc2_percent$PC_region <- factor(
  pc2_percent$PC_region,
  levels = c(
    "Lowest Range PC2",
    "Low-mid Range PC2",
    "High-mid Range PC2",
    "Highest Range PC2"
  )
)

PC2 <- ggplot(
  pc2_percent,
  aes(x = PC_region, y = percent, fill = cell_type)
) +
  scale_fill_manual(values = cols, drop = FALSE) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(
    title = "Cell-type composition across PC2 coordinate bins Up reg genes",
    x = "Coordinate bin",
    y = "Percent of cells"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

PC2

ggsave(
  file.path(dir_figures, "Up_PC2_celltype_composition.png"),
  PC2,
  width = 10,
  height = 6,
  dpi = 300
)
# -----------------------------------------------------------------------------
# 16. Upregulated correlation heatmap
# -----------------------------------------------------------------------------
Upregulated_matrix <- as.matrix(
  GetAssayData(sc_data_Upregulated, layer = "data")[Upregulated_genes, ]
)

Upregulated_variance <- apply(Upregulated_matrix, 1, var)

Upregulated_valid_genes <- names(Upregulated_variance[Upregulated_variance > 0])

Upregulated_matrix_valid <- as.matrix(
  GetAssayData(sc_data_Upregulated, layer = "data")[Upregulated_valid_genes, ]
)

Upregulated_correlations <- cor(
  (t(Upregulated_matrix_valid) + 1),
  method = "pearson"
)

write.csv(
  Upregulated_correlations,
  file.path(dir_tables, "Upregulated_correlations.csv")
)

# -----------------------------------------------------------------------------
# Hierarchical clustering
# -----------------------------------------------------------------------------

gene_distance <- as.dist(1 - Upregulated_correlations)

gene_clustering <- hclust(
  gene_distance,
  method = "average"
)

# -----------------------------------------------------------------------------
# Dynamic tree cutting
# -----------------------------------------------------------------------------

dynamic_modules <- cutreeDynamic(
  dendro = gene_clustering,
  distM = as.matrix(gene_distance),
  
  deepSplit = 1,
  minClusterSize = 5
)

# -----------------------------------------------------------------------------
# Inspect module sizes
# -----------------------------------------------------------------------------

table(dynamic_modules)

# -----------------------------------------------------------------------------
# Create gene-module table
# -----------------------------------------------------------------------------

gene_module_table <- data.frame(
  Gene = rownames(Upregulated_correlations),
  Module = dynamic_modules
)

# -----------------------------------------------------------------------------
# Save module assignments
# -----------------------------------------------------------------------------

write.csv(gene_module_table, file.path(dir_tables, "Upregulated_gene_modules.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# Create module annotation for heatmap in the oreder of the modules + order correlation matrix
# -----------------------------------------------------------------------------

annotation_df <- data.frame(
  Module = factor(dynamic_modules)
)

rownames(annotation_df) <- rownames(Upregulated_correlations)


# -----------------------------------------------------------------------------
# Plot heatmap
# -----------------------------------------------------------------------------
# 1. Define the scale limits and center
min_val <- min(
  Upregulated_correlations,
  na.rm = TRUE
)

max_val <- max(
  Upregulated_correlations,
  na.rm = TRUE
)

center_val <- 0

# Inspect values
min_val
max_val

# 2. Create asymmetric breaks
# Ensure the number of colors on each side of the center is proportional to the range
palette_length <- 100
my_breaks <- c(seq(min_val, center_val, length.out = ceiling(palette_length/2) + 1), 
               seq(max_val/palette_length, max_val, length.out = floor(palette_length/2)))

# 3. Define colors
my_colors <- colorRampPalette(c("#0571b0", "#f7f7f7", "#ca0020"))(palette_length)


p_Up<-pheatmap(
  Upregulated_correlations,
  
  cluster_rows = gene_clustering,
  cluster_cols = gene_clustering,
  
  annotation_row = annotation_df,
  annotation_col = annotation_df,
  
  clustering_method = "average",
  
  color=my_colors,
  breaks = my_breaks,
  
  border_color = NA,

  fontsize_row = 4,
  fontsize_col = 4,

  # pheatmap writes the file itself (ggsave() does not work on a pheatmap object).
  filename = file.path(dir_figures, "Up_heatmap_gene_gene_correlations.png"),
  width = 20,
  height = 20
)
# =============================================================================
# FINAL SUMMARY
# =============================================================================

message("Analysis complete.")
message("Figures saved in: ", dir_figures)
message("Tables saved in:  ", dir_tables)