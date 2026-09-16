# ============================================================
# Title: Single-cell transcriptomic analysis of human retinal
#        pericytes, UROD virtual knockout and cell-cell communication
# R version: >= 4.3.0
# ============================================================

rm(list = ls())
set.seed(123)
options(stringsAsFactors = FALSE)

project_dir <- "."
sam.name    <- "multi"
out_dir     <- file.path(project_dir, sam.name)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------
required_pkgs <- c(
  "Seurat", "ggplot2", "cowplot", "Matrix", "dplyr", "ggsci",
  "harmony", "CellChat", "patchwork", "ggpubr", "scales",
  "clusterProfiler", "org.Hs.eg.db", "fgsea", "msigdbr",
  "tidytext", "forcats", "igraph", "ggrepel"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
}

library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(harmony)
library(CellChat)
library(patchwork)
library(ggpubr)
library(scales)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)
library(msigdbr)
library(tidytext)
library(forcats)
library(igraph)
library(ggrepel)

# scTenifoldKnk is not on CRAN; install locally if needed.
# install.packages("C:/R_packages/scTenifoldKnk-1.0.3", repos = NULL, type = "source")
if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
  stop("Package 'scTenifoldKnk' is required. Please install it first.")
}
library(scTenifoldKnk)

# ------------------------------------------------------------
# 1. Colors
# ------------------------------------------------------------
cluster_cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"
)

# ============================================================
# 2. Read data and preprocessing
# ============================================================

# 2.1 Read 10x data
DR1 <- Read10X(file.path(project_dir, "DR1"))
DR2 <- Read10X(file.path(project_dir, "DR2"))
DR3 <- Read10X(file.path(project_dir, "DR3"))
DR4 <- Read10X(file.path(project_dir, "DR4"))

# If Read10X returns a list, extract gene expression matrix
if (is.list(DR1)) DR1 <- DR1[["Gene Expression"]]
if (is.list(DR2)) DR2 <- DR2[["Gene Expression"]]
if (is.list(DR3)) DR3 <- DR3[["Gene Expression"]]
if (is.list(DR4)) DR4 <- DR4[["Gene Expression"]]

# Add sample suffix to barcodes
colnames(DR1) <- paste(colnames(DR1), "DR1", sep = "_")
colnames(DR2) <- paste(colnames(DR2), "DR2", sep = "_")
colnames(DR3) <- paste(colnames(DR3), "DR3", sep = "_")
colnames(DR4) <- paste(colnames(DR4), "DR4", sep = "_")

# Merge all samples
experiment.data <- cbind(DR1, DR2, DR3, DR4)

# Create Seurat object
experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project       = "multi",
  min.cells     = 10,
  min.features  = 200
)

# Set sample identity from barcode suffix
experiment.aggregate$orig.ident <- sub(".*_", "", colnames(experiment.aggregate))
Idents(experiment.aggregate) <- experiment.aggregate$orig.ident

# 2.2 Mitochondrial percentage
experiment.aggregate[["percent.MT"]] <- PercentageFeatureSet(
  experiment.aggregate,
  pattern = "^MT"
)

# QC violin plot
pdf(file.path(out_dir, "QC-VlnPlot.pdf"), width = 15, height = 5)
print(
  VlnPlot(
    experiment.aggregate,
    features = c("nFeature_RNA", "nCount_RNA", "percent.MT"),
    ncol = 3
  )
)
dev.off()

# 2.3 QC filtering
cat("Before filter:", nrow(experiment.aggregate@meta.data), "cells\n")

experiment.aggregate <- subset(
  experiment.aggregate,
  subset =
    nFeature_RNA > 50 &
    nFeature_RNA < 6000 &
    nCount_RNA > 10 &
    nCount_RNA < 20000 &
    percent.MT < 10
)

cat("After filter:", nrow(experiment.aggregate@meta.data), "cells\n")

# 2.4 Normalization
experiment.aggregate <- NormalizeData(
  experiment.aggregate,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# 2.5 Variable features
experiment.aggregate <- FindVariableFeatures(
  experiment.aggregate,
  selection.method = "vst",
  nfeatures = 2000
)

# 2.6 Scaling
experiment.aggregate <- ScaleData(
  experiment.aggregate,
  features = VariableFeatures(experiment.aggregate)
)

# 2.7 PCA
experiment.aggregate <- RunPCA(
  experiment.aggregate,
  features = VariableFeatures(experiment.aggregate),
  verbose = FALSE,
  npcs = 50
)

# 2.8 Harmony batch correction
experiment.aggregate <- RunHarmony(
  experiment.aggregate,
  reduction = "pca",
  group.by.vars = "orig.ident",
  reduction.save = "harmony"
)

# 2.9 tSNE and UMAP
dim.use <- 1:30

experiment.aggregate <- RunTSNE(
  experiment.aggregate,
  reduction = "harmony",
  dims = dim.use,
  reduction.name = "tsne"
)

experiment.aggregate <- RunUMAP(
  experiment.aggregate,
  reduction = "harmony",
  dims = 1:30
)

# 2.10 Clustering
experiment.aggregate <- FindNeighbors(
  experiment.aggregate,
  reduction = "harmony",
  dims = 1:30
)

experiment.aggregate <- FindClusters(
  experiment.aggregate,
  resolution = 0.5
)

# 2.11 Marker genes
all.markers <- FindAllMarkers(
  experiment.aggregate,
  only.pos = TRUE,
  min.pct = 0.3,
  logfc.threshold = 0.25
)

top10.markers <- all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  ungroup()

write.csv(
  top10.markers,
  file = file.path(out_dir, "Supplementary_Table_Top10_Markers.csv"),
  row.names = FALSE
)

write.table(
  all.markers,
  file = file.path(out_dir, "total_marker_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 2.12 Cell type annotation
celltype <- data.frame(
  ClusterID = 0:16,
  celltype  = "unknown"
)

celltype[celltype$ClusterID %in% c(0, 7, 13), 2] <- "Macrophage"
celltype[celltype$ClusterID %in% c(1, 10), 2]    <- "Mural cell"
celltype[celltype$ClusterID %in% c(2, 6), 2]     <- "Endothelial cell"
celltype[celltype$ClusterID %in% c(3, 9), 2]     <- "T cell"
celltype[celltype$ClusterID %in% c(4, 5, 14), 2] <- "Fibroblast"
celltype[celltype$ClusterID %in% c(8), 2]        <- "Monocyte"
celltype[celltype$ClusterID %in% c(11), 2]       <- "cDC"
celltype[celltype$ClusterID %in% c(12, 15), 2]   <- "B cell"

experiment.aggregate@meta.data$celltype <- "NA"

for (i in seq_len(nrow(celltype))) {
  idx <- which(experiment.aggregate@meta.data$seurat_clusters == celltype$ClusterID[i])
  experiment.aggregate@meta.data$celltype[idx] <- celltype$celltype[i]
}

new.cluster.ids <- c(
  "Macrophage", "Mural cell", "Endothelial cell",
  "T cell", "Fibroblast", "Fibroblast",
  "Endothelial cell", "Macrophage", "Monocyte",
  "T cell", "Mural cell", "cDC",
  "B cell", "Macrophage", "Fibroblast",
  "B cell", "unknown"
)

names(new.cluster.ids) <- levels(experiment.aggregate)
experiment.aggregate <- RenameIdents(experiment.aggregate, new.cluster.ids)

scobj <- experiment.aggregate

all.markers <- FindAllMarkers(
  scobj,
  only.pos = TRUE,
  min.pct = 0.3,
  logfc.threshold = 0.25
)

top_markers <- all.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup()

write.csv(
  top_markers,
  file = file.path(out_dir, "top_markers_per_celltype.csv"),
  row.names = FALSE
)

# ============================================================
# 3. Mural cell subcluster identification
# ============================================================

mural_cells <- subset(scobj, subset = celltype == "Mural cell")

mural_cells <- NormalizeData(mural_cells) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony(
    group.by.vars = "orig.ident",
    reduction = "pca",
    reduction.save = "harmony"
  ) %>%
  RunTSNE(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.5)

# Marker genes for pericytes, VSMC and endothelial contamination
pericyte_pos  <- c("PDGFRB", "CSPG4", "RGS5", "DLK1", "ANPEP")
smooth_muscle <- c("MYH11", "CNN1")
endothelial   <- c("PECAM1", "CDH5")

pdf(file.path(out_dir, "Mural_subcluster_markers.pdf"), width = 14, height = 10)
print(
  FeaturePlot(
    mural_cells,
    features = c(pericyte_pos, smooth_muscle, endothelial),
    reduction = "tsne",
    ncol = 4
  )
)
dev.off()

print(table(Idents(mural_cells)))

pdf(file.path(out_dir, "Mural_subcluster_tsne_markers.pdf"), width = 12, height = 6)
print(DimPlot(mural_cells, reduction = "tsne", label = TRUE))
print(
  DotPlot(
    mural_cells,
    features = c("PDGFRB", "CSPG4", "RGS5", "MYH11", "CNN1")
  ) + RotatedAxis()
)
dev.off()

# ------------------------------------------------------------
# IMPORTANT: manually check cluster numbers from DotPlot
# ------------------------------------------------------------
pericyte_clusters          <- c(0, 1, 4)
vsmc_clusters              <- c(3)
transitional_clusters      <- c(2)

pericyte_barcodes <- colnames(subset(mural_cells, idents = pericyte_clusters))
vsmc_barcodes     <- colnames(subset(mural_cells, idents = vsmc_clusters))
trans_barcodes    <- colnames(subset(mural_cells, idents = transitional_clusters))
mural_other_barcodes <- colnames(subset(mural_cells, idents = mural_unclassified_clusters))

scobj$celltype_fine <- as.character(scobj$celltype)
scobj$celltype_fine[pericyte_barcodes] <- "Pericyte"
scobj$celltype_fine[vsmc_barcodes]     <- "VSMC"
scobj$celltype_fine[trans_barcodes]    <- "Transitional"
scobj$celltype_fine[mural_other_barcodes] <- "Mural_Unclassified"

Idents(scobj) <- factor(scobj$celltype_fine)

print(table(scobj$celltype_fine))

pdf(file.path(out_dir, "Pericyte_VSMC_validation.pdf"), width = 8, height = 5)
print(
  DotPlot(
    subset(scobj, subset = celltype_fine %in% c("Pericyte", "VSMC")),
    features = c("PDGFRB", "CSPG4", "RGS5", "MYH11", "CNN1")
  ) + RotatedAxis()
)
dev.off()

# ------------------------------------------------------------
# Remove low-quality / unclassified cells
# ------------------------------------------------------------
scobj_clean <- subset(
  scobj,
  subset = celltype_fine %in% c(
    "Pericyte", "VSMC", "Transitional",
    "B cell", "cDC", "Endothelial cell",
    "Fibroblast", "Macrophage", "Monocyte", "T cell"
  )
)

Idents(scobj) <- "celltype_fine"

pdf(file.path(out_dir, "Mural_subtypes_tsne.pdf"), width = 8, height = 6)
print(
  DimPlot(
    subset(scobj, idents = c("Pericyte", "VSMC", "Transitional")),
    reduction = "tsne",
    label = TRUE,
    label.size = 4
  ) + ggtitle("Mural Subtypes")
)
dev.off()

scobj_sub <- subset(
  scobj,
  subset = celltype_fine != "Mural_Unclassified" & celltype_fine != "unknown"
)

scobj_sub$celltype_fine <- factor(
  scobj_sub$celltype_fine,
  levels = c(
    "Pericyte", "VSMC", "Transitional",
    "B cell", "cDC", "Endothelial cell",
    "Fibroblast", "Macrophage", "Monocyte", "T cell"
  )
)

Idents(scobj_sub) <- "celltype_fine"

# ============================================================
# 4. UROD expression visualization
# ============================================================

p1 <- DimPlot(
  scobj_sub,
  reduction = "umap",
  label = TRUE,
  label.size = 4,
  pt.size = 0.8,
  cols = cluster_cols
) +
  NoLegend() +
  ggtitle("Global Cell Clusters")

p2 <- FeaturePlot(
  scobj_sub,
  features = "UROD",
  reduction = "umap",
  min.cutoff = "q10",
  max.cutoff = "q90",
  pt.size = 0.8
) +
  ggtitle("UROD Expression")

p3_highlight <- DotPlot(
  scobj_sub,
  features = "UROD",
  group.by = "celltype_fine"
) +
  RotatedAxis() +
  scale_color_gradientn(colors = c("lightgrey", "blue")) +
  ggtitle("UROD: Positive Cell Proportion") +
  theme(axis.text.x = element_text(size = 10, face = "bold"))

p4_fixed <- VlnPlot(
  scobj_sub,
  features = "UROD",
  group.by = "celltype_fine",
  pt.size = 0,
  cols = cluster_cols
) +
  scale_y_continuous(limits = c(0, 1.5), oob = scales::squish) +
  stat_compare_means(
    comparisons = list(
      c("Pericyte", "Fibroblast"),
      c("Pericyte", "Macrophage"),
      c("Pericyte", "T cell")
    ),
    label = "p.signif",
    method = "wilcox.test"
  ) +
  theme(legend.position = "none") +
  ggtitle("UROD Expression (Scale optimized)")

p5 <- plot_grid(
  p1, p2, p3_highlight, p4_fixed,
  labels = c("A", "B", "C", "D"),
  ncol = 2
)

ggsave(
  file.path(out_dir, "Final_UROD_4Panel_Plots.pdf"),
  p5,
  width = 16,
  height = 12
)

# ============================================================
# 5. Virtual knockout with scTenifoldKnk
# ============================================================

mural_cells <- subset(scobj, subset = celltype == "Mural cell")

mural_cells <- NormalizeData(mural_cells) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony(
    group.by.vars = "orig.ident",
    reduction = "pca",
    reduction.save = "harmony"
  ) %>%
  RunTSNE(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.5)

print(table(Idents(mural_cells)))

pdf(file.path(out_dir, "Mural_subcluster_identification.pdf"), width = 14, height = 10)
print(DimPlot(mural_cells, reduction = "tsne", label = TRUE))
print(
  DotPlot(
    mural_cells,
    features = c(
      "PDGFRB", "CSPG4", "RGS5", "DLK1", "ANPEP",
      "MYH11", "CNN1",
      "PECAM1", "CDH5"
    )
  )
)
dev.off()

# ------------------------------------------------------------
# IMPORTANT: adjust cluster numbers if necessary
# ------------------------------------------------------------
mural_cells$celltype_raw <- ifelse(
  mural_cells$seurat_clusters %in% c(0, 1, 4), "Pericyte",
  ifelse(mural_cells$seurat_clusters %in% c(2, 3), "VSMC", "others")
)

pure_pericyte_cells <- colnames(mural_cells)[
  mural_cells$celltype_raw == "Pericyte" &
    FetchData(mural_cells, vars = "PDGFRB")[, 1] > 1.5 &
    FetchData(mural_cells, vars = "MYH11")[, 1] == 0 &
    FetchData(mural_cells, vars = "CNN1")[, 1] == 0
]

pure_vsmc_cells <- colnames(mural_cells)[mural_cells$celltype_raw == "VSMC"]

mural_cells$celltype_fine <- "others"
mural_cells$celltype_fine[pure_pericyte_cells] <- "Pericyte"
mural_cells$celltype_fine[pure_vsmc_cells] <- "VSMC"

print(table(mural_cells$celltype_fine))

# Write back to scobj
if ("celltype_fine" %in% colnames(scobj@meta.data)) {
  scobj$celltype_fine <- NULL
}
scobj$celltype_fine <- "others"
scobj$celltype_fine[colnames(mural_cells)] <- mural_cells$celltype_fine
print(table(scobj$celltype_fine))

# Ensure UMAP exists
if (!"umap" %in% names(scobj@reductions)) {
  message("UMAP not found in scobj. Computing UMAP...")
  if (!"harmony" %in% names(scobj@reductions)) {
    scobj <- NormalizeData(scobj) %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA() %>%
      RunHarmony(group.by.vars = "orig.ident", reduction = "pca") %>%
      RunUMAP(reduction = "harmony", dims = 1:15)
  } else {
    scobj <- RunUMAP(scobj, reduction = "harmony", dims = 1:15)
  }
}

cluster_cols_fine <- c(
  "Pericyte"     = "#E41A1C",
  "VSMC"         = "#377EB8",
  "Transitional" = "#4DAF4A",
  "others"       = "grey80"
)

# ------------------------------------------------------------
# Extract pure pericytes
# ------------------------------------------------------------
pbmc_pure <- subset(
  scobj,
  subset = celltype_fine == "Pericyte" & MYH11 == 0 & CNN1 == 0
)

n_pericyte_clean <- ncol(pbmc_pure)
message("Final pure pericyte number: ", n_pericyte_clean)

pdf(file.path(out_dir, "Pericyte_purification_check.pdf"), width = 12, height = 10)

print(
  DimPlot(
    scobj,
    reduction = "umap",
    group.by = "celltype_fine",
    cols = cluster_cols_fine,
    order = TRUE
  ) + ggtitle("Fine cell subtypes UMAP")
)

print(
  DimPlot(
    pbmc_pure,
    reduction = "umap",
    cols = cluster_cols_fine["Pericyte"]
  ) + ggtitle("Pure pericytes")
)

print(
  VlnPlot(
    scobj,
    features = c("MYH11", "CNN1"),
    group.by = "celltype_fine",
    pt.size = 0,
    cols = cluster_cols_fine
  ) + plot_annotation(title = "Smooth muscle markers before purification")
)

print(
  VlnPlot(
    pbmc_pure,
    features = c("MYH11", "CNN1"),
    group.by = "celltype_fine",
    pt.size = 0,
    cols = cluster_cols_fine["Pericyte"]
  ) + plot_annotation(title = "Smooth muscle markers after purification")
)

print(
  VlnPlot(
    scobj,
    features = c("PDGFRB", "RGS5", "ACTA2", "TAGLN"),
    group.by = "celltype_fine",
    pt.size = 0,
    cols = cluster_cols_fine
  ) + plot_annotation(title = "Pericyte / activation markers")
)

print(
  VlnPlot(
    pbmc_pure,
    features = c("PDGFRB", "RGS5", "ACTA2", "TAGLN"),
    group.by = "celltype_fine",
    pt.size = 0,
    cols = cluster_cols_fine["Pericyte"]
  ) + plot_annotation(title = "Pure pericyte markers")
)

print(
  VlnPlot(
    scobj,
    features = "DES",
    group.by = "celltype_fine",
    pt.size = 0,
    cols = cluster_cols_fine
  ) + ggtitle("DES expression across fine subtypes")
)

print(
  VlnPlot(pbmc_pure, features = "DES", pt.size = 0) +
    ggtitle("DES expression in pure pericytes")
)

dev.off()

# ------------------------------------------------------------
# Prepare count matrix for scTenifoldKnk
# ------------------------------------------------------------
countMat <- tryCatch(
  GetAssayData(pbmc_pure, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(pbmc_pure, assay = "RNA", slot = "counts")
)

pbmc_pure <- FindVariableFeatures(
  pbmc_pure,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 10000
)

hvgs <- VariableFeatures(pbmc_pure)
target_genes <- unique(c("UROD", hvgs))
count_matrix <- as.matrix(countMat[target_genes, ])

message("Final pure pericyte number: ", ncol(pbmc_pure))

print(
  VlnPlot(
    pbmc_pure,
    features = c("PDGFRB", "RGS5", "MYH11", "CNN1", "ACTA2", "TAGLN"),
    pt.size = 0,
    ncol = 3
  )
)

print(table(FetchData(pbmc_pure, vars = "MYH11")[, 1] > 0))

# ------------------------------------------------------------
# Run virtual knockout
# ------------------------------------------------------------
result_knockout <- scTenifoldKnk(
  countMatrix     = count_matrix,
  gKO             = "UROD",
  qc              = TRUE,
  qc_mtThreshold  = 0.15,
  qc_minLSize     = 500,
  nc_nNet         = 3,
  nc_nCells       = min(250, ncol(count_matrix))
)

saveRDS(
  result_knockout,
  file = file.path(out_dir, "scTenifoldKnk_UROD_results.rds")
)

print(names(result_knockout))

urod_exp <- count_matrix["UROD", ]
print(summary(urod_exp))
print(sum(urod_exp > 0) / length(urod_exp))

print(class(result_knockout$diffRegulation))
print(dim(result_knockout$diffRegulation))
print(head(result_knockout$diffRegulation))

sig_genes <- result_knockout$diffRegulation[
  result_knockout$diffRegulation$Z > 2,
]
print(head(sig_genes, 20))

sig_genes_adj <- subset(
  result_knockout$diffRegulation,
  p.adj < 0.05
)
print(nrow(sig_genes_adj))

df_pure <- result_knockout$diffRegulation
df_pure <- df_pure[df_pure$gene != "UROD", ]
df_pure_sorted <- df_pure[order(df_pure$p.adj), ]
print(head(df_pure_sorted, 20))

df <- result_knockout$diffRegulation
df_noUROD <- df[df$gene != "UROD", ]
candidate <- head(df_noUROD[order(-df_noUROD$Z), ], 20)

df$log10p <- -log10(df$p.value + 1e-300)
top_genes <- candidate$gene
df$label_me <- df$gene %in% top_genes

p_differential <- ggplot(df, aes(x = Z, y = log10p)) +
  geom_point(aes(color = label_me), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_text_repel(
    data = subset(df, label_me),
    aes(label = gene),
    size = 3.5,
    max.overlaps = 50
  ) +
  labs(x = "Z-score", y = "-log10(p-value)") +
  theme_minimal() +
  ggtitle("Differential regulation upon UROD knockout") +
  theme(legend.position = "none")

ggsave(
  file.path(out_dir, "Differential_regulation_upon_UROD_knockout.pdf"),
  p_differential,
  width = 12,
  height = 10
)

top15 <- head(df_noUROD[order(-df_noUROD$Z), ], 15)

p_top <- ggplot(top15, aes(x = reorder(gene, Z), y = Z)) +
  geom_col(aes(fill = Z), show.legend = FALSE) +
  scale_fill_gradient2(
    low = "#D6EAF8",
    mid = "#85C1E9",
    high = "#2E86C1",
    midpoint = median(top15$Z)
  ) +
  coord_flip() +
  labs(
    x = "",
    y = "Z-score",
    title = "Top 15 perturbed genes (excluding UROD)"
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir, "Top15_perturbed_genes.pdf"),
  p_top,
  width = 12,
  height = 10
)

# ============================================================
# 6. Enrichment analysis
# ============================================================

df <- result_knockout$diffRegulation
genes_ora <- df$gene[df$gene != "UROD" & abs(df$Z) > 1.5]

geneList <- df$Z
names(geneList) <- df$gene
geneList <- sort(geneList, decreasing = TRUE)

# GO ORA
ego_bp <- enrichGO(
  gene         = genes_ora,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2
)

ego_mf <- enrichGO(
  gene         = genes_ora,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "MF",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2
)

ego_cc <- enrichGO(
  gene         = genes_ora,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "CC",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2
)

bp_df <- as.data.frame(ego_bp)
mf_df <- as.data.frame(ego_mf)
cc_df <- as.data.frame(ego_cc)

bp_df$ONTOLOGY <- "BP"
mf_df$ONTOLOGY <- "MF"
cc_df$ONTOLOGY <- "CC"

all_go <- rbind(bp_df, mf_df, cc_df)

top_go <- all_go %>%
  group_by(ONTOLOGY) %>%
  slice_min(order_by = p.adjust, n = 5) %>%
  ungroup()

p_go <- ggplot(
  top_go,
  aes(
    x = -log10(p.adjust),
    y = reorder_within(Description, p.adjust, ONTOLOGY),
    fill = p.adjust
  )
) +
  geom_col(width = 0.7) +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
  scale_y_reordered() +
  scale_fill_gradient2(
    low = "#D96B7B",
    mid = "#D96B7B",
    high = "#2166AC",
    midpoint = median(top_go$p.adjust),
    name = "p.adjust"
  ) +
  labs(
    x = "-log10(adjusted p-value)",
    y = NULL,
    title = "GO Enrichment (BP, MF, CC)"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black"),
    legend.position = "right"
  )

ggsave(
  file.path(out_dir, "GO_ORA2.pdf"),
  p_go,
  width = 8,
  height = 10
)

write.csv(
  all_go,
  file = file.path(out_dir, "GO_ORA_results.csv"),
  row.names = FALSE
)

# GO GSEA
ego_gsea <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  pvalueCutoff = 0.2,
  verbose      = FALSE
)

write.csv(
  as.data.frame(ego_gsea),
  file = file.path(out_dir, "GO_BP_GSEA_results.csv"),
  row.names = FALSE
)

pdf(file.path(out_dir, "ego_gsea1.pdf"), width = 8, height = 10)
print(ridgeplot(ego_gsea, showCategory = 10) + ggtitle("GO BP (GSEA)"))
dev.off()

# KEGG / WikiPathways fgsea
msig_kegg <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CP:KEGG"
)

msig_wiki <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CP:WIKIPATHWAYS"
)

pathway_list <- c(
  split(msig_kegg$gene_symbol, msig_kegg$gs_name),
  split(msig_wiki$gene_symbol, msig_wiki$gs_name)
)

gene_list_sym <- df$Z
names(gene_list_sym) <- df$gene
gene_list_sym <- sort(gene_list_sym, decreasing = TRUE)

gsea_result <- fgsea(
  pathways = pathway_list,
  stats    = gene_list_sym,
  minSize  = 10,
  maxSize  = 500,
  nperm    = 10000
)

kegg_fgsea <- gsea_result[grepl("^KEGG_", pathway), ]
wiki_fgsea <- gsea_result[grepl("^WP_", pathway), ]

print(head(kegg_fgsea[order(pval), ], 20))
print(head(wiki_fgsea[order(pval), ], 20))

top <- gsea_result %>%
  filter(padj < 0.25) %>%
  arrange(padj) %>%
  head(20)

P5 <- ggplot(
  top,
  aes(
    x = NES,
    y = fct_reorder(pathway, NES),
    size = -log10(padj),
    color = NES
  )
) +
  geom_point() +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(
    x = "Normalized Enrichment Score",
    y = "",
    size = expression(-log[10](padj)),
    color = "NES"
  ) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank()) +
  ggtitle("GSEA (KEGG & WikiPathways)")

ggsave(
  file.path(out_dir, "ego_kegg.pdf"),
  P5,
  width = 8,
  height = 10
)

# ============================================================
# 7. Cell-cell communication analysis with CellChat
# ============================================================

source_cell <- "Pericyte"
target_t    <- "T cell"
target_ec   <- "Endothelial cell"

pericyte_cells <- colnames(pbmc_pure)

scobj$celltype <- as.character(scobj$celltype)
scobj$celltype[pericyte_cells] <- "Pericyte"

mural_all <- WhichCells(scobj, expression = celltype == "Mural cell")
non_pericyte_mural <- setdiff(mural_all, pericyte_cells)
scobj$celltype[non_pericyte_mural] <- "mural_other"

print(table(scobj$celltype))

cellchat <- createCellChat(
  object = scobj,
  meta = scobj@meta.data,
  group.by = "celltype"
)

cellchat <- setIdent(cellchat, ident.use = "celltype")
print(levels(cellchat@idents))

groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(
  cellchat,
  file = file.path(out_dir, "cellchat_object.rds")
)

# Global interaction overview
pdf(file.path(out_dir, "Global_Interaction_Overview.pdf"), width = 12, height = 5)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions"
)
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction weights/strength"
)
dev.off()

# Each cell type as sender
pdf(file.path(out_dir, "Each_Cell_as_Sender_Circle.pdf"), width = 12, height = 10)
mat <- cellchat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}
dev.off()

# Bubble plots
p_bubble_T <- netVisual_bubble(
  cellchat,
  sources.use = "Pericyte",
  targets.use = "T cell",
  remove.isolate = FALSE,
  title.name = "Pericyte -> T cell"
)

ggsave(
  file.path(out_dir, "Pericyte_to_T_bubble.pdf"),
  p_bubble_T,
  width = 8,
  height = 6
)

p_bubble_T2 <- netVisual_bubble(
  cellchat,
  sources.use = "Pericyte",
  targets.use = "Endothelial cell",
  remove.isolate = FALSE,
  title.name = "Pericyte -> Endothelial cell"
)

ggsave(
  file.path(out_dir, "Pericyte_to_EC_bubble.pdf"),
  p_bubble_T2,
  width = 8,
  height = 6
)

# Chord diagrams: Pericyte -> T cell
pathways_T <- c("MIF", "MK")

for (pathway in pathways_T) {
  if (!pathway %in% cellchat@netP$pathways) {
    cat("Pathway", pathway, "not found. Skipping.\n")
    next
  }
  pdf(
    file.path(out_dir, paste0("Pathway_", pathway, "_Pericyte_to_T_Chord.pdf")),
    width = 8,
    height = 7
  )
  netVisual_chord_cell(
    cellchat,
    signaling = pathway,
    sources.use = "Pericyte",
    targets.use = "T cell",
    lab.cex = 0.8,
    title.name = paste0(pathway, " (Pericyte -> T cell)")
  )
  dev.off()
}

# Chord diagrams: Pericyte -> Endothelial cell
pathways_EC <- c("VEGF", "CypA", "MK")

for (pathway in pathways_EC) {
  if (!pathway %in% cellchat@netP$pathways) {
    cat("Pathway", pathway, "not found. Skipping.\n")
    next
  }
  pdf(
    file.path(out_dir, paste0("Pathway_", pathway, "_Pericyte_to_EC_Chord.pdf")),
    width = 8,
    height = 7
  )
  netVisual_chord_cell(
    cellchat,
    signaling = pathway,
    sources.use = source_cell,
    targets.use = target_ec,
    lab.cex = 0.8,
    title.name = paste0(pathway, " (", source_cell, " -> ", target_ec, ")")
  )
  dev.off()
}

# ============================================================
# 8. Validation of ligand-receptor gene expression
# ============================================================

genes_all <- c(
  "MIF",
  "CD74", "CXCR4", "CD44",
  "MDK",
  "NCL", "ITGA4", "ITGB1",
  "VEGFA", "VEGFB", "PGF",
  "KDR", "FLT1", "NRP1", "NRP2",
  "PPIA",
  "BSG", "PTPRC"
)

p_dot <- DotPlot(
  scobj,
  features = genes_all,
  group.by = "celltype"
) +
  scale_color_gradientn(colours = c("lightgrey", "#1E90FF", "#D22B2B")) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11, color = "black")
  ) +
  labs(
    x = "Ligand/Receptor Genes",
    y = "Cell Types",
    title = "Expression of key pathway genes (MIF, MK, VEGF, CypA)"
  )

ggsave(
  file.path(out_dir, "DotPlot_AllPathways_Validation.pdf"),
  p_dot,
  width = 16,
  height = 7,
  dpi = 300
)

ligands_all <- c("MIF", "MDK", "VEGFA", "VEGFB", "PGF", "PPIA")
receptors_all <- c(
  "CD74", "CXCR4", "CD44", "NCL", "ITGA4", "ITGB1",
  "KDR", "FLT1", "NRP1", "NRP2", "BSG", "PTPRC"
)

p_lig <- VlnPlot(
  scobj,
  features = ligands_all,
  group.by = "celltype",
  pt.size = 0
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Ligands expression")

p_rec <- VlnPlot(
  scobj,
  features = receptors_all,
  group.by = "celltype",
  pt.size = 0
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Receptors expression")

ggsave(
  file.path(out_dir, "VlnPlot_Ligands_All.pdf"),
  p_lig,
  width = 14,
  height = 6
)

ggsave(
  file.path(out_dir, "VlnPlot_Receptors_All.pdf"),
  p_rec,
  width = 18,
  height = 6
)

p_feat_lig <- FeaturePlot(
  scobj,
  features = ligands_all,
  reduction = "tsne",
  cols = c("lightgrey", "red"),
  pt.size = 0.5,
  order = TRUE,
  label = TRUE,
  label.size = 4,
  ncol = 3
) +
  theme(plot.title = element_text(size = 12, face = "bold"))

ggsave(
  file.path(out_dir, "FeaturePlot_Ligands_All_tsne.pdf"),
  p_feat_lig,
  width = 14,
  height = 8,
  dpi = 300
)

receptors_focus <- c(
  "CXCR4", "CD44", "ITGA4", "ITGB1",
  "KDR", "FLT1", "NRP1",
  "BSG"
)

p_feat_rec <- FeaturePlot(
  scobj,
  features = receptors_focus,
  reduction = "tsne",
  cols = c("lightgrey", "blue"),
  pt.size = 0.5,
  order = TRUE,
  label = TRUE,
  label.size = 4,
  ncol = 3
) +
  theme(plot.title = element_text(size = 12, face = "bold"))

ggsave(
  file.path(out_dir, "FeaturePlot_Receptors_All_tsne.pdf"),
  p_feat_rec,
  width = 14,
  height = 10,
  dpi = 300
)

# Source and target cells only
cells_of_interest <- c(source_cell, target_t, target_ec)
scobj_sub_cc <- subset(scobj, subset = celltype %in% cells_of_interest)
Idents(scobj_sub_cc) <- "celltype"

genes_focus <- c(
  "MIF", "MDK", "VEGFA", "PPIA",
  "CD74", "CXCR4", "CD44", "NCL",
  "KDR", "FLT1", "BSG"
)

p_dot_sub <- DotPlot(
  scobj_sub_cc,
  features = genes_focus,
  group.by = "celltype"
) +
  scale_color_gradientn(colours = c("lightgrey", "#1E90FF", "#D22B2B")) +
  RotatedAxis() +
  theme_classic() +
  labs(
    x = "Ligand/Receptor Genes",
    y = "Cell Types",
    title = "Key pathway genes in source & target cells"
  )

ggsave(
  file.path(out_dir, "DotPlot_SourceTarget_All.pdf"),
  p_dot_sub,
  width = 12,
  height = 4,
  dpi = 300
)

# ============================================================
# 9. UROD-high vs UROD-low pericyte validation
# ============================================================

urod_exp <- FetchData(pbmc_pure, vars = "UROD")
threshold <- median(urod_exp$UROD)

pbmc_pure$UROD_status <- ifelse(
  urod_exp$UROD > threshold,
  "UROD_high",
  "UROD_low"
)

Idents(pbmc_pure) <- "UROD_status"

p_urod_lig <- VlnPlot(
  pbmc_pure,
  features = c("PGF", "PPIA", "MDK"),
  group.by = "UROD_status",
  pt.size = 0
) +
  stat_compare_means(
    comparisons = list(c("UROD_low", "UROD_high")),
    method = "wilcox.test"
  ) +
  scale_y_continuous(limits = c(0, 2.5), oob = scales::squish) +
  ggtitle("Endothelial Attacking Ligands in UROD KO Pericytes")

ggsave(
  file.path(out_dir, "UROD_high_low_ligands.pdf"),
  p_urod_lig,
  width = 8,
  height = 5
)

endo_sub <- subset(scobj_sub_cc, idents = "Endothelial cell")

p_endo <- VlnPlot(
  endo_sub,
  features = c("FLT1", "BSG", "ITGA4"),
  pt.size = 0
) +
  ggtitle("Endothelial receptors")

ggsave(
  file.path(out_dir, "Endothelial_receptors.pdf"),
  p_endo,
  width = 8,
  height = 5
)

# ============================================================
# 10. Session information
# ============================================================

writeLines(
  capture.output(sessionInfo()),
  con = file.path(out_dir, "sessionInfo.txt")
)

message("Analysis completed successfully.")