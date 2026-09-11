#######################################################################
# 第一部分：加载包
#######################################################################
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(harmony)
library(CellChat)   # 细胞通讯
library(patchwork)  # 拼图用

#######################################################################
# 第二部分：读取数据 & 预处理（基于你原来的代码）
#######################################################################

# 颜色配置（保留你原来的）
cluster_cols <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

# 1. 读取4个样本
DR1 <- Read10X("./DR1/")
DR2 <- Read10X("./DR2/")
DR3 <- Read10X("./DR3/")
DR4 <- Read10X("./DR4/")

# 为每个样本的细胞barcode添加样本标签
colnames(DR1) <- paste(colnames(DR1), "DR1", sep = "_")
colnames(DR2) <- paste(colnames(DR2), "DR2", sep = "_")
colnames(DR3) <- paste(colnames(DR3), "DR3", sep = "_")
colnames(DR4) <- paste(colnames(DR4), "DR4", sep = "_")

# 合并所有样本
experiment.data <- cbind(DR1, DR2, DR3, DR4)
sam.name <- "multi"

if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

# 创建Seurat对象
experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "-"
)

# 修正样本名称（去掉前缀）
experiment.aggregate$orig.ident <- gsub("^1_", "", experiment.aggregate$orig.ident)
cell_names <- colnames(experiment.aggregate)
split_parts <- strsplit(cell_names, "-")
updated_names <- sapply(split_parts, function(x) {
  paste0(x[1], "-", gsub("^1_", "", x[2]))
})
colnames(experiment.aggregate) <- updated_names

Idents(experiment.aggregate) <- experiment.aggregate$orig.ident

# 2. 计算线粒体基因比例
experiment.aggregate[["percent.MT"]] <- PercentageFeatureSet(experiment.aggregate, 
                                                             pattern = "^MT")

# QC可视化（可省略，保留以防需要）
pdf(paste0("./", sam.name, "/QC-VlnPlot.pdf"), width = 15, height = 5)
VlnPlot(experiment.aggregate, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), 
        ncol = 3)
dev.off()

# 3. 质控过滤
cat("Before filter :", nrow(experiment.aggregate@meta.data), "cells\n")
experiment.aggregate <- subset(experiment.aggregate, 
                               subset = 
                                 nFeature_RNA > 50 & 
                                 nFeature_RNA < 6000 & 
                                 nCount_RNA > 10 & 
                                 nCount_RNA < 20000 &
                                 percent.MT < 10)
cat("After filter :", nrow(experiment.aggregate@meta.data), "cells\n")

# 4. 标准化
experiment.aggregate <- NormalizeData(experiment.aggregate, 
                                      normalization.method = "LogNormalize",
                                      scale.factor = 10000)

# 5. 找高变基因
experiment.aggregate <- FindVariableFeatures(experiment.aggregate, 
                                             selection.method = "vst",
                                             nfeatures = 2000)

# 6. 缩放数据（只缩放高变基因，节省内存）
experiment.aggregate <- ScaleData(experiment.aggregate, 
                                  features = VariableFeatures(experiment.aggregate))

# 7. PCA降维
experiment.aggregate <- RunPCA(object = experiment.aggregate, 
                               features = VariableFeatures(experiment.aggregate),
                               verbose = F, npcs = 50)

# 8. Harmony去除批次效应（按 orig.ident）
experiment.aggregate <- RunHarmony(experiment.aggregate, 
                                   reduction = "pca",
                                   group.by.vars = "orig.ident", 
                                   reduction.save = "harmony")

# 9. tSNE可视化
dim.use <- 1:30
experiment.aggregate <- RunTSNE(experiment.aggregate, 
                                reduction = "harmony", 
                                dims = dim.use, 
                                reduction.name = "tsne")
experiment.aggregate <- RunUMAP(experiment.aggregate, reduction = "harmony", dims = 1:30)
# 10. 聚类
experiment.aggregate <- FindNeighbors(experiment.aggregate, 
                                      reduction = "harmony", 
                                      dims = 1:30)
experiment.aggregate <- FindClusters(experiment.aggregate, resolution = 0.5)

# 11. 寻找各cluster的标记基因
library(dplyr)
all.markers <- FindAllMarkers(experiment.aggregate, only.pos = TRUE, 
                              min.pct = 0.3, logfc.threshold = 0.25)

top10.markers <- all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  ungroup()

write.csv(top10.markers,
          file = paste0("./", sam.name, "/Supplementary_Table_Top10_Markers.csv"),
          row.names = FALSE)

write.table(all.markers,
            file = paste0("./", sam.name, "/total_marker_genes.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 12. 细胞类型注释（）
celltype <- data.frame(ClusterID = 0:16, 
                       celltype = 'unknown')
celltype[celltype$ClusterID %in% c(0, 7, 13), 2] <- 'Macrophage'
celltype[celltype$ClusterID %in% c(1, 10), 2] <- 'Mural cell'     # 周细胞
celltype[celltype$ClusterID %in% c(2, 6), 2] <- 'Endothelial cell'
celltype[celltype$ClusterID %in% c(3, 9), 2] <- 'T cell'
celltype[celltype$ClusterID %in% c(4, 5, 14), 2] <- 'Fibroblast'
celltype[celltype$ClusterID %in% c(8), 2] <- 'Monocyte'
celltype[celltype$ClusterID %in% c(11), 2] <- 'cDC'
celltype[celltype$ClusterID %in% c(12, 15), 2] <- 'B cell'

# 将注释添加到Seurat对象
experiment.aggregate@meta.data$celltype <- "NA"
for(i in 1:nrow(celltype)){
  experiment.aggregate@meta.data[which(experiment.aggregate@meta.data$seurat_clusters == celltype$ClusterID[i]),
                                 'celltype'] <- celltype$celltype[i]
}

# 重命名Ident为celltype（方便后续画图）
new.cluster.ids <- c("Macrophage", "Mural cell", "Endothelial cell",
                     "T cell", "Fibroblast", "Fibroblast",
                     "Endothelial cell", "Macrophage", "Monocyte",
                     "T cell", "Mural cell", "cDC",
                     "B cell", "Macrophage", "Fibroblast",
                     "B cell","unknown")
names(new.cluster.ids) <- levels(experiment.aggregate)
experiment.aggregate <- RenameIdents(experiment.aggregate, new.cluster.ids)

# 将最终的注释结果存入 scobj
scobj <- experiment.aggregate
all.markers <- FindAllMarkers(scobj, only.pos = TRUE, 
                              min.pct = 0.3, logfc.threshold = 0.25)

top_markers <- all.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup()

write.csv(top_markers, "top_markers_per_celltype.csv", row.names = FALSE)

# ==================== 周细胞亚群精准鉴定 ====================
# 步骤：从 scobj 中提取 Mural cell，再聚类，用标志物鉴定周细胞 vs 平滑肌细胞

# 1. 提取所有 Mural cells
mural_cells <- subset(scobj, subset = celltype == "Mural cell")

# 2. 对 Mural cells 重新处理（标准化、降维、聚类）
mural_cells <- NormalizeData(mural_cells) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony(group.by.vars = "orig.ident", reduction = "pca", reduction.save = "harmony") %>%
  RunTSNE(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.5)

# 3. 定义周细胞阳性标志物与平滑肌阴性标志物（人源）
pericyte_pos  <- c("PDGFRB", "CSPG4", "RGS5", "DLK1", "ANPEP")
smooth_muscle <- c("MYH11", "CNN1")
endothelial   <- c("PECAM1", "CDH5")   # 排除内皮污染

# 4. 可视化标志物以辅助判断（PDF 输出）
pdf("Mural_subcluster_markers.pdf", width = 14, height = 10)
FeaturePlot(mural_cells, features = c(pericyte_pos, smooth_muscle, endothelial),
            reduction = "tsne", ncol = 4)
dev.off()

# 检查 mural_cells 里有几个 cluster
table(Idents(mural_cells))

# 画出分群图和标志物分布，目测找出哪几个数字对应 MYH11 阴性/阳性
DimPlot(mural_cells, reduction = "tsne", label = T)
DotPlot(mural_cells, features = c("PDGFRB", "CSPG4", "RGS5", "MYH11", "CNN1")) + RotatedAxis()
# 5. 【手动修改】根据 DotPlot 结果确认 Cluster 编号
pericyte_clusters <- c(0, 1,4)  # 周细胞（PDGFRB+, CSPG4+, RGS5+, MYH11-）
vsmc_clusters     <- c(3)     # 平滑肌细胞（MYH11+, CNN1+, CSPG4低）

# 6. 获取对应的 Barcode
pericyte_barcodes <- colnames(subset(mural_cells, idents = pericyte_clusters))
vsmc_barcodes     <- colnames(subset(mural_cells, idents = vsmc_clusters))

# 7. 在 scobj 中创建精细注释列 celltype_fine
scobj$celltype_fine <- as.character(scobj$celltype) 
scobj$celltype_fine[pericyte_barcodes] <- "Pericyte"
scobj$celltype_fine[vsmc_barcodes]     <- "VSMC"

# 【可选】如果你想把 C2 单独标为 "Transitional_Mural"，取消下一行注释
 scobj$celltype_fine[colnames(subset(mural_cells, idents = c(2)))] <- "Transitional"

# 【关键】对于 C1（未知细胞），如果不归入 VSMC，直接写个逻辑让它保持原样或者改成 "Mural_Other"
# 最稳妥的做法：不要把它覆盖成 VSMC。
mural_other_barcodes <- colnames(subset(mural_cells, idents = c(1))) 
# 可以将其命名为 "Mural_unclassified" 或者检查原始注释，如果原始就是 Mural，保持即可。
 scobj$celltype_fine[mural_other_barcodes] <- "Mural_Unclassified" 

# 更新 Idents
Idents(scobj) <- factor(scobj$celltype_fine)

# 8. 查看分布并保存精细注释图
table(scobj$celltype_fine)
# 验证一下你标记的细胞是否正确：
DotPlot(subset(scobj, subset = celltype_fine %in% c("Pericyte", "VSMC")), 
        features = c("PDGFRB", "CSPG4", "RGS5", "MYH11", "CNN1")) + RotatedAxis()
####分离杂质
# 1. 创建新对象（推荐），保留高质量细胞
# 将 Mural_Unclassified, unknown, 及原有的 cDC, unknown 等都可能包含低质量杂质
scobj_clean <- subset(scobj, subset = celltype_fine %in% c("Pericyte", "VSMC", "Transitional", 
                                                           "B cell", "cDC", "Endothelial cell", 
                                                           "Fibroblast", "Macrophage", "Monocyte", "T cell"))

# 2. 如果不想创建新对象，只想在画图时过滤（最快的方法）
Idents(scobj) <- "celltype_fine"
DimPlot(subset(scobj, idents = c("Pericyte", "VSMC", "Transitional")), 
        reduction = "tsne", label = TRUE, label.size = 4) + ggtitle("Mural Subtypes")

# 3. 如果想在整体 tSNE 图中隐藏不相关的细胞，建议用“黑名单”排除
scobj_sub <- subset(scobj, subset = celltype_fine != "Mural_Unclassified" & celltype_fine != "unknown")
# 调整因子顺序，把新的精细亚群排在前面
scobj_sub$celltype_fine <- factor(scobj_sub$celltype_fine, 
                                  levels = c("Pericyte", "VSMC", "Transitional", 
                                             "B cell", "cDC", "Endothelial cell", 
                                             "Fibroblast", "Macrophage", "Monocyte", "T cell"))
# 准备工作：确保颜色是正确的（用你前面定义的 cluster_cols）
Idents(scobj_sub) <- "celltype_fine"

# ================= 全局细胞分群图 (左上) =================
p1 <- DimPlot(scobj_sub, reduction = "umap", label = TRUE, label.size = 4, 
              pt.size = 0.8, cols = cluster_cols) + 
  NoLegend() + ggtitle("Global Cell Clusters")
print(p1)
# ================= UROD 特征图 (右上) =================
# 注意：min.cutoff = 'q10' 是用来过滤掉灰色的背景噪声，只突出表达量前10%的细胞
# 使用 min.cutoff 和 max.cutoff 强制拉伸颜色
# 'q10' 忽略最低的10%背景噪声，max.cutoff 设为 'q90' 或具体数值 1.5，防止被极值拖累
p2 <- FeaturePlot(scobj_sub, features = "UROD", reduction = "umap", 
                  min.cutoff = 'q10', max.cutoff = 'q90', 
                  pt.size = 0.8) + 
  ggtitle("UROD Expression ")
print(p2)
# 专门画一张大的 UROD DotPlot 放在你的图里
p3_highlight <- DotPlot(scobj_sub, features = "UROD", group.by = "celltype_fine") + 
  RotatedAxis() +
  scale_color_gradientn(colors = c("lightgrey", "blue")) + # 确保点颜色有对比
  ggtitle("UROD: Positive Cell Proportion") +
  theme(axis.text.x = element_text(size = 10, face = "bold"))
# 打印这张图，你会发现周细胞那个点的面积非常巨大，远大于其他细胞
print(p3_highlight)
library(ggpubr)
# 把分组切好
library(ggpubr) 
# 【关键修复】：加上 scale_y_continuous + oob = scales::squish
# 这个参数的意思是：将 y 轴锁定在 0 到 1.5 之间，如果有超过 1.5 的极值，把它们压缩（Squish）到 1.5 显示，而不是拉长 Y 轴。
p4_fixed <- VlnPlot(scobj_sub, features = "UROD", group.by = "celltype_fine", pt.size = 0, cols = cluster_cols) + 
  scale_y_continuous(limits = c(0, 1.5), oob = scales::squish) + 
  stat_compare_means(comparisons = list(c("Pericyte", "Fibroblast"), 
                                        c("Pericyte", "Macrophage"), 
                                        c("Pericyte", "T cell")),
                     label = "p.signif", # 改用星号（*, **, ***），这样图不会太拥挤
                     method = "wilcox.test") +
  theme(legend.position = "none") + # 去除图例让图更宽
  ggtitle("UROD Expression (Scale optimized)")

print(p4_fixed)

# ================= 5. 将四张图组合并保存为 PDF =================
library(cowplot) # 确保加载了 cowplot
p5<-plot_grid(p1, p2, p3_highlight, p4_fixed, labels = c("A", "B", "C", "D"), ncol = 2)
ggsave("Final_UROD_4Panel_Plots.pdf", p5, width = 16, height = 12)
# 打印提示
print("四宫格图已成功生成并保存为 PDF！")
#########虚拟敲除
install.packages("C:/R_packages/scTenifoldKnk-1.0.3", repos = NULL, type = "source")
library(scTenifoldKnk)
library(Seurat)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggrepel)

# ============================
# 步骤1：提取Mural cells并重新聚类
# ============================
mural_cells <- subset(scobj, subset = celltype == "Mural cell")

mural_cells <- NormalizeData(mural_cells) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony(group.by.vars = "orig.ident", reduction = "pca", reduction.save = "harmony") %>%
  RunTSNE(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.5)

# 查看有几个cluster
table(Idents(mural_cells))

# 可视化：决定哪些cluster是周细胞、平滑肌、过渡细胞
pdf("Mural_subcluster_identification.pdf", width = 14, height = 10)
DimPlot(mural_cells, reduction = "tsne", label = TRUE)
DotPlot(mural_cells, 
        features = c("PDGFRB", "CSPG4", "RGS5", "DLK1", "ANPEP",   # 周细胞阳性标志
                     "MYH11", "CNN1",                                # 平滑肌标志
                     "PECAM1", "CDH5"))                              # 内皮排除标志
dev.off()

# ============================
# 步骤2：⚠️ 请你根据上面输出的PDF和DotPlot，手动修改这里的簇号
# ============================
# 1. 先用你之前分的簇打个粗标签
mural_cells$celltype_raw <- ifelse(
  mural_cells$seurat_clusters %in% c(0, 1, 4), "Pericyte",
  ifelse(mural_cells$seurat_clusters %in% c(2, 3), "VSMC", "others")
)

# 2. 周细胞提纯：只要 PDGFRB > 1.5 并且 MYH11 == 0 并且 CNN1 == 0
pure_pericyte_cells <- colnames(mural_cells)[
  mural_cells$celltype_raw == "Pericyte" &
    FetchData(mural_cells, vars = "PDGFRB")[,1] > 1.5 &
    FetchData(mural_cells, vars = "MYH11")[,1] == 0 &
    FetchData(mural_cells, vars = "CNN1")[,1] == 0
]

# 3. 平滑肌保持原标签（或也过滤一下）
pure_vsmc_cells <- colnames(mural_cells)[mural_cells$celltype_raw == "VSMC"]

# 4. 覆盖精细标签
mural_cells$celltype_fine <- "others"
mural_cells$celltype_fine[pure_pericyte_cells] <- "Pericyte"
mural_cells$celltype_fine[pure_vsmc_cells] <- "VSMC"

table(mural_cells$celltype_fine)

# ============================
# 步骤4：把精细注释回写到原始scobj（保持一致）
# ============================
# 如果scobj里已经有celltype_fine列，建议先删除再写，避免冲突
if("celltype_fine" %in% colnames(scobj@meta.data)) {
  scobj$celltype_fine <- NULL
}
scobj$celltype_fine <- "others"
scobj$celltype_fine[colnames(mural_cells)] <- mural_cells$celltype_fine
table(scobj$celltype_fine)

# ============================
# 步骤5：确保scobj有UMAP（后面画图需要）
# ============================
if(!"umap" %in% names(scobj@reductions)) {
  message("原始scobj中没有UMAP，正在计算...")
  if(!"harmony" %in% names(scobj@reductions)) {
    # 如果连harmony都没有，先跑一次harmony（通常scobj已经有了）
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

# ============================
# 步骤6：定义颜色（供画图使用）
# ============================
cluster_cols <- c("Pericyte" = "#E41A1C", 
                  "VSMC"     = "#377EB8", 
                  "Transitional" = "#4DAF4A",
                  "others"   = "grey80")
# 如果需要更多颜色，可以自行用RColorBrewer扩展

# ============================
# 步骤7：提取纯净周细胞（去除残留的MYH11/CNN1表达细胞）
# ============================
# 步骤7 修改后
pbmc_pure <- subset(scobj, 
                    subset = celltype_fine == "Pericyte" & 
                      MYH11 == 0 & 
                      CNN1 == 0)   
# 注释可写成：已进行log归一化，表达值为0表示无UMI，直接剔除 
# 查看最终细胞数
n_pericyte_clean <- ncol(pbmc_pure)
message("最终纯净周细胞数量: ", n_pericyte_clean)

# ============================
# 步骤8：对比去污染前后的关键基因表达（四宫格图）
# ============================
pdf("Pericyte_purification_check.pdf", width = 12, height = 10)

# 8.1 UMAP上标记纯净周细胞的位置
DimPlot(scobj, reduction = "umap", group.by = "celltype_fine", cols = cluster_cols, order = TRUE) + 
  ggtitle("所有细胞的精细亚型UMAP")

DimPlot(pbmc_pure, reduction = "umap", cols = cluster_cols["Pericyte"]) + 
  ggtitle("纯净周细胞（去除MYH11/CNN1表达细胞）")

# 8.2 小提琴图：MYH11 和 CNN1 的去除效果
VlnPlot(scobj, features = c("MYH11", "CNN1"), 
        group.by = "celltype_fine", pt.size = 0, cols = cluster_cols) + 
  plot_annotation(title = "原始精细分组中的平滑肌标志表达")

VlnPlot(pbmc_pure, features = c("MYH11", "CNN1"), 
        group.by = "celltype_fine", pt.size = 0, cols = cluster_cols["Pericyte"]) + 
  plot_annotation(title = "纯净周细胞中的平滑肌标志（应全部为0）")

# 8.3 小提琴图：周细胞阳性标志物确认
VlnPlot(scobj, features = c("PDGFRB", "RGS5", "ACTA2", "TAGLN"), 
        group.by = "celltype_fine", pt.size = 0, cols = cluster_cols) + 
  plot_annotation(title = "周细胞/活化标志在各组中的表达")

VlnPlot(pbmc_pure, features = c("PDGFRB", "RGS5", "ACTA2", "TAGLN"), 
        group.by = "celltype_fine", pt.size = 0, cols = cluster_cols["Pericyte"]) + 
  plot_annotation(title = "纯净周细胞的标志物表达")

# 8.4 最后确认一下DES表达（你关心的）
VlnPlot(scobj, features = "DES", group.by = "celltype_fine", pt.size = 0, cols = cluster_cols) + 
  ggtitle("DES在各组中的表达")
VlnPlot(pbmc_pure, features = "DES", pt.size = 0) + ggtitle("纯净周细胞DES表达")

dev.off()

# ============================
# 完成！纯净周细胞保存在 pbmc_pure 对象中
# ============================
message("分析完成。最终周细胞对象: pbmc_pure，含 ", ncol(pbmc_pure), "个细胞")
#虚拟敲除
# ==================== 1. 提取纯净周细胞 ====================
pbmc_pure <- subset(scobj_sub, 
                    subset = celltype_fine == "Pericyte" & 
                      MYH11 == 0 & 
                      CNN1 == 0)

# ==================== 2. 准备输入矩阵 ====================
# 提取原始 counts（注意：scTenifoldKnk 需要原始 counts，不要 log 归一化）
countMat <- GetAssayData(pbmc_pure, assay = "RNA", layer = "counts")

# 选高变基因以降低计算量（加 UROD 防止它被漏掉）
pbmc_pure <- FindVariableFeatures(pbmc_pure, 
                                  assay = "RNA",        # 指定 assay
                                  selection.method = "vst", 
                                  nfeatures = 10000)    # 你可以改成 2000 节约时间
hvgs <- VariableFeatures(pbmc_pure)
target_genes <- unique(c("UROD", hvgs))

# 构建最终用的 count matrix（行是基因，列是细胞）
count_matrix <- as.matrix(countMat[target_genes, ])
n_pericyte_clean <- ncol(pbmc_pure)
message("最终纯净周细胞数量: ", n_pericyte_clean)
#再次检查污染
VlnPlot(pbmc_pure, features = c("PDGFRB", "RGS5", "MYH11", "CNN1", "ACTA2", "TAGLN"), 
        pt.size = 0, ncol = 3)


# 检查 MYH11 是否真的全为 0
table(FetchData(pbmc_pure, vars = "MYH11")[,1] > 0)
# 如果输出是 FALSE n（n是你的周细胞数），就说明 0 表达

# ==================== 3. 运行虚拟敲除 ====================
result_knockout <- scTenifoldKnk(
  countMatrix = count_matrix,      # 你的纯净 count 矩阵
  gKO = "UROD",
  qc = TRUE,
  qc_mtThreshold = 0.15,           # 适当放宽，避免误杀
  qc_minLSize = 500,               # 改为 500，保留更多细胞
  nc_nNet = 3,                     # 网络数不变（已是最低）
  nc_nCells = min(250, ncol(count_matrix))  # 如果细胞数>250，就抽250个，否则全用
)
# ==================== 4. 查看结果 ====================

# 1. 检查整个结果对象里到底有什么
names(result_knockout)

# 2. 看 UROD 基因的表达量分布
urod_exp <- count_matrix["UROD", ]
summary(urod_exp)
# 计算表达细胞比例
sum(urod_exp > 0) / length(urod_exp)
# 看这个对象的类型和结构
class(result_knockout$diffRegulation)
dim(result_knockout$diffRegulation)
head(result_knockout$diffRegulation)
sig_genes <- result_knockout$diffRegulation[result_knockout$diffRegulation$Z > 2, ]
head(sig_genes, 20)
sig_genes_adj <- subset(result_knockout$diffRegulation, p.adj < 0.05)
nrow(sig_genes_adj)
# 5. 提取结果并查看网络扰动最剧烈的基因
df_pure = result_knockout$diffRegulation
df_pure <- df_pure[df_pure$gene != "UROD", ]
# 按 p.adj 从小到大排序
df_pure_sorted <- df_pure[order(df_pure$p.adj), ]

# 打印前20个受敲除影响最大的基因
head(df_pure_sorted, 20)
library(ggplot2)
library(ggrepel)

df <- result_knockout$diffRegulation
# 去除 UROD 本身
df_noUROD <- df[df$gene != "UROD", ]

# 方法1：Z > 1.5 且按 Z 降序，取前 20 个
candidate <- head(df_noUROD[order(-df_noUROD$Z), ], 20)

#####
df$log10p <- -log10(df$p.value + 1e-300)
# 标记这些候选基因（包括血红蛋白和非血红蛋白）
top_genes <- candidate$gene
df$label_me <- df$gene %in% top_genes
pdf("Differential regulation upon UROD knockout.pdf", width = 12, height = 10)
p_differential<-ggplot(df, aes(x = Z, y = log10p)) +
  geom_point(aes(color = label_me), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_text_repel(data = subset(df, label_me), 
                  aes(label = gene), 
                  size = 3.5, 
                  max.overlaps = 50) +
  labs(x = "Z-score", y = "-log10(p-value)") +
  theme_minimal() +
  ggtitle("Differential regulation upon UROD knockout") +
  theme(legend.position = "none")
print(p_differential)
dev.off()
top15 <- head(df_noUROD[order(-df_noUROD$Z), ], 15)
# 可保留血红蛋白在内
pdf("Top 15 perturbed genes.pdf", width = 12, height = 10)
p_top <- ggplot(top15, aes(x = reorder(gene, Z), y = Z)) +
  geom_col(aes(fill = Z), show.legend = FALSE) +
  scale_fill_gradient2(low = "#D6EAF8", mid = "#85C1E9", high = "#2E86C1", midpoint = median(top15$Z)) +  
  coord_flip() +
  labs(x = "", y = "Z-score", 
       title = "Top 15 perturbed genes (excluding UROD)") +
  theme_minimal()
print(p_top)
dev.off()
####富集分析
# 读取差异表（假设变量名为 df，你自己按实际名称替换）
df <- result_knockout$diffRegulation
genes_ora <- df$gene[df$gene != "UROD"  & abs(df$Z) > 1.5]
# 创建排序向量（用 Z 值或符号化的 FC）
geneList <- df$Z
names(geneList) <- df$gene
geneList <- sort(geneList, decreasing = TRUE)   # 高 Z 的排在前面
#GO
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tidytext)   # 用于 reorder_within

ego_bp <- enrichGO(gene         = genes_ora,
                   OrgDb        = org.Hs.eg.db,
                   keyType      = "SYMBOL",
                   ont          = "BP",
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.2)

ego_mf <- enrichGO(gene         = genes_ora,
                   OrgDb        = org.Hs.eg.db,
                   keyType      = "SYMBOL",
                   ont          = "MF",
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.2)

ego_cc <- enrichGO(gene         = genes_ora,
                   OrgDb        = org.Hs.eg.db,
                   keyType      = "SYMBOL",
                   ont          = "CC",
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.2)

bp_df <- as.data.frame(ego_bp)
mf_df <- as.data.frame(ego_mf)
cc_df <- as.data.frame(ego_cc)
bp_df$ONTOLOGY <- "BP"
mf_df$ONTOLOGY <- "MF"
cc_df$ONTOLOGY <- "CC"
all_go <- rbind(bp_df, mf_df, cc_df)
library(dplyr)
top_go <- all_go %>%
  group_by(ONTOLOGY) %>%
  top_n(-5, wt = p.adjust) %>%   # 按调整后 p 值最小的取
  ungroup()

library(ggplot2)

library(ggplot2)

p_go <- ggplot(top_go, aes(x = -log10(p.adjust), 
                           y = reorder_within(Description, p.adjust, ONTOLOGY),
                           fill = p.adjust)) +  # 👉 关键点：将颜色映射改为原始的 p.adjust
  
  geom_col(width = 0.7) +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
  scale_y_reordered() +
  
  # 👉 设置红-紫-蓝三色渐变（图例自动显示为 p.adjust 的数值）
  scale_fill_gradient2(
    mid = "#D96B7B",        # 小 p 值（极显著）：三文鱼色/浅红（对应高柱子）
    high = "#2166AC",       # 大 p 值（不显著）：蓝色（对应低柱子）
    midpoint = median(top_go$p.adjust), # 以原始 p.adjust 的中位数作为紫色的分界点
    name = "p.adjust"       # 强制让图例名字叫 p.adjust
  ) +
  
  labs(x = "-log10(adjusted p-value)", y = NULL, title = "GO Enrichment (BP, MF, CC)") +
  
  # 👉 保持你喜欢的 Nature 风格干净主题
  theme_classic() +   
  theme(
    strip.background = element_blank(),           
    strip.text = element_text(face = "bold", size = 12, color = "black"), 
    axis.text = element_text(color = "black", size = 10),    
    axis.line = element_line(color = "black"),               
    legend.position = "right"
  )
ggsave("GO_ORA2.pdf", p_go,width = 8, height = 10)
write.csv(all_go, "GO_ORA_results.csv")
#GSEA
library(org.Hs.eg.db)
ego_gsea <- gseGO(geneList     = geneList,
                  OrgDb        = org.Hs.eg.db,
                  keyType      = "SYMBOL",
                  ont          = "BP",
                  pvalueCutoff = 0.2,
                  verbose      = FALSE)

head(ego_gsea)
write.csv(ego_gsea, "GO_BP_GSEA_results.csv")

# 可视化（山脊图或气泡图）
ridgeplot(ego_gsea, showCategory = 10) + ggtitle("GO BP (GSEA)")
ggsave("ego_gsea1.pdf", width = 8, height = 10)
# 转换函数
gene.df <- bitr(genes_ora, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)



###GSEA AND KEGG
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)

# 1. 排序向量（Symbol 版）
gene_list_sym <- df$Z
names(gene_list_sym) <- df$gene
gene_list_sym <- sort(gene_list_sym, decreasing = TRUE)

# 2. 用你已有的 pathway_list（就是之前 split 生成的列表）
#    它已经是 list，每条通路对应 Symbol 基因向量
gsea_result <- fgsea(
  pathways = pathway_list,
  stats    = gene_list_sym,
  minSize  = 10,
  maxSize  = 500,
  nperm    = 10000
)

# 3. 分开 KEGG 和 WikiPathways
kegg_fgsea <- gsea_result[grepl("^KEGG_", pathway), ]
wiki_fgsea <- gsea_result[grepl("^WP_", pathway), ]

# 4. 查看显著结果
head(kegg_fgsea[order(pval), ], 20)
head(wiki_fgsea[order(pval), ], 20)
####
library(ggplot2)
library(forcats)   # 用于排序

# 取显著前 20 条（可根据需要调整）
top <- gsea_result %>%
  filter(padj < 0.25) %>%
  arrange(padj) %>%
  head(20)

P5<-ggplot(top, aes(x = NES, 
                y = fct_reorder(pathway, NES),
                size = -log10(padj),
                color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Normalized Enrichment Score", y = "",
       size = expression(-log[10](padj)),
       color = "NES") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank()) +
  ggtitle("GSEA (KEGG & WikiPathways)")
ggsave("ego_kegg.pdf", P5,width = 8, height = 10)
#######################################################################
# 第三部分：细胞通讯分析（CellChat）
#######################################################################
source_cell <- "Pericyte"            # 周细胞
target_t <- "T cell"                # T 细胞
target_ec <- "Endothelial cell"           # 视网膜内皮细胞
# ==============================================================

library(CellChat)
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
cat("\n===== 开始细胞通讯分析 =====\n")

# ---------- 直接复用前面已纯化好的周细胞，不重新聚类 ----------
# pbmc_pure 里的细胞就是经过 PDGFRB>1.5 & MYH11==0 & CNN1==0 纯化的周细胞
pericyte_cells <- colnames(pbmc_pure)

# 在完整 scobj 上更新 celltype：把这批细胞标为 Pericyte
scobj$celltype <- as.character(scobj$celltype)
scobj$celltype[pericyte_cells] <- "Pericyte"

# 原来的 Mural cell 里，除了纯化出来的 Pericyte，其余标记为 mural_other
mural_all <- WhichCells(scobj, expression = celltype == "Mural cell")
non_pericyte_mural <- setdiff(mural_all, pericyte_cells)
scobj$celltype[non_pericyte_mural] <- "mural_other"

table(scobj$celltype)

# ---------- 1. 创建 CellChat 对象 ----------
cellchat <- createCellChat(object = scobj,
                           meta = scobj@meta.data,
                           group.by = "celltype")
cellchat <- setIdent(cellchat, ident.use = "celltype")
cat("检测到的细胞类型：\n")
print(levels(cellchat@idents))
groupSize <- as.numeric(table(cellchat@idents))

# ---------- 2. 加载数据库 ----------
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

# ---------- 3. 预处理与通讯计算 ----------
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = "cellchat_object.rds")

# ============================================================
# 可视化（重点：周细胞 → T 细胞，周细胞 → 视网膜内皮细胞）
# ============================================================

# ---------- 全局互作 ----------
pdf("Global_Interaction_Overview.pdf", width = 12, height = 5)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F,
                 title.name = "Interaction weights/strength")
dev.off()

# ---------- 发送者圆图 ----------
pdf("Each_Cell_as_Sender_Circle.pdf", width = 12, height = 10)
mat <- cellchat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T,
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# ---------- 气泡图：周细胞 → T 细胞 ----------
cat("\n===== 周细胞 -> T 细胞 气泡图 =====\n")
unique(cellchat@idents)
# 或者
levels(cellchat@idents)
p_bubble_T <- netVisual_bubble(cellchat,
                               sources.use = "Pericyte", 
                               targets.use = "T cell",
                               remove.isolate = FALSE,
                               title.name = "Pericyte -> T cell")
ggsave("Pericyte_to_T_bubble.pdf", p_bubble_T, width = 8, height = 6)
cat("✅ 气泡图已保存: Pericyte_to_T_bubble.pdf\n")

# ---------- 气泡图：周细胞 → 视网膜内皮细胞 ----------
p_bubble_T2 <- netVisual_bubble(cellchat,
                               sources.use = "Pericyte", 
                               targets.use = "Endothelial cell",
                               remove.isolate = FALSE,
                               title.name = "Pericyte -> Endothelial cell")
ggsave("Pericyte_to_T2_bubble.pdf", p_bubble_T2, width = 8, height = 6)
cat("✅ 气泡图已保存: Pericyte_to_T_bubble.pdf\n")

# ─── 周细胞 → T 细胞：只关注 MIF, MK ───
pathways_T <- c("MIF", "MK")

# 检查必须的包
library(CellChat)
library(ggplot2)

for (pathway in pathways_T) {
  if (!pathway %in% cellchat@netP$pathways) {
    cat("⚠️ 通路", pathway, "不在数据中，跳过\n")
    next
  }
  cat("\n===== 处理通路:", pathway, " (周细胞 -> T 细胞) =====\n")
  
  # --- 1. 绘制和弦图 ---
  pdf(paste0("Pathway_", pathway, "_Pericyte_to_T_Chord.pdf"), width = 8, height = 7)
  netVisual_chord_cell(cellchat, signaling = pathway,
                       sources.use = "Pericyte",
                       targets.use = "T cell",
                       lab.cex = 0.8,
                       # 【修改点1】把变量 Pericyte 改成了引号括起来的字符串 "Pericyte"
                       title.name = paste0(pathway, " (Pericyte -> T cell)")) 
  dev.off()
}
graphics.off()
# ─── 周细胞 → 视网膜内皮细胞：VEGF, CypA, MK ───
# ─── 步骤 0：务必先定义细胞群变量 ───
source_cell <- "Pericyte"
target_ec <- "Endothelial cell"

# 简单确认是否定义成功
cat("✅ 来源细胞 (source_cell) =", source_cell, "\n")
cat("✅ 目标细胞 (target_ec) =", target_ec, "\n")

pathways_EC <- c("VEGF", "CypA", "MK")
cat("🔄 即将处理通路:", paste(pathways_EC, collapse=", "), "\n")


# ─── 步骤 1：循环绘制和弦图 ───
for (pathway in pathways_EC) {
  # 检测通路是否存在，防止报错中止
  if (!pathway %in% cellchat@netP$pathways) {
    cat("⚠️ 通路", pathway, "不在当前 cellchat 对象中，跳过\n")
    next
  }
  cat("\n===== 处理通路:", pathway, " (周细胞 -> 内皮细胞) =====\n")
  
  # 保存和弦图
  pdf(paste0("Pathway_", pathway, "_Pericyte_to_EC_Chord.pdf"), width = 8, height = 7)
  netVisual_chord_cell(cellchat, signaling = pathway,
                       sources.use = source_cell,
                       targets.use = target_ec,
                       lab.cex = 0.8,
                       title.name = paste0(pathway, " (", source_cell, " -> ", target_ec, ")"))
  dev.off()
  cat("✅ 和弦图已保存\n")
}

# ============================================================
# 基因表达验证（气泡图、小提琴图、特征图）
# ============================================================
cat("\n===== 开始整合通路基因表达验证 =====\n")

# ---------- 定义所有需验证的基因 ----------
# 四条通路的核心配体与受体（根据 CellChat 数据库和常见文献）
genes_all <- c(
  "MIF",          # MIF 通路配体
  "CD74", "CXCR4", "CD44",        # MIF 受体
  "MDK",          # MK 通路配体 (MK = MDK)
  "NCL", "ITGA4", "ITGB1",        # MK 受体
  "VEGFA", "VEGFB", "PGF",        # VEGF 通路配体（可只选 VEGFA 作为代表，但多写无妨）
  "KDR", "FLT1", "NRP1", "NRP2",  # VEGF 受体
  "PPIA",                         # CypA 通路配体 (CypA = PPIA)
  "BSG", "PTPRC"                  # CypA 受体 (CD147=BSG, CD45=PTPRC)
)

# 如果某些基因你的 scobj 里没有，R 会跳过并警告，不影响其他图。

# ---------- 气泡图（全细胞类型） ----------
p_dot <- DotPlot(scobj, features = genes_all, group.by = "celltype") +
  scale_color_gradientn(colours = c("lightgrey", "#1E90FF", "#D22B2B")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11, color = "black")) +
  labs(x = "Ligand/Receptor Genes", y = "Cell Types",
       title = "Expression of all key pathway genes (MIF, MK, VEGF, CypA)")
ggsave("DotPlot_AllPathways_Validation.pdf", plot = p_dot, width = 16, height = 7, dpi = 300)

# ---------- 小提琴图 ----------
# 为了可读性，将基因分成两组：配体一组，受体一组
ligands_all <- c("MIF", "MDK", "VEGFA", "VEGFB", "PGF", "PPIA")
receptors_all <- c("CD74", "CXCR4", "CD44", "NCL", "ITGA4", "ITGB1",
                   "KDR", "FLT1", "NRP1", "NRP2", "BSG", "PTPRC")

p_lig <- VlnPlot(scobj, features = ligands_all, group.by = "celltype", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Ligands expression")
p_rec <- VlnPlot(scobj, features = receptors_all, group.by = "celltype", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Receptors expression")
ggsave("VlnPlot_Ligands_All.pdf", plot = p_lig, width = 14, height = 6)
ggsave("VlnPlot_Receptors_All.pdf", plot = p_rec, width = 18, height = 6)

# ---------- 特征图（t-SNE） ----------
# 配体特征图
p_feat_lig <- FeaturePlot(scobj, features = ligands_all, reduction = "tsne",
                          cols = c("lightgrey", "red"), pt.size = 0.5,
                          order = TRUE, label = TRUE, label.size = 4, ncol = 3) +
  theme(plot.title = element_text(size = 12, face = "bold"))
ggsave("FeaturePlot_Ligands_All_tsne.pdf", plot = p_feat_lig, width = 14, height = 8, dpi = 300)

# 受体特征图（挑选最重要的几个，避免太挤）
receptors_focus <- c("CXCR4", "CD44", "ITGA4", "ITGB1",        # MK/MIF
                     "KDR", "FLT1", "NRP1",                    # VEGF
                     "BSG")                                      # CypA (CD147)
p_feat_rec <- FeaturePlot(scobj, features = receptors_focus, reduction = "tsne",
                          cols = c("lightgrey", "blue"), pt.size = 0.5,
                          order = TRUE, label = TRUE, label.size = 4, ncol = 3) +
  theme(plot.title = element_text(size = 12, face = "bold"))
ggsave("FeaturePlot_Receptors_All_tsne.pdf", plot = p_feat_rec, width = 14, height = 10, dpi = 300)

# ---------- 精简气泡图：仅展示来源和目标细胞 ----------
# 1. 确保变量存在（防止之前没跑通导致 source_cell 等为空）
if(!exists("source_cell")) source_cell <- "Pericyte"
if(!exists("target_t")) target_t <- "T cell"
if(!exists("target_ec")) target_ec <- "Endothelial cell"

cells_of_interest <- c(source_cell, target_t, target_ec)

# ⚠️ 重点提醒：请检查您的 Seurat 对象 scobj 的 metadata 中，细胞注释的列名是否真的叫 "celltype"？
# 如果您的元数据列名叫 "cell_type"，"cluster"，或者 "cell_annotation"，请把下面的 "celltype" 替换成您的真实列名。
# 如果您的注释是直接作为对象的活动 ident，可以写成：scobj_sub <- subset(scobj, idents = cells_of_interest)

scobj_sub <- subset(scobj, subset = celltype %in% cells_of_interest)

# 2. 整理目标基因列表
genes_focus <- c("MIF", "MDK", "VEGFA", "PPIA",      # 配体
                 "CD74", "CXCR4", "CD44", "NCL",      # MIF/MK 受体
                 "KDR", "FLT1", "BSG")               # VEGF/CypA 核心受体

# 3. 绘制气泡图
p_dot_sub <- DotPlot(scobj_sub, features = genes_focus, group.by = "celltype") +
  scale_color_gradientn(colours = c("lightgrey", "#1E90FF", "#D22B2B")) +
  RotatedAxis() + theme_classic() +
  labs(x = "Ligand/Receptor Genes", y = "Cell Types", 
       title = "Key pathway genes in source & target cells")

# 4. 保存图片
ggsave("DotPlot_SourceTarget_All.pdf", plot = p_dot_sub, width = 12, height = 4, dpi = 300)

cat("\n✅ 细胞互作来源与目标细胞的关键基因气泡图已保存！\n")
###=================
# 1. 重新计算中位数，给 pbmc_pure 打上 UROD 高/低 的分组标签
urod_exp <- FetchData(pbmc_pure, vars = "UROD")
threshold <- median(urod_exp$UROD)

# 添加 UROD_status 列
pbmc_pure$UROD_status <- ifelse(urod_exp$UROD > threshold, "UROD_high "," UROD_low")

# 必须把分组设为当前对象的默认 Idents，否则 VlnPlot 依然找不到
Idents(pbmc_pure) <- "UROD_status"

# 2. 再次运行你刚才的小提琴图代码（现在已经彻底通了！）
library(ggpubr)
VlnPlot(pbmc_pure, features = c("PGF", "PPIA", "MDK"), group.by = "UROD_status", pt.size = 0) + 
  stat_compare_means(comparisons = list(c("UROD_WT", "UROD_KO")), method = "wilcox.test") +
  scale_y_continuous(limits = c(0, 2.5), oob = scales::squish) +
  ggtitle("Endothelial Attacking Ligands in UROD KO Pericytes")
# 提取你 scobj_sub 里的内皮细胞
endo_sub <- subset(scobj_sub, idents = "Endothelial cell")

# 把内皮细胞根据你的分组（DR 组 vs Control 组）画图
# 看内皮细胞上的受体在 DR 组是不是表达更高！
VlnPlot(endo_sub, features = c("FLT1", "BSG", "ITGA4"), pt.size = 0) 