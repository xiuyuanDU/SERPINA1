###1.Unclassified cells
setwd('F:/scRNA/workpath/0.allcell')
thyroid.combined <- readRDS("thyroid.combined_after_annotation.rds")
setwd('F:/scRNA/workpath/1.second_annotation')

Subset_Unclassified <- subset(thyroid.combined, idents = "Unclassified cells")
Subset_Unclassified <- ScaleData(Subset_Unclassified)

# Find differentially expressed genes
Subset_Unclassified <- FindVariableFeatures(Subset_Unclassified, selection.method = "vst", nfeatures = 2000)

# Dimensionality reduction
Subset_Unclassified <- RunPCA(Subset_Unclassified, features = VariableFeatures(object = Subset_Unclassified))
Subset_Unclassified <- RunTSNE(Subset_Unclassified, dims = 1:20)
Subset_Unclassified <- RunHarmony(Subset_Unclassified, group.by.vars = "stim")

# Elbow plot
pdf(file = "5.ElbowPlot.pdf", width = 5, height = 4)
ElbowPlot(Subset_Unclassified, ndims = 30)
dev.off()

Subset_Unclassified <- JackStraw(Subset_Unclassified, num.replicate = 100)
Subset_Unclassified <- ScoreJackStraw(Subset_Unclassified, dims = 1:20)

pdf(file = "5.jackstrawplot.pdf", width = 7.5, height = 5.5)
JackStrawPlot(Subset_Unclassified, dims = 1:20)
dev.off()

Subset_UnclassifiedPC <- 18
Subset_Unclassified <- FindNeighbors(Subset_Unclassified, dims = 1:Subset_UnclassifiedPC, reduction = "harmony")

for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.2, 1.5, 2, 2.5, 3)) {
  Subset_Unclassified <- FindClusters(Subset_Unclassified, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}

apply(Subset_Unclassified@meta.data[, grep("RNA_snn_res", colnames(Subset_Unclassified@meta.data))], 2, table)

p2_tree <- clustree(Subset_Unclassified@meta.data, prefix = "RNA_snn_res.")
pdf(file = "5.clustertree_UnclassifiedPC.pdf", width = 12, height = 10)
p2_tree
dev.off()

Subset_Unclassified <- FindNeighbors(Subset_Unclassified, dims = 1:Subset_UnclassifiedPC, reduction = "harmony")
Subset_Unclassified <- FindClusters(Subset_Unclassified, resolution = 2)
Subset_Unclassified <- RunUMAP(Subset_Unclassified, reduction = "harmony", dims = 1:Subset_UnclassifiedPC)
Subset_Unclassified <- RunTSNE(Subset_Unclassified, reduction = "harmony", dims = 1:Subset_UnclassifiedPC)

# Plot by cell cluster
pdf("5.Unclassified_cluster_plot.pdf", width = 8, height = 6, pointsize = 0.5)
DimPlot(Subset_Unclassified, reduction = "umap", label = TRUE, label.size = 3.5, pt.size = 0.5) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

# Keep only positively differentially expressed genes
Subset_Unclassified.markers <- FindAllMarkers(Subset_Unclassified, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Subset_Unclassified.markers, file = "05.cluster_markers2.csv")
head(Subset_Unclassified.markers)

Subset_Unclassified.markers <- read.csv("05.cluster_markers2.csv")

# Top 5 genes
top5Subset_Unclassified.markers <- Subset_Unclassified.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Color palette
library(viridis)
# Extract colors from the Viridis palette
col <- viridis(50)
# Output palette length
length(col)

pdf(file = "05-cluster.hetmap.pdf", width = 50, height = 30)
DoHeatmap(Subset_Unclassified, 
          features = top5Subset_Unclassified.markers$gene,
          slot = "data",
          group.by = "seurat_clusters",
          group.colors = col) +
  scale_fill_viridis(option = "plasma") +  
  scale_fill_gradient2(low = '#0099CC', mid = 'white', high = '#CC0033', name = 'Z-score')
dev.off()

## Annotation
Subset_Unclassified <- FindClusters(Subset_Unclassified, resolution = 2)

genes <- list("T＆NK cells" = c("CD3E", "CD3D", "CD3G", "CD247", "IL7R", "IL32", "TRAC"),
              "B cells" = c("CD79B", "IGHM", "IGHD", "CD79A", "MS4A1", "IGHG1", "IGLC2", "IGHG4", "IGHG2"),
              "thyrocytes" = c("TG", "EPCAM", "KRT18", "KRT19", "CLU", "FN1", "MGST1", "S100A13"),
              "myeloid cells" = c("LYZ", "S100A8", "S100A9", "CD14", "FCER1G", "TYROBP"),          
              "endothelial cells" = c("PECAM1", "CD34", "CDH5", "VWF", "TIMP3", "RAMP2", "CLDN5", "TFPI", "MGP"),
              "myofibroblasts" = c("ACTA2", "CNN1", "TAGLN"),
              "Fibroblasts" = c("COL1A1", "COL1A2", "COL3A1", "RGS5", "IGFBP7"))

pdf(file = "05.ann_cluster_marker2.pdf", width = 40, height = 19)
do_DotPlot(sample = Subset_Unclassified, features = genes, group.by = "seurat_clusters", dot.scale = 12, colors.use = c("yellow", "red"), legend.length = 50, 
           legend.framewidth = 2, font.size = 12)
dev.off()

ann.ids <- c("myeloid cells",
             "B cells",
             "B cells",
             "myeloid cells",
             "T＆NK cells",
             "myeloid cells",
             "thyrocytes",
             "myeloid cells",
             "T＆NK cells",
             "9",
             "T＆NK cells",
             "myeloid cells",
             "B cells",
             "13",
             "14",
             "15",
             "T＆NK cells",
             "B cells",
             "B cells",
             "19",
             "B cells",
             "21",
             "myeloid cells",
             "23",
             "24",
             "25")

Subset_Unclassifiedidens <- mapvalues(Idents(Subset_Unclassified), from = levels(Idents(Subset_Unclassified)), to = ann.ids)
Idents(Subset_Unclassified) <- Subset_Unclassifiedidens
Subset_Unclassified$cellType <- Idents(Subset_Unclassified)

# Visualize UMAP/tSNE
pdf(file = "5-scRNA.UMAP2.pdf", width = 7.5, height = 5.5)
DimPlot(Subset_Unclassified, reduction = "umap", label = TRUE, label.size = 2.5, pt.size = 0.5) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

pdf(file = "5-scRNA.TSEN2.pdf", width = 7.5, height = 5.5)
DimPlot(Subset_Unclassified, reduction = "tsne", label = TRUE, label.size = 3.5, pt.size = 0.5) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

Idents(Subset_Unclassified) <- Subset_Unclassified@active.ident
Subset_Unclassified.markers <- FindAllMarkers(Subset_Unclassified, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox", slot = "data") 
write.csv(Subset_Unclassified.markers, file = "5.cell_markers3.csv")

top5Subset_Unclassified.markers <- Subset_Unclassified.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Save as PDF
library(viridis)
palette <- viridis(n = 50)
# Plot
pdf(file = "5-cell_marker.hetmap3.pdf", width = 15, height = 10)
DoHeatmap(
  Subset_Unclassified,
  features = top5Subset_Unclassified.markers$gene,
  slot = "data",
  group.colors = palette
) +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_gradient2(low = '#0099CC', mid = 'white', high = '#CC0033', name = 'Z-score')
dev.off()

# Assuming you have completed the re-clustering and the results are stored in Subset_Unclassified@active.ident
# Extract the re-clustered classification results and convert them to character vector
new_cluster_labels <- as.character(Subset_Unclassified@active.ident)

# Get the original data indices that need to be replaced
unclassified_indices <- which(thyroid.combined@meta.data[["cellType"]] == "Unclassified cells")

# Convert the cellType column to a character vector for replacement
thyroid.combined@meta.data[["cellType"]] <- as.character(thyroid.combined@meta.data[["cellType"]])

# Replace the corresponding part in the original data with the new classification results
thyroid.combined@meta.data[["cellType"]][unclassified_indices] <- new_cluster_labels

# Check the replacement results
table(thyroid.combined@meta.data[["cellType"]])
Idents(thyroid.combined) <- thyroid.combined@meta.data[["cellType"]]

pdf(file = "5-ann.scRNA.UMAP3.pdf", width = 10, height = 6)
DimPlot(thyroid.combined, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        legend.position = "right",
        legend.text = element_text(size = 14))
dev.off()

saveRDS(thyroid.combined, file = "thyroid.combined_after_second_annotation.rds")
