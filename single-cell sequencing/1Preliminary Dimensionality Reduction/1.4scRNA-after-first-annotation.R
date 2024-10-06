# Read in the integrated data
setwd('F:/scRNA/workpath/0.allcell')
thyroid.combined = readRDS("thyroid.combined.rds")
Idents(thyroid.combined) <- thyroid.combined@meta.data[["stim"]]

# Correlation between metrics
plot2 <- FeatureScatter(thyroid.combined, feature1 = "nCount_RNA", feature2 = "percent.mt") + RotatedAxis()
plot3 <- FeatureScatter(thyroid.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + RotatedAxis()

# Combined plots
library(patchwork)
pdf(file = "01.corqc.pdf", width = 12, height = 10)
plot2 + plot3 + plot_layout(ncol = 2)
dev.off()

# Violin plot
pdf(file = "01.vlnplot.pdf", width = 20, height = 7)
VlnPlot(thyroid.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Normalization using the LogNormalize method
thyroid.combined <- NormalizeData(thyroid.combined, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify 2000 highly variable genes
thyroid.combined <- FindVariableFeatures(thyroid.combined, selection.method = "vst", nfeatures = 2000)

# Extract the top 10 highly variable genes
top10 <- head(VariableFeatures(thyroid.combined), 10)

# Display highly variable genes
plot1 <- VariableFeaturePlot(thyroid.combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(file = "01.topgene.pdf", width = 7, height = 6)
plot2
dev.off()

# Normalize highly variable genes
thyroid.combined <- ScaleData(thyroid.combined)

### Part 2. Preliminary Dimensionality Reduction
# PCA dimensionality reduction using the top 2000 highly variable genes; you can use 'features' to change the gene set used for dimensionality reduction
thyroid.combined <- Seurat::RunPCA(thyroid.combined, features = VariableFeatures(object = thyroid.combined))
thyroid.combined <- Seurat::RunTSNE(thyroid.combined, dims = 1:20)

### Images
library(ggplot2)
pdf(file = "02.rawtsne.pdf", width = 7.5, height = 5.5)
DimPlot(thyroid.combined, reduction = "tsne", pt.size = 0.5) + 
  theme_classic() + 
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

pdf(file = "02.rawpca.pdf", width = 7.5, height = 5.5)
DimPlot(thyroid.combined, reduction = "pca", pt.size = 0.5) + 
  theme_classic() + 
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

# Samples were grouped according to the "stim" variable and differentiated in the t-SNE plot
colaa = distinctColorPalette(100)
pdf(file = "02.raw.tsne.split.pdf", width = 12, height = 7.5)
do_DimPlot(sample = thyroid.combined,
           plot.title = "",
           reduction = "tsne",
           legend.position = "bottom",
           dims = c(1, 2), split.by = "stim", pt.size = 2)
dev.off()

### Part 3. Optional Harmony for Batch Correction and Dimensionality Reduction
thyroid.combined <- RunHarmony(thyroid.combined, group.by.vars = "stim")
thyroid.combined <- FindVariableFeatures(thyroid.combined, selection.method = "vst", nfeatures = 2000)

# Extract the top 10 highly variable genes
top10 <- head(VariableFeatures(thyroid.combined), 10)

# Display highly variable genes
plot1 <- VariableFeaturePlot(thyroid.combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Output images
pdf(file = "03.topgene.pdf", width = 7, height = 6)
plot2                 
dev.off()

# PCA dimensionality reduction using the top 2000 highly variable genes; you can use 'features' to change the gene set used for dimensionality reduction
pdf(file = "03.harmony.pdf", width = 7.5, height = 5.5)
DimPlot(thyroid.combined, reduction = "harmony", pt.size = 1) + 
  theme_classic() + 
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

thyroid.combined <- Seurat::RunTSNE(thyroid.combined, dims = 1:20, reduction = 'harmony')
pdf(file = "03.tsne.pdf", width = 7.5, height = 5.5)
DimPlot(thyroid.combined, reduction = "tsne", pt.size = 3) + 
  theme_classic() + 
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

# Load color palette
library(colorspace)
collist <- rainbow_hcl(20)
names(collist) = names(table(thyroid.combined$stim))

pdf(file = "03.tsne.split.pdf", width = 12, height = 7.5)
do_DimPlot(sample = thyroid.combined,
           plot.title = "",
           reduction = "tsne",
           legend.position = "bottom",
           dims = c(1, 2), split.by = "stim", pt.size = 2)
dev.off()

# Visualize the first two PC feature genes
VizDimLoadings(thyroid.combined, dims = 1:2, reduction = 'harmony')

# Heatmap visualization of the first 20 PCs
pdf(file = "03.pc_heatmap.pdf", width = 7.5, height = 9)
DimHeatmap(thyroid.combined, dims = 1:20, cells = 1000, balanced = TRUE)
dev.off()

### Determine the number of PCs to use
# Each principal component represents a pattern in the original data
thyroid.combined <- JackStraw(thyroid.combined, num.replicate = 100, dims = 50)
thyroid.combined <- ScoreJackStraw(thyroid.combined, dims = 1:50)

pdf(file = "03.jackstrawplot.pdf", width = 7.5, height = 5.5)
JackStrawPlot(thyroid.combined, dims = 1:50)
dev.off()

pdf(file = "03.ElbowPlot.pdf", width = 5, height = 4)
ElbowPlot(thyroid.combined, ndims = 50)
dev.off()

### Select PCs
thyroid.combinedPC = 30

# Cell Clustering
# Set different resolutions to observe the clustering effect; dims is the number of principal components selected in PCA
thyroid.combined = FindNeighbors(thyroid.combined, dims = 1:thyroid.combinedPC, reduction = "harmony")
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.2, 1.5, 2, 2.5, 3)) {
  thyroid.combined = FindClusters(thyroid.combined, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}

apply(thyroid.combined@meta.data[, grep("RNA_snn_res", colnames(thyroid.combined@meta.data))], 2, table)

p2_tree = clustree(thyroid.combined@meta.data, prefix = "RNA_snn_res.")
pdf(file = "03.clustertree.pdf", width = 12, height = 10)
p2_tree
dev.off()

### Part 4: Clustering and Optimization
# Select resolution for dimensionality reduction
thyroid.combined = FindNeighbors(thyroid.combined, dims = 1:thyroid.combinedPC, reduction = "harmony")
thyroid.combined <- FindClusters(thyroid.combined, resolution = 2)

# View cluster IDs
head(Idents(thyroid.combined), 5)

# Check how many cells are in each category
head(thyroid.combined@meta.data)
table(thyroid.combined@meta.data$seurat_clusters)

# only.pos: only retain upregulated differentially expressed genes
thyroid.combined.markers <- FindAllMarkers(thyroid.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(thyroid.combined.markers, file = "04.cluster_markers.csv")
head(thyroid.combined.markers)
thyroid.combined.markers <- read.csv("04.cluster_markers.csv")

# Top 5 genes
top5thyroid.combined.markers <- thyroid.combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Color palette
library(viridis)
# Extract colors from the Viridis palette
col <- viridis(50)
# Output length of the color palette
length(col)

pdf(file = "04-cluster.hetmap.pdf", width = 50, height = 30)
DoHeatmap(thyroid.combined, 
          features = top5thyroid.combined.markers$gene,
          slot = "data",
          group.by = "seurat_clusters",
          group.colors = col) +
  scale_fill_viridis(option = "plasma") +  
  scale_fill_gradient2(low = '#0099CC', mid = 'white', high = '#CC0033',
                       name = 'Z-score')
dev.off()

### Part 5: Visualizing Cells in Low-Dimensional Space with UMAP/tSNE
thyroid.combined <- RunUMAP(thyroid.combined, dims = 1:thyroid.combinedPC, reduction = "harmony")
thyroid.combined <- RunTSNE(thyroid.combined, dims = 1:thyroid.combinedPC, reduction = "harmony")

# Visualize UMAP/tSNE
pdf(file = "05-cluster.UMAP.pdf", width = 12.5, height = 10.5)
DimPlot(thyroid.combined, reduction = "umap", label = TRUE, label.size = 3.5, pt.size = 0.5) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

pdf(file = "05-cluster.TSNE.pdf", width = 6.5, height = 5.5)
DimPlot(thyroid.combined, reduction = "tsne", label = TRUE, label.size = 3.5, pt.size = 0.5) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

# saveRDS(thyroid.combined, file = "thyroid.combined_UMAP.rds")
# thyroid.combined <- JoinLayers(thyroid.combined)
saveRDS(thyroid.combined, file = "thyroid.combined_before.annotation.rds")
thyroid.combined = readRDS("thyroid.combined_before.annotation.rds")

### Cell Annotation ###
### SingleR Cell Annotation
refdata = celldex::HumanPrimaryCellAtlasData()   
refdata
head(colnames(refdata))
head(rownames(refdata))


# Check how many types of cells there are
unique(refdata@colData@listData[["label.main"]])

# The data used is normalized data
# testdata <- GetAssayData(thyroid.combined, slot="data") 
# In Seurat version 5.0.0, the GetAssayData function no longer supports using the slot parameter to retrieve multi-layer data,
# instead, the layer parameter should be used.
testdata <- GetAssayData(thyroid.combined, layer = "data")

dim(testdata)
testdata[1:30, 1:4]
clusters <- thyroid.combined@meta.data$seurat_clusters
table(clusters)

cellpred <- SingleR(test = testdata,  
                    ref = refdata, 
                    labels = refdata$label.main,
                    method = "cluster", 
                    clusters = clusters,
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")

str(cellpred, max.level = 3)
metadata <- cellpred@metadata
head(metadata)

celltype = data.frame(ClusterID = rownames(cellpred), 
                      celltype = cellpred$labels, 
                      stringsAsFactors = FALSE)
celltype
write.csv(celltype, "05.singleR.celltype_anno_SingleR.csv")

# The annotation results on the score heatmap need to be corrected
pdf(file = "05-singleR.pdf", width = 10, height = 8)
p = plotScoreHeatmap(cellpred, clusters = rownames(cellpred), order.by = "cluster")
p
dev.off()


## Visualization of results after SingleR annotation
newLabels = cellpred$labels
names(newLabels) = levels(thyroid.combined)
thyroid.combined = RenameIdents(thyroid.combined, newLabels)
thyroid.combined$cellType = Idents(thyroid.combined)

# Visualize UMAP/tSNE
pdf(file = "05-scRNA.UMAP.pdf", width = 7.5, height = 5.5)
DimPlot(thyroid.combined, reduction = "umap", label = TRUE, label.size = 2.5, pt.size = 0.5) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

pdf(file = "05-scRNA.TSNE.pdf", width = 7.5, height = 5.5)
DimPlot(thyroid.combined, reduction = "tsne", label = TRUE, label.size = 3.5, pt.size = 0.5) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()


### Manual Annotation
thyroid.combined <- FindClusters(thyroid.combined, resolution = 2)

# List of marker genes for different cell types
genes <- list("T＆NK cells" = c("CD3E", "CD3D", "CD3G", "CD247", "IL7R", "IL32", "TRAC"),
              "B cells" = c("CD79B", "IGHM", "IGHD", "CD79A", "MS4A1", "IGLC2", "JCHAIN"),
              "thyrocytes" = c("TG", "EPCAM", "KRT18", "KRT19", "CLU", "FN1", "MGST1", "S100A13"),
              "myeloid cells" = c("LYZ", "S100A8", "S100A9", "CD14", "FCER1G", "TYROBP"),          
              "endothelial cells" = c("PECAM1", "CD34", "CDH5", "VWF", "TIMP3", "RAMP2", "CLDN5", "TFPI", "MGP"),
              "myofibroblasts" = c("ACTA2", "CNN1", "TAGLN"),
              "Fibroblasts" = c("COL1A1", "COL1A2", "COL3A1", "RGS5", "IGFBP7"))

# Create a dot plot for the marker genes
pdf(file = "05.ann_cluster_marker.pdf", width = 40, height = 19)
do_DotPlot(sample = thyroid.combined, features = genes, group.by = "seurat_clusters", 
           dot.scale = 12, colors.use = c("yellow", "red"), legend.length = 50,
           legend.framewidth = 2, font.size = 12)
dev.off()

# Define the annotations for each cluster
ann.ids <- c("T＆NK cells",
             "T＆NK cells",
             "thyrocytes",
             "T＆NK cells",
             "B cells",
             "B cells",
             "thyrocytes",
             "myeloid cells",
             "T＆NK cells",
             "T＆NK cells",
             "thyrocytes",
             "myofibroblasts",
             "thyrocytes",  # 12
             "Unclassified cells",  # 13
             "thyrocytes",
             "myeloid cells",  # 15
             "T＆NK cells",
             "T＆NK cells",
             "T＆NK cells",
             "myeloid cells",
             "endothelial cells",  # 20
             "T＆NK cells",
             "thyrocytes",
             "Unclassified cells",  # 23
             "thyrocytes",
             "T＆NK cells",  # 25
             "thyrocytes",
             "myeloid cells",
             "Unclassified cells",  # 28
             "thyrocytes",
             "T＆NK cells",
             "T＆NK cells",
             "Unclassified cells",  # 32
             "Fibroblasts",
             "myeloid cells",  # 34
             "myeloid cells",
             "T＆NK cells",
             "endothelial cells",
             "Unclassified cells",  # 38
             "Unclassified cells",
             "T＆NK cells",
             "endothelial cells",
             "thyrocytes",
             "Unclassified cells",  # 43
             "thyrocytes",  # 44
             "thyrocytes",
             "Unclassified cells"  # 46
)

# Map the annotations to the clusters
thyroid.combinedidens = mapvalues(Idents(thyroid.combined), from = levels(Idents(thyroid.combined)), to = ann.ids)
Idents(thyroid.combined) = thyroid.combinedidens
thyroid.combined$cellType = Idents(thyroid.combined)

### Visualization of results after manual annotation
# Visualize UMAP/tSNE
pdf(file = "05-ann.scRNA.UMAP.pdf", width = 7.5, height = 5.5)
DimPlot(thyroid.combined, reduction = "umap", label = TRUE, label.size = 3.5, pt.size = 0.1) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

pdf(file = "05-ann.scRNA.TSNE.pdf", width = 7.5, height = 5.5)
DimPlot(thyroid.combined, reduction = "tsne", label = TRUE, label.size = 3.5, pt.size = 0.1) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

# Save the combined object after annotation
saveRDS(thyroid.combined, file = "thyroid.combined_after_annotation.rds")

# Load the annotated data
thyroid.combined = readRDS("thyroid.combined_after_annotation.rds")
