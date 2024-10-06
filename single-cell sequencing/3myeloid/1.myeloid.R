# Save the processed data for future use
saveRDS(Subset_myeloid , file = "Subset_myeloid.rds")

# Change the working directory and load the combined data
setwd('F:/scRNA/workpath/0.allcell')
thyroid.combined = readRDS("thyroid.combined.rds")

# Subset the myeloid cells from the combined data
setwd('F:/scRNA/workpath/3.myeloid')
Subset_myeloid <- subset(thyroid.combined, idents = "myeloid cells")

# Normalize the data and find variable features
Subset_myeloid <- ScaleData(Subset_myeloid)
Subset_myeloid <- FindVariableFeatures(Subset_myeloid, selection.method = "vst", nfeatures = 2000)

# Perform dimensionality reduction with PCA, t-SNE, and Harmony
Subset_myeloid <- RunPCA(Subset_myeloid, features = VariableFeatures(object = Subset_myeloid))
Subset_myeloid <- RunTSNE(Subset_myeloid, dims = 1:20)
Subset_myeloid <- RunHarmony(Subset_myeloid, group.by.vars = "stim")

# Create an Elbow plot to determine the optimal number of dimensions
pdf(file = "10.ElbowPlot_myeloid_ElbowPlot.pdf", width = 5, height = 4)
ElbowPlot(Subset_myeloid, ndims = 30)
dev.off()

Subset_myeloidPC = 20
Subset_myeloid = FindNeighbors(Subset_myeloid, dims = 1:Subset_myeloidPC, reduction = "harmony")

# Perform clustering with different resolutions and visualize clustering
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.2, 1.5, 2, 2.5, 3)) {
  Subset_myeloid = FindClusters(Subset_myeloid, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
apply(Subset_myeloid@meta.data[, grep("RNA_snn_res", colnames(Subset_myeloid@meta.data))], 2, table)

p2_tree = clustree(Subset_myeloid@meta.data, prefix = "RNA_snn_res.")
pdf(file = "10.clustertree_myeloid.pdf", width = 12, height = 10)
p2_tree
dev.off()

Subset_myeloid = FindNeighbors(Subset_myeloid, dims = 1:Subset_myeloidPC, reduction = "harmony")
Subset_myeloid <- FindClusters(Subset_myeloid, resolution = 1)
Subset_myeloid <- RunUMAP(Subset_myeloid, reduction = "harmony", dims = 1:Subset_myeloidPC)
Subset_myeloid <- RunTSNE(Subset_myeloid, reduction = "harmony", dims = 1:Subset_myeloidPC)

Idents(Subset_myeloid) = Subset_myeloid@meta.data[["RNA_snn_res.1"]]

# Plot cluster visualization using UMAP
pdf("10.myeloid_cluster_plot.pdf", width = 6.5, height = 4.5, pointsize = 0.5)
palette("Dark2")
DimPlot(Subset_myeloid, reduction = "umap", label = T, label.size = 3.5, pt.size = 0.5) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

# Plot UMAP visualization by groups
pdf("10.myeloid_umap_plot.pdf", width = 8, height = 6, pointsize = 0.5)
palette("Dark2")
DimPlot(Subset_myeloid, reduction = "umap", group.by = "stim")
dev.off()

# Load color palette and plot UMAP split by "stim"
library(colorspace)
collist <- rainbow_hcl(20)
names(collist) = names(table(Subset_myeloid$stim))
pdf(file = "10.stim.split.pdf", width = 12, height = 7.5)
do_DimPlot(sample = Subset_myeloid, plot.title = "", reduction = "tsne",
           legend.position = "bottom", dims = c(1, 2), split.by = "stim", pt.size = 2)
dev.off()

# Extract and preprocess the "stim" column to create a new indicator "tissue"
stim <- Subset_myeloid$stim
tissue <- gsub("T_WHT\\d+", "T_WHT", stim)
tissue <- gsub("pT_WHT\\d+", "pT_WHT", stim)
tissue <- gsub("pT_WHT\\d+", "pT_WOHT", stim)
tissue <- gsub("T_WOHT\\d+", "T_WOHT", tissue)
tissue <- gsub("\\d+", "", tissue)
Subset_myeloid$tissue <- tissue

# Plot UMAP visualization by "tissue"
pdf("10Subset_myeloid_tissue_plot.pdf", width = 8, height = 6, pointsize = 0.5)
palette("Dark2")
DimPlot(Subset_myeloid, reduction = "umap", group.by = "tissue")
dev.off()

# Save the Cluster Markers information
write.csv(Subset_myeloid.markers, file = "10.cluster_markers.csv")
head(Subset_myeloid)

# Identify and visualize the top 5 genes for each cluster
top5Subset_myeloid.markers <- Subset_myeloid.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Set the color palette for Heatmap
library(viridis)
col <- viridis(50)

# Create a heatmap visualization
pdf(file = "10-cluster.hetmap.pdf", width = 18, height = 10)
DoHeatmap(Subset_myeloid, 
          features = top5Subset_myeloid.markers$gene,
          slot = "data",
          group.by = "seurat_clusters",
          group.colors = col) +
  scale_fill_viridis(option = "plasma") +  
  scale_fill_gradient2(low = '#0099CC', mid = 'white', high = '#CC0033', name = 'Z-score')
dev.off()

# Save the final processed data and workspace
saveRDS(Subset_myeloid , file = "Subset_myeloid.rds")
save.image("my_workspace_myeloid.RData")

