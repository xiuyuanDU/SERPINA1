Idents(Subset_thyrocytes) = Subset_thyrocytes@meta.data[["orig.ident"]]
Subset_thyrocytes <- subset(Subset_thyrocytes, idents = c("NT", "T_WHT3", "T_WOHT5", "T_WOHT6"))
# Create a vector containing the desired samples
subset_samples <- c("NT", "T_WHT3", "T_WOHT5", "T_WOHT6")
# Use the match function to determine which samples in meta.data match subset_samples
match_samples <- match(Subset_thyrocytes@meta.data[["orig.ident"]], subset_samples, nomatch = 0)
# Assign matching sample names directly to the sample column
Subset_thyrocytes@meta.data[["orig.ident"]] <- subset_samples[match_samples]

# Normalization
Subset_thyrocytes <- ScaleData(Subset_thyrocytes)
# Find differentially expressed genes
Subset_thyrocytes <- FindVariableFeatures(Subset_thyrocytes, selection.method = "vst", nfeatures = 2000)
# Dimension reduction
Subset_thyrocytes <- RunPCA(Subset_thyrocytes, features = VariableFeatures(object = Subset_thyrocytes))
Subset_thyrocytes <- RunTSNE(Subset_thyrocytes, dims = 1:20)
Subset_thyrocytes <- RunHarmony(Subset_thyrocytes, group.by.vars = "stim")

# Save Elbow Plot
pdf(file = "07.ElbowPlot_thyrocytes.pdf", width = 5, height = 4)
ElbowPlot(Subset_thyrocytes, ndims = 30)
dev.off()

# JackStraw analysis
Subset_thyrocytes <- JackStraw(Subset_thyrocytes, num.replicate = 100)
Subset_thyrocytes <- ScoreJackStraw(Subset_thyrocytes, dims = 1:20)
pdf(file = "07.jackstrawplot.pdf", width = 7.5, height = 5.5)
JackStrawPlot(Subset_thyrocytes, dims = 1:20)
dev.off()

Subset_thyrocytesPC=12
Subset_thyrocytes <- FindNeighbors(Subset_thyrocytes, dims = 1:Subset_thyrocytesPC, reduction = "harmony")
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.2, 1.5, 2, 2.5, 3)) {
  Subset_thyrocytes <- FindClusters(Subset_thyrocytes, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
apply(Subset_thyrocytes@meta.data[, grep("RNA_snn_res", colnames(Subset_thyrocytes@meta.data))], 2, table)

p2_tree <- clustree(Subset_thyrocytes@meta.data, prefix = "RNA_snn_res.")
pdf(file = "07.clustertree_thyrocytes.pdf", width = 12, height = 10)
p2_tree
dev.off()

Subset_thyrocytes <- FindNeighbors(Subset_thyrocytes, dims = 1:Subset_thyrocytesPC, reduction = "harmony")
Subset_thyrocytes <- FindClusters(Subset_thyrocytes, resolution = 0.8)
Subset_thyrocytes <- RunUMAP(Subset_thyrocytes, reduction = "harmony", dims = 1:Subset_thyrocytesPC)
Subset_thyrocytes <- RunTSNE(Subset_thyrocytes, reduction = "harmony", dims = 1:Subset_thyrocytesPC)

# Cluster-based plot
pdf("07thyrocytes_cluster_plot.pdf", width = 6.5, height = 4.5, pointsize = 0.5)
DimPlot(Subset_thyrocytes, reduction = "umap", label = TRUE, label.size = 3.5, pt.size = 0.5) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

# Differential gene expression
Subset_thyrocytes.markers <- FindAllMarkers(Subset_thyrocytes, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Subset_thyrocytes.markers, file = "07.thyrocytes_cluster_markers.csv")

# Top 5 genes
top5Subset_thyrocytes.markers <- Subset_thyrocytes.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
top5Subset_thyrocytes.markers[, c("gene", "cluster")]

# Color palette
library(viridis)
col <- viridis(50)
# Output palette length
length(col)
pdf(file = "07-cluster.hetmap.pdf", width = 18, height = 10)
DoHeatmap(Subset_thyrocytes, features = top5Subset_thyrocytes.markers$gene,
          slot = "data", group.by = "seurat_clusters", group.colors = col) +
  scale_fill_viridis(option = "plasma") +
  scale_fill_gradient2(low = '#0099CC', mid = 'white', high = '#CC0033', name = 'Z-score')
dev.off()

# Frequency heatmap of thyroid cell clusters by tissue type
table(Subset_thyrocytes$orig.ident, Subset_thyrocytes$seurat_clusters)
# Data for the table
data <- matrix(c(
  "0", "0", "0", "0", "0", "0", "0", "0", "0", "295", "0", "0", "0",
  "239", "135", "938", "12", "928", "580", "143", "409", "26", "0", "121", "20", "46",
  "1307", "46", "363", "1238", "125", "177", "30", "121", "31", "0", "9", "43", "51",
  "27", "1362", "13", "2", "101", "31", "526", "110", "285", "0", "26", "51", "5"
), nrow = 4, byrow = TRUE)

# Set row names and column names
rownames(data) <- c("NT", "T_WHT3", "T_WOHT5", "T_WOHT6")
colnames(data) <- 0:12

# Convert to a data frame
df <- as.data.frame.table(data)
df$Var1 <- factor(df$Var1, levels = c("NT", "T_WHT3", "T_WOHT5", "T_WOHT6"))
df$Var2 <- factor(df$Var2)

# Heatmap plot
p1 <- ggplot(data = df, aes(x = Var2, y = Var1, fill = as.numeric(Freq))) +
  geom_tile() +
  scale_fill_gradientn(colours = c("white", "blue")) +
  labs(x = "Seurat Clusters", y = "Sample", title = "Heatmap of Cell Counts by Sample and Seurat Cluster") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), text = element_text(size = 16))

# Save as PDF file
pdf(file = "07.thyrocytes.matrix.pdf", width = 11.2, height = 5)
print(p1)
dev.off()

# Specify the gene list to extract
selected_genes <- c("TG", "TPO", "SLC26A4", "DIO2", "TSHR", "PAX8", "DUOX1", "DUOX2", "NKX2-1", "GLIS3", "FOXE1", "TFF3", "FHL1")
# Extract gene expression data
# Check which genes in selected_genes do not exist in Subset_thyrocytes
genes_in_data <- intersect(selected_genes, rownames(Subset_thyrocytes))

# Extract gene expression data present in Subset_thyrocytes
gene_expression_matrix <- GetAssayData(object = Subset_thyrocytes, slot = "data")[genes_in_data, ]
gene_expression_matrix <- t(gene_expression_matrix)
# Compute TDS scores
tds_scores <- rowSums(log2(gene_expression_matrix + 1))
# Extract cell cluster labels
seurat_clusters <- Subset_thyrocytes$seurat_clusters
# Create a data frame
data_df <- data.frame(TDS_score = tds_scores, seurat_clusters = seurat_clusters)
# Save TDS scores to Subset_Cells dataset
Subset_thyrocytes$TDS_score <- tds_scores
# Plot a boxplot
pdf(file = "07heatmap.pdf", width = 6.2, height = 4)
boxplot(TDS_score ~ seurat_clusters, data = data_df, 
        xlab = "thyrocytes clusters", ylab = "TDS score", 
        main = "Boxplot of TDS score for different thyrocytes clusters",
        notch = FALSE, outline = TRUE, col = "lightblue")
##############################################################################################




Subset_thyrocytes <- FindClusters(Subset_thyrocytes, resolution = 0.8)

# Define annotation IDs
ann.ids <- c("Moderately",
             "Poorly",
             "Moderately",
             "Moderately",
             "Moderately",
             "Moderately",
             "Poorly",
             "Moderately",
             "Poorly",
             "Highly",
             "Poorly",
             "stem",
             "Highly")

# Map cluster identities to annotations
Subset_thyrocytesidens = mapvalues(Idents(Subset_thyrocytes), from = levels(Idents(Subset_thyrocytes)), to = ann.ids)
Idents(Subset_thyrocytes) = Subset_thyrocytesidens

# Assign the new identities to the cellType column in metadata
Subset_thyrocytes$cellType = Idents(Subset_thyrocytes)
Subset_thyrocytes@meta.data[["cellType"]] <- factor(Subset_thyrocytes@meta.data[["cellType"]], levels = c("Highly", "Moderately", "Poorly", "stem"))

# Update identities with cell types
Idents(Subset_thyrocytes) <- Subset_thyrocytes@meta.data[["cellType"]]

### Visualization after manual annotation
# Visualize UMAP/tSNE
pdf(file = "07-ann.scRNA.SERPINA1.UMAP.pdf", width = 6.5, height = 4.5)
DimPlot(Subset_thyrocytes, reduction = "umap", label = TRUE, label.size = 3.5, pt.size = 0.5) + 
  theme_classic() + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

pdf(file = "07-ann.scRNA.SERPINA1.TSNE.pdf", width = 6.5, height = 4.5)
DimPlot(Subset_thyrocytes, reduction = "tsne", label = TRUE, label.size = 3.5, pt.size = 0.5) + 
  theme_classic() + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

########################################################################
### Part 6: Retain only upregulated differentially expressed genes
Idents(Subset_thyrocytes) <- Subset_thyrocytes@active.ident
Subset_thyrocytes.markers <- FindAllMarkers(Subset_thyrocytes, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox", slot = "data") 
write.csv(Subset_thyrocytes.markers, file = "08.Subset_thyrocytescell_markers.csv")

# Extract top 10 highly expressed genes for each category
top5Subset_thyrocytes.markers <- Subset_thyrocytes.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Plot heatmap
pdf(file = "08-Subset_thyrocytescell_marker.heatmap.pdf", width = 12, height = 8)
DoHeatmap(
  Subset_thyrocytes,
  size = 5.5,
  lines.width = 10,
  features = top5Subset_thyrocytes.markers$gene,
  slot = "data"
) +
  scale_fill_gradient2(low = '#0099CC', mid = 'white', high = '#CC0033', name = 'Z-score')
dev.off()

# Save the object and workspace
saveRDS(Subset_thyrocytes, file = "Subset_thyrocytes.rds")
save.image("my_workspace_thyrocytes.RData")

# Load saved workspace
load("my_workspace_thyrocytes.RData")
# Subset_thyrocytes = readRDS("Subset_thyrocytes.rds")

                                                                                                          
                                                                                                         
