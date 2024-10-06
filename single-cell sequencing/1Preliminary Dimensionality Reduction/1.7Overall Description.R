setwd('F:/scRNA/workpath/1.second_annotation')
thyroid.combined <- readRDS("thyroid.combined.delete.rds")

# UMAP plot
setwd('F:/scRNA/workpath/0.allcell')
pdf(file = "ann.scRNA.UMAP.pdf", width = 7.5, height = 5.5)
DimPlot(thyroid.combined, reduction = "umap", label = TRUE, label.size = 3.5, pt.size = 0.1) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

# tSNE plot
pdf(file = "ann.scRNA.TSEN.pdf", width = 7.5, height = 5.5)
DimPlot(thyroid.combined, reduction = "tsne", label = TRUE, label.size = 3.5, pt.size = 0.1) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

# Only keep upregulated differentially expressed genes
Idents(thyroid.combined) <- thyroid.combined@active.ident
thyroid.combined.markers <- FindAllMarkers(thyroid.combined, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox", slot = "data")
write.csv(thyroid.combined.markers, file = "06.cell_markers.csv")

# Extract top 5 highly expressed genes for each cluster
top5thyroid.combined.markers <- thyroid.combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Heatmap plot
pdf(file = "06-cell_marker.hetmap.pdf", width = 12, height = 8)
DoHeatmap(
  thyroid.combined,
  size = 5.5,
  lines.width = 10,
  features = top5thyroid.combined.markers$gene,
  slot = "data"
) +
  scale_fill_gradient2(low = '#0099CC', mid = 'white', high = '#CC0033', name = 'Z-score')
dev.off()

# Extract gene and cluster columns into a new data frame
df <- as.data.frame(thyroid.combined.markers[, c("gene", "cluster")])
table(df$cluster)

# Create a GMT file from the gene sets
dirnb <- paste0(getwd(), "/D2.afgmt.gmt")
gmtlist <- list()
for (afk in as.character(names(table(df$cluster)))) {
  ddgene <- rownames(df)[df$cluster == afk]
  aaname <- afk
  gmtlist[[aaname]] <- ddgene
}

output_gmt <- function(geneset, file) {
  sink(file)
  lapply(names(geneset), function(i) {
    cat(paste(c(i, 'NA', geneset[[i]]), collapse = '\t'))
    cat('\n')
  })
  sink()
}

output_gmt(gmtlist, dirnb)

# VlnPlot for marked genes
afgenes <- c("CD3D", "CD79A", "EPCAM", "PECAM1", "COL1A1","LYZ", "SERPINA1")
pdf(file = "06-cell_FeaturePlot.pdf", width = 12, height = 6)
FeaturePlot(thyroid.combined, features = afgenes, cols = c("grey", "red"), min.cutoff = 1, max.cutoff = 3, ncol = 4, pt.size = 0.5)
dev.off()

pdf(file = "06-cell_VlnPlot.pdf", width = 8, height = 6)
VlnPlot(thyroid.combined, features = afgenes, group.by = "cellType", stack = TRUE) + theme(legend.position = "none")
dev.off()

# Calculate the relative cell frequency for each sample
thyroid.combined@meta.data[["orig.ident"]] <- factor(thyroid.combined@meta.data[["orig.ident"]], levels = c("NT", "pT_WHT1", "pT_WOHT4", "T_WHT1", "T_WHT2", "T_WHT3", "T_WOHT4", "T_WOHT5", "T_WOHT6"))
celltype_frequency <- table(thyroid.combined@meta.data[["orig.ident"]], thyroid.combined$cellType)
write.csv(celltype_frequency, "celltypefrequency.csv")
# Calculate cell type frequency for each sample
celltype_frequency_df <- as.data.frame(celltype_frequency)
names(celltype_frequency_df) <- c("Sample", "Cell_Type", "Frequency")
celltype_frequency_df <- aggregate(Frequency ~ Sample + Cell_Type, data = celltype_frequency_df, FUN = sum)
total_cell_counts <- aggregate(Frequency ~ Sample, data = celltype_frequency_df, FUN = sum)
celltype_frequency_df <- merge(celltype_frequency_df, total_cell_counts, by = "Sample")
celltype_frequency_df$Relative_Frequency <- celltype_frequency_df$Frequency.x / celltype_frequency_df$Frequency.y

# Merge frequency and total cell count data
total_cell_counts <- data.frame(Sample = c("NT", "pT_WHT1", "pT_WOHT4", "T_WHT1", "T_WHT2", "T_WHT3","T_WOHT4", "T_WOHT5", "T_WOHT6"),
                                Frequency = c(2743, 8425, 8972, 5602, 4924, 8740, 6857, 8935, 11155))

celltype_frequency_df <- merge(celltype_frequency_df, total_cell_counts, by = "Sample")
names(celltype_frequency_df) <- c("Sample", "Cell_Type", "Frequency", "Total_Frequency")
celltype_frequency_df$Relative_Frequency <- celltype_frequency_df$Frequency / celltype_frequency_df$Total_Frequency

# Create a barplot using ggplot2
barplot <- ggplot(celltype_frequency_df, aes(x = Sample, y = Relative_Frequency, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Relative Cell Type Frequency in Samples", x = "Sample", y = "Relative Frequency") +
  geom_text(aes(label = Total_Frequency, y = 1.05), vjust = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display and save the barplot
print(barplot)
ggsave("F:/scRNA/workpath/0.allcell/6.barplot(thyroid.combined).pdf", plot = barplot, width = 8, height = 4, dpi = 300)

# Save the workspace
save.image("my_workspace.RData")
saveRDS(thyroid.combined, file = "thyroid.combined.rds")
# Load the workspace back using:
# load("my_workspace.RData")
