setwd('F:/scRNA/workpath/3.myeloid')
Subset_myeloid = readRDS("Subset_myeloid.rds")
Idents(Subset_myeloid) = Subset_myeloid@meta.data[["cellType"]]
Subset_myeloid <- subset(Subset_myeloid, idents = c("Ms-C0", "Ms-C1", "Ms-C3", "Ms-C4", "Ms-C5", "Ms-C6", 
                                                   "Ms-C7", "Ms-C9", "PMs-C11", "Ms-C12", "Ms-C14", "Ms-C15"))
Idents(Subset_myeloid) <- Subset_myeloid@meta.data[["stim"]]
# Filter the data, keeping only the target cell clusters
Subset_DGEs <- subset(Subset_myeloid, idents = c("tumor_Without_HT4", "tumor_Without_HT5", "tumor_Without_HT6", "tumor_With_HT1", "tumor_With_HT3"))
# Calculate differentially expressed genes for each cell cluster
Subset_DGEs <- FindMarkers(Subset_myeloid, ident.1 = c("tumor_Without_HT6"), ident.2 = c("tumor_Without_HT4", "tumor_With_HT3", "tumor_Without_HT5", "tumor_With_HT1"), test.use = 'wilcox', min.pct = 0.1)
write.csv(Subset_DGEs, file = "11.Subset_myeloid.csv", row.names = TRUE)
# Manually add the missing gene column name "gene_symbol" in the csv file, and read the table again
Subset_DGEs <- read.csv("11.Subset_myeloid.csv")
# Get differentially expressed genes
DE_genes <- Subset_DGEs %>%
  filter(p_val_adj < 0.05) %>%
  pull(gene_symbol)
# Output the list of differentially expressed genes
print(DE_genes)
library(ggrepel)
cluster1.markers <- Subset_DGEs %>%
  mutate(Difference = pct.1 - pct.2) %>% 
  rownames_to_column("gene")
# Custom threshold
log2FC = 1
padj = 0.05 
cluster1.markers$threshold = "ns"
cluster1.markers[which(cluster1.markers$avg_log2FC > log2FC & cluster1.markers$p_val_adj < padj),]$threshold = "up"
cluster1.markers[which(cluster1.markers$avg_log2FC < (-log2FC) & cluster1.markers$p_val_adj < padj),]$threshold = "down"
cluster1.markers$threshold = factor(cluster1.markers$threshold, levels = c('down', 'ns', 'up'))
pdf("11.volcano_plot_myeloid.pdf", width = 8, height = 6) 
ggplot(cluster1.markers, aes(x = Difference, y = avg_log2FC, color = threshold)) + 
  geom_point(size = 0.5) + 
  scale_color_manual(values = c("blue", "grey", "red")) + 
  geom_label_repel(data = subset(cluster1.markers, avg_log2FC >= 1 & Difference >= 0.2 & p_val_adj <= 0.05), 
                   aes(label = gene_symbol),  # Add label
                   color = "black", # Set the color of the labels in the label
                   segment.colour = "black", # Set the color of the label segment
                   label.padding = 0.1, 
                   # max.overlaps = 200,
                   segment.size = 0.3,  # Size of the segment
                   size = 4) +
  geom_label_repel(data = subset(cluster1.markers, avg_log2FC <= -1 & Difference <= -0.2 & p_val_adj <= 0.05), 
                   aes(label = gene_symbol), label.padding = 0.1, 
                   color = "black",
                   segment.colour = "black",
                   segment.size = 0.3, size = 4) +
  geom_vline(xintercept = 0.0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Difference",
       y = "avg_log2FC", 
       title = "Volcano Plot of DE Genes"
  )
theme_classic()
dev.off()

