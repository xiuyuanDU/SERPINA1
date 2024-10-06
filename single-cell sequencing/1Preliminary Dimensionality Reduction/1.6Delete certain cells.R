setwd('F:/scRNA/workpath/1.second_annotation')

# Read the data after the second annotation
thyroid.combined <- readRDS("thyroid.combined_after_second_annotation.rds")

# View the distribution of cell types in the original data
table(thyroid.combined@meta.data[["cellType"]])

# Remove rows with "doublet" in the cellType column
thyroid.combined <- subset(thyroid.combined, cellType != "9")
thyroid.combined <- subset(thyroid.combined, cellType != "13")
thyroid.combined <- subset(thyroid.combined, cellType != "14")
thyroid.combined <- subset(thyroid.combined, cellType != "15")
thyroid.combined <- subset(thyroid.combined, cellType != "19")
thyroid.combined <- subset(thyroid.combined, cellType != "21")
thyroid.combined <- subset(thyroid.combined, cellType != "23")
thyroid.combined <- subset(thyroid.combined, cellType != "24")
thyroid.combined <- subset(thyroid.combined, cellType != "25")

# Update the cellType column with active.ident values
thyroid.combined@meta.data[["cellType"]] <- thyroid.combined@active.ident

# View the distribution of cell types after removal
table(thyroid.combined@meta.data[["cellType"]])

# Visualize the updated UMAP
pdf(file = "5-ann.scRNA.UMAP4.pdf", width = 8, height = 5)
DimPlot(thyroid.combined, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        legend.position = "right",
        legend.text = element_text(size = 14))
dev.off()

# Save the updated data to a new file
saveRDS(thyroid.combined, file = "thyroid.combined.delete.rds")