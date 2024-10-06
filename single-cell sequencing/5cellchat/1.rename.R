setwd('F:/scRNA/workpath/5.cellchat')
thyroid.combined = readRDS("thyroid.combined.rds")
Subset_thyrocytes = readRDS("Subset_thyrocytes.rds")
Subset_lymphocytes = readRDS("Subset_lymphocytes.rds")
Subset_myeloid = readRDS("Subset_myeloid.rds")
Subset_CD = readRDS("Subset_DC.rds")
scRNAsub = readRDS("scRNAsub.rds")

##########################################################
# Re-cluster and store the results in Subset_Unclassified@active.ident
# First, extract the new clustering results and convert them to a character vector
new_cluster_labels <- as.character(Subset_thyrocytes@active.ident)

# Get the indices of the original data that need to be replaced
thyrocytes_indices <- which(thyroid.combined@meta.data[["cellType"]] == "thyrocytes")

# Convert the cellType column to a character vector for replacement
thyroid.combined@meta.data[["cellType"]] <- as.character(thyroid.combined@meta.data[["cellType"]])

# Use the new classification results to replace the corresponding parts of the original data
thyroid.combined@meta.data[["cellType"]][thyrocytes_indices] <- new_cluster_labels[1:length(thyrocytes_indices)]

# Check the replacement results
table(thyroid.combined@meta.data[["cellType"]])
Idents(thyroid.combined) <- thyroid.combined@meta.data[["cellType"]]
pdf(file = "13.1-ann.scRNA.UMAP3.pdf", width = 10, height = 6)
DimPlot(thyroid.combined, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        legend.position = "right",
        legend.text = element_text(size = 14))
dev.off()

##########################################################
# Assuming you have completed re-clustering and stored the results in Subset_Unclassified@active.ident
# First, extract the new clustering results and convert them to a character vector
new_cluster_labels <- as.character(Subset_lymphocytes@active.ident)

# Get the indices of the original data that need to be replaced
T_indices <- which(thyroid.combined@meta.data[["cellType"]] == "T＆NK cells")
B_indices <- which(thyroid.combined@meta.data[["cellType"]] == "B cells")
lymphocytes_indices <- c(T_indices, B_indices)
# Convert the cellType column to a character vector for replacement
thyroid.combined@meta.data[["cellType"]] <- as.character(thyroid.combined@meta.data[["cellType"]])

# Use the new classification results to replace the corresponding parts of the original data
thyroid.combined@meta.data[["cellType"]][lymphocytes_indices] <- new_cluster_labels

# Check the replacement results
table(thyroid.combined@meta.data[["cellType"]])
Idents(thyroid.combined) <- thyroid.combined@meta.data[["cellType"]]
pdf(file = "13.2-ann.scRNA.UMAP3.pdf", width = 10, height = 6)
DimPlot(thyroid.combined, reduction = "umap", label = FALSE, label.size = 4, pt.size = 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        legend.position = "right",
        legend.text = element_text(size = 14))
dev.off()

##########################################################
# Assuming you have completed re-clustering and stored the results in Subset_myeloid@active.ident
# First, extract the new clustering results and convert them to a character vector
new_cluster_labels <- as.character(Subset_myeloid@active.ident)

# Get the indices of the original data that need to be replaced
myeloid_indices <- which(thyroid.combined@meta.data[["cellType"]] == "myeloid cells")

# Convert the cellType column to a character vector for replacement
thyroid.combined@meta.data[["cellType"]] <- as.character(thyroid.combined@meta.data[["cellType"]])

# Use the new classification results to replace the corresponding parts of the original data
thyroid.combined@meta.data[["cellType"]][myeloid_indices] <- new_cluster_labels

# Check the replacement results
table(thyroid.combined@meta.data[["cellType"]])
Idents(thyroid.combined) <- thyroid.combined@meta.data[["cellType"]]
pdf(file = "13.2-ann.scRNA.UMAP4.pdf", width = 10, height = 6)
DimPlot(thyroid.combined, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        legend.position = "right",
        legend.text = element_text(size = 14))
dev.off()

##########################################################
# Re-cluster and store the results in Subset_Unclassified@active.ident
# First, extract the new clustering results and convert them to a character vector
new_cluster_labels <- as.character(scRNAsub@meta.data[["State"]])

# Get indices that meet both conditions
scRNAsub_indices <- which(thyroid.combined@meta.data[["stim"]] %in% c("tumor_Without_HT5", "tumor_Without_HT6", "tumor_With_HT1", "tumor_With_HT3") &
                            thyroid.combined@meta.data[["cellType"]] %in% c("Ms-C0", "Ms-C1", "Ms-C3", "Ms-C4", "Ms-C5", "Ms-C6", "Ms-C7", "Ms-C9", "Ms-C12", "PMs-C11", "Ms-C14", "Ms-C15"))

# Check index lengths
cat("scRNAsub_indices length:", length(scRNAsub_indices), "\n")
cat("new_cluster_labels length:", length(new_cluster_labels), "\n")

# Convert the cellType column to a character vector for replacement
thyroid.combined@meta.data[["cellType"]] <- as.character(thyroid.combined@meta.data[["cellType"]])

# Use the new classification results to replace the corresponding parts of the original data
if (length(new_cluster_labels) >= length(scRNAsub_indices)) {
  thyroid.combined@meta.data[["cellType"]][scRNAsub_indices] <- new_cluster_labels[1:length(scRNAsub_indices)]
} else {
  stop("Error: The length of new_cluster_labels is less than the length of scRNAsub_indices.")
}

# Check the replacement results
print(table(thyroid.combined@meta.data[["cellType"]]))

# Update the object's identifiers
Idents(thyroid.combined) <- thyroid.combined@meta.data[["cellType"]]

# Plot and save
pdf(file = "13.m-ann.scRNA.UMAP3.pdf", width = 10, height = 6)
DimPlot(thyroid.combined, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.01) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        legend.position = "right",
        legend.text = element_text(size = 14))
dev.off()

thyroid.combined <- subset(thyroid.combined, idents = c("endothelial cells", "myofibroblasts", "Highly", "Treg", "B-C1", 
                                                       "Fibroblasts", "CD8-C2", "Tm-C3", "CD8-C3", "Tm-C1", "NK", 
                                                       "Plasma cells", "Moderately", "Poorly", "stem", "Tm-C7", "Tm-C4", 
                                                       "Tm-C2", "B-C2", "CD8-C1", "Tm-C5", "B-C3", "pDCs-C10", "Tm-C6", 
                                                       "CD8-C5", "CD8-C4", "B-C4", "mDCs-C2", "LAMP_mDCs-C13",  
                                                       "mDCs-C8","state1","state2","state3"))
print(table(thyroid.combined@meta.data[["cellType"]]))
saveRDS(thyroid.combined , file = "cellchat.rds")
