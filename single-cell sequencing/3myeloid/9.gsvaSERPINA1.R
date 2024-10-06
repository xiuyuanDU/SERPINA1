# Extract the state information from the mycds object
state_info <- pData(mycds)$State

# Replace the state information
state_info <- recode(state_info, "1" = "state1", "2" = "state2", "3" = "state3")

# Add this state information to the meta.data of scRNAsub
scRNAsub@meta.data$State <- state_info

# Confirm successful addition
head(scRNAsub@meta.data)

# Save the updated scRNAsub object
saveRDS(scRNAsub, file = "scRNAsub.rds")

# Create Subset for Macrophage for GSVA
Subset_Macrophage_GSVA <- scRNAsub
Idents(Subset_Macrophage_GSVA) <- Subset_Macrophage_GSVA@meta.data[["State"]]
table(Idents(Subset_Macrophage_GSVA))

# Load necessary libraries
library(GSVA)
library(pheatmap)
library(GSEABase)
library(limma)

# Load and preprocess the data
exp <- AggregateExpression(Subset_Macrophage_GSVA)[[1]]
exp <- as.matrix(exp)
exp <- exp[rowSums(exp) > 0, ]
print(exp[1:2, 1:2])

# Load the gene set
h_df <- read.gmt("c5.go.v2023.2.Hs.symbols.gmt")[, c(2, 1)]
h_list <- unstack(h_df)

# Filter pathways containing SERPINA1
paths_with_SERPINA1 <- sapply(h_list, function(genes) "SERPINA1" %in% genes)
filtered_h_list <- h_list[which(paths_with_SERPINA1)]

# Check the filtering results
if (length(filtered_h_list) > 0) {
  cat("Paths containing SERPINA1:\n")
  print(names(filtered_h_list))
} else {
  cat("No paths containing SERPINA1 found.\n")
}

# Perform GSVA analysis
ES <- gsva(exp, filtered_h_list)
print(ES[1:41, 1:2])

# Save the GSVA results to a CSV file
write.csv(ES, "11.allDiffGO_SERPINA1.csv")

# Create a heatmap using pheatmap
p <- pheatmap(
  ES,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  breaks = seq(-0.4, 0.4, length.out = 101),
  main = "Pathways involving SERPINA1 in Macrophage",
  border_color = NA,  # Remove clustering lines
  fontsize = 12  # Set font size
)

# Save the heatmap as PNG and PDF files
ggsave("F:/scRNA/workpath/3.myeloid/11SERPINA1_gsva.png", plot = p, width = 13, height = 12, dpi = 300)
ggsave("F:/scRNA/workpath/3.myeloid/11SERPINA1_gsva.pdf", plot = p, width = 13, height = 12, dpi = 300)