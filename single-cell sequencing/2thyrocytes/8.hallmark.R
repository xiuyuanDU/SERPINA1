setwd('F:/scRNA/workpath/2.thyrocytes')

# Load the required libraries
library(GSVA)
library(GSEABase)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load the dataset
Subset_thyrocytes <- readRDS("Subset_thyrocytes.rds")

# Set identities
Idents(Subset_thyrocytes) <- Subset_thyrocytes@meta.data[["cellType"]]

# Subset the data based on conditions
Subset_thyrocytes_GSVA <- Subset_thyrocytes 

# Display the identity table
print(table(Idents(Subset_thyrocytes_GSVA)))

# Aggregate expression data for GSVA
exp <- AggregateExpression(Subset_thyrocytes_GSVA)[[1]]
exp <- as.matrix(exp)
exp <- exp[rowSums(exp) > 0,]
print(exp[1:2, 1:2])

# Load gene sets
h_df <- read.gmt("Hallmarkgenesets.v2023.2.Hs.gmt")[, c(2, 1)]
h_list <- unstack(h_df)

# Filter pathways that contain SERPINA1
paths_with_SERPINA1 <- sapply(h_list, function(genes) "SERPINA1" %in% genes)
filtered_h_list <- h_list[which(paths_with_SERPINA1)]

# Check filtering results
if (length(filtered_h_list) > 0) {
  cat("Paths containing SERPINA1:\n")
  print(names(filtered_h_list))
} else {
  cat("No paths containing SERPINA1 found.\n")
}

# Perform GSVA analysis
ES <- gsva(exp, h_list)
print(ES[1:50, 1:2])

# Save GSVA results to a CSV file
write.csv(ES, "HALLMARK.allDiffGO_SERPINA1.csv")
library(pheatmap)
p <- pheatmap(
  ES,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  breaks = seq(-0.4, 0.4, length.out = 101),
  main = "Hallmark pathways enriched in PTC",
  border_color = NA,
  fontsize = 12
)
ggsave("F:/scRNA/workpath/2.thyrocytes/HALLMARKSERPINA1_gsva.png", plot = p, width = 10, height = 10, dpi = 300)
ggsave("F:/scRNA/workpath/2.thyrocytes/HALLMARKSERPINA1_gsva.pdf", plot = p, width = 10, height = 10, dpi = 300)
