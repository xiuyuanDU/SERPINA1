# Load required libraries
library(randomcoloR)
library(VennDiagram)

# Read in the differentially expressed genes (DEGs) data for moderately differentiated thyrocytes
Subset_DGEs <- read.csv("8.Subset_thyrocytes(Moderately).csv")
library(dplyr)

# Filter for upregulated DEGs with avg_log2FC > 1 and adjusted p-value < 0.05
UP01_DE_genes <- Subset_DGEs %>%
  filter(avg_log2FC > 1 & p_val_adj < 0.05) %>%
  dplyr::select(gene_symbol, avg_log2FC)

# Convert the data frame to a matrix
UP01_DE_genes_matrix <- as.matrix(UP01_DE_genes)

# Write the filtered DEGs to a CSV file
write.csv(UP01_DE_genes, "UP_Moderately_DE_genes.csv", row.names = FALSE)

#################################################################

# Read in the differentially expressed genes (DEGs) data for poorly differentiated thyrocytes
Subset_DGEs <- read.csv("8.Subset_thyrocytes(Poorly).csv")
library(dplyr)

# Filter for upregulated DEGs with avg_log2FC > 1 and adjusted p-value < 0.05
UP01_DE_genes <- Subset_DGEs %>%
  filter(avg_log2FC > 1 & p_val_adj < 0.05) %>%
  dplyr::select(gene_symbol, avg_log2FC)

# Convert the data frame to a matrix
UP01_DE_genes_matrix <- as.matrix(UP01_DE_genes)

# Write the filtered DEGs to a CSV file
write.csv(UP01_DE_genes, "UP_Poorly_DE_genes.csv", row.names = FALSE)

#################################################################

# Save in CSV format for one-way analysis, to be processed with Cytoscape 3.10 (64bit)