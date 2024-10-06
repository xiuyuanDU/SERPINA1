library(dplyr)
library(patchwork)
library(Matrix)
library(stringr)
library(ggplot2)
library(Signac) 
library(tidyverse)
library(clusterProfiler)
library(psych)
library(qgraph)
library(igraph)
library(irGSEA)
library(plyr)
library(biomaRt)
library(enrichplot)
library(tibble)
# Set working directory
dir <- setwd('F:/TCGA/workpath/01serpina1expression')
library(openxlsx)

# Read sample names file
sample_data <- read.xlsx("sample.xlsx", sheet = 1)

# Read expression matrix file and keep original column names
expression_matrix <- read.table("HiSeqV2", header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

# Extract expression data for the SERPINA1 gene
SERPINA1_expression <- as.numeric(expression_matrix["SERPINA1", ])

# Ensure sample names are in column names
selected_samples <- sample_data$sampleID
matching_indices <- match(selected_samples, colnames(expression_matrix))

# Create a new column, initialized to NA
sample_data$SERPINA1 <- NA

# Find samples with matching expression and assign values
matched_samples <- !is.na(matching_indices)
sample_data$SERPINA1[matched_samples] <- SERPINA1_expression[matching_indices[matched_samples]]

# Identify samples without corresponding expression values
unmatched_samples <- sample_data$sampleID[!matched_samples]

# Display samples without corresponding SERPINA1 expression values
print("The following samples do not have corresponding SERPINA1 expression values:")
print(unmatched_samples)

# Save the results back to an Excel file
write.xlsx(sample_data, "sample_with_SERPINA1.xlsx")






















