# Extract the state information from the mycds object
state_info <- pData(mycds)$State

# Replace the state information
state_info <- recode(state_info, "1" = "state1", "2" = "state2", "3" = "state3")

# Add this state information to the meta.data of scRNAsub
scRNAsub@meta.data$State <- state_info

# Confirm successful addition
head(scRNAsub@meta.data)

# Create a subset for Macrophage for GSVA
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

# Read the gene set file and create a gene set list
h_df <- read.gmt("newgenesets.v2023.2.Hs.gmt")[, c(2, 1)]
h_list <- unstack(h_df)

# Perform GSVA gene set enrichment analysis
ES <- gsva(exp, h_list)
print(ES[1:3, 1:3])

# Save the results to a CSV file
write.csv(ES, "11.allDiffGO_Macrophage_GSVA.csv")

# Extract data for different states
State1_ES <- ES[,"state1"]
State2_ES <- ES[,"state2"]
State3_ES <- ES[,"state3"]

# Extract the top enriched pathways for each state
top_pathways_State1 <- names(sort(State1_ES, decreasing = TRUE)[1:15])
top_pathways_State2 <- names(sort(State2_ES, decreasing = TRUE)[1:15])
top_pathways_State3 <- names(sort(State3_ES, decreasing = TRUE)[1:15])

# Create an empty data frame to store the enrichment scores of the top pathways for State 1
top_pathways_scores_State1 <- data.frame(matrix(ncol = ncol(ES), nrow = 15))
rownames(top_pathways_scores_State1) <- top_pathways_State1
colnames(top_pathways_scores_State1) <- colnames(ES)

# Create an empty data frame to store the enrichment scores of the top pathways for State 2
top_pathways_scores_State2 <- data.frame(matrix(ncol = ncol(ES), nrow = 15))
rownames(top_pathways_scores_State2) <- top_pathways_State2
colnames(top_pathways_scores_State2) <- colnames(ES)

# Create an empty data frame to store the enrichment scores of the top pathways for State 3
top_pathways_scores_State3 <- data.frame(matrix(ncol = ncol(ES), nrow = 15))
rownames(top_pathways_scores_State3) <- top_pathways_State3
colnames(top_pathways_scores_State3) <- colnames(ES)

# Fill in the data frames with the respective enrichment scores
for (i in 1:ncol(ES)) {
  top_pathways_scores_State1[, i] <- ES[rownames(top_pathways_scores_State1), i]
}
for (i in 1:ncol(ES)) {
  top_pathways_scores_State2[, i] <- ES[rownames(top_pathways_scores_State2), i]
}
for (i in 1:ncol(ES)) {
  top_pathways_scores_State3[, i] <- ES[rownames(top_pathways_scores_State3), i]
}

# Assuming the necessary libraries are already installed
library(pheatmap)

# Extract the top pathways and their enrichment scores
top_pathways_State1 <- names(sort(State1_ES, decreasing = TRUE)[1:15])
top_pathways_State2 <- names(sort(State2_ES, decreasing = TRUE)[1:15])
top_pathways_State3 <- names(sort(State3_ES, decreasing = TRUE)[1:15])

# Extract the data for these pathways from the ES data frame
top_pathways_data <- ES[c(top_pathways_State1, top_pathways_State2, top_pathways_State3), ]

# Create a heatmap
pdf("11.top_pathways_heatmap.pdf", width = 15, height = 12)  # Set PDF file width and height
pheatmap(
  top_pathways_data,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  breaks = seq(-0.4, 0.4, length.out = 101),
  main = "Top Enriched Pathways",
  border_color = NA,  # Remove clustering lines
  fontsize = 12  # Set font size
)
dev.off()