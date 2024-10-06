# Set the working directory and load the processed data
setwd('F:/scRNA/workpath/3.myeloid')
Subset_myeloid <- readRDS("Subset_myeloid.rds")

# Set the cell type identities and subset for specific cell types
Idents(Subset_myeloid) <- Subset_myeloid@meta.data[["cellType"]]
Subset_DC <- subset(Subset_myeloid, idents = c("mDCs-C2", "mDCs-C8", "pDCs-C10", "LAMP_mDCs-C13"))

# Save the subset for specific cell types
saveRDS(Subset_DC, file = "Subset_DC.rds")
############################################################################

# Load required libraries
library(ggplot2)
library(dplyr)

# Read the processed data again for the following analysis
setwd('F:/scRNA/workpath/3.myeloid')
Subset_myeloid <- readRDS("Subset_myeloid.rds")
Idents(Subset_myeloid) <- Subset_myeloid@meta.data[["cellType"]]

# Subset for specific cell types
Subset_myeloid <- subset(Subset_myeloid, idents = c("mDCs-C2", "mDCs-C8", "pDCs-C10", "LAMP_mDCs-C13"))
Idents(Subset_myeloid) <- Subset_myeloid@meta.data[["orig.ident"]]

# Convert sample names to factor and set levels
Subset_myeloid@meta.data[["orig.ident"]] <- factor(Subset_myeloid@meta.data[["orig.ident"]], 
                                                   levels = c("T_WHT1", "T_WHT3", "T_WOHT5", "T_WOHT6"))

# Compute the frequency of each cell type in each sample
celltype_frequency <- table(Subset_myeloid@meta.data[["orig.ident"]], Subset_myeloid$cellType)

# Select specific cell types to retain
selected_cell_types <- c("mDCs-C2", "mDCs-C8", "pDCs-C10", "LAMP_mDCs-C13")
celltype_frequency <- celltype_frequency[, selected_cell_types]

# Print the modified data table and export to CSV
print(celltype_frequency)
write.csv(celltype_frequency, "10.CDScelltypefrequency.csv")

# Convert to a data frame
celltype_frequency_df <- as.data.frame(celltype_frequency)
names(celltype_frequency_df) <- c("Sample", "Cell_Type", "Frequency")

# Compute total cell counts for each sample
total_cell_counts <- aggregate(Frequency ~ Sample, data = celltype_frequency_df, FUN = sum)
names(total_cell_counts) <- c("Sample", "Total_Frequency")

# Check the total cell counts data frame
print(total_cell_counts)

# Merge frequency data with total cell counts
celltype_frequency_df <- merge(celltype_frequency_df, total_cell_counts, by = "Sample", all.x = TRUE)

# Check the merged data frame
print(celltype_frequency_df)

# Ensure Total_Frequency column has no missing values
if (any(is.na(celltype_frequency_df$Total_Frequency))) {
  stop("Total_Frequency column has missing values, unable to calculate relative frequency")
}

# Compute relative frequency
celltype_frequency_df$Relative_Frequency <- celltype_frequency_df$Frequency / celltype_frequency_df$Total_Frequency

# Convert Sample column to ordered factor with specified order
celltype_frequency_df$Sample <- factor(celltype_frequency_df$Sample, levels = c("T_WHT1", "T_WHT3", "T_WOHT5", "T_WOHT6"))

# Convert Cell_Type column to ordered factor with specified order
celltype_frequency_df$Cell_Type <- factor(celltype_frequency_df$Cell_Type, levels = c("mDCs-C2", "mDCs-C8", "pDCs-C10", "LAMP_mDCs-C13"))

# Create labels with relative frequency and cell count
celltype_frequency_df$Label <- paste0("Freq: ", round(celltype_frequency_df$Relative_Frequency, 2), "\nCount: ", celltype_frequency_df$Frequency)

# Create a heatmap using ggplot2
heatmap_plot <- ggplot(celltype_frequency_df, aes(x = Sample, y = Cell_Type, fill = Relative_Frequency)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(title = "Heatmap of Relative Cell Type Frequencies in Samples",
       x = "Sample",
       y = "Cell Type",
       fill = "Relative Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, vjust = 0.1, size = 14, face = "bold")) +
  geom_text(aes(label = Label), color = "black", size = 3, vjust = 0.1)

# Save the heatmap plot to a PDF file
ggsave("CDSheatmap_relative_celltype_frequencies.pdf", plot = heatmap_plot, width = 8.6, height = 6, dpi = 300)