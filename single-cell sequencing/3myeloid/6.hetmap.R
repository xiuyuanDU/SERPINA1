# Load necessary libraries
library(ggplot2)
library(dplyr)

# Set working directory
setwd('F:/scRNA/workpath/3.myeloid')

# Load myeloid subset data
Subset_myeloid = readRDS("Subset_myeloid.rds")

# Set identities based on cell types
Idents(Subset_myeloid) = Subset_myeloid@meta.data[["cellType"]]

# Subset data to retain specific myeloid cell types
Subset_myeloid <- subset(Subset_myeloid, idents = c("Ms-C0", "Ms-C1", "Ms-C3", "Ms-C4", "Ms-C5", "Ms-C6", 
                                                   "Ms-C7", "Ms-C9", "PMs-C11", "Ms-C12", "Ms-C14", "Ms-C15"))

# Set identities based on original sample identifiers
Idents(Subset_myeloid) <- Subset_myeloid@meta.data[["orig.ident"]]

# Convert sample names to factors and specify levels
Subset_myeloid@meta.data[["orig.ident"]] <- factor(Subset_myeloid@meta.data[["orig.ident"]], 
                                                       levels = c("T_WHT1", "T_WHT3", "T_WOHT5", "T_WOHT6"))

# Calculate frequency of each cell type per sample
celltype_frequency <- table(Subset_myeloid@meta.data[["orig.ident"]], Subset_myeloid$cellType)

# Specify selected cell types to keep
selected_cell_types <- c("Ms-C0", "Ms-C1", "Ms-C3", "Ms-C4", "Ms-C5", "Ms-C6", 
                         "Ms-C7", "Ms-C9", "PMs-C11", "Ms-C12", "Ms-C14", "Ms-C15")

# Filter frequency table to keep only selected cell types
celltype_frequency <- celltype_frequency[, selected_cell_types]

# Print modified frequency table
print(celltype_frequency)

# Save frequency table as CSV
write.csv(celltype_frequency, "10.celltypefrequency.csv")

# Convert frequency table to a data frame
celltype_frequency_df <- as.data.frame(celltype_frequency)
names(celltype_frequency_df) <- c("Sample", "Cell_Type", "Frequency")

# Calculate total cell counts for each sample
total_cell_counts <- aggregate(Frequency ~ Sample, data = celltype_frequency_df, FUN = sum)
names(total_cell_counts) <- c("Sample", "Total_Frequency")

# Merge frequency data with total counts
celltype_frequency_df <- merge(celltype_frequency_df, total_cell_counts, by = "Sample")

# Calculate relative frequency
celltype_frequency_df$Relative_Frequency <- celltype_frequency_df$Frequency / celltype_frequency_df$Total_Frequency

# Convert Sample column to an ordered factor
celltype_frequency_df$Sample <- factor(celltype_frequency_df$Sample, levels = c("T_WHT1", "T_WHT3", "T_WOHT5", "T_WOHT6"))

# Convert Cell_Type column to an ordered factor
celltype_frequency_df$Cell_Type <- factor(celltype_frequency_df$Cell_Type, 
                                           levels = c("Ms-C0", "Ms-C1", "Ms-C3", "Ms-C4", "Ms-C5", "Ms-C6", 
                                                      "Ms-C7", "Ms-C9", "PMs-C11", "Ms-C12", "Ms-C14", "Ms-C15"))

# Create heatmap using ggplot2
heatmap_plot <- ggplot(celltype_frequency_df, aes(x = Sample, y = Cell_Type, fill = Relative_Frequency)) +
  geom_tile(color = "black") +  # Add black borders around tiles
  scale_fill_gradient(low = "white", high = "steelblue") +  # Color gradient for fill
  theme_minimal() +
  labs(title = "Heatmap of Relative Cell Type Frequencies in Samples",
       x = "Sample",
       y = "Cell Type",
       fill = "Relative Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  geom_text(aes(label = round(Relative_Frequency, 2)), color = "black", size = 3)  # Add text labels

# Save the heatmap to a PDF file
ggsave("heatmap_relative_celltype_frequencies.pdf", plot = heatmap_plot, width = 7, height = 5, dpi = 300)
