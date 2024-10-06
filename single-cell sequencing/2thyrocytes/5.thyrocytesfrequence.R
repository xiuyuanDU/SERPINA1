# Load the subset of thyrocytes
Subset_thyrocytes = readRDS("Subset_thyrocytes.rds")

# Relative cell frequency for each sample
Subset_thyrocytes@meta.data[["orig.ident"]] <- factor(Subset_thyrocytes@meta.data[["orig.ident"]], levels = c("NT", "T_WHT3", "T_WOHT5", "T_WOHT6"))
celltype_frequency <- table(Subset_thyrocytes@meta.data[["orig.ident"]], Subset_thyrocytes$cellType)

# Write cell type frequency to CSV
write.csv(table(Subset_thyrocytes@meta.data[["orig.ident"]], Subset_thyrocytes$cellType), "celltypefrequency(thyrocytes).csv")

# Calculate the frequency of each cell type in each sample
celltype_frequency_df <- as.data.frame(celltype_frequency)
names(celltype_frequency_df) <- c("Sample", "Cell_Type", "Frequency")
celltype_frequency_df <- aggregate(Frequency ~ Sample + Cell_Type, data = celltype_frequency_df, FUN = sum)

# Calculate total cell counts for each sample
total_cell_counts <- aggregate(Frequency ~ Sample, data = celltype_frequency_df, FUN = sum)

# Merge frequency data with total cell count data
celltype_frequency_df <- merge(celltype_frequency_df, total_cell_counts, by = "Sample")

# Calculate relative frequency
celltype_frequency_df$Relative_Frequency <- celltype_frequency_df$Frequency.x / celltype_frequency_df$Frequency.y

# Provide total cell count data
total_cell_counts <- data.frame(Sample = c("NT", "T_WHT3", "T_WOHT5", "T_WOHT6"),
                                Frequency = c(295, 3597, 3541, 2539))

# Merge frequency data with total cell count data
celltype_frequency_df <- merge(celltype_frequency_df, total_cell_counts, by = "Sample", suffixes = c(".x", ".y"))
names(celltype_frequency_df) <- c("Sample", "Cell_Type", "Frequency", "Total_Frequency")

# Calculate relative frequency
celltype_frequency_df$Relative_Frequency <- celltype_frequency_df$Frequency / celltype_frequency_df$Total_Frequency

# Convert Sample column to ordered factor and arrange in specified order
celltype_frequency_df$Sample <- factor(celltype_frequency_df$Sample, levels = c("NT", "T_WHT3", "T_WOHT5", "T_WOHT6"))

# Create stacked bar plot using ggplot2
barplot <- ggplot(celltype_frequency_df, aes(x = Sample, y = Relative_Frequency, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Relative Cell Type Frequency in Samples", x = "Sample", y = "Relative Frequency") +
  geom_text(aes(label = Total_Frequency, y = 1.05), vjust = 0) +  # Add labels above each bar
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the bar plot
print(barplot)

# Save the plot to the specified path
ggsave("F:/scRNA/workpath/2.thyrocytes/barplot.pdf", plot = barplot, width = 5.8, height = 5.5, dpi = 300)

