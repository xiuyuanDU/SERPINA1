## Relative cell frequency for each sample
Subset_myeloid@meta.data[["orig.ident"]] <- factor(Subset_myeloid@meta.data[["orig.ident"]], levels = c("NT", "pT_WHT1", "pT_WOHT4", "T_WHT1", "T_WHT2", "T_WHT3", "T_WOHT4", "T_WOHT5", "T_WOHT6"))
celltype_frequency <- table(Subset_myeloid@meta.data[["orig.ident"]], Subset_myeloid$cellType)
write.csv(celltype_frequency, "celltypefrequency.csv")

# Calculate the frequency of each cell type in each sample   
celltype_frequency_df <- as.data.frame(celltype_frequency)
names(celltype_frequency_df) <- c("Sample", "Cell_Type", "Frequency")
celltype_frequency_df <- aggregate(Frequency ~ Sample + Cell_Type, data = celltype_frequency_df, FUN = sum)

# Calculate total cell counts for each sample
total_cell_counts <- aggregate(Frequency ~ Sample, data = celltype_frequency_df, FUN = sum)

# Merge frequency data with total cell counts data
celltype_frequency_df <- merge(celltype_frequency_df, total_cell_counts, by = "Sample")

# Calculate relative frequency
celltype_frequency_df$Relative_Frequency <- celltype_frequency_df$Frequency.x / celltype_frequency_df$Frequency.y

# Provided total cell count data
total_cell_counts <- data.frame(Sample = c("NT", "pT_WHT1", "pT_WOHT4", "T_WHT1", "T_WHT2", "T_WHT3", "T_WOHT4", "T_WOHT5", "T_WOHT6"),
                                Frequency = c(37, 88, 134, 1734, 9, 2402, 326, 1456, 2482))

# Merge frequency data with total cell counts data
celltype_frequency_df <- merge(celltype_frequency_df, total_cell_counts, by = "Sample", suffixes = c(".x", ".y"))
names(celltype_frequency_df) <- c("Sample", "Cell_Type", "Frequency", "Total_Frequency")

# Calculate relative frequency
celltype_frequency_df$Relative_Frequency <- celltype_frequency_df$Frequency / celltype_frequency_df$Total_Frequency

# Convert Sample column to ordered factor and arrange in specified order
celltype_frequency_df$Sample <- factor(celltype_frequency_df$Sample, levels = c("NT", "pT_WHT1", "pT_WOHT4", "T_WHT1", "T_WHT2", "T_WHT3", "T_WOHT4", "T_WOHT5", "T_WOHT6"))

# Use ggplot2 to create a stacked bar plot and add total cell count labels above
barplot <- ggplot(celltype_frequency_df, aes(x = Sample, y = Relative_Frequency, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Relative Cell Type Frequency in Samples", x = "Sample", y = "Relative Frequency") +
  geom_text(aes(label = Total_Frequency, y = 1.05), vjust = 0) +  # Add labels above each bar
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the bar plot
print(barplot)

# Save the image to the specified path
ggsave("F:/scRNA/workpath/3.myeloid/10.barplot_myeloid.pdf", plot = barplot, width = 8, height = 6, dpi = 300)

##################################################################################################
## Relative cell frequency for each sample
Subset_myeloid@meta.data[["orig.ident"]] <- factor(Subset_myeloid@meta.data[["orig.ident"]], levels = c("T_WHT1", "T_WHT3", "T_WOHT5", "T_WOHT6"))
celltype_frequency <- table(Subset_myeloid@meta.data[["orig.ident"]], Subset_myeloid$cellType)
write.csv(celltype_frequency, "celltypefrequency2.csv")

# Calculate the frequency of each cell type in each sample   
celltype_frequency_df <- as.data.frame(celltype_frequency)
names(celltype_frequency_df) <- c("Sample", "Cell_Type", "Frequency")
celltype_frequency_df <- aggregate(Frequency ~ Sample + Cell_Type, data = celltype_frequency_df, FUN = sum)

# Calculate total cell counts for each sample
total_cell_counts <- aggregate(Frequency ~ Sample, data = celltype_frequency_df, FUN = sum)

# Merge frequency data with total cell counts data
celltype_frequency_df <- merge(celltype_frequency_df, total_cell_counts, by = "Sample")

# Calculate relative frequency
celltype_frequency_df$Relative_Frequency <- celltype_frequency_df$Frequency.x / celltype_frequency_df$Frequency.y

# Provided total cell count data
total_cell_counts <- data.frame(Sample = c("T_WHT1", "T_WHT3", "T_WOHT5", "T_WOHT6"),
                                Frequency = c(1734, 2402, 1456, 2482))

# Merge frequency data with total cell counts data
celltype_frequency_df <- merge(celltype_frequency_df, total_cell_counts, by = "Sample", suffixes = c(".x", ".y"))
names(celltype_frequency_df) <- c("Sample", "Cell_Type", "Frequency", "Total_Frequency")

# Calculate relative frequency
celltype_frequency_df$Relative_Frequency <- celltype_frequency_df$Frequency / celltype_frequency_df$Total_Frequency

# Convert Sample column to ordered factor and arrange in specified order
celltype_frequency_df$Sample <- factor(celltype_frequency_df$Sample, levels = c("T_WHT1", "T_WHT3", "T_WOHT5", "T_WOHT6"))


