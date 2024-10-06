library(ggplot2)
# Set working directory
setwd('F:/scRNA/workpath/0.allcell')
thyroid.combined <- readRDS("thyroid.combined.rds")

setwd('F:/scRNA/workpath/2.thyrocytes')
Subset_thyrocytes <- subset(thyroid.combined, idents = "thyrocytes")

# Remove cells expressing LYZ
gene_expression <- FetchData(Subset_thyrocytes, vars = c("LYZ"), cells = Cells(Subset_thyrocytes))
Subset_thyrocytes@meta.data[["LYZ"]] <- gene_expression
Subset_thyrocytes@meta.data[["LYZ"]] <- Subset_thyrocytes@meta.data[["LYZ"]][["LYZ"]]
Idents(Subset_thyrocytes) = Subset_thyrocytes@meta.data[["LYZ"]]
Subset_thyrocytes_filtered <- subset(Subset_thyrocytes, subset = LYZ == 0)
Subset_thyrocytes <- Subset_thyrocytes_filtered

# Remove cells expressing CD3D
gene_expression <- FetchData(Subset_thyrocytes, vars = c("CD3D"), cells = Cells(Subset_thyrocytes))
Subset_thyrocytes@meta.data[["CD3D"]] <- gene_expression
Subset_thyrocytes@meta.data[["CD3D"]] <- Subset_thyrocytes@meta.data[["CD3D"]][["CD3D"]]
Idents(Subset_thyrocytes) = Subset_thyrocytes@meta.data[["CD3D"]]
Subset_thyrocytes_filtered <- subset(Subset_thyrocytes, subset = CD3D == 0)
Subset_thyrocytes <- Subset_thyrocytes_filtered

# Extract gene expression levels of SERPINA1
gene_expression <- FetchData(Subset_thyrocytes, vars = c("SERPINA1"), cells = Cells(Subset_thyrocytes))
Subset_thyrocytes@meta.data[["SERPINA1"]] <- gene_expression
Subset_thyrocytes@meta.data[["SERPINA1"]] <- Subset_thyrocytes@meta.data[["SERPINA1"]][["SERPINA1"]]
Idents(Subset_thyrocytes) = Subset_thyrocytes@meta.data[["SERPINA1"]]

# Assume SERPINA1 is a numeric vector
SERPINA1 <- Subset_thyrocytes@meta.data[["SERPINA1"]]

# Check if SERPINA1 is 0 or missing (NA)
SERPINA1_expression <- SERPINA1 != 0 & !is.na(SERPINA1)
# Assume that Subset_thyrocytes@meta.data contains orig.ident identifying the samples
orig_ident <- Subset_thyrocytes@meta.data[["orig.ident"]]

# Calculate the number of cells expressing SERPINA1 in each sample
counts_by_sample <- table(orig_ident, SERPINA1_expression)
counts_by_sample
# Calculate the total number of cells in each sample
total_cells_by_sample <- as.vector(table(orig_ident))

# Calculate proportions
proportion_expression <- counts_by_sample[, "TRUE"] / total_cells_by_sample

# Print results
proportion_expression
samples_of_interest <- c("NT", "T_WHT3", "T_WOHT5", "T_WOHT6")

# Create a data frame to store this data
plot_data <- data.frame(
  Sample = samples_of_interest,
  Proportion = proportion_expression[samples_of_interest]
)

# Plot bar chart
p <- ggplot(plot_data, aes(x = Sample, y = Proportion)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.5) +
  labs(
    title = "Proportion of Cells Expressing SERPINA1",
    x = "Sample",
    y = "Proportion"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("serpina1_expression.pdf", plot = p, width = 4, height = 3, units = "in", dpi = 300)

Subset_thyrocytes_zero <- subset(Subset_thyrocytes_filtered, subset = SERPINA1 == 0)
Subset_thyrocytes_filtered <- subset(Subset_thyrocytes_filtered, subset = SERPINA1 > 0)
Subset_thyrocytes <- Subset_thyrocytes_filtered
save.image("1.SERPINA1_not_zero.RData")
#load("1.SERPINA1_not_zero.RData")
saveRDS(Subset_thyrocytes_zero, file = "Subset_thyrocytes_zero.rds")
