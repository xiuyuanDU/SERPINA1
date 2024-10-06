# Set the working directory
setwd('F:/scRNA/workpath/2.thyrocytes')

# Load the subset of thyrocytes data
Subset_thyrocytes = readRDS("Subset_thyrocytes.rds")

# Set the identity of the subset based on the cell type metadata
Idents(Subset_thyrocytes) <- Subset_thyrocytes@meta.data[["cellType"]]

# Subset the data for poorly and highly expressed cell types
Subset_DGEs <- subset(Subset_thyrocytes, idents = c("Poorly", "Highly"))

# Find differentially expressed genes using Wilcoxon test
Subset_DGEs <- FindMarkers(Subset_thyrocytes, ident.1 = c("Poorly"), ident.2 = c("Highly"), test.use = 'wilcox', min.pct = 0.1)

# Write the results to a CSV file
write.csv(Subset_DGEs, file = "8.Subset_thyrocytes(Poorly).csv", row.names = TRUE)

# Manually add missing gene symbol column in the CSV format, then read the table again
Subset_DGEs <- read.csv("8.Subset_thyrocytes(Poorly).csv")

# Get the differentially expressed genes
DE_genes <- Subset_DGEs %>%
  filter(p_val_adj < 0.05) %>%
  pull(gene_symbol)

# Output the list of differentially expressed genes
print(DE_genes)

# Load the ggrepel library for better label placement
library(ggrepel)

# Prepare the data for plotting
cluster1.markers <- Subset_DGEs %>%
  mutate(Difference = pct.1 - pct.2) %>%
  rownames_to_column("gene")

# Define custom thresholds
log2FC = 1
padj = 0.05

# Categorize the markers based on thresholds
cluster1.markers$threshold = "ns"
cluster1.markers[which(cluster1.markers$avg_log2FC > log2FC & cluster1.markers$Difference >= 0.1 & cluster1.markers$p_val_adj < padj),]$threshold = "up"
cluster1.markers[which(cluster1.markers$avg_log2FC < (-log2FC) & cluster1.markers$Difference <= -0.1 & cluster1.markers$p_val_adj < padj),]$threshold = "down"
cluster1.markers$threshold = factor(cluster1.markers$threshold, levels = c('down', 'ns', 'up'))

# Create a PDF for the volcano plot
pdf("8.volcano_plot_thyrocytes(Poorly).pdf", width = 8, height = 6)

# Plot the volcano plot
ggplot(cluster1.markers, aes(x = Difference, y = avg_log2FC, color = threshold)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_label_repel(data = subset(cluster1.markers, avg_log2FC >= 1 & Difference >= 0.2 & p_val_adj <= 0.05),
                   aes(label = gene_symbol),  # Add labels
                   color = "black", # Set label text color
                   segment.colour = "black", # Set label segment color
                   label.padding = 0.1, 
                   segment.size = 0.3,  # Size of the label box
                   size = 4) +
  geom_label_repel(data = subset(cluster1.markers, avg_log2FC <= -1 & Difference <= -0.2 & p_val_adj <= 0.05),
                   aes(label = gene_symbol), label.padding = 0.1, 
                   color = "black",
                   segment.colour = "black",
                   segment.size = 0.3, size = 4) +
  geom_vline(xintercept = 0.0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Difference",
       y = "avg_log2FC", 
       title = "Volcano Plot of DE Genes (Poorly)"
  ) +
  theme_classic()

# Close the PDF device
dev.off()
##############################################################################################
# Set working directory
setwd('F:/scRNA/workpath/2.thyrocytes')

# Load the subset of thyrocytes data
Subset_thyrocytes = readRDS("Subset_thyrocytes.rds")

# Set the identity of the cells based on the 'cellType' metadata
Idents(Subset_thyrocytes) <- Subset_thyrocytes@meta.data[["cellType"]]

# Subset the data to include only "Moderately" and "Highly" differentiated cells
Subset_DGEs <- subset(Subset_thyrocytes, idents = c("Moderately", "Highly"))

# Perform differential expression analysis between the two groups
Subset_DGEs <- FindMarkers(Subset_thyrocytes, ident.1 = c("Moderately"), ident.2 = c("Highly"), test.use = 'wilcox', min.pct = 0.1)

# Write the results to a CSV file
write.csv(Subset_DGEs, file = "8.Subset_thyrocytes(Moderately).csv", row.names = TRUE)

# Manually add the missing gene name column "gene_symbol" in the CSV and read it back
Subset_DGEs <- read.csv("8.Subset_thyrocytes(Moderately).csv")

# Extract differentially expressed genes with adjusted p-value < 0.05
DE_genes <- Subset_DGEs %>%
  filter(p_val_adj < 0.05) %>%
  pull(gene_symbol)

# Print the list of differentially expressed genes
print(DE_genes)

# Load necessary library for plotting
library(ggrepel)

# Create a dataframe with the differences and genes
cluster1.markers <- Subset_DGEs %>%
  mutate(Difference = pct.1 - pct.2) %>% 
  rownames_to_column("gene")

# Set custom thresholds for log2 fold change and adjusted p-value
log2FC = 1
padj = 0.05 

# Assign significance thresholds
cluster1.markers$threshold = "ns"
cluster1.markers[which(cluster1.markers$avg_log2FC > log2FC & 
                        cluster1.markers$Difference >= 0.1 & 
                        cluster1.markers$p_val_adj < padj),]$threshold = "up"

cluster1.markers[which(cluster1.markers$avg_log2FC < (-log2FC) & 
                        cluster1.markers$Difference <= -0.1 & 
                        cluster1.markers$p_val_adj < padj),]$threshold = "down"

# Set threshold as a factor for plotting
cluster1.markers$threshold = factor(cluster1.markers$threshold, levels = c('down', 'ns', 'up'))

# Create a PDF for the volcano plot
pdf("8.volcano_plot_thyrocytes(Moderately).pdf", width = 8, height = 6) 

# Generate volcano plot
ggplot(cluster1.markers, aes(x = Difference, y = avg_log2FC, color = threshold)) + 
  geom_point(size = 0.5) + 
  scale_color_manual(values = c("blue", "grey", "red")) + 
  geom_label_repel(data = subset(cluster1.markers, avg_log2FC >= 1 & 
                                  Difference >= 0.2 & p_val_adj <= 0.05), 
                   aes(label = gene_symbol),  # Add gene labels
                   color = "black", 
                   segment.colour = "black",
                   label.padding = 0.1, 
                   segment.size = 0.3,  
                   size = 4) +
  geom_label_repel(data = subset(cluster1.markers, avg_log2FC <= -1 & 
                                  Difference <= -0.2 & p_val_adj <= 0.05), 
                   aes(label = gene_symbol), label.padding = 0.1, 
                   color = "black",
                   segment.colour = "black",
                   segment.size = 0.3, size = 4) +
  geom_vline(xintercept = 0.0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Difference",
       y = "avg_log2FC", 
       title = "Volcano Plot of DE Genes (Moderately)"
  )
theme_classic()
dev.off()
