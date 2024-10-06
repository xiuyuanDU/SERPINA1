setwd('F:/scRNA/workpath/0.allcell')
thyroid.combined = readRDS("thyroid.combined.rds")
setwd('F:/scRNA/workpath/4.lymphocytes')
Idents(thyroid.combined) <- thyroid.combined@meta.data[["cellType"]]
Subset_lymphocytes <- subset(thyroid.combined, idents = c("T＆NK cells", "B cells"))
Subset_lymphocytes <- ScaleData(Subset_lymphocytes)

# Find differentially expressed genes
Subset_lymphocytes <- FindVariableFeatures(Subset_lymphocytes, selection.method = "vst", nfeatures = 2000)

# Dimensionality reduction
Subset_lymphocytes <- RunPCA(Subset_lymphocytes, features = VariableFeatures(object = Subset_lymphocytes))
Subset_lymphocytes <- RunTSNE(Subset_lymphocytes, dims = 1:20)
Subset_lymphocytes <- RunHarmony(Subset_lymphocytes, group.by.vars = "stim")

pdf(file = "12.ElbowPlot_lymphocytes.pdf", width = 5, height = 4)
ElbowPlot(Subset_lymphocytes, ndims = 30)
dev.off()

Subset_lymphocytes <- JackStraw(Subset_lymphocytes, num.replicate = 100)
Subset_lymphocytes <- ScoreJackStraw(Subset_lymphocytes, dims = 1:20)
pdf(file = "12.jackstrawplot_lymphocytes.pdf", width = 7.5, height = 5.5)
JackStrawPlot(Subset_lymphocytes, dims = 1:20)
dev.off()

Subset_lymphocytesPC = 18
Subset_lymphocytes = FindNeighbors(Subset_lymphocytes, dims = 1:Subset_lymphocytesPC, reduction = "harmony")

for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.2, 1.5, 2, 2.5, 3)) {
  Subset_lymphocytes = FindClusters(Subset_lymphocytes, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}

apply(Subset_lymphocytes@meta.data[, grep("RNA_snn_res", colnames(Subset_lymphocytes@meta.data))], 2, table)

p2_tree = clustree(Subset_lymphocytes@meta.data, prefix = "RNA_snn_res.")
pdf(file = "12.clustertree_Subset_lymphocytesPC.pdf", width = 12, height = 10)
p2_tree
dev.off()

Subset_lymphocytes = FindNeighbors(Subset_lymphocytes, dims = 1:Subset_lymphocytesPC, reduction = "harmony")
Subset_lymphocytes <- FindClusters(Subset_lymphocytes, resolution = 1)
Subset_lymphocytes <- RunUMAP(Subset_lymphocytes, reduction = "harmony", dims = 1:Subset_lymphocytesPC)
Subset_lymphocytes <- RunTSNE(Subset_lymphocytes, reduction = "harmony", dims = 1:Subset_lymphocytesPC)

# Plot by cell clusters
pdf("12.lymphocytes_cluster_plot.pdf", width = 8, height = 6, pointsize = 0.5)
DimPlot(Subset_lymphocytes, reduction = "umap", label = TRUE, label.size = 3.5, pt.size = 0.5) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

# Plot by samples
pdf("12.lymphocytes_umap_plot.pdf", width = 8, height = 6, pointsize = 0.5)
palette("Dark2")
DimPlot(Subset_lymphocytes, reduction = "umap", group.by = "stim")
dev.off()
##################################################################################################
Idents(Subset_lymphocytes) <- Subset_lymphocytes@meta.data[["RNA_snn_res.1"]]
# only.pos: Retain only upregulated differentially expressed genes
Subset_lymphocytes.markers <- FindAllMarkers(Subset_lymphocytes, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Subset_lymphocytes.markers, file = "12.cluster_markers.csv")
head(Subset_lymphocytes.markers)
Subset_lymphocytes.markers <- read.csv("12.cluster_markers.csv")

# Top 5 genes
top5Subset_lymphocytes.markers <- Subset_lymphocytes.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
top5Subset_lymphocytes.markers[, c("gene", "cluster")]

# Color palette
library(viridis)
col <- viridis(50)
# Output the length of the color palette
length(col)

pdf(file = "12.lymphocytes-cluster.hetmap.pdf", width = 15, height = 10)
DoHeatmap(Subset_lymphocytes, 
          features = top5Subset_lymphocytes.markers$gene,
          slot = "data",
          group.by = "seurat_clusters",
          group.colors = col) +
  scale_fill_viridis(option = "plasma") +  
  scale_fill_gradient2(low = '#0099CC', mid = 'white', high = '#CC0033',
                       name = 'Z-score')
dev.off()

##########################################################################################
Idents(Subset_lymphocytes) <- Subset_lymphocytes@meta.data[["RNA_snn_res.1"]]
Subset_lymphocytes <- FindClusters(Subset_lymphocytes, resolution = 1)

# Defining gene lists for various cell types
genes <- list("T cells" = c("CD3E", "PTPRC"),  
              "NK" = c("NCAM1", "NCR1", "FCGR3A", "KLRC1", "KLRD1"),  
              "CD8+ T" = c("CD8A", "CD8B"),  # ,"CD69"
              "CD4+" = c("CD4", "CD40LG"),
              "cytotoxic CD8+ T" = c("PRF1", "GZMA", "NKG7", "ID2"),  # ,"PRDM1"
              "Treg" = c("IL2RA", "FOXP3", "CTLA4"),
              "T Memory Cells" = c("SELL", "CCR7", "LEF1"),
              "CD20+ B cells" = c("MS4A1"),
              "naive B cells" = c("IGHD", "FCER2", "TCL1A", "IL4R"),
              "Plasma cells" = c("SDC1"),
              "IgA+ Plasma" = c("IGHA1", "IGHA2"),
              "IgG+ Plasma" = c("IGHG1", "IGHG2", "IGHG3"))

pdf(file = "12.ann_cluster_marker_lymphocytes.pdf", width = 50, height = 15)
do_DotPlot(sample = Subset_lymphocytes, features = genes, group.by = "seurat_clusters", dot.scale = 12, colors.use = c("yellow", "red"), legend.length = 50,
           legend.framewidth = 2, font.size = 12)
dev.off()

ann.ids <- c("Tm-C1",
             "B-C1",
             "CD8-C1",
             "Tm-C2",
             "CD8-C2",
             "Treg",
             "Tm-C3",
             "B-C2",
             "CD8-C3",
             "NK",
             "Tm-C4",
             "Tm-C5",
             "Tm-C6",
             "CD8-C4",
             "Plasma cells",
             "B-C3",
             "Tm-C7",
             "CD8-C5",
             "B-C4")   

Subset_lymphocytesidens = mapvalues(Idents(Subset_lymphocytes), from = levels(Idents(Subset_lymphocytes)), to = ann.ids)
Idents(Subset_lymphocytes) = Subset_lymphocytesidens
Subset_lymphocytes$cellType = Idents(Subset_lymphocytes)

# Visualize UMAP/tSNE
pdf(file = "12-scRNA.UMAP.pdf", width = 7.5, height = 5.5)
DimPlot(Subset_lymphocytes, reduction = "umap", label = TRUE, label.size = 2.5, pt.size = 0.5) + 
  theme_classic() + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

pdf(file = "12-scRNA.TSEN.pdf", width = 7.5, height = 5.5)
DimPlot(Subset_lymphocytes, reduction = "tsne", label = TRUE, label.size = 3.5, pt.size = 0.5) + 
  theme_classic() + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), legend.position = "right")
dev.off()

Idents(Subset_lymphocytes) <- Subset_lymphocytes@active.ident
Subset_lymphocytes.markers <- FindAllMarkers(Subset_lymphocytes, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox", slot = "data") 
write.csv(Subset_lymphocytes.markers, file = "12.cell_markers.csv")

top5Subset_lymphocytes.markers <- Subset_lymphocytes.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Save as PDF
library(viridis)
palette <- viridis(n = 50)

# Plot
pdf(file = "12-cell_marker.hetmap.pdf", width = 15, height = 10)
DoHeatmap(
  Subset_lymphocytes,
  features = top5Subset_lymphocytes.markers$gene,
  slot = "data",
  group.colors = palette
) +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_gradient2(low = '#0099CC', mid = 'white', high = '#CC0033',
                       name = 'Z-score')
dev.off()

saveRDS(Subset_lymphocytes, file = "Subset_lymphocytes.rds")
save.image("my_workspacelymphocytes.RData")
setwd('F:/scRNA/workpath/4.lymphocytes')
load("my_workspacelymphocytes.RData")











