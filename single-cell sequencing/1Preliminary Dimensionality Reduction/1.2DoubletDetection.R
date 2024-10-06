#Load R packages
library(DoubletFinder)
library(tidyverse)
library(Seurat)
library(patchwork)
library(glmGamPoi)
#Set working directory
setwd('F:/scRNA/workpath')
dir.create("DoulbletFinder")
setwd("DoulbletFinder")

## NT
thyroid.single<-NT
thyroid.single<- SCTransform(thyroid.single)
thyroid.single <- RunPCA(thyroid.single, verbose = F)
ElbowPlot(thyroid.single)
pc.num=1:10
thyroid.single <- FindNeighbors(thyroid.single, dims = pc.num) %>% FindClusters(resolution = 0.3)
thyroid.single<- RunUMAP(thyroid.single, dims=pc.num)
## Find the optimal pK value
sweep.res.list <- paramSweep(thyroid.single, PCs = pc.num, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = ncol(NT)*8*1e-6 
homotypic.prop <- modelHomotypic(thyroid.single$seurat_clusters)  
nExp_poi <- round(DoubletRate*ncol(thyroid.single)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Use the determined parameters to identify doublets
thyroid.single<- doubletFinder(thyroid.single, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                               nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
## Results display, classification results are in pbmc@meta.data
DimPlot(thyroid.single, reduction = "umap", group.by = "DF.classifications_0.25_0.11_75")
## Create a subset excluding "Doublet"
NT<- subset(thyroid.single, subset = DF.classifications_0.25_0.11_75 != "Doublet")
# Plot a UMAP using the subsetted data
DimPlot(NT, reduction = "umap", group.by = "DF.classifications_0.25_0.11_75")
saveRDS(NT, file = "NT.rds")


## pT_WHT1
thyroid.single<-pT_WHT1
thyroid.single<- SCTransform(thyroid.single)
thyroid.single <- RunPCA(thyroid.single, verbose = F)
ElbowPlot(thyroid.single)
pc.num=1:10
thyroid.single <- FindNeighbors(thyroid.single, dims = pc.num) %>% FindClusters(resolution = 0.3)
thyroid.single<- RunUMAP(thyroid.single, dims=pc.num)
## Find the optimal pK value
sweep.res.list <- paramSweep(thyroid.single, PCs = pc.num, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = ncol(pT_WHT1)*8*1e-6 
homotypic.prop <- modelHomotypic(thyroid.single$seurat_clusters)   # It is best to provide cell type
nExp_poi <- round(DoubletRate*ncol(thyroid.single)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Use the determined parameters to identify doublets
thyroid.single<- doubletFinder(thyroid.single, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                               nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
## Results display, classification results are in pbmc@meta.data
DimPlot(thyroid.single, reduction = "umap", group.by = "DF.classifications_0.25_0.14_541")
# Create a subset excluding "Doublet"
pT_WHT1<- subset(thyroid.single, subset = DF.classifications_0.25_0.14_541 != "Doublet")
# Plot a UMAP using the subsetted data
DimPlot(pT_WHT1, reduction = "umap", group.by = "DF.classifications_0.25_0.14_541")
saveRDS(pT_WHT1, file = "pT_WHT1.rds")


## pT_WOHT4
thyroid.single<-pT_WOHT4
thyroid.single<- SCTransform(thyroid.single)
thyroid.single <- RunPCA(thyroid.single, verbose = F)
ElbowPlot(thyroid.single)
pc.num=1:10
thyroid.single <- FindNeighbors(thyroid.single, dims = pc.num) %>% FindClusters(resolution = 0.3)
thyroid.single<- RunUMAP(thyroid.single, dims=pc.num)
## Find the optimal pK value
sweep.res.list <- paramSweep(thyroid.single, PCs = pc.num, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = ncol(pT_WOHT4)*8*1e-6
homotypic.prop <- modelHomotypic(thyroid.single$seurat_clusters)   # It is best to provide cell type
nExp_poi <- round(DoubletRate*ncol(thyroid.single)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Use the determined parameters to identify doublets
thyroid.single<- doubletFinder(thyroid.single, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                               nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
## Results display, classification results are in pbmc@meta.data
DimPlot(thyroid.single, reduction = "umap", group.by = "DF.classifications_0.25_0.09_582")
# Create a subset excluding "Doublet"
pT_WOHT4<- subset(thyroid.single, subset = DF.classifications_0.25_0.09_582 != "Doublet")
# Plot a UMAP using the subsetted data
DimPlot(pT_WOHT4, reduction = "umap", group.by = "DF.classifications_0.25_0.09_582")
saveRDS(pT_WOHT4, file = "pT_WOHT4.rds")


## T_WHT1
thyroid.single<-T_WHT1
thyroid.single<- SCTransform(thyroid.single)
thyroid.single <- RunPCA(thyroid.single, verbose = F)
ElbowPlot(thyroid.single)
pc.num=1:10
thyroid.single <- FindNeighbors(thyroid.single, dims = pc.num) %>% FindClusters(resolution = 0.3)
thyroid.single<- RunUMAP(thyroid.single, dims=pc.num)
## Find the optimal pK value
sweep.res.list <- paramSweep(thyroid.single, PCs = pc.num, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = ncol(T_WHT1)*8*1e-6
homotypic.prop <- modelHomotypic(thyroid.single$seurat_clusters)   # It is best to provide cell type
nExp_poi <- round(DoubletRate*ncol(thyroid.single)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Use the determined parameters to identify doublets
thyroid.single<- doubletFinder(thyroid.single, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                               nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
## Results display, classification results are in pbmc@meta.data
DimPlot(thyroid.single, reduction = "umap", group.by = "DF.classifications_0.25_0.14_234")
# Create a subset excluding "Doublet"
T_WHT1<- subset(thyroid.single, subset = DF.classifications_0.25_0.14_234 != "Doublet")
# Plot a UMAP using the subsetted data
DimPlot(T_WHT1, reduction = "umap", group.by = "DF.classifications_0.25_0.14_234")
saveRDS(T_WHT1, file = "T_WHT1.rds")

## T_WHT2
thyroid.single <- T_WHT2
thyroid.single <- SCTransform(thyroid.single)
thyroid.single <- RunPCA(thyroid.single, verbose = F)
ElbowPlot(thyroid.single)
pc.num = 1:9
thyroid.single <- FindNeighbors(thyroid.single, dims = pc.num) %>% FindClusters(resolution = 0.3)
thyroid.single <- RunUMAP(thyroid.single, dims = pc.num)
## Find the optimal pK value
sweep.res.list <- paramSweep(thyroid.single, PCs = pc.num, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = ncol(T_WHT2) * 8 * 1e-6
homotypic.prop <- modelHomotypic(thyroid.single$seurat_clusters)   # It is best to provide cell type
nExp_poi <- round(DoubletRate * ncol(thyroid.single)) 
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
## Use the determined parameters to identify doublets
thyroid.single <- doubletFinder(thyroid.single, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                               nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
## Results display, classification results are in pbmc@meta.data
DimPlot(thyroid.single, reduction = "umap", group.by = "DF.classifications_0.25_0.05_174")
# Create a subset excluding "Doublet"
T_WHT2 <- subset(thyroid.single, subset = DF.classifications_0.25_0.05_174 != "Doublet")
# Plot a UMAP using the subsetted data
DimPlot(T_WHT2, reduction = "umap", group.by = "DF.classifications_0.25_0.05_174")
saveRDS(T_WHT2, file = "T_WHT2.rds")

## T_WHT3
thyroid.single <- T_WHT3
thyroid.single <- SCTransform(thyroid.single)
thyroid.single <- RunPCA(thyroid.single, verbose = F)
ElbowPlot(thyroid.single)
pc.num = 1:10
thyroid.single <- FindNeighbors(thyroid.single, dims = pc.num) %>% FindClusters(resolution = 0.3)
thyroid.single <- RunUMAP(thyroid.single, dims = pc.num)
## Find the optimal pK value
sweep.res.list <- paramSweep(thyroid.single, PCs = pc.num, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = ncol(T_WHT3) * 8 * 1e-6
homotypic.prop <- modelHomotypic(thyroid.single$seurat_clusters)   # It is best to provide cell type
nExp_poi <- round(DoubletRate * ncol(thyroid.single)) 
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
## Use the determined parameters to identify doublets
thyroid.single <- doubletFinder(thyroid.single, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                               nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
## Results display, classification results are in pbmc@meta.data
DimPlot(thyroid.single, reduction = "umap", group.by = "DF.classifications_0.25_0.17_604")
# Create a subset excluding "Doublet"
T_WHT3 <- subset(thyroid.single, subset = DF.classifications_0.25_0.17_604 != "Doublet")
# Plot a UMAP using the subsetted data
DimPlot(T_WHT3, reduction = "umap", group.by = "DF.classifications_0.25_0.17_604")
saveRDS(T_WHT3, file = "T_WHT3.rds")

## T_WOHT4
thyroid.single <- T_WOHT4
thyroid.single <- SCTransform(thyroid.single)
thyroid.single <- RunPCA(thyroid.single, verbose = F)
ElbowPlot(thyroid.single)
pc.num = 1:10
thyroid.single <- FindNeighbors(thyroid.single, dims = pc.num) %>% FindClusters(resolution = 0.3)
thyroid.single <- RunUMAP(thyroid.single, dims = pc.num)
## Find the optimal pK value
sweep.res.list <- paramSweep(thyroid.single, PCs = pc.num, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = ncol(T_WOHT4) * 8 * 1e-6
homotypic.prop <- modelHomotypic(thyroid.single$seurat_clusters)   # It is best to provide cell type
nExp_poi <- round(DoubletRate * ncol(thyroid.single)) 
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
## Use the determined parameters to identify doublets
thyroid.single <- doubletFinder(thyroid.single, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                               nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
## Results display, classification results are in pbmc@meta.data
DimPlot(thyroid.single, reduction = "umap", group.by = "DF.classifications_0.25_0.08_354")
# Create a subset excluding "Doublet"
T_WOHT4 <- subset(thyroid.single, subset = DF.classifications_0.25_0.08_354 != "Doublet")
# Plot a UMAP using the subsetted data
DimPlot(T_WOHT4, reduction = "umap", group.by = "DF.classifications_0.25_0.08_354")
saveRDS(T_WOHT4, file = "T_WOHT4.rds")

## T_WOHT5
thyroid.single <- T_WOHT5
thyroid.single <- SCTransform(thyroid.single)
thyroid.single <- RunPCA(thyroid.single, verbose = F)
ElbowPlot(thyroid.single)
pc.num = 1:10
thyroid.single <- FindNeighbors(thyroid.single, dims = pc.num) %>% FindClusters(resolution = 0.3)
thyroid.single <- RunUMAP(thyroid.single, dims = pc.num)
## Find the optimal pK value
sweep.res.list <- paramSweep(thyroid.single, PCs = pc.num, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = ncol(T_WOHT5) * 8 * 1e-6
homotypic.prop <- modelHomotypic(thyroid.single$seurat_clusters)   # It is best to provide cell type
nExp_poi <- round(DoubletRate * ncol(thyroid.single)) 
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
## Use the determined parameters to identify doublets
thyroid.single <- doubletFinder(thyroid.single, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                               nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
## Results display, classification results are in pbmc@meta.data
DimPlot(thyroid.single, reduction = "umap", group.by = "DF.classifications_0.25_0.29_638")
# Create a subset excluding "Doublet"
T_WOHT5 <- subset(thyroid.single, subset = DF.classifications_0.25_0.29_638 != "Doublet")
# Plot a UMAP using the subsetted data
DimPlot(T_WOHT5, reduction = "umap", group.by = "DF.classifications_0.25_0.29_638")
saveRDS(T_WOHT5, file = "T_WOHT5.rds")

## T_WOHT6
thyroid.single <- T_WOHT6
thyroid.single <- SCTransform(thyroid.single)
thyroid.single <- RunPCA(thyroid.single, verbose = F)
ElbowPlot(thyroid.single)
pc.num = 1:6
thyroid.single <- FindNeighbors(thyroid.single, dims = pc.num) %>% FindClusters(resolution = 0.3)
thyroid.single <- RunUMAP(thyroid.single, dims = pc.num)
## Find the optimal pK value
sweep.res.list <- paramSweep(thyroid.single, PCs = pc.num, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = ncol(T_WOHT6) * 8 * 1e-6
homotypic.prop <- modelHomotypic(thyroid.single$seurat_clusters)   # It is best to provide cell type
nExp_poi <- round(DoubletRate * ncol(thyroid.single)) 
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
## Use the determined parameters to identify doublets
thyroid.single <- doubletFinder(thyroid.single, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                               nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
## Results display, classification results are in pbmc@meta.data
DimPlot(thyroid.single, reduction = "umap", group.by = "DF.classifications_0.25_0.02_1082")
# Create a subset excluding "Doublet"
T_WOHT6 <- subset(thyroid.single, subset = DF.classifications_0.25_0.02_1082 != "Doublet")
# Plot a UMAP using the subsetted data
DimPlot(T_WOHT6, reduction = "umap", group.by = "DF.classifications_0.25_0.02_1082")
saveRDS(T_WOHT6, file = "T_WOHT6.rds")


