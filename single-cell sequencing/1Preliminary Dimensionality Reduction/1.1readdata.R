#Load R packages
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(stringr)
library(cowplot)
library(ggplot2)
library(Signac) 
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(psych)
library(qgraph)
library(igraph)
library(GSVA)
library(GSEABase)
library(limma)
library(hdf5r)
library(SingleR)
library(CCA)
library(clustree)
library(SCpubr)
library(UCell)
library(irGSEA)
library(harmony)
library(plyr)
library(randomcoloR)
library(CellChat)
library(monocle)
library(biomaRt)
library(enrichplot)
library(tibble)
# Set working directory
dir=setwd('F:/scRNA/data')
#Read in data
library(Seurat)
library(Matrix)
library(hdf5r)
NT.data <- Read10X_h5("NT.h5", use.names = T) 
pT_WHT1.data <- Read10X(data.dir = 'F:/scRNA/data/pT_WHT1') 
pT_WOHT4.data <- Read10X(data.dir = 'F:/scRNA/data/pT_WOHT4') 
T_WHT1.data <- Read10X(data.dir = 'F:/scRNA/data/T_WHT1') 
T_WHT2.data <- Read10X(data.dir = 'F:/scRNA/data/T_WHT2') 
T_WHT3.data <- Read10X(data.dir = 'F:/scRNA/data/T_WHT3') 
T_WOHT4.data <- Read10X(data.dir = 'F:/scRNA/data/T_WOHT4') 
T_WOHT5.data <- Read10X(data.dir = 'F:/scRNA/data/T_WOHT5')
T_WOHT6.data <- Read10X(data.dir = 'F:/scRNA/data/T_WOHT6')

#Create Seurat data while performing quality control
#NT
NT <- CreateSeuratObject(counts =NT.data, project = "NT", min.cells = 3, min.features = 200)
NT$stim <- "normal_thyroid"
NT[["percent.mt"]] <- PercentageFeatureSet(NT, pattern = "^MT-")
NT[["percent.rb"]] <- PercentageFeatureSet(NT, pattern = "^RP")
VlnPlot(NT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
NT <- subset(NT, subset = nCount_RNA >= 1000 & nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)

#pT_WHT1
pT_WHT1 <- CreateSeuratObject(counts = pT_WHT1.data, project = "pT_WHT1", min.cells = 3, min.features = 200)
pT_WHT1$stim <- "paratumor_tissue_With_HT1"
pT_WHT1[["percent.mt"]] <- PercentageFeatureSet(pT_WHT1, pattern = "^MT-")
pT_WHT1[["percent.rb"]] <- PercentageFeatureSet(pT_WHT1, pattern = "^RP")
VlnPlot(pT_WHT1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pT_WHT1 <- subset(pT_WHT1, subset = nCount_RNA >= 1000 & nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)

#pT_WOHT4
pT_WOHT4<- CreateSeuratObject(counts =pT_WOHT4.data, project = "pT_WOHT4", min.cells = 3, min.features = 200)
pT_WOHT4$stim <- "paratumor_tissue_Without_HT4"
pT_WOHT4[["percent.mt"]] <- PercentageFeatureSet(pT_WOHT4, pattern = "^MT-")
pT_WOHT4[["percent.rb"]] <- PercentageFeatureSet(pT_WOHT4, pattern = "^RP")
VlnPlot(pT_WOHT4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pT_WOHT4 <- subset(pT_WOHT4, subset = nCount_RNA >= 1000 & nFeature_RNA > 200 & nFeature_RNA < 7500& percent.mt < 10)

#T_WHT1
T_WHT1 <- CreateSeuratObject(counts = T_WHT1.data, project = "T_WHT1", min.cells = 3, min.features = 200)
T_WHT1$stim <- "tumor_With_HT1"
T_WHT1[["percent.mt"]] <- PercentageFeatureSet(T_WHT1, pattern = "^MT-")
T_WHT1[["percent.rb"]] <- PercentageFeatureSet(T_WHT1, pattern = "^RP")
VlnPlot(T_WHT1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T_WHT1 <- subset(T_WHT1, subset = nCount_RNA >= 1000 & nFeature_RNA > 200 & nFeature_RNA < 7500& percent.mt < 10)

#T_WHT2
T_WHT2 <- CreateSeuratObject(counts = T_WHT2.data, project = "T_WHT2", min.cells = 3, min.features = 200)
T_WHT2$stim <- "tumor_With_HT2"
T_WHT2[["percent.mt"]] <- PercentageFeatureSet(T_WHT2, pattern = "^MT-")
T_WHT2[["percent.rb"]] <- PercentageFeatureSet(T_WHT2, pattern = "^RP")
VlnPlot(T_WHT2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T_WHT2 <- subset(T_WHT2, subset = nCount_RNA >= 1000 & nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt <10)

#T_WHT3
T_WHT3 <- CreateSeuratObject(counts = T_WHT3.data, project = "T_WHT3", min.cells = 3, min.features = 200)
T_WHT3$stim <- "tumor_With_HT3"
T_WHT3[["percent.mt"]] <- PercentageFeatureSet(T_WHT3, pattern = "^MT-")
T_WHT3[["percent.rb"]] <- PercentageFeatureSet(T_WHT3, pattern = "^RP")
VlnPlot(T_WHT3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T_WHT3 <- subset(T_WHT3, subset = nCount_RNA >= 1000 & nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)


#T_WOHT4
T_WOHT4<- CreateSeuratObject(counts =T_WOHT4.data, project = "T_WOHT4", min.cells = 3, min.features = 200)
T_WOHT4$stim <- "tumor_Without_HT4"
T_WOHT4[["percent.mt"]] <- PercentageFeatureSet(T_WOHT4, pattern = "^MT-")
T_WOHT4[["percent.rb"]] <- PercentageFeatureSet(T_WOHT4, pattern = "^RP")
VlnPlot(T_WOHT4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T_WOHT4<- subset(T_WOHT4, subset = nCount_RNA >= 1000 & nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)

#T_WOHT5
T_WOHT5<- CreateSeuratObject(counts =T_WOHT5.data, project = "T_WOHT5", min.cells = 3, min.features = 200)
T_WOHT5$stim <- "tumor_Without_HT5"
T_WOHT5[["percent.mt"]] <- PercentageFeatureSet(T_WOHT5, pattern = "^MT-")
T_WOHT5[["percent.rb"]] <- PercentageFeatureSet(T_WOHT5, pattern = "^RP")
VlnPlot(T_WOHT5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T_WOHT5<- subset(T_WOHT5, subset = nCount_RNA >= 1000 & nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)

#T_WOHT6
T_WOHT6<- CreateSeuratObject(counts =T_WOHT6.data, project = "T_WOHT6", min.cells = 3, min.features = 200)
T_WOHT6$stim <- "tumor_Without_HT6"
T_WOHT6[["percent.mt"]] <- PercentageFeatureSet(T_WOHT6, pattern = "^MT-")
T_WOHT6[["percent.rb"]] <- PercentageFeatureSet(T_WOHT6, pattern = "^RP")
VlnPlot(T_WOHT6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T_WOHT6<- subset(T_WOHT6, subset = nCount_RNA >= 1000 & nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)

