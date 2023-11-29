setwd("/data/NCBR/rawdata/Greten_10X_data")

library(Seurat)
library(gtools)
library(dplyr)
library(patchwork)
library(scran)

library(scDblFinder)

library(Matrix)
library(writexl)
library(readxl)

library(scater)
library(loomR)
#library(destiny)
library(ggplot2)
library(ggthemes)
#library(tidyverse)
library(SingleCellExperiment)
#library("ggbeeswarm", lib.loc="~/Library/R/3.5/library")
library(ggpubr)

##########
##########
library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize= 40000*1024^2)
##########
##########


sc_object_list<-readRDS("Greten_10X_sc_object_list_14samples_apr29.rds")

total_obj_count<-14

features <- SelectIntegrationFeatures(sc_object_list)
sc_object_list <- lapply(X = sc_object_list, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE)
  })

 Greten_10X_sc_objects.anchors <- FindIntegrationAnchors(object.list = sc_object_list[c(1:total_obj_count)], dims = 1:20,verbose=TRUE, anchor.features = features, reduction = "rpca",k.anchor = 40, reference = c(5, 6, 13, 14))

 saveRDS(Greten_10X_sc_objects.anchors,"Greten_10X_sc_objects.anchors.integrated_14samples_rpca_fewrefs_apr30.rds")

#Greten_10X_sc_objects.anchors<-readRDS("Greten_10X_sc_objects.anchors.integrated_14samples_rpca_fewrefs_apr30.rds")


Greten_10X_sc_objects.integrated <- IntegrateData(anchorset = Greten_10X_sc_objects.anchors, dims = 1:20)

DefaultAssay(Greten_10X_sc_objects.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
Greten_10X_sc_objects.integrated <- ScaleData(Greten_10X_sc_objects.integrated, verbose = TRUE)
Greten_10X_sc_objects.integrated <- RunPCA(Greten_10X_sc_objects.integrated, npcs = 20, verbose = FALSE)






#saveRDS(Greten_10X_sc_objects.integrated, file = #"Greten_10X_sc_objects.integrated_14samples_pre_umap_apr29.rd#s")

# UMAP and Clustering
Greten_10X_sc_objects.integrated <- RunUMAP(Greten_10X_sc_objects.integrated, reduction = "pca", dims = 1:20)
Greten_10X_sc_objects.integrated <- FindNeighbors(Greten_10X_sc_objects.integrated, reduction = "pca", dims = 1:20)
Greten_10X_sc_objects.integrated <- FindClusters(Greten_10X_sc_objects.integrated, resolution = 0.3)

#####################
saveRDS(Greten_10X_sc_objects.integrated, file = "Greten_10X_sc_objects.integrated_14samples_rpca_fewrefs_apr30.rds")
