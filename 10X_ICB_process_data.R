setwd("/Users/oguzc/Downloads/Greten_10X_data")

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

library(future)
#plan("multiprocess", workers = 40)
#options(future.globals.maxSize= 25091289600)
#options(future.globals.maxSize= 400000*1024^2)

#plan("multiprocess", workers = 4)
#options(future.globals.maxSize= 40000*1024^2)

#######
#######
files_greten_GEX<-read.table("/Users/oguzc/Downloads/Greten_10X_data/GEX_list_greten_h5_files.txt")
#######
#######
samples_greten_GEX<-read.table("/Users/oguzc/Downloads/Greten_10X_data/samples_list_greten_h5_files.txt")
#######
#######
sc_object_list<-vector(mode = "list", length = length(samples_greten_GEX$V1))



#/Users/oguzc/Downloads/NCBR-191/ncbr191_integrate_dec8.R
#, use.names = TRUE, unique.features = TRUE


total_obj_count<-length(samples_greten_GEX$V1)

for (objno in c(1:total_obj_count)) {

print(objno)

datadir<-paste0("/Users/oguzc/Downloads/Greten_10X_data/03_FilteredMatricesH5/",files_greten_GEX$V1[objno])

data<- Read10X_h5(filename =datadir, use.names = TRUE, unique.features = TRUE)

# Initialize the Seurat object with the raw (non-normalized data).
sc_object_list[[objno]]<-CreateSeuratObject(counts = data, project = samples_greten_GEX$V1[objno], min.cells = 100, min.features = 200)


sc_object_list[[objno]][["percent.mt"]] <- PercentageFeatureSet(sc_object_list[[objno]], pattern = "^MT-")

data_object_sc <- as.SingleCellExperiment(sc_object_list[[objno]])
data_object_sc <- scDblFinder(data_object_sc)

ind_doublets_data_object<-which(data_object_sc$scDblFinder.class=="doublet")
data_object_sc <- data_object_sc[,-ind_doublets_data_object]

high_mitoexp_data_object<-isOutlier(data_object_sc@colData@listData$percent.mt, nmads=4, type="higher")
high_ncount_data_object<-isOutlier(data_object_sc@colData@listData$nCount_RNA, nmads=2.5, type="higher")
low_ncount_data_object<-isOutlier(data_object_sc@colData@listData$nCount_RNA, nmads=3, type="lower")
high_nfeature_data_object<-isOutlier(data_object_sc@colData@listData$nFeature_RNA, nmads=2.5, type="higher")
low_nfeature_data_object<-isOutlier(data_object_sc@colData@listData$nFeature_RNA, nmads=3, type="lower")

ind_outliers_data_object<-Reduce("|", list(high_mitoexp_data_object,high_ncount_data_object,low_ncount_data_object,high_nfeature_data_object,low_nfeature_data_object))
table(ind_outliers_data_object)

length(which(high_mitoexp_data_object=="TRUE")) #
length(which(low_nfeature_data_object=="TRUE"))

ind_outliers_data_object<-which(ind_outliers_data_object=="TRUE")
#ind_outliers_data_object<-union(which(ind_outliers_data_object=="TRUE"),which(data_object_sc@colData@listData$nFeature_RNA<1000))

data_object_sc <- data_object_sc[,-ind_outliers_data_object]

length(data_object_sc@colData@rownames)

data_object_filtered <- subset(sc_object_list[[objno]], cells=data_object_sc@colData@rownames)

##normalize, variablefeatures,scale,runpca
##Normalize RNA data with log normalization
data_object_filtered <- NormalizeData(data_object_filtered)
# Find and scale variable features
data_object_filtered <- FindVariableFeatures(data_object_filtered,selection.method = "vst",mean.function = "FastExpMean")

#data_object_filtered <- ScaleData(data_object_filtered)
#data_object_filtered <- RunPCA(data_object_filtered, features = VariableFeatures(object =data_object_filtered))

sc_object_list[[objno]]<-data_object_filtered

    }


#saveRDS(sc_object_list, file = "Greten_10X_sc_object_list_14samples_apr29.rds")


features <- SelectIntegrationFeatures(sc_object_list)
sc_object_list <- lapply(X = sc_object_list, FUN = function(x) {
     x <- ScaleData(x, features = features, verbose = FALSE)
     x <- RunPCA(x, features = features, verbose = FALSE)
 })

saveRDS(sc_object_list, file = "Greten_10X_sc_object_list_14samples_apr29.rds")


sc_object_list<-readRDS("Greten_10X_sc_object_list_14samples_apr29.rds")

Greten_10X_sc_objects.anchors <- FindIntegrationAnchors(object.list = sc_object_list[c(1:total_obj_count)], dims = 1:20,verbose=TRUE, anchor.features = features, reduction = "rpca",k.anchor = 40, reference = c(5, 6))

saveRDS(Greten_10X_sc_objects.anchors,"Greten_10X_sc_objects.anchors.integrated_14samples_apr29.rds")

Greten_10X_sc_objects.anchors<-readRDS("Greten_10X_sc_objects.anchors.integrated_14samples_apr29.rds")


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
saveRDS(Greten_10X_sc_objects.integrated, file = "Greten_10X_sc_objects.integrated_14samples_apr29.rds")
#####################
###################
