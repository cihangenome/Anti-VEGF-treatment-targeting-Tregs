#setwd("/Users/oguzc/Downloads/Greten_10X_data/Greten_followup_recluster_pseudo_sep19")


setwd("/Users/oguzc/Downloads/Greten_10X_data/Greten_followup_visualize_recluster_pseudo_nov15")

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

library(CellChat)

Greten_10X_sc_objects.integrated<-readRDS("/Users/oguzc/Downloads/Greten_10X_data/Greten_10X_sc_objects.integrated_14samples_rpca_fewrefs_apr30.rds")


##########
##########

Greten_10X_sc_objects.integrated@meta.data$Time<-"Time"


Greten_10X_sc_objects.integrated@meta.data$Time[grep("_P",Greten_10X_sc_objects.integrated@meta.data[["orig.ident"]])]<-"pre"


Greten_10X_sc_objects.integrated@meta.data$Time[grep("_3W",Greten_10X_sc_objects.integrated@meta.data[["orig.ident"]])]<-"post"


sample_names<-unique(Greten_10X_sc_objects.integrated@meta.data[["orig.ident"]])

library(readr)

Greten_azimuth_pred<-read_tsv("/Users/oguzc/Downloads/Greten_10X_data/Greten_few_refs_azimuth_cell_type_preds/azimuth_pred.tsv")

df_cell_names<-data.frame(d1=Greten_azimuth_pred$cell,d2=colnames(Greten_10X_sc_objects.integrated))

#table(df_cell_names$d1==df_cell_names$d2)
# TRUE
#53842

#head(colnames(Greten_10X_sc_objects.integrated))


####add the predicted cell type into metadata
predicted.celltype<-unlist(Greten_azimuth_pred[,c(2)])

Greten_10X_sc_objects.integrated <- AddMetaData(object = Greten_10X_sc_objects.integrated,metadata = predicted.celltype,col.name = 'predicted.celltype'
)



Greten_10X_sc_objects.integrated@meta.data$predicted.celltype<-gsub(" ","_",Greten_10X_sc_objects.integrated@meta.data$predicted.celltype)

Greten_10X_sc_objects.integrated@meta.data$predicted.celltype_Time<-paste0(Greten_10X_sc_objects.integrated@meta.data[["predicted.celltype"]],"_",Greten_10X_sc_objects.integrated@meta.data[["Time"]])

labels_new<-as.matrix(levels(as.factor(Greten_10X_sc_objects.integrated@meta.data$predicted.celltype_Time)))

Greten_10X_sc_objects.integrated@meta.data$predicted.celltype_Time<-factor(Greten_10X_sc_objects.integrated@meta.data$predicted.celltype_Time,levels=as.character(rev(mixedsort(unique(labels_new)))))

Greten_10X_sc_objects.integrated@meta.data$predicted.celltype_Time<-as.factor(Greten_10X_sc_objects.integrated@meta.data$predicted.celltype_Time)

levels(Greten_10X_sc_objects.integrated@meta.data$predicted.celltype_Time)


labels_new<-as.matrix(levels(as.factor(Greten_10X_sc_objects.integrated@meta.data$predicted.celltype_Time)))

Greten_10X_sc_objects.integrated@meta.data$predicted.celltype_Time<-factor(Greten_10X_sc_objects.integrated@meta.data$predicted.celltype_Time,levels=as.character(rev(mixedsort(unique(labels_new)))))


#Granulocytes, monocytes, macrophages, and dendritic cells (DCs) represent a subgroup of leukocytes, collectively called myeloid cells.


Idents(Greten_10X_sc_objects.integrated) <- "Time"
# Extract the pre and post subsets
# Greten_10X_sc_objects.integrated_pre <- subset(Greten_10X_sc_objects.integrated, idents = "pre")
# Greten_10X_sc_objects.integrated_post <- subset(Greten_10X_sc_objects.integrated, idents = "post")

df_predicted.celltype<-as.data.frame(table(Greten_10X_sc_objects.integrated@meta.data$predicted.celltype))

selec_predicted.celltype<-df_predicted.celltype$Var1[which(df_predicted.celltype$Freq>500)]


############
############
# Idents(Greten_10X_sc_objects.integrated_pre) <- "predicted.celltype"
#
# Greten_10X_sc_objects.integrated_pre <- subset(Greten_10X_sc_objects.integrated_pre, idents = selec_predicted.celltype)
#
# Idents(Greten_10X_sc_objects.integrated_post) <- "predicted.celltype"
#
# Greten_10X_sc_objects.integrated_post <- subset(Greten_10X_sc_objects.integrated_post , idents = selec_predicted.celltype)

############
############

#CD14_Mono:5777
#CD16_Mono:1203
#Treg:1356
#DC: too few

###################
###################
Idents(Greten_10X_sc_objects.integrated)<-"predicted.celltype"
###################
###################
#Greten_10X_sc_objects.integrated_Treg<-subset(Greten_10X_sc_objects.integrated,idents="Treg")

#saveRDS(Greten_10X_sc_objects.integrated_Treg,"Greten_10X_sc_objects.integrated_14samples_Treg_sep19_2022.rds")

Greten_10X_sc_objects.integrated_Treg<-readRDS("/Users/oguzc/Downloads/Greten_10X_data/Greten_followup_recluster_pseudo_sep19/Greten_10X_sc_objects.integrated_14samples_Treg_sep19_2022.rds")
###################
###################
# Greten_10X_sc_objects.integrated_Mono<-subset(Greten_10X_sc_objects.integrated,idents=c("CD14_Mono","CD16_Mono"))
#
# saveRDS(Greten_10X_sc_objects.integrated_Mono,"Greten_10X_sc_objects.integrated_14samples_Mono_sep19_2022.rds")

Greten_10X_sc_objects.integrated_Mono<-readRDS("/Users/oguzc/Downloads/Greten_10X_data/Greten_followup_recluster_pseudo_sep19/Greten_10X_sc_objects.integrated_14samples_Mono_sep19_2022.rds")
###################
###################
# Greten_10X_sc_objects.integrated_DC<-subset(Greten_10X_sc_objects.integrated,idents=c("cDC1","cDC2","pDC"))
#
# saveRDS(Greten_10X_sc_objects.integrated_DC,"Greten_10X_sc_objects.integrated_14samples_DC_sep19_2022.rds")
###################
###################

# plural noun: granulocytes
# a white blood cell with secretory granules in its cytoplasm, i.e. a neutrophil, basophil, or eosinophil.

# Greten_10X_sc_objects.integrated_DC<-subset(Greten_10X_sc_objects.integrated,idents=c("cDC1","cDC2","pDC"))
#
# saveRDS(Greten_10X_sc_objects.integrated_DC,"Greten_10X_sc_objects.integrated_14samples_DC_sep19_2022.rds")

Greten_10X_sc_objects.integrated_DC<-readRDS("/Users/oguzc/Downloads/Greten_10X_data/Greten_followup_recluster_pseudo_sep19/Greten_10X_sc_objects.integrated_14samples_DC_sep19_2022.rds")


# Greten_10X_sc_objects.integrated_Treg<-readRDS("Greten_10X_sc_objects.integrated_14samples_Treg_sep19_2022.rds")

DefaultAssay(Greten_10X_sc_objects.integrated_Treg)<-"RNA"

#DefaultAssay(Greten_10X_sc_objects.integrated_Treg)<-"integrated"

Greten_10X_sc_objects.integrated_Treg <- FindVariableFeatures(Greten_10X_sc_objects.integrated_Treg)

dim_count=20

Greten_10X_sc_objects.integrated_Treg <- ScaleData(Greten_10X_sc_objects.integrated_Treg, verbose = TRUE)
Greten_10X_sc_objects.integrated_Treg <- RunPCA(Greten_10X_sc_objects.integrated_Treg, verbose = FALSE, npcs = dim_count)
# UMAP and Clustering
Greten_10X_sc_objects.integrated_Treg <- RunUMAP(Greten_10X_sc_objects.integrated_Treg, reduction = "pca", dims = 1:dim_count)
Greten_10X_sc_objects.integrated_Treg<- FindNeighbors(Greten_10X_sc_objects.integrated_Treg, reduction = "pca", dims = 1:dim_count)


Greten_10X_sc_objects.integrated_Treg <- FindClusters(Greten_10X_sc_objects.integrated_Treg, verbose =  TRUE,graph.name="RNA_nn", resolution = 0.2)

Greten_10X_sc_objects.integrated_Treg <- FindClusters(Greten_10X_sc_objects.integrated_Treg, verbose = TRUE,graph.name="RNA_snn", resolution = 0.2)

Idents(Greten_10X_sc_objects.integrated_Treg)<-"RNA_snn_res.0.2"

#p_Treg_clusters<-DimPlot(Greten_10X_sc_objects.integrated_Treg, reduction = "umap", group.by="RNA_snn_res.0.2", label=TRUE, label.size = 5,pt.size=1,raster=FALSE)
#+ plot_annotation(title = 'LP CD4+ (clusters)',subtitle = '')


pdf(file = "umap_Greten_10X_sc_objects.integrated_clusters_Tregs_nov15_2022.pdf",width=9, height=5)
DimPlot(Greten_10X_sc_objects.integrated_Treg, reduction = "umap", group.by="RNA_snn_res.0.2", label=TRUE, label.size = 5,pt.size=1,raster=FALSE)+ plot_annotation(title = 'Treg clusters after reclustering',subtitle = '')
dev.off()

pdf(file = "umap_Greten_10X_sc_objects.integrated_clusters_split_Time_Tregs_nov15_2022.pdf",width=9, height=5)
DimPlot(Greten_10X_sc_objects.integrated_Treg, reduction = "umap", split.by="predicted.celltype_Time", label=FALSE, label.size = 5, pt.size=1, ncol=2,raster=FALSE)+ plot_annotation(title = 'Treg clusters after reclustering',subtitle = '')
dev.off()

###########
###########

##########
##########

# Greten_10X_sc_objects.integrated_Mono<-readRDS("Greten_10X_sc_objects.integrated_14samples_Mono_sep19_2022.rds")

DefaultAssay(Greten_10X_sc_objects.integrated_Mono)<-"RNA"

#DefaultAssay(Greten_10X_sc_objects.integrated_Mono)<-"integrated"

Greten_10X_sc_objects.integrated_Mono <- FindVariableFeatures(Greten_10X_sc_objects.integrated_Mono)

dim_count=20

Greten_10X_sc_objects.integrated_Mono <- ScaleData(Greten_10X_sc_objects.integrated_Mono, verbose = TRUE)
Greten_10X_sc_objects.integrated_Mono <- RunPCA(Greten_10X_sc_objects.integrated_Mono, verbose = FALSE, npcs = dim_count)
# UMAP and Clustering
Greten_10X_sc_objects.integrated_Mono <- RunUMAP(Greten_10X_sc_objects.integrated_Mono, reduction = "pca", dims = 1:dim_count)
Greten_10X_sc_objects.integrated_Mono<- FindNeighbors(Greten_10X_sc_objects.integrated_Mono, reduction = "pca", dims = 1:dim_count)


Greten_10X_sc_objects.integrated_Mono <- FindClusters(Greten_10X_sc_objects.integrated_Mono, verbose =  TRUE,graph.name="RNA_nn", resolution = 0.2)

Greten_10X_sc_objects.integrated_Mono <- FindClusters(Greten_10X_sc_objects.integrated_Mono, verbose = TRUE,graph.name="RNA_snn", resolution = 0.2)

Idents(Greten_10X_sc_objects.integrated_Mono)<-"RNA_snn_res.0.2"



pdf(file = "umap_Greten_10X_sc_objects.integrated_clusters_Monos_nov15_2022.pdf",width=9, height=5)
DimPlot(Greten_10X_sc_objects.integrated_Mono, reduction = "umap", group.by="RNA_snn_res.0.2", label=TRUE, label.size = 5,pt.size=1,raster=FALSE)+ plot_annotation(title = 'Monocyte clusters after reclustering',subtitle = '')
dev.off()

pdf(file = "umap_Greten_10X_sc_objects.integrated_clusters_split_Time_Monos_nov15_2022.pdf",width=9, height=5)
DimPlot(Greten_10X_sc_objects.integrated_Mono, reduction = "umap", split.by="predicted.celltype_Time", label=FALSE, label.size = 5, pt.size=1, ncol=2,raster=FALSE)+ plot_annotation(title = 'Monocyte clusters after reclustering',subtitle = '')
dev.off()

##########
##########

##########
##########

# Greten_10X_sc_objects.integrated_DC<-readRDS("Greten_10X_sc_objects.integrated_14samples_DC_sep19_2022.rds")

DefaultAssay(Greten_10X_sc_objects.integrated_DC)<-"RNA"

#DefaultAssay(Greten_10X_sc_objects.integrated_DC)<-"integrated"

Greten_10X_sc_objects.integrated_DC <- FindVariableFeatures(Greten_10X_sc_objects.integrated_DC)

dim_count=20

Greten_10X_sc_objects.integrated_DC <- ScaleData(Greten_10X_sc_objects.integrated_DC, verbose = TRUE)
Greten_10X_sc_objects.integrated_DC <- RunPCA(Greten_10X_sc_objects.integrated_DC, verbose = FALSE, npcs = dim_count)
# UMAP and Clustering
Greten_10X_sc_objects.integrated_DC <- RunUMAP(Greten_10X_sc_objects.integrated_DC, reduction = "pca", dims = 1:dim_count)
Greten_10X_sc_objects.integrated_DC<- FindNeighbors(Greten_10X_sc_objects.integrated_DC, reduction = "pca", dims = 1:dim_count)


Greten_10X_sc_objects.integrated_DC <- FindClusters(Greten_10X_sc_objects.integrated_DC, verbose =  TRUE,graph.name="RNA_nn", resolution = 0.2)

Greten_10X_sc_objects.integrated_DC <- FindClusters(Greten_10X_sc_objects.integrated_DC, verbose = TRUE,graph.name="RNA_snn", resolution = 0.2)

Idents(Greten_10X_sc_objects.integrated_DC)<-"RNA_snn_res.0.2"



pdf(file = "umap_Greten_10X_sc_objects.integrated_clusters_DCs_nov15_2022.pdf",width=9, height=5)
DimPlot(Greten_10X_sc_objects.integrated_DC, reduction = "umap", group.by="RNA_snn_res.0.2", label=TRUE, label.size = 5,pt.size=1,raster=FALSE)+ plot_annotation(title = 'DC clusters after reclustering',subtitle = '')
dev.off()

pdf(file = "umap_Greten_10X_sc_objects.integrated_clusters_split_Time_DCs_nov15_2022.pdf",width=9, height=5)
DimPlot(Greten_10X_sc_objects.integrated_DC, reduction = "umap", split.by="predicted.celltype_Time", label=FALSE, label.size = 5, pt.size=1, ncol=2,raster=FALSE)+ plot_annotation(title = 'DC clusters after reclustering',subtitle = '')
dev.off()

####CIHAN HERE: NOV 15, 2022


##########
##########
###now report contribution of pre and post to each cluster (per myeloid subpopulation), find cluster-specific markers, perform pseudotime, find markers associated with pseudotime, and find the overlap between such markers and highly variable genes/cluster-specific markers.
##########
##########


##############
##############
df_Treg_metadata<-Greten_10X_sc_objects.integrated_Treg@meta.data
##############
##############

df_Treg_metadata_proportions_clusters_Time<- df_Treg_metadata %>%
  group_by(RNA_snn_res.0.2,Time) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  mutate(proportion_overall = count/sum(count))

df_Treg_metadata_proportions_clusters_Time<-df_Treg_metadata_proportions_clusters_Time[order(-df_Treg_metadata_proportions_clusters_Time$proportion_overall),]



df2_metadata_proportions_clusters_Time<- df_Treg_metadata %>%
  group_by(RNA_snn_res.0.2,Time) %>%
  summarize(count = n()) %>%
  mutate(proportion_in_cluster = count/sum(count))

df2_metadata_proportions_clusters_Time<-df2_metadata_proportions_clusters_Time[order(-df2_metadata_proportions_clusters_Time$proportion_in_cluster),]

##############
##############

df_Treg_proportions_clusters_Time<-merge(df2_metadata_proportions_clusters_Time,df_Treg_metadata_proportions_clusters_Time,by=c("RNA_snn_res.0.2","Time","count"))

df_Treg_proportions_clusters_Time<-df_Treg_proportions_clusters_Time[order(-df_Treg_proportions_clusters_Time$proportion_overall),]
###############
###############


##############
##############
df_DC_metadata<-Greten_10X_sc_objects.integrated_DC@meta.data
##############
##############

df_DC_metadata_proportions_clusters_Time<- df_DC_metadata %>%
  group_by(RNA_snn_res.0.2,Time) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  mutate(proportion_overall = count/sum(count))

df_DC_metadata_proportions_clusters_Time<-df_DC_metadata_proportions_clusters_Time[order(-df_DC_metadata_proportions_clusters_Time$proportion_overall),]



df2_metadata_proportions_clusters_Time<- df_DC_metadata %>%
  group_by(RNA_snn_res.0.2,Time) %>%
  summarize(count = n()) %>%
  mutate(proportion_in_cluster = count/sum(count))

df2_metadata_proportions_clusters_Time<-df2_metadata_proportions_clusters_Time[order(-df2_metadata_proportions_clusters_Time$proportion_in_cluster),]

##############
##############

df_DC_proportions_clusters_Time<-merge(df2_metadata_proportions_clusters_Time,df_DC_metadata_proportions_clusters_Time,by=c("RNA_snn_res.0.2","Time","count"))

df_DC_proportions_clusters_Time<-df_DC_proportions_clusters_Time[order(-df_DC_proportions_clusters_Time$proportion_overall),]
###############
###############

##############
##############
df_Mono_metadata<-Greten_10X_sc_objects.integrated_Mono@meta.data
##############
##############

df_Mono_metadata_proportions_clusters_Time<- df_Mono_metadata %>%
  group_by(RNA_snn_res.0.2,Time) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  mutate(proportion_overall = count/sum(count))

df_Mono_metadata_proportions_clusters_Time<-df_Mono_metadata_proportions_clusters_Time[order(-df_Mono_metadata_proportions_clusters_Time$proportion_overall),]



df2_metadata_proportions_clusters_Time<- df_Mono_metadata %>%
  group_by(RNA_snn_res.0.2,Time) %>%
  summarize(count = n()) %>%
  mutate(proportion_in_cluster = count/sum(count))

df2_metadata_proportions_clusters_Time<-df2_metadata_proportions_clusters_Time[order(-df2_metadata_proportions_clusters_Time$proportion_in_cluster),]

##############
##############

df_Mono_proportions_clusters_Time<-merge(df2_metadata_proportions_clusters_Time,df_Mono_metadata_proportions_clusters_Time,by=c("RNA_snn_res.0.2","Time","count"))

df_Mono_proportions_clusters_Time<-df_Mono_proportions_clusters_Time[order(-df_Mono_proportions_clusters_Time$proportion_overall),]
###############
###############

#slingshot(data,clusterLabels,reducedDim = NULL,start.clus = NULL,end.clus = NULL)

#_Treg start with cluster 0 (cluster_Time: 0_pre)
#_Mono start with cluster 0
#_DC start with cluster 0


#unique(Greten_10X_sc_objects.integrated_Mono@meta.data[["Time"]])
#[1] "pre"  "post"


library(slingshot)


sds_Treg_new<-readRDS("/Users/oguzc/Downloads/Greten_10X_data/Greten_followup_recluster_pseudo_sep19/sds_Treg_new_Greten_10X_sc_objects.integrated_sep19_2022.rds")

sds_DC_new<-readRDS("/Users/oguzc/Downloads/Greten_10X_data/Greten_followup_recluster_pseudo_sep19/sds_DC_new_Greten_10X_sc_objects.integrated_sep19_2022.rds")

sds_Mono_new<-readRDS("/Users/oguzc/Downloads/Greten_10X_data/Greten_followup_recluster_pseudo_sep19/sds_Mono_new_Greten_10X_sc_objects.integrated_sep19_2022.rds")

lnes_Treg<-readRDS("/Users/oguzc/Downloads/Greten_10X_data/Greten_followup_recluster_pseudo_sep19/sling_lineages_Treg_new_Greten_10X_sc_objects.integrated_sep19_2022.rds")

lnes_DC<-readRDS("/Users/oguzc/Downloads/Greten_10X_data/Greten_followup_recluster_pseudo_sep19/sling_lineages_DC_new_Greten_10X_sc_objects.integrated_sep19_2022.rds")

lnes_Mono<-readRDS("/Users/oguzc/Downloads/Greten_10X_data/Greten_followup_recluster_pseudo_sep19/sling_lineages_Mono_new_Greten_10X_sc_objects.integrated_sep19_2022.rds")

library(viridis)
library(chroma)




library(SummarizedExperiment)
sds_Treg_new_pseudo<-slingPseudotime(sds_Treg_new)
sds_Treg_new_pseudo_df <- data.frame(sds_Treg_new_pseudo)


clusters_Treg_df<-data.frame(cluster=Greten_10X_sc_objects.integrated_Treg@meta.data$RNA_snn_res.0.2)
sds_Treg_new_pseudo_df<-cbind(sds_Treg_new_pseudo_df,clusters_Treg_df)


library(gam)

library(Polychrome)

uniq_clusters_Treg<-unique(as.character(Greten_10X_sc_objects.integrated_Treg@meta.data$RNA_snn_res.0.2))

genes_Treg_traj<-read.table("/Users/oguzc/Downloads/Greten_10X_data/Greten_followup_visualize_recluster_pseudo_nov15/genes_Treg_reda_traj_nov_2022.txt",header = T)


###################
###################

sds_Treg_new_pseudo_df$cellbarcode<-rownames(sds_Treg_new_pseudo_df)

###################
###################

ind_relev_gene<-which(rownames(Greten_10X_sc_objects.integrated_Treg@assays[["RNA"]]) %in% genes_Treg_traj$Gene)

# rownames(Greten_10X_sc_objects.integrated_Treg@assays[["RNA"]]@data)[ind_relev_gene]
#  [1] "LCK"    "GNLY"   "IKZF2"  "IL7R"   "PSMB8"  "IL10RA" "LYZ"    "CRIP1"
#  [9] "IL32"   "TGFB1"  "NKG7"   "IL2RA"  "FOXP3"

df_exp<-as.matrix(Greten_10X_sc_objects.integrated_Treg@assays[["RNA"]]@data[ind_relev_gene,])

df_exp<-as.data.frame(df_exp)

df_exp<-t(df_exp)

df_exp<-as.data.frame(df_exp)

df_exp$cellbarcode<-rownames(df_exp)

cellbarcode_check<-data.frame(d1=df_exp$cellbarcode,d2=sds_Treg_new_pseudo_df$cellbarcode)

#table(cellbarcode_check$d1==cellbarcode_check$d2)
#TRUE
#1356


df_Treg_metadata<-Greten_10X_sc_objects.integrated_Treg@meta.data
df_Treg_metadata$cellbarcode<-rownames(df_Treg_metadata)



sds_Treg_new_pseudo_df_expanded<-merge(sds_Treg_new_pseudo_df,df_exp,by="cellbarcode")

sds_Treg_new_pseudo_df_expanded<-merge(sds_Treg_new_pseudo_df_expanded,df_Treg_metadata[,c(7,12)],by="cellbarcode")

sds_Treg_new_pseudo_df_expanded<-sds_Treg_new_pseudo_df_expanded[-which(sds_Treg_new_pseudo_df_expanded$cluster=="3"),]


sds_Treg_new_pseudo_df_expanded_Lineage1<-sds_Treg_new_pseudo_df_expanded[-which(is.na(sds_Treg_new_pseudo_df_expanded$Lineage1)),]

sds_Treg_new_pseudo_df_expanded_Lineage1<-sds_Treg_new_pseudo_df_expanded_Lineage1[order(sds_Treg_new_pseudo_df_expanded_Lineage1$Lineage1),]

####Lineage1

#######
#######
sds_Treg_new_pseudo_df_expanded_Lineage2<-sds_Treg_new_pseudo_df_expanded[-which(is.na(sds_Treg_new_pseudo_df_expanded$Lineage2)),]

sds_Treg_new_pseudo_df_expanded_Lineage2<-sds_Treg_new_pseudo_df_expanded_Lineage2[order(sds_Treg_new_pseudo_df_expanded_Lineage2$Lineage2),]
########
#######
#/Users/oguzc/Documents/DESKTOP_PREVIOUS_MACBOOK/NCBR_67_68/NCBR_68_FOLLOWUP_2020/correl_plots_reference/plo_tec_many_panels_nocorrels.R

#/Users/oguzc/Documents/DESKTOP_PREVIOUS_MACBOOK/NCBR_67_68/NCBR_68_FOLLOWUP_2020/correl_plots_reference/plo_tec_many_panels.R

########
#######

library(ggpubr)

###################
###################
Idents(object = Greten_10X_sc_objects.integrated_Treg) <- "RNA_snn_res.0.2"
Treg_genes_aver<-AverageExpression(Greten_10X_sc_objects.integrated_Treg, slot="data")

Treg_genes_aver<-as.data.frame(Treg_genes_aver[["RNA"]])
Treg_genes_aver$Gene<-rownames(Treg_genes_aver)

Treg_genes_aver_selec<-Treg_genes_aver[which(Treg_genes_aver$Gene %in% genes_Treg_traj$Genes),]

Treg_genes_aver_selec<-Treg_genes_aver_selec[,c(5,1:4)]

Treg_genes_aver<-Treg_genes_aver[,c(5,1:4)]

##################
##################
sds_Treg_new_pseudo_df_expanded_Lineage1_c02<-sds_Treg_new_pseudo_df_expanded_Lineage1[which(sds_Treg_new_pseudo_df_expanded_Lineage1$cluster=="0"| sds_Treg_new_pseudo_df_expanded_Lineage1$cluster=="2"),]
##################
##################

#genes_Group1<-c("IL32","CRIP1")
#genes_Treg_traj$Gene[c(1:8,9,13)]
#c("IL32","CRIP1","IL7R","IKZF2")
#genes_Group1<-setdiff(genes_Treg_traj$Gene[c(1:8,9,13)],c("IL32","CRIP1","IL7R","IKZF2"))

int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5]
}


############
############

p1<-ggscatter(sds_Treg_new_pseudo_df_expanded_Lineage1_c02, x = "Lineage1", y = c("IL32","CRIP1","IL10RA","IL7R","IKZF2"),legend.title="",ylab = "Expression level",xlab = "Pseudotime (Lineage1)",ylim = c(0, 4), show.legend = F, rug = FALSE, add = "loess",merge=TRUE,point=FALSE,title = "Treg clusters 0 & 2")+scale_y_continuous(breaks = function(x) int_breaks(x, n = 2))

p2<-ggscatter(sds_Treg_new_pseudo_df_expanded_Lineage1_c02, x = "Lineage1", y = c("LCK","TGFB1","PSMB8","FOXP3","IL2RA"),legend.title="",ylab = "Expression level",xlab = "Pseudotime (Lineage1)",ylim = c(0, 2), show.legend = F, rug = FALSE, add = "loess",merge=TRUE,point=FALSE,title = "Treg clusters 0 & 2")+scale_y_continuous(breaks = function(x) int_breaks(x, n = 2))

#setdiff(genes_Treg_traj$Gene[c(1:8,9,13)],c("IL32","CRIP1","IL7R","IKZF2","IL10RA"))

p3<-ggscatter(sds_Treg_new_pseudo_df_expanded_Lineage1_c02, x = "Lineage1", y = c("GNLY","LYZ","NKG7"),legend.title="",ylab = "Expression level",xlab = "Pseudotime (Lineage1)",ylim = c(0, 0.6), show.legend = F, rug = FALSE, add = "loess",merge=TRUE,point=FALSE,title = "Treg clusters 0 & 2")
#+scale_y_continuous(breaks = function(x) int_breaks(x, n = 2))
#setdiff(genes_Treg_traj$Gene[c(10,11,12)],c("IL32","CRIP1","IL7R","IKZF2","IL10RA"))

sds_Treg_new_pseudo_df_expanded_Lineage1_c02$Time<-as.factor(sds_Treg_new_pseudo_df_expanded_Lineage1_c02$Time)

sds_Treg_new_pseudo_df_expanded_Lineage1_c02$Time<-factor(sds_Treg_new_pseudo_df_expanded_Lineage1_c02$Time,levels=rev(mixedsort(levels(sds_Treg_new_pseudo_df_expanded_Lineage1_c02$Time))))

p4<-ggscatter(sds_Treg_new_pseudo_df_expanded_Lineage1_c02, x = "Lineage1", y = "cluster",legend.title="",ylab = "cluster",xlab = "Pseudotime (Lineage1)",ylim = c(0, 2), show.legend = F,point=TRUE,color="Time",palette="grey",title = "Treg clusters 0 & 2")  #,size=1

pdf(file = "scatter_plot_all_selec_genes_Treg_Lineage1_nov15_2022.pdf",width=9, height=6)
ggarrange(ggpar(p1,legend=c("right")),ggpar(p2,legend=c("right")),ggpar(p3,legend=c("right")),ggpar(p4,legend=c("right")),ncol = 2,nrow=2)
dev.off()

############
############

#####################
#####################
#####################

sds_Treg_new_pseudo_df_expanded_Lineage1_c01<-sds_Treg_new_pseudo_df_expanded_Lineage1[which(sds_Treg_new_pseudo_df_expanded_Lineage1$cluster=="0"| sds_Treg_new_pseudo_df_expanded_Lineage1$cluster=="1"),]

#####################
#####################
#####################

##################
##################
sds_Treg_new_pseudo_df_expanded_Lineage2_c01<-sds_Treg_new_pseudo_df_expanded_Lineage2[which(sds_Treg_new_pseudo_df_expanded_Lineage2$cluster=="0"| sds_Treg_new_pseudo_df_expanded_Lineage2$cluster=="1"),]
##################
##################

p1<-ggscatter(sds_Treg_new_pseudo_df_expanded_Lineage2_c01, x = "Lineage2", y = c("IL32","CRIP1","IL10RA","IL7R","IKZF2"),legend.title="",ylab = "Expression level",xlab = "Pseudotime (Lineage2)",ylim = c(0, 4), show.legend = F, rug = FALSE, add = "loess",merge=TRUE,point=FALSE,title = "Treg clusters 0 & 1")+scale_y_continuous(breaks = function(x) int_breaks(x, n = 2))

p2<-ggscatter(sds_Treg_new_pseudo_df_expanded_Lineage2_c01, x = "Lineage2", y = c("LCK","TGFB1","PSMB8","FOXP3","IL2RA"),legend.title="",ylab = "Expression level",xlab = "Pseudotime (Lineage2)",ylim = c(0, 2.5), show.legend = F, rug = FALSE, add = "loess",merge=TRUE,point=FALSE,title = "Treg clusters 0 & 1")+scale_y_continuous(breaks = function(x) int_breaks(x, n = 2))

p3<-ggscatter(sds_Treg_new_pseudo_df_expanded_Lineage2_c01, x = "Lineage2", y = c("GNLY","LYZ","NKG7"),legend.title="",ylab = "Expression level",xlab = "Pseudotime (Lineage2)",ylim = c(0, 3), show.legend = F, rug = FALSE, add = "loess",merge=TRUE,point=FALSE,title = "Treg clusters 0 & 1")


sds_Treg_new_pseudo_df_expanded_Lineage2_c01$Time<-as.factor(sds_Treg_new_pseudo_df_expanded_Lineage2_c01$Time)

sds_Treg_new_pseudo_df_expanded_Lineage2_c01$Time<-factor(sds_Treg_new_pseudo_df_expanded_Lineage2_c01$Time,levels=rev(mixedsort(levels(sds_Treg_new_pseudo_df_expanded_Lineage2_c01$Time))))

p4<-ggscatter(sds_Treg_new_pseudo_df_expanded_Lineage2_c01, x = "Lineage2", y = "cluster",legend.title="",ylab = "cluster",xlab = "Pseudotime (Lineage2)",ylim = c(0, 2), show.legend = F,point=TRUE,color="Time",palette="grey",title = "Treg clusters 0 & 1")

pdf(file = "scatter_plot_all_selec_genes_Treg_Lineage2_nov15_2022.pdf",width=9, height=6)
ggarrange(ggpar(p1,legend=c("right")),ggpar(p2,legend=c("right")),ggpar(p3,legend=c("right")),ggpar(p4,legend=c("right")),ncol = 2,nrow=2)
dev.off()

#10X_ICB_study_reclustering_trajectory_analysis_results_sep20_2022.pdf
#10X_ICB_cluster_specific_markers_after_reclustering_sep20_2022.xlsx

#write_xlsx(list(selected_genes=Treg_genes_aver_selec,all_genes=Treg_genes_aver), "10X_ICB_average_exp_by_cluster_Tregs_after_reclustering_nov15_2022.xlsx")

#Treg clusters after reclustering
#Monocyte clusters after reclustering
#DC clusters after reclustering

#/Users/oguzc/Documents/DESKTOP_PREVIOUS_MACBOOK/NCBR_67_68/NCBR_68_FOLLOWUP_2020


#df_predicted.celltype

##############
##############
#CD4_TCM
#CD8_TEM
#NK
#B_naive
#B_intermediate
#CD4_TEM
#dnT--
#Platelet--
#CD8_TCM
#CD8_Naive
#NK_CD56bright
#MAIT--
#CD4_Naive

#CD4_ ,CD8_, B_, NK_


#NK
#NK_CD56bright

#B_naive
#B_intermediate

#CD8_Naive
#CD8_TEM
#CD8_TCM

#CD4_Naive
#CD4_TEM
#CD4_TCM

#https://www.frontiersin.org/articles/10.3389/fimmu.2020.01014/full
#MAIT cell TCR provides an innate capacity to respond to a specific set of ligands without the need for expansion.
