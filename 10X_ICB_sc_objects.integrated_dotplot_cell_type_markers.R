setwd("/Users/oguzc/Downloads/Greten_10X_data/reda_dotplot_cell_type_markers_final_oct31_2023")

###############
###############

###find the expression upregulated among the CDR3s with >10 abundance

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

Idents(object = Greten_10X_sc_objects.integrated) <- "seurat_clusters"

Greten_10X_sc_objects.integrated@meta.data$Time<-"Time"

Greten_10X_sc_objects.integrated@meta.data$Time[grep("_P",Greten_10X_sc_objects.integrated@meta.data[["orig.ident"]])]<-"pre"

Greten_10X_sc_objects.integrated@meta.data$Time[grep("_3W",Greten_10X_sc_objects.integrated@meta.data[["orig.ident"]])]<-"post"

DefaultAssay(Greten_10X_sc_objects.integrated) <- "RNA"


Greten_azimuth_umap<-readRDS("/Users/oguzc/Downloads/Greten_10X_data/Greten_few_refs_azimuth_cell_type_preds/azimuth_umap.Rds")

library(readr)

Greten_azimuth_pred<-read_tsv("/Users/oguzc/Downloads/Greten_10X_data/Greten_few_refs_azimuth_cell_type_preds/azimuth_pred.tsv")

Greten_azimuth_impADT<-readRDS("/Users/oguzc/Downloads/Greten_10X_data/Greten_few_refs_azimuth_cell_type_preds/azimuth_impADT.Rds")

Greten_azimuth_human_pbmc_ref<-readRDS("/Users/oguzc/Downloads/Greten_10X_data/Greten_few_refs_azimuth_cell_type_preds/azimuth_human_pbmc_ref.Rds")

df_cell_names<-data.frame(d1=Greten_azimuth_pred$cell,d2=colnames(Greten_10X_sc_objects.integrated))

#table(df_cell_names$d1==df_cell_names$d2)
# TRUE
#53842

#head(colnames(Greten_10X_sc_objects.integrated))


####add the predicted cell type into metadata
predicted.celltype<-unlist(Greten_azimuth_pred[,c(2)])

Greten_10X_sc_objects.integrated <- AddMetaData(object = Greten_10X_sc_objects.integrated,metadata = predicted.celltype,col.name = 'predicted.celltype'
)


Greten_10X_sc_objects.integrated@reductions$umap2<-Greten_azimuth_umap

Idents(object = Greten_10X_sc_objects.integrated) <- "Time"

Idents(object = Greten_10X_sc_objects.integrated) <- "predicted.celltype"

df_cell_counts<-as.data.frame(table(Greten_10X_sc_objects.integrated@meta.data[["predicted.celltype"]]))

predicted.celltype_selec<-df_cell_counts$Var1[which(df_cell_counts$Freq>500)]


#cells_selec<-lapply(Greten_10X_sc_objects.integrated@meta.data$predicted.celltype, function(x) which(as.character(predicted.celltype_selec)==x))

#cells_selec<-unlist(cells_selec)

#cells_selec<-colnames(Greten_10X_sc_objects.integrated)[cells_selec]

Idents(Greten_10X_sc_objects.integrated) <- "predicted.celltype"

Greten_10X_sc_objects.integrated_selec <- subset(Greten_10X_sc_objects.integrated, idents =predicted.celltype_selec )



Idents(object = Greten_10X_sc_objects.integrated_selec) <- "Time"

df_metadata_integrated_14samples<-Greten_10X_sc_objects.integrated@meta.data

df_metadata_integrated_14samples$cell_sub_type<-predicted.celltype

ind_3W<-grep("_3W",df_metadata_integrated_14samples$orig.ident)
ind_P<-grep("_P",df_metadata_integrated_14samples$orig.ident)
df_metadata_integrated_14samples$Group<-"New"

df_metadata_integrated_14samples$Group[ind_3W]<-"post"
df_metadata_integrated_14samples$Group[ind_P]<-"pre"

table(df_metadata_integrated_14samples$Group)
#post   pre
#25380 28462

colnames(df_metadata_integrated_14samples)[1]<-"Sample"


df_metadata_integrated_14samples$Cell_name<-rownames(df_metadata_integrated_14samples)

df_metadata_integrated_14samples<-df_metadata_integrated_14samples[,c(ncol(df_metadata_integrated_14samples),1:ncol(df_metadata_integrated_14samples)-1)]

###############
###############

#unique(Greten_10X_sc_objects.integrated_selec@meta.data[["predicted.celltype"]])
#  [1] "NK" "CD8 TEM""CD14 Mono"  "B naive"
#  [5] "CD16 Mono"  "CD8 Naive"  "CD4 TCM""CD8 TCM"
#  [9] "NK_CD56bright"  "Treg"   "CD4 TEM""CD4 Naive"
# [13] "B intermediate" "Platelet"   "dnT""MAIT"

###############
###############

# unique(Greten_10X_sc_objects.integrated@meta.data[["predicted.celltype"]])
#  [1] "NK""CD8 TEM"   "CD14 Mono"
#  [4] "B naive"   "CD16 Mono" "CD8 Naive"
#  [7] "CD4 TCM"   "CD8 TCM"   "NK_CD56bright"
# [10] "Treg"  "CD8 Proliferating" "B memory"
# [13] "CD4 TEM"   "CD4 Proliferating" "NK Proliferating"
# [16] "CD4 CTL"   "pDC"   "CD4 Naive"
# [19] "gdT"   "B intermediate""ILC"
# [22] "Platelet"  "cDC1"  "cDC2"
# [25] "Plasmablast"   "dnT"   "HSPC"
# [28] "Eryth" "MAIT"

azimuth_human_pbmc_markers<-read_excel("/Users/oguzc/Downloads/Greten_10X_data/reda_dotplot_cell_type_markers_final_oct31_2023/azimuth_human_pbmc_markers_oct31_2023.xlsx",sheet="level2")


azimuth_human_pbmc_markers$markers<-gsub(" ","",azimuth_human_pbmc_markers$markers)

library(tidyverse)

azimuth_human_pbmc_markers_concise<-azimuth_human_pbmc_markers[,c(2,5)] %>% separate_rows(markers, sep = ',')


azimuth_human_pbmc_markers_concise_list<-unique(azimuth_human_pbmc_markers_concise$markers)

#head(rownames(Greten_10X_sc_objects.integrated))
#[1] "LINC01409" "LINC01128" "NOC2L" "ISG15" "TNFRSF18"  "TNFRSF4"

azimuth_human_pbmc_markers_concise_list_detected<-intersect(azimuth_human_pbmc_markers_concise_list,rownames(Greten_10X_sc_objects.integrated))


#length(unique(Idents(Greten_10X_sc_objects.integrated)))
#[1] 29


############
############
DefaultAssay(Greten_10X_sc_objects.integrated)<-"RNA"
############
############
pbmc_markers_aver<-AverageExpression(Greten_10X_sc_objects.integrated, slot="data",features=azimuth_human_pbmc_markers_concise_list_detected)
############
############
pbmc_markers_aver<-as.data.frame(pbmc_markers_aver[["RNA"]])
pbmc_markers_aver$Gene<-rownames(pbmc_markers_aver)
##pbmc_markers_aver$sumexp<-
pbmc_markers_aver<-pbmc_markers_aver[,c(ncol(pbmc_markers_aver),1:(ncol(pbmc_markers_aver)-1))]
pbmc_markers_aver$sum<-rowSums(pbmc_markers_aver[,c(2:ncol(pbmc_markers_aver))])
############
############

############
############
top2_aver_exp_levels_ratio <- matrix(nrow = nrow(pbmc_markers_aver), ncol = 1)

celltype_count<-length(c(2:30))

for (geneno in c(1:nrow(pbmc_markers_aver))) {

exp_levels_relev<-pbmc_markers_aver[geneno,c(2:30)]

exp_levels_relev<-unname(unlist(exp_levels_relev[which(exp_levels_relev>0)]))

n <- length(exp_levels_relev)

top2_aver_exp_levels_ratio[geneno,1]<-max(exp_levels_relev)/sort(exp_levels_relev,partial=n-1)[n-1]

}
############
############


azimuth_human_pbmc_markers_concise_tomerge<-azimuth_human_pbmc_markers_concise

colnames(azimuth_human_pbmc_markers_concise_tomerge)[2]<-"Gene"


pbmc_markers_aver_extended<-merge(pbmc_markers_aver,azimuth_human_pbmc_markers_concise_tomerge,by="Gene")

############
############
top_aver_exp_level_matching_label <- matrix(nrow = nrow(pbmc_markers_aver_extended), ncol = 1)

celltype_count<-length(c(2:30))

for (geneno in c(1:nrow(pbmc_markers_aver_extended))) {

exp_levels_relev<-pbmc_markers_aver_extended[geneno,c(2:30)]

exp_levels_relev_names<-names(exp_levels_relev[which(exp_levels_relev>0)])

exp_levels_relev<-unname(unlist(exp_levels_relev[which(exp_levels_relev>0)]))



n <- length(exp_levels_relev)

d1<-exp_levels_relev_names[which(exp_levels_relev==max(exp_levels_relev))]

d2<-pbmc_markers_aver_extended$label[geneno]

print(paste0(d1,"_",d2))

if (d1==d2) {
top_aver_exp_level_matching_label[geneno,1]<-"Yes"
}

}
############
############
pbmc_markers_aver_extended_selected<-pbmc_markers_aver_extended[which(top_aver_exp_level_matching_label=="Yes"),]
############
############
length(unique(pbmc_markers_aver_extended_selected$Gene))
############
############
df_label_counts_selected_with_matching_markers<-as.data.frame(table(pbmc_markers_aver_extended_selected$label))
############
############


###pick up to top 2 genes based on $Sum from ###pbmc_markers_aver_extended_selected$Gene per label
###


pbmc_markers_aver_extended_selected %>%
  group_by(label) %>%
top_n(n = 2, wt = sum) -> pbmc_markers_aver_extended_selected_top2_sum

################
################
Greten_10X_sc_objects.integrated$predicted.celltype<-as.factor(Greten_10X_sc_objects.integrated$predicted.celltype)


ind_reorder_test<-pbmc_markers_aver_extended_selected_top2_sum$label[match(levels(Greten_10X_sc_objects.integrated$predicted.celltype),pbmc_markers_aver_extended_selected_top2_sum$label)]


ind_reorder<-match(levels(Greten_10X_sc_objects.integrated$predicted.celltype),pbmc_markers_aver_extended_selected_top2_sum$label)

ind_reorder<-ind_reorder[!is.na(ind_reorder)]

#levels(Greten_10X_sc_objects.integrated$predicted.celltype)
################
################

Idents(Greten_10X_sc_objects.integrated)<-"predicted.celltype"

features_DotPlot<-unique(pbmc_markers_aver_extended_selected_top2_sum$Gene[ind_reorder])

p1<-Seurat::DotPlot(Greten_10X_sc_objects.integrated,features=features_DotPlot, assay="RNA" )+scale_color_viridis_c()+ RotatedAxis()

p1<-p1+theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7),legend.title=element_text(size=8),axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10),legend.text=element_text(size=8))+labs(title="Cell type-specific markers that overlap with the \n canonical markers in the Azimuth Human PBMC reference")+ theme(legend.position="top")

#p2<-#Seurat::DotPlot(Greten_10X_sc_objects.integrated,features=selec_genes[25:48#], assay="RNA" )+scale_color_viridis_c()+ RotatedAxis()

#p2<-p2+theme(axis.text.x = element_text(size = 8), axis.text.y = #element_text(size = 10),legend.title=element_text(size=8),axis.title.x #= element_text(size = 10),axis.title.y = element_text(size = #10),legend.text=element_text(size=8))
#+ theme(legend.position="top")


pdf(file = "ICB_scRNA-Seq_dotplot_cell_type_specific_Azimuth_Human_PBMC_ref_canonical_markers_oct31_2023.pdf", width=9, height=5) #width=10, height=11
#patchwork::wrap_plots(p1,p2, ncol = 1)
print(p1)
dev.off()
