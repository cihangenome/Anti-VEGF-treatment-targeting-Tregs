setwd("/Users/oguzc/Downloads/Greten_10X_data/reda_figs_collec_aug10_2023/epitope_preds")

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


#alternatively spliced Insulin

###CDR3 (epitopes) are rows, columns are antigens, row annotations are cell types (cell type-time point pairs), and the color is proportional to the frequency of the clonotype?


df_CDR3_beta_epitopes<-read_excel("/Users/oguzc/Downloads/Greten_10X_data/reda_figs_collec_aug10_2023/epitope_preds/10X_ICB_epitope_predictions_CDR3_abundance_gt10_mar31_2023.xlsx",sheet = "CDR3_beta_epitopes")

df1_CDR3_beta_epitopes<-df_CDR3_beta_epitopes[grep("Homo s", df_CDR3_beta_epitopes$antigen),]

df2_CDR3_beta_epitopes<-df_CDR3_beta_epitopes[-grep("Homo s", df_CDR3_beta_epitopes$antigen),]


library(tidyverse)
#unique(df_CDR3_beta_epitopes[,-c(3,7)]) %>% separate_rows(epitope, sep = ',') -> df_CDR3_beta_epitopes_reformatted #[,-c(3,4,7)]

unique(df1_CDR3_beta_epitopes[,-c(3,7)]) %>% separate_rows(epitope,antigen, sep = ',') -> df1_CDR3_beta_epitopes_reformatted

unique(df2_CDR3_beta_epitopes[,-c(3,7)]) %>% separate_rows(epitope, sep = ',') -> df2_CDR3_beta_epitopes_reformatted

df_CDR3_beta_epitopes_reformatted<-rbind(df1_CDR3_beta_epitopes_reformatted,df2_CDR3_beta_epitopes_reformatted)

#df_CDR3_beta_epitopes_reformatted %>% separate_rows(antigen, sep = ',') -> df_CDR3_beta_epitopes_reformatted

df_CDR3_beta_epitopes_reformatted<-unique(df_CDR3_beta_epitopes_reformatted)


#length(unique(df_CDR3_beta_epitopes_reformatted$epitope)) #71
#length(unique(df_CDR3_beta_epitopes_reformatted$antigen)) #52
#length(unique(df_CDR3_beta_epitopes_reformatted$cdr3_beta)) #9
#length(unique(df_CDR3_beta_epitopes_reformatted$cell_type)) #3

#unique(df_CDR3_beta_epitopes_reformatted$cell_type)
#[1] "CD8 TEM_pre"  "CD8 TEM_post" "MAIT_pre"

#df_CDR3_beta_epitopes_reformatted<-df_CDR3_beta_epitopes_reformatted[grep("CD8",df_CDR3_beta_epitopes_reformatted$cell_type),]


#this was good, but just moving forward
df_CDR3_beta_epitopes_reformatted%>%
  separate("cell_type", c("cell_type", "time"), sep = "_")-> df_CDR3_beta_epitopes_reformatted


length(unique(df_CDR3_beta_epitopes_reformatted$epitope)) #67
length(unique(df_CDR3_beta_epitopes_reformatted$antigen)) #50
length(unique(df_CDR3_beta_epitopes_reformatted$cdr3_beta)) #7
length(unique(df_CDR3_beta_epitopes_reformatted$cell_type)) #1

df_CDR3_beta_epitopes_reformatted$antigen<-gsub("NUF2R","Kinetochore protein Nuf2",df_CDR3_beta_epitopes_reformatted$antigen)

#df_CDR3_beta_epitopes_reformatted$ "[Homo sapiens]",""
###CDR3 (epitopes) are rows, columns are antigens, row annotations are , and the color is proportional to the frequency of the clonotype?

#cell types (cell type-time point pairs)

df_liver_cancer_targets<-read_excel("/Users/oguzc/Downloads/Greten_10X_data/reda_figs_collec_aug10_2023/epitope_preds/liver_cancer_targets_mar31_2023.xlsx",sheet = "Sheet4")

library(Hmisc)



df_liver_cancer_targets$Name<-tolower(df_liver_cancer_targets$Name)
df_CDR3_beta_epitopes_reformatted$antigen<-tolower(df_CDR3_beta_epitopes_reformatted$antigen)

df_liver_cancer_targets$Name<-str_to_title(df_liver_cancer_targets$Name)
df_CDR3_beta_epitopes_reformatted$antigen<-str_to_title(df_CDR3_beta_epitopes_reformatted$antigen)

#capitalize str_to_title



inters_beta_antigens<-intersect(unique(df_liver_cancer_targets$Name),unique(df_CDR3_beta_epitopes_reformatted$antigen))

df_CDR3_beta_epitopes_reformatted_filtered<-df_CDR3_beta_epitopes_reformatted[which(df_CDR3_beta_epitopes_reformatted$antigen %in% inters_beta_antigens),]


df_CDR3_beta_epitopes_reformatted_filtered<-rbind(df_CDR3_beta_epitopes_reformatted_filtered,df_CDR3_beta_epitopes_reformatted[grep("liver",df_CDR3_beta_epitopes_reformatted$antigen, ignore.case =TRUE),])

df_CDR3_beta_epitopes_reformatted_filtered<-df_CDR3_beta_epitopes_reformatted_filtered[,-c(3)]

#length(unique(df_CDR3_beta_epitopes_reformatted_filtered$epitope))  #32
#length(unique(df_CDR3_beta_epitopes_reformatted_filtered$antigen))  #6
#length(unique(df_CDR3_beta_epitopes_reformatted_filtered$cdr3_beta))  #4
#length(unique(df_CDR3_beta_epitopes_reformatted_filtered$cell_type))  #2


###CDR3 (epitopes) are rows, columns are CDR3s, row annotations are cell types & antigens (no time point info needed), and the color is proportional to the prediction score?

# unique(df_CDR3_beta_epitopes_reformatted_filtered$antigen)
# [1] "Kinetochore Protein Nuf2"
# [2] "C-C Motif Chemokine 3"
# [3] "Kinesin-Like Protein Kif23"
# [4] "Nlr Family Card Domain-Containing Protein 4"
# [5] "Cyclin-Dependent Kinase 4"
# [6] "Soluble Liver Antigen/Liver Pancreas Antigen"


###Epitopes are rows, columns are CDR3s, row annotations are cell types & antigens (no time point info needed), and the color is proportional to the prediction score?

library(pheatmap)
library(ComplexHeatmap)


df_CDR3_beta_epitopes_reformatted_filtered_toplot<-df_CDR3_beta_epitopes_reformatted_filtered

df_CDR3_beta_epitopes_reformatted_filtered_toplot$epitope_antigen_cdr3_beta<- paste0(df_CDR3_beta_epitopes_reformatted_filtered_toplot$epitope," (", df_CDR3_beta_epitopes_reformatted_filtered_toplot$antigen,"): ", df_CDR3_beta_epitopes_reformatted_filtered_toplot$cdr3_beta)

#, "_", df_CDR3_beta_epitopes_reformatted_filtered_toplot$cdr3_beta


df_CDR3_beta_epitopes_reformatted_filtered_toplot<-unique(df_CDR3_beta_epitopes_reformatted_filtered_toplot[,-c(1,4,5)])


df_CDR3_beta_epitopes_reformatted_filtered_toplot<-df_CDR3_beta_epitopes_reformatted_filtered_toplot[order(-df_CDR3_beta_epitopes_reformatted_filtered_toplot$score),]

df_CDR3_beta_epitopes_reformatted_filtered_toplot<-df_CDR3_beta_epitopes_reformatted_filtered_toplot[!duplicated(df_CDR3_beta_epitopes_reformatted_filtered_toplot[,c(1,3)]),]

p_dotplot_epitopes<-ggballoonplot(df_CDR3_beta_epitopes_reformatted_filtered_toplot, y = "epitope_antigen_cdr3_beta",x ="cell_type",size = 5, fill = "score", ggtheme = theme_bw()) + scale_fill_gradientn(colors = my_cols)+ggtitle("Liver cancer-specific antigens and their epitopes \n predicted to be recognized by expanded TRB clonotypes with >10 abundance") + labs(x = "Cell type")

p_dotplot_epitopes<-ggpar(p_dotplot_epitopes, font.title =12, xlab="Cell type",ylab="Epitope (Antigen): TRB CDR3 sequence",legend.title = "prediction score",caption="y-axis shows Epitope (Antigen): TRB CDR3 sequence, x-axis shows the cell type in which the TRB clonotype is detected. \n All clonotypes except the MAIT hit are detected both at the pre- and post- treatment stages \n Antigens that are liver cancer-specific targets are derived from Pharos (Illuminating the Druggable Genome) database https://pharos.nih.gov")

pdf(file = "NCBR358_dotplot_10X_ICB_predicted_liver_cancer_specific_epitopes_sep14_2023.pdf", width=12, height=4)
print(p_dotplot_epitopes)
dev.off()
