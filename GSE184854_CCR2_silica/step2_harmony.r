library(dplyr)
library(cowplot)
library(Seurat)
library(harmony)



getwd()
path="G:/silicosis/geo/GSE184854_scRNA-seq_mouse_lung_ccr2/GSE184854_RAW_all_merged/step_2_harmony/"
dir.create(path)
setwd(path)
getwd()



#load("G:/silicosis/geo/GSE184854_scRNA-seq_mouse_lung_ccr2/GSE184854_RAW/mydata.rds")
load("G:/silicosis/geo/GSE184854_scRNA-seq_mouse_lung_ccr2/GSE184854_RAW_all_merged/merged.rds")
table(Idents(All.merge))
All=subset(All.merge,idents = c('A0301',  'A0302',   'A0303', 
                                'A0304',   'A0305' ,'A0306', 
                                'A0307' , 'A0308',   'A0309'))
dim(All)
table(Idents(All))

#All$percent.mt=PercentageFeatureSet(All,pattern = "^mt-")
pdf("1_contorl_QC.pdf")
VlnPlot(All, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
All <- subset(All, subset = nFeature_RNA > 200 & percent.mt < 20)## 15499,16133;

All = All%>%Seurat::NormalizeData(verbose = FALSE) %>%  
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE)
All = RunPCA(All, npcs = 50, verbose = FALSE)

pdf("2_ElbowPlot.pdf")
ElbowPlot(All, ndims = 50)
dev.off()


table(Idents(All))
All$stim=Idents(All)
#All@meta.data$stim <- c(rep("case", length(grep("1$", All@assays$RNA@counts@Dimnames[[2]]))), rep("ctrl", length(grep("2$", All@assays$RNA@counts@Dimnames[[2]])))) ## 8186,7947;
pdf("2_pre_harmony_harmony_plot.pdf")
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = All, reduction = "pca", pt.size = .1, group.by = "stim")
p2 <- VlnPlot(object = All, features = "PC_1", group.by = "stim", pt.size = .1)
plot_grid(p1, p2)
dev.off()
##########################run harmony

All <- All %>% RunHarmony("stim", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(All, 'harmony') 
pdf("1_contorl_QC_.pdf")
VlnPlot(All, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf("2_after_harmony_harmony_plot.pdf")
options(repr.plot.height = 5, repr.plot.width = 12)
p3 <- DimPlot(object = All, reduction = "harmony", pt.size = .1, group.by = "stim")
p4 <- VlnPlot(object = All, features = "harmony_1", group.by = "stim", pt.size = .1)
plot_grid(p3, p4)
dev.off()

#######################cluster
All <- All %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30)
All<-All%>% FindClusters(resolution = 3) %>% identity()

options(repr.plot.height = 4, repr.plot.width = 10)
pdf("3_after_harmony_umap_two_group.pdf")
DimPlot(All, reduction = "umap", group.by = "stim", pt.size = .1)
dev.off()
pdf("3_after_harmony_cluster_UMAP.pdf")
DimPlot(All, reduction = "umap", label = TRUE, pt.size = .1)
dev.off()
pdf("3_umap_samples_split.pdf")
DimPlot(All, reduction = "umap", pt.size = .1, split.by = "stim", label = T)
dev.off()
pdf("3_after_harmony_tsne_two_group.pdf")
DimPlot(All, reduction = "tsne", group.by = "stim", pt.size = .1)
dev.off()
pdf("3_after_harmony_cluster_tSNE.pdf")
DimPlot(All, reduction = "tsne", label = TRUE, pt.size = .1)
dev.off()
pdf("3_tSNE_samples_split.pdf")
DimPlot(All, reduction = "tsne", pt.size = .1, split.by = "stim", label = T)
dev.off()

#######################################################################
stat = as.matrix(table(Idents(All), All$stim))
write.table(stat, "cluster_stat.txt", sep = "\t", quote = F, col.names = T, row.names = T)
################################################################
All<-All%>% FindClusters(resolution = 2) %>% identity()

Disease.markers <- FindAllMarkers(All, min.pct = 0.25, logfc.threshold = 0.35, only.pos = T)
head(Disease.markers)
top100markers <- Disease.markers %>% group_by(cluster) %>% slice_max(avg_log2FC,n=100)  
openxlsx::write.xlsx(top100markers,'top100markers_for_all_res2.xlsx')
#write.table(top20markers, "top20_markers.txt", sep = "\t", quote = F, col.names = T, row.names = F)
#save(All, file = "G:/silicosis/geo/GSE184854_scRNA-seq_mouse_lung_ccr2/GSE184854_RAW_all_merged/step_2_harmony/all_harmony.rds")

load("G:/silicosis/geo/GSE184854_scRNA-seq_mouse_lung_ccr2/GSE184854_RAW_all_merged/step_2_harmony/all_harmony.rds")

getwd()
DimPlot(All,label = T)
DimPlot(All.merge,label = T)
FeaturePlot(All,features = c('Cd68'))
FeaturePlot(All,features = c('Lyz2'))
FeaturePlot(All,features = c('Mrc1'))
FeaturePlot(All,features = c('Adgre1'))
FeaturePlot(All,features = c('Axl'))

FeaturePlot(All,features = c('Csf1r'))
FeaturePlot(All,features = c('Ly6c2'))
FeaturePlot(All,features = c('Il7r'))
FeaturePlot(All,features = c('Car4'))
FeaturePlot(All,features = c('Siglecf'))
FeaturePlot(All,features = c('Mki67'))
FeaturePlot(All,features = c('Jchain','Tnfrsf17'),split.by = 'mystim')
DotPlot(All,features = c('Jchain','Tnfrsf17'),split.by = 'mystim')+RotatedAxis()+ggplot2::coord_flip()

colnames(All)
grepl(colnames(All),pattern = 'WT')
table(grepl(colnames(All),pattern = 'WT'))
table(grepl(colnames(All),pattern = 'CCR'))

All$orig_stim=Idents(All)
if(1==1){
  All$mystim_validate=ifelse(grepl(colnames(All),pattern = 'WT') &
                               (grepl(colnames(All),pattern = 'A0301')|grepl(colnames(All),pattern = 'A0302')|grepl(colnames(All),pattern = 'A0303')),
                             paste0(colnames(All),'_Day_21_WT'),
                             ifelse(grepl(colnames(All),pattern = 'WT') &
                                      (grepl(colnames(All),pattern = 'A0304')|grepl(colnames(All),pattern = 'A0305')|grepl(colnames(All),pattern = 'A0306'))
                                    ,paste0(colnames(All),'_Day_7_WT'),
                                    ifelse(grepl(colnames(All),pattern = 'WT') &
                                             (grepl(colnames(All),pattern = 'A0307')|grepl(colnames(All),pattern = 'A0308')|grepl(colnames(All),pattern = 'A0309'))
                                           ,paste0(colnames(All),'_Day_3_WT'),
                                           ifelse(grepl(colnames(All),pattern = 'CCR') &
                                                    (grepl(colnames(All),pattern = 'A0301')|grepl(colnames(All),pattern = 'A0302')|grepl(colnames(All),pattern = 'A0303'))
                                                  ,paste0(colnames(All),'_Day_21_CCR'),
                                                  ifelse(grepl(colnames(All),pattern = 'CCR') &
                                                           (grepl(colnames(All),pattern = 'A0304')|grepl(colnames(All),pattern = 'A0305')|grepl(colnames(All),pattern = 'A0306'))
                                                         ,paste0(colnames(All),'_Day_7_CCR'),
                                                         ifelse(grepl(colnames(All),pattern = 'CCR') &
                                                                  (grepl(colnames(All),pattern = 'A0307'))
                                                                ,paste0(colnames(All),'_Day_3_CCR'),paste0(colnames(All),'_Day_3_CCR'))
                                                  )
                                           ))
                             )
  )
  
  All$mystim=ifelse(grepl(colnames(All),pattern = 'WT') &
                      (grepl(colnames(All),pattern = 'A0301')|grepl(colnames(All),pattern = 'A0302')|grepl(colnames(All),pattern = 'A0303')),
                    paste0('_Day_21_WT'),
                    ifelse(grepl(colnames(All),pattern = 'WT') &
                             (grepl(colnames(All),pattern = 'A0304')|grepl(colnames(All),pattern = 'A0305')|grepl(colnames(All),pattern = 'A0306'))
                           ,paste0('_Day_7_WT'),
                           ifelse(grepl(colnames(All),pattern = 'WT') &
                                    (grepl(colnames(All),pattern = 'A0307')|grepl(colnames(All),pattern = 'A0308')|grepl(colnames(All),pattern = 'A0309'))
                                  ,paste0('_Day_3_WT'),
                                  ifelse(grepl(colnames(All),pattern = 'CCR') &
                                           (grepl(colnames(All),pattern = 'A0301')|grepl(colnames(All),pattern = 'A0302')|grepl(colnames(All),pattern = 'A0303'))
                                         ,paste0('_Day_21_CCR'),
                                         ifelse(grepl(colnames(All),pattern = 'CCR') &
                                                  (grepl(colnames(All),pattern = 'A0304')|grepl(colnames(All),pattern = 'A0305')|grepl(colnames(All),pattern = 'A0306'))
                                                ,paste0('_Day_7_CCR'),
                                                ifelse(grepl(colnames(All),pattern = 'CCR') &
                                                         (grepl(colnames(All),pattern = 'A0307'))
                                                       ,paste0('_Day_3_CCR'),paste0('_Day_3_CCR'))
                                         )
                                  ))
                    )
  )
  }
table(All$mystim)
table(All$orig_stim)
table(All$orig.ident)
table(grep(colnames(All),pattern = 'WT') & (grep(colnames(All),pattern = 'A0301')|
                                              grep(colnames(All),pattern = 'A0302')|
                                              grep(colnames(All),pattern = 'A0303'))
      )

All$mystim_2=ifelse(grepl(colnames(All),pattern = 'WT'),'WT','CCR')
table(All$mystim)
table(All$mystim_2)
table(ifelse(grepl(colnames(All),pattern = 'WT'),'WT','CCR'))
length(grep(colnames(All),pattern = 'WT'))
table(All$orig.ident)
Idents(All)=All$mystim_2
table(grepl(colnames(All),pattern = 'WT') )
grepl(c('af','faf','cd'),pattern = 'a')

getwd()

path="G:/silicosis/geo/GSE184854_scRNA-seq_mouse_lung_ccr2/step_2_harmony/"
setwd(path)
save(All,file ="ALL_with_stim_harmonized.rds" )







