


getwd()
path="G:/silicosis/geo/GSE104154_scRNA-seq_fibrotic MC_bleomycin/normalized" #
dir.create(path)
setwd(path)
getwd()
list.files()



raw_counts=read.csv("G:\\silicosis\\geo\\GSE104154_scRNA-seq_fibrotic MC_bleomycin\\GSE104154_d0_d21_sma_tm_Expr_nor\\GSE104154_d0_d21_sma_tm_Expr_norm.csv")

head(raw_counts)[1:4,1:4]
counts=raw_counts[,-1]
head(counts)[1:4,1:4]
rownames(counts)=counts$symbol

head(raw_counts)[1:4,1:4]
counts=raw_counts[,-2]
head(counts)[1:4,1:4]
rownames(counts)=counts$id
counts=counts[,-1]
head(counts)[1:4,1:4]

library(Seurat)
#https://zhuanlan.zhihu.com/p/385206713
rawdata=CreateSeuratObject(counts = counts,project = "blem",assay = "RNA")

ids=raw_counts[,1:2]
head(ids)
colnames(ids)= c('ENSEMBL','SYMBOL')
head(ids)
dim(ids) # [1] 16428 
ids=na.omit(ids)
dim(ids) # [1] 15504 
length(unique(ids$SYMBOL)) # [1] 15494 
# 
# 
ids=ids[!duplicated(ids$SYMBOL),]
ids=ids[!duplicated(ids$ENSEMBL),]
dim(ids)
pos=match(ids$ENSEMBL,rownames(rawdata) )
hp_sce=rawdata[pos,]
hp_sce
#rownames(hp_sce) = ids$SYMBOL
# RenameGenesSeurat  -----------------------------------------------
#
RenameGenesSeurat <- function(obj , 
                              newnames ) { 
  # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. 
  # It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

hp_sce=RenameGenesSeurat(obj = hp_sce, 
                         newnames = ids$SYMBOL)
getwd()
#save(hp_sce,file = 'first_sce.Rdata')



hp_sce
rownames(hp_sce)[grepl('^mt-',rownames(hp_sce))]
rownames(hp_sce)[grepl('^Rp[sl]',rownames(hp_sce))]

hp_sce[["percent.mt"]] <- PercentageFeatureSet(hp_sce, pattern = "^mt-")
fivenum(hp_sce[["percent.mt"]][,1])
rb.genes <- rownames(hp_sce)[grep("^Rp[sl]",rownames(hp_sce))]
C<-GetAssayData(object = hp_sce, slot = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
hp_sce <- AddMetaData(hp_sce, percent.ribo, col.name = "percent.ribo")

getwd()




plot1 <- FeatureScatter(hp_sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hp_sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

VlnPlot(hp_sce, features = c("percent.ribo", "percent.mt"), ncol = 2)
VlnPlot(hp_sce, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(hp_sce, features = c("percent.ribo", "nCount_RNA"), ncol = 2)
hp_sce

hp_sce1 <- subset(hp_sce, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 20)
hp_sce1


sce=hp_sce1
sce
colnames(sce)
grep(colnames(sce),pattern = ".1")
grep(colnames(sce),pattern = ".2")
sce@meta.data$stim <-c(rep("PBS", length(grep("1$", sce@assays$RNA@counts@Dimnames[[2]]))), 
                       rep("PBS", length(grep("2$", sce@assays$RNA@counts@Dimnames[[2]]))),
                       rep("PBS", length(grep("3$", sce@assays$RNA@counts@Dimnames[[2]]))),
                       rep("Bleomycin", length(grep("4$", sce@assays$RNA@counts@Dimnames[[2]]))),
                       rep("Bleomycin", length(grep("5$", sce@assays$RNA@counts@Dimnames[[2]]))),
                       rep("Bleomycin", length(grep("6$", sce@assays$RNA@counts@Dimnames[[2]])))
) ## 8186,7947;

table(sce$stim)

library(dplyr)
sce[["RNA"]]@meta.features <- data.frame(row.names = rownames(sce[["RNA"]]))

All = sce%>%Seurat::NormalizeData(verbose = FALSE) %>%  
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE)
All = RunPCA(All, npcs = 50, verbose = FALSE)

pdf("2_ElbowPlot.pdf")
ElbowPlot(All, ndims = 50)
dev.off()

library(cowplot)
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


pdf("2_after_harmony_harmony_plot.pdf")
options(repr.plot.height = 5, repr.plot.width = 12)
p3 <- DimPlot(object = All, reduction = "harmony", pt.size = .1, group.by = "stim")
p4 <- VlnPlot(object = All, features = "harmony_1", group.by = "stim", pt.size = .1)
plot_grid(p3, p4)
dev.off()
#############cluster
#library(harmony)
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
getwd()
#save(All,file ="G:/silicosis/geo/GSE104154_scRNA-seq_fibrotic MC_bleomycin/normalized/All_normolized_for_clustering.rds" )
load("G:/silicosis/geo/GSE104154_scRNA-seq_fibrotic MC_bleomycin/normalized/All_normolized_for_clustering.rds")



###FIGURE.s4I s4J############################################################
##################################################################################################################################################################
DimPlot(All,label = F,reduction = 'umap',group.by = 'stim')
FeaturePlot(All,features = 'C1qa',reduction = 'tsne')


DimPlot(All,label = T,reduction = 'tsne')
getwd()
Disease.markers <- FindAllMarkers(All, min.pct = 0.35, logfc.threshold = 0.35, only.pos = T)
openxlsx::write.xlsx(Disease.markers,file ="G:/silicosis/geo/GSE104154_scRNA-seq_fibrotic MC_bleomycin/normalized/markers_normolized_for_all.xlsx" )

top20markers <- Disease.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)  
#write.table(top20markers, "top20_markers.txt", sep = "\t", quote = F, col.names = T, row.names = F)
#save(All, file = "sepsis_harmony.rds")


FeaturePlot(All,features = 'Cd68',reduction = 'tsne')
FeaturePlot(All,features = 'Mrc1',reduction = 'tsne')
FeaturePlot(All,features = 'Mrc1',reduction = 'umap')
FeaturePlot(All,features = 'C1qa',reduction = 'tsne')
FeaturePlot(All,features = 'C1qb',reduction = 'tsne')
FeaturePlot(All,features = 'C1qc',reduction = 'tsne')
FeaturePlot(All,features = 'Spp1',reduction = 'tsne')
FeaturePlot(All,features = 'Ear2',reduction = 'tsne')
FeaturePlot(All,features = 'Ear1',reduction = 'tsne')
FeaturePlot(All,features = 'Mmp12',reduction = 'tsne')
FeaturePlot(All,features = 'Mmp14',reduction = 'tsne')
FeaturePlot(All,features = 'Gpnmb',reduction = 'tsne')

subset_data=subset(All,idents = c('0','17','18','22','29','39'))
DimPlot(subset_data,label = T,reduction = 'tsne')
subset_data$orig_cluster_from_all=Idents(subset_data)

subset_data=subset_data %>% RunHarmony("stim", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(All, 'harmony') 
dim(subset_data)
subset_data <- subset_data %>% 
  RunUMAP(reduction = "harmony", dims = 1:22) %>% 
  RunTSNE(reduction = "harmony", dims = 1:22) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:22)
subset_data<-subset_data%>% FindClusters() %>% identity()


DimPlot(subset_data,reduction = 'tsne')
DimPlot(subset_data)
DimPlot(subset_data,reduction = 'tsne',split.by = 'stim')

#
for (res in seq(0.2,1,0.1)) {
  subset_data=FindClusters(subset_data,graph.name = 'RNA_snn',resolution = res,algorithm = 1)
}
apply(subset_data@meta.data[,grep('RNA_snn_res',colnames(subset_data@meta.data))], 2, table)

library(clustree)
p5_tree=clustree::clustree(subset_data@meta.data,prefix='RNA_snn_res.')
p5_tree

#
ggplot(subset_data@meta.data, aes(x=RNA_snn_res.0.2, fill=stim)) + geom_bar(position = "fill")
ggplot(subset_data@meta.data, aes(x=orig_cluster_from_all, fill=stim)) + geom_bar(position = "fill")


Idents(subset_data)=subset_data$orig_cluster_from_all
markers=FindAllMarkers(subset_data,logfc.threshold = 0.5,only.pos = T,min.pct = 0.3)
DimPlot(subset_data,label = T,reduction = 'tsne')


subset_data=subset(All,idents = c('0','18','22','29'))
DimPlot(subset_data,label = T,reduction = 'tsne')
markers=FindAllMarkers(subset_data,logfc.threshold = 0.5,only.pos = T,min.pct = 0.3)
DimPlot(subset_data,label = T,reduction = 'tsne')

subset_data$orig_cluster_from_all=Idents(subset_data)
ggplot(subset_data@meta.data, aes(x=orig_cluster_from_all, fill=stim)) + geom_bar(position = "fill")

subset_data=RenameIdents(subset_data,'0'='AM1',
                         '22'='AM2','18'='AM3',
                         '29'='IM')

DimPlot(subset_data,label = T,reduction = 'tsne')
#
ggplot(subset_data@meta.data, aes(x=Idents(subset_data), fill=stim)) + geom_bar(position = "fill")


DimPlot(subset_data,label = T,reduction = 'tsne')
DimPlot(subset_data,group.by = 'stim')
getwd()
#save(subset_data,file = "G:/silicosis/geo/GSE104154_scRNA-seq_fibrotic MC_bleomycin/normalized/IM_AMs.rds")
load("G:/silicosis/geo/GSE104154_scRNA-seq_fibrotic MC_bleomycin/normalized/IM_AMs.rds")
library(Seurat)

DimPlot(subset_data,label = T)
FeaturePlot(subset_data,features = "Pdgfa")


subset_data=subset(subset_data,idents = c('AM1','AM2','AM3'))
DimPlot(subset_data,label = T,reduction = 'tsne')
#姣斾緥鍥?
ggplot(subset_data@meta.data, aes(x=Idents(subset_data), fill=stim)) + geom_bar(position = "fill")


DimPlot(subset_data,label = T,reduction = 'tsne')
DimPlot(subset_data,group.by = 'stim')
getwd()

subset_data$orig_cluster_from_IM_AMs=Idents(subset_data)

subset_data=subset_data %>% RunHarmony("stim", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(All, 'harmony') 
dim(subset_data)
subset_data <- subset_data %>% 
  RunUMAP(reduction = "harmony", dims = 1:22) %>% 
  RunTSNE(reduction = "harmony", dims = 1:22) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:22)
subset_data<-subset_data%>% FindClusters() %>% identity()


DimPlot(subset_data,reduction = 'tsne',label = T)
DimPlot(subset_data,label = T,label.size = 5)
DimPlot(subset_data,reduction = 'tsne',split.by = 'stim')
DimPlot(subset_data,group.by = 'stim')
DimPlot(subset_data,split.by = 'stim',label = T,label.size = 5)

#
getwd()
ggplot(subset_data@meta.data, aes(x=Idents(subset_data), fill=stim)) + geom_bar(position = "fill")
markers=FindAllMarkers(subset_data,logfc.threshold = 0.5,only.pos = T,min.pct = 0.3)
openxlsx::write.xlsx(markers,file ="G:/silicosis/geo/GSE104154_scRNA-seq_fibrotic MC_bleomycin/normalized/makers_for_AM1_AM2_AM3.xlsx" )
Idents(subset_data)=subset_data$RNA_snn_res.0.8

subset_data=RenameIdents(subset_data,'1'='AM1','4'='AM1',
                         '0'='AM1',
                         '2'='AM2',
                         '3'='AM3','5'='AM3')

getwd()
#save(subset_data,file = "G:/silicosis/geo/GSE104154_scRNA-seq_fibrotic MC_bleomycin/normalized/only_AMs.rds")
load("G:/silicosis/geo/GSE104154_scRNA-seq_fibrotic MC_bleomycin/normalized/only_AMs.rds")

library(Seurat)

##############################FIGURE s4k
DimPlot(subset_data,split.by = 'stim')

####################FIGURE s4L################################################3
############################################################################################
ggplot(subset_data@meta.data, aes(x=Idents(subset_data), fill=stim)) + geom_bar(position = "fill")

DimPlot(subset_data)
DimPlot(subset_data,reduction = "tsne")

FeaturePlot(subset_data,features = "Csf1",pt.size = 1.5)



# Adjust the contrast in the plot
FeaturePlot(subset_data, features = "Csf1",pt.size = 1.5,
            min.cutoff = 0.11)
DotPlot(subset_data,features = c("Csf2rb","Csf1","Csf1r"))
getwd()
setwd("G:/silicosis/geo/GSE104154_scRNA-seq_fibrotic MC_bleomycin/normalized")
jpeg("csf1.png",res = 300,units = "in",height = 7,width = 7)
p=FeaturePlot(subset_data, features = "Csf1",pt.size = 2,
            min.cutoff = 0.11)
print(p)
dev.off()


####
load("G:/silicosis/geo/GSE104154_scRNA-seq_fibrotic MC_bleomycin/normalized/only_AMs.rds")

markers=FindAllMarkers(subset_data,only.pos = T)
#openxlsx::write.xlsx(markers,file = "G:/silicosis/geo/GSE104154_scRNA-seq_fibrotic MC_bleomycin/normalized/makers_for_AM1_AM2_AM3.xlsx")
markers=openxlsx::read.xlsx("G:/silicosis/geo/GSE104154_scRNA-seq_fibrotic MC_bleomycin/normalized/makers_for_AM1_AM2_AM3.xlsx")
head(markers)

mymarkers=list()
for (eachcluster_markers in unique(markers$cluster)) {
  #eachcluster_markers="AM1"
  mymarkers[[paste0(eachcluster_markers)]]=markers[markers$cluster==eachcluster_markers,"gene"] 
} 
head(mymarkers)

fibrosis_gene=openxlsx::read.xlsx("G:/silicosis/fibrosis_geness/harmonizome_gfibrosis.xlsx")
head(fibrosis_gene)
fibrosis_gene=fibrosis_gene$Symbol



#intersect_fibrosis_gene_with_AM
intersected_fibrosis_genes=list()
for (eachcluster_AM in names(mymarkers)) {
  intersected_fibrosis_genes[[eachcluster_AM]]=intersect(mymarkers[[eachcluster_AM]],fibrosis_gene)
}

intersected_fibrosis_genes

lapply(intersected_fibrosis_genes,length)


intersected_fibrosis_genes_as_a_vector=as.vector(unlist(intersected_fibrosis_genes))

################################FIGRUE.S4M
marker_exp=AverageExpression(subset_data, features = intersected_fibrosis_genes_as_a_vector, return.seurat = TRUE)
DoHeatmap(marker_exp, features = intersected_fibrosis_genes_as_a_vector,
          label=TRUE, group.bar = TRUE, draw.lines = FALSE)

getwd()
pdf('heatmap.pdf',width = 7,height = 14)
p=DoHeatmap(marker_exp, features = intersected_fibrosis_genes_as_a_vector,
            label=TRUE, group.bar = TRUE, draw.lines = FALSE)
print(p)
dev.off()



#supplementary table S7

PF_silicosis_bleo=openxlsx::read.xlsx("G:\\silicosis\\sicosis\\Supplementary Table.xlsx",
                                      sheet = "S7_hPF_Silica_Bleomycin")
head(PF_silicosis_bleo)

library(Hmisc)
library(dplyr)
disease_markers_for_AM3=list()
for (eachdisease in colnames(PF_silicosis_bleo)) {
  #eachdisease="Silicosis,.Cluster.AM3"
  disease_markers_for_AM3[[paste0(eachdisease)]]=PF_silicosis_bleo[,eachdisease] %>%
    capitalize()%>% tolower() %>% na.omit() %>% unique()  %>%capitalize()
}
head(disease_markers_for_AM3)

#https://www.jianshu.com/p/58f5429b5402

genes <- Reduce(intersect, list(disease_markers_for_AM3[[1]],
                                disease_markers_for_AM3[[2]],
                                disease_markers_for_AM3[[3]])) %>% capitalize()
length(genes)
getwd()


library(VennDiagram)#https://www.jianshu.com/p/b5a4c40c3a33
venn_list <- list(Human_PF = disease_markers_for_AM3[[1]],
                  Silicosis = disease_markers_for_AM3[[2]],
                  Bleomycin=disease_markers_for_AM3[[3]])

###################FIGRUE S4N################3333
venn.diagram(venn_list, filename = 'venn2_for_humanPF-Silicosi-Bleomycin.png', imagetype = 'png', 
             fill = c('red', 'blue','orange'), alpha = 0.50, cat.col = rep('black', 3), 
             col = 'black', cex = 1.3, fontfamily = 'serif', 
             cat.cex = 1.3, cat.fontfamily = 'serif')

getwd()
#缁х画浠ヤ笂杩?4涓垎缁勪负渚嬶紝缁勯棿浜ら泦鍏冪礌鑾峰緱
inter <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], 'venn4_inter.txt', row.names = FALSE, sep = '\t', quote = FALSE)


for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
openxlsx::write.xlsx(inter[-c(5, 6)], 'venn4_inter.xlsx', row.names = FALSE, sep = '\t', quote = FALSE)




