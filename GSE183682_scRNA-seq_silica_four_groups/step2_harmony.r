rm(list = ls())
library(dplyr)
library(cowplot)
library(Seurat)
library(harmony)
Feature_RNA = 200
setwd("/data/home/yll/code_test/silicosis_scRNAseq/Step2_harmony_200")##	鏀规垚鎯冲瓨鏀炬枃浠剁殑璺緞
case_sparse <- read.csv("/data/home/yll/code_test/silicosis_scRNAseq/Step1_doubletfinder/case_doubletfinder.csv") ## 32285,8327;	鏀硅矾寰?
ctrl_sparse <- read.csv("/data/home/yll/code_test/silicosis_scRNAseq/Step1_doubletfinder/ctrl_doubletfinder.csv") ## 32285,8306;	鏀硅矾寰?


colnames(ctrl_sparse) <- sub(pattern = "1", replacement = "2", x = colnames(ctrl_sparse))
All_counts = cbind(case_sparse,ctrl_sparse) ## 32285 16633


All <- CreateSeuratObject(counts = All_counts, project = "silicosis", min.cells = 10)## 15499,16633;
All[["percent.mt"]] <- PercentageFeatureSet(All, pattern = "^mt-")

pdf("1_contorl_QC.pdf")
VlnPlot(All, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

All <- subset(All, subset = nFeature_RNA > Feature_RNA & percent.mt < 20)## 15499,16133;

All = All%>%Seurat::NormalizeData(verbose = FALSE) %>%  
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE)
All = RunPCA(All, npcs = 50, verbose = FALSE)

pdf("2_ElbowPlot.pdf")
ElbowPlot(All, ndims = 50)
dev.off()

All@meta.data$stim <- c(rep("case", length(grep("1$", All@assays$RNA@counts@Dimnames[[2]]))), rep("ctrl", length(grep("2$", All@assays$RNA@counts@Dimnames[[2]])))) ## 8186,7947;
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

Disease.markers <- FindAllMarkers(All, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
top20markers <- Disease.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)  
write.table(top20markers, "top20_markers.txt", sep = "\t", quote = F, col.names = T, row.names = F)
save(All, file = "silicosis_harmony.rds")


version




# Data in two numeric vectors
women_weight <- c(38.9, 61.2, 73.3, 21.8, 63.4, 64.6, 48.4, 48.8, 48.5)
men_weight <- c(67.8, 60, 63.4, 76, 89.4, 73.3, 67.3, 61.3, 62.4) 

wilcox.test(women_weight,men_weight)
wilcox.test(women_weight,men_weight,alternative = 'two.side')
wilcox.test(women_weight,men_weight,alternative = 'less')
wilcox.test(women_weight,men_weight,alternative = 'greater')

sessionInfo(package = 'harmony')
sessionInfo(package = "DoubletFinder")
sessionInfo(package = "Seurat")


getwd()
setEPS()#瀵煎嚭鐭㈤噺鍥? 
postscript("whatever.eps")
plot(rnorm(100), main="Hey Some Data")
dev.off()

