#here we show how to find doublets using two groups, as we did in four groups
rm(list = ls())
library(dplyr)
library(cowplot)
library(Seurat)
library(DoubletFinder)
##################https://github.com/ddiez/DoubletFinder
setwd("/data/home/yll/code_test/silicosis_scRNAseq/Step1_doubletfinder")##	鏀规垚鎯冲瓨鏀炬枃浠剁殑璺緞
########################################################
case.data <- Read10X(data.dir = "/data/biomath/Spatial_Transcriptome/silicosis/scRNA_seq/no_aggr/exp-1-5/2.2.filtered_feature_bc_matrix")## 閲岄潰鏈変笁涓枃浠?32285  8859	鏀?
ctrl.data <- Read10X(data.dir = "/data/biomath/Spatial_Transcriptome/silicosis/scRNA_seq/no_aggr/con-1-5/2.2.filtered_feature_bc_matrix")## 32285  8836	鏀?
pre <- CreateSeuratObject(counts = cbind(case.data, ctrl.data), project = "case")
pre[["percent.mt"]] <- PercentageFeatureSet(pre, pattern = "^mt-")
pdf("1_pre_doubletfinder.pdf")
VlnPlot(pre, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()
#######################################################
case <- CreateSeuratObject(counts = case.data, project = "case")
load_number <- dim(case)
case <- NormalizeData(case, normalization.method = "LogNormalize", scale.factor = 10000)
case <- FindVariableFeatures(case, selection.method = "vst", nfeatures = 2000)
case <- ScaleData(case)
case <- RunPCA(case)
case <- RunUMAP(case, dims = 1:30)
#####pK identification####################
sweep.res.list <- paramSweep_v3(case, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
nExp_poi <- round(0.06*ncol(case))        ## 8000~0.061 
case <- doubletFinder_v3(case, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
paste("DF.classifications_", "0.25_", mpK, "_", nExp_poi, sep="")
pdf("1_doubletFinder_case.pdf")
DimPlot(case, pt.size = 1, label = TRUE, label.size = 5, reduction = "umap", group.by = "DF.classifications_0.25_0.005_532")
dev.off()
case_filter <- subset(case, DF.classifications_0.25_0.005_532 == "Singlet" )## 32285  8327
DoubletFinder_number <- dim(case_filter)
##############
write.table(as.matrix(GetAssayData(object = case_filter, slot = "counts")),'case_doubletfinder.csv', sep = ',', row.names = T, col.names = T, quote = F)
write.table(rbind(load_number,DoubletFinder_number), 'case_filter_cell_number.csv', sep = ',', row.names = T, col.names = F, quote = F)
#########################################################
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "control")
load_number<-dim(ctrl)
ctrl <- NormalizeData(ctrl, normalization.method = "LogNormalize", scale.factor = 10000)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)
ctrl <- ScaleData(ctrl)
ctrl <- RunPCA(ctrl)
ctrl <- RunUMAP(ctrl, dims = 1:30)
#####pK identification####################
sweep.res.list <- paramSweep_v3(ctrl, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
nExp_poi <- round(0.06*ncol(ctrl))        ## 8000~0.061 
ctrl <- doubletFinder_v3(ctrl, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
paste("DF.classifications_", "0.25_", mpK, "_", nExp_poi, sep="")
pdf("1_doubletFinder_control.pdf")
DimPlot(ctrl, pt.size = 1, label=TRUE, label.size = 5, reduction = "umap", group.by = "DF.classifications_0.25_0.005_530")
dev.off()
ctrl_filter <- subset(ctrl, DF.classifications_0.25_0.005_530 == "Singlet")  ##32285  8306
DoubletFinder_number <- dim(ctrl_filter)
##############
write.table(as.matrix(GetAssayData(object = ctrl_filter, slot = "counts")),'ctrl_doubletfinder.csv', sep = ',', row.names = T, col.names = T, quote = F)
write.table(rbind(load_number,DoubletFinder_number),'ctrl_filter_cell_number.csv', sep = ',', row.names = T, col.names = F, quote = F)
###########################################################
after.data <- cbind(as.matrix(GetAssayData(object = case_filter, slot = "counts")), 
                    as.matrix(GetAssayData(object = ctrl_filter, slot = "counts")))
after <- CreateSeuratObject(counts = after.data, project = "case")
after[["percent.mt"]] <- PercentageFeatureSet(after, pattern = "^mt-")
pdf("1_after_doubletfinder.pdf")
VlnPlot(after, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()
