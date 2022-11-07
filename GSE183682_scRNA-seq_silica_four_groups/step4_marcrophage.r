
#################macrophage-1  鍚堝苟macrophage-2
load("G:/silicosis/sicosis/zxb/1207/Macrophgage/silicosis_Macrophgage.rds")## 17226 4965
table(Idents(subset_data))
as.matrix(table(Idents(subset_data), subset_data$stim))
levels(Idents(subset_data))
table(subset_data$stim)
unique(subset_data$stim)
Idents(subset_data)
subset_data$stim
DimPlot(subset_data,label = TRUE)

for (res in seq(0.1,1,0.2)) {
  subset_data=FindClusters(subset_data,graph.name = "RNA_snn",resolution = res,algorithm = 1)
}
apply(subset_data@meta.data[,grep("RNA_snn_res",colnames(subset_data@meta.data))], 2, table)
getwd()
library(clustree)
p5_tree=clustree::clustree(subset_data@meta.data,prefix="RNA_snn_res.")
p5_tree
getwd()

save(subset_data,file = "macrophage-various-resolution.rds")

Idents(subset_data)<-subset_data$RNA_snn_res.0.3
Idents(subset_data)
DimPlot(subset_data,reduction = "umap",label = TRUE)
library(Seurat)
mysubset=subset(subset_data, idents = c(2), invert = T)
DimPlot(mysubset,reduction = "umap",label = TRUE)

getwd()
path="G:/silicosis/sicosis/yll/macrophage/no cluster2"
dir.create(path)
setwd(path)
subset_data=mysubset

subset_data = subset_data %>%
  Seurat::NormalizeData(verbose = FALSE) %>%  
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 50, verbose = FALSE)
subset_data@meta.data$stim <- c(rep("SiO2_7", length(grep("1$",colnames(subset_data)))),rep("SiO2_56", length(grep("2$",colnames(subset_data)))), rep("NS_7", length(grep("3$",colnames(subset_data)))),rep("NS_56", length(grep("4$",colnames(subset_data)))))
subset_data <- subset_data %>% RunHarmony("stim", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(subset_data, 'harmony') 
dims = 1:30
subset_data <- subset_data %>% 
  RunUMAP(reduction = "harmony", dims = dims) %>% 
  RunTSNE(reduction = "harmony", dims = dims) %>% 
  FindNeighbors(reduction = "harmony", dims = dims)
subset_data$cell.type = as.character(Idents(All.merge)[colnames(subset_data)])
Idents(subset_data) = factor(subset_data$cell.type)
Idents(subset_data)
DefaultAssay(subset_data)
subset_data@graphs

getwd()
save(subset_data,file = "macrophage_clean.rds")
load(file = "macrophage_clean.rds")



##
dev.off()
file="G:/silicosis/sicosis/yll/macrophage/no cluster2"
resolution = c(0.05,0.1,0.3,0.4,0.5,0.7,0.9)
for(r in resolution){
  load(paste(file, "macrophage_clean.rds", sep="/"))
  
  dir.create(paste(file,r,sep="/"))
  setwd(paste(file,r,sep="/"))
  getwd()
  
  
  subset_data <- FindClusters(subset_data, resolution = r)
  save(subset_data, file=paste0("silicosis_Macrophgage_r", r, ".rds"))
  
  pdf("1_Macrophgage_cluster_TSNE.pdf")
  p = DimPlot(subset_data, reduction = "tsne", label = TRUE, pt.size = 1,label.size = 6)
  print(p)
  dev.off()
  
  pdf("1_Macrophgage_cluster_UMAP.pdf")
  p = DimPlot(subset_data, reduction = "umap", label = TRUE, pt.size = 1,label.size = 6)
  print(p)
  dev.off()
  
  stat = cbind(table(Idents(subset_data), subset_data$stim), All = rowSums(table(Idents(subset_data), subset_data$stim)))
  stat = rbind(stat, All = colSums(stat))
  write.xlsx(stat, "2_Macrophgage_cluster_stat.xlsx", col.names=T, row.names=T)
  markers <- FindAllMarkers(subset_data,  min.pct = 0.1, logfc.threshold = 0.25, only.pos=T)
  write.xlsx(markers,"3_Macrophgage_cluster_markers.xlsx", col.names=T,row.names=F)
}
