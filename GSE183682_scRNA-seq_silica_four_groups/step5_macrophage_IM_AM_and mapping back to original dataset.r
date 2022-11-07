library(CellChat)
library(patchwork)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(openxlsx)
library(harmony)
library(dplyr)

getwd()
path="G:/silicosis/sicosis/yll/macrophage/no cluster2/0.3/pure_cluster013_in_allmerge"
dir.create(path)
setwd(path)
#https://www.jianshu.com/p/cef5663888ff
load("G:/silicosis/sicosis/silicosis-1122-merge/silicosis_cluster_merge.rds")##	鏀硅矾寰?
library(Seurat)
table(All.merge$cell.type)
DotPlot(All.merge,features = c('Csf1r','Mafb','Cd38'))


load("G:/silicosis/sicosis/yll/macrophage/no cluster2/macrophage_clean.rds")
table(subset_data$orig.ident)
table(Idents(subset_data))
DimPlot(subset_data,label = TRUE)
getwd()


#鎻愬彇macrophage_clean鐨処M
IMcluster=subset(subset_data,idents = c('3'))
table(IMcluster$orig.ident)

#macrophage_clean  AM鍒嗙兢鎯呭喌
if(TRUE){
  DimPlot(subset_data,label = TRUE)
  FeaturePlot(subset_data,features = c('Fabp4','Hmox1','Ctsk','Fabp5',
                                       'Chil3','S100a1','Wfdc21'))
  DotPlot(subset_data,features = c('Fabp4','Hmox1','Ctsk',
                                   'Chil3','S100a1','Wfdc21')
  )
  
  my015subsetdata=subset(subset_data,idents = c('0','1','4','5'))
  subset_data=my015subsetdata
  #閲嶆柊鍒嗙粍 鍘婚櫎鎵规鏁堝簲
  if(1==1){#閲嶆柊鍒嗙粍 鍘婚櫎鎵规鏁堝簲
    subset_data[["percent.mt"]] <- PercentageFeatureSet(subset_data, pattern = "^mt-")
    VlnPlot(subset_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    subset_data = subset_data %>%
      Seurat::NormalizeData(verbose = FALSE) %>%  
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
      ScaleData(verbose = FALSE) %>%
      RunPCA(npcs = 50, verbose = FALSE)
    ElbowPlot(subset_data, ndims = 50)
    VlnPlot(subset_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    subset_data <- subset_data %>% RunHarmony("stim", plot_convergence = TRUE)
    harmony_embeddings <- Embeddings(subset_data, 'harmony') 
    #######################cluster
    dims = 1:30
    subset_data <- subset_data %>% 
      RunUMAP(reduction = "harmony", dims = dims) %>% 
      RunTSNE(reduction = "harmony", dims = dims) %>% 
      FindNeighbors(reduction = "harmony", dims = dims)
  }
  
  getwd()
  DimPlot(subset_data,label = TRUE)
  table(Idents(subset_data))
  
  
  
  #R涓璼eurat鍖匘imPlot濡備綍鏍规嵁涓嶅悓鐨剅esolution鐢诲浘 seurat鏌ョ湅鏌愪竴鍒嗚鲸鐜囦笅鐨刣implot
  #鏍戠姸鍥?
  library(clustree)
  for (res in seq(0,0.8,0.1)) {
    subset_data=FindClusters(subset_data,graph.name = "RNA_snn",resolution = res,algorithm = 1)
  }
  subset_data=FindClusters(subset_data,resolution = 0.25)
  apply(subset_data@meta.data[,grep("RNA_snn_res",colnames(subset_data@meta.data))],2,table)
  
  mytree=clustree(subset_data@meta.data,prefix="RNA_snn_res.")
  mytree
  
  while (dev.off()) {
    dev.off()
  }
  DimPlot(subset_data,label = TRUE)
  DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.5")
  DotPlot(subset_data,group.by="RNA_snn_res.0.5",features = c('Car4','Ctsk','Chil3','S100a1','Wfdc21',
                                                              'Fabp5','Fabp4','Hmox1'))
  DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.5")
  DimPlot(subset_data,label = TRUE,split.by = 'stim',group.by="RNA_snn_res.0.4")
  
  
  
  #鍘绘帀绾跨矑浣撳惈閲忛珮鐨刢6
  subset_data=subset(subset_data,idents = c('0','1','2','3','4','5','7')) #鍘绘帀绾跨矑浣撳惈閲忛珮鐨刢6
  #閲嶆柊鍒嗙粍 鍘婚櫎鎵规鏁堝簲
  if(1==1){subset_data = subset_data %>%
    Seurat::NormalizeData(verbose = FALSE) %>%  
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE)
  ElbowPlot(subset_data, ndims = 50)
  VlnPlot(subset_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  subset_data <- subset_data %>% RunHarmony("stim", plot_convergence = TRUE)
  harmony_embeddings <- Embeddings(subset_data, 'harmony') 
  #######################cluster
  dims = 1:30
  subset_data <- subset_data %>% 
    RunUMAP(reduction = "harmony", dims = dims) %>% 
    RunTSNE(reduction = "harmony", dims = dims) %>% 
    FindNeighbors(reduction = "harmony", dims = dims)}
  
  
  subset_data=FindClusters(subset_data)
  DimPlot(subset_data,label = TRUE)
  DimPlot(subset_data,label = TRUE,split.by = 'stim')
  
  
  
  
  
  #0.5鍒嗚鲸鐜囦笅 鍚堝苟04   鍚堝苟23 
  if(1==1){table(Idents(subset_data))
    Idents(subset_data)=subset_data@meta.data$RNA_snn_res.0.5
    table(Idents(subset_data))
    subset_data=RenameIdents(subset_data,'4'='0','3'='2')
    table(Idents(subset_data))
    AM_markers012=FindAllMarkers(subset_data,min.pct = 0.1,logfc.threshold = 0.5,
                                 only.pos = TRUE)
    getwd()
    write.xlsx(AM_markers012,file = "G:/silicosis/sicosis/yll/macrophage/no cluster2/0.3/pure_cluster013_in_allmerge/AM_Marker4-0=3-2.xlsx")
  }
  
  
  #0.3鍒嗚鲸鐜囦笅锛屽悎骞?3 0
  if(1==1){table(Idents(subset_data))
    Idents(subset_data)=subset_data@meta.data$RNA_snn_res.0.3
    table(Idents(subset_data))
    DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.3")
    subset_data=RenameIdents(subset_data,'3'='0') #鍚堝苟涓夌兢鍜?0缇? 鍗? 鎶?3鍙樹负0
    table(Idents(subset_data))}
  
  #0.3鍒嗚鲸鐜囦笅锛屽悎骞?2 0
  if(1==1){table(Idents(subset_data))
    Idents(subset_data)=subset_data@meta.data$RNA_snn_res.0.3
    table(Idents(subset_data))
    DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.3")
    DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.5")
    DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.4")
    subset_data=RenameIdents(subset_data,'2'='0') #鍚堝苟2缇ゅ拰0缇? 鍗? 鎶?2鍙樹负0  鎶?
    table(Idents(subset_data))
    subset_data=RenameIdents(subset_data,'3'='2')
    DimPlot(subset_data,label = TRUE)
    DimPlot(subset_data,label = TRUE,split.by='stim')
    subset_data=RenameIdents(subset_data,'0'='AM1',
                             '2'='AM2','1'='AM3')
    #姣斾緥鍥?
    ggplot(subset_data@meta.data, aes(x=Idents(subset_data), fill=stim)) + geom_bar(position = 'fill')
    #save(subset_data,file ="G:/silicosis/sicosis/yll/macrophage/no cluster2/0.3/pure_cluster013_in_allmerge/macrophage0.3combine2-0.rds" )
    load(file = "G:/silicosis/sicosis/yll/macrophage/no cluster2/0.3/pure_cluster013_in_allmerge/macrophage0.3combine2-0.rds")
    table(subset_data$orig.ident)
    table(Idents(subset_data))
    
  }
  
  getwd()
  #0.7鍒嗚鲸鐜囦笅锛屽悎骞?3 2 0   4鍙樹负2
  if(1==1){table(Idents(subset_data))
    Idents(subset_data)=subset_data@meta.data$RNA_snn_res.0.7
    table(Idents(subset_data))
    DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.3")
    DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.5")
    DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.4")
    DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.6")
    DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.7")
    subset_data=RenameIdents(subset_data,'2'='0','3'='0','4'='2') #鍚堝苟2缇ゅ拰0缇? 鍗? 鎶?2鍙樹负0
    table(Idents(subset_data))
  }
  
  
  #0.5鍒嗚鲸鐜囦笅锛屽悎骞?3 2 0   4鍙樹负2
  if(1==1){table(Idents(subset_data))
    Idents(subset_data)=subset_data@meta.data$RNA_snn_res.0.7
    table(Idents(subset_data))
    DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.2")
    DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.25")
    DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.3")
    DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.5")
    DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.4")
    DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.6")
    DimPlot(subset_data, reduction = "umap",label = TRUE,group.by="RNA_snn_res.0.7")
    subset_data=RenameIdents(subset_data,'2'='0','3'='0','4'='2') #鍚堝苟2缇ゅ拰0缇? 鍗? 鎶?2鍙樹负0
    table(Idents(subset_data))
  }
  #save(subset_data,file="G:/silicosis/sicosis/yll/macrophage/no cluster2/0.3/pure_cluster013_in_allmerge/macrophage_0145_noMT.rds")
  load("G:/silicosis/sicosis/yll/macrophage/no cluster2/0.3/pure_cluster013_in_allmerge/macrophage_0145_noMT.rds")
  
  
  #姣斾緥鍥?
  ggplot(subset_data@meta.data, aes(x=RNA_snn_res.0.3, fill=stim)) + geom_bar(position = 'fill')
  
  ggplot(subset_data@meta.data, aes(x=Idents(subset_data), fill=stim)) + geom_bar(position = 'fill')
  
  FeaturePlot(subset_data,features = c('C1qa','C1qb','C1qc'))
  
  #鏌ョ湅markers
  #Idents(subset_data)=subset_data@meta.data$RNA_snn_res.0.3
  DimPlot(subset_data,label = TRUE)
  AM_markers012=FindAllMarkers(subset_data,min.pct = 0.1,logfc.threshold = 0.5,
                               only.pos = TRUE)
  getwd()
  write.xlsx(AM_markers012,file = "G:/silicosis/sicosis/yll/macrophage/no cluster2/0.3/pure_cluster013_in_allmerge/AM_Markero.3r_combine2-0.xlsx")
  
  
}


#鍚堝苟IM AM
mysubsetdata=merge(subset_data,IMcluster)
table(Idents(mysubsetdata))
mysubsetdata=RenameIdents(mysubsetdata,'3'='IM')
table(Idents(mysubsetdata))

subset_data=mysubsetdata
table(Idents(subset_data))

subset_data$new_cluster_id=Idents(subset_data)
table(Idents(All.merge))

Idents(subset_data)=subset_data$new_cluster_id #clusterIM  AM1 AM2 AM3
#閲嶆柊鍒嗙粍 鍘婚櫎鎵规鏁堝簲
if(1==1){subset_data = subset_data %>%
  Seurat::NormalizeData(verbose = FALSE) %>%  
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 50, verbose = FALSE)
ElbowPlot(subset_data, ndims = 50)
VlnPlot(subset_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
subset_data <- subset_data %>% RunHarmony("stim", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(subset_data, 'harmony') 
#######################cluster
dims = 1:30
subset_data <- subset_data %>% 
  RunUMAP(reduction = "harmony", dims = dims) %>% 
  RunTSNE(reduction = "harmony", dims = dims) %>% 
  FindNeighbors(reduction = "harmony", dims = dims)
}

getwd()
DimPlot(subset_data,label = TRUE)
table(subset_data$new_cluster_id)
#save(subset_data,file ="G:/silicosis/sicosis/yll/macrophage/no cluster2/0.3/pure_cluster3_in_allmerge-IM/macrophage_im_am.rds" )
load(file = "G:/silicosis/sicosis/yll/macrophage/no cluster2/0.3/pure_cluster3_in_allmerge-IM/macrophage_im_am.rds")

FeaturePlot(subset_data,features = "Pdgfa")
#NS_7缁勭殑AM3缁勫彧鏈変竴涓粏鑳烇紝鍐嶆敼涓€涓粏鑳炵殑idents  浣夸箣鎴愪负NS_7 AM3
if(1==1){#NS_7缁勭殑AM3缁勫彧鏈変竴涓粏鑳烇紝鍐嶆敼涓€涓粏鑳炵殑idents  浣夸箣鎴愪负NS_7 AM3
  library(Seurat)
  table(subset_data$stim,subset_data$new_cluster_id)
  DimPlot(subset_data,group.by = 'new_cluster_id',split.by = 'stim',pt.size = 3)
  
  #https://zhuanlan.zhihu.com/p/146236794
  plot=DimPlot(subset_data,group.by = 'new_cluster_id',split.by = 'stim',pt.size = 3)
  HoverLocator(plot = plot, information = FetchData(subset_data, vars = c("ident", "stim"))) 
  select.cells <- CellSelector(plot = plot) 
  subset_data <- CellSelector(plot = plot, object = subset_data, ident = "AM3") 
  
  
  DimPlot(subset_data)
  DimPlot(subset_data,
          cells.highlight ="TGCTGCTGTACTTCTT.3",pt.size = 2 )
  
  colnames(subset_data@assays$RNA@data[,Idents(subset_data)=='AM3'])
  colnames(subset_data@assays$RNA@data[,subset_data@meta.data$stim=="NS_7"])
  intersect(colnames(subset_data@assays$RNA@data[,Idents(subset_data)=='AM3']),
            colnames(subset_data@assays$RNA@data[,subset_data@meta.data$stim=="NS_7"]))
  
  Idents(subset_data)
  subset_data <- SetIdent(object = subset_data, #鍚屼竴涓猻tim鍐呯殑鎵嶅彲浠ユ敼
                          cells = c(
                            "TGCTGCTGTACTTCTT.3",
                            "AATCGGTTCTACCTGC.3"
                          ), value = 'AM3')
  table(Idents(subset_data),subset_data$stim)
  table(subset_data$new_cluster_id,subset_data$stim)
  getwd()
 # save(subset_data,file ="G:/silicosis/sicosis/yll/macrophage/no cluster2/0.3/pure_cluster3_in_allmerge-IM/macrophage_im_am_setidnets_ns7d.rds" )
load(file ="G:/silicosis/sicosis/yll/macrophage/no cluster2/0.3/pure_cluster3_in_allmerge-IM/macrophage_im_am_setidnets_ns7d.rds")
  }





cellname_mycluster=list()
for (cluster in levels(Idents(subset_data))) {
  #cluster=0
  mycluster=cluster
  cellname_mycluster[[paste(cluster)]]=colnames(subset_data[,subset_data$new_cluster_id==mycluster])
  print(paste('cluster' ,"_",cluster,"====="));print(length(cellname_mycluster[[paste(cluster)]]))
}

names(cellname_mycluster)

DimPlot(All.merge,label = TRUE, cells.highlight=cellname_mycluster[["0"]])+ggtitle(paste("cluster_","0"))
DimPlot(All.merge,label = TRUE, cells.highlight=cellname_mycluster[["1"]])+ggtitle(paste("cluster_","1"))
DimPlot(All.merge,label = TRUE, cells.highlight=cellname_mycluster[["3"]])+ggtitle(paste("cluster_","3"))
DimPlot(All.merge,label = TRUE, cells.highlight=cellname_mycluster[["4"]])+ggtitle(paste("cluster_","4"))
DimPlot(All.merge,label = TRUE, cells.highlight=cellname_mycluster[["5"]])+ggtitle(paste("cluster_","5"))
length(cellname_mycluster[["0"]])

Idents(All.merge)=All.merge$cell.type
levels(All.merge)
rownames(All.merge@meta.data[Idents(All.merge)=="Monocyte",])
All.merge@meta.data$new.cluster.idents=ifelse( Idents(All.merge)=="Monocyte","Monocyte",
                                               ifelse( Idents(All.merge)=="T cell","T cell",
                                                       ifelse( Idents(All.merge)=="B cell","B cell",
                                                               ifelse( Idents(All.merge)=="Ig-producing B cell","Ig-producing B cell",
                                                                       ifelse( Idents(All.merge)=="Dendritic cell","Dendritic cell",
                                                                               ifelse( Idents(All.merge)=="Neutrophil","Neutrophil",
                                                                                       ifelse( Idents(All.merge)=="NK cell","NK cell",
                                                                                               ifelse(Idents(All.merge)=="Endothelial cell-1","Endothelial cell-1",   
                                                                                                      ifelse( Idents(All.merge)=="Endothelial cell-2","Endothelial cell-2",
                                                                                                              ifelse( Idents(All.merge)=="Fibroblast","Fibroblast",
                                                                                                                      ifelse(Idents(All.merge)=="Myofibroblast/vascular smooth muscle cell","Myofibroblast/vascular smooth muscle cell",
                                                                                                                             ifelse(Idents(All.merge)=="Cycling basal cell","Cycling basal cell",
                                                                                                                                    ifelse(Idents(All.merge)=="Ciliated cell","Ciliated cell",
                                                                                                                                           ifelse(Idents(All.merge)=="Clara cell","Clara cell",
                                                                                                                                                  ifelse(Idents(All.merge)=="AT1 cell","AT1 cell",
                                                                                                                                                         ifelse(Idents(All.merge)=="AT2 cell-1","AT2 cell-1",
                                                                                                                                                                ifelse(Idents(All.merge)=="AT2 cell-2","AT2 cell-2",
                                                                                                                                                                       ifelse(Idents(All.merge)=="Igha+ AT2  cell","Igha+ AT2  cell",
ifelse(  rownames(All.merge@meta.data) %in% union(cellname_mycluster[["IM"]],seq(1,ncol(All.merge)-length(cellname_mycluster[["IM"]]),1) ),"IM",
 ifelse(rownames(All.merge@meta.data) %in% union(cellname_mycluster[["AM1"]],seq(1,ncol(All.merge)-length(cellname_mycluster[["AM1"]]),1) ),"AM1",
 ifelse(rownames(All.merge@meta.data) %in% union(cellname_mycluster[["AM2"]],seq(1,ncol(All.merge)-length(cellname_mycluster[["AM2"]]),1) ),"AM2",
 ifelse(rownames(All.merge@meta.data) %in% union(cellname_mycluster[["AM3"]],seq(1,ncol(All.merge)-length(cellname_mycluster[["AM3"]]),1) ),"AM3","unknow_Macro")
                                                                                                                                                                                       ))   
   )))))))))))))))))))

table(All.merge$new.cluster.idents)
table(All.merge$cell.type,All.merge$orig.ident)
Idents(All.merge)=All.merge$new.cluster.idents
summary(All.merge$cell.type)

DimPlot(All.merge,label = TRUE)
levels(Idents(All.merge))
setdiff(levels(Idents(All.merge)),'unknow_Macro') #鍘婚櫎'unknow_Macro'
All.merge=subset(All.merge,idents = setdiff(levels(Idents(All.merge)),'unknow_Macro'))

getwd()
#save(All.merge,file = "G:/silicosis/sicosis/yll/macrophage/no cluster2/0.3/pure_cluster3_in_allmerge-IM/silicosis_cluster_merge.rds")
getwd()
load("G:/silicosis/sicosis/yll/macrophage/no cluster2/0.3/pure_cluster3_in_allmerge-IM/silicosis_cluster_merge.rds")

FeaturePlot(All.merge,features = c('Cd68','Mrc1','Lyz2','Adgre1','Ax1'))
FeaturePlot(All.merge,features = c('Mrc1'))


getwd()
#鐢诲浘
jpeg(paste0("Mrc1" ,".jpeg"),height = 9, width = 9, units = 'in', res=300)
p=FeaturePlot(All.merge,features = c('Mrc1'),pt.size = 1)
print(p)
dev.off()

table(Idents(All.merge))
table(All.merge$new.cluster.idents)
All.merge=RenameIdents(All.merge,'AM3'='AM','AM2'='AM','AM1'='AM'
                       )
markers=FindAllMarkers(All.merge,min.pct = 0.25,logfc.threshold = 0.6,only.pos = TRUE)

write.xlsx(markers,'markers.xlsx')

DimPlot(All.merge, label = T, 
        cells.highlight = WhichCells(All.merge, 
                                     idents = c("IM", "AM")))

DimPlot(All.merge, 
        cells.highlight = WhichCells(All.merge, 
                                     idents = c("IM", "AM")))

DimPlot(All.merge, label = T, 
        cells.highlight = WhichCells(All.merge, 
                                     idents = c("IM", "AM")), 
        cols.highlight = c("darkblue", "darkred"), cols = "grey")
    

    
jpeg(paste0("allmerge", ".jpeg"),height = 9, width = 18, units = 'in', res=300)
p=DimPlot(All.merge,label = TRUE)    
print(p)
dev.off()

jpeg(paste0("mrc1", ".jpeg"),height = 9, width = 9, units = 'in', res=300)
p=FeaturePlot(All.merge,features = 'Mrc1',pt.size = 1)    
print(p)
dev.off()


DimPlot(subset_data,label = TRUE)
jpeg(paste0("macro_all", ".jpeg"),height = 9, width = 9, units = 'in', res=300)
p=DimPlot(subset_data,label = TRUE,label.size = 8)    
print(p)
dev.off()


#姣斾緥鍥?
ggplot(subset_data@meta.data, aes(x=Idents(subset_data), fill=stim)) + geom_bar(position = 'fill')

DimPlot(subset_data,group.by = 'stim')
DimPlot(subset_data,split.by = 'stim',label = TRUE,pt.size = 1,label.size = 8)




##
load(file = "G:/silicosis/sicosis/yll/macrophage/no cluster2/0.3/pure_cluster3_in_allmerge-IM/macrophage_im_am.rds")

#姣斾緥鍥?
ggplot(subset_data@meta.data, aes(x=stim, fill=Idents(subset_data))) + geom_bar(position = 'fill')


subset_data$stim = factor(subset_data$stim, levels=c("NS_7","SiO2_7","NS_56","SiO2_56"))#鏀瑰彉椤哄簭

table(subset_data$stim,Idents(subset_data))


