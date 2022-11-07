
library(patchwork)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(openxlsx)
library(Hmisc)
#https://www.jianshu.com/p/cef5663888ff
getwd()
path="G:/silicosis/sicosis/silicosis_ST/overlapped_map/addmodule_allmarkers_from_findallmarkers_allmerge_ams_im"
dir.create(path)
setwd(path)
getwd()

load("G:/silicosis/sicosis/yll/macrophage/no cluster2/0.3/pure_cluster3_in_allmerge-IM/silicosis_cluster_merge.rds")
table(All.merge$new.cluster.idents)
#姣斾緥鍥?
markers=FindAllMarkers(All.merge,min.pct = 0.75,logfc.threshold = 0.8,only.pos = T)
head(markers)
unique(markers$cluster)
library(stringr)
Myselectedmarekrs=markers %>% filter(str_detect(markers$cluster,"AM"))
DotPlot(All.merge,features=Myselectedmarekrs$gene[1:30])+RotatedAxis()
#openxlsx::write.xlsx(markers,file = "G:/silicosis/sicosis/silicosis_ST/overlapped_map/addmodule/markers_forallmerge_ams_im.xlsx")

library(openxlsx)
load("G:/silicosis/sicosis/silicosis_ST/yll/0214/harmony_cluster/d_all/silicosis_ST_harmony_SCT_r0.6.rds")
load("G:/silicosis/sicosis/silicosis_ST/yll/0214/harmony_cluster/d_all/silicosis_ST_harmony_SCT_r0.6.rds")
#marker = read.xlsx("G:/silicosis/sicosis/silicosis_ST/overlapped_map/Rigional and cell markers.xlsx", 
 #                  sheet = "SingleCell_markers")

#markers=read.xlsx('G:/silicosis/sicosis/yll/macrophage/no cluster2/0.3/findmarkers_1and2/30cluster_markers.xlsx')
#markers=read.xlsx("G:/silicosis/sicosis/silicosis-1122-merge/silicosis_cluster_merge_markers.xlsx")
markers=read.xlsx("G:/silicosis/sicosis/silicosis_ST/overlapped_map/addmodule/markers_forallmerge_ams_im.xlsx")

head(markers)
library(dplyr)

markers=markers %>% group_by(cluster)  %>% slice_head(n=50) %>%select(cluster,gene)
head(markers)
library(reshape2)
markers2=dcast(markers,gene~cluster) 
head(markers2)
markers2[is.na(markers2)]<-0
head(markers2)
markers2=markers2[,-1]

marker=markers2
head(marker)
colnames(marker)[19]="Myofibroblast-vascular smooth muscle cell"
cellnames=colnames(marker)                ##number=length(marker[,cellname])

library(Hmisc)
getwd()
#path="G:/silicosis/sicosis/silicosis_ST/overlapped_map/addmodule_allmarkers_from_findallmarkers"
path="G:/silicosis/sicosis/silicosis_ST/overlapped_map/addmodule_allmarkers_from_findallmarkers_allmerge_ams_im-2.0"
dir.create(path)
setwd(path)
getwd()

library(Seurat)
for (each in cellnames) {
  #each='Myofibroblast-vascular smooth muscle cell'
  cellname=each
  
  mymarker=marker[,paste0(cellname)]   %>% na.exclude() %>% unique() %>%
     list()         #capitalize()    %>%
  number=length(mymarker[[1]])
  unlist(mymarker)
  
  #瀵圭粰瀹氱殑鍩哄洜闆嗗悎杩涜鎵撳垎  骞剁敾鍥?
  if(1==1){
    d.all=AddModuleScore(d.all,
                                features = mymarker,
                                name = paste0(cellname))
  #缁撴灉淇濆瓨鍦ㄨ繖閲?
  colnames(d.all@meta.data)
  head(d.all@meta.data)
  colnames(d.all@meta.data)[[9]]=paste0(cellname)
  
  ###
  
  p1=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image")+
    ggtitle(paste0("SiO2_7d")) +scale_fill_gradientn(guide = "colourbar",colours=c("blue","yellow","red"), limits=c(0,3)) #sio27d
  p2=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image.1")+ggtitle(paste0("NS_7d")) +scale_fill_gradientn(guide = "colourbar",colours=c("blue","yellow","red"), limits=c(0,3))
  p3=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image.2")+ ggtitle(paste0("SiO2_56d")) +scale_fill_gradientn(guide = "colourbar",colours=c("blue","yellow","red"), limits=c(0,3))
  p4=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image.3")+ggtitle(paste0(("NS_56d"))) +scale_fill_gradientn(guide = "colourbar",colours=c("blue","yellow","red"), limits=c(0,3))
  
  jpeg(paste0(ifelse(grep(paste0(cellname),pattern = "/"),"Myofibroblast-vascular smooth muscle cell",paste0(cellname))
              ,paste0(cellname),"_","total_",length(unlist(mymarker)),"_",paste0(min(number),"-",max(number)), 
              paste(unlist(mymarker)[1:15],collapse = "_"),"_.jpeg"), #鍙彇鍓?15涓?
       height = 12, width = 12, units = 'in', res=600)
  p=ggpubr::ggarrange(p2,p1,p4,p3,ncol = 2,nrow =2)
  print(p)
  dev.off()
  d.all@meta.data=d.all@meta.data[,1:8] 
  }
}


#scale.bar 0-2
for (each in cellnames) {
  #each='Myofibroblast-vascular smooth muscle cell'
  cellname=each
  
  mymarker=marker[,paste0(cellname)]   %>% na.exclude() %>% unique() %>%
    list()         #capitalize()    %>%
  number=length(mymarker[[1]])
  unlist(mymarker)
  
  #瀵圭粰瀹氱殑鍩哄洜闆嗗悎杩涜鎵撳垎  骞剁敾鍥?
  if(1==1){
    d.all=AddModuleScore(d.all,
                         features = mymarker,
                         name = paste0(cellname))
    #缁撴灉淇濆瓨鍦ㄨ繖閲?
    colnames(d.all@meta.data)
    head(d.all@meta.data)
    colnames(d.all@meta.data)[[9]]=paste0(cellname)
    
    ###
    
    p1=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image")+
      ggtitle(paste0("SiO2_7d")) +scale_fill_gradientn(guide = "colourbar",colours=c("blue","yellow","red"), limits=c(0,2)) #sio27d
    p2=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image.1")+ggtitle(paste0("NS_7d")) +scale_fill_gradientn(guide = "colourbar",colours=c("blue","yellow","red"), limits=c(0,2))
    p3=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image.2")+ ggtitle(paste0("SiO2_56d")) +scale_fill_gradientn(guide = "colourbar",colours=c("blue","yellow","red"), limits=c(0,2))
    p4=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image.3")+ggtitle(paste0(("NS_56d"))) +scale_fill_gradientn(guide = "colourbar",colours=c("blue","yellow","red"), limits=c(0,2))
    
    jpeg(paste0(ifelse(grep(paste0(cellname),pattern = "/"),"Myofibroblast-vascular smooth muscle cell",paste0(cellname))
                ,paste0(cellname),"_","total_",length(unlist(mymarker)),"_",paste0(min(number),"-",max(number)), 
                paste(unlist(mymarker)[1:15],collapse = "_"),"_.jpeg"), #鍙彇鍓?15涓?
         height = 12, width = 12, units = 'in', res=600)
    p=ggpubr::ggarrange(p2,p1,p4,p3,ncol = 2,nrow =2)
    print(p)
    dev.off()
    d.all@meta.data=d.all@meta.data[,1:8] 
  }
}


#scale.bar custome
for (each in cellnames) {
  #each="Fibroblast"
  cellname=each
  
  mymarker=marker[,paste0(cellname)]   %>% na.exclude() %>% unique() %>%
    list()         #capitalize()    %>%
  number=length(mymarker[[1]])
  unlist(mymarker)
  
  #瀵圭粰瀹氱殑鍩哄洜闆嗗悎杩涜鎵撳垎  骞剁敾鍥?
  if(1==1){
    d.all=AddModuleScore(d.all,
                         features = mymarker,
                         name = paste0(cellname))
    #缁撴灉淇濆瓨鍦ㄨ繖閲?
    colnames(d.all@meta.data)
    head(d.all@meta.data)
    colnames(d.all@meta.data)[[9]]=paste0(cellname)
    
    ###
    
    p1=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image")+
      ggtitle(paste0("SiO2_7d")) +scale_fill_gradientn(guide = "colourbar",colours=c("blue","yellow","red"), limits=c(0.9,2.1)) #sio27d
    p2=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image.1")+ggtitle(paste0("NS_7d")) +scale_fill_gradientn(guide = "colourbar",colours=c("blue","yellow","red"), limits=c(0.9,2.1))
    p3=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image.2")+ ggtitle(paste0("SiO2_56d")) +scale_fill_gradientn(guide = "colourbar",colours=c("blue","yellow","red"), limits=c(0.9,2.1))
    p4=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image.3")+ggtitle(paste0(("NS_56d"))) +scale_fill_gradientn(guide = "colourbar",colours=c("blue","yellow","red"), limits=c(0.9,2.1))
    
    jpeg(paste0(ifelse(grep(paste0(cellname),pattern = "/"),"Myofibroblast---vascular smooth muscle cell",paste0(cellname))
                ,paste0(cellname),"_","total_",length(unlist(mymarker)),"_",paste0(min(number),"-",max(number)), 
                paste(unlist(mymarker)[1:15],collapse = "_"),"_.jpeg"), #鍙彇鍓?15涓?
         height = 12, width = 12, units = 'in', res=600)
    p=ggpubr::ggarrange(p2,p1,p4,p3,ncol = 2,nrow =2)
    print(p)
    dev.off()
    d.all@meta.data=d.all@meta.data[,1:8] 
  }
}

for (each in c("Neutrophil","NK cell","T cell")) {
  #each='Myofibroblast/vascular smooth muscle cell'
  cellname=each
  
  mymarker=marker[,paste0(cellname)]   %>% na.exclude() %>% unique() %>%
    list()         #capitalize()    %>%
  number=length(mymarker[[1]])
  unlist(mymarker)
  
  #瀵圭粰瀹氱殑鍩哄洜闆嗗悎杩涜鎵撳垎  骞剁敾鍥?
  if(1==1){d.all=AddModuleScore(d.all,
                                features = mymarker,
                                name = paste0(cellname))
  #缁撴灉淇濆瓨鍦ㄨ繖閲?
  colnames(d.all@meta.data)
  head(d.all@meta.data)
  colnames(d.all@meta.data)[[9]]=paste0(cellname)
  
  ###
  
  p1=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image")+ ggtitle(paste0("SiO2_7d"))  #sio27d
  p2=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image.1")+ggtitle(paste0("NS_7d"))
  p3=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image.2")+ ggtitle(paste0("SiO2_56d"))
  p4=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image.3")+ggtitle(paste0(("NS_56d")))
  
  jpeg(paste0(ifelse(grep(paste0(cellname),pattern = "/"),"Myofibroblast-vascular smooth muscle cell",paste0(cellname))
              ,paste0(cellname),"_","total_",length(unlist(mymarker)),"_",paste0(min(number),"-",max(number)), 
              paste(unlist(mymarker)[1:15],collapse = "_"),"_.jpeg"), #鍙彇鍓?15涓?
       height = 12, width = 12, units = 'in', res=600)
  p=ggpubr::ggarrange(p2,p1,p4,p3,ncol = 2,nrow =2)
  print(p)
  dev.off()
  d.all@meta.data=d.all@meta.data[,1:8] }
}


#鍙ns56 鍜宻io2_56d
for (each in cellnames) {
  #each='Myofibroblast/vascular smooth muscle cell'
  cellname=each
  
  mymarker=marker[,paste0(cellname)]   %>% na.exclude() %>% unique() %>%
    list()         #capitalize()    %>%
  number=length(mymarker[[1]])
  unlist(mymarker)
  
  #瀵圭粰瀹氱殑鍩哄洜闆嗗悎杩涜鎵撳垎  骞剁敾鍥?
  if(1==1){d.all=AddModuleScore(d.all,
                                features = mymarker,
                                name = paste0(cellname))
  #缁撴灉淇濆瓨鍦ㄨ繖閲?
  colnames(d.all@meta.data)
  colnames(d.all@meta.data)[[9]]=paste0(cellname)
  ###
  
 # p1=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image")+ ggtitle(paste0("SiO2_7d"))  #sio27d
 # p2=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image.1")+ggtitle(paste0("NS_7d"))
  p3=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image.2")+ ggtitle(paste0("SiO2_56d"))
  p4=SpatialFeaturePlot(d.all, features = paste0(cellname), slot = "scale.data",images = "image.3")+ggtitle(paste0(("NS_56d")))
  
  jpeg(paste0(paste0(cellname),"_","total_",length(unlist(mymarker)),"_",paste0(min(number),"-",max(number)), 
              paste(unlist(mymarker)[1:15],collapse = "_"),"_.jpeg"), #鍙彇鍓?15涓?
       height = 12, width = 12, units = 'in', res=600)
  p=ggpubr::ggarrange(p4,p3,ncol = 1,nrow =2)
  print(p)
  dev.off()}
}



