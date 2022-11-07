rm(list=ls())
library(dplyr)
library(cowplot)
library(Seurat)
library(ggplot2)
library(openxlsx)
require(Matrix)
require(magrittr)
library("reshape")
library("RColorBrewer")

#浠庡幓闄よ繃鍙岀粏鑳烇紝骞朵笖缁忚繃harmony澶勭悊鐨剆epsis_harmony.rds鑾峰緱鏈€缁堢殑鍒嗙兢silicosis_cluster_merge.rds


getwd()
file = "D:/ARDS_scripts_1012/ARDS/Step2_harmony_f200_R3/0805/cluster_merge"##	鏀规垚鎯冲瓨鏀炬枃浠剁殑璺緞
dir.create(file)
setwd(file)
load("D:/ARDS_scripts_1012/silicosis_harmony.rds")##	鏀硅矾寰?
table(Idents(All))
as.matrix(table(Idents(All), All$stim))

DimPlot(All, reduction = "tsne", label = T, pt.size = .1, label.size = 6)
table(All$stim)
All$stim[All$stim=="case"] = "Exp"
All$stim[All$stim=="ctrl"] = "Con"
table(All$stim)
##	1.merge
All.merge = All

cluster = read.xlsx("D:/ARDS_scripts_1012/ARDS_cluster_anno_0805.xlsx")
cluster_sort = cluster[order(cluster$new), ]
new.cluster.ids = cluster$new
names(new.cluster.ids) <- levels(All.merge)
All.merge <- RenameIdents(All.merge, new.cluster.ids)
Idents(All.merge) = factor(Idents(All.merge),levels=1:length(levels(All.merge)))
table(Idents(All.merge))
as.matrix(table(Idents(All.merge), All.merge$stim))
All.merge <- subset(x = All.merge, idents = 16, invert = TRUE)##	16133-5025=11108
All.merge$cluster.sort = Idents(All.merge)

cell.type = rep(NA, ncol(All.merge))

for(i in 1:length(levels(All.merge)))
{
    cell.type[which(Idents(All.merge) == levels(All.merge)[i])] = unique(cluster_sort$cell.type)[i]
}

All.merge$cell.type = factor(cell.type,levels=unique(cluster_sort$cell.type)[-16])
Idents(All.merge) = All.merge$cell.type
table(All.merge$cell.type)
as.matrix(table(All.merge$cell.type, All.merge$stim))
stat = cbind(table(All.merge$cell.type, All.merge$stim), All = rowSums(table(All.merge$cell.type, All.merge$stim)))
stat = rbind(stat, All = colSums(stat))
write.xlsx(stat, "1_silicosis_cluster_merge_stat.xlsx", col.names=T, row.names=T)
legend_labels = paste(All.merge$cluster.sort, All.merge$cell.type, sep=" ")
All.merge$legend_labels = factor(legend_labels, levels=paste(levels(All.merge$cluster.sort), levels(All.merge$cell.type)))

save(All.merge, file="silicosis_cluster_merge.rds")



#load("silicosis_cluster_merge.rds")

pdf("1_cluster_merge_UMAP.pdf", width=10)
DimPlot(All.merge, reduction = "umap", label = T, pt.size = .1, label.size = 6, group.by = "cluster.sort") +
scale_color_hue(labels = levels(All.merge$legend_labels)) +
guides(color = guide_legend(ncol = 1,override.aes = list(size = 4))) +
theme(legend.key.size = unit(0.28, "inches"), legend.text=element_text(size = 17))
dev.off()
pdf("1_cluster_merge_NoLabel_UMAP.pdf", width=10)
DimPlot(All.merge, reduction = "umap", label = F, pt.size = .1, group.by = "cluster.sort") +
scale_color_hue(labels = levels(All.merge$legend_labels)) +
guides(color = guide_legend(ncol = 1,override.aes = list(size = 4))) +
theme(legend.key.size = unit(0.28, "inches"), legend.text=element_text(size = 17))
dev.off()
##	2
pdf("2_cluster_merge_stim_UMAP.pdf", width=8)
DimPlot(All.merge, reduction = "umap", label = F, pt.size = .1, group.by = "stim") +
guides(color = guide_legend(ncol = 1,override.aes = list(size = 4))) +
theme(legend.key.size = unit(0.28, "inches"), legend.text=element_text(size = 17))
dev.off()
-----------------
##added by youngleel 2021-11-3
  ##
  ###

#缁嗚優绫诲瀷姣斾緥鍙鍖?
  #https://www.data-to-viz.com/graph/barplot.html
getwd()
  file = "D:/ARDS_scripts_1012/ARDS/Step2_harmony_f200_R3/0805/cluster_merge/" 
  
#涓€瀹氳娉ㄦ剰rds閫夋嫨鐨勬纭?
load('D:/ARDS_scripts_1012/ARDS/Step2_harmony_f200_R3/0805/cluster_merge/silicosis_cluster_merge.rds')
all<-All.merge

table( Idents(all), all$orig.ident )
table(Idents(all),all$stim)

# PRO@active.ident
  
  num_tab <- table( Idents(all), all$stim )
num_tab
   
  
freq_tab <- prop.table(x= num_tab , margin=2)
freq_tab
rownames(freq_tab)

barplot(height=freq_tab, width=1, xlim=c(0,5), col=c(1:15), legend= rownames(freq_tab), xlab="")







  



