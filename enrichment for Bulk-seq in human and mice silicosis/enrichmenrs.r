#################FIGURE S1A#####################################HUMAN
#########################################################################################################################################################

getwd()
path="G:/silicosis/sicosis/yll/clusterprolifer-frompapers-10 silicosis patients and 7 human donors/thesameparameter_as_paper"
dir.create(path)
setwd(path)

mymarkers=openxlsx:: read.xlsx("G:/silicosis/sicosis/yll/clusterprolifer-frompapers-10 silicosis patients and 7 human donors/papers -human 10 silicosis and 7 donors.xlsx",
                    colNames = TRUE)
DEG=mymarkers
#去重复值 方法1.随机去重  ids = ids[!duplicated(ids$symbol),]#谁排第一个就取谁
table(!duplicated(DEG$Symbol))
DEG=DEG[!duplicated(DEG$Symbol),]
table(!duplicated(DEG$Symbol))
head(DEG)

DEG$p_val=as.numeric(DEG$`p-val`)
DEG$avg_log2FC=as.numeric(DEG$"log2FoldChange")
rownames(DEG)=DEG$"Symbol"    ###################3此处行名有重复
DEG$gene=DEG$Symbol

head(DEG)
dim(DEG)
#https://blog.csdn.net/weixin_33915154/article/details/112069336
DEG$change<-ifelse(as.numeric(DEG$p_val)<1 & abs(as.numeric(DEG$avg_log2FC))>=1.5,  #自己改变筛选条件
                   ifelse(as.numeric(DEG$avg_log2FC)>=1.5,'UP','DOWN'),'NOT')

head(DEG)
table(DEG$"avg_log2FC">0) ;table(DEG$change)
getwd()
#NOT  UP 
#743 426
library("org.Hs.eg.db")
keytypes(org.Hs.eg.db)
df <- bitr(DEG$gene, fromType = "SYMBOL", ###注意输入的基因类型是symbol还是entrezid
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db) #########注意人还是小鼠
head(df)
DEG$SYMBOL =DEG$gene              ###注意输入的基因类型是symbol还是entrezid
DEG=merge(DEG,df,by='SYMBOL')
head(DEG)

gene_up= DEG[DEG$change == 'UP','ENTREZID']  #[1] 417
gene_down=DEG[DEG$change == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] )
#自建函数   输入go富集条目  和文件的名称
clusterProfilerOut = function(enrichRes,outfile){
  enrichRes$function.gene.num = unlist(strsplit(enrichRes[,"BgRatio"],"/"))[seq(1,length(unlist(strsplit(enrichRes[,"BgRatio"],"/"))),2)]
  enrichRes$GeneRatio = enrichRes$Count / as.numeric(enrichRes$function.gene.num) #计算generatio
  write.xlsx(enrichRes,paste0(outfile, ".xlsx"),col.names=T,row.names=F)
  
  if (dim(enrichRes)[1]>=10) {enrichRes.use = enrichRes[1:10,]} else{enrichRes.use = enrichRes[1:dim(enrichRes)[1],]}#显示的富集条目不超过10条
  enrichRes.use = enrichRes.use[order(enrichRes.use$GeneRatio, decreasing=T),] #按照generatio降序显示
  enrichRes.use$Description = factor(enrichRes.use$Description,levels=rev(enrichRes.use$Description))
  gap = (max(enrichRes.use$GeneRatio) - min(enrichRes.use$GeneRatio))/5
  
  #制作并输出pdf
  pdf(paste0(outfile, ".pdf"))
  p = ggplot(enrichRes.use,mapping = aes(x=GeneRatio,y=Description))+
    geom_point(aes(size=Count,color=p.adjust))+theme_bw()+ 
    scale_color_gradient(low = "blue", high = "red")+
    scale_x_continuous(expand = c(gap, gap))+ylab(NULL)+ 
    scale_y_discrete(labels=function(x) str_wrap(x, width=60))
  print(p)
  dev.off()
}  
### GO database analysis 

g_list=list(gene_up=gene_up,
            gene_down=gene_down,
            gene_diff=gene_diff)
length(g_list)
for (each in names(g_list)) {
  print(each)
}

if(T){
  go_enrich_results <- lapply( g_list , function(gene) {
    
    lapply( c('BP','MF','CC') , function(ont) {
      cat(paste('Now process ',ont ))
      ego <- enrichGO(gene          = gene,
                      universe      = gene_all,
                      OrgDb         = org.Hs.eg.db,  #包含人注释信息的数据库OrgDb   = org.Mm.eg.db,
                      ont           = ont ,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.99,
                      qvalueCutoff  = 0.99,
                      readable      = TRUE)
      print( head(ego) )
      return(ego)
    })
  })
  save(go_enrich_results,file = 'go_enrich_results.Rdata')
  
}



################FIGURE S1B#######################################MICE
#########################################################################################################################################################
getwd()
path="G:/silicosis/sicosis/yll/clusterprolifer-frompapers-20 mouse silicosis model and 20 silicosis sham/60d-duplication"
dir.create(path)
setwd(path)

mymarkers=openxlsx:: read.xlsx("G:/silicosis/sicosis/yll/clusterprolifer-frompapers-20 mouse silicosis model and 20 silicosis sham/clusterprolifer-frompapers-20 mouse silicosis model and 20 silicosis sham.xlsx",
                   colNames = TRUE)

head(mymarkers)
colnames(mymarkers)=mymarkers[1,]
mymarkers=mymarkers[-1,]
mymarkers$pval=0.00000001
head(mymarkers)
#View(markers)
DEG<-mymarkers
#去重复值 方法1.随机去重  ids = ids[!duplicated(ids$symbol),]#谁排第一个就取谁
table(!duplicated(DEG$`Gene Name`))
DEG=DEG[!duplicated(DEG$`Gene Name`),]
DEG$p_val=DEG$pval
DEG$avg_log2FC=DEG$"Fold Change"
rownames(DEG)=DEG$"Gene Name"    ###################3此处行名有重复
DEG$gene=DEG$`Gene Name`
head(DEG)
dim(DEG)#[1] 818  11
#https://blog.csdn.net/weixin_33915154/article/details/112069336
DEG$change<-ifelse(DEG$p_val<0.05 & abs(as.numeric(DEG$avg_log2FC))>=1,  #自己改变筛选条件
                   ifelse(DEG$avg_log2FC>=1,'UP','DOWN'),'NOT')

head(DEG)
dim(DEG)#[1] 818  12
#rownames(DEG)<-DEG[,1]
table(DEG$change)
getwd()

library(org.Mm.eg.db)
keytypes(org.Mm.eg.db)
library(clusterProfiler)
df <- bitr(DEG$gene, fromType = "SYMBOL", ###注意输入的基因类型是symbol还是entrezid
           toType = c( "ENTREZID"),
           OrgDb = org.Mm.eg.db)
head(df)
DEG$SYMBOL =DEG$gene              ###注意输入的基因类型是symbol还是entrezid
DEG=merge(DEG,df,by='SYMBOL')
head(DEG)


gene_up= DEG[DEG$change == 'UP','ENTREZID'] 
length(gene_up)#[1] 718
gene_down=DEG[DEG$change == 'DOWN','ENTREZID'] 
length(gene_down) #[1] 64
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] )





### GO database analysis 

g_list=list(gene_up=gene_up,
            gene_down=gene_down,
            gene_diff=gene_diff)
length(g_list)
for (each in names(g_list)) {
  print(each)
}

if(T){
  go_enrich_results <- lapply( g_list , function(gene) {
    
    lapply( c('BP','MF','CC') , function(ont) {
      cat(paste('Now process ',ont ))
      ego <- enrichGO(gene          = gene,
                      universe      = gene_all,
                      OrgDb         = org.Mm.eg.db,
                      ont           = ont ,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.99,
                      qvalueCutoff  = 0.99,
                      readable      = TRUE)
      print( head(ego) )
      return(ego)
    })
  })
  save(go_enrich_results,file = 'go_enrich_results.Rdata')
  
}


#load(file = 'go_enrich_results.Rdata')

n1= c('gene_up','gene_down','gene_diff')
n2= c('BP','MF','CC') 
for (i in 1:3){
  #i=1
  for (j in 1:3){
    #j=1
    fn=paste0('dotplot_',n1[i],'_',n2[j],'.png')
    cat(paste0(fn,'\n'))
    png(fn,res=150,width = 1080,height = 1080)
    print( dotplot(go_enrich_results[[i]][[j]] ))
    dev.off()
    
    ego = go_enrich_results[[i]][[j]]
    res = ego[ego$p.adjust<0.05,]
    print(dim(res))
    print(table(res$ONTOLOGY))
    if(nrow(res) != 0){#go，使用自定义的函数画图 ,注意命名时候把i的 文件属性名去掉
      clusterProfilerOut(res, fn)
    }
    
  }
}


© 2022 GitHub, Inc.
