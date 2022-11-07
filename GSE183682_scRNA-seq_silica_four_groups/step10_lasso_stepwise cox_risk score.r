load("G:/r/duqiang_IPF/GSE70866_metainformation_4_platforms/3_ipf_combined_cox_univariate_Adjuste_for_age_sex.RData")
https://mp.weixin.qq.com/s/ajiNpTDQdfNGY_WgNBVlBQ



#lasso回归筛选变量 得到21个
if(1==1){
  head(cox_results)
  dim(cox_results)
  rownames(cox_results)
  cox_results2=cox_results %>% as.data.frame() %>% filter(p<0.05)
  
  getElement(cox_results,"p")
  cox_results['p']
  head(cox_results2)
  dim(cox_results)
  
  head(exprSet)
  dim(exprSet)
  
  dim(phe)
  head(phe)
  
  identical(colnames(exprSet),rownames(phe))
  
  gene_interested=c("SPP1","LGMN","EMP1","MARCKS","PLA2G7", "MAFB", "BASP1","LHFPL2","FABP5",
                    "CCL3" , "ITGAM" ,"NINJ1", "ENO1","FAM20C","CSTB","PLEKHO1","IER3","CREG1"  
                    ,"PMP22","GADD45B","IFITM2")
  cox_results2=cox_results2[rownames(cox_results2)==gene_interested,]
  cox_results2=cox_results[rownames(cox_results) %in% gene_interested,]
  
  x=exprSet[rownames(exprSet) %in% rownames(cox_results2),]
  x=t(x)
  dim(x)
  y=phe %>%select('time','event')
  head(y)[1:4,1:2]
  head(x)[1:4,1:4]
  
  
  
  table(y$time==0)
  
  
  #OS单位从天转换为年： 是否转换成年不影响结果
  if(1==1){
    y$time <- round(y$time/365,5) #单位年，保留5位小数  time不可以有0
    head(y)
  }
  
str(y$time)
str(y$event) 
  ###################################################lasso#####################
###  #构建模型
library("survival")
library("survminer")
library(glmnet)
  y=data.matrix(Surv(time=y$time,
                     event= y$event))
  head(y)
  head(x)[1:4,1:5]
  fit <- glmnet(x, y, family = 'cox', type.measure = "deviance", nfolds = 10)
  plot(fit,xvar = 'lambda',label = T) #候选DEHGs的lasso系数
  head(coef(fit))
  
  

  #十折交叉检验筛选最佳lambda：
  set.seed(007)
  lasso_fit <- cv.glmnet(x, y, family = 'cox', type.measure = 'deviance', nfolds = 10)
  
  plot(lasso_fit)
  lasso_fit
  head(coef(lasso_fit))
  rownames(lasso_fit$beta)[as.numeric(lasso_fit$beta)>0]
  
  
  #提取最佳λ值(这里选择1se对应lambda)：
  getwd()
  dir.create("G:/r/duqiang_IPF/risk_score_for_mypapers/lse")
  setwd("G:/r/duqiang_IPF/risk_score_for_mypapers/lse")
  if(1==1){
    lambda.1se <- lasso_fit$lambda.1se
    lasso_fit$lambda.min
    lambda.1se  #[1] 0.1925584
    
    
    #使用1se的lambda重新建模：
    model_lasso_1se <- glmnet(x, y, family = 'cox',
                              type.measure = 'deviance', nfolds = 10,
                              lambda = lambda.1se)
    head(model_lasso_1se)
    head(coef(model_lasso_1se))
    
    
    #拎出建模使用基因：
    gene_1se <- rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]#as.numeric后"."会转化为0
    gene_1se #筛选出5个  #"EMP1"   "SPP1"   "IER3"   "LGMN"   "PLA2G7"
    
    
  }
  
  ##########___________________-----------------------------------
  #提取最佳λ值(这里选择min对应lambda)：
  lambda.1se <- lasso_fit$lambda.1se
  lambda.min<-lasso_fit$lambda.min
  lambda.1se 
  lambda.min #[1] 0.06305419
  
  #使用min的lambda重新建模：
  model_lasso_min <- glmnet(x, y, family = 'cox',
                            type.measure = 'deviance', nfolds = 10,
                            lambda = lambda.min)
  head(model_lasso_1se)
  head(coef(model_lasso_1se))
  head(model_lasso_min)
  
  #拎出建模使用基因：
  gene_min <- rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]#as.numeric后"."会转化为0
  gene_min #筛选出10个 "NINJ1"  "EMP1"   "IFITM2" "LHFPL2" "SPP1"   "IER3"   "BASP1"  "MARCKS" "LGMN"   "PLA2G7"
  
  
  #使用1se的lambda重新建模：
  model_lasso_1se <- glmnet(x, y, family = 'cox', type.measure = 'deviance', nfolds = 10,lambda = lambda.1se)
  #拎出建模使用基因：
  gene_1se <- rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]#as.numeric后"."会转化为0
  
  
}
gene_1se #筛选出5个


library(org.Mm.eg.db)

k=keys(org.Mm.eg.db,keytype = "ENSEMBL")
head(k,5)
list=select(org.Mm.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
head(list)
grep(x=list$SYMBOL,"Mt-",value = T)


1. 逐步回归法筛选最优多因素cox回归模型
https://mp.weixin.qq.com/s/lAUGZj3_R874f46dxTXKmg
#相关R包载入：
library(survival)
library(survminer)
#BiocManager::install('My.stepwise')
library(My.stepwise)

#######stepwise cox#############################################
逐步回归得到4个基因
######
if(1==1){
  #表达矩阵保留筛选基因（使用TCGA训练集）：
  exp <- t(exprSet[gene_min,]) #转置，矩阵仅保留lasso-min筛选基因；行为样本，列为基因
  View(exp)
  
  identical(rownames(phe),rownames(exp)) #true
  #将表达矩阵合并到临床信息中：
  dt <- cbind(phe,exp)
  head(dt)
  #逐步回归筛选(拎最后)：
  My.stepwise.coxph(Time = "time",
                    Status = "event",
                    variable.list = gene_min,
                    data = dt)
#  IER3     EMP1     LGMN   LHFPL2 
#  1.181874 1.401688 1.212237 1.395558 
  
  
  #最优cox模型筛选出的4个基因作为最终筛选结果
  final_genes <- c('IER3','EMP1','LGMN','LHFPL2')
  
  
  exp <- as.data.frame(exp[,final_genes])
  head(exp)
  dt <- cbind(phe,exp)
  #构建 最优 多(4)因素 cox比例风险 回归模型：
  colnames(dt)
  train_cox <- coxph(Surv(time, event = event)
                     ~ IER3+EMP1+LGMN+LHFPL2,
                     data = dt)
  train_cox
  #提取回归系数用于不同队列风险评分计算：
  coef <- coef(train_cox)
#     IER3      EMP1      LGMN    LHFPL2 
 # 0.3395320 0.3897767 0.2362561 0.2516435  
}
coef

################risk score###########33

2. 风险评分计算
#首先来看一下大部分文章中计算风险评分所采用的公式（包括本文）：
#Risk score = coefficient1 * 基因1表达量 + ...+ coefficientN * 基因N表达量

######
#2.1.训练集队列风险评分计算：
coef
head(exp)

if(1==1){
  str(exp) #attr(*, "dimnames")=List of 2
  exp=exp %>% as.data.frame()
  #Step1：先将每个基因表达量*对应系数（注意顺序，基因和对应系数不能乱）：
  x <- data.frame(exp$IER3*coef[1],
                  exp$EMP1*coef[2],
                  exp$LGMN*coef[3],
                  exp$LHFPL2*coef[4])
  
  ,
  exp$PLA2G7*coef[5],
                  exp$C10orf107*coef[6],
                  exp$CXCR7*coef[7],
                  exp$SOD3*coef[8],
                  exp$PPP1R14C*coef[9])
  head(x)
  colnames(x) <- names(coef)
  head(x)
  #Step2:将每行相加即为每个样本的风险评分：
  dt$score <- apply(x,1,sum) #相加，并将风险评分列添加到训练集矩阵备用
  #重命名(训练集)：
  train <- dt
  head(train)
  
}


#保存：
getwd()
file="G:/r/duqiang_IPF/risk_score_for_mypapers/lse/gene_interested_stepwise"
#file="G:/r/duqiang_IPF/GSE70866—true—_BAL_IPF_donors_RNA-seq/risk_score"
dir.create(file)
setwd(file)
getwd()
#save(coef,train,file = c('Risk_Score.Rdata')) #此风险评分为9个基因的风险评分
load("G:/r/duqiang_IPF/risk_score_for_mypapers/lse/min_score/Risk_Score.Rdata")







#1. 风险因子联动图绘制
colnames(train)

dt <- train #提取所需列：生存时间、生死、基因表达量、风险评分
dt <- dt[order(dt$score,decreasing = F),] #按风险评分从低到高排序
dt$id <- c(1:length(dt$score)) #根据调整后顺序建立编号id
head(dt)
dt$status <- ifelse(dt$event==0,'alive','death') #0/1转换为字符生死
dt$status <- factor(dt$status,levels = c("death","alive")) #指定因子调整顺序
dt$Risk_Group <- ifelse(dt$score<median(dt$score),'Low Risk','High Risk') #将风险评分按中位数拆分为高/低风险两组
dt$Risk_Group <- factor(dt$Risk_Group,levels = c('Low Risk','High Risk')) #指定顺序
head(dt)

coef
head(dt)
if(1==1){
  exp <- dt[,which(colnames(dt)=='HS3ST1'):
              which(colnames(dt)=="PPP1R14C")] %>%  #as.data.frame.matrix() %>% as.numeric() %>%
    select(everything(),'CXCR7','CDH23') #提取表达矩阵,并调整一下顺序(按风险因子和保护因子排列)
  head(exp)
  
  head(dt)
  #先三个图分别绘制，最后拼接起来；
  library(ggplot2)
  #1.1.风险评分散点图绘制：
  dev.off()
  p1 <- ggplot(dt,aes(x = id,y = score)) +
    xlab(label = 'Case number(ordered by increasing risk score)')+
    ylab(label = 'Risk score')+
    ggtitle(label = 'IPF patients(N=176)')+
    geom_point(aes(col = Risk_Group)) +
    scale_colour_manual(values = c("blue","red")) +
    geom_hline(yintercept = median(dt$score), colour="grey", linetype="dashed", size=0.8) +
    geom_vline(xintercept = sum(dt$Risk_Group == "Low Risk"), colour="grey", linetype = "dashed", size = 0.8) +
    theme_bw()
  p1
  
  
  #1.2.患者生存时间散点图绘制：
  p2 <- ggplot(dt,aes(x = id,y = time)) +
    xlab(label = 'Case number(ordered by increasing risk score)')+
    ylab(label = 'Survival time(days)')+
  ggtitle(label = 'IPF patients(N=176)')+
    geom_point(aes(col = status)) +
    scale_colour_manual(values = c("red","blue")) +
    geom_vline(xintercept = sum(dt$Risk_Group == "Low Risk"), colour = "grey", linetype = "dashed", size = 0.8) +
    theme_bw()
  p2
  
  
  #1.3.表达量热图绘制：
  mycol <- colorRampPalette(c("blue","yellow","red"))(100) #自定义颜色
  exp2 <- t(scale(exp)) #矩阵标准化：
  exp2[1:5,1:4]
  
  
  #添加分组信息：
  head(dt)
  annotation <- data.frame(Type = as.vector(dt[,ncol(dt)]))
  head(annotation)
  rownames(annotation) <- colnames(exp2)
  annotation$Type <- factor(annotation$Type,levels = c('Low Risk','High Risk'))
  head(annotation)
  
  
  ann_colors <- list(Type = c('Low Risk' = "blue",
                              'High Risk' = "red")) #添加分组颜色信息
  #绘图：
  library(pheatmap)
  pheatmap(exp2,
           col= mycol,
           cluster_rows = F,
           cluster_cols = F,
           show_colnames = F,
           annotation_col = annotation,
           annotation_colors = ann_colors,
           annotation_legend = F
  )
  
  
  #将热图转化为ggplot2对象：
  library(ggpubr)
  library(ggplotify )
  p3 <- as.ggplot(as.grob(pheatmap(exp2,
                                   col= mycol,
                                   cluster_rows = F,
                                   cluster_cols = F,
                                   show_colnames = F,
                                   annotation_col = annotation,
                                   annotation_colors = ann_colors,
                                   annotation_legend = F
  )))
  #拼图：
  library(cowplot)
  plot_grid(p1,p2,p3, nrow = 3, align = "v", axis = "tlbr") #可以看到直接拼在一起热图是对不齐的
  
  
  
  
  
}



2. 生存分析绘制    高低风险组的生存分析
#save(coef,train,file = c('06_Risk_Score.Rdata'))
load("G:/r/duqiang_IPF/GSE70866—true—_BAL_IPF_donors_RNA-seq/risk_score/06_Risk_Score.Rdata")


#相关R包载入：
library(survminer)
library(survival)
#重新载入数据：
rm(list = ls()) #先清空工作环境

dt <- train #TCGA训练集
head(dt)
dt$risk <- ifelse(dt$score > median(dt$score),"High","Low") #将风险评分按中位数拆分为高/低风险两组
#2.生存分析绘制：
fit <- survfit(Surv(time, event) ~ risk, data = dt)
fit

p=ggsurvplot(
  fit,
  data = dt,
  censor = T, #是否绘制删失点
  censor.shape = "|", censor.size = 4,
  conf.int = TRUE, #是否显示置信区间
  conf.int.style = "ribbon",
  conf.int.alpha = 0.3,
  pval = TRUE, #是否显示P值
  pval.size = 5,
  legend = "top", #图例位置
  legend.title = 'Risk Score',
  legend.labs = c("High Risk","Low Risk"),
  xlab = "Days",
  ylab = "Survival probablity",
  title = "IER3     EMP1      LGMN   LHFPL2",
  palette = c('#ed0000','#00468b'), #调用色板or自行创建
  ggtheme = theme_bw(), #主题修改
  risk.table = TRUE, #是否风险表添加
  risk.table.col = "strata", #颜色跟随曲线
  risk.table.title = 'Number at risk',
  fontsize = 4,
  risk.table.y.text = FALSE, #是否显示风险表y轴标签
  risk.table.height = 0.2
)
pdf(paste0('risk1', "_surv.pdf"),width = 5, height = 6)
print(p, newpage = FALSE)
dev.off()
getwd()
head(exp)


一、独立预后指标筛选：多因素+单因素cox
##多因素/单因素cox建模与风险森林图绘制

if(1==1){
  #save(coef,train,file = c('06_Risk_Score.Rdata'))
  load("G:/r/duqiang_IPF/GSE70866—true—_BAL_IPF_donors_RNA-seq/risk_score/06_Risk_Score.Rdata")
  
  dt <- train #TCGA总队列
  
  #1.1 数据清洗
  #修改列名：
  #colnames(dt) <- c("submitter_id","status","OS","Age","Gender","Clinical_stage","T_stage","M_stage","N_stage","PHKG1","PLAUR","NDRG1","GAPDH","KIF5A","Risk_score")
  head(dt)
  
  
  #将所有需要观测的变量转换为数值型：
  #年龄转换：
  dt$Age <- as.numeric(dt$age)
  
  #性别转换：
  table(dt$sex) #转换前
  dt$Gender <- ifelse(dt$sex =='0',0,1)
  table(dt$Gender) #转换后
  str(dt$Gender)
  str(dt$sex)
  
  dt$Gender=as.character(dt$Gender)
  #临床分期转换：
  if(1==1){
    table(dt$Clinical_stage) #转换前：初始为罗马数字
    dt$Clinical_stage <- ifelse(dt$Clinical_stage == 'I',1,
                                ifelse(dt$Clinical_stage == 'II',2,
                                       ifelse(dt$Clinical_stage == 'III',3,4)))
    table(dt$Clinical_stage) #转换后：对应数值1234
    
    #T stage转换：
    table(dt$T_stage) #转换前
    dt$T_stage <- ifelse(dt$T_stage == 'T1',1,
                         ifelse(dt$T_stage == 'T2',2,
                                ifelse(dt$T_stage == 'T3',3,4)))
    table(dt$T_stage) #转换后
    
    #M stage转换：
    table(dt$M_stage) #转换前：注意如果存在MX，需转换为NA
    dt$M_stage <- ifelse(dt$M_stage == 'M0',0,
                         ifelse(dt$M_stage == 'M1',1,NA))
    table(dt$M_stage) #转换后
    
  }
  
  class(dt$score) #风险评分不用再转换（已是数值型）
  #保存一下清洗后的数据：
  
  1.2  多因素/单因素cox建模与风险森林图绘制
  #下面我们通过结合单因素与多因素cox，筛选出独立预后指标。
  #多因素cox回归模型构建：
  cox1 <- coxph(Surv(time, event) ~ Age +Gender+ score, data = dt)
  cox1
  summary(cox1)
  if(1==1){
    
    #可以看到年龄、性别、风险评分这三个变量可能是队列患者的独立预后影响因素
    
    #单因素cox回归模型构建：
    #以性别为例：
    cox2 <- coxph(Surv(time, event) ~ score, data = dt)
    cox2
    #通过summary查看模型详细数据(p值、HR值、HR95%CI等)，统计不同队列结果：
    summary(cox2)
    
    论文中将不同队列的单因素、多因素cox结果绘制成了三线表，
    #大家summary统计完数据后Excel绘制即可，这里就不过多赘述，论文表格如下
    
    
    #风险森林图可视化多因素cox：
    options(scipen=1)
    ggforest(cox1,
             data = dt,
             main = "Hazard ratio", #标题
             cpositions = c(0.02, 0.22, 0.4), #前三列距离
             fontsize = 1, #字体相对大小
             refLabel = "1", #相对变量的数值标签
             noDigits = 2) #95%CI、p值、HR值等保留小数点后位数
    
    
    
  }
  
  
  
}


二、timeROC验证预后模型（风险评分）准确度
#save(coef,train,file = c('06_Risk_Score.Rdata'))
load("G:/r/duqiang_IPF/GSE70866—true—_BAL_IPF_donors_RNA-seq/risk_score/06_Risk_Score.Rdata")
head(train)

if(1==1){
  #提取目标列：
  dt <- train[,c(1:4,ncol(train))] #event time sex age    score
  head(dt)
  
  library(timeROC)
  #时间依赖性ROC曲线绘制：
  boxplot(dt$time)
  timeROC <- with(dt, #with限制数据框;
                  timeROC(T = time,
                          delta = event,
                          marker = score, #预测的生死（默认较大的预测值和高事件风险相关，如果较大预测值和低风险事件相关，需要添加一个“-”反转关联）
                          cause = 1, #阳性事件结局（这里阳性事件为死亡，对应1）
                          times = c(365,2*365,3*365), #时间点划分：c(1,3,5)1年、3年、5年
                          ROC = TRUE, #是否保存TP和FP值
                          iid = TRUE, #是否计算ROC曲线下面积
                          weighting = "marginal")) #权重计算方法，默认
  print(timeROC)
 #获取p值
   se <- sqrt(var())  # 获取AUC的SE（标准误）
  b <- r1$auc - .5
  z <- (b / se)  # 计算Z值
  2 * pt(-abs(z), df=Inf)  ## two-sided test  
 
   library(dplyr)
  #BiocManager::install('pROC')
  library(pROC)
  ci.auc(timeROC) 
  #使用ggplot2绘图：
  #提取各生存时间点的TPR、FPR值：
  df <- data.frame(FPR = as.numeric(timeROC$FP),
                   TPR = as.numeric(timeROC$TP),
                   time = rep(as.factor(c(1,2,3)), each = nrow(timeROC$TP)))
  head(df)
  #自定义主题：
  mytheme <- theme(axis.title = element_text(size = 13),
                   axis.text = element_text(size = 11),
                   plot.title = element_text(size = 14,
                                             hjust = 0.5,
                                             face = "bold"),
                   legend.title = element_text(size = 13),
                   legend.text = element_text(size = 11),
                   legend.background = element_rect(linetype = 1, size = 0.25, colour = "black"),
                   legend.position = c(1, 0),
                   legend.justification = c(1, 0))
  
  #timeROC绘制：
  p <- ggplot() +
    geom_line(data = df,aes(x = FPR, y = TPR, color = time), size = 1) +
    geom_line(aes(x = c(0,1),y = c(0,1)),color = "grey") +
    scale_color_manual(name = NULL, values = c("#8AD879","#FA9F42", "#F3533A"),
                       labels = c("1 Years(AUC = 82.25)" ,
                                  "2 Years(AUC = 73.02)",
                                  "3 Years(AUC = 89.71)"))+
    theme_bw() +
    mytheme +
    labs(x = "1 - Specificity (FPR)",
         y = "Sensitivity (TPR)")
  p
  
  
}








