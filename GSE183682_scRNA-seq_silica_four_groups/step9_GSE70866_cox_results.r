load("G:/r/duqiang_IPF/GSE70866_metainformation_4_platforms/3_ipf_combined_cox_univariate_Adjuste_for_age_sex.RData")
#多个基因批量cox 多因素回归分析
mySurv=with(phe,Surv(time, event))

cox_results <-apply(exprSet , 1 , function(gene){ #把每个基因 分为高低表达两组 进行多因素回归分析
  group=ifelse(gene>median(gene),'high','low')
  survival_dat <- data.frame(group=group,  
                             sex=phe$sex,
                             age=phe$age,stringsAsFactors = F)
  m=coxph(mySurv ~ sex + age + group, #  代码核心句 感兴趣的因素放在最后，这里和deseq2的design正好位置相反
          data =  survival_dat)
  
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)  #提取风险因子
  HRse <- HR * se
  
  #signif(summary(m)$conf.int[,"lower .95"], 5) #lower .95 下限值HR (95% CI for HR)
  #signif(summary(m)$conf.int[,"upper .95"], 5) #lower .95 上限值HR (95% CI for HR)
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se),
                     lower_95=signif(summary(m)$conf.int[,"lower .95"], 5),
                     upper_95=signif(summary(m)$conf.int[,"upper .95"], 5)   ), 3)
  return(tmp['grouplow',])
  
}
)
cox_results=t(cox_results)
head(cox_results)
table(cox_results[,4]<0.05,cox_results[,5]>1)
cox_results[cox_results[,4]<0.05,]
head(cox_results)
head(cox_results[cox_results[,4]<0.05,])
summary(m)

#save(phe,phe_final_3,exprSet,cox_results,file = "G:/r/duqiang_IPF/GSE70866_metainformation_4_platforms/3_ipf_combined_cox_univariate_Adjuste_for_age_sex.RData")
load("G:/r/duqiang_IPF/GSE70866_metainformation_4_platforms/3_ipf_combined_cox_univariate_Adjuste_for_age_sex.RData")

