rm(list=ls())

#install.packages("survival")
 
library(survival)
library(tidyverse)
clincal <- read.table("01rawData/TCGA-LUAD_filtered_clinical.xls",sep = "\t",header = T,check.names = F)

expFile="01rawData/TCGA-LUAD_TPM_mRNA.xls"
expr=read.table(expFile,sep="\t",header=T,check.names=F,row.names = 1)

DEGs_venn="02DEGs/venn_intersectGenes.xls"
DEGs_venn=read.table(DEGs_venn,sep="\t",header=F,check.names=F,row.names = 1)

TumorId <- colnames(expr)[substr(colnames(expr),14,15)!="11" & substr(colnames(expr),14,15)!="02"  ]

exprSet_TEX <- expr[rownames(expr) %in% rownames(DEGs_venn),TumorId]

exprSet_TEX <- as.data.frame( t(exprSet_TEX))
exprSet_TEX$id <- substr(rownames(exprSet_TEX),1,12)
exprSet_TEX_Time <- merge(clincal,exprSet_TEX,by="id")
exprSet_TEX_Time[1:3,1:4]

write.table(exprSet_TEX_Time,file = "03Cox/TCGA-LUAD_TEX_DEGs_exprSetTime.xls",row.names = F,quote = F,sep = "\t")

# Single-factor Cox regression analysis
library(survival)
rt=read.table("03Cox/TCGA-LUAD_TEX_DEGs_exprSetTime.xls",header=T,check.names=F,row.names=1,sep = "\t") 
rt1 <- rt
a= apply(rt[,-(1:9)],2,function(x) log2(x+1))
rt =as.data.frame(cbind(rt[,(1:9)],a))

# Select variable names for single-factor analysis
(variable_names<-colnames(rt)[-c(1:9)])

sur<-Surv(time=rt$Overall_Survival, event = rt$Vital_Status)

pFilter=0.05
outTab=data.frame()
sigGenes = colnames(rt1)[1:9]
for(gene in variable_names){
  if(sd(rt[,gene])<0.001){next}
  if(grepl("-", gene)){next}
  cox=coxph(as.formula(paste0('sur~',gene))  ,data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    group=ifelse(rt[,gene]>median(rt[,gene]),"high","low")
    diff=survdiff(sur ~group,data = rt)
    pValue=1-pchisq(diff$chisq,df=1)
    if(pValue<1){
      sigGenes=c(sigGenes,gene)
      outTab=rbind(outTab,
                   cbind(gene=gene,
                         KM=pValue,
                         HR=coxSummary$conf.int[,"exp(coef)"],
                         HR.95L=coxSummary$conf.int[,"lower .95"],
                         HR.95H=coxSummary$conf.int[,"upper .95"],
                         coxPvalue=coxP) )
    }
  }
}

write.table(outTab,file="03Cox/TCGA-LUAD_uniCox.xls",sep="\t",row.names=F,quote=F)   

surSigExp=rt1[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="03Cox/TCGA-LUAD_uniSigExp.xls",sep="\t",row.names=F,quote=F)

# Function to plot forest plot
bioForest=function(coxFile=null,forestFile=null){
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <-  gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\", rownames(rt)) 
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr," (",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$coxPvalue<0.001, "<0.001", sprintf("%.3f", rt$coxPvalue))

  pdf(file=forestFile, width = 5,height = 7)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(5,2))
  
  xlim = c(0,3)
  par(mar=c(4,3,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'Pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard Ratio',cex=text.cex,font=2,adj=1,)
  
  # forest plot
  par(mar=c(4,0,2,1),mgp=c(2,0.5,0))
  xlim = c(0 ,max(as.numeric(hrLow),as.numeric(hrHigh)+0.5))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard Ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=1.5)
  abline(v=1,col="black",lty=2,lwd=1.3)
  boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.8)
  axis(1)
  dev.off()
}

bioForest(coxFile="03Cox/TCGA-LUAD_uniCox.xls",forestFile="03Cox/TCGA-LUAD_uniForest.pdf")

# Performing a Cox regression model with LASSO
options(stringAsFactors=F)
library(glmnet)

rt=read.table("03Cox/TCGA-LUAD_uniSigExp.xls",header=T,sep="\t",check.names=F,row.names = 1)  
rt1 = rt
rt <- rt[rt$Overall_Survival != 0,]

str(rt)
set.seed(8)

x = as.matrix(rt[,-c(1:9)])
y = Surv(time=rt$Overall_Survival, event = rt$Vital_Status)
cvfit = cv.glmnet(x,y, family = "cox", nfold = 10)
pdf("03Cox/cvfit.pdf",width = 4.6,height = 5)
plot(cvfit)
text(x = log(cvfit$lambda.min),y = 13.48,
     paste('Lambda.min\n',round(cvfit$lambda.min,4)),cex=1,adj=0.9)
text(x = log(cvfit$lambda.1se)+0.01,y = 14.56,
     paste('Lambda.lse\n',round(cvfit$lambda.1se,4)),cex=1,adj=0.9)
dev.off()

fit <- glmnet(x, y, family = "cox",nfold = 10)
library("reshape")
library("ggsci")
library("ggplot2")
x <- coef(fit)  
tmp <- as.data.frame(as.matrix(x))
tmp$coef <- row.names(tmp)
tmp <- reshape::melt(tmp, id = "coef")
tmp$variable <- as.numeric(gsub("s", "", tmp$variable))
tmp$coef <- gsub('_','-',tmp$coef)
tmp$lambda <- fit$lambda[tmp$variable+1] # extract the lambda values
tmp$norm <- apply(abs(x[-1,]), 2, sum)[tmp$variable+1] # compute L1 norm  
tmp$coef2 <- ifelse(tmp$norm==max(tmp$norm),tmp$coef,NA)

ggplot(tmp,aes(log(lambda),value,color = coef)) + 
  geom_vline(xintercept = log(cvfit$lambda.min),size=0.8,color='grey60',alpha=0.8,linetype=2)+
  geom_line(size=0.8) + 
  xlab("Log Lambda") + 
  #xlab("L1 norm")+
  ylab('Coefficients')+
  theme_bw()+ 
  scale_color_manual(values = c(pal_npg()(10),pal_d3()(10),pal_jco()(7)))+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.01))+
  theme(panel.grid = element_blank(),
        # axis.title = element_text(size=15,color='black'),
        # axis.text = element_text(size=12,color='black'),
        legend.key.size = unit(10, "pt"),
        legend.title = element_blank(),
        legend.text = element_text(size=8,color='black'),
        legend.position = 'right')+
  annotate('text',x = -3.5,y=0.03,label='lambda.min = 0.068',color='black')+
  guides(col=guide_legend(ncol = 1))

ggsave("03Cox/fit.pdf",width = 5,height = 4)

myCoefs <- coef(cvfit, s=cvfit$lambda.min)
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
lasso_fea
sigGenes= colnames(rt1)[1:9]
rt <- rt1[,c(sigGenes,lasso_fea)]
rt$riskScore <-  apply(rt[,lasso_fea], 1, function(x) {x %*% myCoefs@x})

# Finding the Best Survival Rate

library(survminer)
rt2 <- rt
# determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(rt2, time = "Overall_Survival", event = "Vital_Status",minprop = 0.3,
                         variables = "riskScore")
rt$risk=as.vector(ifelse(rt$riskScore> res.cut$cutpoint$cutpoint,"High Risk","Low Risk"))

write.table(cbind(id=rownames(rt),rt),
            file="03Cox/TCGA-LUAD_riskScore.xls",sep="\t",quote=F,row.names=F)

lassoGene <-  cbind(Gene = lasso_fea,Coef = myCoefs[which(myCoefs != 0 )])
write.table(lassoGene,file="03Cox/TCGA-LUAD_lasso_coef.xls",row.names = F,quote = F,sep = "\t")

# GSE31210 GSE30219 validset
validSet1 <-  read.table("01rawData/GSE31210_exprSet.xls",check.names = F,sep="\t",quote = "\"",header = T,row.names = 1)
clinical1 <- read.table("01rawData/GSE31210_clinical_clean.xls",check.names = F,sep="\t",quote = "\"",header = T)
validSet1 <-  as.data.frame(t( validSet1[lasso_fea,]))  %>% 
                rownames_to_column("id")
clinical1 <- clinical1[clinical1$tissue== "primary lung tumor",c("id","Overall_Survival", "Vital_Status")]
validSet1_time <- merge(clinical1,validSet1,by="id")

validSet1_time$riskScore <- apply(validSet1_time[,lasso_fea], 1, function(x) {x %*% myCoefs@x})
rt2 <- validSet1_time
rt2$Vital_Status <- ifelse(rt2$Vital_Status=="alive",0,1) 
res.cut <- surv_cutpoint(rt2, time = "Overall_Survival", event = "Vital_Status",minprop = 0.3,
                         variables = "riskScore")
validSet1_time$risk=as.vector(ifelse(validSet1_time$riskScore> res.cut$cutpoint$cutpoint,"High Risk","Low Risk"))

write.table( validSet1_time,
             file="03Cox/GSE31210_validSet_riskScore.xls",
             sep="\t",row.names = F,
             quote=F)

validSet2 <-  read.table("01rawData/GSE30219_exprSet.xls",check.names = F,sep="\t",quote = "\"",header = T,row.names = 1)
clinical2 <- read.table("01rawData/GSE30219_clinical_clean.xls",check.names = F,sep="\t",quote = "\"",header = T)

validSet2 <-  as.data.frame(t( validSet2[lasso_fea,]))  %>% 
  rownames_to_column("id")
clinical2 <- clinical2[clinical2$source_name_ch1=="Lung Tumour" ,c("id","Overall_Survival", "Vital_Status")]
clinical2$Overall_Survival <- as.numeric(clinical2$Overall_Survival)
validSet2_time <- merge(clinical2,validSet2,by="id")

validSet2_time$riskScore <- apply(validSet2_time[,lasso_fea], 1, function(x) {x %*% myCoefs@x})
rt2 <- validSet2_time
rt2$Vital_Status <- ifelse(rt2$Vital_Status=="ALIVE",0,1) 
res.cut <- surv_cutpoint(rt2, time = "Overall_Survival", event = "Vital_Status",minprop = 0.3,
                         variables = "riskScore")
validSet2_time$risk=as.vector(ifelse(validSet2_time$riskScore> res.cut$cutpoint$cutpoint,"High Risk","Low Risk"))

write.table( validSet2_time,
             file="03Cox/GSE30219_validSet_riskScore.xls",
             sep="\t",row.names = F,
             quote=F)
