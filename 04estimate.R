# library(utils)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos=rforge, dependencies=TRUE)

rm(list=ls())

library(estimate)

inputFile="01rawData/TCGA-LUAD_TPM_mRNA.xls" 
# estimate
filterCommonGenes(input.f= inputFile,
                  output.f="05immune/TCGA-LUAD_commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "05immune/TCGA-LUAD_commonGenes.gct",
              output.ds="05immune/TCGA-LUAD_estimateScore.gct", 
              platform="illumina")


scores=read.table("05immune/TCGA-LUAD_estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
scores <- scores[substr(rownames(scores),14,15) !=11,]
out=rbind(id=colnames(scores),scores)
write.table(out,file="05immune/TCGA-LUAD_estimateScores.xls",sep="\t",quote=F,col.names=F)


# Plotting boxplot of immune scores in high and low-risk groups
library(tidyr)
library(ggpubr)
library(cowplot)
mycol <- c("#028846","red")

risk <- read.table("03Cox/TCGA-LUAD_riskScore.xls",header = T,sep = "\t",check.names = F)
scores =read.table("05immune/TCGA-LUAD_estimateScores.xls",header=T,sep="\t",check.names=F)
scores$id <- substr(scores$id ,1,12)
Data <- merge(risk,scores,by="id")
write.table(Data,file="05immune/TCGA-LUAD_estimateScores_risk.xls",sep="\t",quote=F,col.names=T,row.names = F)

plotData <- Data[,c("risk","StromalScore","ImmuneScore","ESTIMATEScore")]
plotData$risk <- factor(plotData$risk,levels = c("Low Risk" ,"High Risk"))

p1 <- ggplot(plotData,aes(risk,StromalScore,fill=risk))+
  geom_boxplot(outlier.colour = NA,notch = F,size = 0.3,width=.3,alpha = 0.4)+
  geom_jitter(shape = 21,size=2,width = 0.2,alpha = 0.4)+
  geom_violin(position = position_dodge(width = .75), 
              size = 0.3,alpha = 0.4,trim = T)+
  
  theme_cowplot()+
  theme(legend.position = 'none',
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=10),
        axis.title.x = element_blank())+
  scale_fill_manual(values=mycol[c(1,2)])+
  stat_compare_means(comparisons = list(c("Low Risk" ,"High Risk")),
                      label = 'p.signif')#+
  #stat_compare_means(label.y = max(plotData$StromalScore)+5.5)
p2 <-  ggplot(plotData,aes(risk,ImmuneScore ,fill=risk))+
  geom_boxplot(outlier.colour = NA,notch = F,size = 0.3,width=.3,alpha = 0.4)+
  geom_jitter(shape = 21,size=2,width = 0.2,alpha = 0.4)+
  geom_violin(position = position_dodge(width = .75), 
              size = 0.3,alpha = 0.4,trim = T)+
  theme_cowplot()+
  theme(legend.position = 'none',
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=10),
        axis.title.x = element_blank())+
  scale_fill_manual(values=mycol[c(1,2)])+
  stat_compare_means(comparisons = list(c("Low Risk" ,"High Risk")),
                     label = 'p.signif')#+
  #stat_compare_means(label.y = max(plotData$ImmuneScore )+5.5)
p3 <- ggplot(plotData,aes(risk,ESTIMATEScore ,fill=risk))+
  geom_boxplot(outlier.colour = NA,notch = F,size = 0.3,width=.3,alpha = 0.4)+
  geom_jitter(shape = 21,size=2,width = 0.2,alpha = 0.4)+
  geom_violin(position = position_dodge(width = .75), 
              size = 0.3,alpha = 0.4,trim = T)+
  theme_cowplot()+
  theme(legend.position = 'none',
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=10),
        axis.title.x = element_blank())+
  scale_fill_manual(values=mycol[c(1,2)])+
  stat_compare_means(comparisons = list(c("Low Risk" ,"High Risk")),
                     label = 'p.signif')#+
ggarrange(p1,p2,p3,ncol = 3, nrow = 1)

ggsave("05immune/TCGA-LUAD_estimateScore_violin.pdf",width = 9,height = 3)

