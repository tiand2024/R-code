rm(list=ls())
## Load R packages
library(tidyverse)

datExpr <- read.table("01rawData/TCGA-LUAD_TPM_mRNA.xls",sep="\t",header=T,check.names=F,row.names = 1)
TumorId <- colnames(datExpr)[substr(colnames(datExpr),14,15)!="11" & substr(colnames(datExpr),14,15)!="02"  ]
datExpr <- datExpr[,TumorId]
riskScore <- read.table("03Cox/TCGA-LUAD_riskScore.xls",header = T,sep = "\t")
low <- riskScore$id[ riskScore$risk== "Low Risk"]
high <- riskScore$id[ riskScore$risk== "High Risk"]
lowId <- TumorId[substr(TumorId,1,12) %in% low]
highId <- TumorId[substr(TumorId,1,12) %in% high]
# 343 low, 164 high samples
GSEA_df <- datExpr[,c(lowId,highId)]

GSEA_df[1:3,1:5]
GSEA_df <- cbind(Name = rownames(GSEA_df) ,DESCRIPTION = "na",GSEA_df)

group <- c(rep("low", 343), rep("high", 164))
group <- paste(group, collapse = " ")
group <- c(paste(c(507, 2, 1), collapse = " "), "# low high", group)

write.table(file = "06GSEA/TCGA-LUAD_risk_exp.txt", GSEA_df, sep = "\t", row.names = F, quote = F)
write.table(file = "06GSEA/risk_group.cls", group, col.names = F, row.names = F, quote = F)



GSEA_df <- datExpr[,c(lowId,highId)]
# 343low   164high 
GSEA_df[1:3,1:5]
GSEA_df <- GSEA_df %>% rownames_to_column("id")

write.table(file = "07TIDE/TCGA-LUAD_risk_exp.xls", GSEA_df, sep = "\t", row.names = F, quote = F)


expr <-  read.table("07TIDE/TCGA-LUAD_risk_exp.xls",header = T,sep = "\t",check.names = F,row.names = 1)
xx <- substr(colnames(expr),14,15)
table(xx)

genes <-  read.table("07TIDE/HLA_genes.txt",header = T,sep = "\t",check.names = F)

HLA_expr <- expr[rownames(expr)%in%genes$Name, c(lowId,highId)]

riskScore <- read.table("03Cox/TCGA-LUAD_riskScore.xls",header = T,sep = "\t")

library(tibble)
# HLA_expr =read.table("GSVA/GSVA_geneSet_result.txt",header=T,sep="\t",check.names=F,row.names = 1)
boxdata <- as.data.frame(t( (HLA_expr) )) %>% rownames_to_column("id")
boxdata$id <- substr(boxdata$id,1,12)
boxdata <- merge(boxdata,riskScore[,c("id","risk")],by="id")

write.table(boxdata, "07TIDE/HLA_genes_exprSet.xls", quote = F, row.names = F,sep = "\t")




boxdata <- as.data.frame(t( log2(HLA_expr+1) )) %>% rownames_to_column("id")

boxdata$id <- substr(boxdata$id,1,12)
boxdata <- merge(boxdata,riskScore[,c("id","risk")],by="id")


boxdata1 <- boxdata[,-1]
library(tidyr)
library(ggpubr)

boxdata1 <- gather(boxdata1,cellType, value,colnames(boxdata1)[-ncol(boxdata1)])
boxdata1$risk <- factor(boxdata1$risk,levels = c("Low Risk","High Risk"))
boxdata1$cellType<- factor(boxdata1$cellType,levels = genes$Name)
p2 <- ggboxplot(boxdata1, x = "cellType", y = "value",
                color = "risk", palette = c("#028846","red"),
                add = "jitter",
                legend.title = "",
                add.params = list(alpha=0.5),)+ rotate_x_text(45) +xlab("")+ylab("Gene expression log2(TPM+1)")
p2 + stat_compare_means(aes(group = risk),method="wilcox.test", label = "p.signif")
ggsave("07TIDE/HLA_genes_type_boxplot.pdf",height=4.2,width=9)

# Immune checkpoint genes

genes <-  read.table("07TIDE/immune-checkpoint_genes.txt",header = T,sep = "\t",check.names = F)
checkpoint_expr <- expr[rownames(expr)%in%genes$gene, c(lowId,highId)]

riskScore <- read.table("03Cox/TCGA-LUAD_riskScore.xls",header = T,sep = "\t")

boxdata <- as.data.frame(t( (checkpoint_expr) )) %>% rownames_to_column("id")
boxdata$id <- substr(boxdata$id,1,12)
boxdata <- merge(boxdata,riskScore[,c("id","risk")],by="id")

write.table(boxdata, "07TIDE/immune-checkpoint_genes_exprSet.xls", quote = F, row.names = F,sep = "\t")




boxdata <- as.data.frame(t( log2(checkpoint_expr+1) )) %>% rownames_to_column("id")

boxdata$id <- substr(boxdata$id,1,12)
boxdata <- merge(boxdata,riskScore[,c("id","risk")],by="id")


boxdata1 <- boxdata[,-1]

boxdata1 <- gather(boxdata1,cellType, value,colnames(boxdata1)[-ncol(boxdata1)])
boxdata1$risk <- factor(boxdata1$risk,levels = c("Low Risk","High Risk"))
boxdata1$cellType<- factor(boxdata1$cellType,levels = genes$Symbol)
p2 <- ggboxplot(boxdata1, x = "cellType", y = "value",
                color = "risk", palette = c("#028846","red"),
                add = "jitter",
                legend.title = "",
                add.params = list(alpha=0.5),)+ rotate_x_text(45) +xlab("")+ylab("Gene expression log2(TPM+1)")
p2 + stat_compare_means(aes(group = risk),method="wilcox.test", label = "p.signif")
ggsave("07TIDE/immune-checkpoint_genes_type_boxplot.pdf",height=4.2,width=13)
