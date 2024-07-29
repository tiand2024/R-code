rm(list=ls())
# Load R packages
library(limma)
library(tidyverse)
source("D:/GEO/script/Filter_TCGA_Replicate_Samples.R")
# Extract clinical information
clinical <- read.table("01rawData/TCGA-LUAD_clinical.xls",check.names = F,sep="\t",quote = "\"",header = T)
clinical$id <- clinical$submitter_id
clinical <- clinical[!is.na( clinical$vital_status ),]
clinical$Overall_Survival <-  ifelse(clinical$vital_status =="Alive",clinical$days_to_last_follow_up,clinical$days_to_death)
clinical  <- clinical[!is.na( clinical$Overall_Surviva),]
clinical$Vital_Status <- ifelse(clinical$vital_status =="Alive",0,1)

clinical$Age <- clinical$age_at_index
clinical$Gender <- clinical$gender
clinical$pathologic_stage <-  gsub("[ABCabc]","", clinical$ajcc_pathologic_stage)
clinical$pathologic_t <-  gsub("[ABCabc]","", clinical$ajcc_pathologic_t)
clinical$pathologic_m <-  gsub("[ABCabc]","", clinical$ajcc_pathologic_m)
clinical$pathologic_n <-  gsub("[ABCabc]","", clinical$ajcc_pathologic_n)

clinical$cigarettes_per_day <- clinical$cigarettes_per_day
clincal <- clinical[c("id","Overall_Survival","Vital_Status","Age","Gender","cigarettes_per_day","pathologic_stage","pathologic_t","pathologic_m","pathologic_n")]

clincal <- unique(clincal)

# Read counts data
expMatrix <- readRDS(file = "01rawData/TCGA-LUAD.rds")
expMatrix <- expMatrix %>% column_to_rownames("gene_id")

# Ensure that the row names of the expression matrix match the names in the vector storing gene lengths. 
gene_length <-read.table("01rawData/gene_length.xls", header = T,sep = "\t")
rownames(gene_length) <- gene_length$gene_id 
feature_ids <- rownames(expMatrix)

if (! all(feature_ids %in% rownames(gene_length))){
  tbl <- table(feature_ids %in% rownames(gene_length))
  msg1 <- sprintf("%i gene is shared, %i gene is specified", tbl[[2]],tbl[[1]])
  warning(msg1)
  
} 

if (! identical(feature_ids, rownames(gene_length))){
  msg2 <- sprintf("Given GTF file only contain %i gene, but experssion matrix has %i gene", nrow(gene_length), nrow(expMatrix))
  warning(msg2)
}

# trim the expression matrix and effetive gene length

mm <- match(rownames(expMatrix), rownames(gene_length))
gene_length <- gene_length[mm, ]

if (identical(rownames(gene_length), rownames(expMatrix))){
  print("GTF and expression matix now have the same gene and gene in same order")
}

x <- expMatrix / gene_length$length
expMatrix_tpm <- t( t(x) / colSums(x) ) * 1e6 
expMatrix_tpm <-  as.data.frame( round(expMatrix_tpm,4))
expMatrix_tpm[1:5,1:7]


tcgaReplicateFilterId <- tcgaReplicateFilter(colnames(expMatrix_tpm),analyte_target="RNA")
expMatrix_tpm <- expMatrix_tpm[,tcgaReplicateFilterId]

xx <- substr(colnames(expMatrix_tpm),14,15)
table(xx)
expMatrix_tpm[1:3,1:3]


NormalId <- colnames(expMatrix_tpm)[substr(colnames(expMatrix_tpm),14,15) =="11"]

TumorId <- colnames(expMatrix_tpm)[substr(colnames(expMatrix_tpm),14,15)!="11" & substr(colnames(expMatrix_tpm),14,15)!="02"  ]

duplicated_sample <- colnames(expMatrix_tpm)[ substr(colnames(expMatrix_tpm),14,15)=="02"  ]
TumorId[substr(TumorId,1,12) %in% substr(duplicated_sample,1,12) ]

coTumorId <- TumorId[substr(TumorId,1,12) %in% clincal$id]
# coTumorId <- TumorId[TumorId %in% clincal$id]
# 507 tumor cases, 59 normal cases
str(clincal)
clincal <- clincal[clincal$id %in%  substr(coTumorId,1,12),]

write.table(clincal,file = "01rawData/TCGA-LUAD_filtered_clinical.xls",row.names = F,quote = F,sep = "\t")

STAR_TPM <- expMatrix_tpm[,c(NormalId,coTumorId)]
STAR_TPM[1:3,1:5]

# Save computed TPM to local file
STAR_TPM <- STAR_TPM %>% rownames_to_column("gene_id")
write.table(STAR_TPM, "01rawData/TCGA-LUAD_TPM.xls", sep="\t", quote=F, row.names=F)


load("D:/TCGA/gencode.v36.annotation.Rdata")
gtf_df[1:3,1:3]
gtf_df <- gtf_df[gtf_df$"gene_type" == "protein_coding",]


exprSet <- STAR_TPM %>%
  dplyr::inner_join(gtf_df,by="gene_id") %>%
  # Remove redundant information
  dplyr::select(-c(gene_id,gene_type))%>%
  # Reorder or rearrange
  dplyr::select(gene_name,everything()) %>%
  # Calculate average expression values for duplicate genes
  dplyr::group_by(gene_name) %>%
  dplyr::summarise_all(mean)
exprSet[1:3,1:4]

write.table(exprSet,file = "01rawData/TCGA-LUAD_TPM_mRNA.xls",row.names = F,quote = F,sep = "\t")

exprSet1 <- exprSet %>% column_to_rownames("gene_name")
group_list=c(rep("Normal",length(NormalId)),
             rep("Tumor",length(coTumorId)))  # Grouping information
group_list=factor(group_list)
# Forcefully enforce order
group_list <- relevel(group_list, ref="Normal")
table(group_list)
# Remove genes with low expression
exprSet1=exprSet1[rowMeans(exprSet1)>0.001,]
exprSet2 <- log2(exprSet1+1)

# Differential analysis

design=model.matrix(~ group_list)
fit=lmFit(exprSet2,design)
fit=eBayes(fit)
res=topTable(fit,adjust='fdr',coef="group_listTumor",number=Inf)
allDiff <- na.omit(res)

logFCCutoff <- 1
pvalueCutoff <- 0.05

outDiff=allDiff[(abs(allDiff$logFC)>logFCCutoff & allDiff$adj.P.Val<pvalueCutoff),]
outDiff <- outDiff %>% rownames_to_column(var = "id")
write.table(outDiff,file="02DEGs/TCGA-LUAD_mRNA_limma_diff.xls",row.names=F,quote=F,sep = "\t")

allDiff <- allDiff  %>% rownames_to_column(var = "id")
write.table(allDiff, '02DEGs/TCGA-LUAD_mRNA_limma_alldiff.xls', sep = '\t', row.names=F, quote = FALSE)


# Plot heatmap of differential genes
library(pheatmap)
geneNum=50
outDiff=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
# diffGeneName=rownames( outDiff )
diffGeneName=outDiff$id 
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(geneNum*2)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp= exprSet1[hmGene,]

hmExp = log2(hmExp+1)

max(hmExp)
min(hmExp)

# Normal  Tumor 
# 59    507 
Type= factor(c(rep("Normal",59),rep("Tumor",507)),levels = c("Normal","Tumor"))
names(Type)=colnames(hmExp)
Type=as.data.frame(Type)

pdf(file="02DEGs/mRNA_heatmap.pdf",height=8,width=9)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8,
         border=FALSE)
dev.off()



library(ggplot2)
library(ggrepel)
res <- allDiff

Significant=ifelse((res$adj.P.Val< 0.05 & abs(res$logFC)> 1), ifelse(res$logFC > 1,"Up","Down"), "Not")
# Plot volcano plot
p = ggplot(res, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#004BFB", "#BFBFBF", "#F91F10"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-1,1), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
  xlab(expression("log"["2"]*"FC"))+
  ylab(expression("-log"["10"]*"adj.P.Val"))+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()

# Save plot as PDF file
pdf("02DEGs/mRNA_vol.pdf",width=5.5,height=5)
print(p)
dev.off()

# Perform lncRNA differential expression analysis using limma
load("D:/TCGA/gencode.v36.annotation.Rdata")
gtf_df[1:3,1:3]
gtf_df <- gtf_df[gtf_df$"gene_type" == "lncRNA",]

exprSet_lncRNA <- STAR_TPM %>%
  dplyr::inner_join(gtf_df,by="gene_id") %>%
  dplyr::select(-c(gene_id,gene_type))%>%
  dplyr::select(gene_name,everything()) %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise_all(mean)
exprSet_lncRNA[1:3,1:4]

write.table(exprSet_lncRNA,file = "01rawData/TCGA-LUAD_TPM_lncRNA.xls",row.names = F,quote = F,sep = "\t")

exprSet_lncRNA1 <- exprSet_lncRNA %>% column_to_rownames("gene_name")
group_list=c(rep("Normal",length(NormalId)),
             rep("Tumor",length(coTumorId))) 
group_list=factor(group_list)
group_list <- relevel(group_list, ref="Normal")
table(group_list)
# group_list
# Normal  Tumor 
# 59    507 
exprSet_lncRNA1=exprSet_lncRNA1[rowMeans(exprSet_lncRNA1)>0.001,]
exprSet_lncRNA2 <- log2(exprSet_lncRNA1+1)

design=model.matrix(~ group_list)
fit=lmFit(exprSet_lncRNA2,design)
fit=eBayes(fit)
res_lncRNA=topTable(fit,adjust='fdr',coef="group_listTumor",number=Inf)
allDiff_lncRNA <- na.omit(res_lncRNA)

logFCCutoff <- 1
pvalueCutoff <- 0.05

outDiff_lncRNA=allDiff_lncRNA[(abs(allDiff_lncRNA$logFC)>logFCCutoff & allDiff_lncRNA$adj.P.Val<pvalueCutoff),]
outDiff_lncRNA <- outDiff_lncRNA %>% rownames_to_column(var = "id")
write.table(outDiff_lncRNA,file="02DEGs/TCGA-LUAD_lncRNA_limma_diff.xls",row.names=F,quote=F,sep = "\t")

allDiff_lncRNA <- allDiff_lncRNA  %>% rownames_to_column(var = "id")
write.table(allDiff_lncRNA, '02DEGs/TCGA-LUAD_lncRNA_limma_alldiff.xls', sep = '\t', row.names=F, quote = FALSE)


# Plot heatmap of differentially expressed genes
library(pheatmap)
geneNum=50
outDiff_lncRNA=outDiff_lncRNA[order(as.numeric(as.vector(outDiff_lncRNA$logFC))),]

diffGeneName=outDiff_lncRNA$id 
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(geneNum*2)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp= exprSet_lncRNA1[hmGene,]

hmExp = log2(hmExp+1)

max(hmExp)
min(hmExp)

# Normal  Tumor 
# 59    507 
Type= factor(c(rep("Normal",59),rep("Tumor",507)),levels = c("Normal","Tumor"))
names(Type)=colnames(hmExp)
Type=as.data.frame(Type)

pdf(file="02DEGs/lncRNA_heatmap.pdf",height=8,width=9)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8,
         border=FALSE)
dev.off()


library(ggplot2)

Significant=ifelse((res$adj.P.Val< 0.05 & abs(res$logFC)> 1), ifelse(res$logFC > 1,"Up","Down"), "Not")
# Plot volcano plot
p = ggplot(res, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#004BFB", "#BFBFBF", "#F91F10"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-1,1), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
  xlab(expression("log"["2"]*"FC"))+
  ylab(expression("-log"["10"]*"adj.P.Val"))+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()

pdf("02DEGs/lncRNA_vol.pdf",width=5.5,height=5)
print(p)
dev.off()



# Use ggvenn to draw a Venn diagram

outDiff <- read.table("02DEGs/TCGA-LUAD_mRNA_limma_diff.xls",check.names = F,sep="\t",header = T)

TEX_genes <- read.table("01rawData/TEX-associated genes_GeneCards.csv",check.names = F,sep=",",header = T)
TEX_genes <- TEX_genes[TEX_genes$`Relevance score` >30,]
x <- list(DEGs=outDiff$id,TEXGs=TEX_genes$`Gene Symbol`)
library(ggvenn)
mycol <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#ADD1E5",
           '#CCEBC5', '#FCCDE5', '#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")
ggvenn(x,c("DEGs","TEXGs"),
       stroke_size = 0.3,
       show_percentage = F,
       fill_color =mycol[4:8])
ggsave("02DEGs/venn.pdf",width = 4.2,height = 4.5)

intersectGenes<- intersect(outDiff$id,TEX_genes$`Gene Symbol`)
write.table(intersectGenes,file = "02DEGs/venn_intersectGenes.xls",row.names = F,col.names = F,quote = F,sep = "\t",)


# miRNA differently expression analysis
expMatrix_miRNA <- readRDS(file = "01rawData/TCGA-LUAD_miRNA_matrue_RPM.rds" )
expMatrix_miRNA <-expMatrix_miRNA[,-1] %>% column_to_rownames("NAME")

NormalId <- colnames(expMatrix_miRNA)[substr(colnames(expMatrix_miRNA),14,15) =="11"]
TumorId <- colnames(expMatrix_miRNA)[substr(colnames(expMatrix_miRNA),14,15)!="11" & substr(colnames(expMatrix_miRNA),14,15)!="02"  ]

duplicated_sample <- colnames(expMatrix_miRNA)[ substr(colnames(expMatrix_miRNA),14,15)=="02"  ]
TumorId[substr(TumorId,1,12) %in% substr(duplicated_sample,1,12) ]
exprSet1 <- expMatrix_miRNA[,c(NormalId,TumorId)]
exprSet1_miRNA <- exprSet1 %>% rownames_to_column("id")
write.table(exprSet1_miRNA,file = "01rawData/TCGA-LUAD_miRNA_matrue_RPM.xls",row.names = F,quote = F,sep = "\t",)

# 519 tumor cases, 46 normal cases
group_list=c(rep("Normal",length(NormalId)),
             rep("Tumor",length(TumorId)))  
group_list=factor(group_list)

group_list <- relevel(group_list, ref="Normal")
table(group_list)
# group_list
# Normal  Tumor 
# 46    519 
# Remove genes with low expression
exprSet1=exprSet1[rowMeans(exprSet1)>0.001,]
exprSet2 <- log2(exprSet1+1)

# Differential miRNA analysis
design=model.matrix(~ group_list)
fit=lmFit(exprSet2,design)
fit=eBayes(fit)
res=topTable(fit,adjust='fdr',coef="group_listTumor",number=Inf)
allDiff <- na.omit(res)

logFCCutoff <- 1
pvalueCutoff <- 0.05

outDiff=allDiff[(abs(allDiff$logFC)>logFCCutoff & allDiff$adj.P.Val<pvalueCutoff),]
outDiff <- outDiff %>% rownames_to_column(var = "id")
write.table(outDiff,file="02DEGs/TCGA-LUAD_miRNA_limma_diff.xls",row.names=F,quote=F,sep = "\t")

allDiff <- allDiff  %>% rownames_to_column(var = "id")
write.table(allDiff, '02DEGs/TCGA-LUAD_miRNA_limma_alldiff.xls', sep = '\t', row.names=F, quote = FALSE)


# Plot heatmap of differentially expressed miRNAs
library(pheatmap)
geneNum=50
outDiff=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
diffGeneName=outDiff$id 
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(geneNum*2)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp= exprSet1[hmGene,]

hmExp = log2(hmExp+1)

max(hmExp)
min(hmExp)

# Normal  Tumor 
# 46    519
Type= factor(c(rep("Normal",46),rep("Tumor",519)),levels = c("Normal","Tumor"))
names(Type)=colnames(hmExp)
Type=as.data.frame(Type)

pdf(file="02DEGs/miRNA_heatmap.pdf",height=8,width=9)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8,
         border=FALSE)
dev.off()


# Define significance
library(ggplot2)
library(ggrepel)
res <- allDiff

Significant=ifelse((res$adj.P.Val< 0.05 & abs(res$logFC)> 1), ifelse(res$logFC > 1,"Up","Down"), "Not")
# Plot volcano plot

p = ggplot(res, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#004BFB", "#BFBFBF", "#F91F10"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-1,1), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
  xlab(expression("log"["2"]*"FC"))+
  ylab(expression("-log"["10"]*"adj.P.Val"))+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()

pdf("02DEGs/miRNA_vol.pdf",width=5.5,height=5)
print(p)
dev.off()

