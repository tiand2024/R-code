rm(list = ls()) 
options(stringsAsFactors = F)
# Drug prediction R package pRRophetic
library(pRRophetic)
library(ggplot2)
library(ggpubr)
library(cowplot)

expFile="07TIDE/TCGA-LUAD_risk_exp.xls"                                                 
pdFile = "03Cox/TCGA-LUAD_riskScore.xls"

# Read immune results file and preprocess the data
testExpr <- read.table(expFile,sep="\t",header=T,row.names=1,check.names=F)

colnames(testExpr) <- substr(colnames(testExpr),1,12 )
pd<- read.table(pdFile,sep="\t",header=T,row.names = 1,check.names=F)

mycol <-   c("#028846","red")
# Drug name
GCP.drug <- read.table("07TIDE/pRRophetic/drug_sig_pvalue.txt") 
GCP.drug <- GCP.drug$V1

exprData <- as.matrix(testExpr[,rownames(pd)])
jco <- c("#028846","red")
### Drug sensitivity prediction ###
GCPinfo <- GCP.IC50 <- GCP.expr <- cvOut <- predictedPtype <- predictedBoxdat <- list()# Initialize list
plotp <- list()

for (drug in GCP.drug) {
  set.seed(4) 
  cat(drug," starts!\n") 
  
 
  predictedPtype[[drug]] <- pRRopheticPredict(testMatrix = exprData ,
                                              drug = drug,
                                              tissueType = "allSolidTumors",
                                              batchCorrect = "eb",
                                              selection = 1, 
                                              dataset = "cgp2014")
  if(!all(names(predictedPtype[[drug]])==rownames(pd))) {stop("Name mismatched!\n")} 
  
  predictedBoxdat[[drug]] <- data.frame("est.ic50"=predictedPtype[[drug]],
                                        "risk"=ifelse(pd$risk == "Low Risk","Low Risk","High Risk"), 
                                        row.names = names(predictedPtype[[drug]])) 
  predictedBoxdat[[drug]]$risk <- factor(predictedBoxdat[[drug]]$risk,levels = c("Low Risk","High Risk"),ordered = F) 
  predictedBoxdat[["Paclitaxel"]]

  p <- ggplot(data = predictedBoxdat[[drug]], aes(x=risk, y=est.ic50))
  p <- p + geom_boxplot(aes(fill = risk),outlier.colour = NA,notch = F,size = 0.2,width = 0.2,alpha = 0.4,)+
    geom_jitter(aes(fill = risk),shape = 21,size=2,width = 0.2,alpha = 0.4,)+
    geom_violin(aes(fill = risk),position = position_dodge(width = .75), 
                size = 0.3,width = 0.6,alpha = 0.4,trim = T)+
    scale_fill_manual(values=jco) + 
    theme_classic()+
    stat_compare_means(comparisons = list(c("Low Risk","High Risk")),method="wilcox.test", label = "p.signif")+
    theme(legend.position="none") + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10),plot.title = element_text(size = 12, hjust = 0.5)) + 
    xlab("") + ylab("Estimated IC50") + 
    ggtitle(drug) 
  
  plotp[[drug]] <- p 
  cat(drug," has been finished!\n") 
}


p2 <- plot_grid(plotlist=plotp, ncol=5)
ggsave("07TIDE/pRRophetic/boxplot of predicted IC50_multiple.pdf", width = 10, height = 8)

# Test inter-group differences
p <- vector()
for (drug in GCP.drug) {
  tmp <- wilcox.test(as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$risk %in% "Low Risk"),"est.ic50"]),
                     as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$risk %in% "High Risk"),"est.ic50"]),alternative = "less")$p.value
  p <- append(p,tmp) 
}
names(p) <- GCP.drug
print(p) #

write.table(p,"07TIDE/pRRophetic/drug_sig_pvalue.txt", quote = F, sep = "\t")



# TIDE compares immune response across different subtypes

library(pheatmap)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
expFile="07TIDE/TCGA-LUAD_risk_exp.xls"  
TIDE <- read.table(expFile,sep="\t",header=T,row.names=1,check.names=F)

# two direction median centered
TIDE <- sweep(TIDE,2, apply(TIDE,2,median,na.rm=T))
TIDE <- sweep(TIDE,1, apply(TIDE,1,median,na.rm=T))
write.table(TIDE,"07TIDE/SubMap/TIDE_input.self_subtract",sep = "\t",row.names = T,col.names = NA,quote = F)

# submap predicts subtype immune therapy responsiveness
# Define custom function to generate data format required for submap

generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

# Create data format for submap
skcm.immunotherapy.logNC <- read.table("07TIDE/SubMap/skcm.immunotherapy.47samples.log2CountsNorm.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F) #原文提供的log2转化的标准化count值
rownames(skcm.immunotherapy.logNC) <- toupper(rownames(skcm.immunotherapy.logNC)) 
skcm.immunotherapy.info <- read.table("07TIDE/SubMap/skcm.immunotherapy.47sampleInfo.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$label),]
skcm.immunotherapy.info$rank <- rep(c(1,2,3,4),times=as.character(table(skcm.immunotherapy.info$label))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R

expFile="07TIDE/TCGA-LUAD_risk_exp.xls"  
tmp <- read.table(expFile,sep="\t",header=T,row.names=1,check.names=F)
TumorId <- colnames(tmp)[substr(colnames(tmp),14,15)!="11" & substr(colnames(tmp),14,15)!="02"  ]
tmp <- tmp[,TumorId]
riskScore <- read.table("03Cox/TCGA-LUAD_riskScore.xls",header = T,sep = "\t")
low <- riskScore$id[ riskScore$risk== "Low Risk"]
high <- riskScore$id[ riskScore$risk== "High Risk"]


# colnames(tmp) <- substr(colnames(tmp),1,12 )
GENELIST <- intersect(rownames(tmp),rownames(skcm.immunotherapy.logNC))

sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST,rownames(skcm.immunotherapy.info)]

# Generate output file name for the output data

gct_file <- "skcm.immunotherapy.for.SubMap.gct"
cls_file <- "skcm.immunotherapy.for.SubMap.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")
# Extract samples of a subtype, ordered sequentially
samples.C1 <- TumorId[substr(TumorId,1,12) %in% low]
samples.C2 <- highId <- TumorId[substr(TumorId,1,12) %in% high]

sam_info <- data.frame("Type"=c(samples.C1,samples.C2),row.names = c(samples.C1,samples.C2))
sam_info$rank <- rep(c(1,2),times=c(length(samples.C1),length(samples.C2))) 

# Generate output file name
gct_file <- "TCGA-LUAD_risk.for.SubMap.gct"
cls_file <- "TCGA-LUAD_risk.SubMap.cls"

in_gct <- log2(tmp[GENELIST,rownames(sam_info)] + 1) 

generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

library(GenePattern)
gctfileA <- '07TIDE/SubMap/TCGA-LUAD_risk.for.SubMap.gct'
clsfileA <- '07TIDE/SubMap/TCGA-LUAD_risk.SubMap.cls'

gctfileB <- '07TIDE/SubMap/skcm.immunotherapy.for.SubMap.gct'
clsfileB <- '07TIDE/SubMap/skcm.immunotherapy.for.SubMap.cls'

source('07TIDE/SubMap/submap.R')
submap.main(  gctfileA,
              gctfileB,
              clsfileA,
              clsfileB,
              output.filename="07TIDE/SubMap/SubMap",
              ntag=100,
              nperm=100,
              nperm.fisher=1000,
              weighted.score.type=1,
              null.dist="pool",
              p.corr="Bonferroni",
              clust.row=1,
              clust.col=1,
              nom.p.mat="T",
              create.legend="T",
              rnd.seed=47365321
              # rnd.seed=321
              ) 



# Fill values from the file LUAD_ceRNA/07TIDE/SubMap/SubMap_SubMapResult.txt into respective positions
# Plot heatmap of nominal and adjusted p-values from input file

library(pheatmap)
heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
cherry    <- "#700353"
lightgrey <- "#dcddde"
submap.data.A <- matrix(c( 0.5114885, 0.54145854, 1.00000000, 0.4365634,
                           0.4815185, 0.03296703, 0.09190809, 0.6673327),nrow = 2,byrow = T)
submap.data.B <- matrix(c(1, 1.0000000, 1.0000000,  1,
                          1, 0.2637363, 0.7352647,  1),nrow = 2,byrow = T)
submap.data <- rbind(submap.data.A,submap.data.B)
row.names(submap.data) <- c('Low Risk_p','High Risk_p','Low Risk_b','High Risk_b')
colnames(submap.data) <- c("CTAL4-noR","CTLA4-R","PD1-noR","PD1-R")
annotation_row <- data.frame(pvalue=rep(c('Nominal p value','Bonferroni corrected'),c(2,2)))
rownames(annotation_row) <- row.names(submap.data)
pdf(file = '07TIDE/SubMap/heatmap_submap.pdf',width= 12,height= 8)
pheatmap(submap.data,
         show_colnames = T,
         color = heatmap.YlGnPe[5:1],
         display_numbers = matrix(ifelse(submap.data < 0.05,'p < .05',''),nrow(submap.data)),number_format = "%.3f",
         cluster_rows = F,cluster_cols = F,
         annotation_row=annotation_row,
         annotation_colors = list(pvalue=c("Nominal p value"=lightgrey,"Bonferroni corrected"=cherry)),
         cellwidth=30,cellheight=30,main = "",fontsize_number = 9,
         #fontsize_row = 20,fontsize_col = 25,
         gaps_row=2,
         fontsize=12)
dev.off()


# Plot boxplot of TIDE response

ann <-read.table("03Cox/TCGA-LUAD_riskScore.xls",header = T,sep = "\t",check.names = F)
head(ann)
print(table(ann$risk))
# Obtain output file TIDE_output.csv following the TIDE usage tutorial in the reference folder

TIDE.res <- read.csv("07TIDE/SubMap/TIDE_output.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
TIDE.res$id <- substr(rownames(TIDE.res),1,12 )


# Test correlation between immune therapy responsiveness and subtype; p < 0.05 indicates correlation

ann$TIDE <- TIDE.res[ann$id,"Responder"]
print(table(ann$TIDE,ann$risk))
print(fisher.test(table(ann$TIDE,ann$risk))) 



boxdata <- merge(TIDE.res,ann,by="id")
boxdata$risk <- factor(boxdata$risk,levels=c("Low Risk","High Risk"))

library(ggpubr)
mycol <- c("#028846","red")
ggplot(boxdata,aes(risk,TIDE.x,fill=risk))+
  geom_boxplot(outlier.colour = NA,notch = F,size = 0.3,width=.3,alpha = 0.4)+
  geom_jitter(shape = 21,size=2,width = 0.2,alpha = 0.4)+
  geom_violin(position = position_dodge(width = .75), 
              size = 0.3,alpha = 0.4,trim = T)+
  #ylim(0, 3)+
  theme_classic()+ ylab('TIDE prediction score')+ggtitle('')+
  theme(legend.position = 'none',
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=10),
        axis.title.x = element_blank())+
  scale_fill_manual(values=mycol)+
  stat_compare_means(comparisons = list(c("Low Risk","High Risk")),
                     label = 'p.signif')#+
#stat_compare_means(label.y = max(boxdata$IPS)+0.5)
ggsave("07TIDE/SubMap/TIDE_boxplot.pdf",width = 3,height = 3.5)


# Plot scatter plot of correlation

p1=  ggplot(boxdata,aes(TIDE.x,riskScore))+
  geom_point(shape=21 ,size=2,alpha=0.7,fill="#6baed6")+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+ 
  theme_cowplot()+  
  theme(axis.title = element_text(size = 12))+
  xlab("TIDE prediction score")+
  ylab("riskScore")

p1 <- p1+stat_cor(method = "pearson",digits = 3)
ggsave("07TIDE/cor_plot.pdf",width = 3,height = 3)
