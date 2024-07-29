# Load R packages
library(tidyverse)

##########################################
# starbase predicts miRNA-lncRNA interactions
##########################################
# Read raw regulatory relationship data between lncRNA and miRNA
mir_lnc=read.table("08ceRNA/starbase/lncRNA_miRNA_interaction.txt", sep="\t", stringsAsFactors=F, header=T)
# Extract only 'miRNAname' and 'geneName' columns and remove duplicate rows
mir_lnc=unique(mir_lnc[,c("miRNAname","geneName")])

head(mir_lnc)

# Convert data frame to list with miRNA names as list names

mir_lnc_list=unstack(mir_lnc,geneName~miRNAname)

# miRNA data is saved in the differential miRNA results file; let's read it into R
miRNA <- read.table("08ceRNA/candidates/TCGA-LUAD_miRNA_limma_diff.xls",sep = "\t",header = T,check.names = F)
mir_candidate <- miRNA$id
# Obtain subset of list to get regulatory relationships between miRNA and lncRNA, then convert to data frame
mir_target_lnc=stack(mir_lnc_list[mir_candidate])
# Rename columns of data frame
names(mir_target_lnc)=c("lncRNA","miRNA")
mir_target_lnc <- mir_target_lnc[c("miRNA","lncRNA")]
# Write miRNA-lncRNA regulatory relationships to a file

write.table(file="08ceRNA/results/starbase_miRNA_to_lnc.txt",mir_target_lnc,quote=F,sep="\t",row.names = F)


#######################
# Prediction results of miRNA target genes
########################
diff_mRNA=read.table("08ceRNA/candidates/TCGA-LUAD_mRNA_limma_diff.xls",header=T,sep="\t",stringsAsFactors = F,check.names = F)
diff_miRNA=read.table("08ceRNA/candidates/TCGA-LUAD_miRNA_limma_diff.xls",header=T,sep="\t",stringsAsFactors = F,check.names = F)
diff_lncRNA=read.table("08ceRNA/candidates/TCGA-LUAD_lncRNA_limma_diff.xls",header=T,sep="\t",stringsAsFactors = F,check.names = F)

miRNA_Target <- read.table("08ceRNA/miRNA_Target/result.xls",sep = "\t",header = T,check.names = F)
miRNA_Target <- miRNA_Target[miRNA_Target$Gene %in% c("CCNA2","SLC2A1"),]

starbase_lncRNA <- read.table("08ceRNA/results/starbase_miRNA_to_lnc.txt",sep = "\t",header = T,check.names = F)
starbase_lncRNA <- starbase_lncRNA[starbase_lncRNA$lncRNA %in% diff_lncRNA$id,]

inter_mir <- intersect(miRNA_Target$miRNA,starbase_lncRNA$miRNA)
miRNA_Target <- miRNA_Target[miRNA_Target$miRNA %in% inter_mir,c("miRNA", "Gene")]
starbase_lncRNA <- starbase_lncRNA[starbase_lncRNA$miRNA %in% inter_mir,]
# Integrate targeting relationships

ceRNA <- merge(starbase_lncRNA, miRNA_Target, by = 'miRNA')
ceRNA$link <- 1
ceRNA <- reshape::melt(ceRNA, id = 'link')

variable <- summary(ceRNA$variable)
ceRNA$flow <- rep(1:variable[1], length(variable))

head(ceRNA) 

mycol <- c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462',
           '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8',
           '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#66C2A5')

#ggalluvial 
library(ggalluvial)

p <- ggplot(ceRNA, aes(x = variable, y = link,
                       stratum = value, alluvium = flow, fill = value)) +
  geom_stratum(width = 3/7) + 
  geom_flow(aes.flow = 'forward',width = 3/7) +  
  scale_fill_manual(values = mycol) +  
  geom_text(stat = 'stratum', infer.label = TRUE, size = 2.5) + 
  scale_x_discrete(limits = c('lncRNA', 'miRNA', 'Gene')) +  
  labs(x = '', y = '') + 
  theme(legend.position = 'none', panel.background = element_blank(),
        line = element_blank(), axis.text.y = element_blank())

ggsave("08ceRNA/ceRNA.pdf",width = 5,height = 5)


## ceRNA network miRNA and mRNA survival analysis
library(survival)
library(survminer)
clincal <- read.table("01rawData/TCGA-LUAD_filtered_clinical.xls",sep = "\t",header = T,check.names = F)

expFile="01rawData/TCGA-LUAD_miRNA_matrue_RPM.xls"
expr=read.table(expFile,sep="\t",header=T,check.names=F,row.names = 1)

ceRNA_miRNA= c("hsa-miR-129-5p","hsa-miR-199b-5p","hsa-miR-218-5p")

TumorId <- colnames(expr)[substr(colnames(expr),14,15)!="11" & substr(colnames(expr),14,15)!="02"  ]

exprSet_TEX <- expr[rownames(expr) %in% ceRNA_miRNA,TumorId]

exprSet_TEX <- as.data.frame( t(exprSet_TEX))
# exprSet$Id <- rownames(exprSet)
exprSet_TEX$id <- substr(rownames(exprSet_TEX),1,12)
exprSet_TEX_Time <- merge(clincal,exprSet_TEX,by="id")
# exprSet_Time <- exprSet_Time[,c(ncol(exprSet_Time),2:(ncol(exprSet_Time)-1))]
exprSet_TEX_Time[1:3,1:4]

write.table(exprSet_TEX_Time,file = "08ceRNA/TCGA-LUAD_ceRNA_miRNA_exprSetTime.xls",row.names = F,quote = F,sep = "\t")

library(survival)
library(survminer)

inputfile <- "08ceRNA/TCGA-LUAD_ceRNA_miRNA_exprSetTime.xls"
outputfile <- "TCGA-LUAD"
gene <- c("hsa-miR-129-5p",  "hsa-miR-199b-5p", "hsa-miR-218-5p")
rt <- read.table(inputfile,header=T,sep="\t",check.names = F)
i = "hsa-miR-218-5p"
res.cut <- surv_cutpoint(rt, time = "Overall_Survival", event = "Vital_Status",minprop = 0.3,
                         variables = i)
group=as.vector(ifelse(rt[,i]> res.cut$cutpoint$cutpoint,"High","Low"))
diff=survdiff(Surv(Overall_Survival, Vital_Status) ~group,data = rt)
pValue=1-pchisq(diff$chisq,length(diff$n) - 1)

median(rt$Overall_Survival)

pval = ifelse(pValue < 0.001, "p < 0.001", 
              paste("p = ",round(pValue,3), sep = ""))
pval

fit <- survfit(Surv(Overall_Survival, Vital_Status) ~ group, data = rt)
ggsurvplot(fit,
           pval = pval,
           linetype = "solid",  
           palette = c("#FF0033","#1B9E77"),
           title = i,
           ylab = "Overall survival probability",
           xlab = " Time (Days)",
           legend = c(0.7,0.8),
           legend.labs =c(" High expression "," Low expression "),
           legend.title="",
           risk.table = F,
           risk.table.title="",
           tables.height = 0.2,
           ggtheme = theme_bw(),
           tables.theme = theme_cleantable()
)
dev.copy2pdf(file = paste0("08ceRNA/",i,"_survival.pdf") , width = 4,height = 4)
dev.off()

# model gene KM survival curve
inputfile <- "03Cox/TCGA-LUAD_riskScore.xls"
outputfile <- "TCGA-LUAD"
gene <- c("CCNA2",  "SLC2A1")
rt <- read.table(inputfile,header=T,sep="\t",check.names = F)
i = "SLC2A1"
res.cut <- surv_cutpoint(rt, time = "Overall_Survival", event = "Vital_Status",minprop = 0.3,
                         variables = i)
group=as.vector(ifelse(rt[,i]> res.cut$cutpoint$cutpoint,"High","Low"))
diff=survdiff(Surv(Overall_Survival, Vital_Status) ~group,data = rt)
pValue=1-pchisq(diff$chisq,length(diff$n) - 1)

median(rt$Overall_Survival)

pval = ifelse(pValue < 0.001, "p < 0.001", 
              paste("p = ",round(pValue,3), sep = ""))
pval

fit <- survfit(Surv(Overall_Survival, Vital_Status) ~ group, data = rt)
ggsurvplot(fit,
           pval = pval,
           linetype = "solid",  
           palette = c("#FF0033","#1B9E77"), 
           title = i,
           ylab = "Overall survival probability",
           xlab = " Time (Days)",
           legend = c(0.7,0.8),
           legend.labs =c(" High expression "," Low expression "),
           legend.title="",
           risk.table = F,
           risk.table.title="",
           tables.height = 0.2,
           ggtheme = theme_bw(),
           tables.theme = theme_cleantable()
)
dev.copy2pdf(file = paste0("08ceRNA/",i,"_survival.pdf") , width = 4,height = 4)
dev.off()
