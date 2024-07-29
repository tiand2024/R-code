#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

                                            
expFile="01rawData/TCGA-LUAD_TPM_mRNA.xls"                                                 

source("05immune/CIBERSORT.R")
# nperm specifies the number of permutations.
# QN is set to T for microarray and F for sequencing. For sequencing data, TPM is preferred.
results=CIBERSORT("05immune/LM22.txt", expFile, perm=1000, QN=FALSE)


# Plotting visually appealing CIBERSORT immune cell fraction bar chart
library(tidyverse)
library(ggplot2)
library(paletteer)
mycol<- colorRampPalette(paletteer_d("ggthemes::Classic_20",n=20))(22)
df <- sample(22)
barplot(df,col = mycol)
# Read immune results file and process the data
pFilter=0.05 
pd<- read.table("03Cox/TCGA-LUAD_riskScore.xls",sep="\t",header=T,check.names=F)
immune=read.table("05immune/CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)
TumorId <- rownames(immune)[substr(rownames(immune),14,15)!="11" & substr(rownames(immune),14,15)!="02"  ]
immune=immune[TumorId,]
# immune=immune[immune[,"P-value"]<pFilter,]
immune=(immune[,1:(ncol(immune)-3)])
sort_cell <- names(sort(apply(immune,2,sum),decreasing = T))


data <- immune %>% as.data.frame() %>%
  rownames_to_column("id") %>% 
  gather(key = immuneCell,value = Proportion,-id)
data$id <- substr(data$id ,1,12)
data <- merge(data,pd,by="id")
data$immuneCell <- factor(data$immuneCell ,levels = sort_cell)
data$risk <- factor(data$risk ,levels = c("Low Risk","High Risk"))
p1 <- ggplot(data,aes(id,Proportion,fill = immuneCell)) + 
      geom_bar(stat = "identity") +
      labs(fill = "Cell Type",x = "",y = "Relative Percent") + 
      theme_bw() +
      facet_grid(.~risk,scales = "free",space="free_x")+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y=element_text(colour = "black",size = 8),
            legend.position = "bottom") + 
      scale_y_continuous(expand = c(0.01,0)) +
      scale_fill_manual(values = mycol)

p2 <- p1 + theme(legend.text = element_text(colour = "black",size = 8)) + 
  theme(legend.title = element_text(size = 10,colour = "black")) 

pdf("05immune/CIBERSORT_Proportion.pdf",height = 6,width = 12)
p2
dev.off()

data$Proportion  <- ifelse(data$Proportion  ==0,0.001,data$Proportion )
library(ggpubr)
p3 <- ggboxplot(data, x = "immuneCell", y = "Proportion",
                color = "risk", palette = c("#028846","red"),
                add = "jitter", 
                legend.title = "",
                add.params = list(alpha=0.6),)+ 
                rotate_x_text(45) +xlab("")+ylab("Fraction")
                      
p3 + stat_compare_means(aes(group = risk),method="wilcox.test", label = "p.signif")
ggsave("05immune/CIBERSORT_boxplot.pdf",height=6,width=10)

# Plotting heatmap of correlation between differential genes and immune cells
# Read in expression data


tcga_expr <- read.table("01rawData/TCGA-LUAD_TPM_mRNA.xls", sep = "\t",row.names = 1,header = T,check.names = F )
tcga_expr[1:4,1:5]
tcga_expr <- tcga_expr[,rownames(immune)]
# Batch compute correlations
genelist <- c("CCNA2","SLC2A1") 

immuscore <- function(gene){
  y <- as.numeric(tcga_expr[gene,])
  colnames <- colnames(immune)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(immune[,x]), y , method="spearman")
    data.frame(gene=gene,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
immuscore("CCNA2")

data <- do.call(rbind,lapply(genelist,immuscore))
head(data)
data <- na.omit(data)
write.table(data, "05immune/gene_immuneCell_correlation.xls",sep = "\t" ,quote = F, row.names = F)
data$pstar <- ifelse(data$p.value < 0.05,
                     ifelse(data$p.value < 0.01,"**","*"),
                     "")
data$pstar[1:20]

ggplot(data, aes(immune_cells, gene)) + 
  geom_tile(aes(fill = cor), colour = "white",size=1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=pstar),col ="black",size = 5)+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))+
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
ggsave("05immune/correlation.pdf",width = 7,height = 2.8)


# Output CIBERSORT results for high and low-risk groups
pd<- read.table("03Cox/TCGA-LUAD_riskScore.xls",sep="\t",header=T,check.names=F)
immune=read.table("05immune/CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)
TumorId <- rownames(immune)[substr(rownames(immune),14,15)!="11" & substr(rownames(immune),14,15)!="02"  ]
immune=immune[TumorId,]
# immune=immune[immune[,"P-value"]<pFilter,]
immune=(immune[,1:(ncol(immune)-3)])
immune$id <- substr(rownames(immune) ,1,12)
out <- merge(immune,pd[c("id","risk")],by="id")
write.table(out, "05immune/immuneCell_risk_summary.xls",sep = "\t" ,quote = F, row.names = F)
