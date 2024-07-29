rm(list = ls())  
options(stringsAsFactors = F)
# Load packages
library(ggplot2)
library(grid)
library(gridExtra)

expFile="07TIDE/TCGA-LUAD_risk_exp.xls"                                                 
pdFile = "03Cox/TCGA-LUAD_riskScore.xls"
# Custom function for calculating Immunophenoscore (IPS)
ipsmap<- function (x) {
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=0)
    }
  }
  return(ips)
}
# Assign colors
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
mapcolors<-function (x) {
  za<-NULL
  if (x>=3) {
    za=1000
  } else {
    if (x<=-3) {
      za=1
    } else {
      za=round(166.5*x+500.5,digits=0)
    }
  }
  return(my_palette[za])
}
my_palette2 <- colorRampPalette(c("black", "white"))(n = 1000)
mapbw<-function (x) {
  za2<-NULL
  if (x>=2) {
    za2=1000
  } else {
    if (x<=-2) {
      za2=1
    } else {
      za2=round(249.75*x+500.5,digits=0)
    }
  }
  return(my_palette2[za2])
}
# Read expression profile file
Expr <- read.table(expFile,sep="\t",header=T,row.names=1,check.names=F)
Expr[1:3,1:3]

tumsam <- colnames(Expr)[substr(colnames(Expr),14,15) == "01"] 
gene_expression <- log2(Expr[,tumsam] + 1) 
sample_names <- tumsam

## 修改几个基因名
## IPS基因权重文件中的基因名跟常用的gene symbol不同，我们把表达矩阵里的gene symbol修改成跟IPS基因列表一致的基因名。
if(is.element("CAVIN2",rownames(gene_expression))) {
  rownames(gene_expression) <- gsub("CAVIN2","SDPR",rownames(gene_expression))
}

if(is.element("CCL3L3",rownames(gene_expression))) {
  rownames(gene_expression) <- gsub("CCL3L3","CCL3L1",rownames(gene_expression))
}
if(is.element("IARS1",rownames(gene_expression))) {
  rownames(gene_expression) <- gsub("IARS1","IARS",rownames(gene_expression))
}
if(is.element("DARS1",rownames(gene_expression))) {
  rownames(gene_expression) <- gsub("DARS1","DARS",rownames(gene_expression))
}
# Read IPS-related genes and their weights
IPSG <- read.table("07TIDE/IPS/easy_input_IPS_genes.txt",header=TRUE, sep="\t",check.names=FALSE,stringsAsFactors = FALSE,row.names = NULL)
head(IPSG)
unique_ips_genes <- as.vector(unique(IPSG$NAME))

# Initialize data
IPS <- MHC <- CP <- EC <- SC <- AZ <- NULL

# Obtain gene names from expression profile
GVEC <- row.names(gene_expression)

# Obtain gene names from IPS gene file
VEC <- as.vector(IPSG$GENE)

# Match genes and identify missing genes
ind <- which(is.na(match(VEC,GVEC)))
MISSING_GENES <- VEC[ind]

dat <- IPSG[ind,]
if (length(MISSING_GENES) > 0) { 
  message(paste0("--differently named or missing genes: ",paste(MISSING_GENES,collapse = ", ")))
  print(IPSG[ind,])
  message("please check data and make sure all genes matches!")
} else {
  message("--all genes matched!") 
}
## Iterate over each sample
outTab <- NULL
for (i in 1:length(sample_names)) { 
  GE <- gene_expression[[i]]
  mGE <- mean(GE)
  sGE <- sd(GE)
  Z1 <- (gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
  W1 <- IPSG$WEIGHT
  WEIGHT <- MIG <- NULL
  k <- 1
  for (gen in unique_ips_genes) {
    MIG[k] <- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
    WEIGHT[k] <- mean(W1[which (as.vector(IPSG$NAME)==gen)])
    k <- k + 1
  }
  WG <- MIG * WEIGHT
  MHC[i] <- mean(WG[1:10])
  CP[i] <- mean(WG[11:20])
  EC[i] <- mean(WG[21:24])
  SC[i] <- mean(WG[25:26])
  AZ[i] <- sum(MHC[i],CP[i],EC[i],SC[i])
  IPS[i] <- ipsmap(AZ[i])
  
  tmp <- as.data.frame(t(data.frame(WG))); colnames(tmp) <- unique_ips_genes; rownames(tmp) <- sample_names[i]
  outTab <- rbind.data.frame(outTab,tmp)

}

# Construct results, including scores and IPS

DF <- data.frame(SAMPLE=sample_names,
                 MHC=MHC,
                 EC=EC,
                 SC=SC,
                 CP=CP,
                 AZ=AZ,
                 IPS=IPS,
                 stringsAsFactors = F)

# Output z-score values to file
write.table(outTab, file = "07TIDE/IPS/output_zscore.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep="\t") 

# Output IPS to file
write.table(DF,file = "07TIDE/IPS/output_IPS.txt", row.names = FALSE, col.names = TRUE, quote=FALSE, sep="\t") 


DF<- read.table("07TIDE/IPS/output_IPS.txt",sep="\t",header=T,check.names=F,row.names = 1)

DF$id <- substr(rownames(DF),1,12 )

pd<- read.table(pdFile,sep="\t",header=T,check.names=F)

boxdata <- merge(DF,pd,by="id")
boxdata$risk <- factor(boxdata$risk,levels=c("Low Risk","High Risk"))

library(ggpubr)
mycol <- c("#028846","red")
ggplot(boxdata,aes(risk,IPS,fill=risk))+
  geom_boxplot(outlier.colour = NA,notch = F,size = 0.3,width=.3,alpha = 0.4)+
  geom_jitter(shape = 21,size=2,width = 0.2,alpha = 0.4)+
  geom_violin(position = position_dodge(width = .75), 
              size = 0.3,alpha = 0.4,trim = T)+
  ylim(2, 10)+
  theme_classic()+ ylab('immunophenoscore')+ggtitle('')+
  theme(legend.position = 'none',
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=10),
        axis.title.x = element_blank())+
  scale_fill_manual(values=mycol)+
  stat_compare_means(comparisons = list(c("Low Risk","High Risk")),
                     label = 'p.signif')
ggsave("07TIDE/IPS/ips_boxplot.pdf",width = 3,height = 3.5)




