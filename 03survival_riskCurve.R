rm(list=ls())
# Kaplan-Meier Survival Analysis
library(survival)
library(survminer)
inputfile<- "03Cox/TCGA-LUAD_riskScore.xls"
outputfile <- "TCGA-LUAD"
rt <- read.table(inputfile,header=T,sep="\t",check.names = F)
diff=survdiff(Surv(Overall_Survival, Vital_Status) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,length(diff$n) - 1)

median(rt$Overall_Survival)

pval = ifelse(pValue < 0.001, "p < 0.001", 
              paste("p = ",round(pValue,3), sep = ""))
pval

fit <- survfit(Surv(Overall_Survival, Vital_Status) ~ risk, data = rt)
ggsurvplot(fit,
           pval = pval,
           linetype = "solid",  
           palette = c("#FF0033","#1B9E77"), 
           title = outputfile,
           ylab = "Overall survival probability",
           xlab = " Time (Days)",
           legend = c(0.8,0.8),
           legend.labs = c("High Risk","Low Risk"),
           legend.title="",
           risk.table = T,
           risk.table.title="",
           tables.height = 0.2,
           ggtheme = theme_bw(),
           tables.theme = theme_cleantable()
)
dev.copy2pdf(file = paste0("04ROC/",outputfile,"_survival.pdf") , width = 4.3,height = 5)
dev.off()


library(survivalROC)
period_time <- 365

rt <- read.table(inputfile,header=T,sep="\t",check.names = F)
rocCol=c("#3cb346","#eeb401", "#ef1828","#942d8d")
aucText <- c()
pdf(file=paste0("04ROC/",outputfile,"_ROC.pdf"),width=4.5,height=4.5)
par(mar=c(4,4,2,1),mgp=c(2,0.5,0))
roc=survivalROC(Stime=rt$Overall_Survival,
                status=rt$Vital_Status,
                marker = rt$riskScore ,
                predict.time =period_time ,
                method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="#1B9E77",
     xlab="1-Specificity (False Positive)", ylab="Sensitivity (True Positive)",
     main=outputfile,
     lwd = 1.5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1,lty= 3)
j <- 0
for(i in c(1,2,3)){
        
        roc=survivalROC(Stime=rt$Overall_Survival, status=rt$Vital_Status, marker = rt$riskScore,
                        predict.time =period_time*i, method="KM")
        
        j=j+1
        aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",roc$AUC),")"))
        lines(roc$FP, roc$TP, type="l",col=rocCol[j],xlim=c(0,1), ylim=c(0,1),lwd = 1.8, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
        
}

legend("bottomright", aucText,lty= 1,lwd=1.8,bty="n",col=rocCol)
dev.off()

# Training set risk curve
rt <- read.table(inputfile,header=T,sep="\t",check.names = F)

rt$Vital_Status <- ifelse(rt$Vital_Status==0,'Alive','Dead')
rt$risk <- factor(rt$risk,levels = c("Low Risk","High Risk"))
rt=rt[order(rt$riskScore),] 
rt$patient = 1:nrow(rt)

p1 <- ggplot(data=rt,mapping = aes(x=patient,y=riskScore))+
        geom_point(aes(color=risk),size=0.5)+
        theme_classic()+
        scale_color_manual(name="Risk Group",values = c("#028846","red"))+
        #x-axis
        theme(axis.line.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())+
        ylab("Risk Score")+
        geom_vline(xintercept=sum(rt$risk=="Low Risk"),colour="black", linetype=2,size=0.5)+
        scale_x_continuous(expand = c(0,0))
p1 


p2 <- ggplot(data = rt,mapping = aes(x=patient,y=Overall_Survival,col=Vital_Status))+
        geom_point(aes(col=Vital_Status),size=0.5)+ scale_color_manual(name="Status",values = c("#028846","red")) +
        theme_classic()+
        #x a-xis
        theme(axis.line.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.x=element_blank())+
        ylab("Survival time(Days)")+
        geom_vline(xintercept = sum(rt$risk=="Low Risk") ,linetype=2,size=0.5)+
        scale_x_continuous(expand = c(0,0))
p2

middle = ggplot(rt, aes(
        x = patient,
        y = 1)) +
        geom_tile(aes(fill = risk))+theme_classic()+
        scale_fill_manual(name="Risk Group",  values = c("#028846","red"))+
        theme(
                legend.position = "none",
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                axis.title = element_blank(),
                plot.margin = unit(c(3,0,-3,3), "cm")
        )+theme(legend.title = element_text(size = 12),
                legend.text = element_text(size=12))+
        scale_x_continuous(expand = c(0,1))
middle

rt1=rt[c(11:(ncol(rt)-3))] 
max(rt1)
min(rt1)
pheatmap::pheatmap(t(rt1),scale = "row")
rt1=rt1 %>% scale() %>%as.data.frame()
rt1$id = 1:nrow(rt1)
rt1= tidyr::gather(rt1,variable,value,-id)
rt1$variable <- factor(rt1$variable ,levels =c("CCNA2","SLC2A1") )
p3 = ggplot(rt1, aes_string(x = 'id',y = 'variable',fill = 'value')) +
        geom_raster() + labs(fill="Expression")+
        theme(
                panel.grid = element_blank(),
                panel.background = element_blank(),
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                axis.title = element_blank(),
                plot.background = element_blank() #the key to avoide legend overlap
        ) +
        scale_fill_gradient2(
                low = "blue",
                mid = "white",
                high = "red"
        ) +
        scale_x_continuous(expand = c(0,0))
p3

relative_heights=c(1,0.06,1)
library(aplot)

p1 %>% 
        insert_bottom(p2,height = relative_heights[1])%>% 
        insert_bottom(middle,height=relative_heights[2]) %>%
        insert_bottom(p3,height=relative_heights[3])

ggsave(file = paste0("04ROC/",outputfile,"_riskCurve.pdf") , width = 5,height = 5)

# GEO validset
inputfile<- "03Cox/GSE30219_validSet_riskScore.xls"
outputfile <- "GSE30219"
rt <- read.table(inputfile,header=T,sep="\t",check.names = F)
rt$Vital_Status <- ifelse(rt$Vital_Status =="ALIVE",0,1)
diff=survdiff(Surv(Overall_Survival, Vital_Status) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,length(diff$n) - 1)

median(rt$Overall_Survival)

pval = ifelse(pValue < 0.001, "p < 0.001", 
              paste("p = ",round(pValue,3), sep = ""))
pval

fit <- survfit(Surv(Overall_Survival, Vital_Status) ~ risk, data = rt)
ggsurvplot(fit,
           pval = pval,
           linetype = "solid",  
           palette = c("#FF0033","#1B9E77"), 
           title = outputfile,
           ylab = "Overall survival probability",
           xlab = " Time (Months)",
           legend = c(0.8,0.8),
           legend.labs = c("High Risk","Low Risk"),
           legend.title="",
           risk.table = T,
           risk.table.title="",
           tables.height = 0.2,
           ggtheme = theme_bw(),
           tables.theme = theme_cleantable()
)
dev.copy2pdf(file = paste0("04ROC/",outputfile,"_survival.pdf") , width = 4.3,height = 5)
dev.off()


library(survivalROC)
period_time <- 12

rt <- read.table(inputfile,header=T,sep="\t",check.names = F)
rt$Vital_Status <- ifelse(rt$Vital_Status =="ALIVE",0,1)
rocCol=c("#3cb346","#eeb401", "#ef1828","#942d8d")
aucText <- c()
pdf(file=paste0("04ROC/",outputfile,"_ROC.pdf"),width=4.5,height=4.5)
par(mar=c(4,4,2,1),mgp=c(2,0.5,0))
roc=survivalROC(Stime=rt$Overall_Survival,
                status=rt$Vital_Status,
                marker = rt$riskScore ,
                predict.time =period_time ,
                method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="#1B9E77",
     xlab="1-Specificity (False Positive)", ylab="Sensitivity (True Positive)",
     main=outputfile,
     lwd = 1.5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1,lty= 3)
j <- 0
for(i in c(1,2,3)){
        
        roc=survivalROC(Stime=rt$Overall_Survival, status=rt$Vital_Status, marker = rt$riskScore,
                        predict.time =period_time*i, method="KM")
        
        j=j+1
        aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",roc$AUC),")"))
        lines(roc$FP, roc$TP, type="l",col=rocCol[j],xlim=c(0,1), ylim=c(0,1),lwd = 1.8, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
        
}

legend("bottomright", aucText,lty= 1,lwd=1.8,bty="n",col=rocCol)
dev.off()


rt <- read.table(inputfile,header=T,sep="\t",check.names = F)
rt$Vital_Status <- ifelse(rt$Vital_Status =="ALIVE",0,1)
rt$Vital_Status <- ifelse(rt$Vital_Status==0,'Alive','Dead')
rt$risk <- factor(rt$risk,levels = c("Low Risk","High Risk"))
rt=rt[order(rt$riskScore),] 
rt$patient = 1:nrow(rt)

p1 <- ggplot(data=rt,mapping = aes(x=patient,y=riskScore))+
        geom_point(aes(color=risk),size=0.5)+
        theme_classic()+
        scale_color_manual(name="Risk Group",values = c("#028846","red"))+
        #x-axis
        theme(axis.line.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())+
        ylab("Risk Score")+
        geom_vline(xintercept=sum(rt$risk=="Low Risk"),colour="black", linetype=2,size=0.5)+
        scale_x_continuous(expand = c(0,0))
p1 


p2 <- ggplot(data = rt,mapping = aes(x=patient,y=Overall_Survival,col=Vital_Status))+
        geom_point(aes(col=Vital_Status),size=0.5)+ scale_color_manual(name="Status",values = c("#028846","red")) +
        theme_classic()+
        #x a-xis
        theme(axis.line.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.x=element_blank())+
        ylab("Survival time(Months)")+
        geom_vline(xintercept = sum(rt$risk=="Low Risk") ,linetype=2,size=0.5)+
        scale_x_continuous(expand = c(0,0))
p2

middle = ggplot(rt, aes(
        x = patient,
        y = 1)) +
        geom_tile(aes(fill = risk))+theme_classic()+
        scale_fill_manual(name="Risk Group",  values = c("#028846","red"))+
        theme(
                legend.position = "none",
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                axis.title = element_blank(),
                plot.margin = unit(c(3,0,-3,3), "cm")
        )+theme(legend.title = element_text(size = 12),
                legend.text = element_text(size=12))+
        scale_x_continuous(expand = c(0,1))
middle

rt1=rt[c(4:(ncol(rt)-3))] 
max(rt1)
min(rt1)
pheatmap::pheatmap(t(rt1),scale = "row")
rt1=rt1 %>% scale() %>%as.data.frame()
rt1$id = 1:nrow(rt1)
rt1= tidyr::gather(rt1,variable,value,-id)
rt1$variable <- factor(rt1$variable ,levels =c("CCNA2","SLC2A1") )
p3 = ggplot(rt1, aes_string(x = 'id',y = 'variable',fill = 'value')) +
        geom_raster() + labs(fill="Expression")+
        theme(
                panel.grid = element_blank(),
                panel.background = element_blank(),
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                axis.title = element_blank(),
                plot.background = element_blank() #the key to avoide legend overlap
        ) +
        scale_fill_gradient2(
                low = "blue",
                mid = "white",
                high = "red"
        ) +
        scale_x_continuous(expand = c(0,0))
p3

relative_heights=c(1,0.06,1)
library(aplot)

p1 %>% 
        insert_bottom(p2,height = relative_heights[1])%>% 
        insert_bottom(middle,height=relative_heights[2]) %>%
        insert_bottom(p3,height=relative_heights[3])

ggsave(file = paste0("04ROC/",outputfile,"_riskCurve.pdf") , width = 5,height = 5)

# GEO validset GSE31210
inputfile<- "03Cox/GSE31210_validSet_riskScore.xls"
outputfile <- "GSE31210"
rt <- read.table(inputfile,header=T,sep="\t",check.names = F)
rt$Vital_Status <- ifelse(rt$Vital_Status =="alive",0,1)
diff=survdiff(Surv(Overall_Survival, Vital_Status) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,length(diff$n) - 1)

median(rt$Overall_Survival)

pval = ifelse(pValue < 0.001, "p < 0.001", 
              paste("p = ",round(pValue,3), sep = ""))
pval

fit <- survfit(Surv(Overall_Survival, Vital_Status) ~ risk, data = rt)
ggsurvplot(fit,
           pval = pval,
           linetype = "solid",  
           palette = c("#FF0033","#1B9E77"), 
           title = outputfile,
           ylab = "Overall survival probability",
           xlab = " Time (Days)",
           legend = c(0.8,0.2),
           legend.labs = c("High Risk","Low Risk"),
           legend.title="",
           risk.table = T,
           risk.table.title="",
           tables.height = 0.2,
           ggtheme = theme_bw(),
           tables.theme = theme_cleantable()
)
dev.copy2pdf(file = paste0("04ROC/",outputfile,"_survival.pdf") , width = 4.3,height = 5)
dev.off()


library(survivalROC)
period_time <- 365

rt <- read.table(inputfile,header=T,sep="\t",check.names = F)
rt$Vital_Status <- ifelse(rt$Vital_Status =="ALIVE",0,1)
rocCol=c("#3cb346","#eeb401", "#ef1828","#942d8d")
aucText <- c()
pdf(file=paste0("04ROC/",outputfile,"_ROC.pdf"),width=4.5,height=4.5)
par(mar=c(4,4,2,1),mgp=c(2,0.5,0))
roc=survivalROC(Stime=rt$Overall_Survival,
                status=rt$Vital_Status,
                marker = rt$riskScore ,
                predict.time =period_time ,
                method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="#1B9E77",
     xlab="1-Specificity (False Positive)", ylab="Sensitivity (True Positive)",
     main=outputfile,
     lwd = 1.5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1,lty= 3)
j <- 0
for(i in c(1,2,3)){
        
        roc=survivalROC(Stime=rt$Overall_Survival, status=rt$Vital_Status, marker = rt$riskScore,
                        predict.time =period_time*i, method="KM")
        
        j=j+1
        aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",roc$AUC),")"))
        lines(roc$FP, roc$TP, type="l",col=rocCol[j],xlim=c(0,1), ylim=c(0,1),lwd = 1.8, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
        
}

legend("bottomright", aucText,lty= 1,lwd=1.8,bty="n",col=rocCol)
dev.off()

rt <- read.table(inputfile,header=T,sep="\t",check.names = F)
rt$Vital_Status <- ifelse(rt$Vital_Status =="alive",0,1)
rt$Vital_Status <- ifelse(rt$Vital_Status==0,'Alive','Dead')
rt$risk <- factor(rt$risk,levels = c("Low Risk","High Risk"))
rt=rt[order(rt$riskScore),] 
rt$patient = 1:nrow(rt)

p1 <- ggplot(data=rt,mapping = aes(x=patient,y=riskScore))+
        geom_point(aes(color=risk),size=0.5)+
        theme_classic()+
        scale_color_manual(name="Risk Group",values = c("#028846","red"))+
        #x-axis
        theme(axis.line.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())+
        ylab("Risk Score")+
        
        geom_vline(xintercept=sum(rt$risk=="Low Risk"),colour="black", linetype=2,size=0.5)+
        scale_x_continuous(expand = c(0,0))
p1 


p2 <- ggplot(data = rt,mapping = aes(x=patient,y=Overall_Survival,col=Vital_Status))+
        geom_point(aes(col=Vital_Status),size=0.5)+ scale_color_manual(name="Status",values = c("#028846","red")) +
        theme_classic()+
        #x a-xis
        theme(axis.line.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.x=element_blank())+
        ylab("Survival time(Days)")+
        
        geom_vline(xintercept = sum(rt$risk=="Low Risk") ,linetype=2,size=0.5)+
        scale_x_continuous(expand = c(0,0))
p2

middle = ggplot(rt, aes(
        x = patient,
        y = 1)) +
        geom_tile(aes(fill = risk))+theme_classic()+
        scale_fill_manual(name="Risk Group",  values = c("#028846","red"))+
        theme(
                legend.position = "none",
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                axis.title = element_blank(),
                plot.margin = unit(c(3,0,-3,3), "cm")
        )+theme(legend.title = element_text(size = 12),
                legend.text = element_text(size=12))+
        scale_x_continuous(expand = c(0,1))
middle

rt1=rt[c(4:(ncol(rt)-3))] 
max(rt1)
min(rt1)
pheatmap::pheatmap(t(rt1),scale = "row")
rt1=rt1 %>% scale() %>%as.data.frame()
rt1$id = 1:nrow(rt1)
rt1= tidyr::gather(rt1,variable,value,-id)
rt1$variable <- factor(rt1$variable ,levels =c("CCNA2","SLC2A1") )
p3 = ggplot(rt1, aes_string(x = 'id',y = 'variable',fill = 'value')) +
        geom_raster() + labs(fill="Expression")+
        theme(
                panel.grid = element_blank(),
                panel.background = element_blank(),
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                axis.title = element_blank(),
                plot.background = element_blank() #the key to avoide legend overlap
        ) +
        scale_fill_gradient2(
                low = "blue",
                mid = "white",
                high = "red"
        ) +
        scale_x_continuous(expand = c(0,0))
p3

relative_heights=c(1,0.06,1)
library(aplot)

p1 %>% 
        insert_bottom(p2,height = relative_heights[1])%>% 
        insert_bottom(middle,height=relative_heights[2]) %>%
        insert_bottom(p3,height=relative_heights[3])

ggsave(file = paste0("04ROC/",outputfile,"_riskCurve.pdf") , width = 5,height = 5)
