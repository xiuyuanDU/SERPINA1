library(IOBR)
library(limma)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsci)
setwd('F:/TCGA/workpath/04.CIBERSORT/Lymphocytic Thyroiditis')
expFile="HiSeqV2.txt"     # Input expression matrix file name
sampleFile="sample(HT).txt"    # Phenotype file

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=normalizeBetweenArrays(mat)
mat=mat[rowMeans(mat)>0,]
afcibersort<-deconvo_tme(mat, method = "cibersort", perm = 1000)
write.table(afcibersort,"13.CIBERSORT.txt", sep="\t", quote = F, row.names = F)

#############################################################################################
### Cell Proportions
sampleFile="sample(HT).txt"    # Phenotype file
cibersort=as.data.frame(afcibersort)
rownames(cibersort)=cibersort[,1]
immune=cibersort[,-1]
immune=immune[immune[,"P-value_CIBERSORT"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
colnames(immune)=gsub("_CIBERSORT"," ",colnames(immune))
colnames(data)=gsub("_CIBERSORT"," ",colnames(data))
colnames(immune)=gsub("_"," ",colnames(immune))
colnames(data)=gsub("_"," ",colnames(data))
cluster=read.table(sampleFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(cluster)="Type"
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$Type),]
gaps=c(1, as.vector(cumsum(table(data$Type))))
xlabels=levels(factor(data$Type))
data$Type=factor(data$Type, levels=c(cluster[2,1],cluster[nrow(cluster),1]))
data$GSM=rownames(data)
data=melt(data, id.vars=c("Type","GSM"))
colnames(data)=c("Type","GSM","Celltype", "Freq")
Cellratio=data
colourCount = length(unique(Cellratio$Celltype))
# Define colors for cell infiltration
colaa=col=c(pal_aaas()(10),pal_jco()(10),pal_nejm()(8),pal_d3()(10),pal_jama()(6))
ggplot(Cellratio) + 
  geom_bar(aes(x =GSM, y= Freq, fill = Celltype), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x='Cell cycle phase', y = 'Ratio') +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"), legend.position = "right") +   # Legend: "left" left, "right" right, "bottom" bottom, "top" top
  scale_fill_manual(values=colaa) +
  theme_bw() +
  xlab(NULL) +
  theme(axis.text.x  = element_blank()) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  facet_grid(. ~ Type, scales="free") +
  theme(strip.text.x = element_text(size = 20, colour = "white")) +       # Facet font color
  theme(strip.background.x = element_rect(fill = c("grey"), colour = "black")) # Facet color
ggsave("4.1.immune.ration.pdf", width = 15, height = 8)        # Output image

### Box plot of cell proportion differences
cibersort=as.data.frame(afcibersort)
rownames(cibersort)=cibersort[,1]
immune=cibersort[,-1]
immune=immune[immune[,"P-value_CIBERSORT"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
colnames(immune)=gsub("_CIBERSORT"," ",colnames(immune))
colnames(data)=gsub("_CIBERSORT"," ",colnames(data))
colnames(immune)=gsub("_"," ",colnames(immune))
colnames(data)=gsub("_"," ",colnames(data))
cluster=read.table(sampleFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(cluster)="Type"
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$Type),]
gaps=c(1, as.vector(cumsum(table(data$Type))))
xlabels=levels(factor(data$Type))
data=melt(data, id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")
group=levels(factor(data$Type))
data$Type=factor(data$Type, levels=c(cluster[2,1],cluster[nrow(cluster),1]))
bioCol=pal_jco()(6)
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Type",
                  xlab="",
                  ylab="CIBERSORT Fraction",
                  legend.title="Type", 
                  width=0.8,
                  palette=bioCol, add.params = list(size=0.1))
boxplot=boxplot +
  stat_compare_means(aes(group=Type), symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif") +
  theme_bw() +
  rotate_x_text(50) +
  labs(title = "Lymphocytic thyroiditis")

pdf(file="4.1.immune.diff.pdf", width=9, height=4.5)
print(boxplot)
dev.off()

#################################################################################################
setwd('F:/TCGA/workpath/04.CIBERSORT/nomal')
sampleFile="sample（Nomal）.txt"    # Change phenotype file (with header, sorted)
### Cell Proportions
cibersort=as.data.frame(afcibersort)
rownames(cibersort)=cibersort[,1]
immune=cibersort[,-1]
immune=immune[immune[,"P-value_CIBERSORT"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
colnames(immune)=gsub("_CIBERSORT"," ",colnames(immune))
colnames(data)=gsub("_CIBERSORT"," ",colnames(data))
colnames(immune)=gsub("_"," ",colnames(immune))
colnames(data)=gsub("_"," ",colnames(data))
cluster=read.table(sampleFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(cluster)="Type"
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$Type),]
gaps=c(1, as.vector(cumsum(table(data$Type))))
xlabels=levels(factor(data$Type))
data$Type=factor(data$Type, levels=c(cluster[2,1],cluster[nrow(cluster),1]))
data$GSM=rownames(data)
data=melt(data, id.vars=c("Type","GSM"))
colnames(data)=c("Type","GSM","Celltype", "Freq")
Cellratio=data
colourCount = length(unique(Cellratio$Celltype))
# Define colors for cell infiltration
colaa=col=c(pal_aaas()(10),pal_jco()(10),pal_nejm()(8),pal_d3()(10),pal_jama()(6))
ggplot(Cellratio) + 
  geom_bar(aes(x =GSM, y= Freq, fill = Celltype), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x='Cell cycle phase', y = 'Ratio') +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"), legend.position = "right") +   # Legend: "left" left, "right" right, "bottom" bottom, "top" top
  scale_fill_manual(values=colaa) +
  theme_bw() +
  xlab(NULL) +
  theme(axis.text.x  = element_blank()) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  facet_grid(. ~ Type, scales="free") +
  theme(strip.text.x = element_text(size = 20, colour = "white")) +       # Facet font color
  theme(strip.background.x = element_rect(fill = c("grey"), colour = "black")) # Facet color
ggsave("4.2.immune.ration.pdf", width = 15, height = 8)        # Output image

### Box plot of cell proportion differences
cibersort=as.data.frame(afcibersort)
rownames(cibersort)=cibersort[,1]
immune=cibersort[,-1]
immune=immune[immune[,"P-value_CIBERSORT"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
colnames(immune)=gsub("_CIBERSORT"," ",colnames(immune))
colnames(data)=gsub("_CIBERSORT"," ",colnames(data))
colnames(immune)=gsub("_"," ",colnames(immune))
colnames(data)=gsub("_"," ",colnames(data))
cluster=read.table(sampleFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(cluster)="Type"
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$Type),]
gaps=c(1, as.vector(cumsum(table(data$Type))))
xlabels=levels(factor(data$Type))
data=melt(data, id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")
group=levels(factor(data$Type))
data$Type=factor(data$Type, levels=c(cluster[2,1],cluster[nrow(cluster),1]))
bioCol=pal_jco()(6)
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Type",
                  xlab="",
                  ylab="CIBERSORT Fraction",
                  legend.title="Type", 
                  width=0.8,
                  palette=bioCol, add.params = list(size=0.1))
boxplot=boxplot +
  stat_compare_means(aes(group=Type), symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif") +
  theme_bw() +
  rotate_x_text(50) +
  labs(title = "No Thyroiditis")

pdf(file="4.2.immune.diff.pdf", width=9, height=4.5)
print(boxplot)
dev.off()

##############################################################################################
setwd('F:/TCGA/workpath/04.CIBERSORT')
### Cell Proportions
sampleFile="sample.txt"    # Phenotype file
### Cell Proportions
cibersort=as.data.frame(afcibersort)
rownames(cibersort)=cibersort[,1]
immune=cibersort[,-1]
immune=immune[immune[,"P-value_CIBERSORT"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
colnames(immune)=gsub("_CIBERSORT"," ",colnames(immune))
colnames(data)=gsub("_CIBERSORT"," ",colnames(data))
colnames(immune)=gsub("_"," ",colnames(immune))
colnames(data)=gsub("_"," ",colnames(data))
cluster=read.table(sampleFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(cluster)="Type"
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$Type),]
gaps=c(1, as.vector(cumsum(table(data$Type))))
xlabels=levels(factor(data$Type))
data$Type=factor(data$Type, levels=c(cluster[2,1],cluster[nrow(cluster),1]))
data$GSM=rownames(data)
data=melt(data, id.vars=c("Type","GSM"))
colnames(data)=c("Type","GSM","Celltype", "Freq")
Cellratio=data
colourCount = length(unique(Cellratio$Celltype))
# Define colors for cell infiltration
colaa=col=c(pal_aaas()(10),pal_jco()(10),pal_nejm()(8),pal_d3()(10),pal_jama()(6))
ggplot(Cellratio) + 
  geom_bar(aes(x =GSM, y= Freq, fill = Celltype), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x='Cell cycle phase', y = 'Ratio') +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"), legend.position = "right") +   # Legend: "left" left, "right" right, "bottom" bottom, "top" top
  scale_fill_manual(values=colaa) +
  theme_bw() +
  xlab(NULL) +
  theme(axis.text.x  = element_blank()) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  facet_grid(. ~ Type, scales="free") +
  theme(strip.text.x = element_text(size = 20, colour = "white")) +       # Facet font color
  theme(strip.background.x = element_rect(fill = c("grey"), colour = "black")) # Facet color
ggsave("4.2.immune.ration.pdf", width = 15, height = 8)        # Output image

### Cell Proportion Difference Boxplot
cibersort = as.data.frame(afcibersort)
rownames(cibersort) = cibersort[, 1]
immune = cibersort[, -1]
immune = immune[immune[, "P-value_CIBERSORT"] < 0.05, ]
data = as.matrix(immune[, 1:(ncol(immune) - 3)])
colnames(immune) = gsub("CIBERSORT", " ", colnames(immune))
colnames(data) = gsub("CIBERSORT", " ", colnames(data))
colnames(immune) = gsub("", " ", colnames(immune))
colnames(data) = gsub("", " ", colnames(data))
cluster = read.table(sampleFile, header = FALSE, sep = "\t", check.names = FALSE, row.names = 1)
colnames(cluster) = "Type"
sameSample = intersect(row.names(data), row.names(cluster))
data = cbind(data[sameSample, , drop = FALSE], cluster[sameSample, , drop = FALSE])
data = data[order(data$Type), ]
gaps = c(1, as.vector(cumsum(table(data$Type))))
xlabels = levels(factor(data$Type))
data = melt(data, id.vars = c("Type"))
colnames(data) = c("Type", "Immune", "Expression")
group = levels(factor(data$Type))
data$Type = factor(data$Type, levels = c(cluster[2, 1], cluster[nrow(cluster), 1]))
bioCol = pal_jco()(6)
bioCol = bioCol[1:length(group)]
boxplot = ggboxplot(data, x = "Immune", y = "Expression", fill = "Type",
                    xlab = "",
                    ylab = "CIBERSORT Fraction",
                    legend.title = "Type",
                    width = 0.8,
                    palette = bioCol, add.params = list(size = 0.1))
boxplot = boxplot +
  stat_compare_means(aes(group = Type), symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("", "", "", "ns")), label = "p.signif") +
  theme_bw() +
  rotate_x_text(50) +
  labs(title = "All Samples")

pdf(file = "4.3.immune.diff.pdf", width = 9, height = 4.5)
print(boxplot)
dev.off()

