### Figure 8
# the following files are provided in OSF https://osf.io/9xys4/
Counts <- read.delim("Counts.txt",h=T)
Diff_ctr <- read.delim("Diff_ctr.txt",h=T)
colData_RNA_Diff <- data.frame("Sample"=colnames(Diff_ctr)[3:31],
                               "Donor"=gsub("_.*","",colnames(Diff_ctr)[3:31]),
                               "Timepoint"=gsub(".*_","",colnames(Diff_ctr)[3:31]))
rownames(colData_RNA_Diff)<-colData_RNA_Diff$Sample


### Figure 8A
library(Seurat)
library(gplots)
# the following files are provided in OSF https://osf.io/9xys4/
scOC <- readRDS("ReadyToUse_scOC.rds")

# Expression of resorption associated genes per cell
Matrix <- data.frame(as.matrix(GetAssayData(scOC)))
# upregulated
tmp <- Diff_ctr[Diff_ctr$padj_Resorption_d0 < 0.01 & Diff_ctr$logFC_Resorption_d0 > 0,'Symbol']
tmp <- tmp[tmp %in% rownames(scOC)]
# add expression sum to seurat object
scOC$d0_res_up <- colSums(Matrix[rownames(Matrix) %in% tmp,])
# downregulated
tmp <- Diff_ctr[Diff_ctr$padj_Resorption_d0 < 0.01 & Diff_ctr$logFC_Resorption_d0 < 0,'Symbol']
tmp <- tmp[tmp %in% rownames(scOC)]
# add expression sum to seurat object
scOC$d0_res_down <- colSums(Matrix[rownames(Matrix) %in% tmp,])
# all
tmp <- Diff_ctr[Diff_ctr$padj_Resorption_d0 < 0.01,'Symbol']
tmp <- tmp[tmp %in% rownames(scOC)]
# add expression sum to seurat object
scOC$d0_res <- colSums(Matrix[rownames(Matrix) %in% tmp,])

FeaturePlot(scOC,c('d0_res_up','d0_res_down','d0_res'))
rm(scOC,tmp,Matrix)


### Figure 8B
# The output of the figure is modified in Illustrator in which the genes for lymphocytes and erythrocytes are moved apart from the rest of the heatmap to highlight that those genes are not part of the genes that are differentially expressed with reorption outcome. 

# Day 0 resorption heatmap to which marker genes for lymphocytes and erythrocytes are added to check if d0 genes reflect contamination of the buffy coats
tmp <- colData_RNA_Diff[colData_RNA_Diff$Timepoint=="d0",]
# the following files are provided in OSF https://osf.io/9xys4/
Activity <- read.delim("OC_activity.txt")
tmp <- merge(tmp,Activity, by="Donor")
#remove Donor 7
tmp <- tmp[tmp$Donor != "Donor7",]
tmp[order(tmp$Resorption),'Sample']

# select genes regulated by resorption at day 0, order donor samples at day 0 according to resorptive activity
y <- rbind(
  Diff_ctr[Diff_ctr$padj_Resorption_d0 < 0.01,c('Symbol','d0',tmp[order(tmp$Resorption),'Sample'])],
  Diff_ctr[Diff_ctr$Symbol %in% c("HBA2","HBA1","HBB","CD3E","CD5","LEF1"),c('Symbol','d0',tmp[order(tmp$Resorption),'Sample'])])

rownames(y) <- y$Symbol
y$Symbol <- NULL

y <- cbind(y[,1],t(scale(t(y[,2:6]))))

library(fields)
library(RColorBrewer)
library(gplots)

# color for scaled expression values
Mycol2 <- rev(designer.colors(n=100,col=brewer.pal(9, "Spectral")))
# color for absolute expression at day 0
Mycol <- designer.colors(n=150,col=brewer.pal(9, "Greens"))

heatp <- heatmap.2(as.matrix(y[,2:6]), col=Mycol2,Colv = F,trace='none',labCol= colnames(y),labRow = rownames(y),
                   RowSideColors = Mycol[round((y[,1]+2)*9)],
                   ColSideColors = designer.colors(100,col=brewer.pal(9,"BuPu"))[round(tmp[order(tmp$Resorption),'Resorption']*10)])

# cellular localization of those factors according to data from Human Protein Atlas
# the following files are provided in OSF https://osf.io/9xys4/
Localization <- read.delim("subcellular_location.tsv", h=T)
# Get the Human protein atals information for the genes of heatmap
y3 <- data.frame("Symbol"=rownames(y))
y3 <- merge(y3, Localization[,c('Gene.name','Reliability','Main.location')], by.x="Symbol", by.y="Gene.name", all.x=T)
rownames(y3) <- y3$Symbol
# Split up if there are several Main locations
y3$Main.2 <- NA
for (i in 1:nrow(y3)){
  if(grepl(";",y3[i,"Main.location"])){
    y3[i,'Main.2'] <- gsub(".*;","",y3[i,"Main.location"])
    y3[i,'Main.location'] <- gsub(";.*","",y3[i,"Main.location"])
  } else{}
}
# remove NA
y3[!complete.cases(y3$Main.location),'Main.location'] <- 'unknown'
y3[!complete.cases(y3$Main.2),'Main.2'] <- 'unknown'


# combine "Nuclear bodies", "Nucleoli", "Nucleoplasm" into "Nuclear"
# combine "Actin filaments", "Intermediate filaments", "Microtubules" into "Cytoskeleton"

y3[y3$Main.location=="Actin filaments",'Main.location'] <- 'Cytoskeleton'
y3[y3$Main.location=="Intermediate filaments",'Main.location'] <- 'Cytoskeleton'
y3[y3$Main.location=="Microtubules",'Main.location'] <- 'Cytoskeleton'
y3[y3$Main.location=="Nucleoli",'Main.location'] <- 'Nuclear'
y3[y3$Main.location=="Nucleoplasm",'Main.location'] <- 'Nuclear'
y3[y3$Main.location=="Nuclear bodies",'Main.location'] <- 'Nuclear'

y3[y3$Main.2=="Actin filaments",'Main.2'] <- 'Cytoskeleton'
y3[y3$Main.2=="Intermediate filaments",'Main.2'] <- 'Cytoskeleton'
y3[y3$Main.2=="Microtubules",'Main.2'] <- 'Cytoskeleton'
y3[y3$Main.2=="Nucleoli",'Main.2'] <- 'Nuclear'
y3[y3$Main.2=="Nucleoplasm",'Main.2'] <- 'Nuclear'
y3[y3$Main.2=="Nuclear bodies",'Main.2'] <- 'Nuclear'

# Make levels for the heatmap
y3$Main.location <- factor(y3$Main.location,levels=c("Cytoskeleton","Cytosol","Endoplasmic reticulum","Golgi apparatus","Nuclear","Plasma membrane","Vesicles","unknown"))
y3$Main.2 <- factor(y3$Main.2,levels=c("Cytoskeleton","Cytosol","Endoplasmic reticulum","Golgi apparatus","Nuclear","Plasma membrane","Vesicles","unknown"))

# order data frame with human protein atlas according to the heatmap
y3 <- y3[rownames(y)[rev(heatp$rowInd)],]

heatmap.2(cbind(as.numeric(as.factor(y3$Main.location)),
                as.numeric(as.factor(y3$Main.2))), col=c(brewer.pal(7,"Set1"),'white'),Rowv = F,Colv = F, trace='none',labCol= colnames(y3),labRow = rownames(y3))

heatmap.2(as.matrix(y[,2:6]), col=Mycol2,Colv = F, trace='none',labCol= colnames(y),labRow = rownames(y),
          RowSideColors = c(brewer.pal(7,"Set1"),'white')[as.numeric(as.factor(y3$Main.2))])

# Helpful plot to link Rowsidecolors above to categories
plot(as.numeric(unique(as.factor(y3$Main.location))),pch=21,bg=c(brewer.pal(7,"Set1"),'white')[as.numeric(unique(as.factor(y3$Main.location)))])
text(1:8,as.numeric(unique(as.factor(y3$Main.location))),unique(y3$Main.location))
rm(tmp,y,y3,Localization,Mycol,Mycol2,heatp,Activity)


### Figure 8C
# the following files are provided in OSF https://osf.io/9xys4/
Activity <- read.delim("OC_activity.txt")

par(mfrow=c(2,2), pty="s")
GOI <- c('TNK2','CLEC5A','FLNB','OLR1')
for (i in GOI){
  # extract expression levels of GOI
  tmp <- cbind(colData_RNA_Diff[colData_RNA_Diff$Timepoint=="d0",],
               t(Counts[Counts$Symbol == i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d0','Sample']]))
  names(tmp)[4] <- i
  tmp <- merge(tmp,Activity, by="Donor")
  
  plot(tmp[,i],tmp$Resorption, main=paste(i,"padj d0:",Diff_ctr[Diff_ctr$Symbol==i,"padj_Resorption_d0"]))
  abline(lm(tmp$Resorption ~ tmp[,i]))
}


### Figure 8F
# the following files are provided in OSF https://osf.io/9xys4/
FACS <- read.delim("OC_Resorption_FACS.txt")

par(mfrow=c(2,4), pty="s")
plot(FACS[,"Percent_ACK1"],FACS[,"Resorption"])
abline(lm(FACS[,"Resorption"]~FACS[,"Percent_ACK1"]))
title(summary(lm(FACS[,"Resorption"]~FACS[,"Percent_ACK1"]))$coefficients[2,4])
plot(FACS[,"Percent_CLEC5a"],FACS[,"Resorption"])
abline(lm(FACS[,"Resorption"]~FACS[,"Percent_CLEC5a"]))
title(summary(lm(FACS[,"Resorption"]~FACS[,"Percent_CLEC5a"]))$coefficients[2,4])
plot(FACS[,"Percent_FLNB"],FACS[,"Resorption"])
abline(lm(FACS[,"Resorption"]~FACS[,"Percent_FLNB"]))
title(summary(lm(FACS[,"Resorption"]~FACS[,"Percent_FLNB"]))$coefficients[2,4])
plot(FACS[,"Percent_LOX1"],FACS[,"Resorption"])
abline(lm(FACS[,"Resorption"]~FACS[,"Percent_LOX1"]))
title(summary(lm(FACS[,"Resorption"]~FACS[,"Percent_LOX1"]))$coefficients[2,4])
plot(FACS[,"MFI_ACK1"],FACS[,"Resorption"])
abline(lm(FACS[,"Resorption"]~FACS[,"MFI_ACK1"]))
title(summary(lm(FACS[,"Resorption"]~FACS[,"MFI_ACK1"]))$coefficients[2,4])
plot(FACS[,"MFI_CLEC5a"],FACS[,"Resorption"])
abline(lm(FACS[,"Resorption"]~FACS[,"MFI_CLEC5a"]))
title(summary(lm(FACS[,"Resorption"]~FACS[,"MFI_CLEC5a"]))$coefficients[2,4])
plot(FACS[,"MFI_FLNB"],FACS[,"Resorption"])
abline(lm(FACS[,"Resorption"]~FACS[,"MFI_FLNB"]))
title(summary(lm(FACS[,"Resorption"]~FACS[,"MFI_FLNB"]))$coefficients[2,4])
plot(FACS[,"MFI_LOX1"],FACS[,"Resorption"])
abline(lm(FACS[,"Resorption"]~FACS[,"MFI_LOX1"]))
title(summary(lm(FACS[,"Resorption"]~FACS[,"MFI_LOX1"]))$coefficients[2,4])

rm(list=ls())
