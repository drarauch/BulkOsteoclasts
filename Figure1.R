### Figure 1
# the following files are provided in OSF https://osf.io/9xys4/
Counts <- read.delim("Counts.txt",h=T)
Diff_ctr <- read.delim("Diff_ctr.txt",h=T)
colData_RNA_Diff <- data.frame("Sample"=colnames(Diff_ctr)[3:31],
                               "Donor"=gsub("_.*","",colnames(RNA_raw)[3:31]),
                               "Timepoint"=gsub(".*_","",colnames(RNA_raw)[3:31]))
rownames(colData_RNA_Diff)<-colData_RNA_Diff$Sample


### Figure 1B
# Graphs were merged in illustrator
data <- read.delim("OC_activity.txt",h=T)
barplot(c(t(data[1,2:3]),
          t(data[2,2:3]),
          t(data[3,2:3]),
          t(data[4,2:3]),
          t(data[5,2:3]),
          t(data[6,2:3]),
          t(data[7,2:3]),
          t(data[8,2:3])),col=c('white','grey'),beside=T)
barplot(data[,4],col=c('black'),beside=T)
rm(data)


### Figure 1C
GOI <- c('TREM1','ADGRE1','SELL','CLEC10A','CTSK','ACP5','MMP9','CA2')
for(i in GOI){
  tmp <-colData_RNA_Diff
  tmp$x <- NA
  tmp[tmp$Timepoint=='d0','x'] <- 1
  tmp[tmp$Timepoint=='d2','x'] <- 2
  tmp[tmp$Timepoint=='d5','x'] <- 3
  tmp[tmp$Timepoint=='d9','x'] <- 4
  d0 <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d0','Sample']])
  d2 <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d2','Sample']])
  d5 <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d5','Sample']])
  d9 <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d9','Sample']])
  boxplot(d0,d2,d5,d9,ylab=paste("Counts",i), xaxt="none")
  axis(1,at=c(1:4),labels = c('D0','D2','D5','D9'))
  for(k in unique(colData_RNA_Diff$Donor)){
    a <- tmp[tmp$Donor==k,'x']
    b <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Donor==k,'Sample']])
    lines(a,b)
  }
}
rm(GOI,tmp,d0,d2,d5,d9,a,b,i,k)


### Figure 1D
Diff_ctr[Diff_ctr$Symbol %in% c('LPCAT2','COX11'),grep("cor_Donor",colnames(Diff_ctr))]

GOI <- c('LPCAT2','COX11')
for(i in GOI){
  tmp <-colData_RNA_Diff
  tmp$x <- NA
  tmp[tmp$Timepoint=='d0','x'] <- 1
  tmp[tmp$Timepoint=='d2','x'] <- 2
  tmp[tmp$Timepoint=='d5','x'] <- 3
  tmp[tmp$Timepoint=='d9','x'] <- 4
  d0 <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d0','Sample']])
  d2 <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d2','Sample']])
  d5 <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d5','Sample']])
  d9 <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d9','Sample']])
  boxplot(d0,d2,d5,d9,ylab=paste("Counts",i), xaxt="none")
  axis(1,at=c(1:4),labels = c('D0','D2','D5','D9'))
  for(k in unique(colData_RNA_Diff$Donor)){
    a <- tmp[tmp$Donor==k,'x']
    b <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Donor==k,'Sample']])
    lines(a,b)
  }
}
rm(GOI,tmp,d0,d2,d5,d9,a,b,i,k)


### Figure 1E
barplot(c(
  nrow(Diff_ctr[Diff_ctr$cor_08==0,]),
  nrow(Diff_ctr[Diff_ctr$cor_08==1,]),
  nrow(Diff_ctr[Diff_ctr$cor_08==2,]),
  nrow(Diff_ctr[Diff_ctr$cor_08==3,]),
  nrow(Diff_ctr[Diff_ctr$cor_08==4,]),
  nrow(Diff_ctr[Diff_ctr$cor_08==5,]),
  nrow(Diff_ctr[Diff_ctr$cor_08==6,]),
  nrow(Diff_ctr[Diff_ctr$cor_08==7,]),
  nrow(Diff_ctr[Diff_ctr$cor_08==8,])),col=c('white','white','white','white','white','white','grey','grey','grey'))


### Figure 1F
# PCA plot based on rlog values in Diff_ctr (not shown in the Figure)
object <- Diff_ctr[,rownames(colData_RNA_Diff)]
intgroup <- c('Timepoint', 'Donor')
intgroup.df <- colData_RNA_Diff
# based on the 5000 most variable genes (arbitrary)
ntop=5000
# rank varying genes and select ntop ones, all below is wrapped in the DEseq2 function plotPCA
library(genefilter)
rv <- rowVars(object)
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(object[select,]))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
if (!all(intgroup %in% names(intgroup.df))) {
  stop("the argument 'intgroup' should specify columns of intgroup.df")
}
rownames(intgroup.df) <- intgroup.df$Sample
group <- if (length(intgroup) > 1) {
  factor(apply( intgroup.df, 1, paste, collapse=" : "))
} else {
  intgroup.df[[intgroup]]
}
dat <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3], PC4=pca$x[,4], group=group, intgroup.df, name=colnames(object))
attr(dat, "percentVar") <- percentVar[1:4]

# scale axis to get a squared graph
ratio.values <- (max(dat$PC1)-min(dat$PC1))/(max(dat$PC2)-min(dat$PC2))
#plot results
library(ggplot2)
ggplot(data=dat, aes_string(x="PC1", y="PC2", color="Timepoint"))+ geom_point(size=3) + coord_fixed(ratio.values) +
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance"))


## this one below is plotted in the Figure and is almost identical to the unbiased above
# PCA on selected genes, reproducible across Donors (at least 6 Donors with a pearson's correlation of more than 0.8 to the average) and significant for at least one timepoint with FDR < 0.0001
object <- Diff_ctr[apply(Diff_ctr[,c('padj_d2_vs_d0','padj_d5_vs_d0','padj_d9_vs_d0','padj_d5_vs_d2','padj_d9_vs_d2','padj_d9_vs_d5')],1,min)<0.0001 & Diff_ctr$cor_08 > 5,rownames(colData_RNA_Diff)]
intgroup <- c('Timepoint', 'Donor')
intgroup.df <- colData_RNA_Diff
# based on the 5000 most variable genes (arbitrary)
ntop=5000

# rank varying genes and select ntop ones, all below is wrapped in the DEseq2 function plotPCA
library(genefilter)
rv <- rowVars(object)
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(object[select,]))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
if (!all(intgroup %in% names(intgroup.df))) {
  stop("the argument 'intgroup' should specify columns of intgroup.df")
}
rownames(intgroup.df) <- intgroup.df$Sample
group <- if (length(intgroup) > 1) {
  factor(apply( intgroup.df, 1, paste, collapse=" : "))
} else {
  intgroup.df[[intgroup]]
}
dat <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3], PC4=pca$x[,4], group=group, intgroup.df, name=colnames(object))
attr(dat, "percentVar") <- percentVar[1:4]

# scale axis to get a squared graph
ratio.values <- (max(dat$PC1)-min(dat$PC1))/(max(dat$PC2)-min(dat$PC2))
#plot results
library(ggplot2)
ggplot(data=dat, aes_string(x="PC1", y="PC2", color="Timepoint"))+ geom_point(size=3) + coord_fixed(ratio.values) +
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance"))

rm(group, intgroup, intgroup.df,ntop,object,dat,ratio.values,percentVar,select,rv,pca)


###Figure 1G
# the following files are provided in OSF https://osf.io/9xys4/
data <- read.delim("ReadyToUse_GSE225974.txt",h=T)
tmp <- merge(Diff_ctr[,c("Symbol","logFC_d2_vs_d0","logFC_d5_vs_d0", "logFC_d9_vs_d0","logFC_d5_vs_d2","logFC_d9_vs_d2","logFC_d9_vs_d5")],data[,c("Symbol","padj_OCvsPBMC","logFC_OCvsPBMC")], by="Symbol")

# select genes based on published data
tmp_down <- tmp[tmp$padj_OCvsPBMC < 0.01 & tmp$logFC_OCvsPBMC <0,]
tmp_up <- tmp[tmp$padj_OCvsPBMC < 0.01 & tmp$logFC_OCvsPBMC >0,]

# calculate Pearson's correlation for the log fold changes of the published study (one comparison) against the various combinations of log fold changes from this study (6 comparisons)
y <- cbind(
  cor(tmp_down[,c("logFC_d2_vs_d0","logFC_d5_vs_d0", "logFC_d9_vs_d0","logFC_d5_vs_d2","logFC_d9_vs_d2","logFC_d9_vs_d5","logFC_OCvsPBMC")])[,7],
  cor(tmp_up[,c("logFC_d2_vs_d0","logFC_d5_vs_d0", "logFC_d9_vs_d0","logFC_d5_vs_d2","logFC_d9_vs_d2","logFC_d9_vs_d5","logFC_OCvsPBMC")])[,7])

# plot results in a heatmap not including the published study against itself
library(fields)
mat_col <- designer.colors(50,col=c('blue','white','red'))
mat_col_breaks <- seq(-1,1,length=51)  
library(gplots)
heatmap.2(y[1:6,],col=mat_col, breaks=mat_col_breaks, trace="none", Colv = F, Rowv = F)
rm(tmp_down, tmp_up,mat_col, mat_col_breaks,tmp, data,y)


### Figure 1H
# the following files are provided in OSF https://osf.io/9xys4/
data <- read.delim("ReadyToUse_GSE225974.txt",h=T)

GOI <- c('SLC6A7')
par(mfrow=c(1,2),pty="s")
for(i in GOI){
  tmp <-colData_RNA_Diff
  tmp$x <- NA
  tmp[tmp$Timepoint=='d0','x'] <- 1
  tmp[tmp$Timepoint=='d2','x'] <- 2
  tmp[tmp$Timepoint=='d5','x'] <- 3
  tmp[tmp$Timepoint=='d9','x'] <- 4
  d0 <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d0','Sample']])
  d2 <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d2','Sample']])
  d5 <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d5','Sample']])
  d9 <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d9','Sample']])
  boxplot(d0,d2,d5,d9,ylab=paste("Counts",i), xaxt="none")
  axis(1,at=c(1:4),labels = c('D0','D2','D5','D9'))
  for(k in unique(colData_RNA_Diff$Donor)){
    a <- tmp[tmp$Donor==k,'x']
    b <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Donor==k,'Sample']])
    lines(a,b)
  }
  tmp <- t(data[data$Symbol==i,grep("counts",colnames(data))])
  tmp <- cbind(tmp[c(1,3,5,7,9,11,13,15),],tmp[c(2,4,6,8,10,12,14,16),])
  boxplot(tmp, ylab=paste("Counts",i), xaxt="none")
  axis(1,at=c(1:2),labels = c('PBMC','OC-like'))
  for(k in 1:dim(tmp)[1]){
    lines(c(1,2),c(tmp[k,1:2]))
  }
}

rm(GOI,tmp,d0,d2,d5,d9,a,b,i,k,data)

### Figure 1I
# the following files are provided in OSF https://osf.io/9xys4/
data <- read.delim("ReadyToUse_GSE225974.txt",h=T)
tmp <- merge(Diff_ctr[,c("Symbol","logFC_d2_vs_d0","logFC_d5_vs_d0", "logFC_d9_vs_d0","logFC_d5_vs_d2","logFC_d9_vs_d2","logFC_d9_vs_d5")],data[,c("Symbol","padj_OCvsPBMC","logFC_OCvsPBMC")], by="Symbol")

# select upregulated genes based on published data
tmp_up <- tmp[tmp$padj_OCvsPBMC < 0.01 & tmp$logFC_OCvsPBMC >0,]
plot(tmp_up$logFC_OCvsPBMC,tmp_up$logFC_d9_vs_d2, col=alpha('black',0.3))
rm(tmp_up, data,tmp)


