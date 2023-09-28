#### Data preparation
 
# For reproducing the individual figures processed data files are available at OSF https://osf.io/9xys4/
# This part is needed if you want to start from scratch to make the processed data files

#### DEseq2 analysis of RNA-seq data
library(DESeq2)
## Input, homer based raw counts, day 0 sample for Donor 5, 7, and 8 had no counts due to low RNA yield
RNA_raw <- read.delim("RNAcounts_Differentiation.txt", h=T)

# Remove very low expressed transcripts
RNA_raw <- RNA_raw[apply(RNA_raw[,9:37],1,max)> 5,]

# Extract gene symbol anotation from homer annotated expression file
RNA_raw$Symbol <- NA
for(i in 1:nrow(RNA_raw)){
  RNA_raw[i,'Symbol'] <- unlist(strsplit(as.character(RNA_raw[i,"Annotation.Divergence"]), '\\|'))[1]
}

# Set up colData, using information from sample names
colData_RNA_Diff <- data.frame("Sample"=colnames(RNA_raw)[9:37],
                          "Donor"=gsub("_.*","",colnames(RNA_raw)[9:37]),
                          "Timepoint"=gsub(".*_","",colnames(RNA_raw)[9:37]))
rownames(colData_RNA_Diff)<-colData_RNA_Diff$Sample

# Set up countData
countData_RNA<-as.matrix(RNA_raw[,c(9:37)])
row.names(countData_RNA)<-RNA_raw$RefSeqID

# Make DEseq 2 object, use Timepoint and Donor in the design matrix
dds_Diff <- DESeqDataSetFromMatrix(countData = countData_RNA , colData = colData_RNA_Diff, design = ~Timepoint + Donor)
dds_Diff <- DESeq(dds_Diff)

# Extract normalized and log-transformed expression levels for PCA plot and heatmaps
rld_RNA_Diff <- rlog(dds_Diff, blind=F)

# Extract normalized counts for plotting of genes
Counts <- counts(dds_Diff, normalize=T)
Counts <- cbind(RNA_raw[,c('Symbol','RefSeqID',"Length")],Counts)
# "Counts.txt" is provided at OSF https://osf.io/9xys4/
write.table(Counts,file="Counts.txt",sep="\t",col.names=T,row.names=F,quote=F)

# Extract p-values and log fold changes over time
Diff_ctr <- cbind(RNA_raw[,c('Symbol','RefSeqID')],assay(rld_RNA_Diff))
tmp <- cbind(as.data.frame (results(dds_Diff, contrast = c('Timepoint','d2','d0')))[6],
             as.data.frame (results(dds_Diff, contrast = c('Timepoint','d5','d0')))[6],
             as.data.frame (results(dds_Diff, contrast = c('Timepoint','d9','d0')))[6],
             as.data.frame (results(dds_Diff, contrast = c('Timepoint','d5','d2')))[6],
             as.data.frame (results(dds_Diff, contrast = c('Timepoint','d9','d2')))[6],
             as.data.frame (results(dds_Diff, contrast = c('Timepoint','d9','d5')))[6],
             as.data.frame (results(dds_Diff, contrast = c('Timepoint','d2','d0')))[6],
             as.data.frame (results(dds_Diff, contrast = c('Timepoint','d5','d0')))[2],
             as.data.frame (results(dds_Diff, contrast = c('Timepoint','d9','d0')))[2],
             as.data.frame (results(dds_Diff, contrast = c('Timepoint','d5','d2')))[2],
             as.data.frame (results(dds_Diff, contrast = c('Timepoint','d9','d2')))[2],
             as.data.frame (results(dds_Diff, contrast = c('Timepoint','d9','d5')))[2])
names(tmp) <- c('padj_d2_vs_d0','padj_d5_vs_d0','padj_d9_vs_d0','padj_d5_vs_d2','padj_d9_vs_d2','padj_d9_vs_d5','logFC_d2_vs_d0','logFC_d5_vs_d0','logFC_d9_vs_d0','logFC_d5_vs_d2','logFC_d9_vs_d2','logFC_d9_vs_d5')

# Replace NA in padj columnns with 1 and NA in logFC columns with 0
tmp[!complete.cases(tmp$padj_d2_vs_d0),'padj_d2_vs_d0'] <- 1
tmp[!complete.cases(tmp$padj_d5_vs_d0),'padj_d5_vs_d0'] <- 1
tmp[!complete.cases(tmp$padj_d9_vs_d0),'padj_d9_vs_d0'] <- 1
tmp[!complete.cases(tmp$padj_d5_vs_d2),'padj_d5_vs_d2'] <- 1
tmp[!complete.cases(tmp$padj_d9_vs_d2),'padj_d9_vs_d2'] <- 1
tmp[!complete.cases(tmp$padj_d9_vs_d5),'padj_d9_vs_d5'] <- 1

tmp[!complete.cases(tmp$logFC_d2_vs_d0),'logFC_d2_vs_d0'] <- 0
tmp[!complete.cases(tmp$logFC_d5_vs_d0),'logFC_d5_vs_d0'] <- 0
tmp[!complete.cases(tmp$logFC_d9_vs_d0),'logFC_d9_vs_d0'] <- 0
tmp[!complete.cases(tmp$logFC_d5_vs_d2),'logFC_d5_vs_d2'] <- 0
tmp[!complete.cases(tmp$logFC_d9_vs_d2),'logFC_d9_vs_d2'] <- 0
tmp[!complete.cases(tmp$logFC_d9_vs_d5),'logFC_d9_vs_d5'] <- 0

Diff_ctr <- cbind(Diff_ctr,tmp)
rm(tmp,rld_RNA_Diff,i,dds_Diff)

# Before clustering remove gene that do not show similar temporal profiles for the individual donors and genes that have high expression levels at single time points for single donors
# Check how good the individual donor correlates with expression of the mean
# Calculate mean for the timepoints
Diff_ctr$d0 <- apply(Diff_ctr[,grep('_d0',colnames(Diff_ctr)[1:31])],1,mean)
Diff_ctr$d2 <- apply(Diff_ctr[,grep('_d2',colnames(Diff_ctr)[1:31])],1,mean)
Diff_ctr$d5 <- apply(Diff_ctr[,grep('_d5',colnames(Diff_ctr)[1:31])],1,mean)
Diff_ctr$d9 <- apply(Diff_ctr[,grep('_d9',colnames(Diff_ctr)[1:31])],1,mean)


Diff_ctr$cor_Donor1 <- apply(cbind(Diff_ctr[,c('d0','d2','d5','d9')],Diff_ctr[,grep("Donor1",colnames(Diff_ctr)[1:31])]),1,function(x){
  a <- x[1:4]
  b <- x[5:8]
  return(cor(a,b))
})
Diff_ctr$cor_Donor2 <- apply(cbind(Diff_ctr[,c('d0','d2','d5','d9')],Diff_ctr[,grep("Donor2",colnames(Diff_ctr)[1:31])]),1,function(x){
  a <- x[1:4]
  b <- x[5:8]
  return(cor(a,b))
})
Diff_ctr$cor_Donor3 <- apply(cbind(Diff_ctr[,c('d0','d2','d5','d9')],Diff_ctr[,grep("Donor3",colnames(Diff_ctr)[1:31])]),1,function(x){
  a <- x[1:4]
  b <- x[5:8]
  return(cor(a,b))
})
Diff_ctr$cor_Donor4 <- apply(cbind(Diff_ctr[,c('d0','d2','d5','d9')],Diff_ctr[,grep("Donor4",colnames(Diff_ctr)[1:31])]),1,function(x){
  a <- x[1:4]
  b <- x[5:8]
  return(cor(a,b))
})
Diff_ctr$cor_Donor5 <- apply(cbind(Diff_ctr[,c('d2','d5','d9')],Diff_ctr[,grep("Donor5",colnames(Diff_ctr)[1:31])]),1,function(x){
  a <- x[1:3]
  b <- x[4:6]
  return(cor(a,b))
})
Diff_ctr$cor_Donor6 <- apply(cbind(Diff_ctr[,c('d0','d2','d5','d9')],Diff_ctr[,grep("Donor6",colnames(Diff_ctr)[1:31])]),1,function(x){
  a <- x[1:4]
  b <- x[5:8]
  return(cor(a,b))
})
Diff_ctr$cor_Donor7 <- apply(cbind(Diff_ctr[,c('d2','d5','d9')],Diff_ctr[,grep("Donor7",colnames(Diff_ctr)[1:31])]),1,function(x){
  a <- x[1:3]
  b <- x[4:6]
  return(cor(a,b))
})
Diff_ctr$cor_Donor8 <- apply(cbind(Diff_ctr[,c('d2','d5','d9')],Diff_ctr[,grep("Donor8",colnames(Diff_ctr)[1:31])]),1,function(x){
  a <- x[1:3]
  b <- x[4:6]
  return(cor(a,b))
})

# Check for how many donors is the correlation greater than 0.8
Diff_ctr$cor_08 <- rowSums(Diff_ctr[,grepl("cor_D",colnames(Diff_ctr))]>0.8)


# Cluster the gene expression profiles - THIS IS NON DETERMINISTIC and can lead to slightly different results when re-running
y <- Diff_ctr[apply(Diff_ctr[,c('padj_d2_vs_d0','padj_d5_vs_d0','padj_d9_vs_d0','padj_d5_vs_d2','padj_d9_vs_d2','padj_d9_vs_d5')],1,min)<0.0001,]
y <- y[y$cor_08 > 5, ]
rownames(y) <- y$RefSeqID

# kmeans clustering with transparant colors for publication graph, make many clusters and group them later manually
addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

# Determine 25 clusters
no.cluster <- 25
y2 <- as.matrix(t(scale(t(y[,c("d0","d2","d5","d9")]))))
rownames(y2) <- y$RefSeqID
km <- kmeans(y2, no.cluster, iter.max=60)

library(RColorBrewer)
library(scales)
par(mfrow=c(5,5),pty='s')

# Plot cluster patterns
for(i in 1:no.cluster) {
  tmp1 <- y2[km$cluster == i,1:4]
  Mycol <- addTrans(rainbow(no.cluster)[i],255*20 / dim(tmp1)[1])
  plot(0,0,pch=' ', xlim=c(1,10), ylim=c(-1.5,1.5), xlab="", ylab="", yaxt="n", xaxt="n", main=paste(dim(tmp1)[1]))
  mtext("scaled Signal", side=2, line=2.5, cex=0.5)
  for(j in 1:dim(tmp1)[1]){
    lines(c(1,3,6,10), tmp1[j,], col=Mycol)
  }
  lines(c(1,3,6,10), apply(tmp1,2, function(x){median(x,na.rm=T)}), lwd=3)
  axis(2, at=c(-1,0,1), lab=c(-1,0,1), cex.axis=1)
  axis(1, at=c(1,3,6,10), lab=c("D1","D3","D6","D10"), las=2)
}


y2 <- data.frame(cbind(y2,km$cluster))
names(y2)[5] <- 'cluster'
y2$RefSeqID <- rownames(y2)

# Order and group the clusters into 8 major patterns

y2$clust <- NA
y2[y2$cluster==1,'clust'] <- 8
y2[y2$cluster==2,'clust'] <- 2
y2[y2$cluster==3,'clust'] <- 4
y2[y2$cluster==4,'clust'] <- 3
y2[y2$cluster==5,'clust'] <- 3
y2[y2$cluster==6,'clust'] <- 3
y2[y2$cluster==7,'clust'] <- 5
y2[y2$cluster==8,'clust'] <- 1
y2[y2$cluster==9,'clust'] <- 7
y2[y2$cluster==10,'clust'] <- 6
y2[y2$cluster==11,'clust'] <- 7
y2[y2$cluster==12,'clust'] <- 8
y2[y2$cluster==13,'clust'] <- 5
y2[y2$cluster==14,'clust'] <- 1
y2[y2$cluster==15,'clust'] <- 1
y2[y2$cluster==16,'clust'] <- 6
y2[y2$cluster==17,'clust'] <- 6
y2[y2$cluster==18,'clust'] <- 2
y2[y2$cluster==19,'clust'] <- 4
y2[y2$cluster==20,'clust'] <- 4
y2[y2$cluster==21,'clust'] <- 6
y2[y2$cluster==22,'clust'] <- 7
y2[y2$cluster==23,'clust'] <- 5
y2[y2$cluster==24,'clust'] <- 2
y2[y2$cluster==25,'clust'] <- 7

# Combine cluster information with genes that were not significantly expressed
# Using merge here changes the order in Diff_ctr compared to the DEseq2 object
Diff_ctr <- merge(Diff_ctr,y2[,c("RefSeqID","clust")], by="RefSeqID", all.x=T)
Diff_ctr[!complete.cases(Diff_ctr$clust),'clust'] <- 0

# Reset order to DEseq object
rownames(Diff_ctr) <- Diff_ctr$RefSeqID
Diff_ctr <- Diff_ctr[RNA_raw$RefSeqID,]
rm(i,Mycol,no.cluster, y,y2,j,km)


# Make a new DEseq2 model for resorption but leave out Donor 7 due to the strange resorption activity
# Make DEseq2 model for each day

## day 0
# extract coldata information and merge with Activity table
x <- which(colData_RNA_Diff$Timepoint =="d0")
colData_RNA_Res <- colData_RNA_Diff[x,]
colData_RNA_Res <- merge(colData_RNA_Res, Activity[,c('Donor','Resorption')], by="Donor")
#remove Donor 7
colData_RNA_Res <- colData_RNA_Res[colData_RNA_Res$Donor != "Donor7",]
colData_RNA_Res <- droplevels(colData_RNA_Res)
rownames(colData_RNA_Res) <- colData_RNA_Res$Sample

dds_RNA_Res <- DESeqDataSetFromMatrix(countData = countData_RNA[,colData_RNA_Res$Sample] , colData = colData_RNA_Res, design = ~Resorption)
dds_RNA_Res <- DESeq(dds_RNA_Res)

Result_Res <- cbind(
  as.data.frame (results(dds_RNA_Res, name = "Resorption"))[6],
  as.data.frame (results(dds_RNA_Res, name = "Resorption"))[2])
names(Result_Res) <- c('padj_Resorption_d0','logFC_Resorption_d0')

## day 2
# extract coldata information and merge with Activity table
x <- which(colData_RNA_Diff$Timepoint =="d2")
colData_RNA_Res <- colData_RNA_Diff[x,]
colData_RNA_Res <- merge(colData_RNA_Res, Activity[,c('Donor','Resorption')], by="Donor")
#remove Donor 7
colData_RNA_Res <- colData_RNA_Res[colData_RNA_Res$Donor != "Donor7",]
colData_RNA_Res <- droplevels(colData_RNA_Res)
rownames(colData_RNA_Res) <- colData_RNA_Res$Sample

dds_RNA_Res <- DESeqDataSetFromMatrix(countData = countData_RNA[,colData_RNA_Res$Sample] , colData = colData_RNA_Res, design = ~Resorption)
dds_RNA_Res <- DESeq(dds_RNA_Res)

Result_Res <- cbind(Result_Res,
                    as.data.frame (results(dds_RNA_Res, name = "Resorption"))[6],
                    as.data.frame (results(dds_RNA_Res, name = "Resorption"))[2])
names(Result_Res)[3:4] <- c('padj_Resorption_d2','logFC_Resorption_d2')

## day 5
# extract coldata information and merge with Activity table
x <- which(colData_RNA_Diff$Timepoint =="d5")
colData_RNA_Res <- colData_RNA_Diff[x,]
colData_RNA_Res <- merge(colData_RNA_Res, Activity[,c('Donor','Resorption')], by="Donor")
#remove Donor 7
colData_RNA_Res <- colData_RNA_Res[colData_RNA_Res$Donor != "Donor7",]
colData_RNA_Res <- droplevels(colData_RNA_Res)
rownames(colData_RNA_Res) <- colData_RNA_Res$Sample

dds_RNA_Res <- DESeqDataSetFromMatrix(countData = countData_RNA[,colData_RNA_Res$Sample] , colData = colData_RNA_Res, design = ~Resorption)
dds_RNA_Res <- DESeq(dds_RNA_Res)

Result_Res <- cbind(Result_Res,
                    as.data.frame (results(dds_RNA_Res, name = "Resorption"))[6],
                    as.data.frame (results(dds_RNA_Res, name = "Resorption"))[2])
names(Result_Res)[5:6] <- c('padj_Resorption_d5','logFC_Resorption_d5')

## day 9
# extract coldata information and merge with Activity table
x <- which(colData_RNA_Diff$Timepoint =="d9")
colData_RNA_Res <- colData_RNA_Diff[x,]
colData_RNA_Res <- merge(colData_RNA_Res, Activity[,c('Donor','Resorption')], by="Donor")
#remove Donor 7
colData_RNA_Res <- colData_RNA_Res[colData_RNA_Res$Donor != "Donor7",]
colData_RNA_Res <- droplevels(colData_RNA_Res)
rownames(colData_RNA_Res) <- colData_RNA_Res$Sample

dds_RNA_Res <- DESeqDataSetFromMatrix(countData = countData_RNA[,colData_RNA_Res$Sample] , colData = colData_RNA_Res, design = ~Resorption)
dds_RNA_Res <- DESeq(dds_RNA_Res)

Result_Res <- cbind(Result_Res,
                    as.data.frame (results(dds_RNA_Res, name = "Resorption"))[6],
                    as.data.frame (results(dds_RNA_Res, name = "Resorption"))[2])
names(Result_Res)[7:8] <- c('padj_Resorption_d9','logFC_Resorption_d9')

Result_Res[!complete.cases(Result_Res$padj_Resorption_d0),'padj_Resorption_d0'] <- 1
Result_Res[!complete.cases(Result_Res$padj_Resorption_d2),'padj_Resorption_d2'] <- 1
Result_Res[!complete.cases(Result_Res$padj_Resorption_d5),'padj_Resorption_d5'] <- 1
Result_Res[!complete.cases(Result_Res$padj_Resorption_d9),'padj_Resorption_d9'] <- 1
Result_Res[!complete.cases(Result_Res$logFC_Resorption_d0),'logFC_Resorption_d0'] <- 0
Result_Res[!complete.cases(Result_Res$logFC_Resorption_d2),'logFC_Resorption_d2'] <- 0
Result_Res[!complete.cases(Result_Res$logFC_Resorption_d5),'logFC_Resorption_d5'] <- 0
Result_Res[!complete.cases(Result_Res$logFC_Resorption_d9),'logFC_Resorption_d9'] <- 0

# order Result_Res like the Diff_ctr data frame and bind them together
Result_Res <- Result_Res[Diff_ctr$RefSeqID,]
Diff_ctr <- cbind(Diff_ctr,Result_Res)

# "Diff_ctr.txt" is provided at OSF https://osf.io/9xys4/
write.table(Diff_ctr, file="Diff_ctr.txt", quote=F, col.names=T, row.names=F, sep="\t")
rm(Result_Res)









#### External data sets

### GSE225974 from 10.3390/genes14040916, download raw counts from processed data "GSE225974_gene_counts_PBMC_OC_150219.txt"
# Identification of Differentially Expressed Genes and Molecular Pathways Involved in Osteoclastogenesis Using RNA-seq

data <- read.delim("GSE225974_gene_counts_PBMC_OC_150219.txt",h=T)
rownames(data) <- data$geneid
data <- data[,2:17]
colnames(data) <- c('PBMC_D1','OC_D1','PBMC_D2','OC_D2','PBMC_D3','OC_D3','PBMC_D4','OC_D4','PBMC_D5','OC_D5','PBMC_D6','OC_D6','PBMC_D7','OC_D7','PBMC_D8','OC_D8')
library(DESeq2)
colData_GSE225974 <- data.frame("Sample"=colnames(data), "Donor"=c('D1','D1','D2','D2','D3','D3','D4','D4','D5','D5','D6','D6','D7','D7','D8','D8'), Group=rep(c('PBMC','OC'),8))
GSE225974_dds <-DESeqDataSetFromMatrix(countData = as.matrix(data) , colData = colData_GSE225974, design = ~Group + Donor)
GSE225974_dds <- DESeq(GSE225974_dds)

GSE225974 <- cbind(counts(GSE225974_dds,normalize=T),
                   as.data.frame (results(GSE225974_dds, contrast = c('Group','OC','PBMC')))[6],
                   as.data.frame (results(GSE225974_dds, contrast = c('Group','OC','PBMC')))[2])
names(GSE225974)[1:16] <- paste("counts",colnames(data),sep="_")
names(GSE225974)[17:18] <- c('padj_OCvsPBMC','logFC_OCvsPBMC')
GSE225974[!complete.cases(GSE225974$padj_OCvsPBMC),"padj_OCvsPBMC"] <- 1
GSE225974[!complete.cases(GSE225974$logFC_OCvsPBMC),"logFC_OCvsPBMC"] <- 0
GSE225974$Symbol <- rownames(GSE225974)

# "ReadyToUse_GSE225974.txt" is provided at OSF https://osf.io/9xys4/
write.table(GSE225974,file="ReadyToUse_GSE225974.txt", sep="\t",quote=F, col.names=T,row.names=F )
rm(GSE225974_dds,data,colData_GSE225974,GSE225974)






### GSE56815 from 10.1371/journal.pone.0116792.
# Attenuated monocyte apoptosis, a new mechanism for osteoporosis suggested by a transcriptome-wide expression study of monocytes 
# Gene expression from PBMCs of postmenopausal women with high and low BMD
# Processed data has been downloaded from GSE56815
GSE56815 <- read.delim("GSE56815.Post_lowvshigh.tsv",h=T)
names(GSE56815)[7] <- 'Symbol'
GSE56815_1 <- GSE56815[!grepl("\\///",GSE56815$Symbol),]
GSE56815_2 <- GSE56815[grepl("\\///",GSE56815$Symbol),]
nrow(GSE56815_2)
library(dplyr)

for(i in 1:1223){
  tmp <- unlist(strsplit(GSE56815_2[i,'Symbol'],"\\///"))
  tmp2 <- data.frame(matrix(NA, ncol=ncol(GSE56815_2), nrow=length(tmp)))
  colnames(tmp2) <- colnames(GSE56815_2)
  tmp2[1:length(tmp),] <- GSE56815_2[i,]
  tmp2$Symbol <- tmp
  GSE56815_2 <- rbind(GSE56815_2,tmp2)
}
GSE56815_2 <- GSE56815_2[!grepl("\\///",GSE56815_2$Symbol),]
GSE56815 <- rbind(GSE56815_1,GSE56815_2)
rm(GSE56815_1, GSE56815_2)
# "ReadyToUse_GSE56815.txt" is provided at OSF https://osf.io/9xys4/
write.table(GSE56815,file="ReadyToUse_GSE56815.txt", sep="\t",quote=F, col.names=T,row.names=F )
rm(GSE56815,i,tmp,tmp2)



### Gene expression from iliac crest biopsies https://asbmr.onlinelibrary.wiley.com/doi/full/10.1002/jbmr.396
#Molecular disease map of bone characterizing the postmenopausal osteoporosis phenotype
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)
library(fields)

# Script inspired by
#https://www.bioconductor.org/packages/devel/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html

# Raw data
# https://www.ebi.ac.uk/arrayexpress/files/E-MEXP-1618/
# ArrayExpress
# A-AFFY-44 - Affymetrix GeneChip Human Genome U133 Plus 2.0 [HG-U133_Plus_2]

# load eSet from data
load("../E-MEXP-1618.eSet.r")
# Get group membership and BMI values, was provided by mail contact to Sjur Reppe <sjur.reppe@medisin.uio.no>
Group <- read.delim("E-MEXP-1618_Group.txt")
# Get assay data and add group membership and BMI to the phenotype data
raw_data <- study
rm(study)

pData(raw_data)$Patient = Group$Patient
pData(raw_data)$Group = Group$Group
pData(raw_data)$BMI = Group$BMI
pData(raw_data)$Age = Group$Age

rm(Group)  

# Quality control using PCA plot
exprs(raw_data)[1:5, 1:5]

exp_raw <- log2(Biobase::exprs(raw_data))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Group = pData(raw_data)$Group,
                     Individual = pData(raw_data)$Patient)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Group)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  # scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("black", "orange","red"))

rm(dataGG, PCA_raw, percentVar,sd_ratio, exp_raw)

# Quality control using boxplot of oligo intensities
oligo::boxplot(raw_data, target = "core", 
               main = "Boxplot of log2-intensitites for the raw data")

# Normalization of data
# Relative Log Expression data quality analysis
# Deviation of exprssion from mean before normalization
palmieri_eset <- affy::rma(raw_data, target = "core", normalize = FALSE)
row_medians_assayData <- 
  Biobase::rowMedians(as.matrix(Biobase::exprs(palmieri_eset)))

RLE_data <- sweep(Biobase::exprs(palmieri_eset), 1, row_medians_assayData)

RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <- 
  tidyr::gather(RLE_data, patient_array, log2_expression_deviation)

ggplot2::ggplot(RLE_data_gathered, aes(patient_array,
                                       log2_expression_deviation)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-2, 2)) + 
  theme(axis.text.x = element_text(colour = "aquamarine4", 
                                   angle = 60, size = 6.5, hjust = 1 ,
                                   face = "bold"))

# Deviation of exprssion from mean after normalization
palmieri_eset <- affy::rma(raw_data, target = "core", normalize = TRUE)
row_medians_assayData <- 
  Biobase::rowMedians(as.matrix(Biobase::exprs(palmieri_eset)))

RLE_data <- sweep(Biobase::exprs(palmieri_eset), 1, row_medians_assayData)

RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <- 
  tidyr::gather(RLE_data, patient_array, log2_expression_deviation)

ggplot2::ggplot(RLE_data_gathered, aes(patient_array,
                                       log2_expression_deviation)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-2, 2)) + 
  theme(axis.text.x = element_text(colour = "aquamarine4", 
                                   angle = 60, size = 6.5, hjust = 1 ,
                                   face = "bold"))
rm(RLE_data, RLE_data_gathered, row_medians_assayData)


exp_palmieri <- Biobase::exprs(palmieri_eset)
PCA <- prcomp(t(exp_palmieri), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Group = 
                       Biobase::pData(palmieri_eset)$Group)


ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Group)) +
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  scale_color_manual(values = c("black", "grey","white"))


# Distance Heatmap
group_names <- pData(palmieri_eset)$Group
annotation_for_heatmap <- data.frame(Group = group_names)
row.names(annotation_for_heatmap) <- row.names(pData(palmieri_eset))

dists <- as.matrix(dist(t(exp_palmieri), method = "manhattan"))

rownames(dists) <- row.names(pData(palmieri_eset))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
colnames(dists) <- NULL
diag(dists) <- NA

ann_colors <- list(
  Group = c(healthy = "black", intermediate = "grey", osteoporotic="white"))

pheatmap(dists, col = (hmcol), 
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE, 
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = TRUE), 
                           max(dists, na.rm = TRUE)), 
         legend_labels = (c("small distance", "large distance")),
         main = "Clustering heatmap for the calibrated samples")

# Filtering lowly expressed genes
palmieri_medians <- rowMedians(Biobase::exprs(palmieri_eset))

hist_res <- hist(palmieri_medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
# Set manual threshold
man_threshold = 4.5
abline(v = man_threshold, col = "coral4", lwd = 2)

# The threshold will be used to filter genes that are not above the threshold in at least as many arrays as the samllest group, here osteoporotic as intermediate is only for correlation
no_of_samples <- table(pData(palmieri_eset)$Group)
rm(no_of_samples)

samples_cutoff <- 27
idx_man_threshold <- apply(Biobase::exprs(palmieri_eset), 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)
palmieri_manfiltered <- subset(palmieri_eset, idx_man_threshold)
rm(PCA,exp_palmieri, group_names, hmcol, idx_man_threshold,man_threshold,palmieri_medians, percentVar, samples_cutoff, sd_ratio, hist_res, dists, dataGG,annotation_for_heatmap, ann_colors)

# Annotation of transcripts
library(pd.hugene.1.0.st.v1)
library(hgu133plus2.db)
#head(ls("package:hugene10sttranscriptcluster.db"))

anno_palmieri <- AnnotationDbi::select(hgu133plus2.db,
                                       keys = (featureNames(palmieri_manfiltered)),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype = "PROBEID")

anno_palmieri <- subset(anno_palmieri, !is.na(SYMBOL))

# Look for multiple mappers, e.g. probids that are linked to several symbols
anno_grouped <- group_by(anno_palmieri, PROBEID)
anno_summarized <- dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))
head(anno_summarized)
anno_filtered <- filter(anno_summarized, no_of_matches > 1)
head(anno_filtered)
probe_stats <- anno_filtered 
nrow(probe_stats)

ids_to_exlude <- (featureNames(palmieri_manfiltered) %in% probe_stats$PROBEID)
table(ids_to_exlude)

palmieri_final <- subset(palmieri_manfiltered, !ids_to_exlude)
validObject(palmieri_final)

fData(palmieri_final)$PROBEID <- rownames(fData(palmieri_final))
fData(palmieri_final) <- left_join(fData(palmieri_final), anno_palmieri)
rownames(fData(palmieri_final)) <- fData(palmieri_final)$PROBEID 

validObject(palmieri_final)
rm(anno_grouped, anno_summarized, anno_filtered, anno_palmieri,probe_stats, ids_to_exlude)

# Making linear models for limma
library(limma)
# Correct for confounders (Age and BMI)
# palmieri_final_batch <-removeBatchEffect(palmieri_final,covariates=pData(palmieri_final)$Factor.Value..AGE.)
# Make the test of all contrasts against each other but estimate confounders
Age <- as.numeric(pData(palmieri_final)$Age)
BMI <- as.numeric(pData(palmieri_final)$BMI)

f <- factor(pData(palmieri_final)$Group, levels=c("healthy","intermediate","osteoporotic"))

design <- model.matrix(~0+f+Age+BMI)
colnames(design) <- c("healthy","intermediate","osteoporotic","Age","BMI")

fit <- lmFit(palmieri_final, design)
contrast.matrix <- makeContrasts(intermediate-healthy, osteoporotic-intermediate, osteoporotic-healthy,levels=coef(fit))

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef="osteoporotic - healthy", adjust="BH")

# Collect statistics, and merge with annotation and gene expression for Excel file

final.Data <-data.frame(
  cbind(
    topTable(fit2, coef="osteoporotic - healthy", adjust="BH",number=Inf )[4],
    topTable(fit2, coef="osteoporotic - healthy", adjust="BH",number=Inf )[7],
    topTable(fit2, coef="osteoporotic - healthy", adjust="BH",number=Inf )[8]
  ))
names(final.Data) <- c("LogFC_Osteoporotic_Healthy","Pval_Osteoporotic_Healthy","Padj_Osteoporotic_Healthy")
final.Data$PROBEID <- rownames(final.Data)

tmp <-data.frame(
  cbind(
    topTable(fit2, coef="intermediate - healthy", adjust="BH",number=Inf )[4],
    topTable(fit2, coef="intermediate - healthy", adjust="BH",number=Inf )[7],
    topTable(fit2, coef="intermediate - healthy", adjust="BH",number=Inf )[8]
  ))
names(tmp) <- c("LogFC_Intermediate_Healthy","Pval_Intermediate_Healthy","Padj_Intermediate_Healthy")
tmp$PROBEID <- rownames(tmp)
final.Data <- merge(final.Data, tmp, by="PROBEID")

tmp <-data.frame(
  cbind(
    topTable(fit2, coef="osteoporotic - intermediate", adjust="BH",number=Inf )[4],
    topTable(fit2, coef="osteoporotic - intermediate", adjust="BH",number=Inf )[7],
    topTable(fit2, coef="osteoporotic - intermediate", adjust="BH",number=Inf )[8]
  ))
names(tmp) <- c("LogFC_Osteoporotic_Intermediate","Pval_Osteoporotic_Intermediate","Padj_Osteoporotic_Intermediate")
tmp$PROBEID <- rownames(tmp)
final.Data <- merge(final.Data, tmp, by="PROBEID")

# include annotation
tmp <- data.frame(fData(palmieri_final))
final.Data <- merge(tmp,final.Data,  by="PROBEID")

# include gene expression values for the individual subjects
tmp <- data.frame(Biobase::exprs(palmieri_final))
# include group membership in the colnames to share the data without sharing the sensitive information on Age and BMI
colnames(tmp) <- paste(pData(palmieri_final)$Group,pData(palmieri_final)$Patient,sep="_")
tmp$PROBEID <- rownames(tmp)
final.Data <- merge(final.Data, tmp, by="PROBEID")
final.Data <- final.Data[complete.cases(final.Data$SYMBOL),]

# "ReadyToUse_E-MEXP-1618.txt" is provided at OSF https://osf.io/9xys4/
write.table(final.Data,file="ReadyToUse_E-MEXP-1618.txt",sep="\t", quote=F, col.names=T,row.names=F)
rm(design,final.Data,fit, fit2, palmieri_eset, palmieri_manfiltered, palmieri_final,raw_data,tmp,Age, BMI,f,contrast.matrix,data)



### ISMARA analysis, bam files of RNA-seq have been submitted to https://ismara.unibas.ch/mara/
# output from ISMARA analysis are provided at OSF https://osf.io/9xys4/ in the ISMARA.zip file, which includes the sampe activities and a folder with the predicted target genes for each motif
ISMARA <- read.delim("ISMARA/activity_table.txt",h=T)
ISMARA <- data.frame(t(ISMARA))

# Average
ISMARA$d0 <- apply(ISMARA[,grep("_d0",colnames(ISMARA))],1,mean)
ISMARA$d2 <- apply(ISMARA[,grep("_d2",colnames(ISMARA))],1,mean)
ISMARA$d5 <- apply(ISMARA[,grep("_d5",colnames(ISMARA))],1,mean)
ISMARA$d9 <- apply(ISMARA[,grep("_d9",colnames(ISMARA))],1,mean)

# T-Tests
# D1 versus D3, D6 og D10; D3 versus D6 og D10; D6 versus D10
ISMARA$TTd2vsd0 <- apply(ISMARA[,c(grep("_d0",colnames(ISMARA)),grep("_d2",colnames(ISMARA)))],1,function(x){return(t.test(x[1:5],x[6:13])$p.value)})
ISMARA$TTd5vsd0 <- apply(ISMARA[,c(grep("_d0",colnames(ISMARA)),grep("_d5",colnames(ISMARA)))],1,function(x){return(t.test(x[1:5],x[6:13])$p.value)})
ISMARA$TTd9vsd0 <- apply(ISMARA[,c(grep("_d0",colnames(ISMARA)),grep("_d9",colnames(ISMARA)))],1,function(x){return(t.test(x[1:5],x[6:13])$p.value)})
ISMARA$TTd5vsd2 <- apply(ISMARA[,c(grep("_d2",colnames(ISMARA)),grep("_d5",colnames(ISMARA)))],1,function(x){return(t.test(x[1:8],x[9:16])$p.value)})
ISMARA$TTd9vsd2 <- apply(ISMARA[,c(grep("_d2",colnames(ISMARA)),grep("_d9",colnames(ISMARA)))],1,function(x){return(t.test(x[1:8],x[9:16])$p.value)})
ISMARA$TTd9vsd5 <- apply(ISMARA[,c(grep("_d5",colnames(ISMARA)),grep("_d9",colnames(ISMARA)))],1,function(x){return(t.test(x[1:8],x[9:16])$p.value)})

# seperate combined motifs into individual lines
library(tidyr)
tmp <- ISMARA[!grepl("_",rownames(ISMARA)),]
tmp2 <- ISMARA[grepl("_",rownames(ISMARA)),]
# make two new columns with motif name to keep the information of the combined motifs in one column and to get then singled out in the other
tmp$Motif <- rownames(tmp)
tmp$Motif_sep <- rownames(tmp)
tmp2$Motif <- rownames(tmp2)
tmp2$Motif_sep <- rownames(tmp2)
ISMARA <- rbind(tmp,separate_rows(tmp2, Motif_sep, sep = "_", convert = FALSE))
rownames(ISMARA) <- ISMARA$Motif_sep
ISMARA$Motif_sep <- NULL

# "ReadyToUse_ISMARA.txt" is provided at OSF https://osf.io/9xys4/
write.table(ISMARA,file="ReadyToUse_ISMARA.txt",quote=F,col.names=T,row.names = T,sep="\t")
rm(tmp,tmp2)

# Get target genes
Targets <- data.frame(matrix(NA,ncol=3,nrow=1))
colnames(Targets) <- c('Motif','Target','Score')
Targets <- Targets[-1,]

# set a vector with the file names containing the targets of each individual motifs
list_of_files <- list.files(path = "ISMARA/targets_0.0", recursive = TRUE,
                            pattern = "\\.txt", 
                            full.names = TRUE)

# load the target files of each motif and rbind them in a data frame
for(k in 1:length(list_of_files)){
  tmp_Motif <- unlist(strsplit(list_of_files[k],"\\/"))[3]
  tmp_Motif <- gsub("\\..*","",as.character(tmp_Motif))
  tmp <- read.delim(list_of_files[k], h=F)
  tmp <- tmp[,c(3,2,4)]
  tmp <- tmp[grepl(tmp_Motif,tmp[,1]),]
  tmp[,4] <- apply(tmp,1,function(x){unlist(strsplit(as.character(x[3]),"\\|"))[2]})
  tmp <- tmp[,c(1,4,2)]
  colnames(tmp) <- c('Motif','Target','Score')
  tmp_Motif <- gsub("-",".",tmp_Motif)
  if(tmp_Motif %in% rownames(ISMARA)){
  } else{
    print(paste(tmp_Motif,"not in rownames(ISMARA)"))
  }
  Targets <- rbind(Targets,tmp)
}

rm(tmp, list_of_files,tmp_Motif,k)
Targets$Score <- as.numeric(Targets$Score)
Targets <- Targets[complete.cases(Targets),]

# seperate combined motifs into individual lines
library(tidyr)

tmp <- Targets[!grepl("_",Targets$Motif),]
tmp2 <- Targets[grepl("_",Targets$Motif),]
# make a second motif column to keep the information of the combined motifs in one column and to get then singled out in the other
tmp$Factor <- tmp$Motif
tmp2$Factor <- tmp2$Motif

Targets <- rbind(tmp,separate_rows(tmp2, Factor, sep = "_", convert = FALSE))

# "ReadyToUse_ISMARA_Targets.txt" is provided at OSF https://osf.io/9xys4/
write.table(Targets,file="ReadyToUse_ISMARA_Targets.txt",quote=F,col.names=T,row.names = F,sep="\t")
rm(tmp,tmp2,Targets)



### scRNA-seq data from doi: 10.1002/jbm4.10631o, raw data of 10x pipeline were provided by Yasunori Omata and Mario M. Zaiss
# Interspecies Single-Cell RNA-Seq Analysis Reveals the Novel Trajectory of Osteoclast Differentiation and Therapeutic Targets
library(Matrix)
library(Seurat)
library(valiDrops)

counts <- Read10X(data.dir = "../E-GEAD-553.processed/")
valid <- valiDrops(counts,rank_barcodes = FALSE)

## SIMPLE: Create a Seurat object with the barcodes that pass quality control
scOC <- CreateSeuratObject(counts[, colnames(counts) %in% valid[ valid$qc.pass == "pass","barcode"]], min.cells = 1, min.features = 1)
dim(scOC)

# run seurat pipeline with regression for cell cycle effects
scOC <- NormalizeData(scOC)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- s.genes[s.genes %in% rownames(scOC)]
g2m.genes <- g2m.genes[g2m.genes %in% rownames(scOC)]
scOC <- CellCycleScoring(scOC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#scOC <- ScaleData(scOC,features = rownames(scOC),vars.to.regress = c("S.Score","G2M.Score","ribosomal_fraction"))
#scOC <- FindVariableFeatures(scOC,nfeatures = 1000)
scOC <- SCTransform(scOC,do.scale=TRUE, verbose = FALSE,return.only.var.genes =FALSE,vars.to.regress = c("S.Score","G2M.Score"))
scOC <- RunPCA(scOC, verbose = FALSE)
ElbowPlot(scOC)
scOC <- RunUMAP(scOC, dims = 1:15, verbose = FALSE)
scOC <- FindNeighbors(scOC, dims = 1:2,reduction = "umap", verbose = FALSE)
scOC <- FindClusters(scOC, verbose = FALSE, resolution = 0.9)
DimPlot(scOC, label = TRUE) + NoLegend()

# "ReadyToUse_scOC.rds" is provided at OSF https://osf.io/9xys4/
saveRDS(scOC,file="ReadyToUse_scOC.rds")
rm(scOC,g2m.genes,s.genes,counts,valid)

