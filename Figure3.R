### Figure 3
# the following files are provided in OSF https://osf.io/9xys4/
Counts <- read.delim("Counts.txt",h=T)
Diff_ctr <- read.delim("Diff_ctr.txt",h=T)
colData_RNA_Diff <- data.frame("Sample"=colnames(Diff_ctr)[3:31],
                               "Donor"=gsub("_.*","",colnames(Diff_ctr)[3:31]),
                               "Timepoint"=gsub(".*_","",colnames(Diff_ctr)[3:31]))
rownames(colData_RNA_Diff)<-colData_RNA_Diff$Sample

### Figure 3A

# the following files are provided in OSF https://osf.io/9xys4/
ISMARA <- read.delim("ReadyToUse_ISMARA.txt",h=T)

GOI <- c('NFATC1','JUN','ALX1')
for(i in GOI){
  tmp <- colData_RNA_Diff
  tmp$x <- NA
  tmp[tmp$Timepoint=='d0','x'] <- 1
  tmp[tmp$Timepoint=='d2','x'] <- 2
  tmp[tmp$Timepoint=='d5','x'] <- 3
  tmp[tmp$Timepoint=='d9','x'] <- 4
  d0 <- unlist(ISMARA[rownames(ISMARA) == i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d0','Sample']])
  d2 <- unlist(ISMARA[rownames(ISMARA) == i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d2','Sample']])
  d5 <- unlist(ISMARA[rownames(ISMARA) == i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d5','Sample']])
  d9 <- unlist(ISMARA[rownames(ISMARA) == i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d9','Sample']])
  boxplot(d0,d2,d5,d9,ylab=paste("Motif activity",i), xaxt="none")
  axis(1,at=c(1:4),labels = c('D0','D2','D5','D9'))
  for(k in unique(colData_RNA_Diff$Donor)){
    a <- tmp[tmp$Donor==k,'x']
    b <- unlist(ISMARA[rownames(ISMARA) == i,colData_RNA_Diff[colData_RNA_Diff$Donor==k,'Sample']])
    lines(a,b)
  }
}

rm(GOI,a,b,d0,d2,d5,d9,tmp,tmp2,i,ISMARA)


### Figure 3B

# the following files are provided in OSF https://osf.io/9xys4/
ISMARA <- read.delim("ReadyToUse_ISMARA.txt",h=T)
# TFs from GO OC terms and DOI: 10.1002/jbmr.2229 - manual collection
OC_TFs <- c("JUN","CREB1","NRF1","SREBF2","NFYA","NFYB","NFYC","NR3C1","FOS","ATF1","ATF2","SMARCC2","ETS1","SIX5","ATF4","ZBTB33","BCL6","NFATC1","RFX5","JUND","MAFB","MITF","TFE3","CEBPB","JUNB","ESRRA","FOXP1")


# plot only individual motifs not factors
tmp <- ISMARA
tmp$Dup <- duplicated(tmp$Motif)
tmp <- tmp[tmp$Dup=="FALSE",]
rownames(tmp) <- tmp$Motif

# extract motif activities for significant motifs at average level and cluster them
y <- tmp[apply(tmp[,c("TTd2vsd0","TTd5vsd0","TTd9vsd0","TTd5vsd2","TTd9vsd2","TTd9vsd5")],1, min)< 0.001,c('d0','d2','d5','d9')]
y <- t(scale(t(y)))
hr <- hclust(as.dist(1-cor(t(y))))

# extract motif activities for significant motifs at donor level and use clustering from above for visualization
y <- tmp[apply(tmp[,c("TTd2vsd0","TTd5vsd0","TTd9vsd0","TTd5vsd2","TTd9vsd2","TTd9vsd5")],1, min)< 0.001,colData_RNA_Diff[order(colData_RNA_Diff$Timepoint),"Sample"]]
y <- t(scale(t(y)))
y[y< -2] <- -2
y[y > 2] <- 2

dev.off()
library(fields)
library(gplots)
library(RColorBrewer)
Mycol <- rev(designer.colors(100, col=brewer.pal(9,"Spectral")))
# As heatmap.2 does not plot all the names of the rows, the information is given in a txt file
heatp <- heatmap.2(y, trace="none", Colv=F, Rowv=as.dendrogram(hr), col=Mycol, labRow = rownames(y), )

# Reduce the information to the factors included in OC_TFs 
# Create list with each motif being split up into TFs that are combined by ISAMRA
label_row <- strsplit(rownames(y),"_")
# get the index of the list elements which contain an exact match with one of the factors from OC_TFs
a <- c()
for (i in OC_TFs){
  a <- c(a,which(lapply(label_row, function(x) grep(paste("^",i,"$", sep=""), x))!=0))
}

# replace rownames labels that are not in vector a with nothing
label <- rownames(y)
label[c(1:381)[!c(1:381) %in% a]] <- ""

# putput the labels in excel and combine it with heatmap in illustrator
write.table(rev(label[heatp$rowInd]), file="Motifs_OC_terms_heatmap.txt", sep="\n", quote=F, col.names=F, row.names=F)
rm(a,label_row,label,heatp,ISMARA,tmp,y,hr,Mycol,OC_TFs,i)


### Figure 3C
# Due that the ranking of the genes within a cluster is random, the start and end points of the lines are none-deterministic

library(circlize)
library(gplots)
library(RColorBrewer)
library(biomaRt)
library(ggplot2)
library(fields)

# the following files are provided in OSF https://osf.io/9xys4/
ISMARA <- read.delim("ReadyToUse_ISMARA.txt",h=T)
Targets <- read.delim("ReadyToUse_ISMARA_Targets.txt",h=T)

# Set up Targets by focusing on the factors that are part of the RNA_seq data 
tmp <- Targets[Targets$Score > 1 & Targets$Factor %in% Diff_ctr$Symbol & Targets$Target %in% Diff_ctr$Symbol ,c("Factor","Target")]
tmp <- unique(tmp)

# Make a vector with all genes ISAMRA network has information for
vec <- unique(c(as.character(unique(tmp[,1])), as.character(unique(tmp[,2]))))

# Reduce Diff_ctr to all the genes that are part of the vector ablove and change cluster 0 to 9, due to visualization
test <- Diff_ctr[Diff_ctr$Symbol %in% vec,c('Symbol','clust',"d0","d2","d5","d9")]
names(test)[2] <- 'Cl'
test[test$Cl == 0,'Cl'] <- 9
rm(vec)


# Order genes based on cluster information
test <- test[order(test[,2], decreasing=F),]

# Make a help data frame in which the genes are ranked for each individual cluster from 1 to the number of genes in each cluster
# This information is needed in the circos plot to identify the rank of the regulating TF (from which the lines start) to the rank of the targets (where the lines end)
test_final <- data.frame(matrix(NA, ncol=ncol(test)+1, nrow=1))
names(test_final) <- c(names(test),'Rank')
test_final <- test_final[-1,]

for (i in 1:9){
  test_1 <- test[test$Cl==i,]
  test_1$Rank <- 1:nrow(test_1)
  test_final <- rbind(test_final, test_1)
}

rm(test_1,i,test) 

# Get the rank information into the network data frame tmp 
tmp_net <- merge(tmp, test_final, by.x="Factor", by.y="Symbol")
tmp_net <- merge(tmp_net, test_final[,c(1:2,7)], by.x="Target", by.y="Symbol")
# Reorder the tmp_net
tmp_net <- tmp_net[,c(2,1,4:7,3,8:10)]
names(tmp_net)[c(7:10)] <- c('Cl_Fac','Rank_Fac','Cl_Tar','Rank_Tar')

# Order the network according to the clusters and ranl of the targets
tmp_net <- tmp_net[order(tmp_net$Cl_Fac, tmp_net$Rank_Fac),]

# Plotting using the circos package
circos.clear()
par(mfrow=c(1,2),mar = c(1, 1, 1, 1), pty="s" )
circos.par("gap.degree" = c(rep(2, 8),20))
circos.par(cell.padding = c(0.01, 0.01, 0.01, 0.01))

# Plot for NFATC1 and JUN
GOI <- c('NFATC1','JUN')

for (j in GOI){
  
  tmp_net4 <- tmp_net[tmp_net$Factor== j,]
  # Set up the individual segments of the circos plot (number of segments by number of clusters and size of segments by the rank)
  circos.initialize(factors = test_final$Cl, x = test_final$Rank )
  # Also change Diff_ctr cluster 0 to 9
  Diff_ctr_tmp <- Diff_ctr
  Diff_ctr_tmp[Diff_ctr_tmp$clust==0,'clust'] <- 9
  # Calculate the enrichment of the individual clusters for the target genes of the factor of interest
  # The numerical value will be transformed to a color with (blue depletion and red enrichement and white no change)
  # Define color scale
  mat_col <- designer.colors(n=50, col=c(rgb(0,0,255,max=255),rgb(255,255,255,max=255),rgb(255,0,0,max=255)))
  # Store color information in a vector for each cluster
  col2 <- c()
  for (i in 1:9){
    # Number of target genes in the particular cluster
    a1 <- nrow(tmp_net4[tmp_net4$Cl_Tar ==i,])
    if(a1==0){
      a <- 0
    } else {
      # Numeric enrichment value if there are target genes in the respective cluster
      a <-   log2((a1/nrow(tmp_net4)/(nrow(Diff_ctr_tmp[Diff_ctr_tmp$clust==i,])/nrow(Diff_ctr_tmp))))
      # If enrichment is above 8 or below -8 set to 8 or -8 (below is log scale) to be able to show differences in color scales
      a[a< -3] <- -3
      a[a> 3] <- 3
    }
    # Transform the value between -3 and 3 to a scale from 1 to 50 according to the 50 color scale mat_col
    b <- round((a+3)/0.125)+1
    col2 <- c(col2,b)
  }
  # Plot enrichment ring
  circos.trackPlotRegion(factors = test_final$Cl, x = test_final$Rank, ylim=c(0,1), bg.col = mat_col[col2],bg.lwd=0.1, bg.border = 'black',track.height = 0.1)
  # Plot cluster membership color
  circos.trackPlotRegion(factors = test_final$Cl ,ylim=c(0,1), bg.col = c(rainbow(8)[1:8],'grey') , bg.border = 'black',bg.lwd=0.1,track.height = 0.05)
  # Add a single black line on the cluster membership circle to indicate the rank of the TF of interest
  circos.lines(x=c(
    test_final[test_final$Symbol==j,'Rank'],
    test_final[test_final$Symbol==j,'Rank']),
    y=c(0,1),sector.index=test_final[test_final$Symbol==j,'Cl'], lwd=4)
  
  # Add the lines connecting the TF of interst with its target genes, add transparency to help indicate where the lines start
  for (i in 1:nrow(tmp_net4)){
    circos.link(tmp_net4[i,'Cl_Fac'], tmp_net4[i,'Rank_Fac'],tmp_net4[i,'Cl_Tar'],tmp_net4[i,'Rank_Tar'], directional = 0, col=alpha('black',0.6), lwd=0.5, arr.width = 0.05, arr.length = 0.1)
  }
  title(paste(j,"Target genes"))
}
rm(a,j,b,a1,GOI,i,tmp,ISMARA,Targets,col2,c, mat_col, tmp_net, tmp_net4, n, Nodes, q, Timepoint, test_final,Diff_ctr_tmp)


### Figure 3D
library(Seurat)
# the following files are provided in OSF https://osf.io/9xys4/
scOC <- readRDS("ReadyToUse_scOC.rds")
# This object has 18 clusters
DimPlot(scOC, label = TRUE, group.by = "SCT_snn_res.0.9") + NoLegend()+ theme(aspect.ratio=1)
rm(scOC)


### Figure 3E, F & G
library(Seurat)
library(gplots)
# the following files are provided in OSF https://osf.io/9xys4/
scOC <- readRDS("ReadyToUse_scOC.rds")
Targets <- read.delim("ReadyToUse_ISMARA_Targets.txt",h=T)

# Use scRNA-seq to check ISMARA network, start with extraction of averaged expression data at cluster level
Matrix <- data.frame(AverageExpression(scOC,group.by = 'SCT_snn_res.0.9', assays = 'SCT')[[1]])
Matrix_scale <- data.frame(t(scale(t(AverageExpression(scOC,group.by = 'SCT_snn_res.0.9', assays = 'SCT')[[1]]))))

# Collect average expression of all TFs of the ISAMARA network in one data frame and make a second for their targets with the same order 
sc_TFs <- data.frame(t(Matrix[rownames(Matrix) %in% unique(Targets$Factor),]))
sc_Targets_scaled <- data.frame(matrix(NA,ncol=ncol(sc_TFs), nrow=nrow(sc_TFs)))

colnames(sc_Targets_scaled) <- colnames(sc_TFs)
rownames(sc_Targets_scaled) <- rownames(sc_TFs)

for(i in colnames(sc_TFs)){
  tmp <- Targets[Targets$Factor ==i & Targets$Score > 2,'Target']
  tmp <- tmp[tmp %in% rownames(scOC)]
  sc_Targets_scaled[,i] <- colSums(Matrix[rownames(Matrix) %in% tmp,])
}

# Spearman's correlation between expression of TF and its targets across the clusters of scOC
cor_scOC_TFs_Targets <- c()
for (i in 1:ncol(sc_TFs)){
  cor_scOC_TFs_Targets <- c(cor_scOC_TFs_Targets,cor(sc_TFs[,i],sc_Targets_scaled[,i],method = "spearman"))
}

# Name the correlation vector like the TFs, order it and include a rank for gene set enrichment analysis
names(cor_scOC_TFs_Targets) <- names(sc_TFs)
cor_scOC_TFs_Targets <- cor_scOC_TFs_Targets[order(cor_scOC_TFs_Targets)]
cor_scOC_TFs_Targets <- data.frame(cbind(cor_scOC_TFs_Targets,1:length(cor_scOC_TFs_Targets)))
names(cor_scOC_TFs_Targets) <- c('Spearman','Rank')
# Remove TFs for which Spearman's correlation cannot be calculated
cor_scOC_TFs_Targets <- cor_scOC_TFs_Targets[complete.cases(cor_scOC_TFs_Targets),]


# Select example showing TF and Target gene relation at single cell level
cor_scOC_TFs_Targets[order(-cor_scOC_TFs_Targets$Spearman),][1:10,]

# Figure 3E
# Expression data per cell
Matrix <- data.frame(as.matrix(GetAssayData(scOC)))
# Select MYBL2 targets from ISMARA network (Score 2 to increase confidence of target genes due to comparing two different data sets from different labs with different protocols)
tmp <- Targets[Targets$Factor =="MYBL2" & Targets$Score > 2,'Target']
tmp <- tmp[tmp %in% rownames(scOC)]
# add expression sum of MYBL2 targets to seurat object
scOC$MYBL2_scaled <- colSums(Matrix[rownames(Matrix) %in% tmp,])

FeaturePlot(scOC,c('MYBL2','MYBL2_scaled'))

# Figure 3F
# Expression data per cluster including spearmans correlation and coloring of the clusters
dev.off()
par(mfrow=c(1,1),pty="s")
library(scales)
plot(sc_TFs$MYBL2, sc_Targets_scaled$MYBL2,col=hue_pal()(19), pch=16)
title(paste("Spearman's MYBL2:",cor_scOC_TFs_Targets["MYBL2","Spearman"]))

# Figure 3G

# Identify TFs that change expression across clusters
scOC@active.ident <- scOC$SCT_snn_res.0.9
TF_marker <- FindAllMarkers(scOC,features=colnames(sc_TFs),logfc.threshold = 0.5,min.pct = 0.25)
TF_marker <- TF_marker[TF_marker$p_val_adj < 0.05,]

# Perform gene set enrichment analysis
library(fgsea)
ranks <- cor_scOC_TFs_Targets$Spearman
names(ranks) <- rownames(cor_scOC_TFs_Targets)
tmp_list <- list()
tmp_list[[1]] <- unique(TF_marker$gene)
names(tmp_list)[1] <- 'Diff_TF'
fgseaRes <- fgsea(pathways=tmp_list, stats=ranks, nperm=1000)
p <- plotEnrichment(tmp_list$Diff_TF,
                    ranks) + labs(title=paste("diff TFs in scOC 19cluster:",length(unique(TF_marker$gene)),"padj:", fgseaRes[fgseaRes$pathway=='Diff_TF','padj']))
print(p)
# Add color scale and values for the correlations in illustrator
heatmap.2(as.matrix(cor_scOC_TFs_Targets[,c('Spearman','Spearman')]), col=designer.colors(n=100,col=c('blue','white','red')), breaks=seq(-1,1,length=101), Colv = F, Rowv = F, dendrogram = "none", labRow = F, labCol = F, trace="none")
unique(TF_marker$gene)

rm(Matrix, Matrix_scale,Targets, tmp, tmp_list, scOC,sc_TFs, sc_Targets_scaled,p,fgseaRes,ranks,i,TF_marker,cor_scOC_TFs_Targets)

### Figure 3H
library(tidyr)
library(neat)
library(qgraph)
# the following files are provided in OSF https://osf.io/9xys4/
Targets <- read.delim("ReadyToUse_ISMARA_Targets.txt",h=T)

# Make two list containing the clustered genes
Aset <- list()
Bset <- list()
for (i in c(1:8)){
  Bset[[paste("Cl_",i,sep="")]] <- Diff_ctr[Diff_ctr$clust ==i,"Symbol"]
  Aset[[paste("Cl_",i,sep="")]] <- Diff_ctr[Diff_ctr$clust ==i,"Symbol"]
}

# Reduce ISMARA network to clustered genes
tmp <- Targets[Targets$Factor %in% Diff_ctr[Diff_ctr$clust !=0,"Symbol"],]
tmp <- tmp[tmp$Target %in% Diff_ctr[Diff_ctr$clust !=0,"Symbol"],]

# Define all nodes in the network
Nodes <- rbind(data.frame("Factor" = unique(tmp$Factor)),data.frame("Factor" = unique(tmp$Target)))
Nodes <- unique(Nodes)
Nodes <- as.vector(Nodes[,1])
# Run NEAT analysis using the ISMARA network (Score above 1 to include more confident genes)
test <- neat(alist = Aset, blist =Bset, network = as.matrix(unique(tmp[tmp$Score>0,c('Factor','Target')])), nettype = 'directed', nodes = Nodes)
test$LogOdds <- test$nab / test$expected_nab
test <- data.frame(test)

# Reduce the network enrichment analysis to significant interactions between the clusters
tmp <- test[test$adjusted_p < 0.01 & test$LogOdds > 1 ,c(1:2,7)]
tmp$LogOdds <- log2(tmp$LogOdds)


# Show spring layout for all clusters - nondeterministic and can result in different visualization
layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
qg <- qgraph(as.matrix(tmp), pty="s", layout="spring",directed=T, color='white', width=7, height=7, posCol='black', negCol='red', vsize=9.7 )
title(paste("ISMARA"), cex.main=0.7, line=3)
col <- data.frame("col"=qg$graphAttributes$Edges$color,
                  "value" <- tmp[tmp[,3]>0,3])
col <- col[order(col[,2],decreasing=F),]
col <- col[col[,2]>0,]
legend_image <- as.raster(matrix(rev(col[,1]), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=1.5, y = c(0,0.5,1), labels = c(round(col[1,2],4),round((col[nrow(col),2]+col[1,2])/2,4),round(col[nrow(col),2],4)))
rasterImage(legend_image, 0, 1, 0.5,0)
# Adjust some arrows in illustrator for better visualization and include cluster labels

rm(legend_image, qg, test, tmp,i, legend_image,Nodes,value,Targets,col, Aset, Bset)


### Figure 3I

# the following files are provided in OSF https://osf.io/9xys4/
# Get GWAS causal genes from  doi:10.1038/s41588-018-0302-x part of supplementary table 2
Targets <- read.delim("ReadyToUse_ISMARA_Targets.txt",h=T)
GWAS <-  read.csv("ebmd_con_indepen.csv",h=T)
names(GWAS)[2] <- "Lead_SNP"
GWAS$Start <- GWAS$BP

# Are the factors regulating causal SNP-eBMD-genes or genes from cluster 4
# Reduce Factors in ISMARA network to TFs that are expressed in the RNA-seq data set Diff_ctr
Targets_tmp <- Targets[Targets$Factor %in% Diff_ctr$Symbol,]

# Set up a matrix
mat <- matrix(NA, ncol=2, nrow=length(unique(Targets_tmp$Factor)))
colnames(mat) <- c("Enrich_GWAS","Enrich_CL4")
rownames(mat) <- unique(Targets_tmp$Factor)
# Collect all genes that are causual eBMD target genes and also predicted targets in the ISAMRA network
tmp <- GWAS[GWAS$P.NI < 5e-8,'C.GENE']
tmp <- tmp[tmp %in% Targets_tmp$Target]
# Collect all genes that are in cluster 4 and also predicted targets in the ISAMRA network
tmp2 <- Diff_ctr[Diff_ctr$Symbol %in% Targets_tmp$Target & Diff_ctr$clust==4,"Symbol"]
# Collect p-values using hypergeometric test
for (i in unique(Targets_tmp$Factor)){
  k <- unique(Targets_tmp[Targets_tmp$Factor == i & Targets_tmp$Score >1,'Target'])
  q <- unique(k[k %in% tmp])
  m <- unique(tmp)
  n <- unique(Targets_tmp$Target[!Targets_tmp$Target %in% tmp])
  mat[i,1] <- phyper(length(q), length(m), length(n), length(k), lower.tail = FALSE)
  
  q <- unique(k[k %in% tmp2])
  m <- tmp2
  n <- unique(Targets_tmp$Target[!Targets_tmp$Target %in% tmp2])
  mat[i,2] <- phyper(length(q), length(m), length(n), length(k), lower.tail = FALSE)
 
}
# Combine p-values for the enrichments with nominal p-value from eBMD GWAS
TF_Enrichment <- data.frame(mat)
TF_Enrichment$Factor <- rownames(TF_Enrichment)
TF_Enrichment <- merge(TF_Enrichment,GWAS[,c("C.GENE","P.NI")],by.x="Factor",by.y="C.GENE", all.x=T)
TF_Enrichment[!complete.cases(TF_Enrichment$P.NI),'P.NI'] <- 1

# barplot
# subgroup1: differentially expressed and a causal eBMD gene
y1 <- TF_Enrichment[TF_Enrichment$Factor %in% Diff_ctr[Diff_ctr$clust!=0,'Symbol'] & TF_Enrichment$P.NI < 5e-8, ]
# subgroup2: differentially expressed and NOT a causal eBMD gene
y2 <- TF_Enrichment[TF_Enrichment$Factor %in% Diff_ctr[Diff_ctr$clust!=0,'Symbol'] & TF_Enrichment$P.NI > 5e-8, ]
# subgroup3: NOT differentially expressed and a causal eBMD gene
y3 <- TF_Enrichment[TF_Enrichment$Factor %in% Diff_ctr[Diff_ctr$clust==0,'Symbol'] & TF_Enrichment$P.NI < 5e-8, ]
# subgroup4: NOT differentially expressed and NOT a causal eBMD gene
y4 <- TF_Enrichment[TF_Enrichment$Factor %in% Diff_ctr[Diff_ctr$clust==0,'Symbol'] & TF_Enrichment$P.NI > 5e-8, ]

# Barplot left panel, which is modified later in Illustrator, using the second bar just for upper alignment of the fourth, thrid and fourth bar are then aligned with the first one and the second is deleted
barplot(
  c(
    # All TFs from the group
    nrow(y1),
    # This number is useful for visualization, align the next at lower end and second next at upper end
    nrow(y1[y1$Enrich_CL4 < 0.01 | y1$Enrich_GWAS < 0.01,]),
    # All TFs enriched for cluster 4 genes
    nrow(y1[y1$Enrich_CL4 < 0.01,]),
    # All TFs enriched for eBMD genes
    nrow(y1[y1$Enrich_GWAS < 0.01,]),
    nrow(y2),
    nrow(y2[y2$Enrich_CL4 < 0.01 | y2$Enrich_GWAS < 0.01,]),
    nrow(y2[y2$Enrich_CL4 < 0.01,]),
    nrow(y2[y2$Enrich_GWAS < 0.01,]),
    nrow(y3),
    nrow(y3[y3$Enrich_CL4 < 0.01 | y3$Enrich_GWAS < 0.01,]),
    nrow(y3[y3$Enrich_CL4 < 0.01,]),
    nrow(y3[y3$Enrich_GWAS < 0.01,]),
    nrow(y4),
    nrow(y4[y4$Enrich_CL4 < 0.01 | y4$Enrich_GWAS < 0.01,]),
    nrow(y4[y4$Enrich_CL4 < 0.01,]),
    nrow(y4[y4$Enrich_GWAS < 0.01,])
  )
)
rm(y1,y2,y3,y4)

# Scatter plot right panel
# Select TFs that are enriched for target genes in cluster4 or causal eBMD genes
TF_Enrichment$Enrich_GWAS_log10 <- -log10(TF_Enrichment$Enrich_GWAS)
TF_Enrichment$Enrich_CL4_log10 <- -log10(TF_Enrichment$Enrich_CL4)
# Set very low p-values for cluster 4 to 100 and for eBMD to 15 (according to the distribution for each value)
TF_Enrichment[is.infinite(TF_Enrichment$Enrich_CL4_log10),'Enrich_CL4_log10'] <- 100
TF_Enrichment[is.infinite(TF_Enrichment$Enrich_GWAS_log10),'Enrich_GWAS_log10'] <- 15
TF_Enrichment[TF_Enrichment$Enrich_CL4_log10 >100,'Enrich_CL4_log10'] <- 100
TF_Enrichment[TF_Enrichment$Enrich_GWAS_log10 >15,'Enrich_GWAS_log10'] <- 15

# Select for significant enrichments
tmp <- TF_Enrichment[apply(TF_Enrichment[,c('Enrich_CL4_log10','Enrich_GWAS_log10')],1,max)>2,]
# Replace non-significant ones with 0
tmp[tmp$Enrich_CL4_log10 < 2,'Enrich_CL4_log10'] <- 0
tmp[tmp$Enrich_GWAS_log10 < 2,'Enrich_GWAS_log10'] <- 0
# Introduce color code for the eBMD GWAS if the TF is listed as causal gene
tmp[tmp$P.NI > 5e-8,'P.NI'] <- 1
tmp$P.NI <- -log10(tmp$P.NI)
# Add color value to data.frame
# first make breaks similar to a heatmap function
mat_col_breaks <- c(0,seq(-log10(5e-8),max(tmp$P.NI),length=50))
# Get the closest lower index of the break vector mat_col_breaks
tmp$col <- NA
for (i in 1:nrow(tmp)){
  tmp[i,'col'] <- which(mat_col_breaks==max(mat_col_breaks[mat_col_breaks <= tmp[i,'P.NI']]))
}
# Replace index with color from color vector
mat_col <- c(alpha('white',0),designer.colors(n=50, col=c('plum1','darkmagenta')))
tmp$col <- mat_col[tmp$col]


# Scatter plot
plot(tmp$Enrich_GWAS_log10,tmp$Enrich_CL4_log10,pch=21, bg=tmp$col)
# Show names of the factors that are themselves eBMD genes or which are enriched for both gene groups
tmp2 <- tmp[tmp$P.NI > -log10(5e-8) | (tmp$Enrich_CL4 < 0.01 & tmp$Enrich_GWAS < 0.01),]
text(tmp2$Enrich_GWAS_log10,tmp2$Enrich_CL4_log10,tmp2$Factor)

# Make a heatmap to get the color scale
library(fields)
library(gplots)
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(5e-8),max(tmp$P.NI),length=51))
heatp <- heatmap.2(as.matrix(tmp[,c('P.NI','P.NI')]), col=mat_col, breaks=mat_col_breaks,trace="none",Colv =F, Rowv = F)
rm(Targets_tmp,Targets,tmp, tmp2,a,b,j,k,mat,heatp,mat_col,mat_col_breaks,i,m,n,q,TF_Enrichment,GWAS)

