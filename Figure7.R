### Figure 7
# the following files are provided in OSF https://osf.io/9xys4/
Counts <- read.delim("Counts.txt",h=T)
Diff_ctr <- read.delim("Diff_ctr.txt",h=T)
colData_RNA_Diff <- data.frame("Sample"=colnames(Diff_ctr)[3:31],
                               "Donor"=gsub("_.*","",colnames(Diff_ctr)[3:31]),
                               "Timepoint"=gsub(".*_","",colnames(Diff_ctr)[3:31]))
rownames(colData_RNA_Diff)<-colData_RNA_Diff$Sample

### Figure 7A
# the following files are provided in OSF https://osf.io/9xys4/
Activity <- read.delim("OC_activity.txt")
# extract expression levels of ACP5 and NFACT1 at day 9 together with Donor information
tmp <- cbind(colData_RNA_Diff[colData_RNA_Diff$Timepoint=="d9",],
             t(Counts[Counts$Symbol == 'CTSK',colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d9','Sample']]),
             t(Counts[Counts$Symbol == 'ACP5',colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d9','Sample']]))
names(tmp)[4:5] <- c('CTSK','ACP5')

tmp <- merge(tmp,Activity, by="Donor")

# Combine the both plots in illustrator
plot(tmp$CTSK,tmp$Resorption, pch=2)
plot(tmp$ACP5,tmp$Resorption, pch=8)


### Figure 7B
# Barplot of differentially expressed genes
barplot(c(
  length(Diff_ctr[Diff_ctr$padj_Resorption_d0 < 0.01,'Symbol']),
  length(Diff_ctr[Diff_ctr$padj_Resorption_d2 < 0.01,'Symbol']),
  length(Diff_ctr[Diff_ctr$padj_Resorption_d5 < 0.01,'Symbol']),
  length(Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01,'Symbol'])))


### Figure 7C
library(goseq)
library(gplots)
# gene ontology results can vary with package updates, here we used goseq_1.42.0 with geneLenDataBase_1.26.0 

# Make a list with gene groups
Gene_groups <- list()
Gene_groups[[1]] <- Diff_ctr[Diff_ctr$padj_Resorption_d0 < 0.01 & Diff_ctr$logFC_Resorption_d0 > 0,'Symbol']
Gene_groups[[2]] <- Diff_ctr[Diff_ctr$padj_Resorption_d0 < 0.01 & Diff_ctr$logFC_Resorption_d0 < 0,'Symbol']
Gene_groups[[3]] <- Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01 & Diff_ctr$logFC_Resorption_d9 > 0,'Symbol']
Gene_groups[[4]] <- Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01 & Diff_ctr$logFC_Resorption_d9 < 0,'Symbol']
names(Gene_groups) <- c('d0_up','d0_down','d9_up','d9_down')

# Define background genes (all detected in RNA-seq that were used as inoput for DEseq analysis)
GO <- Diff_ctr$Symbol
tmp_All <- make.names(GO, unique = TRUE)

# Run goseq over all gene groups
for(i in 1:length(Gene_groups)){
  tmp <- data.frame("Symbol"=Gene_groups[[i]])
  tmp_Cl <- make.names(tmp$Symbol, unique = TRUE)
  if(length(tmp_Cl)>0){
    Cl <- as.integer(tmp_All %in% tmp_Cl)
    names(Cl) <- tmp_All
    pwf_Cl <- nullp(Cl, "hg19", "geneSymbol")
    GO.BP_Cl <- goseq(pwf_Cl, "hg19", "geneSymbol", use_genes_without_cat = TRUE)
    GO.BP_Cl$over_represented_p.adjust <- p.adjust(GO.BP_Cl$over_represented_pvalue, method = "BH")
    names(GO.BP_Cl)[8] <- names(Gene_groups)[i]
    if(i ==1){
      GO_cluster <- GO.BP_Cl[,c(1,6,7,8)]
    } else {
      GO_cluster <- merge(GO_cluster,GO.BP_Cl[,c(1,8)],by="category")
    }
  } else {}
}

# Select the long list for valuable terms
tmp <- GO_cluster[GO_cluster$ontology=="BP",]
# Select examples
BP <- unique(c("GO:0045087", "GO:0001816", "GO:0006952", "GO:0008283", "GO:0071674", "GO:0070372", "GO:0097529",
               "GO:0048870", "GO:0007155", "GO:0048771", "GO:0032940", "GO:0008360", "GO:0060326","GO:0006950",
               "GO:0007154", "GO:0009888","GO:0016477", "GO:0046903","GO:0070371","GO:0046849","GO:0045453"))

# Order the matrix for the heatmap
p <- tmp[tmp$category %in% BP,]
p2 <- cbind(p[,1:2],-log10(p[,c(4:7)]))
p2 <- p2[order(-p2[,3]),]
p2 <- p2[p2[,3]> -log10(0.1),]
for (i in 4:6){
  p3 <- cbind(p[,1:2],-log10(p[,c(4:7)]))
  p3 <- p3[order(-p3[,i]),]
  p3 <- p3[p3[,i]> -log10(0.1),]
  p2 <- rbind(p2,p3)
}
# Make heatmap for visualization
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.1),max(p2[,3:6]),length=51))

heatmap.2(as.matrix(p2[,3:6]),main="GO cluster Resorption",Rowv= F,dendrogram = 'none',  Colv=F, scale='none', col=mat_col,breaks=mat_col_breaks, trace='none', labRow=p2$term,labCol=colnames(p2)[3:6] )

rm(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,p2,p,mat_col,mat_col_breaks,hg19.geneSymbol.LENGTH,heatp,GO.BP_Cl,GO_cluster,GO,Gene_groups,pwf_Cl,y,y2,BP,diff_help,tmp_Cl,tmp_All,Cl,tmp)


### Figure 7D
library(org.Hs.eg.db)
library(clusterProfiler)
library(goseq)
library(reactome.db)

# Make a list with gene groups
Gene_groups <- list()
Gene_groups[[1]] <- Diff_ctr[Diff_ctr$padj_Resorption_d0 < 0.01 & Diff_ctr$logFC_Resorption_d0 > 0,'Symbol']
Gene_groups[[2]] <- Diff_ctr[Diff_ctr$padj_Resorption_d0 < 0.01 & Diff_ctr$logFC_Resorption_d0 < 0,'Symbol']
Gene_groups[[3]] <- Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01 & Diff_ctr$logFC_Resorption_d9 > 0,'Symbol']
Gene_groups[[4]] <- Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01 & Diff_ctr$logFC_Resorption_d9 < 0,'Symbol']
names(Gene_groups) <- c('d0_up','d0_down','d9_up','d9_down')

# Setup the reactome database and select all pathways involved in Metabolism
Reactome <- as.data.frame(reactomeEXTID2PATHID)
# the following files are provided in OSF https://osf.io/9xys4/
Relation <- read.delim("ReactomePathwaysRelation.txt", header=FALSE)
Pathways <- read.delim("ReactomePathways.txt", header=FALSE)
Pathways <- Pathways[ Pathways$V3 == "Homo sapiens",]
Pathways$DB_ID <- substr(Pathways$V1, 7, nchar(as.character(Pathways$V1)))

Metabolism <- Relation[ Relation[,1] %in% as.character(Pathways[ Pathways$V2 == "Metabolism",1]),]  
Metabolism <- rbind(Metabolism, Relation[ Relation[,1] %in% Metabolism[,2],])
Metabolism <- Metabolism[ duplicated(Metabolism[,2])==FALSE,]

Current <- nrow(Metabolism)
New <- Current + 1
while (New > Current) {
  Current <- nrow(Metabolism)
  Metabolism <- rbind(Metabolism, Relation[ Relation[,1] %in% Metabolism[,2],])
  Metabolism <- Metabolism[ duplicated(Metabolism[,2])==FALSE,]
  New <- nrow(Metabolism)
}

Pathways <- Pathways[ Pathways$V1 %in% Metabolism[,1] | Pathways$V1 %in% Metabolism[,2],]
Pathways <- Pathways[ duplicated(Pathways$V1)==FALSE,]
Reactome <- Reactome[ Reactome$DB_ID %in% Pathways$V1,]
Reactome <- split(Reactome$gene_id, f = Reactome$DB_ID, drop=T)

# Define background genes (all detected in RNA-seq that were used as inoput for DEseq analysis)
Convert <- bitr(Diff_ctr[, "Symbol"], "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)
Convert <- Convert[ duplicated(Convert$ENTREZID)==FALSE,]
Convert <- Convert[ duplicated(Convert$SYMBOL)==FALSE,]

All <- merge(Convert, Diff_ctr[,c("RefSeqID","Symbol")], by.y="Symbol", by.x="SYMBOL")
All <- unique(All)

# Prepare data frame to collect results
Enrichment <- Pathways[,c(1,2)]
colnames(Enrichment) <- c("category","name") 

for (i in 1:length(Gene_groups)){
  Cl1 <- as.integer(All[,"ENTREZID"] %in% All[ All$SYMBOL %in% Gene_groups[[i]], "ENTREZID"])
  names(Cl1) <- All[,"ENTREZID"]
  Cl1.nullp <- nullp(Cl1, genome="hg19", id="refGene")
  Cl1.EA <- goseq(Cl1.nullp, genome="hg19",id="refGene", gene2cat=Reactome, use_genes_without_cat = T)
  Cl1.EA$Cl1 <- p.adjust(Cl1.EA$over_represented_pvalue, method="fdr")
  colnames(Cl1.EA)[2] <- names(Gene_groups)[i]
  Enrichment <- merge( Enrichment, Cl1.EA[,c(1,2)], by="category", all=T)
}
Enrichment <- Enrichment[complete.cases(Enrichment),]
p <- Enrichment[apply(Enrichment[,3:(2+length(Gene_groups))],1,min)<0.05,]

#Sort the results by cluster for nice visual representation
p2 <- cbind(p[,1:2],-log10(p[,c(3:6)]))
p2 <- p2[order(-p2[,3]),]
p2 <- p2[p2[,3]> -log10(0.05),]

for (i in 4:10){
  p3 <- cbind(p[,1:2],-log10(p[,c(3:6)]))
  p3 <- p3[order(-p3[,i]),]
  p3 <- p3[p3[,i]> -log10(0.05),]
  p2 <- rbind(p2,p3)
}
# Screen for redundancy
p2$rep <- duplicated(p2$name)
p2 <- p2[p2$rep == "FALSE",]

# make a heatmap with the manual order from above
library(gplots)
library(fields)
library(scales)
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.05),max(p2[,3:6]),length=51))
heatp <- heatmap.2(as.matrix(p2[,3:6]),main="GO cluster", Rowv = F, Colv=F, dendrogram='none', scale='none', col=mat_col,breaks=mat_col_breaks, trace='none', labRow =p2$name )

rm(heatp,p2,p3,p,Enrichment,hg19.refGene.LENGTH,Metabolism,All,Cl1.EA,Cl1,Cl1.nullp,Convert, Gene_groups,Pathways,Reactome,Relation,Current,i, mat_col, mat_col_breaks,New)


### Figure 7E

#hypergeometric testing of the differntiaiton and resorption regulated genes
# Make a list with gene groups related to resorption
Gene_groups <- list()
Gene_groups[[1]] <- Diff_ctr[Diff_ctr$padj_Resorption_d0 < 0.01 & Diff_ctr$logFC_Resorption_d0 > 0,'Symbol']
Gene_groups[[2]] <- Diff_ctr[Diff_ctr$padj_Resorption_d0 < 0.01 & Diff_ctr$logFC_Resorption_d0 < 0,'Symbol']
Gene_groups[[3]] <- Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01 & Diff_ctr$logFC_Resorption_d9 > 0,'Symbol']
Gene_groups[[4]] <- Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01 & Diff_ctr$logFC_Resorption_d9 < 0,'Symbol']
names(Gene_groups) <- c('d0_up','d0_down','d9_up','d9_down')

# Make a list with gene groups related to differentiation
Gene_groups2 <- list()
for(i in 1:8){
  Gene_groups2[[i]] <- Diff_ctr[Diff_ctr$clust==i,'Symbol']
  names(Gene_groups2)[i] <- paste("clust_",i,sep="")
}

mat <- matrix(NA,nrow=length(Gene_groups2), ncol=length(Gene_groups))
rownames(mat) <- names(Gene_groups2)
colnames(mat) <- names(Gene_groups)

for (i in 1:length(Gene_groups2)){
  for (k in 1:length(Gene_groups)){
    tmp_i <- Gene_groups2[[i]]
    tmp_k <- Gene_groups[[k]]
    mat[i,k] <- phyper(length(tmp_i[tmp_i %in% tmp_k]), length(tmp_k), length(Diff_ctr[!Diff_ctr$Symbol %in% tmp_k,'Symbol']), length(tmp_i), lower.tail = F)
    
  }
}
mat <- -log10(mat)

#plot as heatmap
library(fields)
library(scales)
library(gplots)
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.05),max(mat),length=51))
heatmap.2(mat,Rowv= F,dendrogram = 'none',  Colv=F, scale='none', col=mat_col,breaks=mat_col_breaks, trace='none', labRow=rownames(mat),labCol=colnames(mat) )
rm(Gene_groups_1, Gene_groups_2,mat,i,k,mat_col,mat_col_breaks,tmp_i,tmp_k)


### Figure 7F

# Day 9 resorption heatmap
tmp <- colData_RNA_Diff[colData_RNA_Diff$Timepoint=="d9",]
# the following files are provided in OSF https://osf.io/9xys4/
Activity <- read.delim("OC_activity.txt")
tmp <- merge(tmp,Activity, by="Donor")
#remove Donor 7
tmp <- tmp[tmp$Donor != "Donor7",]
tmp[order(tmp$Resorption),'Sample']

# select genes regulated by resorption at day 9, order rows according to log fold changes and donor samples at day 9 according to resorptive activity
y <- Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01,tmp[order(tmp$Resorption),'Sample']][order(-Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01,"logFC_Resorption_d9"]),]
y <- t(scale(t(y)))
y[y < -2] <- -2
y[y > 2] <- 2

library(fields)
library(RColorBrewer)
library(gplots)
Mycol <- rev(designer.colors(n=100, col=brewer.pal(9,"Spectral")))
heatmap.2(y,trace="none",Colv = F,Rowv = F,col=Mycol,
          ColSideColors = designer.colors(100,col=brewer.pal(9,"BuPu"))[round(tmp[order(tmp$Resorption),'Resorption']*10)])

# heatmap for the log fold change, will be added to the scaled heatmap in illustrator
Mycol <- designer.colors(n=100, col=c('blue','white','red'))
Mybreaks <- seq(-1,1,length=101)
heatmap.2(cbind(
  Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01,'logFC_Resorption_d9'][order(-Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01,'logFC_Resorption_d9',])],
  Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01,'logFC_Resorption_d9'][order(-Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01,'logFC_Resorption_d9',])]),
  trace="none",Colv = F,Rowv = F,col=Mycol, breaks=Mybreaks)

# heatmap to get the colour scale for the resorption
Mycol <- designer.colors(100,col=brewer.pal(9,"BuPu"))
Mybreaks <- seq(0,10,length=101)
heatmap.2(cbind(tmp$Resorption,tmp$Resorption),trace="none",Colv = F,Rowv = F,col=Mycol, breaks=Mybreaks)
rm(Mybreaks,Mycol,y,tmp)


### Figure 7G
library(Seurat)
library(gplots)
# the following files are provided in OSF https://osf.io/9xys4/
scOC <- readRDS("ReadyToUse_scOC.rds")

# Expression of resorption associated genes per cell
Matrix <- data.frame(as.matrix(GetAssayData(scOC)))
# upregulated
tmp <- Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01 & Diff_ctr$logFC_Resorption_d9 > 0,'Symbol']
tmp <- tmp[tmp %in% rownames(scOC)]
# add expression sum to seurat object
scOC$d9_res_up <- colSums(Matrix[rownames(Matrix) %in% tmp,])
# downregulated
tmp <- Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01 & Diff_ctr$logFC_Resorption_d9 < 0,'Symbol']
tmp <- tmp[tmp %in% rownames(scOC)]
# add expression sum to seurat object
scOC$d9_res_down <- colSums(Matrix[rownames(Matrix) %in% tmp,])
# all
tmp <- Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01,'Symbol']
tmp <- tmp[tmp %in% rownames(scOC)]
# add expression sum to seurat object
scOC$d9_res <- colSums(Matrix[rownames(Matrix) %in% tmp,])

FeaturePlot(scOC,c('d9_res_up','d9_res_down','d9_res'))
rm(scOC,tmp,Matrix)


### Enrichment of resorption associated genes with gene expression dyanmics during fracture repair
Human_Mouse <- read.delim("Ensemble_SYMBOL_Mouse_Human.txt",h=T)

# Read each individual sheet of the processed data that are provided as Excel file
library("readxl")
Fracture_list <- list()
x <- 1
names <- c("full_4h","full_1d","full_3d","full_7d","full_14d","stress_4h","stress_1d","stress_3d","stress_5d","stress_7d")
for (i in c("fullfracture_DEG_DESeq2_full_4h","fullfracture_DEG_DESeq2_full_1d","fullfracture_DEG_DESeq2_full_3d","fullfracture_DEG_DESeq2_full_7d","fullfracture_DEG_DESeq2_full_14",
            "stressfracture_DESeq2_full_4hr","stressfracture_DESeq2_full_1day","stressfracture_DESeq2_full_3day","stressfracture_DESeq2_full_5day","stressfracture_DESeq2_full_7day")){
  data <- read_excel("GSE152677_DEG_DESeq2.xlsx", sheet = i)
  if(x < 6){
    #extract logFC and adjusted p-values
    data <- data[,c(1,5,9)]
  } else {
    #extract logFC and adjusted p-values
    data <- data[,c(1,3,7)]  
  }
  names(data) <- c('SYMBOL_Mouse','logFC','padj')
  # get Human Symbols
  data <- merge(data,Human_Mouse[,c("SYMBOL_Mouse","SYMBOL_Human")], by="SYMBOL_Mouse")
  data$padj <- as.numeric(data$padj)
  data[!complete.cases(data$padj),'padj'] <- 1
  # seperate into upregulated and downregulated genes for each sheet
  Fracture_list[[(2*(x-1)+1)]] <- unique(data[data$padj < 0.01 & data$logFC < 0, "SYMBOL_Human"])
  Fracture_list[[(2*(x-1)+2)]] <- unique(data[data$padj < 0.01 & data$logFC > 0, "SYMBOL_Human"])
  names(Fracture_list)[(2*(x-1)+1)] <- paste(names[x],"_down",sep="")
  names(Fracture_list)[(2*(x-1)+2)] <- paste(names[x],"_up",sep="")
  x <- x+1
}

Gene_groups <- list()
Gene_groups[[1]] <- Diff_ctr[Diff_ctr$padj_Resorption_d0 < 0.01 & Diff_ctr$logFC_Resorption_d0 > 0,'Symbol']
Gene_groups[[2]] <- Diff_ctr[Diff_ctr$padj_Resorption_d0 < 0.01 & Diff_ctr$logFC_Resorption_d0 < 0,'Symbol']
Gene_groups[[3]] <- Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01 & Diff_ctr$logFC_Resorption_d9 > 0,'Symbol']
Gene_groups[[4]] <- Diff_ctr[Diff_ctr$padj_Resorption_d9 < 0.01 & Diff_ctr$logFC_Resorption_d9 < 0,'Symbol']
names(Gene_groups) <- c('d0_up','d0_down','d9_up','d9_down')


# Calculate enrichment based on a hypergeometric test
mat <- matrix(NA,ncol=length(Fracture_list),nrow=4)
colnames(mat) <- names(Fracture_list)
rownames(mat) <- names(Gene_groups)
for (i in 1:length(Fracture_list)){
  for (k in 1:length(Gene_groups)){
    a <- Fracture_list[[i]]
    b <- Gene_groups[[k]]
    mat[k,i] <- phyper(length(a[a %in% b]),length(b),length(Diff_ctr$Symbol)-length(b),length(a),lower.tail=FALSE)
  }
}

mat[mat > 0.01] <- 1
mat <- -log10(mat)
# Order into all upregulated from full fracture, all down from full fracture, all up from stress fracture, and all down from stress fracture
mat <- t(mat[,c(1,3,5,7,9,2,4,6,8,10,11,13,15,17,19,12,14,16,18,20)])

# Plot result as heat map
library(fields)
library(gplots)
mat_col <- c('white',designer.colors(n=49, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(2,max(mat),length=50))
heatmap.2(mat,trace="none",Colv = F,Rowv = F,col=mat_col)

rm(list=ls())
