### Figure 4
# the following files are provided in OSF https://osf.io/9xys4/
Counts <- read.delim("Counts.txt",h=T)
Diff_ctr <- read.delim("Diff_ctr.txt",h=T)
colData_RNA_Diff <- data.frame("Sample"=colnames(Diff_ctr)[3:31],
                               "Donor"=gsub("_.*","",colnames(Diff_ctr)[3:31]),
                               "Timepoint"=gsub(".*_","",colnames(Diff_ctr)[3:31]))
rownames(colData_RNA_Diff)<-colData_RNA_Diff$Sample

### Figure 4A, 4B, 4C and 4D
library(Seurat)
# the following files are provided in OSF https://osf.io/9xys4/
Human_Mouse <- read.delim("Ensemble_SYMBOL_Mouse_Human.txt",h=T)
scOC <- readRDS("ReadyToUse_scOC.rds")
# Osteomorph gene list from doi: 10.1016/j.cell.2021.02.002
# Osteoclasts recycle via osteomorphs during RANKL-stimulated bone resorption
Osteomorphs <- read.delim("Osteomorphs.txt",h=T)
# In doi: 10.1016/j.cell.2021.02.002 they report 151 genes that are only uprgeulated in osteomorphs but not osteoclasts when comparing to macrophages

boxplot(Osteomorphs[Osteomorphs$Common.in.osteoclasts.and.osteomorph!="Common",2],
        Osteomorphs[Osteomorphs$Common.in.osteoclasts.and.osteomorph!="Common",3],
        Osteomorphs[Osteomorphs$Common.in.osteoclasts.and.osteomorph!="Common",4],
        # Delta between Osteomorphs an Osteoclasts
        Osteomorphs[Osteomorphs$Common.in.osteoclasts.and.osteomorph!="Common",2]-Osteomorphs[Osteomorphs$Common.in.osteoclasts.and.osteomorph!="Common",3],
        Osteomorphs[Osteomorphs$Common.in.osteoclasts.and.osteomorph=="Common",2],
        Osteomorphs[Osteomorphs$Common.in.osteoclasts.and.osteomorph=="Common",3],
        Osteomorphs[Osteomorphs$Common.in.osteoclasts.and.osteomorph=="Common",4],
        # Delta between Osteomorphs an Osteoclasts
        Osteomorphs[Osteomorphs$Common.in.osteoclasts.and.osteomorph=="Common",2]-Osteomorphs[Osteomorphs$Common.in.osteoclasts.and.osteomorph=="Common",3],
        names=c('OM','OC','Mac','OM-OC'))
# Selection must be based on p-values only, as some genes are nominal higher expressed in osteoclasts, but there they might not have reached statistical significant higher expression as compared to macrophages 

# Get Human Symbols for osteomorph genes leads to 131 genes
Osteomorphs <- merge(Osteomorphs,Human_Mouse[,c("SYMBOL_Mouse","SYMBOL_Human")], by.x="Symbol",by.y="SYMBOL_Mouse") 

# Get all marker genes of the single cell cluster
AllMarkers <- FindAllMarkers(scOC)

# Barplot for osteomorph genes and their overlap with dynamic RNA-seq and scRNA-seq Marker genes
Gene_groups <- list()
Gene_groups$RNA_seq <- Diff_ctr[Diff_ctr$clust != 0,'Symbol']
Gene_groups$scMarker <- unique(AllMarkers[AllMarkers$p_val_adj < 0.01 & AllMarkers$avg_log2FC >0,'gene'])
Gene_groups$Osteomorph_only <- Osteomorphs[Osteomorphs$Common.in.osteoclasts.and.osteomorph != "Common","SYMBOL_Human"]
Gene_groups$Osteomorph_Osteoclast <- Osteomorphs[Osteomorphs$Common.in.osteoclasts.and.osteomorph == "Common","SYMBOL_Human"]

## Figure 4A, the second and sixth bar is for upper alignment of the 4th and 8th in illustrator.
barplot(c(
  length(Gene_groups[[3]]),
  length(Gene_groups[[3]][Gene_groups[[3]] %in% Gene_groups[[1]] | Gene_groups[[3]] %in% Gene_groups[[2]]]),
  length(Gene_groups[[3]][Gene_groups[[3]] %in% Gene_groups[[1]]]),
  length(Gene_groups[[3]][Gene_groups[[3]] %in% Gene_groups[[2]]]),
  length(Gene_groups[[4]]),
  length(Gene_groups[[4]][Gene_groups[[4]] %in% Gene_groups[[1]] | Gene_groups[[4]] %in% Gene_groups[[2]]]),
  length(Gene_groups[[4]][Gene_groups[[4]] %in% Gene_groups[[1]]]),
  length(Gene_groups[[4]][Gene_groups[[4]] %in% Gene_groups[[2]]])))


# Enrichment of osteomorph genes for clustered genes using hypergeometric test and calculate numeric enrichment
mat <- matrix(NA,nrow=8, ncol=4)
# Left with 118 osteomorph genes that are in the RNA-seq data set
for (i in 1:8){
  tmp1 <- Diff_ctr[Diff_ctr$clust ==i,'Symbol']
  tmp_osteomorph <- Osteomorphs[Osteomorphs$Common.in.osteoclasts.and.osteomorph!="Common" & Osteomorphs$SYMBOL_Human %in% Diff_ctr$Symbol,]
  tmp2 <- tmp_osteomorph[tmp_osteomorph$SYMBOL_Human %in% tmp1,'SYMBOL_Human']
  mat[i,1] <- -log10(phyper(length(tmp2),length(tmp1),length(Diff_ctr$Symbol)-length(tmp1),length(tmp_osteomorph$SYMBOL_Human),lower.tail = F))
  mat[i,2] <- log2((length(tmp2)/length(tmp_osteomorph$SYMBOL_Human))/(length(tmp1)/length(Diff_ctr$Symbol)))
  tmp_osteomorph <- Osteomorphs[Osteomorphs$Common.in.osteoclasts.and.osteomorph=="Common",]
  tmp2 <- tmp_osteomorph[tmp_osteomorph$SYMBOL_Human %in% tmp1,'SYMBOL_Human']
  mat[i,3] <- -log10(phyper(length(tmp2),length(tmp1),length(Diff_ctr$Symbol)-length(tmp1),length(tmp_osteomorph$SYMBOL_Human),lower.tail = F))
  mat[i,4] <- log2((length(tmp2)/length(tmp_osteomorph$SYMBOL_Human))/(length(tmp1)/length(Diff_ctr$Symbol)))
  
}
## Figure 4B Visualize enrichment in a heatmap
library(fields)
library(gplots)
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.01),max(mat[,c(1,3)]),length=51))
heatmap.2(as.matrix(mat[,c(1,3)]),main="GO cluster", Rowv = F, Colv=F, dendrogram='none', scale='none', col=mat_col,breaks=mat_col_breaks, trace='none' )


# Enrichment of osteomorph genes for marker genes of the scRNA-seq clusters using hypergeometric test and calculate numeric enrichment
# Collect the scRNA-seq marker genes in a list
Gene_groups <- list()
x <- 1
for (i in levels(scOC)){
  Gene_groups[[x]] <- AllMarkers[AllMarkers$cluster==i & AllMarkers$p_val_adj < 0.01 & AllMarkers$avg_log2FC >0,'gene']
  names(Gene_groups)[x] <- paste('Cl_',i,sep="")
  x <- x+1
}
# collect the p-value for the enrichment in a matrix
mat <- matrix(NA,nrow=length(Gene_groups), ncol=4)
rownames(mat) <- names(Gene_groups)
# Collect the overlapping genes in a list
mat_list_NOTcommon <- list()
mat_list_common <- list()
for (i in 1:length(Gene_groups)){
  tmp1 <- Gene_groups[[i]]
  tmp_osteomorph <- Osteomorphs[Osteomorphs$Common.in.osteoclasts.and.osteomorph!="Common" & Osteomorphs$SYMBOL_Human %in% rownames(scOC),]
  tmp2 <- tmp_osteomorph[tmp_osteomorph$SYMBOL_Human %in% tmp1,'SYMBOL_Human']
  if(length(tmp2)>0){
    mat[i,1] <- -log10(phyper(length(tmp2),length(tmp1),length(rownames(scOC))-length(tmp1),length(tmp_osteomorph$SYMBOL_Human),lower.tail = F))
    mat[i,2] <- log2((length(tmp2)/length(tmp_osteomorph$SYMBOL_Human))/(length(tmp1)/length(rownames(scOC))))
  } else{
    mat[i,1] <- 0
    mat[i,2] <- 0
  }
  if(mat[i,1] > 2){ mat_list_NOTcommon[[i]] <- tmp2} else{mat_list_NOTcommon[[i]] <- ""}
  
  tmp_osteomorph <- Osteomorphs[Osteomorphs$Common.in.osteoclasts.and.osteomorph=="Common" & Osteomorphs$SYMBOL_Human %in% rownames(scOC),]
  tmp2 <- tmp_osteomorph[tmp_osteomorph$SYMBOL_Human %in% tmp1,'SYMBOL_Human']
  if(length(tmp2)>0){
    mat[i,3] <- -log10(phyper(length(tmp2),length(tmp1),length(rownames(scOC))-length(tmp1),length(tmp_osteomorph$SYMBOL_Human),lower.tail = F))
    mat[i,4] <- log2((length(tmp2)/length(tmp_osteomorph$SYMBOL_Human))/(length(tmp1)/length(rownames(scOC))))
  } else{
    mat[i,3] <- 0
    mat[i,4] <- 0
  }
  if(mat[i,3] > 2){ mat_list_common[[i]] <- tmp2} else{mat_list_common[[i]] <- ""}
  
}
names(mat_list_NOTcommon) <- names(Gene_groups)
names(mat_list_common) <- names(Gene_groups)
lapply(mat_list,length)

## Figure 4C show enrichment in heatmap
library(fields)
library(gplots)
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.01),max(mat[,c(1,3)]),length=51))
heatmap.2(as.matrix(mat[,c(1,3)]),main="GO cluster", Rowv = F, Colv=F, dendrogram='none', scale='none', col=mat_col,breaks=mat_col_breaks, trace='none' )

## Figure 4D
# Show dotplot for the Osteomorph only genes that are among the enriched cluster (only cluster 6)
genes <- unique(unlist(mat_list_NOTcommon))
genes <- genes[!genes %in% ""]
DotPlot(scOC,features=unique(genes))
# Get cluster information of the genes (Add this in illustrator)
Diff_ctr[Diff_ctr$Symbol %in% genes,c('Symbol','clust')]
rm(mat_list_common,mat_list_NOTcommon,tmp_osteomorph,Osteomorphs,scOC,genes, i,mat,Human_Mouse,Gene_groups,AllMarkers,mat_col, mat_col_breaks,k,tmp1,tmp2,x)

### Figure 4E
library(Seurat)
# the following files are provided in OSF https://osf.io/9xys4/
scOC <- readRDS("ReadyToUse_scOC.rds")
# Show cell cycle scoring for the clusters
CCS <- table(scOC$SCT_snn_res.0.9,scOC$Phase)
for(i in 1:nrow(CCS)){
  CCS[i,] <- CCS[i,]/unlist(sum(CCS[i,]))
}
barplot(t(CCS),legend=T)
rm(scOC,CCS,i)

### Figure 4F
library(Seurat)
# the following files are provided in OSF https://osf.io/9xys4/
scOC <- readRDS("ReadyToUse_scOC.rds")

FeaturePlot(scOC,c('CD14'))+ theme(aspect.ratio=1)
FeaturePlot(scOC,c('CTSK'))+ theme(aspect.ratio=1)
rm(scOC)

### Figure 4G, 4H, 4I, 4J, and 4H
library(Seurat)
# the following files are provided in OSF https://osf.io/9xys4/
scOC <- readRDS("ReadyToUse_scOC.rds")
# This object has 18 clusters, reduce to 3 (left, middel and right)
DimPlot(scOC, label = TRUE, group.by = "SCT_snn_res.0.9") + NoLegend()+ theme(aspect.ratio=1)

# set identities into left and right for new DimPlot
new.identities <- c('right','right','right','right','right','right','right','left','right','right','right','right','middle','leftdiff','left','left','right','rightdiff','rightdiff')
names(new.identities) <- levels(scOC)
scOC <- RenameIdents(scOC, new.identities)

## Figure 4G
DimPlot(scOC, label = TRUE) + NoLegend()+ theme(aspect.ratio=1)

# Find signature genes for left and right and for the differentiated ones
tmp <- FindMarkers(scOC,ident.1 = c("left","leftdiff"), ident.2 = c("right","rightdiff"), min.pct = 0, logfc.threshold = 0.01)
tmp$Symbol <- rownames(tmp)
# Include strong changes here which are not only due to a single cluster on each side
tmp <- tmp[abs(tmp$avg_log2FC)>0.9 | tmp$p_val_adj < 1e-30,]

tmp_left <- FindMarkers(scOC,ident.1 = "leftdiff",ident.2 = "middle", min.pct = 0, logfc.threshold = 0.01)
tmp_left$Symbol <- rownames(tmp_left)
tmp_left <- tmp_left[tmp_left$avg_log2FC>0 & tmp_left$p_val_adj < 1e-4,]

tmp_right <- FindMarkers(scOC,ident.1 = "rightdiff",ident.2 = "middle", min.pct = 0., logfc.threshold = 0.01)
tmp_right$Symbol <- rownames(tmp_right)
tmp_right <- tmp_right[tmp_right$avg_log2FC>0 & tmp_right$p_val_adj < 1e-4,]

table(tmp_left$Symbol %in% tmp_right$Symbol)

# Make Gene groups for left, right and differentiated (both arms), used for 4H and the following enrichment analyses
Gene_groups <- list()
Gene_groups[[1]] <- rownames(tmp)[tmp$avg_log2FC >0]
Gene_groups[[2]] <- rownames(tmp)[tmp$avg_log2FC <0]
Gene_groups[[3]] <- rownames(tmp_left)[tmp_left$Symbol %in% tmp_right$Symbol]
Gene_groups[[1]] <- Gene_groups[[1]][!Gene_groups[[1]] %in% Gene_groups[[3]]]
Gene_groups[[2]] <- Gene_groups[[2]][!Gene_groups[[2]] %in% Gene_groups[[3]]]
names(Gene_groups) <- c('left','right','both')
lapply(Gene_groups,length)

## Figure 4H Do Heatmap for signature genes, individually and combine them in illustrator
# Order the cluster for the heatmap
levels(scOC) <- c('leftdiff','left','middle','right','rightdiff')
DoHeatmap(scOC,Gene_groups[[1]],lines.width = 50)
DoHeatmap(scOC,Gene_groups[[2]],lines.width = 50)
DoHeatmap(scOC,Gene_groups[[3]],lines.width = 50)

# Keep Gene_groups for the enrichment and network analysis
rm(genes, tmp, tmp_left, tmp_right,tmp_osteomorph,new.identities,scOC)

# Gene ontology analysis of the signature genes for left, right and differentiated cells
library(goseq)
library(gplots)
# the following files are provided in OSF https://osf.io/9xys4/
scOC <- readRDS("ReadyToUse_scOC.rds")

# Define background genes (all detected in RNA-seq that were used as inoput for DEseq analysis)
GO <- rownames(scOC)
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
# Reduce to biological processes and screen the list for useful terms that show insight into biology
GO_cluster <- GO_cluster[GO_cluster$ontology=='BP',]
# Select terms
GO_cluster <- GO_cluster[GO_cluster$term %in% c('tRNA import into mitochondrion','bone resorption','pH reduction','syncytium formation by plasma membrane fusion','nuclear-transcribed mRNA catabolic process','translation','intracellular protein transport','NADH regeneration','canonical glycolysis','phagosome maturation','transferrin transport','ossification','ron ion transport','mitochondrial unfolded protein response','import into nucleus','negative regulation of ubiquitin protein ligase activity','ATP transport','immune system process','signal transduction','vesicle-mediated transport','myeloid cell activation involved in immune response','cell migration','chemotaxis','cytokine production','cell adhesion','lipid localization','endocytosis','extracellular matrix organization','cell-cell fusion','foam cell differentiation','fatty acid transport','lipid storage'),]

# Collect the p-values in -log10 form and order the matrix according to enrichemnt in the three groups
tmp2 <- cbind(GO_cluster[,1:2], -log10(GO_cluster[,4:6]))
# Some processes have a p-value of 0 giving infinite value after log transformation,change to 25 (higher as all other values) 
tmp2[tmp2=="Inf"] <- 25
p2 <- tmp2
p2 <- p2[order(-p2[,3]),]
p2 <- p2[p2[,3]> -log10(0.01),]

for (i in 4:5){
  p3 <- tmp2
  p3 <- p3[order(-p3[,i]),]
  p3 <- p3[p3[,i]> -log10(0.01),]
  p2 <- rbind(p2,p3)
}
## Figure 4I Show result as heatmap
library(fields)
library(gplots)
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.01),max(p2[,3:5]),length=51))
heatmap.2(as.matrix(p2[,3:5]),Colv=F,Rowv = F, scale='none', col=mat_col,breaks=mat_col_breaks, trace='none', labRow=p2$term,labCol=colnames(p2)[3:5] )
rm(mat_col, p3,p2,GO,Cl,tmp, tmp2, mat_col_breaks,i,tmp_All,tmp_Cl,pwf_Cl,hg19.geneSymbol.LENGTH,GO.BP_Cl,GO_cluster)

# Make second gene list with the clustered genes from the RNA-seq data
Gene_groups_2 <- list()
for (i in 1:8){
  Gene_groups_2[[i]] <- Diff_ctr[Diff_ctr$clust==i,'Symbol']
}
names(Gene_groups_2) <- paste('Cl',1:8,sep="_")
# Reduce the scRNA-seq signature genes to those that are also found in the RNA-seq data
Gene_groups[[1]] <- Gene_groups[[1]][Gene_groups[[1]] %in% Diff_ctr$Symbol]
Gene_groups[[2]] <- Gene_groups[[2]][Gene_groups[[2]] %in% Diff_ctr$Symbol]
Gene_groups[[3]] <- Gene_groups[[3]][Gene_groups[[3]] %in% Diff_ctr$Symbol]

# Enrichment of scRNA-seq signature (left,right, both) genes for clustered genes using hypergeometric test and calculate numeric enrichment
mat <- matrix(NA,nrow=length(Gene_groups), ncol=length(Gene_groups_2))
rownames(mat) <- names(Gene_groups)
colnames(mat) <- names(Gene_groups_2)
for (i in 1:length(Gene_groups)){
  tmp1 <- Gene_groups[[i]]
  for (k in 1:length(Gene_groups_2)){
    tmp2 <- Gene_groups_2[[k]]
    mat[i,k] <- -log10(phyper(length(tmp2[tmp2 %in% tmp1]),length(tmp1),length(Diff_ctr$Symbol)-length(tmp1),length(tmp2),lower.tail = F))
  }
}
## Figure 4J Show result as heatmap¨
library(fields)
library(gplots)
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.01),max(mat),length=51))
heatmap.2(as.matrix(mat),main="GO cluster", Rowv = F, Colv=F, dendrogram='none', scale='none', col=mat_col,breaks=mat_col_breaks, trace='none' )
rm(i,k,mat_col_breaks,mat_col,tmp1,tmp2,mat,Gene_groups_2)


# Check if the ISAMRA network indicates left and right specific interactions of a core network that is important for the regulation of genes in the differentiated cells

# the following files are provided in OSF https://osf.io/9xys4/
Targets <- read.delim("ReadyToUse_ISMARA_Targets.txt",h=T)
# TFs from GO OC terms and DOI: 10.1002/jbmr.2229 - manual collection
OC_TFs <- c("JUN","CREB1","NRF1","SREBF2","NFYA","NFYB","NFYC","NR3C1","FOS","ATF1","ATF2","SMARCC2","ETS1","SIX5","ATF4","ZBTB33","BCL6","NFATC1","RFX5","JUND","MAFB","MITF","TFE3","CEBPB","JUNB","ESRRA","FOXP1")

# Reduce Regulators and Targets to the genes that are detected in the scRNA-seq data
tmp <- Targets[Targets$Factor %in% rownames(scOC) & Targets$Target %in% rownames(scOC),]

# Make enrichment using hypergeometric test
mat <- matrix(NA, ncol=length(Gene_groups), nrow=length(unique(tmp$Factor)))
colnames(mat) <- paste('Pval_',names(Gene_groups),sep="")
rownames(mat) <- unique(tmp$Factor)
for (i in unique(tmp$Factor)){
  for (l in 1:length(Gene_groups)){
    k <- unique(tmp[tmp$Factor == i & tmp$Score >1,'Target'])
    q <- unique(k[k %in% Gene_groups[[l]]])
    m <- Gene_groups[[l]]
    n <- rownames(scOC)[!rownames(scOC) %in% Gene_groups[[l]]]
    if(length(q)>0){
      mat[i,l] <- phyper(length(q), length(m), length(n), length(k), lower.tail = FALSE)
    } else {
      mat[i,l] <- 1
    }
  }
  
}
TF_Enrichment <- data.frame(mat)

# Do adjustment for multiple testing
TF_Enrichment$Padj_left <- p.adjust(TF_Enrichment$Pval_left,method = "BH")
TF_Enrichment$Padj_right <- p.adjust(TF_Enrichment$Pval_right,method = "BH")
TF_Enrichment$Padj_both <- p.adjust(TF_Enrichment$Pval_both,method = "BH")

# Select factors that score at least in one group significance
tmp1 <- -log10(TF_Enrichment[apply(TF_Enrichment[,4:6],1,min) < 0.01,4:6])


mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(2,max(tmp1),length=51))
# Cluster the result, binary based on significant or not with a helping dataframe
tmp2 <- tmp1
# not significant gets 0
tmp2[tmp2< 2] <- 0
# significant gets 1
tmp2[tmp2> 2] <- 1
hr <- hclust(dist(tmp2,method="binary"))

## Figure K show enrichment in a heat map
p <- heatmap.2(as.matrix(tmp1),Colv=F,Rowv = as.dendrogram(hr), scale='none', col=mat_col,breaks=mat_col_breaks, trace='none', labRow=rownames(tmp1),labCol=colnames(tmp1))
# get the index of the factors which match factors from OC_TFs
tmp2 <- tmp1[p$rowInd,]
tmp2$Symbol <- rownames(tmp2)
tmp2[!tmp2$Symbol %in% OC_TFs, 'Symbol'] <- ""

# Make a grid for the alignment of the OC_TF factors in the heatmap. Lower and upper value reflects the rows of the heatmap and can be adjusted in illustrator.
plot(0,0,pch="", xlim=c(0,4),ylim=c(0,nrow(tmp2)+1))
# lower line --> lower end of heatmap
lines(c(1,2), c(1,1))
# upper line --> upper end of heatmap
lines(c(1,2), c(nrow(tmp2),nrow(tmp2)))
# All new lines correspond to the rows of the OC_TFs in the heatmap
for (i in 1:nrow(tmp2)){
  if(tmp2[i,'Symbol']!=""){
    lines(c(1,2), c(i,i))
    text(2,i,tmp2[i,'Symbol'])
  } else {}
}

rm(scOC,p,TF_Enrichment,tmp, tmp1, tmp2,Targets,i,k,l,m,mat_col, mat_col_breaks,n,OC_TFs,hr,mat,q)


