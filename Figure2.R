### Figure 2
# the following files are provided in OSF https://osf.io/9xys4/
Counts <- read.delim("Counts.txt",h=T)
Diff_ctr <- read.delim("Diff_ctr.txt",h=T)
colData_RNA_Diff <- data.frame("Sample"=colnames(Diff_ctr)[3:31],
                               "Donor"=gsub("_.*","",colnames(Diff_ctr)[3:31]),
                               "Timepoint"=gsub(".*_","",colnames(Diff_ctr)[3:31]))
rownames(colData_RNA_Diff)<-colData_RNA_Diff$Sample

### Figure 2A
# Graphs were put together in illustrator and examples were added based on personal selection

# Line plots
library(RColorBrewer)
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

par(mfrow=c(3,3),pty='s')
for(i in 1:8) {
  tmp1 <- t(scale(t(Diff_ctr[Diff_ctr$clust == i,c('d0','d2','d5','d9')])))
  Mycol <- addTrans(rainbow(8)[i],255*40 / dim(tmp1)[1])
  plot(0,0,pch=' ', xlim=c(1,10), ylim=c(-1.5,1.5), xlab="", ylab="", yaxt="n", xaxt="n", main=paste(dim(tmp1)[1]))
  mtext("scaled Signal", side=2, line=2.5, cex=0.5)
  for(j in 1:dim(tmp1)[1]){
    lines(c(1,3,6,10), tmp1[j,], col=Mycol)
  }
  lines(c(1,3,6,10), apply(tmp1,2, function(x){median(x,na.rm=T)}), lwd=3)
  axis(2, at=c(-1,0,1), lab=c(-1,0,1), cex.axis=1)
  axis(1, at=c(1,3,6,10), lab=c("D1","D3","D6","D10"), las=2)
}

# Heatmap takes a lot of time with heatmap.2, as alternative try pheatmap
dev.off()
library(pheatmap)
library(scales)
library(fields)
library(pheatmap)
library(gplots)

y <- Diff_ctr[!Diff_ctr$clust==0,]
y <- y[order(y$clust),]

Mycol <- rainbow(length(unique(y$clust)))
Mycol2 <- designer.colors(n=100,col=brewer.pal(9, "Spectral"))
diff_help <- as.character(colData_RNA_Diff[order(colData_RNA_Diff$Timepoint,colData_RNA_Diff$Donor),'Sample'])
y2 <-t(scale(t(y[,diff_help])))
y2[y2 > 2] <- 2
y2[y2 < -2] <- -2
rownames(y2) <- y$Symbol
heatmap.2(as.matrix(y2), col=rev(Mycol2), RowSideColors = Mycol[y$clust], Colv = F, Rowv = F, trace='none')
rm(diff_help,y2,y,Mycol,Mycol2,i,j)


### Figure 2B
library(goseq)
library(gplots)
# gene ontology results can vary with package updates, here we used goseq_1.42.0 with geneLenDataBase_1.26.0 

# Make a list with gene groups
Gene_groups <- list()
for(i in 1:8){
  Gene_groups[[i]] <- Diff_ctr[Diff_ctr$clust==i,'Symbol']
  names(Gene_groups)[i] <- paste("clust_",i,sep="")
}

# Define background genes (all detected in RNA-seq that were used as input for DEseq analysis)
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

# Select terms of interest for plotting and arrange them according to clusters and if there are common terms across clusters (manually)
BP <- unique(c(
  'GO:0055114','GO:0007005','GO:0010257','GO:0042773','GO:0042407','GO:0019395','GO:0034470',
  'GO:0016126','GO:0045540','GO:0042254','GO:0006364','GO:0010467','GO:0008033','GO:0042773',
  'GO:0051301','GO:0048285','GO:0007049','GO:0007010','GO:0007059','GO:0033108','GO:0006119',
  'GO:0022008', 'GO:0032543','GO:0070085','GO:0007156','GO:0017004','GO:0044782','GO:0007010',
  'GO:0016071','GO:0006397','GO:0006351','GO:0008380','GO:0006325','GO:0006366','GO:0016570',
  'GO:0045087','GO:0030097','GO:0046649','GO:0032635','GO:0002224',
  'GO:0043299','GO:0032606','GO:0036230','GO:0007296','GO:0000278','GO:0007059','GO:0000281',
  'GO:0016043','GO:0007154','GO:0007155','GO:0048812','GO:0030029','GO:0030865','GO:0016477',
  'GO:0090383','GO:0030198','GO:0045851','GO:0046849','GO:0045047','GO:0006612','GO:0000187',
  'GO:0006955','GO:0010467','GO:0006397','GO:0008380','GO:0045321','GO:0006915','GO:0045087',
  'GO:0007249','GO:0000165','GO:0032635','GO:0032612','GO:0032615','GO:0045444','GO:0032623',
  'GO:0060070'))

# Extract p values for categories of interest and log transform them
p <- cbind(GO_cluster[GO_cluster$category %in% BP,1:2],-log10(GO_cluster[GO_cluster$category %in% BP,c(4:ncol(GO_cluster))]))

# Order the enrichments based on the individual clusters for nice representation 
p2 <- -log10(GO_cluster[GO_cluster$category %in% BP,c(4:ncol(GO_cluster))])
rownames(p2) <- p$category
rownames(p) <- p$category
p2[p2< -log10(0.01)] <- 0
p2[p2 >= -log10(0.01)] <- 1
tmp1 <- p2[p2[,1]==1,]
tmp1 <- tmp1[order(tmp1[,2],tmp1[,3],tmp1[,4],tmp1[,5],tmp1[,6],tmp1[,7],tmp1[,8]),]
tmp2 <- p2[p2[,1]==0 & p2[,2]==1,]
tmp2 <- tmp2[order(tmp2[,3],tmp2[,4],tmp2[,5],tmp2[,6],tmp2[,7],tmp2[,8]),]
tmp3 <- p2[p2[,1]==0 & p2[,2]==0 & p2[,3]==1,]
tmp3 <- tmp3[order(tmp3[,4],tmp3[,5],tmp3[,6],tmp3[,7],tmp3[,8]),]
tmp4 <- p2[p2[,1]==0 & p2[,2]==0 & p2[,3]==0 & p2[,4]==1,]
tmp4 <- tmp4[order(tmp4[,5],tmp4[,6],tmp4[,7],tmp4[,8]),]
tmp5 <- p2[p2[,1]==0 & p2[,2]==0 & p2[,3]==0 & p2[,4]==0 & p2[,5]==1,]
tmp5 <- tmp5[order(tmp5[,6],tmp5[,7],tmp5[,8]),]
tmp6 <- p2[p2[,1]==0 & p2[,2]==0 & p2[,3]==0 & p2[,4]==0 & p2[,5]==0 & p2[,6]==1,]
tmp6 <- tmp6[order(tmp6[,7],tmp6[,8]),]
tmp7 <- p2[p2[,1]==0 & p2[,2]==0 & p2[,3]==0 & p2[,4]==0 & p2[,5]==0 & p2[,6]==0 & p2[,7]==1,]
tmp7 <- tmp7[order(tmp7[,8]),]
p2 <- rbind(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7)

# make a heatmap with the manual order from above, since R cannot plot all names on rows, extract the names in a txt file to manually add them to the heatmap in illustrator
library(gplots)
library(fields)
library(scales)
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.01),max(p[,3:ncol(p)]),length=51))
heatp <- heatmap.2(as.matrix(p[rownames(p2),3:ncol(p)]),main="GO cluster BP",Rowv= F,dendrogram = 'none',  Colv=F, scale='none', col=mat_col,breaks=mat_col_breaks, trace='none', labRow=p[rownames(p2),2],labCol=colnames(p)[3:ncol(p)] )
write.table(p[rownames(p2),'term'][rev(heatp$rowInd)],file="GO_BP_names.txt",quote=F, col.names = F, sep="\n",row.names = F)

rm(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,p2,p,mat_col,mat_col_breaks,hg19.geneSymbol.LENGTH,heatp,GO.BP_Cl,GO_cluster,GO,Gene_groups,pwf_Cl,y,y2,BP,diff_help,tmp_Cl,tmp_All,Cl,tmp)


### Figure 1C
library(org.Hs.eg.db)
library(clusterProfiler)
library(goseq)
library(reactome.db)

# Make a list with gene groups
Gene_groups <- list()
for(i in 1:8){
  Gene_groups[[i]] <- Diff_ctr[Diff_ctr$clust==i,'Symbol']
  names(Gene_groups)[i] <- paste("clust_",i,sep="")
}

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
p <- Enrichment[apply(Enrichment[,3:(2+length(Gene_groups))],1,min)<0.001,]

#Sort the results by cluster for nice visual representation
p2 <- cbind(p[,1:2],-log10(p[,c(3:10)]))
p2 <- p2[order(-p2[,3]),]
p2 <- p2[p2[,3]> -log10(0.001),]

for (i in 4:10){
  p3 <- cbind(p[,1:2],-log10(p[,c(3:10)]))
  p3 <- p3[order(-p3[,i]),]
  p3 <- p3[p3[,i]> -log10(0.001),]
  p2 <- rbind(p2,p3)
}
# Screen for redundancy
p2$rep <- duplicated(p2$name)
p2 <- p2[p2$rep == "FALSE",]

# make a heatmap with the manual order from above, since R cannot plot all names on rows, extract the names in a txt file to manually add them to the heatmap in illustrator
library(gplots)
library(fields)
library(scales)
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.01),max(p2[,3:10]),length=51))
heatp <- heatmap.2(as.matrix(p2[,3:10]),main="GO cluster", Rowv = F, Colv=F, dendrogram='none', scale='none', col=mat_col,breaks=mat_col_breaks, trace='none',labRow = p2$name )
write.table(p2$name,file="Txt_Reactome_Cluster.txt", quote=F, col.names=F, row.names=F,sep="\n")

rm(heatp,p2,p3,p,Enrichment,hg19.refGene.LENGTH,Metabolism,All,Cl1.EA,Cl1,Cl1.nullp,Convert, Gene_groups,Pathways,Reactome,Relation,Current,i, mat_col, mat_col_breaks,New)

### Figure 2D

# the following files are provided in OSF https://osf.io/9xys4/
# File to convert Mouse Symbols into Human
Human_Mouse <- read.delim("Ensemble_SYMBOL_Mouse_Human.txt",h=T)
# Mouse phenotype data have been downloaded from IMPC (International Mouse phenotyping consortium) http://ftp.ebi.ac.uk/pub/databases/impc/
Genes_tested <- read.csv("procedureCompletenessAndPhenotypeHits.csv", h=T)
Mouse_Phenotypes <- read.csv("genotype-phenotype-assertions-ALL.csv")

# Identify interesting phenotypes
tmp <- Mouse_Phenotypes[Mouse_Phenotypes$top_level_mp_term_name == "skeleton phenotype",]
tmp2 <- unique(tmp$mp_term_name)
tmp2

Phenotypes <- list(
  c('increased bone mineral content'),
  c('decreased bone mineral content'),
  c('increased bone mineral density'),
  c('decreased bone mineral density'),
  c('abnormal bone structure'),
  c('abnormal bone mineralization'),
  c('abnormal clavicle morphology','abnormal joint morphology','abnormal lumbar vertebrae morphology','abnormal pelvic girdle bone morphology','abnormal rib morphology',
    'abnormal sacral vertebrae morphology','abnormal scapula morphology','abnormal sternum morphology','abnormal vertebrae morphology','abnormal vertebral arch morphology'))

names(Phenotypes) <- c('BMC_up','BMC_down','BMD_up','BMD_down','AbnormalBoneStructure','AbnormalBoneMineralization','AbnormalXMorphology')

Phenotypes_Parameters <- list()
for(i in 1:length(Phenotypes)){
  print(unique(Mouse_Phenotypes[Mouse_Phenotypes$mp_term_name %in% Phenotypes[[i]],'parameter_name']))
}
Phenotypes_Parameters[[1]] <- c("BMC/Body weight","Bone Mineral Content \\(excluding skull\\)","Bone Mineral Content")                       
Phenotypes_Parameters[[2]] <- c("BMC/Body weight","Bone Mineral Content \\(excluding skull\\)","Bone Mineral Content")                       
Phenotypes_Parameters[[3]] <- c("Bone Mineral Density \\(excluding skull\\)")
Phenotypes_Parameters[[4]] <- c("Bone Mineral Density \\(excluding skull\\)")
Phenotypes_Parameters[[5]] <- c("Bone Area","Bone area \\(BMC/BMD\\)","Bone")
Phenotypes_Parameters[[6]] <- c("BMC/Body weight")
Phenotypes_Parameters[[7]] <- c("Shape of vertebrae","Processes on vertebrae","Joints","Sternum","Scapulae","Shape of ribs","Clavicle","Pelvis","Thoracic processes","Lumbar processes","Cervical processes","Number of pelvic vertebrae","Sacral processes","Caudal processes","Number of lumbar vertebrae")

names(Phenotypes_Parameters) <- names(Phenotypes)

# Collect Human symbols of all hit genes from each Phenotype in a list 
Gene_groups <- list()
for (i in 1:length(Phenotypes)){
  Gene_groups[[i]] <- unique(Human_Mouse[Human_Mouse$SYMBOL_Mouse %in% unique(Mouse_Phenotypes[Mouse_Phenotypes$mp_term_name %in% Phenotypes[[i]] & Mouse_Phenotypes$zygosity == "homozygote",c("marker_symbol")]),'SYMBOL_Human'])
}
names(Gene_groups) <- names(Phenotypes)

# Collect Human symbols of all tested genes from each Phenotype
Gene_groups_tested <- list()
for (i in 1:length(Phenotypes)){
  Gene_groups_tested[[i]] <- unique(Human_Mouse[Human_Mouse$SYMBOL_Mouse %in% unique(Genes_tested[grepl(paste(Phenotypes_Parameters[[i]],collapse="|"), Genes_tested$Parameter.Name..Successful) & Genes_tested$Zygosity=="homozygote",'Gene.Symbol']),'SYMBOL_Human'])
}

# Calculate enrichment between genes that affect the phenotypes and genes from the clusters using a hypergeometric test
mat <- matrix(NA, ncol=9, nrow=length(Gene_groups))
rownames(mat) <- names(Gene_groups)
colnames(mat) <- paste("clust_",c(0,1:8),sep="")

for(i in c(0,1:8)){
  for(k in 1:length(Gene_groups)){
    tmp_all <- Diff_ctr[Diff_ctr$Symbol %in% Gene_groups_tested[[k]],]
    tmp_k <- Gene_groups[[k]][Gene_groups[[k]] %in% tmp_all$Symbol]
    tmp_i <- tmp_all[tmp_all$clust == i,"Symbol"]
    mat[k,i+1] <- phyper(length(tmp_k[tmp_k %in% tmp_i])+1, length(tmp_k), length(tmp_all[!tmp_all$Symbol %in% tmp_k,'Symbol']), length(tmp_i), lower.tail = F)
    
  }
}
#plot result as heatmap
library(gplots)
library(fields)
library(scales)
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat_col_breaks <- c(0,seq(-log10(0.05),max(-log10(mat)),length=51))
# Do not plot cluster 0
heatmap.2(-log10(mat[,2:9]),Rowv= F,dendrogram = 'none',  Colv=F, scale='none', col=mat_col,breaks=mat_col_breaks, trace='none' )
rm(i,k,mat_col,mat_col_breaks,tmp_i,tmp_k,tmp2,tmp_all,tmp,Phenotypes_Parameters,Phenotypes,Mouse_Phenotypes,mat,Human_Mouse,Gene_groups, Gene_groups_tested, Genes_tested)


### Figure 2E
# the following files are provided in OSF https://osf.io/9xys4/
data_Iliac <- read.delim("ReadyToUse_E-MEXP-1618.txt",h=T)

GOI <- c("ACTN2")
par(mfrow=c(1,2))
for (i in GOI){
  if(i %in% data_Iliac$SYMBOL){
    a <- unlist(data_Iliac[data_Iliac$SYMBOL==i,grep("healthy_",colnames(data_Iliac))])
    b <- unlist(data_Iliac[data_Iliac$SYMBOL==i,grep("osteoporotic_",colnames(data_Iliac))])
    boxplot(a,b,main=paste(i,"in Iliac Crest"), names=c('healthy','osteoporotic'),las=2)
    title(data_Iliac[data_Iliac$SYMBOL==i,'Pval_Osteoporotic_Healthy'], line=0.5)
    
  } else{
    plot(0,0,pch="", xaxt="none", yaxt="none", ylab="", xlab="", main=paste(i,"not in data_Iliac"))
  }
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
rm(data_Iliac,a,b,d0,d2,d5,d9,GOI,i,k,tmp)

### Figure 2F

#PBMCs were not used for the final figure
# the following files are provided in OSF https://osf.io/9xys4/
data_Iliac <- read.delim("ReadyToUse_E-MEXP-1618.txt",h=T)
data_PBMC <- read.delim("ReadyToUse_GSE56815.txt",h=T)

Gene_groups_1 <- list()
Gene_groups_2 <- list()

# Define gene groups from PBMC and Iliac crest studies, there are very few genes in the Iliac crest study that score an adjusted p-value below 0.05 for which nominal p-value of 0.005 was chosen
Gene_groups_1[[1]] <- unique(data_PBMC[data_PBMC$adj.P.Val < 0.05 & data_PBMC$logFC < 0 & data_PBMC$Symbol %in% Diff_ctr$Symbol,'Symbol'])
Gene_groups_1[[2]] <- unique(data_PBMC[data_PBMC$adj.P.Val < 0.05 & data_PBMC$logFC > 0 & data_PBMC$Symbol %in% Diff_ctr$Symbol,'Symbol'])
Gene_groups_1[[3]] <- unique(data_Iliac[data_Iliac$Pval_Osteoporotic_Healthy < 0.005 & data_Iliac$LogFC_Osteoporotic_Healthy < 0 & data_Iliac$SYMBOL %in% Diff_ctr$Symbol,"SYMBOL"])
Gene_groups_1[[4]] <- unique(data_Iliac[data_Iliac$Pval_Osteoporotic_Healthy < 0.005 & data_Iliac$LogFC_Osteoporotic_Healthy > 0 & data_Iliac$SYMBOL %in% Diff_ctr$Symbol,"SYMBOL"])
names(Gene_groups_1) <- c('PBMC_down','PBMC_up','Iliac_down','Iliac_up')

# Define gene groups for clusters
Gene_groups_2[[1]] <- Diff_ctr[Diff_ctr$clust == 1,'Symbol']
Gene_groups_2[[2]] <- Diff_ctr[Diff_ctr$clust == 2,'Symbol']
Gene_groups_2[[3]] <- Diff_ctr[Diff_ctr$clust == 3,'Symbol']
Gene_groups_2[[4]] <- Diff_ctr[Diff_ctr$clust == 4,'Symbol']
Gene_groups_2[[5]] <- Diff_ctr[Diff_ctr$clust == 5,'Symbol']
Gene_groups_2[[6]] <- Diff_ctr[Diff_ctr$clust == 6,'Symbol']
Gene_groups_2[[7]] <- Diff_ctr[Diff_ctr$clust == 7,'Symbol']
Gene_groups_2[[8]] <- Diff_ctr[Diff_ctr$clust == 8,'Symbol']
names(Gene_groups_2) <- c('Cl1','Cl2','Cl3','Cl4','Cl5','Cl6','Cl7','Cl8')

# Calculate enrichment based on a hypergeometric test
mat <- matrix(NA, ncol=length(Gene_groups_1), nrow=length(Gene_groups_2))
colnames(mat) <-names(Gene_groups_1)
rownames(mat) <-names(Gene_groups_2)

for(i in 1:length(Gene_groups_1)){
  for(k in 1:length(Gene_groups_2)){
    tmp_i <- Gene_groups_1[[i]]
    tmp_k <- Gene_groups_2[[k]]
    mat[k,i] <- phyper(length(tmp_i[tmp_i %in% tmp_k]), length(tmp_k), length(Diff_ctr[!Diff_ctr$Symbol %in% tmp_k,'Symbol']), length(tmp_i), lower.tail = F)
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
rm(Gene_groups_1, Gene_groups_2,mat,i,k,mat_col,mat_col_breaks,tmp_i,tmp_k, data_Iliac,data_PBMC)

### Figure 2G

# Define TSS of genes from RNA-seq on hg19 coordinates
library(data.table)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# Get transcripts for hg19
GR <- data.frame(transcripts(txdb))
# Get symbols to UCSC IDs
GR2 <- select(Homo.sapiens,GR$tx_name, "SYMBOL","TXNAME")
GR <- merge(GR,GR2, by.x="tx_name", by.y="TXNAME")
rm(GR2)
GR <- GR[GR$seqnames %in% paste('chr',c(1:23,'X','Y','M'),sep=""),]
GR <- unique(GR[GR$SYMBOL %in% Diff_ctr$Symbol,c("SYMBOL","seqnames","start","strand")])
# Get lowest TSS for "+"-stranded
GR1 <- GR[GR$strand =="+",] 
GR1 <- GR1[order(GR1$SYMBOL,GR1$start),]
GR1$dup <- duplicated(GR1$SYMBOL)
GR1 <- GR1[GR1$dup =="FALSE",]
GR1$dup <- NULL
names(GR1) <- c("Symbol",'Chr',"TSS","Strand")
# Get highest TSS for "-"-stranded
GR2 <- GR[GR$strand =="-",] 
GR2 <- GR2[order(GR2$SYMBOL,-GR2$start),]
GR2$dup <- duplicated(GR2$SYMBOL)
GR2 <- GR2[GR2$dup =="FALSE",]
GR2$dup <- NULL
names(GR2) <- c("Symbol",'Chr',"TSS","Strand")
GR <- rbind(GR1,GR2)

rm(GR1,GR2, txdb)
# Merge GR and Diff_ctr to get the cluster information
GR <- merge(GR,Diff_ctr[,c('Symbol','clust')], by="Symbol")
write.table(GR,file="GR.txt", quote=F, col.names=T, row.names=F, sep="\t")

## Run overlap on server (normal PC has not enough power)
library(data.table)
library(GenomicRanges)
GR <- read.delim("GR.txt",h=T)

# the following files are provided in OSF https://osf.io/9xys4/
# Get GWAS summary statistics from http://www.gefos.org/?q=content/data-release-2018 doi:10.1038/s41588-018-0302-x
GWAS_summary <- fread("Biobank2-British-Bmd-As-C-Gwas-SumStats.txt",h=T)
GWAS_summary <- GWAS_summary[,c('RSID','CHR','BP','P.NI','BETA')]
names(GWAS_summary)[1] <- 'SNP'
GWAS_summary$Chr <- paste("chr",GWAS_summary$CHR, sep="")
GWAS_summary$Start <- as.numeric(GWAS_summary$BP)
GWAS_summary <- GWAS_summary[complete.cases(GWAS_summary$Start),c('SNP','P.NI','Chr','Start','BETA')]
GWAS_summary <- data.frame(GWAS_summary)
GWAS_summary$Pval <- as.numeric(GWAS_summary$P.NI)

# Do overlap of SNPs in window of 5 Mb
grGenes <- with(unique(GR[,c("Symbol","Chr","TSS")]) , GRanges(Chr, IRanges(start=TSS - 5000000, end=TSS + 5000000, names=Symbol)))
grSummary <- with(unique(GWAS_summary[,c("SNP","Chr","Start")]) , GRanges(Chr, IRanges(start=Start, end=Start, names=SNP)))

hits = findOverlaps(grGenes,grSummary)
tmp2 <- cbind(data.frame(ranges(grGenes)[queryHits(hits)]),data.frame(ranges(grSummary)[subjectHits(hits)]))
colnames(tmp2) <- c('TSSminus','TSSplus','Window','Symbol','SNP_Start','SNP_End','SNP_Length','SNP')
head(tmp2)

# Calculate distance bewteen TSS and SNP
tmp2$Distance <- tmp2$TSSminus+5000000 - tmp2$SNP_Start
tmp2 <- tmp2[,c('Symbol','SNP','Distance')]
tmp2 <- merge(tmp2[,c('Symbol','SNP','Distance')], GWAS_summary[,c("SNP","Pval","BETA")], by="SNP")
tmp2 <- merge( GR,tmp2[,c('Symbol','SNP','Pval','Distance')], by="Symbol")


#Chisquare test for enrichment with different distances from the TSS
mat2 <- matrix(NA, ncol=9, nrow=13)
x <- 1
for (k in c(50000,100000,250000,seq(500000,5000000, length=10))){
  for (i in c(0,1:8)){
    a <- length(tmp2[tmp2$clust==i & tmp2$Pval < 5E-8 & abs(tmp2$Distance) < k ,'SNP'])
    b <- length(tmp2[tmp2$clust==i & tmp2$Pval > 5E-8 & abs(tmp2$Distance) < k,'SNP'])
    c <- length(tmp2[tmp2$Pval < 5E-8 & abs(tmp2$Distance) < k ,'SNP'])
    d <- length(tmp2[tmp2$Pval > 5E-8 & abs(tmp2$Distance) < k ,'SNP'])
    M <- as.table(rbind(c(a, b), c(c, d)))
    dimnames(M) <- list(group = c("in_cluster", "in_gwascat"),
                        categ = c("HBMD_assoc","not_HBMD_assoc"))
    mat2[x,i+1] <- chisq.test(M)$p.value
  }
  x <- x+1
}
rownames(mat2) <- paste('Dist_',c(50000,100000,250000,seq(500000,5000000, length=10)),sep="")
colnames(mat2) <- paste('Clust_',c(0,1:8),sep="")
saveRDS(mat2,file="mat2.rds")

# Run on local machine again
mat <- -log10(readRDS("mat2.rds"))
library(fields)
library(gplots)
library(RColorBrewer)
mat_col <- c('white',designer.colors(n=50, col=c('plum1','darkmagenta')))
mat[mat>150] <- 150
mat_col_breaks <- c(0,seq(5,max(mat),length=51))
heatmap.2(as.matrix(mat),Rowv= F,dendrogram = 'none',  Colv=F, scale='none', col=mat_col,breaks=mat_col_breaks, trace='none', labRow=rownames(mat),labCol=colnames(mat) )

rm(GR, mat_col, mat_col_breaks, txdb)

### Figure 2H

# GSE152677 has processed data "GSE152677_DEG_DESeq2.xls" which is available at OSF https://osf.io/9xys4/
# 10.1016/j.bone.2019.07.022: Transcriptional Profiling of Intramembranous and Endochondral Ossification after Fracture in Mice (10.1016/j.bone.2019.07.022)
# the following files are provided in OSF https://osf.io/9xys4/
# File to convert Mouse Symbols into Human
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

# Calculate enrichment based on a hypergeometric test
mat <- matrix(NA,ncol=length(Fracture_list),nrow=9)
colnames(mat) <- names(Fracture_list)
rownames(mat) <- paste("clust",c(0,1,2,3,4,5,6,7,8),sep="_")
for (i in 1:length(Fracture_list)){
  for (k in c(0,1,2,3,4,5,6,7,8)){
    a <- Fracture_list[[i]]
    b <- Diff_ctr[Diff_ctr$clust==k, 'Symbol']
    mat[k+1,i] <- phyper(length(a[a %in% b]),length(b),length(Diff_ctr$Symbol)-length(b),length(a),lower.tail=FALSE)
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

rm(mat_col, mat_col_breaks,mat,a,b,i,k,names,x,Fracture_list,Human_Mouse,data)

