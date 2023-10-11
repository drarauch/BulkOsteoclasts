### Figure 5
# the following files are provided in OSF https://osf.io/9xys4/
Counts <- read.delim("Counts.txt",h=T)
Diff_ctr <- read.delim("Diff_ctr.txt",h=T)
colData_RNA_Diff <- data.frame("Sample"=colnames(Diff_ctr)[3:31],
                               "Donor"=gsub("_.*","",colnames(Diff_ctr)[3:31]),
                               "Timepoint"=gsub(".*_","",colnames(Diff_ctr)[3:31]))
rownames(colData_RNA_Diff)<-colData_RNA_Diff$Sample

### Figure 5A, 5B, 5C, 5D and 5E

# the following files are provided in OSF https://osf.io/9xys4/
# Download from the https://gpcrdb.org/
GPCR <- read.csv("GPCRTargets.csv")
Ligand <- read.csv("gpcrdb_endogenousligands.csv",h=T)
GWAS <-  read.csv("ebmd_con_indepen.csv",h=T)

GPCR <- unique(GPCR$HGNC.symbol)
Ligand <- unique(Ligand[,c("Gene","name","pKi_value","pKd_value","LigandPrGene")])
names(Ligand) <- c('Symbol','Ligand',"Ligand_pKi_value","Ligand_pKd_value","Ligand_Symbol")

# Expression bias in GPCRs compared to all other genes (dynamic ones are due to statistics not lowly expressed)
rownames(Counts) <- Counts$RefSeqID
Counts <- Counts[Diff_ctr$RefSeqID,]

# log transform average counts per day that are normalized to transcript length
Counts$d0 <- log2(apply(Counts[,grep("_d0",colnames(Counts))],1,mean)*1000/Counts$Length)
Counts$d2 <- log2(apply(Counts[,grep("_d2",colnames(Counts))],1,mean)*1000/Counts$Length)
Counts$d5 <- log2(apply(Counts[,grep("_d5",colnames(Counts))],1,mean)*1000/Counts$Length)
Counts$d9 <- log2(apply(Counts[,grep("_d9",colnames(Counts))],1,mean)*1000/Counts$Length)

# Keep the highest expression value for each gene
tmp <- data.frame('Symbol' =Counts$Symbol,'max' =apply(Counts[,c('d0','d2','d5','d9')],1,max))

# Plot distribution of maximum gene expression
plot(density(tmp[,2]))
lines(density(tmp[tmp$Symbol %in% GPCR,2]),lty=2,col="grey")
lines(density(tmp[tmp$Symbol %in% Ligand$Ligand_Symbol,2]),lty=2,col="grey")
# Make a threshold for the GPCR data to distangel expressed and not-expressed, with local minimum for both lines at 3
thres <- 3
abline(v=thres)

# Subset all genes to that threshold
Diff_ctr_thres <- Diff_ctr[Diff_ctr$Symbol %in% tmp[tmp$max>thres,'Symbol'],]

# Connect Ligand and Receptor
Ligand_Receptor <- merge(unique(Ligand[,c("Symbol","Ligand_Symbol")]),unique(data.frame("Symbol"=GPCR)),by="Symbol", all.y=T, all.x=T)
Ligand_Receptor <- unique(Ligand_Receptor[complete.cases(Ligand_Receptor),])
Ligand_Receptor <- Ligand_Receptor[Ligand_Receptor$Ligand_Symbol != "",]

# expressed GPCR, Ligand and Pairs
expGPCR <- GPCR[GPCR %in% Diff_ctr_thres$Symbol]
expLig <- unique(Ligand$Ligand_Symbol[Ligand$Ligand_Symbol %in% Diff_ctr_thres$Symbol])
expLR <- Ligand_Receptor[Ligand_Receptor$Symbol %in% Diff_ctr_thres$Symbol & Ligand_Receptor$Ligand_Symbol %in% Diff_ctr_thres$Symbol,]

## Figure 5A
barplot(c(
  length(expGPCR),
  length(expLig),
  nrow(expLR)))

# differential expressed GPCR, Ligand and Pairs, as enrichment over all differentially expressed ones
tmp <- Diff_ctr_thres[Diff_ctr_thres$clust != 0,'Symbol']
# Background probability regulated
reg <- length(tmp)/nrow(Diff_ctr_thres)
# Background probability non-regulated
nonreg <- (nrow(Diff_ctr_thres)-length(tmp))/nrow(Diff_ctr_thres)

## Figure 5B
barplot(c(
  # Single factors GPCR
  log2((length(expGPCR[expGPCR %in% tmp])/length(expGPCR))/reg),
  # Single factors Ligand
  log2((length(expLig[expLig %in% tmp])/length(expLig))/reg),
  # Pairs with one (relate to probability of having one regulated and one non-regulated in the background)
  log2((nrow(expLR[(expLR$Symbol %in% tmp | expLR$Ligand_Symbol %in% tmp) & !(expLR$Symbol %in% tmp & expLR$Ligand_Symbol %in% tmp),])/nrow(expLR))/(reg*nonreg)),
  # Pairs with both (relate to probability of having two regulated in the background)
  log2((nrow(expLR[expLR$Symbol %in% tmp & expLR$Ligand_Symbol %in% tmp,])/nrow(expLR))/(reg*reg))
))


## Figure 5C
# heatmap of expressed GPCRs and Ligands, order by cluster and save row number of each entry to connect both heatmaps GPCR
y <- rbind(
  Diff_ctr_thres[Diff_ctr_thres$Symbol %in% expGPCR,][order(Diff_ctr_thres[Diff_ctr_thres$Symbol %in% expGPCR,]$clust,apply(Diff_ctr_thres[Diff_ctr_thres$Symbol %in% expGPCR,][,c('d0','d2','d5','d9')],1,max)),],
  Diff_ctr_thres[Diff_ctr_thres$Symbol %in% expLig,][order(Diff_ctr_thres[Diff_ctr_thres$Symbol %in% expLig,]$clust,apply(Diff_ctr_thres[Diff_ctr_thres$Symbol %in% expLig,][,c('d0','d2','d5','d9')],1,max)),])
rownames(y) <- y$Symbol

# Add black rowsidecolor for GPCRs that are in cluster 0
library(fields)
library(RColorBrewer)
Mycol <- c('black',rainbow(8))
Mycol2 <- rev(designer.colors(n=100,col=brewer.pal(9, "Greens")))

# order samples according to day
y2 <-y[,colData_RNA_Diff[order(colData_RNA_Diff$Timepoint),'Sample']]
# Plot absolute expression
mat_col_breaks <- seq(0,max(y2),length=101)

length(expGPCR)
length(expLig)

library(gplots)
# heatmap receptor (GPCR)
heatp_R <- heatmap.2(as.matrix(y2[1:length(expGPCR),]), col=rev(Mycol2), RowSideColors = Mycol[y[1:length(expGPCR),'clust']+1],breaks=mat_col_breaks, Colv = F, Rowv = F, trace='none')
# heatmap ligand
heatp_L <- heatmap.2(as.matrix(y2[(length(expGPCR)+1):(length(expGPCR)+length(expLig)),]), col=rev(Mycol2), RowSideColors = Mycol[y[(length(expGPCR)+1):(length(expGPCR)+length(expLig)),'clust']+1],breaks=mat_col_breaks, Colv = F, Rowv = F, trace='none')

# combine receptor-ligand pairs
R <- data.frame(cbind(y[1:length(expGPCR),'Symbol'], rev(heatp_R$rowInd)))
names(R) <- c('Symbol','Receptor_index')
L <- data.frame(cbind(y[(length(expGPCR)+1):(length(expGPCR)+length(expLig)),'Symbol'], rev(heatp_L$rowInd)))
names(L) <- c('Ligand_Symbol','Ligand_index')
# Merge ligand index with Ligand data frame which contains the information of the receptor
L <- merge(L, unique(Ligand[,c("Ligand_Symbol","Symbol")]), by="Ligand_Symbol")
LR <- merge(L,R, by="Symbol", all.x=T, all.y=T)  
LR$Receptor_index <- as.numeric(LR$Receptor_index)
LR$Ligand_index <- as.numeric(LR$Ligand_index)
LR <- LR[order(LR$Receptor_index,LR$Ligand_index),]  

plot(0,0,pch="", xlim=c(1,4), ylim=c(-max(LR$Receptor_index, na.rm = T),0))
LR_tmp <- unique(LR[,c("Ligand_Symbol","Ligand_index")])
LR_tmp <- LR_tmp[complete.cases(LR_tmp),]
text(rep(1,nrow(LR_tmp)), -LR_tmp$Ligand_index, LR_tmp$Ligand_Symbol)
LR_tmp <- unique(LR[,c("Symbol","Receptor_index")])
LR_tmp <- LR_tmp[complete.cases(LR_tmp),]
text(rep(4,nrow(LR_tmp)), -LR_tmp$Receptor_index, LR_tmp$Symbol)     
LR_tmp <- LR[complete.cases(LR),]
for(i in 1:nrow(LR_tmp)){
  lines(c(1.5,2,3,3.5), c(-LR_tmp[i,"Ligand_index"],-LR_tmp[i,"Ligand_index"],-LR_tmp[i,"Receptor_index"],-LR_tmp[i,"Receptor_index"]))
}


## Figure 5D
# enrichment of GPCR among causal eBMD genes
# the following files are provided in OSF https://osf.io/9xys4/
GWAS <-  read.csv("ebmd_con_indepen.csv",h=T)
# Background 20595 genes in genome according to doi: 10.1093/nar/gkaa1080 related to GPCRdb
eBMDgene <- unique(GWAS[GWAS$P.NI < 5e-8,'C.GENE'])
tmp <- Diff_ctr_thres[Diff_ctr_thres$clust !=0,'Symbol']

barplot(c(
  log2((length(GPCR[GPCR %in% eBMDgene])/length(GPCR))/(length(eBMDgene)/20595)),
  log2((length(expGPCR[expGPCR %in% eBMDgene])/length(expGPCR))/(length(eBMDgene[eBMDgene %in% Diff_ctr_thres$Symbol])/length(Diff_ctr_thres$Symbol))),
  log2((length(expGPCR[expGPCR %in% tmp][expGPCR[expGPCR %in% tmp] %in% eBMDgene])/length(expGPCR[expGPCR %in% tmp]))/(length(eBMDgene[eBMDgene %in% tmp])/length(tmp)))))


## Figure 5E
# Sort GPCRs by signaling
# the following files are provided in OSF https://osf.io/9xys4/
GPCR_couplingsub <- read.delim("subtypes_coupling.txt",h=T)
GPCR_couplingsub <- GPCR_couplingsub[GPCR_couplingsub$Biosensor=="Mean",]

# Data is strange to work with in R, partly formated and not all rows have values is similar range (sometimes binart sometimes continious) - work with continious data
GPCR_couplingsub_con <- GPCR_couplingsub[apply(GPCR_couplingsub[,8:21],1,max)>3,]
GPCR_couplingsub_con <- unique(GPCR_couplingsub_con)

GPCR_couplingsub_con[,8:21][GPCR_couplingsub_con[,8:21]=="1'"] <- -1
GPCR_couplingsub_con[,8:21][GPCR_couplingsub_con[,8:21]=="2'"] <- -1

# Here GPCRdb has used uniprot data and not gene names
# the following files are provided in OSF https://osf.io/9xys4/
GPCR_Uniprot_Symbol <-read.delim("GPCR_coupling_Uniprot_Symbol.txt")
GPCR_couplingsub_con <- merge(GPCR_Uniprot_Symbol,GPCR_couplingsub_con,by="Uniprot")

GPCR_couplingsub_con <- GPCR_couplingsub_con[GPCR_couplingsub_con$Symbol %in% expGPCR[expGPCR %in% Diff_ctr_thres[Diff_ctr_thres$clust!=0,'Symbol']],]

GPCR_couplingsub_con <- merge(GPCR_couplingsub_con,Diff_ctr_thres[,c('Symbol','clust')],by="Symbol")
x=1
for (i in unique(GPCR_couplingsub_con$Uniprot)){
  tmp <- GPCR_couplingsub_con[GPCR_couplingsub_con$Uniprot==i,]
  tmp <- tmp[which.min(rowSums(tmp[,9:22]=="-")),]
  if(x==1){
    GPCR_couplingsub_con2 <- tmp
  } else{
    GPCR_couplingsub_con2 <- rbind(GPCR_couplingsub_con2,tmp)
  }
  x <- x+1
}

GPCR_couplingsub_con <- GPCR_couplingsub_con2
GPCR_couplingsub_con[,9:22][GPCR_couplingsub_con[,9:22]=="-"] <- -1
rownames(GPCR_couplingsub_con) <- GPCR_couplingsub_con$Symbol

GPCR_couplingsub_con[,9:22] <- apply(GPCR_couplingsub_con[,9:22],2,as.numeric)
# Show the sensor data in a heatmap
heatmap.2(as.matrix(GPCR_couplingsub_con[,9:22]),trace="none", col=c('black','white',designer.colors(n=99,col=brewer.pal(9,"Greens")[2:9])),breaks=c(-1,0,seq(5,10,length=100)), Colv = F)

rm(mat_col_breaks, Mycol,Mycol2,reg,nonreg,x,thres,i,k,GPCR,expLig,expGPCR,eBMDgene,y2,y,tmp,R,LR_tmp,LR,Ligand_Receptor,Ligand,heatp_L,heatp_R, GWAS,GPCR_couplingsub,GPCR_couplingsub_con,GPCR_couplingsub_con2,GPCR_Uniprot_Symbol,Gene_groups,expLR,Diff_ctr_thres,L)
