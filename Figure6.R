### Figure 6
# the following files are provided in OSF https://osf.io/9xys4/
Counts <- read.delim("Counts.txt",h=T)
Diff_ctr <- read.delim("Diff_ctr.txt",h=T)
colData_RNA_Diff <- data.frame("Sample"=colnames(Diff_ctr)[3:31],
                               "Donor"=gsub("_.*","",colnames(Diff_ctr)[3:31]),
                               "Timepoint"=gsub(".*_","",colnames(Diff_ctr)[3:31]))
rownames(colData_RNA_Diff)<-colData_RNA_Diff$Sample

### Figure 6A

GOI <- c("C5AR1")
par(pty="s")
for (i in GOI){
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
rm(a,b,d0,d2,d5,d9,GOI,i,k,tmp)


### Figure 6B
library(Seurat)
library(gplots)
# the following files are provided in OSF https://osf.io/9xys4/
scOC <- readRDS("ReadyToUse_scOC.rds")

FeaturePlot(scOC,'C5AR1')
rm(scOC)

### Figure 6C
# the following files are provided in OSF https://osf.io/9xys4/
C5AR1 <- read.delim("C5AR1.txt",h=T)

# left panel, osteoclast number
bp <- boxplot(C5AR1[,grep("Number_",colnames(C5AR1))], names=c("DMSO","Agonist","Ago+Antago","Antagonist"), ylab="OC numbers")
for (i in 1:nrow(C5AR1)){
  lines(1:4,unlist(C5AR1[i,grep("Number_",colnames(C5AR1))]))
}

t.test(C5AR1[,grep("Number_",colnames(C5AR1))][,1],C5AR1[,grep("Number_",colnames(C5AR1))][,2],paired = T)
t.test(C5AR1[,grep("Number_",colnames(C5AR1))][,2],C5AR1[,grep("Number_",colnames(C5AR1))][,3],paired = T)
t.test(C5AR1[,grep("Number_",colnames(C5AR1))][,1],C5AR1[,grep("Number_",colnames(C5AR1))][,4],paired = F)
rm(bp)

# right panel, nuclei per osteoclast
bp <- boxplot(C5AR1[,grep("NucOc_",colnames(C5AR1))], names=c("DMSO","Agonist","Ago+Antago","Antagonist"), ylab="Nuclei per OC")
for (i in 1:nrow(C5AR1)){
  lines(1:4,unlist(C5AR1[i,grep("NucOc_",colnames(C5AR1))]))
}

t.test(C5AR1[,grep("NucOc_",colnames(C5AR1))][,1],C5AR1[,grep("NucOc_",colnames(C5AR1))][,2],paired = T)
t.test(C5AR1[,grep("NucOc_",colnames(C5AR1))][,2],C5AR1[,grep("NucOc_",colnames(C5AR1))][,3],paired = T)
t.test(C5AR1[,grep("NucOc_",colnames(C5AR1))][,1],C5AR1[,grep("NucOc_",colnames(C5AR1))][,4],paired = F)
rm(C5AR1,bp,i)


### Figure 6D
# the following files are provided in OSF https://osf.io/9xys4/
C5AR1 <- read.delim("C5AR1.txt",h=T)
bp <- boxplot(C5AR1[,grep("Res_",colnames(C5AR1))], names=c("DMSO","Agonist","Ago+Antago","Antagonist"), ylab="% eroded surface/bone surface")
for (i in 1:nrow(C5AR1)){
  lines(1:4,unlist(C5AR1[i,grep("Res_",colnames(C5AR1))]))
}

t.test(C5AR1[,grep("Res_",colnames(C5AR1))][,1],C5AR1[,grep("Res_",colnames(C5AR1))][,2],paired = T)
t.test(C5AR1[,grep("Res_",colnames(C5AR1))][,2],C5AR1[,grep("Res_",colnames(C5AR1))][,3],paired = T)
t.test(C5AR1[,grep("Res_",colnames(C5AR1))][,1],C5AR1[,grep("Res_",colnames(C5AR1))][,4],paired = F)
rm(C5AR1,bp)


### Figure 6E
GOI <- c("FFAR4")
par(pty="s")
for (i in GOI){
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
rm(a,b,d0,d2,d5,d9,GOI,i,k,tmp)


### Figure 6F
# the following files are provided in OSF https://osf.io/9xys4/
FFAR4 <- read.delim("FFAR4_Signaling_Resorption.txt",h=T)
bp <- boxplot(FFAR4[,grep("IP1_",colnames(FFAR4))], names=c("0.1nMTUG","1nMTUG","10nMTUG"), ylab="IP1 levels")
for (i in 1:nrow(FFAR4)){
  lines(1:3,unlist(FFAR4[i,grep("IP1_",colnames(FFAR4))]))
}
rm(FFAR4,bp)


### Figure 6G
# the following files are provided in OSF https://osf.io/9xys4/
FFAR4 <- read.delim("FFAR4_Diff.txt",h=T)

# left panel Nuclei numbers
bp <- boxplot(FFAR4[,grep("Number_",colnames(FFAR4))], names=c("Veh","1nMTUG"), ylab="Osteoclast number")
for (i in 1:nrow(FFAR4)){
  lines(1:2,unlist(FFAR4[i,grep("Number_",colnames(FFAR4))]))
}
title(paste("T.Test TUG Veh:",t.test(FFAR4[,"Number_Veh"],FFAR4[,"Number_TUG_1nM"], paired=T)$p.value))
rm(bp)

# right panel Nuclei per OC
bp <- boxplot(FFAR4[,grep("NucOc_",colnames(FFAR4))], names=c("Veh","1nMTUG"), ylab="Nuclei per OC")
for (i in 1:nrow(FFAR4)){
  lines(1:2,unlist(FFAR4[i,grep("NucOc_",colnames(FFAR4))]))
}
rm(FFAR4,bp)


### Figure 6H
# the following files are provided in OSF https://osf.io/9xys4/
FFAR4 <- read.delim("FFAR4_Signaling_Resorption.txt",h=T)
bp <- boxplot(FFAR4[,grep("Res_",colnames(FFAR4))], names=c("Veh","1nMTUG"), ylab="% eroded surface/bone surface")
for (i in 1:nrow(FFAR4)){
  lines(1:2,unlist(FFAR4[i,grep("Res_",colnames(FFAR4))]))
}

title(paste("T.Test TUG Veh:",t.test(FFAR4[,"Res_Veh"],FFAR4[,"Res_TUG891_1nM"], paired=T)$p.value))
rm(FFAR4,bp)


### Figure 6I

GOI <- c("SSTR2")
par(pty="s")
for (i in GOI){
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

# Focus on day 5 and 9 only
GOI <- c("SSTR2")
par(pty="s")
for (i in GOI){
  tmp <-colData_RNA_Diff
  tmp$x <- NA
  tmp[tmp$Timepoint=='d5','x'] <- 1
  tmp[tmp$Timepoint=='d9','x'] <- 2
  d5 <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d5','Sample']])
  d9 <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Timepoint=='d9','Sample']])
  boxplot(d5,d9,ylab=paste("Counts",i), xaxt="none")
  axis(1,at=c(1:2),labels = c('D5','D9'))
  for(k in unique(colData_RNA_Diff$Donor)){
    a <- tmp[tmp$Donor==k,'x']
    b <- unlist(Counts[Counts$Symbol==i,colData_RNA_Diff[colData_RNA_Diff$Donor==k,'Sample']])
    lines(a,b)
  }
}
rm(a,b,d0,d2,d5,d9,GOI,i,k,tmp)


### Figure 6J
# the following files are provided in OSF https://osf.io/9xys4/
SSTR2 <- read.delim("SSTR2_Signaling_Resorption.txt",h=T)
bp <- boxplot(SSTR2[,grep("cAMP_",colnames(SSTR2))], names=c("1nMSST","10nMSST","100nMSST","1nMSST_100nMBIM","10nMSST_100nMBIM","100nMSST_100nMBIM","1nMOCT","10nMOCT","1nMOCT_100nMBIM","10nMOCT_100nMBIM"), ylab="cAMP inhibition")
for (i in 1:nrow(SSTR2)){
  lines(1:6,unlist(SSTR2[i,grep("cAMP_",colnames(SSTR2))][1:6]))
  lines(7:10,unlist(SSTR2[i,grep("cAMP_",colnames(SSTR2))][7:10]))
}
rm(SSTR2,bp)


### Figure 6K
# the following files are provided in OSF https://osf.io/9xys4/
SSTR2 <- read.delim("SSTR2_Diff.txt",h=T)

# left panel Nuclei numbers
bp <- boxplot(SSTR2[,grep("Number_",colnames(SSTR2))], names=c("Veh","100nMSST","10nMOCT","10nMOCT_100nMBIM"), ylab="Osteoclast number")
for (i in 1:nrow(SSTR2)){
  lines(1:4,unlist(SSTR2[i,grep("Number_",colnames(SSTR2))]))
}
rm(bp)

# right panel Nuclei per OC
bp <- boxplot(SSTR2[,grep("NucOc_",colnames(SSTR2))], names=c("Veh","100nMSST","10nMOCT","10nMOCT_100nMBIM"), ylab="Nuclei per OC")
for (i in 1:nrow(SSTR2)){
  lines(1:4,unlist(SSTR2[i,grep("NucOc_",colnames(SSTR2))]))
}
rm(SSTR2,bp)


### Figure 6L
# the following files are provided in OSF https://osf.io/9xys4/
SSTR2 <- read.delim("SSTR2_Signaling_Resorption.txt",h=T)
bp <- boxplot(SSTR2[,grep("Res_",colnames(SSTR2))], names=c("Veh","100nMSST","10nMOCT","10nMOCT_100nMBIM"), ylab="% eroded surface/bone surface")
for (i in 1:nrow(SSTR2)){
  lines(1:4,unlist(SSTR2[i,grep("Res_",colnames(SSTR2))]))
}

title(paste("T.Test SST Veh:",t.test(SSTR2[,"Res_Veh"],SSTR2[,"Res_SST_100nM"], paired=T)$p.value))
title(paste("T.Test Oct Veh:",t.test(SSTR2[,"Res_Veh"],SSTR2[,"Res_OCT_10nM"], paired=T)$p.value), line = 0)
title(paste("T.Test Oct Anta:",t.test(SSTR2[,"Res_OCT_10nM"],SSTR2[,"Res_OCT_10nM_Anta_100nM"], paired=T)$p.value), line = 0.8)
rm(SSTR2,bp)


### Figure 6M, 6N and 6O 
# the following files are provided in OSF https://osf.io/9xys4/
SSTR2_KD_RNA <- read.delim("SSTR2_KD_RNA.txt",h=T)
SSTR2_KD_resorption <- read.delim("SSTR2_KD_resorption.txt",h=T)

## Figure 6M knockdown efficiency
rownames(SSTR2_KD_RNA) <- SSTR2_KD_RNA$Donor
SSTR2_KD_RNA$Donor <- NULL
boxplot(SSTR2_KD_RNA,ylim=c(0,max(SSTR2_KD_RNA)), xlab="",ylab="Expression SSTR2")
for (i in 1:nrow(SSTR2_KD_RNA)){
  lines(1:2,unlist(SSTR2_KD_RNA[i,]))
}
axis(1,at=c(1,2),labels=c('siCTR','siSSTR2'),las=2)
title(paste("T.Test SST Veh:",t.test(SSTR2_KD_RNA[,"Negative"],SSTR2_KD_RNA[,"siSSTR2"], paired=T)$p.value))

## Figure 6N effect of knockdown on resorption
SSTR2_KD_resorption <- read.delim("SSTR2_KD_resorption.txt",h=T)
rownames(SSTR2_KD_resorption) <- SSTR2_KD_resorption$Donor
SSTR2_KD_resorption$Donor <- NULL
boxplot(SSTR2_KD_resorption,ylim=c(0,max(SSTR2_KD_resorption)), xlab="",ylab="Resorption")
for (i in 1:nrow(SSTR2_KD_resorption)){
  lines(1:4,unlist(SSTR2_KD_resorption[i,1:4]))
}
axis(1,at=c(1,2,3,4),labels=c('siCTR Veh','siCTR SST','siSSTR2 Veh','siSSTR2 SST'),las=2)
title(paste("T.Test SST Veh in sictr:",t.test(SSTR2_KD_resorption[,"Negative_wo"],SSTR2_KD_resorption[,"Negative_SST"], paired=T)$p.value))
title(paste("T.Test si ctr si SSTR2:",t.test(SSTR2_KD_resorption[,"Negative_wo"],SSTR2_KD_resorption[,"siSSTR2_wo"], paired=T)$p.value), line=0)
title(paste("T.Test SST Veh in siSSTR2:",t.test(SSTR2_KD_resorption[,"siSSTR2_wo"],SSTR2_KD_resorption[,"siSSTR2_SST"], paired=T)$p.value), line=0.8)



# Calculate delta values for expression and resorption, and align the two data frames
SSTR2_KD_RNA <- SSTR2_KD_RNA[rownames(SSTR2_KD_resorption),]
SSTR2_KD_resorption$Delta_Negative <- (SSTR2_KD_resorption[,1]-SSTR2_KD_resorption[,2])/SSTR2_KD_resorption[,1]
SSTR2_KD_resorption$Delta_siSSTR2 <- (SSTR2_KD_resorption[,3]-SSTR2_KD_resorption[,4])/SSTR2_KD_resorption[,3]


## Figure 6O Plot Effect of knockdown against changes in resorption
plot(SSTR2_KD_RNA[,1]-SSTR2_KD_RNA[,2],
     (SSTR2_KD_resorption$Negative_wo - SSTR2_KD_resorption$siSSTR2_wo)/SSTR2_KD_resorption$Negative_wo,
     xlab="SSTR2 knockdown",ylab="Change in resorption \n by knockdown (no treatment)")


## Figure 6P Plot Response to SST as function of SSTR2 expression 
plot(SSTR2_KD_RNA[,1],SSTR2_KD_resorption$Delta_Negative,
     xlim=c(min(SSTR2_KD_RNA[,1:2]),max(SSTR2_KD_RNA[,1:2])),
     ylim=c(min(SSTR2_KD_resorption[,5:6]),max(SSTR2_KD_resorption[,5:6])),
     xlab="SSTR2 expression levels", ylab="Anti-resorptive effect \n of SST treatment")
points(SSTR2_KD_RNA[,2],SSTR2_KD_resorption$Delta_siSSTR2,pch=2)
for(i in 1:nrow(SSTR2_KD_RNA)){
  lines(c(SSTR2_KD_RNA[i,1],SSTR2_KD_RNA[i,2]),c(SSTR2_KD_resorption[i,'Delta_Negative'],SSTR2_KD_resorption[i,'Delta_siSSTR2']),lty=2)
}

# Plot Delta of SST response in ctr and siSSTR2 versus knockdown effect
plot((SSTR2_KD_RNA[,1]-SSTR2_KD_RNA[,2]),
     (SSTR2_KD_resorption$Delta_Negative-SSTR2_KD_resorption$Delta_siSSTR2),
     xlab="Delta in SSTR2 expression \n upon knockdown", ylab="Loss of Anti-resorptive \n treatment effect")
abline(lm((SSTR2_KD_resorption$Delta_Negative-SSTR2_KD_resorption$Delta_siSSTR2)~I(SSTR2_KD_RNA[,1]-SSTR2_KD_RNA[,2])),lty=2)
summary(lm((SSTR2_KD_resorption$Delta_Negative-SSTR2_KD_resorption$Delta_siSSTR2)~(SSTR2_KD_RNA[,1]-SSTR2_KD_RNA[,2])))
