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


### Figure 6F
# the following files are provided in OSF https://osf.io/9xys4/
SSTR2 <- read.delim("SSTR2_Signaling_Resorption.txt",h=T)
bp <- boxplot(SSTR2[,grep("cAMP_",colnames(SSTR2))], names=c("1nMSST","10nMSST","100nMSST","1nMSST_100nMBIM","10nMSST_100nMBIM","100nMSST_100nMBIM","1nMOCT","10nMOCT","1nMOCT_100nMBIM","10nMOCT_100nMBIM"), ylab="cAMP inhibition")
for (i in 1:nrow(SSTR2)){
  lines(1:6,unlist(SSTR2[i,grep("cAMP_",colnames(SSTR2))][1:6]))
  lines(7:10,unlist(SSTR2[i,grep("cAMP_",colnames(SSTR2))][7:10]))
}
rm(SSTR2,bp)


### Figure 6G
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

### Figure 6H
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


### Figure 6I
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

### Figure 6J
# the following files are provided in OSF https://osf.io/9xys4/
FFAR4 <- read.delim("FFAR4_Signaling_Resorption.txt",h=T)
bp <- boxplot(FFAR4[,grep("IP1_",colnames(FFAR4))], names=c("0.1nMTUG","1nMTUG","10nMTUG"), ylab="IP1 levels")
for (i in 1:nrow(FFAR4)){
  lines(1:3,unlist(FFAR4[i,grep("IP1_",colnames(FFAR4))]))
}
rm(FFAR4,bp)


### Figure 6K
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

### Figure 6L
# the following files are provided in OSF https://osf.io/9xys4/
FFAR4 <- read.delim("FFAR4_Signaling_Resorption.txt",h=T)
bp <- boxplot(FFAR4[,grep("Res_",colnames(FFAR4))], names=c("Veh","1nMTUG"), ylab="% eroded surface/bone surface")
for (i in 1:nrow(FFAR4)){
  lines(1:2,unlist(FFAR4[i,grep("Res_",colnames(FFAR4))]))
}

title(paste("T.Test TUG Veh:",t.test(FFAR4[,"Res_Veh"],FFAR4[,"Res_TUG891_1nM"], paired=T)$p.value))
rm(FFAR4,bp)
