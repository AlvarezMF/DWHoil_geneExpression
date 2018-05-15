## Overlaps
# Dec 19 2017

setwd("/Volumes/Analysis1/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/")

library(data.table)
library(dplyr)

Armand<-fread("SuppFile4.txt",data.table = FALSE,header = FALSE)
colnames(Armand)<-"locus"
MDP<-fread("overlays/pass5/MDPOverlaptable.csv",data.table = FALSE)

comb<-merge(Armand,MDP,by="locus",all=FALSE)

# lookup functions of overlapping genes

ov<-fread("SuppFile4.txt",data.table = FALSE,header = FALSE)

library(AnnotationDbi)
library(org.At.tair.db)
toget<-ov[,1]
toget<-toget[!is.na(toget)]
#fwrite(data.frame(toget),file="Genes.txt",sep="\t",quote=FALSE,row.names = FALSE)
mapped<-mapIds(org.At.tair.db, keys=toget, column=c("GENENAME"), keytype="TAIR")
mapped<-data.frame(names(mapped),mapped)
colnames(mapped)<-c("ATnumber","Description")
mapped<-mapped[!duplicated(mapped),]
#mapped<-filter(mapped,!is.na(Description))
fwrite(mapped,file = "FileS4.txt",sep="\t",quote = FALSE,row.names = FALSE)

# Dumas and Weisman enrichment

# GO (VPlant)
W<-fread('/Volumes/Analysis1/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/Overlays/Pass5/WeismanDiffExp.txt',data.table = FALSE)
colnames(W)<-c("locus","dir")
D<-fread('/Volumes/Analysis1/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/Overlays/Pass5/DumasDiffExp.txt',data.table = FALSE)
colnames(D)<-c("locus","dir")
ov<-rbind(W,D)
rm(W,D)
fwrite(data.frame(ov[,1]),file="DumasWeismanOverlap.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)

# Pfam
