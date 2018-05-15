#Script for condensing ATH_GO_SLIM annotations into a single row for JMP analysis
#Written by MFA
#July 5, 2016

setwd('~/documents/code/r/2015BPoil/Annotation/')
library(plyr)
library(biomaRt)

listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)

#read in data frame
df <- read.csv("ATH_GO_GOSLIM.csv", header=FALSE)


#biological process

#subset biological process
BioP <- subset(df, V8=="P")

#remove duplicates
#First subset the two columns that we need anyway
BioP2<- BioP[c("V1","V6")]
#then output
UniQ<-unique(BioP2)

#Merge based on AT numbers and GO numbers
GOmerge <- ddply(UniQ, .(V1), summarize, V6=paste(V6, collapse="//"))

#change GOmerge names
names(GOmerge)<- c("AT","GO")

#Load SpAlt annotation file
SpAlt.ann <- read.csv("BPma_annotation.csv", header=TRUE)

#Merge GoSLIM annotations with 
Merged <- merge(SpAlt.ann, GOmerge, by="AT")

#Write out the file
write.csv(Merged, "SpAlt_BioProc_annotation.csv", row.names=FALSE)


#####
