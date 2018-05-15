# Diff Exp script for BP oil data
# Mariano Alvarez
# Aug 4 2016

#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
#biocLite("edgeR")
library(edgeR)
library(Biobase)
library(pnmath0)

setwd('/Volumes/Rdrive/R/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/EdgeR/')

#load in expression data
exprs <- as.matrix(read.csv("medsummclean.csv", header=TRUE, row.names=1))

#impute missing values
#library(impute)
#any(is.na(exprs))
#if (any(is.na(exprs)))
#  exprsImp <- impute.knn(exprs)$data
#any(is.na(exprsImp))

#create metadata
metadata <- data.frame(labelDescription=c("Oil exposure, 1", "Field plot", "Which State"),
                       rownames=c("Treatment", "Site", "State"))

#load phenotype data
pData = read.csv("pheno.csv")

#reassign first row of phenotype data as rownames
pData2 <- pData[,-1]
rownames(pData) <- pData[,1]
pData <- pData[,-c(1)]
remove(pData2)

#do the rownames of the phenotype data match the column names of expression data?
all(rownames(pData)==colnames(exprs))

#merge metadata and phenotype data
phenoData = new("AnnotatedDataFrame", data=pData, varMetadata=metadata)

#make ExpressionSet
CADset <- ExpressionSet(assayData=exprs, phenoData=phenoData)




#PCA
#source("https://bioconductor.org/biocLite.R")
#biocLite("pcaMethods")
library(pcaMethods)
PC<-pca(CADset, nPcs=5)
#slplot(PC)
df <- merge(scores(PC), pData(CADset), by=0)
library(ggplot2)
ggplot(df, aes(PC1, PC2, shape=Site, color=Treatment)) +
  geom_point() +
  xlab(paste("PC1", PC@R2[1] * 100, "% of variance")) +
  ylab(paste("PC2", PC@R2[2] * 100, "% of variance"))
