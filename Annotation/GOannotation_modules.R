## SpAlt enrichment
## Feb 1 2017

setwd('/Volumes/Analysis1/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/Overlays/Pass5')
library(dplyr)
library(AnnotationDbi)
library(org.At.tair.db)
library(topGO)
library(data.table)
library(Rgraphviz)
library(edgeR)
library(qvalue)
data(geneList)

background<-read.delim("pctlnorm_amr.txt",header=TRUE, sep = "\t")
selGenes<-dplyr::filter(background, !is.na(AT))
selGenes<-as.character(selGenes$AT)

WGCNA<-read.csv("/Volumes/Analysis1/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/Final_ClusterPVCA/MergedDataset.csv", header=TRUE)
colnames(WGCNA)[20]<-"contig"
WGCNA<-merge(WGCNA,background,by="contig",all.x = TRUE)

WGCNA<-WGCNA[,c(21:22,45)]
colnames(WGCNA)[3]<-"Pval"
WGCNA$Pval<-sub("<","",WGCNA$Pval)
WGCNA$Pval<-sub("\\..","0.",WGCNA$Pval)
WGCNA$Pval<-as.numeric(WGCNA$Pval)
colors<-as.list(levels(as.factor(WGCNA$moduleColor)))


WGCNA$LocusLinkID<-as.character(WGCNA$LocusLinkID)
WGCNA<-dplyr::filter(WGCNA, !is.na(LocusLinkID))
colnames(WGCNA)[1]<-"AtNumber"
#WGCNA<-merge(WGCNA,background,by="AtNumber",all.x = TRUE)

setwd('/Volumes/Analysis1/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/Annotation')

MakeTopGO<-function(MOD) {
  x<-WGCNA
  y<-x
  x<-filter(x, moduleColor == MOD)
  #x$AtNumber<-sub("\\..","",x$AtNumber)
  #x<-x[!duplicated(x$AtNumber),]
  TAIRs<-x$Pval
  names(TAIRs)<-x$AtNumber
  #selGenes <- sample(ls(org.At.tairGO))
  gene2GO <- lapply(mget(selGenes, envir = org.At.tairGO), names)
  gene2GO[sapply(gene2GO, is.null)] <- NA
  topGOobj <- new("topGOdata",
                  ontology = "BP",
                  allGenes = TAIRs, geneSel = topDiffGenes,
                  nodeSize = 10,
                  annot = annFUN.gene2GO, gene2GO = gene2GO)
  print("topGO object created")
  return(topGOobj)
}
GetGOTable<-function(topGOobj) {
  resultKS <- runTest(topGOobj, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(topGOobj, algorithm = "elim", statistic = "ks")
  resultFisher <- runTest(topGOobj, algorithm = "classic", statistic = "fisher")
  allRes<-tryCatch(GenTable(topGOobj, classicFisher = resultFisher,
                            classicKS = resultKS, 
                            elimKS = resultKS.elim,
                            orderBy = "classicKS", ranksOf = "classicFisher", topNodes = 500),
                   error = function(e) print("No results found"))
  return(allRes)
}
GetGOChart<-function(topGOobj) {
  resultKS <- runTest(topGOobj, algorithm = "classic", statistic = "ks")
  #resultKS.elim <- runTest(topGOobj, algorithm = "elim", statistic = "ks")
  #resultFisher <- runTest(topGOobj, algorithm = "classic", statistic = "fisher")
  # if you feel like plotting the relationship between the two algorithms
  #pValue.classic <- score(resultKS)
  #pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
  #gstat <- termStat(topGOobj, names(pValue.classic))
  #gSize <- gstat$Annotated / max(gstat$Annotated) * 4
  #colMap <- function(x) {
  #.col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  #return(.col[match(1:length(x), order(x))])
  #}
  #gCol <- colMap(gstat$Significant)
  #plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize, col = gCol)
  showSigOfNodes(topGOobj, score(resultKS), firstSigNodes = 5, useInfo = 'all')
}
GetGOTable2<-function(topGOobj) {
  #resultKS <- runTest(topGOobj, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(topGOobj, algorithm = "elim", statistic = "ks")
  #resultFisher <- runTest(topGOobj, algorithm = "classic", statistic = "fisher")
  allRes<-tryCatch(score(resultKS.elim),
                   error = function(e) NA)
  return(allRes)
}
GetKSstats<-function(input) {
  topGOobj<-MakeTopGO(input)
  resultKS <- runTest(topGOobj, algorithm = "classic", statistic = "ks")
  #x<-qvalue(score(resultKS))
  x<-data.frame(score(resultKS))
  x$module<-paste(input)
  x$GOterm<-rownames(x)
  return(x)
}
KStoQ<-function(input) {
  x<-tryCatch(qvalue(input[,1]), error = function(e) NA)
  input[,1]<-x["qvalues"]
  return(input)
}



brown<-GetKSstats("brown")
x<-qvalue(brown[,1])
x["qvalues"]
#x<-qvalue(score(resultKS))
#x$qvalues
#KS<-data.table(GetKSscores(brown))
#sum(KS<0.05)
q<-data.frame(GetQ.KSstats(brown))


brown<-MakeTopGO("brown")
resultKS <- runTest(brown, algorithm = "classic", statistic = "ks")
resultsDF<-score(resultKS)
resultsDF<-qvalue(resultsDF)$qvalues
#resultsDF<-p.adjust(resultsDF, method = "BH")
resultsDF[resultsDF<0.05]


#resultsDF<-sort(resultsDF)
#sig<-resultsDF[1:5]
pdf("brownGraph_KS.pdf",height=5,width = 6.5)
showSigOfNodes(brown, score(resultKS), firstSigNodes = 5, useInfo = 'all')
dev.off()
#printGraph(brown, resultKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
brownTable<-data.frame(GetGOTable(brown))
resultsDF_df<-data.frame(resultsDF)
resultsDF_df$GO.ID<-rownames(resultsDF_df)
brownTable2<-merge(brownTable,resultsDF_df,by="GO.ID", all.x=TRUE)
brownTable2<-brownTable2[,-c(6:9)]
colnames(brownTable2)[6]<-"Q.value"
brownTable2<-arrange(brownTable2,Q.value)
write.table(brownTable2,"brownTable.txt",sep = "\t",quote = FALSE,row.names = FALSE)


browntable<-GetGOTable(brown)
pdf("brownGOchart.pdf",width = 6, height = 6)
GetGOChart(brown)
dev.off()




#### use this compartmentalized version

Table1<-lapply(colors,GetKSstats)
Table2<-lapply(Table1,KStoQ)

Test1<-data.frame(Table1[[1]])
Test2<-data.frame(Table2[[1]])

Table3<-rbindlist(Table2)
colnames(Table3)[1]<-"KS_qvalue"
colnames(Table3)

Table4<-dplyr::filter(Table3, KS_qvalue <= 0.05)
Table4[,1]<-as.numeric(Table4[,1])
Table4[,1]<-round(Table4[,1], digits = 4)

library(GO.db)
Table4$GOname<-Term(Table4[,3])
#org.At.tair.db

setwd('/Volumes/External1/R/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/Annotation')
write.table(Table4,"ModuleGO_supp.txt", quote = FALSE, row.names = FALSE)



#### EXTRA CODE

GOtable<-function(input) {
  report <- data.frame(Module=character(),TopGOcat=character(),KS_Pval=numeric(),stringsAsFactors=FALSE)
  Obj<-MakeTopGO(WGCNA,input)
  Objtable<-tryCatch(GetGOTable(Obj),error = function(e) NA)
  result<-tryCatch(data.frame(paste(input),Objtable[1,2],Objtable[1,9]), error= function(e) data.frame(c(NA,NA)))
  colnames(result)<-tryCatch(colnames(report), error= function(e) NULL)
  report<-rbind(report,result)
  return(report)
}
GOtable2<-function(input) {
  #report <- data.frame(Module=character(),TopGOcat=character(),KS_Pval=numeric(),stringsAsFactors=FALSE)
  Obj<-MakeTopGO(WGCNA,input)
  x<-GetQ.KSstats(Obj)
  y<-data.frame(x,paste(input))
  return(y)
}
GOtable2("brown")

library(parallel)
# start cluster
parallelCluster <- parallel::makeCluster(detectCores()-1, type = "FORK")
print(parallelCluster)
Table1<-parLapply(parallelCluster,colors,)
# Shutdown cluster neatly
if(!is.null(parallelCluster)) {
  parallel::stopCluster(parallelCluster)
  parallelCluster <- c()
}
colnames<-c("Module","topGOcat")
#lapply(Table1, setNames, colnames)

Table2<-rbindlist(Table1)

setwd('/Volumes/External1/R/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/Annotation')
write.table(Table2,"ModuleGO.txt", quote = FALSE)

