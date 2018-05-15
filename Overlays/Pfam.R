## PFAM enrichment
# Nov 30, 2017

setwd("/Volumes/Analysis1/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/")

library(data.table)
library(dplyr)

GetArray<-function(){
  #Array<-fread("Spart_array_1_info.csv",data.table=FALSE,na.strings = "",stringsAsFactors = FALSE)
  
  Array<-Array[,c(1,3,8)]
  library(stringr)
  Array.spl<-data.frame(str_split_fixed(Array$Functionnal_Annotation_IDs,",",10))
  Array.spl[Array.spl==""]<-NA
  library(tidyr)
  Array.spl<-gather_(Array.spl,"var","pfams",colnames(Array.spl))
  
  Array<-data.frame(cbind(Array[,1],Array.spl[,2]))
  colnames(Array)<-c("contig","pfam")
  Array<-mutate_all(Array,as.character)
  rm(Array.spl)
  return(Array)
}
#Array<-GetArray()

Array<-fread("pfam2.txt",data.table=FALSE,na.strings = "",skip=3,stringsAsFactors = FALSE)

Array<-Array[,c(4,2)]
Array$V4<-as.character(Array$V4)
Array$V2<-as.character(Array$V2)

library(stringr)

Array.temp<-data.frame(str_split_fixed(Array[,2],fixed("."),2))
Array$V2<-Array.temp$X1
rm(Array.temp)
colnames(Array)<-c("contig","pfam")
Array<-mutate_all(Array,as.character)


GetEnrichment<-function(Module=fread("Overlays/Pass5/all_oil_responsive2.csv",data.table = FALSE,stringsAsFactors = FALSE)){
  #Module<-fread("Final_ClusterPVCA/brown.csv",stringsAsFactors = FALSE)
  #Module<-Module[,c(21:23)]
  #colnames(Module)[1]<-"contig"
  #Module<-fread("Overlays/Pass5/all_oil_responsive2.csv",data.table = FALSE,stringsAsFactors = FALSE)
  colnames(Module)[2]<-"contig"
  Module<-merge(Module,Array,by="contig",all.x=TRUE)
  
  # Remove NA
  Module<-data.frame(Module$pfam)
  Module<-filter(Module,!is.na(Module[,1]))
  Array<-data.frame(Array$pfam)
  Array<-filter(Array,!is.na(Array[,1]))
  
  ## Fisher's exact test
  
  Enrichment<-function(mod=Module,universe=Array){
    PFams<-unique(as.character(mod[,1]))
    Report<- data.frame(Pfam=character(),Pvalue=numeric(),stringsAsFactors=FALSE)
    
    Iteration<-function(i){
      ToTest<-PFams[i]
      GetTable<-function(x){
        Table<-data.frame()
        Table[1,1]<-sum(str_count(Module[,1], x))
        Table[1,2]<-nrow(Module)
        Table[2,1]<-sum(str_count(Array[,1], x))
        Table[2,2]<-nrow(Array)
        return(data.matrix(Table))
      }
      Table<-GetTable(ToTest)
      Outcome<-fisher.test(Table)$p.value
      ToReport<-data.frame(ToTest,Outcome)
      colnames(ToReport)<-colnames(Report)
      ToReport<-rbind(Report,ToReport)
      return(ToReport)
    }
    
    library(parallel)
    # start cluster
    parallelCluster <- parallel::makeCluster(detectCores()-1, type = "FORK")
    print(parallelCluster)
    
    templist<-parLapply(parallelCluster,1:length(PFams), Iteration)
    
    # Shutdown cluster neatly
    if(!is.null(parallelCluster)) {
      parallel::stopCluster(parallelCluster)
      parallelCluster <- c()
    }
    
    report2<-rbindlist(templist)
    return(report2)
    
  }
  
  Outcomes<-Enrichment()
  
  library(multtest)
  Corrections<-mt.rawp2adjp(Outcomes$Pvalue)
  Corrections$h0.TSBH
  Corrections<-data.frame(Corrections$index,Corrections$adjp)
  rows<-as.numeric(as.character(filter(Corrections, BH <= 0.05)[,1]))
  print(Module[rows,])
  
}
GetEnrichment()

Mod<-fread("Final_ClusterPVCA/brown.csv",stringsAsFactors = FALSE)
Mod<-Mod[,c(22,21,23)]
GetEnrichment(Module=Mod)
