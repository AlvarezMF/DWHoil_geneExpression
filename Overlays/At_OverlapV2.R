#2016 At Phenanthrene Overlay with BP oil experiment
#Written by Mariano Alvarez
#June 17, 2016

setwd('~/documents/code/R/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/Overlays/Pass5')
library(limma)
library(VennDiagram)
library(dplyr)

# First, use only the genes that appear on both the ATH1 array and the Spartina array
ath1 <- read.delim("affy_ATH1_array_elements-2010-12-20.txt", sep='\t', header=TRUE)
SpAlt1 <- read.csv("Spart_array_annotation.csv", header=TRUE)
colnames(SpAlt1)[8]<-"locus"
SpAlt1$locus<-gsub("\\..*","",SpAlt1$locus)
ConsensusGenes<-merge(ath1,SpAlt1,by="locus", all = FALSE)
ConsensusGenes<-as.data.frame(ConsensusGenes$locus)
#remove duplicates
ConsensusGenes<-as.data.frame(ConsensusGenes[!duplicated(ConsensusGenes[,1]),])
colnames(ConsensusGenes)<-"locus"


#read in Weismann data and BP oil data
At.data <- read.csv("12870_2009_553_MOESM2_ESM.csv")
colnames(At.data)[2]<-"locus"
At.origL <- nrow(At.data)
Sa.dataDUP <- read.csv("all_oil_responsive2.csv")
colnames(Sa.dataDUP)[1]<-"locus"
#remove duplicate observations in S. alterniflora
Sa.data<-Sa.dataDUP[!duplicated(Sa.dataDUP[,1]),]
Sa.origL <- nrow(Sa.data)
print(c(At.origL, "original observations in Arabidopsis"))
print(c(Sa.origL, "original observations in Spartina"))
print(c(nrow(Sa.dataDUP)-Sa.origL, "duplicate observations removed in Spartina"))

# filter out observations of each that don't appear on the consensus gene list
At.consensus<-merge(At.data,ConsensusGenes,by="locus", all = FALSE)
Sa.consensus<-merge(Sa.data,ConsensusGenes,by="locus", all = FALSE)
print(c(nrow(At.consensus), "consensus observations in Arabidopsis", At.origL-nrow(At.consensus), "observations removed"))
print(c(nrow(Sa.consensus), "consensus observations in Spartina",Sa.origL-nrow(Sa.consensus), "observations removed"))

### Parse out overlapping genes
#make ID column
Sa.consensus$Spartina<-1
At.consensus$Arabidopsis<-1
#merge
Merged<-merge(Sa.consensus, At.consensus, by="locus", all=TRUE)
Merged<-Merged[,c(1,3,9)]
#OverlapTable<-OverlapTable[!is.na(OverlapTable$AGI),]
write.csv(Merged, "overlaptable.csv")

Merged[is.na(Merged)] <- 0
nrow(dplyr::filter(Merged,Arabidopsis==1 & Spartina==1))
nrow(dplyr::filter(Merged,Arabidopsis==1 & Spartina==0))
nrow(dplyr::filter(Merged,Arabidopsis==0 & Spartina==1))


rm("At.data","ath1","Merged","OverlapTable","Sa.data","Sa.merge")

















##### log2 ratios ----- CAN START FROM HERE

#setup data
Ov<-read.csv("overlappingAT.csv",header=FALSE)
colnames(Ov)<-"AT"
SAexp<-read.csv("medsummclean.csv")
ATsig<-read.csv("12870_2009_553_MOESM2_ESM.csv")
colnames(ATsig)[2]<-"AT"
PctNorm<-read.csv("pctlnorm_amr.csv")
SA<-merge(SAexp,PctNorm,by="ProbeName")

#remove non-overlapping genes
AT2<-merge(ATsig,Ov,by="AT",all=FALSE)
AT2<-AT2[!duplicated(AT2$AT),]
SA2<-merge(SA,Ov,by="AT",all=FALSE)
SA2<-SA2[!duplicated(SA2$AT),]

#Create SA log2ratio
SA2$OilMean <- rowMeans(subset(SA2, select = c(3:10,19)), na.rm = TRUE)
SA2$NoOilMean <- rowMeans(subset(SA2, select = c(11:18,20)), na.rm = TRUE)
SA2$SAlog2ratio<-(log2(as.matrix(SA2$NoOilMean)))-(log2(as.matrix(SA2$OilMean)))

#rename log2ratio for AT 
colnames(AT2)[4]<-"ATlog2ratio"

#merge
Ov2<-merge(SA2,AT2,by="AT", all=FALSE)
keeps <- c("AT","SAlog2ratio","ATlog2ratio")
Ov2<-Ov2[keeps]
#colnames(Ov2)<-c("AT","SAlog2ratio","ATlog2ratio")
#Ov2<-data.matrix(Ov2)
hist(as.numeric(Ov2$SAlog2ratio))
hist(as.numeric(Ov2$ATlog2ratio))
write.csv(Ov2,"Ov2.csv")





##########
Ov2<-read.csv("Ov2.csv")
Ov2<-Ov2[,-c(1)]
Ov3<-na.omit(Ov2)
colnames(Ov3)<-c("TAIR_number", "S_alterniflora", "A_thaliana")

library(ggplot2)
library(reshape2)
dfm <- melt(Ov3[,c("TAIR_number", "S_alterniflora", "A_thaliana")],id.vars = 1)
colnames(dfm)<-c("ID","Species","Log2ratio")

nrow(dfm)
dfm<-dfm[with(dfm, order(ID)), ]
subset1<-dfm[1:50,]
subset2<-dfm[51:100,]
subset3<-dfm[101:150,]

ggplot(subset1,aes(x = ID, y= Log2ratio)) + geom_bar(stat="identity",aes(fill = Species),position = "dodge")
ggplot(subset2,aes(x = ID, y= Log2ratio)) + geom_bar(stat="identity",aes(fill = Species),position = "dodge")
ggplot(subset3,aes(x = ID, y= Log2ratio)) + geom_bar(stat="identity",aes(fill = Species),position = "dodge")


#Ov3<-na.omit(Ov2)
#subset <- t(data.frame(Ov3$SAlog2ratio, Ov3$ATlog2ratio))
#barplot(subset, legend = c("SAlog2ratio", "ATlog2ratio"), names.arg=Ov3$AT, beside=TRUE)




#2016 MDP stress annotation overlap
#Written by Mariano Alvarez
#Oct 17, 2016

setwd('~/documents/code/R/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/Overlays/Pass5')
library(limma)
library(VennDiagram)
library(dplyr)

# First, use only the genes that appear on both the ATH1 array and the Spartina array
ath1 <- read.delim("affy_ATH1_array_elements-2010-12-20.txt", sep='\t', header=TRUE)
SpAlt1 <- read.csv("Spart_array_annotation.csv", header=TRUE)
colnames(SpAlt1)[8]<-"locus"
SpAlt1$locus<-gsub("\\..*","",SpAlt1$locus)
ConsensusGenes<-merge(ath1,SpAlt1,by="locus", all = FALSE)
ConsensusGenes<-as.data.frame(ConsensusGenes$locus)
#remove duplicates
ConsensusGenes<-as.data.frame(ConsensusGenes[!duplicated(ConsensusGenes[,1]),])
colnames(ConsensusGenes)<-"locus"

# Then, subset genes from the MDP stress annotation that appear in the list of consensus genes
MDP<-read.csv("MDPStress_GO.csv", header=TRUE)
MDP.init<-nrow(MDP)
MDP$ATnumber<-toupper(MDP$ATnumber)
colnames(MDP)[1]<-"locus"
MDP<-merge(ConsensusGenes,MDP,by="locus")
paste("Rows removed:",MDP.init-nrow(MDP))
paste("Rows left:", nrow(MDP))

# Load in oil-responsive Spartina genes
Sa.dataDUP <- read.csv("all_oil_responsive2.csv")
colnames(Sa.dataDUP)[1]<-"locus"
#remove duplicate observations in S. alterniflora
Sa.data<-Sa.dataDUP[!duplicated(Sa.dataDUP[,1]),]
Sa.origL <- nrow(Sa.data)
print(c(Sa.origL, "original observations in Spartina"))
print(c(nrow(Sa.dataDUP)-Sa.origL, "duplicate observations removed in Spartina"))

# Prep Spartina data for overlap comparison
Sa.data$SpartinaOil<-1
Sa.data<-Sa.data[,-c(2)]

#Make sure to only compare consensus genes that appear in new overlap
Sa.origL <- nrow(Sa.data)
Sa.data<-merge(ConsensusGenes,Sa.data,by="locus")
paste("Rows removed:",Sa.origL-nrow(Sa.data))
#paste("Rows left:", nrow(Sa.data)

# Merge Spartina and MDP
StressOverlap<-merge(Sa.data,MDP,by="locus", all=TRUE)
# Why do I have extra genes? I guess I'll remove them?
StressOverlap<-filter(StressOverlap, !is.na(BacterialPathogen))
# Recode Spartina NA as 0
StressOverlap[is.na(StressOverlap)] <- 0

# Look for super responders (genes with 1s for more than half of stresses, excluding SpartinaOil)
StressOverlap$BaseResponseSum<-rowSums(StressOverlap[,3:11])
SuperResponders<-filter(StressOverlap, BaseResponseSum>4)
SpartinaResponders<-filter(SuperResponders, SpartinaOil==1)

# make SpartinaResponders table
SRtable<-merge(SpartinaResponders,ath1,by="locus", all=FALSE)
write.csv(SRtable, "SRtable.csv")




























raw<-StressOverlap[,-c(14,15)]
factors<-colnames(raw)[2:13]

ggplot(raw, aes(locus,factors)) + 
  geom_bar(stat="identity", position = "dodge")



# Plotting a different way - don't use!
raw<-StressOverlap[,2:13]
freq=table(col(raw), as.matrix(raw))
Names=colnames(raw)     # create list of names
data=data.frame(cbind(freq),Names)   # combine them into a data frame

library(reshape2)
# melt the data frame for plotting
data.m <- melt(data, id.vars='Names')

library(ggplot2)
# plot everything
ggplot(data.m, aes(Names, value)) +   
  geom_bar(aes(fill = variable), position = "dodge", stat="identity")

