#2016 At Phenanthrene Overlay with BP oil experiment
#Written by Mariano Alvarez
#Nov 7, 2016

##### log2 ratios ----- CAN START FROM HERE
## this has been changed to just look at directionality, since our normalized
## Spartina data is going to mess up fold change calculations
setwd('/Volumes/Analysis1/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/Overlays/Pass5')
library(limma)
library(VennDiagram)
library(dplyr)


#setup data
SAexp<-read.csv("medsummclean.csv")
PctNorm<-read.csv("pctlnorm_amr.csv")
SAexp<-merge(SAexp,PctNorm,by="ProbeName")
SAexp<-SAexp[,c(20,2:19)]
colnames(SAexp)[1]<-"locus"

#remove non-overlapping genes
#SAexp<-merge(SAexp,ConsensusGenes,by="locus",all=FALSE)
#SAexp<-as.data.frame(SAexp[!duplicated(SAexp[,1]),])

#Create SA log2ratio
SAexp$OilMean <- (rowMeans(subset(SAexp, select = c(2:9,18)), na.rm = TRUE))
SAexp$NoOilMean <- (rowMeans(subset(SAexp, select = c(10:17,19)), na.rm = TRUE))
SAexp$SAlog2ratio<-SAexp$NoOilMean > SAexp$OilMean
#SAexp$SAlog2ratio<-log2(SAexp$SAlog2ratio)
SAexp$SAlog2ratio[SAexp$SAlog2ratio=="TRUE"] <- "Downregulated"
SAexp$SAlog2ratio[SAexp$SAlog2ratio=="FALSE"] <- "Upregulated"

SAexp<-SAexp[,c("locus","SAlog2ratio")]
#SAexp<-merge(SAexp,Sa.consensus,by="locus",all=FALSE)
#SAexp<-SAexp[,1:2]
SAexp<-as.data.frame(SAexp[!duplicated(SAexp[,1]),])

# Create Dumas average log2
DumasExp<-read.csv("DumasLog2.csv")
colnames(DumasExp)<-c("locus","30m","2h","4h","8h","24h")
DumasExp$DumasLog2ratio<-rowMeans(DumasExp[,2:ncol(DumasExp)])
DumasExp$DumasLog2ratio<-DumasExp$DumasLog2ratio < 0
DumasExp$DumasLog2ratio[DumasExp$DumasLog2ratio=="TRUE"] <- "Downregulated"
DumasExp$DumasLog2ratio[DumasExp$DumasLog2ratio=="FALSE"] <- "Upregulated"


DumasExp$locus<-toupper(DumasExp$locus)
DumasExp<-DumasExp[,c("locus","DumasLog2ratio")]
library(data.table)
fwrite(DumasExp,file="DumasDiffExp.txt",sep="\t",row.names=FALSE,quote = FALSE)
#DumasExp<-merge(DumasExp,Dumas.consensus,by="locus",all=FALSE)
#DumasExp<-DumasExp[,1:2]
#colnames(DumasExp)[2]<-"log2ratio"

# Create Weisman log2
Weisman <- read.csv("Weisman.csv")
colnames(Weisman)[2]<-"locus"
WeismanExp<-Weisman[,c("locus","log2ratio")]
colnames(WeismanExp)<-c("locus","WeisLog2ratio")
WeismanExp$WeisLog2ratio<-WeismanExp$WeisLog2ratio < 0
WeismanExp$WeisLog2ratio[WeismanExp$WeisLog2ratio=="TRUE"] <- "Downregulated"
WeismanExp$WeisLog2ratio[WeismanExp$WeisLog2ratio=="FALSE"] <- "Upregulated"


WeismanExp$locus<-as.character(WeismanExp$locus)
s <- strsplit(WeismanExp$locus, split = ",")
WeismanExp<-data.frame(V1 = rep(WeismanExp$WeisLog2ratio, sapply(s, length)), V2 = unlist(s))
WeismanExp<-WeismanExp[,2:1]
colnames(WeismanExp)<-c("locus","WeisLog2ratio")
#colnames(WeismanExp)<-c("locus","log2ratio")
fwrite(WeismanExp,file="WeismanDiffExp.txt",sep="\t",row.names=FALSE,quote = FALSE)


#ConsensusCheck<-merge(WeismanExp,DumasExp,by="locus",all=TRUE)
#ConsensusCheck2<-dplyr::filter(ConsensusCheck, !is.na(log2ratio.x) & !is.na(log2ratio.y))
#Error<-dplyr::filter(ConsensusCheck2, log2ratio.x != log2ratio.y)

#AtExp<-rbind(WeismanExp,DumasExp)
#AtExp<-as.data.frame(AtExp[!duplicated(AtExp[,1]),])

#sum(is.na(AtExp$log2ratio))

#WeismanExp<-merge(WeismanExp,Weisman.consensus,by="locus",all=FALSE)
#WeismanExp<-WeismanExp[,1:2]

#load 187 overlap
Xenome <- read.table("DE_At_common_loci.txt", sep = "\t", colClasses = "character")
Xenome<-data.frame(toupper(Xenome[,1]))
colnames(Xenome)<-"locus"

#merge
Log2Overlap<-merge(SAexp,Xenome,by="locus",all=FALSE)
#Log2Overlap<-merge(Log2Overlap,AtExp,by="locus",all.x = TRUE)

Log2Overlap<-merge(Log2Overlap,WeismanExp,by="locus",all.x=TRUE)
Log2Overlap<-merge(Log2Overlap,DumasExp,by="locus",all.x=TRUE)

#Log2Overlap<-dplyr::filter(Log2Overlap, SAlog2ratio != "NaN")
Log2Overlap<-Log2Overlap[with(Log2Overlap, order(locus)), ]

write.csv(Log2Overlap,"log2overlap4.csv")

WSAoverlap<-dplyr::filter(Log2Overlap, SAlog2ratio == WeisLog2ratio)
DSAoverlap<-dplyr::filter(Log2Overlap, SAlog2ratio == DumasLog2ratio)
sum(WSAoverlap$locus == DSAoverlap$locus)
sum(DSAoverlap$locus == WSAoverlap$locus)
WDSAoverlap<-rbind(WSAoverlap,DSAoverlap)

187-45
(142/187)*100

# graphs
library(ggplot2)
library(reshape2)
dfm <- melt(Log2Overlap,id.vars = 1)
colnames(dfm)<-c("ID","Species","Log2ratio")
#dfm$Species<-sub("SAlog2ratio","Spartina",dfm$Species)
dfm$Species<-sub("SAlog2ratio","1",dfm$Species)
#dfm$Species<-sub("WeisLog2ratio","Weisman",dfm$Species)
dfm$Species<-sub("WeisLog2ratio","2",dfm$Species)
#dfm$Species<-sub("DumasLog2ratio","Dumas",dfm$Species)
dfm$Species<-sub("DumasLog2ratio","3",dfm$Species)

nrow(dfm)
dfm<-dfm[with(dfm, order(ID, Species)), ]

dfm$Species<-sub("1","Spartina",dfm$Species)
dfm$Species<-sub("2","Weisman",dfm$Species)
dfm$Species<-sub("3","Dumas",dfm$Species)
dfm$Species<-as.factor(dfm$Species)
#dfm$Species<-factor(dfm$Species, levels = dfm$Species[order(dfm$Species)])
dfm$Species<-factor(dfm$Species,levels(dfm$Species)[c(2,3,1)])
dfm$Log2ratio<-as.factor(dfm$Log2ratio)
dfm$Log2ratio<-factor(dfm$Log2ratio,levels(dfm$Log2ratio)[c(1,2)])

# how much overlap?
nrow(filter(Log2Overlap, SAlog2ratio == WeisLog2ratio | SAlog2ratio == DumasLog2ratio))
nrow(Log2Overlap)
42/188
100-22


write.table(dfm,"heatmap.txt",row.names = FALSE,quote = FALSE,sep="\t")

subset1<-dfm[1:273,]
subset2<-dfm[274:nrow(dfm),]
#subset3<-dfm[101:150,]
#subset4<-dfm[151:nrow(dfm),]

#ggplot(dfm,aes(x = ID, y= Log2ratio)) + geom_bar(stat="identity",aes(fill = Species),position = "dodge") +
#  labs(x="Locus",y="Log 2 Ratio") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Elaboration of heatmap (white - steelblue)

pdf("heatmap.pdf",width=30,height=9)
ggplot(dfm, aes(x = ID, y= Species)) +
  geom_tile(aes(fill = Log2ratio), color = "white") +
  scale_fill_discrete() +
  ylab("List of genes ") +
  xlab("List of patients") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Expression level")
dev.off()

# subsetted
#ggplot(subset1,aes(x = ID, y= Log2ratio)) + geom_bar(stat="identity",aes(fill = Species),position = "dodge")+ coord_flip() +
#  labs(x="Locus",y="Log 2 Ratio")
#ggplot(subset2,aes(x = ID, y= Log2ratio)) + geom_bar(stat="identity",aes(fill = Species),position = "dodge")+ coord_flip() +
#  labs(x="Locus",y="Log 2 Ratio")
#ggplot(subset3,aes(x = ID, y= Log2ratio)) + geom_bar(stat="identity",aes(fill = Species),position = "dodge")
#ggplot(subset4,aes(x = ID, y= Log2ratio)) + geom_bar(stat="identity",aes(fill = Species),position = "dodge")

# Elaboration of heatmap (white - steelblue)
pdf("heatmap_subset1.pdf",width=3,height=13)
ggplot(subset1, aes(x = ID, y= Species)) +
  coord_flip() +
  geom_tile(aes(fill = Log2ratio), color = "white") +
  guides(fill=FALSE) +
  scale_fill_discrete() +
  ylab("Experiment ") +
  xlab("Genes") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Expression level")
dev.off()

# Elaboration of heatmap (white - steelblue)

pdf("heatmap_subset2.pdf",width=3,height=13)
ggplot(subset2, aes(x = ID, y= Species)) +
  coord_flip() +
  geom_tile(aes(fill = Log2ratio), color = "white") +
  guides(fill=FALSE) +
  scale_fill_discrete() +
  ylab("Experiment") +
  xlab("Genes") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Expression level")
dev.off()






#Ov3<-na.omit(Ov2)
#subset <- t(data.frame(Ov3$SAlog2ratio, Ov3$ATlog2ratio))
#barplot(subset, legend = c("SAlog2ratio", "ATlog2ratio"), names.arg=Ov3$AT, beside=TRUE)




#2016 MDP stress annotation overlap
#Written by Mariano Alvarez
#Oct 17, 2016

setwd('/volumes/External1/R/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/Overlays/Pass5')
library(limma)
library(VennDiagram)
library(ggplot2)
library(reshape2)
### REDO OF MDP OVERLAP

library(dplyr)
library(stringr)

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
#MDP<-read.csv("MDPStress_GO.csv", header=TRUE)
MDP<-read.delim("MDPStress_GO.txt", header=TRUE, sep="\t")
MDP<-MDP[,c(1,3)]
MDP$GOStress<-as.character(MDP$GOStress)
#split<-strsplit(MDP$GOStress, " /// ", fixed = TRUE)
#split<-data.frame(matrix(unlist(split), nrow=22746, byrow=T),stringsAsFactors=FALSE)
#split<-tidyr::separate(MDP, GOStress, c(as.character(1:11)), sep = " /// ", remove = TRUE, fill = "right")
#colnames(split)[3:13]<-c("")
MDP$GOStress<-gsub("///", "",MDP$GOStress, fixed = TRUE)


stresses<-as.character(c("Herbivore","COLD","HEAT","OSMOTIC",
            "OXIDATIVE","UVB","VIRUS","FungalPathogen",
            "WOUND","SALT","GENOTOXIC","BacterialPathogen","HL","DROUGHT"))

for (x in (stresses)) {
  Stress<-filter(MDP, grepl(paste(x),GOStress))
  #Stress<-MDP %> filter(str_detect(GOStress, paste(x)))
  #Stress<- MDP %>% filter(grepl(paste(x), GOStress))
  #Stress<- MDP %>% filter(grepl("DROUGHT", GOStress))
  #Stress <- MDP[grepl(paste(x),]
  Stress$New<-1
  Stress<-Stress[,c(1,ncol(Stress))]
  colnames(Stress)[2]<-paste(x)
  MDP<-merge(MDP,Stress,by="Transcript_ID_Array_Design_", all.x=TRUE)
}

#MDP<-MDP[,-c(2)]
MDP[is.na(MDP)] <- 0
MDP$sum<-rowSums(MDP[,c(3:ncol(MDP))])
MDP<-filter(MDP, sum != 0)
MDP<-MDP[,-c(2)]

colnames(MDP)[1]<-"ATnumber"
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
StressOverlap<-merge(MDP,Sa.data,by="locus", all.x=TRUE)
# Why do I have extra genes? I guess I'll remove them?
#StressOverlap<-filter(StressOverlap, !is.na(BacterialPathogen))
# Recode Spartina NA as 0
StressOverlap[is.na(StressOverlap)] <- 0

# Look for super responders (genes with 1s for more than half of stresses, excluding SpartinaOil)
hist(StressOverlap$sum)
#StressOverlap$BaseResponseSum<-rowSums(StressOverlap[,3:11])
SuperResponders<-filter(StressOverlap, sum>4)
SpartinaResponders<-filter(SuperResponders, SpartinaOil==1)

# make SpartinaResponders table
SRtable<-merge(SpartinaResponders,ath1,by="locus", all=FALSE)
write.csv(SRtable, "SRtable.csv")

# overlap with Xenome
Xenome <- read.table("xenomic_AT_loci_and_related_spartina_probes.txt", sep = "\t", colClasses = "character")
#Xenome<-data.frame(toupper(Xenome[,1]))
colnames(Xenome)<-c(Xenome[1,1],Xenome[1,2],Xenome[1,3])
Xenome<-Xenome[-c(1),]
Xenome<-data.frame(Xenome[,2])
Xenome$AtXenome<-1
colnames(Xenome)[1]<-"locus"

XenomeOverlap<-merge(StressOverlap,Xenome, by="locus", all.x = TRUE)
#XenomeOverlap<-XenomeOverlap[,c(1:2,16,3:15)]
XenomeOverlap[is.na(XenomeOverlap)] <- 0

write.csv(XenomeOverlap,"MDPOverlaptable.csv")
nrow(filter(XenomeOverlap, AtXenome == 1 & SpartinaOil == 1 & sum > 4))






#### Final Edits

setwd('/Volumes/External1/R/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/Overlays/LastMinute')
library(limma)
library(VennDiagram)
library(dplyr)
library(progress)

Supp2 <- read.csv("SuppFile2.csv")
Supp3 <- read.csv("SuppFile3.csv")

Yellow<-filter(Supp2, moduleColor == "yellow")
Yellow<-Yellow[,21:22]
VP20<-Supp3
VP20<-VP20[1:20,]
VP20$Gene<-toupper(VP20$Gene)
colnames(VP20)[1]<-"LocusLinkID"

VP20<-distinct(VP20)
Yellow<-distinct(Yellow)
overlap<-merge(Yellow, VP20, by = "LocusLinkID", all = FALSE)


ColrNum<-data.frame(levels(Supp2$moduleColor))
ColrNum$Num<-1:nrow(ColrNum)
colnames(ColrNum)[1]<-"color"
Supp2$moduleColor<-as.character(Supp2$moduleColor)

pb<-txtProgressBar(title = "Working...", max = nrow(Supp2))
for (x in (1:nrow(Supp2))){
  getTxtProgressBar(pb)
  value<-as.character(Supp2[x,22])
  selected<-filter(ColrNum, color == paste(value))
  Supp2[x,22]<-as.character(selected[1,2])
  setTxtProgressBar(pb, x)
}
close(pb)

write.csv(Supp2, "SuppFile2_2.csv")



#################################################

### UNUSED CODE

#################################################













setwd('/Volumes/External1/R/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/Overlays/Pass5')
library(limma)
library(VennDiagram)
library(dplyr)

# First, use only the genes that appear on both the ATH1 array and the Spartina array
ath1 <- read.delim("affy_ATH1_array_elements-2010-12-20.txt", sep='\t', header=TRUE)
SpAlt1 <- read.csv("Spart_array_annotation.csv", header=TRUE)
colnames(SpAlt1)[8]<-"locus"
SpAlt1$locus<-gsub("\\..*","",SpAlt1$locus)
ConsensusGenes<-merge(ath1,SpAlt1,by="locus", all = FALSE)
write.csv(ConsensusGenes, "ATH1_SpAlt_consensus.csv",row.names = FALSE)
ConsensusGenes<-as.data.frame(ConsensusGenes$locus)
#remove duplicates
before<-nrow(as.data.frame(ConsensusGenes))
ConsensusGenes<-as.data.frame(ConsensusGenes[!duplicated(ConsensusGenes[,1]),])
print(paste("Loci removed:",before-nrow(ConsensusGenes)))
print(paste("Loci remaining:",nrow(ConsensusGenes)))
colnames(ConsensusGenes)<-"locus"


#read in Weismann data, Dumas data, and BP oil data
Weisman <- read.csv("Weisman.csv")
Dumas <- read.csv("Dumas.csv")
Dumas<-data.frame(toupper(Dumas[,1]))
colnames(Weisman)[2]<-"locus"
colnames(Dumas)<-"locus"
Weisman.origL <- nrow(Weisman)
Dumas.origL <- nrow(Dumas)
Sa.dataDUP <- read.csv("all_oil_responsive2.csv")
colnames(Sa.dataDUP)[1]<-"locus"
#remove duplicate observations in S. alterniflora
Sa.data<-Sa.dataDUP[!duplicated(Sa.dataDUP[,1]),]
Sa.origL <- nrow(Sa.data)
print(c(Weisman.origL, "original observations in Weisman"))
print(c(Dumas.origL, "original observations in Dumas"))
print(c(Sa.origL, "original observations in Spartina"))
print(c(nrow(Sa.dataDUP)-Sa.origL, "duplicate observations removed in Spartina"))

# filter out observations of each that don't appear on the consensus gene list
Weisman.consensus<-merge(Weisman,ConsensusGenes,by="locus", all = FALSE)
Dumas.consensus<-merge(Dumas,ConsensusGenes,by="locus", all = FALSE)
Sa.consensus<-merge(Sa.data,ConsensusGenes,by="locus", all = FALSE)
paste(nrow(Weisman.consensus), "consensus observations in Weisman ||", Weisman.origL-nrow(Weisman.consensus), "observations removed")
paste(nrow(Dumas.consensus), "consensus observations in Dumas ||", Dumas.origL-nrow(Dumas.consensus), "observations removed")
paste(nrow(Sa.consensus), "consensus observations in Spartina ||",Sa.origL-nrow(Sa.consensus), "observations removed")
paste("Total observations:",nrow(Weisman.consensus)+nrow(Dumas.consensus)+nrow(Sa.consensus))

### Parse out overlapping genes
#make ID column
Sa.consensus$Spartina<-1
Weisman.consensus$Weisman<-1
Dumas.consensus$Dumas<-1
#merge
Merged<-merge(Sa.consensus, Weisman.consensus, by="locus", all=TRUE)
Merged<-merge(Merged, Dumas.consensus, by="locus", all=TRUE)
Merged<-Merged[,c("Spartina","Weisman","Dumas")]
#OverlapTable<-OverlapTable[!is.na(OverlapTable$AGI),]
write.csv(Merged, "overlaptable.csv")

Merged[is.na(Merged)] <- 0
nrow(dplyr::filter(Merged,Weisman==1 & Spartina==1))
nrow(dplyr::filter(Merged,Weisman==1 & Spartina==0))
nrow(dplyr::filter(Merged,Weisman==0 & Spartina==1))
nrow(dplyr::filter(Merged,Dumas==1 & Spartina==1))
nrow(dplyr::filter(Merged,Dumas==1 & Spartina==0))
nrow(dplyr::filter(Merged,Dumas==0 & Spartina==1))
nrow(dplyr::filter(Merged,Dumas==1 & Spartina==1 & Weisman==1))
nrow(dplyr::filter(Merged,Weisman==1 & Dumas==1))

# Venn
library(VennDiagram)

draw.triple.venn(area1 = nrow(subset(Merged, Weisman == 1)), area2 = nrow(subset(Merged, Dumas == 1)),
                 area3 = nrow(subset(Merged, Spartina == 1)), n12 = nrow(subset(Merged, Weisman == 1 & Dumas == 1)),
                 n23 = nrow(subset(Merged, Dumas == 1 & Spartina == 1)), n13 = nrow(subset(Merged,Weisman == 1 & Spartina == 1)), 
                 n123 = nrow(subset(Merged, Weisman == 1 & Dumas == 1 & Spartina ==1)), 
                 category = c("Weisman et al. 2010", "Dumas et al. 2016", "Spartina"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))




#rm(list=setdiff(ls(), c("")))










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

