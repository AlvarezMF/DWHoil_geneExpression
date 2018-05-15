#Final script for WGCNA clustering, followed by PCVA on each cluster
#Written by MFA
#WGCNA code adapted from tutorial
#July 21, 2016

setwd('/Volumes/Analysis1/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/Final_ClusterPVCA_2/')
library(cluster)
#library(pnmath0)
#library(pvclust)
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)


#This is normalized expression data from the BPoil microarray
#All genes
BPdata <- read.csv("medsummclean.csv")
#BPdata<-BPdata[,-c(18:19)]

#Tutorial code to get it into a dataframe
datExpr0 = as.data.frame(t(BPdata[, -c(1:1)]));
names(datExpr0) = BPdata$ProbeName;
rownames(datExpr0) = names(BPdata)[-c(1:1)];

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

#=====================================================================================
#
#  Code chunk 4: remove bad genes
#
#=====================================================================================


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
ncol(datExpr0)

datTemp<-data.frame(t(datExpr0))
datTemp$geneName<-rownames(datTemp)
badGenes<-names(datExpr0)[!gsg$goodGenes]
datTemp<-dplyr::filter(datTemp, geneName %in% badGenes)

apply(datTemp, 1, function(x) sum(is.na(x)))
rowMeans(data.matrix(datTemp),na.rm = TRUE)
apply(datTemp, 1, function(x) sd(x,na.rm = TRUE))





#=====================================================================================
#
#  Code chunk 5: heirarchal clustering
#
#=====================================================================================


sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

#=====================================================================================
#
#  Code chunk 7: trait data
#
#=====================================================================================

datExpr = datExpr0

traitData = read.csv("BPoilTraits.csv");
dim(traitData)
names(traitData)
allTraits = traitData
# Form a data frame analogous to expression data that will hold the traits. 
femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$SamplePool);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage();


#=====================================================================================
#
#  Code chunk 8: cluster with colors
#
#=====================================================================================


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Hierarchical clustering of sample pools")

#=====================================================================================
#
# Intermediate: plot PCA
#
#=====================================================================================

library(vegan)

datExpr2<-data.frame(rownames(datExpr),datExpr)
colnames(datExpr2)[1]<-"Sample"
Pheno<-read.csv("Pheno.csv")
Pheno$Sample<-paste(Pheno$Sample,"_Median",sep="")
Pheno<-dplyr::filter(Pheno,Sample %in% datExpr2$Sample)
Pheno2<-data.table::fread("/Volumes/Analysis1/2016SpAlt_epiGBS/depth_10/UPheno.txt",data.table=FALSE)
Pheno$Site<-paste(Pheno$Site,Pheno$Treatment,sep="_")
Pheno$Site<-gsub("GIA_oil","GIO1",Pheno$Site)
Pheno$Site<-gsub("GIB_oil","GIO2",Pheno$Site)
Pheno$Site<-gsub("BSL_oil","MSO",Pheno$Site)
Pheno$Site<-gsub("BSL_no_oil","MSN",Pheno$Site)
Pheno$Site<-gsub("GIA_no_oil","GIN1",Pheno$Site)
Pheno$Site<-gsub("GIB_no_oil","GIN2",Pheno$Site)


datExpr2<-dist(datExpr)

datExpr2<-metaMDS(datExpr2,k=5,trymax = 50)

datExpr2<-data.frame(datExpr2$points,Pheno)

library(ggplot2)
library(ggthemes)
pdf("MDSplot.pdf",height = 5,width = 6)
print(
ggplot(datExpr2,aes(x=MDS1,y=MDS2,colour=Treatment,shape=Site)) + geom_point() + labs(x="MDS1",y="MDS2") + 
  theme_minimal() + scale_color_fivethirtyeight() 
)
dev.off()

library(GGally)

ggparcoord(datExpr2, scale = "uniminmax",columns = 1:5, groupColumn = 7) + 
  labs(x="MDS dimension",y="") +
  theme_minimal() + scale_color_fivethirtyeight() 


#=====================================================================================
#
#  Code chunk 9: save as Rdata
#
#=====================================================================================


save(datExpr, datTraits, file = "AllGenes_InitialClustering_Cleanup.RData")



#=====================================================================================
#
#  Part II
#
#=====================================================================================

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

setwd('/Volumes/Analysis1/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/Final_ClusterPVCA_2/')
library(cluster)
#library(pnmath0)
library(pvclust)
library(WGCNA)

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
#workingDir = ".";
#setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
#enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "AllGenes_InitialClustering_Cleanup.RData");
#The variable lnames contains the names of loaded variables.
lnames


#=====================================================================================
#
#  Code chunk 2: Building/plotting the network
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, corFnc = "bicor", blockSize = 5000)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#=====================================================================================
#
#  Code chunk 3: auto construction of network
#
#=====================================================================================

#maxBlockSize adjusted for 16gb RAM

net = blockwiseModules(data.matrix(datExpr), power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       checkMissingData = FALSE,
                       replaceMissingAdjacencies = TRUE,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, #pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       #saveTOMFileBase = "femaleMouseTOM", 
                       verbose = 3,
                       maxBlockSize=5000,
                       corType = "bicor")


#saveRDS(net, "net.rds")
net<-readRDS("net.rds")

#=====================================================================================
#
#  Code chunk 4: Plot cluster dendrogram
#
#=====================================================================================


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#=====================================================================================
#
#  Code chunk 5: 
#
#=====================================================================================


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
#save(MEs, moduleLabels, moduleColors, geneTree,file = "BPoil_autocluster.RData")
load("BPoil_autocluster.RData")

#=====================================================================================
#=====================================================================================
#     Phenotypic correlations
#=====================================================================================
#=====================================================================================







#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
saveRDS(MEs, "MEs.rds")
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#=====================================================================================
#
#  Code chunk 3 - not needed
#
#=====================================================================================


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Define variable OIL containing the weight column of datTrait
weight = as.data.frame(datTraits$Treatment);
names(weight) = "Treatment"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");



#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


names(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================

#Changed to turq
names(datExpr)[moduleColors=="turquoise"]


#=====================================================================================
#
#  Code chunk 8: loading annotations
#
#=====================================================================================


annot = read.csv(file = "GeneAnnotation.csv");
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$ProbeName)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.
#Jk not in our data! We have a few.


#=====================================================================================
#
#  Code chunk 9: Ordering data and creating a giant table
#
#=====================================================================================


# Create the starting data frame
geneInfo0 = data.frame(ProbeName = probes,
                       geneSymbol = annot$TargetID[probes2annot],
                       LocusLinkID = annot$AT[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight OR WHATEVER TRAIT YOU CHOSE
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Treatment));
geneInfo = geneInfo0[geneOrder, ]


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================

#change filename as needed
write.csv(geneInfo, file = "geneInfo_AllGenes.csv")



#=====================================================================================
#=====================================================================================
#=====================================================================================
#        PVCA
#=====================================================================================
#=====================================================================================
#=====================================================================================

#impute missing values
BPdata <- read.csv("medsummclean.csv")
library(impute)
ProbeNames <- BPdata$ProbeName
dfclean<-as.matrix(BPdata[,2:19])
if (any(is.na(dfclean)))
  df_imp <- impute.knn(dfclean[])$data
df_imp<-data.frame(df_imp)
df_imp$ProbeName<-ProbeNames


#merge BPdata and mod memberships

MergedDataset <- merge(df_imp, geneInfo, "ProbeName")
MergedDataset <- MergedDataset[,1:22]
MergedDataset$moduleColor<-as.numeric(as.factor(MergedDataset$moduleColor))


write.csv(MergedDataset, "MergedDataSet.csv", row.names=FALSE)



nrow(dplyr::filter(MergedDataset, moduleColor==3))








###### Can start here



setwd('~/documents/code/R/2015BPoil/Deepwater-Horizon-Oil-Spill-analyses/Final_ClusterPVCA/')





MergedDataset<- read.csv("MergedDataSet.csv", header=TRUE)

#create color names string
ModuleNames<-unique(MergedDataset$moduleColor)

#load libraries
library(Biobase)
library(pvca)

pData = read.table("BPoilTraits.csv",row.names=1, header=TRUE, sep=",")
sapply(pData,class)

metadata <- data.frame(labelDescription=c("Oil exposure, 1", "Field plot", "Which State"),
                       rownames=c("Treatment", "Population", "State"))
phenoData = new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
phenoData

getPVCA <- function(x) {
  #Get each color from the input
  x <- MergedDataset[ which(MergedDataset$moduleColor==x), ]
  x <- x[,2:19]
  exprs <- as.matrix(x)
  CADset <- ExpressionSet(assayData=exprs, phenoData=phenoData)
  pct_threshold = 0.7
  batch.factors <- c("Population", "State", "Treatment")
  pvcaObj <- pvcaBatchAssess(CADset, batch.factors, pct_threshold)
  return(pvcaObj)
}

ModuleNames

plotPVCA <- function(x) {
  thiscluster<-x
  pvcaObj<-getPVCA(thiscluster)
  pdf(paste(x, "PVCA.pdf", sep=""),width=10,height=10)
  bp <- barplot(pvcaObj$dat,
                ylab = "Weighted average proportion variance", ylim= c(0,1.1),
                col = c("blue"), las=2, main=paste(x,"PVCA estimation bar chart", sep=" "))
  
  axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 0.5, las=2)
  values = pvcaObj$dat
  new_values = round(values , 2)
  text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.8)
  dev.off()
}



length(ModuleNames)

lapply(ModuleNames, plotPVCA)








#=====================================================================================
#=====================================================================================
#=====================================================================================
#       Generalized linear models for eigengenes
#=====================================================================================
#=====================================================================================
#=====================================================================================

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
#enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "AllGenes_InitialClustering_Cleanup.RData");
#The variable lnames contains the names of loaded variables.
lnames


Eigens<-readRDS("MEs.rds")
Names<-colnames(Eigens)
Eigens$SamplePool<-rownames(datTraits)
datTraits$SamplePool<-rownames(datTraits)
Eigens<-merge(Eigens, datTraits, by="SamplePool")

## save point
saveRDS(Eigens, "eigens.rds")
## load point
Eigens<-readRDS("eigens.rds")


library(lme4)
library(nlme)
library(car)
library(MASS)
library(lmerTest)
library(dplyr)
library(data.table)

#hist(Eigens$MEdarkgreen)

list<-colnames(Eigens[,2:(ncol(Eigens)-3)])

system("rm PsTreat.txt")
getGLMs <- function(x) {
  #system(paste("rm ","EigenFit_", x, ".txt",sep=""))
  formula = "~ Treatment + State + (1|Treatment:Population)"
  #ME.fit=lmer(Eigens[,c("x")] ~ Treatment*State + Treatment:(1|Population), data=Eigens)
  #ME.fit=lmer(formula(paste(x, formula)), data=Eigens)
  #formula2 = "~ Treatment"
  ME.fit=lmer(formula(paste(x, formula)), data=Eigens)
  #ME.fit2=lm(formula(paste(x, formula2)), data=Eigens)
  A=data.frame(x,summary(ME.fit)$coefficients[2,5])
  #cat("\n\n", file = paste("EigenFit_", x, ".txt", sep=""), append = TRUE)
  # export anova test output
  #cat(c(paste(x),"Summary\n"), file = paste("EigenFit_", x, ".txt", sep=""), append = TRUE)
  #capture.output(A, file = paste("EigenFit_", x, ".txt", sep=""), append = TRUE)
  #B=car::Anova(ME.fit)
  #cat("\n\n", file = paste("EigenFit_", x, ".txt", sep=""), append = TRUE)
  #cat(c(paste(x),"ANOVA\n"), file = paste("EigenFit_", x, ".txt", sep=""), append = TRUE)
  #capture.output(B, file = paste("EigenFit_", x, ".txt", sep=""), append = TRUE)
  cat("\n\n", file = "PsTreat.txt", append = TRUE)
  fwrite(A, file = "PsTreat.txt", append = TRUE,col.names=FALSE,row.names = FALSE,quote = FALSE,sep="\t")
  #cat("\n\n", file = "PsPop.txt", append = TRUE)
  #cat(paste(x, A$coefficients[3,4], sep="\t"), file = "PsPop.txt", append = TRUE)
  #C<-anova(ME.fit,ME.fit2)
}

# start cluster
parallelCluster <- parallel::makeCluster(detectCores()-1, type = "FORK")
print(parallelCluster)

parLapply(parallelCluster, list, getGLMs)

# Shutdown cluster neatly
if(!is.null(parallelCluster)) {
  parallel::stopCluster(parallelCluster)
  parallelCluster <- c()
}

gc()


ForFDR.treat<-read.delim("PsTreat.txt", sep = "\t", header=FALSE)
ForFDR.treat$FDR<-p.adjust(ForFDR.treat$V2, method = "holm")
#ForFDR.pop<-read.delim("PsPop.txt", sep = "\t", header=FALSE)
#ForFDR.pop$FDR<-p.adjust(ForFDR.pop$V2, method = "holm")


filter(ForFDR.treat, FDR <=0.05)

#Sig fits for treatment: yellow**,brown*,blue*,pink**,
SigModTreat<-c("brown")
#Sig fits for state: green*,yellow*,white*,red*
#SigModState<-c("green","yellow","white","red")

MergedDataset<- read.csv("MergedDataSet.csv", header=TRUE)
x<-dplyr::filter(MergedDataset, moduleColor== SigModTreat)
write.csv(x, "brown.csv")
nrow(x)
Ats<-data.frame(x$LocusLinkID)
fwrite(Ats,file="brown.ats.txt",sep="\t",quote=FALSE,row.names = FALSE)

VPlant<-fread("/Users/marianoalvarez/Dropbox (Personal)/Alvarez/Microarray_MS/Final microarray ms/BPOilMolEcol/SuppFile3.csv",data.table = FALSE)

colnames(VPlant)[1]<-"LocusLinkID"
x2<-merge(VPlant,x,all.x=TRUE,by="LocusLinkID")
attach(x2)
x2<-x2[order(-Connections),]
detach(x2)
rownames(x2)<-1:nrow(x2)

Sequences<-fread("/Users/marianoalvarez/Library/Mobile Documents/com~apple~CloudDocs/Research/Projects/2012BPOilMicroarray/sequence_array.csv")
colnames(x)[20]<-"TargetID"

x.seq<-merge(x,Sequences,by="TargetID",all.x=TRUE)
x.seq<-x.seq[,c(1,24)]
library(seqinr)

for (i in 1:nrow(x.seq)){
  x.seq[i,2]<-c2s(translate(s2c(x.seq[i,2])))
}

fwrite(x.seq,file="protein_sequences.txt",sep="\t",col.names = FALSE,row.names = FALSE,quote = FALSE)

library(dplyr)
MergedCount<-1
for (i in SigModTreat) {
  x<-dplyr::filter(MergedDataset, moduleColor== i)
  write.csv(x, paste(SigModTreat[MergedCount], "Merged.csv", sep=""))
  assign(i, x)
  MergedCount<-MergedCount+1
}

read.csv("geneInfo_AllGenes.csv", header=TRUE)

MergedDataset$moduleNumber<-as.numeric(as.factor(MergedDataset$moduleColor))

#=====================================================================================
#
#  Code chunk 5: making mod membership plots
#
#=====================================================================================

###CHANGE MODULE HERE AS NEEDED

#module = "purple"
#column = match(module, modNames);
#moduleGenes = moduleColors==module;

#sizeGrWindow(7, 7);
#par(mfrow = c(1,1));
#verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                   abs(geneTraitSignificance[moduleGenes, 1]),
#                   xlab = paste("Module Membership in", module, "module"),
#                   ylab = "Gene significance for treatment",
#                   main = paste("Module membership vs. gene significance\n"),
#                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

plotModMem <- function(x) {
  module = x
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  pdf(paste(x, "ModMem.pdf", sep=""),width=10,height=10)
  sizeGrWindow(7, 7);
  plot<-verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for treatment",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
}

lapply(SigModTreat, plotModMem)
lapply(SigModState, plotModMem)




#### make a 4 panel plot
sizeGrWindow(7, 7);
par(mfrow = c(2,2));

module = "pink"
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for treatment",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


#=====================================================================================
#=====================================================================================
#=====================================================================================
#        Network visualization using WGCNA functions
#=====================================================================================
#=====================================================================================
#=====================================================================================

#COMPUTATIONALLY INTENSIVE - erase all other variable besides datExpr before attempting
rm(list= ls()[!(ls() %in% c('datExpr','geneTree', 'moduleColors', 'datTraits'))])

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

library(pnmath0)
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
#dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
#dissTOM = ""
#load("femaleMouseTOM-block.1.RData")
#TOM<-data.matrix(TOM)
#plotTOM<-TOM
#rm(TOM)
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
#plotTOM = dissTOM^7;
#save(plotTOM, file = "plotTom.RData")
#load("plotTom.RData")
load("femaleMouseTOM-block.1.RData")
plotTOM<-TOM
rm(TOM)
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

##### EIGENGENE NETWORK

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate treatment from the  traits
treatment = as.data.frame(datTraits$Treatment);
names(treatment) = "treatment"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, treatment))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)






#### CYTOSCAPE EXPORT

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

#Sig fits for treatment: lightyellow*, tan**, magenta**, purple**, brown**, red**, yellow**, 
#Sig fits for state: turquoise*, green***, 

load("femaleMouseTOM-block.1.RData")
lnames = load(file = "AllGenes_InitialClustering_Cleanup.RData");
cluster = load(file = "BPoil_autocluster.RData");
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Read in the annotation file
annot = read.csv(file = "GeneAnnotation.csv");
# Select modules
modules = c("yellow");
# Select module probes

probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);


