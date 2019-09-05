#####  Results with Age at exam
##########  WGCNA
source("https://bioconductor.org/biocLite.R")
biocLite("impute")
biocLite("preprocessCore")
install.packages("WGCNA")

##########
library(WGCNA)
library(DESeq2)
library(flashClust)
library(readr)
library(GSEABase)
library(genefilter)
library(org.Hs.eg.db)

options(stringsAsFactors = FALSE)     # The following setting is important, do not omit.
allowWGCNAThreads(nThreads = 4)      # disableWGCNAThreads()

#ALLOW_WGCNA_THREADS=12


###########################################################################################3
##### method from https://www.bioconductor.org/packages/release/bioc/vignettes/CVE/inst/doc/WGCNA_from_TCGA_RNAseq.html
##### and https://wikis.utexas.edu/display/bioiteam/Clustering+using+WGCNA
#load("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-AllCell_MS-all and HC_withAge.RData")
load("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-AllCells_MSvsHC_withAge_updated.RData")
setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge-updated")
rm(dds.cd14, dds.cd4, dds.cd8, geneCounts)
#setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge")


## CD4
#load("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-AllCell_MS-all and HC_withAge.RData")
setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge-updated/CD4_signed")
cds <- cds.cd4
vst <- vst.cd4
res_deseq <- res3.cd4
celltype <- "CD4"
sTable.sub <- subset(sTable.cd4, Subset=="CD4" & (Last_Known_Treat_Stat=="TreatmentNaive" | Last_Known_Treat_Stat=="Healthy"))
datTraits <- data.frame(MSvsHC=as.factor(sTable.sub$Last_Known_Treat_Stat), DiseaseCourse=sTable.sub$DiseaseCourse,
                        Sex=sTable.sub$Sex, DiseaseDuration=sTable.sub$DiseaseDuration, Age=sTable.sub$AgeAtExam,
                        EDSS=sTable.sub$EDSS, MSSS=sTable.sub$MSSS)
row.names(datTraits) <- sTable.sub$Sample_ID
# disease course
datTraits$DiseaseCourse <- as.character(datTraits$DiseaseCourse)
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "Healthy"] <- NA
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "Unknown"] <- NA
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "RIS"] <- NA
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "CIS"] <- 1
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "RR"] <- 2
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "PP"] <- 3
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "SP"] <- 3
datTraits$DiseaseCourse <- as.numeric(datTraits$DiseaseCourse)
#levels(datTraits$DiseaseCourse) <- 1:length(levels(datTraits$DiseaseCourse))
#datTraits$DiseaseCourse <- as.numeric(datTraits$DiseaseCourse)
# age
datTraits[datTraits$MSvsHC == "Healthy",]$Age <- NA
# sex
datTraits[datTraits$MSvsHC == "Healthy",]$Sex <- NA
levels(datTraits$Sex) <- 1:length(levels(datTraits$Sex))
datTraits$Sex <- as.numeric(datTraits$Sex)
#
datTraits$MSvsHC <- as.character(datTraits$MSvsHC)
datTraits$MSvsHC[datTraits$MSvsHC == "TreatmentNaive"] <- 1
#datTraits$MSvsHC[datTraits$MSvsHC == "Treated"] <- 2
datTraits$MSvsHC[datTraits$MSvsHC == "Healthy"] <- 2
#levels(datTraits$MSvsHC) <- 1:length(levels(datTraits$MSvsHC))
datTraits$MSvsHC <- as.numeric(datTraits$MSvsHC)
head(datTraits)



## CD8
#load("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-AllCell_MS-all and HC_withAge.RData")
setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge-updated/CD8_signed")
cds <- cds.cd8
vst <- vst.cd8
res_deseq <- res3.cd8
celltype <- "CD8"
sTable.sub <- subset(sTable.cd8, Subset=="CD8" & (Last_Known_Treat_Stat=="TreatmentNaive" | Last_Known_Treat_Stat=="Healthy"))
datTraits <- data.frame(MSvsHC=as.factor(sTable.sub$Last_Known_Treat_Stat), DiseaseCourse=sTable.sub$DiseaseCourse,
                        Sex=sTable.sub$Sex, DiseaseDuration=sTable.sub$DiseaseDuration, Age=sTable.sub$AgeAtExam,
                        EDSS=sTable.sub$EDSS, MSSS=sTable.sub$MSSS)
row.names(datTraits) <- sTable.sub$Sample_ID
# disease course
datTraits$DiseaseCourse <- as.character(datTraits$DiseaseCourse)
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "Healthy"] <- NA
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "Unknown"] <- NA
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "RIS"] <- NA
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "CIS"] <- 1
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "RR"] <- 2
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "PP"] <- 3
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "SP"] <- 3
datTraits$DiseaseCourse <- as.numeric(datTraits$DiseaseCourse)
#levels(datTraits$DiseaseCourse) <- 1:length(levels(datTraits$DiseaseCourse))
#datTraits$DiseaseCourse <- as.numeric(datTraits$DiseaseCourse)
# age
datTraits[datTraits$MSvsHC == "Healthy",]$Age <- NA
# sex
datTraits[datTraits$MSvsHC == "Healthy",]$Sex <- NA
levels(datTraits$Sex) <- 1:length(levels(datTraits$Sex))
datTraits$Sex <- as.numeric(datTraits$Sex)
#
datTraits$MSvsHC <- as.character(datTraits$MSvsHC)
datTraits$MSvsHC[datTraits$MSvsHC == "TreatmentNaive"] <- 1
#datTraits$MSvsHC[datTraits$MSvsHC == "Treated"] <- 2
datTraits$MSvsHC[datTraits$MSvsHC == "Healthy"] <- 2
#levels(datTraits$MSvsHC) <- 1:length(levels(datTraits$MSvsHC))
datTraits$MSvsHC <- as.numeric(datTraits$MSvsHC)
head(datTraits)



## CD14
#load("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-AllCell_MS-all and HC_withAge.RData")
setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge-updated/CD14_signed")
cds <- cds.cd14
vst <- vst.cd14
res_deseq <- res3.cd14
celltype <- "CD14"
sTable.sub <- subset(sTable.cd14, Subset=="CD14" & (Last_Known_Treat_Stat=="TreatmentNaive" | Last_Known_Treat_Stat=="Healthy"))
datTraits <- data.frame(MSvsHC=as.factor(sTable.sub$Last_Known_Treat_Stat), DiseaseCourse=sTable.sub$DiseaseCourse,
                        Sex=sTable.sub$Sex, DiseaseDuration=sTable.sub$DiseaseDuration, Age=sTable.sub$AgeAtExam,
                        EDSS=sTable.sub$EDSS, MSSS=sTable.sub$MSSS)
row.names(datTraits) <- sTable.sub$Sample_ID
# disease course
datTraits$DiseaseCourse <- as.character(datTraits$DiseaseCourse)
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "Healthy"] <- NA
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "Unknown"] <- NA
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "RIS"] <- NA
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "CIS"] <- 1
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "RR"] <- 2
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "PP"] <- 3
datTraits$DiseaseCourse[datTraits$DiseaseCourse == "SP"] <- 3
datTraits$DiseaseCourse <- as.numeric(datTraits$DiseaseCourse)
#levels(datTraits$DiseaseCourse) <- 1:length(levels(datTraits$DiseaseCourse))
#datTraits$DiseaseCourse <- as.numeric(datTraits$DiseaseCourse)
# age
datTraits[datTraits$MSvsHC == "Healthy",]$Age <- NA
# sex
datTraits[datTraits$MSvsHC == "Healthy",]$Sex <- NA
levels(datTraits$Sex) <- 1:length(levels(datTraits$Sex))
datTraits$Sex <- as.numeric(datTraits$Sex)
#
datTraits$MSvsHC <- as.character(datTraits$MSvsHC)
datTraits$MSvsHC[datTraits$MSvsHC == "TreatmentNaive"] <- 1
#datTraits$MSvsHC[datTraits$MSvsHC == "Treated"] <- 2
datTraits$MSvsHC[datTraits$MSvsHC == "Healthy"] <- 2
#levels(datTraits$MSvsHC) <- 1:length(levels(datTraits$MSvsHC))
datTraits$MSvsHC <- as.numeric(datTraits$MSvsHC)
head(datTraits)




##########  prepare matrix  ########## ########## ########## ########## ########## ########## ########## ##########
#setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC/CD4")
#setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC/CD8")
#setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC/CD14")

setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge-updated/CD4_signed/withALL")
setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge-updated/CD4_signed")

setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge-updated/CD8_signed/withALL")
setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge-updated/CD8_signed")

setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge-updated/CD14_signed/withALL")
setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge-updated/CD14_signed")


### filter by variance
#exp_matrix <- counts(cds, normalized=TRUE)
#exp_matrix <- exp_matrix[,colnames(exp_matrix) %in% sTable.sub$Sample_ID]
#exp_matrix_filt <- exp_matrix[,sTable.sub$Sample_ID]     # select samples
#exp_matrix_filt <- varFilter(exp_matrix_filt, var.func=IQR, var.cutoff=0.6, filterByQuantile=TRUE)
#nrow(exp_matrix_filt)          # 30% genes filtered by variance (CD4:11,553, CD8:11,787, CD14: 11,074)
#  40% genes filtered for non-coding RNAs (because significant ncRNAs are dropped when I filtered by 30%) 
#  40% genes filtered: CD4, ; CD8, ; CD14, 
#tmp <- data.frame(row.names(exp_matrix_filt))


### filtering
exp_matrix <- counts(cds, normalized=TRUE)
nrow(exp_matrix)
exp_matrix_filt <- exp_matrix[subset(res_deseq, padj < 0.2)$gene, sTable.sub$Sample_ID]     # select samples   # padj < 0.2 & 
nrow(exp_matrix_filt)


#####
vst_matrix <- assay(vst)
vst_matrix_filt <- vst_matrix[,sTable.sub$Sample_ID]     # select samples
vst_matrix_filt <- vst_matrix_filt[row.names(vst_matrix_filt) %in% row.names(exp_matrix_filt),]
WGCNA_matrix = as.data.frame(t(vst_matrix_filt))
#WGCNA_matrix = as.data.frame(t(vst_matrix))
#WGCNA_matrix = as.data.frame(t(vst_matrix[order(apply(vst_matrix,1,mad), decreasing = T)[1:10000],]))   # median absolute devision (MAD) was used as a robust measure of variability.
head(WGCNA_matrix)[1:5]




##########  network construction and module detection - manual  ##########
##### select power #####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 11, to=25, by=1))
# Call the network topology analysis function
sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5)
## Plot the results:
sizeGrWindow(9, 5); par(mfrow = c(1,2)); cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.9,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# setting threshold
#softPower = 18     # CD4=18; CD8=; CD14= 14     # padj < 0.1
softPower = 6     # CD4=7; CD8=12; CD14= 6       # padj < 0.2 


##### from WGCNA vignette & Wiki of UT Austin #####
#calculation a adjacency, "correlation" matrix

#adjacency = adjacency(WGCNA_matrix, power = softPower, type = "unsigned")
adjacency = adjacency(WGCNA_matrix, power = softPower, type = "signed")


#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
#TOM = TOMsimilarity(adjacency, TOMType="unsigned") # specify network type

dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 10         # 25 previous size
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
table(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors); length(unique(dynamicColors))
write.csv(data.frame(table(dynamicColors)), paste0("MS-HC_",celltype,"_moduleColor_numbers.csv"))   # number of genes in each color

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
# title: MS-HC_CD4_plot-dendrogram_colors
png(filename=paste0("MS-HC_",celltype,"_plot-dendrogram_colors.png"), width=800, height=600)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()


### Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList = moduleEigengenes(WGCNA_matrix, colors = dynamicColors, softPower = softPower)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
png(filename=paste0("MS-HC_",celltype,"_plot-Clustering_module_eigengenes.png"), width=1000, height=800)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()


#
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
table(dynamicColors); length(unique(dynamicColors))
# choose a height cut (0.25), corresponding to correlation (0.75), to merge
#MEDissThres = 0.0

MEDissThres = 0.4     # CD4, 0.3; CD8 0.4; CD14 0.4
abline(h=MEDissThres, col = "red")





### If do NOT merge                ###############
moduleColors <- dynamicColors 
table(moduleColors); length(unique(dynamicColors))
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1





### If DO Merge                     ###############
# Call an automatic merging function
merge = mergeCloseModules(WGCNA_matrix, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs
sizeGrWindow(12, 9)

png(filename=paste0("MS-HC_",celltype,"_plot-dendrogram_colors_merged.png"), width=800, height=600)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()
png(filename=paste0("MS-HC_",celltype,"_plot-dendrogram_colors_merged-only.png"), width=800, height=600)
plotDendroAndColors(geneTree, cbind(mergedColors), c("Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()
# Rename to moduleColors
moduleColors = mergedColors
table(moduleColors)
write.csv(data.frame(table(moduleColors)), paste0("MS-HC_",celltype,"_moduleColor_numbers-merged.csv"))

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
# Save module colors and labels for use in subsequent parts
#save(MEs, moduleLabels, moduleColors, geneTree, file = "02-networkConstruction-stepByStep.RData")







##########  Relating modules to external information (phenotype) and identifying important genes  ##########
### Relating modules to external clinical traits (phenotype) - Quantifying module{trait associations
# Define numbers of genes and samples
nGenes = ncol(WGCNA_matrix); nGenes
nSamples = nrow(WGCNA_matrix); nSamples

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(WGCNA_matrix, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 10, 3, 2))
# Display the correlation values within a heatmap plot
# title: MS-HC_CD4_heatmap_Module-trait associations
#png("MS-HC_CD4_heatmap_Module-trait associations.png", width=1200, height=800)
pdf(paste0("MS-HC_",celltype,"_Module_heatmap-trait associations.pdf"), paper="USr", width = 10, height = 9)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


table(moduleColors)


### Gene relationship to trait and important modules: Gene Significance and Module Membership
### quantify associations of individual genes with our trait of interest (weight) by defining Gene Significance GS as (the absolute value of) the correlation between the gene and the trait.
### quantify the similarity of all genes on the array to every module.
# Define variable weight containing the weight column of datTrait

# for Treatment Status:  c("green","turquoise","blue","grey60")  
# for Disease Course:  c("pink","grey60")
# for EDSS: c(black","lightyellow","red")

phenotype = as.data.frame(datTraits$MSvsHC)
names(phenotype) = "MSvsHC"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(WGCNA_matrix, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(WGCNA_matrix, phenotype, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(phenotype), sep="");
names(GSPvalue) = paste("p.GS.", names(phenotype), sep="");




### intramodular connectivity
# whole network connectivity (denoted by kTotal) or the intramodular connectivity (kWithin)
IM = intramodularConnectivity(adjacency, moduleColors, scaleByMax = FALSE) 
genes <- names(WGCNA_matrix)    # anotation

#names(WGCNA_matrix)
#names(WGCNA_matrix)[moduleColors==module]
length(names(WGCNA_matrix)[moduleColors==intModules[1]])

#library(readr)
#annot <- read_csv("~/epic.neb/analysis/info_table/geneList.csv", col_types = cols(X1 = col_skip()))
#dim(annot)
#names(annot)
#probes2annot <- unlist(lapply(strsplit(genes, "_"), "[[",1))
# The following is the number or probes without annotation:
#sum(is.na(probes2annot))        # Should return 0.


## create a data frame holding the following information for all probes: probe ID, gene symbol, Locus Link ID (Entrez code),
## module color, gene significance for selected phenotype "TreatmentStatus", and module membership and p-values in all modules.
# Create the starting data frame
geneInfo0 = data.frame(ID=genes,
                       geneID = unlist(lapply(strsplit(genes, "_"), "[[",1)),
                       geneSymbol = unlist(lapply(strsplit(genes, "_"), "[[",2)),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
write.csv(geneInfo0, file = paste0("MS-HC_",celltype,"_geneInfo.csv"))
# Order modules by their significance for selected phenotype "TreatmentStatus"
modOrder = order(-abs(cor(MEs, phenotype, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.MSvsHC))      # change column name for phenotype
geneInfo = geneInfo0[geneOrder, ]
head(geneInfo)[1:7]

geneInfo$entrez <- mapIds(org.Hs.eg.db,
                          keys=substr(geneInfo$geneID, 1,15),
                          column="ENTREZID",
                          keytype="ENSEMBL",
                          multiVals="first")
geneInfo$alias <- mapIds(org.Hs.eg.db,
                         keys=substr(geneInfo$geneID, 1,15),
                         column="ALIAS",
                         keytype="ENSEMBL",
                         multiVals="CharacterList")
geneInfo$alias <- do.call(rbind, lapply(geneInfo$alias,paste0,collapse="/"))

geneInfo = merge(geneInfo, res_deseq[!duplicated(res_deseq),], by.x="ID", by.y="gene",  sort = FALSE, all.x = TRUE)
write.csv(geneInfo, file = paste0("MS-HC_",celltype,"_geneInfo_moduleInfo.csv"))




#####
##########  Network visualization using WGCNA functions  ##########

### Visualizing the gene network
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
#dissTOM = 1-TOMsimilarityFromExpr(WGCNA_matrix, power = softPower)
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
# title: MS-HC_CD4_plot-Network_heatmap_all
#pdf(paste0("MS-HC_",celltype,"_plot-Network_heatmap_all.pdf"), paper="USr", width = 0, height = 0)
png(paste0("MS-HC_",celltype,"_plot-Network_heatmap_all.png"), width=2400, height=2000)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

nSelect = 500   # select number of genes
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
# title: MS-HC_CD4_plot-Network_heatmap_selected
png(paste0("MS-HC_",celltype,"_plot-Network_heatmap_selected.png"), width=2400, height=2000)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()


### Visualizing the network of eigengenes
# Recalculate module eigengenes
MEs = moduleEigengenes(WGCNA_matrix, moduleColors)$eigengenes
# Isolate trait (phenotype) from the clinical traits
#phenotype = as.data.frame(datTraits$MSvsHC);
#names(phenotype) = "MSvsHC"
# Add the trait (phenotype) to existing module eigengenes
MET = orderMEs(cbind(MEs, phenotype))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
# title: MS-HC_CD4_plot-dendrogram_heatmap
png(paste0("MS-HC_",celltype,"_plot-dendrogram_heatmap.png"), width=1200, height=1600)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)
dev.off()


## To split the dendrogram and heatmap plots
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
# title: MS-HC_CD4_plot-dendrogram
png(paste0("MS-HC_",celltype,"_plot-dendrogram.png"), width=800, height=600)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
# title: MS-HC_CD4_plot-heatmap
png(paste0("MS-HC_",celltype,"_plot-heatmap.png"), width=800, height=600)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), 
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()




##### Intramodular analysis: identifying genes with high GS and MM

# Select significant modules  :  royalblue
#intModules <- c("darkred","pink","royalblue")   # CD8
#intModules <- c("darkturquoise","brown","turquoise","yellow")   # CD4
#intModules <- c("cyan","turquoise")    # CD14


#intModules <- c("darkgreen","grey","blue","lightcyan", "lightgreen","green","darkgrey","yellow")      # CD4 unsigned with ncRNA
#intModules <- c("turquoise","darkgrey","lightyellow","blue", "brown")      # CD4 signed with ncRNA
intModules <- unique(moduleColors)


for (i in 1:length(intModules)){
  modules = intModules[i]
  column = match(modules, modNames)
  moduleGenes = moduleColors==modules
  #  sizeGrWindow(7, 7)
  #  par(mfrow = c(1,1))
  png(filename=paste0("MS-HC_",celltype,"_modulePlot_",names(phenotype),"_",modules,".png"), width=700, height=600)
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", modules, "module"), ylab = "Gene significance for treatment status",
                     main = paste("Module membership vs. gene significance\n"), cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = modules)
  dev.off()
}



#####
##########  Interfacing network analysis - functional annotation and gene ontology  ##########
### Output gene lists for use with online software and services
# Read in the probe annotation
#annot = read_csv("~/epic.neb/analysis/info_table/geneList.csv", col_types = cols(X1 = col_skip()))
# intramodular connectivity
IM = intramodularConnectivity(adjacency, moduleColors, scaleByMax = FALSE) 

# $ Choose interesting modules
genes = names(WGCNA_matrix)
intModules

for (module in intModules){
  # Select module probes
  modGenes = genes[moduleColors==module]
  # Get their entrez ID codes
  modIDs = geneInfo[geneInfo$ID %in% modGenes,c(1:4,(ncol(geneInfo)-1):(ncol(geneInfo)))]
  # Write them into a file
  fileName = paste("MS-HC_",celltype,"_moduleIDs-", module, ".csv", sep="");
  write.table(as.data.frame(modIDs), file = fileName, sep=",",
              row.names = FALSE, col.names = TRUE)
}

# additional information
for (module in intModules){
  # Select module probes
  modGenes = genes[moduleColors==module]
  # Get their entrez ID codes
  modIDs = geneInfo[geneInfo$ID %in% modGenes, c(1:4,(ncol(geneInfo)-1):(ncol(geneInfo)))]
  # add p-value and log2fc
  modIDs = data.frame(modIDs, IM[row.names(IM) %in% modGenes,],
                      geneInfo[geneInfo$ID %in% modGenes, c(paste("MM.", module, sep=""),paste("p.MM.", module, sep=""))])
  modIDs = merge(modIDs, res_deseq, by.x="ID", by.y="gene", all.x=TRUE)
  # Write them into a file
  fileName = paste("MS-HC_",celltype,"_moduleIDs-", module, "+.csv", sep="");
  write.table(as.data.frame(modIDs), file = fileName, sep=",",
              row.names = FALSE, col.names = TRUE)
}

# As background in the enrichment analysis, we will use all probes in the analysis.
#fileName = paste("LocusLinkIDs-all.txt", sep="");
#write.table(as.data.frame(allLLIDs), file = fileName,
#            row.names = FALSE, col.names = FALSE)






### Enrichment analysis directly within R
# required packages: GO.db, AnnotationDBI, org.Hs.eg.db
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)

GOenr = GOenrichmentAnalysis(moduleColors, geneInfo$entrez, organism = "human", nBestP = 10)
tab = GOenr$bestPTerms[[4]]$enrichment
module.num = length(unique(tab$module))
for(m in 1:module.num){
  for(i in 1:10){
    entrez = GOenr$bestPTerms[[4]]$forModule[[m]][[i]]$geneCodes
    genesymbol = geneInfo[match(entrez, geneInfo$entrez), "geneSymbol"]
    tab$geneSymbol[(m-1)*10 + i] = paste0(genesymbol,collapse="|")
  }
}
tab <- tab[,c(1:5,7,6,15,8:14)]
write.table(tab, file = paste0("MS-HC_",celltype,"_GOEnrichmentTable.csv"), sep = ",", quote = TRUE, row.names = FALSE)



#####  find hub genes
annot <- read_csv("~/epic.neb/analysis/info_table/geneList.csv", col_types = cols(X1 = col_skip()))
annot$idName <- paste0(annot$gene_id,"_",annot$gene_name)

hubGenes <- chooseTopHubInEachModule(WGCNA_matrix, moduleColors, type="signed") 
hubGenes <- as.data.frame(hubGenes); hubGenes$colorMod <- row.names(hubGenes)
hubGenes <- merge(hubGenes, annot, by.x="hubGenes", by.y="idName")
write.csv(hubGenes, paste0("MS-HC_",celltype,"_topHubGenes_signed.csv"))

hubGenes <- chooseTopHubInEachModule(WGCNA_matrix, moduleColors, type="unsigned") 
hubGenes <- as.data.frame(hubGenes); hubGenes$colorMod <- row.names(hubGenes)
hubGenes <- merge(hubGenes, annot, by.x="hubGenes", by.y="idName")
write.csv(hubGenes, paste0("MS-HC_",celltype,"_topHubGenes_unsigned.csv"))






##########  Exporting a gene network to external visualization software  ##########
### Exporting to Cytoscape
# Recalculate topological overlap if needed
#TOM = TOMsimilarityFromExpr(WGCNA_matrix, power = softPower)

# Read in the annotation file
#annot = read.csv(file = "GeneAnnotation.csv")
#annot <- data.frame(ID=names(WGCNA_matrix), geneID=unlist(lapply(strsplit(genes, "_"), "[[",1)), 
#                    geneName=unlist(lapply(strsplit(genes, "_"), "[[",2)))    # anotation


setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge-updated/CD4_signed/withALL")

setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge-updated/CD4_signed/cytoscape")
setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge-updated/CD8_signed/cytoscape")
setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge-updated/CD14_signed/cytoscape")

#setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge/CD4_signed/cytoscape")
#setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge/CD8_signed/cytoscape")
#setwd("~/epic.neb/analysis/results/coexp_wgcna/MS_untreat-HC_ResultWithAge/CD14_signed/cytoscape")


## Select modules
intModules  # CD4: "turquoise", "black"; CD8: "5 turquoise", "1 black", "4 magenta"; CD14: "1 turquoise" "2 blue" "4 green"
modules <- intModules[c(1,3)]; 

modules <- intModules[4]; 
modules

##
annot <- read_csv("~/epic.neb/analysis/info_table/geneList.csv", col_types = cols(X1 = col_skip()))
annot$idName <- paste0(annot$gene_id,"_",annot$gene_name)

probes = names(WGCNA_matrix)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]; length(modProbes)
modGenes = annot$gene_name[match(modProbes, annot$idName)]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

## Export the network into edge and node list files Cytoscape can read
threshold = 0.3
threshold = 0.25
threshold = 0.2
threshold = 0.15
threshold = 0.1
threshold = 0.09
threshold = 0.05
threshold = 0.01


cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("Cyto-",paste(modules, collapse="-"),"-edges","(",threshold, ").txt", sep=""),
                               nodeFile = paste("Cyto-",paste(modules, collapse="-"),"-nodes","(",threshold, ").txt", sep=""),
                               weighted = TRUE,
                               threshold = threshold,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
nrow(cyt$nodeData); nrow(cyt$edgeData)
#tmp <- cyt$nodeData
#tmp
#


## select top intramodular hub genes
modInfo <- read_csv(paste("../MS-HC_",celltype,"_moduleIDs-", modules, "+.csv", sep=""))

probes = names(WGCNA_matrix)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]; length(modProbes)

#nTop = length(modProbes)*0.1; nTop
nTop = 30
kIN = softConnectivity(WGCNA_matrix[, modProbes], power=softPower) 
selectHubs = (rank (-kIN) <= nTop)
modProbes[selectHubs]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
#threshold = 0.15
cyt = exportNetworkToCytoscape(modTOM[selectHubs,selectHubs],
                               edgeFile = paste(celltype,"_Cyto-",paste(modules, collapse="-"),"-edges","-30HubGene_",threshold,".txt", sep=""),
                               nodeFile = paste(celltype,"_Cyto-",paste(modules, collapse="-"),"-nodes","-30HubGene_",threshold,".txt", sep=""),
                               weighted = TRUE,
                               threshold = threshold,
                               nodeNames = modProbes[selectHubs],
                               altNodeNames = annot$gene_name[match(modProbes[selectHubs], annot$idName)],
                               nodeAttr = moduleColors[inModule][1:nTop])
write.table(modInfo[modInfo$ID %in% modProbes[selectHubs],], 
            file=paste(celltype,"_Cyto-",paste(modules, collapse="-"),"-geneInfo","-30HubGene_",threshold,".txt", sep=""), quote = F, row.names = F, sep="\t")
nrow(cyt$nodeData); nrow(cyt$edgeData)
#tmp <- cyt$nodeData
#tmp
#




## select genes significant in DEG analysis (padj < 0.1)
modInfo <- read_csv(paste("../MS-HC_",celltype,"_moduleIDs-", modules, "+.csv", sep=""))

modInfo_sig <- as.data.frame(subset(modInfo, padj < 0.05)); nrow(modInfo_sig)
write.table(modInfo_sig, file=paste(celltype,"_Cyto-",paste(modules, collapse="-"),"-geneInfo", "-sigGenes.txt", sep=""), quote = F, row.names = F, sep="\t")

probes = names(WGCNA_matrix)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]; length(modProbes)

selectHubs <- is.finite(match(modProbes,modInfo_sig$ID))
length(selectHubs[selectHubs==TRUE])
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

threshold = 0.02
cyt = exportNetworkToCytoscape(modTOM[selectHubs,selectHubs],
                               edgeFile = paste(celltype,"_Cyto-",paste(modules, collapse="-"),"-edges","-sigGenes_",threshold,".txt", sep=""),
                               nodeFile = paste(celltype,"_Cyto-",paste(modules, collapse="-"),"-nodes","-sigGenes_",threshold,".txt", sep=""),
                               weighted = TRUE,
                               threshold = threshold,
                               nodeNames = modProbes[selectHubs],
                               altNodeNames = annot$gene_name[match(modProbes[selectHubs], annot$idName)],
                               nodeAttr = moduleColors[inModule][1:length(selectHubs[selectHubs==TRUE])])
nrow(cyt$nodeData); nrow(cyt$edgeData)
#
tmp <- cyt$nodeData
#################


