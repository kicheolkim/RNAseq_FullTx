source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
install.packages(c("rgl", "car"))

setwd("~/epic.neb/analysis")

##### prepare sample information (metadata) #####
library(readr)
libInfo <- read_csv("~/epic.neb/analysis/info_table/library_info_original.csv")
metaInfo <- read_csv("~/epic.neb/analysis/info_table/EPIC_HCvB_metadata.csv",
                     col_types = cols(DMTsAtVisit = col_factor(levels = c("Copaxone", "Gilenya", "Rebif", "Aubagio", "Tecfidera", "Tysabri", "Rituximab", "none")), 
                                      DiseaseCourse = col_factor(levels = c("CIS", "PP", "RR", "RIS", "SP", "Healthy","Unknown")), 
                                      DiseaseStatus = col_factor(levels = c("MS", "Healthy", "NotMS", "Unknown")), 
                                      Last_Known_Treat_Stat = col_factor(levels = c("Treated", "TreatmentNaive", "Healthy", "NA")), 
                                      Sex = col_factor(levels = c("F", "M", "NA"))))
sTable <- merge(libInfo, metaInfo, by.x="HCVB_ID", by.y="HCVB_ID")

### remove samples in analysis
sTable <- sTable[sTable$Sample_ID != "78217c-CD14",]
sTable <- sTable[sTable$Sample_ID != "64615a-CD8",]



#####################################################################################################################333
#####
setwd("~/epic.neb/analysis/results/qc-gene_level-1")

### prepare sampleTable for deseq2 input
sTable <- sTable[sTable$Sample_ID %in% colnames(rsem.gene$counts),]
rownames(sTable) <- sTable$Sample_ID
sTable <- sTable[colnames(rsem.gene$counts),]
sTable$Subset <- as.factor(sTable$Subset)
sTable$DiseaseStatus <- as.factor(sTable$DiseaseStatus)


### distrigution of TPM values
tropical= c('darkorange', 'dodgerblue','hotpink', 'limegreen', 'yellow')
palette(tropical)

logTPM <- log2(rsem.gene$abundance + 1)
pdf("boxplot_all.pdf", paper="USr", width = 0, height = 0)
boxplot(logTPM, col=as.numeric(sTable$Subset), ylab='log2(TPM+1)')
dev.off()
boxplot(logTPM, col=as.numeric(sTable$Subset),las=2,ylab='log2(TPM+1)')



########=====  QC - using DESeq2  =====########
### input DESeq2
library(DESeq2)
setwd("~/epic.neb/analysis/results/qc-gene_level-1")

# RSEM - gene
gene.count <- rsem.gene
list <- gsub("-","_",rownames(sTable))
colnames(gene.count$abundance) <- list
colnames(gene.count$counts) <- list
colnames(gene.count$length) <- list
row.names(sTable) <- list
rm(list)

gene.count$length[gene.count$length == 0] <- 1
cds.gene <- DESeqDataSetFromTximport(gene.count, sTable, ~1)
cds.gene <- cds.gene[ rowSums(counts(cds.gene)) > 1, ]
cds.gene <- DESeq(cds.gene)
rm(gene.count)


##### PCA plot
##### loading selected samples from RSEM #####
setwd("~/epic.neb/analysis/results/qc-gene_level-2(edited label)")
# subset metadata
sTable.sub <- subset(sTable, Subset=="CD4"); sTable.sub$Subset <- droplevels(sTable.sub$Subset)
sTable.sub <- subset(sTable, Subset=="CD8"); sTable.sub$Subset <- droplevels(sTable.sub$Subset)
sTable.sub <- subset(sTable, Subset=="CD14"); sTable.sub$Subset <- droplevels(sTable.sub$Subset)

## PCA using TPM, all genes
pcaMatrix <- rsem.gene$abundance; pcaMatrix <- pcaMatrix[ rowSums(pcaMatrix) > 1, ]
pcaMatrix <- log2(pcaMatrix+1)
pcaNorm <- prcomp(t(pcaMatrix), scale. = FALSE)
#barplot(((pcaNorm$sdev)^2 / sum((pcaNorm$sdev)^2) * 100)[1:15], main = "Variance plot, scaled PCA")
#write.csv(pcaNorm$x[, 1:2], "TPM_pc-table_pc1-2.csv")
save(sTable, pcaMatrix, pcaNorm, file="~/epic.neb/analysis/R_neb/RData/TPM_PCA_data.RData")
rm(pcaMatrix, pcaNorm)

## PCA using vst transformed value from DESeq2, all genes
pcaMatrix <- assay(vst.cds)
pcaNorm <- prcomp(t(pcaMatrix), scale. = FALSE)
#barplot(((pcaNorm$sdev)^2 / sum((pcaNorm$sdev)^2) * 100)[1:15], main = "Variance plot, scaled PCA")
#write.csv(pcaNorm$x[, 1:2], "vst_pc-table_pc1-2.csv")
save(sTable, pcaMatrix, pcaNorm, file="~/epic.neb/analysis/R_neb/RData/vst_PCA_data.RData")
rm(pcaMatrix, pcaNorm)
##


###
library("car")
library("rgl")
PC1 <- pcaNorm$x[,1]
PC2 <- pcaNorm$x[,2]
PC3 <- pcaNorm$x[,3]
group <- sTable$Subset
group[group=="CD8minus"] <- "CD4"
group[group=="CD4_CD14"] <- "CD4"
group <- factor(group)
scatter3d(x=PC2, y=PC1, z=PC3, 
          point.col = RColorBrewer::brewer.pal(5, "Dark2")[sTable$Subset], sphere.size = 1.2,
          group = group, axis.col = c("black", "black", "black"),
          surface=FALSE, grid = FALSE, ellipsoid = TRUE, surface.col = RColorBrewer::brewer.pal(3, "Dark2")[group])
rgl.snapshot(filename = "PCA123_plot.png")

rm(PC1, PC2, PC3, group)


##
plot(pcaNorm$x[, 1:2], col = RColorBrewer::brewer.pal(5, "Dark2")[sTable$Subset], pch = 20, main = "PCs plot (PC1-PC2)")
legend("top", c("CD14","CD4","CD4+CD14","CD8","CD8-"), fill= RColorBrewer::brewer.pal(5, "Dark2"), inset=.0, bty='n', cex=0.7)

plot(pcaNorm$x[, 1:2], col = RColorBrewer::brewer.pal(5, "Dark2")[sTable$Subset], pch = 26, main = "PCs plot (PC1-PC2)")
text(pcaNorm$x[, 1:2], col = RColorBrewer::brewer.pal(5, "Dark2")[sTable$Subset], label=sTable$Sample_ID, cex = 0.8)
legend("top", c("CD14","CD4","CD4+CD14","CD8","CD8-"), fill= RColorBrewer::brewer.pal(5, "Dark2"), inset=.00, bty='n', cex=0.7)

#
plot(pcaNorm$x[, 1:2], col = RColorBrewer::brewer.pal(3, "Set1")[sTable$prep], pch = 20, main = "PCs plot (PC1-PC2)")
legend("top", c("Kicheol","Ryan"), fill= RColorBrewer::brewer.pal(3, "Set1"), inset=.02, bty='n', cex=0.8)


##
plot(pcaNorm$x[, 3:4], col = RColorBrewer::brewer.pal(5, "Dark2")[sTable$Subset], pch = 20, main = "PCs plot (PC3-PC4)")
legend("bottomright", c("CD14","CD4","CD4+CD14","CD8","CD8-"), fill= RColorBrewer::brewer.pal(5, "Dark2"), inset=.02, bty='n', cex=0.8)

plot(pcaNorm$x[, 3:4], col = RColorBrewer::brewer.pal(5, "Dark2")[sTable$Subset], pch = 26, main = "PCs plot (PC3-PC4)")
text(pcaNorm$x[, 3:4], col = RColorBrewer::brewer.pal(5, "Dark2")[sTable$Subset], label=sTable$Sample_ID, cex = 0.8)
legend("bottomright", c("CD14","CD4","CD4+CD14","CD8","CD8-"), fill= RColorBrewer::brewer.pal(5, "Dark2"), inset=.02, bty='n', cex=0.8)

##
plot(pcaNorm$x[, c(1,3)], col = RColorBrewer::brewer.pal(5, "Dark2")[sTable$Subset], pch = 20, main = "PCs plot (PC1-PC3)")
legend("top", c("CD14","CD4","CD4+CD14","CD8","CD8-"), fill= RColorBrewer::brewer.pal(5, "Dark2"), inset=.00, bty='n', cex=0.7)

plot(pcaNorm$x[, c(1,3)], col = RColorBrewer::brewer.pal(5, "Dark2")[sTable$Subset], pch = 26, main = "PCs plot (PC1-PC3)")
text(pcaNorm$x[, c(1,3)], col = RColorBrewer::brewer.pal(5, "Dark2")[sTable$Subset], label=sTable$Sample_ID, cex = 0.5)
legend("top", c("CD14","CD4","CD4+CD14","CD8","CD8-"), fill= RColorBrewer::brewer.pal(5, "Dark2"), inset=.00, bty='n', cex=0.7)



### heatmap : sample-to-sample distance 
library("RColorBrewer")
library("pheatmap")
sampleDists <- dist(t(assay(vst.cds)))

sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- paste(sTable$Sample_ID)
rownames(sampleDistMatrix) <- paste(sTable$Subset)
df <- data.frame(Cell=sTable$Subset, Gender=sTable$Sex, Stauts=sTable$DiseaseStatus)
colors <- colorRampPalette( rev(brewer.pal(9, 'YlOrRd')) )(255)
rownames(df) <- (sTable$Sample_ID)

pdf("Heatmap_sample_distances_all.pdf", paper="USr", width = 0, height = 0)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, 
         col=colors, fontsize_row = 0.8, fontsize_col = 1, annotation_col = df, show_colnames = T, fontsize=3)  #cellwidth = 30, cellheight = 30, 
dev.off()

rm(sampleDists, sampleDistMatrix, df, colors)
