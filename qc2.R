##########=====  DEG using DESeq2  =====##########
library(DESeq2)
library(tximport)
library(readr)
library("RColorBrewer")
library("pheatmap")
setwd("~/epic.neb/analysis/results")

load(file="~/epic.neb/analysis/R_neb/RData/infoTable_filtered_updated.RData")

geneList <- read_csv("~/epic.neb/analysis/info_table/geneList.csv", col_types = cols(X1 = col_skip()))
geneList$gene <- paste0(geneList$gene_id,"_",geneList$gene_name)
write.csv(geneList, file="~/epic.neb/analysis/info_table/geneList.csv")


##### loading metadata
libInfo <- read_csv("~/epic.neb/analysis/info_table/library_info.csv")
metaInfo <- read_csv("~/epic.neb/analysis/info_table/EPIC_HCvB_metadata.csv",
                     col_types = cols(DMTsAtVisit = col_factor(levels = c("Copaxone", "Gilenya", "Rebif", "Aubagio", "Tecfidera", "Tysabri", "Rituximab", "none")), 
                                      DiseaseCourse = col_factor(levels = c("CIS", "PP", "RR", "RIS", "SP", "Healthy","Unknown")), 
                                      DiseaseStatus = col_factor(levels = c("MS", "Healthy", "NotMS", "Unknown")), 
                                      Last_Known_Treat_Stat = col_factor(levels = c("Treated", "TreatmentNaive", "Healthy", "NA")), 
                                      Sex = col_factor(levels = c("F", "M", "NA"))))
metaInfo_base <- read_csv("~/epic.neb/analysis/info_table/EPIC_HCvB_metadata_baseline.csv",
                     col_types = cols(DMTsAtVisit = col_factor(levels = c("Copaxone", "Gilenya", "Rebif", "Aubagio", "Tecfidera", "Tysabri", "Rituximab", "none")), 
                                      DiseaseCourse = col_factor(levels = c("CIS", "PP", "RR", "RIS", "SP", "Healthy","Unknown")), 
                                      DiseaseStatus = col_factor(levels = c("MS", "Healthy", "NotMS", "Unknown")), 
                                      Last_Known_Treat_Stat = col_factor(levels = c("Treated", "TreatmentNaive", "Healthy", "NA")), 
                                      Sex = col_factor(levels = c("F", "M", "NA"))))

sTable <- merge(libInfo, metaInfo, by.x="HCVB_ID", by.y="HCVB_ID")
sTable$Subset <- factor(sTable$Subset)
sTable$prep <- factor(sTable$prep)
sTable$Lane <- factor(sTable$Lane)
sTable$DiseaseStatus <- factor(sTable$DiseaseStatus)
sTable$Sex <- factor(sTable$Sex)
sTable$Last_Known_Treat_Stat <- factor(sTable$Last_Known_Treat_Stat)
sTable$DiseaseCourse <- factor(sTable$DiseaseCourse)

save(libInfo, metaInfo, metaInfo_base, sTable, file="~/epic.neb/analysis/R_neb/RData/infoTable.RData")


##### loading selected samples from RSEM #####
# subset metadata
sTable.sub <- subset(sTable, Subset=="CD4"); sTable.sub$Subset <- droplevels(sTable.sub$Subset)
sTable.sub <- subset(sTable, Subset=="CD8"); sTable.sub$Subset <- droplevels(sTable.sub$Subset)
sTable.sub <- subset(sTable, Subset=="CD14"); sTable.sub$Subset <- droplevels(sTable.sub$Subset)

setwd("~/epic.neb/analysis/results/qc-gene_level-2(edited label)")

# file list for rsem
setwd("~/epic.neb/counts/rsem")
files <- list.files(path=".", pattern="*.genes.results", recursive = TRUE)
tmp <- unlist(lapply(strsplit(files, "/"), "[[",2))
names(files) <- sub(".genes.results", "", tmp)
files <- files[names(files) %in% sTable.sub$Sample_ID]

row.names(sTable.sub) <- sTable.sub$Sample_ID
sTable.sub <- sTable.sub[names(files),]

# loading gene counts - subset
rsem.gene.sub <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
rm(files, tmp)


##### input and run DESeq2 for subset#####
rsem.gene.sub$length[rsem.gene.sub$length == 0] <- 1
cds.sub <- DESeqDataSetFromTximport(rsem.gene.sub, sTable.sub, ~DiseaseStatus)
cds.sub <- cds.sub[ rowSums(counts(cds.sub)) > 1, ]
cds.sub <- DESeq(cds.sub)
vst.cds <- vst(cds.sub)

cds.cd4 <- cds.sub
sTable.cd4 <- sTable.sub
save(cds.cd4, sTable.cd4, file="~/epic.neb/analysis/R_neb/RData/qc-DESeq2_object-CD4.RData")
rm(cds.cd4, sTable.cd4)



##########

### loading DESeq2 object for PCA of each cell subset
load(file="~/epic.neb/analysis/R_neb/RData/qc-DESeq2_object-CD4.RData")
load(file="~/epic.neb/analysis/R_neb/RData/qc-DESeq2_object-CD8.RData")
load(file="~/epic.neb/analysis/R_neb/RData/qc-DESeq2_object-CD14.RData")

cds.sub <- cds.cd4
sTable.sub <- sTable.cd4
vst.cds <- vst(cds.sub)

cds.sub <- cds.cd8
sTable.sub <- sTable.cd8
vst.cds <- vst(cds.sub)

cds.sub <- cds.cd14
sTable.sub <- sTable.cd14
vst.cds <- vst(cds.sub)

cds.sub <- cds.edited
sTable.sub <- sTable.edited
vst.cds <- vst(cds.sub)

rm(cds.edited, sTable.edited)


###  PCA plot using deseq2 function (plotPCA)  ###
setwd("~/epic.neb/analysis/results/qc-gene_level-2(edited label)")
plotPCA(vst.cds, intgroup = c("Subset"), ntop = 60000)
plotPCA(vst.cds, intgroup = c("DiseaseStatus"), ntop = 60000)

plotPCA(vst.cds, intgroup = c("Sex", "DiseaseStatus"), ntop = 10000)
plotPCA(vst.cds, intgroup = c("DiseaseDuration"), ntop = 1000)


# check gender (sex)
plotPCA(vst.cds, intgroup = c("Sex"), ntop = 500)
tmpPCA <- plotPCA(vst.cds, intgroup = c("Sex"), ntop = 500, returnData = TRUE)
write.csv(tmpPCA, file="PCAtable_CD14-Sex.csv")

plot(tmpPCA$PC1, tmpPCA$PC2, col = RColorBrewer::brewer.pal(3, "Dark2")[tmpPCA$Sex], pch = 26, main = "PCs plot, scaled PCA (PC1-PC2)")
text(tmpPCA$PC1, tmpPCA$PC2, col = RColorBrewer::brewer.pal(3, "Dark2")[tmpPCA$Sex], label=tmpPCA$name, cex = 0.8)
legend("bottom", levels(tmpPCA$Sex), fill= RColorBrewer::brewer.pal(3, "Dark2"), inset=.02, bty='n', cex=0.8)



###  PCA plot using plot  ###
setwd("~/epic.neb/analysis/results/qc-gene_level-2(edited label)")
pcaMatrix <- assay(vst.cds)
pcaMatrix <- pcaMatrix[rowVars(pcaMatrix) > 0.1,]
pcaNorm <- prcomp(t(pcaMatrix), scale. = FALSE)
barplot(((pcaNorm$sdev)^2 / sum((pcaNorm$sdev)^2) * 100)[1:15], main = "Variance plot, scaled PCA")

#
plot(pcaNorm$x[, c(1,2)], col = RColorBrewer::brewer.pal(3, "Dark2")[sTable.sub$Subset], pch = 20, main = "PCs plot, scaled PCA (PC1-PC2)")
legend("bottomleft", levels(sTable.sub$Subset), fill= RColorBrewer::brewer.pal(3, "Dark2"), inset=.02, bty='n', cex=0.8)

plot(pcaNorm$x[, c(1,2)], col = RColorBrewer::brewer.pal(4, "Set1")[sTable.sub$DiseaseStatus], pch = 20, main = "PCs plot, scaled PCA (PC1-PC2)")
legend("bottomleft", levels(sTable.sub$DiseaseStatus), fill= RColorBrewer::brewer.pal(4, "Set1"), inset=.00, bty='n', cex=0.8)

plot(pcaNorm$x[, c(1,2)], col = RColorBrewer::brewer.pal(3, "Dark2")[sTable.sub$Subset], pch = 26, main = "PCs plot, scaled PCA (PC1-PC2)")
text(pcaNorm$x[, c(1,2)], col = RColorBrewer::brewer.pal(3, "Dark2")[sTable.sub$Subset], label=sTable.sub$Sample_ID, cex = 0.8)
legend("bottomleft", levels(sTable.sub$Subset), fill= RColorBrewer::brewer.pal(3, "Dark2"), inset=.02, bty='n', cex=0.8)
#


### heatmap : sample-to-sample distance ###
sampleDists <- dist(t(assay(vst.cds)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- paste(sTable.sub$Sample_ID)
rownames(sampleDistMatrix) <- paste(sTable.sub$Subset)
df <- data.frame(Cell=sTable.sub$Subset, Gender=sTable.sub$Sex, Stauts=sTable.sub$DiseaseStatus)
colors <- colorRampPalette( rev(brewer.pal(9, 'YlOrRd')) )(255)
rownames(df) <- (sTable.sub$Sample_ID)

dev.off()
pdf("Heatmap_sample_distances_CD4.pdf", paper="USr", width = 0, height = 0)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, 
         col=colors, fontsize_row = 2, fontsize_col = 4, annotation_col = df, show_colnames = T, fontsize=3)  #cellwidth = 30, cellheight = 30, 
dev.off()

rm(sampleDists, sampleDistMatrix, df, colors)
#



##### input and run DESeq2#####
gene.count$length[gene.count$length == 0] <- 1
cds <- DESeqDataSetFromTximport(gene.count, sTable, ~DiseaseStatus)
cds <- cds[ rowSums(counts(cds)) > 1, ]
cds <- DESeq(cds)

# test for disease course
cds <- DESeqDataSetFromTximport(gene.count, sTable, ~DiseaseCourse)
cds <- cds[ rowSums(counts(cds)) > 1, ]
cds <- DESeq(cds)
count.table <- counts(cds, normalized=TRUE)


#####  PCA plot  #####
pcaMatrix <- assay(vst.cds)
pcaMatrix <- pcaMatrix[rowVars(pcaMatrix) > 0.1,]
pcaNorm <- prcomp(t(pcaMatrix), scale. = FALSE)
barplot(((pcaNorm$sdev)^2 / sum((pcaNorm$sdev)^2) * 100)[1:15], main = "Variance plot, scaled PCA")

#
plot(pcaNorm$x[, c(1,2)], col = RColorBrewer::brewer.pal(5, "Set1")[sTable.sub$DiseaseStatus], pch = 20, main = "PCs plot, scaled PCA (PC1-PC2)")
legend("bottomright", c("MS","Healthy"), fill= RColorBrewer::brewer.pal(3, "Set1"), inset=.02, bty='n', cex=0.8)

plot(pcaNorm$x[, c(1,2)], col = RColorBrewer::brewer.pal(3, "Set1")[sTable.sub$DiseaseStatus], pch = 26, main = "PCs plot, scaled PCA (PC1-PC2)")
text(pcaNorm$x[, c(1,2)], col = RColorBrewer::brewer.pal(3, "Set1")[sTable$DiseaseStatus], label=sTable$DiseaseStatus, cex = 0.8)
legend("topleft", c("MS","Healthy"), fill= RColorBrewer::brewer.pal(3, "Set1"), inset=.02, bty='n', cex=0.8)

plot(pcaNorm$x[, c(3,4)], col = RColorBrewer::brewer.pal(3, "Set1")[sTable.sub$DiseaseStatus], pch = 20, main = "PCs plot, scaled PCA (PC1-PC2)")
legend("bottomright", c("MS","Healthy"), fill= RColorBrewer::brewer.pal(3, "Set1"), inset=.02, bty='n', cex=0.8)

#
plot(pcaNorm$x[, c(1,2)], col = RColorBrewer::brewer.pal(7, "Dark2")[sTable.sub$DiseaseCourse], pch = 20, main = "PCs plot, scaled PCA (PC1-PC2)")
legend("topleft", c("CIS","PP","RR","RIS","SP","Healthy","Unknown"), fill= RColorBrewer::brewer.pal(7, "Dark2"), inset=.02, bty='n', cex=0.8)

plot(pcaNorm$x[, c(1,2)], col = RColorBrewer::brewer.pal(7, "Dark2")[sTable.sub$DiseaseCourse], pch = 26, main = "PCs plot, scaled PCA (PC1-PC2)")
text(pcaNorm$x[, c(1,2)], col = RColorBrewer::brewer.pal(7, "Dark2")[sTable.sub$DiseaseCourse], label=sTable$DiseaseCourse, cex = 0.8)
legend("bottomright", c("CIS","PP","RR","RIS","SP","Healthy","Unknown"), fill= RColorBrewer::brewer.pal(7, "Dark2"), inset=.02, bty='n', cex=0.8)


#
rm(pcaMatrix, pcaNorm)



##### PCA using plotPCA (using DESeq2 PCA function)  #####
#
plotPCA(vst(cds, blind=FALSE), intgroup = c("Subset", "DiseaseStatus"))

#
pcaData <- plotPCA(vsd, intgroup = c( "dex", "cell"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()



##### MDS plot  #####
# using VST data
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed()

# using Poisson distance
mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed()


###############################################################################################################################################333
###############################################################################################################################################333
###############################################################################################################################################333


##### Quick check: DEG results #####
res_deseq <- results(cds, alpha=0.05)
summary(res_deseq)
sum(res_deseq$padj < 0.05, na.rm=TRUE)
res_deseq <- as.data.frame(res_deseq)
res_deseq <- res_deseq[order(res_deseq$padj),]
write.csv(as.data.frame(res_deseq), file="CD8_result_MS-HC.csv")



##### Heatmap #####
library("pheatmap"); library("RColorBrewer")
heatmap.expr <- assay(vst(cds, blind=FALSE))

select <- row.names(subset(res_deseq, (log2FoldChange > 1.5 | log2FoldChange < -1.5) & padj<0.05)); length(select)
#select <- select[!select %in% c("ENSG00000188536.13_HBA2","ENSG00000244734.4_HBB","ENSG00000206172.8_HBA1","ENSG00000143546.9_S100A8")]
select <- select[!select %in% c("ENSG00000113070.7_HBEGF","ENSG00000169508.6_GPR183","ENSG00000123689.5_G0S2","ENSG00000167680.15_SEMA6B")]

df <- data.frame(sTable$DiseaseStatus)
row.names(df) <- sTable$filename3
colnames(df) <- c("Status")
pheatmap(heatmap.expr[select,], color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100),   # colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)
         fontsize_row = 8, fontsize_col = 8, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, 
         annotation_col=df)  # heatmap.expr[select,!colnames(heatmap.expr[select,])==c("48413b_CD4_L020")]
#
rm(select, df, heatmap.expr)



##### Plotting results #####
# quick way
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("dex"))

# custom plot using ggplot2
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex","cell"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = cell)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

ggplot(geneCounts, aes(x = dex, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()



##### single gene plot #####
library(ggplot2)
### using normalized counts
count.table <- counts(cds, normalized=TRUE)
gene <- "ENST00000368655.4_GATAD2B"
gene <- select[9]
genename <- strsplit(gene,"_")[[1]][2]

# plot
plot_table <- data.frame("expr"=as.numeric(count.table[row.names(count.table)==gene,]), "CellType"=sTable$Subset, 
                         "Status"=sTable$DiseaseStatus, "Course"=sTable$DiseaseCourse)
ggplot(plot_table, aes(Status, expr, col=Course)) + geom_boxplot(outlier.shape = 1) + 
  geom_jitter(aes(col=Course),width=0.3, shape = 20, size=1.8) + 
  ggtitle(paste0(plot_table$CellType[1],", expression for ",genename)) + xlab("") + 
  ylab("Expression level (DESeq2 normalization)") + theme_bw()

ggplot(plot_table, aes(Status, expr)) + geom_boxplot(outlier.colour = "Black", outlier.shape = 1) + 
  geom_jitter(width=0.3, shape = 20, size=2) + 
  ggtitle(paste0(plot_table$CellType[1],", expression for ",genename)) + xlab("") + 
  ylab("Expression level (DESeq2 normalization)") + theme_bw()


### using TPM
count.table <- rsem.tpm[,colnames(rsem.tpm) %in% sTable$names]
count.table <- count.table[,sTable$names]   # matching order
gene <- "ENST00000368655.4_GATAD2B"
gene <- select[9]
genename <- strsplit(gene,"_")[[1]][2]

# plot
plot_table <- data.frame("expr"=as.numeric(count.table[row.names(count.table)==gene,]), "CellType"=sTable$Subset, 
                         "Status"=sTable$DiseaseStatus, "Course"=sTable$DiseaseCourse)
ggplot(plot_table, aes(Status, expr, col=Course)) + geom_boxplot(outlier.shape = 1) + 
  geom_jitter(aes(col=Course),width=0.3, shape = 20, size=1.8) + 
  ggtitle(paste0(plot_table$CellType[1],", expression for ",genename)) + xlab("") + 
  ylab("Expression level (DESeq2 normalization)") + theme_bw()

ggplot(plot_table, aes(Status, expr)) + geom_boxplot(outlier.colour = "Black", outlier.shape = 1) + 
  geom_jitter(width=0.3, shape = 20, size=2) + 
  ggtitle(paste0(plot_table$CellType[1],", expression for ",genename)) + xlab("") + 
  ylab("Expression level (DESeq2 normalization)") + theme_bw()

#
rm(gene, genename, plot_table, count.table)



##### MA-plot #####
library("apeglm")
resultsNames(dds)

res <- lfcShrink(dds, coef="dex_trt_vs_untrt", type="apeglm")
plotMA(res, ylim = c(-5, 5))

res.noshr <- results(dds, name="dex_trt_vs_untrt")
plotMA(res.noshr, ylim = c(-5, 5))

# label individual points on the MA-plot
plotMA(res, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

# histogram of the p values
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")




##### Independent filtering #####
# The ratio of small p values for genes binned by mean normalized count
qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(resLFC1$pvalue, bins, function(p)
  mean(p < .05, na.rm = TRUE))
barplot(fractionSig, xlab = "mean normalized count",
        ylab = "fraction of small p values")




##### Annotating and exporting results #####
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

# add the gene symbol and Entrez ID to result table
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

resOrdered <- res[order(res$pvalue),]
head(resOrdered)


### Exporting results #####
resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
write.csv(resOrderedDF, file = "results.csv")

# automatically generate dynamic HTML documents
library("ReportingTools")
htmlRep <- HTMLReport(shortName="report", title="My report",
                      reportDirectory="./report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)


### Plotting fold changes in genomic space
resGR <- results(dds, name="dex_trt_vs_untrt", format="GRanges")
resGR$log2FoldChange <- res$log2FoldChange
resGR
resGR$symbol <- mapIds(org.Hs.eg.db, names(resGR), "SYMBOL", "ENSEMBL")

# plotting the GRanges and associated metadata
library("Gviz")
window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)
status <- factor(ifelse(resGRsub$padj < 0.1 & !is.na(resGRsub$padj),
                        "sig", "notsig"))

options(ucscChromosomeNames = FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name = "gene ranges", feature = status)
d <- DataTrack(resGRsub, data = "log2FoldChange", baseline = 0,
               type = "h", name = "log2 fold change", strand = "+")
plotTracks(list(g, d, a), groupAnnotation = "group",
           notsig = "grey", sig = "hotpink")




##########
rm(cds, res, res.sig)
sessionInfo()
