# -----
# title: "DEG analysis (DESeq2) - MS-all patients (treated+untreated) vs Healthy controls"
# -----
##########=====  DEG analysis using DESeq2  =====##########
library(DESeq2)
library(tximport)
library(readr)
library(RColorBrewer)
library(pheatmap)


##### load all dataset
load("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-AllCells_MSvsHC_withAge_updated.RData")



##########  CD4 vs CD8 vs CD14
##########  Venn diagram for subsets - comparison among different cell types  ##########
library("VennDiagram")
library("RColorBrewer")
library("gplots")

# genes (l2fc > 1 or < -1) and (adjusted p-value < 0.05)
venn <- list("CD4"=subset(res.cd4, (log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.05)$gene, 
             "CD8"=subset(res.cd8, (log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.05)$gene, 
             "CD14"=subset(res.cd14, (log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.05)$gene)
grid.newpage()
venn.plot <- venn.diagram(venn, NULL, resolution=600, fill=RColorBrewer::brewer.pal(3, "Dark2"), col="transparent", alpha=c(0.5,0.5,0.5), 
                          main.cex = 1.5, sub.cex=3, cat.cex=1.5,
                          category=c("CD4","CD8","CD14"), main="")
grid.draw(venn.plot)



##### comparison between cell types: treated vs hc - for paper figure
venn <- list("CD4"=subset(res2.cd4, (baseMean > 3) & (padj < 0.05))$gene, 
             "CD8"=subset(res2.cd8, (baseMean > 3) & (padj < 0.05))$gene, 
             "CD14"=subset(res2.cd14, (baseMean > 3) & (padj < 0.05))$gene)
grid.newpage()
venn.plot <- venn.diagram(venn, NULL, resolution=600, fill=RColorBrewer::brewer.pal(3, "Dark2"), col="transparent", alpha=c(0.5,0.5,0.5), 
                          main.cex = 1.5, sub.cex=3, cat.cex=1.5,
                          category=c("CD4","CD8","CD14"), main="")
grid.draw(venn.plot)


##### comparison between cell types: untreated vs hc - for paper figure
venn <- list("CD4"=subset(res3.cd4, (baseMean > 3) & (padj < 0.05))$gene, 
             "CD8"=subset(res3.cd8, (baseMean > 3) & (padj < 0.05))$gene, 
             "CD14"=subset(res3.cd14, (baseMean > 3) & (padj < 0.05))$gene)
grid.newpage()
venn.plot <- venn.diagram(venn, NULL, resolution=600, fill=RColorBrewer::brewer.pal(3, "Dark2"), col="transparent", alpha=c(0.5,0.5,0.5), 
                          main.cex = 1.5, sub.cex=3, cat.cex=1.5,
                          category=c("CD4","CD8","CD14"), main="")
grid.draw(venn.plot)



##### comparison between each comparison conditions - for paper figure
# CD4
venn <- list("CD4:Untreated vs HC"=subset(res3.cd4, (baseMean > 3) & (padj < 0.05))$gene, 
             "CD4:Treated vs HC"=subset(res2.cd4, (baseMean > 3) & (padj < 0.05))$gene)
venn.plot <- venn.diagram(venn, NULL, resolution=600, fill=c(RColorBrewer::brewer.pal(4, "Dark2")[1],RColorBrewer::brewer.pal(4, "Dark2")[2]), 
                          col="transparent", alpha=c(0.5,0.5), 
                          sub.cex=3, cat.cex=1.5,
                          category=c("CD4:Untreated vs HC","CD4:Treated vs HC"))
grid.newpage()
grid.draw(venn.plot)

# CD8
venn <- list("CD8:Untreated vs HC"=subset(res3.cd8, (baseMean > 3) & (padj < 0.05))$gene, 
             "CD8:Treated vs HC"=subset(res2.cd8, (baseMean > 3) & (padj < 0.05))$gene)
venn.plot <- venn.diagram(venn, NULL, resolution=600, fill=c(RColorBrewer::brewer.pal(4, "Dark2")[1],RColorBrewer::brewer.pal(4, "Dark2")[2]), 
                          col="transparent", alpha=c(0.5,0.5), 
                          sub.cex=3, cat.cex=1.5,
                          category=c("CD8:Untreated vs HC","CD8:Treated vs HC"))
grid.newpage()
grid.draw(venn.plot)

# CD14
venn <- list("CD14:Untreated vs HC"=subset(res3.cd14, (baseMean > 3) & (padj < 0.05))$gene, 
             "CD14:Treated vs HC"=subset(res2.cd14, (baseMean > 3) & (padj < 0.05))$gene)
venn.plot <- venn.diagram(venn, NULL, resolution=600, fill=c(RColorBrewer::brewer.pal(4, "Dark2")[1],RColorBrewer::brewer.pal(4, "Dark2")[2]), 
                          col="transparent", alpha=c(0.5,0.5), 
                          sub.cex=3, cat.cex=1.5,
                          category=c("CD14:Untreated vs HC","CD14:Treated vs HC"))
grid.newpage()
grid.draw(venn.plot)



### extract list from venn diagram
tmp.list <- venn(venn, show.plot=FALSE)
str(tmp.list)
tmp.inters <- attr(tmp.list,"intersections")
lapply(tmp.inters, head, n=50)
#ven.list1 <- data.frame("CD4"=tmp.inters$CD14)
#ven.list2 <- data.frame("CD8"=tmp.inters$CD8)
tmp.inters$`CD4:CD8:CD14`
tmp.inters$`CD4:CD8`

#
do.call(rbind, strsplit(as.character(ven.list1[,1]), '_'))
#




###############################   boxplot by Sex
#####  single gene boxplot  ##### selected genes
target <- "ENSG00000164867.10_NOS3"

##### CD4
cds <- cds.cd4
res <- res3.cd4   #res <- res.cd4
cell <- "CD4"

### select gene for boxplot     # ENSG00000184557.4_SOCS3, ENSG00000126861.4_OMG, ENSG00000163930.9_BAP1
sig_gene_list <- res
#sig_gene_list <- subset(res, padj < 0.1); nrow(sig_gene_list)

gene <- sig_gene_list[sig_gene_list$gene==target,]$gene
genename <- strsplit(gene,"_")[[1]][2]; genename

# plot
geneCounts <- plotCounts(cds, gene = gene, intgroup = c("Sex"),
                         returnData = TRUE)
#geneCounts <- subset(geneCounts, geneCounts$Last_Known_Treat_Stat != "Treated")
colnames(geneCounts) <- c("count","Status")
geneCounts$Status <- as.character(geneCounts$Status)
#geneCounts[geneCounts$Status == "TreatmentNaive","Status"] <- "MS-untreated"
#geneCounts[geneCounts$Status == "Treated","Status"] <- "MS-treated"
#geneCounts$Status <- factor(geneCounts$Status, levels=c("MS-untreated","MS-treated","Healthy"))
ggplot(geneCounts, aes(x = Status, y = count, color = Status)) + geom_boxplot() + 
  scale_y_log10() +  geom_beeswarm(cex = 3) + 
  ggtitle(paste0(genename," in ", cell)) + xlab("") + 
  ylab("log10(DESeq2 normalized count)") + theme_bw()
# save plot
ggsave(paste0(cell,"_singlePlot_",genename,".png"), plot = last_plot(), device="png",
       scale = 1, width = 4.66, height = 5.34, units = "in")






############################################################
############################################################
############################################################
##### loop for single gene boxplot
sig_gene_list <- subset(res, padj < 0.1)$gene; length(sig_gene_list)
#
for (i in 1:length(sig_gene_list)){
  gene <- sig_gene_list[i]
  genename <- strsplit(gene,"_")[[1]][2]; genename
  geneCounts <- plotCounts(cds, gene = gene, intgroup = c("Last_Known_Treat_Stat"),
                           returnData = TRUE)
  geneCounts <- subset(geneCounts, geneCounts$Last_Known_Treat_Stat != "Treated")
  ggplot(geneCounts, aes(x = Last_Known_Treat_Stat, y = count, color = Last_Known_Treat_Stat)) + geom_boxplot() + 
    scale_y_log10() +  geom_beeswarm(cex = 1.5) + 
    ggtitle(paste0(genename, " (adj. p-val = ",round(res[res$gene==gene,"padj"],digit=9),")")) + xlab("") + 
    ylab("log10(DESeq2 normalized count)") + theme_bw()
  
  ggsave(paste0(cell,"_singlePlot_",i,"_",genename,".png"), plot = last_plot(), device="png",
         scale = 1, width = 4.66, height = 5.34, units = "in")
}
#



setwd("~/epic.neb/analysis/results/gene-level/MS_all-HC/Trt-Untreat-HC_withAge/single_plots_2_PRE")
PREgene_T <- read_csv("~/epic.neb/analysis/results/gene-level/MS_all-HC/Trt-Untreat-HC_withAge/Gene_PRE_network-Tcell.csv")
PREgene_M <- read_csv("~/epic.neb/analysis/results/gene-level/MS_all-HC/Trt-Untreat-HC_withAge/Gene_PRE_network-monocyte.csv")

sig_gene_list <- res[res$gene_name %in% PREgene_T$gene,]$gene   # CD4 and CD8
sig_gene_list <- res[res$gene_name %in% PREgene_M$gene,]$gene   # CD14

for (i in 1:length(sig_gene_list)){
  gene <- sig_gene_list[i]
  genename <- strsplit(gene,"_")[[1]][2]; genename
  geneCounts <- plotCounts(cds, gene = gene, intgroup = c("Last_Known_Treat_Stat"),
                           returnData = TRUE)
  geneCounts <- subset(geneCounts, geneCounts$Last_Known_Treat_Stat != "Treated")
  ggplot(geneCounts, aes(x = Last_Known_Treat_Stat, y = count, color = Last_Known_Treat_Stat)) + geom_boxplot() + 
    scale_y_log10() +  geom_beeswarm(cex = 1.5) + 
    ggtitle(paste0(genename, " (adj. p-val = ",round(res[res$gene==gene,"padj"],digit=9),")")) + xlab("") + 
    ylab("log10(DESeq2 normalized count)") + theme_bw()
  
  ggsave(paste0(cell,"_PRE-T_",i,"_",genename,".png"), plot = last_plot(), device="png",
         scale = 1, width = 4.66, height = 5.34, units = "in")
}

write.csv(res[res$gene_name %in% PREgene_T$gene,], "MS-HC_CD4_result-Naive.vs.HC_PRE-Tcell.csv")
write.csv(res[res$gene_name %in% PREgene_T$gene,], "MS-HC_CD8_result-Naive.vs.HC_PRE-Tcell.csv")
write.csv(res[res$gene_name %in% PREgene_M$gene,], "MS-HC_CD14_result-Naive.vs.HC_PRE-Tcell.csv")
#####





########################################################################################################################
########## custom plot using ggplot2 - disease courses with healthy controls ##########
library(ggplot2)
library(ggbeeswarm)
# CD4
cds <- cds.cd4
res <- res.cd4
cell <- "CD4"
# CD8
cds <- cds.cd8
res <- res.cd8
cell <- "CD8"
# CD14
cds <- cds.cd14
res <- res.cd14
cell <- "CD14"
#


###   Disease Course
gene <- "ENSG00000200913.1_SNORD46"
genename <- strsplit(gene,"_")[[1]][2]; genename
### using normalized count
geneCounts <- plotCounts(cds, gene = gene, intgroup = c("DiseaseCourse"),
                         returnData = TRUE)
#geneCounts$count <- geneCounts$count + 1
ggplot(geneCounts, aes(x = DiseaseCourse, y = count, color = DiseaseCourse)) + geom_boxplot() + 
  scale_y_log10() +  geom_beeswarm(cex = 1) + 
  ggtitle(paste0(gene, ", ", genename)) + xlab("") + 
  ylab("log10(DESeq2 normalized count)") + theme_bw()

#
ggsave(paste0(cell,"_singlePlot_",genename,"_withHC.png"), plot = last_plot(), device="png",
       scale = 1, width = 4.66, height = 5.34, units = "in")
#
rm(gene, genename)



### using TPM
#rsem.tpm <- rsem.gene$abundance
#count.table <- rsem.tpm[,colnames(rsem.tpm) %in% sTable.sub$Sample_ID]
#count.table <- count.table[,sTable.sub$Sample_ID]   # matching order
#plot_table <- data.frame("count"=as.numeric(count.table[row.names(count.table)==gene,]), 
#                         "TreatStatus"=sTable.sub$Last_Known_Treat_Stat)
#plot_table$count <- plot_table$count +2
#ggplot(plot_table, aes(x = TreatStatus, y = count, color = TreatStatus)) + geom_boxplot() + 
#  scale_y_log10() +   geom_beeswarm(cex = 0.1) + 
#  ggtitle(paste0(genename, " (adj. p-value = ",round(res[res$gene==gene,"padj"],digit=6),")")) + xlab("") + 
#  ylab("log10(TPM)") + theme_bw()

#
rm(select, gene, genename, plot_table, count.table)



### additional check
geneCounts <- plotCounts(cds, gene = gene, intgroup = c("DiseaseStatus","ActiveFlareAtVisit"),
                         returnData = TRUE)
# plot
ggplot(geneCounts, aes(x = DiseaseStatus, y = count, color = ActiveFlareAtVisit)) + geom_boxplot() + 
  scale_y_log10() +  geom_beeswarm(cex = 2) + 
  ggtitle(paste0(genename, " (adj. p-value = ",round(res[res$gene==gene,"padj"],digit=6),")")) + xlab("") + 
  ylab("Expression level (DESeq2 normalized count)") + theme_bw()
#

## with disease courses
ggplot(plot_table, aes(Status, expr, col=Course)) + geom_boxplot(outlier.shape = 1) + 
  geom_jitter(aes(col=Course),width=0.3, shape = 20, size=1.8) + 
  ggtitle(paste0(plot_table$CellType[1],", expression for ",genename)) + xlab("") + 
  ylab("Expression level (DESeq2 normalization)") + theme_bw()
#




##########  Annotating and exporting results  ##########
# CD4
res <- res.cd4
# CD8
res <- res.cd8
# CD14
res <- res.cd14
###
library("AnnotationDbi")
library("org.Hs.eg.db")
#columns(org.Hs.eg.db)

### add the gene symbol and Entrez ID to result table
keys <- as.character(lapply(strsplit(res$gene, "_"), "[[",1))
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=substr(keys,1,15),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
resOrdered <- res[order(res$padj),]
head(resOrdered)


##### Exporting results #####
resOrderedDF <- as.data.frame(resOrdered)
resOrderedDF <- subset(resOrderedDF, padj<0.05)
#write.csv(resOrderedDF, file = "results.csv")

### automatically generate dynamic HTML documents
library("ReportingTools")

htmlRep <- HTMLReport(shortName="CD14_report", title="CD14, Treat+Untreat vs Healthy",
                      reportDirectory="./report_CD14")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)
#
rm(keys, resOrdered, resOrderedDF, htmlRep, url)
#



##### heatmap : sample-to-sample distance #####
vst <- vst.cd4
cds <- cds.cd4
res <- res.cd4
sTable.sub <- sTable.cd4
#
vst <- vst.cd8
cds <- cds.cd8
res <- res.cd8
sTable.sub <- sTable.cd8
#
vst <- vst.cd14
cds <- cds.cd14
res <- res.cd14
sTable.sub <- sTable.cd14

##
heatmap.expr <- assay(normTransform(cds))
heatmap.expr <- heatmap.expr - rowMeans(heatmap.expr)

res <- res[order(res$padj),]
select <- subset(res, (log2FoldChange > 0.5 | log2FoldChange < -0.5) & padj < 0.05)$gene; length(select)
selected.expr <- heatmap.expr[select,]
#select2 <- subset(sTable.sub, Last_Known_Treat_Stat!="Healthy")$Sample_ID
#selected.expr <- heatmap.expr[select,select2]

row.names(sTable.sub) <- sTable.sub$Sample_ID
sTable.sub <- sTable.sub[colnames(selected.expr),]
sTable.sub <- sTable.sub[order(sTable.sub$Last_Known_Treat_Stat),]
selected.expr <- selected.expr[,sTable.sub$Sample_ID]
row.names(selected.expr) <- as.character(lapply(strsplit(row.names(selected.expr), "_"), "[[",2))

df <- data.frame(sTable.sub$Last_Known_Treat_Stat)
row.names(df) <- sTable.sub$Sample_ID
colnames(df) <- c("Status")

pheatmap(selected.expr, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(165),   # RdYlBu
         fontsize_row = 1, fontsize_col = 4, scale = "row", breaks=seq(-4,4, by=0.05), 
         cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, 
         annotation_col=df)

##
heatmap.expr <- assay(vst)
res <- res[order(res$padj),]
select <- subset(res, (log2FoldChange > 0.5 | log2FoldChange < -0.5) & padj < 0.05)$gene; length(select)
selected.expr <- heatmap.expr[select,]
#select2 <- subset(sTable.sub, Last_Known_Treat_Stat!="Healthy")$Sample_ID
#selected.expr <- heatmap.expr[unique(select),select2]
##
sampleDists <- dist(t(selected.expr))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- paste(sTable.sub$Sample_ID)
rownames(sampleDistMatrix) <- paste(sTable.sub$DiseaseCourse)
df <- data.frame(Status=sTable.sub$DiseaseStatus, Course=sTable.sub$DiseaseCourse, Treat=sTable.sub$Last_Known_Treat_Stat)
colors <- colorRampPalette( rev(brewer.pal(9, 'YlOrRd')) )(255)
rownames(df) <- (sTable.sub$Sample_ID)

#dev.off()
#pdf("Heatmap_sample_distances_CD4.pdf", paper="USr", width = 0, height = 0)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, 
         col=colors, fontsize_row = 5, fontsize_col = 4, annotation_col = df, show_colnames = T, fontsize=5)  #cellwidth = 30, cellheight = 30, 
#dev.off()

rm(sampleDists, sampleDistMatrix, df, colors)
#




#####
rm(cds, res, res.sig)
sessionInfo()
#####
