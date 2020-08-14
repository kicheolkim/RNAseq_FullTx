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



##########  Heatmap for DEGs  ##########
setwd("~/epic.neb/analysis/results/gene-level/MS_all-HC_updated/plots")
library("pheatmap"); library("RColorBrewer")

##### heatmap for genes from wgcna network module
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = -sv, agglo.FUN=mean)
  as.hclust(dend)
}
ann_colors = list(Status = c(MS="orangered2", HC="deepskyblue1"))
#ann_colors = list(Status = c(TreatmentNaive="orangered2", Healthy="deepskyblue1", Treated="purple2"))


# CD4
cds <- cds.cd4
#sTable.sub <- sTable.cd4
sTable.sub <- subset(sTable.cd4, Subset=="CD4" & (Last_Known_Treat_Stat=="TreatmentNaive" | Last_Known_Treat_Stat=="Healthy"))
celltype <- "CD4"
res <- res.cd4     # both - healthy
res1 <- res1.cd4
res2 <- res2.cd4   # treat - healthy
res3 <- res3.cd4   # untreat - healthy

# CD8
cds <- cds.cd8
#sTable.sub <- sTable.cd8
sTable.sub <- subset(sTable.cd8, Subset=="CD8" & (Last_Known_Treat_Stat=="TreatmentNaive" | Last_Known_Treat_Stat=="Healthy"))
celltype <- "CD8"
res <- res.cd8     # both - healthy
res1 <- res1.cd8
res2 <- res2.cd8   # treat - healthy
res3 <- res3.cd8   # untreat - healthy

# CD14
cds <- cds.cd14
#sTable.sub <- sTable.cd14
sTable.sub <- subset(sTable.cd14, Subset=="CD14" & (Last_Known_Treat_Stat=="TreatmentNaive" | Last_Known_Treat_Stat=="Healthy"))
celltype <- "CD14"
res <- res.cd14     # both - healthy
res1 <- res1.cd14
res2 <- res2.cd14   # treat - healthy
res3 <- res3.cd14   # untreat - healthy


#####
#heatmap.expr <- assay(vst(cds, blind=FALSE))
heatmap.expr <- assay(normTransform(cds))
heatmap.expr <- heatmap.expr - rowMeans(heatmap.expr)
heatmap.expr <- heatmap.expr[,colnames(heatmap.expr) %in% sTable.sub$Sample_ID]

### untreated vs healthy by adjusted p-value
select <- subset(res3, baseMean > 3 & padj < 0.05 & gene_type == "protein_coding"); nrow(select)
elected.expr <- heatmap.expr[select$gene,]

# make heatmap
row.names(sTable.sub) <- sTable.sub$Sample_ID
sTable.sub <- sTable.sub[colnames(selected.expr),]
sTable.sub <- sTable.sub[order(sTable.sub$Last_Known_Treat_Stat),]
selected.expr <- selected.expr[,sTable.sub$Sample_ID]
row.names(selected.expr) <- as.character(lapply(strsplit(row.names(selected.expr), "_"), "[[",2))
sTable.sub$Last_Known_Treat_Stat <- droplevels(sTable.sub$Last_Known_Treat_Stat)

df <- data.frame(sTable.sub$Last_Known_Treat_Stat)
row.names(df) <- sTable.sub$Sample_ID
colnames(df) <- c("Status")

df$Status <- as.character(df$Status)
df[df$Status == "TreatmentNaive",] <- "MS"
df[df$Status == "Healthy",] <- "HC"


# for CD4 result
pheatmap(selected.expr, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(165),   # RdYlBu
         fontsize_row = 1, fontsize_col = 3, scale = "row", breaks=seq(-4,4, by=0.05), 
         cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, show_colnames=FALSE, 
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D2", 
         cutree_rows=3, cutree_cols=2, main=paste0(celltype," (FDR < 0.05 & baseMean > 3)"), clustering_callback=callback,
         annotation_col=df, annotation_colors=ann_colors)

# for CD8 result
pheatmap(selected.expr, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(165),   # RdYlBu
         fontsize_row = 7, fontsize_col = 3, scale = "row", breaks=seq(-4,4, by=0.05), 
         cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, show_colnames=FALSE, 
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D2", 
         cutree_rows=2, cutree_cols=2, main=paste0(celltype," (FDR < 0.05 & baseMean > 3)"), clustering_callback=callback,
         annotation_col=df, annotation_colors=ann_colors)

# for CD14 result
pheatmap(selected.expr, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(165),   # RdYlBu
         fontsize_row = 1, fontsize_col = 3, scale = "row", breaks=seq(-4,4, by=0.05), 
         cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, show_colnames=FALSE, 
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D2", 
         cutree_rows=3, cutree_cols=2, main=paste0(celltype," (FDR < 0.05 & baseMean > 3)"), clustering_callback=callback,
         annotation_col=df, annotation_colors=ann_colors)
#
#




#####
rm(cds, res, res.sig)
sessionInfo()
#####
