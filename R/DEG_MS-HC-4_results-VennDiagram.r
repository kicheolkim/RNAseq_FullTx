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


### CD4
cds <- cds.cd4
sTable.sub <- sTable.cd4
res <- res.cd4     # both - healthy
res1 <- res1.cd4   # treat - untreat
res2 <- res2.cd4   # treat - healthy
res3 <- res3.cd4   # untreat - healthy
celltype <- "CD4"

### CD8
cds <- cds.cd8
sTable.sub <- sTable.cd8
res <- res.cd8     # both - healthy
res1 <- res1.cd8
res2 <- res2.cd8   # treat - healthy
res3 <- res3.cd8   # untreat - healthy
celltype <- "CD8"

### CD14
cds <- cds.cd14
sTable.sub <- sTable.cd14
res <- res.cd14     # both - healthy
res1 <- res1.cd14
res2 <- res2.cd14   # treat - healthy
res3 <- res3.cd14   # untreat - healthy
celltype <- "CD14"
#


##### Venn Diagram for all  #####
library("VennDiagram")
library("RColorBrewer")
library("gplots")

# genes (baseMean > 3) and (adjusted p-value < 0.05)
venn <- list("res1"=subset(res1, baseMean > 3 & padj < 0.05)$gene,           # (log2FoldChange > 0.5 | log2FoldChange < -0.5) & 
             "res2"=subset(res2, baseMean > 3 & padj < 0.05)$gene, 
             "res3"=subset(res3, baseMean > 3 & padj < 0.05)$gene)
venn.plot <- venn.diagram(venn, NULL, resolution=600, fill=RColorBrewer::brewer.pal(3, "Dark2"), col="transparent", alpha=c(0.7,0.7,0.7), 
                          sub.cex=1, cat.cex=1.2,
                          category=c("Trt-Untrt","Trt-HC","Untrt-HC"))
grid.newpage(); grid.draw(venn.plot)
grid.draw(venn.plot)

### extract list from venn diagram
tmp.list <- venn(venn, show.plot=FALSE); str(tmp.list)
tmp.inters <- attr(tmp.list,"intersections")
lapply(tmp.inters, head, n=10)
#ven.list1 <- data.frame("CD4"=tmp.inters$CD14)
#
res_overlap <- data.frame(geneid=tmp.inters$res2, overlap="Trt-HC")
res_overlap <- rbind(res_overlap, data.frame(geneid=tmp.inters$res3, overlap="Untrt-HC"))
res_overlap <- rbind(res_overlap, data.frame(geneid=tmp.inters$`res2:res3`, overlap="Trt-HC:Untrt-HC"))
write.csv(res_overlap, paste0("Venn_",celltype,"_DEG_overlap-list.csv"))

###
res3$venn <- NA
res3[res3$gene %in% tmp.inters$`res2:res3`,]$venn <- "Trt-HC:Untrt-HC"
res3[res3$gene %in% tmp.inters$res3,]$venn <- "Untrt-HC"
write.csv(res3, paste0("MS-HC_",celltype,"_result-TNaive.vs.HC_venn.csv"))
#
res2$venn <- NA
res2[res2$gene %in% tmp.inters$`res2:res3`,]$venn <- "Trt-HC:Untrt-HC"
res2[res2$gene %in% tmp.inters$res2,]$venn <- "Trt-HC"
write.csv(res2, paste0("MS-HC_",celltype,"_result-Trt.vs.HC_venn.csv"))
#
rm(venn, venn.plot, tmp.list, tmp.inters)





#####
rm(cds, res, res.sig)
sessionInfo()
#####
