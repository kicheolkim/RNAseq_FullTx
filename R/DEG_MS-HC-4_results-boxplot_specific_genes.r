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



##########  Plotting results for selected gene    ##########

## custom single gene boxplot using ggplot2
library(DESeq2)
library(ggplot2)
library(ggbeeswarm)

##### load all dataset
load("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-AllCells_MSvsHC_withAge_updated.RData")


#####  single gene boxplot  ##### selected genes
target <- "ENSG00000159593.14_NAE1"
sig_gene_list <- subset(res, res$gene==target)
i=1
ann_colors = list(Status = c(CIS="#fe9929", RMS="#e34a33", HC="#31a354"))


setwd("~/epic.neb/analysis/results/gene-level/MS_all-HC_updated/plots/boxplots_untreated")
setwd("~/epic.neb/analysis/results/gene-level/MS_all-HC_updated/plots/boxplots_treated")

##### CD4 
cds <- cds.cd4
#res <- res2.cd4   # treat - healthy
res <- res3.cd4   # untreat - healthy
cell <- "CD4"

### select gene for boxplot     # ENSG00000184557.4_SOCS3, ENSG00000126861.4_OMG, ENSG00000163930.9_BAP1
sig_gene_list <- subset(res, baseMean > 3 & padj < 0.05); nrow(sig_gene_list)
sig_gene_list <- sig_gene_list[order(sig_gene_list$padj),]

for (i in 1:nrow(sig_gene_list)){
  gene <- sig_gene_list$gene[i]
  genename <- strsplit(gene,"_")[[1]][2]; genename
  # plot
  geneCounts <- plotCounts(cds, gene = gene, intgroup = c("Last_Known_Treat_Stat"),
                           returnData = TRUE)
  #geneCounts <- subset(geneCounts, geneCounts$Last_Known_Treat_Stat != "Treated")
  colnames(geneCounts) <- c("count","Status")
  geneCounts$Status <- as.character(geneCounts$Status)
  geneCounts[geneCounts$Status == "TreatmentNaive","Status"] <- "MS-untreated"
  geneCounts[geneCounts$Status == "Treated","Status"] <- "MS-treated"
  geneCounts$Status <- factor(geneCounts$Status, levels=c("MS-untreated","MS-treated","Healthy"))
  
  geneCounts <- geneCounts[geneCounts$Status != "MS-treated",]   # remove treated patients
  
  ggplot(geneCounts, aes(x = Status, y = count, color = Status)) + geom_boxplot() + 
    scale_y_log10() +  geom_beeswarm(cex = 2.5) + scale_color_manual(values=c("#fe9929", "#31a354", "#e34a33")) +
    ggtitle(paste0(genename, " (adj. p-val = ",round(res[res$gene==gene,"padj"],digit=6),") in ", cell)) + xlab("") + 
    ylab("log10(DESeq2 normalized count)") + theme_bw()
  # save plot
  ggsave(paste0(cell,"_boxPlot_",i,"_",genename,".pdf"), plot = last_plot(), device="pdf",
         scale = 1, width = 4.66, height = 5.34, units = "in")
  ggsave(paste0(cell,"_boxPlot_",i,"_",genename,".png"), plot = last_plot(), device="png",
         scale = 1, width = 4.66, height = 5.34, units = "in")
}
rm(sig_gene_list, cds, res, cell)



##### CD8
cds <- cds.cd8
res <- res2.cd8   # treat - healthy
#res <- res3.cd8   # untreat - healthy
cell <- "CD8"

### select gene for boxplot     # ENSG00000184557.4_SOCS3, ENSG00000126861.4_OMG, ENSG00000163930.9_BAP1
sig_gene_list <- subset(res, baseMean > 3 & padj < 0.05); nrow(sig_gene_list)
sig_gene_list <- sig_gene_list[order(sig_gene_list$padj),]

for (i in 1:nrow(sig_gene_list)){
  gene <- sig_gene_list$gene[i]
  genename <- strsplit(gene,"_")[[1]][2]; genename
  # plot
  geneCounts <- plotCounts(cds, gene = gene, intgroup = c("Last_Known_Treat_Stat"),
                           returnData = TRUE)
  #geneCounts <- subset(geneCounts, geneCounts$Last_Known_Treat_Stat != "Treated")
  colnames(geneCounts) <- c("count","Status")
  geneCounts$Status <- as.character(geneCounts$Status)
  geneCounts[geneCounts$Status == "TreatmentNaive","Status"] <- "MS-untreated"
  geneCounts[geneCounts$Status == "Treated","Status"] <- "MS-treated"
  geneCounts$Status <- factor(geneCounts$Status, levels=c("MS-untreated","MS-treated","Healthy"))
  
  geneCounts <- geneCounts[geneCounts$Status != "MS-treated",]   # remove treated patients
  
  ggplot(geneCounts, aes(x = Status, y = count, color = Status)) + geom_boxplot() + 
    scale_y_log10() +  geom_beeswarm(cex = 2.5) + scale_color_manual(values=c("#fe9929", "#31a354", "#e34a33")) +
    ggtitle(paste0(genename, " (adj. p-val = ",round(res[res$gene==gene,"padj"],digit=6),") in ", cell)) + xlab("") + 
    ylab("log10(DESeq2 normalized count)") + theme_bw()
  # save plot
  ggsave(paste0(cell,"_boxPlot_",i,"_",genename,".pdf"), plot = last_plot(), device="pdf",
         scale = 1, width = 4.66, height = 5.34, units = "in")
  ggsave(paste0(cell,"_boxPlot_",i,"_",genename,".png"), plot = last_plot(), device="png",
         scale = 1, width = 4.66, height = 5.34, units = "in")
}
rm(sig_gene_list, cds, res, cell)


##### CD14
cds <- cds.cd14
#res <- res2.cd14   # treat - healthy
res <- res3.cd14   # untreat - healthy
cell <- "CD14"

### select gene for boxplot     # ENSG00000184557.4_SOCS3, ENSG00000126861.4_OMG, ENSG00000163930.9_BAP1
sig_gene_list <- subset(res, baseMean > 3 & padj < 0.05); nrow(sig_gene_list)
sig_gene_list <- sig_gene_list[order(sig_gene_list$padj),]

for (i in 1:nrow(sig_gene_list)){
  gene <- sig_gene_list$gene[i]
  genename <- strsplit(gene,"_")[[1]][2]; genename
  # plot
  geneCounts <- plotCounts(cds, gene = gene, intgroup = c("Last_Known_Treat_Stat"),
                           returnData = TRUE)
  #geneCounts <- subset(geneCounts, geneCounts$Last_Known_Treat_Stat != "Treated")
  colnames(geneCounts) <- c("count","Status")
  geneCounts$Status <- as.character(geneCounts$Status)
  geneCounts[geneCounts$Status == "TreatmentNaive","Status"] <- "MS-untreated"
  geneCounts[geneCounts$Status == "Treated","Status"] <- "MS-treated"
  geneCounts$Status <- factor(geneCounts$Status, levels=c("MS-untreated","MS-treated","Healthy"))
  
  geneCounts <- geneCounts[geneCounts$Status != "MS-treated",]   # remove treated patients
  
  ggplot(geneCounts, aes(x = Status, y = count, color = Status)) + geom_boxplot() + 
    scale_y_log10() +  geom_beeswarm(cex = 2.5) + scale_color_manual(values=c("#fe9929", "#31a354", "#e34a33")) +
    ggtitle(paste0(genename, " (adj. p-val = ",round(res[res$gene==gene,"padj"],digit=6),") in ", cell)) + xlab("") + 
    ylab("log10(DESeq2 normalized count)") + theme_bw()
  # save plot
  ggsave(paste0(cell,"_boxPlot_",i,"_",genename,".pdf"), plot = last_plot(), device="pdf",
         scale = 1, width = 4.66, height = 5.34, units = "in")
  ggsave(paste0(cell,"_boxPlot_",i,"_",genename,".png"), plot = last_plot(), device="png",
         scale = 1, width = 4.66, height = 5.34, units = "in")
  #
}
rm(sig_gene_list, cds, res, cell)
#####




#####
rm(cds, res, res.sig)
sessionInfo()
#####
