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


##### summary of overall DEGs
setwd("~/epic.neb/analysis/results/gene-level/MS_all-HC_updated")
res.cd4 <- read_csv("MS-HC_CD4_result-Trt+Untrt.vs.HC.csv", col_types = cols(X1 = col_skip()))
res1.cd4 <- read_csv("MS-HC_CD4_result-Trt.vs.Untrt.csv", col_types = cols(X1 = col_skip()))
res2.cd4 <- read_csv("MS-HC_CD4_result-Trt.vs.HC.csv", col_types = cols(X1 = col_skip()))
res3.cd4 <- read_csv("MS-HC_CD4_result-Untrt.vs.HC.csv", col_types = cols(X1 = col_skip()))
res_sex.cd4 <- read_csv("MS-HC_CD4_result-covar_Sex.csv", col_types = cols(X1 = col_skip()))
res_age.cd4 <- read_csv("MS-HC_CD4_result-covar_Age.csv", col_types = cols(X1 = col_skip()))

res.cd8 <- read_csv("MS-HC_CD8_result-Trt+Untrt.vs.HC.csv", col_types = cols(X1 = col_skip()))
res1.cd8 <- read_csv("MS-HC_CD8_result-Trt.vs.Untrt.csv", col_types = cols(X1 = col_skip()))
res2.cd8 <- read_csv("MS-HC_CD8_result-Trt.vs.HC.csv", col_types = cols(X1 = col_skip()))
res3.cd8 <- read_csv("MS-HC_CD8_result-Untrt.vs.HC.csv", col_types = cols(X1 = col_skip()))
res_sex.cd8 <- read_csv("MS-HC_CD8_result-covar_Sex.csv", col_types = cols(X1 = col_skip()))
res_age.cd8 <- read_csv("MS-HC_CD8_result-covar_Age.csv", col_types = cols(X1 = col_skip()))

res.cd14 <- read_csv("MS-HC_CD14_result-Trt+Untrt.vs.HC.csv", col_types = cols(X1 = col_skip()))
res1.cd14 <- read_csv("MS-HC_CD14_result-Trt.vs.Untrt.csv", col_types = cols(X1 = col_skip()))
res2.cd14 <- read_csv("MS-HC_CD14_result-Trt.vs.HC.csv", col_types = cols(X1 = col_skip()))
res3.cd14 <- read_csv("MS-HC_CD14_result-Untrt.vs.HC.csv", col_types = cols(X1 = col_skip()))
res_sex.cd14 <- read_csv("MS-HC_CD14_result-covar_Sex.csv", col_types = cols(X1 = col_skip()))
res_age.cd14 <- read_csv("MS-HC_CD14_result-covar_Age.csv", col_types = cols(X1 = col_skip()))


summary_deg(res3.cd4)
summary_deg(res2.cd4)
summary_deg(res1.cd4)
summary_deg(res_sex.cd4)
res_sex.cd4.sub <- subset(res_sex.cd4, seqnames=="chrX" | seqnames=="chrY")
summary_deg(res_sex.cd4.sub)
summary_deg(res_age.cd4)

summary_deg(res3.cd8)
summary_deg(res2.cd8)
result_table <- res2.cd8
rm(result_table)
summary_deg(res1.cd8)
summary_deg(res_sex.cd8)
res_sex.cd8.sub <- subset(res_sex.cd8, seqnames=="chrX" | seqnames=="chrY")
summary_deg(res_sex.cd8.sub)
summary_deg(res_age.cd8)




##### barplot of DEGs group #####
library(ggplot2)
library(gridExtra)
library(plyr)

## treated
res2.cd4_sub <- subset(res2.cd4, baseMean > 3 & padj < 0.05)
res2.cd4_sub$cell <- "CD4"
res2.cd8_sub <- subset(res2.cd8, baseMean > 3 & padj < 0.05)
res2.cd8_sub$cell <- "CD8"
res2.cd14_sub <- subset(res2.cd14, baseMean > 3 & padj < 0.05)
res2.cd14_sub$cell <- "CD14"

res2_sub <- rbind(res2.cd4_sub, res2.cd8_sub)
res2_sub <- rbind(res2_sub, res2.cd14_sub)
nrow(res2_sub)

unique(res2_sub$gene_type)
res2_sub$type <- "Non-coding"
res2_sub[res2_sub$gene_type == "protein_coding",]$type <- "Protein-coding"
#res2_sub[res2_sub$gene_type == "TEC",]$type <- "Protein-coding"
res2_sub[res2_sub$gene_type == "IG_C_gene",]$type <- "Protein-coding"
res2_sub[res2_sub$gene_type == "TR_V_gene",]$type <- "Protein-coding"
res2_sub[res2_sub$gene_type == "TR_C_gene",]$type <- "Protein-coding"

res2_sub$l2fc <- "down"
res2_sub[res2_sub$log2FoldChange > 0,]$l2fc <- "up"
res2_sub$treat <- "Treated vs HC"


## untreated
res3.cd4_sub <- subset(res3.cd4, baseMean > 3 & padj < 0.05)
res3.cd4_sub$cell <- "CD4"
res3.cd8_sub <- subset(res3.cd8, baseMean > 3 & padj < 0.05)
res3.cd8_sub$cell <- "CD8"
res3.cd14_sub <- subset(res3.cd14, baseMean > 3 & padj < 0.05)
res3.cd14_sub$cell <- "CD14"

res3_sub <- rbind(res3.cd4_sub, res3.cd8_sub)
res3_sub <- rbind(res3_sub, res3.cd14_sub)
nrow(res3_sub)

unique(res3_sub$gene_type)
res3_sub$type <- "Non-coding"
res3_sub[res3_sub$gene_type == "protein_coding",]$type <- "Protein-coding"
#res3_sub[res3_sub$gene_type == "TEC",]$type <- "Protein-coding"
res3_sub[res3_sub$gene_type == "IG_C_gene",]$type <- "Protein-coding"
res3_sub[res3_sub$gene_type == "TR_V_gene",]$type <- "Protein-coding"

res3_sub$l2fc <- "down"
res3_sub[res3_sub$log2FoldChange > 0,]$l2fc <- "up"
res3_sub$treat <- "Untreated vs HC"


## plot
res_sub <- rbind(res2_sub, res3_sub)
res_sub$l2fc <- factor(res_sub$l2fc, levels=c("up","down"))
res_sub$type <- factor(res_sub$type, levels=c("Protein-coding","Non-coding"))

res_sub_cd4 <- subset(res_sub, cell=="CD4")
res_sub_cd8 <- subset(res_sub, cell=="CD8")
res_sub_cd14 <- subset(res_sub, cell=="CD14")

p1 <- ggplot(res_sub_cd4, aes(x=type, fill=treat)) + ylim(0,430) +
  geom_bar(aes(fill=l2fc), stat="count") + scale_colour_brewer(palette="BuPu") + 
  facet_wrap(~treat, nrow=1) + theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1))
p2 <- ggplot(res_sub_cd8, aes(x=type, fill=treat)) + ylim(0,430) +
  geom_bar(aes(fill=l2fc), stat="count") + scale_colour_brewer(palette="BuPu") + 
  facet_wrap(~treat, nrow=1) + theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1))
p3 <- ggplot(res_sub_cd14, aes(x=type, fill=treat)) + ylim(0,430) +
  geom_bar(aes(fill=l2fc), stat="count") + scale_colour_brewer(palette="BuPu") + 
  facet_wrap(~treat, nrow=1) + theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1))

grid.arrange(p1, p2, p3, nrow = 1)



count(res2_sub, c("cell","treat","type","l2fc"))
count(res3_sub, c("cell","treat","type","l2fc"))





#####
rm(cds, res, res.sig)
sessionInfo()
#####
