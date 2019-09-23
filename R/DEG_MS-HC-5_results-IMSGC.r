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


#####  check gene name from DEG with IMSGC (comparison with IMSGC genes)  #####
setwd("~/epic.neb/analysis/results/gene-level/MS_all-HC_updated")


###
IMSGC_genes_GW_effects <- read_csv("~/epic.neb/analysis/IMSGC_genes_GW_effects.csv")
IMSGC_genes_cis_eQTL <- read_csv("~/epic.neb/analysis/IMSGC_genes_cis_eQTL.csv")
IMSGC_genes_priotirize <- read_csv("~/epic.neb/analysis/IMSGC_genes_priotirize.csv")

# CD4
cds <- cds.cd4
sTable.sub <- sTable.cd4
celltype <- "CD4"
res <- res.cd4     # both - healthy
res1 <- res1.cd4
res2 <- res2.cd4   # treat - healthy
res3 <- res3.cd4   # untreat - healthy

# CD8
cds <- cds.cd8
sTable.sub <- sTable.cd8
celltype <- "CD8"
res <- res.cd8     # both - healthy
res1 <- res1.cd8
res2 <- res2.cd8   # treat - healthy
res3 <- res3.cd8   # untreat - healthy

# CD14
cds <- cds.cd14
sTable.sub <- sTable.cd14
celltype <- "CD14"
res <- res.cd14     # both - healthy
res1 <- res1.cd14
res2 <- res2.cd14   # treat - healthy
res3 <- res3.cd14   # untreat - healthy


###
select <- subset(res3, baseMean > 3 & padj < 0.05); nrow(select)
#select <- subset(res3, baseMean > 3 & padj < 0.1); nrow(select)

###
expr = paste(IMSGC_genes_GW_effects$`Proximal Gene(s)`[!is.na(IMSGC_genes_GW_effects$`Proximal Gene(s)`)], collapse = "|")
overlap_GW <- select[grepl(expr, select$gene_name),]
#tmp_imsgc = IMSGC_genes_GW_effects[IMSGC_genes_GW_effects$`Proximal Gene(s)`%in% overlap_GW$gene_name,]
overlap_GW$IMSGC <- "IMSGC_GW_effects"

expr = paste(IMSGC_genes_cis_eQTL$gene[!is.na(IMSGC_genes_cis_eQTL$gene)], collapse = "|")
overlap_eQTL <- select[grepl(expr, select$gene_name),]
overlap_eQTL$IMSGC <- "IMSGC_eQTL"

expr = paste(IMSGC_genes_priotirize$`Exonic genes`[!is.na(IMSGC_genes_priotirize$`Exonic genes`)], collapse = "|")
overlap_pri_exonic <- select[grepl(expr, select$gene_name),]
overlap_pri_exonic$IMSGC <- "pri_exonic"

expr = paste(IMSGC_genes_priotirize$`eQTL genes`[!is.na(IMSGC_genes_priotirize$`eQTL genes`)], collapse = "|")
overlap_pri_eQTL <- select[grepl(expr, select$gene_name),]
overlap_pri_eQTL$IMSGC <- "pri_eqtl"

expr = paste(IMSGC_genes_priotirize$`Regulatory network`[!is.na(IMSGC_genes_priotirize$`Regulatory network`)], collapse = "|")
overlap_pri_regnet <- select[grepl(expr, select$gene_name),]
overlap_pri_regnet$IMSGC <- "pri_regnet"

###
overlap_all <- rbind(overlap_GW, overlap_eQTL)
overlap_all <- rbind(overlap_all, overlap_pri_exonic)
overlap_all <- rbind(overlap_all, overlap_pri_eQTL)
overlap_all <- rbind(overlap_all, overlap_pri_regnet)


write.csv(overlap_all, file=paste0("IMSGC_overlap_genes_",celltype,"_p0.05.csv"))
#write.csv(overlap_all, file=paste0("~/epic.neb/analysis/results/IMSGC_overlap_genes_",celltype,"_p0.1.csv"))
#



#####
rm(cds, res, res.sig)
sessionInfo()
#####
