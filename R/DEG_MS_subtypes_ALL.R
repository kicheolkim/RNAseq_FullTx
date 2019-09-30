# -----
# title: "DEG analysis (DESeq2) - MS-subtype (Treatment Naive patients) vs Healthy controls"
# -----
##########=====  DEG using DESeq2  =====##########
library(DESeq2)
library(tximport)
library(readr)
library(RColorBrewer)
library(pheatmap)

##### loading metadata
setwd("~/epic.neb/analysis/results/gene-level/MS_disease_subset_updated")
# loading base information tables
load(file="~/epic.neb/analysis/R_neb/RData/infoTable_filtered_updated.RData")


##### load all dataset - saved DESeq2 results
load("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-MS_naive-subset_AllCells.RData")



########## ====================================================================================================================3
# reference for model matrix (consider treatment status) => https://support.bioconductor.org/p/64844/#64861
# check number of samples
sTable.tmp <- subset(sTable, Subset=="CD14" & (Last_Known_Treat_Stat!="Treated") &
                       (DiseaseCourse=="CIS"|DiseaseCourse=="PP"|DiseaseCourse=="RR"|DiseaseCourse=="Healthy") &
                       !is.na(Last_Known_Treat_Stat))
table(sTable.tmp$Last_Known_Treat_Stat,sTable.tmp$DiseaseCourse)
rm(sTable.tmp)



########## ====================================================================================================================3
##### CD4: loading selected samples from RSEM #####
# subset metadata
sTable.cd4 <- subset(sTable, Subset=="CD4" & (Last_Known_Treat_Stat=="TreatmentNaive") &
                       (DiseaseCourse=="CIS"|DiseaseCourse=="PP"|DiseaseCourse=="SP"|DiseaseCourse=="RR"))
sTable.cd4$DiseaseCourse2 <- sTable.cd4$DiseaseCourse
sTable.cd4[sTable.cd4$DiseaseCourse=="SP",]$DiseaseCourse <- "PP"

factor_cols <- c("Subset","DiseaseCourse","Last_Known_Treat_Stat","DiseaseStatus")
sTable.cd4[factor_cols] <- lapply(sTable.cd4[factor_cols], droplevels)


# file list for rsem
sTable.sub <- sTable.cd4
setwd("~/epic.neb/counts/rsem")
files <- list.files(path=".", pattern="*.genes.results", recursive = TRUE)
tmp <- unlist(lapply(strsplit(files, "/"), "[[",2))
names(files) <- sub(".genes.results", "", tmp)
files <- files[names(files) %in% sTable.sub$Sample_ID]

# matching list order
row.names(sTable.sub) <- sTable.sub$Sample_ID
sTable.sub <- sTable.sub[names(files),]
### loading gene counts - subset
rsem.gene.sub <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
rm(files, tmp)

### input and run DESeq2 for subset  ###
rsem.gene.sub$length[rsem.gene.sub$length == 0] <- 1
cds.sub <- DESeqDataSetFromTximport(rsem.gene.sub, sTable.sub, ~DiseaseCourse+Sex+AgeAtExamGrp-1)
cds.sub <- cds.sub[ rowSums(counts(cds.sub)) > 1, ]
cds.sub <- DESeq(cds.sub)

cds.cd4 <- cds.sub
vst.cd4 <- vst(cds.cd4)
sTable.cd4 <- sTable.sub

rm(cds.sub, sTable.sub, rsem.gene.sub)



##### CD8: loading selected samples from RSEM #####
# subset metadata
sTable.cd8 <- subset(sTable, Subset=="CD8" & (Last_Known_Treat_Stat=="TreatmentNaive") &
                       (DiseaseCourse=="CIS"|DiseaseCourse=="PP"|DiseaseCourse=="SP"|DiseaseCourse=="RR"))
sTable.cd8$DiseaseCourse2 <- sTable.cd8$DiseaseCourse
sTable.cd8[sTable.cd8$DiseaseCourse=="SP",]$DiseaseCourse <- "PP"
factor_cols <- c("Subset","DiseaseCourse","Last_Known_Treat_Stat","DiseaseStatus")
sTable.cd8[factor_cols] <- lapply(sTable.cd8[factor_cols], droplevels)

# file list for rsem
sTable.sub <- sTable.cd8
setwd("~/epic.neb/counts/rsem")
files <- list.files(path=".", pattern="*.genes.results", recursive = TRUE)
tmp <- unlist(lapply(strsplit(files, "/"), "[[",2))
names(files) <- sub(".genes.results", "", tmp)
files <- files[names(files) %in% sTable.sub$Sample_ID]
# matching list order
row.names(sTable.sub) <- sTable.sub$Sample_ID
sTable.sub <- sTable.sub[names(files),]
### loading gene counts - subset
rsem.gene.sub <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
rm(files, tmp)

### input and run DESeq2 for subset  ###
rsem.gene.sub$length[rsem.gene.sub$length == 0] <- 1
cds.sub <- DESeqDataSetFromTximport(rsem.gene.sub, sTable.sub, ~DiseaseCourse+Sex+AgeAtExamGrp-1)
cds.sub <- cds.sub[ rowSums(counts(cds.sub)) > 1, ]
cds.sub <- DESeq(cds.sub)

cds.cd8 <- cds.sub
vst.cd8 <- vst(cds.cd8)
sTable.cd8 <- sTable.sub

rm(cds.sub, sTable.sub, rsem.gene.sub)



##### CD14: loading selected samples from RSEM #####
# subset metadata
sTable.cd14 <- subset(sTable, Subset=="CD14" & (Last_Known_Treat_Stat=="TreatmentNaive") &
                        (DiseaseCourse=="CIS"|DiseaseCourse=="PP"|DiseaseCourse=="SP"|DiseaseCourse=="RR"))
sTable.cd14$DiseaseCourse2 <- sTable.cd14$DiseaseCourse
sTable.cd14[sTable.cd14$DiseaseCourse=="SP",]$DiseaseCourse <- "PP"
factor_cols <- c("Subset","DiseaseCourse","Last_Known_Treat_Stat","DiseaseStatus")
sTable.cd14[factor_cols] <- lapply(sTable.cd14[factor_cols], droplevels)


# file list for rsem
sTable.sub <- sTable.cd14
setwd("~/epic.neb/counts/rsem")
files <- list.files(path=".", pattern="*.genes.results", recursive = TRUE)
tmp <- unlist(lapply(strsplit(files, "/"), "[[",2))
names(files) <- sub(".genes.results", "", tmp)
files <- files[names(files) %in% sTable.sub$Sample_ID]
# matching list order
row.names(sTable.sub) <- sTable.sub$Sample_ID
sTable.sub <- sTable.sub[names(files),]
### loading gene counts - subset
rsem.gene.sub <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
rm(files, tmp)

### input and run DESeq2 for subset  ###
rsem.gene.sub$length[rsem.gene.sub$length == 0] <- 1
cds.sub <- DESeqDataSetFromTximport(rsem.gene.sub, sTable.sub, ~DiseaseCourse+Sex+AgeAtExamGrp-1)
cds.sub <- cds.sub[ rowSums(counts(cds.sub)) > 1, ]
cds.sub <- DESeq(cds.sub)

cds.cd14 <- cds.sub
vst.cd14 <- vst(cds.cd14)
sTable.cd14 <- sTable.sub

rm(cds.sub, sTable.sub, rsem.gene.sub)


### SAVE all
save.image("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-MS_naive-subset_AllCells.RData")



##########  CHECK Results  ##########  ##########  ##########  ##########  ##########  ##########
##### CD4: DEG results  ###
cds <- cds.cd4

# res1: CIS_TreatmentNaive vs PP_TreatmentNaive
res1 <- results(cds, alpha=0.05, contrast=c("DiseaseCourse","CIS","PP"))
summary(res1)   # CIS_TreatmentNaive vs PP_TreatmentNaive
# res2: CIS_TreatmentNaive vs RR_TreatmentNaive
res2 <- results(cds, alpha=0.05, contrast=c("DiseaseCourse","CIS","RR"))
summary(res2)   # CIS_TreatmentNaive vs RR_TreatmentNaive
# res3: PP_TreatmentNaive vs RR_TreatmentNaive
res3 <- results(cds, alpha=0.05, contrast=c("DiseaseCourse","PP","RR"))
summary(res3)   # PP_TreatmentNaive vs RR_TreatmentNaive
# res4: CIS+RR vs PP
res4 <- results(cds, alpha=0.05, 
                contrast=list(c("DiseaseCourseCIS","DiseaseCourseRR"),"DiseaseCoursePP"),
                listValues = c(1/2, -1))
summary(res4)   # PP_TreatmentNaive vs RR_TreatmentNaive

# convert to data frame: CIS_TreatmentNaive vs PP_TreatmentNaive
res1 <- as.data.frame(res1)
res1$gene <- row.names(res1)
# convert to data frame: CIS_TreatmentNaive vs RR_TreatmentNaive
res2 <- as.data.frame(res2)
res2$gene <- row.names(res2)
# convert to data frame: PP_TreatmentNaive vs RR_TreatmentNaive
res3 <- as.data.frame(res3)
res3$gene <- row.names(res3)
# convert to data frame: CIS+RR vs PP
res4 <- as.data.frame(res4)
res4$gene <- row.names(res4)


## find gene related to Sex
dds <- DESeq(cds, test="LRT", reduced=~DiseaseCourse+AgeAtExamGrp-1)
res.sex <- results(dds, alpha=0.05)
summary(res.sex)   ### sex influenced gene
res.sex <- as.data.frame(res.sex)
res.sex$gene <- row.names(res.sex)

## find gene related to age at exam
dds <- DESeq(cds, test="LRT", reduced=~DiseaseCourse+Sex-1)
res.age <- results(dds, alpha=0.05)
summary(res.age)   ### sex influenced gene
res.age <- as.data.frame(res.age)
res.age$gene <- row.names(res.age)


# res1: CIS_TreatmentNaive vs PP_TreatmentNaive
res1 <- merge(res1, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res1 <- res1[order(res1$padj),]
write.csv(res1, file="CD4_result1-MS_untreat-CISvsPP.csv")
# res2: CIS_TreatmentNaive vs RR_TreatmentNaive
res2 <- merge(res2, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res2 <- res2[order(res2$padj),]
write.csv(res2, file="CD4_result2-MS_untreat-CISvsRR.csv")
# res3: PP_TreatmentNaive vs RR_TreatmentNaive
res3 <- merge(res3, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res3 <- res3[order(res3$padj),]
write.csv(res3, file="CD4_result3-MS_untreat-PPvsRR.csv")
# res4: CIS+RR vs PP
res4 <- merge(res4, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res4 <- res4[order(res4$padj),]
write.csv(res4, file="CD4_result4-MS_untreat-CIS+RRvsPP.csv")

res.sex <- merge(res.sex, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res.sex <- res.sex[order(res.sex$padj),]
write.csv(res.sex, file="CD4_result-MS_untreat-covar_Sex.csv")
res.age <- merge(res.age, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res.age <- res.sex[order(res.age$padj),]
write.csv(res.age, file="CD4_result-MS_untreat-covar_Age.csv")

res1.cd4 <- res1
res2.cd4 <- res2
res3.cd4 <- res3
res4.cd4 <- res4

rm(cds, dds, res.sex, res.age, res1, res2, res3, res4)



##### CD8: DEG results  ###
cds <- cds.cd8

# res1: CIS_TreatmentNaive vs PP_TreatmentNaive
res1 <- results(cds, alpha=0.05, contrast=c("DiseaseCourse","CIS","PP"))
summary(res1)   # CIS_TreatmentNaive vs PP_TreatmentNaive
# res2: CIS_TreatmentNaive vs RR_TreatmentNaive
res2 <- results(cds, alpha=0.05, contrast=c("DiseaseCourse","CIS","RR"))
summary(res2)   # CIS_TreatmentNaive vs RR_TreatmentNaive
# res3: PP_TreatmentNaive vs RR_TreatmentNaive
res3 <- results(cds, alpha=0.05, contrast=c("DiseaseCourse","PP","RR"))
summary(res3)   # PP_TreatmentNaive vs RR_TreatmentNaive
# res4: CIS+RR vs PP
res4 <- results(cds, alpha=0.05, 
                contrast=list(c("DiseaseCourseCIS","DiseaseCourseRR"),"DiseaseCoursePP"),
                listValues = c(1/2, -1))
summary(res4)   # PP_TreatmentNaive vs RR_TreatmentNaive

# convert to data frame: CIS_TreatmentNaive vs PP_TreatmentNaive
res1 <- as.data.frame(res1)
res1$gene <- row.names(res1)
# convert to data frame: CIS_TreatmentNaive vs RR_TreatmentNaive
res2 <- as.data.frame(res2)
res2 <- res2[order(res2$padj),]
res2$gene <- row.names(res2)
# convert to data frame: PP_TreatmentNaive vs RR_TreatmentNaive
res3 <- as.data.frame(res3)
res3$gene <- row.names(res3)
# convert to data frame: CIS+RR vs PP
res4 <- as.data.frame(res4)
res4$gene <- row.names(res4)

## find gene related to Sex
dds <- DESeq(cds, test="LRT", reduced=~DiseaseCourse+AgeAtExamGrp-1)
res.sex <- results(dds, alpha=0.05)
summary(res.sex)   ### sex influenced gene
res.sex <- as.data.frame(res.sex)
res.sex$gene <- row.names(res.sex)

## find gene related to age at exam
dds <- DESeq(cds, test="LRT", reduced=~DiseaseCourse+Sex-1)
res.age <- results(dds, alpha=0.05)
summary(res.age)   ### sex influenced gene
res.age <- as.data.frame(res.age)
res.age$gene <- row.names(res.age)


# res1: CIS_TreatmentNaive vs PP_TreatmentNaive
res1 <- merge(res1, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res1 <- res1[order(res1$padj),]
write.csv(res1, file="CD8_result1-MS_untreat-CISvsPP.csv")
# res2: CIS_TreatmentNaive vs RR_TreatmentNaive
res2 <- merge(res2, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res2 <- res2[order(res2$padj),]
write.csv(res2, file="CD8_result2-MS_untreat-CISvsRR.csv")
# res3: PP_TreatmentNaive vs RR_TreatmentNaive
res3 <- merge(res3, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res3 <- res3[order(res3$padj),]
write.csv(res3, file="CD8_result3-MS_untreat-PPvsRR.csv")
# res4: CIS+RR vs PP
res4 <- merge(res4, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res4 <- res4[order(res4$padj),]
write.csv(res4, file="CD8_result4-MS_untreat-CIS+RRvsPP.csv")

res.sex <- merge(res.sex, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res.sex <- res.sex[order(res.sex$padj),]
write.csv(res.sex, file="CD8_result-MS_untreat-covar_Sex.csv")
res.age <- merge(res.age, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res.age <- res.sex[order(res.age$padj),]
write.csv(res.age, file="CD8_result-MS_untreat-covar_Age.csv")

res1.cd8 <- res1
res2.cd8 <- res2
res3.cd8 <- res3
res4.cd8 <- res4

rm(cds, dds, res.sex, res.age, res1, res2, res3, res4)



##### CD14: DEG results  ###
cds <- cds.cd14

# res1: CIS_TreatmentNaive vs PP_TreatmentNaive
res1 <- results(cds, alpha=0.05, contrast=c("DiseaseCourse","CIS","PP"))
summary(res1)   # CIS_TreatmentNaive vs PP_TreatmentNaive
# res2: CIS_TreatmentNaive vs RR_TreatmentNaive
res2 <- results(cds, alpha=0.05, contrast=c("DiseaseCourse","CIS","RR"))
summary(res2)   # CIS_TreatmentNaive vs RR_TreatmentNaive
# res3: PP_TreatmentNaive vs RR_TreatmentNaive
res3 <- results(cds, alpha=0.05, contrast=c("DiseaseCourse","PP","RR"))
summary(res3)   # PP_TreatmentNaive vs RR_TreatmentNaive
# res4: CIS+RR vs PP
res4 <- results(cds, alpha=0.05, 
                contrast=list(c("DiseaseCourseCIS","DiseaseCourseRR"),"DiseaseCoursePP"),
                listValues = c(1/2, -1))
summary(res4)   # PP_TreatmentNaive vs RR_TreatmentNaive

# convert to data frame: CIS_TreatmentNaive vs PP_TreatmentNaive
res1 <- as.data.frame(res1)
res1$gene <- row.names(res1)
# convert to data frame: CIS_TreatmentNaive vs RR_TreatmentNaive
res2 <- as.data.frame(res2)
res2$gene <- row.names(res2)
# convert to data frame: PP_TreatmentNaive vs RR_TreatmentNaive
res3 <- as.data.frame(res3)
res3$gene <- row.names(res3)
# convert to data frame: CIS+RR vs PP
res4 <- as.data.frame(res4)
res4$gene <- row.names(res4)

## find gene related to Sex
dds <- DESeq(cds, test="LRT", reduced=~DiseaseCourse+AgeAtExamGrp-1)
res.sex <- results(dds, alpha=0.05)
summary(res.sex)   ### sex influenced gene
res.sex <- as.data.frame(res.sex)
res.sex$gene <- row.names(res.sex)

## find gene related to age at exam
dds <- DESeq(cds, test="LRT", reduced=~DiseaseCourse+Sex-1)
res.age <- results(dds, alpha=0.05)
summary(res.age)   ### sex influenced gene
res.age <- as.data.frame(res.age)
res.age$gene <- row.names(res.age)


# res1: CIS_TreatmentNaive vs PP_TreatmentNaive
res1 <- merge(res1, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res1 <- res1[order(res1$padj),]
write.csv(res1, file="CD14_result1-MS_untreat-CISvsPP.csv")
# res2: CIS_TreatmentNaive vs RR_TreatmentNaive
res2 <- merge(res2, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res2 <- res2[order(res2$padj),]
write.csv(res2, file="CD14_result2-MS_untreat-CISvsRR.csv")
# res3: PP_TreatmentNaive vs RR_TreatmentNaive
res3 <- merge(res3, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res3 <- res3[order(res3$padj),]
write.csv(res3, file="CD14_result3-MS_untreat-PPvsRR.csv")
# res4: CIS+RR vs PP
res4 <- merge(res4, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res4 <- res4[order(res4$padj),]
write.csv(res4, file="CD14_result4-MS_untreat-CIS+RRvsPP.csv")

res.sex <- merge(res.sex, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res.sex <- res.sex[order(res.sex$padj),]
write.csv(res.sex, file="CD14_result-MS_untreat-covar_Sex.csv")
res.age <- merge(res.age, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res.age <- res.sex[order(res.age$padj),]
write.csv(res.age, file="CD14_result-MS_untreat-covar_Age.csv")

res1.cd14 <- res1
res2.cd14 <- res2
res3.cd14 <- res3
res4.cd14 <- res4

rm(cds, dds, res.sex, res.age, res1, res2, res3, res4)



### SAVE all
save.image("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-MS_naive-subset_AllCells.RData")






##########  CHECK Results  ########## ########## ########## ########## ########## ########## ##########333
### Load all dataset
load("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-MS_naive-subset_AllCells.RData")
setwd("~/epic.neb/analysis/results/gene-level/MS_disease_subset_updated")


##########  compare between comparison in each cell types  ##########
cds <- cds.cd4
res1 <- res1.cd4   # CISvsPP
res2 <- res2.cd4   # CISvsRR
res3 <- res3.cd4   # PPvsRR
res4 <- res4.cd4   # CIS+RRvsPP
celltype <- "CD4"

#
cds <- cds.cd8
res1 <- res1.cd8
res2 <- res2.cd8
res3 <- res3.cd8
res4 <- res4.cd8
celltype <- "CD8"

#
cds <- cds.cd14
res1 <- res1.cd14
res2 <- res2.cd14
res3 <- res3.cd14
res4 <- res4.cd14
celltype <- "CD14"
#

# res1: CIS_TreatmentNaive vs PP_TreatmentNaive
# res2: CIS_TreatmentNaive vs RR_TreatmentNaive
# res3: PP_TreatmentNaive vs RR_TreatmentNaive
# res4: RR_TreatmentNaive + CIS_TreatmentNaive vs PP_TreatmentNaive

subset(res1, padj < 0.05 & (log2FoldChange > 1))[1:8]
subset(res1, padj < 0.05 & (log2FoldChange < -1))[1:8]
subset(res2, padj < 0.05 & (log2FoldChange < -1 | log2FoldChange > 1))[1:8]
nrow(subset(res3, padj < 0.05 & (log2FoldChange < -1 | log2FoldChange > 1))[1:8])
nrow(subset(res4, padj < 0.05 & (log2FoldChange > 1))[1:8])
nrow(subset(res4, padj < 0.05 & (log2FoldChange < -1))[1:8])



##### Venn Diagram #####
library("VennDiagram")
library("RColorBrewer")
library("gplots")
# genes (baseMean > 3) and (adjusted p-value < 0.05)
venn <- list("res1"=subset(res1, baseMean > 3 & padj < 0.05)$gene,         # CIS-PP
             "res2"=subset(res2, baseMean > 3 & padj < 0.05)$gene,         # CIS-RR
             "res4"=subset(res4, baseMean > 3 & padj < 0.05)$gene,         # (CIS+RR)-PP
             "res3"=subset(res3, baseMean > 3 & padj < 0.05)$gene)         # PP-RR
venn.plot <- venn.diagram(venn, NULL, resolution=600, 
                          fill=RColorBrewer::brewer.pal(4, "Dark2"), col="transparent", alpha=c(0.5,0.5,0.5,0.5), 
                          sub.cex=1, cat.cex=1.2,
                          category=c("CIS-PP","CIS-RR","(CIS+RR)-PP","PP-RR"))
grid.newpage(); grid.draw(venn.plot)

### extract list from venn diagram
tmp.list <- venn(venn, show.plot=FALSE); str(tmp.list)
tmp.inters <- attr(tmp.list,"intersections")
lapply(tmp.inters, head, n=100)
#ven.list1 <- data.frame("CD4"=tmp.inters$CD14)

res_overlap <- data.frame(geneid=tmp.inters$res2, overlap="CIS-RR")
res_overlap <- rbind(res_overlap, data.frame(geneid=tmp.inters$`res1:res4:res3`, overlap="CIS-PP:(CIS+RR)-PP:PP-RR"))
write.csv(res_overlap, paste0("Venn_CD4_DEG_overlap-list.csv"))

res_overlap <- data.frame(geneid=tmp.inters$res2, overlap="CIS-RR")
res_overlap <- rbind(res_overlap, data.frame(geneid=tmp.inters$res3, overlap="PP-RR"))
res_overlap <- rbind(res_overlap, data.frame(geneid=tmp.inters$`res4:res3`, overlap="(CIS+RR)-PP:PP-RR"))
write.csv(res_overlap, paste0("Venn_CD8_DEG_overlap-list.csv"))

res_overlap <- data.frame(geneid=tmp.inters$res1, overlap="CIS-PP")
res_overlap <- rbind(res_overlap, data.frame(geneid=tmp.inters$res2, overlap="CIS-RR"))
res_overlap <- rbind(res_overlap, data.frame(geneid=tmp.inters$res3, overlap="PP-RR"))
write.csv(res_overlap, paste0("Venn_CD14_DEG_overlap-list.csv"))


rm(venn, venn.plot, tmp.list, tmp.inters)
#


#
setwd("~/epic.neb/analysis")
res3$venn <- NA
res3[res3$gene %in% tmp.inters$`res:res2:res3`,]$venn <- "overlap3"
res3[res3$gene %in% tmp.inters$`res:res3`,]$venn <- "overlap2"
res3[res3$gene %in% tmp.inters$res3,]$venn <- "overlap0"

write.csv(res3, "MS-HC_CD4_.csv")
#



##########  Heatmap  ##########
library("pheatmap"); library("RColorBrewer")
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = -sv, agglo.FUN=mean)
  as.hclust(dend)
}

# CD4
cds <- cds.cd4
sTable.sub <- sTable.cd4
res1 <- res1.cd4
res1$comparison <- "CIS-PP"
res2 <- res2.cd4
res2$comparison <- "CIS-RR"
res3 <- res3.cd4
res3$comparison <- "PP-RR"
res4 <- res4.cd4
res4$comparison <- "(CIS+RR)-PP"
res <- rbind(subset(res1, padj<0.05), subset(res2, padj<0.05), subset(res3, padj<0.05), subset(res4, padj<0.05))

# CD8
cds <- cds.cd8
sTable.sub <- sTable.cd8
res1 <- res1.cd8
res1$comparison <- "CIS-PP"
res2 <- res2.cd8
res2$comparison <- "CIS-RR"
res3 <- res3.cd8
res3$comparison <- "PP-RR"
res4 <- res4.cd8
res4$comparison <- "(CIS+RR)-PP"
res <- rbind(subset(res1, padj<0.05), subset(res2, padj<0.05), subset(res3, padj<0.05), subset(res4, padj<0.05))

# CD14
cds <- cds.cd14
sTable.sub <- sTable.cd14
res1 <- res1.cd14
res1$comparison <- "CIS-PP"
res2 <- res2.cd14
res2$comparison <- "CIS-RR"
res3 <- res3.cd14
res3$comparison <- "PP-RR"
res4 <- res4.cd14
res4$comparison <- "(CIS+RR)-PP"
res <- rbind(subset(res1, padj<0.05), subset(res2, padj<0.05), subset(res3, padj<0.05), subset(res4, padj<0.05))
#
View(res)
# res1: CIS_TreatmentNaive vs PP_TreatmentNaive
# res2: CIS_TreatmentNaive vs RR_TreatmentNaive
# res3: PP_TreatmentNaive vs RR_TreatmentNaive
# res4: (CIS_TreatmentNaive+RR_TreatmentNaive) vs PP_TreatmentNaive


###
#heatmap.expr <- assay(vst(cds, blind=FALSE))
heatmap.expr <- assay(normTransform(cds))
heatmap.expr <- heatmap.expr - rowMeans(heatmap.expr)


# padj < 0.1 
select <- subset(res2, padj<0.1 & baseMean > 10); nrow(select)
select <- select[order(select$padj),]
selected.expr <- heatmap.expr[select$gene,]


# padj < 0.1 from (CIS+RR)-PP
select <- subset(res4, padj<0.1); nrow(select)
select <- select[order(select$padj),]
selected.expr <- heatmap.expr[select$gene,]


# table for heatmap
row.names(sTable.sub) <- sTable.sub$Sample_ID
sTable.sub <- sTable.sub[colnames(selected.expr),]
sTable.sub <- sTable.sub[order(sTable.sub$DiseaseCourse),]
selected.expr <- selected.expr[,sTable.sub$Sample_ID]
row.names(selected.expr) <- as.character(lapply(strsplit(row.names(selected.expr), "_"), "[[",2))

df <- data.frame(sTable.sub$DiseaseCourse)
row.names(df) <- sTable.sub$Sample_ID
colnames(df) <- c("Status")
#
pheatmap(selected.expr, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(165),   # RdYlBu
         fontsize_row = 8, fontsize_col = 4, scale = "row", breaks=seq(-4,4, by=0.05), 
         cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, 
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "ward.D2", 
         cutree_rows=1, cutree_cols=2, main="CD14: CIS vs RR (FDR < 0.1 & baseMean > 10)",
         annotation_col=df, clustering_callback=callback)
#


###
select <- subset(res4, padj<0.1 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)); nrow(select)
selected.expr <- heatmap.expr[select$gene,]

row.names(sTable.sub) <- sTable.sub$Sample_ID
sTable.sub <- sTable.sub[colnames(selected.expr),]
sTable.sub <- sTable.sub[order(sTable.sub$DiseaseCourse),]
selected.expr <- selected.expr[,sTable.sub$Sample_ID]
row.names(selected.expr) <- as.character(lapply(strsplit(row.names(selected.expr), "_"), "[[",2))

df <- data.frame(sTable.sub$DiseaseCourse)
row.names(df) <- sTable.sub$Sample_ID
colnames(df) <- c("Status")
#
pheatmap(selected.expr, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(165),   # RdYlBu
         fontsize_row = 8, fontsize_col = 4, scale = "row", breaks=seq(-4,4, by=0.05), 
         cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, 
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "ward.D2", 
         cutree_rows=1, cutree_cols=1, main="CD14: (CIS+RR) vs PP (FDR < 0.1 & (log2FC > 0.5 | log2FC < -0.5))",
         annotation_col=df, clustering_callback=callback)
#

rm(df, heatmap.expr, selected.expr)
#




##########  volcano plot  ##########
library(dplyr)
library(ggplot2)
library(ggrepel)

### CD4
vlcn_plot <- mutate(res4.cd4, sig=ifelse(res4.cd4$padj<0.1 & (log2FoldChange > 0.5 | log2FoldChange < -0.5), "Sig", "NoSig"))
row.names(vlcn_plot) <- row.names(res4.cd4$gene_name)
ggplot(vlcn_plot, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(aes(col=sig), show.legend=FALSE) + scale_color_manual(values=c("black","red")) +
  geom_text_repel(data=subset(vlcn_plot, padj<0.1), aes(label=subset(vlcn_plot, padj<0.1)$gene_name), colour="blue", segment.colour="black", size=3) +
  xlab("log2FC") + ylab("-log10(adj. p-value)") + ggtitle("Volcano plot for CD4: (CIS+RRMS) vs PPMS") +
  geom_hline(yintercept=-log10(0.1), linetype="dashed", color="orange") + 
  geom_vline(xintercept=c(0.5,-0.5), linetype="dashed", color="orange") +
  theme_bw() 
rm(vlcn_plot)



### CD8
vlcn_plot <- mutate(res4.cd8, sig=ifelse(res4.cd8$padj<0.1 & (log2FoldChange > 0.5 | log2FoldChange < -0.5), "Sig", "NoSig"))
row.names(vlcn_plot) <- row.names(res4.cd8$gene_name)
ggplot(vlcn_plot, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(aes(col=sig), show.legend=FALSE) + scale_color_manual(values=c("black","red")) +
  geom_text_repel(data=subset(vlcn_plot, padj<0.1), aes(label=subset(vlcn_plot, padj<0.1)$gene_name), colour="blue", segment.colour="black", size=3) +
  xlab("log2FC") + ylab("-log10(adj. p-value)") + ggtitle("Volcano plot for CD8: (CIS+RRMS) vs PPMS") +
  geom_hline(yintercept=-log10(0.1), linetype="dashed", color="orange") + xlim(c(-9, 9)) +
  geom_vline(xintercept=c(0.5,-0.5), linetype="dashed", color="orange") +
  theme_bw() 
rm(vlcn_plot)



### CD14
vlcn_plot <- mutate(res4.cd14, sig=ifelse(res4.cd14$padj<0.1 & (log2FoldChange > 0.5 | log2FoldChange < -0.5), "Sig", "NoSig"))
row.names(vlcn_plot) <- row.names(res4.cd14$gene_name)
ggplot(vlcn_plot, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(aes(col=sig), show.legend=FALSE) + scale_color_manual(values=c("black","red")) +
  geom_text_repel(data=subset(vlcn_plot, padj<0.1), aes(label=subset(vlcn_plot, padj<0.1)$gene_name), colour="blue", segment.colour="black", size=3) +
  xlab("log2FC") + ylab("-log10(adj. p-value)") + ggtitle("Volcano plot for CD14: (CIS+RRMS) vs PPMS") +
  geom_hline(yintercept=-log10(0.1), linetype="dashed", color="orange") + 
  geom_vline(xintercept=c(0.5,-0.5), linetype="dashed", color="orange") + xlim(c(-8, 12)) +
  theme_bw() 
rm(vlcn_plot)
#





##########  CD4 vs CD8 vs CD14
##########  comparison among different cell types  ##########
library("VennDiagram")
library("RColorBrewer")
library("gplots")

# genes (l2fc > 1 or < -1) and (adjusted p-value < 0.05)
venn <- list("CD4"=subset(res4.cd4, (log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.05)$gene, 
             "CD8"=subset(res4.cd8, (log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.05)$gene, 
             "CD14"=subset(res4.cd14, (log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.05)$gene)
grid.newpage()
venn.plot <- venn.diagram(venn, NULL, resolution=600, fill=RColorBrewer::brewer.pal(3, "Dark2"), col="transparent", alpha=c(0.5,0.5,0.5), 
                          main.cex = 1.5, sub.cex=3, cat.cex=1.5,
                          category=c("CD4","CD8","CD14"), main="")
grid.draw(venn.plot)

### extract list from venn diagram
tmp.list <- venn(venn, show.plot=FALSE)
str(tmp.list)
tmp.inters <- attr(tmp.list,"intersections")
lapply(tmp.inters, head, n=50)
#ven.list1 <- data.frame("CD4"=tmp.inters$CD14)
#ven.list2 <- data.frame("CD8"=tmp.inters$CD8)

#
do.call(rbind, strsplit(as.character(ven.list1[,1]), '_'))
#




##########  Plotting results  ##########
# CD4
cds <- cds.cd4
sTable.sub <- sTable.cd4
res1 <- res1.cd4
res1$comparison <- "CIS-PP"
res2 <- res2.cd4
res2$comparison <- "CIS-RR"
res3 <- res3.cd4
res3$comparison <- "PP-RR"
res4 <- res4.cd4
res4$comparison <- "(CIS+RR)-PP"
res <- rbind(subset(res1, padj<0.05), subset(res2, padj<0.05), subset(res3, padj<0.05), subset(res4, padj<0.05))
# CD8
cds <- cds.cd8
sTable.sub <- sTable.cd8
res1 <- res1.cd8
res1$comparison <- "CIS-PP"
res2 <- res2.cd8
res2$comparison <- "CIS-RR"
res3 <- res3.cd8
res3$comparison <- "PP-RR"
res4 <- res4.cd8
res4$comparison <- "(CIS+RR)-PP"
res <- rbind(subset(res1, padj<0.05), subset(res2, padj<0.05), subset(res3, padj<0.05), subset(res4, padj<0.05))
# CD14
cds <- cds.cd14
sTable.sub <- sTable.cd14
res1 <- res1.cd14
res1$comparison <- "CIS-PP"
res2 <- res2.cd14
res2$comparison <- "CIS-RR"
res3 <- res3.cd14
res3$comparison <- "PP-RR"
res4 <- res4.cd14
res4$comparison <- "(CIS+RR)-PP"
res <- rbind(subset(res1, padj<0.05), subset(res2, padj<0.05), subset(res3, padj<0.05), subset(res4, padj<0.05))
#
head(res3[order(res3$padj),])
# res1: CIS_TreatmentNaive vs PP_TreatmentNaive
# res2: CIS_TreatmentNaive vs RR_TreatmentNaive
# res3: PP_TreatmentNaive vs RR_TreatmentNaive


#### custom plot using ggplot2 (selected gene)
library(ggplot2)
library(ggbeeswarm)
setwd("~/epic.neb/analysis/results/gene-level/MS_disease_subset_updated/plots")
#res <- res[order(res$padj),]
#select <- subset(res, (log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.05); nrow(select)

#res2 <- res2[order(res2$padj),]
#select <- subset(res2, padj < 0.05); nrow(select)

select <- res[!duplicated(res$gene),]

# CD20
gene <- "ENSG00000156738.17_MS4A1"
#
gene <- subset(select, comparison=="CIS-PP" | comparison=="PP-RR")$gene[2]

#
gene <- res$gene[1]
genename <- strsplit(gene,"_")[[1]][2]; genename
### using normalized count
geneCounts <- plotCounts(cds, gene = gene, intgroup = c("DiseaseCourse2"),
                         returnData = TRUE)
#geneCounts$count <- geneCounts$count + 1
ggplot(geneCounts, aes(x = DiseaseCourse, y = count, color = DiseaseCourse)) + geom_boxplot() + 
  scale_y_log10() +  geom_beeswarm(cex = 1, groupOnX = TRUE) + 
  ggtitle(paste0(genename, " (adj. p-val = ",round(res[res$gene==gene,"padj"],digit=6),")")) + xlab("") + 
  ylab("log10(DESeq2 normalized count)") + theme_bw()
#
ggsave(paste0("DiseaseCourse_CD8_boxPlot_",genename,".png"), plot = last_plot(), device="png",
       scale = 1, width = 4.66, height = 5.34, units = "in")

rm(gene, genename)


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




#####  loop ploting for single gene boxplot #####
library(ggplot2)
library(ggbeeswarm)
#setwd("~/epic.neb/analysis/results/gene-level/MS_disease_subset_updated/plots")
setwd("~/epic.neb/analysis/results/gene-level/MS_disease_subset_updated/plots/boxplot")

# CD4
cds <- cds.cd4
cell <- "CD4"
sTable.sub <- sTable.cd4
res1 <- res1.cd4     # CIS-PP
res2 <- res2.cd4     # CIS-RR
res3 <- res3.cd4     # PP-RR
res4 <- res4.cd4     # (CIS+RR) vs PP
res <- rbind(subset(res1, padj<0.1), subset(res2, padj<0.1), subset(res3, padj<0.1), subset(res4, padj<0.1))

# CD8
cds <- cds.cd8
cell <- "CD8"
sTable.sub <- sTable.cd8
res1 <- res1.cd8     # CIS-PP
res2 <- res2.cd8     # CIS-RR
res3 <- res3.cd8     # PP-RR
res4 <- res4.cd8     # (CIS+RR) vs PP
res <- rbind(subset(res1, padj<0.1), subset(res2, padj<0.1), subset(res3, padj<0.1), subset(res4, padj<0.1))

# CD14
cds <- cds.cd14
cell <- "CD14"
sTable.sub <- sTable.cd14
res1 <- res1.cd14     # CIS-PP
res2 <- res2.cd14     # CIS-RR
res3 <- res3.cd14     # PP-RR
res4 <- res4.cd14     # (CIS+RR) vs PP
res <- rbind(subset(res1, padj<0.1), subset(res2, padj<0.1), subset(res3, padj<0.1), subset(res4, padj<0.1))

##### repeat
select <- unique(res$gene); length(unique(select))
#
for (i in 1:length(select)){
  gene <- select[i]
  genename <- strsplit(gene,"_")[[1]][2]; genename
  geneCounts <- plotCounts(cds, gene = gene, intgroup = c("DiseaseCourse"),
                           returnData = TRUE)
  geneCounts$DiseaseCourse <- factor(geneCounts$DiseaseCourse, levels=c("CIS","RR","PP"))
  ggplot(geneCounts, aes(x = DiseaseCourse, y = count, color = DiseaseCourse)) + geom_boxplot() + 
    scale_y_log10() +  geom_beeswarm(cex = 1.5) + 
    ggtitle(paste0(genename)) + xlab("") + 
    ylab("log10(DESeq2 normalized count)") + theme_bw(base_size=10)
  
  ggsave(paste0(cell,"_singlePlot_",i,"_",genename,".png"), plot = last_plot(), device="png",
         scale = 1, width = 4.66, height = 5.34, units = "in")
}
###




##########  Annotating and exporting results  ##########
# CD4
cds <- cds.cd4
sTable.sub <- sTable.cd4
res1 <- res1.cd4
res1$comparison <- "CIS-PP"
res2 <- res2.cd4
res2$comparison <- "CIS-RR"
res3 <- res3.cd4
res3$comparison <- "PP-RR"
res4 <- res4.cd4
res4$comparison <- "(CIS+RR)-PP"
res <- rbind(subset(res1, padj<0.05), subset(res2, padj<0.05), subset(res3, padj<0.05), subset(res4, padj<0.05))
# CD8
cds <- cds.cd8
sTable.sub <- sTable.cd8
res1 <- res1.cd8
res1$comparison <- "CIS-PP"
res2 <- res2.cd8
res2$comparison <- "CIS-RR"
res3 <- res3.cd8
res3$comparison <- "PP-RR"
res4 <- res4.cd8
res4$comparison <- "(CIS+RR)-PP"
res <- rbind(subset(res1, padj<0.05), subset(res2, padj<0.05), subset(res3, padj<0.05), subset(res4, padj<0.05))
# CD14
cds <- cds.cd14
sTable.sub <- sTable.cd14
res1 <- res1.cd14
res1$comparison <- "CIS-PP"
res2 <- res2.cd14
res2$comparison <- "CIS-RR"
res3 <- res3.cd14
res3$comparison <- "PP-RR"
res4 <- res4.cd14
res4$comparison <- "(CIS+RR)-PP"
res <- rbind(subset(res1, padj<0.05), subset(res2, padj<0.05), subset(res3, padj<0.05), subset(res4, padj<0.05))
#


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



#####
rm(cds, res, res.sig)
sessionInfo()
#####
