# -----
# title: "DEG analysis (DESeq2) - MS-all patients (treated+untreated) vs Healthy controls"
# -----
##########=====  DEG analysis using DESeq2  =====##########
library(DESeq2)
library(tximport)
library(readr)
library(RColorBrewer)
library(pheatmap)


##### loading metadata - from updated file
libInfo <- read_csv("~/epic.neb/analysis/info_table/library_info.csv")
metaInfo <- read_csv("~/epic.neb/analysis/info_table/EPIC_HCvB_metadata_updated.csv",
                     col_types = cols(DMTsAtVisit = col_factor(levels = c("Copaxone", "Gilenya", "Rebif", "Aubagio", "Tecfidera", "Tysabri", "Rituximab", "none")), 
                                      DiseaseCourse = col_factor(levels = c("CIS", "PP", "RR", "RIS", "SP", "Healthy","Unknown")), 
                                      DiseaseStatus = col_factor(levels = c("MS", "Healthy", "NotMS", "Unknown")), 
                                      Last_Known_Treat_Stat = col_factor(levels = c("Treated", "TreatmentNaive", "Healthy", "NA")), 
                                      Sex = col_factor(levels = c("F", "M", "NA"))))
metaInfo_base <- read_csv("~/epic.neb/analysis/info_table/EPIC_HCvB_metadata_baseline_updated.csv",
                          col_types = cols(DMTsAtVisit = col_factor(levels = c("Copaxone", "Gilenya", "Rebif", "Aubagio", "Tecfidera", "Tysabri", "Rituximab", "none")), 
                                           DiseaseCourse = col_factor(levels = c("CIS", "PP", "RR", "RIS", "SP", "Healthy","Unknown")), 
                                           DiseaseStatus = col_factor(levels = c("MS", "Healthy", "NotMS", "Unknown")), 
                                           Last_Known_Treat_Stat = col_factor(levels = c("Treated", "TreatmentNaive", "Healthy", "NA")), 
                                           Sex = col_factor(levels = c("F", "M", "NA"))))

sTable <- merge(libInfo, metaInfo, by.x="HCVB_ID", by.y="HCVB_ID")
factor_cols <- c("Subset","prep","Lane","DiseaseStatus","Sex","Last_Known_Treat_Stat","DiseaseCourse")
sTable[factor_cols] <- lapply(sTable[factor_cols], factor)


save(libInfo, metaInfo, metaInfo_base, sTable, file="~/epic.neb/analysis/R_neb/RData/infoTable_updated.RData")


##### loading metadata
load(file="~/epic.neb/analysis/R_neb/RData/infoTable_updated.RData")





########## ====================================================================================================================3
setwd("~/epic.neb/analysis/results/gene-level/MS_all-HC_updated")


### Sample Table for baseline
sTable <- merge(libInfo, metaInfo_base, by.x="HCVB_ID", by.y="HCVB_ID")

sTable$AgeAtExamGrp <- cut(sTable$AgeAtExam, breaks=c(10, 30, 40, 50, 60, 90), labels = c("10_30","30_40","40_50","50_60","60_90"))
sTable[is.na(sTable$AgeAtExamGrp),]$AgeAtExamGrp <- "30_40"  # missing age considered as average age

factor_cols <- c("Subset","prep","Lane","DiseaseStatus","Sex","Last_Known_Treat_Stat","DiseaseCourse","AgeAtExamGrp")
sTable[factor_cols] <- lapply(sTable[factor_cols], factor)


### remove outlier samples
# 49013b-CD4, 75216a_CD4 (CD14 subset), 80517a-CD14, 80617a-CD14
sTable <- sTable[sTable$Sample_ID != "49013b-CD4",]
sTable <- sTable[sTable$Sample_ID != "75216a-CD4",]
sTable <- sTable[sTable$Sample_ID != "80517a-CD14",]
sTable <- sTable[sTable$Sample_ID != "80617a-CD14",]

write.csv(sTable, "sampleTable.csv")
# Memo: total 518 baseline samples -> 3 samples removed (75216a-CD4, 80517a-CD14, 80617a-CD14)

geneList <- read_csv("~/epic.neb/analysis/info_table/geneList.csv", col_types = cols(X1 = col_skip()))

save(libInfo, metaInfo, metaInfo_base, sTable, geneList, file="~/epic.neb/analysis/R_neb/RData/infoTable_filtered_updated.RData")


########## ====================================================================================================================3
########## ====================================================================================================================3
# reference for model matrix (consider treatment status) => https://support.bioconductor.org/p/64844/#64861

load(file="~/epic.neb/analysis/R_neb/RData/infoTable_filtered_updated.RData")
setwd("~/epic.neb/analysis/results/gene-level/MS_all-HC_updated")


##### CD4: loading selected samples from RSEM #####
# subset metadata
sTable.cd4 <- subset(sTable, Subset=="CD4")
sTable.cd4$Subset <- droplevels(sTable.cd4$Subset)
sTable.cd4$Last_Known_Treat_Stat <- droplevels(sTable.cd4$Last_Known_Treat_Stat)
sTable.cd4$DiseaseStatus <- droplevels(sTable.cd4$DiseaseStatus)

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

### input and run DESeq2 for subset - with Age at exam as a covariate
rsem.gene.sub$length[rsem.gene.sub$length == 0] <- 1
cds.sub <- DESeqDataSetFromTximport(rsem.gene.sub, sTable.sub, ~Last_Known_Treat_Stat+Sex+AgeAtExamGrp-1)   # 1 NA -> add average age (36)
cds.sub <- cds.sub[ rowSums(counts(cds.sub)) > 1, ]
cds.sub <- DESeq(cds.sub)

cds.cd4 <- cds.sub
vst.cd4 <- vst(cds.cd4)
sTable.cd4 <- sTable.sub

rm(cds.sub, sTable.sub, rsem.gene.sub)



##### CD8: loading selected samples from RSEM #####
# subset metadata
sTable.cd8 <- subset(sTable, Subset=="CD8")
sTable.cd8$Subset <- droplevels(sTable.cd8$Subset)
sTable.cd8$Last_Known_Treat_Stat <- droplevels(sTable.cd8$Last_Known_Treat_Stat)
sTable.cd8$DiseaseStatus <- droplevels(sTable.cd8$DiseaseStatus)

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
cds.sub <- DESeqDataSetFromTximport(rsem.gene.sub, sTable.sub, ~Last_Known_Treat_Stat+Sex+AgeAtExamGrp-1)   # 1 NA -> add average age (36)
cds.sub <- cds.sub[ rowSums(counts(cds.sub)) > 1, ]
cds.sub <- DESeq(cds.sub)

cds.cd8 <- cds.sub
vst.cd8 <- vst(cds.cd8)
sTable.cd8 <- sTable.sub

rm(cds.sub, sTable.sub, rsem.gene.sub)



##### CD14: loading selected samples from RSEM #####
# subset metadata
sTable.cd14 <- subset(sTable, Subset=="CD14")
sTable.cd14$Subset <- droplevels(sTable.cd14$Subset)
sTable.cd14$Last_Known_Treat_Stat <- droplevels(sTable.cd14$Last_Known_Treat_Stat)
sTable.cd14$DiseaseStatus <- droplevels(sTable.cd14$DiseaseStatus)

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
cds.sub <- DESeqDataSetFromTximport(rsem.gene.sub, sTable.sub, ~Last_Known_Treat_Stat+Sex+AgeAtExamGrp-1)   # 1 NA -> add average age (36)
cds.sub <- cds.sub[ rowSums(counts(cds.sub)) > 1, ]
cds.sub <- DESeq(cds.sub)

cds.cd14 <- cds.sub
vst.cd14 <- vst(cds.cd14)
sTable.cd14 <- sTable.sub

rm(cds.sub, sTable.sub, rsem.gene.sub)


### SAVE all
save.image("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-AllCells_MSvsHC_withAge_updated.RData")




#################### ===========================================================================33
#################### 
load("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-AllCells_MSvsHC_withAge_updated.RData")
setwd("~/epic.neb/analysis/results/gene-level/MS_all-HC_updated")
#################### 
#################### ===========================================================================33

##### for Table summary
library(tableone)

# CD4
sub <- subset(metaInfo_base, Last_Known_Treat_Stat=="TreatmentNaive")
CreateTableOne(vars=c("Sex","AgeAtExam","DiseaseCourse2","DiseaseDuration","EDSS"), data=sub)
min(sub$AgeAtExam, na.rm=TRUE); max(sub$AgeAtExam, na.rm=TRUE)
median(sub$DiseaseDuration, na.rm=TRUE); min(sub$DiseaseDuration, na.rm=TRUE); max(sub$DiseaseDuration, na.rm=TRUE)
median(sub$EDSS, na.rm=TRUE); min(sub$EDSS, na.rm=TRUE); max(sub$EDSS, na.rm=TRUE)

# CD8
sub <- subset(metaInfo_base, Last_Known_Treat_Stat=="Treated")
CreateTableOne(vars=c("Sex","AgeAtExam","DiseaseCourse2","DiseaseDuration","EDSS"), data=sub)
min(sub$AgeAtExam, na.rm=TRUE); max(sub$AgeAtExam, na.rm=TRUE)
median(sub$DiseaseDuration, na.rm=TRUE); min(sub$DiseaseDuration, na.rm=TRUE); max(sub$DiseaseDuration, na.rm=TRUE)
median(sub$EDSS, na.rm=TRUE); min(sub$EDSS, na.rm=TRUE); max(sub$EDSS, na.rm=TRUE)

# CD14
sub <- subset(metaInfo_base, Last_Known_Treat_Stat=="Healthy")
CreateTableOne(vars=c("Sex","AgeAtExam"), data=sub)
min(sub$AgeAtExam, na.rm=TRUE); max(sub$AgeAtExam, na.rm=TRUE)

# disease subset
sub <- subset(sTable, Subset=="CD14" & Last_Known_Treat_Stat=="TreatmentNaive")
sub <- subset(sTable, Subset=="CD14" & Last_Known_Treat_Stat=="Treated")
CreateTableOne(vars=c("DiseaseCourse"), data=sub)



##### PCA plot
plotPCA(vst.cd4, intgroup = c("Sex","Last_Known_Treat_Stat"), ntop = 5000)
plotPCA(vst.cd4, intgroup = c("AgeAtExamGrp"), ntop = 10000)
plotPCA(vst.cd8, intgroup = c("Sex","Last_Known_Treat_Stat"), ntop = 5000)
plotPCA(vst.cd14, intgroup = c("Sex","Last_Known_Treat_Stat"), ntop = 5000)

# Sex difference is less in CD14 than CD4 and CD8



########## Comparison results ##########
# DESeq2 objects: cds.cd4, cds.cd8, cds.cd14
# metadata tables: sTable.cd4, sTable.cd8, sTable.cd14
## 
geneList <- read_csv("~/epic.neb/analysis/info_table/geneList.csv", col_types = cols(X1 = col_skip()))

select <- geneList
select$alias <- mapIds(org.Hs.eg.db,
                       keys=substr(select$gene_id, 1,15),
                       column="ALIAS",
                       keytype="ENSEMBL",
                       multiVals="CharacterList")
select$alias <- do.call(rbind, lapply(select$alias,paste0,collapse="|"))
select$alias <- as.character(select$alias)
select$alias <- paste0(select$gene_name,"|",select$alias)

select$GeneFullName <- mapIds(org.Hs.eg.db,
                              keys=substr(select$gene_id, 1,15),
                              column="GENENAME",
                              keytype="ENSEMBL",
                              multiVals="CharacterList")
select$GeneFullName <- do.call(rbind, lapply(select$GeneFullName,paste0,collapse="|"))
select$GeneFullName <- as.character(select$GeneFullName)
geneList <- select; rm(select)





########## CD4: DEG results  ##########
cds <- cds.cd4
vst <- vst.cd4

res1 <- results(cds, alpha=0.05, contrast=c("Last_Known_Treat_Stat","Treated","TreatmentNaive"))
summary(res1)   ### Treated vs Untreat
res2 <- results(cds, alpha=0.05, contrast=c("Last_Known_Treat_Stat","Treated","Healthy"))
summary(res2)   ### Treated vs Healthy
res3 <- results(cds, alpha=0.05, contrast=c("Last_Known_Treat_Stat","TreatmentNaive","Healthy"))
summary(res3)   ### Untreated vs Healthy
res <- results(cds, alpha=0.05, 
               contrast=list(c("Last_Known_Treat_StatTreated","Last_Known_Treat_StatTreatmentNaive"),"Last_Known_Treat_StatHealthy"),
               listValues = c(1/2, -1))
summary(res)   ### (Treat + Untreat) vs Healthy

## find gene related to Sex
dds <- DESeq(cds, test="LRT", reduced=~Last_Known_Treat_Stat+AgeAtExamGrp-1)
res.sex <- results(dds, alpha=0.05)
summary(res.sex)

## find gene related to Age
dds <- DESeq(cds, test="LRT", reduced=~Last_Known_Treat_Stat+Sex-1)
res.age <- results(dds, alpha=0.05)
summary(res.age)

## change order by adjusted p-value
res <- as.data.frame(res)
res <- res[order(res$padj),]
res1 <- as.data.frame(res1)
res1 <- res1[order(res1$padj),]
res2 <- as.data.frame(res2)
res2 <- res2[order(res2$padj),]
res3 <- as.data.frame(res3)
res3 <- res3[order(res3$padj),]
res.sex <- as.data.frame(res.sex)
res.sex <- res.sex[order(res.sex$padj),]
res.age <- as.data.frame(res.age)
res.age <- res.age[order(res.age$padj),]

res.sex$gene <- row.names(res.sex)
res.sex <- merge(res.sex, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
write.csv(res.sex, file="MS-HC_CD4_result-covar_Sex.csv")
res.age$gene <- row.names(res.age)
res.age <- merge(res.age, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
write.csv(res.age, file="MS-HC_CD4_result-covar_Age.csv")

res$gene <- row.names(res)
res <- merge(res, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
write.csv(res, file="MS-HC_CD4_result-Trt+Untrt.vs.HC.csv")

res1$gene <- row.names(res1)  ### Treated vs Untreat
res1 <- merge(res1, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res1.cd4 <- res1
write.csv(res1, file="MS-HC_CD4_result-Trt.vs.Untrt.csv")
res2$gene <- row.names(res2)  ### Treated vs Healthy
res2 <- merge(res2, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res2.cd4 <- res2
write.csv(res2, file="MS-HC_CD4_result-Trt.vs.HC.csv")
res3$gene <- row.names(res3)  ### Untreated vs Healthy
res3 <- merge(res3, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res3.cd4 <- res3
write.csv(res3, file="MS-HC_CD4_result-Untrt.vs.HC.csv")


res_coding <- subset(res3, gene_type=="protein_coding"|gene_type=="IG_C_gene"|gene_type=="IG_D_gene"|
                       gene_type=="IG_J_gene"|gene_type=="IG_LV_gene"|gene_type=="IG_V_gene"|
                       gene_type=="TR_C_gene"|gene_type=="TR_J_gene"|gene_type=="TR_V_gene"|gene_type=="TR_D_gene")
res_other <- subset(res3, !(gene_type=="protein_coding"|gene_type=="IG_C_gene"|gene_type=="IG_D_gene"|
                             gene_type=="IG_J_gene"|gene_type=="IG_LV_gene"|gene_type=="IG_V_gene"|
                             gene_type=="TR_C_gene"|gene_type=="TR_J_gene"|gene_type=="TR_V_gene"|gene_type=="TR_D_gene"))
write.csv(res_coding, file="MS-HC_CD4_result-Untrt.vs.HC-coding.csv")
write.csv(res_other, file="MS-HC_CD4_result-Untrt.vs.HC-noncoding.csv")

dds.cd4 <- dds
res.cd4 <- res
rm(cds, dds, vst, res, res.sex, res.age, res1, res2, res3, res_coding, res_other)





########## CD8: DEG results  ##########
cds <- cds.cd8
vst <- vst.cd8

res1 <- results(cds, alpha=0.05, contrast=c("Last_Known_Treat_Stat","Treated","TreatmentNaive"))
summary(res1)   ### Treated vs Untreat
res2 <- results(cds, alpha=0.05, contrast=c("Last_Known_Treat_Stat","Treated","Healthy"))
summary(res2)   ### Treated vs Healthy
res3 <- results(cds, alpha=0.05, contrast=c("Last_Known_Treat_Stat","TreatmentNaive","Healthy"))
summary(res3)   ### Untreated vs Healthy
res <- results(cds, alpha=0.05, 
               contrast=list(c("Last_Known_Treat_StatTreated","Last_Known_Treat_StatTreatmentNaive"),"Last_Known_Treat_StatHealthy"),
               listValues = c(1/2, -1))
summary(res)   ### (Treat + Untreat) vs Healthy

## find gene related to Sex
dds <- DESeq(cds, test="LRT", reduced=~Last_Known_Treat_Stat+AgeAtExamGrp-1)
res.sex <- results(dds, alpha=0.05)
summary(res.sex)

## find gene related to Age
dds <- DESeq(cds, test="LRT", reduced=~Last_Known_Treat_Stat+Sex-1)
res.age <- results(dds, alpha=0.05)
summary(res.age)

## change order by adjusted p-value
res <- as.data.frame(res)
res <- res[order(res$padj),]
res1 <- as.data.frame(res1)
res1 <- res1[order(res1$padj),]
res2 <- as.data.frame(res2)
res2 <- res2[order(res2$padj),]
res3 <- as.data.frame(res3)
res3 <- res3[order(res3$padj),]
res.sex <- as.data.frame(res.sex)
res.sex <- res.sex[order(res.sex$padj),]
res.age <- as.data.frame(res.age)
res.age <- res.age[order(res.age$padj),]

res.sex$gene <- row.names(res.sex)
res.sex <- merge(res.sex, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
write.csv(res.sex, file="MS-HC_CD8_result-covar_Sex.csv")
res.age$gene <- row.names(res.age)
res.age <- merge(res.age, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
write.csv(res.age, file="MS-HC_CD8_result-covar_Age.csv")
res$gene <- row.names(res)
res <- merge(res, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
write.csv(res, file="MS-HC_CD8_result-Trt+Untrt.vs.HC.csv")

res1$gene <- row.names(res1)  ### Treated vs Untreat
res1 <- merge(res1, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res1.cd8 <- res1
write.csv(res1, file="MS-HC_CD8_result-Trt.vs.Untrt.csv")
res2$gene <- row.names(res2)  ### Treated vs Healthy
res2 <- merge(res2, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res2.cd8 <- res2
write.csv(res2, file="MS-HC_CD8_result-Trt.vs.HC.csv")
res3$gene <- row.names(res3)  ### Untreated vs Healthy
res3 <- merge(res3, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res3.cd8 <- res3
write.csv(res3, file="MS-HC_CD8_result-Untrt.vs.HC.csv")


res_coding <- subset(res3, gene_type=="protein_coding"|gene_type=="IG_C_gene"|gene_type=="IG_D_gene"|
                       gene_type=="IG_J_gene"|gene_type=="IG_LV_gene"|gene_type=="IG_V_gene"|
                       gene_type=="TR_C_gene"|gene_type=="TR_J_gene"|gene_type=="TR_V_gene"|gene_type=="TR_D_gene")
res_other <- subset(res3, !(gene_type=="protein_coding"|gene_type=="IG_C_gene"|gene_type=="IG_D_gene"|
                             gene_type=="IG_J_gene"|gene_type=="IG_LV_gene"|gene_type=="IG_V_gene"|
                             gene_type=="TR_C_gene"|gene_type=="TR_J_gene"|gene_type=="TR_V_gene"|gene_type=="TR_D_gene"))
write.csv(res_coding, file="MS-HC_CD8_result-Untrt.vs.HC-coding.csv")
write.csv(res_other, file="MS-HC_CD8_result-Untrt.vs.HC-noncoding.csv")

dds.cd8 <- dds
res.cd8 <- res
rm(cds, dds, vst, res, res.sex, res.age, res1, res2, res3, res_coding, res_other)





########## CD14: DEG results  ##########
cds <- cds.cd14
vst <- vst.cd14

res1 <- results(cds, alpha=0.05, contrast=c("Last_Known_Treat_Stat","Treated","TreatmentNaive"))
summary(res1)   ### Treated vs Untreat
res2 <- results(cds, alpha=0.05, contrast=c("Last_Known_Treat_Stat","Treated","Healthy"))
summary(res2)   ### Treated vs Healthy
res3 <- results(cds, alpha=0.05, contrast=c("Last_Known_Treat_Stat","TreatmentNaive","Healthy"))
summary(res3)   ### Untreated vs Healthy
res <- results(cds, alpha=0.05, 
               contrast=list(c("Last_Known_Treat_StatTreated","Last_Known_Treat_StatTreatmentNaive"),"Last_Known_Treat_StatHealthy"),
               listValues = c(1/2, -1))
summary(res)   ### (Treat + Untreat) vs Healthy

## find gene related to Sex
dds <- DESeq(cds, test="LRT", reduced=~Last_Known_Treat_Stat+AgeAtExamGrp-1)
res.sex <- results(dds, alpha=0.05)
summary(res.sex)

## find gene related to Age
dds <- DESeq(cds, test="LRT", reduced=~Last_Known_Treat_Stat+Sex-1)
res.age <- results(dds, alpha=0.05)
summary(res.age)


## change order by adjusted p-value
res <- as.data.frame(res)
res <- res[order(res$padj),]
res1 <- as.data.frame(res1)
res1 <- res1[order(res1$padj),]
res2 <- as.data.frame(res2)
res2 <- res2[order(res2$padj),]
res3 <- as.data.frame(res3)
res3 <- res3[order(res3$padj),]
res.sex <- as.data.frame(res.sex)
res.sex <- res.sex[order(res.sex$padj),]
res.age <- as.data.frame(res.age)
res.age <- res.age[order(res.age$padj),]

res.sex$gene <- row.names(res.sex)
res.sex <- merge(res.sex, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
write.csv(res.sex, file="MS-HC_CD14_result-covar_Sex.csv")
res.age$gene <- row.names(res.age)
res.age <- merge(res.age, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
write.csv(res.age, file="MS-HC_CD14_result-covar_Age.csv")
res$gene <- row.names(res)
res <- merge(res, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
write.csv(res, file="MS-HC_CD14_result-Trt+Untrt.vs.HC.csv")

res1$gene <- row.names(res1)  ### Treated vs Untreat
res1 <- merge(res1, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res1.cd14 <- res1
write.csv(res1, file="MS-HC_CD14_result-Trt.vs.Untrt.csv")
res2$gene <- row.names(res2)  ### Treated vs Healthy
res2 <- merge(res2, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res2.cd14 <- res2
write.csv(res2, file="MS-HC_CD14_result-Trt.vs.HC.csv")
res3$gene <- row.names(res3)  ### Untreated vs Healthy
res3 <- merge(res3, geneList[,c("gene","gene_type","gene_name","seqnames","start","end","width","strand","gene_id","GeneFullName","alias")], by.x="gene", by.y="gene")
res3.cd14 <- res3
write.csv(res3, file="MS-HC_CD14_result-Untrt.vs.HC.csv")


res_coding <- subset(res3, gene_type=="protein_coding"|gene_type=="IG_C_gene"|gene_type=="IG_D_gene"|
                       gene_type=="IG_J_gene"|gene_type=="IG_LV_gene"|gene_type=="IG_V_gene"|
                       gene_type=="TR_C_gene"|gene_type=="TR_J_gene"|gene_type=="TR_V_gene"|gene_type=="TR_D_gene")
res_other <- subset(res3, !(gene_type=="protein_coding"|gene_type=="IG_C_gene"|gene_type=="IG_D_gene"|
                             gene_type=="IG_J_gene"|gene_type=="IG_LV_gene"|gene_type=="IG_V_gene"|
                             gene_type=="TR_C_gene"|gene_type=="TR_J_gene"|gene_type=="TR_V_gene"|gene_type=="TR_D_gene"))
write.csv(res_coding, file="MS-HC_CD14_result-Untrt.vs.HC-coding.csv")
write.csv(res_other, file="MS-HC_CD14_result-Untrt.vs.HC-noncoding.csv")

dds.cd14 <- dds
res.cd14 <- res
rm(cds, dds, vst, res, res.sex, res.age, res1, res2, res3, res_coding, res_other)



### SAVE all  ########## 
save.image("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-AllCells_MSvsHC_withAge_updated.RData")






########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 333
########## CHECK the results ########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 333


##### load all dataset
load("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-AllCells_MSvsHC_withAge_updated.RData")


##
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


## check number of significant genes - protein coding
summary_deg <- function(result_table){
  sub <- subset(result_table, (gene_type=="protein_coding") | (gene_type=="IG_C_gene") | (gene_type=="IG_V_gene") | 
                  (gene_type=="TR_C_gene") | (gene_type=="TR_V_gene") | (gene_type=="TR_J_gene"))
  sub <- subset(sub, baseMean > 3 & padj < 0.05)
  sub$l2fc <- "down"
  sub[sub$log2FoldChange > 0,]$l2fc <- "up"
  coding <- paste0("Protein: total= ",nrow(sub)," / up-regulated= ",nrow(sub[sub$l2fc=="up",])," / down-regulated= ",nrow(sub[sub$l2fc=="down",]))
  
  
  sub <- subset(result_table, (gene_type!="protein_coding") & (gene_type!="IG_C_gene") & (gene_type!="IG_V_gene") & 
                  (gene_type!="TR_C_gene") & (gene_type!="TR_V_gene") & (gene_type!="TR_J_gene"))
  sub <- subset(sub, baseMean > 3 & padj < 0.05)
  sub$l2fc <- "down"
  sub[sub$log2FoldChange > 0,]$l2fc <- "up"
  noncoding <- paste0("ncRNA/pseudogene: total= ",nrow(sub)," / up-regulated= ",nrow(sub[sub$l2fc=="up",])," / down-regulated= ",nrow(sub[sub$l2fc=="down",]))
  
  return(list(coding, noncoding))
}


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




##########  CHECK Results  ########## ########## ########## ########## ########## ########## ##########333
##### load all dataset
load("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-AllCells_MSvsHC_withAge_updated.RData")

setwd("~/epic.neb/analysis/results/gene-level/MS_all-HC_updated/plots")


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


#table(c(res2_sub$cell, res2_sub$treat, res2_sub$type, res2_sub$l2fc))

count(res2_sub, c("cell","treat","type","l2fc"))
count(res3_sub, c("cell","treat","type","l2fc"))





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
select <- subset(res3, baseMean > 3 & padj < 0.05); nrow(select)
select <- select[order(select$padj),]
selected.expr <- heatmap.expr[select$gene,]

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




##########  volcano plot  ##########
library(dplyr)
library(ggplot2)
library(ggrepel)

### CD4
vlcn_plot <- mutate(res3.cd4, sig=ifelse(res3.cd4$padj<0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5), "Sig", "NoSig"))
vlcn_plot <- mutate(vlcn_plot, sigName=ifelse(res3.cd4$padj<0.0001 & (log2FoldChange > 0.5 | log2FoldChange < -0.5), "Sig", "NoSig"))
row.names(vlcn_plot) <- row.names(res3.cd4$gene_name)
ggplot(vlcn_plot, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(aes(col=sig), show.legend=FALSE) + scale_color_manual(values=c("black","red")) +
  geom_text_repel(data=subset(vlcn_plot, sigName=="Sig"), aes(label=subset(vlcn_plot, sigName=="Sig")$gene_name), colour="blue", segment.colour="black", size=3) +
  xlab("log2FC") + ylab("-log10(FDR-adjusted p-value)") + ggtitle("Volcano plot, CD4: MS-Untreated v.s Healthy") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="orange") + 
  geom_vline(xintercept=c(0.5,-0.5), linetype="dashed", color="orange") + xlim(-5,5) + 
  theme_bw() 
ggsave(paste0("CD4_volcano-MS_untreat-HC.pdf"), plot = last_plot(), device="pdf",
       scale = 1, width = 6, height = 5.4, units = "in", useDingbats=F)
rm(vlcn_plot)



### CD8
vlcn_plot <- mutate(res3.cd8, sig=ifelse(res3.cd8$padj<0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5), "Sig", "NoSig"))
vlcn_plot <- mutate(vlcn_plot, sigName=ifelse(res3.cd8$padj<0.02 & (log2FoldChange > 0.5 | log2FoldChange < -0.5), "Sig", "NoSig"))
row.names(vlcn_plot) <- row.names(res3.cd8$gene_name)
ggplot(vlcn_plot, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(aes(col=sig), show.legend=FALSE) + scale_color_manual(values=c("black","red")) +
  geom_text_repel(data=subset(vlcn_plot, sigName=="Sig"), 
                  aes(label=subset(vlcn_plot, sigName=="Sig")$gene_name), colour="blue", segment.colour="black", size=2.5) +
  xlab("log2FC") + ylab("-log10(FDR-adjusted p-value)") + ggtitle("Volcano plot, CD8: MS-Untreated v.s Healthy") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="orange") + 
  geom_vline(xintercept=c(0.5,-0.5), linetype="dashed", color="orange") + 
  xlim(-6,6) + ylim(0,3) + theme_bw() 
ggsave(paste0("CD8_volcano-MS_untreat-HC.pdf"), plot = last_plot(), device="pdf",
       scale = 1, width = 6, height = 5.4, units = "in", useDingbats=F)



### CD14
vlcn_plot <- mutate(res3.cd14, sig=ifelse(res3.cd14$padj<0.05 & (log2FoldChange > 0.5 | log2FoldChange < -0.5), "Sig", "NoSig"))
vlcn_plot <- mutate(vlcn_plot, sigName=ifelse(res3.cd14$padj<0.005 & (log2FoldChange > 0.5 | log2FoldChange < -0.5), "Sig", "NoSig"))
row.names(vlcn_plot) <- row.names(res3.cd14$gene_name)
ggplot(vlcn_plot, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(aes(col=sig), show.legend=FALSE) + scale_color_manual(values=c("black","red")) +
  geom_text_repel(data=subset(vlcn_plot, sigName=="Sig"), 
                  aes(label=subset(vlcn_plot, sigName=="Sig")$gene_name), colour="blue", segment.colour="black", size=3) +
  xlab("log2FC") + ylab("-log10(FDR-adjusted p-value)") + ggtitle("Volcano plot, CD14: MS-Untreated v.s Healthy") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="orange") + 
  geom_vline(xintercept=c(0.5,-0.5), linetype="dashed", color="orange") + 
  xlim(-4,4) + theme_bw() 
ggsave(paste0("CD14_volcano-MS_untreat-HC.pdf"), plot = last_plot(), device="pdf",
       scale = 1, width = 6, height = 5.4, units = "in", useDingbats=F)

#




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





############################################################33
##########  Plotting results for selected gene    ##########
############################################################33

## custom single gene boxplot using ggplot2
library(DESeq2)
library(ggplot2)
library(ggbeeswarm)

##### load all dataset
load("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-AllCells_MSvsHC_withAge_updated.RData")


#####  single gene boxplot  ##### selected genes
target <- "ENSG00000179388.8_EGR3"
target <- "ENSG00000181418.7_DDN"
target <- "ENSG00000233913.7_CTC-575D19.1"

target <- "ENSG00000158985.13_CDC42SE2"
target <- "ENSG00000144744.16_UBA3"
target <- "ENSG00000183813.6_CCR4"
sig_gene_list <- subset(res, res$gene==target)
i=1


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
    scale_y_log10() +  geom_beeswarm(cex = 2.5) + 
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
    scale_y_log10() +  geom_beeswarm(cex = 2.5) + 
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
    scale_y_log10() +  geom_beeswarm(cex = 2.5) + 
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




### count table
geneCounts$Sample_ID <- row.names(geneCounts)
geneCounts_table <- merge(geneCounts, sTable.cd14, by.x="Sample_ID", by.y="Sample_ID")
geneCounts_table$HCVB <- gsub("[a-z]", "", geneCounts_table$HCVB_ID)

HLA_allele <- read_csv("~/epic.neb/analysis/results/HLA/HLA_allele.csv")
HLA_allele_sub <- HLA_allele[HLA_allele$HCVB_ID %in% unique(geneCounts_table$HCVB_ID),]
geneCounts_table <- merge(geneCounts_table, HLA_allele_sub, by.x="HCVB_ID", by.y="HCVB_ID")
write.csv(geneCounts_table, paste0(cell,"_geneCounts_",genename,".csv"))
#

# plots to check with phenotypes
plot(geneCounts_table$DiseaseDuration, geneCounts_table$count)
plot(geneCounts_table$AgeAtExam.x, geneCounts_table$count)
ggplot(geneCounts_table, aes(x = HLA_DRB1_1, y = count, color = HLA_DRB1_1)) + geom_boxplot() + 
  scale_y_log10() +  geom_beeswarm(cex = 1.5) + ggtitle(paste0(genename)) + xlab("") + ylab("log10(DESeq2 normalized count)") + theme_bw()
ggplot(geneCounts_table, aes(x = HLA_DRB1_2, y = count, color = HLA_DRB1_2)) + geom_boxplot() + 
  scale_y_log10() +  geom_beeswarm(cex = 1.5) + ggtitle(paste0(genename)) + xlab("") + ylab("log10(DESeq2 normalized count)") + theme_bw()
ggplot(geneCounts_table, aes(x = Sex, y = count, color = Sex)) + geom_boxplot() + 
  scale_y_log10() +  geom_beeswarm(cex = 1.5) + ggtitle(paste0(genename)) + xlab("") + ylab("log10(DESeq2 normalized count)") + theme_bw()
ggplot(geneCounts_table, aes(x = AgeAtExam.x, y = count, color = DiseaseDuration)) + geom_density() + 
  scale_y_log10() +  ggtitle(paste0(genename)) + xlab("") + ylab("log10(DESeq2 normalized count)") + theme_bw()
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
##### for loop for single gene boxplot
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





#############################################################################################3
#####  check gene name from DEG with IMSGC (comparison with IMSGC genes)  #####
########################################3
##### load all dataset
load("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-AllCells_MSvsHC_withAge_updated.RData")
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
