# -----
# title: "DEG analysis (DESeq2) - MS-all patients (treated+untreated) vs Healthy controls"
# -----
##########=====  DEG analysis using DESeq2  =====##########
library(DESeq2)
library(tximport)
library(readr)
library(RColorBrewer)
library(pheatmap)


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



#####
rm(cds, res, res.sig)
sessionInfo()
#####
