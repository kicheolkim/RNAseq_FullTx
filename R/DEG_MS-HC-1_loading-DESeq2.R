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


##### loading previously saved metadata
load(file="~/epic.neb/analysis/R_neb/RData/infoTable_updated.RData")



########## ====================================================================================================================3
setwd("~/epic.neb/analysis/results/gene-level/MS_all-HC_updated")


### Sample Table for baseline
sTable <- merge(libInfo, metaInfo_base, by.x="HCVB_ID", by.y="HCVB_ID")
sTable$Subset <- factor(sTable$Subset)
sTable$prep <- factor(sTable$prep)
sTable$Lane <- factor(sTable$Lane)
sTable$DiseaseStatus <- factor(sTable$DiseaseStatus)
sTable$Sex <- factor(sTable$Sex)
sTable$Last_Known_Treat_Stat <- factor(sTable$Last_Known_Treat_Stat)
sTable$DiseaseCourse <- factor(sTable$DiseaseCourse)

sTable$AgeAtExamGrp <- cut(sTable$AgeAtExam, breaks=c(10, 30, 40, 50, 60, 90), labels = c("10_30","30_40","40_50","50_60","60_90"))
sTable[is.na(sTable$AgeAtExamGrp),]$AgeAtExamGrp <- "30_40"  # missing age considered as average age

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





#####
sessionInfo()
#####
