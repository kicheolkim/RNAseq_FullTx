################################################
### FastQTL using variants called from RNA-Seq dataset ###
################################################
load(file="~/epic.neb/analysis/R_neb/RData/fastqtl_input_tables.RData")



########  make phenotype (gene/transcript expression) table for FastQTL  ########
library(readr)
library(tximport)
library(gtools)
#load(file="~/epic.neb/analysis/R_neb/RData/infoTable_filtered.RData")
#geneList <- read_csv("~/epic.neb/analysis/info_table/geneList.csv", col_types = cols(X1 = col_skip()))
load(file="~/epic.neb/analysis/R_neb/RData/infoTable_filtered_updated.RData")
MSchip_RNAseq_samples <- read_csv("~/epic.neb/qtl/GWAS_mschip/MSchip-RNAseq.samples_ed_baseline.csv")
vcf_list_all <- read_csv("~/epic.neb/analysis/all_vcf_list.csv")


#chrOrder <- c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
#chrOrder_gl <- unique(gene.expr_table_merge$seqnames)
#chrOrder_gl <- chrOrder_gl[mixedorder(chrOrder_gl)]
#chrOrder <- c(chrOrder,chrOrder_gl[26:length(chrOrder_gl)])



##### CD4: loading selected samples from RSEM #####
# subset metadata
sTable.cd4 <- subset(sTable, Subset=="CD4")
#sTable.cd4 <- subset(sTable, Subset=="CD4" & (DiseaseStatus=="MS" | DiseaseStatus=="Healthy") & !is.na(Last_Known_Treat_Stat))
sTable.cd4$Subset <- droplevels(sTable.cd4$Subset)
sTable.cd4$Last_Known_Treat_Stat <- droplevels(sTable.cd4$Last_Known_Treat_Stat)
sTable.cd4$DiseaseStatus <- droplevels(sTable.cd4$DiseaseStatus)
sTable.cd4[is.na(sTable.cd4$AgeAtExam),]$AgeAtExam <- 36
sTable.cd4$Sex <- as.character(sTable.cd4$Sex)
sTable.cd4[sTable.cd4$Sex == "M",]$Sex <- 0
sTable.cd4[sTable.cd4$Sex == "F",]$Sex <- 1

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

geneList$gene_id2 <- paste0(geneList$gene_id,"_",geneList$gene_name)
geneList <- as.data.frame(geneList)
geneList$gene_id2 <- make.unique(geneList$gene_id2)
row.names(geneList) <- geneList$gene_id2
geneList <- geneList[row.names(rsem.gene.sub$abundance),]

gene.expr_table <- as.data.frame(rsem.gene.sub$abundance)
colnames(gene.expr_table) <- colnames(rsem.gene.sub$abundance)
gene.expr_table$gene_id <- row.names(gene.expr_table)
gene.expr_table_merge <- merge(gene.expr_table, geneList[,c("seqnames","start","end","gene_id2")], by.x="gene_id", by.y="gene_id2")
expr_ncol <- ncol(gene.expr_table_merge)
gene.expr_table_merge <- gene.expr_table_merge[,c((expr_ncol-2):expr_ncol,1,2:(expr_ncol-3))]
chrOrder <- c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
chrOrder_gl <- unique(gene.expr_table_merge$seqnames)
chrOrder_gl <- chrOrder_gl[mixedorder(chrOrder_gl)]
chrOrder <- c(chrOrder,chrOrder_gl[26:length(chrOrder_gl)])
gene.expr_table_merge$seqnames <- factor(gene.expr_table_merge$seqnames, levels=chrOrder)
gene.expr_table_merge <- gene.expr_table_merge[order(gene.expr_table_merge$seqnames, gene.expr_table_merge$start),]
colnames(gene.expr_table_merge) <- c("#Chr","start","end","ID",colnames(gene.expr_table_merge)[5:length(gene.expr_table_merge)])

gene.expr_table_merge2 <- gene.expr_table_merge
#gene.expr_table_merge2 <- subset(gene.expr_table_merge, substr(gene.expr_table_merge$`#Chr`, 1, 2) != "KI")
#gene.expr_table_merge2 <- subset(gene.expr_table_merge2, substr(gene.expr_table_merge2$`#Chr`, 1, 2) != "GL")
#gene.expr_table_merge2$`#Chr` <- substr(gene.expr_table_merge2$`#Chr`, 4, 7)



### For GWAS: select samples that have GWAS data
# modify code for splicing
MSchip_RNAseq_samples_sub <- subset(MSchip_RNAseq_samples, Subset=="CD4")
row.names(MSchip_RNAseq_samples_sub) <- MSchip_RNAseq_samples_sub$Sample_ID2

expr_overlapGWAS <- gene.expr_table_merge2[, colnames(gene.expr_table_merge2) %in% MSchip_RNAseq_samples_sub$Sample_ID2]
expr_overlapGWAS <- cbind(gene.expr_table_merge2[,c(1:4)], expr_overlapGWAS)
MSchip_RNAseq_samples_sub <- MSchip_RNAseq_samples_sub[colnames(expr_overlapGWAS)[5:ncol(expr_overlapGWAS)],]
colnames(expr_overlapGWAS)[5:ncol(expr_overlapGWAS)] <- MSchip_RNAseq_samples_sub$vcfID

nrow(expr_overlapGWAS); length(unique(row.names(expr_overlapGWAS)))
gene_expr_CD4 <- expr_overlapGWAS
#write.table(expr_overlapGWAS, "~/epic.neb/qtl/pheno_expr/gene_expr_CD4.txt", sep="\t", row.names=FALSE, quote = FALSE)
rm(expr_overlapGWAS, MSchip_RNAseq_samples_sub, sTable.sub)


### For RNA-seq variants: CD4 - random select 20 patients
sTable.sub <- subset(sTable.cd4, Last_Known_Treat_Stat=="TreatmentNaive")
sTable.sub <- sTable.sub[sTable.sub$Sample_ID %in% sample(subset(sTable.cd4, Last_Known_Treat_Stat=="TreatmentNaive")$Sample_ID, 20),]
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
expr <- gene.expr_table_merge2[, colnames(gene.expr_table_merge2) %in% sTable.sub$Sample_ID]
expr <- expr[,vcf_list_sub$sample]
expr <- cbind(gene.expr_table_merge2[,c(1:4)], expr)
nrow(expr); length(unique(row.names(expr)))
#
expr <- expr_transform_new(expr)
write.table(expr, "~/epic.neb/qtl/test_20samp/CD4_test2.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(vcf_list_sub$path, file="~/epic.neb/qtl/test_20samp/CD4_test2_vcf.list", quote=FALSE, row.names=FALSE, col.names=FALSE)
#
covar <- sTable.sub[,c("Sex","AgeAtExam")]
row.names(covar) <- sTable.sub$Sample_ID
covar$id <- sTable.sub$Sample_ID
covar <- covar[vcf_list_sub$sample,]
covar <- t(covar[,c(3,1:2)])
write.table(covar, "~/epic.neb/qtl/test_20samp/CD4_test2_covar.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)




### For RNA-seq variants: CD4: select samples for eQTL test     #####
sTable.sub <- subset(sTable.cd4, Last_Known_Treat_Stat=="TreatmentNaive")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
expr <- gene.expr_table_merge2[, colnames(gene.expr_table_merge2) %in% sTable.sub$Sample_ID]
expr <- expr[,vcf_list_sub$sample]
expr <- cbind(gene.expr_table_merge2[,c(1:4)], expr)
nrow(expr); length(unique(row.names(expr)))
expr_CD4_MSuntreat <- expr
rm(expr, sTable.sub, vcf_list_sub)

sTable.sub <- subset(sTable.cd4, Last_Known_Treat_Stat=="Healthy")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
expr <- gene.expr_table_merge2[, colnames(gene.expr_table_merge2) %in% sTable.sub$Sample_ID]
expr <- expr[,vcf_list_sub$sample]
expr <- cbind(gene.expr_table_merge2[,c(1:4)], expr)
nrow(expr); length(unique(row.names(expr)))
expr_CD4_HC <- expr
rm(expr, sTable.sub, vcf_list_sub)

sTable.sub <- subset(sTable.cd4, DiseaseStatus=="MS")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
expr <- gene.expr_table_merge2[, colnames(gene.expr_table_merge2) %in% sTable.sub$Sample_ID]
expr <- expr[,vcf_list_sub$sample]
expr <- cbind(gene.expr_table_merge2[,c(1:4)], expr)
nrow(expr); length(unique(row.names(expr)))
expr_CD4_MSall <- expr
rm(expr, sTable.sub, vcf_list_sub)

rm(gene.expr_table_merge, gene.expr_table_merge2, expr_ncol, chrOrder_gl)





##### CD8: loading selected samples from RSEM #####
# subset metadata
sTable.cd8 <- subset(sTable, Subset=="CD8")
#sTable.cd8 <- subset(sTable, Subset=="CD8" & (DiseaseStatus=="MS" | DiseaseStatus=="Healthy") & !is.na(Last_Known_Treat_Stat)) # & !is.na(Last_Known_Treat_Stat)
sTable.cd8$Subset <- droplevels(sTable.cd8$Subset)
sTable.cd8$Last_Known_Treat_Stat <- droplevels(sTable.cd8$Last_Known_Treat_Stat)
sTable.cd8$DiseaseStatus <- droplevels(sTable.cd8$DiseaseStatus)
sTable.cd8[is.na(sTable.cd8$AgeAtExam),]$AgeAtExam <- 36
sTable.cd8$Sex <- as.character(sTable.cd8$Sex)
sTable.cd8[sTable.cd8$Sex == "M",]$Sex <- 0
sTable.cd8[sTable.cd8$Sex == "F",]$Sex <- 1

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

geneList$gene_id2 <- paste0(geneList$gene_id,"_",geneList$gene_name)
geneList <- as.data.frame(geneList)
geneList$gene_id2 <- make.unique(geneList$gene_id2)
row.names(geneList) <- geneList$gene_id2
geneList <- geneList[row.names(rsem.gene.sub$abundance),]

gene.expr_table <- as.data.frame(rsem.gene.sub$abundance)
colnames(gene.expr_table) <- colnames(rsem.gene.sub$abundance)
gene.expr_table$gene_id <- row.names(gene.expr_table)
gene.expr_table_merge <- merge(gene.expr_table, geneList[,c("seqnames","start","end","gene_id2")], by.x="gene_id", by.y="gene_id2")
expr_ncol <- ncol(gene.expr_table_merge)
gene.expr_table_merge <- gene.expr_table_merge[,c((expr_ncol-2):expr_ncol,1,2:(expr_ncol-3))]
chrOrder <- c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
chrOrder_gl <- unique(gene.expr_table_merge$seqnames)
chrOrder_gl <- chrOrder_gl[mixedorder(chrOrder_gl)]
chrOrder <- c(chrOrder,chrOrder_gl[26:length(chrOrder_gl)])
gene.expr_table_merge$seqnames <- factor(gene.expr_table_merge$seqnames, levels=chrOrder)
gene.expr_table_merge <- gene.expr_table_merge[order(gene.expr_table_merge$seqnames, gene.expr_table_merge$start),]
colnames(gene.expr_table_merge) <- c("#Chr","start","end","ID",colnames(gene.expr_table_merge)[5:length(gene.expr_table_merge)])

gene.expr_table_merge2 <- gene.expr_table_merge
#gene.expr_table_merge2 <- subset(gene.expr_table_merge, substr(gene.expr_table_merge$`#Chr`, 1, 2) != "KI")
#gene.expr_table_merge2 <- subset(gene.expr_table_merge2, substr(gene.expr_table_merge2$`#Chr`, 1, 2) != "GL")
#gene.expr_table_merge2$`#Chr` <- substr(gene.expr_table_merge2$`#Chr`, 4, 7)



### For GWAS: select samples that have GWAS data
# modify code for splicing
MSchip_RNAseq_samples_sub <- subset(MSchip_RNAseq_samples, Subset=="CD8")
row.names(MSchip_RNAseq_samples_sub) <- MSchip_RNAseq_samples_sub$Sample_ID2

expr_overlapGWAS <- gene.expr_table_merge2[, colnames(gene.expr_table_merge2) %in% MSchip_RNAseq_samples_sub$Sample_ID2]
expr_overlapGWAS <- cbind(gene.expr_table_merge2[,c(1:4)], expr_overlapGWAS)
MSchip_RNAseq_samples_sub <- MSchip_RNAseq_samples_sub[colnames(expr_overlapGWAS)[5:ncol(expr_overlapGWAS)],]
colnames(expr_overlapGWAS)[5:ncol(expr_overlapGWAS)] <- MSchip_RNAseq_samples_sub$vcfID

nrow(expr_overlapGWAS); length(unique(row.names(expr_overlapGWAS)))
gene_expr_CD8 <- expr_overlapGWAS
#write.table(expr_overlapGWAS, "~/epic.neb/qtl/pheno_expr/gene_expr_CD8.txt", sep="\t", row.names=FALSE, quote = FALSE)
rm(expr_overlapGWAS, MSchip_RNAseq_samples_sub, sTable.sub)


### For RNA-seq variants: select samples for eQTL test     #####
sTable.sub <- subset(sTable.cd8, Last_Known_Treat_Stat=="TreatmentNaive")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
expr <- gene.expr_table_merge2[, colnames(gene.expr_table_merge2) %in% sTable.sub$Sample_ID]
expr <- expr[,vcf_list_sub$sample]
expr <- cbind(gene.expr_table_merge2[,c(1:4)], expr)
nrow(expr); length(unique(row.names(expr)))
expr_CD8_MSuntreat <- expr
rm(expr, sTable.sub, vcf_list_sub)

sTable.sub <- subset(sTable.cd8, Last_Known_Treat_Stat=="Healthy")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
expr <- gene.expr_table_merge2[, colnames(gene.expr_table_merge2) %in% sTable.sub$Sample_ID]
expr <- expr[,vcf_list_sub$sample]
expr <- cbind(gene.expr_table_merge2[,c(1:4)], expr)
nrow(expr); length(unique(row.names(expr)))
expr_CD8_HC <- expr
rm(expr, sTable.sub, vcf_list_sub)

sTable.sub <- subset(sTable.cd8, DiseaseStatus=="MS")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
expr <- gene.expr_table_merge2[, colnames(gene.expr_table_merge2) %in% sTable.sub$Sample_ID]
expr <- expr[,vcf_list_sub$sample]
expr <- cbind(gene.expr_table_merge2[,c(1:4)], expr)
nrow(expr); length(unique(row.names(expr)))
expr_CD8_MSall <- expr
rm(expr, sTable.sub, vcf_list_sub)

rm(gene.expr_table_merge, gene.expr_table_merge2, expr_ncol, chrOrder_gl)




##### CD14: loading selected samples from RSEM ##### 
# subset metadata
sTable.cd14 <- subset(sTable, Subset=="CD14")
#sTable.cd14 <- subset(sTable, Subset=="CD14" & (DiseaseStatus=="MS" | DiseaseStatus=="Healthy") & !is.na(Last_Known_Treat_Stat))  # & !is.na(Last_Known_Treat_Stat)
sTable.cd14$Subset <- droplevels(sTable.cd14$Subset)
sTable.cd14$Last_Known_Treat_Stat <- droplevels(sTable.cd14$Last_Known_Treat_Stat)
sTable.cd14$DiseaseStatus <- droplevels(sTable.cd14$DiseaseStatus)
sTable.cd14[is.na(sTable.cd14$AgeAtExam),]$AgeAtExam <- 36
sTable.cd14$Sex <- as.character(sTable.cd14$Sex)
sTable.cd14[sTable.cd14$Sex == "M",]$Sex <- 0
sTable.cd14[sTable.cd14$Sex == "F",]$Sex <- 1

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

geneList$gene_id2 <- paste0(geneList$gene_id,"_",geneList$gene_name)
geneList <- as.data.frame(geneList)
geneList$gene_id2 <- make.unique(geneList$gene_id2)
row.names(geneList) <- geneList$gene_id2
geneList <- geneList[row.names(rsem.gene.sub$abundance),]

gene.expr_table <- as.data.frame(rsem.gene.sub$abundance)
colnames(gene.expr_table) <- colnames(rsem.gene.sub$abundance)
gene.expr_table$gene_id <- row.names(gene.expr_table)
gene.expr_table_merge <- merge(gene.expr_table, geneList[,c("seqnames","start","end","gene_id2")], by.x="gene_id", by.y="gene_id2")
expr_ncol <- ncol(gene.expr_table_merge)
gene.expr_table_merge <- gene.expr_table_merge[,c((expr_ncol-2):expr_ncol,1,2:(expr_ncol-3))]
chrOrder <- c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
chrOrder_gl <- unique(gene.expr_table_merge$seqnames)
chrOrder_gl <- chrOrder_gl[mixedorder(chrOrder_gl)]
chrOrder <- c(chrOrder,chrOrder_gl[26:length(chrOrder_gl)])
gene.expr_table_merge$seqnames <- factor(gene.expr_table_merge$seqnames, levels=chrOrder)
gene.expr_table_merge <- gene.expr_table_merge[order(gene.expr_table_merge$seqnames, gene.expr_table_merge$start),]
colnames(gene.expr_table_merge) <- c("#Chr","start","end","ID",colnames(gene.expr_table_merge)[5:length(gene.expr_table_merge)])

gene.expr_table_merge2 <- gene.expr_table_merge
#gene.expr_table_merge2 <- subset(gene.expr_table_merge, substr(gene.expr_table_merge$`#Chr`, 1, 2) != "KI")
#gene.expr_table_merge2 <- subset(gene.expr_table_merge2, substr(gene.expr_table_merge2$`#Chr`, 1, 2) != "GL")
#gene.expr_table_merge2$`#Chr` <- substr(gene.expr_table_merge2$`#Chr`, 4, 7)



### For GWAS: select samples that have GWAS data
# modify code for splicing
MSchip_RNAseq_samples_sub <- subset(MSchip_RNAseq_samples, Subset=="CD14")
row.names(MSchip_RNAseq_samples_sub) <- MSchip_RNAseq_samples_sub$Sample_ID2

expr_overlapGWAS <- gene.expr_table_merge2[, colnames(gene.expr_table_merge2) %in% MSchip_RNAseq_samples_sub$Sample_ID2]
expr_overlapGWAS <- cbind(gene.expr_table_merge2[,c(1:4)], expr_overlapGWAS)
MSchip_RNAseq_samples_sub <- MSchip_RNAseq_samples_sub[colnames(expr_overlapGWAS)[5:ncol(expr_overlapGWAS)],]
colnames(expr_overlapGWAS)[5:ncol(expr_overlapGWAS)] <- MSchip_RNAseq_samples_sub$vcfID

nrow(expr_overlapGWAS); length(unique(row.names(expr_overlapGWAS)))
gene_expr_CD14 <- expr_overlapGWAS
#write.table(expr_overlapGWAS, "~/epic.neb/qtl/pheno_expr/gene_expr_CD14.txt", sep="\t", row.names=FALSE, quote = FALSE)
rm(expr_overlapGWAS, MSchip_RNAseq_samples_sub, sTable.sub)


### For RNA-seq variants: CD14 - random select 20 patients
sTable.sub <- subset(sTable.cd14, DiseaseStatus=="MS")
sTable.sub <- sTable.sub[sTable.sub$Sample_ID %in% sample(subset(sTable.cd14, DiseaseStatus=="MS")$Sample_ID, 22),]
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
expr <- gene.expr_table_merge2[, colnames(gene.expr_table_merge2) %in% sTable.sub$Sample_ID]
expr <- expr[,vcf_list_sub$sample]
expr <- cbind(gene.expr_table_merge2[,c(1:4)], expr)
nrow(expr); length(unique(row.names(expr)))
#
expr <- expr_transform_new(expr)
write.table(expr, "~/epic.neb/qtl/test_20samp/CD14_test2_expr.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(vcf_list_sub$path, file="~/epic.neb/qtl/test_20samp/CD14_test2_vcf.list", quote=FALSE, row.names=FALSE, col.names=FALSE)
#
covar <- sTable.sub[,c("Sex","AgeAtExam")]
row.names(covar) <- sTable.sub$Sample_ID
covar$id <- sTable.sub$Sample_ID
covar <- covar[vcf_list_sub$sample,]
covar <- t(covar[,c(3,1:2)])
write.table(covar, "~/epic.neb/qtl/test_20samp/CD14_test2_covar.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)



### For RNA-seq variants: select samples for eQTL test     #####
sTable.sub <- subset(sTable.cd14, Last_Known_Treat_Stat=="TreatmentNaive")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
expr <- gene.expr_table_merge2[, colnames(gene.expr_table_merge2) %in% sTable.sub$Sample_ID]
expr <- expr[,vcf_list_sub$sample]
expr <- cbind(gene.expr_table_merge2[,c(1:4)], expr)
nrow(expr); length(unique(row.names(expr)))
expr_CD14_MSuntreat <- expr
rm(expr, sTable.sub, vcf_list_sub)

sTable.sub <- subset(sTable.cd14, Last_Known_Treat_Stat=="Healthy")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
expr <- gene.expr_table_merge2[, colnames(gene.expr_table_merge2) %in% sTable.sub$Sample_ID]
expr <- expr[,vcf_list_sub$sample]
expr <- cbind(gene.expr_table_merge2[,c(1:4)], expr)
nrow(expr); length(unique(row.names(expr)))
expr_CD14_HC <- expr
rm(expr, sTable.sub, vcf_list_sub)

sTable.sub <- subset(sTable.cd14, DiseaseStatus=="MS")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
expr <- gene.expr_table_merge2[, colnames(gene.expr_table_merge2) %in% sTable.sub$Sample_ID]
expr <- expr[,vcf_list_sub$sample]
expr <- cbind(gene.expr_table_merge2[,c(1:4)], expr)
nrow(expr); length(unique(row.names(expr)))
expr_CD14_MSall <- expr
rm(expr, sTable.sub, vcf_list_sub)

rm(gene.expr_table_merge, gene.expr_table_merge2, expr_ncol, chrOrder_gl)





##### eQTL: filtering and normalization of gene expression from RSEM abundance
# TPM >= 0.05 in 20% of samples (filt) -> quantile normalization (qn) -> inverse-normal transformation (int)
# TPM >= 0.05 in 20% of samples (filt) -> inverse-normal transformation (int) (it's same as 'int' after 'qn')
gene_expr_CD4[,5:ncol(gene_expr_CD4)]
hist(t(gene_expr_CD4[,5:ncol(gene_expr_CD4)]))

### filtering and normalization
### function for filtering
expr_transform <- function(x){
  # remove low expressed genes
  expr_filt <- data.frame()
  for (i in 1:nrow(x)){
    tmp <- x[i,5:ncol(x)]
    if (ncol(tmp[which(tmp>= 0.05)]) >= ((ncol(x)-4)*0.2)){
      expr_filt <- rbind(expr_filt, x[i,1:ncol(x)])
    } 
  }
  # inverse-normal tansformation 
  mat <- expr_filt[,5:ncol(expr_filt)]
  expr_int <- data.frame()
  for (i in 1:nrow(mat)){
    mat[i,] <- qnorm((rank(mat[i,],na.last="keep")-0.5)/sum(!is.na(mat[i,])))  # inverse-normal transformation: y <‐ qnorm((rank(x, na.last="keep") ‐ 0.5) / sum(!is.na(x))
  }
  mat <- cbind(expr_filt[,1:4], mat)
  return(mat)
}

expr_transform_new <- function(x){
  # remove low expressed genes
  colsum = rowSums(x[,5:ncol(x)] >= 0.05)
  expr_filt <- x[colsum >= ((ncol(x)-4)*0.2),]
  # inverse-normal tansformation
  mat = apply(expr_filt[,5:ncol(expr_filt)], 1, function(x){
    y = qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
    y
  })
  mat = cbind(expr_filt[,1:4], t(mat))
  mat
}


# quantile normalization (from MatrixEQTL: http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/faq.html#qnorm)
#mat = expr_filt[,5:ncol(expr_filt)]
#mat = t(apply(mat, 1, rank, ties.method = "average"));
#mat = qnorm(mat / (ncol(mat)+1));
#expr_qn = mat;
#expr_qn_int <- qnorm((rank(expr_qn[2,],na.last="keep")-0.5)/sum(!is.na(expr_qn[2,])))  # inverse-normal transformation: y <‐ qnorm((rank(x, na.last="keep") ‐ 0.5) / sum(!is.na(x))

#hist(t(expr_filt[,5:ncol(expr_filt)]))
#hist(t(expr_int[3,5:ncol(expr_int)]))


#gene_expr_CD14_int <- expr_int; 
#rm(expr, expr_filt, expr_int)

#gene_expr_CD8_int <- expr_int; 
#rm(expr, expr_filt, expr_int)

### For GWAS: filtering genes for GWAS-MSchip
gene_expr_CD4_int <- expr_transform_new(gene_expr_CD4)
head(gene_expr_CD4_int)[1:8]
gene_expr_CD8_int <- expr_transform_new(gene_expr_CD8)
head(gene_expr_CD8_int)[1:8]
gene_expr_CD14_int <- expr_transform_new(gene_expr_CD14)
head(gene_expr_CD14_int)[1:8]

head(gene_expr_CD4_int)
tail(gene_expr_CD4_int)
hist(t(gene_expr_CD4_int[3,5:ncol(gene_expr_CD4_int)]))

write.table(gene_expr_CD4_int, "~/epic.neb/qtl/pheno_expr/with_GWAS_MSchip/gene_expr_filt_CD4.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(gene_expr_CD8_int, "~/epic.neb/qtl/pheno_expr/with_GWAS_MSchip/gene_expr_filt_CD8.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(gene_expr_CD14_int, "~/epic.neb/qtl/pheno_expr/with_GWAS_MSchip/gene_expr_filt_CD14.txt", sep="\t", row.names=FALSE, quote = FALSE)


### For RNA-seq variants: filtering genes for RNA-Seq variants (int = inverse-normal tansformed)
expr_CD4_MSuntreat_int <- expr_transform_new(expr_CD4_MSuntreat); dim(expr_CD4_MSuntreat_int)
expr_CD4_HC_int <- expr_transform_new(expr_CD4_HC); dim(expr_CD4_HC_int)
expr_CD4_MSall_int <- expr_transform_new(expr_CD4_MSall); dim(expr_CD4_MSall_int)

expr_CD8_MSuntreat_int <- expr_transform_new(expr_CD8_MSuntreat); dim(expr_CD8_MSuntreat_int)
expr_CD8_HC_int <- expr_transform_new(expr_CD8_HC); dim(expr_CD8_HC_int)
expr_CD8_MSall_int <- expr_transform_new(expr_CD8_MSall); dim(expr_CD8_MSall_int)

expr_CD14_MSuntreat_int <- expr_transform_new(expr_CD14_MSuntreat); dim(expr_CD14_MSuntreat_int)
expr_CD14_HC_int <- expr_transform_new(expr_CD14_HC); dim(expr_CD14_HC_int)
expr_CD14_MSall_int <- expr_transform_new(expr_CD14_MSall); dim(expr_CD14_MSall_int)


write.table(expr_CD4_MSuntreat_int, "~/epic.neb/qtl/pheno_expr/CD4_expr_MS-untreat.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(expr_CD4_HC_int, "~/epic.neb/qtl/pheno_expr/CD4_expr_HC.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(expr_CD4_MSall_int, "~/epic.neb/qtl/pheno_expr/CD4_expr_MSall.txt", sep="\t", row.names=FALSE, quote = FALSE)

write.table(expr_CD8_MSuntreat_int, "~/epic.neb/qtl/pheno_expr/CD8_expr_MS-untreat.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(expr_CD8_HC_int, "~/epic.neb/qtl/pheno_expr/CD8_expr_HC.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(expr_CD8_MSall_int, "~/epic.neb/qtl/pheno_expr/CD8_expr_MSall.txt", sep="\t", row.names=FALSE, quote = FALSE)

write.table(expr_CD14_MSuntreat_int, "~/epic.neb/qtl/pheno_expr/CD14_expr_MS-untreat.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(expr_CD14_HC_int, "~/epic.neb/qtl/pheno_expr/CD14_expr_HC.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(expr_CD14_MSall_int, "~/epic.neb/qtl/pheno_expr/CD14_expr_MSall.txt", sep="\t", row.names=FALSE, quote = FALSE)


save.image("~/epic.neb/analysis/R_neb/RData/fastqtl_input_tables_indivBAM.RData")


save.image("~/epic.neb/analysis/R_neb/RData/fastqtl_input_tables.RData")    # merged bam

load(file="~/epic.neb/analysis/R_neb/RData/fastqtl_input_tables.RData")





#####  covariates  #####
## CD4
sTable.sub <- subset(sTable.cd4, Last_Known_Treat_Stat=="TreatmentNaive")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
covar <- sTable.sub[,c("Sex","AgeAtExam")]
row.names(covar) <- sTable.sub$Sample_ID
covar$id <- sTable.sub$Sample_ID
covar <- covar[vcf_list_sub$sample,]
covar <- t(covar[,c(3,1:2)])
write.table(covar, "~/epic.neb/qtl/pheno_expr/covar_CD4_MS-untreat.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)
covar_cd4_MSuntreat <- covar; rm(covar, sTable.sub, vcf_list_sub)

sTable.sub <- subset(sTable.cd4, Last_Known_Treat_Stat=="Healthy")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
covar <- sTable.sub[,c("Sex","AgeAtExam")]
row.names(covar) <- sTable.sub$Sample_ID
covar$id <- sTable.sub$Sample_ID
covar <- covar[vcf_list_sub$sample,]
covar <- t(covar[,c(3,1:2)])
write.table(covar, "~/epic.neb/qtl/pheno_expr/covar_CD4_HC.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)
covar_cd4_HC <- covar; rm(covar)

sTable.sub <- subset(sTable.cd4, DiseaseStatus=="MS")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
covar <- sTable.sub[,c("Sex","AgeAtExam")]
row.names(covar) <- sTable.sub$Sample_ID
covar$id <- sTable.sub$Sample_ID
covar <- covar[vcf_list_sub$sample,]
covar <- t(covar[,c(3,1:2)])
write.table(covar, "~/epic.neb/qtl/pheno_expr/covar_CD4_MSall.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)
covar_cd4_MSall <- covar; rm(covar)



## CD8
sTable.sub <- subset(sTable.cd8, Last_Known_Treat_Stat=="TreatmentNaive")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
covar <- sTable.sub[,c("Sex","AgeAtExam")]
row.names(covar) <- sTable.sub$Sample_ID
covar$id <- sTable.sub$Sample_ID
covar <- covar[vcf_list_sub$sample,]
covar <- t(covar[,c(3,1:2)])
write.table(covar, "~/epic.neb/qtl/pheno_expr/covar_CD8_MS-untreat.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)
covar_cd8_MSuntreat <- covar; rm(covar)

sTable.sub <- subset(sTable.cd8, Last_Known_Treat_Stat=="Healthy")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
covar <- sTable.sub[,c("Sex","AgeAtExam")]
row.names(covar) <- sTable.sub$Sample_ID
covar$id <- sTable.sub$Sample_ID
covar <- covar[vcf_list_sub$sample,]
covar <- t(covar[,c(3,1:2)])
write.table(covar, "~/epic.neb/qtl/pheno_expr/covar_CD8_HC.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)
covar_cd8_HC <- covar; rm(covar)

sTable.sub <- subset(sTable.cd8, DiseaseStatus=="MS")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
covar <- sTable.sub[,c("Sex","AgeAtExam")]
row.names(covar) <- sTable.sub$Sample_ID
covar$id <- sTable.sub$Sample_ID
covar <- covar[vcf_list_sub$sample,]
covar <- t(covar[,c(3,1:2)])
write.table(covar, "~/epic.neb/qtl/pheno_expr/covar_CD8_MSall.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)
covar_cd8_MSall <- covar; rm(covar)


## CD14
sTable.sub <- subset(sTable.cd14, Last_Known_Treat_Stat=="TreatmentNaive")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
covar <- sTable.sub[,c("Sex","AgeAtExam")]
row.names(covar) <- sTable.sub$Sample_ID
covar$id <- sTable.sub$Sample_ID
covar <- covar[vcf_list_sub$sample,]
covar <- t(covar[,c(3,1:2)])
write.table(covar, "~/epic.neb/qtl/pheno_expr/covar_CD14_MS-untreat.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)
covar_cd14_MSuntreat <- covar; rm(covar)

sTable.sub <- subset(sTable.cd14, Last_Known_Treat_Stat=="Healthy")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
covar <- sTable.sub[,c("Sex","AgeAtExam")]
row.names(covar) <- sTable.sub$Sample_ID
covar$id <- sTable.sub$Sample_ID
covar <- covar[vcf_list_sub$sample,]
covar <- t(covar[,c(3,1:2)])
write.table(covar, "~/epic.neb/qtl/pheno_expr/covar_CD14_HC.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)
covar_cd14_HC <- covar; rm(covar)

sTable.sub <- subset(sTable.cd14, DiseaseStatus=="MS")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]
covar <- sTable.sub[,c("Sex","AgeAtExam")]
row.names(covar) <- sTable.sub$Sample_ID
covar$id <- sTable.sub$Sample_ID
covar <- covar[vcf_list_sub$sample,]
covar <- t(covar[,c(3,1:2)])
write.table(covar, "~/epic.neb/qtl/pheno_expr/covar_CD14_MSall.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)
covar_cd14_MSall <- covar; rm(covar)




### covariates for GWAS eqtl
# CD8
MSchip_RNAseq_samples_sub <- subset(MSchip_RNAseq_samples, Subset=="CD8")
row.names(MSchip_RNAseq_samples_sub) <- MSchip_RNAseq_samples_sub$Sample_ID2

sTable.cd8.sub <- sTable.cd8[sTable.cd8$Sample_ID %in% MSchip_RNAseq_samples$Sample_ID2,]
row.names(sTable.cd8.sub) <- sTable.cd8.sub$Sample_ID
sTable.cd8.sub <- sTable.cd8.sub[MSchip_RNAseq_samples_sub$Sample_ID2,]
covar_cd8 <- sTable.cd8.sub[,c("Sex","AgeAtExam")]
covar_cd8$id <- MSchip_RNAseq_samples_sub$vcfID
covar_cd8 <- t(covar_cd8[,c(3,1:2)])

write.table(covar_cd8, "~/epic.neb/qtl/pheno_expr/gene_covar_CD8.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)

# CD14
MSchip_RNAseq_samples_sub <- subset(MSchip_RNAseq_samples, Subset=="CD14")
row.names(MSchip_RNAseq_samples_sub) <- MSchip_RNAseq_samples_sub$Sample_ID2

sTable.cd14.sub <- sTable.cd14[sTable.cd14$Sample_ID %in% MSchip_RNAseq_samples$Sample_ID2,]
row.names(sTable.cd14.sub) <- sTable.cd14.sub$Sample_ID
sTable.cd14.sub <- sTable.cd14.sub[MSchip_RNAseq_samples_sub$Sample_ID2,]
covar_cd14 <- sTable.cd14.sub[,c("Sex","AgeAtExam")]
covar_cd14$id <- MSchip_RNAseq_samples_sub$vcfID
covar_cd14 <- t(covar_cd14[,c(3,1:2)])

write.table(covar_cd14, "~/epic.neb/qtl/pheno_expr/gene_covar_CD14.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)



##### 



#####  change sample name (remove cell type in the sample ID)   #####
tmp <- expr_CD4_MSuntreat_int
tmp2 <- changeName(expr_CD4_MSuntreat_int)
rm(tmp, tmp2)

changeName <- function(table){
  colnames(table)[5:ncol(table)] <- as.character(lapply(strsplit(colnames(table)[5:ncol(table)], "-"), "[[",1))
  return(table)
}

# gene expression
expr_CD4_MSuntreat_int2 <- changeName(expr_CD4_MSuntreat_int)
expr_CD4_HC_int2 <- changeName(expr_CD4_HC_int)
expr_CD4_MSall_int2 <- changeName(expr_CD4_MSall_int)
write.table(expr_CD4_MSuntreat_int2, "~/epic.neb/qtl/pheno_expr/CD4_expr_MS-untreat.bed", sep="\t", row.names=FALSE, quote = FALSE)
write.table(expr_CD4_HC_int2, "~/epic.neb/qtl/pheno_expr/CD4_expr_HC.bed", sep="\t", row.names=FALSE, quote = FALSE)
write.table(expr_CD4_MSall_int2, "~/epic.neb/qtl/pheno_expr/CD4_expr_MSall.bed", sep="\t", row.names=FALSE, quote = FALSE)

expr_CD8_MSuntreat_int2 <- changeName(expr_CD8_MSuntreat_int)
expr_CD8_HC_int2 <- changeName(expr_CD8_HC_int)
expr_CD8_MSall_int2 <- changeName(expr_CD8_MSall_int)
write.table(expr_CD8_MSuntreat_int2, "~/epic.neb/qtl/pheno_expr/CD8_expr_MS-untreat.bed", sep="\t", row.names=FALSE, quote = FALSE)
write.table(expr_CD8_HC_int2, "~/epic.neb/qtl/pheno_expr/CD8_expr_HC.bed", sep="\t", row.names=FALSE, quote = FALSE)
write.table(expr_CD8_MSall_int2, "~/epic.neb/qtl/pheno_expr/CD8_expr_MSall.bed", sep="\t", row.names=FALSE, quote = FALSE)

expr_CD14_MSuntreat_int2 <- changeName(expr_CD14_MSuntreat_int)
expr_CD14_HC_int2 <- changeName(expr_CD14_HC_int)
expr_CD14_MSall_int2 <- changeName(expr_CD14_MSall_int)
write.table(expr_CD14_MSuntreat_int2, "~/epic.neb/qtl/pheno_expr/CD14_expr_MS-untreat.bed", sep="\t", row.names=FALSE, quote = FALSE)
write.table(expr_CD14_HC_int2, "~/epic.neb/qtl/pheno_expr/CD14_expr_HC.bed", sep="\t", row.names=FALSE, quote = FALSE)
write.table(expr_CD14_MSall_int2, "~/epic.neb/qtl/pheno_expr/CD14_expr_MSall.bed", sep="\t", row.names=FALSE, quote = FALSE)


###
changeNameCovar <- function(table){
  table <- as.data.frame(t(table))
  table$id <- as.character(table$id)
  table$id <- as.character(lapply(strsplit(table$id, "-"), "[[",1))
  table <- t(table)
  return(table)
}
changeNameCovar()

# covariates
covar_cd4_MSuntreat2 <- changeNameCovar(covar_cd4_MSuntreat)
write.table(covar_cd4_MSuntreat2, "~/epic.neb/qtl/pheno_expr/covar_CD4_MS-untreat.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)

covar_cd4_HC2 <- changeNameCovar(covar_cd4_HC)
write.table(covar_cd4_HC2, "~/epic.neb/qtl/pheno_expr/covar_CD4_HC.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)

covar_cd4_MSall2 <- changeNameCovar(covar_cd4_MSall)
write.table(covar_cd4_MSall2, "~/epic.neb/qtl/pheno_expr/covar_CD4_MSall.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)

covar_cd8_MSuntreat2 <- changeNameCovar(covar_cd8_MSuntreat)
write.table(covar_cd8_MSuntreat2, "~/epic.neb/qtl/pheno_expr/covar_CD8_MS-untreat.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)

covar_cd8_HC2 <- changeNameCovar(covar_cd8_HC)
write.table(covar_cd8_HC2, "~/epic.neb/qtl/pheno_expr/covar_CD8_HC.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)

covar_cd8_MSall2 <- changeNameCovar(covar_cd8_MSall)
write.table(covar_cd8_MSall2, "~/epic.neb/qtl/pheno_expr/covar_CD8_MSall.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)

covar_cd14_MSuntreat2 <- changeNameCovar(covar_cd14_MSuntreat)
write.table(covar_cd14_MSuntreat2, "~/epic.neb/qtl/pheno_expr/covar_CD14_MS-untreat.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)

covar_cd14_HC2 <- changeNameCovar(covar_cd14_HC)
write.table(covar_cd14_HC2, "~/epic.neb/qtl/pheno_expr/covar_CD14_HC.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)

covar_cd14_MSall2 <- changeNameCovar(covar_cd14_MSall)
write.table(covar_cd14_MSall2, "~/epic.neb/qtl/pheno_expr/covar_CD14_MSall.txt", sep="\t", row.names=TRUE, col.names = FALSE, quote = FALSE)







##### VCF file list for merge individual variants #####
###  loading metadata
load(file="~/epic.neb/analysis/R_neb/RData/infoTable.RData")
load("~/epic.neb/analysis/R_neb/RData/DESeq2-gene_level-AllCell_MS-all and HC_withAge.RData")

vcf_list_all <- read_csv("~/epic.neb/analysis/all_vcf_list.csv")


###############################
# select specific subjects from vcf file
setwd("~/epic.neb/variantCall")

sTable.sub <- subset(sTable.cd4, Last_Known_Treat_Stat=="TreatmentNaive")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]; rm(sTable.sub)
vcf_list_sub$vcf.file <- paste0("/pool0/home/kicheol/epic.neb/variantCall/indiv_vcf/from_merged_bam/",vcf_list_sub$subject,".vcf.gz")   # /wynton/scratch/kkim/vcf_filtered/
write.table(vcf_list_sub$vcf.file, file="vcf_CD4_TrtNaive.list", quote=FALSE, row.names=FALSE, col.names=FALSE)
rm(vcf_list_sub)

sTable.sub <- subset(sTable.cd4, Last_Known_Treat_Stat=="Healthy")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]; rm(sTable.sub)
vcf_list_sub$vcf.file <- paste0("/wynton/scratch/kkim/vcf_filtered/",vcf_list_sub$subject,".vcf.gz")
write.table(vcf_list_sub$vcf.file, file="vcf_CD4_HC.list", quote=FALSE, row.names=FALSE, col.names=FALSE)
rm(vcf_list_sub)

sTable.sub <- subset(sTable.cd4, DiseaseStatus=="MS")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]; rm(sTable.sub)
vcf_list_sub$vcf.file <- paste0("/wynton/scratch/kkim/vcf_filtered/",vcf_list_sub$subject,".vcf.gz")
write.table(vcf_list_sub$vcf.file, file="vcf_CD4_all-MS.list", quote=FALSE, row.names=FALSE, col.names=FALSE)
rm(vcf_list_sub)

sTable.sub <- subset(sTable.cd8, Last_Known_Treat_Stat=="TreatmentNaive")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]; rm(sTable.sub)
vcf_list_sub$vcf.file <- paste0("/pool0/home/kicheol/epic.neb/variantCall/indiv_vcf/from_merged_bam/",vcf_list_sub$subject,".vcf.gz")
write.table(vcf_list_sub$vcf.file, file="vcf_CD8_TrtNaive.list", quote=FALSE, row.names=FALSE, col.names=FALSE)
rm(vcf_list_sub)

sTable.sub <- subset(sTable.cd8, Last_Known_Treat_Stat=="Healthy")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]; rm(sTable.sub)
vcf_list_sub$vcf.file <- paste0("/wynton/scratch/kkim/vcf_filtered/",vcf_list_sub$subject,".vcf.gz")
write.table(vcf_list_sub$vcf.file, file="vcf_CD8_HC.list", quote=FALSE, row.names=FALSE, col.names=FALSE)
rm(vcf_list_sub)

sTable.sub <- subset(sTable.cd8, DiseaseStatus=="MS")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]; rm(sTable.sub)
vcf_list_sub$vcf.file <- paste0("/wynton/scratch/kkim/vcf_filtered/",vcf_list_sub$subject,".vcf.gz")
write.table(vcf_list_sub$vcf.file, file="vcf_CD8_all-MS.list", quote=FALSE, row.names=FALSE, col.names=FALSE)
rm(vcf_list_sub)


sTable.sub <- subset(sTable.cd14, Last_Known_Treat_Stat=="TreatmentNaive")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]; rm(sTable.sub)
vcf_list_sub$vcf.file <- paste0("/pool0/home/kicheol/epic.neb/variantCall/indiv_vcf/from_merged_bam/",vcf_list_sub$subject,".vcf.gz")
write.table(vcf_list_sub$vcf.file, file="vcf_CD14_TrtNaive.list", quote=FALSE, row.names=FALSE, col.names=FALSE)
rm(vcf_list_sub)

sTable.sub <- subset(sTable.cd14, Last_Known_Treat_Stat=="Healthy")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]; rm(sTable.sub)
vcf_list_sub$vcf.file <- paste0("/wynton/scratch/kkim/vcf_filtered/",vcf_list_sub$subject,".vcf.gz")
write.table(vcf_list_sub$vcf.file, file="vcf_CD14_HC.list", quote=FALSE, row.names=FALSE, col.names=FALSE)
rm(vcf_list_sub)

sTable.sub <- subset(sTable.cd14, DiseaseStatus=="MS")
vcf_list_sub <- vcf_list_all[vcf_list_all$sample %in% sTable.sub$Sample_ID,]; rm(sTable.sub)
vcf_list_sub$vcf.file <- paste0("/wynton/scratch/kkim/vcf_filtered/",vcf_list_sub$subject,".vcf.gz")
write.table(vcf_list_sub$vcf.file, file="vcf_CD14_all-MS.list", quote=FALSE, row.names=FALSE, col.names=FALSE)
rm(vcf_list_sub)

