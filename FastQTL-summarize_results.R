load("~/epic.neb/qtl/results/tmp_eQTL_dataset.RData")



#####
library(readr)

SNPs_info <- read_csv("~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/SNPs_info_grch38.csv")
mschip_hg19 <- read_csv("~/epic.neb/qtl/results/GWAS_mschip/hg19_gwas/MSCHIP_FullHG19BP.csv", col_types = cols(Chr = col_character()))
#imsgc_snps_list <- read_csv("~/epic.neb/qtl/results/GWAS_mschip/imsgc_snps_list.csv")
geneList <- read_csv("~/epic.neb/analysis/info_table/geneList.csv", col_types = cols(X1 = col_skip()))

compare_snps <- merge(SNPs_info, mschip_hg19, by.x="ID", by.y="Name", all.x=TRUE, sort=FALSE, no.dups=FALSE)
length(intersect(paste(compare_snps$Chr,compare_snps$HG19), paste(imsgc_snps_list$Chr, imsgc_snps_list$pos_hg19)))   # n=126 : # of snps included in QTL analysis
imsgc_snps <- compare_snps[paste(compare_snps$Chr,compare_snps$HG19) %in% intersect(paste(compare_snps$Chr,compare_snps$HG19), paste(imsgc_snps_list$Chr, imsgc_snps_list$pos_hg19)),]
imsgc_snps <- unique(paste(imsgc_snps$Chr,imsgc_snps$HG19))

library(biomaRt)
snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
snp_ids = c("rs16828074", "rs17232800")
snp_attributes = c("refsnp_id", "chr_name", "chrom_start","chrom_end","allele","ensembl_type","associated_gene","clinical_significance","phenotype_name")
snp_locations = getBM(attributes=snp_attributes, filters="snp_filter", values=snp_ids, mart=snp_mart)
snp_locations

# variants list
variant_input = "~/epic.neb/qtl/variant_rnaseq/CD4_HC_filt_addID_split.recode.vcf.gz"
vcf <- read.vcfR(variant_input)
variant_table <- data.frame(vcf@fix[,c(1:3)])
colnames(variant_table) <- c("varCHROM","varPOS","varID")





##### Checking that the experiment went well (http://fastqtl.sourceforge.net/pages/cis_permutation.html)
#setwd("~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam")
setwd("~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated")


d = read.table("CD14_expr_HC_permute.results", hea=F, stringsAsFactors=F)
colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
plot(d$ppval, d$bpval, xlab="Direct method", ylab="Beta approximation", main="Check plot")
abline(0, 1, col="red")

d = read.table("CD14_expr_HC_permute.results", hea=F, stringsAsFactors=F)
colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
plot(d$ppval, d$bpval, xlab="Direct method", ylab="Beta approximation", main="Check plot")
abline(0, 1, col="red")

d = read.table("CD4_expr_MSall_permute.results", hea=F, stringsAsFactors=F)
colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
plot(d$ppval, d$bpval, xlab="Direct method", ylab="Beta approximation", main="Check plot")
abline(0, 1, col="red")

d = read.table("CD4_expr_MS-untreat_permute.results", hea=F, stringsAsFactors=F)
colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
plot(d$ppval, d$bpval, xlab="Direct method", ylab="Beta approximation", main="Check plot")
abline(0, 1, col="red")

d = read.table("CD8_expr_MS-untreat_permute.results", hea=F, stringsAsFactors=F)
colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
plot(d$ppval, d$bpval, xlab="Direct method", ylab="Beta approximation", main="Check plot")
abline(0, 1, col="red")
rm(d)


#########################################################################
##########  eQTL results - results with variants from RNA-Seq  ##########
#########################################################################
#load(file="~/epic.neb/analysis/R_neb/RData/fastqtl_results.RData")
load(file="~/epic.neb/analysis/R_neb/RData/fastqtl_results_variant-rnaseq_filtRNAedit-updated.RData")


library(readr)
setwd("~/epic.neb/qtl/results/eQTL_RNAseqVariant")
setwd("~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam")
setwd("~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated")
#setwd("~/epic.neb/qtl/results/eQTL_RNAseqVariant/RNAvar_eqtl_HC")


res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/CD8_expr_MS-untreat_permute.results"
out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD8_MSuntreat_permute_corrected.csv"
variant_input = "~/epic.neb/variantCall/merged_vcf-updated/bcftools-merge/CD8_TrtNaive_filtered_final.recode.vcf.gz"


### function for permutation results
require(vcfR)
fastqtl_summary_permute <- function(res_file, out_file, variant_input){
  require(vcfR)
  vcf <- read.vcfR(variant_input)
  variant_table <- data.frame(vcf@fix[,c(1:3)])
  colnames(variant_table) <- c("varCHROM","varPOS","varID")
  
  eqtl_res_permute <- read_table2(res_file, col_names = FALSE)
  colnames(eqtl_res_permute) <- c("gene","nvar","shape1","shape2","dummy","sid","dist","npval","slope","ppval","bpval")
  eqtl_res_permute$bonferroni = p.adjust(eqtl_res_permute$bpval, method="bonferroni")
  eqtl_res_permute$bh = p.adjust(eqtl_res_permute$bpval, method="fdr")
  eqtl_res_permute$geneID <- unlist(lapply(strsplit(eqtl_res_permute$gene, "_"), "[[", 1))
  eqtl_res_permute$geneName <- unlist(lapply(strsplit(eqtl_res_permute$gene, "_"), "[[", 2))
  
  tmp_geneList <- unique(geneList[,c("gene_id","seqnames","start","end","gene_type","gene_name")])
  tmp_geneList <- tmp_geneList[tmp_geneList$gene_id %in% eqtl_res_permute$geneID,]
  eqtl_res_permute_add <- merge(eqtl_res_permute, tmp_geneList, by.x="geneID", by.y="gene_id", sort=FALSE)
  
  tmp_variants <- variant_table[variant_table$varID %in% eqtl_res_permute_add$sid,]
  tmp_variants <- tmp_variants[!duplicated(paste0(tmp_variants$varCHROM, tmp_variants$varPOS)),]
  eqtl_res_permute_add <- merge(eqtl_res_permute_add, tmp_variants, by.x="sid", by.y="varID", sort=FALSE)
  
  eqtl_res_permute_add <- eqtl_res_permute_add[,c("gene","nvar","shape1","shape2","dummy","sid","dist","npval",
                                                  "slope","ppval","bpval","bonferroni","bh",
                                                  "geneName","seqnames","start","end","varCHROM","varPOS","gene_type")]
  eqtl_res_permute_add <- eqtl_res_permute_add[order(eqtl_res_permute_add$bpval),]
  write.csv(eqtl_res_permute_add, file=out_file)
  eqtl_res_permute_add
}  


### function for nominal results
fastqtl_summary_nominal <- function(res_file, out_file, variant_input){
  require(vcfR)
  vcf <- read.vcfR(variant_input)
  variant_table <- data.frame(vcf@fix[,c(1:3)])
  colnames(variant_table) <- c("varCHROM","varPOS","varID")
  
  eqtl_res_nominal <- read_table2(res_file, col_names = FALSE)
  colnames(eqtl_res_nominal) <- c("gene","sid","dist","npval","slope")
  eqtl_res_nominal$bonferroni = p.adjust(eqtl_res_nominal$npval, method="bonferroni")
  eqtl_res_nominal$bh = p.adjust(eqtl_res_nominal$npval, method="fdr")
  eqtl_res_nominal$geneID <- unlist(lapply(strsplit(eqtl_res_nominal$gene, "_"), "[[", 1))

  tmp_geneList <- unique(geneList[,c("gene_id","seqnames","start","end","gene_type","gene_name")])
  tmp_geneList <- tmp_geneList[tmp_geneList$gene_id %in% eqtl_res_nominal$geneID,]
  eqtl_res_nominal_add <- merge(eqtl_res_nominal, tmp_geneList, by.x="geneID", by.y="gene_id", sort=FALSE)
  
  tmp_variants <- variant_table[variant_table$varID %in% eqtl_res_nominal_add$sid,]
  tmp_variants <- tmp_variants[!duplicated(paste0(tmp_variants$varCHROM,tmp_variants$varPOS)),]
  eqtl_res_nominal_add <- merge(eqtl_res_nominal_add, tmp_variants, by.x="sid", by.y="varID", sort=FALSE, all.x=TRUE)
  
  eqtl_res_nominal_add <- eqtl_res_nominal_add[,c("gene","sid","dist","npval","bonferroni","bh","slope",
                                                  "gene_name","seqnames","start","end","varCHROM","varPOS","gene_type")]
  eqtl_res_nominal_add <- eqtl_res_nominal_add[order(eqtl_res_nominal_add$npval),]
  write.csv(eqtl_res_nominal_add, file=out_file)
  eqtl_res_nominal_add
}  



# gene list
geneList <- read_csv("~/epic.neb/analysis/info_table/geneList.csv", col_types = cols(X1 = col_skip()))

# run function - permutation results
setwd("~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/")
#
res_CD4_MSuntreat_permute <- fastqtl_summary_permute(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/CD4_expr_MS-untreat_permute.results", 
                                                     out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD4_MS-untreat_permute_corrected.csv",
                                                     variant_input = "~/epic.neb/variantCall/merged_vcf-updated/bcftools-merge/CD4_TrtNaive_filtered_final.recode.vcf.gz")
res_CD8_MSuntreat_permute <- fastqtl_summary_permute(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/CD8_expr_MS-untreat_permute.results", 
                                                     out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD8_MS-untreat_permute_corrected.csv",
                                                     variant_input = "~/epic.neb/variantCall/merged_vcf-updated/bcftools-merge/CD8_TrtNaive_filtered_final.recode.vcf.gz")
res_CD14_MSuntreat_permute <- fastqtl_summary_permute(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/CD14_expr_MS-untreat_permute.results", 
                                                      out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD14_MS-untreat_permute_corrected.csv",
                                                      variant_input = "~/epic.neb/variantCall/merged_vcf-updated/bcftools-merge/CD14_TrtNaive_filtered_final.recode.vcf.gz")
#
res_CD4_HC_permute <- fastqtl_summary_permute(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit/CD4_expr_HC_permute.results", 
                                              out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD4_HC_permute_corrected.csv",
                                              variant_input = "~/epic.neb/qtl/variant_rnaseq/merged_bam/filter_RNAedit/_old_files/CD4_HC_filtered_final.recode.vcf.gz")                                              
res_CD8_HC_permute <- fastqtl_summary_permute(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit/CD8_expr_HC_permute.results", 
                                              out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD8_HC_permute_corrected.csv",
                                              variant_input = "~/epic.neb/qtl/variant_rnaseq/merged_bam/filter_RNAedit/_old_files/CD8_HC_filtered_final.recode.vcf.gz")
res_CD14_HC_permute <- fastqtl_summary_permute(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit/CD14_expr_HC_permute.results", 
                                              out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD14_HC_permute_corrected.csv",
                                              variant_input = "~/epic.neb/qtl/variant_rnaseq/merged_bam/filter_RNAedit/_old_files/CD14_HC_filtered_final.recode.vcf.gz")

#
#res_CD4_MSall_permute <- fastqtl_summary_permute(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/CD4_expr_MSall_permute.results", 
#                                                 out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD4_MSall_permute_corrected.csv",
#                                                 variant_input = "~/epic.neb/qtl/variant_rnaseq/merged_bam/filter_RNAedit/CD4_MSall_filtered_final.recode.vcf.gz")
#res_CD8_MSall_permute <- fastqtl_summary_permute(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/CD8_expr_MSall_permute.results", 
#                                                 out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD8_MSall_permute_corrected.csv",
#                                                 variant_input = "~/epic.neb/qtl/variant_rnaseq/merged_bam/filter_RNAedit/CD8_MSall_filtered_final.recode.vcf.gz")
#res_CD14_MSall_permute <- fastqtl_summary_permute(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/CD14_expr_MSall_permute.results", 
#                                                 out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD14_MSall_permute_corrected.csv",
#                                                 variant_input = "~/epic.neb/qtl/variant_rnaseq/merged_bam/filter_RNAedit/CD14_MSall_filtered_final.recode.vcf.gz")
#


# run function - nominal results
res_CD4_MSuntreat_nominal <- fastqtl_summary_nominal(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/CD4_expr_MS-untreat_nominal.results", 
                                                     out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD4_MS-untreat_nominal_corrected.csv",
                                                     variant_input = "~/epic.neb/variantCall/merged_vcf-updated/bcftools-merge/CD4_TrtNaive_filtered_final.recode.vcf.gz")
res_CD8_MSuntreat_nominal <- fastqtl_summary_nominal(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/CD8_expr_MS-untreat_nominal.results", 
                                                     out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD8_MS-untreat_nominal_corrected.csv",
                                                     variant_input = "~/epic.neb/variantCall/merged_vcf-updated/bcftools-merge/CD8_TrtNaive_filtered_final.recode.vcf.gz")
res_CD14_MSuntreat_nominal <- fastqtl_summary_nominal(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/CD14_expr_MS-untreat_nominal.results", 
                                                      out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD14_MS-untreat_nominal_corrected.csv",
                                                      variant_input = "~/epic.neb/variantCall/merged_vcf-updated/bcftools-merge/CD14_TrtNaive_filtered_final.recode.vcf.gz")

res_CD4_HC_nominal <- fastqtl_summary_nominal(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit/CD4_expr_HC_nominal.results", 
                                              out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD4_HC_nominal_corrected.csv",
                                              variant_input = "~/epic.neb/qtl/variant_rnaseq/merged_bam/filter_RNAedit/_old_files/CD4_HC_filtered_final.recode.vcf.gz")
res_CD8_HC_nominal <- fastqtl_summary_nominal(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit/CD8_expr_HC_nominal.results", 
                                              out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD8_HC_nominal_corrected.csv",
                                              variant_input = "~/epic.neb/qtl/variant_rnaseq/merged_bam/filter_RNAedit/_old_files/CD8_HC_filtered_final.recode.vcf.gz")
res_CD14_HC_nominal <- fastqtl_summary_nominal(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit/CD14_expr_HC_nominal.results", 
                                              out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD14_HC_nominal_corrected.csv",
                                              variant_input = "~/epic.neb/qtl/variant_rnaseq/merged_bam/filter_RNAedit/_old_files/CD14_HC_filtered_final.recode.vcf.gz")


#res_CD4_MSall_nominal <- fastqtl_summary_nominal(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/CD4_expr_MSall_nominal.results", 
#                                                 out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD4_MSall_nominal_corrected.csv",
#                                                 variant_input = "~/epic.neb/qtl/variant_rnaseq/merged_bam/filter_RNAedit/CD4_MSall_filtered_final.recode.vcf.gz")
#res_CD8_MSall_nominal <- fastqtl_summary_nominal(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/CD8_expr_MSall_nominal.results", 
#                                                 out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD8_MSall_nominal_corrected.csv",
#                                                 variant_input = "~/epic.neb/qtl/variant_rnaseq/merged_bam/filter_RNAedit/CD8_MSall_filtered_final.recode.vcf.gz")
#res_CD14_MSall_nominal <- fastqtl_summary_nominal(res_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/CD14_expr_MSall_nominal.results", 
#                                                 out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD14_MSall_nominal_corrected.csv",
#                                                 variant_input = "~/epic.neb/qtl/variant_rnaseq/merged_bam/filter_RNAedit/CD14_MSall_filtered_final.recode.vcf.gz")



#####  results with GWAS (new)
#res_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD4_expr_permute.results"
#out_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD4_expr_corrected.csv"
#variant_input = "~/epic.neb/qtl/GWAS_mschip/UCSF-MSchip_29overlap_remapGRCh38_rehead_sorted_filt.recode.vcf.gz"


res_CD4_gwas_permute <- fastqtl_summary_permute(res_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD4_expr_permute.results", 
                                                out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD4-GWAS_expr_permute_corrected.csv",
                                                variant_input = "~/epic.neb/qtl/GWAS_mschip/UCSF-MSchip_29overlap_remapGRCh38_rehead_sorted_filt.recode.vcf.gz")
res_CD8_gwas_permute <- fastqtl_summary_permute(res_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD8_expr_permute.results", 
                                                out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD8-GWAS_expr_permute_corrected.csv",
                                                variant_input = "~/epic.neb/qtl/GWAS_mschip/UCSF-MSchip_29overlap_remapGRCh38_rehead_sorted_filt.recode.vcf.gz")
res_CD14_gwas_permute <- fastqtl_summary_permute(res_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD14_expr_permute.results", 
                                                 out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD14-GWAS_expr_permute_corrected.csv",
                                                 variant_input = "~/epic.neb/qtl/GWAS_mschip/UCSF-MSchip_29overlap_remapGRCh38_rehead_sorted_filt.recode.vcf.gz")

res_CD4_gwas_nominal <- fastqtl_summary_nominal(res_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD4_expr_nominal.results", 
                                                out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD4-GWAS_expr_nominal_corrected.csv",
                                                variant_input = "~/epic.neb/qtl/GWAS_mschip/UCSF-MSchip_29overlap_remapGRCh38_rehead_sorted_filt.recode.vcf.gz")
res_CD8_gwas_nominal <- fastqtl_summary_nominal(res_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD8_expr_nominal.results", 
                                                out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD8-GWAS_expr_nominal_corrected.csv",
                                                variant_input = "~/epic.neb/qtl/GWAS_mschip/UCSF-MSchip_29overlap_remapGRCh38_rehead_sorted_filt.recode.vcf.gz")
res_CD14_gwas_nominal <- fastqtl_summary_nominal(res_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD14_expr_nominal.results", 
                                                 out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD14-GWAS_expr_nominal_corrected.csv",
                                                 variant_input = "~/epic.neb/qtl/GWAS_mschip/UCSF-MSchip_29overlap_remapGRCh38_rehead_sorted_filt.recode.vcf.gz")


##
#res_CD4_gwas_permute <- fastqtl_summary_permute(res_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD4_expr_permute.results", 
#                        out_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD4_expr_permute_corrected.csv",
#                        variant_input = "~/epic.neb/qtl/GWAS_mschip/UCSF-MSchip_29overlap_remapGRCh38_rehead_sorted_filt.recode.vcf.gz")
#res_CD8_gwas_permute <- fastqtl_summary_permute(res_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD8_expr_permute.results", 
#                        out_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD8_expr_permute_corrected.csv",
#                        variant_input = "~/epic.neb/qtl/GWAS_mschip/UCSF-MSchip_29overlap_remapGRCh38_rehead_sorted_filt.recode.vcf.gz")
#res_CD14_gwas_permute <- fastqtl_summary_permute(res_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD14_expr_permute.results", 
#                        out_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD14_expr_permute_corrected.csv",
#                        variant_input = "~/epic.neb/qtl/GWAS_mschip/UCSF-MSchip_29overlap_remapGRCh38_rehead_sorted_filt.recode.vcf.gz")

#res_CD4_gwas_nominal <- fastqtl_summary_nominal(res_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD4_expr_nominal.results", 
#                        out_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD4_expr_nominal_corrected.csv",
#                        variant_input = "~/epic.neb/qtl/GWAS_mschip/UCSF-MSchip_29overlap_remapGRCh38_rehead_sorted_filt.recode.vcf.gz")
#res_CD8_gwas_nominal <- fastqtl_summary_nominal(res_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD8_expr_nominal.results", 
#                        out_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD8_expr_nominal_corrected.csv",
#                        variant_input = "~/epic.neb/qtl/GWAS_mschip/UCSF-MSchip_29overlap_remapGRCh38_rehead_sorted_filt.recode.vcf.gz")
#res_CD14_gwas_nominal <- fastqtl_summary_nominal(res_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD14_expr_nominal.results", 
#                        out_file="~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr3/CD14_expr_nominal_corrected.csv",
#                        variant_input = "~/epic.neb/qtl/GWAS_mschip/UCSF-MSchip_29overlap_remapGRCh38_rehead_sorted_filt.recode.vcf.gz")

###
save.image("~/epic.neb/analysis/R_neb/RData/fastqtl_results_variant-rnaseq_filtRNAedit-updated.RData")

#save.image("~/epic.neb/analysis/R_neb/RData/fastqtl_results_variant-rnaseq.RData")


###
#load(file="~/epic.neb/analysis/R_neb/RData/fastqtl_results_variant-rnaseq.RData")
#load(file="~/epic.neb/analysis/R_neb/RData/fastqtl_results_variant-rnaseq_filtRNAedit.RData")
load(file="~/epic.neb/analysis/R_neb/RData/fastqtl_results_variant-rnaseq_filtRNAedit-updated.RData")



#####  test run - random sampling 
setwd("~/epic.neb/qtl/test_20samp_random")
res_CD4_20random <- test_fastqtl_summary(res_file="results/CD4_test1_permute.results",
                        out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD4_test1_permute_corrected.csv",
                        variant_input = "CD4_test1_vcf.list.vcf.gz_filt.recode.vcf.gz")
res_CD14_20random <- test_fastqtl_summary(res_file="results/CD14_test1_permute.results",
                                          out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD14_test1_permute_corrected.csv",
                                          variant_input = "CD14_test1_vcf.list.vcf.gz_filt.recode.vcf.gz")

test_fastqtl_summary(res_file="results/CD4_test2_permute.results",
                        out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD4_test2_permute_corrected.csv",
                        variant_input = "CD4_test2_vcf.list.vcf.gz_filt.recode.vcf.gz")

test_fastqtl_summary(res_file="results/CD14_test2_permute.results",
                        out_file="~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD14_test2_permute_corrected.csv",
                        variant_input = "CD14_test2_vcf.list.vcf.gz_filt.recode.vcf.gz")


fastqtl_summary_nominal(res_file="results/CD4_test1_nominal.results", 
                        out_file="results/CD4_test1_nominal_corrected.csv",
                        variant_input = "CD4_test1_vcf.list.vcf.gz_filt.recode.vcf.gz")
fastqtl_summary_nominal(res_file="results/CD14_test1_nominal.results", 
                        out_file="results/CD14_test1_nominal_corrected.csv",
                        variant_input = "CD14_test1_vcf.list.vcf.gz_filt.recode.vcf.gz")


res_file="results/CD4_test1_permute.results"
out_file="results/CD4_test1_permute_corrected.csv"
variant_input = "CD4_test1_vcf.list.vcf.gz_filt.recode.vcf.gz"


##
test_fastqtl_summary <- function(res_file, out_file, variant_input){
vcf <- read.vcfR(variant_input)
variant_table <- data.frame(vcf@fix[,c(1:3)])
colnames(variant_table) <- c("varCHROM","varPOS","varID")

eqtl_res_permute <- read_table2(res_file, col_names = FALSE)
colnames(eqtl_res_permute) <- c("gene","nvar","shape1","shape2","dummy","sid","dist","npval","slope","ppval","bpval")
eqtl_res_permute$bonferroni = p.adjust(eqtl_res_permute$bpval, method="bonferroni")
eqtl_res_permute$bh = p.adjust(eqtl_res_permute$bpval, method="fdr")
eqtl_res_permute$geneID <- unlist(lapply(strsplit(eqtl_res_permute$gene, "_"), "[[", 1))
eqtl_res_permute$geneName <- unlist(lapply(strsplit(eqtl_res_permute$gene, "_"), "[[", 2))

tmp_geneList <- unique(geneList[,c("gene_id","seqnames","start","end","gene_type","gene_name")])
tmp_geneList <- tmp_geneList[tmp_geneList$gene_id %in% eqtl_res_permute$geneID,]
eqtl_res_permute_add <- merge(eqtl_res_permute, tmp_geneList, by.x="geneID", by.y="gene_id", sort=FALSE)

tmp_variants <- variant_table[variant_table$varID %in% eqtl_res_permute_add$sid,]
tmp_variants <- unique(tmp_variants)
#eqtl_res_permute_add <- merge(eqtl_res_permute_add, tmp_variants, by.x="sid", by.y="varID", sort=FALSE)

eqtl_res_permute_add <- eqtl_res_permute_add[,c("gene","nvar","shape1","shape2","dummy","sid","dist","npval",
                                                "slope","ppval","bpval","bonferroni","bh",
                                                "geneName","seqnames","start","end","gene_type")]
eqtl_res_permute_add <- eqtl_res_permute_add[order(eqtl_res_permute_add$bpval),]
write.csv(eqtl_res_permute_add, file=out_file)
eqtl_res_permute_add
}




##########################################################################################
###
#load(file="~/epic.neb/analysis/R_neb/RData/fastqtl_results_variant-rnaseq.RData")
#load(file="~/epic.neb/analysis/R_neb/RData/fastqtl_results_variant-rnaseq_filtRNAedit.RData")
load(file="~/epic.neb/analysis/R_neb/RData/fastqtl_results_variant-rnaseq_filtRNAedit-updated.RData")



##### QQ plot
library(RColorBrewer)

head(res_CD4_HC_permute)
head(res_CD4_MSuntreat_permute)
#res_CD4_MSall_permute


### QQplot - final
my.pvalue.list <- list("CD4_MS-naive"=res_CD4_MSuntreat_permute$bpval[!is.na(res_CD4_MSuntreat_permute$bpval)], 
                       "CD8_MS-naive"=res_CD8_MSuntreat_permute$bpval[!is.na(res_CD8_MSuntreat_permute$bpval)], 
                       "CD14_MS-naive"=res_CD14_MSuntreat_permute$bpval[!is.na(res_CD14_MSuntreat_permute$bpval)],
                       "CD4_genotyped"=res_CD4_gwas_permute$bpval[!is.na(res_CD4_gwas_permute$bpval)], 
                       "CD8_genotyped"=res_CD8_gwas_permute$bpval[!is.na(res_CD8_gwas_permute$bpval)], 
                       "CD14_genotyped"=res_CD14_gwas_permute$bpval[!is.na(res_CD14_gwas_permute$bpval)],
                       "CD4_HC"=res_CD4_HC_permute$bpval[!is.na(res_CD4_HC_permute$bpval)], 
                       "CD8_HC"=res_CD8_HC_permute$bpval[!is.na(res_CD8_HC_permute$bpval)], 
                       "CD14_HC"=res_CD14_HC_permute$bpval[!is.na(res_CD14_HC_permute$bpval)],
                       "CD4_20random"=res_CD4_20random$bpval[!is.na(res_CD4_20random$bpval)],
                       "CD14_20random"=res_CD14_20random$bpval[!is.na(res_CD14_20random$bpval)])

qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.02)), xlim=c(0,5.5), draw.conf=TRUE,
            par.settings=list(superpose.symbol=list(pch=20, 
                                                    col = c(brewer.pal(9, "YlOrRd")[5], brewer.pal(9, "YlOrRd")[6], brewer.pal(9, "YlOrRd")[7], 
                                                            brewer.pal(9, "Blues")[5], brewer.pal(9, "Blues")[6], brewer.pal(9, "Blues")[7], 
                                                            brewer.pal(9, "Greens")[5], brewer.pal(9, "Greens")[6], brewer.pal(9, "Greens")[7], 
                                                            brewer.pal(9, "Greys")[7], brewer.pal(9, "Greys")[8]))))



###
table_plot <- res_CD4_MSuntreat_permute
qqunif.plot(table_plot$bh[!is.na(table_plot$bh)])

my.pvalue.list<-list("CD4_HC"=res_CD4_HC_permute$bpval[!is.na(res_CD4_HC_permute$bpval)], 
                     "CD8_HC"=res_CD8_HC_permute$bpval[!is.na(res_CD8_HC_permute$bpval)], 
                     "CD14_HC"=res_CD14_HC_permute$bpval[!is.na(res_CD14_HC_permute$bpval)])
qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))

my.pvalue.list<-list("CD4_MSuntrt"=res_CD4_MSuntreat_permute$bpval[!is.na(res_CD4_MSuntreat_permute$bpval)], 
                     "CD8_MSuntrt"=res_CD8_MSuntreat_permute$bpval[!is.na(res_CD8_MSuntreat_permute$bpval)], 
                     "CD14_MSuntrt"=res_CD14_MSuntreat_permute$bpval[!is.na(res_CD14_MSuntreat_permute$bpval)])
qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))

my.pvalue.list<-list("CD4_MSall"=res_CD4_MSall_permute$bpval[!is.na(res_CD4_MSall_permute$bpval)], 
                     "CD8_MSall"=res_CD8_MSall_permute$bpval[!is.na(res_CD8_MSall_permute$bpval)], 
                     "CD14_MSall"=res_CD14_MSall_permute$bpval[!is.na(res_CD14_MSall_permute$bpval)])
qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))


res_CD4_MSall_permute <- read_csv("~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD4_MSall_permute_corrected.csv")
res_CD8_MSall_permute <- read_csv("~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD8_MSall_permute_corrected.csv")
res_CD14_MSall_permute <- read_csv("~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD14_MSall_permute_corrected.csv")

res_CD4_HC_permute <- read_csv("~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD4_HC_permute_corrected.csv")
res_CD8_HC_permute <- read_csv("~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD8_HC_permute_corrected.csv")
res_CD14_HC_permute <- read_csv("~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD14_HC_permute_corrected.csv")

res_CD4_MSuntreat_permute <- read_csv("~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD4_MSuntreat_permute_corrected.csv")
res_CD8_MSuntreat_permute <- read_csv("~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD8_MSuntreat_permute_corrected.csv")
res_CD14_MSuntreat_permute <- read_csv("~/epic.neb/qtl/results/eQTL_RNAseqVariant/from_merged_bam/filt0.6-filtRNAedit-updated/corrected/CD14_MSuntreat_permute_corrected.csv")



my.pvalue.list<-list("CD4_HC"=res_CD4_HC_permute$bpval[!is.na(res_CD4_HC_permute$bpval)], 
                     "CD8_HC"=res_CD8_HC_permute$bpval[!is.na(res_CD8_HC_permute$bpval)], 
                     "CD14_HC"=res_CD14_HC_permute$bpval[!is.na(res_CD14_HC_permute$bpval)],
                     "CD4_MSuntrt"=res_CD4_MSuntreat_permute$bpval[!is.na(res_CD4_MSuntreat_permute$bpval)], 
                     "CD8_MSuntrt"=res_CD8_MSuntreat_permute$bpval[!is.na(res_CD8_MSuntreat_permute$bpval)], 
                     "CD14_MSuntrt"=res_CD14_MSuntreat_permute$bpval[!is.na(res_CD14_MSuntreat_permute$bpval)],
                     "CD4_MSall"=res_CD4_MSall_permute$bpval[!is.na(res_CD4_MSall_permute$bpval)], 
                     "CD8_MSall"=res_CD8_MSall_permute$bpval[!is.na(res_CD8_MSall_permute$bpval)], 
                     "CD14_MSall"=res_CD14_MSall_permute$bpval[!is.na(res_CD14_MSall_permute$bpval)],
                     "CD4_gwas"=res_CD4_gwas_permute$bpval[!is.na(res_CD4_gwas_permute$bpval)], 
                     "CD8_gwas"=res_CD8_gwas_permute$bpval[!is.na(res_CD8_gwas_permute$bpval)], 
                     "CD14_gwas"=res_CD14_gwas_permute$bpval[!is.na(res_CD14_gwas_permute$bpval)])
qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))




my.pvalue.list<-list("CD4_HC"=res_CD4_HC_permute$bpval[!is.na(res_CD4_HC_permute$bpval)], 
                     "CD4_MSuntrt"=res_CD4_MSuntreat_permute$bpval[!is.na(res_CD4_MSuntreat_permute$bpval)], 
                     "CD4_MSall"=res_CD4_MSall_permute$bpval[!is.na(res_CD4_MSall_permute$bpval)], 
                     "CD4_gwas"=CD4_expr$bpval[!is.na(CD4_expr$bpval)])
qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))

my.pvalue.list<-list("CD8_HC"=res_CD8_HC_permute$bpval[!is.na(res_CD8_HC_permute$bpval)], 
                     "CD8_MSuntrt"=res_CD8_MSuntreat_permute$bpval[!is.na(res_CD8_MSuntreat_permute$bpval)], 
                     "CD8_MSall"=res_CD8_MSall_permute$bpval[!is.na(res_CD8_MSall_permute$bpval)], 
                     "CD8_gwas"=CD8_expr$bpval[!is.na(CD8_expr$bpval)])
qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))

my.pvalue.list<-list("CD14_HC"=res_CD14_HC_permute$bpval[!is.na(res_CD14_HC_permute$bpval)],
                     "CD14_MSuntrt"=res_CD14_MSuntreat_permute$bpval[!is.na(res_CD14_MSuntreat_permute$bpval)],
                     "CD14_MSall"=res_CD14_MSall_permute$bpval[!is.na(res_CD14_MSall_permute$bpval)],
                     "CD14_gwas"=CD14_expr$bpval[!is.na(CD14_expr$bpval)])
qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))


#






##### eQTL results - results with GWAS
library(readr)
setwd("~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted")

CD4_expr <- read_table2("~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr2/CD4_expr.results", col_names = FALSE)
colnames(CD4_expr) <- c("gene","nvar","shape1","shape2","dummy","sid","dist","npval","slope","ppval","bpval")
CD4_expr$bonferroni = p.adjust(CD4_expr$bpval, method="bonferroni")
CD4_expr$bh = p.adjust(CD4_expr$bpval, method="fdr")

CD8_expr <- read_table2("~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr2/CD8_expr.results", col_names = FALSE)
colnames(CD8_expr) <- c("gene","nvar","shape1","shape2","dummy","sid","dist","npval","slope","ppval","bpval")
CD8_expr$bonferroni = p.adjust(CD8_expr$bpval, method="bonferroni")
CD8_expr$bh = p.adjust(CD8_expr$bpval, method="fdr")

CD14_expr <- read_table2("~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/results_expr2/CD14_expr.results", col_names = FALSE)
colnames(CD14_expr) <- c("gene","nvar","shape1","shape2","dummy","sid","dist","npval","slope","ppval","bpval")
CD14_expr$bonferroni = p.adjust(CD14_expr$bpval, method="bonferroni")
CD14_expr$bh = p.adjust(CD14_expr$bpval, method="fdr")


SNPs_info <- read_csv("~/epic.neb/qtl/results/GWAS_mschip/grch38_gwas_lifted/SNPs_info_grch38.csv")
CD4_expr <- merge(CD4_expr, SNPs_info, by.x="sid", by.y="ID", all.x=TRUE, sort = FALSE, no.dups = FALSE)
CD8_expr <- merge(CD8_expr, SNPs_info, by.x="sid", by.y="ID", all.x=TRUE, sort = FALSE, no.dups = FALSE)
CD14_expr <- merge(CD14_expr, SNPs_info, by.x="sid", by.y="ID", all.x=TRUE, sort = FALSE, no.dups = FALSE)

mschip_hg19 <- read_csv("../hg19_gwas/MSCHIP_FullHG19BP.csv", col_types = cols(Chr = col_character()))
CD4_expr <- merge(CD4_expr, mschip_hg19[,c(1,3)], by.x="sid", by.y="Name", all.x=TRUE, sort = FALSE, no.dups = FALSE)
CD8_expr <- merge(CD8_expr, mschip_hg19[,c(1,3)], by.x="sid", by.y="Name", all.x=TRUE, sort = FALSE, no.dups = FALSE)
CD14_expr <- merge(CD14_expr, mschip_hg19[,c(1,3)], by.x="sid", by.y="Name", all.x=TRUE, sort = FALSE, no.dups = FALSE)


CD4_expr <- CD4_expr[order(CD4_expr$bh),]
CD4_expr$geneID <- unlist(lapply(strsplit(CD4_expr$gene, "_"), "[[", 1))
CD4_expr$geneName <- unlist(lapply(strsplit(CD4_expr$gene, "_"), "[[", 2))
CD8_expr <- CD8_expr[order(CD8_expr$bh),]
CD8_expr$geneID <- unlist(lapply(strsplit(CD8_expr$gene, "_"), "[[", 1))
CD8_expr$geneName <- unlist(lapply(strsplit(CD8_expr$gene, "_"), "[[", 2))
CD14_expr <- CD14_expr[order(CD14_expr$bh),]
CD14_expr$geneID <- unlist(lapply(strsplit(CD14_expr$gene, "_"), "[[", 1))
CD14_expr$geneName <- unlist(lapply(strsplit(CD14_expr$gene, "_"), "[[", 2))

CD4_expr <- CD4_expr[,c(2:13, 19:20, 1, 14:18)]
CD8_expr <- CD8_expr[,c(2:13, 19:20, 1, 14:18)]
CD14_expr <- CD14_expr[,c(2:13, 19:20, 1, 14:18)]

setwd("~/epic.neb/qtl/results/grch38_gwas_lifted/results_expr2")
write.csv(CD4_expr, "CD4_expr_corrected.csv")
write.csv(CD8_expr, "CD8_expr_corrected.csv")
write.csv(CD14_expr, "CD14_expr_corrected.csv")



my.pvalue.list<-list("CD4"=CD4_expr$bpval[!is.na(CD4_expr$bpval)], 
                     "CD8"=CD8_expr$bpval[!is.na(CD8_expr$bpval)], 
                     "CD14"=CD14_expr$bpval[!is.na(CD14_expr$bpval)])
qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))




##### sQTL results
library(readr)
setwd("~/epic.neb/qtl/results/grch38_gwas_lifted/splicing_with_PCs")

CD4_splice <- read_table2("results_cd4/CD4.results", col_names = FALSE)
colnames(CD4_splice) <- c("intron_clu","nvar","shape1","shape2","dummy","sid","dist","npval","slope","ppval","bpval")
CD8_splice <- read_table2("results_cd8/CD8.results", col_names = FALSE)
colnames(CD8_splice) <- c("intron_clu","nvar","shape1","shape2","dummy","sid","dist","npval","slope","ppval","bpval")
CD14_splice <- read_table2("results_cd14/CD14.results", col_names = FALSE)
colnames(CD14_splice) <- c("intron_clu","nvar","shape1","shape2","dummy","sid","dist","npval","slope","ppval","bpval")

SNPs_info <- read_csv("~/epic.neb/qtl/results/grch38_gwas_lifted/SNPs_info_grch38.csv")
CD4_splice <- merge(CD4_splice, SNPs_info, by.x="sid", by.y="ID", all.x=TRUE, sort = FALSE, no.dups = FALSE)
CD8_splice <- merge(CD8_splice, SNPs_info, by.x="sid", by.y="ID", all.x=TRUE, sort = FALSE, no.dups = FALSE)
CD14_splice <- merge(CD14_splice, SNPs_info, by.x="sid", by.y="ID", all.x=TRUE, sort = FALSE, no.dups = FALSE)

mschip_hg19 <- read_csv("~/epic.neb/qtl/results/hg19_gwas/MSCHIP_FullHG19BP.csv", col_types = cols(Chr = col_character()))
CD4_splice <- merge(CD4_splice, mschip_hg19[,c(1,3)], by.x="sid", by.y="Name", all.x=TRUE, sort = FALSE, no.dups = FALSE)
CD8_splice <- merge(CD8_splice, mschip_hg19[,c(1,3)], by.x="sid", by.y="Name", all.x=TRUE, sort = FALSE, no.dups = FALSE)
CD14_splice <- merge(CD14_splice, mschip_hg19[,c(1,3)], by.x="sid", by.y="Name", all.x=TRUE, sort = FALSE, no.dups = FALSE)

CD4_splice <- CD4_splice[,c(2:11, 1, 12:16)]
CD4_splice$bonferroni = p.adjust(CD4_splice$bpval, method="bonferroni")
CD4_splice$bh = p.adjust(CD4_splice$bpval, method="fdr")
CD4_splice <- CD4_splice[order(CD4_splice$bh),]

CD8_splice <- CD8_splice[,c(2:11, 1, 12:16)]
CD8_splice$bonferroni = p.adjust(CD8_splice$bpval, method="bonferroni")
CD8_splice$bh = p.adjust(CD8_splice$bpval, method="fdr")
CD8_splice <- CD8_splice[order(CD8_splice$bh),]

CD14_splice <- CD14_splice[,c(2:11, 1, 12:16)]
CD14_splice$bonferroni = p.adjust(CD14_splice$bpval, method="bonferroni")
CD14_splice$bh = p.adjust(CD14_splice$bpval, method="fdr")
CD14_splice <- CD14_splice[order(CD14_splice$bh),]


write.csv(CD4_splice, "CD4_splice_corrected.csv")
write.csv(CD8_splice, "CD8_splice_corrected.csv")
write.csv(CD14_splice, "CD14_splice_corrected.csv")

nrow(subset(CD4_splice, CD4_splice$bh < 0.05))
nrow(subset(CD8_splice, CD8_splice$bh < 0.05))
nrow(subset(CD14_splice, CD14_splice$bh < 0.05))


library(ggplot2)
qplot(sample = bpval, data = CD4_splice)
ggplot(CD4_splice, aes(sample=bpval)) + stat_qq() + stat_qq_line()
ggplot(CD14_splice, aes(sample=bpval)) + stat_qq() + stat_qq_line()

nq <- 80
p <- (1 : nq) / nq - 0.5 / nq
ggplot() + geom_point(aes(x = qnorm(p), y = quantile(CD4_splice$bpval, p, na.rm = TRUE)))
ggplot() + geom_point(aes(x = qexp(p), y = quantile(CD4_splice$bpval, p, na.rm = TRUE)))
ggplot() + geom_point(aes(x = qexp(p), y = quantile(CD8_splice$bpval, p, na.rm = TRUE)))
ggplot() + geom_point(aes(x = qexp(p), y = quantile(CD14_splice$bpval, p, na.rm = TRUE)))



#############################################################

### check overlap with 200 IMSGC SNPs
imsgc_snps <- compare_snps[paste(compare_snps$Chr,compare_snps$HG19) %in% intersect(paste(compare_snps$Chr,compare_snps$HG19), paste(imsgc_snps_list$Chr, imsgc_snps_list$pos_hg19)),]
imsgc_snps <- unique(paste(imsgc_snps$Chr,imsgc_snps$HG19)); length(imsgc_snps)

# CD4 expression
CD4_expr_sig <- subset(CD4_expr, bh < 0.05); nrow(CD4_expr_sig)
length(unique(CD4_expr$gene)); length(unique(paste(CD4_expr$CHROM,CD4_expr$POS))) # 24400 gene; 13295 snp
length(unique(CD4_expr_sig$gene)); length(unique(paste(CD4_expr_sig$CHROM,CD4_expr_sig$POS))) # 53 gene; 45 snp

length(intersect(paste(CD4_expr$CHROM, CD4_expr$HG19), imsgc_snps))  # 21 SNP included in CD4 eQTL test among 126 IMSGC SNPs
length(intersect(paste(CD4_expr_sig$CHROM, CD4_expr_sig$HG19), imsgc_snps))  # 0 SNP found in CD4 among 126 IMSGC SNPs
CD4_expr_sig[paste(CD4_expr_sig$CHROM, CD4_expr_sig$HG19) %in% intersect(paste(CD4_expr_sig$CHROM, CD4_expr_sig$HG19), imsgc_snps),]


# CD8 expression
CD8_expr_sig <- subset(CD8_expr, bh < 0.05); nrow(CD8_expr_sig)
length(unique(CD8_expr$gene)); length(unique(paste(CD8_expr$CHROM,CD8_expr$POS))) # 24400 gene; 13295 snp
length(unique(CD8_expr_sig$gene)); length(unique(paste(CD8_expr_sig$CHROM,CD8_expr_sig$POS))) # 36 gene; 28 snp

length(intersect(paste(CD8_expr$CHROM, CD8_expr$HG19), imsgc_snps))  # 21 SNP tested in CD4 among 126 IMSGC SNPs
length(intersect(paste(CD8_expr_sig$CHROM, CD8_expr_sig$HG19), imsgc_snps))  # 0 SNP found in CD4 among 126 IMSGC SNPs
CD8_expr_sig[paste(CD8_expr_sig$CHROM, CD8_expr_sig$HG19) %in% intersect(paste(CD8_expr_sig$CHROM, CD8_expr_sig$HG19), imsgc_snps),]


# CD14 expression
CD14_expr_sig <- subset(CD14_expr, bh < 0.05); nrow(CD14_expr_sig)
length(unique(CD14_expr$gene)); length(unique(paste(CD14_expr$CHROM,CD14_expr$POS))) # 45930 gene; 19699 snp
length(unique(CD14_expr_sig$gene)); length(unique(paste(CD14_expr_sig$CHROM,CD14_expr_sig$POS))) # 1937 gene; 1584 snp

length(intersect(paste(CD14_expr$CHROM, CD14_expr$HG19), imsgc_snps))  # 25 SNP tested in CD4 among 126 IMSGC SNPs
length(intersect(paste(CD14_expr_sig$CHROM, CD14_expr_sig$HG19), imsgc_snps))  # 2 SNP found in CD4 among IMSGC 126 SNPs
CD14_expr_sig[paste(CD14_expr_sig$CHROM, CD14_expr_sig$HG19) %in% intersect(paste(CD14_expr_sig$CHROM, CD14_expr_sig$HG19), imsgc_snps),]



# CD4 splicing
CD4_splice_sig <- subset(CD4_splice, bpval < 0.05)
length(unique(CD4_splice$intron_clu)); length(unique(paste(CD4_splice$CHROM,CD4_splice$POS))) # 150516 intron cluster; 34727 snp
length(unique(CD4_splice_sig$intron_clu)); length(unique(paste(CD4_splice_sig$CHROM,CD4_splice_sig$POS))) # 9453 intron cluster; 6939 snp

length(intersect(paste(CD4_splice$CHROM, CD4_splice$HG19), imsgc_snps))  # 63 SNP tested in CD4 among 126 SNPs
length(intersect(paste(CD4_splice_sig$CHROM, CD4_splice_sig$HG19), imsgc_snps))  # 10 SNP found in CD4 among 126 SNPs
CD4_splice_sig[paste(CD4_splice_sig$CHROM, CD4_splice_sig$HG19) %in% intersect(paste(CD4_splice_sig$CHROM, CD4_splice_sig$HG19), imsgc_snps),]


# CD8 splicing
CD8_splice_sig <- subset(CD8_splice, bpval < 0.05)
length(unique(CD8_splice$intron_clu)); length(unique(paste(CD8_splice$CHROM,CD8_splice$POS))) # 149207 intron cluster; 34524 snp
length(unique(CD8_splice_sig$intron_clu)); length(unique(paste(CD8_splice_sig$CHROM,CD8_splice_sig$POS))) # 9098 intron cluster; 6800 snp

length(intersect(paste(CD8_splice$CHROM, CD8_splice$HG19), imsgc_snps))  # 69 SNP tested in CD4 among 126 SNPs
length(intersect(paste(CD8_splice_sig$CHROM, CD8_splice_sig$HG19), imsgc_snps))  # 13 SNP found in CD4 among 126 SNPs
CD8_splice_sig[paste(CD8_splice_sig$CHROM, CD8_splice_sig$HG19) %in% intersect(paste(CD8_splice_sig$CHROM, CD8_splice_sig$HG19), imsgc_snps),]


# CD14 splicing
CD14_splice_sig <- subset(CD14_splice, bpval < 0.05)
length(unique(CD14_splice$intron_clu)); length(unique(paste(CD14_splice$CHROM,CD14_splice$POS))) # 138075 intron cluster; 33738 snp
length(unique(CD14_splice_sig$intron_clu)); length(unique(paste(CD14_splice_sig$CHROM,CD14_splice_sig$POS))) # 8872 intron cluster; 6500 snp

length(intersect(paste(CD14_splice$CHROM, CD14_splice$HG19), imsgc_snps))  # 63 SNP tested in CD4 among 126 SNPs
length(intersect(paste(CD14_splice_sig$CHROM, CD14_splice_sig$HG19), imsgc_snps))  # 12 SNP found in CD4 among 126 SNPs
CD14_splice_sig[paste(CD14_splice_sig$CHROM, CD14_splice_sig$HG19) %in% intersect(paste(CD14_splice_sig$CHROM, CD14_splice_sig$HG19), imsgc_snps),]



# merge significant overlap expression results
expr_all <- data.frame(CD4_expr_sig[paste(CD4_expr_sig$CHROM, CD4_expr_sig$HG19) %in% intersect(paste(CD4_expr_sig$CHROM, CD4_expr_sig$HG19), imsgc_snps),], "cell"="CD4")
expr_all <- rbind(expr_all, data.frame(CD8_expr_sig[paste(CD8_expr_sig$CHROM, CD8_expr_sig$HG19) %in% intersect(paste(CD8_expr_sig$CHROM, CD8_expr_sig$HG19), imsgc_snps),], "cell"="CD8"))
expr_all <- rbind(expr_all, data.frame(CD14_expr_sig[paste(CD14_expr_sig$CHROM, CD14_expr_sig$HG19) %in% intersect(paste(CD14_expr_sig$CHROM, CD14_expr_sig$HG19), imsgc_snps),], "cell"="CD14"))
write.csv(expr_all, "eQTL_IMSGC-overlap.csv")

# merge significant overlap splicing results
CD4_introns_all <- read_csv("C:/Users/kicheol/Desktop/Projects/EPIC_RNASeq/Exp2_NEB_FullTx/Analysis/results/leafcutter/CD4_introns_all.csv")
CD4_introns_all$area <- paste0(CD4_introns_all$chr, ":", CD4_introns_all$start, ":", CD4_introns_all$end)
CD8_introns_all <- read_csv("C:/Users/kicheol/Desktop/Projects/EPIC_RNASeq/Exp2_NEB_FullTx/Analysis/results/leafcutter/CD8_introns_all.csv")
CD8_introns_all$area <- paste0(CD8_introns_all$chr, ":", CD8_introns_all$start, ":", CD8_introns_all$end)
CD14_introns_all <- read_csv("C:/Users/kicheol/Desktop/Projects/EPIC_RNASeq/Exp2_NEB_FullTx/Analysis/results/leafcutter/CD14_introns_all.csv")
CD14_introns_all$area <- paste0(CD14_introns_all$chr, ":", CD14_introns_all$start, ":", CD14_introns_all$end)


splice_all <- data.frame(CD4_splice_sig[paste(CD4_splice_sig$CHROM, CD4_splice_sig$HG19) %in% intersect(paste(CD4_splice_sig$CHROM, CD4_splice_sig$HG19), imsgc_snps),], "cell"="CD4")
splice_all <- rbind(splice_all, data.frame(CD8_splice_sig[paste(CD8_splice_sig$CHROM, CD8_splice_sig$HG19) %in% intersect(paste(CD8_splice_sig$CHROM, CD8_splice_sig$HG19), imsgc_snps),], "cell"="CD8"))
splice_all <- rbind(splice_all, data.frame(CD14_splice_sig[paste(CD14_splice_sig$CHROM, CD14_splice_sig$HG19) %in% intersect(paste(CD14_splice_sig$CHROM, CD14_splice_sig$HG19), imsgc_snps),], "cell"="CD14"))


splice_cd4 <- data.frame(CD4_splice_sig[paste(CD4_splice_sig$CHROM, CD4_splice_sig$HG19) %in% intersect(paste(CD4_splice_sig$CHROM, CD4_splice_sig$HG19), imsgc_snps),], "cell"="CD4")
splice_cd4$area <- paste0(paste0("chr",splice_cd4$CHROM), ":", unlist(lapply(strsplit(splice_cd4$intron_clu, ":"), "[[", 2)), ":", unlist(lapply(strsplit(splice_cd4$intron_clu, ":"), "[[", 3)))
CD4_introns_all_sub <- CD4_introns_all[CD4_introns_all$area %in% splice_cd4$area,]
splice_cd4 <- merge(splice_cd4, CD4_introns_all_sub[,c(2:ncol(CD4_introns_all_sub))], by.x="area", by.y="area")

splice_cd8 <- data.frame(CD8_splice_sig[paste(CD8_splice_sig$CHROM, CD8_splice_sig$HG19) %in% intersect(paste(CD8_splice_sig$CHROM, CD8_splice_sig$HG19), imsgc_snps),], "cell"="CD8")
splice_cd8$area <- paste0(paste0("chr",splice_cd8$CHROM), ":", unlist(lapply(strsplit(splice_cd8$intron_clu, ":"), "[[", 2)), ":", unlist(lapply(strsplit(splice_cd8$intron_clu, ":"), "[[", 3)))
CD8_introns_all_sub <- CD8_introns_all[CD8_introns_all$area %in% splice_cd8$area,]
splice_cd8 <- merge(splice_cd8, CD8_introns_all_sub[,c(2:ncol(CD8_introns_all_sub))], by.x="area", by.y="area")

splice_cd14 <- data.frame(CD14_splice_sig[paste(CD14_splice_sig$CHROM, CD14_splice_sig$HG19) %in% intersect(paste(CD14_splice_sig$CHROM, CD14_splice_sig$HG19), imsgc_snps),], "cell"="CD14")
splice_cd14$area <- paste0(paste0("chr",splice_cd14$CHROM), ":", unlist(lapply(strsplit(splice_cd14$intron_clu, ":"), "[[", 2)), ":", unlist(lapply(strsplit(splice_cd14$intron_clu, ":"), "[[", 3)))
CD14_introns_all_sub <- CD14_introns_all[CD14_introns_all$area %in% splice_cd14$area,]
splice_cd14 <- merge(splice_cd14, CD14_introns_all_sub[,c(2:ncol(CD14_introns_all_sub))], by.x="area", by.y="area")

splice_all_plus <- rbind(splice_cd4, splice_cd8)
splice_all_plus <- rbind(splice_all_plus, splice_cd14)


write.csv(splice_all, "sQTL_IMSGC-overlap.csv")
write.csv(splice_all_plus, "sQTL_IMSGC-overlap_plus.csv")



### check overlapped SNPs between eQTL and sQTL

overlap_cd4 <- intersect(CD4_expr_sig$sid, CD4_splice_sig$sid)
overlap_cd8 <- intersect(CD8_expr_sig$sid, CD8_splice_sig$sid)
overlap_cd14 <- intersect(CD14_expr_sig$sid, CD14_splice_sig$sid)

# Venn Diagram
library("VennDiagram"); library("RColorBrewer"); library("gplots")
venn <- list("cd4"=overlap_cd4,           
             "cd8"=overlap_cd8, 
             "cd14"=overlap_cd14)
venn.plot <- venn.diagram(venn, NULL, resolution=600, fill=RColorBrewer::brewer.pal(3, "Dark2"), col="transparent", alpha=c(0.5,0.5,0.5), 
                          sub.cex=1, cat.cex=1.2,
                          category=c("CD4","CD8","CD14"))
grid.newpage(); grid.draw(venn.plot)

### extract list from venn diagram
tmp.list <- venn(venn, show.plot=FALSE); str(tmp.list)
tmp.inters <- attr(tmp.list,"intersections")
lapply(tmp.inters, head, n=10)
#ven.list1 <- data.frame("CD4"=tmp.inters$CD14)
overlap_cd4cd8 <- tmp.inters$`cd4:cd8`
overlap_cd4cd14 <- tmp.inters$`cd4:cd14`
overlap_cd8cd14 <- tmp.inters$`cd8:cd14`
overlap_cd4cd8cd14 <- tmp.inters$`cd4:cd8:cd14`

overlap_list <- data.frame(snp=tmp.inters$`cd4:cd8:cd14`, overlap="CD4-CD8-CD14")
overlap_list <- rbind(overlap_list, data.frame(snp=tmp.inters$`cd4:cd8`, overlap="CD4-CD8"))
overlap_list <- rbind(overlap_list, data.frame(snp=tmp.inters$`cd4:cd14`, overlap="CD4-CD14"))
overlap_list <- rbind(overlap_list, data.frame(snp=tmp.inters$`cd8:cd14`, overlap="CD8-CD14"))
overlap_list <- rbind(overlap_list, data.frame(snp=tmp.inters$cd4, overlap="CD4"))
overlap_list <- rbind(overlap_list, data.frame(snp=tmp.inters$cd8, overlap="CD8"))
overlap_list <- rbind(overlap_list, data.frame(snp=tmp.inters$cd14, overlap="CD14"))
write.csv(overlap_list, "eQTL-sQTL_overlap_snp_list.csv")


CD4_expr_overlap <- CD4_expr[CD4_expr$sid %in% overlap_cd4,]
CD4_splice_overlap <- CD4_splice[CD4_splice$sid %in% overlap_cd4,]
write.csv(CD4_expr_overlap, "CD4_expr_overlap.csv")
write.csv(CD4_splice_overlap, "CD4_splice_overlap.csv")

CD8_expr_overlap <- CD8_expr[CD8_expr$sid %in% overlap_cd8,]
CD8_splice_overlap <- CD8_splice[CD8_splice$sid %in% overlap_cd8,]
write.csv(CD8_expr_overlap, "CD8_expr_overlap.csv")
write.csv(CD8_splice_overlap, "CD8_splice_overlap.csv")

CD14_expr_overlap <- CD14_expr[CD14_expr$sid %in% overlap_cd14,]
CD14_splice_overlap <- CD14_splice[CD14_splice$sid %in% overlap_cd14,]
write.csv(CD14_expr_overlap, "CD14_expr_overlap.csv")
write.csv(CD14_splice_overlap, "CD14_splice_overlap.csv")




#####  covariates (sex and age at exam) for splicing (sQTL)
perind_overlapGWAS_gz <- read_table2("C:/Users/kicheol/Desktop/CD4_perind_overlapGWAS.gz.PCs")
perind_overlapGWAS_gz <- read_table2("C:/Users/kicheol/Desktop/CD8_perind_overlapGWAS.gz.PCs")
perind_overlapGWAS_gz <- read_table2("C:/Users/kicheol/Desktop/CD14_perind_overlapGWAS.gz.PCs")
perind_overlapGWAS_gz <- data.frame(t(perind_overlapGWAS_gz))

gene_covar <- read_table2("C:/Users/kicheol/Desktop/gene_covar_CD4.txt")
gene_covar <- read_table2("C:/Users/kicheol/Desktop/gene_covar_CD8.txt")
gene_covar <- read_table2("C:/Users/kicheol/Desktop/gene_covar_CD14.txt")
gene_covar <- as.data.frame(t(gene_covar))
gene_covar$id <- row.names(gene_covar)
colnames(gene_covar) <- c("Sex","AgeAtExam","id")
gene_covar <- gene_covar[2:30,c(3,1,2)]

gene_covar <- gene_covar[row.names(perind_overlapGWAS_gz[2:30,]),]
gene_covar <- t(gene_covar)
write.table(gene_covar, file="CD14_perind_overlapGWAS.covar.txt", quote=FALSE, sep="\t", col.names=FALSE)
write.table(gene_covar, file="CD8_perind_overlapGWAS.covar.txt", quote=FALSE, sep="\t", col.names=FALSE)
write.table(gene_covar, file="CD4_perind_overlapGWAS.covar.txt", quote=FALSE, sep="\t", col.names=FALSE)






#############################################################################################################################
##########  QQplot function by Matthew Flickinger (https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R)
#############################################################################################################################
library(lattice)
qqunif.plot<-function(pvalues, 
                      should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                      xlab=expression(paste("Expected (",-log[10], " p-value)")),
                      ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                      draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                      already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                      par.settings=list(superpose.symbol=list(pch=pch)), ...) {
  
  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed==FALSE) {
    if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }
  
  
  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }
  
  
  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    require(grid)
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
    }
    grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
  }
  
  #reduce number of points to plot
  if (should.thin==T) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()
  
  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }
  
  #draw the plot
  xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
         prepanel=prepanel, scales=list(axs="i"), pch=pch,
         panel = function(x, y, ...) {
           if (draw.conf) {
             panel.qqconf(n, conf.points=conf.points, 
                          conf.col=conf.col, conf.alpha=conf.alpha)
           };
           panel.xyplot(x,y, ...);
           panel.abline(0,1);
         }, par.settings=par.settings, ...
  )
}
