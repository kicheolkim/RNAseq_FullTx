### function for permutation results
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
