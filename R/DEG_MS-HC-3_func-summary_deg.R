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
