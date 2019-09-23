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



#####
rm(cds, res, res.sig)
sessionInfo()
#####
