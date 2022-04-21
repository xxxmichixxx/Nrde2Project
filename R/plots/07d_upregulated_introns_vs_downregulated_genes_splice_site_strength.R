rm(list=ls())

library(tidyverse)
library(Rsubread)
library(GenomicFeatures)
library(GenomicRanges)
library(MiniChip)
library(patchwork)

setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")

#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
#   #load upregulated introns results 
#----------------------------------------------------------------------------------------------------------------

#load intronic read differential expression
#normal
res2 <- read.table(file="FigureData/Nrde2_Mtr4_Ccdc174_batch1and2_intronicRNAseq_DEseq2_results_0.05_0.5FC_5readspersample_seperateWTs.txt",sep="\t",header=TRUE)
#plus Cyclohexamide (CHX)
res2C <- read.table(file="FigureData/CHX_Nrde2_Ccdc174_intronicRNAseq_DEseq2_results_0.05_0.5FC_5readspersample_seperateWTs.txt",sep="\t",header=TRUE)
#RIBOseq
res2R <- read.table(file="FigureData/Nrde2_intronicRIBOseq_DEseq2_results_0.05_0.5FC_5readspersample.txt",sep="\t",header=TRUE)
res2R$Contrast <- paste(res2R$Contrast,"Ribo",sep="_")

#CHX WT vs WT
res2W <- read.table(file="FigureData/CHX_WT_intronicRNAseq_DEseq2_results_0.05_0.5FC_5readspersample.txt",sep="\t",header=TRUE)

#combine all contrasts
res2 <- rbind(res2,res2C,res2R,res2W)

#add GeneID column
load("splice_sites_expressed_and_notexpressed_transcripts_log2TPMabove0_and_genes_FPKMabove0.RData")
names(spliceDonors) <- paste(names(spliceDonors),start(spliceDonors),sep="_")
spliceDonors.df <- data.frame(spliceDonors)
spliceDonors.df$ID <- names(spliceDonors)

res2 <- left_join(res2,spliceDonors.df,by="ID")

#----------------------------------------------------------------------------------------------------------------
#   #add the gene expression DEseq data 
#-----------------------------------------------------------------------------------------------------------------
#normal
#res1 <- read.table("FigureData/Nrde2_Mtr4_Ccdc174_RNAseq_DEseq2_results.txt",sep="\t",header=TRUE)
res1 <- read.table("FigureData/Nrde2_Mtr4_Ccdc174_RNAseq_DEseq2_results_v2_cutoff200.txt",sep="\t",header=TRUE)

#plus Cyclohexamide (CHX)
res1C <- read.table("FigureData/CHX_Nrde2_Ccdc174_RNAseq_DEseq2_results_cutoff80.txt",sep="\t",header=TRUE)

#Riboseq
res1R <- read.table("FigureData/Nrde2_RIBOseq_DEseq2_results_v2_cutoff10.txt",sep="\t",header=TRUE)
res1R$Contrast <- paste(res1R$Contrast,"Ribo",sep="_")

#CHX WT vs WT
res1W <- read.table(file="FigureData/CHX_WT_RNAseq_DEseq2_results_cutoff80.txt",sep="\t",header=TRUE)

#combine all contrasts
res1 <- rbind(res1,res1C,res1R,res1W)




res2$GeneID <- matrix(unlist(strsplit(as.character(res2$GeneID),".",fixed=TRUE)),ncol=2,byrow = TRUE)[,1]
res2$GeneConID <- paste(res2$GeneID,res2$Contrast,sep="_")
res1$GeneConID <- paste(res1$GeneID,res1$Contrast,sep="_")

res3 <- left_join(res2,res1,by="GeneConID")

res3[res3$GeneID.x=="ENSMUSG00000025358" & res3$intron_regulated=="up",]
res3[res3$GeneID.x=="ENSMUSG00000004099" & res3$intron_regulated=="up",]

#x <- filter(res3,Contrast.x==contrasts[i],intron_regulated=="up",log2FoldChange.y < 0)
#----------------------------------------------------------------------------------------------------------------
#   #add the peaks and MATT results
#-----------------------------------------------------------------------------------------------------------------
matt.res <- read.table("FigureData/MATT_scores_on_all_introns_with_peak_overlap_info.txt",header=TRUE,sep="\t")
res4 <- left_join(res3,matt.res,by="ID")

res4[res4$GeneID.x=="ENSMUSG00000025358" & res4$intron_regulated=="up",]
res4[res4$GeneID.x=="ENSMUSG00000004099" & res3$intron_regulated=="up",]

write.table(res4,file="FigureData/Introns_up_data_with_gene_expression_and_peaks_and_MATTresults_With_WTCHX.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE, append=FALSE)
#----------------------------------------------------------------------------------------------------------------
#   #select introns with nrde2 peak
#-----------------------------------------------------------------------------------------------------------------
contrasts <- unique(res4$Contrast.x)[c(1,2,5:9)]

intron_up_list <- list()
intron_notup_list <- list()

for (i in seq_along(contrasts)){
#select single contrast RNAseq data with Nrde2 peaks
#res5 <- dplyr::filter(res4,Nrde2.peak=="yes",Contrast.x==contrasts[i])
#OR: select single contrast RNAseq data with Ccdc174 peaks
res5 <- dplyr::filter(res4,Ccdc174.peak=="yes",Contrast.x==contrasts[i])

#select the not upregulated genes (all genes with Fc< 0)
res6 <- dplyr::filter(res5,log2FoldChange.y < 0)
table(res6$intron_regulated)
res6$intron_regulated2 <- ifelse(res6$intron_regulated=="up","up","not_up")
intron_up_list[[i]] <- res6$ID[res6$intron_regulated2=="up"]
intron_notup_list[[i]] <- res6$ID[res6$intron_regulated2=="not_up"]

plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A")
#plotcols2 <- c("darkgrey",plotcols[1])
plotcols2 <- c("darkgrey",plotcols[5])

pval1 <- wilcox.test(res6$MAXENTSCR_HSAMODEL_5SS[res6$intron_regulated2=="up"],res6$MAXENTSCR_HSAMODEL_5SS[res6$intron_regulated2=="not_up"])
p1 <- ggplot(res6,aes(x=MAXENTSCR_HSAMODEL_5SS,fill=intron_regulated2,color=intron_regulated2)) +geom_density(alpha=0.05) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2) + xlim(c(-5,15)) + ggtitle(pval1$p.value)

pval2 <- wilcox.test(res6$MAXENTSCR_HSAMODEL_3SS[res6$intron_regulated2=="up"],res6$MAXENTSCR_HSAMODEL_3SS[res6$intron_regulated2=="not_up"])
p2 <- ggplot(res6,aes(x=MAXENTSCR_HSAMODEL_3SS,fill=intron_regulated2,color=intron_regulated2)) +geom_density(alpha=0.05) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2) + xlim(c(-5,15)) + ggtitle(pval2$p.value)

pval3 <- wilcox.test(res6$BPSCORE_MAXBP[res6$intron_regulated2=="up"],res6$BPSCORE_MAXBP[res6$intron_regulated2=="not_up"])
p3 <- ggplot(res6,aes(x=BPSCORE_MAXBP,fill=intron_regulated2,color=intron_regulated2)) +geom_density(alpha=0.05) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2) + ggtitle(pval3$p.value)

p4 <- ggplot(res6,aes(x=INTRON_GCC,fill=intron_regulated2,color=intron_regulated2)) +geom_density(alpha=0.05) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)

p5 <- ggplot(res6,aes(x=log2(INTRON_LENGTH),fill=intron_regulated2,color=intron_regulated2)) +geom_density(alpha=0.05) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)

p6 <- ggplot(res6,aes(x=RATIO_UPEXON_INTRON_GCC,fill=intron_regulated2,color=intron_regulated2)) +geom_density(alpha=0.05) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)

p7 <- ggplot(res6,aes(x=RATIO_UPEXON_INTRON_LENGTH,fill=intron_regulated2,color=intron_regulated2)) +geom_density(alpha=0.05) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)

p1 + p2 + p3 #+p4 + p5 + p6 + p7
#ggsave(sprintf("Figures/MATT_scores_in_upregulated_vs_not_upregulated_introns_in_contrast_%s_within_Nrde2_peak_introns_within_downregulated_genes_210427.pdf",contrasts[i]),device="pdf",height=4,width=16)
ggsave(sprintf("Figures/MATT_scores_in_upregulated_vs_not_upregulated_introns_in_contrast_%s_within_Ccdc174_peak_introns_within_downregulated_genes_210614.pdf",contrasts[i]),device="pdf",height=4,width=16)

}

names(intron_up_list) <- contrasts
names(intron_notup_list) <- contrasts

#overkaps between upregulated introns
sapply(intron_up_list,length)

#Nrde2 and Ccdc174
length(which(intron_up_list[[2]] %in% intron_up_list[[4]]))
#Nrde2 and Mtr4
length(which(intron_up_list[[2]] %in% intron_up_list[[1]]))
#Ccdc174 and Mtr4
length(which(intron_up_list[[4]] %in% intron_up_list[[1]]))
#Nrde2 CHX vs Nrde2 normal
length(which(intron_up_list[[2]] %in% intron_up_list[[5]]))
#Ccdc174 CHX vs Nrde2 normal
length(which(intron_up_list[[4]] %in% intron_up_list[[6]]))
#Nrde2 CHX vsCcdc174 CHX
length(which(intron_up_list[[5]] %in% intron_up_list[[6]]))
#Nrde2 CHX vs Mtr4 normal
length(which(intron_up_list[[5]] %in% intron_up_list[[1]]))
#ccdc174 CHX vs Mtr4 normal
length(which(intron_up_list[[6]] %in% intron_up_list[[1]]))
#Nrde2 CHX vs Ccdc174 normal
length(which(intron_up_list[[5]] %in% intron_up_list[[4]]))
#Ccdc174 CHX vs Nrde2 normal
length(which(intron_up_list[[6]] %in% intron_up_list[[2]]))

#save the intron identifiers as table
for (i in seq_along(intron_up_list)){
  intron_up_list[[i]] <- data.frame(ID=intron_up_list[[i]],contrast=names(intron_up_list)[i])
}
intron_up_df <- do.call("rbind",intron_up_list)
#write.table(intron_up_df,"FigureData/intron_IDs_of_upregulated_introns_in_downregulated_genes_per_contrast_Nrde2_peaks_210427.txt",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
write.table(intron_up_df,"FigureData/intron_IDs_of_upregulated_introns_in_downregulated_genes_per_contrast_Ccdc174_peaks_210427.txt",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

for (i in seq_along(intron_notup_list)){
  intron_notup_list[[i]] <- data.frame(ID=intron_notup_list[[i]],contrast=names(intron_notup_list)[i])
}
intron_notup_df <- do.call("rbind",intron_notup_list)
#write.table(intron_notup_df,"FigureData/intron_IDs_of_not_upregulated_introns_in_downregulated_genes_per_contrast_Nrde2_peaks_210427.txt",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
write.table(intron_notup_df,"FigureData/intron_IDs_of_not_upregulated_introns_in_downregulated_genes_per_contrast_Ccdc174_peaks_210427.txt",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)



#----------------------------------------------------------------------------------------------------------------
#   #compare single vs double contrast up introns
#-----------------------------------------------------------------------------------------------------------------
res5 <- dplyr::filter(res4,Nrde2.peak=="yes",Contrast.x=="Ccdc174KD_vs_WT")

#select the not upregulated genes (all genes with Fc< 0)
res6 <- dplyr::filter(res5,log2FoldChange.y < 0)
table(res6$intron_regulated)
res6$intron_regulated2 <- ifelse(res6$intron_regulated=="up","up","not_up")

#select the Nrde2 up introns
res5N <- dplyr::filter(res4,Nrde2.peak=="yes",Contrast.x=="Nrde2KO_vs_WT")
res6N <- dplyr::filter(res5N,log2FoldChange.y < 0)
res6N$intron_regulated2 <- ifelse(res6N$intron_regulated=="up","up","not_up")

#split the Ccdc174 up introns in Nrde2 up and not up
res6$both_up_introns <- ifelse(res6$intron_regulated2=="up" & res6$ID %in% res6N$ID[res6N$intron_regulated2=="up"]==TRUE,"both_up",
                          ifelse(res6$intron_regulated2=="up" & res6$ID %in% res6N$ID[res6N$intron_regulated2=="up"]==FALSE,"Ccdc174_up", "other"))   
table(res6$both_up_introns)


plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A")
plotcols2 <- c("darkgrey",plotcols[1])

ggplot(res6,aes(x=MAXENTSCR_HSAMODEL_5SS,fill=both_up_introns,color=both_up_introns)) +geom_density(alpha=0.05) + 
  theme_classic()  + xlim(c(-15,15))


#----------------------------------------------------------------------------------------------------------------
#   #consolidate to genes and plot the logFC of gene expression
#-----------------------------------------------------------------------------------------------------------------
contrasts <- unique(res4$Contrast.x)[c(1,2,5:9)]
plotlist <- list()

for (i in seq_along(contrasts)){
  #select single contrast RNAseq data with Nrde2 peaks
  #res5 <- dplyr::filter(res4,Nrde2.peak=="yes",Contrast.x==contrasts[i])
  #select single contrast RNAseq data with Ccdc174 peaks
  res5 <- dplyr::filter(res4,Ccdc174.peak=="yes",Contrast.x==contrasts[i])
  
  #select the not upregulated genes (all genes with Fc< 0)
  res6 <- dplyr::filter(res5,log2FoldChange.y < 0)
  table(res6$intron_regulated)
  res6$intron_regulated2 <- ifelse(res6$intron_regulated=="up","up","not_up")

  genes_with_upregulated_introns <- unique(res6$GeneID.x[res6$intron_regulated=="up"])
  genes_with_not_upregulated_introns <- unique(res6$GeneID.x[res6$intron_regulated!="up"])
  genes_with_not_upregulated_introns <- genes_with_not_upregulated_introns[genes_with_not_upregulated_introns %in% genes_with_upregulated_introns == FALSE]

  res1.Nrde2 <- dplyr::filter(res1,Contrast==contrasts[i])
  res1.Nrde2$introns_reg <- ifelse(res1.Nrde2$GeneID %in% genes_with_upregulated_introns,"up",
                                 ifelse(res1.Nrde2$GeneID %in% genes_with_not_upregulated_introns,"not.up","other"))
table(res1.Nrde2$introns_reg)
res1.Nrde2 <- res1.Nrde2[res1.Nrde2$introns_reg !="other",]

plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A")
#plotcols2 <- c("darkgrey",plotcols[1])
plotcols2 <- c("darkgrey",plotcols[5])

#plotlist[[i]] <- ggplot(res1.Nrde2,aes(x=log2FoldChange,col=introns_reg,fill=introns_reg)) +geom_density(alpha=0.05) + 
 # theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2) + ggtitle(contrasts[i])

pval <- wilcox.test(res1.Nrde2$log2FoldChange[res1.Nrde2$introns_reg=="up"],res1.Nrde2$log2FoldChange[res1.Nrde2$introns_reg=="not.up"])
sample_size = res1.Nrde2 %>% group_by(introns_reg) %>% summarize(num=n())

#plotlist[[i]] <- ggplot(res1.Nrde2,aes(y=log2FoldChange,x=introns_reg,fill=introns_reg)) +geom_violin(alpha=0.25) + geom_boxplot(width=0.1, color="black", alpha=0.8) +
#  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2) + ggtitle(contrasts[i],subtitle = pval$p.value) 

#plotlist[[i]] <- res1.Nrde2 %>%
#  left_join(sample_size) %>%
#  mutate(introns_reg2 = paste0(introns_reg, "\n", "n=", num)) %>% 
 # ggplot(aes(y=log2FoldChange,x=introns_reg2,fill=introns_reg)) +geom_violin(alpha=0.25) + geom_boxplot(width=0.1, color="black", alpha=0.8) +
#  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2) + ggtitle(contrasts[i],subtitle = pval$p.value) 

#boxplot
plotlist[[i]] <- res1.Nrde2 %>%
  left_join(sample_size) %>%
  mutate(introns_reg2 = paste0(introns_reg, "\n", "n=", num)) %>% 
  ggplot(aes(y=log2FoldChange,x=introns_reg2,fill=introns_reg))  + geom_boxplot(color="black",outlier.shape=NA) +
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2) + ggtitle(contrasts[i],subtitle = pval$p.value) +
  scale_y_continuous(limits=c(-1,0))


}

plotlist[[1]] + plotlist[[2]] + plotlist[[3]] + plotlist[[4]] + plotlist[[5]] + plotlist[[6]] 
#ggsave("Figures/Exon_logFC_in_upregulated_vs_not_upregulated_introns_within_Nrde2_peak_introns_within_downregulated_genes_boxplot_210615.pdf",device="pdf",height=10,width=10)
ggsave("Figures/Exon_logFC_in_upregulated_vs_not_upregulated_introns_within_Ccdc174_peak_introns_within_downregulated_genes_boxplot_210615.pdf",device="pdf",height=10,width=10)

#----------------------------------------------------------------------------------------------------------------
#   #consolidate to genes,select the upregulated introns absed on RNAseq, and plot the logFC of Riboseq
#-----------------------------------------------------------------------------------------------------------------
contrasts <- unique(res4$Contrast.x)[c(1,2,5:9)]
plotlist <- list()

for (i in seq_along(contrasts)){
  #select single contrast RNAseq data with Nrde2 peaks
   res5 <- dplyr::filter(res4,Nrde2.peak=="yes",Contrast.x==contrasts[i])
  #select single contrast RNAseq data with Ccdc174 peaks
  res5 <- dplyr::filter(res4,Ccdc174.peak=="yes",Contrast.x==contrasts[i])
  
  #select the not upregulated genes (all genes with Fc< 0)
  res6 <- dplyr::filter(res5,log2FoldChange.y < 0)
  table(res6$intron_regulated)
  res6$intron_regulated2 <- ifelse(res6$intron_regulated=="up","up","not_up")
  
  genes_with_upregulated_introns <- unique(res6$GeneID.x[res6$intron_regulated=="up"])
  genes_with_not_upregulated_introns <- unique(res6$GeneID.x[res6$intron_regulated!="up"])
  genes_with_not_upregulated_introns <- genes_with_not_upregulated_introns[genes_with_not_upregulated_introns %in% genes_with_upregulated_introns == FALSE]
  
  res1.Nrde2 <- dplyr::filter(res1,Contrast==contrasts[i])
  res1.Nrde2$introns_reg <- ifelse(res1.Nrde2$GeneID %in% genes_with_upregulated_introns,"up",
                                   ifelse(res1.Nrde2$GeneID %in% genes_with_not_upregulated_introns,"not.up","other"))
  table(res1.Nrde2$introns_reg)
  res1.Nrde2 <- res1.Nrde2[res1.Nrde2$introns_reg !="other",]
  #add RIBOseq FCs
  RIBO_res <- dplyr::filter(res1,Contrast=="Nrde2KO_vs_WT_Ribo") %>% dplyr::select("GeneID","log2FoldChange")
  res1.Nrde2 <- left_join(res1.Nrde2,RIBO_res,by="GeneID")
  
  plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A")
 plotcols2 <- c("darkgrey",plotcols[1])
 plotcols2 <- c("darkgrey",plotcols[5])
  
 
  
  pval <- wilcox.test(res1.Nrde2$log2FoldChange.y[res1.Nrde2$introns_reg=="up"],res1.Nrde2$log2FoldChange.y[res1.Nrde2$introns_reg=="not.up"])
  sample_size = res1.Nrde2 %>% group_by(introns_reg) %>% summarize(num=n())
  
  #plotlist[[i]] <- ggplot(res1.Nrde2,aes(y=log2FoldChange,x=introns_reg,fill=introns_reg)) +geom_violin(alpha=0.25) + geom_boxplot(width=0.1, color="black", alpha=0.8) +
  #  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2) + ggtitle(contrasts[i],subtitle = pval$p.value) 
  
 #violin
  #plotlist[[i]] <- res1.Nrde2 %>%
   # left_join(sample_size) %>%
   # mutate(introns_reg2 = paste0(introns_reg, "\n", "n=", num)) %>% 
   # ggplot(aes(y=log2FoldChange.y,x=introns_reg2,fill=introns_reg)) +geom_violin(alpha=0.25) + geom_boxplot(width=0.1, color="black", alpha=0.8) +
   # theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2) + ggtitle(contrasts[i],subtitle = pval$p.value) 
  
 #boxplot
  plotlist[[i]] <- res1.Nrde2 %>%
    left_join(sample_size) %>%
    mutate(introns_reg2 = paste0(introns_reg, "\n", "n=", num)) %>% 
    ggplot(aes(y=log2FoldChange.y,x=introns_reg2,fill=introns_reg))  + geom_boxplot(color="black") +
    theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2) + ggtitle(contrasts[i],subtitle = pval$p.value) 
  
  
}

plotlist[[1]] + plotlist[[2]] + plotlist[[3]] + plotlist[[4]] + plotlist[[5]] + plotlist[[6]] 
#ggsave("Figures/Exon_logFC_in_upregulated_vs_not_upregulated_introns_within_Nrde2_peak_introns_within_downregulated_genes_boxplot_with_Ribo_210427.pdf",device="pdf",height=10,width=10)
ggsave("Figures/Exon_logFC_in_upregulated_vs_not_upregulated_introns_within_Ccdc174_peak_introns_within_downregulated_genes_boxplot_with_Ribo_210427.pdf",device="pdf",height=10,width=10)

#----------------------------------------------------------------------------------------------------------------
#   #RMATS analysis of upregulated Nrde2 bound introns
#-----------------------------------------------------------------------------------------------------------------

comp <- c("Nrde2KO_vs_WT","Nrde2d200_vs_WT","Nrde2D174R_vs_WT","Nrde2KO_vs_WT_CHX","Ccdc174KO_vs_WT","Ccdc174KO_vs_WT_CHX","Mtr4KO_vs_WT","WT_CHX_vs_WT")
RIlist <- list()
A5SSlist <- list()
A3SSlist <- list()

for (i in seq_along(comp)){
  #retained introns
  RI.MATS.JCEC <- read.table(sprintf("rmats_analysis/%s/RI.MATS.JCEC.txt",comp[i]),header=TRUE)
  RI.MATS.JCECsig <- RI.MATS.JCEC %>% dplyr::filter(FDR < 0.01,IncLevelDifference > 0)
  #  RI.MATS.JCECsig <- RI.MATS.JCEC %>% dplyr::filter(FDR < 0.01)
  RI.MATS.JCECsig$comp <- comp[i]
  RIlist[[i]] <- RI.MATS.JCECsig
  #5'ASS
  A5SS.MATS.JCEC <- read.table(sprintf("rmats_analysis/%s/A5SS.MATS.JCEC.txt",comp[i]),header=TRUE)
  A5SS.MATS.JCECsig <- A5SS.MATS.JCEC %>% dplyr::filter(FDR < 0.01,IncLevelDifference > 0)
  # A5SS.MATS.JCECsig <- A5SS.MATS.JCEC %>% dplyr::filter(FDR < 0.01)
  A5SS.MATS.JCECsig$comp <- comp[i]
  A5SSlist[[i]] <- A5SS.MATS.JCECsig
  #3'ASS
  A3SS.MATS.JCEC <- read.table(sprintf("rmats_analysis/%s/A3SS.MATS.JCEC.txt",comp[i]),header=TRUE)
  A3SS.MATS.JCECsig <- A3SS.MATS.JCEC %>% dplyr::filter(FDR < 0.01,IncLevelDifference > 0)
  # A3SS.MATS.JCECsig <- A3SS.MATS.JCEC %>% dplyr::filter(FDR < 0.01)
  A3SS.MATS.JCECsig$comp <- comp[i]
  A3SSlist[[i]] <- A3SS.MATS.JCECsig
}

RIres <- do.call(rbind,RIlist)
A5SSres <- do.call(rbind,A5SSlist)
A3SSres <- do.call(rbind,A3SSlist)

RIres$chr_5SS_start_ID <- ifelse(RIres$strand=="+",paste(RIres$chr,(RIres$upstreamEE+1),sep="_"),
                                 paste(RIres$chr,(RIres$downstreamES),sep="_"))
A5SSres$chr_5SS_start_ID <- ifelse(A5SSres$strand=="+",paste(A5SSres$chr,(A5SSres$shortEE+1),sep="_"),
                                 paste(A5SSres$chr,(A5SSres$shortES),sep="_"))
A3SSres$chr_5SS_start_ID <- ifelse(A3SSres$strand=="+",paste(A3SSres$chr,(A3SSres$flankingEE+1),sep="_"),
                                   paste(A3SSres$chr,(A3SSres$flankingES),sep="_"))

#compare to total number of splice sites
spliceDonors$chr_5SS_start_ID <- paste(seqnames(spliceDonors),start(spliceDonors),sep="_")
spliceDonors$intron_retention <- ifelse(spliceDonors$chr_5SS_start_ID %in% RIres$chr_5SS_start_ID,"RMATS_Intron_retention","no_RMATS_Intron_retention")
table(spliceDonors$intron_retention)
spliceDonors$SS5 <- ifelse(spliceDonors$chr_5SS_start_ID %in% A5SSres$chr_5SS_start_ID,"RMATS_Intron_retention","no_RMATS_Intron_retention")
table(spliceDonors$SS5)
spliceDonors$SS3 <- ifelse(spliceDonors$chr_5SS_start_ID %in% A3SSres$chr_5SS_start_ID,"RMATS_Intron_retention","no_RMATS_Intron_retention")
table(spliceDonors$SS3)



#select the nrde2 bound Nrde2 upreulated introns

contrasts <- unique(res4$Contrast.x)[c(1,2,5:8,10)]


for (i in seq_along(contrasts)){
  #select single contrast RNAseq data with Nrde2 peaks
  res5 <- dplyr::filter(res4,Nrde2.peak=="yes",Contrast.x==contrasts[i])
  #OR: select single contrast RNAseq data with Ccdc174 peaks
 # res5 <- dplyr::filter(res4,Ccdc174.peak=="yes",Contrast.x==contrasts[i])
  
  #select the not upregulated genes (all genes with Fc< 0)
  res6 <- dplyr::filter(res5,log2FoldChange.y < 0)
 # table(res6$intron_regulated)
  res6$intron_regulated2 <- ifelse(res6$intron_regulated=="up","up","not_up")
  #make an ID to match the RMATS data (start and end are the same in res6)
  res6$chr_5SS_start_ID <- paste(res6$seqnames.x,res6$start.x,sep="_")
  
  intron_retention <- matrix(ncol=length(comp),nrow=nrow(res6))
  ASS5 <- matrix(ncol=length(comp),nrow=nrow(res6))
  ASS3 <- matrix(ncol=length(comp),nrow=nrow(res6))
  
  for (o in seq_along(comp)){
  #add intron retention info
  intron_retention[,o] <- ifelse(res6$chr_5SS_start_ID %in% RIres$chr_5SS_start_ID[RIres$comp==comp[o]],"RMATS_Intron_retention","no_RMATS_Intron_retention")
  #add new 5'SS info
  ASS5[,o] <- ifelse(res6$chr_5SS_start_ID %in% A5SSres$chr_5SS_start_ID[A5SSres$comp==comp[o]],"RMATS_new5SS","no_RMATS_new5SS")
  #add new 5'SS info
  ASS3[,o] <- ifelse(res6$chr_5SS_start_ID %in% A3SSres$chr_5SS_start_ID[A3SSres$comp==comp[o]],"RMATS_new3SS","no_RMATS_new3SS")
  }
  colnames(intron_retention) <- paste("IR",comp,sep=".")
  colnames(ASS5) <- paste("ASS5",comp,sep=".")
  colnames(ASS3) <- paste("ASS3",comp,sep=".")
  res7 <- cbind(res6,intron_retention,ASS5,ASS3)

  #save table for current contrast
  
write.table(res7,file=sprintf("FigureData/selected_introns_with_Nrde2_peak_geneFCbelow0_with_intron_up_info_for_%s_and_allRMATS_info.txt",contrasts[i]),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE ) 
#write.table(res7,file=sprintf("FigureData/selected_introns_with_Ccdc174_peak_geneFCbelow0_with_intron_up_info_for_%s_and_allRMATS_info.txt",contrasts[i]),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE ) 

  
}
  

#plot (optional)
ASS.counts <- data.frame(RI=table(RIres$comp),
                         A5SS=table(A5SSres$comp),
                         A3SS=table(A3SSres$comp)) %>% pivot_longer(cols = c("RI.Freq","A5SS.Freq","A3SS.Freq"))

ggplot(ASS.counts,aes(x=RI.Var1, y=value,fill=name)) + geom_bar(stat="identity",position = "dodge") + theme_classic()

#----------------------------------------------------------------------------------------------------------------
#   # compare +/- RMATS  FPKMs for upregulated Nrde2 bound introns
#-----------------------------------------------------------------------------------------------------------------
#find RNAseq bam files
all.bamFiles <- list.files("/work2/gbuehler/deepSeqRepos/bam/", full.names=TRUE,pattern="*bam$")
#bamFiles <- c(grep("2191F9",all.bamFiles,value=TRUE),grep("2191F10",all.bamFiles,value=TRUE),grep("2191F11",all.bamFiles,value=TRUE),grep("2191F12",all.bamFiles,value=TRUE))
bamFilesR <- c(grep("2191F",all.bamFiles,value=TRUE))
bamFilesR <- c(grep("spliced",bamFilesR,value=TRUE))
bamFilesR <- c(grep("Nxf1",bamFilesR,value=TRUE,invert=TRUE))
bamFilesR2 <- grep("2447F",all.bamFiles,value=TRUE)
bamFilesR <- c(bamFilesR,bamFilesR2)


#select Nrde2 KO
bamFiles <- c(grep("Nrde2KO",bamFilesR,value=TRUE),grep("Nrde2-KO",bamFilesR,value=TRUE))
bamFiles <- c(grep("CHX",bamFiles,value=TRUE,invert=FALSE))
bamNames <- gsub("/work2/gbuehler/deepSeqRepos/bam//","",bamFiles)
bamNames <- gsub("_Aligned.sortedByCoord.out.bam","",bamNames)

res7 <- read.table("FigureData/selected_introns_with_Nrde2_peak_geneFCbelow0_with_intron_up_info_for_Nrde2KOCHX_vs_WTCHX_and_allRMATS_info.txt",sep="\t",header=TRUE)
res8 <- filter(res7,intron_regulated2=="up")

names(introns3g) <- names(spliceDonors)
introns3gsel <- introns3g[names(introns3g) %in% res8$ID]
introns.saf <- data.frame(GeneID= names(introns3gsel), Chr=seqnames(introns3gsel),
                          Start=start(introns3gsel), End=end(introns3gsel),Strand=strand(introns3gsel))

#calculate the number of reads overlapping introns
counts.introns <- featureCounts(bamFiles,annot.ext=introns.saf,
                                useMetaFeatures=FALSE,allowMultiOverlap=TRUE,
                                minOverlap=10,countMultiMappingReads=FALSE,fraction=FALSE,
                                minMQS=255,strandSpecific=2,nthreads=10,verbose=FALSE,isPairedEnd=FALSE)

counts.introns2 <- data.frame(counts.introns$counts)
colnames(counts.introns2) <- bamNames

#calculate CPMs and FPKMs
fpkm <- counts.introns$counts
cpm <- counts.introns$counts

for (i in seq_along(bamFiles)){
  scaling_factor <-  sum(counts.introns$stat[counts.introns$stat$Status=="Assigned" | counts.introns$stat$Status=="Unassigned_NoFeatures",i+1])/1e9
  cpm[,i] <- counts.introns$counts[,i]/scaling_factor
  fpkm[,i] <- cpm[,i]/(counts.introns$annotation$Length/1000)
}

cpm <- data.frame(cpm)
colnames(cpm) <- bamNames

fpkm <- data.frame(fpkm)
colnames(fpkm) <- bamNames
fpkm$Nrde2KO.FPKM <- rowSums(fpkm)
fpkm$ID <- row.names(fpkm)

res9 <- left_join(res8,fpkm[,3:4],by="ID")

plot(log2(res9$INTRON_LENGTH[res9$ASS5.Nrde2KO_vs_WT_CHX=="no_RMATS_new5SS"]),log2(res9$Nrde2KO.FPKM[res9$ASS5.Nrde2KO_vs_WT_CHX=="no_RMATS_new5SS"]+1),pch=20,col="grey")
points(log2(res9$INTRON_LENGTH[res9$ASS5.Nrde2KO_vs_WT_CHX=="RMATS_new5SS"]),log2(res9$Nrde2KO.FPKM[res9$ASS5.Nrde2KO_vs_WT_CHX=="RMATS_new5SS"]+1),pch=20,col="red")

