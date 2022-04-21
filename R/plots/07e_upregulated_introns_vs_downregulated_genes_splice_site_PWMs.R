rm(list=ls())

library(tidyverse)
library(Rsubread)
library(GenomicFeatures)
library(GenomicRanges)
library(MiniChip)
library(patchwork)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TFBSTools)

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

#----------------------------------------------------------------------------------------------------------------
#   #add the peaks and MATT results
#-----------------------------------------------------------------------------------------------------------------
matt.res <- read.table("FigureData/MATT_scores_on_all_introns_with_peak_overlap_info.txt",header=TRUE,sep="\t")
res4 <- left_join(res3,matt.res,by="ID")

#----------------------------------------------------------------------------------------------------------------
#   #resize the spliceDonor sites to et  -3 to +6 bp around junctions (donors)
#----------------------------------------------------------------------------------------------------------------
spliceDonors.plus <- spliceDonors[strand(spliceDonors)=="+"]
start(spliceDonors.plus) <- start(spliceDonors.plus) -3
end(spliceDonors.plus) <- start(spliceDonors.plus) +8
getSeq(BSgenome.Mmusculus.UCSC.mm10,spliceDonors.plus)

spliceDonors.minus <- spliceDonors[strand(spliceDonors)=="-"]
start(spliceDonors.minus) <- start(spliceDonors.minus) -5
end(spliceDonors.minus) <- start(spliceDonors.minus) +8
getSeq(BSgenome.Mmusculus.UCSC.mm10,spliceDonors.minus)

spliceDonors8n <- c(spliceDonors.plus,spliceDonors.minus)


#----------------------------------------------------------------------------------------------------------------
#   #select introns with nrde2 peak
#-----------------------------------------------------------------------------------------------------------------
contrasts <- unique(res4$Contrast.x)[c(1,2,5:9)]



for (i in seq_along(contrasts)){
#select single contrast RNAseq data with Nrde2 peaks
#res5 <- dplyr::filter(res4,Nrde2.peak=="yes",Contrast.x==contrasts[i])
#OR: select single contrast RNAseq data with Ccdc174 peaks
res5 <- dplyr::filter(res4,Ccdc174.peak=="yes",Contrast.x==contrasts[i])

#select the not upregulated genes (all genes with Fc< 0)
res6 <- dplyr::filter(res5,log2FoldChange.y < 0)
table(res6$intron_regulated)
res6$intron_regulated2 <- ifelse(res6$intron_regulated=="up","up","not_up")
i.reg <- c("up","not_up")

#plot the 2 seq logos

for (r in seq_along(i.reg)){
#get the sequecne  and make logo of upregulated 5'SS
  up.seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10,spliceDonors8n[names(spliceDonors8n) %in% res6$ID[res6$intron_regulated2==i.reg[r]]])
  up.matrix <- consensusMatrix(up.seqs,as.prob=FALSE) [1:4,]
  up.pfm <- PFMatrix(ID="Nrde2_up", name="Nrde2_up", 
                               matrixClass="Zipper-Type", strand="+",
                               bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                               tags=list(family="Helix-Loop-Helix", species="10090",
                                         tax_group="vertebrates",medline="7592839", 
                                         type="SELEX",ACC="P53762", pazar_tf_id="TF0000003",
                                         TFBSshape_ID="11", TFencyclopedia_ID="580"),
                               profileMatrix=up.matrix
  )
  Nrde2_up.icm <- toICM(up.pfm)
  png(sprintf("Figures/seqlogos_in_Ccdc174bound_introns_%s_in_%s.png",i.reg[r],contrasts[i]),height=400,width=400)
  seqLogo(Nrde2_up.icm)
  dev.off()
  
}
  
}


