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
#   #load bam files
#----------------------------------------------------------------------------------------------------------------
#new BCLIP data
new_files <- c(list.files("/tungstenfs/scratch/gbuehler/flemrmat/data/2743/raw_data/", full.names=TRUE,pattern="*LCstripped.prinseq_Aligned.sortedByCoord.out.bam$"))
new_names <- gsub("/tungstenfs/scratch/gbuehler/flemrmat/data/2743/raw_data//2743_BCLIP_","",new_files)
new_names <- gsub(".LCstripped.prinseq_Aligned.sortedByCoord.out.bam","",new_names)

bamFiles <- c(list.files("/tungstenfs/scratch/gbuehler/flemrmat/data/2456/raw_data/", full.names=TRUE,pattern="*LCstripped.prinseq_Aligned.sortedByCoord.out.bam$"))[1:12]
bamNames <- c("Ccdc174_rep1","Ccdc174_rep2","Eif4a3_rep1","Eif4a3_rep2","Mtr4_rep1","Mtr4_rep2","Nrde2WT_rep1","Nrde2WT_rep2",
              "Nrde2D174R_rep1","Nrde2D174R_rep2","Nrde2d200_rep1","Nrde2d200_rep2")
bamFiles <- c(bamFiles,new_files)
bamNames <- c(bamNames,new_names)

#find RNAseq bam files
all.bamFiles <- list.files("/tungstenfs/scratch/gbuehler/deepSeqRepos/bam/", full.names=TRUE,pattern="*bam$")
#bamFiles <- c(grep("2191F9",all.bamFiles,value=TRUE),grep("2191F10",all.bamFiles,value=TRUE),grep("2191F11",all.bamFiles,value=TRUE),grep("2191F12",all.bamFiles,value=TRUE))
bamFilesR <- c(grep("2191F",all.bamFiles,value=TRUE))
bamFilesR <- c(grep("spliced",bamFilesR,value=TRUE))
bamFilesR <- c(grep("Nxf1",bamFilesR,value=TRUE,invert=TRUE))
bamFilesR2 <- grep("2447F",all.bamFiles,value=TRUE)
bamFilesR3 <- c(grep("2701F",all.bamFiles,value=TRUE))[5:8]

bamFilesR <- c(bamFilesR,bamFilesR2,bamFilesR3)

bamNamesR <- gsub("/tungstenfs/scratch/gbuehler/deepSeqRepos/bam//","",bamFilesR)
bamNamesR <- gsub("_Aligned.sortedByCoord.out.bam","",bamNamesR)

#find RIBO bam files
bamFilesRi <- list.files("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/Riboseq/experiment_1/raw_data", full.names=TRUE,pattern="*LCstripped.prinseq_Aligned.sortedByCoord.out.bam$")
bamNamesRi <- gsub("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/Riboseq/experiment_1/raw_data/riboseq_","",bamFilesRi)
bamNamesRi <- gsub(".LCstripped.prinseq_Aligned.sortedByCoord.out.bam","",bamNamesRi)

#----------------------------------------------------------------------------------------------------------------
#   #load upregulated introns results 
#----------------------------------------------------------------------------------------------------------------

#load intronic read differential expression
#normal
res2 <- read.table(file="FigureData/Nrde2_Mtr4_Ccdc174_batch1and2_intronicRNAseq_DEseq2_results_0.05_0.5FC_5readspersample_seperateWTs.txt",sep="\t",header=TRUE)
#plus Cyclohexamide (CHX)
#res2 <- read.table(file="FigureData/CHX_Nrde2_Ccdc174_intronicRNAseq_DEseq2_results_0.05_0.5FC_5readspersample_seperateWTs.txt",sep="\t",header=TRUE)

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
#res1 <- read.table("FigureData/CHX_Nrde2_Ccdc174_RNAseq_DEseq2_results_cutoff80.txt",sep="\t",header=TRUE)

res2$GeneID <- matrix(unlist(strsplit(as.character(res2$GeneID),".",fixed=TRUE)),ncol=2,byrow = TRUE)[,1]
res2$GeneConID <- paste(res2$GeneID,res2$Contrast,sep="_")
res1$GeneConID <- paste(res1$GeneID,res1$Contrast,sep="_")

res3 <- left_join(res2,res1,by="GeneConID")
#----------------------------------------------------------------------------------------------------------------
#   # select upregulated introns in downregulated genes and plot BCLIP levels in them
#-----------------------------------------------------------------------------------------------------------------
contrasts <- unique(res3$Contrast.x)
#res3 <- res3[is.na(res3$intron_regulated)==FALSE,]

res3$intron_gene_regulated <- ifelse(res3$intron_regulated =="up" & res3$log2FoldChange.y < 0,"Intron up, gene down",
                                     ifelse(res3$log2FoldChange.x < 0 & res3$log2FoldChange.y < 0,"Intron down, gene down","other"))
table(res3$intron_gene_regulated)
#which genes have increased splcied reads in introns
unique(res3$gene_symbol[res3$Contrast.x=="Nrde2KO_vs_WT" & res3$intron_gene_regulated =="Intron up, gene down"])

#turn it into a GRanges object,selecting one contrast
i=6
res3.gr <- makeGRangesFromDataFrame(res3[res3$Contrast.x==contrasts[i],],
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames"),
                         start.field=c("start"),
                         end.field=c("end"),
                         strand.field=c("strand"),
                         starts.in.df.are.0based=FALSE)
names(res3.gr) <- res3.gr$ID
table(res3.gr$intron_gene_regulated)
res3.gr$intron.up <- ifelse(res3.gr$intron_gene_regulated=="Intron up, gene down","yes","no")
table(res3.gr$intron.up)
length(unique( names(res3.gr))) == length( names(res3.gr))

#remove all genes with FC >=0
inx0 <- which(res3.gr$log2FoldChange.y < 0)
res3.gr <- res3.gr[inx0]
table(res3.gr$intron_gene_regulated)

  #----------------------------------------------------------------------------------------------------------------
  #   #calculate heatmaps
  #----------------------------------------------------------------------------------------------------------------
  
  #define parameters for heatmaps
  span <- 200.5
  step <- 1
  
  #calculate heatmaps for splice donor sites
  counts.donors <- SummitHeatmap(res3.gr,bamFiles,bamNames,span,step,useCPM=TRUE,strand=1,read2pos=0,readShiftSize=0)
  counts.RNA <- SummitHeatmap(res3.gr,bamFilesR,bamNamesR,span,step,useCPM=TRUE,strand=2,read2pos=0,readShiftSize=0)
  counts.Ribo <- SummitHeatmap(res3.gr,bamFilesRi,bamNamesRi,span,step,useCPM=TRUE,strand=1,read2pos=0,readShiftSize=0)
  
  
  #----------------------------------------------------------------------------------------------------------------
  #   #combine replicates for heatmaps
  #----------------------------------------------------------------------------------------------------------------
  sampleList <- list(
    Nrde2=c("Nrde2WT_rep1","Nrde2WT_rep2"),
    Ccdc174=c("Ccdc174_rep1","Ccdc174_rep2"),
    Eif4a3=c("Eif4a3_rep1","Eif4a3_rep2")
  )
  counts.donors.means <- SummarizeHeatmaps(counts.donors,sampleList)
  
  sampleList <- list(
    Nrde2WT=c("Nrde2-WT_1_spliced_2191F9","Nrde2-WT_2_spliced_2191F10","Nrde2WT_1_2447F1","Nrde2WT_2_2447F2"),
    Nrde2KO=c("Nrde2-KO_1_spliced_2191F11","Nrde2-KO_2_spliced_2191F12","Nrde2KO_1_2447F5","Nrde2KO_2_2447F6" ),
    Nrde2WT_CHX=c("Nrde2WT_CHX_1_2447F3","Nrde2WT_CHX_2_2447F4"),
    Nrde2KO_CHX=c("Nrde2KO_CHX_1_2447F7","Nrde2KO_CHX_2_2447F8"),
    
    Cdc174WT=c("Ccdc174WT_1_2447F9","Ccdc174WT_2_2447F10"),
    Cdc174KO=c("Ccdc174KD_1_2447F13","Ccdc174KD_2_2447F14"),
    Cdc174WT_CHX=c("Ccdc174WT_CHX_1_2447F11", "Ccdc174WT_CHX_2_2447F12"),
    Cdc174KO_CHX=c("Ccdc174KD_CHX_1_2447F15","Ccdc174KD_CHX_2_2447F16"),
    Mtr4WT=c("Mtr4_WT_1_spliced_2191F19","Mtr4_WT_2_spliced_2191F20"),
    Mtr4KO=c("Mtr4_KO_1_spliced_2191F21","Mtr4_KO_2_spliced_2191F22"),
    SmEWT=c("SmE_untreated_rep1_2701F1","SmE_untreated_rep2_2701F2"),
    SmEKO=  c("SmE_dTAG_rep1_2701F3","SmE_dTAG_rep2_2701F4")
  )
  counts.RNA.means <- SummarizeHeatmaps(counts.RNA,sampleList)
  
  sampleList <- list(
    Nrde2WTribo=c("Nrde2WT_rep1","Nrde2WT_rep2"),
    Nrde2KOribo=c("Nrde2KO_rep1","Nrde2KO_rep2")
  )
  counts.ribo.means <- SummarizeHeatmaps(counts.Ribo,sampleList)
  
  save(counts.donors.means,counts.RNA.means,counts.ribo.means,file="FigureData/heatmap_counts_at_5pSS_Ccdc174KOvsWT_wholeUniqueReads_reps_combined_210429_only_genesdown.RData")
  
  load("FigureData/heatmap_counts_at_5pSS_Nrde2KOvsWT_wholeUniqueReads_reps_combined_210429_only_genesdown.RData")
  load("FigureData/heatmap_counts_at_5pSS_Ccdc174KOvsWT_wholeUniqueReads_reps_combined_210429_only_genesdown.RData")
  
  #----------------------------------------------------------------------------------------------------------------
  #   #make cumulative plots Eif4a3, Nrde2 
  #----------------------------------------------------------------------------------------------------------------
  plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A")
  
  cumu.counts <- CumulativePlots(
    counts.donors.means[c(1,3)],
    bamNames=names(counts.donors.means[c(1,3)]),
    span = span,
    step = step,
    summarizing = "mean",
    plot = FALSE,
    confInterval = 0.95,
    overlapNames = names(res3.gr)[res3.gr$intron.up=="yes"]
  )

  #option 1: including CHX
  cumu.counts.RNA <- CumulativePlots(
    counts.RNA.means[1:4],
    bamNames=names(counts.RNA.means[1:4]),
    span = span,
    step = step,
    summarizing = "mean",
    plot = FALSE,
    confInterval = 0.95,
    overlapNames = names(res3.gr)[res3.gr$intron.up=="yes"]
  )
  
  #option 2: incuding Riboseq
  cumu.counts.RNA <- CumulativePlots(
    c(counts.RNA.means[1:2],counts.ribo.means),
    bamNames=c(names(counts.RNA.means[1:2]),names(counts.ribo.means)),
    span = span,
    step = step,
    summarizing = "mean",
    plot = FALSE,
    confInterval = 0.95,
    overlapNames = names(res3.gr)[res3.gr$intron.up=="yes"]
  )
  
  

  cumu.counts2 <- cumu.counts %>% pivot_longer(contains("overlap"),names_to = c(".value","overlap"),names_sep="_")
  cumu.counts2$name <- factor(cumu.counts2$name,levels=c("Nrde2","Eif4a3"),labels=c("Nrde2","Eif4a3"))
  cumu.counts2$introns <- factor(ifelse(cumu.counts2$overlap=="overlap1",sprintf("%d upregulated",table(res3.gr$intron.up)[2]),
                                     sprintf("%d not upregulated",table(res3.gr$intron.up)[1])),
                                 levels=c(sprintf("%d upregulated",table(res3.gr$intron.up)[2]),sprintf("%d not upregulated",table(res3.gr$intron.up)[1])),
  labels=c(sprintf("%d upregulated",table(res3.gr$intron.up)[2]),sprintf("%d not upregulated",table(res3.gr$intron.up)[1]))
)

  
  cumu.counts.RNA2 <- cumu.counts.RNA %>% pivot_longer(contains("overlap"),names_to = c(".value","overlap"),names_sep="_")
  cumu.counts.RNA2$introns <- factor(ifelse(cumu.counts.RNA2$overlap=="overlap1",sprintf("%d upregulated",table(res3.gr$intron.up)[2]),
                                         sprintf("%d not upregulated",table(res3.gr$intron.up)[1])),levels=c(sprintf("%d upregulated",table(res3.gr$intron.up)[2]),sprintf("%d not upregulated",table(res3.gr$intron.up)[1])),
                                    labels=c(sprintf("%d upregulated",table(res3.gr$intron.up)[2]),sprintf("%d not upregulated",table(res3.gr$intron.up)[1]))
  )
  
  p <- ggplot(cumu.counts2,aes(x=position, y=mean,col=name,linetype=introns)) + geom_smooth(aes(ymin=ci.lower,ymax=ci.upper,fill=name,color=name),stat="identity") +
    facet_wrap(vars(name))
  p <- p + theme_classic() + ylab("log2(cpm)") + ylim(c(0,2.5)) +
    scale_color_manual(values=(plotcols[c(1,6)]),labels = names(counts.donors.means)[c(1,3)]) + 
    scale_fill_manual(values=(plotcols[c(1,6)]),labels = names(counts.donors.means)[c(1,3)]) 
  p
  
  p1 <- ggplot(cumu.counts.RNA2,aes(x=position, y=mean,col=name,linetype=introns)) + geom_smooth(aes(ymin=ci.lower,ymax=ci.upper,fill=name,color=name),stat="identity") +
    facet_wrap(vars(name))
  p1 <- p1 + theme_classic() + ylab("log2(cpm)") + 
  #  scale_color_manual(values=(plotcols[c(1,1,6,6)]),labels = names(counts.RNA.means[c(2,4,1,3)])) + 
   # scale_fill_manual(values=(plotcols[c(1,1,6,6)]),labels = names(counts.RNA.means[c(2,4,1,3)]))
    scale_color_manual(values=(plotcols[c(1,1,6,6)]),labels = names(counts.RNA.means[c(2,2,1,1)])) + 
    scale_fill_manual(values=(plotcols[c(1,1,6,6)]),labels = names(counts.RNA.means[c(2,2,1,1)]))
  p1
  p1+p
  
  
  ggsave("Figures/metaplots_at_Nrde2KO_vs_WT_CHX_introns_upregulated_5p-splice_junctions_wholeUniqueReads_reps_combined_Nrde2_data_210430_genesdownonly_ylim2.5.pdf",device="pdf",height=5,width=12)
 # ggsave("Figures/metaplots_at_Nrde2KO_vs_WT_RIBO_introns_upregulated_5p-splice_junctions_wholeUniqueReads_reps_combined_Nrde2_data_210430_genesdownonly.pdf",device="pdf",height=5,width=12)
  ggsave("Figures/metaplots_at_Ccdc174KO_vs_WT_CHX_introns_upregulated_5p-splice_junctions_wholeUniqueReads_reps_combined_Nrde2_data_210430_genesdownonly_ylim2.5.pdf",device="pdf",height=5,width=12)
  
  
  #----------------------------------------------------------------------------------------------------------------
  #   #make cumulative plots Eif4a3, Ccdc174 
  #----------------------------------------------------------------------------------------------------------------
  plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A")
  
  cumu.counts <- CumulativePlots(
    counts.donors.means[c(2,3)],
    bamNames=names(counts.donors.means[c(2,3)]),
    span = span,
    step = step,
    summarizing = "mean",
    plot = FALSE,
    confInterval = 0.95,
    overlapNames = names(res3.gr)[res3.gr$intron.up=="yes"]
  )
  
  cumu.counts.RNA <- CumulativePlots(
    counts.RNA.means[5:8],
    bamNames=names(counts.RNA.means[5:8]),
    span = span,
    step = step,
    summarizing = "mean",
    plot = FALSE,
    confInterval = 0.95,
    overlapNames = names(res3.gr)[res3.gr$intron.up=="yes"]
  )
  
  cumu.counts2 <- cumu.counts %>% pivot_longer(contains("overlap"),names_to = c(".value","overlap"),names_sep="_")
  cumu.counts2$name <- factor(cumu.counts2$name,levels=c("Ccdc174","Eif4a3"),labels=c("Ccdc174","Eif4a3"))
  cumu.counts2$introns <- factor(ifelse(cumu.counts2$overlap=="overlap1",sprintf("%d upregulated",table(res3.gr$intron.up)[2]),
                                        sprintf("%d not upregulated",table(res3.gr$intron.up)[1])),
                                 levels=c(sprintf("%d upregulated",table(res3.gr$intron.up)[2]),sprintf("%d not upregulated",table(res3.gr$intron.up)[1])),
                                 labels=c(sprintf("%d upregulated",table(res3.gr$intron.up)[2]),sprintf("%d not upregulated",table(res3.gr$intron.up)[1]))
  )
  
  
  cumu.counts.RNA2 <- cumu.counts.RNA %>% pivot_longer(contains("overlap"),names_to = c(".value","overlap"),names_sep="_")
  cumu.counts.RNA2$introns <- factor(ifelse(cumu.counts.RNA2$overlap=="overlap1",sprintf("%d upregulated",table(res3.gr$intron.up)[2]),
                                            sprintf("%d not upregulated",table(res3.gr$intron.up)[1])),levels=c(sprintf("%d upregulated",table(res3.gr$intron.up)[2]),sprintf("%d not upregulated",table(res3.gr$intron.up)[1])),
                                     labels=c(sprintf("%d upregulated",table(res3.gr$intron.up)[2]),sprintf("%d not upregulated",table(res3.gr$intron.up)[1]))
  )
  
  p <- ggplot(cumu.counts2,aes(x=position, y=mean,col=name,linetype=introns)) + geom_smooth(aes(ymin=ci.lower,ymax=ci.upper,fill=name,color=name),stat="identity") +
    facet_wrap(vars(name))
  p <- p + theme_classic() + ylab("log2(cpm)") + ylim(c(0,2.5)) +
    scale_color_manual(values=(plotcols[c(5,6)]),labels = names(counts.donors.means)[c(2,3)]) + 
    scale_fill_manual(values=(plotcols[c(5,6)]),labels = names(counts.donors.means)[c(2,3)])
  p
  
  p1 <- ggplot(cumu.counts.RNA2,aes(x=position, y=mean,col=name,linetype=introns)) + geom_smooth(aes(ymin=ci.lower,ymax=ci.upper,fill=name,color=name),stat="identity") +
    facet_wrap(vars(name))
  p1 <- p1 + theme_classic() + ylab("log2(cpm)") + 
    scale_color_manual(values=(plotcols[c(5,5,6,6)]),labels = names(counts.RNA.means[c(6,8,5,7)])) + 
    scale_fill_manual(values=(plotcols[c(5,5,6,6)]),labels = names(counts.RNA.means[c(6,8,5,7)]))
  p1
  p1+p
  
  
  ggsave("Figures/metaplots_at_Ccdc174KO_vs_WT_introns_upregulated_5p-splice_junctions_wholeUniqueReads_reps_combined_Ccdc174data_210430_genesdownonly_ylim2.5.pdf",device="pdf",height=5,width=12)
  
  
  #----------------------------------------------------------------------------------------------------------------
  #   #make cumulative plots RNA Nrde2,Ccdc174,Mtr4,SmE
  #----------------------------------------------------------------------------------------------------------------
  plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A")
  
  
  
  cumu.counts.RNA <- CumulativePlots(
    counts.RNA.means,
    bamNames=names(counts.RNA.means),
    span = span,
    step = step,
    summarizing = "mean",
    plot = FALSE,
    confInterval = 0.95,
    overlapNames = names(res3.gr)[res3.gr$intron.up=="yes"]
  )
  
  
  
cumu.counts.RNA2 <- cumu.counts.RNA %>% pivot_longer(contains("overlap"),names_to = c(".value","overlap"),names_sep="_")
cumu.counts.RNA2$introns <- factor(ifelse(cumu.counts.RNA2$overlap=="overlap1",sprintf("%d upregulated",table(res3.gr$intron.up)[2]),
                                            sprintf("%d not upregulated",table(res3.gr$intron.up)[1])),levels=c(sprintf("%d upregulated",table(res3.gr$intron.up)[2]),sprintf("%d not upregulated",table(res3.gr$intron.up)[1])),
                                     labels=c(sprintf("%d upregulated",table(res3.gr$intron.up)[2]),sprintf("%d not upregulated",table(res3.gr$intron.up)[1]))
  )
  

  
  p1 <- ggplot(cumu.counts.RNA2,aes(x=position, y=mean,col=name,linetype=introns)) + geom_smooth(aes(ymin=ci.lower,ymax=ci.upper,fill=name,color=name),stat="identity") +
    facet_wrap(vars(name))
  p1 <- p1 + ylab("log2(cpm)") 
  p1
  
  ggsave("Figures/metaplots_at_Nrde2KO_vs_WT_introns_upregulated_5p-splice_junctions_wholeUniqueReads_reps_combined_all_RNAseq_data.pdf",device="pdf",height=6,width=12)
  