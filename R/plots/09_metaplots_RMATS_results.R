rm(list=ls())

library(tidyverse)
library(Rsubread)
library(GenomicFeatures)
library(GenomicRanges)
library(MiniChip)

setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")

#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#   #load bam files
#----------------------------------------------------------------------------------------------------------------
bamFiles <- c(list.files("/tungstenfs/scratch/gbuehler/flemrmat/data/2456/raw_data/", full.names=TRUE,pattern="*LCstripped.prinseq_Aligned.sortedByCoord.out.bam$"))[1:12]
bamNames <- c("Ccdc174_rep1","Ccdc174_rep2","Eif4a3_rep1","Eif4a3_rep2","Mtr4_rep1","Mtr4_rep2","Nrde2WT_rep1","Nrde2WT_rep2",
              "Nrde2D174R_rep1","Nrde2D174R_rep2","Nrde2d200_rep1","Nrde2d200_rep2")

#find RNAseq bam files
all.bamFiles <- list.files("/tungstenfs/scratch/gbuehler/deepSeqRepos/bam/", full.names=TRUE,pattern="*bam$")
#bamFiles <- c(grep("2191F9",all.bamFiles,value=TRUE),grep("2191F10",all.bamFiles,value=TRUE),grep("2191F11",all.bamFiles,value=TRUE),grep("2191F12",all.bamFiles,value=TRUE))
bamFilesR <- c(grep("2191F",all.bamFiles,value=TRUE))
bamFilesR <- c(grep("spliced",bamFilesR,value=TRUE))
bamFilesR <- c(grep("Nxf1",bamFilesR,value=TRUE,invert=TRUE))
bamFilesR2 <- grep("2447F",all.bamFiles,value=TRUE)
bamFilesR <- c(bamFilesR,bamFilesR2)

bamNamesR <- gsub("/tungstenfs/scratch/gbuehler/deepSeqRepos/bam//","",bamFilesR)
bamNamesR <- gsub("_Aligned.sortedByCoord.out.bam","",bamNamesR)

#----------------------------------------------------------------------------------------------------------------
#   #load RMATS results of alternative 5' SS
#----------------------------------------------------------------------------------------------------------------
comp <- c("Nrde2KO_vs_WT","Nrde2d200_vs_WT","Nrde2D174R_vs_WT","Nrde2KO_vs_WT_CHX","Ccdc174KO_vs_WT","Ccdc174KO_vs_WT_CHX","Mtr4KO_vs_WT","WT_CHX_vs_WT")
i=8
#for (i in seq_along(comp)){
  
  A5SS.MATS.JCEC <- read.table(sprintf("rmats_analysis/%s/A5SS.MATS.JCEC.txt",comp[i]),header=TRUE)
  A5SS.gr <- makeGRangesFromDataFrame(A5SS.MATS.JCEC,
                                      keep.extra.columns=TRUE,
                                      seqinfo=NULL,
                                      seqnames.field=c("chr"),
                                      start.field=c("shortES"),
                                      end.field=c("shortEE"),
                                      strand.field=c("strand"),
                                      ignore.strand=FALSE,
                                      starts.in.df.are.0based=FALSE)
  A5SS.gr$upregulated <- ifelse(A5SS.gr$FDR < 0.05 & A5SS.gr$IncLevelDifference > 0, "yes","no")
  table(A5SS.gr$upregulated)
  
  A5SS.gr$regulated <- ifelse(A5SS.gr$FDR < 0.05, "yes","no")
  table(A5SS.gr$regulated)
  
  A5SS.gr2 <- resize(A5SS.gr,1L,fix="end")
  names(A5SS.gr2) <- paste(seqnames(A5SS.gr2),start(A5SS.gr2),A5SS.gr2$ID,sep="_")
  length(unique( names(A5SS.gr2))) == length( names(A5SS.gr2))
  #----------------------------------------------------------------------------------------------------------------
  #   #calculate heatmaps
  #----------------------------------------------------------------------------------------------------------------
  
  #define parameters for heatmaps
  span <- 200.5
  step <- 1
  
  #calculate heatmaps for splice donor sites
  counts.donors <- SummitHeatmap(A5SS.gr2,bamFiles,bamNames,span,step,useCPM=TRUE,strand=1,read2pos=0,readShiftSize=0)
  counts.RNA <- SummitHeatmap(A5SS.gr2,bamFilesR,bamNamesR,span,step,useCPM=TRUE,strand=2,read2pos=0,readShiftSize=0)
  
  #----------------------------------------------------------------------------------------------------------------
  #   #combine replicates for heatmaps
  #----------------------------------------------------------------------------------------------------------------
  sampleList <- list(
    Nrde2=c("Nrde2WT_rep1","Nrde2WT_rep2"),
    Cdc174=c("Ccdc174_rep1","Ccdc174_rep2"),
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
    Cdc174KO_CHX=c("Ccdc174KD_CHX_1_2447F15","Ccdc174KD_CHX_2_2447F16")
    
  )
  counts.RNA.means <- SummarizeHeatmaps(counts.RNA,sampleList)
  
  
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
    overlapNames = names(A5SS.gr2)[A5SS.gr2$upregulated=="yes"]
  )

  cumu.counts.RNA <- CumulativePlots(
    counts.RNA.means[1:4],
    bamNames=names(counts.RNA.means[1:4]),
    span = span,
    step = step,
    summarizing = "mean",
    plot = FALSE,
    confInterval = 0.95,
    overlapNames = names(A5SS.gr2)[A5SS.gr2$upregulated=="yes"]
  )
  
  cumu.counts2 <- cumu.counts %>% pivot_longer(contains("overlap"),names_to = c(".value","overlap"),names_sep="_")
  cumu.counts2$name <- factor(cumu.counts2$name,levels=c("Nrde2","Eif4a3"),labels=c("Nrde2","Eif4a3"))
  cumu.counts2$FivePrimeSS <- ifelse(cumu.counts2$overlap=="overlap1",sprintf("%d upregulated",table(A5SS.gr$upregulated)[2]),
                                     sprintf("%d not regulated",table(A5SS.gr$upregulated)[1]))
  
  cumu.counts.RNA2 <- cumu.counts.RNA %>% pivot_longer(contains("overlap"),names_to = c(".value","overlap"),names_sep="_")
  cumu.counts.RNA2$FivePrimeSS <- ifelse(cumu.counts.RNA2$overlap=="overlap1",sprintf("%d upregulated",table(A5SS.gr$upregulated)[2]),
                                         sprintf("%d not regulated",table(A5SS.gr$upregulated)[1]))
  
  p <- ggplot(cumu.counts2,aes(x=position, y=mean,col=name,linetype=FivePrimeSS)) + geom_smooth(aes(ymin=ci.lower,ymax=ci.upper,fill=name,color=name),stat="identity") +
    facet_wrap(vars(name))
  p <- p + theme_classic() + ylab("log2(cpm)") + 
    scale_color_manual(values=(plotcols[c(1,6)]),labels = names(counts.donors.means)[c(1,3)]) + 
    scale_fill_manual(values=(plotcols[c(1,6)]),labels = names(counts.donors.means)[c(1,3)]) 
  p
  
  p1 <- ggplot(cumu.counts.RNA2,aes(x=position, y=mean,col=name,linetype=FivePrimeSS)) + geom_smooth(aes(ymin=ci.lower,ymax=ci.upper,fill=name,color=name),stat="identity") +
    facet_wrap(vars(name))
  p1 <- p1 + theme_classic() + ylab("log2(cpm)") + 
    scale_color_manual(values=(plotcols[c(1,1,6,6)]),labels = names(counts.RNA.means[c(2,4,1,3)])) + 
    scale_fill_manual(values=(plotcols[c(1,1,6,6)]),labels = names(counts.RNA.means[c(2,4,1,3)]))
  p1
  p1+p
  
  
  ggsave("Figures/metaplots_at_CHX_vs_WT_RMATS_upregulated_FDR0.05_5p-splice_junctions_wholeUniqueReads_reps_combined_Nrde2_data.pdf",device="pdf",height=5,width=12)
  
  
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
    overlapNames = names(A5SS.gr2)[A5SS.gr2$upregulated=="yes"]
  )
  
  cumu.counts.RNA <- CumulativePlots(
    counts.RNA.means[5:8],
    bamNames=names(counts.RNA.means[5:8]),
    span = span,
    step = step,
    summarizing = "mean",
    plot = FALSE,
    confInterval = 0.95,
    overlapNames = names(A5SS.gr2)[A5SS.gr2$upregulated=="yes"]
  )
  
  cumu.counts2 <- cumu.counts %>% pivot_longer(contains("overlap"),names_to = c(".value","overlap"),names_sep="_")
 # cumu.counts2$name <- factor(cumu.counts2$name,levels=c("Ccdc174","Eif4a3"),labels=c("Ccdc174","Eif4a3"))
  cumu.counts2$FivePrimeSS <- ifelse(cumu.counts2$overlap=="overlap1",sprintf("%d upregulated",table(A5SS.gr$upregulated)[2]),
                                     sprintf("%d not regulated",table(A5SS.gr$upregulated)[1]))
  
  cumu.counts.RNA2 <- cumu.counts.RNA %>% pivot_longer(contains("overlap"),names_to = c(".value","overlap"),names_sep="_")
  cumu.counts.RNA2$FivePrimeSS <- ifelse(cumu.counts.RNA2$overlap=="overlap1",sprintf("%d upregulated",table(A5SS.gr$upregulated)[2]),
                                         sprintf("%d not regulated",table(A5SS.gr$upregulated)[1]))
  
  p <- ggplot(cumu.counts2,aes(x=position, y=mean,col=name,linetype=FivePrimeSS)) + geom_smooth(aes(ymin=ci.lower,ymax=ci.upper,fill=name,color=name),stat="identity") +
    facet_wrap(vars(name))
  p <- p + theme_classic() + ylab("log2(cpm)") + 
    scale_color_manual(values=(plotcols[c(5,6)]),labels = names(counts.donors.means)[c(2,3)]) + 
    scale_fill_manual(values=(plotcols[c(5,6)]),labels = names(counts.donors.means)[c(2,3)])
  p
  
  p1 <- ggplot(cumu.counts.RNA2,aes(x=position, y=mean,col=name,linetype=FivePrimeSS)) + geom_smooth(aes(ymin=ci.lower,ymax=ci.upper,fill=name,color=name),stat="identity") +
    facet_wrap(vars(name))
  p1 <- p1 + theme_classic() + ylab("log2(cpm)") + 
    scale_color_manual(values=(plotcols[c(5,5,6,6)]),labels = names(counts.RNA.means[c(6,8,5,7)])) + 
    scale_fill_manual(values=(plotcols[c(5,5,6,6)]),labels = names(counts.RNA.means[c(6,8,5,7)]))
  p1
  p1+p
  
  
  ggsave("Figures/metaplots_at_CHX_vs_WT_RMATS_upregulated_FDR0.05_5p-splice_junctions_wholeUniqueReads_reps_combined_Ccdc174data.pdf",device="pdf",height=5,width=12)
  
  
#}