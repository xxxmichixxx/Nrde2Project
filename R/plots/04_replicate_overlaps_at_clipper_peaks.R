rm(list=ls())

library(tidyverse)
library(GenomicFeatures)
library(GenomicRanges)
library(MiniChip)
library(ComplexHeatmap)
library(wesanderson)
library(viridis)
setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")

#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#   #load peaks
#----------------------------------------------------------------------------------------------------------------
print(load("clipper_peaks_5readsmin_200202.RData"))

#set the summit start and ed to start and end
start(clipper.peaks) <- clipper.peaks$summit.start
end(clipper.peaks) <- clipper.peaks$summit.end

clipper.peaks$ID <- names(clipper.peaks)
#----------------------------------------------------------------------------------------------------------------
#   #load bam files
#----------------------------------------------------------------------------------------------------------------

#new BCLIP data
#new_files <- c(list.files("/tungstenfs/scratch/gbuehler/flemrmat/data/2456/raw_data/", full.names=TRUE,pattern="*LCstripped.prinseq_Aligned.sortedByCoord.out.bam$"))


bamFiles <- c(list.files("/tungstenfs/scratch/gbuehler/flemrmat/data/2456/raw_data/", full.names=TRUE,pattern="*LCstripped.prinseq_Aligned.sortedByCoord.out.bam$"))[1:12]
bamNames <- c("Ccdc174_rep1","Ccdc174_rep2","Eif4a3_rep1","Eif4a3_rep2","Mtr4_rep1","Mtr4_rep2","Nrde2WT_rep1","Nrde2WT_rep2",
              "Nrde2D174R_rep1","Nrde2D174R_rep2","Nrde2d200_rep1","Nrde2d200_rep2")

#----------------------------------------------------------------------------------------------------------------
#   #calculate cpms per peak
#----------------------------------------------------------------------------------------------------------------

counts <- CountPeakReads(
  clipper.peaks,
  bamFiles,
  bamNames = bamNames,
  chips = bamNames,
  width = 60,
  minOverlap = 1,
  PairedEnd = FALSE,
  minMQS = 255,
  strand = 1,
  read2pos = 0)

counts2 <- data.frame(counts[[1]])
counts2$ID <- row.names(counts2)
counts2 <- left_join(counts2,data.frame(mcols(clipper.peaks)),by="ID")

samples <- unique(counts2$sample)
reps1 <- seq(1,12,2)
reps2 <- seq(2,12,2)

pdf(file="Figures/logCPM_at_peaks_pairwise_between_replicates_v2.pdf",height=9,width=6)
par(mfrow=c(3,2))
for (s in seq_along(samples)){
counts3 <- dplyr::filter(counts2,sample==samples[s])
smoothScatter(log2(counts3[,reps1[s]]+1),log2(counts3[,reps2[s]]+1),xlim=c(0,14),ylim=c(0,14),
              xlab="replicate 1",ylab="replicate 2",main=sprintf("%d %s peaks",nrow(counts3),samples[s]))
}
dev.off()