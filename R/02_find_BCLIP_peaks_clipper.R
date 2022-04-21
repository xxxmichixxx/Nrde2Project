conda activate clipper3

cd /tungstenfs/scratch/gbuehler/flemrmat/data/2456/raw_data/
  
  for sample in Ccdc174_1 Ccdc174_2 Eif4a3_1 Eif4a3_2 Mtr4_1 Mtr4_2 Nrde2_1 Nrde2_2 Nrde2D174R_1 Nrde2D174R_2 Nrde2d200_1 Nrde2d200_2
do
clipper -b 2456_BCLIP_${sample}.LCstripped.prinseq_Aligned.sortedByCoord.out.bam -o ${sample}.clipper_peaks.bed -s mm10
done

#ENCODE version of calling clipper (doesnt find option bonferroni)
for sample in Ccdc174_1 Ccdc174_2 Eif4a3_1 Eif4a3_2 Mtr4_1 Mtr4_2 Nrde2_1 Nrde2_2 Nrde2D174R_1 Nrde2D174R_2 Nrde2d200_1 Nrde2d200_2
do
clipper -b 2456_BCLIP_${sample}.LCstripped.prinseq_Aligned.sortedByCoord.out.bam -s mm10 -o ${sample}.clipper_peaks_vE.bed --bonferroni -- superlocal --threshold-method binomial --save-pickle 
done


#-----------------------------------------------------------
# start R script####
#-----------------------------------------------------------
rm(list=ls())


library(tidyverse)
library(Rsubread)
library(GenomicFeatures)
library(GenomicRanges)
library(MiniChip)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(EnsDb.Mmusculus.v79)
library(ChIPseeker)
library(ComplexHeatmap)
library(UpSetR)
library(wesanderson)


setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")


#-----------------------------------------------------------
# combine peaks####
#-----------------------------------------------------------
#load all peaks
bedFiles <- c(list.files("/tungstenfs/scratch/gbuehler/flemrmat/data/2456/raw_data/", full.names=TRUE,pattern="*.clipper_peaks.bed$"))
bedNames <- gsub("/tungstenfs/scratch/gbuehler/flemrmat/data/2456/raw_data//","",bedFiles)
bedNames <- gsub(".clipper_peaks.bed","",bedNames)


#make a GRanges object for each sample
peaks.list <- GenomicRanges::GRangesList()
for (i in seq_along(bedFiles)){
  peaks.df <- read.table(bedFiles[i],sep="\t",header=FALSE)
  colnames(peaks.df) <- c("chr","start","end","name","pval","strand","summit.start","summit.end")
  peaks.list[[i]] <- makeGRangesFromDataFrame(peaks.df,
                                              keep.extra.columns=TRUE,
                                              seqinfo=NULL,
                                              seqnames.field=c("chr"),
                                              start.field=c("start"),
                                              end.field=c("end"),
                                              strand.field=c("strand"),
                                              ignore.strand=FALSE,
                                              starts.in.df.are.0based=FALSE)
  names(peaks.list[[i]]) <- peaks.list[[i]]$name
}

#combine replicates by overlaps
groups <- c("Ccdc174","Eif4a3","Mtr4","Nrde2","Nrde2D174R","Nrde2d200")
peaks.list2 <- peaks.list
for (g in 1:length(peaks.list)) {
  if (g %% 2 ==1) {
    peaks.list2[[g]] <- subsetByOverlaps(peaks.list[[g]],peaks.list[[g+1]])
    peaks.list2[[g]] <- peaks.list2[[g]][width(peaks.list2[[g]]) > 35]
  }
}
peaks.list3 <- peaks.list2[1:length(peaks.list) %% 2 ==1]
names(peaks.list3) <- groups

for (g in seq_along(peaks.list3)){
  plot(density(peaks.list3[[g]]$pval),main="pval")
  plot(density(width(peaks.list3[[g]])),main="width")
  
}

#remove peaks with < 5 reads in 50bp windnow around peak middle
bamFiles <- c(list.files("/tungstenfs/scratch/gbuehler/flemrmat/data/2456/raw_data/", full.names=TRUE,pattern="*LCstripped.prinseq_Aligned.sortedByCoord.out.bam$"))[1:12]
bamNames <- c("Ccdc174_rep1","Ccdc174_rep2","Eif4a3_rep1","Eif4a3_rep2","Mtr4_rep1","Mtr4_rep2","Nrde2WT_rep1","Nrde2WT_rep2",
              "Nrde2D174R_rep1","Nrde2D174R_rep2","Nrde2d200_rep1","Nrde2d200_rep2","Nrde2plusIGG_rep1","Nrde2plusIGG_rep2",
              "Nrde2plusThA_rep1","Nrde2plusThA_rep2")[1:12]
cpm <- list()
read.counts <- list(0)
peaks.list4 <- peaks.list3
for (g in seq_along(peaks.list3)){
  peaks <- peaks.list3[[g]]
  bamFilesG <- bamFiles[grep(paste(groups[g],"_",sep=""),bamFiles)]
  bamNamesG <- bamNames[grep(paste(groups[g],"_",sep=""),bamFiles)]
  
  counts <-  CountPeakReads(
    peaks,
    bamFilesG,
    bamNamesG,
    chips=bamNamesG,
    width = 50,
    minMQS = 0,
    strand = 1,
    read2pos = 0
  )
  cpm[[g]] <- counts[[1]]
  read.counts[[g]] <- counts[[2]]
  peaks.list4[[g]] <- peaks[rowSums(read.counts[[g]]) > 5]
}




#save bed files of combined peaks
for (i in seq_along(peaks.list4)){
  peaks.df2 <- data.frame(chrom=seqnames(peaks.list4[[i]]),
                          start=start(peaks.list4[[i]]),
                          stop=end(peaks.list4[[i]]),
                          gene=names(peaks.list4[[i]]),
                          score=-log10(peaks.list4[[i]]$pval),
                          strand=ifelse(strand(peaks.list4[[i]])=="-",-1,1), 
                          type=rep("exon",length(peaks.list4[[i]]))
  )
  write.table(peaks.df2,file=sprintf("clipper_peaks_combined_reps_5readsmin_%s.bed",names(peaks.list4)[i]),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
}

#add BCLIP name column to GRanges lists
for (i in seq_along(peaks.list4)){
  peaks.list4[[i]]$sample <- names(peaks.list4)[i]
}

clipper.peaks <- unlist(peaks.list4)
table(clipper.peaks$sample)
save(clipper.peaks,file="clipper_peaks_5readsmin_200202.RData")
