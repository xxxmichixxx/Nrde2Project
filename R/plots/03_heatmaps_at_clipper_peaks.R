rm(list=ls())

library(tidyverse)
library(GenomicFeatures)
library(GenomicRanges)
library(MiniChip)
library(ComplexHeatmap)
library(wesanderson)

setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")

#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#   #load peaks
#----------------------------------------------------------------------------------------------------------------
print(load("clipper_peaks_5readsmin_200202.RData"))

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
#----------------------------------------------------------------------------------------------------------------
#   #calculate heatmaps
#----------------------------------------------------------------------------------------------------------------

#define parameters for heatmaps
span <- 200.5
step <- 1

#set the summit start and ed to start and end
start(clipper.peaks) <- clipper.peaks$summit.start
end(clipper.peaks) <- clipper.peaks$summit.end

#calculate heatmaps for splice donor sites
counts <- SummitHeatmap(clipper.peaks,bamFiles,bamNames,span,step,useCPM=TRUE,strand=1,read2pos=0,readShiftSize=0)

#re-order counts list
counts2 <- counts[c(7,8,11,12,9,10,15:18,5,6,1,2,13,14,3,4,19,20)]
#----------------------------------------------------------------------------------------------------------------
#   #draw  heatmaps
#----------------------------------------------------------------------------------------------------------------


plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A")


plotcols <- c(plotcols[1],plotcols[1],plotcols[2],plotcols[2],plotcols[3],plotcols[3],plotcols[4],plotcols[4],plotcols[5],plotcols[5],plotcols[6],plotcols[6])

medianCpm <- rep(1,length(names(counts2)))

heatmap_list <- DrawSummitHeatmaps(counts2, names(counts2), 
                                   plotcols = plotcols, use.log = TRUE,
                                   orderSample = 1, medianCpm = medianCpm,
                                   summarizing = "mean", show_axis = FALSE,
                                   splitHM = clipper.peaks$sample)

pdf("Figures/heatmaps_at_all_clipper_peaks_orderedbyNrde2_wholereads_unique.pdf",
    width=30,height=30)
draw(heatmap_list, padding = unit(c(3, 8, 8, 2), "mm"),show_heatmap_legend=TRUE)
#for(i in 1:length(heatmap_list)){
#  decorate_heatmap_body(names(heatmap_list)[i], {
#    grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lty = 2, lwd = 1))
#  })
#}
dev.off()


#----------------------------------------------------------------------------------------------------------------
#   #draw  heatmaps without Eif4a3 peaks
#----------------------------------------------------------------------------------------------------------------
counts3 <- counts2
for(i in seq_along(counts2)){
  counts3[[i]] <- counts2[[i]][clipper.peaks$sample != "Eif4a3",]
}



medianCpm <- rep(10,length(names(counts3)))

heatmap_list <- DrawSummitHeatmaps(counts3, names(counts3), 
                                   plotcols = plotcols, use.log = TRUE,
                                   orderSample = 1, medianCpm = medianCpm,
                                   summarizing = "mean", show_axis = FALSE,
                                   splitHM = clipper.peaks$sample[clipper.peaks$sample != "Eif4a3"])

pdf("Figures/heatmaps_at_all_clipper_peaks_except_Eif4a3_orderedbyNrde2_wholereads_unique.pdf",
    width=30,height=30)
draw(heatmap_list, padding = unit(c(3, 8, 8, 2), "mm"),show_heatmap_legend=TRUE)
#for(i in 1:length(heatmap_list)){
#  decorate_heatmap_body(names(heatmap_list)[i], {
#    grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lty = 2, lwd = 1))
#  })
#}
dev.off()

#----------------------------------------------------------------------------------------------------------------
#   #draw  heatmaps for each peak set
#----------------------------------------------------------------------------------------------------------------
peak.names <- unique(clipper.peaks$sample)[c(4,6,5,3,1,2)]

medianCpm <- rep(2,length(names(counts2)))

plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A")
plotcols <- c(plotcols[1],plotcols[1],plotcols[2],plotcols[2],plotcols[3],plotcols[3],
              "red","red","red","red",
              plotcols[4],plotcols[4],plotcols[5],plotcols[5],
              "red","red",
              plotcols[6],plotcols[6],
              "red","red")


#ordersamples <- seq(1,12,by = 2)
ordersamples <- c(1,2,3,11,13,17)
  
for (p in seq_along(peak.names)){
counts3 <- counts2
for(i in seq_along(counts2)){
  counts3[[i]] <- counts2[[i]][clipper.peaks$sample == peak.names[p],]
}


heatmap_list <- DrawSummitHeatmaps(counts3, names(counts3), 
                                   plotcols = plotcols, use.log = TRUE,TargetHeight = 300,
                                   orderSample = ordersamples[p], medianCpm = medianCpm,
                                   summarizing = "mean", show_axis = FALSE)

pdf(sprintf("Figures/heatmaps_%s_all_clipper_peaks_ordered_wholereads_unique_TH300_with_new_Bclip.pdf",peak.names[p]),
    width=20,height=10)
draw(heatmap_list, padding = unit(c(3, 8, 8, 2), "mm"),show_heatmap_legend=TRUE)
#for(i in 1:length(heatmap_list)){
#  decorate_heatmap_body(names(heatmap_list)[i], {
#    grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lty = 2, lwd = 1))
#  })
#}
dev.off()
}

#----------------------------------------------------------------------------------------------------------------
#   #draw  heatmaps for Ccdc174 peaks in Nrd2 WT and KO
#----------------------------------------------------------------------------------------------------------------
p=5
counts2C <- counts2[13:16]
medianCpm <- rep(2,length(names(counts2C)))

plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A")
plotcols <- c(plotcols[5],plotcols[5],
              "red","red")

  counts3 <- counts2C
  for(i in seq_along(counts2C)){
    counts3[[i]] <- counts2C[[i]][clipper.peaks$sample == peak.names[p],]
  }
  
  
  heatmap_list <- DrawSummitHeatmaps(counts3, names(counts3), 
                                     plotcols = plotcols, use.log = TRUE,TargetHeight = 500,
                                     orderSample = 1, medianCpm = medianCpm,
                                     summarizing = "mean", show_axis = FALSE)
  
  pdf("Figures/heatmaps_Ccdc174_clipper_peaks_ordered_wholereads_unique_TH500_Ccdc174_in_Nrde2_WT_and_KO.pdf",
      width=10,height=10)
  draw(heatmap_list, padding = unit(c(3, 8, 8, 2), "mm"),show_heatmap_legend=TRUE)
  dev.off()


