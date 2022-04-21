rm(list=ls())

#install.packages("/tungstenfs/scratch/gbuehler/bioinfo/Rpackages/MiniChip_0.0.0.9000.tar.gz",repos=NULL)
library(tidyverse)
#library(Rsubread)
library(GenomicFeatures)
library(GenomicRanges)
library(MiniChip)
library(ComplexHeatmap)
library(wesanderson)
library(job)

#setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")
setwd("/tungstenfs/nobackup/gbuehler/michi/Projects/Nrde2/")

#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#   #load exon junctions
#----------------------------------------------------------------------------------------------------------------

load("splice_sites_expressed_and_notexpressed_transcripts_log2TPMabove0_and_genes_FPKMabove0.RData")

#spliceDonors <- unique(spliceDonors[seqnames(spliceDonors)=="chr1"])
names(spliceDonors) <- paste(names(spliceDonors),start(spliceDonors),sep="_")

#spliceDonors.ne <- unique(spliceDonors.ne[seqnames(spliceDonors.ne)=="chr1"])
names(spliceDonors.ne) <- paste(names(spliceDonors.ne),start(spliceDonors.ne),sep="_")

#spliceAcceptors <- unique(spliceAcceptors[seqnames(spliceAcceptors)=="chr1"])
names(spliceAcceptors) <- paste(names(spliceAcceptors),start(spliceAcceptors),sep="_")

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

#calculate heatmaps for splice donor sites (whole reads)
job::job(
{counts.donors <- SummitHeatmap(spliceDonors,bamFiles,bamNames,span,step,useCPM=TRUE,strand=1,read2pos=0,readShiftSize=0)}
)

#calculate heatmaps for splice donor sites (3' ends only)
job::job(
  {counts.donors <- SummitHeatmap(spliceDonors,bamFiles,bamNames,span,step,useCPM=TRUE,strand=1,read2pos=3,readShiftSize=0,nonSplitOnly=TRUE)}
)

#calculate heatmaps for splice donor sites (3' ends only), spliced reads only
job::job(
  {counts.donors <- SummitHeatmap(spliceDonors,bamFiles,bamNames,span,step,useCPM=TRUE,strand=1,read2pos=3,readShiftSize=0,splitOnly=TRUE)}
)

#----------------------------------------------------------------------------------------------------------------
#   #combine replicates for heatmaps
#----------------------------------------------------------------------------------------------------------------
#old data
sampleList <- list(
                   Nrde2WT=c("Nrde2WT_rep1","Nrde2WT_rep2"),
                   Nrde2d200= c("Nrde2d200_rep1","Nrde2d200_rep2"),
                   Nrde2D174R= c("Nrde2D174R_rep1","Nrde2D174R_rep2"),
                   Mtr4=c("Mtr4_rep1","Mtr4_rep2"),
                   Cdc174=c("Ccdc174_rep1","Ccdc174_rep2"),
                   Eif4a3=c("Eif4a3_rep1","Eif4a3_rep2")
                   )
counts.donors.means <- SummarizeHeatmaps(counts.donors,sampleList)

#save(counts.donors.means,file="FigureData/mean_cpms_over_2reps_in_all_active_gene_junctions_wholeReads_unique_210428.RData")
save(counts.donors.means,file="FigureData/mean_cpms_over_2reps_in_all_active_gene_junctions_wholeReads_unique_3primeEnds_210924.RData")

#new data
sampleList <- list(
  Nrde2dTSme=c("Nrde2_dTSmE_1","Nrde2_dTSmE_2"),
  Nrde2OE= c("Nrde2OE_1","Nrde2OE_2"),
  Ccdc174_Nrde2KO= c("Ccdc174_Nrde2KO_1","Ccdc174_Nrde2KO_2"),
  SmE=c("SmE_1","SmE_2")
  
)
counts.donors.means.new <- SummarizeHeatmaps(counts.donors,sampleList)
#save(counts.donors.means.new,file="FigureData/mean_cpms_over_2reps_in_all_active_gene_junctions_wholeReads_unique_new_BCLIP_210428.RData")
save(counts.donors.means.new,file="FigureData/mean_cpms_over_2reps_in_all_active_gene_junctions_wholeReads_unique_new_BCLIP_3primeEnds_210924.RData")

#load heatmap data
#print(load("FigureData/mean_cpms_over_2reps_in_all_active_gene_junctions_wholeReads_unique_210428.RData"))
#print(load("FigureData/mean_cpms_over_2reps_in_all_active_gene_junctions_wholeReads_unique_new_BCLIP_210428.RData"))
counts.donors.means <- c(counts.donors.means,counts.donors.means.new)

#----------------------------------------------------------------------------------------------------------------
#   #draw combined heatmaps without new data
#----------------------------------------------------------------------------------------------------------------

#plotcols1 <- wes_palette("Darjeeling1")
#plotcols2 <- wes_palette("Darjeeling2")
#plotcols <- c(plotcols2[c(2,4,3)],plotcols1[c(1,5,2)])

plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A")

#plot(1:6,1:6,col=plotcols,pch=20)


medianCpm <- c(rep(0.9,5),0.9)
topCpm <- rep(10,6)

heatmap_list <- DrawSummitHeatmaps(counts.donors.means, names(counts.donors.means), 
                                   plotcols = plotcols, use.log = TRUE,
                                   orderSample = 1, orderWindows = 200, medianCpm = medianCpm,
                                   summarizing = "mean", show_axis = FALSE,topCpm=topCpm)

pdf("Figures/heatmaps_at_5p-splice_junctions_orderedbyNrde2_60windows__whole_UniqueReads_reps_combined_adjusted_scale_210428.pdf",
    width=10,height=10)
draw(heatmap_list, padding = unit(c(3, 8, 8, 2), "mm"),show_heatmap_legend=TRUE)
for(i in 1:length(heatmap_list)){
  decorate_heatmap_body(names(heatmap_list)[i], {
    grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lty = 2, lwd = 1))
  })
}
dev.off()


#----------------------------------------------------------------------------------------------------------------
#   #draw combined heatmaps with new data
#----------------------------------------------------------------------------------------------------------------

#plotcols1 <- wes_palette("Darjeeling1")
#plotcols2 <- wes_palette("Darjeeling2")
#plotcols <- c(plotcols2[c(2,4,3)],plotcols1[c(1,5,2)])

plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A",rep("red",4))

#plot(1:6,1:6,col=plotcols,pch=20)


medianCpm <- c(rep(0.9,10))

heatmap_list <- DrawSummitHeatmaps(counts.donors.means, names(counts.donors.means), 
                                   plotcols = plotcols, use.log = TRUE,
                                   orderSample = 1, orderWindows = 200, medianCpm = medianCpm,
                                   summarizing = "mean", show_axis = FALSE)

pdf("Figures/heatmaps_at_5p-splice_junctions_orderedbyNrde2_60windows__whole_UniqueReads_reps_combined_including_new_BCLIPs.pdf",
    width=15,height=10)
draw(heatmap_list, padding = unit(c(3, 8, 8, 2), "mm"),show_heatmap_legend=TRUE)
for(i in 1:length(heatmap_list)){
  decorate_heatmap_body(names(heatmap_list)[i], {
    grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lty = 2, lwd = 1))
  })
}
dev.off()

#----------------------------------------------------------------------------------------------------------------
#   #draw only Ccdc174 heatmaps
#----------------------------------------------------------------------------------------------------------------

medianCpm <- c(rep(0.9,2))
topCpm <- rep(10,2)
plotcols <- c(  "#F98400","#F98400")

heatmap_list <- DrawSummitHeatmaps(counts.donors.means[c(5,9)], names(counts.donors.means[c(5,9)]), 
                                   plotcols = plotcols, use.log = TRUE,
                                   orderSample = 1, orderWindows = 200, medianCpm = medianCpm, topCpm=topCpm,
                                   summarizing = "mean", show_axis = FALSE)

pdf("Figures/heatmaps_at_5p-splice_junctions_Ccdc174_in_WT_and_Nrde2KO_210428.pdf",
    width=5,height=10)
draw(heatmap_list, padding = unit(c(3, 8, 8, 2), "mm"),show_heatmap_legend=TRUE)
for(i in 1:length(heatmap_list)){
  decorate_heatmap_body(names(heatmap_list)[i], {
    grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lty = 2, lwd = 1))
  })
}
dev.off()

#----------------------------------------------------------------------------------------------------------------
#   #draw only Nrde2 heatmaps
#----------------------------------------------------------------------------------------------------------------

medianCpm <- c(rep(0.9,2))
topCpm <- rep(10,2)
plotcols <- c(  "#046C9A","#046C9A")

heatmap_list <- DrawSummitHeatmaps(counts.donors.means[c(1,7)], names(counts.donors.means[c(1,7)]), 
                                   plotcols = plotcols, use.log = TRUE,
                                   orderSample = 1, orderWindows = 200, medianCpm = medianCpm, topCpm=topCpm,
                                   summarizing = "mean", show_axis = FALSE)

pdf("Figures/heatmaps_at_5p-splice_junctions_Nrde2_in_WT_and_SmEKO_210428.pdf",
    width=5,height=10)
draw(heatmap_list, padding = unit(c(3, 8, 8, 2), "mm"),show_heatmap_legend=TRUE)
for(i in 1:length(heatmap_list)){
  decorate_heatmap_body(names(heatmap_list)[i], {
    grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lty = 2, lwd = 1))
  })
}
dev.off()

#----------------------------------------------------------------------------------------------------------------
#   #make cumulative plots (60bp) Eif4a3, Nrde2WT, Ccdc174
#----------------------------------------------------------------------------------------------------------------
plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A")

cumu.counts <- CumulativePlots(
  counts.donors.means[c(1,5,6)],
  bamNames=names(counts.donors.means)[c(1,5,6)],
  span = span,
  step = step,
  summarizing = "mean",
  plot = FALSE,
  confInterval = 0.95
)
cumu.counts2 <- filter(cumu.counts,abs(position)<=100)

p <- ggplot(cumu.counts2,aes(x=position, y=mean_overlap2)) + geom_smooth(aes(ymin=ci.lower_overlap2,ymax=ci.upper_overlap2,fill=name,color=name),stat="identity") 
p <- p + theme_classic() + ylab("log2(cpm)") + 
  scale_color_manual(values=(plotcols[c(5,6,1)]),labels = names(counts.donors.means)[c(5,6,1)]) + 
  scale_fill_manual(values=(plotcols[c(5,6,1)]),labels = names(counts.donors.means)[c(5,6,1)])
p
#ggsave("Figures/metaplots_at_5p-splice_junctions_wholeUniqueReads_reps_combined_210428.pdf",device="pdf",height=4,width=5)
ggsave("Figures/metaplots_at_5p-splice_junctions_wholeUniqueReads_reps_combined_3primeEnds_210924_spliced_reads_only.pdf",device="pdf",height=4,width=5)


#----------------------------------------------------------------------------------------------------------------
#   #make cumulative plots split by gene expression (Nrde2 vs Mtr4) Eif4a3, Nrde2WT, Ccdc174, Mtr4
#----------------------------------------------------------------------------------------------------------------
res2 <- read.table("FigureData/Nrde2_Mtr4_Ccdc174_RNAseq_DEseq2_results_v2_cutoff200.txt",sep="\t",header=TRUE)
res.reg1 <- dplyr::filter(res2,padj < 0.01, abs(log2FoldChange)>1 & Contrast== "Nrde2KO_vs_WT")
res.reg2 <- dplyr::filter(res2,padj < 0.01, abs(log2FoldChange)>1 & Contrast== "Mtr4KO_vs_WT")
res.reg3 <- dplyr::filter(res2,padj < 0.01, abs(log2FoldChange)>1 & Contrast== "Nrde2D174R_vs_WT")
res.reg4 <- dplyr::filter(res2,padj < 0.01, abs(log2FoldChange)>1 & Contrast== "Nrde2d200_vs_WT")
res.reg5 <- dplyr::filter(res2,padj < 0.01, abs(log2FoldChange)>1 & Contrast== "Ccdc174KD_vs_WT")
regulated.genes <- unique(c(res.reg1$GeneID,res.reg2$GeneID,res.reg3$GeneID,res.reg4$GeneID,res.reg5$GeneID))

res3 <- pivot_wider(res2[,c("log2FoldChange","Contrast","GeneID")],names_from = Contrast,values_from = log2FoldChange)
res4 <- res3[res3$GeneID %in% regulated.genes,c(3,5,4,6,2)]
res5 <- res3[res3$GeneID %in% unique(c(res.reg1$GeneID)),c(3,5,4,6,2)]
res6 <- res3[res3$GeneID %in% unique(c(res.reg2$GeneID)),c(3,5,4,6,2)]

#group the genes into Mtr4 and Nrde2 regulated
Nrde2.up <- dplyr::filter(res2,padj < 0.01, log2FoldChange > 1 & Contrast== "Nrde2KO_vs_WT")
Nrde2.down <- dplyr::filter(res2,padj < 0.01, log2FoldChange < -1 & Contrast== "Nrde2KO_vs_WT")
Mtr4.up <- dplyr::filter(res2,padj < 0.01, log2FoldChange > 1 & Contrast== "Mtr4KO_vs_WT")
Mtr4.down <- dplyr::filter(res2,padj < 0.01, log2FoldChange < -1 & Contrast== "Mtr4KO_vs_WT")

NMreg <- ifelse(res3$GeneID %in% Nrde2.up$GeneID & res3$GeneID %in% Mtr4.up$GeneID, "Nrde2&Mtr4_up",
                ifelse(res3$GeneID %in% Nrde2.down$GeneID & res3$GeneID %in% Mtr4.down$GeneID, "Nrde2&Mtr4_down", 
                       ifelse(res3$GeneID %in% Nrde2.down$GeneID & res3$GeneID %in% Mtr4.down$GeneID==FALSE, "Nrde2_down",
                              ifelse(res3$GeneID %in% Nrde2.up$GeneID & res3$GeneID %in% Mtr4.up$GeneID==FALSE, "Nrde2_up",
                                     ifelse(res3$GeneID %in% Nrde2.down$GeneID==FALSE & res3$GeneID %in% Mtr4.down$GeneID, "Mtr4_down",
                                            ifelse(res3$GeneID %in% Nrde2.up$GeneID==FALSE & res3$GeneID %in% Mtr4.up$GeneID, "Mtr4_up","other"))))))
table(NMreg)
res3$NMreg <- NMreg
Nrde2_up.genes <- res3$GeneID[res3$NMreg=="Nrde2_up"]
Nrde2_down.genes <- res3$GeneID[res3$NMreg=="Nrde2_down"]
Nrde2Mtr4_up.genes <- res3$GeneID[res3$NMreg=="Nrde2&Mtr4_up"]
Nrde2Mtr4_down.genes <- res3$GeneID[res3$NMreg=="Nrde2&Mtr4_down"]
reg.genes <- list(Nrde2_up.genes,Nrde2_down.genes,Nrde2Mtr4_up.genes,Nrde2Mtr4_down.genes)
names(reg.genes) <- c("Nrde2_up","Nrde2_down","Nrde2Mtr4_up","Nrde2Mtr4_down")
#select the IDs from the GRanges object of 5'SS based on DE gene IDs
spliceDonors$GeneID <- matrix(unlist(strsplit(spliceDonors$GeneID,split = ".",fixed=TRUE)),ncol=2,byrow=TRUE)[,1]

cumu.counts <- list()
for (i in seq_along(reg.genes)){
  
cumu.counts[[i]] <- CumulativePlots(
  counts.donors.means[c(1,5,4,6)],
  bamNames=names(counts.donors.means)[c(1,5,4,6)],
  span = span,
  step = step,
  summarizing = "mean",
  plot = FALSE,
  confInterval = 0.95,overlapNames = names(spliceDonors)[spliceDonors$GeneID %in% reg.genes[[i]]]
)
cumu.counts[[i]]$DE <- names(reg.genes)[i]
}

cumu.counts2 <- do.call("rbind",cumu.counts)

p <- ggplot(cumu.counts2,aes(x=position, y=mean_overlap1)) + geom_smooth(aes(ymin=ci.lower_overlap1,ymax=ci.upper_overlap1,fill=DE,color=DE),stat="identity") +
  facet_wrap(vars(name))
p

p2 <- ggplot(cumu.counts2,aes(x=position, y=mean_overlap2)) + geom_smooth(aes(ymin=ci.lower_overlap2,ymax=ci.upper_overlap2,fill=DE,color=DE),stat="identity") +
  facet_wrap(vars(name))
p2

p <- p + theme_classic() + ylab("log2(cpm)") 
p
ggsave("Figures/metaplots_at_5p-splice_junctions_wholeUniqueReads_reps_combined_split_by_Nrde2_Mtr4_geneDE_210428.pdf",device="pdf",height=4,width=5)


#----------------------------------------------------------------------------------------------------------------
#   #make cumulative plots  Ccdc174 in Wt and Nrde2 KO
#----------------------------------------------------------------------------------------------------------------

cumu.counts <- CumulativePlots(
  counts.donors.means[c(5,9)],
  bamNames=names(counts.donors.means)[c(5,9)],
  span = span,
  step = step,
  summarizing = "mean",
  plot = FALSE,
  confInterval = 0.95
)
cumu.counts2 <- cumu.counts
plotcols <- c(  "#F98400","#F98400")

p <- ggplot(cumu.counts2,aes(x=position, y=mean_overlap2)) + geom_smooth(aes(ymin=ci.lower_overlap2,ymax=ci.upper_overlap2,fill=name,color=name,linetype=name),stat="identity") 
p <- p + theme_classic() + ylab("log2(cpm)") + 
  scale_color_manual(values=plotcols,labels = names(counts.donors.means)[c(5,9)]) + 
  scale_fill_manual(values=plotcols,labels = names(counts.donors.means)[c(5,9)]) +
  scale_linetype_manual(values=c("dashed", "solid"))
p
ggsave("Figures/metaplots_at_5p-splice_junctions_wholeUniqueReads_reps_combined_Ccdc174_210428.pdf",device="pdf",height=4,width=5)


#----------------------------------------------------------------------------------------------------------------
#   #make cumulative plots  Nrde2 in Wt and SmE KO
#----------------------------------------------------------------------------------------------------------------

cumu.counts <- CumulativePlots(
  counts.donors.means[c(1,7)],
  bamNames=names(counts.donors.means)[c(1,7)],
  span = span,
  step = step,
  summarizing = "mean",
  plot = FALSE,
  confInterval = 0.95
)
cumu.counts2 <- cumu.counts
plotcols <- c(  "#046C9A","#046C9A")

p <- ggplot(cumu.counts2,aes(x=position, y=mean_overlap2)) + geom_smooth(aes(ymin=ci.lower_overlap2,ymax=ci.upper_overlap2,fill=name,color=name,linetype=name),stat="identity") 
p <- p + theme_classic() + ylab("log2(cpm)") + 
  scale_color_manual(values=plotcols,labels = names(counts.donors.means)[c(1,7)]) + 
  scale_fill_manual(values=plotcols,labels = names(counts.donors.means)[c(1,7)]) +
  scale_linetype_manual(values=c("dashed", "solid"))
p
ggsave("Figures/metaplots_at_5p-splice_junctions_wholeUniqueReads_reps_combined_Nrde2_210428.pdf",device="pdf",height=4,width=5)


