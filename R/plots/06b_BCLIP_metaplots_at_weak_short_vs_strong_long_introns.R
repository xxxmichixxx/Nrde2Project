
rm(list=ls())


#library(csaw)
library(tidyverse)
library(Rsubread)
library(GenomicFeatures)
library(GenomicRanges)
library(MiniChip)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(EnsDb.Mmusculus.v79)
#library(ChIPseeker)
#library(ComplexHeatmap)
#library(UpSetR)
library(patchwork)
library(ggpubr)

setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")



#load data
matt.res  <- read.table("FigureData/MATT_scores_on_all_introns_with_peak_overlap_info.txt",header=TRUE,sep="\t")

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

#find RIBO bam files
bamFilesRi <- list.files("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/Riboseq/experiment_1/raw_data", full.names=TRUE,pattern="*LCstripped.prinseq_Aligned.sortedByCoord.out.bam$")
bamNamesRi <- gsub("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/Riboseq/experiment_1/raw_data/riboseq_","",bamFilesRi)
bamNamesRi <- gsub(".LCstripped.prinseq_Aligned.sortedByCoord.out.bam","",bamNamesRi)


#----------------------------------------------------------------------------------------------------------------
#   #define intron features####
#----------------------------------------------------------------------------------------------------------------

matt.results2 <- dplyr::filter(matt.res,!is.na(BPSCORE_MAXBP),!is.na(NUM_PREDICTED_BPS),!is.na(SF1_HIGHESTSCORE_3SS),!is.na(MAXENTSCR_HSAMODEL_5SS),!is.na(MAXENTSCR_HSAMODEL_3SS))

MAXENTSCR_HSAMODEL_5SS.q <- quantile(matt.results2$MAXENTSCR_HSAMODEL_5SS)
MAXENTSCR_HSAMODEL_3SS.q <- quantile(matt.results2$MAXENTSCR_HSAMODEL_3SS)
BPSCORE_MAXBP.q <- quantile(matt.results2$BPSCORE_MAXBP)
NUM_PREDICTED_BPS.q <- quantile(matt.results2$NUM_PREDICTED_BPS)
SF1_HIGHESTSCORE_3SS.q <- quantile(matt.results2$SF1_HIGHESTSCORE_3SS)

#version1
matt.results2$splicing <- ifelse(matt.results2$MAXENTSCR_HSAMODEL_5SS < MAXENTSCR_HSAMODEL_5SS.q[3] &
                                   matt.results2$BPSCORE_MAXBP < BPSCORE_MAXBP.q[3] &
                                   matt.results2$NUM_PREDICTED_BPS < NUM_PREDICTED_BPS.q[3] &
                                   matt.results2$SF1_HIGHESTSCORE_3SS < SF1_HIGHESTSCORE_3SS.q[3],"weak",
                                 
                                 ifelse(matt.results2$MAXENTSCR_HSAMODEL_5SS > MAXENTSCR_HSAMODEL_5SS.q[3] &
                                          matt.results2$BPSCORE_MAXBP > BPSCORE_MAXBP.q[3] &
                                          matt.results2$NUM_PREDICTED_BPS > NUM_PREDICTED_BPS.q[3] &
                                          matt.results2$SF1_HIGHESTSCORE_3SS > SF1_HIGHESTSCORE_3SS.q[3],"strong",
                                        "mixed"))

#version2
plot(density(matt.results2$MAXENTSCR_HSAMODEL_5SS))
matt.results2$splicing <- ifelse(matt.results2$MAXENTSCR_HSAMODEL_5SS < MAXENTSCR_HSAMODEL_5SS.q[2] |
                                   matt.results2$MAXENTSCR_HSAMODEL_3SS < MAXENTSCR_HSAMODEL_3SS.q[2] |
                                   matt.results2$BPSCORE_MAXBP < BPSCORE_MAXBP.q[2] |
                                   matt.results2$NUM_PREDICTED_BPS < NUM_PREDICTED_BPS.q[2] |
                                   matt.results2$SF1_HIGHESTSCORE_3SS < SF1_HIGHESTSCORE_3SS.q[2],"weak",
                                 
                                 ifelse(matt.results2$MAXENTSCR_HSAMODEL_5SS > MAXENTSCR_HSAMODEL_5SS.q[4] |
                                          matt.results2$MAXENTSCR_HSAMODEL_3SS > MAXENTSCR_HSAMODEL_3SS.q[4] |
                                          matt.results2$BPSCORE_MAXBP > BPSCORE_MAXBP.q[4] |
                                          matt.results2$NUM_PREDICTED_BPS > NUM_PREDICTED_BPS.q[4] |
                                          matt.results2$SF1_HIGHESTSCORE_3SS > SF1_HIGHESTSCORE_3SS.q[4],"strong",
                                        "mixed"))

#version3
matt.results2$splicing <- ifelse(matt.results2$MAXENTSCR_HSAMODEL_5SS < MAXENTSCR_HSAMODEL_5SS.q[3] &
                                   matt.results2$MAXENTSCR_HSAMODEL_3SS < MAXENTSCR_HSAMODEL_3SS.q[3] &
                                   matt.results2$BPSCORE_MAXBP < BPSCORE_MAXBP.q[3] &
                                   matt.results2$SF1_HIGHESTSCORE_3SS < SF1_HIGHESTSCORE_3SS.q[3],"weak",
                                 
                                 ifelse(matt.results2$MAXENTSCR_HSAMODEL_5SS > MAXENTSCR_HSAMODEL_5SS.q[3] &
                                          matt.results2$MAXENTSCR_HSAMODEL_3SS > MAXENTSCR_HSAMODEL_3SS.q[3] &
                                          matt.results2$BPSCORE_MAXBP > BPSCORE_MAXBP.q[3] &
                                          matt.results2$SF1_HIGHESTSCORE_3SS > SF1_HIGHESTSCORE_3SS.q[3],"strong",
                                        "mixed"))

#version4
matt.results2$splicing <- ifelse(matt.results2$MAXENTSCR_HSAMODEL_5SS < MAXENTSCR_HSAMODEL_5SS.q[3] &
                                   matt.results2$BPSCORE_MAXBP < BPSCORE_MAXBP.q[3] &
                                   matt.results2$SF1_HIGHESTSCORE_3SS < SF1_HIGHESTSCORE_3SS.q[3],"weak",
                                 
                                 ifelse(matt.results2$MAXENTSCR_HSAMODEL_5SS > MAXENTSCR_HSAMODEL_5SS.q[3] &
                                          matt.results2$BPSCORE_MAXBP > BPSCORE_MAXBP.q[3] &
                                          matt.results2$SF1_HIGHESTSCORE_3SS > SF1_HIGHESTSCORE_3SS.q[3],"strong",
                                        "mixed"))

#version5
matt.results2$splicing <- ifelse(matt.results2$MAXENTSCR_HSAMODEL_5SS < MAXENTSCR_HSAMODEL_5SS.q[3] &
                                   matt.results2$BPSCORE_MAXBP < BPSCORE_MAXBP.q[3] &
                                   matt.results2$NUM_PREDICTED_BPS < NUM_PREDICTED_BPS.q[3],"weak",
                                 
                                 ifelse(matt.results2$MAXENTSCR_HSAMODEL_5SS > MAXENTSCR_HSAMODEL_5SS.q[3] &
                                          matt.results2$BPSCORE_MAXBP > BPSCORE_MAXBP.q[3] &
                                          matt.results2$NUM_PREDICTED_BPS > NUM_PREDICTED_BPS.q[3],"strong",
                                        "mixed"))

#version6
matt.results2$splicing <- ifelse(matt.results2$MAXENTSCR_HSAMODEL_5SS < MAXENTSCR_HSAMODEL_5SS.q[2] |
                                   matt.results2$BPSCORE_MAXBP < BPSCORE_MAXBP.q[2] |
                                   matt.results2$NUM_PREDICTED_BPS < NUM_PREDICTED_BPS.q[2],"weak",
                                 
                                 ifelse(matt.results2$MAXENTSCR_HSAMODEL_5SS > MAXENTSCR_HSAMODEL_5SS.q[4] |
                                          matt.results2$BPSCORE_MAXBP > BPSCORE_MAXBP.q[4] |
                                          matt.results2$NUM_PREDICTED_BPS > NUM_PREDICTED_BPS.q[4],"strong",
                                        "mixed"))


table(matt.results2$splicing)
matt.results2 <- dplyr::filter(matt.results2,splicing != "mixed")



#----------------------------------------------------------------------------------------------------------------
#   #turn MATT results into a GRAnges object
#----------------------------------------------------------------------------------------------------------------
matt.gr <- makeGRangesFromDataFrame(matt.results2,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=FALSE,
                                    seqinfo=NULL,
                                    seqnames.field=c("seqnames"),
                                    start.field=c("start"),
                                    end.field=c("end"),
                                    strand.field=c("strand"),
                                    starts.in.df.are.0based=FALSE)
names(matt.gr) <- matt.gr$INTRON_ID
matt.gr <- resize(matt.gr,width=1L,fix="end")
table(matt.gr$splicing)
length(unique( names(matt.gr))) == length( names(matt.gr))


#----------------------------------------------------------------------------------------------------------------
#   #calculate heatmaps
#----------------------------------------------------------------------------------------------------------------

#define parameters for heatmaps
span <- 200.5
step <- 1

#calculate heatmaps for splice donor sites
counts.donors <- SummitHeatmap(matt.gr,bamFiles,bamNames,span,step,useCPM=TRUE,strand=1,read2pos=0,readShiftSize=0)
#counts.RNA <- SummitHeatmap(matt.gr,bamFilesR,bamNamesR,span,step,useCPM=TRUE,strand=2,read2pos=0,readShiftSize=0)
#counts.Ribo <- SummitHeatmap(matt.gr,bamFilesRi,bamNamesRi,span,step,useCPM=TRUE,strand=1,read2pos=0,readShiftSize=0)


#----------------------------------------------------------------------------------------------------------------
#   #combine replicates for heatmaps
#----------------------------------------------------------------------------------------------------------------
sampleList <- list(
  Nrde2=c("Nrde2WT_rep1","Nrde2WT_rep2"),
  Ccdc174=c("Ccdc174_rep1","Ccdc174_rep2"),
  Eif4a3=c("Eif4a3_rep1","Eif4a3_rep2")
)
counts.donors.means <- SummarizeHeatmaps(counts.donors,sampleList)

#----------------------------------------------------------------------------------------------------------------
#   #make cumulative plots Eif4a3, Nrde2, Ccdc174
#----------------------------------------------------------------------------------------------------------------
plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A")

cumu.counts <- CumulativePlots(
  counts.donors.means[c(1,2,3)],
  bamNames=names(counts.donors.means[c(1,2,3)]),
  span = span,
  step = step,
  summarizing = "mean",
  plot = FALSE,
  confInterval = 0.95,
  overlapNames = names(matt.gr)[matt.gr$splicing=="weak"]
)


cumu.counts2 <- cumu.counts %>% pivot_longer(contains("overlap"),names_to = c(".value","overlap"),names_sep="_")
cumu.counts2$name <- factor(cumu.counts2$name,levels=c("Nrde2","Ccdc174","Eif4a3"),labels=c("Nrde2","Ccdc174","Eif4a3"))
cumu.counts2$splicing <- ifelse(cumu.counts2$overlap=="overlap1",sprintf("%d weak",table(matt.gr$splicing)[2]),
                               sprintf("%d strong",table(matt.gr$splicing)[1]))

p <- ggplot(cumu.counts2,aes(x=position, y=mean,col=name,linetype=splicing)) + geom_smooth(aes(ymin=ci.lower,ymax=ci.upper,fill=name,color=name),stat="identity") +
  facet_wrap(vars(name))
p <- p + theme_classic() + ylab("log2(cpm)") + 
  scale_color_manual(values=(plotcols[c(1,5,6)]),labels = names(counts.donors.means)[c(1,2,3)]) + 
  scale_fill_manual(values=(plotcols[c(1,5,6)]),labels = names(counts.donors.means)[c(1,2,3)]) 
p

ggsave("Figures/metaplots_at_strong_vs_weak_5p-splice_junctions_wholeUniqueReads_reps_combined_Nrde2_Ccdc174_Eif4a3_v5_210429.pdf",device="pdf",height=4,width=12)


#----------------------------------------------------------------------------------------------------------------
#   #define intron length and GC content features####
#----------------------------------------------------------------------------------------------------------------

matt.results2 <- dplyr::filter(matt.res,!is.na(INTRON_GCC),!is.na(INTRON_LENGTH))

INTRON_GCC.q <- quantile(matt.results2$INTRON_GCC)
INTRON_LENGTH.q <- quantile(matt.results2$INTRON_LENGTH)


#version1
matt.results2$intron.type <- ifelse(matt.results2$INTRON_LENGTH < INTRON_LENGTH.q[3] &
                                   matt.results2$INTRON_GCC > INTRON_GCC.q[3],"short & high GC",
                                 
                                 ifelse(matt.results2$INTRON_LENGTH > INTRON_LENGTH.q[3] &
                                          matt.results2$INTRON_GCC < INTRON_GCC.q[3],"long and low GC",
                                        "mixed"))


table(matt.results2$intron.type)
matt.results2 <- dplyr::filter(matt.results2,intron.type != "mixed")



#----------------------------------------------------------------------------------------------------------------
#   #turn MATT results into a GRAnges object
#----------------------------------------------------------------------------------------------------------------
matt.gr <- makeGRangesFromDataFrame(matt.results2,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=FALSE,
                                    seqinfo=NULL,
                                    seqnames.field=c("seqnames"),
                                    start.field=c("start"),
                                    end.field=c("end"),
                                    strand.field=c("strand"),
                                    starts.in.df.are.0based=FALSE)
names(matt.gr) <- matt.gr$INTRON_ID
matt.gr <- resize(matt.gr,width=1L,fix="end")
table(matt.gr$intron.type)
length(unique( names(matt.gr))) == length( names(matt.gr))


#----------------------------------------------------------------------------------------------------------------
#   #calculate heatmaps
#----------------------------------------------------------------------------------------------------------------

#define parameters for heatmaps
span <- 200.5
step <- 1

#calculate heatmaps for splice donor sites
counts.donors <- SummitHeatmap(matt.gr,bamFiles,bamNames,span,step,useCPM=TRUE,strand=1,read2pos=0,readShiftSize=0)
#counts.RNA <- SummitHeatmap(matt.gr,bamFilesR,bamNamesR,span,step,useCPM=TRUE,strand=2,read2pos=0,readShiftSize=0)
#counts.Ribo <- SummitHeatmap(matt.gr,bamFilesRi,bamNamesRi,span,step,useCPM=TRUE,strand=1,read2pos=0,readShiftSize=0)


#----------------------------------------------------------------------------------------------------------------
#   #combine replicates for heatmaps
#----------------------------------------------------------------------------------------------------------------
sampleList <- list(
  Nrde2=c("Nrde2WT_rep1","Nrde2WT_rep2"),
  Ccdc174=c("Ccdc174_rep1","Ccdc174_rep2"),
  Eif4a3=c("Eif4a3_rep1","Eif4a3_rep2")
)
counts.donors.means <- SummarizeHeatmaps(counts.donors,sampleList)

#----------------------------------------------------------------------------------------------------------------
#   #make cumulative plots Eif4a3, Nrde2, Ccdc174
#----------------------------------------------------------------------------------------------------------------
plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A")

cumu.counts <- CumulativePlots(
  counts.donors.means[c(1,2,3)],
  bamNames=names(counts.donors.means[c(1,2,3)]),
  span = span,
  step = step,
  summarizing = "mean",
  plot = FALSE,
  confInterval = 0.95,
  overlapNames = names(matt.gr)[matt.gr$intron.type=="short & high GC"]
)


cumu.counts2 <- cumu.counts %>% pivot_longer(contains("overlap"),names_to = c(".value","overlap"),names_sep="_")
cumu.counts2$name <- factor(cumu.counts2$name,levels=c("Nrde2","Ccdc174","Eif4a3"),labels=c("Nrde2","Ccdc174","Eif4a3"))
cumu.counts2$splicing <- ifelse(cumu.counts2$overlap=="overlap1",sprintf("%d short & high GC",table(matt.gr$intron.type)[2]),
                                sprintf("%d long and low GC",table(matt.gr$intron.type)[1]))

p <- ggplot(cumu.counts2,aes(x=position, y=mean,col=name,linetype=splicing)) + geom_smooth(aes(ymin=ci.lower,ymax=ci.upper,fill=name,color=name),stat="identity") +
  facet_wrap(vars(name))
p <- p + theme_classic() + ylab("log2(cpm)") + 
  scale_color_manual(values=(plotcols[c(1,5,6)]),labels = names(counts.donors.means)[c(1,2,3)]) + 
  scale_fill_manual(values=(plotcols[c(1,5,6)]),labels = names(counts.donors.means)[c(1,2,3)]) 
p

ggsave("Figures/metaplots_at_shortGCrich_vs_longGCpoor_introns_wholeUniqueReads_reps_combined_Nrde2_Ccdc174_Eif4a3_v1_210429.pdf",device="pdf",height=4,width=12)

