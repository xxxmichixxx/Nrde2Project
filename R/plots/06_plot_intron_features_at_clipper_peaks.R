
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
#   #plot intron features####
#----------------------------------------------------------------------------------------------------------------
matt.res1 <- matt.res[,c(11:13,31,32,40,27,17,22)]

matt.res1$INTRON_LENGTH <- log2(matt.res1$INTRON_LENGTH)
matt.res2 <- matt.res1[matt.res1$Nrde2.peak=="yes",]
matt.res2$peak <- "Nrde2"
matt.res3 <- matt.res1[matt.res1$Ccdc174.peak=="yes",]
matt.res3$peak <- "Ccdc174"
matt.res4 <- matt.res1[matt.res1$Eif4a3.peak=="yes",]
matt.res4$peak <- "Eif4a3"
matt.res5 <- rbind(matt.res2,matt.res3,matt.res4)


#define colors
plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A")

plotcols2 <- plotcols[c(5,6,1)]

#density plots
p1 <- ggplot(matt.res5,aes(x=MAXENTSCR_HSAMODEL_5SS,fill=peak,color=peak)) +geom_density(alpha=0.05) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2) + xlim(c(-15,15))

p2 <- ggplot(matt.res5,aes(x=MAXENTSCR_HSAMODEL_3SS,fill=peak,color=peak)) +geom_density(alpha=0.05) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2) + xlim(c(-15,15))

p3 <- ggplot(matt.res5,aes(x=BPSCORE_MAXBP,fill=peak,color=peak)) +geom_density(alpha=0.05) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2) + xlim(c(-3,3))

p4 <- ggplot(matt.res5,aes(x=SF1_HIGHESTSCORE_3SS,fill=peak,color=peak)) +geom_density(alpha=0.05) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)

p5 <- ggplot(matt.res5,aes(x=INTRON_LENGTH,fill=peak,color=peak)) +geom_density(alpha=0.05) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)

p6 <- ggplot(matt.res5,aes(x=INTRON_GCC,fill=peak,color=peak)) +geom_density(alpha=0.05) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)

p1 + p2 + p3 + p4 + p5 + p6

ggsave("Figures/intron_features_at_peaks_density_v3_210429.pdf",device="pdf",height=8,width=14)


#ridgeline plots
library(ggridges)

p1 <- ggplot(matt.res5,aes(x=MAXENTSCR_HSAMODEL_5SS,fill=peak,y=peak)) + geom_density_ridges(scale = 6,alpha=0.8) + 
  theme_minimal()  + scale_fill_manual(values=plotcols2) 
p1 <- p1 + scale_y_discrete(limits=c("Nrde2", "Ccdc174", "Eif4a3"))

p2 <- ggplot(matt.res5,aes(x=MAXENTSCR_HSAMODEL_3SS,fill=peak,y=peak)) + geom_density_ridges(scale = 6,alpha=0.8) + 
  theme_minimal()  + scale_fill_manual(values=plotcols2) 
p2 <- p2 + scale_y_discrete(limits=c("Nrde2", "Ccdc174", "Eif4a3"))

p3 <- ggplot(matt.res5,aes(x=BPSCORE_MAXBP,fill=peak,y=peak)) + geom_density_ridges(scale = 6,alpha=0.8) + 
  theme_minimal()  + scale_fill_manual(values=plotcols2) 
p3 <- p3 + scale_y_discrete(limits=c("Nrde2", "Ccdc174", "Eif4a3"))

p4 <- ggplot(matt.res5,aes(x=SF1_HIGHESTSCORE_3SS,fill=peak,y=peak)) + geom_density_ridges(scale = 6,alpha=0.8) + 
  theme_minimal()  + scale_fill_manual(values=plotcols2) 
p4 <- p4 + scale_y_discrete(limits=c("Nrde2", "Ccdc174", "Eif4a3"))

p5 <- ggplot(matt.res5,aes(x=INTRON_LENGTH,fill=peak,y=peak)) + geom_density_ridges(scale = 6,alpha=0.8) + 
  theme_minimal()  + scale_fill_manual(values=plotcols2) 
p5 <- p5 + scale_y_discrete(limits=c("Nrde2", "Ccdc174", "Eif4a3"))

p6 <- ggplot(matt.res5,aes(x=INTRON_GCC,fill=peak,y=peak)) + geom_density_ridges(scale = 6,alpha=0.8) + 
  theme_minimal()  + scale_fill_manual(values=plotcols2) 
p6 <- p6 + scale_y_discrete(limits=c("Nrde2", "Ccdc174", "Eif4a3"))

p1 + p2 + p3 + p4 + p5 + p6

ggsave("Figures/intron_features_at_peaks_ridgeline.pdf",device="pdf",height=8,width=14)



#violin plots
p1 <- ggplot(matt.res5,aes(y=MAXENTSCR_HSAMODEL_5SS,x=peak,fill=peak,color=peak)) + geom_violin(alpha=0.25) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)
p1 <- p1 + scale_x_discrete(limits=c("Nrde2", "Ccdc174", "Eif4a3"))

p2 <- ggplot(matt.res5,aes(y=MAXENTSCR_HSAMODEL_3SS,x=peak,fill=peak,color=peak)) + geom_violin(alpha=0.25) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)
p2 <- p2 + scale_x_discrete(limits=c("Nrde2", "Ccdc174", "Eif4a3"))

p3 <- ggplot(matt.res5,aes(y=BPSCORE_MAXBP,x=peak,fill=peak,color=peak)) + geom_violin(alpha=0.25) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)
p3 <- p3 + scale_x_discrete(limits=c("Nrde2", "Ccdc174", "Eif4a3"))

p4 <- ggplot(matt.res5,aes(y=INTRON_LENGTH,x=peak,fill=peak,color=peak)) + geom_violin(alpha=0.25) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)
p4 <- p4 + scale_x_discrete(limits=c("Nrde2", "Ccdc174", "Eif4a3"))

p5 <- ggplot(matt.res5,aes(y=INTRON_GCC,x=peak,fill=peak,color=peak)) + geom_violin(alpha=0.25) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)
p5 <- p5 + scale_x_discrete(limits=c("Nrde2", "Ccdc174", "Eif4a3"))

p1 + p2 + p3 + p4 + p5

ggsave("Figures/intron_features_at_peaks_violin.pdf",device="pdf",height=8,width=14)

#letter-value plots
library(lvplot)
p1 <- ggplot(matt.res5,aes(y=MAXENTSCR_HSAMODEL_5SS,x=peak,fill=peak,color=peak)) + geom_lv(alpha=0.25) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)
p1 <- p1 + scale_x_discrete(limits=c("Nrde2", "Ccdc174", "Eif4a3"))

p2 <- ggplot(matt.res5,aes(y=MAXENTSCR_HSAMODEL_3SS,x=peak,fill=peak,color=peak)) + geom_lv(alpha=0.25) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)
p2 <- p2 + scale_x_discrete(limits=c("Nrde2", "Ccdc174", "Eif4a3"))

p3 <- ggplot(matt.res5,aes(y=BPSCORE_MAXBP,x=peak,fill=peak,color=peak)) + geom_lv(alpha=0.25) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)
p3 <- p3 + scale_x_discrete(limits=c("Nrde2", "Ccdc174", "Eif4a3"))

p4 <- ggplot(matt.res5,aes(y=INTRON_LENGTH,x=peak,fill=peak,color=peak)) + geom_lv(alpha=0.25) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)
p4 <- p4 + scale_x_discrete(limits=c("Nrde2", "Ccdc174", "Eif4a3"))

p5 <- ggplot(matt.res5,aes(y=INTRON_GCC,x=peak,fill=peak,color=peak)) + geom_lv(alpha=0.25) + 
  theme_classic()  + scale_color_manual(values=plotcols2) + scale_fill_manual(values=plotcols2)
p5 <- p5 + scale_x_discrete(limits=c("Nrde2", "Ccdc174", "Eif4a3"))

p1 + p2 + p3 + p4 + p5

#TSNE
library(Rtsne)
row.names(matt.res1) <- matt.res$INTRON_ID


TDAT <- unique(matt.res1[complete.cases(matt.res1[,4:8]),4:8])


matt.tsne <- Rtsne(TDAT, dims = 2, initial_dims = 50,
      perplexity = 500, theta = 0.5, 
      pca = FALSE, partial_pca = FALSE, max_iter = 1000,
      verbose = getOption("verbose", FALSE), is_distance = FALSE,
      Y_init = NULL, pca_center = TRUE, pca_scale = FALSE,
      normalize = TRUE, num_threads = 2)

matt.tsne.res <- data.frame(matt.tsne$Y)
colnames(matt.tsne.res) <- c("dim1","dim2")

#add peak info
matt.tsne.res$INTRON_ID <- row.names(TDAT)
matt.tsne.res <- left_join(matt.tsne.res,matt.res[,c(11:14)],by="INTRON_ID")

#plot
ggplot(matt.tsne.res,aes(x=dim1,y=dim2,col=Nrde2.peak)) + geom_point(alpha=0.1)

#PCA
pca <- prcomp(TDAT[,1:3])
PCA <- data.frame(pca$x)
PCA$INTRON_ID <- row.names(PCA)
#add peak info
PCA <- left_join(PCA,matt.res[,c(11:14)],by="INTRON_ID")

ggplot(PCA,aes(x=PC1,y=PC2,col=Nrde2.peak)) + geom_point(alpha=0.5)

ggplot(matt.res1,aes(x=INTRON_LENGTH,y=INTRON_GCC,col=Nrde2.peak)) + geom_point(alpha=0.5) + facet_wrap(vars(Nrde2.peak))



