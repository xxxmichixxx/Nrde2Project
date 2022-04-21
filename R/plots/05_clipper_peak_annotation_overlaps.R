rm(list=ls())

library(tidyverse)
library(GenomicFeatures)
library(GenomicRanges)
library(MiniChip)
library(ComplexHeatmap)
library(wesanderson)
library(viridis)
library(UpSetR)
library(EnsDb.Mmusculus.v79)

setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")

#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#   #load peaks
#----------------------------------------------------------------------------------------------------------------
print(load("clipper_peaks_5readsmin_200202.RData"))


clipper.peaks$ID <- names(clipper.peaks)
samples <- unique(clipper.peaks$sample)

#-----------------------------------------------------------
# intersect peaks with Gencode and repeat annos####
#-----------------------------------------------------------
#exon
#intron
#5'SS
#3'SS
#snRNA
#snoRNA
#rRNA

#annotate peaks with genes
edb <- EnsDb.Mmusculus.v79
edbg <- genes(edb)
#txdb <- makeTxDbFromGFF("RNAseq/stringtie_all_lib643_and_lib2191and2447_groupwise_transcripts.annotated.gtf",
#                     format=c("gtf"))
txdb <- makeTxDbFromGFF("/work/gbioinfo/DB/GENCODE/Mouse/release_M23/gencode.vM23.annotation.gtf",
                        format=c("gtf"))
genes <- genes(txdb)
genes$gene_id <- matrix(unlist(strsplit(genes$gene_id,".",fixed=TRUE)),ncol=2,byrow = TRUE)[,1]
values(genes) <- dplyr::left_join(data.frame(values(genes)),data.frame(values(edbg)),by="gene_id")
genes2 <- genes[is.na(genes$gene_biotype)==FALSE]


protein_coding <- genes2[genes2$gene_biotype=="protein_coding"]
lincRNA <- genes2[genes2$gene_biotype=="lincRNA"]
snRNA <- genes2[genes2$gene_biotype=="snRNA"]
snoRNA <- genes2[genes2$gene_biotype=="snoRNA"]
rRNA <- genes2[genes2$gene_biotype=="rRNA"]

exons <- unlist(exonsBy(txdb,by = "tx"))
introns <- unlist(intronsByTranscript(txdb))
SS5 <- resize(introns,1L,fix="start")
SS3 <- resize(introns,1L,fix="end")

#annotate each peak set seperately and make upsetplot
for (i in seq_along(samples)){
  peaks2 <- clipper.peaks[clipper.peaks$sample==samples[i]]
  
  #annotate peaks
  peaks2$protein_coding <- ifelse(overlapsAny(peaks2,protein_coding)==TRUE,1,0)
  peaks2$lincRNA <- ifelse(overlapsAny(peaks2,lincRNA)==TRUE,1,0)
  peaks2$snRNA <- ifelse(overlapsAny(peaks2,snRNA)==TRUE,1,0)
  peaks2$snoRNA <- ifelse(overlapsAny(peaks2,snoRNA)==TRUE,1,0)
  peaks2$rRNA <- ifelse(overlapsAny(peaks2,rRNA)==TRUE,1,0)
  peaks2$exon <- ifelse(overlapsAny(peaks2,exons)==TRUE,1,0)
  peaks2$intron <- ifelse(overlapsAny(peaks2,introns)==TRUE,1,0)
  peaks2$FivePSS <- ifelse(overlapsAny(peaks2,SS5)==TRUE,1,0)
  peaks2$ThreePSS <- ifelse(overlapsAny(peaks2,SS3)==TRUE,1,0)
  
  
  peaks2upset <- data.frame(mcols(peaks2))
  
  pdf(sprintf("Figures/upset_plot_clipper_peaks_5readsmin_%s_peak_and_gene_overlaps.pdf",samples[i]),height=6,width=9)
  upset(peaks2upset[,12:ncol(peaks2upset)], nsets = 4, number.angles = 30, point.size = 3.5, line.size = 2, 
        mainbar.y.label = samples[i] , sets.x.label = "number of genes",
        text.scale = c(1.5, 1.5, 1.5, 1.5, 2, 1.5))
  dev.off()
}

#annotate all peaks and make barplogts
peaks2 <- clipper.peaks
  #annotate peaks
  peaks2$protein_coding <- ifelse(overlapsAny(peaks2,protein_coding)==TRUE,1,0)
  peaks2$lincRNA <- ifelse(overlapsAny(peaks2,lincRNA)==TRUE,1,0)
  peaks2$snRNA <- ifelse(overlapsAny(peaks2,snRNA)==TRUE,1,0)
  peaks2$snoRNA <- ifelse(overlapsAny(peaks2,snoRNA)==TRUE,1,0)
 # peaks2$rRNA <- ifelse(overlapsAny(peaks2,rRNA)==TRUE,1,0)
  peaks2$exon <- ifelse(overlapsAny(peaks2,exons)==TRUE,1,0)
  peaks2$intron <- ifelse(overlapsAny(peaks2,introns)==TRUE,1,0)
  peaks2$FivePSS <- ifelse(overlapsAny(peaks2,SS5)==TRUE,1,0)
  peaks2$ThreePSS <- ifelse(overlapsAny(peaks2,SS3)==TRUE,1,0)
  
  
  peaks2upset <- data.frame(mcols(peaks2))[,c(5,7:14)]
  peaks2upset2 <- pivot_longer(peaks2upset,cols = colnames(peaks2upset)[2:9],names_to = "type",values_to="overlap")
  
  #define colors
  plotcols1 <- wes_palette("Darjeeling1")
  plotcols2 <- wes_palette("Darjeeling2")
  
  plotcols <- c(plotcols2[c(2,4,3)],plotcols1[c(1,5,2)])
  plotcols2 <- plotcols[c(5,6,4,1,3,2)]
  
  #with Eif4a3
  ggplot(peaks2upset2[peaks2upset2$overlap==1,],aes(y=type,fill=sample)) +geom_bar() + facet_wrap(vars(sample),scales="free") +
    theme_classic() + scale_fill_manual(values=plotcols2)
  ggsave("Figures/barplots_of_peak_counts_over_gene_annotations_withoutrrna.pdf",device="pdf",height=7,width=10)
  #without Eif4a3
  #ggplot(peaks2upset2[peaks2upset2$overlap==1 & peaks2upset2$sample !="Eif4a3",],aes(y=type,fill=sample)) + geom_bar(position='dodge')
  ggplot(peaks2upset2[peaks2upset2$overlap==1,],aes(y=sample,fill=type)) + geom_bar(stat="count",position="fill")  +
    theme_classic() + scale_fill_manual(values=plotcols2)
  
  
  
#annotate expressed genes with peaks
print(load("splice_sites_expressed_and_notexpressed_transcripts_log2TPMabove0.1_and_genes_FPKMabove8below0.RData"))
genes2 <- exons3g
names(genes2) <- genes2$GeneID

plotlist <- list()
peak.genes.list <- list()

for (i in seq_along(samples)){
  peaks2 <- clipper.peaks[clipper.peaks$sample==samples[i]]
    peaks2genes <- findOverlaps(genes2,peaks2)
  peaks3 <- peaks2[to(peaks2genes)]
  genes3 <- genes2[from(peaks2genes)]
  values(genes3) <- c(values(genes3),values(peaks3))
  genes3 <- genes3[unique(genes3$GeneID)]
  
  genes2 <- genes2[unique(genes2$GeneID)]
  
  all.genes <- data.frame(protein_coding=length(na.omit(genes2$gene_biotype[genes2$gene_biotype=="protein_coding"])),
                          lincRNA=length(na.omit(genes2$gene_biotype[genes2$gene_biotype=="lincRNA"])),
                          snRNA=length(na.omit(genes2$gene_biotype[genes2$gene_biotype=="snRNA"])),
                          snoRNA=length(na.omit(genes2$gene_biotype[genes2$gene_biotype=="snoRNA"])),
                          rRNA=length(na.omit(genes2$gene_biotype[genes2$gene_biotype=="rRNA"])))
  
  peak.genes <- data.frame(protein_coding=length(na.omit(genes3$gene_biotype[genes3$gene_biotype=="protein_coding"])),
                           lincRNA=length(na.omit(genes3$gene_biotype[genes3$gene_biotype=="lincRNA"])),
                           snRNA=length(na.omit(genes3$gene_biotype[genes3$gene_biotype=="snRNA"])),
                           snoRNA=length(na.omit(genes3$gene_biotype[genes3$gene_biotype=="snoRNA"])),
                           rRNA=length(na.omit(genes3$gene_biotype[genes3$gene_biotype=="rRNA"])))
  
  
  
  peak.genes.df <- data.frame(t(data.frame(no = all.genes[1,]-peak.genes[1,],
                                           yes=peak.genes[1,]) ))
  peak.genes.df$gene_biotype <- matrix(unlist(strsplit(row.names(peak.genes.df),".",fixed=TRUE)),ncol=2,byrow = TRUE)[,2]
  peak.genes.df$peak <- matrix(unlist(strsplit(row.names(peak.genes.df),".",fixed=TRUE)),ncol=2,byrow = TRUE)[,1]
  peak.genes.df$number_of_genes <- peak.genes.df$X1
  peak.genes.df$sample <- samples[i]
  peak.genes.list[[i]] <- peak.genes.df[,2:5]
  
  #plotlist[[i]] <- ggplot(peak.genes.df,aes(x=gene_biotype,y=number_of_genes,fill=peak)) + geom_bar(stat="identity") + 
  #  coord_flip() + theme_classic() + ggtitle(names(peaks.list4)[i])
  plotlist[[i]] <- ggplot(peak.genes.df[peak.genes.df$gene_biotype != "protein_coding",],aes(x=gene_biotype,y=number_of_genes,fill=peak)) + geom_bar(stat="identity") + 
    coord_flip() + theme_classic() + ggtitle(samples[i])
}

#plotlist[[1]] / plotlist[[2]] / plotlist[[3]] / plotlist[[4]] / plotlist[[5]] /plotlist[[6]]
#ggsave("plots/number_of_expressed_genes_with_peaks_barplot_zoomin.png",height=10,width=15)

peak.genes <- do.call(rbind,peak.genes.list)
peak.genes <- dplyr::filter(peak.genes,gene_biotype != "rRNA")
peak.genes$sample <- as.factor(peak.genes$sample)
peak.genes <- mutate(peak.genes, sample = fct_relevel(sample, 
                                                      "Nrde2d200","Nrde2D174R","Nrde2","Ccdc174","Mtr4","Eif4a3")) 

ggplot(peak.genes,aes(x=sample,y=number_of_genes,fill=peak)) + geom_bar(stat="identity") + 
  coord_flip() + theme_classic() + facet_wrap(vars(gene_biotype),scales = "free") + scale_fill_manual(values=wes_palette("Darjeeling2",2, type="discrete"))
ggsave("Figures/number_of_expressed_genes_with_peaks_barplot_by_gene_biotype.png",height=7,width=8)


