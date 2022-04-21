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

#load peaks
print(load("clipper_peaks_5readsmin_200202.RData"))
clipper.peaks$ID <- names(clipper.peaks)
samples <- unique(clipper.peaks$sample)

#load introns
load("splice_sites_expressed_and_notexpressed_transcripts_log2TPMabove0_and_genes_FPKMabove0.RData")
names(spliceDonors) <- paste(names(spliceDonors),start(spliceDonors),sep="_")
names(spliceAcceptors) <- names(spliceDonors)
spliceDonors2 <- GenomicRanges::promoters(spliceDonors,upstream=100,downstream=1)


#----------------------------------------------------------------------------------------------------------------
#   #get intron features with MATT####
#----------------------------------------------------------------------------------------------------------------
#With table introns.tab describing human introns with columns START, END, SCAFF, STRAND, ENSEMBL_GID
introns3g.matt <- data.frame(START=start(introns3g),END=end(introns3g),SCAFF=seqnames(introns3g),
                             STRAND=strand(introns3g),ENSEMBL_GID=introns3g$GeneID
)
write.table(introns3g.matt,file="introns3g.matt.txt",sep="\t",col.names=TRUE,row.names=FALSE, append=FALSE, quote=FALSE)

#save introns as tab delimited file

#run matt
cd /tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/
  module load R-BioC/current

ln -s /tungstenfs/scratch/gbuehler/michi/Tools/matt/matt
./matt get_ifeatures explain

genomepath="/work/gbioinfo/DB/GENCODE/Mouse/release_M23"
#cat ${genomepath}/gencode.vM23.annotation.gtf | grep "protein_coding" > gencode.vM23.proteincoding.annotation.gtf

./matt get_ifeatures introns3g.matt.txt START END SCAFF STRAND ENSEMBL_GID ${genomepath}/gencode.vM23.annotation.gtf ${genomepath}/GRCm38.primary_assembly.genome.fa Mmus -notrbts > introns3g.matt.ifeatures.tab

introns3g.matt.res <- read.delim("introns3g.matt.ifeatures.tab",sep="\t",header=TRUE)

#----------------------------------------------------------------------------------------------------------------
#   #which 5' splice sites have Nrde2 at 5' splice site ####
#----------------------------------------------------------------------------------------------------------------

#Nrde2
peaks <- clipper.peaks[clipper.peaks$sample==samples[4]]
spliceDonors$Nrde2.peak <- ifelse(overlapsAny(spliceDonors2,peaks,
                                        minoverlap = 10),"yes" ,"no")
table(spliceDonors$Nrde2.peak)

#Ccdc174
peaks <- clipper.peaks[clipper.peaks$sample==samples[1]]
spliceDonors$Ccdc174.peak <- ifelse(overlapsAny(spliceDonors2,peaks,
                                              minoverlap = 10),"yes" ,"no")
table(spliceDonors$Ccdc174.peak)

#Eif4a3
peaks <- clipper.peaks[clipper.peaks$sample==samples[2]]
spliceDonors$Eif4a3.peak <- ifelse(overlapsAny(spliceDonors2,peaks,
                                                minoverlap = 10),"yes" ,"no")
table(spliceDonors$Eif4a3.peak)

matt.res <- cbind(data.frame(spliceDonors),introns3g.matt.res)
matt.res$ID <- names(spliceDonors)

#----------------------------------------------------------------------------------------------------------------
#   #add BCLIP read counts at 5' splice site ####
#----------------------------------------------------------------------------------------------------------------
#bamFiles <- c(list.files("/tungstenfs/scratch/gbuehler/flemrmat/data/2456/raw_data/", full.names=TRUE,pattern="*LCstripped.prinseq_Aligned.sortedByCoord.out.bam$"))[1:12]
#bamNames <- c("Ccdc174_rep1","Ccdc174_rep2","Eif4a3_rep1","Eif4a3_rep2","Mtr4_rep1","Mtr4_rep2","Nrde2WT_rep1","Nrde2WT_rep2",
              "Nrde2D174R_rep1","Nrde2D174R_rep2","Nrde2d200_rep1","Nrde2d200_rep2")
#spliceDonors.saf <- data.frame(GeneID= names(spliceDonors), Chr=seqnames(spliceDonors),
#                                Start=start(spliceDonors), End=end(spliceDonors),Strand=strand(spliceDonors))
#calculate bclip read counts overlapping donor sites 
#bclip.donors <- featureCounts(bamFiles,annot.ext=spliceDonors.saf,
 #                             useMetaFeatures=FALSE,allowMultiOverlap=TRUE,
#                              minOverlap=1,countMultiMappingReads=FALSE,fraction=TRUE,
#                              minMQS=255,strandSpecific=1,nthreads=10,verbose=FALSE,isPairedEnd=FALSE)
#bclip.donors2 <- data.frame(bclip.donors$counts)
#colnames(bclip.donors2) <- bamNames


write.table(matt.res,"FigureData/MATT_scores_on_all_introns_with_peak_overlap_info.txt",col.names=TRUE,sep="\t",row.names=FALSE,quote=FALSE,append=FALSE)
