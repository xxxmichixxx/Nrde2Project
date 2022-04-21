rm(list=ls())

library(tidyverse)
library(Rsubread)
library(GenomicFeatures)
library(GenomicRanges)
library(MiniChip)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)

setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")

#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#   #load RMATS results of alternative 5' SS
#----------------------------------------------------------------------------------------------------------------

comp <- c("Nrde2KO_vs_WT","Nrde2d200_vs_WT","Nrde2D174R_vs_WT","Nrde2KO_vs_WT_CHX","Ccdc174KO_vs_WT","Ccdc174KO_vs_WT_CHX","Mtr4KO_vs_WT","WT_CHX_vs_WT")
i=1
#for (i in seq_along(comp)){

A5SS.MATS.JCEC <- read.table(sprintf("rmats_analysis/%s/A5SS.MATS.JCEC.txt",comp[i]),header=TRUE)
#short exon 5'SS
A5SS.gr <- makeGRangesFromDataFrame(A5SS.MATS.JCEC,
                                    keep.extra.columns=TRUE,
                                    seqinfo=NULL,
                                    seqnames.field=c("chr"),
                                    start.field=c("shortES"),
                                    end.field=c("shortEE"),
                                    strand.field=c("strand"),
                                    ignore.strand=FALSE,
                                    starts.in.df.are.0based=TRUE)
#long exon 5'SS
A5SS.gr.long <- makeGRangesFromDataFrame(A5SS.MATS.JCEC,
                                    keep.extra.columns=TRUE,
                                    seqinfo=NULL,
                                    seqnames.field=c("chr"),
                                    start.field=c("longExonStart_0base"),
                                    end.field=c("longExonEnd"),
                                    strand.field=c("strand"),
                                    ignore.strand=FALSE,
                                    starts.in.df.are.0based=TRUE)

A5SS.gr$upregulated <- ifelse(A5SS.gr$FDR < 0.05 & A5SS.gr$IncLevelDifference > 0, "yes","no")
table(A5SS.gr$upregulated)


A5SS.gr2 <- resize(A5SS.gr,1L,fix="end")
A5SS.gr.long2 <- resize(A5SS.gr.long,1L,fix="end")

names(A5SS.gr2) <- paste(seqnames(A5SS.gr2),start(A5SS.gr2),A5SS.gr2$ID,sep="_")
names(A5SS.gr.long2) <- paste(seqnames(A5SS.gr.long2),start(A5SS.gr.long2),A5SS.gr.long2$ID,sep="_")

length(unique( names(A5SS.gr2))) == length( names(A5SS.gr2))

#----------------------------------------------------------------------------------------------------------------
  #   #resize the spliceDonor sites to et  -3 to +6 bp around junctions (donors)
  #----------------------------------------------------------------------------------------------------------------
#short exons
spliceDonors <- A5SS.gr2[seqnames(A5SS.gr2) != "chrJH584304.1"]
spliceDonors.plus <- spliceDonors[strand(spliceDonors)=="+"]
start(spliceDonors.plus) <- start(spliceDonors.plus) -2
end(spliceDonors.plus) <- start(spliceDonors.plus) +8
getSeq(BSgenome.Mmusculus.UCSC.mm10,spliceDonors.plus)

spliceDonors.minus <- spliceDonors[strand(spliceDonors)=="-"]
start(spliceDonors.minus) <- start(spliceDonors.minus) -6
end(spliceDonors.minus) <- start(spliceDonors.minus) +8
getSeq(BSgenome.Mmusculus.UCSC.mm10,spliceDonors.minus)

spliceDonors8n <- c(spliceDonors.plus,spliceDonors.minus)

#long exons
spliceDonors <- A5SS.gr.long2[seqnames(A5SS.gr.long2) != "chrJH584304.1"]
spliceDonors.plus <- spliceDonors[strand(spliceDonors)=="+"]
start(spliceDonors.plus) <- start(spliceDonors.plus) -2
end(spliceDonors.plus) <- start(spliceDonors.plus) +8
getSeq(BSgenome.Mmusculus.UCSC.mm10,spliceDonors.plus)

spliceDonors.minus <- spliceDonors[strand(spliceDonors)=="-"]
start(spliceDonors.minus) <- start(spliceDonors.minus) -6
end(spliceDonors.minus) <- start(spliceDonors.minus) +8
getSeq(BSgenome.Mmusculus.UCSC.mm10,spliceDonors.minus)

spliceDonors8n.long <- c(spliceDonors.plus,spliceDonors.minus)
#----------------------------------------------------------------------------------------------------------------
#   # save the sequecnes as fasta file and run MaxEntScan, 
#----------------------------------------------------------------------------------------------------------------
spliceDonors8n.seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10,spliceDonors8n)
spliceDonors8n.long.seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10,spliceDonors8n.long)

#remove Ns and save fasta:short
inxN <- which(letterFrequency(spliceDonors8n.seq,"N")== 0)
writeXStringSet(spliceDonors8n.seq[inxN], 'spliceDonors9n.seq_RMATSshortExons.fa')
spliceDonors9n <- spliceDonors8n[inxN]

#remove Ns and save fasta:long
inxN <- which(letterFrequency(spliceDonors8n.long.seq,"N")== 0)
writeXStringSet(spliceDonors8n.long.seq[inxN], 'spliceDonors9n.seq_RMATSlongExons.fa')
spliceDonors9n.long <- spliceDonors8n.long[inxN]

#run MaxEntScan
cd /tungstenfs/scratch/gbuehler/michi/Tools/MaxEntScan/
perl score5.pl /tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/spliceDonors9n.seq_RMATSshortExons.fa > /tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/spliceDonors9n.seq_RMATSshortExons.maxentscan.txt
perl score5.pl /tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/spliceDonors9n.seq_RMATSlongExons.fa > /tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/spliceDonors9n.seq_RMATSlongExons.maxentscan.txt

spliceDonors9n.maxentscan <- read.table("spliceDonors9n.seq_RMATSshortExons.maxentscan.txt",sep="\t")
spliceDonors9n.long.maxentscan <- read.table("spliceDonors9n.seq_RMATSlongExons.maxentscan.txt",sep="\t")

plot(density(spliceDonors9n.maxentscan$V2))
lines(density(spliceDonors9n.long.maxentscan$V2),col="red")

spliceDonors9n$MaxEntScan <- spliceDonors9n.maxentscan$V2
spliceDonors9n.long$MaxEntScan.long <- spliceDonors9n.long.maxentscan$V2

res <- data.frame(ID=names(spliceDonors9n),upregulated=spliceDonors9n$upregulated,
                  short.score=spliceDonors9n$MaxEntScan,long.score=spliceDonors9n.long$MaxEntScan.long)
ggplot(res,aes(x=short.score,y=long.score,col=upregulated)) + geom_point()

res$FC <- res$long.score - res$short.score
ggplot(res,aes(x=upregulated,y=FC,col=upregulated)) + geom_boxplot()

#cpompare to long intron length
spliceDonors9n.plus <- spliceDonors9n[strand(spliceDonors9n)=="+"]
spliceDonors9n.plus$intron.length <- spliceDonors9n.plus$flankingES - end(spliceDonors9n.plus)

spliceDonors9n.minus <- spliceDonors9n[strand(spliceDonors9n)=="-"]
spliceDonors9n.minus$intron.length <- start(spliceDonors9n.minus) - spliceDonors9n.minus$flankingEE

res$intron.length <- log2(c(spliceDonors9n.plus$intron.length,spliceDonors9n.minus$intron.length))
ggplot(res,aes(x=upregulated,y=intron.length,col=upregulated)) + geom_violin()
