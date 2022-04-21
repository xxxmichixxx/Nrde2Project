rm(list=ls())

setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")

#----------------------------------------------------------------------------------------------------------------
#   #use GViz to make Sashimi plots
#----------------------------------------------------------------------------------------------------------------
library(Gviz)
library(GenomicFeatures)
library(EnsDb.Mmusculus.v79)
options(ucscChromosomeNames=FALSE)
plotcols <- c("#046C9A","#5BBCD6","#C6CDF7","#FD6464",  "#F98400","#00A08A")
plotcols2 <- plotcols[c(5,6,1)]
#----------------------------------------------------------------------------------------------------------------
#   #load bam files
#----------------------------------------------------------------------------------------------------------------

#load BCLIP files
bamFiles <- c(list.files("/tungstenfs/scratch/gbuehler/flemrmat/data/2456/raw_data/", full.names=TRUE,pattern="*LCstripped.prinseq_Aligned.sortedByCoord.out.bam$"))[1:12]
bamNames <- c("Ccdc174_rep1","Ccdc174_rep2","Eif4a3_rep1","Eif4a3_rep2","Mtr4_rep1","Mtr4_rep2","Nrde2WT_rep1","Nrde2WT_rep2",
              "Nrde2D174R_rep1","Nrde2D174R_rep2","Nrde2d200_rep1","Nrde2d200_rep2")

#######merge bam files############
cd /tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/BCLIP/experiment3/merged_bam
bamdir="/tungstenfs/scratch/gbuehler/flemrmat/data/2456/raw_data/"
module load SAMtools/1.10-foss-2019b

samtools merge -f Eif4a3.bam $bamdir/2456_BCLIP_Eif4a3_1.LCstripped.prinseq_Aligned.sortedByCoord.out.bam $bamdir/2456_BCLIP_Eif4a3_2.LCstripped.prinseq_Aligned.sortedByCoord.out.bam
samtools merge -f Nrde2.bam $bamdir/2456_BCLIP_Nrde2_1.LCstripped.prinseq_Aligned.sortedByCoord.out.bam $bamdir/2456_BCLIP_Nrde2_2.LCstripped.prinseq_Aligned.sortedByCoord.out.bam
samtools merge -f Ccdc174.bam $bamdir/2456_BCLIP_Ccdc174_1.LCstripped.prinseq_Aligned.sortedByCoord.out.bam $bamdir/2456_BCLIP_Ccdc174_2.LCstripped.prinseq_Aligned.sortedByCoord.out.bam

samtools index Eif4a3.bam
samtools index Nrde2.bam
samtools index Ccdc174.bam
############################################

BCLIP_files <- c(list.files("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/BCLIP/experiment3/merged_bam", full.names=TRUE,pattern="*.bam$"))
BCLIP_files_names <- gsub("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/BCLIP/experiment3/merged_bam/","",BCLIP_files)
BCLIP_files_names <- gsub(".bam","",BCLIP_files_names)


#load hisat2 mapped RNAseq files
all.bamFiles <- list.files("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/RNAseq/hisat_groupwise",full.names = TRUE,pattern="bam$")               
bamFiles.Nrde2KO <- all.bamFiles[c(8,7)]
bamFiles.Ccdc174KO <- all.bamFiles[c(2,1)]
bamFiles.WT <- all.bamFiles[c(10,9)]

#the dTAG-SmE RNAseq is experiment 2701
all.bamFiles <- list.files("/tungstenfs/scratch/gbuehler/deepSeqRepos/bam/", full.names=TRUE,pattern="*bam$")
bamFilesR <- c(grep("2701F",all.bamFiles,value=TRUE))[5:8]

#----------------------------------------------------------------------------------------------------------------
#   #load genes and peaks
#----------------------------------------------------------------------------------------------------------------

#load txdb
txdb <- makeTxDbFromGFF("/work/gbioinfo/DB/GENCODE/Mouse/release_M23/gencode.vM23.annotation.gtf",
                        format=c("gtf"))
genes.txdb <- genes(txdb)

#to select gene to plot
EnsDb <- EnsDb.Mmusculus.v79
gene_names <- genes(EnsDb, columns=c("gene_name", "gene_biotype"))

#load BCLIP peaks
print(load("clipper_peaks_5readsmin_200202.RData"))

#----------------------------------------------------------------------------------------------------------------
#   #select plot region (by gene symbol)
#----------------------------------------------------------------------------------------------------------------

gene <- "Cdk2"
#gene <- "Tti1"

#axis limits for each gene
#Cdk2:
lim.RNA <- c(0,1500)
lim.BCLIP <- c(0,300)
#Tti1
lim.RNA <- c(0,600)
lim.BCLIP <- c(0,60)

#by gene name
geneinx <- which(gene_names$gene_name==gene)
gene.gr <- gene_names[geneinx]


afrom <- start(gene.gr) - 100
ato <- end(gene.gr) + 100
#chr <- seqnames(gene.gr)
chr <- paste("chr",seqnames(gene.gr),sep="")

strand <-  ifelse(strand(gene.gr) == "-",TRUE, FALSE) # True for - strand, FALSE for + strand

#----------------------------------------------------------------------------------------------------------------
#   #generate tracks
#----------------------------------------------------------------------------------------------------------------

#trasncripts track
txTr <- GeneRegionTrack(txdb, chromosome = chr,cex.legend=80,
                        start = afrom,  end = ato, size=8,fill = "darkgrey", col = NULL, 
                        transcriptAnnotation = "transcript",showId=FALSE)
#axis track
axisTrack <- GenomeAxisTrack(from=afrom, to =ato, fontsize=10,labelPos = "above")

#RNAseq tracks
Nrde2KO.TrackList <- list()
for (i in seq_along(bamFiles.Nrde2KO)){
  Nrde2KO.TrackList[[i]] <- AlignmentsTrack(
    bamFiles.Nrde2KO[i],type=c("coverage", "sashimi"),
    isPaired = FALSE,chromosome =chr,name="Nrde2KO",cex.legend=80,
    col=plotcols[1],col.sashimi=plotcols[1],coverageHeight=0.2,fill.coverage=plotcols[1],
    sashimiScore=5,mapq=255,ylim=lim.RNA)
}

Ccdc174KO.TrackList <- list()
for (i in seq_along(bamFiles.Ccdc174KO)){
  Ccdc174KO.TrackList[[i]] <- AlignmentsTrack(
    bamFiles.Ccdc174KO[i],type=c("coverage", "sashimi"),
    isPaired = FALSE,chromosome =chr,name="Ccdc174KO",cex.legend=80,
    col=plotcols[5],col.sashimi=plotcols[5],coverageHeight=0.2,fill.coverage=plotcols[5],
    sashimiScore=5,mapq=255,ylim=lim.RNA)
}

WT.TrackList <- list()
for (i in seq_along(bamFiles.WT)){
  WT.TrackList[[i]] <- AlignmentsTrack(
    bamFiles.WT[i],type=c("coverage", "sashimi"),
    isPaired = FALSE,chromosome =chr,name="WT",cex.legend=80,
    col="black",col.sashimi="black",coverageHeight=0.2,fill.coverage="black",
    sashimiScore=5,mapq=255,ylim=lim.RNA)
}

#BCLIP tracks
covH <- c(1,1,1)
BCLIP.TrackList <- list()
for (i in seq_along(BCLIP_files)){
  BCLIP.TrackList[[i]] <- AlignmentsTrack(
    BCLIP_files[i],type="coverage",
    isPaired = FALSE,chromosome =chr,name=BCLIP_files_names[i],cex.legend=80,
    col=plotcols2[i],coverageHeight=covH,minCoverageHeight=60,
    fill.coverage=plotcols2[i],mapq=255,ylim=lim.BCLIP)
}


#bclip peaks tracks
peak.TrackList <- list()
for (i in seq_along(BCLIP_files_names)){
  peaks <- clipper.peaks[seqnames(clipper.peaks)==chr  & start(clipper.peaks) > afrom &
                                 end(clipper.peaks) < ato & clipper.peaks$sample==BCLIP_files_names[i]]
  peak.TrackList[[i]] <- AnnotationTrack(peaks,chromosome = chr,cex.legend=80,start = afrom,  end = ato, size=18,
                                         name=BCLIP_files_names[i],showFeatureId=FALSE,fill = plotcols2[i], col = "#252525")
}
#----------------------------------------------------------------------------------------------------------------
#   #plots
#----------------------------------------------------------------------------------------------------------------

pdf(sprintf("Figures/%s_Nrde2-hisat_RNAseq_sashimi_plot_fixed_ylim.pdf",gene),height=20,width=30)
plotTracks(c(peak.TrackList[c(3)],BCLIP.TrackList[c(3)],peak.TrackList[c(2)],BCLIP.TrackList[c(2)],Nrde2KO.TrackList[1:2],WT.TrackList[1:2],txTr,axisTrack), from = afrom, to = ato, 
           chromosome = chr,reverseStrand = strand,sizes=c(0.5,1,0.5,1,2,2,2,2,1,1))
dev.off()

pdf(sprintf("Figures/%s_Ccdc174-hisat_RNAseq_sashimi_plot_fixed_ylim.pdf",gene),height=20,width=30)
plotTracks(c(peak.TrackList[c(1)],BCLIP.TrackList[c(1)],peak.TrackList[c(2)],BCLIP.TrackList[c(2)],Ccdc174KO.TrackList[1:2],WT.TrackList[1:2],txTr,axisTrack), from = afrom, to = ato, 
           chromosome = chr, reverseStrand = strand,sizes=c(0.5,1,0.5,1,2,2,2,2,1,1))
dev.off()


#----------------------------------------------------------------------------------------------------------------
#   #SmE plots
#----------------------------------------------------------------------------------------------------------------
#RNAseq tracks
SmEnames <- c("SmE","SmE","WT","WT")
SmE.TrackList <- list()
for (i in seq_along(bamFilesR)){
  SmE.TrackList[[i]] <- AlignmentsTrack(
    bamFilesR[i],type=c("coverage", "sashimi"),
    isPaired = FALSE,chromosome =chr,name=SmEnames[i],cex.legend=80,
    col=plotcols[1],col.sashimi=plotcols[1],coverageHeight=0.2,fill.coverage=plotcols[1],
    sashimiScore=1,mapq=255,ylim=c(0,50))
}

pdf(sprintf("Figures/%s_Nrde2-hisat_RNAseq_sashimi_plot_fixed_ylim_with_SmE.pdf",gene),height=24,width=30)
plotTracks(c(peak.TrackList[c(3)],BCLIP.TrackList[c(3)],peak.TrackList[c(2)],BCLIP.TrackList[c(2)],Nrde2KO.TrackList[1:2],WT.TrackList[1:2],SmE.TrackList,txTr,axisTrack), from = afrom, to = ato, 
           chromosome = chr,reverseStrand = strand,sizes=c(0.5,1,0.5,1,2,2,2,2,2,2,2,2,1,1))
dev.off()
