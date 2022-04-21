rm(list=ls())

library(tidyverse)
library(Rsubread)
library(GenomicFeatures)
library(GenomicRanges)

setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")

#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#   #find expressed genes
#----------------------------------------------------------------------------------------------------------------

all.bamFiles <- list.files("/tungstenfs/scratch/gbuehler/deepSeqRepos/bam/", full.names=TRUE,pattern="*bam$")
#bamFiles <- c(grep("2191F9",all.bamFiles,value=TRUE),grep("2191F10",all.bamFiles,value=TRUE),grep("2191F11",all.bamFiles,value=TRUE),grep("2191F12",all.bamFiles,value=TRUE))
bamFilesR <- c(grep("2191F",all.bamFiles,value=TRUE))
bamFilesR <- c(grep("spliced",bamFilesR,value=TRUE))
bamFilesR <- c(grep("Nxf1",bamFilesR,value=TRUE,invert=TRUE))
bamFilesR <- c(grep("WT",bamFilesR,value=TRUE))

bamFilesR2 <- grep("2447F",all.bamFiles,value=TRUE)
bamFilesR2 <- c(grep("WT",bamFilesR2,value=TRUE))
bamFilesR2 <- c(grep("CHX",bamFilesR2,value=TRUE,invert=TRUE))

bamFilesR <- c(bamFilesR,bamFilesR2)

bamNamesR <- gsub("/tungstenfs/scratch/gbuehler/deepSeqRepos/bam//","",bamFilesR)
bamNamesR <- gsub("_Aligned.sortedByCoord.out.bam","",bamNamesR)

#calculate the number of reads per gene 
f_counts <- featureCounts(bamFilesR,annot.ext="/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/RNAseq/stringtie_groupwise_dta/WT.annotated.gtf",isGTFAnnotationFile = TRUE,
                          useMetaFeatures=TRUE,allowMultiOverlap=FALSE,
                          minOverlap=5,countMultiMappingReads=FALSE,fraction=FALSE,
                          minMQS=255,strandSpecific=2,nthreads=20,verbose=FALSE,isPairedEnd=FALSE)


counts <- data.frame(f_counts$counts)
colnames(counts) <- bamNamesR
counts$GeneID <- matrix(unlist(strsplit(rownames(counts),".",fixed=TRUE)),ncol=2,byrow = TRUE)[,1]
rownames(counts) <- counts$GeneID

#calculate CPMs and FPKMs
fpkm <- f_counts$counts
cpm <- f_counts$counts

for (i in seq_along(bamFilesR)){
  scaling_factor <-  f_counts$stat[f_counts$stat$Status=="Assigned",i+1]/1e6
  cpm[,i] <- f_counts$counts[,i]/scaling_factor
  fpkm[,i] <- cpm[,i]/(f_counts$annotation$Length/1000)
}

cpm <- data.frame(cpm)
colnames(cpm) <- bamNamesR
cpm$GeneID <- counts$GeneID

fpkm <- data.frame(fpkm)
colnames(fpkm) <- bamNamesR
fpkm$GeneID <- counts$GeneID

expr.genes <- row.names(fpkm)[rowSums(fpkm[,1:8]) > 0]
non.expr.genes <- row.names(fpkm)[rowSums(fpkm[,1:8]) == 0]

plot(density(log2(fpkm[,1]+1)))
plot(density(fpkm[,2]))

#----------------------------------------------------------------------------------------------------------------
#   #find expressed transcripts
#----------------------------------------------------------------------------------------------------------------
cd /tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/RNAseq/stringtie_groupwise_dta/SalmonIndex

#make fasta
/work/gbioinfo/Appz/cufflinks/cufflinks-current/gffread ../WT.annotated.gtf \
-g /tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M23/GRCm38.primary_assembly.genome.fa \
-w WT.fa
#conatenate trasncriptome and genome sequecnes (for using genome as decoy)
cat WT.fa \
/tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M23/GRCm38.primary_assembly.genome.fa > \
GRCm38.primary_assembly.genome.gencode.WT.fa

module load Salmon/1.2.0

#generate salmon reference based on Stringtie gtf (using genome as decoy)
salmon index -t GRCm38.primary_assembly.genome.gencode.WT.fa \
-k 23 -i salmon-1.2.0_GRCm38.primary_assembly.genome.gencode.WT_spliced_k23_gentrome.sidx --type puff -p 16 \
-d /tungstenfs/groups/gbioinfo/DB/GENCODE/Mouse/release_M23/GRCm38.primary_assembly.chromosome_names.txt

#map reads to transcripts using salmon on xenon6
cd /tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/RNAseq/salmon_transcript_expression_WT_gtf

salmon quant -i /tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/RNAseq/stringtie_groupwise_dta/SalmonIndex/salmon-1.2.0_GRCm38.primary_assembly.genome.gencode.WT_spliced_k23_gentrome.sidx -l A -r \
/tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/flemrmat_121119_143452_2191F1-2191F22_/Nrde2-WT_1_2191F9/Nrde2-WT_1_2191F9.fastq.gz \
/tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/flemrmat_121119_143452_2191F1-2191F22_/Mtr4_WT_1_2191F19/Mtr4_WT_1_2191F19.fastq.gz \
/tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/michi_300620_125154_2447F1-2447F16_Nrde2RNA/fastq_files/Nrde2WT_1_2447F1.fastq.gz \
/tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/michi_300620_125154_2447F1-2447F16_Nrde2RNA/fastq_files/Ccdc174WT_1_2447F9.fastq.gz \
--validateMappings -o all_WT_1_transcripts_quant

salmon quant -i /tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/RNAseq/stringtie_groupwise_dta/SalmonIndex/salmon-1.2.0_GRCm38.primary_assembly.genome.gencode.WT_spliced_k23_gentrome.sidx -l A -r \
/tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/flemrmat_121119_143452_2191F1-2191F22_/Nrde2-WT_2_2191F10/Nrde2-WT_2_2191F10.fastq.gz \
/tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/flemrmat_121119_143452_2191F1-2191F22_/Mtr4_WT_2_2191F20/Mtr4_WT_2_2191F20.fastq.gz \
/tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/michi_300620_125154_2447F1-2447F16_Nrde2RNA/fastq_files/Nrde2WT_2_2447F2.fastq.gz \
/tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/michi_300620_125154_2447F1-2447F16_Nrde2RNA/fastq_files/Ccdc174WT_2_2447F10.fastq.gz \
--validateMappings -o all_WT_2_transcripts_quant



Nrde2_WT_1 <- read.table("RNAseq/salmon_transcript_expression_WT_gtf/all_WT_1_transcripts_quant/quant.sf",header=TRUE)
Nrde2_WT_2 <- read.table("RNAseq/salmon_transcript_expression_WT_gtf/all_WT_2_transcripts_quant/quant.sf",header=TRUE)

Nrde2_WT_1$logTPM <- log2(Nrde2_WT_1$TPM + 1)
Nrde2_WT_2$logTPM <- log2(Nrde2_WT_2$TPM + 1)

plot(density(Nrde2_WT_1$logTPM),xlim=c(0,1))
abline(v=0.1)

expr.transcripts <- Nrde2_WT_1$Name[Nrde2_WT_1$logTPM > 0.1 & Nrde2_WT_2$logTPM > 0.1]
nonexpr.transcripts <- Nrde2_WT_1$Name[Nrde2_WT_1$logTPM < 0.1 & Nrde2_WT_2$logTPM < 0.1]


#----------------------------------------------------------------------------------------------------------------
#   #find exon junctions in expressed transcripts
#----------------------------------------------------------------------------------------------------------------
/tungstenfs/scratch/gbuehler/michi/Tools/gffcompare-0.11.6.Linux_x86_64/gffcompare -Q -M -r /work/gbioinfo/DB/GENCODE/Mouse/release_M23/gencode.vM23.annotation.gtf -G -o RNAseq/stringtie_groupwise_dta/WT RNAseq/stringtie_groupwise_dta/WT.gtf

txdb23 <- makeTxDbFromGFF("RNAseq/stringtie_groupwise_dta/WT.annotated.gtf","gtf")                                         
saveDb(txdb23, file="RNAseq/stringtie_groupwise_dta/WT.annotated.txdb.sqlite")


txdb=loadDb("RNAseq/stringtie_groupwise_dta/WT.annotated.txdb.sqlite")

exons <- exonsBy(txdb,by="tx",use.names=TRUE)
exons2 <- exons[names(exons) %in% expr.transcripts]
exons3 <- unique(unlist(exons2))

introns <- intronsByTranscript(txdb,use.names=TRUE)

#split into expressed and not expressed transcripts
introns2 <- introns[names(introns) %in% expr.transcripts]
introns2n <- introns[names(introns) %in% nonexpr.transcripts]

introns3 <- unique(unlist(introns2))
introns3n <- unique(unlist(introns2n))

#add GeneiDs to introns3 and exons3 GRanges objects
genes <- unique(unlist(transcriptsBy(txdb,by="gene")))
genes$GeneID <- names(genes)
genes.df <- data.frame(genes)

introns3$tx_name <- names(introns3)
introns3.df <- data.frame(introns3)
introns3.df2 <- left_join(introns3.df,genes.df,by=c("tx_name"="tx_name"))
introns3$GeneID <- introns3.df2$GeneID

introns3n$tx_name <- names(introns3n)
introns3n.df <- data.frame(introns3n)
introns3n.df2 <- left_join(introns3n.df,genes.df,by=c("tx_name"="tx_name"))
introns3n$GeneID <- introns3n.df2$GeneID

exons3$tx_name <- names(exons3)
exons3.df <- data.frame(exons3)
exons3.df2 <- left_join(exons3.df,genes.df,by=c("tx_name"="tx_name"))
exons3$GeneID <- exons3.df2$GeneID

#split into expressed and not expressed genes
introns3g <- introns3[introns3$GeneID %in% expr.genes]
introns3ng <- introns3n[introns3n$GeneID %in% non.expr.genes]
exons3g <- exons3[exons3$GeneID %in% expr.genes]

spliceDonors <- resize(introns3g, width = 1L, fix = "start")
spliceAcceptors <- resize(introns3g, width = 1L, fix = "end")

spliceDonors.ne <- resize(introns3ng, width = 1L, fix = "start")
spliceAcceptors.ne <- resize(introns3ng, width = 1L, fix = "end")


save(spliceDonors,spliceAcceptors,spliceDonors.ne,spliceAcceptors.ne,introns3g,introns3ng,exons3g,file="WT_strtie_gtf_splice_sites_expressed_and_notexpressed_transcripts_log2TPMabove0.1_and_genes_FPKMabove8below0.RData")
