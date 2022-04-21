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
f_counts <- featureCounts(bamFilesR,annot.ext="/work/gbioinfo/DB/GENCODE/Mouse/release_M23/gencode.vM23.annotation.gtf",isGTFAnnotationFile = TRUE,
                          useMetaFeatures=TRUE,allowMultiOverlap=FALSE,
                          minOverlap=1,countMultiMappingReads=FALSE,fraction=FALSE,
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
#map reads to transcripts using salmon on xenon6
cd /tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/RNAseq/salmon_transcript_expression
module load Salmon/1.2.0
salmon quant -i /work/gbioinfo/DB/GENCODE/Mouse/release_M23/refsalmon_spliced_k23/salmon-1.0.0_GRCm38.primary_assembly_gencodeM23_spliced_k23.sidx -l A -r /tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/flemrmat_121119_143452_2191F1-2191F22_/Nrde2-WT_1_2191F9/Nrde2-WT_1_2191F9.fastq.gz --validateMappings --gcBias -o Nrde2_WT_1_transcripts_quant
salmon quant -i /work/gbioinfo/DB/GENCODE/Mouse/release_M23/refsalmon_spliced_k23/salmon-1.0.0_GRCm38.primary_assembly_gencodeM23_spliced_k23.sidx -l A -r /tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/flemrmat_121119_143452_2191F1-2191F22_/Nrde2-WT_2_2191F10/Nrde2-WT_2_2191F10.fastq.gz --validateMappings --gcBias -o Nrde2_WT_2_transcripts_quant
salmon quant -i /work/gbioinfo/DB/GENCODE/Mouse/release_M23/refsalmon_spliced_k23/salmon-1.0.0_GRCm38.primary_assembly_gencodeM23_spliced_k23.sidx -l A -r /tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/flemrmat_121119_143452_2191F1-2191F22_/Mtr4_WT_1_2191F19/Mtr4_WT_1_2191F19.fastq.gz --validateMappings --gcBias -o Mtr4_WT_1_transcripts_quant
salmon quant -i /work/gbioinfo/DB/GENCODE/Mouse/release_M23/refsalmon_spliced_k23/salmon-1.0.0_GRCm38.primary_assembly_gencodeM23_spliced_k23.sidx -l A -r /tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/flemrmat_121119_143452_2191F1-2191F22_/Mtr4_WT_2_2191F20/Mtr4_WT_2_2191F20.fastq.gz --validateMappings --gcBias -o Mtr4_WT_2_transcripts_quant

salmon quant -i /work/gbioinfo/DB/GENCODE/Mouse/release_M23/refsalmon_spliced_k23/salmon-1.0.0_GRCm38.primary_assembly_gencodeM23_spliced_k23.sidx -l A -r /tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/michi_300620_125154_2447F1-2447F16_Nrde2RNA/fastq_files/Nrde2WT_1_2447F1.fastq.gz --validateMappings --gcBias -o Nrde2_WT_3_transcripts_quant
salmon quant -i /work/gbioinfo/DB/GENCODE/Mouse/release_M23/refsalmon_spliced_k23/salmon-1.0.0_GRCm38.primary_assembly_gencodeM23_spliced_k23.sidx -l A -r /tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/michi_300620_125154_2447F1-2447F16_Nrde2RNA/fastq_files/Nrde2WT_2_2447F2.fastq.gz --validateMappings --gcBias -o Nrde2_WT_4_transcripts_quant
salmon quant -i /work/gbioinfo/DB/GENCODE/Mouse/release_M23/refsalmon_spliced_k23/salmon-1.0.0_GRCm38.primary_assembly_gencodeM23_spliced_k23.sidx -l A -r /tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/michi_300620_125154_2447F1-2447F16_Nrde2RNA/fastq_files/Ccdc174WT_1_2447F9.fastq.gz --validateMappings --gcBias -o Ccdc174_WT_1_transcripts_quant
salmon quant -i /work/gbioinfo/DB/GENCODE/Mouse/release_M23/refsalmon_spliced_k23/salmon-1.0.0_GRCm38.primary_assembly_gencodeM23_spliced_k23.sidx -l A -r /tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/michi_300620_125154_2447F1-2447F16_Nrde2RNA/fastq_files/Ccdc174WT_2_2447F10.fastq.gz --validateMappings --gcBias -o Ccdc174_WT_2_transcripts_quant

salmon quant -i /work/gbioinfo/DB/GENCODE/Mouse/release_M23/refsalmon_spliced_k23/salmon-1.0.0_GRCm38.primary_assembly_gencodeM23_spliced_k23.sidx -l A -r \
/tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/flemrmat_121119_143452_2191F1-2191F22_/Nrde2-WT_1_2191F9/Nrde2-WT_1_2191F9.fastq.gz \
/tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/flemrmat_121119_143452_2191F1-2191F22_/Mtr4_WT_1_2191F19/Mtr4_WT_1_2191F19.fastq.gz \
/tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/michi_300620_125154_2447F1-2447F16_Nrde2RNA/fastq_files/Nrde2WT_1_2447F1.fastq.gz \
/tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/michi_300620_125154_2447F1-2447F16_Nrde2RNA/fastq_files/Ccdc174WT_1_2447F9.fastq.gz \
--validateMappings --gcBias -o all_WT_1_transcripts_quant

salmon quant -i /work/gbioinfo/DB/GENCODE/Mouse/release_M23/refsalmon_spliced_k23/salmon-1.0.0_GRCm38.primary_assembly_gencodeM23_spliced_k23.sidx -l A -r \
/tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/flemrmat_121119_143452_2191F1-2191F22_/Nrde2-WT_2_2191F10/Nrde2-WT_2_2191F10.fastq.gz \
/tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/flemrmat_121119_143452_2191F1-2191F22_/Mtr4_WT_2_2191F20/Mtr4_WT_2_2191F20.fastq.gz \
/tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/michi_300620_125154_2447F1-2447F16_Nrde2RNA/fastq_files/Nrde2WT_2_2447F2.fastq.gz \
/tungstenfs/scratch/gbuehler/deepSeqRepos/RNAseq/michi_300620_125154_2447F1-2447F16_Nrde2RNA/fastq_files/Ccdc174WT_2_2447F10.fastq.gz \
--validateMappings --gcBias -o all_WT_2_transcripts_quant



Nrde2_WT_1 <- read.table("RNAseq/salmon_transcript_expression/all_WT_1_transcripts_quant/quant.sf",header=TRUE)
Nrde2_WT_2 <- read.table("RNAseq/salmon_transcript_expression/all_WT_2_transcripts_quant/quant.sf",header=TRUE)

Nrde2_WT_1$logTPM <- log2(Nrde2_WT_1$TPM + 1)
Nrde2_WT_2$logTPM <- log2(Nrde2_WT_2$TPM + 1)

plot(density(Nrde2_WT_1$logTPM),xlim=c(0,1))
abline(v=0.1)

expr.transcripts <- Nrde2_WT_1$Name[Nrde2_WT_1$logTPM > 0.1 & Nrde2_WT_2$logTPM > 0]
nonexpr.transcripts <- Nrde2_WT_1$Name[Nrde2_WT_1$logTPM < 0.1 & Nrde2_WT_2$logTPM ==0]


#----------------------------------------------------------------------------------------------------------------
#   #find exon junctions in expressed transcripts
#----------------------------------------------------------------------------------------------------------------
#txdb23 <- makeTxDbFromGFF("/work/gbioinfo/DB/GENCODE/Mouse/release_M23/gencode.vM23.annotation.gtf","gtf","GENCODE","Mus musculus",10090)                                         
#saveDb(txdb23, file="/tungstenfs/scratch/gbuehler/michi/Annotations/GENCODE/Mouse/release_M23/gencode.vM23.annotation.txdb.sqlite")


txdb=loadDb("/tungstenfs/scratch/gbuehler/michi/Annotations/GENCODE/Mouse/release_M23/gencode.vM23.annotation.txdb.sqlite")

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


save(spliceDonors,spliceAcceptors,spliceDonors.ne,spliceAcceptors.ne,introns3g,introns3ng,exons3g,file="splice_sites_expressed_and_notexpressed_transcripts_log2TPMabove0_and_genes_FPKMabove0.RData")

#----------------------------------------------------------------------------------------------------------------
#   #add transcript biotypes
#----------------------------------------------------------------------------------------------------------------
#
print(load("splice_sites_expressed_and_notexpressed_transcripts_log2TPMabove0_and_genes_FPKMabove0.RData"))

#get biotype info
library(EnsDb.Mmusculus.v79)
edb <- EnsDb.Mmusculus.v79
txtypes <- genes(edb, columns=c("gene_name", "gene_biotype", "tx_biotype", "tx_id"))
head(txtypes)

#join introns and splice sites with annotations
tx_id <- matrix(unlist(strsplit(spliceDonors$tx_name,".",fixed=TRUE)),ncol=2,byrow = TRUE)[,1]

spliceDonors2 <- left_join(data.frame(tx_id=tx_id,GeneID=spliceDonors$GeneID),data.frame(txtypes),by="tx_id")
head(spliceDonors2)
dim(spliceDonors2)

spliceDonors$gene_name <- spliceDonors2$gene_name
spliceDonors$gene_biotype <- spliceDonors2$gene_biotype
spliceDonors$tx_biotype <- spliceDonors2$tx_biotype
table(spliceDonors$tx_biotype)

spliceAcceptors$gene_name <- spliceDonors2$gene_name
spliceAcceptors$gene_biotype <- spliceDonors2$gene_biotype
spliceAcceptors$tx_biotype <- spliceDonors2$tx_biotype

introns3g$gene_name <- spliceDonors2$gene_name
introns3g$gene_biotype <- spliceDonors2$gene_biotype
introns3g$tx_biotype <- spliceDonors2$tx_biotype

#join exons with annotations
tx_id <- matrix(unlist(strsplit(exons3g$tx_name,".",fixed=TRUE)),ncol=2,byrow = TRUE)[,1]

exons3g2 <- left_join(data.frame(tx_id=tx_id,GeneID=exons3g$GeneID),data.frame(txtypes),by="tx_id")
head(exons3g2)
dim(exons3g2)

exons3g$gene_name <- exons3g2$gene_name
exons3g$gene_biotype <- exons3g2$gene_biotype
exons3g$tx_biotype <- exons3g2$tx_biotype
table(exons3g$tx_biotype)

save(spliceDonors,spliceAcceptors,spliceDonors.ne,spliceAcceptors.ne,introns3g,introns3ng,exons3g,file="splice_sites_expressed_and_notexpressed_transcripts_log2TPMabove0_and_genes_FPKMabove0.RData")

