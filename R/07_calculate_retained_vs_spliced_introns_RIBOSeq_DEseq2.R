
rm(list=ls())

#install.packages("/tungstenfs/scratch/gbuehler/bioinfo/Rpackages/MiniChip_0.0.0.9000.tar.gz",repos=NULL)
library(tidyverse)
library(Rsubread)
library(GenomicFeatures)
library(GenomicRanges)
library(BSgenome)
library(DESeq2)

setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")


#----------------------------------------------------------------------------------------------------------------
#   #load RNAseq bam files and junctions
#----------------------------------------------------------------------------------------------------------------
#find RIBO bam files
bamFilesR <- list.files("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/Riboseq/experiment_1/raw_data", full.names=TRUE,pattern="*LCstripped.prinseq_Aligned.sortedByCoord.out.bam$")
bamNamesR <- gsub("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/Riboseq/experiment_1/raw_data/riboseq_","",bamFilesR)
bamNamesR <- gsub(".LCstripped.prinseq_Aligned.sortedByCoord.out.bam","",bamNamesR)


#load junctions
load("splice_sites_expressed_and_notexpressed_transcripts_log2TPMabove0_and_genes_FPKMabove0.RData")
names(spliceDonors) <- paste(names(spliceDonors),start(spliceDonors),sep="_")

#----------------------------------------------------------------------------------------------------------------
#   #generate intron RNAseq read counts
#-----------------------------------------------------------------------------------------------------------------

names(introns3g) <- names(spliceDonors)
#generate  saf format data frame
introns.saf <- data.frame(GeneID= names(introns3g), Chr=seqnames(introns3g),
                          Start=start(introns3g), End=end(introns3g),Strand=strand(introns3g))


#calculate the number of reads overlapping introns
counts.introns <- featureCounts(bamFilesR,annot.ext=introns.saf,
                                useMetaFeatures=FALSE,allowMultiOverlap=TRUE,
                                minOverlap=10,countMultiMappingReads=FALSE,fraction=FALSE,
                                minMQS=255,strandSpecific=1,nthreads=10,verbose=FALSE,isPairedEnd=FALSE)

counts.introns2 <- data.frame(counts.introns$counts)
colnames(counts.introns2) <- bamNamesR
#plot(density(counts.introns2$`Nrde2-WT_1_spliced_2191F9`))
#counts.introns2$ID <- row.names(counts.introns2)
#counts.introns2$intron_length <- width(introns3g)


#----------------------------------------------------------------------------------------------------------------
#   #perform differential exression analysis on introns
#-----------------------------------------------------------------------------------------------------------------


# prepare annotations
Group <- c(rep("Nrde2KO",2),rep("WT",2))
replicate <- rep(c("rep1","rep2"),2)
annots <- data.frame(sample=colnames(counts.introns2),Group=Group,replicate=replicate)
rownames(annots) <- annots$sample

# check that everything is fine
#counts <- counts[,rownames(annots)] # reorder raw.counts after annots
identical(rownames(annots),colnames(counts.introns2))


#generate DEseq dataset
dds <- DESeqDataSetFromMatrix(countData = counts.introns2,
                              colData = annots,
                              design = ~ Group)
# add feature annotation
#mcols(dds) <- DataFrame(mcols(dds), features)

# filter out genes with less than 5 counts per sample combined in all samples
keep <- rowSums(counts(dds)) >= 20
table(keep)
dds <- dds[keep,]


# differential expression analysis
dds <- DESeq(dds)

#PCA 
norm_counts <- counts(dds,normalized=TRUE)
pca <- prcomp(t(norm_counts))
PCA <- data.frame(pca$x)
PCA$sample <- rownames(PCA)
PCA <- inner_join(PCA,annots,by="sample")
ggplot(PCA,aes(x=PC1,y=PC2)) + geom_point(aes(color=Group))  + theme_classic() +
  scale_color_manual(values=c("darkblue","darkgreen"))
ggsave("Figures/intronic_RIBOseq_reads_Nrde2_PCA_5readspersample.pdf",device="pdf",height=5,width=7)



#get DE results
contrasts <- list(
                  Nrde2KO_vs_WT=c("Group","Nrde2KO","WT")
                 
)
res2 <- list()
for (i in seq_along(contrasts)){
  res <- data.frame(results(dds, contrast=contrasts[[i]], independentFiltering=FALSE))
  res$Contrast <- names(contrasts)[i]
  res$ID <- as.character(rownames(res))
  res2[[i]] <- res
}

res2 <- do.call("rbind",res2)
res2$log.padj <- -log10(res2$padj)
#res2 <- left_join(res2,intron.RNA.bclip2[,c(1:16,23)],by="ID")
res2$intron_regulated <- ifelse(res2$log2FoldChange > 0.5 & res2$padj < 0.05,"up",ifelse(res2$log2FoldChange < -0.5 & res2$padj < 0.05,"down","no"))
table(res2$intron_regulated)
write.table(res2,file="FigureData/Nrde2_intronicRIBOseq_DEseq2_results_0.05_0.5FC_5readspersample.txt",sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

#Volcanoplot
vpa <- ggplot(res2,aes(x=log2FoldChange,y=log.padj)) +geom_point() + theme_classic() +facet_wrap(vars(Contrast))
vpa <- vpa + geom_vline(xintercept=c(-1,1),alpha=0.4) + geom_hline(yintercept = 2,alpha=0.4)
vpa
ggsave("Figures/Intron_RIBOseq_Nrde2_volcano_plots_5readspersample.pdf",device="pdf",height=15,width=15)

