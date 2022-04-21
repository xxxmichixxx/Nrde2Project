
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
#find RNAseq bam files
all.bamFiles <- list.files("/work2/gbuehler/deepSeqRepos/bam/", full.names=TRUE,pattern="*bam$")
#bamFiles <- c(grep("2191F9",all.bamFiles,value=TRUE),grep("2191F10",all.bamFiles,value=TRUE),grep("2191F11",all.bamFiles,value=TRUE),grep("2191F12",all.bamFiles,value=TRUE))
bamFilesR <- c(grep("2191F",all.bamFiles,value=TRUE))
bamFilesR <- c(grep("spliced",bamFilesR,value=TRUE))
bamFilesR <- c(grep("Nxf1",bamFilesR,value=TRUE,invert=TRUE))
bamFilesR2 <- grep("2447F",all.bamFiles,value=TRUE)
bamFilesR <- c(bamFilesR,bamFilesR2)

bamNamesR <- gsub("/work2/gbuehler/deepSeqRepos/bam//","",bamFilesR)
bamNamesR <- gsub("_Aligned.sortedByCoord.out.bam","",bamNamesR)

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
                                minMQS=255,strandSpecific=2,nthreads=10,verbose=FALSE,isPairedEnd=FALSE)

counts.introns2 <- data.frame(counts.introns$counts)
colnames(counts.introns2) <- bamNamesR
#plot(density(counts.introns2$`Nrde2-WT_1_spliced_2191F9`))
#counts.introns2$ID <- row.names(counts.introns2)
#counts.introns2$intron_length <- width(introns3g)


#----------------------------------------------------------------------------------------------------------------
#   #perform differential exression analysis on introns
#-----------------------------------------------------------------------------------------------------------------


# prepare annotations
Group <- c(rep("Mtr4_KO",2),rep("Mtr4_WT",2),rep("Nrde2_d100",2),rep("Nrde2_D174R",2),
           rep("Nrde2_d200",2),rep("Nrde2_KO",2),rep("Nrde2_WT",2),
           rep("Ccdc174_KD",2),rep("Ccdc174_KD_CHX",2), rep("Ccdc174_WT",2),rep("Ccdc174_WT_CHX",2),
           rep("Nrde2_KO",2),rep("Nrde2_KO_CHX",2), rep("Nrde2_WT",2),rep("Nrde2_WT_CHX",2))
replicate <- rep(c("rep1","rep2"),15)
batch <- c(rep("batch1",14),rep("batch2",16))
annots <- data.frame(sample=colnames(counts.introns2),Group=Group,replicate=replicate,batch=batch)
rownames(annots) <- annots$sample

# check that everything is fine
#counts <- counts[,rownames(annots)] # reorder raw.counts after annots
identical(rownames(annots),colnames(counts.introns2))

#split data into CHX and non CHX
inx.CHX <- grep(pattern = "CHX",colnames(counts.introns2))
counts.introns <- counts.introns2[,-inx.CHX]
counts.introns.CHX <- counts.introns2[,inx.CHX]
annots1 <- annots[-inx.CHX,]
annots.CHX <- annots[inx.CHX,]

#keep only WT samples from counts and annots 
inx.WT <- grepl(pattern = "WT",colnames(counts.introns2))
inx.WT2 <- grepl(pattern = "2447F",colnames(counts.introns2))

counts.introns.WT <- counts.introns2[,inx.WT & inx.WT2]
annots.WT <- annots[inx.WT & inx.WT2,]
annots.WT$Group <- ifelse(grepl("CHX",annots.WT$Group),"WTCHX","WT")

#generate DEseq dataset
dds <- DESeqDataSetFromMatrix(countData = counts.introns,
                              colData = annots1,
                              design = ~ Group)
# add feature annotation
#mcols(dds) <- DataFrame(mcols(dds), features)

# filter out genes with less than x counts combined in all samples (x = 5 * number of samples)
keep <- rowSums(counts(dds)) >= 110
table(keep)
dds <- dds[keep,]


# differential expression analysis
dds <- DESeq(dds)

#PCA non CHX
norm_counts <- counts(dds,normalized=TRUE)
pca <- prcomp(t(norm_counts))
PCA <- data.frame(pca$x)
PCA$sample <- rownames(PCA)
PCA <- inner_join(PCA,annots,by="sample")
ggplot(PCA,aes(x=PC1,y=PC2)) + geom_point(aes(color=Group,shape=batch))  + theme_classic() +
  scale_color_manual(values=c("darkgreen","darkblue","red","orange","darkred","purple","black"))
ggsave("Figures/intronic_reads_Nrde2_batch1and2_PCA_5persamplereads.pdf",device="pdf",height=5,width=7)

#generate DEseq dataset CHX
ddsCHX <- DESeqDataSetFromMatrix(countData = counts.introns.CHX,
                              colData = annots.CHX,
                              design = ~ Group)
# add feature annotation
#mcols(dds) <- DataFrame(mcols(dds), features)

# filter out genes with less than x counts combined in all samples (x = 5 * number of samples)
keep <- rowSums(counts(ddsCHX)) >= 40
table(keep)
ddsCHX <- ddsCHX[keep,]


# differential expression analysis
ddsCHX <- DESeq(ddsCHX)

#PCA  CHX
norm_counts <- counts(ddsCHX,normalized=TRUE)
pca <- prcomp(t(norm_counts))
PCA <- data.frame(pca$x)
PCA$sample <- rownames(PCA)
PCA <- inner_join(PCA,annots,by="sample")
ggplot(PCA,aes(x=PC1,y=PC2)) + geom_point(aes(color=Group,shape=batch))  + theme_classic() #+
 # scale_color_manual(values=c("green","pink","darkgrey"))
ggsave("Figures/intronic_reads_Nrde2_batch1and2_PCA_CHX_5readspersample.pdf",device="pdf",height=5,width=7)



#generate DEseq dataset WT
ddsWT <- DESeqDataSetFromMatrix(countData = counts.introns.WT,
                                colData = annots.WT,
                                design = ~ Group)
# add feature annotation
#mcols(dds) <- DataFrame(mcols(dds), features)

# filter out genes with less than x counts combined in all samples (x = 5 * number of samples)
keep <- rowSums(counts(ddsWT)) >= 40
table(keep)
ddsWT <- ddsWT[keep,]
# differential expression analysis
ddsWT <- DESeq(ddsWT)



#get DE results
contrasts <- list(Mtr4KO_vs_WT=c("Group","Mtr4_KO","Mtr4_WT"),
                  Nrde2KO_vs_WT=c("Group","Nrde2_KO","Nrde2_WT"),
                  Nrde2D174R_vs_WT=c("Group","Nrde2_D174R","Nrde2_WT"),
                  Nrde2d100_vs_WT=c("Group","Nrde2_d100","Nrde2_WT"),
                  Nrde2d200_vs_WT=c("Group","Nrde2_d200","Nrde2_WT"),
                  Ccdc174KD_vs_WT=c("Group","Ccdc174_KD","Ccdc174_WT")
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
write.table(res2,file="FigureData/Nrde2_Mtr4_Ccdc174_batch1and2_intronicRNAseq_DEseq2_results_0.05_0.5FC_5readspersample_seperateWTs.txt",sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

#get DE results CHX
contrasts <- list(
                  Nrde2KOCHX_vs_WTCHX=c("Group","Nrde2_KO_CHX","Nrde2_WT_CHX"),
                  Ccdc174KDCHX_vs_WTCHX=c("Group","Ccdc174_KD_CHX","Ccdc174_WT_CHX")
)
res2.CHX <- list()
for (i in seq_along(contrasts)){
  res <- data.frame(results(ddsCHX, contrast=contrasts[[i]], independentFiltering=FALSE))
  res$Contrast <- names(contrasts)[i]
  res$ID <- as.character(rownames(res))
  res2.CHX[[i]] <- res
}


res2.CHX <- do.call("rbind",res2.CHX)
res2.CHX$log.padj <- -log10(res2.CHX$padj)
#res2 <- left_join(res2,intron.RNA.bclip2[,c(1:16,23)],by="ID")
res2.CHX$intron_regulated <- ifelse(res2.CHX$log2FoldChange > 0.5 & res2.CHX$padj < 0.05,"up",ifelse(res2.CHX$log2FoldChange < -0.5 & res2.CHX$padj < 0.05,"down","no"))
table(res2.CHX$intron_regulated)
write.table(res2.CHX,file="FigureData/CHX_Nrde2_Ccdc174_intronicRNAseq_DEseq2_results_0.05_0.5FC_5readspersample_seperateWTs.txt",sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)


#get DE results WT
contrasts <- list(WTCHX_vs_WT=c("Group","WTCHX","WT"))
                  
res2.WT <- list()
for (i in seq_along(contrasts)){
  res <- data.frame(results(ddsWT, contrast=contrasts[[i]], independentFiltering=FALSE))
  res$Contrast <- names(contrasts)[i]
  res$ID <- as.character(rownames(res))
  res2.WT[[i]] <- res
}


res2.WT <- do.call("rbind",res2.WT)
res2.WT$log.padj <- -log10(res2.WT$padj)
#res2 <- left_join(res2,intron.RNA.bclip2[,c(1:16,23)],by="ID")
res2.WT$intron_regulated <- ifelse(res2.WT$log2FoldChange > 0.5 & res2.WT$padj < 0.05,"up",ifelse(res2.WT$log2FoldChange < -0.5 & res2.WT$padj < 0.05,"down","no"))
table(res2.WT$intron_regulated)
write.table(res2.WT,file="FigureData/CHX_WT_intronicRNAseq_DEseq2_results_0.05_0.5FC_5readspersample.txt",sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)




#Volcanoplot
vpa <- ggplot(res2,aes(x=log2FoldChange,y=log.padj)) +geom_point() + theme_classic() +facet_wrap(vars(Contrast))
vpa <- vpa + geom_vline(xintercept=c(-1,1),alpha=0.4) + geom_hline(yintercept = 2,alpha=0.4)
vpa
ggsave("Figures/Intron_RNAseq_Nrde2_batch1and2_volcano_plots_5readspersample_seperateWTs.pdf",device="pdf",height=15,width=15)

#Volcanoplot CHX
vpa <- ggplot(res2.CHX,aes(x=log2FoldChange,y=log.padj)) +geom_point() + theme_classic() +facet_wrap(vars(Contrast))
vpa <- vpa + geom_vline(xintercept=c(-1,1),alpha=0.4) + geom_hline(yintercept = 2,alpha=0.4)
vpa
ggsave("Figures/CHX_Intron_RNAseq_Nrde2_batch1and2_volcano_plots_5readspersample_seperateWTs.pdf",device="pdf",height=15,width=15)



