rm(list=ls())

library(tidyverse)
library(Rsubread)
library(GenomicFeatures)
library(GenomicRanges)
library(DESeq2)
library(limma)
library(biomaRt)
library(pheatmap)
library(ComplexHeatmap)


setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")

#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#   #calculate gene expression
#----------------------------------------------------------------------------------------------------------------

bamFilesR <- list.files("/tungstenfs/nobackup/gbuehler/michi/Projects/Nrde2/Riboseq/experiment_1/raw_data", full.names=TRUE,pattern="*LCstripped.prinseq_Aligned.sortedByCoord.out.bam$")
bamNamesR <- gsub("/tungstenfs/nobackup/gbuehler/michi/Projects/Nrde2/Riboseq/experiment_1/raw_data/riboseq_","",bamFilesR)
bamNamesR <- gsub(".LCstripped.prinseq_Aligned.sortedByCoord.out.bam","",bamNamesR)


#calculate the number of reads per gene 
f_counts <- featureCounts(bamFilesR,annot.ext="/work/gbioinfo/DB/GENCODE/Mouse/release_M23/gencode.vM23.annotation.gtf",isGTFAnnotationFile = TRUE,
                           useMetaFeatures=TRUE,allowMultiOverlap=FALSE,
                           minOverlap=5,countMultiMappingReads=FALSE,fraction=FALSE,
                           minMQS=255,strandSpecific=1,nthreads=20,verbose=FALSE,isPairedEnd=FALSE)


counts <- data.frame(f_counts$counts)
colnames(counts) <- bamNamesR
counts$GeneID <- matrix(unlist(strsplit(rownames(counts),".",fixed=TRUE)),ncol=2,byrow = TRUE)[,1]
rownames(counts) <- counts$GeneID
write.table(counts,"FigureData/feature_counts_per_gene_Riboseq.txt",sep="\t",col.names = TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

#calculate CPMs and FPKMs
fpkm <- f_counts$counts
cpm <- f_counts$counts

for (i in seq_along(bamFilesR)){
  scaling_factor <-  sum(f_counts$stat[f_counts$stat$Status=="Assigned"| f_counts$stat$Status=="Unassigned_NoFeatures",i+1])/1e6
  cpm[,i] <- f_counts$counts[,i]/scaling_factor
  fpkm[,i] <- cpm[,i]/(f_counts$annotation$Length/1000)
}

cpm <- data.frame(cpm)
colnames(cpm) <- bamNamesR
cpm$GeneID <- counts$GeneID
write.table(cpm,"FigureData/cpm_per_gene_Riboseq.txt",sep="\t",col.names = TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

fpkm <- data.frame(fpkm)
colnames(fpkm) <- bamNamesR
fpkm$GeneID <- counts$GeneID
write.table(fpkm,"FigureData/fpkm_per_gene_Riboseq.txt",sep="\t",col.names = TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

#----------------------------------------------------------------------------------------------------------------
#   #calculate differential gene expression
#----------------------------------------------------------------------------------------------------------------
#counts <- read.table("RNAseq/feature_counts_per_gene.txt",sep="\t",header=TRUE)

# prepare annotations
Group <- c(rep("Nrde2_KO",2),rep("Nrde2_WT",2))

annots <- data.frame(sample=colnames(counts[,!colnames(counts) %in% c("GeneID")]),Group=Group,
                     replicate=rep(c("rep1","rep2"),2))
rownames(annots) <- annots$sample

#extract species name
species <- "Mus musculus"


# check that everything is fine
#counts <- counts[,rownames(annots)] # reorder raw.counts after annots
identical(rownames(annots),colnames(counts[,!colnames(counts) %in% c("GeneID")]))

# prepare features
features <- f_counts$annotation[,c(1,6)]
features$GeneID <- matrix(unlist(strsplit(features$GeneID,".",fixed=TRUE)),ncol=2,byrow = TRUE)[,1]
rownames(features) <- features$GeneID

#________________________________________________________
#  Gene mapping and description with BioMart 
#________________________________________________________
# Load the organism-specific biomart
ensembl <- biomaRt::useEnsembl(
  biomart = 'ensembl', 
 # mirror="useast",
  dataset = paste0('mmusculus', '_gene_ensembl')
)
#
geneName <- biomaRt::getBM(attributes = c('ensembl_gene_id','mgi_symbol', 
                                          "entrezgene_id", "description","chromosome_name"), 
                           filters = 'ensembl_gene_id',
                           values = rownames(features), 
                           mart = ensembl)

geneName<-geneName[!duplicated(geneName$ensembl_gene_id),]
rownames(geneName)<-geneName$ensembl_gene_id
#geneName<-geneName[,-1]
description<-lapply(seq(length(geneName$description)),function(i){
  strsplit(geneName$description[i],"[[]")[[1]][1]
})
description<-(unlist(description))
geneName$description<-description
colnames(geneName) <- c("GeneID","gene_symbol", "entrezgene" ,"description" ,"chromosome_name")
#----------------------------------------------------------------------

features <- left_join(features,geneName[,c(1,2,3)],by="GeneID")
#save(features,file="RNAseq/features.RData")
#load("RNAseq/features.RData")
#rownames(features) <- features$GeneID



#generate DEseq dataset
dds <- DESeqDataSetFromMatrix(countData = counts[,!colnames(counts) %in% c("GeneID")],
                              colData = annots,
                              design = ~ Group)

# add feature annotation
mcols(dds) <- DataFrame(mcols(dds), features)

# filter out genes with less than 100 counts combined in all samples = 10 counts/sample on average
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]



# make sure the control group is the first level
#dds$Group <- factor(dds$Group, levels = c("Nrde2_WT","Nrde2_KO","Nrde2_D174R","Nrde2_d100","Nrde2_d200","Mtr4_WT","Mtr4_KO"))

# differential expression analysis
dds <- DESeq(dds)


#PCA
norm_counts <- counts(dds,normalized=TRUE)
pca <- prcomp(t(norm_counts))
PCA <- data.frame(pca$x)
PCA$sample <- rownames(PCA)
PCA <- inner_join(PCA,annots,by="sample")
ggplot(PCA,aes(x=PC1,y=PC2)) + geom_point(aes(color=Group))  + theme_classic() 
ggsave("Figures/Nrde2_RIBOseq_PCA_cutoff1.pdf",device="pdf",height= 5,width=7)




#get DE results
contrasts <- list(
                  Nrde2KO_vs_WT=c("Group","Nrde2_KO","Nrde2_WT")
                  )
res2 <- list()
for (i in seq_along(contrasts)){
res <- data.frame(results(dds, contrast=contrasts[[i]], independentFiltering=FALSE))
res$Contrast <- names(contrasts)[i]
res$GeneID <- as.character(rownames(res))
res2[[i]] <- res
}

res2 <- do.call("rbind",res2)
res2$log.padj <- -log10(res2$padj)
res2 <- left_join(res2,features,by="GeneID")
res2$regulated <- ifelse(res2$log2FoldChange > 0.5 & res2$padj < 0.05,"up",ifelse(res2$log2FoldChange < -0.5 & res2$padj < 0.05,"down","no"))
table(res2$regulated)
write.table(res2,file="FigureData/Nrde2_RIBOseq_DEseq2_results_v2_cutoff10.txt",sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

#Volcanoplot
vpa <- ggplot(res2,aes(x=log2FoldChange,y=log.padj)) +geom_point() + theme_classic() +facet_wrap(vars(Contrast))
vpa <- vpa + geom_vline(xintercept=c(-1,1),alpha=0.4) + geom_hline(yintercept = 2,alpha=0.4)
vpa
ggsave("Figures/Nrde2RIBOseq_volcano_plots_cutoff10.png",device="png",height=10,width=10)
