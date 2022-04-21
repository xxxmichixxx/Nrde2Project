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


setwd("/tungstenfs/nobackup/gbuehler/michi/Projects/Nrde2/")

#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#   #calculate gene expression
#----------------------------------------------------------------------------------------------------------------

#find RNAseq bam files
all.bamFiles <- list.files("/tungstenfs/scratch/gbuehler/deepSeqRepos/bam/", full.names=TRUE,pattern="*bam$")
#bamFiles <- c(grep("2191F9",all.bamFiles,value=TRUE),grep("2191F10",all.bamFiles,value=TRUE),grep("2191F11",all.bamFiles,value=TRUE),grep("2191F12",all.bamFiles,value=TRUE))
bamFilesR <- c(grep("2191F",all.bamFiles,value=TRUE))
bamFilesR <- c(grep("spliced",bamFilesR,value=TRUE))
bamFilesR <- c(grep("Nxf1",bamFilesR,value=TRUE,invert=TRUE))
bamFilesR <- c(grep("d100",bamFilesR,value=TRUE,invert=TRUE))

bamFilesR2 <- grep("2447F",all.bamFiles,value=TRUE)
bamFilesR <- c(bamFilesR,bamFilesR2)

bamFilesUT <- c(grep("1441F",all.bamFiles,value=TRUE))
bamFilesR <- c(bamFilesR,bamFilesUT)

bamNamesR <- gsub("/tungstenfs/scratch/gbuehler/deepSeqRepos/bam//","",bamFilesR)
bamNamesR <- gsub("_Aligned.sortedByCoord.out.bam","",bamNamesR)

#calculate the number of reads per gene 
f_counts <- featureCounts(bamFilesR,annot.ext="/work/gbioinfo/DB/GENCODE/Mouse/release_M23/gencode.vM23.annotation.gtf",isGTFAnnotationFile = TRUE,
                           useMetaFeatures=TRUE,allowMultiOverlap=FALSE,
                           minOverlap=5,countMultiMappingReads=FALSE,fraction=FALSE,
                           minMQS=255,strandSpecific=2,nthreads=20,verbose=FALSE,isPairedEnd=FALSE)


counts <- data.frame(f_counts$counts)
colnames(counts) <- bamNamesR
counts$GeneID <- matrix(unlist(strsplit(rownames(counts),".",fixed=TRUE)),ncol=2,byrow = TRUE)[,1]
rownames(counts) <- counts$GeneID
write.table(counts,"FigureData/feature_counts_per_gene.txt",sep="\t",col.names = TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

#calculate CPMs and FPKMs
fpkm <- f_counts$counts
cpm <- f_counts$counts

for (i in seq_along(bamFilesR)){
  scaling_factor <-  sum(f_counts$stat[f_counts$stat$Status=="Assigned" | f_counts$stat$Status=="Unassigned_NoFeatures",i+1])/1e6
  cpm[,i] <- f_counts$counts[,i]/scaling_factor
  fpkm[,i] <- cpm[,i]/(f_counts$annotation$Length/1000)
}

cpm <- data.frame(cpm)
colnames(cpm) <- bamNamesR
cpm$GeneID <- counts$GeneID
write.table(cpm,"FigureData/cpm_per_gene.txt",sep="\t",col.names = TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

fpkm <- data.frame(fpkm)
colnames(fpkm) <- bamNamesR
fpkm$GeneID <- counts$GeneID
write.table(fpkm,"FigureData/fpkm_per_gene.txt",sep="\t",col.names = TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

#----------------------------------------------------------------------------------------------------------------
#   #calculate differential gene expression
#----------------------------------------------------------------------------------------------------------------
#counts <- read.table("RNAseq/feature_counts_per_gene.txt",sep="\t",header=TRUE)

# prepare annotations
Group <- c(rep("Mtr4_KO",2),rep("Mtr4_WT",2),rep("Nrde2_d100",2),rep("Nrde2_D174R",2),
           rep("Nrde2_d200",2),rep("Nrde2_KO",2),rep("Nrde2_WT",2),
           rep("Ccdc174_KD",2),rep("Ccdc174_KD_CHX",2), rep("Ccdc174_WT",2),rep("Ccdc174_WT_CHX",2),
           rep("Nrde2_KO",2),rep("Nrde2_KO_CHX",2), rep("Nrde2_WT",2),rep("Nrde2_WT_CHX",2))

annots <- data.frame(sample=colnames(counts[,!colnames(counts) %in% c("GeneID")]),Group=Group,
                     replicate=rep(c("rep1","rep2"),15),batch=c(rep("batch1",14),rep("batch2",16)))
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

#####split into CHX and non CHX data####
#generate index for CHX samples
inx.CHX <- grep(pattern = "CHX",colnames(counts))
#remove CHX samples from counts and annots 
counts1 <- counts[,-inx.CHX]
annots1 <- annots[-inx.CHX,]
#also remove d100 samples
inx.100 <- grep(pattern = "d100",colnames(counts1))
counts1 <- counts1[,-inx.100]
annots1 <- annots1[-inx.100,]

#keep only CHX samples from counts and annots 
counts.CHX <- counts[,inx.CHX]
annots.CHX <- annots[inx.CHX,]

#keep only WT samples from counts and annots 
inx.WT <- grepl(pattern = "WT",colnames(counts))
inx.WT2 <- grepl(pattern = "2447F",colnames(counts))

counts.WT <- counts[,inx.WT & inx.WT2]
annots.WT <- annots[inx.WT & inx.WT2,]
annots.WT$Group <- ifelse(grepl("CHX",annots.WT$Group),"WTCHX","WT")

#generate DEseq dataset
dds <- DESeqDataSetFromMatrix(countData = counts1[,!colnames(counts1) %in% c("GeneID")],
                              colData = annots1,
                              design = ~ Group)
#generate DEseq dataset fro CHX
ddsCHX <- DESeqDataSetFromMatrix(countData = counts.CHX[,!colnames(counts.CHX) %in% c("GeneID")],
                              colData = annots.CHX,
                              design = ~ Group)
#generate DEseq dataset fro WT
ddsWT <- DESeqDataSetFromMatrix(countData = counts.WT[,!colnames(counts.WT) %in% c("GeneID")],
                                 colData = annots.WT,
                                 design = ~ Group)
# add feature annotation
mcols(dds) <- DataFrame(mcols(dds), features)
mcols(ddsCHX) <- DataFrame(mcols(ddsCHX), features)
mcols(ddsWT) <- DataFrame(mcols(ddsWT), features)

# filter out genes with less than 100 counts combined in all samples = 10 counts/sample on average
keep <- rowSums(counts(dds)) >= 200
table(keep)
dds <- dds[keep,]

# CHX: filter out genes with less than 40 counts combined in all samples =5counts/sample on average
keep <- rowSums(counts(ddsCHX)) >= 80
table(keep)
ddsCHX <- ddsCHX[keep,]

# WT: filter out genes with less than 80 counts combined in all samples =10counts/sample on average
keep <- rowSums(counts(ddsWT)) >= 80
table(keep)
ddsWT <- ddsWT[keep,]

# make sure the control group is the first level
#dds$Group <- factor(dds$Group, levels = c("Nrde2_WT","Nrde2_KO","Nrde2_D174R","Nrde2_d100","Nrde2_d200","Mtr4_WT","Mtr4_KO"))

# differential expression analysis
dds <- DESeq(dds)
ddsCHX <- DESeq(ddsCHX)
ddsWT <- DESeq(ddsWT)


#PCA
norm_counts <- counts(dds,normalized=TRUE)
pca <- prcomp(t(norm_counts))
PCA <- data.frame(pca$x)
PCA$sample <- rownames(PCA)
PCA <- inner_join(PCA,annots,by="sample")
ggplot(PCA,aes(x=PC1,y=PC2)) + geom_point(aes(color=Group,shape=batch))  + theme_classic() +
scale_color_manual(values=c("orange","black","red","darkgrey","pink","violet","purple","grey"))
#+ geom_text(aes(label=sample),hjust=0, vjust=0)
ggsave("Figures/Nrde2_Mtr4_Ccdc174_RNAseq_PCA_cutoff200.pdf",device="pdf",height= 5,width=7)

#PCA CHX
norm_counts <- counts(ddsCHX,normalized=TRUE)
pca <- prcomp(t(norm_counts))
PCA <- data.frame(pca$x)
PCA$sample <- rownames(PCA)
PCA <- inner_join(PCA,annots,by="sample")
ggplot(PCA,aes(x=PC1,y=PC2)) + geom_point(aes(color=Group,shape=batch))  + theme_classic() +
  scale_color_manual(values=c("orange","black","purple","grey"))
#+ geom_text(aes(label=sample),hjust=0, vjust=0)
ggsave("Figures/Nrde2_Mtr4_Ccdc174_RNAseq_PCA_CHX_cutoff80.pdf",device="pdf",height= 5,width=7)


#get DE results
contrasts <- list(Mtr4KO_vs_WT=c("Group","Mtr4_KO","Mtr4_WT"),
                  Nrde2KO_vs_WT=c("Group","Nrde2_KO","Nrde2_WT"),
                  Nrde2D174R_vs_WT=c("Group","Nrde2_D174R","Nrde2_WT"),
                #  Nrde2d100_vs_WT=c("Group","Nrde2_d100","Nrde2_WT"),
                  Nrde2d200_vs_WT=c("Group","Nrde2_d200","Nrde2_WT"),
                  Ccdc174KD_vs_WT=c("Group","Ccdc174_KD","Ccdc174_WT")
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
write.table(res2,file="FigureData/Nrde2_Mtr4_Ccdc174_RNAseq_DEseq2_results_v2_cutoff200.txt",sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
res2 <- read.table("FigureData/Nrde2_Mtr4_Ccdc174_RNAseq_DEseq2_results_v2_cutoff200.txt",sep="\t",header=TRUE)


#get DE results WT
contrasts <- list(WTCHX_vs_WT=c("Group","WTCHX","WT")
)
res2.WT <- list()
for (i in seq_along(contrasts)){
  res <- data.frame(results(ddsWT, contrast=contrasts[[i]], independentFiltering=FALSE))
  res$Contrast <- names(contrasts)[i]
  res$GeneID <- as.character(rownames(res))
  res2.WT[[i]] <- res
}

res2.WT <- do.call("rbind",res2.WT)
res2.WT$log.padj <- -log10(res2.WT$padj)
res2.WT <- left_join(res2.WT,features,by="GeneID")
res2.WT$regulated <- ifelse(res2.WT$log2FoldChange > 0.5 & res2.WT$padj < 0.05,"up",ifelse(res2.WT$log2FoldChange < -0.5 & res2.WT$padj < 0.05,"down","no"))
table(res2.WT$regulated)
write.table(res2.WT,file="FigureData/CHX_WT_RNAseq_DEseq2_results_cutoff80.txt",sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

#get DE results CHX
contrasts <- list(Nrde2KOCHX_vs_WTCHX=c("Group","Nrde2_KO_CHX","Nrde2_WT_CHX"),
                  Ccdc174KDCHX_vs_WTCHX=c("Group","Ccdc174_KD_CHX","Ccdc174_WT_CHX")
)
res2.CHX <- list()
for (i in seq_along(contrasts)){
  res <- data.frame(results(ddsCHX, contrast=contrasts[[i]], independentFiltering=FALSE))
  res$Contrast <- names(contrasts)[i]
  res$GeneID <- as.character(rownames(res))
  res2.CHX[[i]] <- res
}

res2.CHX <- do.call("rbind",res2.CHX)
res2.CHX$log.padj <- -log10(res2.CHX$padj)
res2.CHX <- left_join(res2.CHX,features,by="GeneID")
res2.CHX$regulated <- ifelse(res2.CHX$log2FoldChange > 0.5 & res2.CHX$padj < 0.05,"up",ifelse(res2.CHX$log2FoldChange < -0.5 & res2.CHX$padj < 0.05,"down","no"))
table(res2.CHX$regulated)
write.table(res2.CHX,file="FigureData/CHX_Nrde2_Ccdc174_RNAseq_DEseq2_results_cutoff80.txt",sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)


#Volcanoplot
vpa <- ggplot(res2,aes(x=log2FoldChange,y=log.padj)) +geom_point() + theme_classic() +facet_wrap(vars(Contrast))
vpa <- vpa + geom_vline(xintercept=c(-1,1),alpha=0.4) + geom_hline(yintercept = 2,alpha=0.4)
vpa
ggsave("Figures/Nrde2_Mtr4_Ccdc174_RNAseq_volcano_plots_cutoff200.png",device="png",height=10,width=10)

#Volcanoplot CHX
vpa <- ggplot(res2.CHX,aes(x=log2FoldChange,y=log.padj)) +geom_point() + theme_classic() +facet_wrap(vars(Contrast))
vpa <- vpa + geom_vline(xintercept=c(-1,1),alpha=0.4) + geom_hline(yintercept = 2,alpha=0.4)
vpa
ggsave("Figures/CHX_Nrde2_Ccdc174_RNAseq_volcano_plots_cutoff80.png",device="png",height=10,width=10)



#select regulated genes
res.reg <- dplyr::filter(res2,padj < 0.05, abs(log2FoldChange) > 1)

#check out the counts data as a heatmap 
#ntd <- normTransform(dds)
rlog_counts <- rlog(dds,blind=FALSE)

#reg.ntd <- assay(ntd)[rownames(assay(ntd)) %in% res.reg$GeneID,]
reg.rlog_counts <- assay(rlog_counts)[rownames(assay(rlog_counts)) %in% res.reg$GeneID,]

#select <- order(rowMeans(counts(dds,normalized=TRUE)),
#                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Group","replicate")])
#png(file="RNAseq/heatmap_rlogtransformed_counts.png")
pheatmap(reg.rlog_counts, cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)
#pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#        cluster_cols=FALSE, annotation_col=df)
#dev.off()

#PCA of rlog regulated counts
pca <- prcomp(t(reg.rlog_counts))
PCA <- data.frame(pca$x)
PCA$sample <- rownames(PCA)
PCA <- inner_join(PCA,annots,by="sample")
ggplot(PCA,aes(x=PC1,y=PC2)) + geom_point(aes(color=Group,shape=batch))  + theme_classic() +
  scale_color_manual(values=c("darkgreen","darkblue","red","orange","darkred","purple","black"))
#+ geom_text(aes(label=sample),hjust=0, vjust=0)
ggsave("Figures/Nrd2_Mtr4_Ccdc174_PCA_regulated_genes_0.05_1FC.pdf",device="pdf",height=5,width=7)


#heatmap of fold changes
res.reg1 <- dplyr::filter(res2,padj < 0.01, abs(log2FoldChange)>1 & Contrast== "Nrde2KO_vs_WT")
res.reg2 <- dplyr::filter(res2,padj < 0.01, abs(log2FoldChange)>1 & Contrast== "Mtr4KO_vs_WT")
res.reg3 <- dplyr::filter(res2,padj < 0.01, abs(log2FoldChange)>1 & Contrast== "Nrde2D174R_vs_WT")
res.reg4 <- dplyr::filter(res2,padj < 0.01, abs(log2FoldChange)>1 & Contrast== "Nrde2d200_vs_WT")
res.reg5 <- dplyr::filter(res2,padj < 0.01, abs(log2FoldChange)>1 & Contrast== "Ccdc174KD_vs_WT")
regulated.genes <- unique(c(res.reg1$GeneID,res.reg2$GeneID,res.reg3$GeneID,res.reg4$GeneID,res.reg5$GeneID))

res3 <- pivot_wider(res2[,c("log2FoldChange","Contrast","GeneID")],names_from = Contrast,values_from = log2FoldChange)
res4 <- res3[res3$GeneID %in% regulated.genes,c(3,5,4,6,2)]
res5 <- res3[res3$GeneID %in% unique(c(res.reg1$GeneID)),c(3,5,4,6,2)]
res6 <- res3[res3$GeneID %in% unique(c(res.reg2$GeneID)),c(3,5,4,6,2)]

#res5 <- res5[order(res5$Nrde2KO_vs_WT),]

#group the genes into Mtr4 and Nrde2 regulated
Nrde2.up <- dplyr::filter(res2,padj < 0.01, log2FoldChange > 1 & Contrast== "Nrde2KO_vs_WT")
Nrde2.down <- dplyr::filter(res2,padj < 0.01, log2FoldChange < -1 & Contrast== "Nrde2KO_vs_WT")
Mtr4.up <- dplyr::filter(res2,padj < 0.01, log2FoldChange > 1 & Contrast== "Mtr4KO_vs_WT")
Mtr4.down <- dplyr::filter(res2,padj < 0.01, log2FoldChange < -1 & Contrast== "Mtr4KO_vs_WT")

NMreg <- ifelse(res3$GeneID %in% Nrde2.up$GeneID & res3$GeneID %in% Mtr4.up$GeneID, "Nrde2&Mtr4_up",
                ifelse(res3$GeneID %in% Nrde2.down$GeneID & res3$GeneID %in% Mtr4.down$GeneID, "Nrde2&Mtr4_down", 
                       ifelse(res3$GeneID %in% Nrde2.down$GeneID & res3$GeneID %in% Mtr4.down$GeneID==FALSE, "Nrde2_down",
                              ifelse(res3$GeneID %in% Nrde2.up$GeneID & res3$GeneID %in% Mtr4.up$GeneID==FALSE, "Nrde2_up",
                                     ifelse(res3$GeneID %in% Nrde2.down$GeneID==FALSE & res3$GeneID %in% Mtr4.down$GeneID, "Mtr4_down",
                                            ifelse(res3$GeneID %in% Nrde2.up$GeneID==FALSE & res3$GeneID %in% Mtr4.up$GeneID, "Mtr4_up","other"))))))
table(NMreg)
NMreg2 <- NMreg[res3$GeneID %in% regulated.genes]
table(NMreg2)

pdf(file="Figures/Nrd2_Mtr4_Ccdc174_heatmap_fold_changes_0.01_1_all_regulated_genes_split_v2.pdf",height=15,width=15)
Heatmap(res4,cluster_rows = TRUE, cluster_columns = FALSE,name = "logFC",show_row_dend = FALSE,
        split = NMreg2)
dev.off()

NMreg3 <- NMreg[res3$GeneID %in% unique(res.reg1$GeneID)]
table(NMreg3)

pdf(file="Figures/Nrd2_Mtr4_Ccdc174_heatmap_fold_changes_0.05_1_Nrde2_regulated_genes_split_v2.pdf",height=10,width=10)
ComplexHeatmap::Heatmap(res5,cluster_rows = FALSE, cluster_columns = FALSE,name = "logFC",show_row_dend = FALSE,
                        split = NMreg3)
dev.off()

NMreg <- ifelse(res3$GeneID %in% Nrde2.up$GeneID & res3$GeneID %in% Mtr4.up$GeneID, "Nrde2&Mtr4_up",
                ifelse(res3$GeneID %in% Nrde2.down$GeneID & res3$GeneID %in% Mtr4.down$GeneID, "Nrde2&Mtr4_down", 
                                     ifelse(res3$GeneID %in% Nrde2.down$GeneID==FALSE & res3$GeneID %in% Mtr4.down$GeneID, "Mtr4_down",
                                            ifelse(res3$GeneID %in% Nrde2.up$GeneID==FALSE & res3$GeneID %in% Mtr4.up$GeneID, "Mtr4_up",
                                                   ifelse(res3$GeneID %in% Nrde2.down$GeneID & res3$GeneID %in% Mtr4.down$GeneID==FALSE, "Nrde2_down",
                                                          ifelse(res3$GeneID %in% Nrde2.up$GeneID & res3$GeneID %in% Mtr4.up$GeneID==FALSE, "Nrde2_up","other"))))))
table(NMreg)
NMreg4 <- NMreg[res3$GeneID %in% unique(res.reg2$GeneID)]
table(NMreg4)

pdf(file="Figures/Nrd2_Mtr4_Ccdc174_heatmap_fold_changes_0.05_1_Nrde2_regulated_genes_splitMtr4_v2.pdf",height=10,width=10)
ComplexHeatmap::Heatmap(res6,cluster_rows = FALSE, cluster_columns = FALSE,name = "logFC",show_row_dend = FALSE,
                        split = NMreg4)
dev.off()


#----------------------------------------------------------------------------------------------------------------
#   #pairwise comparison of Nrde2 and Mtr4 for 2C genes
#----------------------------------------------------------------------------------------------------------------


cpm2c <- read.table("public_data/Macfarlan_2012/cpm_per_gene.txt",sep="\t",header = TRUE)
inx <- which(rowSums(cpm2c[,1:2]) > 10)
log2FC <- log2((cpm2c[inx,2]+1)/(cpm2c[inx,1]+1))
plot(density(log2FC))
cpm2c2 <- cpm2c[inx,]

up2Cgenes <- cpm2c2$GeneID[log2FC>8]
res3$genes2C <- ifelse(res3$GeneID %in% up2Cgenes,"2C_upregulated"," ")

#using ggplot
ggplot(res3,aes(x=Mtr4KO_vs_WT,y=Nrde2KO_vs_WT,col=genes2C)) + geom_point(alpha=0.5) + theme_classic() + 
  scale_color_manual(values=c("lightgrey","red")) 
ggsave("Figures/Nrde2_vs_Mtr4_gene_expression_vs_2Cgenes.pdf",device="pdf",height=4,width=8)

#plot red dots second
pdf("Figures/Nrde2_vs_Mtr4_gene_expression_vs_2Cgenes_red_dots_second.pdf",height=5,width=8)
smoothScatter(res3$Mtr4KO_vs_WT,res3$Nrde2KO_vs_WT)
points(res3$Mtr4KO_vs_WT[res3$genes2C=="2C_upregulated"],res3$Nrde2KO_vs_WT[res3$genes2C=="2C_upregulated"],col="red",pch=20)
dev.off()

#plot red dots second, simple
pdf("Figures/Nrde2_vs_Mtr4_gene_expression_vs_2Cgenes_red_dots_second_simple.pdf",height=5,width=8)
plot(res3$Mtr4KO_vs_WT,res3$Nrde2KO_vs_WT,col="lightgrey",pch=20)
points(res3$Mtr4KO_vs_WT[res3$genes2C=="2C_upregulated"],res3$Nrde2KO_vs_WT[res3$genes2C=="2C_upregulated"],col="red",pch=20)
dev.off()

p1 <- ggplot(res3,aes(y=Nrde2KO_vs_WT,x=genes2C,fill=genes2C)) + geom_violin() + theme_classic() + scale_fill_manual(values=c("lightgrey","red")) + ylim(c(-10,10))
p2 <- ggplot(res3,aes(y=Mtr4KO_vs_WT,x=genes2C,fill=genes2C)) + geom_violin() + theme_classic() + scale_fill_manual(values=c("lightgrey","red")) + ylim(c(-10,10))
p1 + p2 
ggsave("Figures/Nrde2_vs_Mtr4_gene_expression_vs_2Cgenes_violin_v2.pdf",device="pdf",height=8,width=8)

#----------------------------------------------------------------------------------------------------------------
#   #pairwise comparison of Nrde2 and Nrde2 mutants
#----------------------------------------------------------------------------------------------------------------

p1 <- ggplot(res3,aes(x=Nrde2KO_vs_WT,y=Nrde2d200_vs_WT,col=genes2C)) + geom_point(alpha=0.5) + theme_classic() + 
  scale_color_manual(values=c("lightgrey","red")) +xlim(range(c(res3$Nrde2KO_vs_WT,res3$Nrde2d200_vs_WT))) +ylim(range(c(res3$Nrde2KO_vs_WT,res3$Nrde2d200_vs_WT)))
p2 <- ggplot(res3,aes(x=Nrde2KO_vs_WT,y=Nrde2D174R_vs_WT,col=genes2C)) + geom_point(alpha=0.5) + theme_classic() + 
  scale_color_manual(values=c("lightgrey","red")) +xlim(range(c(res3$Nrde2KO_vs_WT,res3$Nrde2D174R_vs_WT))) +ylim(range(c(res3$Nrde2KO_vs_WT,res3$Nrde2D174R_vs_WT)))
p1 + p2
ggsave("Figures/Nrde2_vs_d200_D174R_gene_expression_vs_2Cgenes.pdf",device="pdf",height=4,width=10)

#plot red dots second, simple
pdf("Figures/Nrde2_vs_d200_gene_expression_vs_2Cgenes_red_dots_second_simple.pdf",height=5,width=8)
plot(res3$Nrde2d200_vs_WT,res3$Nrde2KO_vs_WT,col="lightgrey",pch=20)
points(res3$Nrde2d200_vs_WT[res3$genes2C=="2C_upregulated"],res3$Nrde2KO_vs_WT[res3$genes2C=="2C_upregulated"],col="red",pch=20)
dev.off()

pdf("Figures/Nrde2_vs_D174R_gene_expression_vs_2Cgenes_red_dots_second_simple.pdf",height=5,width=8)
plot(res3$Nrde2D174R_vs_WT,res3$Nrde2KO_vs_WT,col="lightgrey",pch=20)
points(res3$Nrde2D174R_vs_WT[res3$genes2C=="2C_upregulated"],res3$Nrde2KO_vs_WT[res3$genes2C=="2C_upregulated"],col="red",pch=20)
dev.off()


#----------------------------------------------------------------------------------------------------------------
#   #pairwise comparison of Nrde2 RNAseq and Riboseq
#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#   #overlap of regulated genes with (Nrde2) peaks
#----------------------------------------------------------------------------------------------------------------
#add Nrde2 vs Mtr4 regulation info to res3 table  
res3$NMreg <- NMreg

print(load("clipper_peaks_5readsmin_200202.RData"))
peak.names <- unique(clipper.peaks$sample)[c(4,5,3,1,2)]

peak.list <- list()
for (i in seq_along(peak.names)){
peaks <- clipper.peaks[clipper.peaks$sample==peak.names[i]]
peaks$GeneID <- matrix(unlist(strsplit(peaks$name,split = "_")),ncol=3,byrow=TRUE)[,1]
  
peak.counts <- data.frame(n.no.peak=NA,n.peak=NA)
#count the number of genes without/with a peak
Nrde2_up.genes <- res3$GeneID[res3$NMreg=="Nrde2_up"]
peak.counts[1,] <- table(Nrde2_up.genes %in% peaks$GeneID)
Nrde2_down.genes <- res3$GeneID[res3$NMreg=="Nrde2_down"]
peak.counts[2,] <- table(Nrde2_down.genes %in% peaks$GeneID)
Nrde2Mtr4_up.genes <- res3$GeneID[res3$NMreg=="Nrde2&Mtr4_up"]
peak.counts[3,] <- table(Nrde2Mtr4_up.genes %in% peaks$GeneID)
Nrde2Mtr4_down.genes <- res3$GeneID[res3$NMreg=="Nrde2&Mtr4_down"]
peak.counts[4,] <- table(Nrde2Mtr4_down.genes %in% peaks$GeneID)
peak.counts$genereg <- c("Nrde2_up","Nrde2_down","Nrde2&Mtr4_up","Nrde2&Mtr4_down")
peak.counts$peakname <- peak.names[i]
peak.list[[i]] <- peak.counts
}
  
peaks.at.reggenes <- do.call("rbind",peak.list)  
  
write.table(peaks.at.reggenes,file="FigureData/number_of_genes_with_and_without_peaks_in_Nrde2_or_Nrde2andMtr4_regulated_genes.txt",
            sep="\t",append=FALSE,quote=FALSE, row.names=FALSE,col.names=TRUE)  
  
  



#pathway analysis 

#extract regulated genes: entrez gene
all_genes <- as.character(res2$entrezgene)
res.reg <- dplyr::filter(res2,padj < 0.05, abs(log2FoldChange)>0.5)
go.genes <- list()
for (i in seq_along(contrasts)){
  go.genes[[i]] <- res.reg$entrezgene[res.reg$Contrast==names(contrasts)[i]]
}
names(go.genes) <- names(contrasts)

#run goana
go.de <- goana(go.genes,species="Mm")

View(topGO(go.de, sort = "Nrde2KO_vs_WT"))
write.table(go.de,file="Nrd2_Mtr4_Ccdc174_GOterms_gene_expression_regulated_genes_p0.05_FC0.5_cutoff30reads.txt",sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

#extract regulated genes: ensembl gene
all_genes <- as.character(res2$GeneID)
res.reg <- dplyr::filter(res2,padj < 0.05, abs(log2FoldChange)>0.5)
go.genes <- list()
for (i in seq_along(contrasts)){
  go.genes[[i]] <- res.reg$GeneID[res.reg$Contrast==names(contrasts)[i]]
}
names(go.genes) <- names(contrasts)

# Run GO enrichment analysis with clusterProfiler
library(clusterProfiler)
library(org.Mm.eg.db)
plotlist <- list()
for (i in seq_along(names(contrasts))){
ego <- enrichGO(gene = go.genes[[i]], universe = all_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Mm.eg.db,
                ont = "BP",
                pvalueCutoff = 0.1,
                pAdjustMethod = "BH",
                qvalueCutoff = 0.2,
                readable = TRUE)
# Output results from GO analysis to a table
#cluster_summary <- data.frame(ego)
#cluster_summary
#dotplot(ego, showCategory=50)
#png(file=sprintf("RNAseq/GOterms_BP_genes_0.05_FC0.5_plot_%s.png",names(contrasts)[i]))
plotlist[[i]] <- emapplot(ego, showCategory=50)
#dev.off()
}

library(cowplot)
cowplot::plot_grid(plotlist=plotlist, nrow=2,ncol=5, labels=names(contrasts))
ggsave("RNAseq/GOterms_BP_genes_0.05_FC0.5_plot_allcontrasts.png",height=20,width=80,limitsize = FALSE)

#GSEA analysis
geneList <- list()
res2u <- res2[!is.na(res2$entrezgene),]
for (i in seq_along(contrasts)){
  log2FoldChange <- res2u$log2FoldChange[res2u$Contrast==names(contrasts)[i]]
  names(log2FoldChange) <- res2u$entrezgene[res2u$Contrast==names(contrasts)[i]]
  log2FoldChange <- log2FoldChange[unique(names(log2FoldChange))]
  geneList[[i]] <- sort(log2FoldChange, decreasing = TRUE)
}
names(geneList) <- names(contrasts)
plotlist <- list()

for (i in seq_along(names(contrasts))){
  ego3 <- gseGO(geneList     = geneList[[i]],
                OrgDb        = org.Mm.eg.db,
                ont          = "BP",
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 0.1,
                verbose      = FALSE)
  
  plotlist[[i]] <- emapplot(ego3, showCategory=50)
}
cowplot::plot_grid(plotlist=plotlist, nrow=2,ncol=5, labels=names(contrasts))
ggsave("RNAseq/GOterms_BP_genes_GSEAonFC_plot_allcontrasts.png",height=20,width=80,limitsize = FALSE)

#Msigdb
library(msigdbr)
m_t2g <- msigdbr(species = "Mus musculus", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

plotlist <- list()
for (i in seq_along(names(contrasts))){
  em2 <- GSEA(geneList[[i]], TERM2GENE = m_t2g,pvalueCutoff = 0.5)
  plotlist[[i]] <- emapplot(em2, showCategory=50)
}
cowplot::plot_grid(plotlist=plotlist, nrow=2,ncol=5, labels=names(contrasts))
ggsave("RNAseq/MSIGDB_C2_genes_GSEAonFC_plot_allcontrasts.png",height=20,width=80,limitsize = FALSE)


# To color genes by log2 fold changes
signif_res <- res.reg[res.reg$Contrast==names(contrasts)[2] & !is.na(res.reg$gene_symbol),]
row.names(signif_res) <- signif_res$gene_symbol
signif_res_lFC <- signif_res$log2FoldChange
names(signif_res_lFC) <- row.names(signif_res)
cnetplot(ego,
         categorySize="pvalue",
         showCategory = 100,
         foldChange= signif_res_lFC, vertex.label.font=6)

heatplot(ego, foldChange=signif_res_lFC,showCategory = 100)
#-----------------
#add gene biotyoe info to find snRNAs
#-----------------
res2 <- read.table("RNAseq/Nrde2_Mtr4_Ccdc174_RNAseq_DEseq2_results_0.05_0.5FC_cutoff30reads.txt",sep="\t",header=TRUE)

#get biotype info
library(EnsDb.Mmusculus.v79)
edb <- EnsDb.Mmusculus.v79
edbg <- data.frame(genes(edb))
res3 <- left_join(res2,edbg,by=c("GeneID"="gene_id"))
inx.sn <- which(res3$gene_biotype=="snRNA")
res3.sn <- res3[inx.sn,]
res3.sn[res3.sn$padj < 0.1 & abs(res3.sn$log2FoldChange) > 0.5 & res3.sn$Contrast=="Nrde2KO_vs_WT",]

#Volcanoplot
vpa <- ggplot(res3.sn,aes(x=log2FoldChange,y=log.padj)) +geom_point(aes()) + 
  theme_classic() +facet_wrap(vars(Contrast))
vpa <- vpa + geom_vline(xintercept=c(-1,1),alpha=0.4) + geom_hline(yintercept = 2,alpha=0.4)
vpa
ggsave("RNAseq/Nrde2_Mtr4_Ccdc174_RNAseq_volcano_plots_cutoff30reads_snRNAs.png",device="png",height=10,width=10)


res3.sn$snRNA.gene <- ifelse(grepl("Gm",res3.sn$symbol),"pseudogene","gene")
table(res3.sn$snRNA.gene)
vpa <- ggplot(res3.sn,aes(x=log2FoldChange,y=log.padj)) +geom_point(aes(col=log2(baseMean))) + 
  theme_classic() +facet_grid(Contrast ~ snRNA.gene)
vpa <- vpa + geom_vline(xintercept=c(-1,1),alpha=0.4) + geom_hline(yintercept = 2,alpha=0.4)
vpa
ggsave("RNAseq/Nrde2_Mtr4_Ccdc174_RNAseq_volcano_plots_cutoff30reads_snRNAs_split_pseudogenes.png",device="png",height=10,width=6)

res3.sn.wide <- pivot_wider(res3.sn,id_cols = c(symbol,snRNA.gene),names_from = Contrast,values_from = log2FoldChange)
library(ComplexHeatmap)
Heatmap(res3.sn.wide[3:length(res3.sn.wide)],split = res3.sn.wide$snRNA.gene)

#numbers of biotypes
contrasts <- unique(res3$Contrast)
type.numbers <- data.frame(table(res3$gene_biotype))

for (contrast in seq_along(contrasts)){
table(res3$gene_biotype[res3$regulated=="up" & res3$Contrast==contrasts[contrast]])
table(res3$gene_biotype[res3$regulated=="down" & res3$Contrast==contrasts[contrast]])
}

#cehck single genes
inx <- which(res2$gene_symbol=="Trp53")
inx <- which(res2$gene_symbol=="Tex19.2")

View(res2[inx,])


#

#----------------------------------------------------------------------------------------------------------------
#   # Sushi plots of interesting regions####
#----------------------------------------------------------------------------------------------------------------
#files
bwFilesPlus <- c(list.files("/tungstenfs/scratch/gbuehler/flemrmat/data/2456/raw_data/star_trackhub/", full.names=TRUE,pattern="*_plus.bw$"))
bwFilesMinus <- c(list.files("/tungstenfs/scratch/gbuehler/flemrmat/data/2456/raw_data/star_trackhub/", full.names=TRUE,pattern="*_minus.bw$"))
bedFiles <- c(list.files(".", full.names=TRUE,pattern="*_q90_size60_stranded.bed$"))
bwNames <- gsub("/tungstenfs/scratch/gbuehler/flemrmat/data/2456/raw_data/star_trackhub//2456_BCLIP_","",bwFilesMinus)
bwNames <- gsub("_STARnorm_uni_minus.bw","",bwNames)

#files
bwFilesPlus <- c(list.files("/tungstenfs/scratch/gbuehler/deepSeqRepos/bigwig/", full.names=TRUE,pattern="*_multi_plus.bw$"))
bwFilesPlus <- grep("spliced_2191F",bwFilesPlus,value=TRUE)
bwFilesPlus <- grep("Nxf1dT",bwFilesPlus,value=TRUE,invert=TRUE)

bwFilesMinus <- c(list.files("/tungstenfs/scratch/gbuehler/deepSeqRepos/bigwig/", full.names=TRUE,pattern="_multi_minus.bw$"))
bwFilesMinus <- grep("spliced_2191F",bwFilesMinus,value=TRUE)
bwFilesMinus <- grep("Nxf1dT",bwFilesMinus,value=TRUE,invert=TRUE)

bwNames <- gsub("/tungstenfs/scratch/gbuehler/deepSeqRepos/bigwig//","",bwFilesMinus)
bwNames <- gsub("_multi_minus.bw","",bwNames)

#regions
res3.sn[res3.sn$log2FoldChange > 2 & res3.sn$Contrast=="Nrde2KO_vs_WT",]

library(MiniChip)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
png("RNAseq/RNA_sushi_plot_of_upregulated_snRNAs_Gm22485.png",width=1200,height=2000)
plotTracks(
  bwFilesPlus,
  bwFilesMinus,
  bwNames=bwNames,
  txdb=TxDb.Mmusculus.UCSC.mm10.ensGene,
  EnsDb=EnsDb.Mmusculus.v79,
  #bedFiles=bedFiles,
  gene="Gm22485",
  plotExtension = 500,
  plotranges = rep(0.5, length(bwNames))
)
dev.off()


