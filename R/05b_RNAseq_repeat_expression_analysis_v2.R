rm(list=ls())


## Libraries
library(tidyverse)
#library(DESeq2)
library(GenomicRanges)
library(Rsubread)
library(limma)
library(edgeR)
library(patchwork)
library(GenomicFeatures)


setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")

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
bamFilesR2 <- grep("2447F",all.bamFiles,value=TRUE)
bamFilesR <- c(bamFilesR,bamFilesR2)

bamNamesR <- gsub("/tungstenfs/scratch/gbuehler/deepSeqRepos/bam//","",bamFilesR)
bamNamesR <- gsub("_Aligned.sortedByCoord.out.bam","",bamNamesR)

bamFiles <- bamFilesR[c(1:4,7:16,19:20,23:24,27:28)]
bamNames <- bamNamesR[c(1:4,7:16,19:20,23:24,27:28)]

#' make repeat saf table
repeats <- read.table("/tungstenfs/scratch/gbuehler/fabiom/manuHisCript/repeats/repeatsMM10_UCSC_1792018.bed")
names(repeats) <- c("chr","start","end","repeat_name","repeat_class","strand","swScore","milliDiv","milliDel","milliIns")
reps <- makeGRangesFromDataFrame(repeats,
                                 keep.extra.columns=TRUE,
                                 ignore.strand=FALSE,
                                 seqinfo=NULL,
                                 seqnames.field=c("chr"),
                                 start.field=c("start"),
                                 end.field=c("end"),
                                 strand.field=c("strand"),
                                 starts.in.df.are.0based=TRUE)
table(reps$repeat_class)
inx <- c(grep("LTR",reps$repeat_class),grep("LINE",reps$repeat_class),grep("SINE",reps$repeat_class),grep("DNA",reps$repeat_class))
reps2 <- reps[inx]
table(reps2$repeat_class)
names(reps2) <- paste(reps2$repeat_name,seqnames(reps2),sep="_")
names(reps2) <- paste(names(reps2),start(reps2),sep="_")
length(names(reps2)) == length(unique(names(reps2)))

#remove repeats that overlap genes
txdb=loadDb("/tungstenfs/scratch/gbuehler/michi/Annotations/GENCODE/Mouse/release_M23/gencode.vM23.annotation.txdb.sqlite")
genes <- genes(txdb)
reps2 <- reps2[overlapsAny(reps2,genes)==FALSE]

#' generate  saf format data frame
saf <- data.frame(GeneID= names(reps2), Chr=seqnames(reps2),
                  Start=start(reps2), End=end(reps2),Strand=strand(reps2))



#calculate the number of reads per gene 
f_counts <- featureCounts(bamFiles,annot.ext=saf,isGTFAnnotationFile = FALSE,
                           useMetaFeatures=TRUE,allowMultiOverlap=TRUE,
                           minOverlap=1,countMultiMappingReads=FALSE,fraction=FALSE,
                           minMQS=255,strandSpecific=2,nthreads=20,verbose=FALSE,isPairedEnd=FALSE)


counts <- data.frame(f_counts$counts)
colnames(counts) <- bamNames

write.table(counts,"FigureData/feature_counts_per_repeat.txt",sep="\t",col.names = TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

#calculate CPMs and FPKMs
fpkm <- f_counts$counts
cpm <- f_counts$counts


for (i in seq_along(bamFiles)){
  scaling_factor <-  f_counts$stat[f_counts$stat$Status=="Assigned" | f_counts$stat$Status=="Unassigned_NoFeatures",i+1]/1e6
  cpm[,i] <- f_counts$counts[,i]/scaling_factor
  fpkm[,i] <- cpm[,i]/(f_counts$annotation$Length/1000)
}

cpm <- data.frame(cpm)
colnames(cpm) <- bamNames
cpm$ID <- row.names(cpm)
reps2$Length <- width(reps2)
cpm2 <- left_join(cpm,data.frame(ID=names(reps2),mcols(reps2)),by="ID")
write.table(cpm2,"FigureData/feature_cpm_per_repeat.txt",sep="\t",col.names = TRUE,row.names=FALSE,append=FALSE,quote=FALSE)


#----------------------------------------------------------------------------------------------------------------
#   #calculate differential gene expression
#----------------------------------------------------------------------------------------------------------------
#counts <- read.table("RNAseq/feature_counts_per_gene.txt",sep="\t",header=TRUE)

# prepare annotations
Group <- c(rep("Mtr4_KO",2),rep("Mtr4_WT",2),rep("Nrde2_D174R",2),
           rep("Nrde2_d200",2),rep("Nrde2_KO",2),rep("Nrde2_WT",2),
           rep("Ccdc174_KD",2), rep("Ccdc174_WT",2),
           rep("Nrde2_KO",2), rep("Nrde2_WT",2))

annots <- data.frame(sample=colnames(counts),Group=Group,
                     replicate=rep(c("rep1","rep2"),10),batch=c(rep("batch1",12),rep("batch2",8)))
rownames(annots) <- annots$sample

#extract species name
species <- "Mus musculus"


# check that everything is fine
#counts <- counts[,rownames(annots)] # reorder raw.counts after annots
identical(rownames(annots),colnames(counts))

# prepare features
features <- data.frame(reps2)[,c(6:11)]
rownames(features) <- names(reps2)



mapped.reads <- apply(f_counts$stat[c(1,12),-1],2,sum)
names(mapped.reads) <- bamNames

library(limma)
library(edgeR)
d0 <- DGEList(counts,lib.size=mapped.reads)
d0 <- calcNormFactors(d0,method="TMMwsp")

#filter peaks with low number of reads
cutoff <- 0.5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
length(drop)
d0b <- d0[-drop,] 

d <- d0b 
dim(d) # number of genes left


plotMDS(d, top = 1000, labels=Group,col = as.numeric(Group))
mm <- model.matrix(~0 + Group)
y <- voom(d, mm, plot = T)

fit <- lmFit(y, mm)
head(coef(fit))




contrasts <- list(makeContrasts(Nrde2KO_vs_WT=GroupNrde2_KO - GroupNrde2_WT, levels = colnames(coef(fit))),
                  makeContrasts(Nrde2d200_vs_WT=GroupNrde2_d200 - GroupNrde2_WT, levels = colnames(coef(fit))),
                  makeContrasts(Nrde2D174R_vs_WT=GroupNrde2_D174R - GroupNrde2_WT, levels = colnames(coef(fit))),
                  makeContrasts(Ccdc174KD_vs_WT=GroupCcdc174_KD - GroupCcdc174_WT, levels = colnames(coef(fit))),
                  makeContrasts(Mtr4KO_vs_WT=GroupMtr4_KO - GroupMtr4_WT, levels = colnames(coef(fit)))
                 )
                  
names(contrasts) <- c("Nrde2KO_vs_WT","Nrde2d200_vs_WT","Nrde2D174R_vs_WT","Ccdc174KD_vs_WT","Mtr4KO_vs_WT")

res <- list()
for (i in seq_along(contrasts)){
  tmp <- contrasts.fit(fit, contrasts[[i]])
  tmp <- eBayes(tmp)
  limma::plotMA(tmp)
  
  top.table <- as.data.frame(topTable(tmp, sort.by = "none", n = Inf))
  top.table$log.adj.P.Val <- -log10(top.table$adj.P.Val)
  top.table$Contrast <- names(contrasts)[i]
  top.table$ID <- row.names(top.table)
  
  res[[i]] <- top.table
}
res2 <- do.call("rbind",res)


res2 <- left_join(res2,data.frame(cbind(names(reps2),data.frame(reps2))),by=c("ID"="names.reps2."))
res2$regulated <- ifelse(res2$logFC > 1 & res2$adj.P.Val < 0.05,"up",ifelse(res2$logFC < -1 & res2$adj.P.Val < 0.05,"down","no"))
table(res2$regulated)
write.table(res2,file="FigureData/Nrde2_Mtr4_Ccdc174_RNAseq_Limma_repeats_results_cutoff0.2.txt",sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
res2 <- read.table("FigureData/Nrde2_Mtr4_Ccdc174_RNAseq_Limma_repeats_results_cutoff0.2.txt",sep="\t",header=TRUE)


#Volcanoplot
res2$repeat.class <- matrix(unlist(strsplit(res2$repeat_class,"/",fixed=TRUE)),ncol=2,byrow=TRUE)[,1]

vpa <- ggplot(res2,aes(x=logFC,y=log.adj.P.Val)) +geom_point(aes(color=repeat.class)) + theme_classic() +facet_wrap(vars(Contrast))
vpa <- vpa + geom_vline(xintercept=c(-1,1),alpha=0.4) + geom_hline(yintercept = 2,alpha=0.4)
vpa
ggsave("Figures/Nrde2_Mtr4_Ccdc174_RNAseq_Repeat_volcano_plots_cutoff0.2.png",device="png",height=10,width=10)



#select upregulated repeats
res.reg1 <- dplyr::filter(res2,adj.P.Val < 0.01, (logFC)>1 & Contrast== "Nrde2KO_vs_WT")
res.reg2 <- dplyr::filter(res2,adj.P.Val < 0.01, (logFC)>1 & Contrast== "Mtr4KO_vs_WT")
res.reg3 <- dplyr::filter(res2,adj.P.Val < 0.01, (logFC)>1 & Contrast== "Nrde2D174R_vs_WT")
res.reg4 <- dplyr::filter(res2,adj.P.Val < 0.01, (logFC)>1 & Contrast== "Nrde2d200_vs_WT")
res.reg5 <- dplyr::filter(res2,adj.P.Val < 0.01, (logFC)>1 & Contrast== "Ccdc174KD_vs_WT")
#regulated.genes <- unique(c(res.reg1$GeneID,res.reg2$GeneID,res.reg3$GeneID,res.reg4$GeneID,res.reg5$GeneID))

#table of all selected repeat names
res2.names <- data.frame(table(res2$repeat_name[res2$Contrast=="Nrde2KO_vs_WT"]))
res2.names <- res2.names[order(res2.names$Freq,decreasing=TRUE),]
res2.names2 <- res2.names[1:10,]

#table of upregulated repeat names
Nrde2.up.names <- data.frame(table(res.reg1$repeat_name))
Nrde2.up.names <- Nrde2.up.names[order(Nrde2.up.names$Freq,decreasing=TRUE),]

d200.up.names <- data.frame(table(res.reg4$repeat_name))
d200.up.names <- d200.up.names[order(d200.up.names$Freq,decreasing=TRUE),]

D174R.up.names <- data.frame(table(res.reg3$repeat_name))
D174R.up.names <- D174R.up.names[order(D174R.up.names$Freq,decreasing=TRUE),]

Mtr4.up.names <- data.frame(table(res.reg2$repeat_name))
Mtr4.up.names <- Mtr4.up.names[order(Mtr4.up.names$Freq,decreasing=TRUE),]


#merge tables
res.names.all <- full_join(Nrde2.up.names,d200.up.names,by="Var1")
res.names.all <- full_join(res.names.all,D174R.up.names,by="Var1")
res.names.all <- full_join(res.names.all,Mtr4.up.names,by="Var1")
res.names.all <- full_join(res.names.all,res2.names,by="Var1")

names(res.names.all) <- c("repeat_name","Nrde2KO_vs_WT","Nrde2d200_vs_WT","Nrde2D174R_vs_WT","Mtr4KO_vs_WT","all_cpm0.5")
#replace nas with 0
res.names.all$Nrde2KO_vs_WT <- replace_na(res.names.all$Nrde2KO_vs_WT,replace=0)
res.names.all$Nrde2d200_vs_WT <- replace_na(res.names.all$Nrde2d200_vs_WT,replace=0)
res.names.all$Nrde2D174R_vs_WT <- replace_na(res.names.all$Nrde2D174R_vs_WT,replace=0)
res.names.all$Mtr4KO_vs_WT <- replace_na(res.names.all$Mtr4KO_vs_WT,replace=0)


#add repeat_class
repname2repclass <- unique(res2[,15:16])
res.names.all <- left_join(res.names.all,repname2repclass,by="repeat_name")
res.names.all$repeat.class <- matrix(unlist(strsplit(res.names.all$repeat_class,"/",fixed=TRUE)),ncol=2,byrow=TRUE)[,1]

Nrde2KO_vs_WT.ratio <- log2((res.names.all$Nrde2KO_vs_WT+1)/(res.names.all$all_cpm0.5+1))
Nrde2d200_vs_WT.ratio <- log2((res.names.all$Nrde2d200_vs_WT+1)/(res.names.all$all_cpm0.5+1))
Nrde2D174R_vs_WT.ratio <- log2((res.names.all$Nrde2D174R_vs_WT+1)/(res.names.all$all_cpm0.5+1))
Mtr4KO_vs_WT.ratio <- log2((res.names.all$Mtr4KO_vs_WT+1)/(res.names.all$all_cpm0.5+1))

#eg plot upregulated vs total in the table or vs total in genome 
p1 <- ggplot(res.names.all,aes(x=all_cpm0.5,y=Nrde2KO_vs_WT,color=repeat.class)) + geom_point() + theme_classic() +
  geom_text(aes(label=ifelse(Nrde2KO_vs_WT > 10 & Nrde2KO_vs_WT.ratio > -4,as.character(repeat_name),'')),hjust=0.5,vjust=-0.7,size=2)
p2 <- ggplot(res.names.all,aes(x=all_cpm0.5,y=Nrde2d200_vs_WT,color=repeat.class)) + geom_point() + theme_classic() +
  geom_text(aes(label=ifelse(Nrde2d200_vs_WT > 5 & Nrde2d200_vs_WT.ratio > -5,as.character(repeat_name),'')),hjust=0.5,vjust=-0.7,size=2)
p3 <- ggplot(res.names.all,aes(x=all_cpm0.5,y=Nrde2D174R_vs_WT,color=repeat.class)) + geom_point() + theme_classic() +
  geom_text(aes(label=ifelse(Nrde2D174R_vs_WT > 1 & Nrde2D174R_vs_WT.ratio > -5,as.character(repeat_name),'')),hjust=0.5,vjust=-0.7,size=2)
p4 <- ggplot(res.names.all,aes(x=all_cpm0.5,y=Mtr4KO_vs_WT,color=repeat.class)) + geom_point() + theme_classic() +
  geom_text(aes(label=ifelse(Mtr4KO_vs_WT > 50 & Mtr4KO_vs_WT.ratio > -3,as.character(repeat_name),'')),hjust=0.5,vjust=-0.7,size=2)

p1 + p2 + p3 + p4          
ggsave(filename="Figures/repeat_names_at_Nrde2_Mtr4_upregulated_numbers_limma0.2.pdf",height=8, width=12,device = "pdf")


#directly compare Mtr4 and Nrde2 repeat numbers
#res.names.all$repeat.class <- gsub("DNA?","DNA",res.names.all$repeat.class)
#res.names.all$repeat.class <- gsub("LTR?","LTR",res.names.all$repeat.class)

p1 <- ggplot(res.names.all,aes(x=Nrde2KO_vs_WT,y=Mtr4KO_vs_WT,shape=repeat.class,color=all_cpm0.5)) + geom_point() + theme_classic() +
  geom_text(aes(label=ifelse(Nrde2KO_vs_WT > 20 & Mtr4KO_vs_WT > 20,as.character(repeat_name),'')),hjust=0.5,vjust=-0.7,size=2)
p2 <- ggplot(res.names.all,aes(x=Nrde2KO_vs_WT,y=Nrde2d200_vs_WT,shape=repeat.class,color=all_cpm0.5)) + geom_point() + theme_classic() +
  geom_text(aes(label=ifelse(Nrde2d200_vs_WT > 5 | Nrde2KO_vs_WT > 20,as.character(repeat_name),'')),hjust=0.5,vjust=-0.7,size=2)
p3 <- ggplot(res.names.all,aes(x=Nrde2KO_vs_WT,y=Nrde2D174R_vs_WT,shape=repeat.class,color=all_cpm0.5)) + geom_point() + theme_classic() +
  geom_text(aes(label=ifelse(Nrde2D174R_vs_WT > 1 | Nrde2KO_vs_WT >20,as.character(repeat_name),'')),hjust=0.5,vjust=-0.7,size=2)

p1 + p2 + p3
ggsave(filename="Figures/repeat_names_at_Nrde2_vs_Mtr4_upregulated_numbers_limma0.2.pdf",height=4, width=14,device = "pdf")

#----------------------------------------------------------------------------------------------------------------
#'   #plot Nrde2 vs Mtr4 logFC####
#----------------------------------------------------------------------------------------------------------------
#select the oens which ahve > 4 occurances in cpm selected repeats
res2b <- res2[res2$repeat_name %in% res2.names$Var1[res2.names$Freq>4],]


#summarize logFC by repeat name
resFC <- res2b %>% group_by(repeat_name,Contrast) %>% summarize(aveLogFC=mean(logFC),repeat_name=repeat_name,
                                                      repeat.class=repeat.class,Contrast=Contrast)

resFC <- unique(resFC) %>% dplyr::filter(Contrast %in% c("Nrde2KO_vs_WT","Mtr4KO_vs_WT")) %>%
  pivot_wider(id_cols = c("repeat_name", "Contrast","repeat.class"),names_from=Contrast,values_from=aveLogFC)

ggplot(resFC,aes(x=Mtr4KO_vs_WT,y=Nrde2KO_vs_WT,col=repeat.class)) + geom_point() + theme_classic() +
  geom_text(aes(label=ifelse(Mtr4KO_vs_WT > 2.5 & Nrde2KO_vs_WT >1,as.character(repeat_name),'')),hjust=0.5,vjust=-0.7,size=3)
ggsave(filename="Figures/repeat_names_at_Nrde2_vs_Mtr4_upregulated_aveLogFC_limma0.2.pdf",height=4, width=6,device = "pdf")




library(viridis)
smoothScatter(res2$logFC[res2$Contrast=="Nrde2KO_vs_WT"],res2$logFC[res2$Contrast=="Mtr4KO_vs_WT"],colramp=viridis)

plot(res2$logFC[res2$Contrast=="Nrde2KO_vs_WT"],res2$logFC[res2$Contrast=="Mtr4KO_vs_WT"])

#----------------------------------------------------------------------------------------------------------------
#'   #boxplots of CPMs####
#----------------------------------------------------------------------------------------------------------------
TPMs <- read.table("FigureData/feature_cpm_per_repeat.txt",sep="\t",header=TRUE)

TPMs2 <- dplyr::filter(TPMs,repeat_name=="MT2_Mm" & Length > 450) 
TPMs2 <- dplyr::filter(TPMs,repeat_name=="MERVL-int" & Length > 4000) 


TPMs2 <- TPMs2[rowSums(TPMs2[,1:20])>0,]
TPMs2 <- TPMs2[rowSums(TPMs2[,1:20])>0 & apply(TPMs2[,1:20],1,max)>10,]

TPMs2 <- pivot_longer(TPMs2,cols=colnames(TPMs)[1:20],names_to="samples",values_to="TPM")
TPMs2$logCPM <- log2(TPMs2$TPM+1)
annots$sample <- gsub("-",".",annots$sample)
TPMs2 <- left_join(TPMs2,annots,by=c("samples"="sample"))
TPMs2$Group2 <- ifelse(grepl("_WT",TPMs2$Group)==TRUE,"WT",as.character(TPMs2$Group))

TPMs2$Group2 <- factor(TPMs2$Group2,
                      levels = c("Nrde2_KO",
                                 "Nrde2_d200","Nrde2_D174R",
                                 "Mtr4_KO","Ccdc174_KD","WT"),ordered = TRUE)

#plotcols <- c("#046C9A","#046C9A","#5BBCD6","#C6CDF7","#FD6464", "#FD6464", "#F98400","#F98400")
plotcols <- c("#046C9A","#5BBCD6","#C6CDF7", "#FD6464", "#F98400","#00A08A")

ggplot(TPMs2,aes(x=Group2,y=logCPM,color=Group2)) + geom_boxplot(notch=TRUE) + theme_classic() +
  geom_jitter(color="black", size=0.4, alpha=0.2) + scale_color_manual(values=plotcols)
ggsave(filename="Figures/MERVL_Mm_CPM_notch_boxplots_per_group_unique_reads_10maxcpm.pdf",device="pdf",height=5,width=8)


