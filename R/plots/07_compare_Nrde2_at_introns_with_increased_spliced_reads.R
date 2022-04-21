
rm(list=ls())

#install.packages("/tungstenfs/scratch/gbuehler/bioinfo/Rpackages/MiniChip_0.0.0.9000.tar.gz",repos=NULL)
library(tidyverse)
library(Rsubread)
library(GenomicFeatures)
library(GenomicRanges)
library(MiniChip)
library(patchwork)
library(ggpubr)
library(ComplexHeatmap)

setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")

#load intronic read differential expression
res2 <- read.table(file="FigureData/Nrde2_Mtr4_Ccdc174_batch1and2_intronicRNAseq_DEseq2_results_0.05_0.5FC.txt",sep="\t",header=TRUE)
#add GeneID column
load("splice_sites_expressed_and_notexpressed_transcripts_log2TPMabove0.1_and_genes_FPKMabove8below0.RData")
names(spliceDonors) <- paste(names(spliceDonors),start(spliceDonors),sep="_")
spliceDonors.df <- data.frame(spliceDonors)
spliceDonors.df$ID <- names(spliceDonors)

res2 <- left_join(res2,spliceDonors.df,by="ID")
#----------------------------------------------------------------------------------------------------------------
#   #add the gene expression DEseq data 
#-----------------------------------------------------------------------------------------------------------------
res1 <- read.table("FigureData/Nrde2_Mtr4_Ccdc174_RNAseq_DEseq2_results.txt",sep="\t",header=TRUE)
res2$GeneID <- matrix(unlist(strsplit(as.character(res2$GeneID),".",fixed=TRUE)),ncol=2,byrow = TRUE)[,1]
res2$GeneConID <- paste(res2$GeneID,res2$Contrast,sep="_")
res1$GeneConID <- paste(res1$GeneID,res1$Contrast,sep="_")

res3 <- left_join(res2,res1,by="GeneConID")
#----------------------------------------------------------------------------------------------------------------
#   # select upregulated introns in downregulated genes and plot BCLIP levels in them
#-----------------------------------------------------------------------------------------------------------------
contrasts <- unique(res3$Contrast.x)
#res3 <- res3[is.na(res3$intron_regulated)==FALSE,]

res3$intron_gene_regulated <- ifelse(res3$intron_regulated =="up" & res3$log2FoldChange.y < 0,"Intron up, gene down",
       ifelse(res3$log2FoldChange.x < 0 & res3$log2FoldChange.y < 0,"Intron down, gene down",NA))
table(res3$intron_gene_regulated)
#which genes have increased splcied reads in introns
unique(res3$gene_symbol[res3$Contrast.x=="Nrde2KO_vs_WT" & res3$intron_gene_regulated =="Intron up, gene down"])

#add repeat information
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
names(introns3g) <- names(spliceDonors)
introns3g.reps <- subsetByOverlaps(introns3g,reps)
res3$reps <- ifelse(res3$ID %in% names(introns3g.reps),"repeats","no repeats")
table(res3$reps)

#----------------------------------------------------------------------------------------------------------------
#   # load BCLIP peaks and intersect them with 5' SS
#-----------------------------------------------------------------------------------------------------------------

load("clipper_peaks_5readsmin_200202.RData")


#intersect them with 5' SS
res3$Nrde2.peaks <- ifelse( res3$ID %in% names(spliceDonors)[overlapsAny(promoters(spliceDonors,upstream=100,downstream=10),clipper.peaks[clipper.peaks$sample=="Nrde2"])],"Nrde2Peak","noPeak")
res3$Ccdc174.peaks <- ifelse( res3$ID %in% names(spliceDonors)[overlapsAny(promoters(spliceDonors,upstream=100,downstream=10),clipper.peaks[clipper.peaks$sample=="Ccdc174"])],"Ccdc174Peak","noPeak")
res3$Mtr4.peaks <- ifelse( res3$ID %in% names(spliceDonors)[overlapsAny(promoters(spliceDonors,upstream=100,downstream=10),clipper.peaks[clipper.peaks$sample=="Mtr4"])],"Mtr4Peak","noPeak")
res3$Eif4a3.peaks <- ifelse( res3$ID %in% names(spliceDonors)[overlapsAny(promoters(spliceDonors,upstream=100,downstream=10),clipper.peaks[clipper.peaks$sample=="Eif4a3"])],"Eif4a3Peak","noPeak")

table(res3$Nrde2.peaks)

#----------------------------------------------------------------------------------------------------------------
#   # plot RNAseq FC at peaks vs no peaks version 1 
#-----------------------------------------------------------------------------------------------------------------

plist <- list(NA,NA,NA,NA,NA,NA)
for (i in 1:length(contrasts)){
  
res4 <- res3[res3$Contrast.x==contrasts[i] & res3$log2FoldChange.y < 0,]
pvalues <- compare_means(log2FoldChange.x ~ Nrde2.peaks, data = res4)

p <- ggplot(res4,aes(x=log2FoldChange.x,group=Nrde2.peaks)) + geom_density(aes(col=Nrde2.peaks)) +
  theme_classic() + scale_color_manual(values=c("darkgrey","darkblue")) + 
  ggtitle(contrasts[i]) + labs(subtitle = sprintf("p=%s",pvalues$p.adj))
p.data <- ggplot_build(p)$data[[1]]
p.text <- lapply(split(p.data, f = p.data$group), function(df){
  df[which.max(df$scaled), ]
})
p.text <- do.call(rbind, p.text)
if(nrow(p.text) > 1){
  if(abs(p.text$x[1]-p.text$x[2]) < 1){
    p.text$x[1] <- p.text$x[1] -1
    p.text$x[2] <- p.text$x[2] +1
  }}
plist[[i]] <- p + annotate('text', x = p.text$x, y = p.text$y,
                           label = sprintf('n = %d', p.text$n), vjust = 0)
}
plist[[1]] + plist[[2]] + plist[[3]] + plist[[4]] + plist[[5]] + plist[[6]]

#----------------------------------------------------------------------------------------------------------------
#   # plot RNAseq FC at peaks vs no peaks version 2
#-----------------------------------------------------------------------------------------------------------------
library(wesanderson)
#define colors
plotcols1 <- wes_palette("Darjeeling1")
plotcols2 <- wes_palette("Darjeeling2")

plotcols <- c(plotcols1[c(5,2)],plotcols2[c(2)])

res4 <- rbind(res3[res3$Nrde2.peaks=="Nrde2Peak",],
                   res3[res3$Ccdc174.peaks=="Ccdc174Peak", ], 
                        res3[res3$Eif4a3.peaks=="Eif4a3Peak",])
res4$peaks <- c(rep("Nrde2",nrow(res3[res3$Nrde2.peaks=="Nrde2Peak",])),
                rep("Ccdc174",nrow(res3[res3$Ccdc174.peaks=="Ccdc174Peak", ])),
                rep("Eif4a3",nrow(res3[res3$Eif4a3.peaks=="Eif4a3Peak",]))
                  )

plist <- list(NA,NA,NA,NA,NA,NA)

for (i in 1:length(contrasts)){
  
  res5 <- res4[res4$Contrast.x==contrasts[i] & res4$log2FoldChange.y < 0,]
  
p <- ggplot(res5,aes(x=log2FoldChange.x,col=peaks)) + geom_density(aes(fill=peaks),alpha=0.2) + theme_classic() +
    scale_color_manual(values=plotcols) + scale_fill_manual(values=plotcols) + ggtitle(contrasts[i]) +
    xlab("log2 Intronic FoldChange")
plist[[i]] <- p
}

plist[[1]] + plist[[2]] + plist[[3]] + plist[[4]] + plist[[5]] + plist[[6]]
ggsave("Figures/Intronic_Read_logFC_in_Nrde2_Ccdc174_Eif4a3_peaks.pdf",height=8,width=14)

#----------------------------------------------------------------------------------------------------------------
#   # plot RNAseq FC at peaks vs no peaks version 3
#-----------------------------------------------------------------------------------------------------------------
plotcols3 <- c("#046C9A","#F98400","#FD6464","#00A08A")

res4 <- dplyr::filter(res3,reps=="no repeats")

p1 <- dplyr::filter(res4,Contrast.x %in% c("Nrde2KO_vs_WT","Ccdc174KD_vs_WT","Mtr4KO_vs_WT")) %>%
ggplot(aes(x=log2FoldChange.x,y=log.padj.x)) +geom_point(aes(color=Nrde2.peaks,size=Nrde2.peaks,alpha=Nrde2.peaks)) + 
  facet_grid(vars(Contrast.x),scales = "free") + scale_color_manual(values=c("grey",plotcols3[1])) + theme_classic() +
  scale_alpha_manual(values=c(0.4,0.8)) + scale_size_manual(values=c(1,2)) 

p2 <- dplyr::filter(res4,Contrast.x %in% c("Nrde2KO_vs_WT","Ccdc174KD_vs_WT","Mtr4KO_vs_WT")) %>%
  ggplot(aes(x=log2FoldChange.x,y=log.padj.x)) +geom_point(aes(color=Ccdc174.peaks,size=Ccdc174.peaks,alpha=Ccdc174.peaks)) + 
  facet_grid(vars(Contrast.x),scales = "free") + scale_color_manual(values=c(plotcols3[2],"grey")) + theme_classic() +
  scale_alpha_manual(values=c(0.8,0.4)) + scale_size_manual(values=c(2,1))

p3 <- dplyr::filter(res4,Contrast.x %in% c("Nrde2KO_vs_WT","Ccdc174KD_vs_WT","Mtr4KO_vs_WT")) %>%
  ggplot(aes(x=log2FoldChange.x,y=log.padj.x)) +geom_point(aes(color=Mtr4.peaks,size=Mtr4.peaks,alpha=Mtr4.peaks)) + 
  facet_grid(vars(Contrast.x),scales = "free") + scale_color_manual(values=c(plotcols3[3],"grey")) + theme_classic() +
  scale_alpha_manual(values=c(0.8,0.4)) + scale_size_manual(values=c(2,1))

p4 <- dplyr::filter(res4,Contrast.x %in% c("Nrde2KO_vs_WT","Ccdc174KD_vs_WT","Mtr4KO_vs_WT")) %>%
  ggplot(aes(x=log2FoldChange.x,y=log.padj.x)) +geom_point(aes(color=Eif4a3.peaks,size=Eif4a3.peaks,alpha=Eif4a3.peaks)) + 
  facet_grid(vars(Contrast.x),scales = "free") + scale_color_manual(values=c(plotcols3[4],"grey")) + theme_classic() +
  scale_alpha_manual(values=c(0.8,0.4)) + scale_size_manual(values=c(2,1))
 
p1 + p2 + p3 + p4
ggsave("Figures/Volcano_plots_intronic_reads_colored_by_peaks_no_reps.png",device="png",height=10, width=14)

#----------------------------------------------------------------------------------------------------------------
#   # plot BCLIP FC at upregulated introns version 1
#-----------------------------------------------------------------------------------------------------------------

#get BCLIP counts at splice donor sites
donor_counts2_log <- read.table("BCLIP_counts_in100bp_5prime_splice_site_region_v3.txt",sep="\t",header=TRUE,row.names=1)
donor_counts2_log$ID <- row.names(donor_counts2_log)



#combine RNAseq with BCLIP data
res4 <- inner_join(res3,donor_counts2_log,by="ID")

#select the ones which do not have exons within them 
#up.introns <- introns3g[names(introns3g) %in% res2$ID[res2$Contrast=="Nrde2_KO_vs_Nrde2_WT" & res2$intron_regulated=="up"]]
#up.introns2 <- subsetByOverlaps(up.introns,exons3g,invert = TRUE,type="within")
#unique(res3$gene_symbol[res3$ID %in% names(up.introns2)])

######Nrde2 vs Eif4a3#####
#calculate all p-values 
pvalues <- compare_means(Nrde2vsEif4a3 ~ intron_gene_regulated, data = res4, group.by="Contrast.x")
#make density plots for all contrasts
plist <- list(NA,NA,NA,NA)
for (i in 1:length(contrasts)){
p <- ggplot(res4[res4$Contrast.x==contrasts[i],],aes(x=Nrde2vsEif4a3,group=intron_gene_regulated)) + geom_density(aes(col=intron_gene_regulated)) +
  theme_classic() + scale_color_manual(values=c("darkgrey","darkblue")) + 
  ggtitle(contrasts[i]) + labs(subtitle = sprintf("p=%s",pvalues$p.format[pvalues$Contrast.x == contrasts[i]]))
p.data <- ggplot_build(p)$data[[1]]
p.text <- lapply(split(p.data, f = p.data$group), function(df){
  df[which.max(df$scaled), ]
})
p.text <- do.call(rbind, p.text)
if(nrow(p.text) > 1){
if(abs(p.text$x[1]-p.text$x[2]) < 1){
p.text$x[1] <- p.text$x[1] -1
p.text$x[2] <- p.text$x[2] +1
}}
plist[[i]] <- p + annotate('text', x = p.text$x, y = p.text$y,
              label = sprintf('n = %d', p.text$n), vjust = 0)
}

plist[[1]] + plist[[2]] + plist[[3]] + plist[[4]] + plist[[5]] + plist[[6]] 

ggsave("Figures/All_reads_Introns_up_gene_downFC0.5_pval0.1_Nrde2vsEIf4a3_density_plots3.png",width=20, height=10)


######Nrde2, Eif4a3,.....#####
#calculate all p-values 
pvalues <- compare_means(Nrde2 ~ intron_gene_regulated, data = res4, group.by="Contrast.x")
#make density plots for all contrasts
plist <- list(NA,NA,NA,NA)
for (i in 1:length(contrasts)){
  p <- ggplot(res4[res4$Contrast.x==contrasts[i],],aes(x=Nrde2,group=intron_gene_regulated)) + geom_density(aes(col=intron_gene_regulated)) +
    theme_classic() + scale_color_manual(values=c("darkgrey","darkblue")) + 
    ggtitle(contrasts[i]) + labs(subtitle = sprintf("p=%s",pvalues$p.format[pvalues$Contrast.x == contrasts[i]]))
  p.data <- ggplot_build(p)$data[[1]]
  p.text <- lapply(split(p.data, f = p.data$group), function(df){
    df[which.max(df$scaled), ]
  })
  p.text <- do.call(rbind, p.text)
  if(nrow(p.text) > 1){
    if(abs(p.text$x[1]-p.text$x[2]) < 1){
      p.text$x[1] <- p.text$x[1] -1
      p.text$x[2] <- p.text$x[2] +1
    }}
  plist[[i]] <- p + annotate('text', x = p.text$x, y = p.text$y,
                             label = sprintf('n = %d', p.text$n), vjust = 0)
}

plist[[1]] + plist[[2]] + plist[[3]] + plist[[4]] + plist[[5]] + plist[[6]] 
ggsave("Figures/All_reads_Introns_up_gene_downFC0.5_pval0.1_Nrde2_density_plots3.png",width=20, height=10)

#----------------------------------------------------------------------------------------------------------------
#   #compare all intron spliced read regulation data by heatmap
#-----------------------------------------------------------------------------------------------------------------

#find all introns that are regulated at all in any contarst
regIDs <- unique(res3$ID[res3$intron_regulated=="up" | res3$intron_regulated=="down"])
#OR: find all introns that are upregulated at all in any contarst and downregulated on the gene level
regIDs <- unique(res3$ID[res3$intron_gene_regulated=="Intron up, gene down"])

res3reg <- res3[res3$ID %in%regIDs,]
#make wide matrix of fold changes
res3regFCs <- res3reg %>% dplyr::select(ID,reps,log2FoldChange.x,Contrast.x) %>%
  pivot_wider(id_cols = c("ID","reps"),values_from=log2FoldChange.x,names_from=Contrast.x)
res3regFCs.m <- as.matrix(res3regFCs[,3:ncol(res3regFCs)])
#row.names(res3regFCs.m) <- res3regFCs$ID
reps.FCs <- res3regFCs$reps

#heatmap of fold chanegs for all contrasts, clustered
pdf("plots/heatmap_of_FCs_spliced_intronic_reads_all_upregulated_introns_splitbyreps.pdf",height=10,width=20)
draw(Heatmap(res3regFCs.m),split=reps.FCs,)
dev.off()

#----------------------------------------------------------------------------------------------------------------
#   # plot BCLIP FC at upregulated introns version 2
#-----------------------------------------------------------------------------------------------------------------

#get BCLIP counts at splice donor sites
donor_counts2_log <- read.table("BCLIP_counts_in100bp_5prime_splice_site_region_v3.txt",sep="\t",header=TRUE,row.names=1)
donor_counts2_log$ID <- row.names(donor_counts2_log)



#combine RNAseq with BCLIP data
res4 <- inner_join(res3,donor_counts2_log,by="ID")
#filter out introns containing repeats
res4 <- dplyr::filter(res4,reps=="no repeats")


#make it long and skinny to plot all violin plots
res5 <- pivot_longer(res4,cols=c("Nrde2","Ccdc174","Mtr4","Eif4a3"),names_to = "BCLIP",values_to = "cpm")
res5 <- dplyr::filter(res5,Contrast.x %in% c("Nrde2KO_vs_WT","Ccdc174KD_vs_WT","Mtr4KO_vs_WT"))
res5$intron_gene_upregulated <- ifelse(res5$intron_gene_regulated=="Intron up, gene down","up","no")
res5$intron_gene_upregulated <- replace_na(res5$intron_gene_upregulated,replace = "no")

ggplot(res5,aes(x=BCLIP,y=cpm,fill=intron_gene_upregulated)) + geom_boxplot() + facet_grid(vars(Contrast.x)) +
  theme_classic() + scale_fill_manual(values=c("grey","red"))
ggsave("Figures/All_reads_Introns_up_gene_downFC0.5_pval0.1_Nrde2_Eif4a3_Ccdc174_Mtr4_BCLIP_boxplots.pdf",device="pdf",width=10, height=10)



