
rm(list=ls())


library(tidyverse)
library(GenomicFeatures)
library(GenomicRanges)




setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")



#----------------------------------------------------------------------------------------------------------------
#   #load exon junctions and introns
#----------------------------------------------------------------------------------------------------------------

#print(load("splice_sites_expressed_and_notexpressed_transcripts_log2TPMabove0.1_and_genes_FPKMabove8below0.RData"))
print(load("WT_strtie_gtf_splice_sites_expressed_and_notexpressed_transcripts_log2TPMabove0.1_and_genes_FPKMabove8below0.RData"))

names(spliceDonors) <- paste(names(spliceDonors),start(spliceDonors),sep="_")
names(introns3g) <- names(spliceDonors)
names(spliceAcceptors) <- paste(names(spliceAcceptors),start(spliceAcceptors),sep="_")

spliceDonors2 <- GenomicRanges::promoters(spliceDonors,upstream=100,downstream=1)

#----------------------------------------------------------------------------------------------------------------
#   #which 5' splice sites have Nrde2 at 5' splice site ####
#----------------------------------------------------------------------------------------------------------------
#load peaks
print(load("clipper_peaks_5readsmin_200202.RData"))
clipper.peaks$ID <- names(clipper.peaks)
samples <- unique(clipper.peaks$sample)

#Nrde2
peaks <- clipper.peaks[clipper.peaks$sample==samples[4]]
spliceDonors$Nrde2.peak <- ifelse(overlapsAny(spliceDonors2,peaks,
                                              minoverlap = 10),1 ,0)
table(spliceDonors$Nrde2.peak)

#Ccdc174
peaks <- clipper.peaks[clipper.peaks$sample==samples[1]]
spliceDonors$Ccdc174.peak <- ifelse(overlapsAny(spliceDonors2,peaks,
                                                minoverlap = 10),1 ,0)
table(spliceDonors$Ccdc174.peak)

#Eif4a3
peaks <- clipper.peaks[clipper.peaks$sample==samples[2]]
spliceDonors$Eif4a3.peak <- ifelse(overlapsAny(spliceDonors2,peaks,
                                               minoverlap = 10),1 ,0)
table(spliceDonors$Eif4a3.peak)

#Nrde2d200
peaks <- clipper.peaks[clipper.peaks$sample==samples[6]]
spliceDonors$Nrde2d200.peak <- ifelse(overlapsAny(spliceDonors2,peaks,
                                              minoverlap = 10),1 ,0)
table(spliceDonors$Nrde2d200.peak)

#Nrde2D174R
peaks <- clipper.peaks[clipper.peaks$sample==samples[5]]
spliceDonors$Nrde2D174R.peak <- ifelse(overlapsAny(spliceDonors2,peaks,
                                                  minoverlap = 10),1 ,0)
table(spliceDonors$Nrde2D174R.peak)

#Mtr4
peaks <- clipper.peaks[clipper.peaks$sample==samples[3]]
spliceDonors$Mtr4.peak <- ifelse(overlapsAny(spliceDonors2,peaks,
                                                   minoverlap = 10),1 ,0)
table(spliceDonors$Mtr4.peak)

#----------------------------------------------------------------------------------------------------------------
#   #find the number of 5'SS per unique 3'SS ####
#----------------------------------------------------------------------------------------------------------------

spliceDonors$SS5 <- paste(seqnames(spliceDonors),start(spliceDonors),sep="_")
spliceDonors$SS3 <- paste(seqnames(spliceAcceptors),start(spliceAcceptors),sep="_")
unique3SS <- unique(spliceDonors$SS3)

res <- data.frame(ncol=8,nrow=length(unique3SS))
for (i in seq_along(unique3SS)){
  ins <- spliceDonors[spliceDonors$SS3==unique3SS[i]]

  res[i,1] <- ins$GeneID[1]
  res[i,2] <- length(unique(ins$SS5))
  res[i,3] <- sum(ins$Nrde2.peak)
  res[i,4] <- sum(ins$Ccdc174.peak)
  res[i,5] <- sum(ins$Eif4a3.peak)
  res[i,6] <- sum(ins$Nrde2d200.peak)
  res[i,7] <- sum(ins$Nrde2D174R.peak)
  res[i,8] <- sum(ins$Mtr4.peak)
}
colnames(res) <- c("GeneID","n5SS","Nrde2.peak","Ccdc174.peak","Eif4a3.peak","Nrde2d200.peak","Nrde2D174R.peak","Mtr4.peak")
res2 <- rbind(data.frame(res[res$Nrde2.peak > 0,],peak="Nrde2"),
              data.frame(res[res$Ccdc174.peak > 0,],peak="Ccdc174"),
              data.frame(res[res$Eif4a3.peak > 0,],peak="Eif4a3"),
              data.frame(res[res$Nrde2d200.peak > 0,],peak="Nrde2d200"),
              data.frame(res[res$Nrde2D174R.peak > 0,],peak="Nrde2D174R"),
              data.frame(res[res$Mtr4.peak > 0,],peak="Mtr4")
                         )
res2$multi5SS <- ifelse(res2$n5SS>1,"yes","no")

res2 <- res2 %>% group_by(peak) %>% dplyr::count(multi5SS) %>%
  mutate(perc=(n/sum(n))*100) 

#calculate p-values
M <- res2 %>% pivot_wider(id_cols = c("peak"),values_from="n",names_from="multi5SS") %>%
  dplyr::filter(peak %in% c("Eif4a3","Nrde2"))
  ct.Nrde2 <- chisq.test(M[,2:3])
M <- res2 %>% pivot_wider(id_cols = c("peak"),values_from="n",names_from="multi5SS") %>%
    dplyr::filter(peak %in% c("Eif4a3","Ccdc174"))
  ct.Ccdc174 <- chisq.test(M[,2:3])

#plot
  ggplot(res2,aes(x=peak,y=perc,fill=multi5SS)) + geom_bar(stat="identity") + theme_classic() +
  scale_fill_manual(values=c("grey","red")) + scale_x_discrete(limits = c("Nrde2","Ccdc174","Eif4a3","Nrde2d200","Nrde2D174R","Mtr4")) +
    ggtitle(sprintf("Nrde2=%s",ct.Nrde2$p.value),subtitle=sprintf("Ccdc174=%s",ct.Ccdc174$p.value))
  ggsave("Figures/percent_of_multi5SS_per_unique_3SS_all_peaks_WTstrtieGTF.pdf",height=6,width=4)


#----------------------------------------------------------------------------------------------------------------
#   #find the number of 3'SS per unique 5'SS ####
#----------------------------------------------------------------------------------------------------------------

spliceDonors$SS5 <- paste(seqnames(spliceDonors),start(spliceDonors),sep="_")
spliceDonors$SS3 <- paste(seqnames(spliceAcceptors),start(spliceAcceptors),sep="_")
unique5SS <- unique(spliceDonors$SS5)

res3 <- data.frame(ncol=5,nrow=length(unique5SS))
for (i in seq_along(unique5SS)){
  ins <- spliceDonors[spliceDonors$SS5==unique5SS[i]]
  
  res3[i,1] <- ins$GeneID[1]
  res3[i,2] <- length(unique(ins$SS3))
  res3[i,3] <- sum(ins$Nrde2.peak)
  res3[i,4] <- sum(ins$Ccdc174.peak)
  res3[i,5] <- sum(ins$Eif4a3.peak)
  res3[i,6] <- sum(ins$Nrde2d200.peak)
  res3[i,7] <- sum(ins$Nrde2D174R.peak)
  res3[i,8] <- sum(ins$Mtr4.peak)
}
colnames(res3) <- c("GeneID","n3SS","Nrde2.peak","Ccdc174.peak","Eif4a3.peak","Nrde2d200.peak","Nrde2D174R.peak","Mtr4.peak")
res4 <- rbind(data.frame(res3[res3$Nrde2.peak > 0,],peak="Nrde2"),
              data.frame(res3[res3$Ccdc174.peak > 0,],peak="Ccdc174"),
              data.frame(res3[res3$Eif4a3.peak > 0,],peak="Eif4a3"),
              data.frame(res3[res3$Nrde2d200.peak > 0,],peak="Nrde2d200"),
              data.frame(res3[res3$Nrde2D174R.peak > 0,],peak="Nrde2D174R"),
              data.frame(res3[res3$Mtr4.peak > 0,],peak="Mtr4")
)
res4$multi3SS <- ifelse(res4$n3SS>1,"yes","no")

res4 <- res4 %>% group_by(peak) %>% dplyr::count(multi3SS) %>%
  mutate(perc=(n/sum(n))*100)
  
#calculate p-values
M <- res4 %>% pivot_wider(id_cols = c("peak"),values_from="n",names_from="multi3SS") %>%
  dplyr::filter(peak %in% c("Eif4a3","Nrde2"))
ct.Nrde2 <- chisq.test(M[,2:3])
M <- res4 %>% pivot_wider(id_cols = c("peak"),values_from="n",names_from="multi3SS") %>%
  dplyr::filter(peak %in% c("Eif4a3","Ccdc174"))
ct.Ccdc174 <- chisq.test(M[,2:3])

#plot
ggplot(res4,aes(x=peak,y=perc,fill=multi3SS)) + geom_bar(stat="identity") + theme_classic() +
  scale_fill_manual(values=c("grey","red")) + scale_x_discrete(limits = c("Nrde2","Ccdc174","Eif4a3")) +
  ggtitle(sprintf("Nrde2=%s",ct.Nrde2$p.value),subtitle=sprintf("Ccdc174=%s",ct.Ccdc174$p.value))
  
ggsave("Figures/percent_of_multi3SS_per_unique_5SS_all_peaks.pdf",height=6,width=4)

