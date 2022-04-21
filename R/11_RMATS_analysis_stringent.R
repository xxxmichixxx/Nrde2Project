library(tidyverse)


#----------------------------------------------------------------------------------------------------------------
#   #prepare bam file lists for each sample for RMATS - single end data-
#----------------------------------------------------------------------------------------------------------------
setwd("/tungstenfs/nobackup/gbuehler/michi/Projects/Nrde2/RNAseq/")

#generate RMATS sample files
all.bamFiles <- list.files("/tungstenfs/scratch/gbuehler/deepSeqRepos/bam/",full.names = TRUE,pattern="bam$")               
bamFiles <- grep("2191F",all.bamFiles,value =TRUE)
bamFiles <- c(grep("spliced",bamFiles,value =TRUE)[1:14],grep("2447F",all.bamFiles,value =TRUE))
bamFiles.WT <- grep("WT",bamFiles,value =TRUE)
bamFiles.WT <- grep("CHX",bamFiles.WT,value =TRUE,invert=TRUE)
bamFiles.Nrde2KO <- c(grep("Nrde2KO",bamFiles,value =TRUE),grep("Nrde2-KO",bamFiles,value =TRUE))
bamFiles.Nrde2KO <- grep("CHX",bamFiles.Nrde2KO,value =TRUE,invert=TRUE)
bamFiles.Mtr4KO <- c(grep("Mtr4_KO",bamFiles,value =TRUE))
bamFiles.Ccdc174KO <- c(grep("Ccdc174KD",bamFiles,value =TRUE))
bamFiles.Ccdc174KO <- grep("CHX",bamFiles.Ccdc174KO,value =TRUE,invert=TRUE)
bamFiles.Nrde2d200 <- c(grep("Nrde2-d200",bamFiles,value =TRUE))
bamFiles.Nrde2D174R <- c(grep("Nrde2-D174R",bamFiles,value =TRUE))

write.table(data.frame(matrix(nrow=1,ncol=length(bamFiles.WT),data=bamFiles.WT)),
            file="../rmats_analysis/bamFiles.WT.txt",quote = FALSE, append = FALSE,
            col.names = FALSE, row.names=FALSE, sep=",")
write.table(data.frame(matrix(nrow=1,ncol=length(bamFiles.Nrde2KO),data=bamFiles.Nrde2KO)),
            file="../rmats_analysis/bamFiles.Nrde2KO.txt",quote = FALSE, append = FALSE,
            col.names = FALSE, row.names=FALSE, sep=",")
write.table(data.frame(matrix(nrow=1,ncol=length(bamFiles.Mtr4KO),data=bamFiles.Mtr4KO)),
            file="../rmats_analysis/bamFiles.Mtr4KO.txt",quote = FALSE, append = FALSE,
            col.names = FALSE, row.names=FALSE, sep=",")
write.table(data.frame(matrix(nrow=1,ncol=length(bamFiles.Ccdc174KO),data=bamFiles.Ccdc174KO)),
            file="../rmats_analysis/bamFiles.Ccdc174KO.txt",quote = FALSE, append = FALSE,
            col.names = FALSE, row.names=FALSE, sep=",")
write.table(data.frame(matrix(nrow=1,ncol=length(bamFiles.Nrde2d200),data=bamFiles.Nrde2d200)),
            file="../rmats_analysis/bamFiles.Nrde2d200.txt",quote = FALSE, append = FALSE,
            col.names = FALSE, row.names=FALSE, sep=",")
write.table(data.frame(matrix(nrow=1,ncol=length(bamFiles.Nrde2D174R),data=bamFiles.Nrde2D174R)),
            file="../rmats_analysis/bamFiles.Nrde2D174R.txt",quote = FALSE, append = FALSE,
            col.names = FALSE, row.names=FALSE, sep=",")

bamFiles.WT <- grep("WT",bamFiles,value =TRUE)
bamFiles.WT <- grep("CHX",bamFiles.WT,value =TRUE,invert=FALSE)
bamFiles.Nrde2KO <- c(grep("Nrde2KO",bamFiles,value =TRUE),grep("Nrde2-KO",bamFiles,value =TRUE))
bamFiles.Nrde2KO <- grep("CHX",bamFiles.Nrde2KO,value =TRUE,invert=FALSE)
bamFiles.Ccdc174KO <- c(grep("Ccdc174KD",bamFiles,value =TRUE))
bamFiles.Ccdc174KO <- grep("CHX",bamFiles.Ccdc174KO,value =TRUE,invert=FALSE)

write.table(data.frame(matrix(nrow=1,ncol=length(bamFiles.WT),data=bamFiles.WT)),
            file="../rmats_analysis/bamFiles.WT_CHX.txt",quote = FALSE, append = FALSE,
            col.names = FALSE, row.names=FALSE, sep=",")
write.table(data.frame(matrix(nrow=1,ncol=length(bamFiles.Nrde2KO),data=bamFiles.Nrde2KO)),
            file="../rmats_analysis/bamFiles.Nrde2KO_CHX.txt",quote = FALSE, append = FALSE,
            col.names = FALSE, row.names=FALSE, sep=",")
write.table(data.frame(matrix(nrow=1,ncol=length(bamFiles.Ccdc174KO),data=bamFiles.Ccdc174KO)),
            file="../rmats_analysis/bamFiles.Ccdc174KO_CHX.txt",quote = FALSE, append = FALSE,
            col.names = FALSE, row.names=FALSE, sep=",")


#----------------------------------------------------------------------------------------------------------------
#   #prepare bam file lists for each sample for RMATS - paired end data-
#----------------------------------------------------------------------------------------------------------------
setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/RNAseq/")

#generate RMATS sample files
all.bamFiles <- list.files("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/RNAseq/hisat2_lib643_dta/",full.names = TRUE,pattern="bam$")               

bamFiles.WT <- grep("0h",all.bamFiles,value =TRUE)
bamFiles.Nrde2KO <- c(grep("6h",all.bamFiles,value =TRUE))

write.table(data.frame(matrix(nrow=1,ncol=length(bamFiles.WT),data=bamFiles.WT)),
            file="../rmats_analysis/bamFiles.WT_PE.txt",quote = FALSE, append = FALSE,
            col.names = FALSE, row.names=FALSE, sep=",")
write.table(data.frame(matrix(nrow=1,ncol=length(bamFiles.Nrde2KO),data=bamFiles.Nrde2KO)),
            file="../rmats_analysis/bamFiles.Nrde2KO_PE.txt",quote = FALSE, append = FALSE,
            col.names = FALSE, row.names=FALSE, sep=",")

#----------------------------------------------------------------------------------------------------------------
#   #run RMATS to find differentially expressed exons and introns based on stringtie gtf
#----------------------------------------------------------------------------------------------------------------
#on xenon 6
cd /tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/rmats_analysis
#install rmats vioa conda
#conda create -n rmats -c bioconda -c conda-forge rmats=4.1.0
#activate rmats conda environment
conda activate  rmats
#rmats.py --help

#run rmats on stringtie gtf
rmats.py --gtf ../RNAseq/stringtie_all_lib643_and_lib2191and2447_groupwise_transcripts.annotated.gtf --b1 bamFiles.Nrde2KO.txt --b2 bamFiles.WT.txt --od Nrde2KO_vs_WT --tmp Nrde2KO_vs_WT.tmp -t single --readLength 51 --libType fr-secondstrand --nthread 6 
rmats.py --gtf ../RNAseq/stringtie_all_lib643_and_lib2191and2447_groupwise_transcripts.annotated.gtf --b1 bamFiles.Nrde2KO_CHX.txt --b2 bamFiles.WT_CHX.txt --od Nrde2KO_vs_WT_CHX --tmp Nrde2KO_vs_WT_CHX.tmp -t single --readLength 51 --libType fr-secondstrand --nthread 6 
rmats.py --gtf ../RNAseq/stringtie_all_lib643_and_lib2191and2447_groupwise_transcripts.annotated.gtf --b1 bamFiles.Mtr4KO.txt --b2 bamFiles.WT.txt --od Mtr4KO_vs_WT --tmp Mtr4KO_vs_WT.tmp -t single --readLength 51 --libType fr-secondstrand --nthread 6 
rmats.py --gtf ../RNAseq/stringtie_all_lib643_and_lib2191and2447_groupwise_transcripts.annotated.gtf --b1 bamFiles.Ccdc174KO.txt --b2 bamFiles.WT.txt --od Ccdc174KO_vs_WT --tmp Ccdc174KO_vs_WT.tmp -t single --readLength 51 --libType fr-secondstrand --nthread 6 
rmats.py --gtf ../RNAseq/stringtie_all_lib643_and_lib2191and2447_groupwise_transcripts.annotated.gtf --b1 bamFiles.Ccdc174KO_CHX.txt --b2 bamFiles.WT_CHX.txt --od Ccdc174KO_vs_WT_CHX --tmp Ccdc174KO_vs_WT_CHX.tmp -t single --readLength 51 --libType fr-secondstrand --nthread 6 
rmats.py --gtf ../RNAseq/stringtie_all_lib643_and_lib2191and2447_groupwise_transcripts.annotated.gtf --b1 bamFiles.WT_CHX.txt --b2 bamFiles.WT.txt --od WT_CHX_vs_WT --tmp WT_CHX_vs_WT.tmp -t single --readLength 51 --libType fr-secondstrand --nthread 6 
rmats.py --gtf ../RNAseq/stringtie_all_lib643_and_lib2191and2447_groupwise_transcripts.annotated.gtf --b1 bamFiles.Nrde2d200.txt --b2 bamFiles.WT.txt --od Nrde2d200_vs_WT --tmp Nrde2d200_vs_WT.tmp -t single --readLength 51 --libType fr-secondstrand --nthread 6 
rmats.py --gtf ../RNAseq/stringtie_all_lib643_and_lib2191and2447_groupwise_transcripts.annotated.gtf --b1 bamFiles.Nrde2D174R.txt --b2 bamFiles.WT.txt --od Nrde2D174R_vs_WT --tmp Nrde2D174R_vs_WT.tmp -t single --readLength 51 --libType fr-secondstrand --nthread 6 

#rmats.py --gtf ../RNAseq/stringtie_all_lib643_and_lib2191and2447_groupwise_transcripts.annotated.gtf --b1 bamFiles.Nrde2KO_PE.txt --b2 bamFiles.WT_PE.txt --od Nrde2KO_vs_WT_PE --tmp Nrde2KO_vs_WT_PE.tmp -t paired --readLength 100 --libType fr-firststrand --nthread 6 


#----------------------------------------------------------------------------------------------------------------
#   #analyze RMATS results: intron retention
#----------------------------------------------------------------------------------------------------------------
setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")
list.files("rmats_analysis/Nrde2KO_vs_WT")

RI.MATS.JCEC <- read.table("rmats_analysis/Nrde2KO_vs_WT_CHX//RI.MATS.JCEC.txt",header=TRUE)
RI.MATS.JCEC$logFDR <- -log10(RI.MATS.JCEC$FDR)
RI.MATS.JCEC$logPVAL <- -log10(RI.MATS.JCEC$PValue)

ggplot(RI.MATS.JCEC,aes(x=IncLevelDifference,y=logFDR)) + geom_point()
ggplot(RI.MATS.JCEC,aes(x=IncLevelDifference,y=logPVAL)) + geom_point()

RI.MATS.JCEC %>% dplyr::filter(FDR < 0.05, IncLevelDifference > 0.1)
RI.MATS.JCEC %>% dplyr::filter(geneSymbol=="Cdk2")
RI.MATS.JCEC %>% dplyr::filter(GeneID=="MSTRG.5268")


#----------------------------------------------------------------------------------------------------------------
#   #analyze RMATS results: new 5'SS
#----------------------------------------------------------------------------------------------------------------
setwd("/tungstenfs/scratch/gbuehler/michi/Projects/Nrde2/")
list.files("rmats_analysis/Nrde2KO_vs_WT")

A5SS.MATS.JCEC <- read.table("rmats_analysis/Nrde2KO_vs_WT_CHX//A5SS.MATS.JCEC.txt",header=TRUE)
A5SS.MATS.JCEC$logFDR <- -log10(A5SS.MATS.JCEC$FDR)
A5SS.MATS.JCEC$logPVAL <- -log10(A5SS.MATS.JCEC$PValue)

ggplot(A5SS.MATS.JCEC,aes(x=IncLevelDifference,y=logFDR)) + geom_point()
ggplot(A5SS.MATS.JCEC,aes(x=IncLevelDifference,y=logPVAL)) + geom_point()

A5SS.MATS.JCEC %>% dplyr::filter(FDR < 0.01, IncLevelDifference > 0)
A5SS.MATS.JCEC %>% dplyr::filter(geneSymbol=="Cdk2")
A5SS.MATS.JCEC.Cdk2 <- A5SS.MATS.JCEC %>% dplyr::filter(GeneID=="MSTRG.5268")
A5SS.MATS.JCEC %>% dplyr::filter(GeneID=="MSTRG.40736")

write.table(A5SS.MATS.JCEC.Cdk2,file="rmats_analysis/A5SS.MATS.JCEC.Cdk2.txt",col.names=TRUE,row.names=FALSE,sep="\t",append=FALSE,quote=FALSE)

#----------------------------------------------------------------------------------------------------------------
#   #analyze RMATS results: all samples
#----------------------------------------------------------------------------------------------------------------
setwd("/tungstenfs/nobackup/gbuehler/michi/Projects/Nrde2/")
list.dirs("rmats_analysis/")

comp <- c("Nrde2KO_vs_WT","Nrde2d200_vs_WT","Nrde2D174R_vs_WT","Nrde2KO_vs_WT_CHX","Ccdc174KO_vs_WT","Ccdc174KO_vs_WT_CHX","Mtr4KO_vs_WT","WT_CHX_vs_WT")
RIlist <- list()
A5SSlist <- list()
A3SSlist <- list()

for (i in seq_along(comp)){
  #retained introns
  RI.MATS.JCEC <- read.table(sprintf("rmats_analysis/%s/RI.MATS.JCEC.txt",comp[i]),header=TRUE)
  RI.MATS.JCECsig <- RI.MATS.JCEC %>% dplyr::filter(FDR < 0.01,IncLevelDifference > 0.1)
  #calculate total number of supporting reads across samples and filter
  RI.MATS.JCECsig$numReads <- apply(RI.MATS.JCECsig[,13:16],1,function(x){sum(as.numeric(unlist(strsplit(x,","))))})
  RI.MATS.JCECsig <- RI.MATS.JCECsig %>% dplyr::filter(numReads > 30)
#  RI.MATS.JCECsig <- RI.MATS.JCEC %>% dplyr::filter(FDR < 0.01)
  RI.MATS.JCECsig$comp <- comp[i]
  RIlist[[i]] <- RI.MATS.JCECsig
 
   #5'ASS
  A5SS.MATS.JCEC <- read.table(sprintf("rmats_analysis/%s/A5SS.MATS.JCEC.txt",comp[i]),header=TRUE)
  A5SS.MATS.JCECsig <- A5SS.MATS.JCEC %>% dplyr::filter(FDR < 0.01,IncLevelDifference > 0.1)
  #calculate total number of supporting reads across samples and filter
  A5SS.MATS.JCECsig$numReads <- apply(A5SS.MATS.JCECsig[,13:16],1,function(x){sum(as.numeric(unlist(strsplit(x,","))))})
  A5SS.MATS.JCECsig <- A5SS.MATS.JCECsig %>% dplyr::filter(numReads > 30)
 # A5SS.MATS.JCECsig <- A5SS.MATS.JCEC %>% dplyr::filter(FDR < 0.01)
  A5SS.MATS.JCECsig$comp <- comp[i]
  A5SSlist[[i]] <- A5SS.MATS.JCECsig

    #3'ASS
  A3SS.MATS.JCEC <- read.table(sprintf("rmats_analysis/%s/A3SS.MATS.JCEC.txt",comp[i]),header=TRUE)
  A3SS.MATS.JCECsig <- A3SS.MATS.JCEC %>% dplyr::filter(FDR < 0.01,IncLevelDifference > 0.1)
  #calculate total number of supporting reads across samples and filter
  A3SS.MATS.JCECsig$numReads <- apply(A3SS.MATS.JCECsig[,13:16],1,function(x){sum(as.numeric(unlist(strsplit(x,","))))})
  A3SS.MATS.JCECsig <- A3SS.MATS.JCECsig %>% dplyr::filter(numReads > 30)
 # A3SS.MATS.JCECsig <- A3SS.MATS.JCEC %>% dplyr::filter(FDR < 0.01)
  A3SS.MATS.JCECsig$comp <- comp[i]
  A3SSlist[[i]] <- A3SS.MATS.JCECsig
}

RIres <- do.call(rbind,RIlist)
A5SSres <- do.call(rbind,A5SSlist)
A3SSres <- do.call(rbind,A3SSlist)

ASS.counts <- data.frame(RI=table(RIres$comp),
A5SS=table(A5SSres$comp),
A3SS=table(A3SSres$comp)) %>% pivot_longer(cols = c("RI.Freq","A5SS.Freq","A3SS.Freq"))

ggplot(ASS.counts,aes(x=RI.Var1, y=value,fill=name)) + geom_bar(stat="identity",position = "dodge") + theme_classic()
ggsave("Figures/RMATS_sig_up_counts2_inclusion_level_diff_0.1_min30reads.png",width=14,height=6)

#----------------------------------------------------------------------------------------------------------------
#   #compare RMATS results to BCLIP peaks
#----------------------------------------------------------------------------------------------------------------
print(load("clipper_peaks_5readsmin_200202.RData"))

txdb <- makeTxDbFromGFF("RNAseq/stringtie_all_lib643_and_lib2191and2447_groupwise_transcripts.annotated.gtf","gtf")    
exons <- unlist(exonsBy(txdb,by="gene",use.names=FALSE))
exons$Nrde2.peak <- overlapsAny(exons,clipper.peaks[clipper.peaks$sample=="Nrde2"])
exons$Eif4a3.peak <- overlapsAny(exons,clipper.peaks[clipper.peaks$sample=="Eif4a3"])

table(exons$Nrde2.peak)
table(exons$Eif4a3.peak)


#select the peaks that overlap an exon upstream of an intron
introns <- unique(unlist(intronsByTranscript(txdb,use.names=TRUE)))
SS5 <- promoters(introns,upstream = 100,downstream = 1)
Nrde2.peaks <- subsetByOverlaps(clipper.peaks[clipper.peaks$sample=="Nrde2"],SS5)
Ccdc174.peaks <- subsetByOverlaps(clipper.peaks[clipper.peaks$sample=="Ccdc174"],SS5)
Eif4a3.peaks <- subsetByOverlaps(clipper.peaks[clipper.peaks$sample=="Eif4a3"],SS5)


#5' SS
#take the long exon for mking a gRanges object
up.peaks <- matrix(ncol=3,nrow=length(comp))
for (i in seq_along(comp)){

A5SS.MATS.JCEC <- read.table(sprintf("rmats_analysis/%s/A5SS.MATS.JCEC.txt",comp[i]),header=TRUE)
A5SS.MATS.JCEC$numReads <- apply(A5SS.MATS.JCEC[,13:16],1,function(x){sum(as.numeric(unlist(strsplit(x,","))))})

A5SS.gr <- makeGRangesFromDataFrame(A5SS.MATS.JCEC,
                         keep.extra.columns=TRUE,
                         seqinfo=NULL,
                         seqnames.field=c("chr"),
                         start.field=c("longExonStart_0base"),
                         end.field=c("longExonEnd"),
                         strand.field=c("strand"),
                         ignore.strand=FALSE,
                         starts.in.df.are.0based=FALSE)


#check if overlaps with peak , yes/no
A5SS.gr$Nrde2.peak <- overlapsAny(A5SS.gr,Nrde2.peaks)
A5SS.gr$Ccdc174.peak <- overlapsAny(A5SS.gr,Ccdc174.peaks)
A5SS.gr$Eif4a3.peak <- overlapsAny(A5SS.gr,Eif4a3.peaks)

A5SS.gr$upregulated <- ifelse(A5SS.gr$FDR < 0.01 & A5SS.gr$IncLevelDifference > 0.1 & A5SS.gr$numReads > 30, "yes","no")


plot(density(A5SS.gr$IncLevelDifference[A5SS.gr$Nrde2.peak==TRUE]),col="red",xlim=c(-0.2,0.2))
lines(density(A5SS.gr$IncLevelDifference[A5SS.gr$Nrde2.peak==FALSE]))

table(A5SS.gr$Nrde2.peak[A5SS.gr$upregulated=="yes"])
table(A5SS.gr$Nrde2.peak[A5SS.gr$upregulated=="no"])

up.peaks[i,1] <- table(A5SS.gr$upregulated[A5SS.gr$Nrde2.peak==TRUE])[2]
up.peaks[i,2] <-table(A5SS.gr$upregulated[A5SS.gr$Ccdc174.peak==TRUE])[2]
up.peaks[i,3] <-table(A5SS.gr$upregulated[A5SS.gr$Eif4a3.peak==TRUE])[2]
}
colnames(up.peaks) <- c("Nrde2",'Ccdc174',"Eif4a3")
row.names(up.peaks) <- comp

up.peaks.perc <- data.frame(up.peaks)
up.peaks.perc$Nrde2 <- up.peaks.perc$Nrde2/length(Nrde2.peaks)*100
up.peaks.perc$Ccdc174 <- up.peaks.perc$Ccdc174/length(Ccdc174.peaks)*100
up.peaks.perc$Eif4a3 <- up.peaks.perc$Eif4a3/length(Eif4a3.peaks)*100
up.peaks.perc$comp <- row.names(up.peaks.perc)
up.peaks.perc2 <- pivot_longer(up.peaks.perc,cols = c("Nrde2",'Ccdc174',"Eif4a3"),names_to = "peak",values_to = "percent")
  
  
ggplot(up.peaks.perc2,aes(x=comp,y=percent,fill=peak)) + geom_bar(stat="identity",position="dodge") + theme_classic()
ggsave("Figures/RMATS_A5SS_sig_up_percent_of_peaks_inclLevel0.1_min30reads.png",width=14,height=6)

#3' SS
#take the long exon for mking a gRanges object
up.peaks <- matrix(ncol=3,nrow=length(comp))
for (i in seq_along(comp)){
  
  A5SS.MATS.JCEC <- read.table(sprintf("rmats_analysis/%s/A3SS.MATS.JCEC.txt",comp[i]),header=TRUE)
  A5SS.MATS.JCEC$numReads <- apply(A5SS.MATS.JCEC[,13:16],1,function(x){sum(as.numeric(unlist(strsplit(x,","))))})
  A5SS.gr <- makeGRangesFromDataFrame(A5SS.MATS.JCEC,
                                      keep.extra.columns=TRUE,
                                      seqinfo=NULL,
                                      seqnames.field=c("chr"),
                                      start.field=c("flankingES"),
                                      end.field=c("flankingEE"),
                                      strand.field=c("strand"),
                                      ignore.strand=FALSE,
                                      starts.in.df.are.0based=FALSE)
  
  
  #check if overlaps with peak , yes/no
  A5SS.gr$Nrde2.peak <- overlapsAny(A5SS.gr,Nrde2.peaks)
  A5SS.gr$Ccdc174.peak <- overlapsAny(A5SS.gr,Ccdc174.peaks)
  A5SS.gr$Eif4a3.peak <- overlapsAny(A5SS.gr,Eif4a3.peaks)
  
  A5SS.gr$upregulated <- ifelse(A5SS.gr$FDR < 0.01 & A5SS.gr$IncLevelDifference > 0 & A5SS.gr$numReads > 30, "yes","no")
  
  
  plot(density(A5SS.gr$IncLevelDifference[A5SS.gr$Nrde2.peak==TRUE]),col="red",xlim=c(-0.2,0.2))
  lines(density(A5SS.gr$IncLevelDifference[A5SS.gr$Nrde2.peak==FALSE]))
  
  table(A5SS.gr$Nrde2.peak[A5SS.gr$upregulated=="yes"])
  table(A5SS.gr$Nrde2.peak[A5SS.gr$upregulated=="no"])
  
  up.peaks[i,1] <- table(A5SS.gr$upregulated[A5SS.gr$Nrde2.peak==TRUE])[2]
  up.peaks[i,2] <-table(A5SS.gr$upregulated[A5SS.gr$Ccdc174.peak==TRUE])[2]
  up.peaks[i,3] <-table(A5SS.gr$upregulated[A5SS.gr$Eif4a3.peak==TRUE])[2]
}
colnames(up.peaks) <- c("Nrde2",'Ccdc174',"Eif4a3")
row.names(up.peaks) <- comp

up.peaks.perc <- data.frame(up.peaks)
up.peaks.perc$Nrde2 <- up.peaks.perc$Nrde2/length(Nrde2.peaks)*100
up.peaks.perc$Ccdc174 <- up.peaks.perc$Ccdc174/length(Ccdc174.peaks)*100
up.peaks.perc$Eif4a3 <- up.peaks.perc$Eif4a3/length(Eif4a3.peaks)*100
up.peaks.perc$comp <- row.names(up.peaks.perc)
up.peaks.perc2 <- pivot_longer(up.peaks.perc,cols = c("Nrde2",'Ccdc174',"Eif4a3"),names_to = "peak",values_to = "percent")


ggplot(up.peaks.perc2,aes(x=comp,y=percent,fill=peak)) + geom_bar(stat="identity",position="dodge") + theme_classic()
ggsave("Figures/RMATS_A3SS_sig_up_percent_of_peaks_inclLevel0.1_min30reads.png",width=14,height=6)


