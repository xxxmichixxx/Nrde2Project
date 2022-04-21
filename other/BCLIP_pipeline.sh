#!/bin/bash

# This script should be run in a separate directory
# The directory must contain the raw data, e.g. Raw_230718.fq, and a samples file, e.g. samples.list:
#
#     NNNNTAAGC_L5Aa     eIF4a3_2i_exp37
#     NNNNATTAGC_L5Ab    eIF4a3_2i_HAR30_exp37
#     NNNNGCGCAGC_L5Ac   Nrde2_2i_exp37
#
# Example usage:
#	nohup /tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/BCLIP_pipeline.sh -S samples.list -F data.fastq -l 51 -b /tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/newbarcodes_bclip.list &

# provide relevant paths
    #PATH=$PATH:/work2/gbuehler/tuckalex/bin
    PATH=$PATH:/tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/bin
    PATH=$PATH:/work/gbioinfo/linux/bin/ncbi_blast+/bin
    PATH=$PATH:/work/gbioinfo/linux/bin
   # PATH=$PATH:/work2/gbuehler/tuckalex/bin/bowtie2-2.3.4.1-linux-x86_64
    PATH=$PATH:/tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/bin/bowtie2-2.3.4.1-linux-x86_64
    export PATH


LINKER3="TGGAATTCTCGGGTGCCAAGGC"        # miRCat-33 linker; full length is TGGAATTCTCGGGTGCCAAGGC
FASTQ_OFFSET=33                         # set to 33 for Sanger, Illumina 1.8+; set to 64 for Solexa, Illumina 1.3+, and Illumina 1.5+;
                                        # -Q 33 use with Sanger or Illumina 1.8+ fastq. Example quality string: BC@FFFFFFHHHHJJGHIJJ################ <--- this is the current version
                                        # -Q 64 use with Illumina 1.3+ or 1.5+ fastq. Example quality string: hhhhhhhhddhhdhfddfeBBBBBBBBBBBBBBB
BARCODE_FILE="/tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/barcodes_bclip.list"      #the barcodes.list file is for CRAC data, the barcodes_bclip.list file is for bclip data

SAMPLE_FILE=""                          # tab separated file listing barcode and sample, e.g. NNNNTAAGC_L5Aa     eIF4a3_2i_exp37
RAW_DATA_FILE=""                        # raw data file, e.g. Raw_230718.fq
SPLIT_LINES_TOGGLE="Y"                  # Split lines for LC stripping into 10? If you enter a number, instead split into this many LINES (e.g. should enter a number like 10,000,000)
                                        # Default behaviour is to divide file into 10
READ_LENGTH=51				# How long are the reads? Important to decide which had adapters clipped
	
# parsing command line options
while getopts "L: Q: S: F: s: l: b: " OPTION ; do
     case $OPTION in
       L) LINKER3=$OPTARG ;;
       Q) FASTQ_OFFSET=$OPTARG ;;
       S) SAMPLE_FILE=$OPTARG ;;
       F) RAW_DATA_FILE=$OPTARG ;;
       s) SPLIT_LINES_TOGGLE=$OPTARG ;;
       l) READ_LENGTH=$OPTARG ;;
       b) BARCODE_FILE=$OPTARG ;;
       ?) echo "incorrect option" ; exit 0 ;;
     esac
done
shift $((OPTIND -1))                    # this makes sure that $@ does not re-use the options already handled by the getopts


## PREPROCESSING - ONLY TRIMS ADAPTERS FROM READS => READ LENGTH DISTINGUISHES WHETHER ADAPTER WAS MAPPED OR NOT


echo "Processing file: " $RAW_DATA_FILE

RAW_DATA_BASENAME=${RAW_DATA_FILE%%.*}

echo "Basename for this experiment is: " $RAW_DATA_BASENAME

echo "Barcodes and samples are: "

cat $SAMPLE_FILE | column -t | awk '{print "\t\t"$0}'

echo "Using Trimmomatic to detect adapters (7 or longer)"
echo "The alignment score is calculated as follows: +0.6 per match, -Q/10 per mismatch. By considering the quality of the base calls, mismatches caused by read errors have less impact. A 7 bp perfect match will score 4.2 (thus with min score set to 4, only adapters 7 or longer will be clipped)"
echo "Also require reads to be at least 27 long (later filter for read >=18, but this is just to exclude really short things), and have all bases with a quality score of at least 25"
echo "Collapse duplicates (identical reads)"
echo "Annotate read names where adapter was clipped (lets you clip more, later on, e.g. by quality, without using the adapter information e.g. >47_5_adapter"
 
java -jar /tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/bin/Trimmomatic-0.32/trimmomatic-0.32.jar SE -threads 7 -phred33 $RAW_DATA_FILE ${RAW_DATA_BASENAME}.adapterStripped.fq ILLUMINACLIP:/tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/References/Solexa_adapters_CRAConly.fa:2:30:4 MINLEN:27
 
cat ${RAW_DATA_BASENAME}.adapterStripped.fq | /work/gbioinfo/linux/bin/fastq_quality_filter -q 10 -p 100 -Q 33 | awk '(NR-2)%4==0{print}' | sort | uniq -c | sort -nr | awk '{printf "%s_%s\t%s\n", NR,$1,$2}' | awk -v read_length=$READ_LENGTH 'BEGIN{FS="\t"; OFS="\t"} length($2)<read_length{$1 = $1 "_adapter"} {print ">"$1; print $2}' > ${RAW_DATA_BASENAME}_comp.fasta


echo "Debarcoding"

module use /tungstenfs/groups/gbioinfo/Appz/easybuild/modules/all	# Now using modules to run pyCRAC
module use /tungstenfs/groups/gbioinfo/Appz/modules
module purge								# Unload all modules, to avoid conflicts/odd behaviour (e.g. pyCRAC module fails if R-BioC/current module loaded
module load pyCRAC/1.3.3						# Need for pyCRAC
module load Kent_tools/20190212-linux.x86_64				# Need for bed2bigwig
#module load R-BioC/current						# Loads the Rstudio version of R - but only load later as it conflicts with pyCRAC

#/tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/venv-pycrac/bin/pyBarcodeFilter.py -f ${RAW_DATA_BASENAME}_comp.fasta -b /tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/barcodes.list --file_type=fasta
#/tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/venv-pycrac/bin/pyBarcodeFilter.py -f ${RAW_DATA_BASENAME}_comp.fasta -b $BARCODE_FILE --file_type=fasta
pyBarcodeFilter.py -f ${RAW_DATA_BASENAME}_comp.fasta -b $BARCODE_FILE --file_type=fasta


echo "Renaming debarcoded files, split, then either run just prinseq, or LCstripper and prinseq. At this point, filter sequences out that are < 18 nt"

if [ "$SPLIT_LINES_TOGGLE" = "Y" ]; then
      echo "Will split files into 10"
      else
      SPLIT_LINES_LC=$SPLIT_LINES_TOGGLE
      echo "Will split files into" $SPLIT_LINES_LC "lines, as you specified on the command line"
fi

while read -r sample
do
      barcode=`echo ${sample//$'\t'*/}`
      sample_name=`echo ${sample//*$'\t'/}`
      echo mv ${RAW_DATA_BASENAME}_comp_${barcode}.fasta ${sample_name}.fasta
      mv ${RAW_DATA_BASENAME}_comp_${barcode}.fasta ${sample_name}.fasta
      if [ "$SPLIT_LINES_TOGGLE" = "Y" ]; then
            SPLIT_LINES_LC=`wc -l ${sample_name}.fasta | awk '{printf "%d", $1/20}' | awk '{printf "%d", 4+2*$1}'`
            echo "Splitting file into 10, so" $SPLIT_LINES_LC "lines"
      else
            echo "Splitting file into" $SPLIT_LINES_LC "lines, as you specified on the command line"
      fi

      split -l $SPLIT_LINES_LC ${sample_name}.fasta split.${sample_name}.fasta.
      rm ${sample_name}.fasta
	ls split.${sample_name}.fasta.* | parallel '/tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/bin/prinseq-lite-0.20.4/prinseq-lite -lc_threshold 20 -lc_method dust -min_len 18 -fasta "{}" --out_good prinseq."{}" --out_bad rejected."{}"'
	ls split.${sample_name}.fasta.* | parallel '/tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/bin/Low_complexity_stripper.awk MIN_FINAL_LENGTH=18 THRESHOLD=80 MIN_TRIM=2 "{}" > LCstripped."{}"'
	ls LCstripped.split.${sample_name}.fasta.* | parallel '/tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/bin/prinseq-lite-0.20.4/prinseq-lite -lc_threshold 20 -lc_method dust -fasta "{}" --out_good prinseq."{}" --out_bad rejected."{}"'
	cat prinseq.split.${sample_name}.fasta.* > ${sample_name}.prinseq.fasta
	cat prinseq.LCstripped.split.${sample_name}.fasta.* > ${sample_name}.LCstripped.prinseq.fasta

done < $SAMPLE_FILE





echo "Taking adapterStripped file, quality filtering (as before), running prinseq LC filter, but this time not collapsing duplicate reads or running my low complexity stripper"
echo "This is used for rRNA alignments (where collapsing reads might lead to compression of real peaks)"

cat ${RAW_DATA_BASENAME}.adapterStripped.fq | /work/gbioinfo/linux/bin/fastq_quality_filter -q 25 -p 100 -Q 33 | awk '(NR-1)%4==0{print ">read_"(1+(NR-1)/4)} (NR-2)%4==0{print}' > ${RAW_DATA_BASENAME}_noComp_noLCstrip.fasta

echo "Debarcoding"

#/tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/venv-pycrac/bin/pyBarcodeFilter.py -f ${RAW_DATA_BASENAME}_noComp_noLCstrip.fasta -b $BARCODE_FILE --file_type=fasta
pyBarcodeFilter.py -f ${RAW_DATA_BASENAME}_noComp_noLCstrip.fasta -b $BARCODE_FILE --file_type=fasta

while read -r sample
do
      barcode=`echo ${sample//$'\t'*/}`
      sample_name=`echo ${sample//*$'\t'/}`
      echo mv ${RAW_DATA_BASENAME}_noComp_noLCstrip_${barcode}.fasta ${sample_name}_noComp_noLCstrip.fasta
      mv ${RAW_DATA_BASENAME}_noComp_noLCstrip_${barcode}.fasta ${sample_name}_noComp_noLCstrip.fasta

      if [ "$SPLIT_LINES_TOGGLE" = "Y" ]; then
            SPLIT_LINES_LC=`wc -l ${sample_name}_noComp_noLCstrip.fasta | awk '{printf "%d", $1/20}' | awk '{printf "%d", 4+2*$1}'`
            echo "Splitting file into 10, so" $SPLIT_LINES_LC "lines"
      else
            echo "Splitting file into" $SPLIT_LINES_LC "lines, as you specified on the command line"
      fi

      split -l $SPLIT_LINES_LC ${sample_name}_noComp_noLCstrip.fasta split.${sample_name}_noComp_noLCstrip.fasta.
      rm ${sample_name}_noComp_noLCstrip.fasta
        ls split.${sample_name}_noComp_noLCstrip.fasta.* | parallel '/tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/bin/prinseq-lite-0.20.4/prinseq-lite -lc_threshold 20 -lc_method dust -min_len 18 -fasta "{}" --out_good prinseq."{}" --out_bad rejected."{}"'
        cat prinseq.split.${sample_name}_noComp_noLCstrip.fasta.* > ${sample_name}_noComp_noLCstrip.prinseq.fasta

done < $SAMPLE_FILE


echo "Some tidying"

rm *split*

#rm ${RAW_DATA_BASENAME}.adapterStripped.fq

rm ${RAW_DATA_BASENAME}_comp_NNNN*
rm ${RAW_DATA_BASENAME}_comp_random_nucleotide_statistics.txt
rm ${RAW_DATA_BASENAME}_comp_others.fasta
#rm ${RAW_DATA_BASENAME}_comp_barcode_statistics.txt

rm ${RAW_DATA_BASENAME}_noComp_noLCstrip_NNNN*
rm ${RAW_DATA_BASENAME}_noComp_noLCstrip_random_nucleotide_statistics.txt
rm ${RAW_DATA_BASENAME}_noComp_noLCstrip_others.fasta
#rm ${RAW_DATA_BASENAME}_noComp_noLCstrip_barcode_statistics.txt



## Make bigwigs from standard preprocessing

for f in *.LCstripped.prinseq.fasta ; do ( /tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/bin/bowtie2-2.2.3/bowtie2 --sensitive -f -x /tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/References/Mus_musculus.GRCm38.75.dna.primary_assembly -S ${f%.fasta}.sam -U $f -p 7 ) ; done

for f in *.LCstripped.prinseq.sam ; do ( cat $f | awk -F "\t" '/^@/{print} !/^@/&&$5>=8{print}' > ${f%.sam}'_unique.sam' ) ; done

for f in *.LCstripped.prinseq_unique.sam ; do ( ( /work/gbioinfo/Appz/SAM_tools/samtools-current/samtools view -b -S $f > ${f%.sam}.bam ; /work/gbioinfo/Appz/SAM_tools/samtools-current/samtools sort ${f%.sam}.bam > ${f%.sam}-sorted.bam ; /work/gbioinfo/Appz/SAM_tools/samtools-current/samtools index ${f%.sam}-sorted.bam ) ) ; done

for f in *.LCstripped.prinseq_unique-sorted.bam ; do ( /work/gbioinfo/Appz/bedtools/bedtools-current/bin/bedtools bamtobed -cigar -i $f  > ${f%.bam}.bed ) ; done

for f in *.LCstripped.prinseq_unique-sorted.bed ; do ( /tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/bed2bigWig.sh -f $f -s ) ; done

echo "Information for trackDb.txt is output into trackDb.info"

for f in `cat ${SAMPLE_FILE}  | cut -f2,2` ; do (echo $f | awk '{print "track "$1"_CRAC_unique\ncontainer multiWig\nshortLabel "$1"_CRAC_unique\nlongLabel Depth of alignments of CRAC-seq reads\ntype bigWig\nconfigurable on\nvisibility full\naggregate transparentOverlay\nshowSubtrackColorOnUi on\nautoScale off\nviewLimits -10:10\nmaxHeightPixels 128:64:8\nwindowingFunction mean+whiskers\nsmoothingWindow off\npriority 1.1\n\ntrack "$1"_CRAC_fw\nparent "$1"_CRAC_unique\nbigDataUrl "$1".LCstripped.prinseq_unique-sorted_fw_UCSC.bigWig\nshortLabel "$1"_fw\nlongLabel Depth of alignments of CRAC-seq reads\ntype bigWig\ncolor 0,102,204\n\ntrack "$1"_CRAC_rev\nparent "$1"_CRAC_unique\nbigDataUrl "$1".LCstripped.prinseq_unique-sorted_rev_UCSC.bigWig\nshortLabel "$1"_rev\nlongLabel Depth of alignments of CRAC-seq reads\ntype bigWig\ncolor 0,102,204\n\n"}' ) ; done > ${RAW_DATA_BASENAME}.trackDb.info




#map and count reads per gene using STAR (splicing aware)
#edited: --outFilterMismatchNmax 3 to --outFilterMismatchNoverLmax 0.05: This means that for reads <20b no mismatches are allowed, 20-39b: 1 mismatch, 40-59b 2 mismatches and so on.
while read -r sample
do
    sample_name=`echo ${sample//*$'\t'/}`
    
    /work/gbioinfo/Appz/STAR/STAR-2.7.0a/source/STAR \
	--runMode alignReads --readFilesIn ${sample_name}.LCstripped.prinseq.fasta --runThreadN 16 \
	--genomeDir /tungstenfs/scratch/gbuehler/michi/Annotations/GENCODE/Mouse/release_M23/refstar_spliced/star2_7_0a_GRCm38.primary_assembly_gencodeM23_spliced_sjdb50 \
	--quantMode TranscriptomeSAM GeneCounts --outSJfilterReads All --outFilterType BySJout \
	--outFilterMultimapNmax 20 --outFilterMismatchNoverLmax 0.05 \
	--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMunmapped Within --outFileNamePrefix ${sample_name}.LCstripped.prinseq_

    /work/gbioinfo/Appz/STAR/STAR-2.7.0a/source/STAR \
	--runMode alignReads --readFilesIn ${sample_name}_noComp_noLCstrip.prinseq.fasta --runThreadN 16 \
	--genomeDir /tungstenfs/scratch/gbuehler/michi/Annotations/GENCODE/Mouse/release_M23/refstar_spliced/star2_7_0a_GRCm38.primary_assembly_gencodeM23_spliced_sjdb50 \
	--quantMode TranscriptomeSAM GeneCounts --outSJfilterReads All --outFilterType BySJout \
	--outFilterMultimapNmax 20 --outFilterMismatchNoverLmax 0.05 \
	--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMunmapped Within --outFileNamePrefix ${sample_name}_noComp_noLCstrip.prinseq_
    
done < $SAMPLE_FILE

module purge

#generate indexes for bam files
module load SAMtools
while read -r sample
do
    sample_name=`echo ${sample//*$'\t'/}`
    samtools index ${sample_name}.LCstripped.prinseq_Aligned.sortedByCoord.out.bam
    samtools index ${sample_name}_noComp_noLCstrip.prinseq_Aligned.sortedByCoord.out.bam

done < $SAMPLE_FILE


#generate bigwig files for tracks for genome browser (split by strand, normalized)
module load BEDTools/2.27.1-foss-2018b #to replace /work/gbioinfo/Appz/bedtools/bedtools-current/bin/

mkdir star_trackhub
while read -r sample
do
    sample_name=`echo ${sample//*$'\t'/}`
    chromSize="/tungstenfs/scratch/gbuehler/michi/Annotations/GENCODE/Mouse/release_M23/refstar_spliced/star2_7_0a_GRCm38.primary_assembly_gencodeM23_spliced_sjdb50/chrNameLength.txt"
    #normalisation afctor for read counts
    norm=$(/work/gbioinfo/Appz/SAM_tools/samtools-current/samtools idxstats ${sample_name}.LCstripped.prinseq_Aligned.sortedByCoord.out.bam | awk '{total+=$3} END{print total}')
    #generate bedgraph files
       samtools view --threads ${t} -b -F 0x100 -q 255 ${sample_name}.LCstripped.prinseq_Aligned.sortedByCoord.out.bam | genomeCoverageBed -split -strand + -bg -ibam stdin -g ${chromSize} | awk -vSIZE=${norm} '{print $1"\t"$2"\t"$3"\t"$4/SIZE*1000000}'> ${sample_name}_STARnorm_plus.bedgraph 
       samtools view --threads ${t} -b -F 0x100 -q 255 ${sample_name}.LCstripped.prinseq_Aligned.sortedByCoord.out.bam | genomeCoverageBed -split -strand - -bg -ibam stdin -g ${chromSize} | awk -vSIZE=${norm} '{print $1"\t"$2"\t"$3"\t"(-$4/SIZE*1000000)}'> ${sample_name}_STARnorm_minus.bedgraph
#convert bedgraph to bigwig
       /work2/gbuehler/fabiom/tools/bedGraphToBigWig ${sample_name}_STARnorm_plus.bedgraph ${chromSize} star_trackhub/${sample_name}_STARnorm_uni_plus.bw 
 /work2/gbuehler/fabiom/tools/bedGraphToBigWig ${sample_name}_STARnorm_minus.bedgraph ${chromSize} star_trackhub/${sample_name}_STARnorm_uni_minus.bw 
done < $SAMPLE_FILE


#map reads to repeat consesus sequences using Bowtie2
module purge
#module load Bowtie2/2.3.5.1-GCC-8.3.0
for f in *_noComp_noLCstrip.prinseq.fasta ; do ( /tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/bin/bowtie2-2.2.3/bowtie2 --sensitive -f -x /tungstenfs/scratch/gbuehler/michi/Annotations/Bowtie2/repeat.consensus/mm.repeat.consensus.seqs -S ${f%.fasta}.repeat.sam -U $f -p 7 ) ; done

for f in *_noComp_noLCstrip.prinseq.repeat.sam ; do ( ( /work/gbioinfo/Appz/SAM_tools/samtools-current/samtools view -b -S $f > ${f%.sam}.bam ; /work/gbioinfo/Appz/SAM_tools/samtools-current/samtools sort ${f%.sam}.bam > ${f%.sam}-sorted.bam ; /work/gbioinfo/Appz/SAM_tools/samtools-current/samtools index ${f%.sam}-sorted.bam ) ) ; done

#generate bedgraph files
while read -r sample
do
    sample_name=`echo ${sample//*$'\t'/}`
    chromSize="/tungstenfs/scratch/gbuehler/michi/Annotations/Bowtie2/mm.repeat.consensus.seqs.sizes.txt"
    #normalisation afctor for read counts
    norm=$(/work/gbioinfo/Appz/SAM_tools/samtools-current/samtools idxstats ${sample_name}_noComp_noLCstrip.prinseq.repeat-sorted.bam | awk '{total+=$3} END{print total}')
    #generate bedgraph files
       /work/gbioinfo/Appz/SAM_tools/samtools-current/samtools view -b ${sample_name}_noComp_noLCstrip.prinseq.repeat-sorted.bam | /work/gbioinfo/Appz/bedtools/bedtools-current/bin/genomeCoverageBed -split -strand + -bg -ibam stdin -g ${chromSize} | awk -vSIZE=${norm} '{print $1"\t"$2"\t"$3"\t"$4/SIZE*1000000}'> ${sample_name}_noComp_noLCstrip.prinseq.repeat_plus.bedgraph 
       /work/gbioinfo/Appz/SAM_tools/samtools-current/samtools view -b ${sample_name}_noComp_noLCstrip.prinseq.repeat-sorted.bam | /work/gbioinfo/Appz/bedtools/bedtools-current/bin/genomeCoverageBed -split -strand - -bg -ibam stdin -g ${chromSize} | awk -vSIZE=${norm} '{print $1"\t"$2"\t"$3"\t"(-$4/SIZE*1000000)}'> ${sample_name}_noComp_noLCstrip.prinseq.repeat_minus.bedgraph
done < $SAMPLE_FILE

rm -f *_noComp_noLCstrip.prinseq.repeat.sam


#load R module
module purge
module load R-BioC/current

## Run read assignment pipeline (splits into categories, makes plots, etc)

echo "Running pipeline to divide reads into repeat, genome and transcript files, and to make plots"

ls *.LCstripped.prinseq.fasta | parallel /tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/bin/categorise_CRAC_reads.sh -f {}


## Run rRNA analysis

for f in *_noComp_noLCstrip.prinseq.fasta ; do ( nice nohup /tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/bin/bowtie2-2.2.3/bowtie2 --sensitive -f -x /tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/alignment_to_rRNA_v2/pymol_rRNA_sequences_corrected_to_mouse -S ${f%_noComp_noLCstrip.prinseq.fasta}_pymol.sam -U $f -p 20 ) ; done

ls *_pymol.sam | parallel '/tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/bin/process_pymol_sam.awk {} > {.}_20ntOrLongerPerfectMatches.sam'

for f in *_pymol_20ntOrLongerPerfectMatches.sam ; do ( ( /work/gbioinfo/Appz/SAM_tools/samtools-current/samtools view -b $f > ${f%.sam}.bam ; /work/gbioinfo/Appz/SAM_tools/samtools-current/samtools sort ${f%.sam}.bam > ${f%.sam}-sorted.bam ; /work/gbioinfo/Appz/SAM_tools/samtools-current/samtools index ${f%.sam}-sorted.bam ) ) ; done

for f in *_pymol_20ntOrLongerPerfectMatches-sorted.bam ; do ( nohup /work/gbioinfo/Appz/SAM_tools/samtools-1.3/samtools depth -d 1000000000 -a $f > ${f%_20ntOrLongerPerfectMatches-sorted.bam}_rRNA_pileup.txt ) ; done

for f in *_pymol_rRNA_pileup.txt ; do ( awk 'BEGIN{FS="\t"; OFS="\t"; maxval=0} NR==FNR&&$3>maxval{maxval=$3} NR!=FNR{$3=999*$3/maxval; printf "%s\t%i\t%6.2f\n", $1, $2, $3}' $f $f > ${f%.txt}_NormAndPadded.txt) ; done

for f in /tungstenfs/scratch/gbuehler/bioinfo/BCLIP_pipelines/alignment_to_rRNA_v2/*S.pdb ; do ( for g in *_pymol_rRNA_pileup_NormAndPadded.txt ; do ( newname=`basename $f`; newname=${newname/.pdb/}; newname=${newname}_${g%_pymol_rRNA_pileup_NormAndPadded.txt}.pdb; awk -v filename=${f##*/} 'BEGIN{FS="\t"; OFS="\t"} NR==FNR{score[$1,$2] = $3} NR!=FNR{match($0, / [SL] *[0-9]+ /); subunit = substr($0, RSTART+1,1); coord = substr($0, RSTART+1, RLENGTH-2); sub(/[LS] */, "", coord); print substr($0,1,60) score[filename,coord] "           " substr($0,78,4)}' $g $f > $newname ) ; done ) ; done




## Tidying

#mv *_pymol* /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/pymol/intermediate_files
#mv *pdb /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/pymol/pdb

#mv ${RAW_DATA_BASENAME}_comp.fasta /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/standard_preprocessing
#rm ${RAW_DATA_BASENAME}_noComp_noLCstrip.fasta # Intermediate file for making the rRNA alignments => this is the pre-debarcoded file that had adapters clipped and then quality filtering (100 % >= 25)

#mv *bam *bai /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/bam
#mv *sam /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/sam
#mv ${RAW_DATA_BASENAME}.trackDb.info  *bigWig /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/bigWig
#mv *LCstripped.prinseq.fasta /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/comp_LCstripped_prinseq_fasta
#mv *_noComp_noLCstrip.prinseq.fasta /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/noComp_noLCstrip_prinseq_fasta
#mv *prinseq.fasta /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/comp_prinseq_fasta

#rm *bedGraph
rm -f *bed

#mv *EntireReadClassOut.table /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/read_classification_output/EntireReadClassOut/
#mv read_classification_output/*start* read_classification_output/*stop* read_classification_output/*CDS_phasing* /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/read_classification_output/CDS_analyses/
#mv read_classification_output/*TSS_metaplot* /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/read_classification_output/TSS_analyses/
#mv read_classification_output/*sam /tung/Users/matyasflemrstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/read_classification_output/sam/
#mv read_classification_output/*bam read_classification_output/*bai /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/read_classification_output/bam/
#mv read_classification_output/*counts.tab read_classification_output/*biotype_totals.tab /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/read_classification_output/counts/
#mv read_classification_output/*tails.pdf read_classification_output/*appris_unique_counts_byNonencodedTail.tab /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/read_classification_output/tail_analyses/
#mv read_classification_output/*E1_E2* read_classification_output/*E2_E3* read_classification_output/*Exon1to5* read_classification_output/*last_EE* /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/read_classification_output/splicing_analyses/
#mv read_classification_output/*read_classifications.tab /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/read_classification_output/read_classifications/
#mv read_classification_output/*.table /tungstenfs/scratch/gbuehler/tuckalex/pipeline_v2_output/read_classification_output/tables/
#rm -r read_classification_output -f

#rm -f *.fq
#rm -f *.fastq
#rm -f *.fasta
