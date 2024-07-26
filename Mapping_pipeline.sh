## Input files

files="E2_rep1.fastq E2_rep2.fastq untr_rep1.fastq untr_rep2.fastq"
in_path="/labs/Mellone/Asna/RNAseq"
out_path="/labs/Mellone/Asna/RNAseq/processed"
log_path="/labs/Mellone/Asna/RNAseq/log"
quality="30"
index='/isg/shared/databases/alignerIndex/animal/hg19_ucsc/hg19_HISAT2/hg19'
trimmed='/labs/Mellone/Asna/RNAseq/Trimmed/'
bam='/labs/Mellone/Asna/RNAseq/bam/'
sorted_bam='/labs/Mellone/Asna/RNAseq/bam/*_sorted.bam'
gtf='/labs/Mellone/Asna/RNAseq/gtf/'
mergelist='/labs/Mellone/Asna/RNAseq/gtf/mergelist.txt'
Hg19_chr='/labs/Mellone/Asna/hg19/hg19_chrInfo.txt'
Hg19Bed='/labs/Mellone/Asna/gencode_hg19.gtf_2021.bed'
SummitsBed='/labs/Mellone/Asna/ER_peaks/ER_E2_summits.bed'
Hg19GTF='/labs/Mellone/Asna/gencode_hg19.gtf'
E2_gtf='/labs/Mellone/Asna/E2_merged.annotated.gtf'
DETableOutput='/labs/Mellone/Asna/E2_gene_matrix.csv'
sampletxt='/labs/Mellone/Asna/samples.txt'
ER peaks='/labs/Mellone/Asna/ER_peaks/'

#!/bin/bash
#SBATCH --job-name=pipeline.sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=20G
#SBATCH --mail-user=asna.amjad@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

###### PIPEINE ######
#####################
module load fastqc
module load fastx
module load hisat2
module load samtools
module load stringtie
module laod htseq
module load gffcompare
module load bedtools

for file in ${files}
do filename=`echo ${file} | cut -d "." -f 1` echo "Processing ${file}..." echo "Performing FastQC of raw data..."
fastqc "${in_path}/${file}"
echo "Removing adapter sequences from raw data..." fastx_clipper -Q 33 -a ${adapter} -i ${in_path}/${file} -o ${out_path}/${filename}_clip.fastq
echo "Trimming reads to Q-score of ${quality} and minimum length of ${min_len}..." fastq_quality_trimmer -Q 33 -t ${quality} -l ${min_len} -i ${in_path}/${file} -o ${out_path}/${filename}_clip_trim.fastq
echo "Performing FastQC on clipped & trimmed data..."
fastqc "${out_path}/${filename}_clip_trim.fastq" echo "Mapping data to reference genome..." hisat2 -p8 --rna-strandness F -x ${index} -U t${in_path}/${file} -S ${out_path}/${filename}.sam
samtools view -@ 16 -bS ${out_path}/${filename}.sam > ${out_path}/${filename}.bam
echo "Sorting data by genome position" samtools sort -@ 16 ${out_path}/${filename}.bam -o ${out_path}/${filename}_sorted.bam

#Stringtie to assemble transcripts stringtie -p 8 -G ${Hg19GTF} -o ${gtf}${filename}.gtf -l Treated1 ${out_path}/${filename}_sorted.bam #merge transcripts stringtie --merge -p 8 -G ${Hg19GTF} -o stringtie_merged.gtf mergelist #gffcompare to compare to annotations
gffcompare -R -r ${Hg19GTF} -o stringtie_merged stringtie_merged.gtf echo stringtie merge complete

# Get counts htseq-count -m union -r pos -f bam -i gene_id -a 10 --stranded=no ${out_path}/${filename}_sorted.bam ${E2_gtf} > ${filename}.counts
stringtie -e -G ${E2_gtf} -p 4 -A ${filename}_geneCounts.tab ${out_path}/${filename}_sorted.bam -o ${filename}_geneCounts.gtf

done 2>&1 | tee -a ${log_path}/${filename}_log.txt

###### Generate BedGraph Files ######
#####################################

### Convert bed files to bed then sorting

bedtools bamtobed -i E2_rep1_sorted.bam > E2_rep1.bed
sortBed -i E2_rep1.bed > E2_rep1_sorted.bed
bedtools bamtobed -i E2_rep2_sorted.bam > E2_rep2.bed
sortBed -i E2_rep2.bed > E2_rep2_sorted.bed
bedtools bamtobed -i untr_rep1_sorted.bam > untr_rep1.bed
sortBed -i untr_rep1.bed > untr_rep1_sorted.bed
bedtools bamtobed -i untr_rep2_sorted.bam > untr_rep2.bed
sortBed -i untr_rep2.bed > untr_rep2_sorted.bed

### Generate bedgraph files

bedtools genomecov -bg -i E2_rep1_sorted.bed -g ${Hg19_chr} > E2_rep1.bedgraphbedtools genomecov -bg -i E2_rep2_sorted.bed -g ${Hg19_chr} > E2_rep2.bedgraphbedtools genomecov -bg -i untr_rep1_sorted.bed -g ${Hg19_chr} > untr_rep1.bedgraphbedtools genomecov -bg -i untr_rep2_sorted.bed -g ${Hg19_chr} > untr_rep2.bedgraph # Generate bedgraph tracklinesawk 'BEGIN { print "browser position chr11:5,289,521-5,291,937"
print "track type=bedGraph name=\"E2_rep1_bedgraph\" description=\"E2_rep1_bedgraph\" visibility=full autoScale=on alwaysZero=on
color=0,125,0 windowingFunction=maximum"} {print $0}' tracks/E2_rep1.bedgraph > tracks/E2_rep1_ucsc_track.bedgraph

### Generate bedgraph tracklines

awk 'BEGIN { print "browser position chr11:5,289,521-5,291,937"
print "track type=bedGraph name=\"E2_rep1_bedgraph\" description=\"E2_rep1_bedgraph\" visibility=full autoScale=on alwaysZero=on
color=0,125,0 windowingFunction=maximum"} {print $0}' tracks/E2_rep1.bedgraph > tracks/E2_rep1_ucsc_track.bedgraph

awk 'BEGIN { print "browser position chr11:5,289,521-5,291,937"
print "track type=bedGraph name=\"E2_rep2_bedgraph\" description=\"E2_rep2_bedgraph\" visibility=full autoScale=on alwaysZero=on
color=0,125,0 windowingFunction=maximum"} {print $0}' tracks/E2_rep2.bedgraph > tracks/E2_rep2_ucsc_track.bedgraph

awk 'BEGIN { print "browser position chr11:5,289,521-5,291,937"
print "track type=bedGraph name=\"Untr_rep1_bedgraph\" description=\"Untr_rep1_bedgraph\" visibility=full autoScale=on alwaysZero=on
color=0,0,205 windowingFunction=maximum"}
{print $0}' tracks/untr_rep1.bedgraph > tracks/untr_rep1_ucsc_track.bedgraph

awk 'BEGIN { print "browser position chr11:5,289,521-5,291,937"
print "track type=bedGraph name=\"Untr_rep2_bedgraph\" description=\"Untr_rep2_bedgraph\" visibility=full autoScale=on alwaysZero=on
color=0,0,205 windowingFunction=maximum"}
{print $0}' tracks/untr_rep2.bedgraph > tracks/untr_rep2_ucsc_track.bedgraph

###### Compare ChIP-seq with Regulated Genes ######
###################################################
bedtools closest -d -a ${SummitsBed} -b ${Hg19Bed} > closest.bed

### Find the distances of up, down, &non reg from peaks
cat E2upreg_genes.txt > all_genes_combine.txt
cat E2downreg_genes.txt >> all_genes_combine.txt

### Grep up, down,&non regulated genes

grep -v -F -f all_genes_combine.txt all_genes.txt > nonregulatedgenes.txt
grep -F -f E2upreg_genes.txt closest.bed > upgenesclosest.txt
grep -F -f E2downreg_genes.txt closest.bed > downgenesclosest.txt
grep -F -f nonregulatedgenes.txt closest.bed > nonregulatedgenesclosest.txt
grep -f E2upreg_genes.txt closest.bed | wc -l
grep -f E2downreg_genes.txt closest.bed | wc -l
grep -f E2nonreg_genes.txt closest.bed | wc -l# Upload these into RStudio to create CDF plot & box plots
