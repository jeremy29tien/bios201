#!/bin/bash

# link to full course located https://github.com/zaczap/bios201/tree/master/Workshop2

#GenePred format
less -S annotation/UCSC_table_browser_chr10.txt
grep "SFTPA2" annotation/UCSC_table_browser_chr10.txt | column -t | less -S

#GTF/GFF format
less -S annotation/gencode.v19.annotation_chr10.gtf
grep "SFTPA2" annotation/gencode.v19.annotation_chr10.gtf | less -S

#Spliced Alignment and read counting

##This builds the genome index
## Do not run this, we did it for you already to save time.
#genomeDir=/afs/ir/users/e/t/etsang/bios201/workshop2/STAR_hg19_chr10
#STAR --runThreadN 20 \
#     --runMode genomeGenerate \
#     --genomeDir $genomeDir \
#     --genomeFastaFiles chr10.fa \
#     --sjdbGTFfile annotation/gencode.v19.annotation_chr10.gtf \
#     --sjdbOverhang 100  # readLength - 1

#This will map the RNA-seq reads (in GTF gene-annotated format) to the genome
#First pass through just the Norm1 read
genomeDir=/users/mgloud/bios201/workshop2/STAR_hg19_chr10 

STAR-2.5.2b --runThreadN 4 \
     --genomeDir $genomeDir \
     --readFilesIn fastq/Norm1_R1.fastq fastq/Norm1_R2.fastq \
     --outFileNamePrefix bam_pass1/Norm1_ \
     --outSAMtype BAM Unsorted

#Second pass using the splice junctions identified by STAR from the first pass
# The first line collects the names of all the junction files from the first pass.
# We can then pass this information to the aligner below.

junctions=`ls bam_pass1/*_SJ.out.tab`

STAR-2.5.2b --runThreadN 4 \
     --genomeDir $genomeDir \
     --readFilesIn fastq/Norm1_R1.fastq fastq/Norm1_R2.fastq \
     --outFileNamePrefix bam_pass2/Norm1_ \
     --outSAMtype BAM Unsorted \
     --sjdbFileChrStartEnd $junctions \
     --quantMode GeneCounts

#Looking at SAM/BAM files
samtools view bam_pass2/Norm1_Aligned.out.bam | head -n1 #use https://broadinstitute.github.io/picard/explain-flags.html to understand flags
samtools view bam_pass2/Norm1_Aligned.out.bam | \
	 grep "HWI-ST689:184:D0YYWACXX:1:2315:14384:11932_1:N:0:CGATGT"

