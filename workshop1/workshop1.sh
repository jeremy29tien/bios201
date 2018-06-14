#!/bin/bash

#Looks at fastq files
head -8 NA12878_R1.fastq
head -8 NA12878_R2.fastq

#Sequencing quality control (fastqc)
fastqc --outdir fastqcOut --format fastq NA12878_R1.fastq NA12878_R2.fastq

# <output_dir> is the name of a folder where you want the results written
# <format> is the format of the sequencing reads (optional)
# <r1_fastq> is the first set of reads
# <r2_fastq> is the second set of reads 


#Removing adapters, trimming reads, filtering
#cutadapt -q <minimum_quality> --minimum-length <minimum_length> -a <adapter on first read> -A <adapter on second read> -o <r1_trimmed_fastq> -p <r2_trimmed_fastq> <r1_fastq> <r2_fastq>
cutadapt -q 10 --minimum-length 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o r1_trimmed_fastq -p r2_trimmed_fastq NA12878_R1.fastq NA12878_R2.fastq
# <minimum_quality> is the minimum base quality to allow (trimmed otherwise)
# <minimum_length> is the minimum length both pairs of reads need to be to be included
# <adapter on first read> for us: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
# <adapter on second read> for us: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
# <r1_trimmed_fastq> the name of the output file for the first read file 
# <r2_trimmed_fastq> the name of the output file for the second read file 
# <r1_fastq> is the first set of reads
# <r2_fastq> is the second set of reads 

# note: the adapter sequence is the reverse complement of the actual adapter sequence

#Align paired-end DNA sequences
#bwa index grch37.fa
bwa mem grch37.fa NA12878_R1.qc.trimmed.fastq NA12878_R2.qc.trimmed.fastq | samtools-1.3.1 sort --output-fmt BAM > NA12878.bam
bwa mem grch37.fa NA12891_R1.qc.trimmed.fastq NA12891_R2.qc.trimmed.fastq | samtools-1.3.1 sort --output-fmt BAM > NA12891.bam
bwa mem grch37.fa NA12892_R1.qc.trimmed.fastq NA12892_R2.qc.trimmed.fastq | samtools-1.3.1 sort --output-fmt BAM > NA12892.bam

#Mark PCR duplicates
java -Xmx2g -jar $PICARD MarkDuplicates \
        INPUT=NA12878.bam \
        OUTPUT=NA12878.markduplicates.bam \
        METRICS_FILE=NA12878.markduplicates.metrics.txt \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        CREATE_INDEX=true \
        TMP_DIR=/tmp

java -Xmx2g -jar $PICARD MarkDuplicates \
        INPUT=NA12891.bam \
        OUTPUT=NA12891.markduplicates.bam \
        METRICS_FILE=NA12891.markduplicates.metrics.txt \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        CREATE_INDEX=true \
        TMP_DIR=/tmp

java -Xmx2g -jar $PICARD MarkDuplicates \
        INPUT=NA12892.bam \
        OUTPUT=NA12892.markduplicates.bam \
        METRICS_FILE=NA12892.markduplicates.metrics.txt \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        CREATE_INDEX=true \
        TMP_DIR=/tmp

#Mark all reads as coming from the same read group
java -jar $PICARD AddOrReplaceReadGroups \
        I=NA12878.markduplicates.bam \
        O=NA12878.markduplicates.rg.bam \
        RGID=1 RGLB=lib1 RGPL=illumina \
        RGPU=unit1 RGSM=NA12878 \
        SORT_ORDER=coordinate CREATE_INDEX=True

java -jar $PICARD AddOrReplaceReadGroups \
        I=NA12891.markduplicates.bam \
        O=NA12891.markduplicates.rg.bam \
        RGID=1 RGLB=lib1 RGPL=illumina \
        RGPU=unit1 RGSM=NA12891 \
        SORT_ORDER=coordinate CREATE_INDEX=True

java -jar $PICARD AddOrReplaceReadGroups \
        I=NA12892.markduplicates.bam \
        O=NA12892.markduplicates.rg.bam \
        RGID=1 RGLB=lib1 RGPL=illumina \
        RGPU=unit1 RGSM=NA12892 \
        SORT_ORDER=coordinate CREATE_INDEX=True

#Perform basequality recalibration aka BQSR
gatk BaseRecalibrator -R grch37.fa \
        -I NA12878.markduplicates.rg.bam \
        -known-sites knownSites.vcf -O NA12878.recal_data.table
gatk BaseRecalibrator -R grch37.fa \
        -I NA12891.markduplicates.rg.bam \
        -known-sites knownSites.vcf -O NA12891.recal_data.table
gatk BaseRecalibrator -R grch37.fa \
        -I NA12892.markduplicates.rg.bam \
        -known-sites knownSites.vcf -O NA12892.recal_data.table

gatk ApplyBQSR -R grch37.fa -I NA12878.markduplicates.rg.bam --bqsr-recal-file NA12878.recal_data.table -O NA12878.markduplicates.rg.bqsr.bam
gatk ApplyBQSR -R grch37.fa -I NA12891.markduplicates.rg.bam --bqsr-recal-file NA12891.recal_data.table -O NA12891.markduplicates.rg.bqsr.bam
gatk ApplyBQSR -R grch37.fa -I NA12892.markduplicates.rg.bam --bqsr-recal-file NA12892.recal_data.table -O NA12892.markduplicates.rg.bqsr.bam

#Calling variants in individual samples
gatk HaplotypeCaller -R grch37.fa -I NA12878.markduplicates.rg.bqsr.bam --emit-ref-confidence GVCF -O NA12878.g.vcf
gatk HaplotypeCaller -R grch37.fa -I NA12891.markduplicates.rg.bqsr.bam --emit-ref-confidence GVCF -O NA12891.g.vcf
gatk HaplotypeCaller -R grch37.fa -I NA12892.markduplicates.rg.bqsr.bam --emit-ref-confidence GVCF -O NA12892.g.vcf

#Joint calling across samples
gatk GenomicsDBImport     -V NA12878.g.vcf     -V NA12891.g.vcf     -V NA12892.g.vcf     --genomicsdb-workspace-path combined-reads --intervals 17
gatk GenotypeGVCFs -R grch37.fa -V gendb://combined-reads/ -O raw_variants.vcf

#Genotype phasing
gatk CalculateGenotypePosteriors    -R grch37.fa    -V raw_variants.vcf   -ped family.ped -O phased_variants.vcf

#Filtering low-quality variants
gatk VariantFiltration -R grch37.fa -V phased_variants.vcf --filter-expression "DP < 10 || QD < 2.0 || FS > 60.0 || MQ < 40.0" --filter-name "BIOS201_FILTER" -O flagged_snps.vcf
gatk SelectVariants -R grch37.fa -V flagged_snps.vcf -O filtered_snps.vcf --select-type-to-include SNP --restrict-alleles-to BIALLELIC

#Compare variants to platinum genomes
vcftools --vcf raw_variants.vcf --diff platinum_trio.vcf --diff-indv-discordance  --out raw_to_platinum_comparison
vcftools --vcf filtered_snps.vcf --diff platinum_trio.vcf --diff-indv-discordance  --out filtered_to_platinum_comparison

less raw_to_platinum_comparison.diff.indv
less filtered_to_platinum_comparison.diff.indv
