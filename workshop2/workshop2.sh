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

#Visualizing alignments with IGV
#Sorting alignments by coordinate and indexing
samtools sort bam_pass2/Norm1_Aligned.out.bam bam_pass2/Norm1_Aligned.out.sorted
samtools index bam_pass2/Norm1_Aligned.out.sorted.bam

#Open new terminal window, ssh with a -X tag, and run `igv`

#Generating mapping metrics with Picard
java -jar $PICARD CollectRnaSeqMetrics \
     REF_FLAT=annotation/refFlat.chr10.txt \
     I=bam_pass2/Norm1_Aligned.out.sorted.bam \
     O=bam_pass2/Norm1_Aligned.out.sorted.metrics.txt \
     STRAND=SECOND_READ_TRANSCRIPTION_STRAND

###NOT ON PERSONAL COMPUTER 
##ON PERSONAL COMPUTER --> IGNORE
#Preliminary inspection of count data
#download files from /users/j29tien/bios201/web.stanford.edu/class/bios201/workshop2
#scp j29tien@durga.stanford.edu:/users/j29tien/bios201/web.stanford.edu/class/bios201/workshop2/bam_pass2/*_ReadsPerGene.out.tab /Users/JeremyTien/Documents/SIMR\ 2018/bios201/workshop2
#scp j29tien@durga.stanford.edu:/users/j29tien/bios201/web.stanford.edu/class/bios201/workshop2/annotation/ensg2hgnc.txt /Users/JeremyTien/Documents/SIMR\ 2018/bios201/workshop2

#Run this before launching R/RStudio:
#export CONDA_BUILD_SYSROOT=$(xcrun --show-sdk-path)

#USING DURGA (ALREADY INSTALLED R, INCLUDING DESEQ2)
ssh j29tien@durga.stanford.edu -X #IMPORTANT that you launch X11
cd bios201/web.stanford.edu/class/bios201/workshop2
R
> getwd()
> library(DESeq2)
> library(pheatmap)
> samples = paste0(rep(c('Norm','IPF'), each = 8), rep(c(1:8), 2))
> 
> ## Read in counts for the first sample and look at the head
> data = read.table(paste0(samples[1], '_ReadsPerGene.out.tab'), 
+   header = FALSE, stringsAsFactors = FALSE, 
+   col.names = c('Gene','Unstranded','First','Second'))
> head(data)

> counts = data[-c(1:4), c('Gene','Second')]
> colnames(counts)[2] = samples[1]
> ## Take a look at what it looks like now
> head(counts)

> for (sample in samples[2:length(samples)]) {
+     sampleData = read.table(paste0(sample,'_ReadsPerGene.out.tab'), 
+     header = FALSE, stringsAsFactors = FALSE, skip = 4)[,c(1,4)]
+     colnames(sampleData) = c('Gene', sample)
+     counts = merge(counts, sampleData, by = 'Gene')
+ }
> genes = read.table('ensg2hgnc.txt', header = FALSE, 
+       stringsAsFactors = FALSE, col.names = c('ENSG', 'HGNC'))
> counts = merge(counts, genes, by.x = 'Gene', by.y = 'ENSG')
> rownames(counts) = paste(counts$HGNC, counts$Gene, sep = '_')
> counts = counts[, samples] # drop gene columns
> group = sapply(samples, function (s) substr(s, 1, nchar(s)-1))
> columnData = as.data.frame(group)
> dataset = DESeqDataSetFromMatrix(countData = counts,
+                                  colData = columnData,
+                                  design = ~ group)

> dataset = estimateSizeFactors(dataset)
> log2normCounts = log2(counts(dataset, normalized = TRUE) + 1)

> top =  order(rowMeans(log2normCounts), decreasing = TRUE)[1:30]
> pheatmap(log2normCounts[top, ], show_rownames = TRUE, annotation_col = columnData)
