#!/bin/bash

#Look at fragment size
# EXAMPLE FORMAT - DO NOT RUN
# bamPEFragmentSize -b <bamfiles> \
#                  -hist <plot output filename> \
#                  --maxFragmentLength <maximum fragment length> \
#                  -T <plot title>

##Don't use below, the Donor7256-NK.chr4.bam file is corrupted
#bamPEFragmentSize -b *.bam -hist fragmentSizes.png --maxFragmentLength 500 \
#                  -T "Fragment sizes of ATAC-seq data"

bamPEFragmentSize -b Donor2596-NK.chr4.bam Donor4983-HSC.chr4.bam Donor5483-NK.chr4.ba -hist fragmentSizes.png --maxFragmentLength 500 -T "Fragment sizes of ATAC-seq data"
#What does the peak at around 200 bp likely represent?
#Answer is 'Reads spanning nucleosomes (DNA-histone structures)', but why 200bp?

cp fragmentSizes.png ~/bios201/workshop4


#Making coverage track to show how many fragments map to a particular region of the genome
# EXAMPLE FORMAT -- DO NOT RUN
# bamCoverage --bam <input.bam> -o <output> --binSize 100 --normalizeUsingRPKM --extendReads
bamCoverage --bam Donor2596-NK.chr4.bam -o Donor2596-NK.chr4.bw --binSize 100 --normalizeUsingRPKM --extendReads
bamCoverage --bam Donor4983-HSC.chr4.bam -o Donor4983-HSC.chr4.bw --binSize 100 --normalizeUsingRPKM --extendReads
bamCoverage --bam Donor5483-NK.chr4.bam -o Donor5483-NK.chr4.bw --binSize 100 --normalizeUsingRPKM --extendReads
#Create per base track of ATAC-seq insertions; precise base pair that represents center of Tn5...?
bamCoverage --bam Donor2596-NK.chr4.bam -o Donor2596-NK.ins.bw  --binSize 1 --Offset 5 5

#Copy data to github to view on UCSC's Gene Browser
cp *.bw ~/bios201/workshop4
#Once uploaded to github, need to use this link for UCSC to access: https://github.com/jeremy29tien/bios201/blob/master/workshop4/Donor2596-NK.chr4.bw?raw=true
#EDIT: the genome browser is unable to uncompress the above link. this is somewhat problematic ... 

#The hg19.refGeneReviewedValidated.tss.chr4.bed file doesn't exist
computeMatrix reference-point -S Donor2596-NK.ins.bw \
  -R hg19.refGeneReviewedValidated.tss.chr4.bed \
  -o Donor2596-NK.matrix.gz --referencePoint TSS \
  --binSize 10 --missingDataAsZero -b 2000 -a 2000

plotHeatmap -m Donor2596-NK.matrix.gz -out Donor2596-NK_heatmap.png \
  --xAxisLabel 'Distance (bp)' --samplesLabel Insertions --zMin 0 -z ATAC

#Use MACS2 to call peaks. Peaks are areas of genome where we have a pileup of signal. In ATAC-seq, this represents accessible regions of genome
#ISSUE: The program 'macs2' is currently not installed. To run 'macs2' please ask your administrator to install the package 'macs'
macs2 callpeak --treatment Donor2596-NK.chr4.bam --name Donor2596-NK --format BAMPE --nomodel --call-summits --nolambda --keep-dup all -q 0.01 -g 1.7e8
macs2 callpeak --treatment Donor4983-HSC.chr4.bam --name Donor4983-HSC --format BAMPE --nomodel --call-summits --nolambda --keep-dup all -q 0.01 -g 1.7e8
macs2 callpeak --treatment Donor5483-NK.chr4.bam --name Donor5483-NK --format BAMPE --nomodel --call-summits --nolambda --keep-dup all -q 0.01 -g 1.7e8

#Visualizing Peaks -- be sure to update the URL
echo 'track type=narrowPeak name="Donor2596 NK Peaks"' | cat - Donor2596-NK_peaks.narrowPeak > \
      ~/afs-home/WWW/Donor2596-NK_peaks.narrowPeak 

#Manipulating peak files
#Finds the intersection between peaks
bedtools intersect -a Donor2596-NK_peaks.narrowPeak -b Donor5483-NK_peaks.narrowPeak > NK_intersection.bed

bedtools intersect -a Donor7256-HSC_peaks.narrowPeak -b Donor4983-HSC_peaks.narrowPeak > HSC_intersection.bed

#Combines peaks and sorts them
cat *.narrowPeak | bedtools sort -i stdin > all_peaks.bed
#Merge overlapping peaks
bedtools merge -i all_peaks.bed

#Find motifs (short, repeated sequences in genome) in peaks
##Slow - do not run during workshop
#findMotifsGenome.pl NK_intersection.bed chr4.fa NK_motifs -size given

#findMotifsGenome.pl HSC_intersection.bed chr4.fa HSC_motifs -size given
