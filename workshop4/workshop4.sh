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


