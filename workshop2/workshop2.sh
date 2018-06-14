#!/bin/bash

#GenePred format
less -S annotation/UCSC_table_browser_chr10.txt
grep "SFTPA2" annotation/UCSC_table_browser_chr10.txt | column -t | less -S

#GTF/GFF format
less -S annotation/gencode.v19.annotation_chr10.gtf
grep "SFTPA2" annotation/gencode.v19.annotation_chr10.gtf | less -S

#Spliced Alignment and read counting
