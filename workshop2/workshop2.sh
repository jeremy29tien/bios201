#!/bin/bash

less -S annotation/UCSC_table_browser_chr10.txt
grep "SFTPA2" annotation/UCSC_table_browser_chr10.txt | column -t | less -S


