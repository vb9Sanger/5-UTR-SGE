#!/bin/bash

#The sequence have the same strand as the primer used.
#Remove Hiseq/novaseq adaptor and revese complement the reverse primer

# Trim the illumina adaptor
#bash ../script/bsub_SGE_demultiplex_and_trim_prep.sh -J demultiplex -I [1-6] --memory 5000 -t demultiplex_UTR.txt -s

#Remove the primer sequence by cutadapt and count by awk
bash ../script/bsub_extract_and_count.sh -J counting -I [1-6] -t demultiplex_UTR.txt -p SE --core 5 --memory 5000 --me cutadapt --mc awk -s


##Manual quick QC

##Quick check the cutadapt output. Divide the file size of un extracted to the extracted fq.gz. It should be small, <5%
#cd extracted;
#ls -l | sort -k 9,9 | grep fq.gz| awk 'NR%2{p=$5;next}{print ($5/p)}' | sort -n

##Check cutadapt output. read with adaptor should be high ,>98% and have enough read, >1M
#find *.stdout.txt | parallel -j 1 "paste <(echo {} | sed 's/_R1_cutadapt_extracted.stdout.txt//g') <(cat {} | grep 'Reads with adapters:')"
