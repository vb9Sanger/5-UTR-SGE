#!/bin/bash

INPUT_FASTQ=$1
OUTDIR=subsampled
mkdir -p $OUTDIR

NB_SUBSAMPLES=10000000

OUTPUT_FASTQ=$OUTDIR/$(basename $INPUT_FASTQ)
echo $OUTPUT_FASTQ
zcat $INPUT_FASTQ | awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t");} }' |
awk -v k=$NB_SUBSAMPLES 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x- 1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' |
awk  -F"\t" -v out=$OUTPUT_FASTQ '{print $1"\n"$2"\n"$3"\n"$4}' | gzip > $OUTPUT_FASTQ
