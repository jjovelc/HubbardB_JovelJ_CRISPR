#!/usr/bin/bash

for FILE in *[1-2]/*_merged.fastq.gz
do
	cutadapt -j 8 --nextseq-trim=20 -e 0.10 -O 16 --discard-untrimmed \
	-a 'ATGGTGAGCAAGGGCGAGG...GCATGGACGAGCTGTACAAGTAA' \
	-o "${FILE/_merged.fastq.gz/_merged.extracted.fastq.gz}" \
	"$FILE"
done
