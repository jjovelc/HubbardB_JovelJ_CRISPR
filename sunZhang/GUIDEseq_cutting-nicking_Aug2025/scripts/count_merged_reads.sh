#!/usr/bin/bash

for DIR in *[1-2]
do
    echo "SAMPLE: $DIR"
    
    # Count reads in original library (using R1 reads)
    READS_ORIGINAL="$(($(zcat ${DIR}/*_1.fq.gz | wc -l) / 4))"
    
    # Count reads in merged library
    READS_MERGED="$(($(zcat ${DIR}/*_merged.fastq.gz | wc -l) / 4))"
    
    # Calculate percentage of reads that were successfully merged
    if [ $READS_ORIGINAL -gt 0 ]; then
        PERC_MERGED=$((($READS_MERGED * 100) / $READS_ORIGINAL))
    else
        PERC_MERGED=0
    fi
    
    echo "# reads original library: $READS_ORIGINAL"
    echo "# reads merged library: $READS_MERGED"
    echo "Percentage of reads merged: $PERC_MERGED%"
    echo "______________________"
done
