#!/bin/bash

echo "Starting CRISPR indel detection analysis..."
echo "=========================================="

# Create output directory
mkdir -p indels_results

# Clean reference file if needed
if [ ! -f "genomic_regions_clean.fa" ]; then
    echo "Cleaning reference file..."
    sed 's/-//g' genomic_regions.fa > genomic_regions_clean.fa
fi

# Process each sample directory
for DIR in *_[12]; do
    if [ -d "$DIR" ]; then
        echo ""
        echo "Processing directory: $DIR"
        
        # Check if merged file exists
        if ls ${DIR}/*_merged.fastq.gz 1> /dev/null 2>&1; then
            echo "Found merged file(s) in $DIR"
            
            # Run analysis for each merged file in the directory
            for FASTQ in ${DIR}/*_merged.fastq.gz; do
                if [ -f "$FASTQ" ]; then
                    BASENAME=$(basename "$FASTQ" _merged.fastq.gz)
                    OUTPUT="indels_results/${DIR}_${BASENAME}_indels.txt"
                    
                    echo "Analyzing: $FASTQ"
                    echo "Output: $OUTPUT"
                    
                    python3 CRISPR_indel_detect.py \
                        "$FASTQ" \
                        genomic_regions_clean.fa \
                        --output "$OUTPUT" \
                        --min-quality 20 \
                        --max-indel 100
                    
                    if [ $? -eq 0 ]; then
                        echo "✓ Successfully processed $FASTQ"
                    else
                        echo "✗ Error processing $FASTQ"
                    fi
                fi
            done
        else
            echo "Warning: No merged files found in $DIR"
        fi
    fi
done

echo ""
echo "Analysis complete! Check indels_results/ directory"
