# CRISPR Indel Detection Pipeline

## Overview
This pipeline detects insertions and deletions (indels) at CRISPR-Cas9 cut sites from paired-end FASTQ sequencing data. It analyzes merged reads to identify non-homologous end joining (NHEJ) repair outcomes at target loci.

## Experimental Design

### Target genes
- **AAVS1**, **CEP290**, **TRAC**

### Guide RNAs
- **AAVS1**: `GTCCCTAGTGGCCCCACTGTNGG`
- **CEP290**: `GGAGTCACATGGGAGTCACANGG` 
- **TRAC**: `GCTGGTACACGGCAGGGTCANGG`

### ODN sequence
`GTTTAATTGAGTTGTCATATGTTAATAACGGTAT`

### Sequencing
Paired-end reads merged into single FASTQ files

## Pipeline Components

### 1. Reference Preparation

Reference sequences can be seen [here](genomic_regions.fa)

```bash
# Remove cut site markers from reference sequences
sed 's/-//g' genomic_regions.fa > genomic_regions_clean.fa
```

**Reference sequences**: ~200bp genomic regions (100bp flanking each side of cut site)
- Cut site positioned 3bp upstream of PAM (between positions 17-18 of protospacer)
- CEP290 and TRAC sequences in reverse complement orientation (genomic strand)

### 2. Indel Detection Algorithm

#### Cut Site Identification
- **Forward strand**: Cut occurs between bases 17|18 of guide sequence
- **Algorithm**: Searches for guide sequences within reference to locate precise cut positions
- **Fallback**: Uses sequence midpoint if guide not found

#### Read Analysis Strategy
1. **Pre-filtering**: Screen reads containing partial target region sequences (25bp flanks)
2. **Bidirectional alignment**: Test both forward and reverse complement orientations
3. **Sliding window approach**: Variable flank matching (15-50bp) for optimal alignment
4. **Indel classification**:
   - **Perfect match**: No editing detected
   - **Insertion**: Extra sequence at cut site
   - **Deletion**: Missing sequence at cut site

#### Alignment Scoring
- **Minimum flank match**: 15bp on each side
- **Score threshold**: ≥25 matching bases required
- **Maximum indel size**: 100bp (configurable)

### 3. Analysis Script Features

#### Input Validation
- File existence and permission checks
- FASTQ format validation
- Reference sequence parsing without BioPython dependency
- File size and content verification

#### Processing Capabilities
- **Gzipped FASTQ support**: Automatic detection and decompression
- **Manual parsing**: No external dependencies beyond Python standard library
- **Progress tracking**: Real-time processing updates
- **Error handling**: Comprehensive exception management

#### Output Metrics
- **Total reads processed**: All reads in FASTQ file
- **Analyzed reads**: Reads containing target region sequences
- **Coverage**: Percentage of reads analyzed per target
- **Editing efficiency**: Percentage of analyzed reads with indels
- **Event breakdown**: Detailed indel size and sequence distribution

### 4. Batch Processing Script

```bash
#!/bin/bash
# Process all sample directories matching pattern *_[12]
for DIR in *_[12]; do
    for FASTQ in ${DIR}/*_merged.fastq.gz; do
        python3 CRISPR_indel_detect.py \
            "$FASTQ" \
            genomic_regions_clean.fa \
            --output "indels_results/${DIR}_${BASENAME}_indels.txt" \
            --min-quality 20 \
            --max-indel 100
    done
done
```

## Sample Categories Analyzed

### Treatment Groups
- **Cas9_only**: Cas9 nuclease without ODN
- **ODN_only**: ODN template without Cas9
- **Cas9_ODN**: Combined Cas9 + ODN treatment
- **Cas9_ODN_10pmol**: High concentration ODN treatment
- **Cas9_ODN_old**: Previous ODN batch
- **Untreated**: Negative control

### Replicates
- Each treatment includes biological replicates (_1, _2)
- Each target gene analyzed separately (AAVS1, CEP290, TRAC)

## Technical Specifications

### Parameters
- **Minimum quality score**: 20 (configurable)
- **Maximum indel size**: 100bp
- **Minimum flank match**: 15bp
- **Alignment threshold**: 25 matching bases
- **Search window**: 25bp partial sequence matching

### File Structure
```
├── genomic_regions.fa          # Reference sequences with cut site markers
├── genomic_regions_clean.fa    # Cleaned reference sequences
├── *_[12]/                     # Sample directories
│   └── *_merged.fastq.gz       # Merged paired-end reads
├── indels_results/             # Output directory
│   └── *_indels.txt           # Analysis results per sample
└── CRISPR_indel_detect.py     # Main analysis script
```

## Output Format

### Summary Statistics
```
Total reads processed: XXX,XXX
Reads analyzed (containing target region): XX,XXX
Analysis coverage: XX.XX%
Editing efficiency: XX.XX%
```

### Event Classification
```
perfect_match: XXX (XX.XX%)
del_3bp: XXX (XX.XX%)
ins_1bp_A: XXX (XX.XX%)
del_15bp: XXX (XX.XX%)
```

## Quality Control Considerations

### Strand Orientation
- CEP290 and TRAC targets on genomic negative strand
- Algorithm handles both orientations automatically
- Reference sequences maintained in genomic coordinates

### Filtering Criteria
- Minimum alignment score prevents false positives
- Bidirectional search captures all orientations
- Quality score filtering removes low-quality reads
- Size limits focus on biologically relevant indels

### Validation Metrics
- Coverage percentage indicates target capture efficiency
- Perfect match rate shows baseline editing levels
- Event distribution reveals repair pathway preferences

---

This pipeline provides comprehensive quantification of CRISPR-induced indels across multiple targets and experimental conditions, enabling detailed analysis of nuclease activity and repair outcomes.
