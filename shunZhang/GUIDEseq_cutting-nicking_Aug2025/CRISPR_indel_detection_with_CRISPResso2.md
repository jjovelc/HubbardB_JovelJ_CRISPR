# CRISPResso2 Analysis Procedure Guide

Many approaches were tried. At the end using [CIRSPResso2](https://crispresso2.pinellolab.org/submission) pipeline was satisfactory.

## Overview

This guide describes the automated analysis procedure for detecting indels in CRISPR experiments using CRISPResso2. **All analyses described in this document were performed using the automated batch processing script `run_all_CRISPResso2.sh`** (see [Analysis Script](#analysis-script) section for complete code). The analysis targets three human genes: **AAVS1**, **CEP290**, and **TRAC** using optimized amplicon sequences determined through BLAST alignment validation.

## Prerequisites

### Software Requirements
- **CRISPResso2** (version 2.3.3 or higher)
- **Conda/Mamba** environment manager
- **Bash** shell environment

### Data Requirements
- **Merged FASTQ files**: Pre-processed paired-end reads that have been merged into single sequences
- **Directory structure**: Each sample should be in its own directory with the merged FASTQ file named as `{DIRECTORY_NAME}_merged.fastq.gz`

## Experimental Setup

### Target Genes and Sequences

### Target Genes and Sequences

**Note: All amplicon sequences were validated and optimized through BLAST alignment analysis to ensure accurate mapping of sequencing reads.**

#### 1. AAVS1 (ADENO-ASSOCIATED VIRUS INTEGRATION SITE 1)
- **Guide RNA**: `GTCCCTAGTGGCCCCACTGT`
- **Amplicon sequence**: `CAGATAAGGAATCTGCCTAACAGGAGGTGGGGGTTAGACCCAATATCAGGAGACTAGGAAGGAGGAGGCCTAAGGATGGGGCTTTTCTGTCACCAATCCTGTCCCTAGTGGCCCCACTGTGGGGTGGAGGGGACAGATAAAAGTACCCAGAACCAGAGCCACATTAACCGGCCCTGGGAATATAAGGTGGTCCCAGCTCGGGGACACAGGATCCCTGGAGGCA`
- **Amplicon length**: 253 bp
- **Function**: Safe harbor locus for gene integration

#### 2. CEP290 (CENTROSOMAL PROTEIN 290)
- **Guide RNA**: `GGAGTCACATGGGAGTCACA`  
- **Amplicon sequence**: `CTCAGCTAGTTCATCTCTTGCTCTAGATGACATGAGGTAAGTAGGGGACTTGACTTTTACCCTTCAGGTAACCGGTAGCTTTTGACAGTTTTTAAGGCGGGGAGTCACATGGGAGTCACAGGGTAGGATTCATGTTTAGAATGATCATTCTTGTGGCAGTAAGGAGGATGTAAGACTGGAGATAGAGACAGGAATAATGGCTGCCACAATAAAGATAATAAAA`
- **Amplicon length**: 256 bp
- **Function**: Centrosomal protein involved in ciliary function

#### 3. TRAC (T-CELL RECEPTOR ALPHA CONSTANT)
- **Guide RNA**: `GCTGGTACACGGCAGGGTCA`
- **Amplicon sequence**: `GTGATATACACATCAGAATCCTTACTTTGTGACACATTTGTTTGAGAATCAAAATCGGTGAATAGGCAGACAGACTTGTCACTGGATTTAGAGTCTCTCAGCTGGTACACGGCAGGGTCAGGGTTCTGGATATCTGTGGGACAAGAGGATCAGGGTTAGGACATGATCTCATTTCCCTCTTTGCCCCAACCCAGGCTGGAGTCCAGATGCCAGTGATGGACAA`
- **Amplicon length**: 244 bp
- **Function**: T-cell receptor constant region

## Analysis Parameters

**All parameters listed below are implemented in the automated analysis script `run_all_CRISPResso2.sh`.**

### Standard Parameters (AAVS1 & CEP290)
```bash
--quantification_window_size 20      # Analysis window around cut site
--exclude_bp_from_left 15           # Exclude 15bp from left end (primers)
--exclude_bp_from_right 15          # Exclude 15bp from right end (primers)
--min_frequency_alleles_around_cut_to_plot 0.05  # Show alleles ≥5% frequency
--plot_histogram_outliers           # Include large indels in plots (-50 to +50 bp range)
```

### Special Parameters (TRAC)
Due to the similar length but different sequence characteristics of the TRAC amplicon:
```bash
--quantification_window_size 10      # Smaller window for optimized analysis
--exclude_bp_from_left 0            # No trimming (optimized for sequence)
--exclude_bp_from_right 0           # No trimming (optimized for sequence)
```

## Analysis Script

**The complete automated analysis was performed using the following bash script (`run_all_CRISPResso2.sh`):**

```bash
#!/usr/bin/env bash
set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate crispresso2

# Define amplicon sequences (validated through BLAST alignment)
AAVS1_amplicon_seq='CAGATAAGGAATCTGCCTAACAGGAGGTGGGGGTTAGACCCAATATCAGGAGACTAGGAAGGAGGAGGCCTAAGGATGGGGCTTTTCTGTCACCAATCCTGTCCCTAGTGGCCCCACTGTGGGGTGGAGGGGACAGATAAAAGTACCCAGAACCAGAGCCACATTAACCGGCCCTGGGAATATAAGGTGGTCCCAGCTCGGGGACACAGGATCCCTGGAGGCA'
CEP290_amplicon_seq='CTCAGCTAGTTCATCTCTTGCTCTAGATGACATGAGGTAAGTAGGGGACTTGACTTTTACCCTTCAGGTAACCGGTAGCTTTTGACAGTTTTTAAGGCGGGGAGTCACATGGGAGTCACAGGGTAGGATTCATGTTTAGAATGATCATTCTTGTGGCAGTAAGGAGGATGTAAGACTGGAGATAGAGACAGGAATAATGGCTGCCACAATAAAGATAATAAAA'
TRAC_amplicon_seq='GTGATATACACATCAGAATCCTTACTTTGTGACACATTTGTTTGAGAATCAAAATCGGTGAATAGGCAGACAGACTTGTCACTGGATTTAGAGTCTCTCAGCTGGTACACGGCAGGGTCAGGGTTCTGGATATCTGTGGGACAAGAGGATCAGGGTTAGGACATGATCTCATTTCCCTCTTTGCCCCAACCCAGGCTGGAGTCCAGATGCCAGTGATGGACAA'

# Check arguments
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 DIRECTORY [DIRECTORY2 DIRECTORY3 ...]" >&2
  echo "Example: $0 sample1_AAVS1 sample2_CEP290 sample3_TRAC" >&2
  exit 1
fi

# Process each directory argument
for DIR in "$@"; do
  echo "Processing directory: $DIR"
  
  # Extract gene name from directory name
  if [[ "$DIR" =~ (AAVS1|CEP290|TRAC) ]]; then
    GENE="${BASH_REMATCH[1]}"
  else
    echo ">> Skipping $DIR (no gene tag found in directory name)."
    continue
  fi
  
  # Set amplicon and guide sequences based on gene
  case "$GENE" in
    AAVS1) 
      AMPLICON_SEQ="$AAVS1_amplicon_seq"
      GUIDE_SEQ="GTCCCTAGTGGCCCCACTGT"
      ;;
    CEP290) 
      AMPLICON_SEQ="$CEP290_amplicon_seq"
      GUIDE_SEQ="GGAGTCACATGGGAGTCACA"
      ;;
    TRAC)  
      AMPLICON_SEQ="$TRAC_amplicon_seq"
      GUIDE_SEQ="GCTGGTACACGGCAGGGTCA"
      ;;
    *)     
      echo "!! Unknown gene '$GENE' in $DIR"
      continue 
      ;;
  esac
      
  # Define input file path
  INFILE="${DIR}/${DIR}_merged.fastq.gz"
  
  # Check if input file exists
  if [[ -f "$INFILE" ]]; then
    echo "Submitting: $INFILE  gene=$GENE"
    echo "Amplicon length: ${#AMPLICON_SEQ} bp"
    echo "Guide sequence: $GUIDE_SEQ"
    
    # Run CRISPResso analysis with gene-specific parameters
    if [[ "$GENE" == "TRAC" ]]; then
      # Special parameters for TRAC amplicon
      CRISPResso --fastq_r1 "$INFILE" \
                --amplicon_seq "$AMPLICON_SEQ" \
                --guide_seq "$GUIDE_SEQ" \
                --output_folder "${DIR}_${GENE}_CRISPResso2_analysis" \
                --quantification_window_size 10 \
                --min_frequency_alleles_around_cut_to_plot 0.05 \
                --plot_histogram_outliers \
                --exclude_bp_from_left 0 \
                --exclude_bp_from_right 0 \
                --name "${DIR}_${GENE}"
    else
      # Standard parameters for AAVS1 and CEP290
      CRISPResso --fastq_r1 "$INFILE" \
                --amplicon_seq "$AMPLICON_SEQ" \
                --guide_seq "$GUIDE_SEQ" \
                --output_folder "${DIR}_${GENE}_CRISPResso2_analysis" \
                --quantification_window_size 20 \
                --min_frequency_alleles_around_cut_to_plot 0.05 \
                --plot_histogram_outliers \
                --exclude_bp_from_left 15 \
                --exclude_bp_from_right 15 \
                --name "${DIR}_${GENE}"
    fi
    
    echo "Completed analysis for $DIR ($GENE)"
    echo "Results saved in: ${DIR}_${GENE}_CRISPResso2_analysis/"
    echo "---"
    
  else
    echo "WARN: Missing FASTQ ${INFILE} — skipping"
  fi
  
done

echo "All analyses completed!"
```

### Script Usage Examples
```bash
# Make script executable
chmod +x run_all_CRISPResso2.sh

# Analyze single sample
./run_all_CRISPResso2.sh Cas9_ODN_old_AAVS1_1

# Analyze multiple samples
./run_all_CRISPResso2.sh Cas9_ODN_old_AAVS1_1 Cas9_ODN_old_CEP290_1 Cas9_ODN_old_TRAC_1

# Use wildcards for batch processing
./run_all_CRISPResso2.sh *AAVS1* *CEP290* *TRAC*
```

### Step 1: Environment Setup
```bash
# Activate the CRISPResso2 conda environment
eval "$(conda shell.bash hook)"
conda activate crispresso2
```

### Step 2: Directory Structure Validation
The script expects the following directory structure:
```
project_directory/
├── run_all_CRISPResso2.sh
├── Sample1_AAVS1/
│   └── Sample1_AAVS1_merged.fastq.gz
├── Sample2_CEP290/
│   └── Sample2_CEP290_merged.fastq.gz
└── Sample3_TRAC/
    └── Sample3_TRAC_merged.fastq.gz
```

### Step 3: Gene Detection and Parameter Assignment
The script automatically:
1. **Parses directory names** to identify target genes (AAVS1, CEP290, or TRAC)
2. **Assigns appropriate amplicon sequences** based on gene identity
3. **Sets gene-specific guide RNA sequences**
4. **Configures analysis parameters** optimized for each gene

### Step 4: Quality Control Checks
Before analysis, the script performs:
- **File existence verification**: Confirms merged FASTQ files are present
- **Parameter validation**: Displays amplicon length and guide sequence
- **Gene identification confirmation**: Shows detected gene for each directory

### Step 5: CRISPResso2 Analysis Execution
For each valid sample, the script automatically runs the appropriate CRISPResso2 command with gene-specific parameters:

**For AAVS1 and CEP290:**
```bash
CRISPResso --fastq_r1 INPUT_FILE \
          --amplicon_seq REFERENCE_SEQUENCE \
          --guide_seq GUIDE_RNA_SEQUENCE \
          --output_folder OUTPUT_DIRECTORY \
          --quantification_window_size 20 \
          --exclude_bp_from_left 15 \
          --exclude_bp_from_right 15 \
          --plot_histogram_outliers \
          --min_frequency_alleles_around_cut_to_plot 0.05
```

**For TRAC (optimized parameters):**
```bash
CRISPResso --fastq_r1 INPUT_FILE \
          --amplicon_seq REFERENCE_SEQUENCE \
          --guide_seq GUIDE_RNA_SEQUENCE \
          --output_folder OUTPUT_DIRECTORY \
          --quantification_window_size 10 \
          --exclude_bp_from_left 0 \
          --exclude_bp_from_right 0 \
          --plot_histogram_outliers \
          --min_frequency_alleles_around_cut_to_plot 0.05
```

### Step 6: Output Generation
Each analysis generates:
- **HTML reports**: Interactive visualization of results
- **Statistical summaries**: Editing efficiency and indel frequencies
- **Plots**: Indel size distributions, nucleotide frequency heatmaps
- **Data tables**: Detailed allele information for downstream analysis

## Usage Instructions

**All usage instructions apply to the `run_all_CRISPResso2.sh` script provided in the [Analysis Script](#analysis-script) section.**

### Basic Usage
```bash
# Analyze single sample
./run_all_CRISPResso2.sh Sample1_AAVS1

# Analyze multiple samples
./run_all_CRISPResso2.sh Sample1_AAVS1 Sample2_CEP290 Sample3_TRAC

# Use wildcards for batch processing
./run_all_CRISPResso2.sh *_AAVS1 *_CEP290 *_TRAC
```

### Make Script Executable
```bash
chmod +x run_all_CRISPResso2.sh
```

## Output Structure

### Generated Directories
For each input sample, the analysis creates:
```
{SAMPLE_NAME}_{GENE}_CRISPResso2_analysis/
├── CRISPResso2_report.html                    # Main interactive report
├── CRISPResso_quantification_of_editing_frequency.txt
├── Alleles_frequency_table.txt                # Detailed allele data
├── Indel_histogram.txt                        # Indel size distribution data
├── 3a.Indel_size_distribution.pdf            # Indel size plots (-50 to +50 bp)
└── [Additional plots and analysis files]
```

### Key Output Files

#### 1. Interactive HTML Report
- **File**: `CRISPResso2_report.html`
- **Contents**: Complete analysis summary with interactive plots
- **Usage**: Open in web browser for comprehensive results review

#### 2. Quantification Summary
- **File**: `CRISPResso_quantification_of_editing_frequency.txt`
- **Contents**: Overall editing statistics
- **Key metrics**: 
  - Total editing efficiency (%)
  - Frameshift frequency (%)
  - In-frame indel frequency (%)

#### 3. Allele Frequency Table
- **File**: `Alleles_frequency_table.txt`
- **Contents**: Detailed information for each detected allele
- **Includes**: Sequence, frequency, indel type, and size

#### 4. Indel Size Distribution
- **File**: `Indel_histogram.txt`
- **Contents**: Raw data for indel size plotting
- **Range**: Extended range (-50 to +50 bp) due to `--plot_histogram_outliers`

## Analysis Interpretation

### Key Metrics to Evaluate

1. **Overall Editing Efficiency**
   - Target: >10% for successful editing
   - Calculation: (Reads with indels / Total aligned reads) × 100

2. **Indel Size Distribution**
   - **Small deletions** (-1 to -10 bp): Most common NHEJ products
   - **Large deletions** (>-10 bp): Less frequent but significant
   - **Insertions** (+1 to +10 bp): Typically less frequent than deletions

3. **Cut Site Precision**
   - Verify indels cluster around expected cut site
   - Check for off-target cutting patterns

4. **Frameshift Analysis**
   - Important for knockout experiments
   - Indels not divisible by 3 cause frameshifts

### Quality Control Indicators

#### Good Results
- Editing efficiency >10%
- Clear indel peak around expected cut site
- Low background in control samples
- Consistent results across replicates

#### Problematic Results
- Very low editing efficiency (<2%)
- Indels distributed randomly across amplicon
- High background in controls
- Poor read alignment rates

## Troubleshooting

### Common Issues and Solutions

1. **"Quantification window excluded" Error**
   - **Cause**: Amplicon too short relative to trimming parameters
   - **Solution**: Adjust `--exclude_bp_from_left/right` parameters (set to 0 for short amplicons)

2. **Low Alignment Rate**
   - **Cause**: Amplicon sequence doesn't match actual PCR product
   - **Solution**: Verify amplicon sequence with BLAST alignment

3. **No Indels Detected**
   - **Cause**: Guide RNA sequence incorrect or experiment failed
   - **Solution**: Verify guide RNA sequence and check experimental conditions

4. **Memory/Performance Issues**
   - **Cause**: Large FASTQ files or insufficient resources
   - **Solution**: Use `--fastq_r1` parameter optimization or increase system resources

## Advanced Options

### HDR Analysis
If using homology-directed repair (HDR) templates:
```bash
--expected_hdr_amplicon_seq HDR_TEMPLATE_SEQUENCE
```

### Custom Quantification Windows
For specific analysis requirements:
```bash
--quantification_window_coordinates START:END
```

### Batch Processing Optimization
For large datasets:
```bash
--n_processes 4  # Parallel processing
--cache_size 1000  # Memory management
```

## Data Export and Downstream Analysis

### Extracting Summary Statistics
```bash
# Combine results from multiple samples
grep "Total" */CRISPResso_quantification_of_editing_frequency.txt > combined_stats.txt
```

### Custom Plotting
Use `Indel_histogram.txt` files for custom visualization:
```python
import pandas as pd
import matplotlib.pyplot as plt

# Load indel data
df = pd.read_csv('Indel_histogram.txt', sep='\t')

# Create custom plots with extended range (-50 to +50 bp)
plt.figure(figsize=(12,8))
plt.bar(df['Size'], df['Percentage'])
plt.xlim(-50, 50)
plt.xlabel('Indel Size (bp)')
plt.ylabel('Frequency (%)')
plt.title('CRISPR Indel Distribution')
plt.show()
```

## Best Practices

1. **Always include control samples** (untreated cells) for background assessment
2. **Use biological replicates** (minimum n=3) for statistical significance
3. **Verify amplicon sequences** against reference genome before analysis
4. **Check FASTQ quality** before running CRISPResso2
5. **Save analysis parameters** for reproducibility
6. **Review HTML reports** for comprehensive result interpretation

## References

- CRISPResso2: [https://github.com/pinellolab/CRISPResso2](https://github.com/pinellolab/CRISPResso2)
- Original publication: Clement et al. Nature Biotechnology (2019)
- Documentation: [https://crispresso.pinellolab.partners.org/](https://crispresso.pinellolab.partners.org/)
