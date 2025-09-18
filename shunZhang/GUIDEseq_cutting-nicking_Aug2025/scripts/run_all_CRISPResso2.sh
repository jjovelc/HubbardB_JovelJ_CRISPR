#!/usr/bin/env bash
set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate crispresso2

# Define amplicon sequences
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
    
    # Run CRISPResso analysis with appropriate parameters for your TRAC sequence
    if [[ "$GENE" == "TRAC" ]]; then
      # Special parameters for shorter TRAC amplicon
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

# Usage examples:
# ./run_all_CRISPResso2.sh Cas9_ODN_old_AAVS1_1
# ./run_all_CRISPResso2.sh Cas9_ODN_old_AAVS1_1 Cas9_ODN_old_CEP290_1 Cas9_ODN_old_TRAC_1
# ./run_all_CRISPResso2.sh *AAVS1* *CEP290* *TRAC*
