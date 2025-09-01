# GUIDE-seq cutting & ODN incorporation rate

## 1. Goal of the Group-1 pipeline

- Quantify, per sample and target (CEP290/AAVS1/TRAC):
- Cutting efficiency: indel % in a fixed window around the SpCas9 cut site (from amplicon reads).
- ODN tag incorporation rate: fraction of reads containing the ODN motif within the extracted insert.
- QC: read quality and successful recovery of the insert between the two assay constant sequences.

## 2. Plain-text, step-by-step pipeline

### Inputs
- Demultiplexed R1/R2 FASTQs (one directory per sample).
- The two constant sequences (left/right, R1 orientation).
- Per-target amplicon sequence and sgRNA (20 nt; no PAM) for CEP290, AAVS1, TRAC.
- A small manifest (TSV) mapping each sample folder → target + paths to R1/R2 (the SLURM script below can auto-build this from your joined_samples_with_barcodes.csv and your FASTQ root).

### QC and trimming

- Run fastqc on R1/R2; keep reports per sample.
- Insert extraction (constant→constant)
- Use cutadapt linked adapters to keep only the sequence between LEFT ... RIGHT, removing the constants.
- Allow modest error rate (e.g., -e 0.10) and minimum overlap (-O 16).
- Discard reads that don’t contain both constants (--discard-untrimmed).

### Output: 
- extracted_R1.fastq.gz (and optionally extracted_R2.fastq.gz).
- ODN tag detection
- Scan extracted_R1.fastq.gz for your ODN IUPAC motif (e.g., NNCC...TCGTCTNN):
- Strict: regex with IUPAC expansion (N → [ACGT]).
- Lenient: allow ≤1 mismatch at fixed (non-N) positions.
- Report: total, strict, lenient, and rates.
- Cutting efficiency (indel %)
- Run CRISPResso2 (amplicon mode) on extracted_R1.fastq.gz with the correct amplicon_seq and guide_seq.
- Use a ±5 bp (or your preference) window around the expected cut site for quantification.
- Parse its quantification table to obtain % indel.

### Outputs & roll-up

- Per sample: FastQC, cutadapt logs, odn_counts.tsv/json, CRISPResso folder.
- Optional: merge per-sample metrics later into a CSV.
