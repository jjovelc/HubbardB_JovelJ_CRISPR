# GUIDE-seq Cutting and Nicking Analysis - August 2025

## Sample Metadata
This project contains **164 samples** with comprehensive barcode and index information. For the complete sample metadata table including all barcodes, I5 sequences, and I7 sequences, see: **[Sample Metadata](sample_metadata.md)**

**Sample Summary:**
- **Control samples:** Untreated, ODN only, Cas9 only
- **Experimental samples:** Cas9+ODN, nickase library (nlib) with XNA treatments
- **Targets:** CEP290, AAVS1, TRAC
- **XNA types:** RNA, FANA, LNA, DNA
- **Concentrations:** 0.5 pmol, 5 pmol, 10 pmol, 100 pmol

## Analysis Protocol
For detailed analysis instructions and protocols, see: <a href="nicking_lib_NGS_analysis_protocol_2024Dec_inline.html" download>
  NGS Analysis Protocol
</a>

**Note:** To view this protocol properly:
1. **Right-click the link above** and select "Save link as..." or "Download linked file"
2. **Save the .html file** to your computer
3. **Open the saved file** in your web browser for proper formatting

Alternatively, you can view the file directly on GitHub and click the "Download" button.

## Overview
This repository contains analysis scripts and documentation for GUIDE-seq experiments investigating Cas9 cutting efficiency, ODN tag incorporation, and Cas9 nickase high-throughput in vitro library analysis.

## Sample Overview
**Total Samples:** 164 samples divided into two experimental groups

## Experimental Design

### Group 1: Cas9 Cutting and ODN Tag Incorporation (Samples 1-30, 153-164)
**Objective:** Measure Cas9 cutting efficiency and ODN tag incorporation rates at three genomic targets:
- **CEP290**
- **AAVS1** 
- **TRAC**

**Sample Details:**
- **Samples 1-30:** Standard ODN concentration (100 pmol)
- **Samples 153-164:** Reduced ODN concentration (10 pmol) - used to investigate concentration effects on tag incorporation

**Experimental Protocol:**
1. Cas9/ODN treatment of target cells
2. Two rounds of PCR amplification
3. Barcoding with unique dual indices (P5 & P7)
4. Company demultiplexing completed
5. Ready for analysis

**Controls:**
- Untreated samples
- ODN-only samples (no Cas9)

**Analysis Focus:**
- Cas9 cutting efficiency across different samples
- ODN incorporation rate assessment
- Troubleshooting ODN incorporation effectiveness

### Group 2: Cas9 Nickase High-Throughput In Vitro Library (Samples 31-152)
**Objective:** Analyze mismatch patterns and XNA treatment effects on Cas9 nickase activity

**Sample Details:**
- **Samples 151-152:** Pre-selection libraries (no XNA treatment) - **control samples**
- **Samples 31-150:** XNA-treated libraries at two dosages:
  - 5 pmol XNA
  - 0.5 pmol XNA

**Analysis Protocol:**
- Same protocol as previous experiments
- Targeted amplicon sequencing
- Barcoding with unique dual indices

## Data Access
**Google Drive Link:** [GUIDE-seq Data](https://drive.google.com/drive/folders/14hpWdR1cHqZE3kIW0SOeMU0pWVYMw35i?usp=sharing)

The drive contains all sample FASTQ files described above, ready for analysis.

## Key Research Questions
1. **ODN Incorporation Efficiency:** How effectively does ODN incorporate into target sites?
2. **Concentration Effects:** Does ODN concentration (100 pmol vs 10 pmol) affect tag incorporation patterns?
3. **XNA Treatment Effects:** How do different XNA treatments and dosages affect Cas9 nickase mismatch patterns?
4. **Cutting Efficiency:** What is the Cas9 cutting efficiency across different samples and targets?

## Notes
- Company has completed demultiplexing of samples
- Analysis protocols from previous experiments should be referenced
- Focus on troubleshooting ODN incorporation effectiveness for assay optimization
