# GUIDE-seq Cutting and Nicking Analysis - August 2025

## Sample Metadata
This project contains 164 samples with the following barcode and index information:

| GUIDE-seq | barcodes (I5-I7) | I5 | I7 (non-reverse) | I7 (reverse) |
| :-------------------------- | :----------------- | :--------- | :--------------- | :----------- |
| Untreated-CEP290-1 | 1-1 | CCTGAACC | GAGATAAC | GTTATCTC |
| Untreated-CEP290-2 | 1-2 | CCTGAACC | CATGGAAT | ATTCCATG |
| Untreated-AAVS1-1 | 1-3 | CCTGAACC | CAGTTCCA | TGGAACTG |
| Untreated-AAVS1-2 | 1-4 | CCTGAACC | TCTCGCCT | AGGCGAGA |
| Untreated-TRAC-1 | 1-5 | CCTGAACC | CGCAAGCT | AGCTTGCG |
| Untreated-TRAC-2 | 1-6 | CCTGAACC | ATAGAGTC | GACTCTAT |
| Cas9+ODN-CEP290-1 | 2-7 | TATCCAGT | GGCTCTTG | CAAGAGCC |
| Cas9+ODN-CEP290-2 | 2-8 | TATCCAGT | TCGCCAGA | TCTGGCGA |
| Cas9+ODN-AAVS1-1 | 2-9 | TATCCAGT | TTACCGGG | CCCGGTAA |
| Cas9+ODN-AAVS1-2 | 2-10 | TATCCAGT | GCAGATAA | TTATCTGC |
| Cas9+ODN-TRAC-1 | 2-11 | TATCCAGT | AGTTGCTG | CAGCAACT |
| Cas9+ODN-TRAC-2 | 2-12 | TATCCAGT | ATCATTGC | GCAATGAT |
| ODN only-CEP290-1 | 3-1 | GAGATAAC | GAGATAAC | GTTATCTC |
| ODN only-CEP290-2 | 3-2 | GAGATAAC | CATGGAAT | ATTCCATG |
| ODN only-AAVS1-1 | 3-3 | GAGATAAC | CAGTTCCA | TGGAACTG |
| ODN only-AAVS1-2 | 3-4 | GAGATAAC | TCTCGCCT | AGGCGAGA |
| ODN only-TRAC-1 | 3-5 | GAGATAAC | CGCAAGCT | AGCTTGCG |
| ODN only-TRAC-2 | 3-6 | GAGATAAC | ATAGAGTC | GACTCTAT |
| Cas9 only-CEP290-1 | 3-7 | GAGATAAC | GGCTCTTG | CAAGAGCC |
| Cas9 only-CEP290-2 | 3-8 | GAGATAAC | TCGCCAGA | TCTGGCGA |
| Cas9 only-AAVS1-1 | 3-9 | GAGATAAC | TTACCGGG | CCCGGTAA |
| Cas9 only-AAVS1-2 | 3-10 | GAGATAAC | GCAGATAA | TTATCTGC |
| Cas9 only-TRAC-1 | 3-11 | GAGATAAC | AGTTGCTG | CAGCAACT |
| Cas9 only-TRAC-2 | 3-12 | GAGATAAC | ATCATTGC | GCAATGAT |
| Cas9+ODN old-CEP290-1 | 3-13 | GAGATAAC | GAGATAAC | GTTATCTC |
| Cas9+ODN old-CEP290-2 | 3-14 | GAGATAAC | CATGGAAT | ATTCCATG |
| Cas9+ODN old-AAVS1-1 | 3-15 | GAGATAAC | CAGTTCCA | TGGAACTG |
| Cas9+ODN old-AAVS1-2 | 3-16 | GAGATAAC | TCTCGCCT | AGGCGAGA |
| Cas9+ODN old-TRAC-1 | 3-17 | GAGATAAC | CGCAAGCT | AGCTTGCG |
| Cas9+ODN old-TRAC-2 | 3-18 | GAGATAAC | ATAGAGTC | GACTCTAT |
| nlib-RNA.1-5pmol-1 | 4-1 | ACACAGGC | GAGATAAC | GTTATCTC |
| nlib-RNA.1-5pmol-2 | 4-2 | ACACAGGC | CATGGAAT | ATTCCATG |
| nlib-RNA.1-0.5pmol-1 | 4-3 | ACACAGGC | CAGTTCCA | TGGAACTG |
| nlib-RNA.1-0.5pmol-2 | 4-4 | ACACAGGC | TCTCGCCT | AGGCGAGA |
| nlib-RNA.2-5pmol-1 | 4-5 | ACACAGGC | CGCAAGCT | AGCTTGCG |
| nlib-RNA.2-5pmol-2 | 4-6 | ACACAGGC | ATAGAGTC | GACTCTAT |
| nlib-RNA.2-0.5pmol-1 | 4-7 | ACACAGGC | GGCTCTTG | CAAGAGCC |
| nlib-RNA.2-0.5pmol-2 | 4-8 | ACACAGGC | TCGCCAGA | TCTGGCGA |
| nlib-RNA.3-5pmol-1 | 4-9 | ACACAGGC | TTACCGGG | CCCGGTAA |
| nlib-RNA.3-5pmol-2 | 4-10 | ACACAGGC | GCAGATAA | TTATCTGC |
| nlib-RNA.3-0.5pmol-1 | 4-11 | ACACAGGC | AGTTGCTG | CAGCAACT |
| nlib-RNA.3-0.5pmol-2 | 4-12 | ACACAGGC | ATCATTGC | GCAATGAT |
| nlib-FANA.1-5pmol-1 | 5-1 | CCGTATAT | GAGATAAC | GTTATCTC |
| nlib-FANA.1-5pmol-2 | 5-2 | CCGTATAT | CATGGAAT | ATTCCATG |
| nlib-FANA.1-0.5pmol-1 | 5-3 | CCGTATAT | CAGTTCCA | TGGAACTG |
| nlib-FANA.1-0.5pmol-2 | 5-4 | CCGTATAT | TCTCGCCT | AGGCGAGA |
| nlib-FANA.2-5pmol-1 | 5-5 | CCGTATAT | CGCAAGCT | AGCTTGCG |
| nlib-FANA.2-5pmol-2 | 5-6 | CCGTATAT | ATAGAGTC | GACTCTAT |
| nlib-FANA.2-0.5pmol-1 | 5-7 | CCGTATAT | GGCTCTTG | CAAGAGCC |
| nlib-FANA.2-0.5pmol-2 | 5-8 | CCGTATAT | TCGCCAGA | TCTGGCGA |
| nlib-FANA.3-5pmol-1 | 5-9 | CCGTATAT | TTACCGGG | CCCGGTAA |
| nlib-FANA.3-5pmol-2 | 5-10 | CCGTATAT | GCAGATAA | TTATCTGC |
| nlib-FANA.3-0.5pmol-1 | 5-11 | CCGTATAT | AGTTGCTG | CAGCAACT |
| nlib-FANA.3-0.5pmol-2 | 5-12 | CCGTATAT | ATCATTGC | GCAATGAT |
| nlib-LNA.1-5pmol-1 | 6-1 | GCTGAAGA | GAGATAAC | GTTATCTC |
| nlib-LNA.1-5pmol-2 | 6-2 | GCTGAAGA | CATGGAAT | ATTCCATG |
| nlib-LNA.1-0.5pmol-1 | 6-3 | GCTGAAGA | CAGTTCCA | TGGAACTG |
| nlib-LNA.1-0.5pmol-2 | 6-4 | GCTGAAGA | TCTCGCCT | AGGCGAGA |
| nlib-LNA.2-5pmol-1 | 6-5 | GCTGAAGA | CGCAAGCT | AGCTTGCG |
| nlib-LNA.2-5pmol-2 | 6-6 | GCTGAAGA | ATAGAGTC | GACTCTAT |
| nlib-LNA.2-0.5pmol-1 | 6-7 | GCTGAAGA | GGCTCTTG | CAAGAGCC |
| nlib-LNA.2-0.5pmol-2 | 6-8 | GCTGAAGA | TCGCCAGA | TCTGGCGA |
| nlib-LNA.3-5pmol-1 | 6-9 | GCTGAAGA | TTACCGGG | CCCGGTAA |
| nlib-LNA.3-5pmol-2 | 6-10 | GCTGAAGA | GCAGATAA | TTATCTGC |
| nlib-LNA.3-0.5pmol-1 | 6-11 | GCTGAAGA | AGTTGCTG | CAGCAACT |
| nlib-LNA.3-0.5pmol-2 | 6-12 | GCTGAAGA | ATCATTGC | GCAATGAT |
| nlib-LNA.4-5pmol-1 | 6-13 | GCTGAAGA | GAGATAAC | GTTATCTC |
| nlib-LNA.4-5pmol-2 | 6-14 | GCTGAAGA | CATGGAAT | ATTCCATG |
| nlib-LNA.4-0.5pmol-1 | 6-15 | GCTGAAGA | CAGTTCCA | TGGAACTG |
| nlib-LNA.4-0.5pmol-2 | 6-16 | GCTGAAGA | TCTCGCCT | AGGCGAGA |
| nlib-LNA.5-5pmol-1 | 6-17 | GCTGAAGA | CGCAAGCT | AGCTTGCG |
| nlib-LNA.5-5pmol-2 | 6-18 | GCTGAAGA | ATAGAGTC | GACTCTAT |
| nlib-LNA.5-0.5pmol-1 | 6-19 | GCTGAAGA | GGCTCTTG | CAAGAGCC |
| nlib-LNA.5-0.5pmol-2 | 6-20 | GCTGAAGA | TCGCCAGA | TCTGGCGA |
| nlib-LNA.6-5pmol-1 | 6-21 | GCTGAAGA | TTACCGGG | CCCGGTAA |
| nlib-LNA.6-5pmol-2 | 6-22 | GCTGAAGA | GCAGATAA | TTATCTGC |
| nlib-LNA.6-0.5pmol-1 | 6-23 | GCTGAAGA | AGTTGCTG | CAGCAACT |
| nlib-LNA.6-0.5pmol-2 | 6-24 | GCTGAAGA | ATCATTGC | GCAATGAT |
| nlib-LNA.7-5pmol-1 | 6-25 | GCTGAAGA | GAGATAAC | GTTATCTC |
| nlib-LNA.7-5pmol-2 | 6-26 | GCTGAAGA | CATGGAAT | ATTCCATG |
| nlib-LNA.7-0.5pmol-1 | 6-27 | GCTGAAGA | CAGTTCCA | TGGAACTG |
| nlib-LNA.7-0.5pmol-2 | 6-28 | GCTGAAGA | TCTCGCCT | AGGCGAGA |
| nlib-LNA.8-5pmol-1 | 6-29 | GCTGAAGA | CGCAAGCT | AGCTTGCG |
| nlib-LNA.8-5pmol-2 | 6-30 | GCTGAAGA | ATAGAGTC | GACTCTAT |
| nlib-LNA.8-0.5pmol-1 | 6-31 | GCTGAAGA | GGCTCTTG | CAAGAGCC |
| nlib-LNA.8-0.5pmol-2 | 6-32 | GCTGAAGA | TCGCCAGA | TCTGGCGA |
| nlib-LNA.9-5pmol-1 | 6-33 | GCTGAAGA | TTACCGGG | CCCGGTAA |
| nlib-LNA.9-5pmol-2 | 6-34 | GCTGAAGA | GCAGATAA | TTATCTGC |
| nlib-LNA.9-0.5pmol-1 | 6-35 | GCTGAAGA | AGTTGCTG | CAGCAACT |
| nlib-LNA.9-0.5pmol-2 | 6-36 | GCTGAAGA | ATCATTGC | GCAATGAT |
| nlib-LNA.10-5pmol-1 | 6-37 | GCTGAAGA | GAGATAAC | GTTATCTC |
| nlib-LNA.10-5pmol-2 | 6-38 | GCTGAAGA | CATGGAAT | ATTCCATG |
| nlib-LNA.10-0.5pmol-1 | 6-39 | GCTGAAGA | CAGTTCCA | TGGAACTG |
| nlib-LNA.10-0.5pmol-2 | 6-40 | GCTGAAGA | TCTCGCCT | AGGCGAGA |
| nlib-LNA.11-5pmol-1 | 6-41 | GCTGAAGA | CGCAAGCT | AGCTTGCG |
| nlib-LNA.11-5pmol-2 | 6-42 | GCTGAAGA | ATAGAGTC | GACTCTAT |
| nlib-LNA.11-0.5pmol-1 | 6-43 | GCTGAAGA | GGCTCTTG | CAAGAGCC |
| nlib-LNA.11-0.5pmol-2 | 6-44 | GCTGAAGA | TCGCCAGA | TCTGGCGA |
| nlib-LNA.12-5pmol-1 | 6-45 | GCTGAAGA | TTACCGGG | CCCGGTAA |
| nlib-LNA.12-5pmol-2 | 6-46 | GCTGAAGA | GCAGATAA | TTATCTGC |
| nlib-LNA.12-0.5pmol-1 | 6-47 | GCTGAAGA | AGTTGCTG | CAGCAACT |
| nlib-LNA.12-0.5pmol-2 | 6-48 | GCTGAAGA | ATCATTGC | GCAATGAT |
| nlib-LNA.13-5pmol-1 | 6-49 | GCTGAAGA | GAGATAAC | GTTATCTC |
| nlib-LNA.13-5pmol-2 | 6-50 | GCTGAAGA | CATGGAAT | ATTCCATG |
| nlib-LNA.13-0.5pmol-1 | 6-51 | GCTGAAGA | CAGTTCCA | TGGAACTG |
| nlib-LNA.13-0.5pmol-2 | 6-52 | GCTGAAGA | TCTCGCCT | AGGCGAGA |
| preselection-lib-1 | 8-1 | CGCAAGCT | GAGATAAC | GTTATCTC |
| preselection-lib-2 | 8-2 | CGCAAGCT | CATGGAAT | ATTCCATG |
| ODN 10pmol-CEP290-1 | 14-1 | ACCTTGAA | GAGATAAC | GTTATCTC |
| ODN 10pmol-CEP290-2 | 14-2 | ACCTTGAA | CATGGAAT | ATTCCATG |
| ODN 10pmol-AAVS1-1 | 14-3 | ACCTTGAA | CAGTTCCA | TGGAACTG |
| ODN 10pmol-AAVS1-2 | 14-4 | ACCTTGAA | TCTCGCCT | AGGCGAGA |
| ODN 10pmol-TRAC-1 | 14-5 | ACCTTGAA | CGCAAGCT | AGCTTGCG |
| ODN 10pmol-TRAC-2 | 14-6 | ACCTTGAA | ATAGAGTC | GACTCTAT |
| Cas9+ODN 10pmol-TRAC-1 | 14-7 | ACCTTGAA | AGTTGCTG | CAGCAACT |
| Cas9+ODN 10pmol-TRAC-2 | 14-8 | ACCTTGAA | TCGCCAGA | TCTGGCGA |

**Note:** The table shows all 164 samples with their corresponding barcodes (I5-I7), I5 sequences, and I7 sequences (both non-reverse and reverse complement).

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
