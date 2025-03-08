# Sun Zhang's page

Sun Zhang has been working of GUIDEseq.

## Description of experiment

### Experiment: GUIDEseq January 2025.

Cells are incubated with Cas9 RNP system, along with a short, double-stranded oligodeoxynucleotide (dsODN) tag. When Cas9 cuts the DNA, creates double-strand breaks (DSBs) at both intended and unintended sites, the cell's natural repair process incorporates the dsODN into these break points. After extracting genomic DNA, PCR is used to recognize the dsODN sequence, enriching the DNA fragments that have incorporated the tag. These fragments are then analyzed by NGS. By mapping where the dsODN has been inserted, we can generate genome-wide off-target effects caused by Cas9. 

The following table explains the barcoding scheme. Sense and anti-sense are not barcoded separately, they have the same barcodes.

| Sample     | i5        | i7        |
|------------|----------|----------|
| CEP290-1 + | TAGATCGC | TAAGGCGA |
| CEP290-2 + | TAGATCGC | CGTACTAG |
| CEP290-3 + | TAGATCGC | AGGCAGAA |
| PDCD1-1 +  | CTCTCTAT | TAAGGCGA |
| PDCD1-2 +  | CTCTCTAT | CGTACTAG |
| PDCD1-3 +  | CTCTCTAT | AGGCAGAA |
| HBB-1 +    | TATCCTCT | TAAGGCGA |
| HBB-2 +    | TATCCTCT | CGTACTAG |
| HBB-3 +    | TATCCTCT | AGGCAGAA |
| ODN-1 +    | AGAGTAGA | TAAGGCGA |
| ODN-2 +    | AGAGTAGA | CGTACTAG |
| ODN-3 +    | AGAGTAGA | AGGCAGAA |
| CEP290-1 - | TAGATCGC | TAAGGCGA |
| CEP290-2 - | TAGATCGC | CGTACTAG |
| CEP290-3 - | TAGATCGC | AGGCAGAA |
| PDCD1-1 -  | CTCTCTAT | TAAGGCGA |
| PDCD1-2 -  | CTCTCTAT | CGTACTAG |
| PDCD1-3 -  | CTCTCTAT | AGGCAGAA |
| HBB-1 -    | TATCCTCT | TAAGGCGA |
| HBB-2 -    | TATCCTCT | CGTACTAG |
| HBB-3 -    | TATCCTCT | AGGCAGAA |
| ODN-1 -    | AGAGTAGA | TAAGGCGA |
| ODN-2 -    | AGAGTAGA | CGTACTAG |
| ODN-3 -    | AGAGTAGA | AGGCAGAA |


The GUIDE-Seq pipeline was cloned from the corresponding Github repo [Github GUIDEseq repo](https://github.com/aryeelab/guideseq).
