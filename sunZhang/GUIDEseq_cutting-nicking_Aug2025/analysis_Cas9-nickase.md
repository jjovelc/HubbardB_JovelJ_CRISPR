# CARBON-G1: Merged-read extraction and nicking tolerance quantification

This guide documents the Group 1 workflow for analyzing Cas9 nickase data:

1. Merge paired reads per library
2. Extract the constant-bounded insert with cutadapt
3. Quantify on-target/off-target nicking tolerance from the 14-nt randomized window using `quantify_libs_eif3d.py`

## Prerequisites

- **BBTools** (`bbmerge.sh`) or vsearch/fastp for merging
- **cutadapt** ≥ 3
- **Python** 3.8+
- **Script**: `quantify_libs_eif3d.py` (with optional `--require-pam`)

## Constants (R1/forward orientation)

- **Forward constant (F)**: `ATGGTGAGCAAGGGCGAGG`
- **Reverse primer (R)**: `GCATGGACGAGCTGTACAAGTAA`

**Reverse complements** (for reference):
- RC(F): `CCTCGCCCTTGCTCACCAT`
- RC(R): `TTACTTGTACAGCTCGTCCATGC`

**EIF3D on-target protospacer** (forward, 5′→3′):
`AGACGACCCTGTCATCCGCA` (PAM on genome: NGG)

## 1. Merge R1/R2 per library

Merging produces a single sequence per fragment in R1 orientation.

```bash
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

THREADS=${THREADS:-8}

for DIR in */; do
  base=${DIR%/}
  [[ -e "$DIR"*_1.fq.gz && -e "$DIR"*_2.fq.gz ]] || continue

  # Concatenate lanes (if present) into temporary files
  cat "$DIR"*"_1.fq.gz" > "$DIR/__tmp_R1.fastq.gz"
  cat "$DIR"*"_2.fq.gz" > "$DIR/__tmp_R2.fastq.gz"

  # Merge
  bbmerge.sh in1="$DIR/__tmp_R1.fastq.gz" in2="$DIR/__tmp_R2.fastq.gz" \
             out="${DIR}${base}_merged.fastq.gz" \
             outu1="${DIR}${base}_unmerged_R1.fastq.gz" \
             outu2="${DIR}${base}_unmerged_R2.fastq.gz" \
             qtrim=r trimq=20 threads="$THREADS" ziplevel=6 overwrite=t

  rm -f "$DIR/__tmp_R1.fastq.gz" "$DIR/__tmp_R2.fastq.gz"
done
```

### Sanity checks

```bash
# Read counts (FASTQ lines/4)
zcat <LIB>/*_1.fq.gz | wc -l
zcat <LIB>/*_2.fq.gz | wc -l
zcat <LIB>_merged.fastq.gz | wc -l

# Spot the constants
zcat <LIB>_merged.fastq.gz | grep -c 'ATGGTGAGCAAGGGCGAGG'    # F
zcat <LIB>_merged.fastq.gz | grep -c 'GCATGGACGAGCTGTACAAGTAA' # R (some reads will include this string too)
```

> **Note**: If merge rate is unexpectedly low, skip merging and use paired-end linked extraction (not covered here since your merge rate is high).

## 2. Extract the insert (between constants) with cutadapt

We keep only the sequence between the forward constant and the reverse primer sequence as it appears in the merged read.

```bash
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

for FILE in *[1-2]/*_merged.fastq.gz; do
  cutadapt -j 8 --nextseq-trim=20 -e 0.10 -O 16 --discard-untrimmed \
    -a 'ATGGTGAGCAAGGGCGAGG...GCATGGACGAGCTGTACAAGTAA' \
    -o "${FILE/_merged.fastq.gz/_merged.extracted.fastq.gz}" \
    "$FILE"
done
```

### Parameters explained
- `...` = linked adapters → require both anchors in order; output is the insert between them (constants removed)
- `--discard-untrimmed` = drop reads that don't contain both anchors
- Tune `-e` (max error rate) and `-O` (min overlap) if needed

### Sanity check (lengths)

```bash
zcat <LIB>_merged.extracted.fastq.gz | awk 'NR%4==2{print length($0)}' | head
```

Expect consistent insert lengths per assay design.

## 3. Quantify nicking tolerance from the 14-nt randomized window

Use `quantify_libs_eif3d.py` to:

- Find the Lib1..Lib7 frame in each read (by fixed flanks)
- Extract the 14-nt randomized block (N14)
- Compare that N14 to the correct 14-nt slice of the reverse-complement of the EIF3D 20-mer, based on the Lib frame
- Optionally enforce PAM: require CCN immediately upstream (5′) of the rc(20-mer) in read orientation (equivalent to forward-strand NGG after the 20-mer)
- Bin reads by mismatch count MM0..MM14 (where MM0 = on-target; MM≥1 = off-target)

### Run on one library

```bash
python quantify_libs_eif3d.py \
  --fastq nlib_DNA_10_0_5pmol_1_merged.extracted.fastq.gz \
  --out-prefix work/nlib_DNA_10_0_5pmol_1/nick \
  --ref-protospacer AGACGACCCTGTCATCCGCA \
  --require-pam
```

### Outputs

- `nick_lib_mm_counts.tsv` — Lib1..Lib7 × MM0..MM14 counts (per-lib histogram)
- `nick_lib_totals.tsv` — per-lib totals + on-target (MM0) and off-target (≥MM1) + pam_filtered
- `nick_diag.txt` — diagnostics (ref used, counts of no-match/multi-match/pam-filtered)
- `nick_lib_unique14.tsv` *(optional)* — top N unique 14-mers with counts and MM (enable with `--dump-unique [--topk N]`)

### Batch all libraries

```bash
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

REF="AGACGACCCTGTCATCCGCA"
OUTROOT="work"

for f in *_[1-2]/*_merged.extracted.fastq.gz; do
  sample="$(basename "$(dirname "$f")")"
  outdir="${OUTROOT}/${sample}"
  mkdir -p "$outdir"

  echo "[INFO] $sample → ${outdir}/nick_*"
  python quantify_libs_eif3d.py \
    --fastq "$f" \
    --out-prefix "${outdir}/nick" \
    --ref-protospacer "$REF" \
    --require-pam
done
```
