# Detecting dsODN Insertions from CRISPR Cuts (Merged PE Libraries)

### NOTE
Also see the CRISPResso2 report [here](CRISPR_indel_detection.md) where ODN are also reported as indels at position +34.

This document explains the end-to-end procedure to detect dsODN (GUIDE-style) insertion junctions in standard Illumina libraries without UMIs, using:

- A small Python script that builds a compact panel of expected junction sequences for each guide and ODN orientation
- A Bash script that maps your merged (or unmerged fallback) reads to that panel and tallies junction hits

It's designed for your directory layout (duplicate libraries ending in `_1` and `_2`) and your `targets.txt` format.

## Biological rationale (why 17 bp?)

SpCas9 cuts 3 bp upstream of the PAM on the target strand. If the target is written 5′→3′ as 20nt + NGG, the cleavage occurs between bases 17|18 of the protospacer:

```
5' [1..........................20] NGG
             ^ cut between 17|18
LEFT17 = bases 1..17
```

If a blunt dsODN integrates at the cut, junction reads will contain:

- A **genomic anchor** = the LEFT17 bases immediately upstream of the cut
- The **ODN** (either orientation)

We therefore search for **four junctions per guide**:
- `LEFT17+ODN`
- `LEFT17+revcomp(ODN)` 
- `revcomp(LEFT17)+ODN`
- `revcomp(LEFT17)+revcomp(ODN)`

## Inputs & expected layout

**One directory per library** (e.g., `Cas9_ODN_10pmol_AAVS1_1`, `..._2`, etc.).

Inside each directory, either:
- A merged reads file (`*_merged*.fastq.gz` or `*.fq.gz`), or
- Unmerged pair (`*_unmerged*_R1*.fastq.gz` and `*_unmerged*_R2*.fastq.gz`)

**A `targets.txt` like:**
```
AAVS1: GTCCCTAGTGGCCCCACTGTNGG
CEP290: GGAGTCACATGGGAGTCACANGG
TRAC: GCTGGTACACGGCAGGGTCANGG

ODN sequence: GTTTAATTGAGTTGTCATATGTTAATAACGGTAT
```

### Notes
- Each line before the blank line is `GENE: PROTOSPACER+PAM` (we use the first 20 nt as the protospacer; PAM can include N)
- The ODN is the double-stranded tag sequence used in integration

## Software requirements

- **Python** 3.7+
- **bwa** (e.g., 0.7.17+), **samtools** (1.10+), **awk**
- **Optional**: seqkit if you want exact deduplication before mapping

## Step 1 — Build the junction panel (Python)

Create `make_junctions.py`:

```python
#!/usr/bin/env python3
import sys, re

def revcomp(s):
    comp = str.maketrans('ACGTNacgtn','TGCANtgcan')
    return s.translate(comp)[::-1]

if len(sys.argv) != 3:
    sys.stderr.write("Usage: make_junctions.py targets.txt junctions.fa\n")
    sys.exit(1)

targets_txt = sys.argv[1]
out_fa      = sys.argv[2]

protospacers = {}
odn = None

with open(targets_txt) as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        m = re.match(r'^([A-Za-z0-9_]+):\s+([ACGTNacgtn]+)$', line)
        if m:
            name, seq = m.group(1), m.group(2).upper()
            protospacer = seq[:20]      # 20 nt before PAM
            left17 = protospacer[:17]   # anchor up to cut
            protospacers[name] = left17
            continue
        m2 = re.match(r'^ODN sequence:\s+([ACGTNacgtn]+)', line)
        if m2:
            odn = m2.group(1).upper()

if not odn:
    raise SystemExit("ERROR: ODN sequence not found in targets file.")

with open(out_fa, 'w') as out:
    for gene, left17 in protospacers.items():
        left17_rc = revcomp(left17)
        odn_rc    = revcomp(odn)
        combos = [
            (f"{gene}__L17__ODN",          left17 + odn),
            (f"{gene}__L17__ODNrc",        left17 + odn_rc),
            (f"{gene}__L17rc__ODN",        left17_rc + odn),
            (f"{gene}__L17rc__ODNrc",      left17_rc + odn_rc),
        ]
        for name, seq in combos:
            out.write(f">{name}\n{seq}\n")
```

**Run:**
```bash
python3 make_junctions.py targets.txt junctions.fa
```

**Output**: `junctions.fa` with 4 sequences × (# guides). This is the only reference we align to.

## Step 2 — Map reads and count junction hits (Bash)

Create `map_merged_reads.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob extglob

REF="junctions.fa"          # built in Step 1
OUTDIR="map_out"
THREADS="${THREADS:-8}"

mkdir -p "$OUTDIR"

# Index once
if [ ! -e "${REF}.bwt" ]; then
  echo "[INFO] Indexing $REF"
  bwa index "$REF"
fi

# Header
echo -e "sample\tjunction\tcount" > "$OUTDIR/junction_counts.tsv"

for DIR in *[12]; do
  [ -d "$DIR" ] || continue
  sample="$DIR"

  # NOTE: underscore before merged so it won't match 'unmerged'
  merged=( "$DIR"/*_merged*.fastq.gz "$DIR"/*_merged*.fq.gz )
  r1=( "$DIR"/*_unmerged*_R1*.fastq.gz "$DIR"/*_unmerged*_R1*.fq.gz )
  r2=( "$DIR"/*_unmerged*_R2*.fastq.gz "$DIR"/*_unmerged*_R2*.fq.gz )

  if (( ${#merged[@]} >= 1 )); then
    in1="${merged[0]}"
    echo "[INFO] ${sample}: using merged: $in1"
    bwa mem -t "$THREADS" -k 15 -T 20 -a "$REF" "$in1" \
    | samtools view -@ "$THREADS" -F 0x904 -q 20 - \
    | awk -v s="$sample" 'BEGIN{OFS="\t"} $3!="*" {c[$3]++} END{for (j in c) print s,j,c[j]}' \
    >> "$OUTDIR/junction_counts.tsv"

  elif (( ${#r1[@]} == 1 && ${#r2[@]} == 1 )); then
    echo "[INFO] ${sample}: using unmerged pair: ${r1[0]} ${r2[0]}"
    bwa mem -t "$THREADS" -k 15 -T 20 -a "$REF" "${r1[0]}" "${r2[0]}" \
    | samtools view -@ "$THREADS" -F 0x904 -q 20 - \
    | awk -v s="$sample" 'BEGIN{OFS="\t"} $3!="*" {c[$3]++} END{for (j in c) print s,j,c[j]}' \
    >> "$OUTDIR/junction_counts.tsv"

  else
    echo "[WARN] ${sample}: no reads found (merged or R1/R2); skipping" >&2
  fi
done

echo "[DONE] Wrote $OUTDIR/junction_counts.tsv"
```

**Run:**
```bash
bash map_merged_reads.sh
```

**Output**: `map_out/junction_counts.tsv` with columns:
- **sample** (directory name)
- **junction** (e.g., `AAVS1__L17__ODNrc`)
- **count** (primary alignments with MAPQ≥20)

## Why `-k 15 -T 20 -F 0x904 -q 20`?

- **`-k 15`**: shorter seeds improve sensitivity for short junctions
- **`-T 20`**: drop weak alignments  
- **`-F 0x904`**: remove secondary (0x100), supplementary (0x800), and QC-fail (0x200)
- **`-q 20`**: MAPQ filter to curb spurious hits
