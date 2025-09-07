#!/bin/bash
#SBATCH -J g1_guide_odn
#SBATCH -o slurm_logs/g1_%A_%a.out
#SBATCH -e slurm_logs/g1_%A_%a.err
#SBATCH -p compute                 # <-- EDIT
#SBATCH -A YOUR_ACCOUNT            # <-- EDIT
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G

set -euo pipefail

############################################
# USER SETTINGS — EDIT THESE
############################################
# Where your demuxed FASTQs live (one dir per sample, as you showed)
FASTQ_DIR="/path/to/demux_fastqs"          # <-- EDIT

# Pattern to find R1/R2 inside each sample directory
R1_GLOB="*_R1*.fastq.gz"
R2_GLOB="*_R2*.fastq.gz"

# The joined CSV we already made (has sample_dir/target/etc.)
JOINED_CSV="joined_samples_with_barcodes.csv"  # <-- put it next to this sbatch file or give full path

# Constants (R1 orientation)
LEFT_CONSTANT="ATGGTGAGCAAGGGCGAGG"        # <-- EDIT (5' constant)
RIGHT_CONSTANT="GCATGGACGAGCTGTACAAGTAA"   # <-- EDIT (3' constant)
ERROR_RATE="0.10"                          # cutadapt -e
MIN_OVERLAP="16"                           # cutadapt -O

# ODN motif (IUPAC); lenient mismatch allowance at fixed positions
ODN_IUPAC="NNCCNNNNNNNNNNNNNNNTCGTCTNN"    # <-- EDIT if needed
LENIENT_MM=1

# Target amplicons/guides (20nt sgRNA; no PAM); edit all three
AMP_CEP290="FILL_ME_AMP_CEP290"
GUIDE_CEP290="FILL_ME_SGRNA_CEP290"
AMP_AAVS1="FILL_ME_AMP_AAVS1"
GUIDE_AAVS1="FILL_ME_SGRNA_AAVS1"
AMP_TRAC="FILL_ME_AMP_TRAC"
GUIDE_TRAC="FILL_ME_SGRNA_TRAC"

# CRISPResso indel window half-size (bp) around cut
CRISP_WIN=5

# Tools (assumes conda env provides these)
CUTADAPT="cutadapt"
FASTQC="fastqc"
CRISPRESSO="CRISPResso"    # CRISPResso2 entrypoint is "CRISPResso"
PYTHON="python"

############################################
# ENV
############################################
mkdir -p slurm_logs
module purge || true
# If you use conda:
if command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  # conda activate guidepipe   # uncomment if you made an env
fi

############################################
# BUILD MANIFEST (once) FROM JOINED CSV
# Produces: manifests/group1_manifest.tsv with columns:
# sample_dir  official_sample_name  category  target  dose_pmol  replicate  R1  R2
############################################
mkdir -p manifests
MANIFEST="manifests/group1_manifest.tsv"

if [ ! -s "$MANIFEST" ]; then
  if [ ! -s "$JOINED_CSV" ]; then
    echo "[ERROR] Missing $MANIFEST and $JOINED_CSV. Provide one of them." >&2
    exit 1
  fi
  echo "[INFO] Creating $MANIFEST from $JOINED_CSV and FASTQ_DIR=$FASTQ_DIR"
  $PYTHON - "$JOINED_CSV" "$FASTQ_DIR" "$R1_GLOB" "$R2_GLOB" "$MANIFEST" <<'PYCODE'
import sys, os, glob, pandas as pd
join_csv, fastq_dir, r1_glob, r2_glob, out_tsv = sys.argv[1:]
df = pd.read_csv(join_csv)
df = df[df["group"]=="cutting"].copy()
def find(globpat):
    hits = sorted(glob.glob(globpat))
    return hits[0] if hits else ""
r1s, r2s = [], []
for s in df["sample_dir"]:
    r1 = find(os.path.join(fastq_dir, s, r1_glob))
    r2 = find(os.path.join(fastq_dir, s, r2_glob))
    r1s.append(r1); r2s.append(r2)
df["R1"] = r1s; df["R2"] = r2s
cols = ["sample_dir","official_sample_name","category","target","dose_pmol","replicate","R1","R2"]
df[cols].to_csv(out_tsv, sep="\t", index=False)
print(f"[OK] wrote {out_tsv} with {len(df)} rows")
PYCODE
fi

############################################
# ARRAY SIZING / DRY RUN
############################################
# N = data rows (skip header)
N=$(($(wc -l < "$MANIFEST") - 1))
if [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 0 ]; then
  echo "[INFO] Submit this as:  sbatch --array=1-${N} $0"
  exit 0
fi

############################################
# PICK ROW
############################################
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$MANIFEST")
IFS=$'\t' read -r SAMPLE OFFNAME CATEGORY TARGET DOSE REPL R1 R2 <<< "$LINE"

if [ ! -s "$R1" ]; then
  echo "[WARN] No R1 found for $SAMPLE — skipping"
  exit 0
fi

WDIR="work/group1/${SAMPLE}"
mkdir -p "$WDIR/fastqc"

echo "[INFO] Sample=$SAMPLE  Target=$TARGET"
echo "[INFO] R1=$R1"
[ -n "${R2:-}" ] && echo "[INFO] R2=$R2"

############################################
# 1) FastQC
############################################
$FASTQC -t ${SLURM_CPUS_PER_TASK} -o "$WDIR/fastqc" "$R1" ${R2:+$R2}

############################################
# 2) Extract insert between constants with cutadapt linked adapters
#    Keeps only the sequence BETWEEN LEFT ... RIGHT
############################################
$CUTADAPT -e ${ERROR_RATE} -O ${MIN_OVERLAP} \
  -g "${LEFT_CONSTANT}...${RIGHT_CONSTANT}" \
  --discard-untrimmed \
  -j ${SLURM_CPUS_PER_TASK} \
  -o "$WDIR/extracted_R1.fastq.gz" "$R1" > "$WDIR/cutadapt_R1.txt"

############################################
# 3) ODN motif scan (strict + lenient)
############################################
$PYTHON - "$WDIR/extracted_R1.fastq.gz" "$ODN_IUPAC" "$LENIENT_MM" "$WDIR/odn_counts.tsv" "$WDIR/odn_counts.json" <<'PYCODE'
import sys, gzip, re, json
fastq, motif, lenient_k, out_tsv, out_json = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
IUPAC = {"A":"A","C":"C","G":"G","T":"T","R":"[AG]","Y":"[CT]","S":"[GC]","W":"[AT]","K":"[GT]","M":"[AC]","B":"[CGT]","D":"[AGT]","H":"[ACT]","V":"[ACG]","N":"[ACGT]"}
rx = re.compile("".join(IUPAC.get(b,b) for b in motif.upper()))
def fq(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            h=fh.readline().rstrip()
            if not h: break
            s=fh.readline().rstrip()
            fh.readline(); fh.readline()
            yield s
strict=lenient=total=0
L = len(motif)
fixed_positions = [i for i,ch in enumerate(motif.upper()) if ch not in ("N")]
for seq in fq(fastq):
    total += 1
    if rx.search(seq):
        strict += 1
        continue
    # lenient: allow <=k mismatches at fixed (non-N) positions in any window of length L
    best=None
    for i in range(0, max(0, len(seq)-L+1)):
        window = seq[i:i+L].upper()
        if len(window)<L: continue
        mm=0
        for j in fixed_positions:
            if window[j] != motif[j]:
                mm += 1
                if mm > lenient_k: break
        if mm <= lenient_k:
            lenient += 1
            break
with open(out_json,"w") as jf:
    json.dump({"total":total,"strict":strict,"lenient":lenient,
               "strict_rate": (strict/total if total else 0.0),
               "lenient_rate": (lenient/total if total else 0.0)}, jf, indent=2)
with open(out_tsv,"w") as tf:
    tf.write("total\tstrict\tlenient\tstrict_rate\tlenient_rate\n")
    tf.write(f"{total}\t{strict}\t{lenient}\t{(strict/total if total else 0.0):.6f}\t{(lenient/total if total else 0.0):.6f}\n")
PYCODE

############################################
# 4) CRISPResso2 (amplicon mode) for indel %
############################################
AMP=""
GUIDE=""
case "$TARGET" in
  CEP290) AMP="$AMP_CEP290"; GUIDE="$GUIDE_CEP290" ;;
  AAVS1)  AMP="$AMP_AAVS1";  GUIDE="$GUIDE_AAVS1"  ;;
  TRAC)   AMP="$AMP_TRAC";   GUIDE="$GUIDE_TRAC"   ;;
  *) echo "[WARN] Unknown target '$TARGET' — skipping CRISPResso";;
esac

if [ -n "$AMP" ] && [ -n "$GUIDE" ] && [[ "$AMP" != FILL_ME_* ]] && [[ "$GUIDE" != FILL_ME_* ]]; then
  $CRISPRESSO \
    --fastq_r1 "$WDIR/extracted_R1.fastq.gz" \
    --amplicon_seq "$AMP" \
    --guide_seq "$GUIDE" \
    --quantification_window_size ${CRISP_WIN} \
    --name "${SAMPLE}" \
    --output_folder "$WDIR" \
    --n_processes ${SLURM_CPUS_PER_TASK} || echo "[WARN] CRISPResso failed for $SAMPLE"
else
  echo "[WARN] Missing amplicon/guide for $TARGET; skipped CRISPResso."
fi

echo "[OK] Done: $SAMPLE"

