#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob   # patterns with no matches expand to empty

THREADS=${THREADS:-8}

for DIR in */; do
  base=${DIR%/}
  # process only sample dirs that end in _1 or _2 (your naming scheme)
  [[ $base =~ _[12]$ ]] || continue

  # collect R1/R2 files (handles 1 or many lanes)
  r1=( "$DIR"*"_1.fq.gz" )
  r2=( "$DIR"*"_2.fq.gz" )
  if (( ${#r1[@]} == 0 || ${#r2[@]} == 0 )); then
    echo "[WARN] $base: no *_1.fq.gz or *_2.fq.gz — skipping"
    continue
  fi

  # if multiple lanes, concatenate to temp files
  tmpR1="$DIR/__tmp_R1.fastq.gz"
  tmpR2="$DIR/__tmp_R2.fastq.gz"
  cat "${r1[@]}" > "$tmpR1"
  cat "${r2[@]}" > "$tmpR2"

  out_merged="$DIR/${base}_merged.fastq.gz"
  outu1="$DIR/${base}_unmerged_R1.fastq.gz"
  outu2="$DIR/${base}_unmerged_R2.fastq.gz"

  echo "[INFO] Merging $base → $out_merged"
  bbmerge.sh in1="$tmpR1" in2="$tmpR2" \
             out="$out_merged" \
             outu1="$outu1" outu2="$outu2" \
             qtrim=r trimq=20 threads="$THREADS" ziplevel=6 overwrite=t

  rm -f "$tmpR1" "$tmpR2"
done
