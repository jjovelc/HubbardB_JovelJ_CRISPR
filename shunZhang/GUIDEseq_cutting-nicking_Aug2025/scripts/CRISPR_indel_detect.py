#!/usr/bin/env python3
"""
CRISPR_indel_detect.py  —  strand-aware indel calling around a Cas9 cut site.

Inputs:
  1) FASTQ(.gz) of merged reads
  2) Genomic regions FASTA (≈200 bp per target around protospacer)
Options:
  --targets FILE        targets.txt (e.g., 'AAVS1: PROTOSPACER+PAM')
  --output FILE         output TSV with per-event counts
  --min-quality N       minimum Phred33 mean quality to keep read (default 20)
  --max-indel N         max absolute indel size to report (default 100)
  --min-flank N         minimal flank size on each side (default 15)
  --max-flank N         maximal flank size to try (default 50)
Notes:
  * Strand is detected automatically by searching both '+' and '-' for the 20 nt
    protospacer (PAM is inferred). SpCas9 cut site is:
      '+' strand: cut = PAM_start - 3
      '-' strand: cut = PAM_end + 3     (PAM on ref is 'CC' for NGG on target)
  * Indel sign convention: insertion = +k, deletion = -k (independent of strand)
"""

import argparse
import gzip
import io
import os
import sys
from collections import Counter, defaultdict
from typing import Dict, Tuple, Iterable

# -----------------------
# Utilities (I/O, FASTA)
# -----------------------

def open_maybe_gzip(path: str) -> io.TextIOBase:
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path, "r")

def read_fasta(path: str) -> Dict[str, str]:
    seqs = {}
    name = None
    buf = []
    with open(path, "r") as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).replace("-", "").replace(" ", "").upper()
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            seqs[name] = "".join(buf).replace("-", "").replace(" ", "").upper()
    return seqs

def read_targets(path: str) -> Dict[str, str]:
    """
    Parse a 'targets.txt' file like:
      AAVS1: GTCCCTAGTGGCCCCACTGTNGG
      CEP290: GGAGTCACATGGGAGTCACANGG
      TRAC:   GCTGGTACACGGCAGGGTCANGG
    Returns dict target -> guide_with_pam (string).
    """
    out = {}
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" in line:
                k, v = line.split(":", 1)
                k = k.strip()
                v = v.strip().split()[0]
                if k and v:
                    out[k] = v
    return out

# -----------------------
# Strand-aware helpers
# -----------------------

def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn",
                          "TGCAYRKMVHDBNtgcayrkmvhdbn")
    return seq.translate(table)[::-1]

def strip_pam(guide_with_pam: str) -> str:
    g = guide_with_pam.upper().replace(" ", "")
    # keep 20 nts for protospacer if something like NGG is appended
    return g[:20] if len(g) >= 21 else g

def find_guide_and_cut(region: str,
                       guide_with_pam: str,
                       pam="GG") -> Tuple[int, str, int, Tuple[int,int]]:
    """
    Locate protospacer on either strand and compute Cas9 cut position.
    Returns:
        cut_pos  (int): 0-based on forward reference coordinates
        strand   (str): '+' or '-'
        gstart   (int): 0-based start of protospacer on forward ref
        pam_span (tuple): (pam_start, pam_end) on forward ref
    """
    s = region.upper()
    prot = strip_pam(guide_with_pam)

    # '+' strand: protospacer appears as-is; PAM (GG) should be 3' (downstream)
    g0 = s.find(prot)
    if g0 != -1:
        pam_start = g0 + len(prot)
        pam_end = pam_start + 2
        # allow small displacement if reference was trimmed oddly
        if not (pam_end <= len(s) and s[pam_start:pam_end] == pam):
            found = False
            for off in range(1, 6):
                if pam_start + off + 2 <= len(s) and s[pam_start+off:pam_start+off+2] == pam:
                    pam_start, pam_end = pam_start + off, pam_start + off + 2
                    found = True
                    break
            if not found:
                pam_start, pam_end = g0 + len(prot), g0 + len(prot) + 2
        cut = pam_start - 3
        return cut, "+", g0, (pam_start, pam_end)

    # '-' strand: protospacer appears as reverse-complement; PAM on forward is 'CC'
    rc_prot = reverse_complement(prot)
    g0 = s.find(rc_prot)
    if g0 != -1:
        pam_start = max(0, g0 - 2)
        pam_end = g0
        if not (pam_end <= len(s) and s[pam_start:pam_end] == "CC"):
            found = False
            for off in range(1, 6):
                left = g0 - off - 2
                if left >= 0 and s[left:left+2] == "CC":
                    pam_start, pam_end = left, left + 2
                    found = True
                    break
            if not found:
                pam_start, pam_end = max(0, g0 - 2), g0
        cut = pam_end + 3
        return cut, "-", g0, (pam_start, pam_end)

    # Fallback: midpoint (keeps old behavior, but you’ll see strand='?')
    return len(s) // 2, "?", -1, (len(s)//2, len(s)//2)

# -----------------------
# FASTQ parsing / QC
# -----------------------

def iter_fastq(path: str) -> Iterable[Tuple[str, str, str]]:
    """
    Yield (name, seq, qual). Supports gzipped input.
    """
    with open_maybe_gzip(path) as fh:
        while True:
            name = fh.readline()
            if not name:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                break
            if not name.startswith("@"):
                continue
            yield name.strip(), seq.strip().upper(), qual.strip()

def mean_phred33(q: str) -> float:
    if not q:
        return 0.0
    return sum(ord(c) - 33 for c in q) / len(q)

# -----------------------
# Indel detection (flank matching)
# -----------------------

def find_indel_by_flanks(read: str,
                         ref: str,
                         cut_pos: int,
                         min_flank: int,
                         max_flank: int,
                         max_indel: int) -> Tuple[bool, int]:
    """
    Detect indel size at cut site by searching flanks around cut in the read.

    Strategy:
      1) Choose a flank size w from max_flank down to min_flank.
      2) Extract left = ref[cut_pos-w:cut_pos], right = ref[cut_pos:cut_pos+w].
      3) Find last occurrence of 'left' in read and first occurrence of 'right' in read.
      4) If both found and left_idx < right_idx, compute:
            ins = (right_idx) - (left_idx + len(left))   # extra bases between flanks in read
            del_len = max(0, (2*w) - (len(left) + ins + len(right)))  # collapse across cut in read
         Because left|right are adjacent in the reference, the "gap" in the read
         is insertion if positive, and if flanks overlap/touch, any missing
         reference bases across the cut is treated as deletion.
    """
    L = len(ref)
    for w in range(max_flank, min_flank - 1, -1):
        l0 = max(0, cut_pos - w)
        r1 = min(L, cut_pos + w)
        if r1 - l0 < 2*w:  # near edges
            continue
        left = ref[cut_pos - w:cut_pos]
        right = ref[cut_pos:cut_pos + w]
        if len(left) != w or len(right) != w:
            continue

        li = read.rfind(left)
        if li == -1:
            continue
        ri = read.find(right, li + len(left) - w//2)  # allow small overlap search
        if ri == -1:
            continue
        if li >= ri:
            continue

        # bases between the two flanks in the read
        ins = max(0, ri - (li + len(left)))
        # If the read glues left and right with little/no gap, treat deficit as deletion
        total_span_in_read = (ri + len(right)) - li
        expected_span = len(left) + len(right)  # flanks only
        del_len = max(0, expected_span - total_span_in_read)
        indel = ins - del_len  # insertion positive; deletion negative

        if abs(indel) <= max_indel:
            return True, indel

    return False, 0

# -----------------------
# Main
# -----------------------

def main():
    ap = argparse.ArgumentParser(description="Strand-aware CRISPR indel detector")
    ap.add_argument("fastq", help="merged FASTQ(.gz)")
    ap.add_argument("genomic_fasta", help="FASTA with ~200bp regions per target")
    ap.add_argument("--targets", default="targets.txt", help="targets file (default: targets.txt)")
    ap.add_argument("--output", required=True, help="output TSV path")
    ap.add_argument("--min-quality", type=int, default=20)
    ap.add_argument("--max-indel", type=int, default=100)
    ap.add_argument("--min-flank", type=int, default=15)
    ap.add_argument("--max-flank", type=int, default=50)

    ap.add_argument("--min-score", type=int, default=25, help=argparse.SUPPRESS)
    ap.add_argument("--prefilter-k", type=int, default=15, help=argparse.SUPPRESS)
    ap.add_argument("--anchor-offsets", type=str, default="0,12,24", help=argparse.SUPPRESS)
    ap.add_argument("--odn", type=str, default=None, help=argparse.SUPPRESS)
    ap.add_argument("--odn-k", type=int, default=20, help=argparse.SUPPRESS)
    args = ap.parse_args()

    if not os.path.exists(args.fastq):
        sys.exit(f"[error] FASTQ not found: {args.fastq}")
    if not os.path.exists(args.genomic_fasta):
        sys.exit(f"[error] FASTA not found: {args.genomic_fasta}")
    if not os.path.exists(args.targets):
        sys.exit(f"[error] targets file not found: {args.targets}")

    targets = read_targets(args.targets)           # e.g., {'AAVS1': '...NGG', ...}
    regions = read_fasta(args.genomic_fasta)      # headers should contain AAVS1, CEP290, TRAC keys

    # map region name to the best target key by prefix match
    name_map = {}
    for rname in regions:
        key = None
        for t in targets:
            if rname.upper().startswith(t.upper()):
                key = t
                break
        name_map[rname] = key

    # Pre-compute cut sites per region/target
    region_info = {}
    for rname, seq in regions.items():
        tkey = name_map[rname]
        if tkey is None:
            print(f"[warn] No matching target for region '{rname}'. Skipping.", file=sys.stderr)
            continue
        guide = targets[tkey]
        cut, strand, gstart, pam_span = find_guide_and_cut(seq, guide, pam="GG")
        region_info[rname] = {
            "target": tkey,
            "guide": guide,
            "strand": strand,
            "cut": cut,
            "pam": pam_span,
            "seq": seq,
        }
        print(f"[info] {rname}: target={tkey} strand={strand} cut={cut} guide_start={gstart}", file=sys.stderr)

    if not region_info:
        sys.exit("[error] No regions with detected guides were found.")

    # Counters
    total_reads = 0
    kept_reads = 0
    per_target_counts = defaultdict(Counter)   # target -> Counter(indel_size)
    per_target_analyzed = Counter()           # target -> reads with both flanks found

    # Process reads
    for name, seq, qual in iter_fastq(args.fastq):
        total_reads += 1
        if args.min_quality > 0 and mean_phred33(qual) < args.min_quality:
            continue
        kept_reads += 1

        # Try each region; stop at first success to avoid double counting
        assigned = False
        for rname, info in region_info.items():
            ref = info["seq"]
            cut = info["cut"]
            ok, indel = find_indel_by_flanks(seq, ref, cut,
                                             args.min_flank, args.max_flank, args.max_indel)
            if ok:
                per_target_analyzed[info["target"]] += 1
                per_target_counts[info["target"]][indel] += 1
                assigned = True
                break

        # If not assigned, we ignore the read (either off-target or low information)

    # Write output
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w") as out:
        print("# CRISPR indel results", file=out)
        print(f"# fastq\t{args.fastq}", file=out)
        print(f"# genomic_fasta\t{args.genomic_fasta}", file=out)
        print(f"# total_reads\t{total_reads}", file=out)
        print(f"# reads_passing_qc\t{kept_reads}", file=out)
        print("# target\tguide\tstrand\tcut_pos\tanalyzed_reads\tindel_size\tcount", file=out)

        for rname, info in region_info.items():
            t = info["target"]
            guide = info["guide"]
            strand = info["strand"]
            cut = info["cut"]
            analyzed = per_target_analyzed[t]
            counts = per_target_counts[t]
            if counts:
                for k, v in sorted(counts.items(), key=lambda kv: kv[0]):
                    out.write(f"{t}\t{guide}\t{strand}\t{cut}\t{analyzed}\t{k}\t{v}\n")
            else:
                out.write(f"{t}\t{guide}\t{strand}\t{cut}\t{analyzed}\tNA\t0\n")

    print(f"[done] Wrote: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()

