#!/usr/bin/env python3
"""
CRISPR_indel_detect.py
Detect insertions and deletions at CRISPR cut sites from merged FASTQ reads.

Key fixes vs. your previous version:
- Proper bidirectional guide matching with correct cut-site logic
- True deletion detection (right flank can start before left flank end)
- Actual use of --min-quality on FASTQ qualities (PHRED+33)
- Clear denominators and true editing efficiency
- Optional ODN-aware insertion tagging
"""

import os
import sys
import re
import gzip
import argparse
from collections import Counter

# ----------------------------
# Utilities
# ----------------------------

def reverse_complement(seq: str) -> str:
    tbl = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(tbl)[::-1]

def mean_phred33(q: str) -> float:
    # PHRED+33
    return sum(ord(c) - 33 for c in q) / max(1, len(q))

def find_all(haystack: str, needle: str, start: int = 0):
    """Yield all positions of needle in haystack from start (inclusive)."""
    i = haystack.find(needle, start)
    while i != -1:
        yield i
        i = haystack.find(needle, i + 1)

# ----------------------------
# Reference & cut-site logic
# ----------------------------

GUIDES_20NT = {
    # Use the canonical 20-nt protospacers (no PAM)
    # The user-supplied guides had an N before PAM; we trim to 20 nt.
    "AAVS1": "GTCCCTAGTGGCCCCACTGT",
    "CEP290": "GGAGTCACATGGGAGTCACA",
    "TRAC":   "GCTGGTACACGGCAGGGTCA",
}

def parse_fasta(path: str):
    if not os.path.exists(path):
        print(f"ERROR: Reference file not found: {path}")
        sys.exit(1)

    with open(path, "r") as fh:
        name, seq = None, []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq).upper()
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name:
            yield name, "".join(seq).upper()

def find_cut_position(seq_id: str, seq: str):
    """
    Discover the cut site by looking for either the forward guide OR its reverse complement.
    If forward protospacer is found at index 'pos' (0-based), cut = pos + 17.
    If reverse protospacer (rc) is found at index 'pos', cut = pos + 3.
    Fallback: midpoint.
    Returns: (cut_pos, matched_gene, strand)
    """
    # Normalize AASV1 typo if present
    uc_id = seq_id.upper().replace("AASV1", "AAVS1")

    for gene, guide in GUIDES_20NT.items():
        if gene in uc_id:
            # Try forward
            pos_f = seq.find(guide)
            if pos_f != -1:
                return pos_f + 17, gene, '+'
            # Try reverse complement
            rc = reverse_complement(guide)
            pos_r = seq.find(rc)
            if pos_r != -1:
                return pos_r + 3, gene, '-'
            # Not found; still this gene id matched
            return len(seq) // 2, gene, '?'

    # Gene name not in header; try all guides anyway
    for gene, guide in GUIDES_20NT.items():
        pos_f = seq.find(guide)
        if pos_f != -1:
            return pos_f + 17, gene, '+'
        rc = reverse_complement(guide)
        pos_r = seq.find(rc)
        if pos_r != -1:
            return pos_r + 3, gene, '-'

    # Last resort
    return len(seq) // 2, "UNKNOWN", '?'

def load_references(ref_fasta_path: str):
    refs = {}
    print(f"Reading reference file: {ref_fasta_path}")
    for sid, sseq in parse_fasta(ref_fasta_path):
        cut, gene, strand = find_cut_position(sid, sseq)
        refs[sid] = {
            "sequence": sseq,
            "cut_position": cut,
            "left_flank": sseq[:cut],
            "right_flank": sseq[cut:],
            "gene": gene,
            "strand": strand,
        }
        print(f"  {sid}: {len(sseq)} bp | cut @ {cut} | gene={gene} | strand={strand}")
    if not refs:
        print("ERROR: No sequences parsed from reference FASTA.")
        sys.exit(1)
    return refs

# ----------------------------
# Alignment around cut
# ----------------------------

def find_best_alignment(
    read_seq: str,
    ref_left: str,
    ref_right: str,
    max_indel_size: int = 100,
    min_flank_match: int = 15,
    min_score: int = 30,
):
    """
    Align read against reference split at cut:
      - Match suffix of left_flank and prefix of right_flank
      - Allow right flank to begin within a window [left_end - max_indel, left_end + max_indel]
        so that negative deltas (deletions) are discoverable.
    Returns dict or None.
    """
    read = read_seq.upper()
    L = ref_left.upper()
    R = ref_right.upper()

    best = None
    best_score = -1

    # Cap flank scans at 50 nt, but not beyond bounds
    max_left = min(len(L), len(read) - min_flank_match, 50)
    if max_left < min_flank_match:
        return None

    for left_len in range(min_flank_match, max_left + 1):
        left_ref = L[-left_len:]  # suffix
        left_match_pos = read.find(left_ref)
        if left_match_pos == -1:
            continue
        left_end = left_match_pos + left_len

        # Right flank length loop
        max_right = min(len(R), len(read), 50)
        if max_right < min_flank_match:
            continue

        for right_len in range(min_flank_match, max_right + 1):
            right_ref = R[:right_len]  # prefix

            # Window where right flank can begin (allow overlap for deletions)
            win_start = max(0, left_end - max_indel_size)
            win_end   = min(len(read), left_end + max_indel_size)

            for rs in find_all(read, right_ref, start=win_start):
                if rs > win_end:
                    break

                delta = rs - left_end  # >0 insertion; <0 deletion; 0 perfect
                if abs(delta) > max_indel_size:
                    continue

                score = left_len + right_len
                if score < min_score:
                    continue

                if score > best_score:
                    if delta > 0:
                        seq_ins = read[left_end:rs]
                        call = {
                            "type": "insertion",
                            "size": delta,
                            "sequence": seq_ins,
                            "score": score,
                            "left_len": left_len,
                            "right_len": right_len,
                        }
                    elif delta < 0:
                        call = {
                            "type": "deletion",
                            "size": -delta,
                            "sequence": "",
                            "score": score,
                            "left_len": left_len,
                            "right_len": right_len,
                        }
                    else:
                        call = {
                            "type": "perfect",
                            "size": 0,
                            "sequence": "",
                            "score": score,
                            "left_len": left_len,
                            "right_len": right_len,
                        }
                    best = call
                    best_score = score

    return best

# ----------------------------
# ODN tagging
# ----------------------------

def contains_odn(insert_seq: str, odn_seq: str, k: int = 20) -> bool:
    """Return True if a k-mer from ODN (or its RC) occurs in the inserted sequence."""
    if not odn_seq or not insert_seq:
        return False
    odn = odn_seq.upper()
    rco = reverse_complement(odn)
    ins = insert_seq.upper()

    if len(odn) < k:
        k = len(odn)
    if k <= 0:
        return False

    for i in range(len(odn) - k + 1):
        sub = odn[i:i+k]
        if sub in ins:
            return True
    for i in range(len(rco) - k + 1):
        sub = rco[i:i+k]
        if sub in ins:
            return True
    return False

# ----------------------------
# Main analysis
# ----------------------------

def analyze_reads(
    fastq_file: str,
    references: dict,
    min_quality: int = 20,
    max_indel_size: int = 100,
    min_flank_match: int = 15,
    min_score: int = 30,
    odn_seq: str = "",
    odn_k: int = 20,
):
    if not os.path.exists(fastq_file):
        print(f"ERROR: FASTQ file not found: {fastq_file}")
        return {}

    print(f"FASTQ size: {os.path.getsize(fastq_file):,} bytes")

    # Results, per reference
    results = {
        ref_name: {
            "total_reads": 0,
            "prefilter_hits": 0,   # contains either left OR right 25-mer
            "both_flanks_hits": 0, # an alignment was found (>= min_flank on both sides)
            "classified": 0,       # perfect + insertion + deletion
            "indel_counts": Counter(),
        }
        for ref_name in references
    }

    # FASTQ reader (4-line records)
    if fastq_file.endswith(".gz"):
        fh = gzip.open(fastq_file, "rt")
    else:
        fh = open(fastq_file, "r")

    try:
        while True:
            hdr = fh.readline()
            if not hdr:
                break
            seq = fh.readline().strip()
            plus = fh.readline()
            qual = fh.readline().strip()

            # Defensive: malformed FASTQ
            if not seq or not qual:
                continue

            # Quality gate
            if mean_phred33(qual) < min_quality:
                # Still count as total for every reference
                for ref_name in results:
                    results[ref_name]["total_reads"] += 1
                continue

            # Process for each reference
            for ref_name, ref in references.items():
                results[ref_name]["total_reads"] += 1

                L = ref["left_flank"]
                R = ref["right_flank"]

                left_partial  = L[-25:] if len(L) >= 25 else L
                right_partial = R[:25]  if len(R) >= 25 else R

                # Prefilter: OR (fast), then real alignment decides
                if (left_partial and left_partial in seq) or (right_partial and right_partial in seq):
                    results[ref_name]["prefilter_hits"] += 1
                else:
                    continue

                # Try forward read
                aln_f = find_best_alignment(
                    seq, L, R, max_indel_size, min_flank_match, min_score
                )
                # Try reverse-complement read
                aln_r = find_best_alignment(
                    reverse_complement(seq), L, R, max_indel_size, min_flank_match, min_score
                )

                aln = None
                if aln_f and aln_r:
                    aln = aln_f if aln_f["score"] >= aln_r["score"] else aln_r
                else:
                    aln = aln_f or aln_r

                if not aln:
                    continue  # prefilter hit but not both flanks with required score

                results[ref_name]["both_flanks_hits"] += 1

                # Classify event
                if aln["type"] == "perfect":
                    results[ref_name]["indel_counts"]["perfect_match"] += 1
                    results[ref_name]["classified"] += 1

                elif aln["type"] == "deletion":
                    key = f"del_{aln['size']}bp"
                    results[ref_name]["indel_counts"][key] += 1
                    results[ref_name]["classified"] += 1

                elif aln["type"] == "insertion":
                    ins_seq = aln["sequence"]
                    key = f"ins_{aln['size']}bp"
                    # annotate short sequences
                    if ins_seq and len(ins_seq) <= 20:
                        key += f"_{ins_seq}"

                    # ODN tagging
                    if odn_seq and contains_odn(ins_seq, odn_seq, k=odn_k):
                        key += "_ODN"

                    results[ref_name]["indel_counts"][key] += 1
                    results[ref_name]["classified"] += 1

        return results

    finally:
        fh.close()

# ----------------------------
# Reporting
# ----------------------------

def print_results(results: dict):
    for ref_name, data in results.items():
        total = data["total_reads"]
        pre  = data["prefilter_hits"]
        both = data["both_flanks_hits"]
        cls  = data["classified"]
        events = sum(data["indel_counts"].values())
        edits = events - data["indel_counts"].get("perfect_match", 0)

        print("\n" + "="*60)
        print(f"Results for {ref_name}")
        print("="*60)
        print(f"Total reads processed:         {total:,}")
        print(f"Prefilter hits (≥1 flank):      {pre:,} ({(pre/total*100 if total else 0):.2f}%)")
        print(f"Both flanks aligned (≥{15}bp): {both:,} ({(both/total*100 if total else 0):.2f}%)")
        print(f"Classified (events called):     {cls:,} ({(cls/total*100 if total else 0):.2f}%)")

        if cls > 0:
            true_eff = edits / cls * 100.0
            print(f"True editing efficiency:        {true_eff:.2f}%  (#edited / classified)")
            print(f"Total events (incl. perfect):   {events:,}")
            print("\nEvent breakdown (as % of CLASSIFIED):")
            for et, c in data["indel_counts"].most_common():
                print(f"  {et}: {c:,} ({(c/cls*100):.2f}%)")
        else:
            print("No classifiable events.")

def save_results(results: dict, out_path: str, fastq: str, ref_fa: str, args):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        f.write(f"CRISPR Indel Detection Results - Sample: {os.path.basename(out_path).replace('_indels.txt','')}\n")
        f.write("="*80 + "\n")
        f.write(f"Input file: {fastq}\n")
        f.write(f"Reference file: {ref_fa}\n")
        f.write(f"Min quality: {args.min_quality}\n")
        f.write(f"Max indel size: {args.max_indel}\n")
        f.write(f"Min flank match: {args.min_flank}\n")
        f.write(f"Min score: {args.min_score}\n")
        if args.odn:
            f.write(f"ODN sequence: {args.odn} (k={args.odn_k})\n")
        f.write("\n")

        for ref_name, data in results.items():
            total = data["total_reads"]
            pre  = data["prefilter_hits"]
            both = data["both_flanks_hits"]
            cls  = data["classified"]
            events = sum(data["indel_counts"].values())
            edits  = events - data["indel_counts"].get("perfect_match", 0)

            f.write(f"\n{ref_name}\n")
            f.write("-"*len(ref_name) + "\n")
            f.write(f"Total reads:                   {total:,}\n")
            f.write(f"Prefilter hits (≥1 flank):     {pre:,} ({(pre/total*100 if total else 0):.2f}%)\n")
            f.write(f"Both flanks aligned:           {both:,} ({(both/total*100 if total else 0):.2f}%)\n")
            f.write(f"Classified (events called):    {cls:,} ({(cls/total*100 if total else 0):.2f}%)\n")

            if cls > 0:
                f.write(f"True editing efficiency:       {(edits/cls*100):.2f}%\n")

            f.write("\nEvents:\n")
            if events > 0:
                for et, c in data["indel_counts"].most_common():
                    pct = (c/cls*100) if cls else 0.0
                    f.write(f"  {et}: {c:,} ({pct:.2f}%)\n")
            else:
                f.write("  No editing events detected\n")

    print(f"\nResults saved to: {out_path}")

# ----------------------------
# CLI
# ----------------------------

def main():
    p = argparse.ArgumentParser(description="Detect CRISPR indels from merged FASTQ reads")
    p.add_argument("fastq", help="Input FASTQ (gz or plain)")
    p.add_argument("reference", help="Reference FASTA (∼200 bp windows)")
    p.add_argument("--min-quality", type=int, default=20, help="Min mean PHRED quality (default 20)")
    p.add_argument("--max-indel", type=int, default=100, help="Max indel size to detect (default 100)")
    p.add_argument("--min-flank", type=int, default=15, help="Min matching bases on EACH flank (default 15)")
    p.add_argument("--min-score", type=int, default=30, help="Min (left_len + right_len) to accept (default 30)")
    p.add_argument("--odn", type=str, default="GTTTAATTGAGTTGTCATATGTTAATAACGGTAT",
                   help="ODN template sequence (default: provided in your doc). Use '' to disable.")
    p.add_argument("--odn-k", type=int, default=20, help="Substring length to tag ODN insertions (default 20)")
    p.add_argument("--output", help="Write a human-readable report here")
    args = p.parse_args()

    print("\n" + "="*80)
    sample_name = (args.output or os.path.basename(args.fastq)).replace("_indels.txt", "")
    print(f"CRISPR Indel Detection Analysis - Sample: {sample_name}")
    print("="*80)

    # Load references and compute cut positions
    refs = load_references(args.reference)

    # Analyze
    results = analyze_reads(
        args.fastq,
        refs,
        min_quality=args.min_quality,
        max_indel_size=args.max_indel,
        min_flank_match=args.min_flank,
        min_score=args.min_score,
        odn_seq=args.odn,
        odn_k=args.odn_k,
    )
    if not results:
        print("ERROR: No results generated")
        sys.exit(1)

    # Print to stdout
    print_results(results)

    # Optional save
    if args.output:
        try:
            save_results(results, args.output, args.fastq, args.reference, args)
        except Exception as e:
            print(f"ERROR: Cannot write output file: {e}")

    print("\nCompleted.")
    print("-"*80)

if __name__ == "__main__":
    main()

