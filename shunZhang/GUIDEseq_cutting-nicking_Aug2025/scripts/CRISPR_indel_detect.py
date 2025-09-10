#!/usr/bin/env python3
"""
CRISPR_indel_detect.py  —  deletion-sensitive, ODN-aware indel caller for merged FASTQ

Key features
- Streams merged FASTQ (gz OK), optional mean-quality filter.
- Flexible, offset-aware anchoring around the Cas9 cut (catches deletions that
  remove proximal bases).
- Event calling: perfect | insertion | deletion with sizes; insertions get
  sequence tags (<=20 nt) and "_ODN" if they contain an ODN k-mer (either strand).
- Clear metrics: Prefilter hits, Both flanks aligned, Classified, True editing efficiency.
- Robust FASTA parsing; tries to locate cut by guide (both strands) else falls back to midpoint.

CLI tips
- Most impactful knobs for deletion sensitivity: --anchor-offsets, --prefilter-k,
  --min-flank, --min-score (try a slightly lower min-score like 28).
"""

import os
import re
import sys
import gzip
import argparse
from collections import Counter, defaultdict

# ------------------------------ FASTA / GUIDES ------------------------------

def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

def parse_reference_sequences(ref_file: str):
    """Manual FASTA parser. Returns dict[id] = {sequence, cut_position, left_flank, right_flank}."""
    if not os.path.exists(ref_file):
        print(f"ERROR: Reference file not found: {ref_file}")
        sys.exit(1)

    with open(ref_file, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read().strip()

    if not content.startswith(">"):
        print("ERROR: Reference file should be FASTA (starts with '>').")
        sys.exit(1)

    seqs = {}
    cur_id, cur_seq = None, []
    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur_id:
                seqs[cur_id] = "".join(cur_seq).upper().replace(" ", "")
            cur_id = line[1:].split()[0]
            cur_seq = []
        else:
            cur_seq.append(line)
    if cur_id:
        seqs[cur_id] = "".join(cur_seq).upper().replace(" ", "")

    print(f"Found {len(seqs)} reference sequences")

    # Known protospacers (20nt) without PAM; include common name variants
    guides = {
        "AAVS1": "GTCCCTAGTGGCCCCACTGT",
        "AASV1": "GTCCCTAGTGGCCCCACTGT",  # handle typo/variant
        "CEP290": "GGAGTCACATGGGAGTCACA",
        # TRAC protospacer (common TCR alpha guide; 20nt)
        "TRAC":   "GCTGGTACACGGCAGGGTCA",
    }

    references = {}
    for seq_id, sequence in seqs.items():
        cut_pos = len(sequence) // 2  # fallback

        # Try to find matching guide by name hint
        target_key = None
        up = seq_id.upper()
        for k in guides:
            if k in up:
                target_key = k
                break

        if target_key:
            g = guides[target_key]
            rcg = reverse_complement(g)
            pos = sequence.find(g)
            if pos != -1:
                # cut between bases 17|18 of protospacer in forward orientation
                cut_pos = pos + 17
                print(f"Guide {target_key} found in {seq_id} at {pos}; cut @ {cut_pos}")
            else:
                pos = sequence.find(rcg)
                if pos != -1:
                    # Treat similarly; this heuristically places the cut
                    cut_pos = pos + 17
                    print(f"Guide {target_key} (RC) found in {seq_id} at {pos}; cut @ {cut_pos}")
                else:
                    print(f"Guide {target_key} not found in {seq_id}; using midpoint @ {cut_pos}")
        else:
            print(f"No guide key matched for {seq_id}; using midpoint @ {cut_pos}")

        cut_pos = max(1, min(cut_pos, len(sequence) - 1))
        left_flank, right_flank = sequence[:cut_pos], sequence[cut_pos:]

        references[seq_id] = {
            "sequence": sequence,
            "cut_position": cut_pos,
            "left_flank": left_flank,
            "right_flank": right_flank,
        }
        print(f"  {seq_id}: {len(sequence)} bp, cut @ {cut_pos}")

    return references

# ------------------------------ ODN detection ------------------------------

class DonorMatcher:
    """Detects whether an inserted fragment contains any k-mer from ODN (either strand)."""
    def __init__(self, odn_seq: str | None, k: int = 20):
        self.enabled = bool(odn_seq)
        self.k = k
        self.kmers = set()
        if self.enabled:
            d = odn_seq.upper().replace(" ", "")
            rc = reverse_complement(d)
            if k <= len(d):
                for i in range(len(d) - k + 1):
                    self.kmers.add(d[i:i+k])
                for i in range(len(rc) - k + 1):
                    self.kmers.add(rc[i:i+k])

    def has_kmer(self, ins: str) -> bool:
        if not self.enabled or not ins:
            return False
        s = ins.upper()
        k = self.k
        if len(s) < k:
            return False
        for i in range(len(s) - k + 1):
            if s[i:i+k] in self.kmers:
                return True
        return False

# ------------------------------ Alignment core ------------------------------

def get_anchor(seq: str, k: int, side: str, offset: int) -> str | None:
    """Return an anchor k-mer from a flank with a given offset away from the cut."""
    if k <= 0 or offset < 0 or len(seq) < k + offset:
        return None
    if side == "L":
        return seq[len(seq) - offset - k : len(seq) - offset]
    else:
        return seq[offset : offset + k]

def find_best_alignment(
    read_seq: str,
    ref_left: str,
    ref_right: str,
    *,
    max_indel_size: int = 100,
    min_flank_match: int = 15,
    min_score: int = 30,
    anchor_offsets: tuple[int, ...] = (0, 12, 24),
    max_flank_len: int = 50,
):
    """
    Offset-aware dual-anchor exact matching across the cut.
    Tries suffix (left) and prefix (right) anchors possibly shifted away from the cut
    to tolerate deletions that remove proximal bases.
    """
    read = read_seq.upper()
    L = ref_left.upper()
    R = ref_right.upper()

    best = None
    best_score = -1

    max_left = min(len(L), len(read) - min_flank_match, max_flank_len)
    max_right = min(len(R), len(read), max_flank_len)
    if max_left < min_flank_match or max_right < min_flank_match:
        return None

    l_off_choices = [o for o in anchor_offsets if o >= 0]
    r_off_choices = [o for o in anchor_offsets if o >= 0]

    for l_off in l_off_choices:
        for left_len in range(min_flank_match, max_left + 1):
            left_ref = get_anchor(L, left_len, "L", l_off)
            if not left_ref:
                continue
            left_match_pos = read.find(left_ref)
            if left_match_pos == -1:
                continue
            left_end = left_match_pos + left_len

            for r_off in r_off_choices:
                for right_len in range(min_flank_match, max_right + 1):
                    right_ref = get_anchor(R, right_len, "R", r_off)
                    if not right_ref:
                        continue

                    expected_right_start = left_end + (l_off + r_off)
                    win_start = max(0, expected_right_start - max_indel_size)
                    win_end = min(len(read), expected_right_start + max_indel_size)

                    rs = read.find(right_ref, win_start)
                    while rs != -1 and rs <= win_end:
                        delta = rs - expected_right_start  # signed indel
                        if abs(delta) <= max_indel_size:
                            score = left_len + right_len
                            if score >= min_score and score > best_score:
                                if delta > 0:
                                    call = {
                                        "type": "insertion",
                                        "size": delta,
                                        "sequence": read[expected_right_start:rs],
                                        "score": score,
                                        "left_len": left_len,
                                        "right_len": right_len,
                                        "l_off": l_off,
                                        "r_off": r_off,
                                    }
                                elif delta < 0:
                                    call = {
                                        "type": "deletion",
                                        "size": -delta,
                                        "sequence": "",
                                        "score": score,
                                        "left_len": left_len,
                                        "right_len": right_len,
                                        "l_off": l_off,
                                        "r_off": r_off,
                                    }
                                else:
                                    call = {
                                        "type": "perfect",
                                        "size": 0,
                                        "sequence": "",
                                        "score": score,
                                        "left_len": left_len,
                                        "right_len": right_len,
                                        "l_off": l_off,
                                        "r_off": r_off,
                                    }
                                best = call
                                best_score = score

                        nxt = read.find(right_ref, rs + 1)
                        if nxt == -1 or nxt > win_end:
                            break
                        rs = nxt

    return best

# ------------------------------ FASTQ streaming ------------------------------

def open_textmaybe_gz(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def mean_phred(line: str) -> float:
    # Sanger / Illumina 1.8+: ASCII 33 offset
    if not line:
        return 0.0
    total = sum((ord(c) - 33) for c in line.rstrip("\n\r"))
    return total / max(1, len(line.strip()))

# ------------------------------ Main analysis ------------------------------

def analyze_reads(
    fastq_file: str,
    references: dict,
    *,
    min_quality: int = 20,
    max_indel_size: int = 100,
    min_flank_match: int = 15,
    min_score: int = 30,
    anchor_offsets: tuple[int, ...] = (0, 12, 24),
    prefilter_k: int = 15,
    donor: DonorMatcher | None = None,
):
    """Return per-target results dict."""
    if not os.path.exists(fastq_file):
        print(f"ERROR: FASTQ file not found: {fastq_file}")
        return {}

    results = {}

    for ref_name, ref_data in references.items():
        L = ref_data["left_flank"]
        R = ref_data["right_flank"]

        print(f"\nAnalyzing {ref_name}...")
        indel_counts = Counter()
        total_reads = 0
        prefilter_hits = 0
        both_flanks_aligned = 0
        classified = 0

        try:
            fh = open_textmaybe_gz(fastq_file)
            line_no = 0
            read_seq = None
            qual = None

            # Precompute anchors used in prefilter
            offsets = tuple(o for o in anchor_offsets if o >= 0)
            prefilter_anchors = []
            for off in offsets:
                la = get_anchor(L, prefilter_k, "L", off)
                ra = get_anchor(R, prefilter_k, "R", off)
                if la:
                    prefilter_anchors.append(la)
                if ra:
                    prefilter_anchors.append(ra)
            prefilter_anchors = [a for a in prefilter_anchors if a]

            for line in fh:
                line_no += 1
                mod = line_no % 4
                if mod == 2:       # sequence
                    read_seq = line.strip().upper()
                elif mod == 0:     # quality line
                    qual = line.rstrip("\n\r")
                    total_reads += 1

                    # Quality filter
                    if min_quality > 0 and mean_phred(qual) < min_quality:
                        continue

                    # Prefilter: OR over several anchors (forward only, fast)
                    if prefilter_anchors:
                        hit = any(a in read_seq for a in prefilter_anchors)
                        if not hit:
                            continue
                        prefilter_hits += 1
                    else:
                        prefilter_hits += 1  # degenerate but consistent

                    # Alignment: forward then reverse
                    aln = find_best_alignment(
                        read_seq, L, R,
                        max_indel_size=max_indel_size,
                        min_flank_match=min_flank_match,
                        min_score=min_score,
                        anchor_offsets=anchor_offsets
                    )

                    if not aln or aln.get("score", 0) < min_score:
                        rc = reverse_complement(read_seq)
                        aln_rc = find_best_alignment(
                            rc, L, R,
                            max_indel_size=max_indel_size,
                            min_flank_match=min_flank_match,
                            min_score=min_score,
                            anchor_offsets=anchor_offsets
                        )
                        if aln_rc and (not aln or aln_rc["score"] > aln["score"]):
                            aln = aln_rc

                    if aln and aln.get("score", 0) >= min_score:
                        both_flanks_aligned += 1
                        classified += 1

                        if aln["type"] == "perfect":
                            indel_counts["perfect_match"] += 1

                        elif aln["type"] == "insertion":
                            size = aln["size"]
                            seq = aln["sequence"]
                            key = f"ins_{size}bp"
                            # include sequence tag if short
                            if seq and len(seq) <= 20:
                                key += f"_{seq}"
                            # add ODN tag if donor k-mer present
                            if donor and donor.has_kmer(seq):
                                key += "_ODN"
                            indel_counts[key] += 1

                        elif aln["type"] == "deletion":
                            size = aln["size"]
                            key = f"del_{size}bp"
                            indel_counts[key] += 1

                    # progress
                    if total_reads % 50000 == 0:
                        print(f"  Processed {total_reads:,} reads  |  prefilter hits {prefilter_hits:,}  |  aligned {both_flanks_aligned:,}")

            fh.close()

        except Exception as e:
            print(f"ERROR while reading FASTQ: {e}")
            return {}

        results[ref_name] = {
            "total_reads": total_reads,
            "prefilter_hits": prefilter_hits,
            "both_flanks_aligned": both_flanks_aligned,
            "classified": classified,
            "indel_counts": indel_counts,
        }

        print(f"  Total reads:            {total_reads:,}")
        print(f"  Prefilter hits (≥1 flank): {prefilter_hits:,} ({(prefilter_hits/max(1,total_reads))*100:.2f}%)")
        print(f"  Both flanks aligned:    {both_flanks_aligned:,} ({(both_flanks_aligned/max(1,total_reads))*100:.2f}%)")
        print(f"  Classified:             {classified:,} ({(classified/max(1,total_reads))*100:.2f}%)")
        print(f"  Reads with events:      {sum(indel_counts.values()):,}")

    return results

# ------------------------------ Reporting ------------------------------

def print_results(results: dict, *, min_flank_match: int, min_score: int, odn_seq: str | None, odn_k: int):
    for ref_name, data in results.items():
        print("\n" + "=" * 80)
        print(f"Results for {ref_name}")
        print("=" * 80)

        total_reads = data["total_reads"]
        prefilter_hits = data["prefilter_hits"]
        both_flanks_aligned = data["both_flanks_aligned"]
        classified = data["classified"]
        counts = data["indel_counts"]

        print(f"Total reads:                   {total_reads:,}")
        print(f"Prefilter hits (≥1 flank):     {prefilter_hits:,} ({(prefilter_hits/max(1,total_reads))*100:.2f}%)")
        print(f"Both flanks aligned:           {both_flanks_aligned:,} ({(both_flanks_aligned/max(1,total_reads))*100:.2f}%)")
        print(f"Classified (events called):    {classified:,} ({(classified/max(1,total_reads))*100:.2f}%)")

        total_events = sum(counts.values())
        if classified > 0:
            perfect = counts.get("perfect_match", 0)
            true_editing = 100.0 * (classified - perfect) / classified
            print(f"True editing efficiency:       {true_editing:.2f}%")

        print("\nEvents:")
        if total_events == 0:
            print("  No editing events detected")
        else:
            for ev, c in counts.most_common():
                pct = (c / max(1, classified)) * 100.0
                print(f"  {ev}: {c:,} ({pct:.2f}%)")

# ------------------------------ Save to file ------------------------------

def save_results(out_path: str, sample_name: str, args, references: dict, results: dict, odn_seq: str | None):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        f.write(f"CRISPR Indel Detection Results - Sample: {sample_name}\n")
        f.write("=" * 80 + "\n")
        f.write(f"Input file: {args.fastq}\n")
        f.write(f"Reference file: {args.reference}\n")
        f.write(f"Min quality: {args.min_quality}\n")
        f.write(f"Max indel size: {args.max_indel}\n")
        f.write(f"Min flank match: {args.min_flank}\n")
        f.write(f"Min score: {args.min_score}\n")
        if odn_seq:
            f.write(f"ODN sequence: {odn_seq} (k={args.odn_k})\n")
        f.write("\n")

        for ref_name, data in results.items():
            f.write(f"{ref_name}\n")
            f.write("-" * len(ref_name) + "\n")
            total_reads = data["total_reads"]
            prefilter_hits = data["prefilter_hits"]
            both_flanks_aligned = data["both_flanks_aligned"]
            classified = data["classified"]
            counts = data["indel_counts"]

            f.write(f"Total reads:                   {total_reads:,}\n")
            f.write(f"Prefilter hits (≥1 flank):     {prefilter_hits:,} ({(prefilter_hits/max(1,total_reads))*100:.2f}%)\n")
            f.write(f"Both flanks aligned:           {both_flanks_aligned:,} ({(both_flanks_aligned/max(1,total_reads))*100:.2f}%)\n")
            f.write(f"Classified (events called):    {classified:,} ({(classified/max(1,total_reads))*100:.2f}%)\n")

            if classified > 0:
                perfect = counts.get("perfect_match", 0)
                true_edit = 100.0 * (classified - perfect) / classified
                f.write(f"True editing efficiency:       {true_edit:.2f}%\n")

            f.write("\nEvents:\n")
            if sum(counts.values()) == 0:
                f.write("  No editing events detected\n\n")
            else:
                for ev, c in counts.most_common():
                    pct = (c / max(1, classified)) * 100.0
                    f.write(f"  {ev}: {c:,} ({pct:.2f}%)\n")
                f.write("\n")

    print(f"\nResults saved to: {out_path}")

# ------------------------------ CLI ------------------------------

def main():
    p = argparse.ArgumentParser(description="Detect CRISPR indels from merged FASTQ reads")
    p.add_argument("fastq", help="Input FASTQ (merged; .gz OK)")
    p.add_argument("reference", help="FASTA with ~200 bp windows around cut (clean; no '-')")
    p.add_argument("--min-quality", type=int, default=20, help="Minimum mean read quality (Phred)")
    p.add_argument("--max-indel", type=int, default=100, help="Maximum indel size to detect")
    p.add_argument("--min-flank", type=int, default=15, help="Minimum matched bases per flank")
    p.add_argument("--min-score", type=int, default=30, help="Minimum left+right match length to accept an alignment")
    p.add_argument("--prefilter-k", type=int, default=15, help="k-mer for prefilter anchors (OR over offsets)")
    p.add_argument("--anchor-offsets", default="0,12,24",
                   help="Comma-separated offsets (bp) away from the cut to place anchors (e.g., 0,8,16,24)")
    p.add_argument("--odn", default="GTTTAATTGAGTTGTCATATGTTAATAACGGTAT",
                   help="ODN sequence (leave empty to disable ODN tagging)")
    p.add_argument("--odn-k", type=int, default=20, help="ODN k-mer length for tagging (default 20)")
    p.add_argument("--output", help="Write a results text file (recommended)")
    args = p.parse_args()

    sample_name = os.path.basename(args.output).replace("_indels.txt", "") if args.output else os.path.basename(args.fastq)

    print("\n" + "=" * 80)
    print(f"CRISPR Indel Detection Analysis - Sample: {sample_name}")
    print("=" * 80)

    print("Loading reference sequences...")
    references = parse_reference_sequences(args.reference)

    # Donor matcher (ODN tagging)
    odn_seq = args.odn.strip() if args.odn else None
    donor = DonorMatcher(odn_seq, k=max(1, args.odn_k)) if odn_seq else None

    # Parse offsets
    try:
        anchor_offsets = tuple(int(x) for x in args.anchor_offsets.split(",") if x.strip() != "")
    except Exception:
        print("ERROR: --anchor-offsets must be a comma-separated list of integers, e.g., 0,12,24")
        sys.exit(1)

    # Analyze
    print(f"\nAnalyzing reads from: {args.fastq}")
    results = analyze_reads(
        args.fastq,
        references,
        min_quality=args.min_quality,
        max_indel_size=args.max_indel,
        min_flank_match=args.min_flank,
        min_score=args.min_score,
        anchor_offsets=anchor_offsets,
        prefilter_k=args.prefilter_k,
        donor=donor,
    )
    if not results:
        print("ERROR: No results generated")
        sys.exit(1)

    # Print to stdout
    print_results(results, min_flank_match=args.min_flank, min_score=args.min_score, odn_seq=odn_seq, odn_k=args.odn_k)

    # Save results
    if args.output:
        save_results(args.output, sample_name, args, references, results, odn_seq)

    print("\nCompleted analysis.")
    print("-" * 80)

if __name__ == "__main__":
    main()

