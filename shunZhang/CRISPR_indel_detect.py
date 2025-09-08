# --- add near top (helpers) ---
def get_anchor(seq, k, side, offset):
    """Return an anchor k-mer from a flank with a given offset away from the cut."""
    if side == 'L':
        if len(seq) < offset + k: 
            return None
        return seq[len(seq) - offset - k : len(seq) - offset]
    else:
        if len(seq) < offset + k:
            return None
        return seq[offset : offset + k]

# --- replace find_best_alignment with this version ---
def find_best_alignment(
    read_seq: str,
    ref_left: str,
    ref_right: str,
    max_indel_size: int = 100,
    min_flank_match: int = 15,
    min_score: int = 30,
    anchor_offsets=(0,),          # try 0 first; add 12,24,... to catch deletions chewing into proximal bases
    max_flank_len: int = 50,
):
    """
    Offset-aware anchor alignment across the cut.
    - We match a suffix of left flank and a prefix of right flank, optionally shifted
      away from the cut by l_off / r_off.
    - Expected spacing between matches with no indel is (l_off + r_off).
    - Insertion: actual_right_start > expected => delta > 0
    - Deletion:  actual_right_start < expected => delta < 0
    """
    read = read_seq.upper()
    L = ref_left.upper()
    R = ref_right.upper()

    best = None
    best_score = -1

    max_left = min(len(L), len(read) - min_flank_match, max_flank_len)
    if max_left < min_flank_match:
        return None
    max_right = min(len(R), len(read), max_flank_len)
    if max_right < min_flank_match:
        return None

    # First, offsets in the relaxed order (0 first for speed)
    l_off_choices = [o for o in anchor_offsets if o >= 0]
    r_off_choices = [o for o in anchor_offsets if o >= 0]

    for l_off in l_off_choices:
        for left_len in range(min_flank_match, max_left + 1):
            left_ref = get_anchor(L, left_len, 'L', l_off)
            if not left_ref:
                continue
            left_match_pos = read.find(left_ref)
            if left_match_pos == -1:
                continue
            left_end = left_match_pos + left_len

            for r_off in r_off_choices:
                # shift right anchor away from cut if needed
                for right_len in range(min_flank_match, max_right + 1):
                    right_ref = get_anchor(R, right_len, 'R', r_off)
                    if not right_ref:
                        continue

                    # window around where we expect the right anchor to land
                    expected_right_start = left_end + (l_off + r_off)
                    win_start = max(0, expected_right_start - max_indel_size)
                    win_end   = min(len(read), expected_right_start + max_indel_size)

                    rs = read.find(right_ref, win_start)
                    while rs != -1 and rs <= win_end:
                        delta = rs - expected_right_start  # signed indel
                        if abs(delta) <= max_indel_size:
                            score = left_len + right_len
                            if score >= min_score and score > best_score:
                                if delta > 0:
                                    call = {"type":"insertion","size": delta,"sequence": read[expected_right_start:rs],
                                            "score": score, "left_len": left_len, "right_len": right_len,
                                            "l_off": l_off, "r_off": r_off}
                                elif delta < 0:
                                    call = {"type":"deletion","size": -delta,"sequence": "",
                                            "score": score, "left_len": left_len, "right_len": right_len,
                                            "l_off": l_off, "r_off": r_off}
                                else:
                                    call = {"type":"perfect","size": 0,"sequence": "",
                                            "score": score, "left_len": left_len, "right_len": right_len,
                                            "l_off": l_off, "r_off": r_off}
                                best = call
                                best_score = score
                        # find next occurrence inside window
                        next_pos = read.find(right_ref, rs + 1)
                        if next_pos == -1 or next_pos > win_end:
                            break
                        rs = next_pos

    return best

# --- in analyze_reads(...), add CLI params and tweak prefilter ---
# Add args in argparse (main): 
#   p.add_argument("--anchor-offsets", default="0,12,24", help="Comma-separated offsets (bp) away from cut to try for flank anchors")
#   p.add_argument("--prefilter-k", type=int, default=15, help="k-mer length for prefilter (default 15)")

# Inside main(), pass anchor_offsets tuple into analyze_reads(...)

# In analyze_reads signature add new params:
#   anchor_offsets=(0,12,24), prefilter_k=15

# Replace the prefilter block with this:
anchor_offsets_tuple = tuple(int(x) for x in anchor_offsets if int(x) >= 0)

# ...
# Prefilter using multiple anchors per flank (OR)
prefilter_hit = False
for off in anchor_offsets_tuple:
    left_anchor  = get_anchor(L, prefilter_k, 'L', off)
    right_anchor = get_anchor(R, prefilter_k, 'R', off)
    if (left_anchor and left_anchor in seq) or (right_anchor and right_anchor in seq):
        prefilter_hit = True
        break
if not prefilter_hit:
    continue

# When calling alignment, pass the offsets:
aln_f = find_best_alignment(seq, L, R, max_indel_size, min_flank_match, min_score,
                            anchor_offsets=anchor_offsets_tuple)
aln_r = find_best_alignment(reverse_complement(seq), L, R, max_indel_size, min_flank_match, min_score,
                            anchor_offsets=anchor_offsets_tuple)

