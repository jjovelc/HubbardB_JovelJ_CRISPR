#!/usr/bin/env python3
"""
plot_indel_sizes.py

Read a CRISPR_indel_detect results file and plot a signed indel size distribution.
- Negative sizes = deletions (from "del_XXbp" bins)
- Positive sizes = insertions (from "ins_XXbp" bins)
- Ignores "perfect_match" rows
- Robust to suffixes (e.g., _ODN, sequence tags)

Outputs:
- PNG barplot of -X..+X (default X=100)
- CSV table of size and count
"""

import re
import os
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd

HEADER_PATTERN = re.compile(r'^\s*([A-Za-z0-9_]+)\s*\n-{3,}\s*$', re.MULTILINE)
EVENTS_HEADER_RE = re.compile(r'^\s*Events:\s*$', re.MULTILINE)
CLASSIFIED_RE = re.compile(r'Classified\s*\(events called\)\s*:\s*([\d,]+)', re.IGNORECASE)

# Lines like:
#   "  ins_34bp_ODN: 550,216 (24.48%)"
#   "  del_15bp: 3,210 (0.14%)"
EVENT_LINE_RE = re.compile(
    r'^\s*(ins|del)_(\d+)bp[^:]*:\s*([\d,]+)',
    re.IGNORECASE | re.MULTILINE
)

# Canonical protospacers (20nt, no PAM); tolerate section typos & rc suffixes
GUIDES = {
    "AAVS1": "GTCCCTAGTGGCCCCACTGT",
    "CEP290": "GGAGTCACATGGGAGTCACA",
    "TRAC":   "GCTGGTACACGGCAGGGTCA",
}
ALIASES = {
    "AASV1": "AAVS1",   # common typo
}

def guide_for_section(section_name: str):
    """
    Map a section name (e.g., 'AASV1', 'CEP290rc') to (display_name, guide_seq, is_rc).
    """
    u = section_name.strip().upper()
    is_rc = False
    if u.endswith("RC"):
        is_rc = True
        u = u[:-2]
    u = ALIASES.get(u, u)
    return u, GUIDES.get(u), is_rc

def read_text(path: str) -> str:
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        return f.read()

def find_sections(text: str):
    """Return list of (name, start_idx, end_idx) for blocks like:
       AAVS1
       -----
       ..."""
    sections = []
    matches = list(HEADER_PATTERN.finditer(text))
    for i, m in enumerate(matches):
        name = m.group(1).strip()
        start = m.end()
        end = matches[i+1].start() if i + 1 < len(matches) else len(text)
        sections.append((name, start, end))
    return sections

def extract_events_block(section_text: str) -> str | None:
    m = EVENTS_HEADER_RE.search(section_text)
    if not m:
        return None
    start = m.end()
    return section_text[start:]  # up to section end

def parse_event_counts(events_block: str):
    """Return:
       counts: dict[int, int] mapping signed size to count
       total:  sum of parsed counts (excludes perfect_match)"""
    counts = defaultdict(int)
    total = 0
    for m in EVENT_LINE_RE.finditer(events_block):
        kind = m.group(1).lower()
        size = int(m.group(2))
        count = int(m.group(3).replace(',', ''))
        signed = -size if kind == 'del' else size
        counts[signed] += count
        total += count
    return counts, total

def _canonical_name(x: str) -> str:
    u = x.strip().upper()
    if u.endswith("RC"):           # allow CEP290 ↔ CEP290rc, TRAC ↔ TRACrc
        u = u[:-2]
    if u == "AASV1":               # tolerate typo in sections
        u = "AAVS1"
    if u == "TRAP":                # (safe guard; won't trigger if you don't have TRAP)
        u = "TRAC"
    return u

def pick_target_section(text: str, preferred: str | None):
    sections = find_sections(text)
    if not sections:
        raise ValueError("No target sections (e.g., AAVS1, CEP290rc) found in file.")
    if preferred:
        pref_can = _canonical_name(preferred)
        for name, s, e in sections:
            if _canonical_name(name) == pref_can:
                return name, text[s:e]
        # If not found, show what *was* present
        avail = ", ".join(n for n, _, _ in sections)
        raise ValueError(f"Target '{preferred}' not found. Available: {avail}")
    # Fallback: auto-pick busiest section
    best = None
    best_total = -1
    for name, s, e in sections:
        block = extract_events_block(text[s:e]) or ""
        _, total = parse_event_counts(block)
        if total > best_total:
            best = (name, text[s:e])
            best_total = total
    return best




def main():
    ap = argparse.ArgumentParser(description="Plot signed indel size distribution from results file.")
    ap.add_argument("results_file", help="Path to a CRISPR indel results text file")
    ap.add_argument("--target", help="Target to plot (e.g., AAVS1, CEP290rc). If omitted, picks the busiest section.")
    ap.add_argument("--xmax", type=int, default=100, help="Half-range to plot: [-X..+X]")
    ap.add_argument("--outdir", default="", help="Output directory (default: alongside input)")
    ap.add_argument("--normalize", choices=["none","classified","events"], default="none",
                    help="none=counts; classified=%% of 'Classified (events called)'; events=%% of parsed events")
    ap.add_argument("--title", default="", help="Custom plot title (if set, overrides filename-based default)")
    ap.add_argument("--basename", default="", help="Base name for output files (PNG/CSV). Defaults to input+target.")
    args = ap.parse_args()

    # Derive the filename-based title once, up front
    input_base_full = os.path.basename(args.results_file)
    if input_base_full.endswith("_indels.txt"):
        input_base = input_base_full[:-len("_indels.txt")]
    else:
        input_base = os.path.splitext(input_base_full)[0]

    text = read_text(args.results_file)
    target_name, section_text = pick_target_section(text, args.target)
    events_block = extract_events_block(section_text)
    if not events_block:
        raise ValueError(f"No 'Events:' block found under target '{target_name}'.")

    counts, total_parsed = parse_event_counts(events_block)

    # Denominator for normalization
    denom = None
    if args.normalize == "classified":
        m = CLASSIFIED_RE.search(section_text)
        denom = int(m.group(1).replace(",", "")) if m else total_parsed
    elif args.normalize == "events":
        denom = total_parsed

    # Build data frame for -X..+X inclusive
    X = max(1, int(args.xmax))
    xs = list(range(-X, X+1))
    ys = []
    for x in xs:
        v = counts.get(x, 0)
        if denom and denom > 0:
            v = v / denom * 100.0
        ys.append(v)
    df = pd.DataFrame({"size_bp": xs, "value": ys})
    if args.normalize == "none":
        df.rename(columns={"value": "count"}, inplace=True)
        value_label = "count"
    else:
        df.rename(columns={"value": "percent"}, inplace=True)
        value_label = "percent"

    # Output paths
    outdir = args.outdir or os.path.dirname(os.path.abspath(args.results_file))
    os.makedirs(outdir, exist_ok=True)
    base = args.basename or f"{input_base}_{target_name}"
    png_path = os.path.join(outdir, f"{base}_indel_size_barplot.png")
    csv_path = os.path.join(outdir, f"{base}_indel_size_counts.csv")

    # Legend labels (target + guide) and ODN marker
    canonical_name, guide_seq, is_rc = guide_for_section(target_name)
    if guide_seq:
        legend_label = f"{canonical_name}{' (rc)' if is_rc else ''} — guide: {guide_seq}"
    else:
        legend_label = f"{target_name}"

    # Plot
    plt.figure(figsize=(12, 5))
    plt.bar(df["size_bp"], df.iloc[:, 1], label=legend_label)  # default matplotlib color
    vline_label = None
    if 34 <= X:
        vline_label = "+34 (ODN)"
        plt.axvline(34, linestyle='--', label=vline_label, color='red')

    # Title: default to the input filename minus _indels.txt (unless user passes --title)
    plt.title(args.title if args.title else input_base)
    plt.xlabel("Indel size (bp; negative = deletions, positive = insertions)")
    plt.ylabel(f"Read {value_label}")

    # Show legend reflecting which guide the counts belong to
    plt.legend()

    plt.xlim([-X, X])
    plt.tight_layout()
    plt.savefig(png_path, dpi=150)

    # Save table
    df.to_csv(csv_path, index=False)

    print(f"[ok] Wrote: {png_path}")
    print(f"[ok] Wrote: {csv_path}")

if __name__ == "__main__":
    main()

