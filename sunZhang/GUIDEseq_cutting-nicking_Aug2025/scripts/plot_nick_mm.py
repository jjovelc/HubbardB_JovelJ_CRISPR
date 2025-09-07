#!/usr/bin/env python3
import os, sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def infer_title(tsv_path: str, explicit: str | None = None) -> str:
    if explicit:
        return explicit
    p = Path(tsv_path)
    # Use parent directory name as the sample title (generic, repo-agnostic)
    return p.parent.name or p.name

def plot_mm_counts(tsv_path: str, out_png: str | None = None, title: str | None = None, ylog: bool = False):
    df = pd.read_csv(tsv_path, sep='\t')
    bins = [c for c in df.columns if c.startswith('MM')]
    title = infer_title(tsv_path, title)

    # 7 distinct colors (user request)
    colors = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2']

    fig, ax = plt.subplots(figsize=(10, 5.5))
    for i, row in df.iterrows():
        lib = row['lib']
        y = [row[b] for b in bins]
        ax.plot(range(len(bins)), y, marker='o', linewidth=1.8, label=lib, color=colors[i % 7])

    ax.set_xticks(range(len(bins)))
    ax.set_xticklabels([b.replace('MM','') for b in bins])
    ax.set_xlabel('Mismatch count (MM)')
    ax.set_ylabel('Reads')
    if ylog:
        ax.set_yscale('log')
    ax.set_title(title)
    ax.legend(title='Library', ncol=4, frameon=False)
    ax.grid(True, linestyle='--', alpha=0.3)
    fig.tight_layout()

    if out_png is None:
        out_png = os.path.join(os.path.dirname(tsv_path), 'nick_mm_lineplot.png')
    fig.savefig(out_png, dpi=200)
    return out_png

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1].startswith('-'):
        print('Usage: plot_nick_mm.py <nick_lib_mm_counts.tsv> [out.png] [--title "Sample"] [--ylog]', file=sys.stderr)
        sys.exit(1)
    tsv = sys.argv[1]
    out = None
    title = None
    ylog = False
    # Parse optional args
    for i, arg in enumerate(sys.argv[2:]):
        if arg == '--ylog':
            ylog = True
        elif arg == '--title' and i+3 <= len(sys.argv)-1:
            title = sys.argv[i+3]
        elif not arg.startswith('--') and out is None:
            out = arg
    out_path = plot_mm_counts(tsv, out, title, ylog=ylog)
    print(out_path)