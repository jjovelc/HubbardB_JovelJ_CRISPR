#!/usr/bin/env python3
import os, sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from typing import List

REQUIRED_COLUMNS = ["lib","total"] + [f"MM{i}" for i in range(0,15)]
LIBS = [f"Lib{i}" for i in range(1,8)]
CONTROLS = ["nlib_RNA_5pmol_1","nlib_RNA_5pmol_2","preselection_lib_1","preselection_lib_2"]

def read_mm_table(tsv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t")
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {tsv_path}: {missing}")
    return df

def mm_series_for_lib(df: pd.DataFrame, lib: str) -> List[int]:
    row = df.loc[df["lib"]==lib]
    if row.empty:
        raise ValueError(f"Library {lib} not found in table {df}")
    bins = [c for c in df.columns if c.startswith("MM")]
    return [int(row.iloc[0][b]) for b in bins]

def plot_one(base: str, root: Path, ylog: bool=False, tight: bool=False) -> Path:
    # Files to load (6 datasets): base_1, base_2 + controls
    samples = [f"{base}_1", f"{base}_2", "nlib_RNA_5pmol_1", "nlib_RNA_5pmol_2", "preselection_lib_1", "preselection_lib_2"]
    labels = [f"{base}_1", f"{base}_2", "RNA_5pmol_1", "RNA_5pmol_2", "preselection_1", "preselection_2"]
    paths = [root / s / "nick_lib_mm_counts.tsv" for s in samples]

    # Read all
    dfs = []
    for p in paths:
        if not p.exists():
            raise FileNotFoundError(f"Missing input: {p}")
        dfs.append(read_mm_table(p))

    outdir = root / f"{base}_plots"
    outdir.mkdir(parents=True, exist_ok=True)

    bins = [f"MM{i}" for i in range(0,15)]
    x = list(range(len(bins)))
    xticklabels = [str(i) for i in range(0,15)]

    for lib in LIBS:
        fig, ax = plt.subplots(figsize=(9.5, 5.3))
        for df, lbl in zip(dfs, labels):
            y = mm_series_for_lib(df, lib)
            ax.plot(x, y, marker="o", linewidth=1.8, label=lbl)
        ax.set_xticks(x)
        ax.set_xticklabels(xticklabels)
        ax.set_xlabel("Mismatch count (MM)")
        ax.set_ylabel("Reads")
        if ylog:
            ax.set_yscale("log")
        ax.set_title(f"{base} — {lib}")
        ax.legend(ncol=3, frameon=False)
        ax.grid(True, linestyle="--", alpha=0.3)
        if tight:
            fig.tight_layout()
        out_png = outdir / f"{base}_{lib}.png"
        fig.savefig(out_png, dpi=200)
        plt.close(fig)
    return outdir

def main():
    import argparse
    ap = argparse.ArgumentParser(description="Plot Lib1..Lib7 MM curves comparing duplicates vs controls for each base sample.")
    ap.add_argument("--root", required=True, help="Parent directory containing sample subdirectories")
    ap.add_argument("--base", action="append", required=True, help="Base sample name (e.g., nlib_LNA_9_5pmol). Repeat for multiple bases.")
    ap.add_argument("--ylog", action="store_true", help="Log-scale y-axis")
    ap.add_argument("--tight", action="store_true", help="Tight layout")
    args = ap.parse_args()

    root = Path(args.root)
    for b in args.base:
        outdir = plot_one(b, root, ylog=args.ylog, tight=args.tight)
        print(outdir)

if __name__ == "__main__":
    main()
