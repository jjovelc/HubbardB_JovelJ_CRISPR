#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from math import isfinite

# Paths
TSV_PATH = "../on_target_sums_42samples_with_reads.tsv"
OUT_PNG = "../pct_cleaved_barplot.png"

def main(tsv_path=TSV_PATH, out_png=OUT_PNG):
    df = pd.read_csv(tsv_path, sep='\t')
    # Ensure the required columns exist
    required = {"sample","sample_gene","total_reads","on_target_total_all4","pct_nicked"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns in TSV: {missing}")

    # If pct_nicked is missing or NaN, recompute
    if df["pct_nicked"].isna().any():
        df["pct_nicked"] = (df["on_target_total_all4"] / df["total_reads"]).where(df["total_reads"]>0, 0) * 100

    # Sort by sample name for a stable order
    df = df.sort_values("sample").reset_index(drop=True)

    # Color palette by gene
    palette = {
        "AAVS1": "#1f77b4",   # blue
        "CEP290": "#2ca02c",  # green
        "TRAC": "#d62728",    # red
    }
    untreated_color = "#d3d3d3"  # light gray

    # Build colors per bar: untreated samples override to gray
    colors = []
    for s, g in zip(df["sample"], df["sample_gene"]):
        if s.startswith("Untreated_"):
            colors.append(untreated_color)
        else:
            colors.append(palette.get(g, "#7f7f7f"))  # fallback gray

    # Plot
    plt.figure(figsize=(18, 7))  # wide for 42 bars
    x = range(len(df))
    bars = plt.bar(x, df["pct_nicked"].values, edgecolor="black", linewidth=0.8, color=colors)

    # X labels
    plt.xticks(ticks=x, labels=df["sample"].tolist(), rotation=45, ha="right")
    plt.ylabel("% reads cleaved (on-target)")
    plt.title("On-target dsODN junction rate per library")

    # Build legend (gene colors + Untreated)
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor=palette["AAVS1"], edgecolor="black", label="AAVS1"),
        Patch(facecolor=palette["CEP290"], edgecolor="black", label="CEP290"),
        Patch(facecolor=palette["TRAC"], edgecolor="black", label="TRAC"),
        Patch(facecolor=untreated_color, edgecolor="black", label="Untreated"),
    ]
    plt.legend(handles=legend_handles, title="Guide / Condition", frameon=False, ncol=4, loc="upper right")

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    print(f"Saved figure to: {out_png}")

if __name__ == "__main__":
    # Allow overriding paths via CLI
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        main()
