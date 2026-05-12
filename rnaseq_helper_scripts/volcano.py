#!/usr/bin/env python3
#volcano.py

import logging
from pathlib import Path
import sys

import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/.config/matplotlib")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
import matplotlib

matplotlib.use("Agg")
# Avoid findfont spam when Arial (or other UI fonts) are missing (e.g. Linux/Docker).
matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Bitstream Vera Sans", "sans-serif"],
    }
)
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid", font="DejaVu Sans")

def generate_volcano_plot(results_dir):
    results_dir = Path(results_dir)
    deseq2_results_path = results_dir / "deseq2" / "deseq2_results.csv"
    output_plot_path = results_dir / "deseq2" / "volcano_plot.png"

    # Read the DESeq2 results file
    df = pd.read_csv(deseq2_results_path, sep=',')

    # Define significance cut-offs
    log2fc_cutoff = 2
    padj_cutoff = 0.05

    # Create a new column for significance
    df['significance'] = 'Not Significant'
    df.loc[(df['log2FoldChange'] >= log2fc_cutoff) & (df['padj'] < padj_cutoff), 'significance'] = 'Upregulated'
    df.loc[(df['log2FoldChange'] <= -log2fc_cutoff) & (df['padj'] < padj_cutoff), 'significance'] = 'Downregulated'

    # Set colors for significance
    colors = {'Not Significant': 'grey', 'Upregulated': 'c', 'Downregulated': 'salmon'}

    # Create the volcano plot
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=df, x='log2FoldChange', y=-np.log10(df['padj']), hue='significance', palette=colors, edgecolor=None)
    plt.axhline(y=-np.log10(padj_cutoff), color='black', linestyle='--')
    plt.axvline(x=log2fc_cutoff, color='black', linestyle='--')
    plt.axvline(x=-log2fc_cutoff, color='black', linestyle='--')
    plt.title('Volcano Plot of DESeq2 Results')
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 Adjusted P-Value')
    plt.legend(title='Significance')
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_plot_path)
    plt.close()
    print(f"Volcano plot saved at '{output_plot_path}'")

def main():
    if len(sys.argv) != 2:
        print("Usage: python volcano.py <results_dir>")
        sys.exit(1)

    results_dir = Path(sys.argv[1])
    generate_volcano_plot(results_dir)

if __name__ == "__main__":
    main()
