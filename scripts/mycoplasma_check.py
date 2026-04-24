#!/usr/bin/env python3
#mycoplasma_check.py

import pandas as pd
from pathlib import Path
import sys

def check_mycoplasma_sequences(results_dir):
    mycoplasma_sequences = [
        'Mycoplasma_hyorhinis_MCLD',
        'Mycoplasma_hominis_ATCC23114',
        'Mycoplasma_fermentans_M64',
        'Acholeplasma_laidlawii_PG-8A'
    ]

    results_dir = Path(results_dir)
    samples = [d.name for d in (results_dir / "kallisto").iterdir() if d.is_dir()]

    report_data = []

    for sample in samples:
        sample_path = results_dir / "kallisto" / sample / f"{sample}_Foreign_sequences" / "abundance.tsv"
        if not sample_path.exists():
            print(f"Error: File '{sample_path}' not found.")
            continue

        df = pd.read_table(sample_path, sep='\t')

        sample_report = {'Sample': sample}
        for sequence in mycoplasma_sequences:
            if sequence in df['target_id'].values:
                est_counts = df.loc[df['target_id'] == sequence, 'est_counts'].values[0]
                sample_report[sequence] = 'positive' if est_counts > 2000 else 'negative'
            else:
                sample_report[sequence] = 'negative'

        report_data.append(sample_report)

    report_df = pd.DataFrame(report_data)
    output_path = results_dir / "mycoplasma_report.tsv"
    report_df.to_csv(output_path, sep='\t', index=False)
    print(f"Mycoplasma report created at '{output_path}'")

def main():
    if len(sys.argv) != 2:
        print("Usage: python mycoplasma_check.py <results_dir>")
        sys.exit(1)

    results_dir = Path(sys.argv[1])
    check_mycoplasma_sequences(results_dir)

if __name__ == "__main__":
    main()
