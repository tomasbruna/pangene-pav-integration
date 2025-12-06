#!/usr/bin/env python3

import pandas as pd
import sys

def classify_pangenome(filepath):
    df = pd.read_csv(filepath, sep='\t')
    genome_counts = df.groupby('globHOG')['genome'].nunique().reset_index()
    genome_counts.columns = ['globHOG', 'nGenome']
    df = df.merge(genome_counts, on='globHOG')
    total_genomes = df['genome'].nunique()
    df['propGenome'] = df['nGenome'] / total_genomes

    def classify(row):
        if row['propGenome'] == 1:
            return 'core'
        elif row['propGenome'] >= 0.9:
            return 'softcore'
        elif row['propGenome'] >= 0.5:
            return 'shell'
        elif row['nGenome'] == 1:
            return 'private'
        else:
            return 'cloud'

    df['cls'] = df.apply(classify, axis=1)

    # Summarize by category
    total_genes = len(df)
    total_ogs = df['globHOG'].nunique()

    summary = df.groupby('cls').agg(
        nGenes=('globHOG', 'size'),
        nOGs=('globHOG', 'nunique')
    ).reset_index()

    summary['percGenes'] = (summary['nGenes'] / total_genes * 100).round(1)
    summary['percOGs'] = (summary['nOGs'] / total_ogs * 100).round(1)

    summary = summary[['cls', 'nGenes', 'nOGs', 'percGenes', 'percOGs']]

    category_order = ['core', 'softcore', 'shell', 'cloud', 'private']
    summary['cls'] = pd.Categorical(summary['cls'], categories=category_order, ordered=True)
    summary = summary.sort_values('cls').reset_index(drop=True)

    return summary

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: computePAV.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    result = classify_pangenome(input_file)
    print(result.to_string(index=False))
