#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
#
# Generate heatmaps from pre-computed Liftoff mappings to assess annotation consistency.
# Analyzes how many genes that map perfectly between genomes are also predicted identically
# and tracks gene presence/absence variation (PAV) across all genome pairs.
# ==============================================================


import argparse
import csv
import re
import sys
import os
import subprocess
import glob
import shutil
import tempfile
import random
import string
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from collections import Counter


class Pair:
    def __init__(self, source, target, sourceFile, targetFile, liftoff):
        self.source = source
        self.target = target
        self.sourceFile = sourceFile
        self.targetFile = targetFile
        self.liftoff = liftoff
        self.comparison = None
        self.gene_identity = 0
        self.gene_presence = 0
        self.lociNum = 0

    def _select(self, mappingThreshold, usemRNA=False):
        selected = False
        outF = tempfile.NamedTemporaryFile(mode='w', delete=False,
                                           suffix='.gff3')
        for row in csv.reader(open(self.liftoff), delimiter='\t'):
            if len(row) != 9:
                continue
            overallFeature = "gene"
            if usemRNA:
                overallFeature = "mRNA"
            if row[2] == overallFeature:
                selected = False
                coverage = float(extractFeatureGff(row[8], "coverage"))
                sequence_ID = float(extractFeatureGff(row[8], "sequence_ID"))
                valid_ORFs = extractFeatureGff(row[8], "valid_ORFs")

                if (valid_ORFs == "1" and coverage >= mappingThreshold and
                sequence_ID >= mappingThreshold):
                    selected = True

            if selected:
                outF.write("\t".join(row) + "\n")
        outF.close()
        selectedLiftoff = outF.name
        return selectedLiftoff

    def makeComparison(self, liftoffFolder, mappingThreshold, usemRNA=False):
        outFname = (f"{liftoffFolder}/{os.path.basename(self.liftoff)}_"
                    f"{mappingThreshold}.gffcmp.tracking")
        if os.path.exists(outFname):
            self.comparison = outFname
            return

        selectedLiftoff = self._select(mappingThreshold, usemRNA)

        prefix = ''.join(random.choice(string.ascii_letters) for _ in range(5))
        systemCall(f"gffcompare --strict-match -e 3 -D -T -o {prefix} -r {self.targetFile} "
                   f"{selectedLiftoff}")

        shutil.move(f"{prefix}.tracking", outFname)
        self.comparison = outFname

        # Cleanup
        for file_path in glob.glob(f"{prefix}*"):
            os.remove(file_path)
        os.unlink(selectedLiftoff)

    def parseComparisonsFromStats(self):
        with open(self.comparison) as file:
            content = file.read()
        gene_identity_match = re.search(r'Locus level:\s+(\d+\.\d+)', content)
        self.gene_identity = float(gene_identity_match.group(1))

        missedLoci_match = re.search(r'Missed loci:\s+(\d+)/(\d+)', content)
        missed = float(missedLoci_match.group(1))
        total = float(missedLoci_match.group(2))

        self.lociNum = int(total)
        self.gene_presence = round(100 * (total - missed) / total, 2)


    def parseComparisons(self):
        matchTypes = []
        with open(self.comparison) as file:
            for l in file:
                matchTypes.append(l.split()[3])

        total = len(matchTypes)
        matchCounts = Counter(matchTypes)
        # Count a subset as an exact match (for lifted truncated genes).
        exactMatches = sum(matchCounts[key] for key in
                         ["=", "c"])
        anyOverlap = sum(matchCounts[key] for key in
                         ["=", "c", "k", "m", "n", "j", "e", "o", "~"])

        self.lociNum = total
        self.gene_identity = round(100 * exactMatches / total, 2)
        self.gene_presence = round(100 * anyOverlap / total, 2)


def read_order_file(order_file):
    """Read the order file and return a list of items in the specified order."""
    order = []
    with open(order_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                order.append(line)
    return order


def reorder_by_file(items, order_file):
    """Reorder items according to the order specified in the file."""
    if not order_file or not os.path.exists(order_file):
        return sorted(items)
    
    file_order = read_order_file(order_file)
    
    order_map = {item: i for i, item in enumerate(file_order)}
    
    # Sort items: first by whether they're in the order file (and their position),
    # then alphabetically for items not in the file
    ordered_items = []
    unordered_items = []
    
    for item in items:
        if item in order_map:
            ordered_items.append((order_map[item], item))
        else:
            unordered_items.append(item)
    
    ordered_items.sort(key=lambda x: x[0])
    
    result = [item for _, item in ordered_items] + sorted(unordered_items)
    
    return result


def create_heatmap(pairs, color_metric, mappingThreshold,
                              count_metric='lociNum',
                              title="", figsize=(10.2, 8.5), cmap="magma_r",
                              noCellText=False,
                              save_csv=False,
                              csv_dir=None,
                              order_file=None,
                              vmin=30,
                              vmax=100):

    color_data = {}
    count_data = {}

    for pair in pairs:
        try:
            color_value = float(getattr(pair, color_metric))
            count_value = int(getattr(pair, count_metric))
            color_data[(pair.source, pair.target)] = color_value
            count_data[(pair.source, pair.target)] = count_value
        except (ValueError, TypeError, AttributeError) as e:
            print(f"Warning: Could not process metrics for {pair.source}-{pair.target}: {e}")
            continue

    sources = list(set(pair.source for pair in pairs))
    targets = list(set(pair.target for pair in pairs))
    
    sources = reorder_by_file(sources, order_file)
    targets = reorder_by_file(targets, order_file)

    color_matrix = []
    count_matrix = []

    for source in sources:
        color_row = []
        count_row = []
        for target in targets:
            color_row.append(color_data.get((source, target), np.nan))
            count_row.append(count_data.get((source, target), np.nan))
        color_matrix.append(color_row)
        count_matrix.append(count_row)

    color_array = np.array(color_matrix, dtype=float)
    count_array = np.array(count_matrix, dtype=object)

    color_df = pd.DataFrame(color_array, index=sources, columns=targets)

    # Fill diagonal with 100
    for i in range(len(sources)):
        if i < len(targets) and sources[i] == targets[i]:
            color_df.iloc[i, i] = 100
            count_array[i, i] = 0  # or any appropriate value

    if save_csv:
        csv_filename = f'{color_metric}_{mappingThreshold}.csv'
        if csv_dir:
            csv_filename = os.path.join(csv_dir, csv_filename)
        color_df.fillna(100).to_csv(csv_filename)
        print(f"Saved matrix to {csv_filename}")

    plt.figure(figsize=figsize)

    mask = np.isnan(color_df.values)

    annot = np.empty_like(color_array, dtype=object)
    for i in range(len(sources)):
        for j in range(len(targets)):
            if i == j:
                annot[i, j] = ""   # No cell text on diagonal
            elif not np.isnan(color_array[i, j]) and not pd.isna(count_array[i, j]):
                if not noCellText:
                    annot[i, j] = f"{color_array[i, j]:.1f}\n(n={int(count_array[i, j])})"
                else:
                    annot[i, j] = ""
            else:
                annot[i, j] = ""

    ax = sns.heatmap(color_df,
                    annot=annot,
                    fmt="",
                    cmap=cmap,
                    linewidths=0.5,
                    mask=mask,
                    vmin=vmin,
                    vmax=vmax,
                    cbar_kws={'label': '%'})

    plt.title(title)
    plt.ylabel('Source')
    plt.xlabel('Target')

    plt.xticks(rotation=90)

    plt.tight_layout()

    plt.savefig(f'{color_metric}_{mappingThreshold}_{cmap}.pdf', dpi=600, bbox_inches='tight')


def extractFeatureGff(text, feature):
    regex = ";" + feature + '=([^;]+)'
    result = re.search(regex, text)
    if result:
        return result.groups()[0]
    else:
        return None


def systemCall(cmd):
    if subprocess.call(["bash", "-c", cmd]) != 0:
        sys.exit('error: Program exited due to an ' +
                 'error in command: ' + cmd)


def extractIdFromFile(filename):
    return os.path.basename(filename).split(".")[0]


def buildPairs(annotationsFolder, liftoffFolder):
    pairs = []
    gff3_files = glob.glob(f"{annotationsFolder}/*.gff3")
    for file1 in gff3_files:
        for file2 in gff3_files:
            if file1 == file2:
                continue
            name1 = extractIdFromFile(file1)
            name2 = extractIdFromFile(file2)
            print(name1, name2)
            liftoffs = glob.glob(f"{liftoffFolder}/*{name1}.*_to_*{name2}.*.gff3")
            if len(liftoffs) != 1:
                sys.exit(f"error: too many or no pair matches between "
                         f"{name1} and {name2}: {liftoffs}")
            liftoff = liftoffs[0]
            pairs.append(Pair(name1, name2, file1, file2, liftoff))
    return pairs


def makeComparisons(pairs, liftoffFolder, mappingThreshold, usemRNA=False):
    for pair in pairs:
        pair.makeComparison(liftoffFolder, mappingThreshold, usemRNA)
        pair.parseComparisons()


def main():
    args = parseCmd()
    if not os.path.exists(args.comparisonFolder):
        os.makedirs(args.comparisonFolder)
    pairs = buildPairs(args.annotationsFolder, args.liftoffFolder)
    makeComparisons(pairs, args.comparisonFolder, args.mappingThreshold, args.usemRNA)

    create_heatmap(pairs, 'gene_identity', args.mappingThreshold,
                           title="Percentage of perfectly mapped genes predicted identically",
                           noCellText=args.noCellText,
                           save_csv=args.save_csv,
                           csv_dir=args.csv_dir,
                           order_file=args.order_file,
                           vmin=args.vmin,
                           vmax=args.vmax,
                           cmap=args.palette)
    create_heatmap(pairs, 'gene_presence', args.mappingThreshold,
                           title="Percentage of perfectly mapped genes present",
                           noCellText=args.noCellText,
                           save_csv=args.save_csv,
                           csv_dir=args.csv_dir,
                           order_file=args.order_file,
                           vmin=args.vmin,
                           vmax=args.vmax,
                           cmap=args.palette)


def parseCmd():

    parser = argparse.ArgumentParser(
        description='Generate heatmaps from pre-computed Liftoff mappings to assess annotation consistency. '
                    'Analyzes how many genes that map perfectly between genomes are also predicted identically, '
                    'and tracks gene presence/absence variation (PAV) across all genome pairs.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('annotationsFolder', type=str)
    parser.add_argument('liftoffFolder', type=str)
    parser.add_argument('comparisonFolder', type=str)
    parser.add_argument('--mappingThreshold', type=float, default=1)
    parser.add_argument('--usemRNA',  default=False, action='store_true')
    parser.add_argument('--noCellText',  default=False, action='store_true')
    parser.add_argument('--save-csv', default=False, action='store_true',
                        help='Save the resulting matrices to CSV files')
    parser.add_argument('--csv-dir', type=str, default=None,
                        help='Directory to save CSV files (default: current directory)')
    parser.add_argument('--order-file', type=str, default=None,
                        help='File containing the desired order for matrix rows/columns (one item per line)')
    parser.add_argument('--vmin', type=float, default=30,
                        help='Minimum value for color scale')
    parser.add_argument('--vmax', type=float, default=100,
                        help='Maximum value for color scale')
    parser.add_argument('--palette', type=str, default="magma_r",
                        help='Color palette.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
