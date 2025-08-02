#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
#
# Estimate pangenome growth with cd-hit from a single clustering
# run of all proteins.
# ==============================================================


import argparse
import csv
import re
import sys
import subprocess
import os
import shutil
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
import tempfile


def systemCall(cmd):
    if subprocess.call(["bash", "-c", cmd]) != 0:
        sys.exit('error: Program exited due to an ' +
                 'error in command: ' + cmd)


def combine_and_assign_unique_names(input_folder):
    temp_file = tempfile.mkstemp(suffix='.fasta')[1]

    with open(temp_file, 'w') as out_handle:
        for filename in os.listdir(input_folder):
            if filename.endswith(('.fasta', '.fa', '.fna')):
                input_path = os.path.join(input_folder, filename)
                base_name = os.path.splitext(filename)[0]

                sequences = SeqIO.parse(input_path, 'fasta')
                for i, record in enumerate(sequences, 1):
                    record.id = f"{base_name}_{i}"
                    record.description = ""
                    original_seq = str(record.seq)
                    if "." in original_seq:
                        print(f"Warning: Dots found in sequence {record.id} and will be removed")
                        record.seq = Seq(original_seq.replace(".", ""))
                    SeqIO.write(record, out_handle, 'fasta')

    return temp_file


def isDNA(sequence):
    sequence = sequence.upper().replace(' ', '')
    nucleotides = set('ATCGN')
    nuc_count = sum(1 for char in sequence if char in nucleotides)
    nuc_percentage = (nuc_count / len(sequence)) * 100
    return True if nuc_percentage > 85 else False


def run_diamond(combinedFile, identity, coverage, threads):
    outf = tempfile.mkstemp(suffix='.clstr')[1]
    database = tempfile.mkstemp(suffix='.dmnd')[1]


    cmd=(f"diamond makedb -d {database} --in {combinedFile} --threads {threads}")
    sys.stderr.write(f'runnig: {cmd}\n')
    systemCall(cmd)

    with open(combinedFile) as f:
        record = str(next(SeqIO.parse(f, 'fasta')).seq)
        if isDNA(record):
            sys.exit("Detected a DNA sequence. DIAMOND only works w/ prots\n")

    cmd=(f"diamond cluster -d {database} -o {outf} --approx-id {identity} --member-cover {coverage} "
         f"--threads {threads}")
    sys.stderr.write(f'runnig: {cmd}\n')
    systemCall(cmd)

    os.remove(database)
    return outf


def parse_cluster_file(clusters):
    line_clusters = {}
    current_cluster = None
    pav_counts = Counter()
    lines = set()
    prev_cluster = ""

    with open(clusters, 'r') as f:
        for line in f:
            current_cluster = line.strip().split()[0]

            if prev_cluster != current_cluster:
                if len(lines) != 0:
                    pav_counts[len(lines)] += 1
                lines = set()

            line_name = ".".join(line.strip().split()[1].split('_')[:-1])
            lines.add(line_name)

            if line_name not in line_clusters:
                line_clusters[line_name] = []
            line_clusters[line_name].append(current_cluster)
            prev_cluster = current_cluster

    return line_clusters, pav_counts

def estimate_from_clusters(clusters, outfh):
    parsed_clusters, pav_counts = parse_cluster_file(clusters)
    seenClusters = set()
    for name, clusterList in sorted(parsed_clusters.items()):
        for c in clusterList:
            seenClusters.add(c)
        outfh.write(f"{name}\t{len(seenClusters)}\n")
    for i in range(1, len(pav_counts) + 1):
        outfh.write(f'{i}: {pav_counts[i]}\t')
    outfh.write('\n')


def main():
    args = parseCmd()
    combinedFile = combine_and_assign_unique_names(args.input_folder)
    clusters = run_diamond(combinedFile,
                           args.identity, args.coverage,
                           args.threads)
    outfh = open(args.outputFile, "w")
    outfh.write(f'#identity={args.identity}\n')
    print(clusters)
    estimate_from_clusters(clusters, outfh)
    outfh.close()
    os.remove(combinedFile)
    os.remove(clusters)

def parseCmd():

    parser = argparse.ArgumentParser(description=' Estimate pangenome growth \
                                    with cd-hit from a single clustering \
                                    run of all proteins',
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_folder', type=str,
                        help='Input folder with proteins or CDS seqs.')

    parser.add_argument('outputFile', type=str,
                        help='OutFile')

    parser.add_argument('--identity', type=int, default=80)
    parser.add_argument('--coverage', type=int, default=90)
    parser.add_argument('--threads', type=int, default=8)

    return parser.parse_args()


if __name__ == '__main__':
    main()
