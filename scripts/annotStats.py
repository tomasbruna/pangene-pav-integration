#!/usr/bin/env python3

import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


def parse_gff3(gff_file):
    cds_by_transcript = defaultdict(list)
    transcript_to_gene = {}
    all_genes = set()
    with open(gff_file) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            chrom, src, ftype, start, end, score, strand, phase, attrs = line.strip().split("\t")
            start, end = int(start), int(end)
            attr_dict = {}
            for part in attrs.split(";"):
                if "=" in part:
                    k, v = part.split("=", 1)
                    attr_dict[k] = v
            if ftype == "gene":
                gid = attr_dict.get("ID")
                if gid:
                    all_genes.add(gid)
            elif ftype in ("mRNA", "transcript"):
                tid = attr_dict.get("ID")
                gid = attr_dict.get("Parent")
                if tid and gid:
                    transcript_to_gene[tid] = gid
                    all_genes.add(gid)
            elif ftype == "CDS":
                parents = attr_dict.get("Parent", "").split(",")
                for tid in parents:
                    cds_by_transcript[tid].append((chrom, start, end, strand))
    return cds_by_transcript, transcript_to_gene, all_genes


def get_splice_sites(seq, intron_start, intron_end, strand):
    """Return donor and acceptor dinucleotides in transcript orientation."""
    if strand == "+":
        donor = seq[intron_start - 1:intron_start + 1].upper()
        acceptor = seq[intron_end - 2:intron_end].upper()
    else:
        donor = str(Seq(seq[intron_end - 2:intron_end]).reverse_complement()).upper()
        acceptor = str(Seq(seq[intron_start - 1:intron_start + 1]).reverse_complement()).upper()
    return donor, acceptor


def collect_gene_stats(gff_file, fasta_file, exon_thresh_small, exon_thresh_ultra, intron_thresh):
    cds_by_transcript, transcript_to_gene, all_genes = parse_gff3(gff_file)
    genome = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}

    # Gene sets
    genes_with_small_exons = set()
    genes_with_ultra_exons = set()
    genes_with_short_introns = set()
    genes_with_gtag = set()
    genes_with_gcag = set()
    genes_with_noncanonical = set()

    # Absolute counts
    small_exon_count = 0
    ultra_exon_count = 0
    short_intron_count = 0
    noncanonical_count = 0

    for tid, cds_list in cds_by_transcript.items():
        if tid not in transcript_to_gene:
            continue
        gid = transcript_to_gene[tid]
        strand = cds_list[0][3]

        # Sort CDS features
        if strand == "+":
            cds_list.sort(key=lambda x: x[1])
        else:
            cds_list.sort(key=lambda x: x[1], reverse=True)

        # Internal exon checks
        if len(cds_list) > 2:
            for chrom, start, end, strand in cds_list[1:-1]:
                exon_len = end - start + 1
                if exon_len < exon_thresh_small:
                    genes_with_small_exons.add(gid)
                    small_exon_count += 1
                if exon_len < exon_thresh_ultra:
                    genes_with_ultra_exons.add(gid)
                    ultra_exon_count += 1

        # Intron checks
        for i in range(len(cds_list) - 1):
            c1 = cds_list[i]
            c2 = cds_list[i + 1]
            chrom = c1[0]
            strand = c1[3]

            if strand == "+":
                intron_start = c1[2] + 1
                intron_end = c2[1] - 1
            else:
                intron_start = c2[2] + 1
                intron_end = c1[1] - 1

            intron_len = intron_end - intron_start + 1
            if intron_len < intron_thresh:
                genes_with_short_introns.add(gid)
                short_intron_count += 1

            if intron_len >= 2:
                donor, acceptor = get_splice_sites(genome[chrom], intron_start, intron_end, strand)
                if donor == "GT" and acceptor == "AG":
                    genes_with_gtag.add(gid)
                elif donor == "GC" and acceptor == "AG":
                    genes_with_gcag.add(gid)
                else:
                    genes_with_noncanonical.add(gid)
                    noncanonical_count += 1

    return {
        "all_genes": all_genes,
        "genes_with_small_exons": genes_with_small_exons,
        "genes_with_ultra_exons": genes_with_ultra_exons,
        "genes_with_short_introns": genes_with_short_introns,
        "genes_with_gtag": genes_with_gtag,
        "genes_with_gcag": genes_with_gcag,
        "genes_with_noncanonical": genes_with_noncanonical,
        "small_exon_count": small_exon_count,
        "ultra_exon_count": ultra_exon_count,
        "short_intron_count": short_intron_count,
        "noncanonical_count": noncanonical_count
    }


def main():
    parser = argparse.ArgumentParser(description="Collect gene statistics from GFF3 and FASTA.")
    parser.add_argument("gff3", help="GFF3 annotation file")
    parser.add_argument("fasta", help="Genome FASTA file")
    parser.add_argument("--exon-threshold-small", type=int, default=10,
                        help="Length threshold in bp for 'small' internal exons (default: 10)")
    parser.add_argument("--exon-threshold-ultra", type=int, default=3,
                        help="Length threshold in bp for 'ultra-tiny' internal exons (default: 3)")
    parser.add_argument("--intron-threshold", type=int, default=10,
                        help="Length threshold in bp for short coding introns")
    parser.add_argument("--list-genes", action="store_true",
                        help="Also print the list of gene IDs for each category")
    args = parser.parse_args()

    stats = collect_gene_stats(
        args.gff3, args.fasta,
        args.exon_threshold_small,
        args.exon_threshold_ultra,
        args.intron_threshold
    )

    total_genes = len(stats['all_genes'])

    # Percentages of genes with each feature
    small_exon_pct = (len(stats['genes_with_small_exons']) / total_genes * 100) if total_genes else 0
    ultra_exon_pct = (len(stats['genes_with_ultra_exons']) / total_genes * 100) if total_genes else 0
    short_intron_pct = (len(stats['genes_with_short_introns']) / total_genes * 100) if total_genes else 0
    gtag_pct = (len(stats['genes_with_gtag']) / total_genes * 100) if total_genes else 0
    gcag_pct = (len(stats['genes_with_gcag']) / total_genes * 100) if total_genes else 0
    noncanonical_pct = (len(stats['genes_with_noncanonical']) / total_genes * 100) if total_genes else 0

    print(f"Genes with small internal exons (<{args.exon_threshold_small} bp): {len(stats['genes_with_small_exons'])} ({small_exon_pct:.2f}%)")
    print(f"Absolute number of small internal exons: {stats['small_exon_count']}")
    if args.list_genes:
        print("\n".join(sorted(stats['genes_with_small_exons'])))

    print(f"Genes with ultra-tiny internal exons (<{args.exon_threshold_ultra} bp): {len(stats['genes_with_ultra_exons'])} ({ultra_exon_pct:.2f}%)")
    print(f"Absolute number of ultra-tiny internal exons: {stats['ultra_exon_count']}")
    if args.list_genes:
        print("\n".join(sorted(stats['genes_with_ultra_exons'])))

    print(f"Genes with short coding introns (<{args.intron_threshold} bp): {len(stats['genes_with_short_introns'])} ({short_intron_pct:.2f}%)")
    print(f"Absolute number of short coding introns: {stats['short_intron_count']}")
    if args.list_genes:
        print("\n".join(sorted(stats['genes_with_short_introns'])))

    print(f"Genes with GT-AG splice sites: {len(stats['genes_with_gtag'])} ({gtag_pct:.2f}%)")
    if args.list_genes:
        print("\n".join(sorted(stats['genes_with_gtag'])))

    print(f"Genes with GC-AG splice sites: {len(stats['genes_with_gcag'])} ({gcag_pct:.2f}%)")
    if args.list_genes:
        print("\n".join(sorted(stats['genes_with_gcag'])))

    print(f"Genes with non-canonical splice sites: {len(stats['genes_with_noncanonical'])} ({noncanonical_pct:.2f}%)")
    if args.list_genes:
        print("\n".join(sorted(stats['genes_with_noncanonical'])))
    print(f"Absolute number of introns with noncanonical splice sites: {stats['noncanonical_count']}")

    # Summary: genes with any problem
    problem_genes = (
        stats['genes_with_small_exons'] |
        stats['genes_with_ultra_exons'] |
        stats['genes_with_short_introns'] |
        stats['genes_with_noncanonical']
    )
    total_genes = len(stats['all_genes'])
    pct_problem = (len(problem_genes) / total_genes * 100) if total_genes > 0 else 0
    print(f"\nSummary: {len(problem_genes)} genes ({pct_problem:.2f}%) "
          f"have at least one problem (short exon/intron or non-canonical splice site) "
          f"out of {total_genes} total genes.")


if __name__ == "__main__":
    main()
