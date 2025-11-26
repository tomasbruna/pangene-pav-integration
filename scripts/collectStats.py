#!/usr/bin/env python3

import re
from tabulate import tabulate
import argparse
import csv
import sys

# grep "Summary:" *.stats | cut -f2 -d \( | cut -f1 -d % | awk 'BEGIN{sum=0;cnt=0}{sum+=$1;cnt+=1}END{print sum/cnt}'

def parse_stats_file(filename):
    data = {}
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            match_small_internal_exons_absolute = re.match(r'Absolute number of small internal exons: (\d+)', line)
            match_small_internal_exons = re.match(r'Genes with small internal exons \(<\d+ bp\): (\d+) \((\d+\.\d+%)\)', line)
            match_ultra_tiny_internal_exons_absolute = re.match(r'Absolute number of ultra-tiny internal exons: (\d+)', line)
            match_ultra_tiny_internal_exons = re.match(r'Genes with ultra-tiny internal exons \(<\d+ bp\): (\d+) \((\d+\.\d+%)\)', line)
            match_non_canonical_splice_sites_absolute = re.match(r'Absolute number of introns with noncanonical splice sites: (\d+)', line)
            match_non_canonical_splice_sites = re.match(r'Genes with non-canonical splice sites: (\d+) \((\d+\.\d+%)\)', line)
            match_genes_with_any_issues = re.match(r'Summary: (\d+) genes \((\d+\.\d+%)\) have at least one problem', line)

            if match_small_internal_exons_absolute:
                if 'Small Internal Exons (<10 bp)' in data:
                    data['Small Internal Exons (<10 bp)']['absolute'] = match_small_internal_exons_absolute.group(1)
                else:
                    data['Small Internal Exons (<10 bp)'] = {'absolute': match_small_internal_exons_absolute.group(1), 'percentage': None}
            if match_small_internal_exons:
                if 'Small Internal Exons (<10 bp)' in data:
                    data['Small Internal Exons (<10 bp)']['percentage'] = match_small_internal_exons.group(2)
                else:
                    data['Small Internal Exons (<10 bp)'] = {'absolute': None, 'percentage': match_small_internal_exons.group(2)}

            if match_ultra_tiny_internal_exons_absolute:
                if 'Ultra-Tiny Internal Exons (<3 bp)' in data:
                    data['Ultra-Tiny Internal Exons (<3 bp)']['absolute'] = match_ultra_tiny_internal_exons_absolute.group(1)
                else:
                    data['Ultra-Tiny Internal Exons (<3 bp)'] = {'absolute': match_ultra_tiny_internal_exons_absolute.group(1), 'percentage': None}
            if match_ultra_tiny_internal_exons:
                if 'Ultra-Tiny Internal Exons (<3 bp)' in data:
                    data['Ultra-Tiny Internal Exons (<3 bp)']['percentage'] = match_ultra_tiny_internal_exons.group(2)
                else:
                    data['Ultra-Tiny Internal Exons (<3 bp)'] = {'absolute': None, 'percentage': match_ultra_tiny_internal_exons.group(2)}

            if match_non_canonical_splice_sites_absolute:
                if 'Non-canonical introns' in data:
                    data['Non-canonical introns']['absolute'] = match_non_canonical_splice_sites_absolute.group(1)
                else:
                    data['Non-canonical introns'] = {'absolute': match_non_canonical_splice_sites_absolute.group(1), 'percentage': None}
            if match_non_canonical_splice_sites:
                if 'Non-canonical introns' in data:
                    data['Non-canonical introns']['percentage'] = match_non_canonical_splice_sites.group(2)
                else:
                    data['Non-canonical introns'] = {'absolute': None, 'percentage': match_non_canonical_splice_sites.group(2)}

            if match_genes_with_any_issues:
                data['Genes with any Issues'] = f"{match_genes_with_any_issues.group(1)} ({match_genes_with_any_issues.group(2)})"

    return data

def main():
    parser = argparse.ArgumentParser(description='Parse stats files and print in a table')
    parser.add_argument('-o', '--order_file', help='File with list of filenames in a specific order')
    args = parser.parse_args()

    if args.order_file:
        try:
            with open(args.order_file, 'r') as file:
                filenames = [line.strip() for line in file.readlines()]
        except FileNotFoundError:
            print(f"Error: File '{args.order_file}' not found.")
            return

    table_data = []

    for filename in filenames:
        full_filename = f"{filename}.stats"
        if full_filename.endswith('.stats'):
            try:
                data = parse_stats_file(full_filename)
                row = [
                    f"{data.get('Small Internal Exons (<10 bp)', {'absolute': 'Not Found', 'percentage': 'Not Found'}).get('absolute', 'Not Found')} ({data.get('Small Internal Exons (<10 bp)', {'absolute': 'Not Found', 'percentage': 'Not Found'}).get('percentage', 'Not Found')})",
                    f"{data.get('Ultra-Tiny Internal Exons (<3 bp)', {'absolute': 'Not Found', 'percentage': 'Not Found'}).get('absolute', 'Not Found')} ({data.get('Ultra-Tiny Internal Exons (<3 bp)', {'absolute': 'Not Found', 'percentage': 'Not Found'}).get('percentage', 'Not Found')})",
                    f"{data.get('Non-canonical introns', {'absolute': 'Not Found', 'percentage': 'Not Found'}).get('absolute', 'Not Found')} ({data.get('Non-canonical introns', {'absolute': 'Not Found', 'percentage': 'Not Found'}).get('percentage', 'Not Found')})",
                    data.get('Genes with any Issues', 'Not Found')
                ]
                table_data.append([full_filename] + row)
            except FileNotFoundError:
                print(f"Warning: File '{full_filename}' not found.")
        else:
            print(f"Warning: Skipping '{filename}' as it does not end with '.stats'")

    headers = ['Filename', 'Small Internal Exons (<10 bp)', 'Ultra-Tiny Internal Exons (<3 bp)', 'Non-canonical introns', 'Genes with any Issues']

    with sys.stdout as csvfile:
        writer = csv.writer(csvfile, delimiter=';')
        writer.writerow(headers)
        for row in table_data:
            new_row = []
            for cell in row:
                match = re.match(r'(\d+) \((\d+\.\d+%)\)', cell)
                if match:
                    absolute, percentage = match.groups()
                    new_row.append(f"{int(absolute):,} ({float(percentage[:-1]):.1f}%)")
                else:
                    new_row.append(str(cell).replace('\n', ' '))
            writer.writerow(new_row)

if __name__ == "__main__":
    main()
