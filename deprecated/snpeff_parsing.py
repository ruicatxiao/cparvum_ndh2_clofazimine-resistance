#!/usr/bin/env python3

import argparse
import glob
import os
import sys

def parse_vcf(file_path, gene_name):
    """
    Parses a VCF file to count the number of impacts per category for a specific gene,
    excluding upstream and downstream gene variants.

    Parameters:
        file_path (str): Path to the VCF file.
        gene_name (str): The gene name to filter annotations.

    Returns:
        dict: A dictionary with counts of each impact category.
    """
    counts = {'HIGH': 0, 'MODERATE': 0, 'LOW': 0, 'MODIFIER': 0}

    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue  # Skip header lines
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue  # Ensure there are enough fields

                info_field = fields[7]  # INFO is the 8th column

                # Extract the ANN field (snpEff annotation)
                ann_list = []
                for info in info_field.split(';'):
                    if info.startswith('ANN='):
                        ann_value = info[4:]  # Remove 'ANN=' prefix
                        ann_list = ann_value.split(',')  # Split multiple annotations
                        break

                # If no ANN annotations found, skip this variant
                if not ann_list:
                    continue

                # Iterate over each annotation
                for ann in ann_list:
                    ann_fields = ann.split('|')
                    if len(ann_fields) < 4:
                        continue  # Ensure there are enough annotation fields

                    allele = ann_fields[0]
                    annotation = ann_fields[1]
                    impact = ann_fields[2].upper()
                    gene = ann_fields[3]

                    # Exclude annotations that are upstream or downstream variants
                    if gene != gene_name:
                        continue
                    if annotation in {'upstream_gene_variant', 'downstream_gene_variant'}:
                        continue

                    # Increment counts if impact category matches
                    if impact in counts:
                        counts[impact] += 1

    except Exception as e:
        print(f"Error processing file {file_path}: {e}", file=sys.stderr)

    return counts

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(
        description='Parse snpEff annotated VCF files for a specific gene, excluding upstream and downstream variants, and count impact categories.'
    )
    parser.add_argument('gene', help='Gene name to parse annotations for (e.g., cpbgf_7001900).')
    args = parser.parse_args()
    gene_name = args.gene

    # Find all *_ann.vcf files recursively
    vcf_files = glob.glob('**/*_ann.vcf', recursive=True)

    if not vcf_files:
        print("No *_ann.vcf files found in the current directory and subdirectories.")
        sys.exit(1)

    total_files = len(vcf_files)
    print(f"Found {total_files} *_ann.vcf files to process.\n")

    # Initialize a list to store output data
    output_data = []

    for idx, vcf_file in enumerate(vcf_files, start=1):
        # Extract SRA name from the filename
        sra_name = os.path.basename(vcf_file).replace('_ann.vcf', '')

        # Log progress
        print(f"Processing file {idx}/{total_files}: {sra_name}")

        # Parse the VCF file for the specified gene
        counts = parse_vcf(vcf_file, gene_name)

        # Append the results to the output data list
        output_data.append([
            sra_name,
            counts.get('HIGH', 0),
            counts.get('MODERATE', 0),
            counts.get('LOW', 0),
            counts.get('MODIFIER', 0)
        ])

    # Sort the output data by SRA number for consistency
    output_data.sort()

    # Define the output filename
    output_file = f"{gene_name}_parsed.txt"

    # Write the results to the output file
    try:
        with open(output_file, 'w') as out:
            # Write the header
            header = ['SRAnumber', 'HIGH', 'MODERATE', 'LOW', 'MODIFIER']
            out.write('\t'.join(header) + '\n')

            # Write each row of data
            for row in output_data:
                out.write('\t'.join(map(str, row)) + '\n')

        print(f"\nParsing complete. Results written to {output_file}")
    except Exception as e:
        print(f"Error writing to output file {output_file}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
