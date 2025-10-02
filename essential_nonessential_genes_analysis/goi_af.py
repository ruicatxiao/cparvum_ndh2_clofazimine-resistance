#!/usr/bin/env python3
import os
import sys
import glob

def print_help():
    print("Usage: python3 goi_af.py <regions_file> <vcf_path_or_pattern>")
    print("Examples:")
    print("  python3 goi_af.py goi.bed /vcf/anno/")
    print("  python3 goi_af.py goi.bed vcf/anno/")
    print("  python3 goi_af.py goi.bed 'vcf/anno/*.vcf'")
    print("")
    print("This script generates a matrix of mean allele frequencies of HIGH-impact variants")
    print("in specified genomic regions. INDEL variants use IMF directly. SNP variants compute")
    print("allele frequency from AD/DP fields for the HIGH-impact allele(s).")
    sys.exit(1)

if len(sys.argv) < 3:
    print_help()

REGION_FILE = sys.argv[1]
VCF_INPUT = sys.argv[2]

if not os.path.exists(REGION_FILE):
    print(f"Error: Regions file {REGION_FILE} does not exist.")
    sys.exit(1)

# Determine output file name
base_name = os.path.splitext(REGION_FILE)[0]
output_file = f"{base_name}_af.tsv"

# If VCF_INPUT is a directory, append *_anno.vcf to form the pattern
if os.path.isdir(VCF_INPUT):
    vcf_pattern = os.path.join(VCF_INPUT, "*_anno.vcf")
else:
    # Treat VCF_INPUT as a glob pattern directly
    vcf_pattern = VCF_INPUT

# Parse regions
regions = []
with open(REGION_FILE) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        chr_, start, end, rname = line.split()
        start = int(start)
        end = int(end)
        regions.append((chr_, start, end, rname))

def parse_info(info_str):
    info_dict = {}
    fields = info_str.split(";")
    for field in fields:
        if "=" in field:
            key, val = field.split("=", 1)
            info_dict[key] = val
        else:
            # flags like "INDEL"
            info_dict[field] = True
    return info_dict, fields

def parse_ann_field(info_dict):
    """
    Parse ANN field and return list of (allele, impact) tuples.
    ANN=alt_allele|effect|impact|...
    """
    annotations = []
    if "ANN" not in info_dict:
        return annotations
    ann_value = info_dict["ANN"]
    ann_entries = ann_value.split(",")
    for ann in ann_entries:
        parts = ann.split("|")
        if len(parts) > 2:
            allele = parts[0]
            impact = parts[2]
            annotations.append((allele, impact))
    return annotations

def compute_variant_frequency(record):
    # record: [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, ...]
    chrom, pos, vid, ref, alt, qual, filt, info = record[:8]
    info_dict, info_fields = parse_info(info)
    alt_alleles = alt.split(",")

    # Check if variant has HIGH impact
    annotations = parse_ann_field(info_dict)
    high_alt_indices = []
    for i, a in enumerate(alt_alleles):
        # Check if this alt allele is HIGH impact
        allele_high = any((ann_allele == a and imp == "HIGH") for (ann_allele, imp) in annotations)
        if allele_high:
            high_alt_indices.append(i)

    if len(high_alt_indices) == 0:
        # No HIGH impact alleles at this site
        return []

    # Determine if INDEL or SNP:
    # If the first field in info_fields is "INDEL", it's an INDEL variant.
    first_field = info_fields[0]
    is_indel = (first_field == "INDEL")

    if is_indel:
        # INDEL variant
        # Just use IMF if present
        if "IMF" in info_dict:
            imf_val = float(info_dict["IMF"])
            # Return one frequency for the variant (IMF)
            return [imf_val]
        else:
            # No IMF? Then skip (no data)
            return []
    else:
        # SNP variant
        if "DP" not in info_dict or "AD" not in info_dict:
            return []
        dp = int(info_dict["DP"])
        if dp == 0:
            return []
        ad_vals = [int(x) for x in info_dict["AD"].split(",")]
        freqs = []
        for i in high_alt_indices:
            if i+1 < len(ad_vals):
                alt_depth = ad_vals[i+1]
                freq = alt_depth / dp
                freqs.append(freq)
            else:
                freqs.append(0.0)

        if len(freqs) == 0:
            return []
        # Average freq for multiple HIGH alleles if any
        return [sum(freqs)/len(freqs)]

vcf_files = glob.glob(vcf_pattern)
vcf_files.sort()

if len(vcf_files) == 0:
    print(f"No VCF files found with pattern: {vcf_pattern}")
    sys.exit(1)

with open(output_file, "w") as out_f:
    # Print header
    out_f.write("vcfName")
    for r in regions:
        out_f.write("\t" + r[3])
    out_f.write("\n")

    total_files = len(vcf_files)
    for idx, vcf in enumerate(vcf_files, start=1):
        if os.path.isdir(vcf):
            # Skip directories if any match
            continue
        sample = os.path.basename(vcf).replace("_anno.vcf", "").replace(".vcf", "")

        # Print progress
        print(f"Processing file {idx}/{total_files}: {vcf}", file=sys.stderr)

        # Read all variants
        variants = []
        with open(vcf) as vf:
            for line in vf:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                variants.append(fields)

        # For each region, compute freq
        region_freqs = []
        for (chr_, start, end, rname) in regions:
            region_variants = [v for v in variants if v[0] == chr_ and start <= int(v[1]) <= end]
            all_variant_freqs = []
            for rv in region_variants:
                vfreqs = compute_variant_frequency(rv)
                all_variant_freqs.extend(vfreqs)

            if len(all_variant_freqs) == 0:
                avg_freq = 0.0
            else:
                avg_freq = sum(all_variant_freqs)/len(all_variant_freqs)

            region_freqs.append(avg_freq)

        out_f.write(sample)
        for af in region_freqs:
            out_f.write("\t" + str(af))
        out_f.write("\n")

print(f"Output written to {output_file}")
