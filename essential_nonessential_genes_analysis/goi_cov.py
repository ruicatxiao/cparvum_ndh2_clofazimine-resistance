#!/usr/bin/env python3
import sys
import os
import glob
import subprocess

def print_help():
    print("Usage: python3 coverage_calc.py <regions_file> <bam_directory>")
    print("")
    print("This script computes average coverage per region for each BAM file in the specified directory.")
    print("It replicates the logic of the provided bash script using Python.")
    print("")
    print("Example:")
    print("  python3 coverage_calc.py goi.bed bam/")
    sys.exit(1)

if len(sys.argv) < 3:
    print_help()

REGION_FILE = sys.argv[1]
BAM_DIR = sys.argv[2]

if not os.path.exists(REGION_FILE):
    print(f"Error: Region file {REGION_FILE} does not exist.")
    sys.exit(1)

if not os.path.isdir(BAM_DIR):
    print(f"Error: {BAM_DIR} is not a directory.")
    sys.exit(1)

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

# Output file name
base_name = os.path.splitext(REGION_FILE)[0]
output_file = f"{base_name}_cov.tsv"

# Extract region names for header
region_names = [r[3] for r in regions]

# Find all BAM files
bam_files = glob.glob(os.path.join(BAM_DIR, "*.bam"))
bam_files.sort()

if len(bam_files) == 0:
    print(f"No BAM files found in {BAM_DIR}")
    sys.exit(1)

with open(output_file, "w") as out_f:
    # Print header
    out_f.write("bamName")
    for rn in region_names:
        out_f.write("\t" + rn)
    out_f.write("\n")

    # Process each BAM file
    for bam in bam_files:
        sample = os.path.basename(bam).replace(".bam", "")
        print(f"Processing {sample}...")
        coverages = []
        for (chr_, start, end, rname) in regions:
            region_str = f"{chr_}:{start}-{end}"
            # Run samtools depth
            cmd = ["samtools", "depth", "-r", region_str, bam]
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                print(f"Error running samtools depth on {bam} for region {region_str}: {result.stderr}")
                avg_cov = 0.0
            else:
                # Compute average coverage from output
                lines = result.stdout.strip().split("\n")
                sum_cov = 0.0
                count = 0
                for line_dep in lines:
                    if line_dep.strip():
                        parts = line_dep.split()
                        # Format: CHR POS DEPTH
                        if len(parts) == 3:
                            depth = float(parts[2])
                            sum_cov += depth
                            count += 1
                if count > 0:
                    avg_cov = sum_cov / count
                else:
                    avg_cov = 0.0
            
            coverages.append(str(avg_cov))

        # Write line for this sample
        out_f.write(sample + "\t" + "\t".join(coverages) + "\n")

print(f"Coverage calculation complete. Results written to {output_file}")
