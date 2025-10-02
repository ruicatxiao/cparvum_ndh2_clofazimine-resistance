#!/usr/bin/env bash

REGION_FILE="goi.bed"
BAM_DIR="bam"

# Extract region names from the region file (column 4)
regions=($(awk '{print $4}' $REGION_FILE))

# Print header: "bamName" followed by all region names separated by tabs
echo -en "bamName"
for region in "${regions[@]}"; do
    echo -en "\t${region}"
done
echo

# Loop through all BAM files
for bam in $(find $BAM_DIR -name '*.bam'); do
    # Extract sample name from bam filename
    sample=$(basename "$bam" .bam)
    
    # Initialize an array to hold coverage values for this sample
    coverages=()
    
    # For each region, run samtools depth and compute average coverage
    while read -r chr start end rname; do
        # Use samtools depth for the specific region
        # If the region might have no coverage, handle gracefully with default "0"
        avg_cov=$(samtools depth -r "${chr}:${start}-${end}" "$bam" | \
                  awk 'BEGIN{sum=0;count=0} {sum+=$3; count++} END{if (count>0) {print sum/count} else {print 0}}')
        
        coverages+=("$avg_cov")
    done < "$REGION_FILE"
    
    # Print the line: sample name followed by coverage values for each region
    echo -en "${sample}"
    for cov in "${coverages[@]}"; do
        echo -en "\t${cov}"
    done
    echo
done