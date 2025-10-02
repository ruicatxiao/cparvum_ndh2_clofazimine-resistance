#!/usr/bin/env bash


VCF_DIR="vcf/anno"
REGIONS="goi.bed"

# Extract all region names into an array
regions=($(awk '{print $4}' $REGIONS))

# Print header line
echo -en "vcfName"
for r in "${regions[@]}"; do
    echo -en "\t${r}"
done
echo

# Loop through VCF files
for vcf in $(find $VCF_DIR -name '*_anno.vcf'); do
    sample=$(basename "$vcf" _anno.vcf)
    
    # Create an array to hold mean allele frequencies for each region
    region_freqs=()
    
    # For each region, we will:
    # - Extract variants in that coordinate range
    # - Filter for '|HIGH|' impact variants
    # - Parse allele frequencies
    # - Compute mean
    
    while read -r chr start endPos rname; do
        # Use awk to filter variants within region and with HIGH impact.
        # Steps:
        # 1. Exclude header lines
        # 2. Match chromosome
        # 3. Position between start and end
        # 4. Check for '|HIGH|' in ANN (INFO field)
        
        # We'll do all extraction in a single awk command for efficiency.
        # After extracting the lines, we will parse the INFO field to compute allele frequencies.
        
        # Explanation of the awk parsing logic:
        # - For each variant line, extract INFO field (column 8).
        # - Check if INDEL is present:
        #     If yes, extract IMF=... value as freq.
        # - If not INDEL:
        #     Extract DP and AD fields. AD=ref,alt1,alt2,... 
        #     Sum all alt alleles and divide by DP.
        #
        # We'll store frequencies, then average them after.
        
        freqs=$(awk -v CHR="$chr" -v START="$start" -v ENDPOS="$endPos" '
            BEGIN {FS="\t"; OFS="\t"}
            $0 !~ /^#/ && $1 == CHR && $2 >= START && $2 <= ENDPOS {
                # Check INFO field for HIGH impact
                # ANN field typically present in INFO, something like ANN=...
                # We just need to see if "|HIGH|" occurs anywhere in the INFO.
                if(index($8, "|HIGH|") > 0) {
                    # We must compute the allele frequency.
                    info=$8
                    # Split INFO by semicolon
                    n=split(info, arr, ";")
                    dp=""; ad=""; imf=""
                    indel=0
                    for(i=1;i<=n;i++){
                        if(arr[i] ~ /^INDEL/){indel=1}
                        else if(arr[i] ~ /^IMF=/){imf=substr(arr[i],5)}
                        else if(arr[i] ~ /^DP=/){dp=substr(arr[i],4)}
                        else if(arr[i] ~ /^AD=/){ad=substr(arr[i],4)}
                    }
                    
                    if(indel==1 && imf != "") {
                        # Use IMF for INDEL
                        print imf
                    } else {
                        # Non-INDEL case:
                        # AD typically: "0,138" for single alt or "0,121,1" for multi-alt
                        # sum alt alleles: everything except first value
                        split(ad,adarr,",")
                        alt_sum=0
                        for(j=2; j<=length(adarr); j++){
                            alt_sum+=adarr[j]
                        }
                        if(dp != "" && dp > 0) {
                            freq=alt_sum/dp
                            print freq
                        }
                    }
                }
            }
        ' "$vcf")
        
        # Compute average frequency for the region
        if [ -z "$freqs" ]; then
            # No variants found in this region => frequency = 0
            avg=0
        else
            # Average the frequencies (freqs could have multiple lines)
            avg=$(echo "$freqs" | awk '{sum+=$1;count++} END{if(count>0) print sum/count; else print 0}')
        fi
        
        region_freqs+=("$avg")
        
    done < "$REGIONS"
    
    # Print line: sample followed by frequencies
    echo -en "$sample"
    for af in "${region_freqs[@]}"; do
        echo -en "\t${af}"
    done
    echo
done