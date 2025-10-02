#!/bin/bash


SOURCE="/data/ruicatxiao/snp_test/vcf"
OUTPUT="/data/ruicatxiao/snp_test/vcf/anno"
SNPEFF="/data/ruicatxiao/snp_test/snpEff"

mkdir -p "$OUTPUT"

find "$SOURCE" -name '*_final.vcf' | while read -r finalvcf; do
    echo "Processing: $finalvcf"
    

    vcfname=$(basename "$finalvcf" _final.vcf)
    echo "Base name: $vcfname"

    # Run the snpEff 
    java -Xmx8g -jar "${SNPEFF}/snpEff.jar" -v \
        CpBGFT2T -csvStats "${OUTPUT}/${vcfname}_stats.csv" \
        -no-downstream -no-upstream \
        "$finalvcf" > "${OUTPUT}/${vcfname}_anno.vcf"
done