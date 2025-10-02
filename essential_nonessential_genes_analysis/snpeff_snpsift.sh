# Building database fpr snpEff with Cparvum Iowa BGF T2T assembly and annotation

java -jar snpEff.jar build -gtf22 -v CpBGFT2T

# Code for running SnpSift

cat vcf/anno_merged_filtered_sorted.vcf \
| java -jar snpEff/SnpSift.jar filter "(ANN[*].IMPACT = 'HIGH')" \
| awk '$1~"^CP14"' | head -2 | awk '{print NF}'
 
cat vcf/anno_merged_filtered_sorted.vcf \
| java -jar snpEff/SnpSift.jar filter "(ANN[*].IMPACT = 'HIGH')" \
| awk '$1~"^CP14"' | head -2 | cut -f 1-9

cat vcf/anno_merged_filtered_sorted.vcf \
| java -jar snpEff/SnpSift.jar filter "(ANN[*].IMPACT = 'HIGH')" \
| awk '$1~"^CP14"' | cut -f 10-89 \
| awk '{
    count = 0
    for (i = 1; i <= NF; i++) {
        if ($i == "./.:.") {
            count++
        }
    }
    print ((80-count)/80)*100
}' > analysis/sampleper_high_impact.txt

cat vcf/anno_merged_filtered_sorted.vcf \
|  awk '$1~"^CP14"' | cut -f 8 | awk -F ';' '{print $1}' \
| awk '$1!="INDEL"' | wc -l

cat vcf/anno_merged_filtered_sorted.vcf \
| java -jar snpEff/SnpSift.jar filter "(ANN[*].IMPACT = 'HIGH')" \
| awk '$1~"^CP14"' |  cut -f 8 | awk -F ';' '{print $1}' \
| awk '{if($1=="INDEL"){print $1} else {print "SNP"}}' | sort -V | uniq -c



cat vcf/anno_merged_filtered_sorted.vcf \
| java -jar snpEff/SnpSift.jar filter "(countVariant() > 60) & ANN[*].IMPACT = 'HIGH' & (QUAL >= 60)" \
| awk '$1~"^CP14"' | cut -f 1,2,8 | awk -F ';' '{print $1}' OFS='\t'



# MultiQC report generation
multiqc --outdir vcf/anno/multiqc_report \
--filename vcf_report \
--clean-up vcf/anno

multiqc --outdir vcf/anno/multiqc_report_wo_up_down_stream \
--filename vcf_report_wo_up_down_stream \
--clean-up vcf/anno

multiqc --outdir vcf/anno/multiqc_report \
--filename final_vcf_snpeff_report \
--clean-up vcf/anno