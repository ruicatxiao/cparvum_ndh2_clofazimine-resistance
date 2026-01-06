# Scripts for checking AA INDEL for newest Cmel genome assembly
# Gene name for cgd7_1900 (cpbgf_7001900) in Cmel is cmbei_7001900 

# UKMEL3 SRR6813720
# UKMEL4 SRR6813896
# CMEL_piglet SRR793561


fasterq-dump --split-3 --skip-technical --progress --temp . --outdir fastq SRR793561

# Chr7 is CM077713.1
# cmbei_7001900 gene range is: 494885  496591

# 494964 is the Cmel new ref the location of AA deletion location


# Setting up sample sheet
sampleType,moleculeType,libraryType,techType,sampleID,read1,read2
SRA,DNA,Paired,Illumina,SRR6813720,,
SRA,DNA,Paired,Illumina,SRR6813896,,
LOCAL,DNA,Paired,Illumina,SRR793561,/data/ruicatxiao/snp_test/fastq/SRR793561_1.fastq,/data/ruicatxiao/snp_test/fastq/SRR793561_2.fastq

# Setting up variant_list.txt
CHR LOCATION
CM077713.1  494964

# Running the program
python sra2genesnv.py \
    --ref_genome cmel.fasta \
    --gtf cmel.gtf \
    --sample cmel_samplesheet.csv \
    --variant_list variant_list.txt \
    --threads 32

# Output
sampleID    moleculeType    libraryType techType    IDV DP  IMF AD  AF
SRR6813720  DNA Paired  Illumina    NA  104 NA  104,0   0.0000
SRR6813896  DNA Paired  Illumina    NA  164 NA  164,0   0.0000
SRR793561   DNA Paired  Illumina    NA  618 NA  618,0   0.0000

# Conclusion
# Cmel for these 3 samples do not have any AA INDEL. They have sufficient coverage at the location for this to be statiscally significant

