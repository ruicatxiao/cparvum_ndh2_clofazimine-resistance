
# cparvum_ndh2_clofazimine-resistance

Scripts used for the manuscript ttiled "Widespread genomic heterogeneity at the NAD(P)H dehydrogenase 2 locus predisposes Cryptosporidium to clofazimine resistance"

Here is a collection of python, shell and R script used for varies analysis

## Manuscript Abstract

```
PLACEHOLDER
```


## Requirements

All tools should be available in your `$PATH`:

| Software | Version | Purpose |
|----------|---------|---------|
| Python | 3.10.11+ | Pipeline execution |
| SRA Toolkit | 3.0.5+ | SRA data download and extraction |
| Trim Galore | 0.6.10+ | Read quality trimming |
| BWA | 0.7.17+ | Short-read DNA alignment |
| SnpEff | 5.2e | Genetic variant annotation, and functional effect prediction |
| SnpSift | 5.2e+ | Annotated variants manipulation |
| samtools | 1.19.2+ | SAM/BAM manipulation |
| bcftools | 1.20+ | Variant calling and filtering |
| R | 4.4.2+ | Result visualization through various pkgs |
| GATK | 4.6+ | Variant software suite |

## Per Folder Scripts and Files Overview
### amplicon_analysis
    amplicon_process.sh ==> For processing amplicon-seq reads 
### bsa_analysis
    bsa_genotyping.filtering.sh ==> For initial bsa pre-processing
    cfz.R ==> For R based segregant analysis and plotting
### deprecated
    goi_af.sh ==> Shell script for analyzing target genes' high impact variants' allele frequency
    goi_cov.sh ==> Shell script for analyzing target genes' high impact variants' coverage
    snpeff_parsing.py ==> Initial python script for processing snpeff output
### essential_nonessential_genes_analysis
    goi_af.py ==> Python script that parse out SnpEff and SnpSift output for gene-of-interest high impact variants' allele frequency
    goi_af.py ==> Python script that parse out SnpEff and SnpSift output for gene-of-interest high impact variants' coverage
    snpeff_loop.sh ==> Bash script used for annotating vairants' effect via SnpEff
    snpeff_snpsift.sh ==> Varies commands used for building SnpEff database with Cparvum genome, and preliminary analysis
    snpeff_snpsift_qc_plotting.R ==> R plotting scripts used for QC plots on initial SnpEff outputs
### cmel_analysis
    cmel_analysis.sh ==> Bash script for commands used to perform AA INDEL analysis with Cmel Illumina reads on newest Cmel ref genome 

## Citation

Please cite our paper:
[Insert citation information here]

## Data availability
All data used for this manuscript can be accessed from PRJNA1336748 (WGS) and PRJNA1337473 (Amplicon-seq)
