#!/bin/bash

set -euo pipefail  

# Help function
show_help() {
    cat << EOF
Usage: $0 [OPTIONS]

Pipeline for mapping, sorting, deduplicating, and variant calling of sequencing data.

OPTIONS:
    -h, --help          Show this help message
    -c, --config        Configuration file (optional, see CONFIGURATION section)

CONFIGURATION:
    Before running this script, you must ensure the following tools are available in your PATH:
    - bwa (tested with version 0.7.18)
    - picard (tested with version 2.27.4)
    - gatk (tested with version 3.8)

    The following files must exist in the current directory:
    - CpBGF_genome_v16.fasta (reference genome)
    - Input FASTQ files in the format: {sample}_R1_001.fastq.gz and {sample}_R2_001.fastq.gz

    The following directories will be created automatically:
    - SAM: for SAM alignment files
    - sorted.bam: for sorted BAM files
    - dedup.sorted.bam: for deduplicated BAM files
    - metrics: for duplicate metrics
    - g.vcf: for gVCF output files

REQUIRED TOOLS:
    This script requires the following tools to be installed and available in your PATH:
    1. bwa (alignment)
    2. picard (sorting, duplicate marking, indexing)
    3. gatk (variant calling)

    Example PATH setup:
    export PATH="/path/to/bwa:\$PATH"
    export PATH="/path/to/picard:\$PATH"
    export PATH="/path/to/gatk:\$PATH"

SAMPLE CONFIGURATION:
    Edit the 'samples' array in the script to include your sample names.
    Sample names should correspond to input files named:
    {sample_name}_R1_001.fastq.gz and {sample_name}_R2_001.fastq.gz

EXAMPLE USAGE:
    $0

EOF
}

# Configuration - UPDATE THESE VALUES BEFORE RUNNING
readonly REFERENCE_GENOME="CpBGF_genome_v16.fasta" 
readonly THREADS=6  # Number of threads to use

# Sample identifiers - UPDATE THIS ARRAY WITH YOUR SAMPLE NAMES
readonly samples=(
    "04_Seb_Cp_Tom_S4"
    "05_Seb_Cp_Isra_S5" 
    "Vehicle_D_6_10_S10"
    "Vehicle_D_13_16_S11"
    "Vehicle_D_26_28_S12"
    "CFZ_D_6_10_S7"
    "CFZ_D_13_16_S8"
    "CFZ_D_26_28_S9"
)

# Logging functionlog() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

error_exit() {
    log "ERROR: $1"
    exit 1
}

# Check if required tools are in PATH
check_tools_in_path() {
    local missing_tools=()
    
    command -v bwa >/dev/null 2>&1 || missing_tools+=("bwa")
    command -v java >/dev/null 2>&1 || missing_tools+=("java")
    
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        log "Missing required tools in PATH:"
        printf '%s\n' "${missing_tools[@]}" >&2
        error_exit "Please install and add the missing tools to your PATH"
    fi
    
    # Check if picard and gatk jars exist by testing if java can find them
    if ! java -jar picard.jar --version >/dev/null 2>&1 && ! java -cp "$(ls -1 picard-*.jar 2>/dev/null | head -n1)" picard.PicardCommandLine --version >/dev/null 2>&1; then
        log "WARNING: Could not verify picard is available. Make sure picard.jar is in PATH or current directory."
    fi
    
    if ! java -jar GenomeAnalysisTK.jar --version >/dev/null 2>&1; then
        log "WARNING: Could not verify GATK is available. Make sure GenomeAnalysisTK.jar is in PATH or current directory."
    fi
}

# Check if required files exist
check_prerequisites() {
    local missing_files=()
    
    # Check if reference genome exists
    [[ ! -f "$REFERENCE_GENOME" ]] && missing_files+=("$REFERENCE_GENOME")
    
    # Check if input files exist for each sample
    for sample in "${samples[@]}"; do
        local r1="${sample}_R1_001.fastq.gz"
        local r2="${sample}_R2_001.fastq.gz"
        [[ ! -f "$r1" ]] && missing_files+=("$r1")
        [[ ! -f "$r2" ]] && missing_files+=("$r2")
    done
    
    if [[ ${#missing_files[@]} -gt 0 ]]; then
        log "Missing required files:"
        printf '%s\n' "${missing_files[@]}" >&2
        error_exit "Prerequisites check failed. Please ensure all required files exist."
    fi
}

# Create output directories
create_directories() {
    local dirs=("SAM" "sorted.bam" "dedup.sorted.bam" "metrics" "g.vcf")
    for dir in "${dirs[@]}"; do
        mkdir -p "$dir" || error_exit "Failed to create directory: $dir"
    done
}

# Align reads with BWA
align_with_bwa() {
    local sample="$1"
    local output_sam="SAM/${sample}.sam"
    
    log "Aligning ${sample}..."
    
    bwa mem \
        -t "$THREADS" \
        -M \
        -R "@RG\tID:${sample}\tLB:${sample}\tPL:ILLUMINA\tPM:HISEQ\tSM:${sample}/" \
        "$REFERENCE_GENOME" \
        "${sample}_R1_001.fastq.gz" \
        "${sample}_R2_001.fastq.gz" \
        > "$output_sam" || error_exit "BWA alignment failed for ${sample}"
}

# Sort SAM file
sort_sam() {
    local sample="$1"
    local input_sam="SAM/${sample}.sam"
    local output_bam="sorted.bam/${sample}.sorted.bam"
    
    log "Sorting ${sample}..."
    
    java -jar picard.jar SortSam \
        -INPUT "$input_sam" \
        -OUTPUT "$output_bam" \
        -SORT_ORDER coordinate \
        || error_exit "Sorting failed for ${sample}"
}

# Mark duplicates
mark_duplicates() {
    local sample="$1"
    local input_bam="sorted.bam/${sample}.sorted.bam"
    local output_bam="dedup.sorted.bam/${sample}.dedup.sorted.bam"
    local metrics_file="metrics/${sample}.metrics.txt"
    
    log "Marking duplicates for ${sample}..."
    
    java -jar picard.jar MarkDuplicates \
        -INPUT "$input_bam" \
        -OUTPUT "$output_bam" \
        -METRICS_FILE "$metrics_file" \
        || error_exit "MarkDuplicates failed for ${sample}"
}

# Build BAM index
build_bam_index() {
    local sample="$1"
    local input_bam="dedup.sorted.bam/${sample}.dedup.sorted.bam"
    
    log "Building index for ${sample}..."
    
    java -jar picard.jar BuildBamIndex \
        -INPUT "$input_bam" \
        || error_exit "BuildBamIndex failed for ${sample}"
}

# Call variants with GATK
call_variants() {
    local sample="$1"
    local input_bam="dedup.sorted.bam/${sample}.dedup.sorted.bam"
    local output_gvcf="g.vcf/${sample}.g.vcf"
    
    log "Calling variants for ${sample}..."
    
    java -jar GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -R "$REFERENCE_GENOME" \
        -I "$input_bam" \
        --emitRefConfidence BP_RESOLUTION \
        -o "$output_gvcf" \
        --useNewAFCalculator \
        --sample_ploidy 1 \
        --dontUseSoftClippedBases \
        --num_cpu_threads_per_data_thread "$THREADS" \
        || error_exit "Variant calling failed for ${sample}"
    
    log "${sample} processing complete!"
}

# Main processing function
process_sample() {
    local sample="$1"
    
    log "Starting processing for ${sample}"
    
    # Validate input files exist
    [[ ! -f "${sample}_R1_001.fastq.gz" ]] && error_exit "R1 file missing for ${sample}: ${sample}_R1_001.fastq.gz"
    [[ ! -f "${sample}_R2_001.fastq.gz" ]] && error_exit "R2 file missing for ${sample}: ${sample}_R2_001.fastq.gz"
    
    # Process the sample
    align_with_bwa "$sample"
    sort_sam "$sample"
    mark_duplicates "$sample"
    build_bam_index "$sample"
    call_variants "$sample"
    
    log "Completed processing for ${sample}"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -c|--config)
            log "Configuration file option not implemented in this version"
            shift
            ;;
        *)
            log "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
    shift
done

# Main execution
main() {
    log "Starting pipeline execution"
    log "Ensure all required tools are in your PATH before continuing."
    log "See help (-h) for requirements and setup instructions."
    
    # Check tools are in PATH
    check_tools_in_path
    
    # Check prerequisites
    check_prerequisites
    
    # Create directories
    create_directories
    
    # Process each sample
    for sample in "${samples[@]}"; do
        process_sample "$sample" || {
            log "Failed to process sample: $sample"
            continue
        }
    done
    
    log "Pipeline completed successfully"
}

# Run main function
main "$@"