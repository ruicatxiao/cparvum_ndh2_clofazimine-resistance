#!/bin/bash

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Help function
show_help() {
    cat << EOF
Usage: $0 [OPTIONS]

Pipeline for filtering and analyzing VCF files from variant calling.

OPTIONS:
    -h, --help          Show this help message
    -m, --memory        Maximum memory for Java (default: 50g)
    -r, --reference     Reference genome file (default: CpBGF_genome_v16.fasta)
    -i, --input-dir     Input directory containing VCF files (default: g.vcf)
    -o, --output-dir    Output directory for results (default: g.vcf)

CONFIGURATION:
    Before running this script, ensure the following:
    - vcflib tools are available in your PATH (vcffilter)
    - Java is available in your PATH
    - GATK jar file is available in your PATH or current directory
    - Reference genome file exists
    - Required input VCF files exist (see INPUT FILES section)

REQUIRED TOOLS:
    This script requires the following tools to be available:
    1. vcflib (vcffilter command)
    2. java
    3. GATK (GenomeAnalysisTK.jar)

INPUT FILES:
    The script expects the following input files in the input directory:
    - CPCP.parents.SNP.vcf (parental variants)
    - CFZ.SNP.vcf (comparison sample variants)

    These files should be the output from a previous variant calling pipeline.

OUTPUT FILES:
    The script will generate:
    - CPCP.parents.SNP.hardfilter.vcf (hard-filtered variants)
    - CPCP.parents.SNP.filtered.vcf (allele frequency filtered variants)
    - CFZ.SNP.BSA.filter.vcf (concordance filtered variants)
    - CFZ.SNP.BSA.filter.table (tabular format of final variants)

EXAMPLE USAGE:
    $0
    $0 -m 30g -r my_genome.fasta -i input_vcf_dir -o output_dir

NOTES:
    - Memory requirements may vary based on input file sizes
    - The script assumes biallelic variants for the comparison steps
    - Allele frequency filter: 0.25 < AF < 1.0

EOF
}

# Default configuration
MEMORY="50g"
REFERENCE_GENOME="CpBGF_genome_v16.fasta"
INPUT_DIR="g.vcf"
OUTPUT_DIR="g.vcf"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -m|--memory)
            MEMORY="$2"
            shift 2
            ;;
        -r|--reference)
            REFERENCE_GENOME="$2"
            shift 2
            ;;
        -i|--input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1" >&2
            show_help
            exit 1
            ;;
    esac
done

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

# Error handling function
error_exit() {
    log "ERROR: $1"
    exit 1
}

# Check if required tools are in PATH
check_tools_in_path() {
    local missing_tools=()
    
    command -v vcffilter >/dev/null 2>&1 || missing_tools+=("vcffilter (from vcflib)")
    command -v java >/dev/null 2>&1 || missing_tools+=("java")
    
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        log "Missing required tools in PATH:"
        printf '%s\n' "${missing_tools[@]}" >&2
        error_exit "Please install and add the missing tools to your PATH"
    fi
    
    # Check if GATK jar is available
    if ! java -jar GenomeAnalysisTK.jar --version >/dev/null 2>&1; then
        error_exit "GenomeAnalysisTK.jar not found in PATH or current directory"
    fi
}

# Check if required files exist
check_prerequisites() {
    local missing_files=()
    
    # Check if reference genome exists
    [[ ! -f "$REFERENCE_GENOME" ]] && missing_files+=("$REFERENCE_GENOME")
    
    # Check if input VCF files exist
    [[ ! -f "$INPUT_DIR/CPCP.parents.SNP.vcf" ]] && missing_files+=("$INPUT_DIR/CPCP.parents.SNP.vcf")
    [[ ! -f "$INPUT_DIR/CFZ.SNP.vcf" ]] && missing_files+=("$INPUT_DIR/CFZ.SNP.vcf")
    
    # Create output directory if it doesn't exist
    mkdir -p "$OUTPUT_DIR" || error_exit "Failed to create output directory: $OUTPUT_DIR"
    
    if [[ ${#missing_files[@]} -gt 0 ]]; then
        log "Missing required files:"
        printf '%s\n' "${missing_files[@]}" >&2
        error_exit "Prerequisites check failed. Please ensure all required files exist."
    fi
}

# Apply hard filters to parental variants
apply_hard_filters() {
    local input_vcf="$INPUT_DIR/CPCP.parents.SNP.vcf"
    local output_vcf="$OUTPUT_DIR/CPCP.parents.SNP.hardfilter.vcf"
    
    log "Applying hard filters to parental variants..."
    log "Filtering: QD > 2.0 & FS < 60.0 & SOR < 3.0, Genotype: DP > 10 & GQ > 90"
    
    vcffilter -f "QD > 2.0 & FS < 60.0 & SOR < 3.0" -g "DP > 10 & GQ > 90" "$input_vcf" > "$output_vcf" || \
        error_exit "Hard filtering failed"
    
    local variant_count=$(grep -c -v "^#" "$output_vcf" 2>/dev/null || echo 0)
    log "Hard filtering complete. Variants remaining: $variant_count"
}

# Filter variants by allele frequency
filter_by_allele_frequency() {
    local input_vcf="$OUTPUT_DIR/CPCP.parents.SNP.hardfilter.vcf"
    local output_vcf="$OUTPUT_DIR/CPCP.parents.SNP.filtered.vcf"
    
    log "Filtering variants by allele frequency (0.25 < AF < 1.0)..."
    
    java -Xmx"$MEMORY" -jar GenomeAnalysisTK.jar \
        -T SelectVariants \
        -R "$REFERENCE_GENOME" \
        -V "$input_vcf" \
        --excludeNonVariants \
        --maxNOCALLnumber 0 \
        --removeUnusedAlternates \
        --restrictAllelesTo BIALLELIC \
        -select "AF > 0.25 && AF < 1.0" \
        -o "$output_vcf" || error_exit "Allele frequency filtering failed"
    
    local variant_count=$(grep -c -v "^#" "$output_vcf" 2>/dev/null || echo 0)
    log "Allele frequency filtering complete. Variants remaining: $variant_count"
}

# Find concordant variants between samples
find_concordant_variants() {
    local input_vcf="$INPUT_DIR/CFZ.SNP.vcf"
    local concordance_vcf="$OUTPUT_DIR/CPCP.parents.SNP.filtered.vcf"
    local output_vcf="$OUTPUT_DIR/CFZ.SNP.BSA.filter.vcf"
    
    log "Finding concordant variants between CFZ and parental samples..."
    
    java -Xmx"$MEMORY" -jar GenomeAnalysisTK.jar \
        -T SelectVariants \
        -R "$REFERENCE_GENOME" \
        -V "$input_vcf" \
        --concordance "$concordance_vcf" \
        --restrictAllelesTo BIALLELIC \
        -o "$output_vcf" || error_exit "Concordance filtering failed"
    
    local variant_count=$(grep -c -v "^#" "$output_vcf" 2>/dev/null || echo 0)
    log "Concordance filtering complete. Variants remaining: $variant_count"
}

# Convert VCF to table format
convert_to_table() {
    local input_vcf="$OUTPUT_DIR/CFZ.SNP.BSA.filter.vcf"
    local output_table="$OUTPUT_DIR/CFZ.SNP.BSA.filter.table"
    
    log "Converting final VCF to tabular format..."
    
    java -Xmx"$MEMORY" -jar GenomeAnalysisTK.jar \
        -T VariantsToTable \
        -R "$REFERENCE_GENOME" \
        -V "$input_vcf" \
        -F CHROM -F POS -F REF -F ALT \
        -GF AD -GF DP -GF GQ -GF PL \
        -o "$output_table" || error_exit "VCF to table conversion failed"
    
    local line_count=$(wc -l < "$output_table" 2>/dev/null || echo 0)
    # Subtract 1 for header line
    local variant_count=$((line_count - 1))
    log "Table conversion complete. Variants in table: $variant_count"
    log "Output table saved to: $output_table"
}

# Main execution
main() {
    log "Starting VCF filtering and analysis pipeline"
    log "Memory: $MEMORY"
    log "Reference genome: $REFERENCE_GENOME"
    log "Input directory: $INPUT_DIR"
    log "Output directory: $OUTPUT_DIR"
    
    # Check tools are in PATH
    check_tools_in_path
    
    # Check prerequisites
    check_prerequisites
    
    # Execute pipeline steps
    apply_hard_filters
    filter_by_allele_frequency
    find_concordant_variants
    convert_to_table
    
    log "Pipeline completed successfully!"
    log "Final results:"
    log "  - Hard-filtered VCF: $OUTPUT_DIR/CPCP.parents.SNP.hardfilter.vcf"
    log "  - Allele frequency filtered VCF: $OUTPUT_DIR/CPCP.parents.SNP.filtered.vcf"
    log "  - Concordance filtered VCF: $OUTPUT_DIR/CFZ.SNP.BSA.filter.vcf"
    log "  - Final table: $OUTPUT_DIR/CFZ.SNP.BSA.filter.table"
}

# Run main function
main "$@"