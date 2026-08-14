#!/bin/bash

# script to perform read trimming for tradis
# and run tradis with optimisation

#SBATCH --job-name=RUNtradis
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --output=RUNtradis_%j.out

START_TIME=$(date +%s)
set -euo pipefail

###############################################################################
# Usage
###############################################################################

usage() {
    cat <<'EOF'
Usage:
    sbatch RUNtradis.sh --input INPUT_FASTQ --output OUTPUT_DIRECTORY [OPTIONS]

Required arguments:
    -i, --input       FASTQ file or directory containing FASTQ files
    -o, --output      Output directory name (must not already exist)

Optional arguments:
    -r, --reference   Reference genome in FASTA format
                      Default: /share/bryant_lab/reference_genomes/GCF_000750555.1_ASM75055v1_genomic.fna

    -a, --annotation  Genome annotation file
                      Default: /share/bryant_lab/reference_genomes/GCF_000750555.1_ASM75055v1_genomic.gff

    -t, --tag         DNA sequence of transposon tag
                      Default: CGAGCTCGAATTCATCGATGATGGTTGAGATGTGTATAAGAGACAG

    -h, --help        Show this help message

Example:
    sbatch RUNtradis.sh \
        --input /path/to/fastq_directory \
        --output results \
        --reference /path/to/reference.fasta \
        --annotation /path/to/annotation.gff \
        --tag CGAGCTCGAATTCATCGATGATGGTTGAGATGTGTATAAGAGACAG
EOF
    exit 1
}

###############################################################################
# Default values
###############################################################################

REFERENCE_GENOME="/share/bryant_lab/reference_genomes/GCF_000750555.1_ASM75055v1_genomic.fna"
GENOME_ANNOTATION="/share/bryant_lab/reference_genomes/GCF_000750555.1_ASM75055v1_genomic.gff"
TRANSPOSON_TAG="CGAGCTCGAATTCATCGATGATGGTTGAGATGTGTATAAGAGACAG"

###############################################################################
# Parse arguments
###############################################################################

TEMP=$(getopt \
    --options i:o:r:a:t:h \
    --longoptions input:,output:,reference:,annotation:,tag:,help \
    --name "$0" \
    -- "$@"
)

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to parse arguments." >&2
    usage
fi

eval set -- "$TEMP"

while true; do
    case "$1" in
        -i|--input)
            INPUT_FASTQ="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIRECTORY="$2"
            shift 2
            ;;
        -r|--reference)
            REFERENCE_GENOME="$2"
            shift 2
            ;;
        -a|--annotation)
            GENOME_ANNOTATION="$2"
            shift 2
            ;;
        -t|--tag)
            TRANSPOSON_TAG="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "ERROR: Unexpected argument: $1" >&2
            usage
            ;;
    esac
done

###############################################################################
# Check required arguments
###############################################################################

if [[ -z "${INPUT_FASTQ:-}" ]]; then
    echo "ERROR: --input is required." >&2
    usage
fi

if [[ -z "${OUTPUT_DIRECTORY:-}" ]]; then
    echo "ERROR: --output is required." >&2
    usage
fi

###############################################################################
# Configuration
###############################################################################

THREADS="${SLURM_CPUS_PER_TASK:-16}"

###############################################################################
# Load environment
###############################################################################

echo "Loading software environment..."

module --force purge 2>/dev/null || module purge

module load fastqc-uoneasy/0.12.1-Java-11
module load multiqc-uoneasy/1.14-foss-2023a
module load fastp-uoneasy/0.23.4-GCC-12.3.0
module load cutadapt-uon/gcc12.3.0/4.6

echo
echo "... software loaded successfully"

###############################################################################
# Validate software and input
###############################################################################

check_command() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "ERROR: $1 was not found in PATH." >&2
        echo "       Please load the appropriate module." >&2
        exit 1
    fi
}

check_command fastqc
check_command multiqc
check_command fastp
check_command cutadapt

if [[ ! -e "$INPUT_FASTQ" ]]; then
    echo "ERROR: Input FASTQ file or directory does not exist:" >&2
    echo "       $INPUT_FASTQ" >&2
    exit 1
fi

if [[ ! -r "$INPUT_FASTQ" ]]; then
    echo "ERROR: Input FASTQ file or directory is not readable:" >&2
    echo "       $INPUT_FASTQ" >&2
    exit 1
fi

###############################################################################
# Prepare output path
###############################################################################

if [[ -e "$OUTPUT_DIRECTORY" ]]; then
    echo "ERROR: Output directory already exists:" >&2
    echo "       $OUTPUT_DIRECTORY" >&2
    echo "Remove it or choose a different output directory name." >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIRECTORY"

OUTPUT_DIRECTORY="$(
    cd "$OUTPUT_DIRECTORY"
    pwd -P
)"

mkdir -p "$OUTPUT_DIRECTORY"/reports/fastqc
mkdir -p "$OUTPUT_DIRECTORY"/reports/multiqc
mkdir -p "$OUTPUT_DIRECTORY"/reports/fastp
mkdir -p "$OUTPUT_DIRECTORY"/reports/cutadapt
mkdir -p "$OUTPUT_DIRECTORY"/trimmed_fastqs
mkdir -p "$OUTPUT_DIRECTORY"/trimmed_fastqs/1_fastp
mkdir -p "$OUTPUT_DIRECTORY"/trimmed_fastqs/2_cutadapt
mkdir -p "$OUTPUT_DIRECTORY"/biotradis

if [[ ! -w "$OUTPUT_DIRECTORY" ]]; then
    echo "ERROR: Output directory is not writable:" >&2
    echo "       $OUTPUT_DIRECTORY" >&2
    exit 1
fi

###############################################################################
# Report configuration
###############################################################################

echo
echo "RUNtradis job"
echo "============="
echo
echo "Job ID:              ${SLURM_JOB_ID:-not-running-under-slurm}"
echo "Job name:            ${SLURM_JOB_NAME:-unknown}"
echo "Compute node:        $(hostname)"
echo "Start time:          $(date)"
echo
echo "FastQC location:     $(command -v fastqc)"
echo "FastQC version:      $(fastqc --version 2>&1 | head -n 1)"
echo "MultiQC location:    $(command -v multiqc)"
echo "MultiQC version:     $(multiqc --version 2>&1 | head -n 1)"
echo "fastp location:      $(command -v fastp)"
echo "fastp version:       $(fastp --version 2>&1 | head -n 1)"
echo "cutadapt location:   $(command -v cutadapt)"
echo "cutadapt version:    $(cutadapt --version 2>&1 | head -n 1)"
echo
echo "Input FASTQ:         $INPUT_FASTQ"
echo "Output directory:    $OUTPUT_DIRECTORY"
echo "Reference genome:    $REFERENCE_GENOME"
echo "Transposon tag:      $TRANSPOSON_TAG"
echo
echo "CPU threads:         $THREADS"
echo

###############################################################################
# Construct commands
###############################################################################

FASTQC_COMMAND=(
    fastqc
    "$INPUT_FASTQ"
    -o "$OUTPUT_DIRECTORY"/reports/fastqc
    -t "$THREADS")

MULTIQC_COMMAND=(
    multiqc
    "$OUTPUT_DIRECTORY"/reports/fastqc
    -o "$OUTPUT_DIRECTORY"/reports/multiqc)

FASTP_COMMAND=(
    fastp
    -i "$INPUT_FASTQ"
    -o "$OUTPUT_DIRECTORY"/trimmed_fastqs/1_fastp/$(basename ${INPUT_FASTQ})
    --disable_quality_filtering
    --disable_adapter_trimming
    --disable_length_filtering
    --trim_poly_g
    --poly_g_min_len 5
    --trim_poly_x
    --poly_x_min_len 4
    --thread $THREADS
    --html "$OUTPUT_DIRECTORY"/reports/fastp/$(basename ${INPUT_FASTQ})_fastp.html
    --json "$OUTPUT_DIRECTORY"/reports/fastp/$(basename ${INPUT_FASTQ})_fastp.json)

CUTADAPT_COMMAND=(
    cutadapt
    --cores $THREADS
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    --revcomp
    --poly-a
    --minimum-length 50
    --info-file "$OUTPUT_DIRECTORY"/reports/cutadapt/$(basename ${INPUT_FASTQ})_info.tsv
    -o "$OUTPUT_DIRECTORY"/trimmed_fastqs/2_cutadapt/$(basename ${INPUT_FASTQ})
    "$OUTPUT_DIRECTORY"/trimmed_fastqs/1_fastp/$(basename ${INPUT_FASTQ}))

TRADIS_COMMAND=(
    bacteria_tradis
    -v
    --smalt
    --smalt_r 0
    --smalt_k 10
    --smalt_s 1
    --smalt_y .90
    -m 0
    -mm 15
    -f "$OUTPUT_DIRECTORY"/biotradis/files.txt
    -t "$TRANSPOSON_TAG"
    -r "$REFERENCE_GENOME")

TRADIS_GIS_COMMAND=(
    tradis_gene_insert_sites
    "$GENOME_ANNOTATION"
    "$OUTPUT_DIRECTORY"/biotradis/*.insert_site_plot.gz
    )

###############################################################################
# Run commands 
###############################################################################

echo "Running fastqc command:"
printf ' %q' "${FASTQC_COMMAND[@]}"
echo
echo
"${FASTQC_COMMAND[@]}"
echo
echo

###################

echo "Running multiqc command:"
printf ' %q' "${MULTIQC_COMMAND[@]}"
echo
echo
"${MULTIQC_COMMAND[@]}"
echo
echo

###################

echo "Running fastp command:"
printf ' %q' "${FASTP_COMMAND[@]}"
echo
echo
"${FASTP_COMMAND[@]}"
echo
echo

###################

echo "Running cutadapt command:"
printf ' %q' "${CUTADAPT_COMMAND[@]}"
echo
echo
"${CUTADAPT_COMMAND[@]}"
echo
echo

###################

echo "Running tradis command:"
echo "$INPUT_FASTQ" > "$OUTPUT_DIRECTORY"/biotradis/files.txt
printf ' %q' "${TRADIS_COMMAND[@]}" 
echo
echo
source $HOME/.bash_profile
conda activate biotradis
cd "$OUTPUT_DIRECTORY"/biotradis
"${TRADIS_COMMAND[@]}"
echo
echo

###################

echo "Running tradis_gene_insert_sites command:"
printf ' %q' "${TRADIS_GIS_COMMAND[@]}" 
echo
echo
"${TRADIS_GIS_COMMAND[@]}"
conda deactivate
echo
echo

###############################################################################
# Cleanup environment
###############################################################################

module unload fastqc-uoneasy/0.12.1-Java-11
module unload multiqc-uoneasy/1.14-foss-2023a
module unload fastp-uoneasy/0.23.4-GCC-12.3.0
module unload cutadapt-uon/gcc12.3.0/4.6

###############################################################################
# Run summary
###############################################################################

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo
echo "========================================"
echo "RUN SUMMARY"
echo "========================================"
echo
echo "RUNtradis.sh script completed."
echo
printf 'Total run time: %02d:%02d:%02d\n' $((ELAPSED / 3600)) $(((ELAPSED % 3600) / 60)) $((ELAPSED % 60))
echo
echo "Completion time: $(date)"

