#!/bin/bash

# script to perform read trimming for tradis
# and run tradis with optimisation

#SBATCH --job-name=RUNtradis
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --time=08:00:00
#SBATCH --output=RUNtradis_%j.out
#SBATCH --error=RUNtradis_%j.err

set -euo pipefail

###############################################################################
# Usage
###############################################################################

usage() {
    cat <<'EOF'
Usage:
    sbatch RUNtradis.sh INPUT_FASTQ OUTPUT_DIRECTORY

Arguments:
    INPUT_FASTQ         FASTQ file or directory containing FASTQ files
    OUTPUT_DIRECTORY    Output directory name (must not already exist)

Example:
    sbatch RUNtradis.sh \
        /path/to/fastq_directory \
        results
EOF

    exit 1
}

###############################################################################
# Parse arguments
###############################################################################

if [[ $# -ne 2 ]]; then
    usage
fi

INPUT_FASTQ="$1"
OUTPUT_DIRECTORY="$2"

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
module load fastp-uoneasy/0.23.4-GCC-12.3.0
module load cutadapt-uon/gcc12.3.0/4.6

###############################################################################
# Validate software and input
###############################################################################

if [[ ! -x "fastqc" ]]; then
    echo "ERROR: fastqc executable was not found or is not executable:" >&2
    echo "       fastqc" >&2
    exit 1
fi

if [[ ! -x "multiqc" ]]; then
    echo "ERROR: multiqc executable was not found or is not executable:" >&2
    echo "       multiqc" >&2
    exit 1
fi

if [[ ! -x "fastp" ]]; then
    echo "ERROR: fastp executable was not found or is not executable:" >&2
    echo "       fastp" >&2
    exit 1
fi

if [[ ! -x "cutadapt" ]]; then
    echo "ERROR: cutadapt executable was not found or is not executable:" >&2
    echo "       cutadapt" >&2
    exit 1
fi

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

mkdir -p "$OUTPUT_DIRECTORY"

OUTPUT_DIRECTORY="$(
    cd "$OUTPUT_DIRECTORY"
    pwd -P
)"

mkdir -p "$OUTPUT_DIRECTORY"/reports/fastqc
mkdir -p "$OUTPUT_DIRECTORY"/reports/multiqc

#OUTPUT_BAM="${OUTPUT_DIRECTORY}/$(basename "$OUTPUT_BAM")"

#if [[ -e "$OUTPUT_BAM" ]]; then
#    echo "ERROR: Output BAM already exists:" >&2
#    echo "       $OUTPUT_BAM" >&2
#    echo "Remove it or choose a different output filename." >&2
#    exit 1
#fi

if [[ ! -w "$OUTPUT_DIRECTORY" ]]; then
    echo "ERROR: Output directory is not writable:" >&2
    echo "       $OUTPUT_DIRECTORY" >&2
    exit 1
fi

###############################################################################
# Temporary-directory configuration
###############################################################################

###  TMP_ROOT="${TMPDIR:-/tmp}"
###  JOB_TMP="${TMP_ROOT}/dorado_${SLURM_JOB_ID:-manual}"
###  
###  mkdir -p "$JOB_TMP"
###  
###  export TMPDIR="$JOB_TMP"
###  export TMP="$JOB_TMP"
###  export TEMP="$JOB_TMP"
###  
###  cleanup() {
###      local status=$?
###  
###      echo
###  
###      if [[ $status -eq 0 ]]; then
###          rm -rf "$JOB_TMP"
###          echo "Temporary directory removed: $JOB_TMP"
###      else
###          echo "Dorado exited with status $status." >&2
###          echo "Temporary directory retained: $JOB_TMP" >&2
###      fi
###  
###      return "$status"
###  }
###  
###  trap cleanup EXIT

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
echo "fastqc version:      $(fastqc --version 2>&1 | head -n 1)"
echo "multiqc version:     $(multiqc --version 2>&1 | head -n 1)"
echo "fastp version:       $(fastp --version 2>&1 | head -n 1)"
echo "cutadapt version:    $(cutadapt --version 2>&1 | head -n 1)"
echo "Input FASTQ:         $INPUT_FASTQ"
echo "Output directory:    $OUTPUT_DIRECTORY"
echo
echo "CPU threads:         $THREADS"
#echo "Temporary directory: $JOB_TMP"
echo

###############################################################################
# Construct commands
###############################################################################

FASTQC_COMMAND=(
    fastqc
    "$INPUT_FASTQ"
    -o "$OUTPUT_DIRECTORY"/reports/fastqc
    -t "$THREADS"
)

MULTIQC_COMMAND=(
    multiqc
    "$OUTPUT_DIRECTORY"/reports/fastqc
    "$OUTPUT_DIRECTORY"/reports/multiqc
)

###############################################################################
# Run commands 
###############################################################################

echo "Running fastqc command:"
printf ' %q' "${FASTQC_COMMAND[@]}"
"${FASTQC_COMMAND[@]}"

###################

echo "Running multiqc command:"
printf ' %q' "${MULTI_COMMAND[@]}"
"${MULTIQC_COMMAND[@]}"

###############################################################################
# Validate output
###############################################################################

###  if [[ ! -s "$OUTPUT_BAM" ]]; then
###      echo "ERROR: Dorado completed but the output BAM is empty:" >&2
###      echo "       $OUTPUT_BAM" >&2
###      exit 1
###  fi

echo
echo "fastqc and multiqc completed successfully."
echo "Completion time: $(date)"
#echo "Output BAM:      $OUTPUT_BAM"
#echo "Output size:     $(du -h "$OUTPUT_BAM" | cut -f1)"

###############################################################################
# Cleanup environment
###############################################################################

module unload fastqc-uoneasy/0.12.1-Java-11
module unload fastp-uoneasy/0.23.4-GCC-12.3.0
module unload cutadapt-uon/gcc12.3.0/4.6

