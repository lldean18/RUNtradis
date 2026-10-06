#!/bin/bash
# Laura Dean
# 1/10/26

# script to perform read trimming for tradis
# and run tradis with optimization
# and draw a circos plot of the insertion sites

#SBATCH --job-name=RUNtradis
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=06:00:00
#SBATCH --output=RUNtradis_%j.out

START_TIME=$(date +%s)
set -euo pipefail

###############################################################################
# Usage
###############################################################################

usage() {
    cat <<'EOF'
Usage:
    sbatch /share/bryant_lab/RUNtradis.sh --input INPUT_FASTQ --output OUTPUT_DIRECTORY [OPTIONS]

Required arguments:
    -i, --input       input FASTQ file (must be compressed with gzip)
    -o, --output      Output directory name (must not already exist)

Optional arguments:
    -r, --reference   Reference genome in FASTA format
                      Default: /share/bryant_lab/reference_genomes/GCF_000750555.1_ASM75055v1_genomic.fna

    -a, --annotation  Genome annotation file (must be in .embl format) NOTE: for Circos plot creation a .gff version of the annotation must also exist in the same location with the same file prefix as the .embl file
                      Default: /share/bryant_lab/reference_genomes/GCF_000750555.1_ASM75055v1_genomic.embl

    -t, --tag         DNA sequence of transposon tag
                      Default: CGAGCTCGAATTCATCGATGATGGTTGAGATGTGTATAAGAGACAG

    -m, --mismatches  Number of mismatches allowed when matching the transposon tag
                      Default: 6

    -q, --quality     Minimum mapping quality score to use a read
                      Default: 0 (must be 0 if you want to retain multi-mapping reads, e.g. so repetitive regions are not incorrectly labelled essential)

    -p, --percentage  Minimum percentage of identical bases between read and reference genome
                      Default: 0.90 (90%)

    -h, --help        Show this help message

Example:
    sbatch /share/bryant_lab/RUNtradis.sh \
        --input /path/to/file.fastq.gz \
        --output ~/results \
        --reference /path/to/reference.fasta \
        --annotation /path/to/annotation.embl \
        --tag CGAGCTCGAATTCATCGATGATGGTTGAGATGTGTATAAGAGACAG \
        --mismatches 6 \
        --quality 0 \
        --percentage 0.90
EOF
    exit 1
}

###############################################################################
# Default values
###############################################################################

REFERENCE_GENOME="/share/bryant_lab/reference_genomes/GCF_000750555.1_ASM75055v1_genomic.fna"
GENOME_ANNOTATION="/share/bryant_lab/reference_genomes/GCF_000750555.1_ASM75055v1_genomic.embl"
TRANSPOSON_TAG="CGAGCTCGAATTCATCGATGATGGTTGAGATGTGTATAAGAGACAG"
MISMATCHES="6"
QUALITY="0"
PERCENTAGE="0.90"

###############################################################################
# Parse arguments
###############################################################################

TEMP=$(getopt \
    --options i:o:r:a:t:m:q:p:h \
    --longoptions input:,output:,reference:,annotation:,tag:,mismatches:,quality:,percentage:,help \
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
        -m|--mismatches)
            MISMATCHES="$2"
            shift 2
            ;;
        -q|--quality)
            QUALITY="$2"
            shift 2
            ;;
        -p|--percentage)
            PERCENTAGE="$2"
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
module load biotradis-uon/1.4.5
source /gpfs01/software/easybuild5-uon/software/Miniforge3/25.3.0-3/etc/profile.d/conda.sh
conda activate /gpfs01/software/conda-extras/biotradis-1.4.5/envs

echo
echo "... initial software loaded successfully"
echo "due to a GCC version conflict samtools and circos software must be loaded separately at the end."

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
check_command bacteria_tradis

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
mkdir -p "$OUTPUT_DIRECTORY"/trimmed_fastqs/3_cutadapt
mkdir -p "$OUTPUT_DIRECTORY"/biotradis
mkdir -p "$OUTPUT_DIRECTORY"/circos

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
echo "Job ID:                      ${SLURM_JOB_ID:-not-running-under-slurm}"
echo "Job name:                    ${SLURM_JOB_NAME:-unknown}"
echo "Compute node:                $(hostname)"
echo "Start time:                  $(date)"
echo
echo "FastQC location:             $(command -v fastqc)"
echo "FastQC version:              $(fastqc --version 2>&1 | head -n 1)"
echo "MultiQC location:            $(command -v multiqc)"
echo "MultiQC version:             $(multiqc --version 2>&1 | head -n 1)"
echo "fastp location:              $(command -v fastp)"
echo "fastp version:               $(fastp --version 2>&1 | head -n 1)"
echo "cutadapt location:           $(command -v cutadapt)"
echo "cutadapt version:            $(cutadapt --version 2>&1 | head -n 1)"
echo "biotradis location:          $(command -v bacteria_tradis)"
echo "biotradis version:           1.4.5"
echo
echo "Input FASTQ:                 $INPUT_FASTQ"
echo "Output directory:            $OUTPUT_DIRECTORY"
echo "Reference genome:            $REFERENCE_GENOME"
echo "Transposon tag:              $TRANSPOSON_TAG"
echo "Mismatches allowed in tag:   $MISMATCHES"
echo "Min read mapping quality:    $QUALITY"
echo "Min read / ref match:        $PERCENTAGE"
echo
echo "CPU threads:                 $THREADS"
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

CUTADAPT_ADAPTER_COMMAND=(
    cutadapt
    --cores $THREADS
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    --poly-a
    --minimum-length 50
    --info-file "$OUTPUT_DIRECTORY"/reports/cutadapt/$(basename ${INPUT_FASTQ})_adapter_info.tsv
    -o "$OUTPUT_DIRECTORY"/trimmed_fastqs/2_cutadapt/$(basename ${INPUT_FASTQ})
    "$OUTPUT_DIRECTORY"/trimmed_fastqs/1_fastp/$(basename ${INPUT_FASTQ}))

CUTADAPT_TAG_COMMAND=(
    cutadapt
    --cores $THREADS
    -g "$TRANSPOSON_TAG"
    -e "$MISMATCHES"
    --overlap 37
    --action=retain
    --info-file "$OUTPUT_DIRECTORY"/reports/cutadapt/$(basename ${INPUT_FASTQ})_tradis_tag_info.tsv
    -o "$OUTPUT_DIRECTORY"/trimmed_fastqs/3_cutadapt/$(basename ${INPUT_FASTQ})
    "$OUTPUT_DIRECTORY"/trimmed_fastqs/2_cutadapt/$(basename ${INPUT_FASTQ}))

# if you wanted to throw away the reads without the tag add --discard-untrimmed to the cutadapt command above

TRADIS_COMMAND=(
    bacteria_tradis
    -v
    --smalt
    --smalt_r 0
    --smalt_k 10
    --smalt_s 1
    --smalt_y "$PERCENTAGE"
    -m "$QUALITY"
    -mm "$MISMATCHES"
    -f "$OUTPUT_DIRECTORY"/biotradis/files.txt
    -t "$TRANSPOSON_TAG"
    -r "$REFERENCE_GENOME")

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

echo "Running cutadapt 3' adapter removal command:"
printf ' %q' "${CUTADAPT_ADAPTER_COMMAND[@]}"
echo
echo
"${CUTADAPT_ADAPTER_COMMAND[@]}"
echo
echo

###################

echo "Running cutadapt trim sequence before transposon tag command:"
printf ' %q' "${CUTADAPT_TAG_COMMAND[@]}"
echo
echo
"${CUTADAPT_TAG_COMMAND[@]}"
echo
echo

###################

echo "Running tradis command:"
echo "$OUTPUT_DIRECTORY/trimmed_fastqs/3_cutadapt/$(basename ${INPUT_FASTQ})" > "$OUTPUT_DIRECTORY"/biotradis/files.txt
printf ' %q' "${TRADIS_COMMAND[@]}" 
echo
echo
cd "$OUTPUT_DIRECTORY"/biotradis
"${TRADIS_COMMAND[@]}"
echo
echo

###################

# loop for this command to run over multiple files if there are multiple contigs
while IFS= read -r file; do
    echo "Running tradis_gene_insert_sites command for file: $file"
    echo
    
    TRADIS_GIS_COMMAND=(
        tradis_gene_insert_sites
        "$GENOME_ANNOTATION"
        "$file"
        )
    
    echo "Running tradis_gene_insert_sites command:"
    printf ' %q' "${TRADIS_GIS_COMMAND[@]}" 
    echo
    echo
    "${TRADIS_GIS_COMMAND[@]}"
    echo
    echo

done < <(find "$OUTPUT_DIRECTORY/biotradis" -maxdepth 1 -type f -name '*.insert_site_plot.gz')

###################

while IFS= read -r file; do
    echo "Running gene_essentiality command for file: $file"
    echo

    GENE_ESSENTIALITY_COMMAND=(
        tradis_essentiality.R
        "$file"
        )

    echo "Running gene_essentiality command:"
    printf ' %q' "${GENE_ESSENTIALITY_COMMAND[@]}"
    echo
    echo
    "${GENE_ESSENTIALITY_COMMAND[@]}"
    echo
    echo

done < <(find "$OUTPUT_DIRECTORY/biotradis" -maxdepth 1 -type f -name '*.tradis_gene_insert_sites.csv')

#################################################################################
# Cleanup environment up to this point as circos requires a different GCC version
#################################################################################

module unload fastqc-uoneasy/0.12.1-Java-11
module unload multiqc-uoneasy/1.14-foss-2023a
module unload fastp-uoneasy/0.23.4-GCC-12.3.0
module unload cutadapt-uon/gcc12.3.0/4.6
module unload biotradis-uon/1.4.5
conda deactivate

###############################################################################
# Draw circos plot
###############################################################################

# load samtools software
module load samtools-uoneasy/1.22.1-GCC-14.2.0
check_command samtools

echo "samtools location:           $(command -v samtools)"
echo "samtools version:            $(samtools --version 2>&1 | head -n 1)"
echo
echo
echo "indexing the reference genome as the index file is needed for circos plotting..."
echo
echo

# index the reference genome
samtools faidx $REFERENCE_GENOME

# unload samtools
module unload samtools-uoneasy/1.22.1-GCC-14.2.0

# load circos software
module load circos-uoneasy/0.69-9-GCCcore-11.3.0
check_command circos

echo "circos location:           $(command -v circos)"
echo "circos version:            $(circos --version 2>&1 | head -n 1)"
echo
echo

echo "Preparing files for the circos plot..."
echo
echo

# concatenate the contig level insert site plots back into one file for the genome
cut -f1 "${REFERENCE_GENOME}.fai" |
while read contig; do
    cat $OUTPUT_DIRECTORY/biotradis/*."${contig}".insert_site_plot.gz
done > $OUTPUT_DIRECTORY/biotradis/combined.insert_site_plot.gz

#####################
### PREP ASSEMBLY ###
#####################

# generate the karyotype file
awk '{print "chr - " $1 " " $1 " 0 " $2 " chr"}' $REFERENCE_GENOME.fai > $OUTPUT_DIRECTORY/circos/karyotype.txt

##############################
### PREP GENOME ANNOTATION ###
##############################

# convert annotation to circos format
awk 'tolower($3) == "cds" || tolower($3) == "gene"' ${GENOME_ANNOTATION%.*}.gff | awk '{print $1, $4, $5}' OFS="\t" > $OUTPUT_DIRECTORY/circos/genes.txt
sort -k1,1 -k2,2n $OUTPUT_DIRECTORY/circos/genes.txt > $OUTPUT_DIRECTORY/circos/genes.txt.tmp
mv $OUTPUT_DIRECTORY/circos/genes.txt.tmp $OUTPUT_DIRECTORY/circos/genes.txt

# make separate annotation files for genes on fwd and rev strands (strand info is 7th column)

# fwd strand
awk '(tolower($3) == "cds" || tolower($3) == "gene") && $7=="+"' ${GENOME_ANNOTATION%.*}.gff |
awk '{print $1, $4, $5}' OFS="\t" > $OUTPUT_DIRECTORY/circos/genes_fwd_strand.txt
sort -k1,1 -k2,2n $OUTPUT_DIRECTORY/circos/genes_fwd_strand.txt > $OUTPUT_DIRECTORY/circos/genes_fwd_strand.txt.tmp
mv $OUTPUT_DIRECTORY/circos/genes_fwd_strand.txt.tmp $OUTPUT_DIRECTORY/circos/genes_fwd_strand.txt

# rev strand
awk '(tolower($3) == "cds" || tolower($3) == "gene") && $7=="-"' ${GENOME_ANNOTATION%.*}.gff |
awk '{print $1, $4, $5}' OFS="\t" > $OUTPUT_DIRECTORY/circos/genes_rev_strand.txt
sort -k1,1 -k2,2n $OUTPUT_DIRECTORY/circos/genes_rev_strand.txt > $OUTPUT_DIRECTORY/circos/genes_rev_strand.txt.tmp
mv $OUTPUT_DIRECTORY/circos/genes_rev_strand.txt.tmp $OUTPUT_DIRECTORY/circos/genes_rev_strand.txt

################################
### PREP INSERTION SITE DATA ###
################################

# convert the tradis insertion site output to bed format
zcat $OUTPUT_DIRECTORY/biotradis/combined.insert_site_plot.gz | awk '
BEGIN {
    while ((getline < "'$REFERENCE_GENOME.fai'") > 0) {
        chr[++n] = $1
        len[n] = $2
    }
    c = 1
    pos = 0
}
{
    print chr[c], pos, pos+1, $1 >> "'$OUTPUT_DIRECTORY/circos/insertions_fwd_strand.bed'"
    print chr[c], pos, pos+1, $2 >> "'$OUTPUT_DIRECTORY/circos/insertions_rev_strand.bed'"
    pos++
    if (pos >= len[c]) {
        c++
        pos = 0
    }
}
' OFS='\t'

# make genome windows to count insertion sites in
bedtools makewindows -g $REFERENCE_GENOME.fai -w 5000 > $OUTPUT_DIRECTORY/circos/windows_5kb.bed

# count the insertions per window
bedtools map -a $OUTPUT_DIRECTORY/circos/windows_5kb.bed \
-b $OUTPUT_DIRECTORY/circos/insertions_fwd_strand.bed \
-c 4 -o sum -null 0 > $OUTPUT_DIRECTORY/circos/insertions_fwd_strand_5kb.bed
bedtools map -a $OUTPUT_DIRECTORY/circos/windows_5kb.bed \
-b $OUTPUT_DIRECTORY/circos/insertions_rev_strand.bed \
-c 4 -o sum -null 0 > $OUTPUT_DIRECTORY/circos/insertions_rev_strand_5kb.bed

##########################
### PREP THE ORIC FILE ###
##########################

awk '$3=="oriC" {
    midpoint=int(($4+$5)/2)
    print $1, midpoint, midpoint, 1
}' OFS="\t" ${GENOME_ANNOTATION%.*}.gff > $OUTPUT_DIRECTORY/circos/oriC.txt

awk '$3=="oriC" {
    midpoint=int(($4+$5)/2)
    print $1, midpoint, midpoint, "oriC"
}' OFS="\t" ${GENOME_ANNOTATION%.*}.gff > $OUTPUT_DIRECTORY/circos/oriC_label.txt

#####################################
### WRITE THE CIRCOS CONFIG FILES ###
#####################################

echo "
<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
<<include ticks.conf>>

<image>
<<include etc/image.conf>>
#file  = circos.png
#radius = 1500p
</image>

karyotype = karyotype.txt

########################################

<ideogram>

show_label       = yes
label_font       = default
label_radius     = 1r + 110p
#label_radius     = 1.17r
label_size       = 40
label_parallel   = yes

<spacing>
default = 0.005r
</spacing>

radius*    = 0.85r
thickness = 20p
fill      = yes
color = black
</ideogram>

########################################

<plots>

########################
# gene annotation ring
########################

<plot>
type = tile
file = genes_fwd_strand.txt
r1   = 0.97r
r0   = 0.92r
color = dgrey
layers = 1
margin      = 0.05u
orientation = center
stroke_thickness = 1
stroke_color     = dgrey
thickness = 50
padding = 8
</plot>

<plot>
type = tile
file = genes_rev_strand.txt
r1   = 0.90r
r0   = 0.85r
color = dgrey
layers = 1
margin = 0.05u
orientation = center
stroke_thickness = 1
stroke_color = dgrey
thickness = 50
padding = 8
</plot>

############################
### INSERTION SITE PLOTS ###
############################

<plot>
type = histogram
file = insertions_fwd_strand_5kb.bed
r1   = 0.80r
r0   = 0.55r
color = vdred
fill_color = vdred
thickness = 0.5
</plot>

<plot>
type = histogram
file = insertions_rev_strand_5kb.bed
r1   = 0.55r
r0   = 0.30r
orientation = in
color = orange
fill_color = orange
thickness = 0.5
</plot>

##################
### oriC marker ###
##################

<plot>
type = scatter
file = oriC.txt
r1 = 0.83r
r0 = 0.81r
color = dblue
stroke_color = blue
stroke_thickness = 2
glyph = circle
glyph_size = 22
</plot>

##################
### oriC label ###
##################

<plot>
type = text
file = oriC_label.txt
r0 = 1.05r
r1 = 1.15r
label_size = 35p
label_font = bold
color = dblue
orientation = out
</plot>

</plots>
" > $OUTPUT_DIRECTORY/circos/circos.conf

echo "
show_ticks          = yes
show_tick_labels    = yes

<ticks>
skip_first_label = no
skip_last_label  = no
radius           = dims(ideogram,radius_outer)
tick_separation  = 3p
label_separation = 1p
multiplier       = 1e-6
color            = black
thickness        = 4p
size             = 20p

<tick>
spacing        = 10000
show_label     = no
thickness      = 3p
color          = dgrey
</tick>

<tick>
spacing        = 500000
show_label     = yes
label_size     = 20p
label_offset   = 10p
suffix         = Mb
format         = %.1f
thickness      = 5p
color          = vdgrey
</tick>

<tick>
spacing        = 1000000
show_label     = yes
label_size     = 26p
label_offset   = 10p
suffix         = Mb
format         = %d
grid           = yes
grid_color     = dgrey
grid_thickness = 1p
grid_start     = 0.5r
grid_end       = 0.999r
</tick>

</ticks>
" > $OUTPUT_DIRECTORY/circos/ticks.conf

##################
### RUN CIRCOS ###
##################

echo
echo "Running circos to generate the circos plot with all contigs in the same plot..."
echo

cd $OUTPUT_DIRECTORY/circos
circos

##################################################################
# Generate input files for making individual contig circos plots #
##################################################################

# set the window size in kb for plotting insertion sites
window_size=2

# make a list of the contigs in the genome
CONTIG_LIST=( $(cut -f1 $REFERENCE_GENOME.fai) )

##########################
# create karyotype files #
##########################

for CONTIG in ${CONTIG_LIST[@]}
do
grep "$CONTIG" $REFERENCE_GENOME.fai |  awk '{print "chr - " $1 " " $1 " 0 " $2 " chr"}' > karyotype_$CONTIG.txt
done

###########################
# create annotation files #
###########################

for CONTIG in ${CONTIG_LIST[@]}
do
# convert annotation to circos format
awk 'match($1, "'$CONTIG'") && (tolower($3) == "cds" || tolower($3) == "gene")' $GENOME_ANNOTATION |
awk '{print $1, $4, $5}' OFS="\t" > genes_$CONTIG.txt
sort -k1,1 -k2,2n genes_$CONTIG.txt > genes_$CONTIG.txt.tmp && mv genes_$CONTIG.txt.tmp genes_$CONTIG.txt

# make separate annotation files for genes on fwd and rev strands (strand info is 7th column)
# fwd strand
awk 'match($1, "'$CONTIG'") && (tolower($3) == "cds" || tolower($3) == "gene") && $7=="+"' $GENOME_ANNOTATION |
awk '{print $1, $4, $5}' OFS="\t" > genes_fwd_strand_$CONTIG.txt
sort -k1,1 -k2,2n genes_fwd_strand_$CONTIG.txt > genes_fwd_strand_$CONTIG.txt.tmp
mv genes_fwd_strand_$CONTIG.txt.tmp genes_fwd_strand_$CONTIG.txt
# rev strand
awk 'match($1, "'$CONTIG'") && (tolower($3) == "cds" || tolower($3) == "gene") && $7=="-"' $GENOME_ANNOTATION |
awk '{print $1, $4, $5}' OFS="\t" > genes_rev_strand_$CONTIG.txt
sort -k1,1 -k2,2n genes_rev_strand_$CONTIG.txt > genes_rev_strand_$CONTIG.txt.tmp
mv genes_rev_strand_$CONTIG.txt.tmp genes_rev_strand_$CONTIG.txt
done

# prep the insertion site files
for CONTIG in ${CONTIG_LIST[@]}; do
# reformat the tradis insertion site data to have the contig name and location at the start
grep "$CONTIG" $REFERENCE_GENOME.fai > $CONTIG.info.txt
zcat $OUTPUT_DIRECTORY/biotradis/*.$CONTIG.insert_site_plot.gz |
awk -v OFS='\t' '
NR == FNR {
    chr = $1
    pos = 0
    next
}
{
    print chr, pos, pos+1, $1 >> "insertions_fwd_strand_'$CONTIG'.bed"
    print chr, pos, pos+1, $2 >> "insertions_rev_strand_'$CONTIG'.bed"
    pos++
}
' $CONTIG.info.txt -
# make genome windows to count insertion sites in
bedtools makewindows -g $CONTIG.info.txt -w ${window_size}000 > windows_${window_size}kb_$CONTIG.bed
# count the insertions per window
bedtools map -a windows_${window_size}kb_$CONTIG.bed -b insertions_fwd_strand_$CONTIG.bed -c 4 -o sum -null 0 > insertions_fwd_strand_${window_size}kb_$CONTIG.bed
bedtools map -a windows_${window_size}kb_$CONTIG.bed -b insertions_rev_strand_$CONTIG.bed -c 4 -o sum -null 0 > insertions_rev_strand_${window_size}kb_$CONTIG.bed
# cleanup
rm $CONTIG.info.txt
done

##########################
### PREP THE ORIC FILE ###
##########################

for CONTIG in ${CONTIG_LIST[@]}; do
awk 'match($1, "'$CONTIG'") && $3=="oriC" {
    midpoint=int(($4+$5)/2)
    print $1, midpoint, midpoint, 1
}' OFS="\t" "$GENOME_ANNOTATION" > oriC_$CONTIG.txt

awk 'match($1, "'$CONTIG'") && $3=="oriC" {
    midpoint=int(($4+$5)/2)
    print $1, midpoint, midpoint, "oriC"
}' OFS="\t" "$GENOME_ANNOTATION" > oriC_label_$CONTIG.txt
done

##############################
# prep the contig label file #
##############################

for CONTIG in ${CONTIG_LIST[@]}; do
cut -f1,2 $REFERENCE_GENOME.fai | grep "$CONTIG" > tmp_$CONTIG
awk -v OFS='\t' -v contig="$CONTIG" '{print $1, int($2/4) -1, int($2/4) +1, contig}' tmp_$CONTIG > ${CONTIG}_label.txt
rm tmp_$CONTIG
done

#############################
# PREP THE CIRCOS CONF FILE #
#############################

for CONTIG in ${CONTIG_LIST[@]}; do
echo "
<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
<<include ticks.conf>>

<image>
background = white
dir   = .
#dir  = conf(configdir)
file  = circos_$CONTIG.png
png   = yes
svg   = yes

# radius of inscribed circle in image
radius         = 1500p
# by default angle=0 is at 3 o'clock position
angle_offset      = -90
#angle_orientation = counterclockwise
auto_alpha_colors = yes
auto_alpha_steps  = 5

</image>

karyotype = karyotype_$CONTIG.txt

########################################

<ideogram>

show_label       = yes
label_font       = default
#label_radius     = 1r + 110p
label_radius     = 1.17r
label_size       = 40
label_parallel   = yes

<spacing>
default = 0.005r
</spacing>

radius*    = 0.85r
thickness = 20p
fill      = yes
color = black
</ideogram>

########################################

<plots>

##############################
# contig label in the centre #
##############################

<plot>
type = text
file = ${CONTIG}_label.txt
r1 = 0.30r
r0 = 0.01r
label_size = 60
label_font = bold
color = black
horizontal_align = left
vertical_align = middle
rpadding = -0.5r
</plot>

########################
# gene annotation ring
########################

<plot>
type = tile
file = genes_fwd_strand_$CONTIG.txt
r1   = 0.97r
r0   = 0.92r
color = dgrey
layers = 1
margin      = 0.05u
orientation = center
stroke_thickness = 1
stroke_color     = dgrey
thickness = 50
padding = 8
</plot>

<plot>
type = tile
file = genes_rev_strand_$CONTIG.txt
r1   = 0.90r
r0   = 0.85r
color = dgrey
layers = 1
margin = 0.05u
orientation = center
stroke_thickness = 1
stroke_color = dgrey
thickness = 50
padding = 8
</plot>

############################
### INSERTION SITE PLOTS ###
############################

<plot>
type = histogram
file = insertions_fwd_strand_${window_size}kb_$CONTIG.bed
r1   = 0.80r
r0   = 0.55r
color = vdred
fill_color = vdred
thickness = 0.5
</plot>

<plot>
type = histogram
file = insertions_rev_strand_${window_size}kb_$CONTIG.bed
r1   = 0.55r
r0   = 0.30r
orientation = in
color = orange
fill_color = orange
thickness = 0.5
</plot>

##################
### oriC marker ###
##################

<plot>
type = scatter
file = oriC_$CONTIG.txt
r1 = 0.83r
r0 = 0.81r
color = dblue
stroke_color = blue
stroke_thickness = 2
glyph = circle
glyph_size = 22
</plot>

##################
### oriC label ###
##################

<plot>
type = text
file = oriC_label_$CONTIG.txt
r0 = 1.05r
r1 = 1.15r
label_size = 35p
label_font = bold
color = dblue
orientation = out
</plot>

</plots>
" > circos_$CONTIG.conf

#############################################
### RUN CIRCOS FOR EACH INDIVIDUAL CONTIG ###
#############################################

circos -conf circos_$CONTIG.conf
done

###############################################################################
# Cleanup environment
###############################################################################

module unload circos-uoneasy/0.69-9-GCCcore-11.3.0

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

