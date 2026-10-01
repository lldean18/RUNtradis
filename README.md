# RUNtradis
Executable for read processing and optimising &amp; running tradis software from a single simple script.

## USEAGE:

Usage:
    sbatch /share/bryant_lab/RUNtradis.sh --input INPUT_FASTQ --output OUTPUT_DIRECTORY [OPTIONS]

Required arguments:
    -i, --input       input FASTQ file (must be compressed with gzip)
    -o, --output      Output directory name (must not already exist)

Optional arguments:
    -r, --reference   Reference genome in FASTA format
                      Default: /share/bryant_lab/reference_genomes/GCF_000750555.1_ASM75055v1_genomic.fna

    -a, --annotation  Genome annotation file (must be in .embl format)
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





