#!/usr/bin/env python3

"""
Calculate insertion-site counts per gene from genome coordinates.

Python:
    3.12.3

Input:
    1. Genome assembly FASTA
    2. Genome annotation in EMBL format
    3. Gene list 1 (one gene name per line)
    4. Gene list 2 (one gene name per line)
    5. Two-column numerical file

The numerical file must contain exactly one row for every base
in the genome, in the same chromosome/contig order as the FASTA.

For each gene, the two numerical columns are summed across all
positions belonging to that gene.

Output:
    - genome_coordinates.txt
    - gene_insertion_site_counts.txt
    - insertion_sites_per_gene.png

Usage:
    python insertion_sites_per_gene.py \
        assembly.fasta \
        annotation.embl \
        genes1.txt \
        genes2.txt \
        insertion_sites.txt
"""

import argparse
import os
import re
import sys
import gzip

import matplotlib.pyplot as plt
from Bio import SeqIO


def read_arguments():
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Calculate insertion-site counts per gene from genome coordinates."
    )

    parser.add_argument(
        "fasta",
        help="Genome assembly in FASTA format"
    )

    parser.add_argument(
        "embl",
        help="Genome annotation in EMBL format"
    )

    parser.add_argument(
        "genes1",
        help="First gene list: one gene name per line"
    )

    parser.add_argument(
        "genes2",
        help="Second gene list: one gene name per line"
    )

    parser.add_argument(
        "numbers",
        help="Two-column numerical file, one row per genome base"
    )

    return parser.parse_args()


def open_text_file(filename, mode="rt"):
    """
    Open a normal text file or a gzip-compressed text file.
    """

    if filename.endswith(".gz"):
        return gzip.open(filename, mode)

    return open(filename, mode)


def read_fasta_lengths(fasta_file):
    """
    Read FASTA and return:

        contigs = [
            (contig_name, length),
            ...
        ]

    Also calculate the total genome length.
    """

    contigs = []

    print("\nReading genome assembly...")

    for record in SeqIO.parse(fasta_file, "fasta"):
        contig_name = record.id
        length = len(record.seq)

        contigs.append((contig_name, length))

        print(f"  {contig_name}: {length:,} bp")

    if not contigs:
        raise ValueError("No FASTA sequences were found.")

    total_length = sum(length for _, length in contigs)

    print(f"\nNumber of contigs/chromosomes: {len(contigs):,}")
    print(f"Total genome length: {total_length:,} bp")

    return contigs, total_length


def create_genome_coordinates(contigs, output_file):
    """
    Create a file containing:

        chromosome    position

    with one row for every base in the genome.
    """

    print("\nCreating genome coordinate file...")

    total_rows = 0

    with open(output_file, "w") as out:

        out.write("chromosome\tposition\n")

        for contig_name, length in contigs:

            for position in range(1, length + 1):
                out.write(f"{contig_name}\t{position}\n")

                total_rows += 1

    print(f"Created: {output_file}")
    print(f"Rows written: {total_rows:,}")

    return total_rows


def count_numeric_rows(numbers_file):
    """
    Count rows in the numerical input file.

    Blank lines are ignored.
    """

    count = 0

    with open_text_file(numbers_file, "rt") as infile:

        for line_number, line in enumerate(infile, start=1):

            line = line.strip()

            if not line:
                continue

            fields = line.split()

            if len(fields) != 2:
                raise ValueError(
                    f"Numerical file line {line_number} does not contain "
                    f"exactly two columns:\n{line}"
                )

            try:
                float(fields[0])
                float(fields[1])
            except ValueError:
                raise ValueError(
                    f"Numerical file line {line_number} contains "
                    f"non-numerical values:\n{line}"
                )

            count += 1

    return count


def create_coordinate_number_file(
    contigs,
    numbers_file,
    output_file
):
    """
    Combine genome coordinates with the two numerical columns.

    Output:

        chromosome    position    number1    number2
    """

    print("\nCombining genome coordinates with numerical data...")

    contig_index = 0
    position = 1

    current_contig = contigs[contig_index][0]
    current_length = contigs[contig_index][1]

    rows_written = 0

    with open_text_file(numbers_file, "rt") as infile, open(output_file, "w") as outfile:

        outfile.write(
            "chromosome\tposition\tvalue1\tvalue2\n"
        )

        for line_number, line in enumerate(infile, start=1):

            line = line.strip()

            if not line:
                continue

            fields = line.split()

            if len(fields) != 2:
                raise ValueError(
                    f"Line {line_number} in numerical file does not "
                    f"contain exactly two columns."
                )

            value1 = float(fields[0])
            value2 = float(fields[1])

            outfile.write(
                f"{current_contig}\t"
                f"{position}\t"
                f"{value1}\t"
                f"{value2}\n"
            )

            rows_written += 1

            position += 1

            # Move to the next contig when the current one is finished.
            if position > current_length:

                contig_index += 1

                if contig_index < len(contigs):

                    current_contig = contigs[contig_index][0]
                    current_length = contigs[contig_index][1]
                    position = 1

    expected_rows = sum(length for _, length in contigs)

    if rows_written != expected_rows:

        raise ValueError(
            "\nThe number of rows in the numerical file does not match "
            "the total genome length.\n\n"
            f"Genome bases expected: {expected_rows:,}\n"
            f"Numerical rows found: {rows_written:,}\n"
            f"Difference: {rows_written - expected_rows:+,}\n"
        )

    print(f"Created: {output_file}")
    print(f"Rows written: {rows_written:,}")

    return rows_written


def read_gene_list(filename):
    """
    Read gene names from a text file.

    One gene name per line.

    Gene names are stored case-insensitively.
    """

    genes = []

    with open(filename, "r") as infile:

        for line in infile:

            gene = line.strip()

            if not gene:
                continue

            genes.append(gene)

    return genes


def get_gene_name_from_feature(feature):
    """
    Try to obtain a useful gene name from an EMBL feature.

    The function checks several common qualifiers.
    """

    qualifiers = feature.qualifiers

    possible_qualifiers = [
        "gene",
        "locus_tag",
        "gene_synonym",
        "Name",
        "name",
        "label",
        "protein_id",
    ]

    for qualifier in possible_qualifiers:

        if qualifier in qualifiers:

            values = qualifiers[qualifier]

            if values:
                return str(values[0])

    return None


def parse_embl_genes(embl_file, contigs):
    """
    Parse EMBL annotations.

    Returns a dictionary:

        gene_name_lowercase -> list of annotations

    Each annotation contains:

        chromosome
        start
        end
        strand
        original_name

    Coordinates are converted to 1-based inclusive coordinates.
    """

    print("\nReading EMBL annotation...")

    genes = {}

    # Get chromosome/contig names from the FASTA.
    fasta_contig_names = [name for name, length in contigs]

    # If there is only one FASTA contig, use its name for all
    # EMBL features. This handles EMBL files with generic IDs
    # such as "XXX".
    if len(fasta_contig_names) == 1:
        embl_chromosome = fasta_contig_names[0]
    else:
        embl_chromosome = None

    for record in SeqIO.parse(embl_file, "embl"):

        if embl_chromosome is not None:
            chromosome = embl_chromosome
        else:
            chromosome = record.id

        for feature in record.features:

            # Use gene features only.
            #
            # This avoids counting the same gene twice when the
            # EMBL file contains both a gene feature and a CDS
            # feature with the same gene name.
            if feature.type != "gene":
                continue

            gene_name = get_gene_name_from_feature(feature)

            if gene_name is None:
                continue

            # Biopython locations are zero-based internally.
            start = int(feature.location.start) + 1
            end = int(feature.location.end)

            if feature.location.strand == -1:
                strand = "-"
            elif feature.location.strand == 1:
                strand = "+"
            else:
                strand = "."

            key = gene_name.lower()

            annotation = {
                "chromosome": chromosome,
                "start": start,
                "end": end,
                "strand": strand,
                "original_name": gene_name,
            }

            if key not in genes:
                genes[key] = []

            genes[key].append(annotation)

    print(f"Gene names found in EMBL: {len(genes):,}")

    return genes


def calculate_gene_counts(
    gene_list,
    embl_genes,
    coordinate_number_file,
):
    """
    Calculate the summed two numerical columns for every gene.

    Returns a list of dictionaries.
    """

    print("\nCalculating gene insertion-site counts...")

    results = []

    # Create a lookup containing only genes we need.
    requested_genes = {
        gene.lower(): gene
        for gene in gene_list
    }

    # Store all matching genomic intervals.
    intervals = []

    not_found = []

    for gene_lower, original_gene_name in requested_genes.items():

        if gene_lower not in embl_genes:

            not_found.append(original_gene_name)
            continue

        for annotation in embl_genes[gene_lower]:

            intervals.append(
                {
                    "requested_name": original_gene_name,
                    "annotation_name": annotation["original_name"],
                    "chromosome": annotation["chromosome"],
                    "start": annotation["start"],
                    "end": annotation["end"],
                    "strand": annotation["strand"],
                    "value1": 0.0,
                    "value2": 0.0,
                }
            )

    if not_found:

        print(
            f"\nWARNING: {len(not_found)} genes from the list "
            f"were not found in the EMBL annotation."
        )

        for gene in not_found:
            print(f"  Not found: {gene}")

    # Group intervals by chromosome so that we only need to
    # examine relevant genes for each coordinate.
    intervals_by_chromosome = {}

    for interval in intervals:

        chromosome = interval["chromosome"]

        if chromosome not in intervals_by_chromosome:
            intervals_by_chromosome[chromosome] = []

        intervals_by_chromosome[chromosome].append(interval)

    # Read the coordinate/value file.
    with open(coordinate_number_file, "r") as infile:

        header = infile.readline()

        for line_number, line in enumerate(infile, start=2):

            line = line.strip()

            if not line:
                continue

            fields = line.split()

            chromosome = fields[0]
            position = int(fields[1])
            value1 = float(fields[2])
            value2 = float(fields[3])

            if chromosome not in intervals_by_chromosome:
                continue

            for interval in intervals_by_chromosome[chromosome]:

                if (
                    interval["start"]
                    <= position
                    <= interval["end"]
                ):

                    interval["value1"] += value1
                    interval["value2"] += value2

    # Convert intervals into final results.
    for interval in intervals:

        total = interval["value1"] + interval["value2"]

        results.append(
            {
                "gene": interval["requested_name"],
                "annotation_name": interval["annotation_name"],
                "chromosome": interval["chromosome"],
                "start": interval["start"],
                "end": interval["end"],
                "strand": interval["strand"],
                "value1_sum": interval["value1"],
                "value2_sum": interval["value2"],
                "total": total,
            }
        )

    return results


def write_results(results, output_file):
    """Write per-gene results."""

    with open(output_file, "w") as outfile:

        outfile.write(
            "gene\tannotation_name\tchromosome\tstart\tend\tstrand\t"
            "value1_sum\tvalue2_sum\ttotal\n"
        )

        for result in results:

            outfile.write(
                f"{result['gene']}\t"
                f"{result['annotation_name']}\t"
                f"{result['chromosome']}\t"
                f"{result['start']}\t"
                f"{result['end']}\t"
                f"{result['strand']}\t"
                f"{result['value1_sum']}\t"
                f"{result['value2_sum']}\t"
                f"{result['total']}\n"
            )

    print(f"\nResults written to: {output_file}")


## NEW VERSION OF THIS FUNCTION TO SPLIT PLOT IN 2
##  def make_plot(results1, results2, output_file):
##      """
##      Plot two separate histograms, one above the other.
##  
##      The two histograms have independent x-axis scales so that
##      distributions with very different ranges can be seen clearly.
##      """
##  
##      values1 = [result["total"] for result in results1]
##      values2 = [result["total"] for result in results2]
##  
##      if not values1:
##          raise ValueError(
##              "No genes with annotations were found in the essential genes list."
##          )
##  
##      if not values2:
##          raise ValueError(
##              "No genes with annotations were found in the non-essential genes list."
##          )
##  
##      # Create two vertically stacked plots with independent x axes.
##      fig, axes = plt.subplots(
##          2,
##          1,
##          figsize=(10, 10)
##      )
##  
##      # ---------------------------------------------------------
##      # Essential genes
##      # ---------------------------------------------------------
##  
##      axes[0].hist(
##          values1,
##          bins=200,
##          alpha=0.7
##      )
##  
##      axes[0].set_xlabel("Insertion sites per gene")
##      axes[0].set_ylabel("Number of genes")
##      axes[0].set_title("Essential genes")
##  
##      # ---------------------------------------------------------
##      # Non-essential genes
##      # ---------------------------------------------------------
##  
##      axes[1].hist(
##          values2,
##          bins=200,
##          alpha=0.7
##      )
##  
##      axes[1].set_xlabel("Insertion sites per gene")
##      axes[1].set_ylabel("Number of genes")
##      axes[1].set_title("Non-essential genes")
##  
##      fig.suptitle(
##          "Distribution of insertion sites per gene",
##          fontsize=16
##      )
##  
##      plt.tight_layout()
##  
##      # Leave room for the overall title.
##      plt.subplots_adjust(
##          top=0.92
##      )
##  
##      plt.savefig(
##          output_file,
##          dpi=300,
##          bbox_inches="tight"
##      )
##  
##      plt.close()
##  
##      print(f"Plot written to: {output_file}")

# newest version of function
def make_plot(results1, results2, output_file):
    """
    Plot the distribution of total insertion sites per gene
    for the two gene lists.
    """

    values1 = [result["total"] for result in results1]
    values2 = [result["total"] for result in results2]

    if not values1:
        raise ValueError(
            "No genes with annotations were found in gene list 1."
        )

    if not values2:
        raise ValueError(
            "No genes with annotations were found in gene list 2."
        )

    plt.figure(figsize=(10, 7))

    # Width of each histogram bin.
    bin_width = 100

    # Create separate bin ranges, but with the same bin width.
    min1 = min(values1)
    max1 = max(values1)

    min2 = min(values2)
    max2 = max(values2)

    bins1 = range(
        int(min1 // bin_width) * bin_width,
        int(max1 // bin_width + 2) * bin_width,
        bin_width
    )

    bins2 = range(
        int(min2 // bin_width) * bin_width,
        int(max2 // bin_width + 2) * bin_width,
        bin_width
    )

    plt.hist(
        values1,
        bins=bins1,
        alpha=0.7,
        label="essential genes",
    )

    plt.hist(
        values2,
        bins=bins2,
        alpha=0.7,
        label="non-essential genes",
    )

    plt.xlabel("Insertion sites per gene")
    plt.ylabel("Number of genes")
    plt.title("Distribution of insertion sites per gene")

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(f"Plot written to: {output_file}")

# OLD VERSION OF THIS FUNCTION
###  def make_plot(results1, results2, output_file):
###      """
###      Plot the distribution of total insertion sites per gene
###      for the two gene lists.
###      """
###  
###      values1 = [result["total"] for result in results1]
###      values2 = [result["total"] for result in results2]
###  
###      if not values1:
###          raise ValueError(
###              "No genes with annotations were found in gene list 1."
###          )
###  
###      if not values2:
###          raise ValueError(
###              "No genes with annotations were found in gene list 2."
###          )
###  
###      plt.figure(figsize=(10, 7))
###  
###      # Use automatically determined bins covering both distributions.
###      all_values = values1 + values2
###  
###      min_value = min(all_values)
###      max_value = max(all_values)
###  
###      if min_value == max_value:
###          bins = 100
###      else:
###          bins = 100
###  
###      plt.hist(
###          values1,
###          bins=bins,
###          alpha=0.7,
###          label="essential genes",
###      )
###  
###      plt.hist(
###          values2,
###          bins=bins,
###          alpha=0.7,
###          label="non-essential genes",
###      )
###  
###      plt.xlabel("Insertion sites per gene")
###      plt.ylabel("Number of genes")
###      plt.title("Distribution of insertion sites per gene")
###  
###      plt.legend()
###  
###      plt.tight_layout()
###  
###      plt.savefig(
###          output_file,
###          dpi=300,
###          bbox_inches="tight"
###      )
###  
###      plt.close()
###  
###      print(f"Plot written to: {output_file}")


def main():

    args = read_arguments()

    # ---------------------------------------------------------
    # Check input files
    # ---------------------------------------------------------

    input_files = [
        args.fasta,
        args.embl,
        args.genes1,
        args.genes2,
        args.numbers,
    ]

    for filename in input_files:

        if not os.path.isfile(filename):

            raise FileNotFoundError(
                f"Input file not found: {filename}"
            )

    # ---------------------------------------------------------
    # 1. Read FASTA
    # ---------------------------------------------------------

    contigs, genome_length = read_fasta_lengths(args.fasta)

    # ---------------------------------------------------------
    # 2. Create genome coordinate file
    # ---------------------------------------------------------

    coordinate_file = "genome_coordinates.txt"

    coordinate_rows = create_genome_coordinates(
        contigs,
        coordinate_file
    )

    # ---------------------------------------------------------
    # 3. Check numerical input file
    # ---------------------------------------------------------

    print("\nChecking numerical input file...")

    number_rows = count_numeric_rows(args.numbers)

    print(f"Genome coordinate rows: {coordinate_rows:,}")
    print(f"Numerical input rows:    {number_rows:,}")

    if coordinate_rows != number_rows:

        raise ValueError(
            "\nERROR: The number of rows in the two files does not match.\n"
            f"Genome coordinates: {coordinate_rows:,}\n"
            f"Numerical file:     {number_rows:,}\n"
        )

    print("Row counts match.")

    # ---------------------------------------------------------
    # 4. Merge coordinates and numerical values
    # ---------------------------------------------------------

    coordinate_number_file = "genome_coordinates_with_values.txt"

    create_coordinate_number_file(
        contigs,
        args.numbers,
        coordinate_number_file
    )

    # ---------------------------------------------------------
    # 5. Read gene lists
    # ---------------------------------------------------------

    genes1 = read_gene_list(args.genes1)
    genes2 = read_gene_list(args.genes2)

    print(f"\nGene list 1: {len(genes1):,} genes")
    print(f"Gene list 2: {len(genes2):,} genes")

    # ---------------------------------------------------------
    # 6. Parse EMBL
    # ---------------------------------------------------------

    embl_genes = parse_embl_genes(args.embl, contigs)

    # ---------------------------------------------------------
    # 7. Calculate counts for gene list 1
    # ---------------------------------------------------------

    print("\n========================================")
    print("Processing gene list 1")
    print("========================================")

    results1 = calculate_gene_counts(
        genes1,
        embl_genes,
        coordinate_number_file
    )

    print(
        f"Genes successfully matched: {len(results1):,}"
    )

    # ---------------------------------------------------------
    # 8. Calculate counts for gene list 2
    # ---------------------------------------------------------

    print("\n========================================")
    print("Processing gene list 2")
    print("========================================")

    results2 = calculate_gene_counts(
        genes2,
        embl_genes,
        coordinate_number_file
    )

    print(
        f"Genes successfully matched: {len(results2):,}"
    )

    # ---------------------------------------------------------
    # 9. Write results
    # ---------------------------------------------------------

    combined_results_file = "gene_insertion_site_counts.txt"

    # Add a group column to the output.
    with open(combined_results_file, "w") as outfile:

        outfile.write(
            "gene_group\tgene\tannotation_name\tchromosome\t"
            "start\tend\tstrand\tinsertions_fwd_strand\tinsertions_rev_strand\ttotal_insertions\n"
        )

        for result in results1:

            outfile.write(
                f"essential_genes\t"
                f"{result['gene']}\t"
                f"{result['annotation_name']}\t"
                f"{result['chromosome']}\t"
                f"{result['start']}\t"
                f"{result['end']}\t"
                f"{result['strand']}\t"
                f"{result['value1_sum']}\t"
                f"{result['value2_sum']}\t"
                f"{result['total']}\n"
            )

        for result in results2:

            outfile.write(
                f"non-essential_genes\t"
                f"{result['gene']}\t"
                f"{result['annotation_name']}\t"
                f"{result['chromosome']}\t"
                f"{result['start']}\t"
                f"{result['end']}\t"
                f"{result['strand']}\t"
                f"{result['value1_sum']}\t"
                f"{result['value2_sum']}\t"
                f"{result['total']}\n"
            )

    print(
        f"\nCombined results written to: "
        f"{combined_results_file}"
    )

    # ---------------------------------------------------------
    # 10. Plot distributions
    # ---------------------------------------------------------

    plot_file = "insertion_sites_per_gene.png"

    make_plot(
        results1,
        results2,
        plot_file
    )

    # ---------------------------------------------------------
    # Finished
    # ---------------------------------------------------------

    print("\n========================================")
    print("Finished")
    print("========================================")

    print(f"\nGenome length: {genome_length:,} bp")
    print(f"Gene list 1 matched: {len(results1):,}")
    print(f"Gene list 2 matched: {len(results2):,}")

    print("\nOutput files:")
    print(f"  {coordinate_file}")
    print(f"  {coordinate_number_file}")
    print(f"  {combined_results_file}")
    print(f"  {plot_file}")


if __name__ == "__main__":
    try:
        main()

    except Exception as error:

        print(
            f"\nERROR: {error}",
            file=sys.stderr
        )

        sys.exit(1)
