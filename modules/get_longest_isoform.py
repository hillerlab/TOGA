#!/usr/bin/env python3
"""Get list of longest isoforms."""
import sys
from collections import defaultdict

# may be variable!
GENE_INDEX_ISOF = 0
TRANSCRIPT_INDEX_IFOS = 1
NAME_IND_BED = 3
SIZES_IND_BED = 10


def read_isoforms(isoforms_file):
    """Create gene to transcripts dictionary."""
    gene_to_trans = defaultdict(list)
    f = open(isoforms_file, "r")
    for line in f:
        line_data = line.rstrip().split("\t")
        gene = line_data[GENE_INDEX_ISOF]
        trans = line_data[TRANSCRIPT_INDEX_IFOS]
        gene_to_trans[gene].append(trans)
    f.close()
    return gene_to_trans


def read_bed_file(bed_file):
    """Read bed file."""
    trans_to_data = {}
    trans_to_line = {}
    f = open(bed_file, "r")
    for line in f:
        line_plain = line.rstrip()
        line_data = line_plain.split("\t")
        trans_name = line_data[NAME_IND_BED]
        trans_block_sizes_raw = line_data[SIZES_IND_BED]
        trans_block_sizes = [int(x) for x in trans_block_sizes_raw.split(",")
                             if x != ""]
        trans_block_sum = sum(trans_block_sizes)
        trans_to_data[trans_name] = (trans_block_sum,
                                     trans_name)
        trans_to_line[trans_name] = line
    f.close()
    return trans_to_data, trans_to_line


def main(bed_file, isoforms_file, output_file, errors_file):
    """Perform longest isoform extraction."""
    # read bed and isoforms
    gene_to_transcripts = read_isoforms(isoforms_file)
    trans_to_data, trans_to_line = read_bed_file(bed_file)
    errors = []
    output = []

    # get the longest for each gene
    # filter and keep errors
    for gene, transcripts in gene_to_transcripts.items():
        transcripts_data = []
        for trans in transcripts:
            trans_data = trans_to_data.get(trans)
            if trans_data:
                transcripts_data.append(trans_data)
            else:
                errors.append(trans)
        if len(transcripts_data) == 0:
            # all transcripts don't match
            continue
        tr_sorted = sorted(transcripts_data, key=lambda x: x[0])
        longest_name = tr_sorted[-1][1]
        result = (gene, longest_name)
        output.append(result)

    # save output
    f = open(output_file, "w")
    for elem in output:
        # f.write(f"{elem[0]}\t{elem[1]}\n")
        trans = elem[1]
        line = trans_to_line[trans]
        f.write(line)
    f.close()

    f = open(errors_file, "w")
    f.write("\n".join(errors) + "\n")
    f.close()


if __name__ == "__main__":
    try:
        bed_file = sys.argv[1]
        isoforms_file = sys.argv[2]
        not_found = sys.argv[3]
        result = sys.argv[4]
    except IndexError:
        sys.exit(f"Usage: {sys.argv[0]} [BED_FILE] [ISOFORMS_FILE]"
                 f" [NOT_FOUND] [RESULT]")
    main(bed_file, isoforms_file, result, not_found)
