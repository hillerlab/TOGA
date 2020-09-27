#!/usr/bin/env python3
"""Merge and process CESAR output files.

After the CESAR part there would appear numerous files.
This script applies parse_cesar_bdb function to each file
and then does the following:
1) Bed annotation track for query.
2) Nucleotide and protein fasta files.
3) Saves exons metadata into a tsv file.
4) Saves a list of problematic projections.
"""
import sys
import argparse
import os
from collections import defaultdict
try:
    from modules.parse_cesar_output import classify_exon
    from modules.common import eprint
    from modules.common import die
except ImportError:
    from parse_cesar_output import classify_exon
    from common import eprint
    from common import die

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]


# constants
MAX_SCORE = 1000
MAX_COLOR = 255
PID_HQ_THR = 65
BLOSUM_HQ_THR = 35
Q_HEADER_FIELDS_NUM = 12
BLACK = "0,0,0"
DEFAULT_SCORE = 1000

# header for exons meta data file
META_HEADER = "\t".join("gene exon_num chain_id act_region exp_region"
                        " in_exp pid blosum gap class paralog q_mark".split())


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("input_dir", help="Directory containing output BDB files")
    app.add_argument("output_bed", help="Save pre_final bed12 file to...")
    app.add_argument("output_fasta", help="Save fasta fasta to...")
    app.add_argument("meta_data", help="Save exons metadata to...")
    app.add_argument("prot_fasta", help="Save protein fasta to...")
    app.add_argument("codon_fasta", help="Save codon alignment fasta to...")
    app.add_argument("skipped", help="Save skipped genes")
    app.add_argument("--output_trash", default=None)
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    # print help if there are no args
    return args


def read_fasta(fasta_line, v=False):
    """Read fasta, return dict and type."""
    fasta_data = fasta_line.split(">")
    eprint(f"fasta_data[0] is:\n{fasta_data[0]}") if v else None
    eprint(f"fasta_data[1] is:\n{fasta_data[1]}") if v else None
    if fasta_data[0] != "":
        # this is a bug
        eprint("ERROR! Cesar output is corrupted")
        eprint(f"Issue detected in the following string:\n{fasta_line}")
        die("Abort")
    del fasta_data[0]  # remove it "" we don't need that
    sequences = {}  # accumulate data here
    order = []  # to have ordered list

    # there is no guarantee that dict will contain elements in the
    # same order as they were added
    for elem in fasta_data:
        raw_lines = elem.split("\n")
        # it must be first ['capHir1', 'ATGCCGCGCCAATTCCCCAAGCTGA... ]
        header = raw_lines[0]
        # separate nucleotide-containing lines
        lines = [x for x in raw_lines[1:] if x != "" and not x.startswith("!")]
        if len(lines) == 0:  # it is a mistake - empty sequence --> get rid of
            continue
        fasta_content = "".join(lines)
        sequences[header] = fasta_content
        order.append(header)
    return sequences, order


def read_region(region):
    """Return convenient region representation."""
    chrom, grange = region.split(":")
    start = int(grange.split("-")[0])
    end = int(grange.split("-")[1])
    return {"chrom": chrom, "start": start, "end": end}


def parse_cesar_bdb(arg_input, v=False):
    """Parse CESAR bdb file core function."""
    in_ = open(arg_input, "r")  # read cesar bdb file
    # two \n\n divide each unit of information
    content = [x for x in in_.read().split("\n\n") if x]
    in_.close()
    # GLP-related data is already filtered out by cesar_runner

    # initiate collectors
    bed_lines = []  # save bed lines here
    skipped = []  # save skipper projections here
    pred_seq_chain = {}  # for nucleotide sequences to fasta
    t_exon_seqs = defaultdict(dict)  # reference exon sequences
    wrong_exons = []  # exons that are predicted but actually deleted/missing
    all_meta_data = [META_HEADER]  # to collect exons meta data
    prot_data = []  # protein sequences
    codon_data = []  # codon sequences

    for elem in content:
        # one elem - one CESAR call (one ref transcript and >=1 chains)
        # now loop gene-by-gene
        gene = elem.split("\n")[0][1:]
        eprint(f"Reading gene {gene}") if v else None
        cesar_out = "\n".join(elem.split("\n")[1:])
        # basically this is a fasta file with headers
        # saturated with different information
        sequences, order = read_fasta(cesar_out, v=v)
        # initiate dicts to fill later
        ranges_chain, chain_dir = defaultdict(dict), {}
        pred_seq_chain[gene] = defaultdict(dict)

        # split fasta headers in different classes
        # query, ref and prot sequence headers are explicitly marked
        query_headers = [h for h in order if h.endswith("query_exon")]
        ref_headers = [h for h in order if h.endswith("reference_exon")]
        prot_ids = [h for h in order if "PROT" in h]
        codon_ids = [h for h in order if "CODON" in h]

        # parse reference exons, quite simple
        for header in ref_headers:
            # one header for one exon
            # fields look like this:
            # FIELD_1 | FIELD_2 | FIELD_3\n
            header_fields = [s.replace(" ", "") for s in header.split("|")]
            exon_num = int(header_fields[1])  # 0-based!
            exon_seq = sequences[header].replace("-", "")  # header is also a key for seq dict
            t_exon_seqs[gene][exon_num] = exon_seq

        # save protein data
        for prot_id in prot_ids:
            prot_seq = sequences[prot_id]
            prot_line = f">{prot_id}\n{prot_seq}\n"
            prot_data.append(prot_line)

        # save codon alignment data
        for codon_id in codon_ids:
            codon_seq = sequences[codon_id]
            codon_line_line = f">{codon_id}\n{codon_seq}\n"
            codon_data.append(codon_line_line)

        # get gene: exons dict to trace deleted exons
        gene_chain_exon_status = defaultdict(dict)

        # parse query headers
        for header in query_headers:
            header_fields = [s.replace(" ", "") for s in header.split("|")]
            if len(header_fields) != Q_HEADER_FIELDS_NUM:
                continue  # ref exon?

            # extract metadata, parse query header
            trans = header_fields[0]
            exon_num = int(header_fields[1])
            chain_id = int(header_fields[2])
            exon_region = read_region(header_fields[3])
            pid = float(header_fields[4])  # nucleotide %ID
            blosum = float(header_fields[5])
            is_gap = header_fields[6]  # asm gap in the expected region
            exon_class = header_fields[7]  # how it aligns to chain
            exp_region_str = header_fields[8]  # expected region
            in_exp = header_fields[9]  # detected in the expected region or not
            in_exp_b = True if in_exp == "INC" else False

            # mark that it's paralogous projection:
            para_annot = True if header_fields[10] == "True" else False
            stat_key = (trans, chain_id)  # projection ID
            # classify exon, check whether it's deleted/missing
            exon_decision, q_mark = classify_exon(exon_class, in_exp_b, pid, blosum)

            if exon_decision is False:
                # exon is deleted/missing
                wrong_exons.append(header)  # save this data
                gene_chain_exon_status[stat_key][exon_num] = False
            else:  # exon is not deleted
                # get/write necessary info
                gene_chain_exon_status[stat_key][exon_num] = True
                chain_dir[chain_id] = exon_region["end"] > exon_region["start"]
                ranges_chain[chain_id][exon_num] = exon_region
                pred_seq_chain[gene][chain_id][exon_num] = sequences[header]
            # collect exon meta-data -> write to file later
            meta_data = "\t".join([gene, header_fields[1], header_fields[2],
                                   header_fields[3], exp_region_str,
                                   in_exp, header_fields[4], header_fields[5], is_gap,
                                   exon_class, str(para_annot), q_mark])
            all_meta_data.append(meta_data)

        # check if there are any exons
        for name, stat in gene_chain_exon_status.items():
            any_exons_left = any(stat.values())
            if any_exons_left:
                continue
            # projection has no exons: log it
            name_ = f"{name[0]}.{name[1]}"
            skipped.append(f"{name_}\tall exons are deleted.")

        # make bed tracks
        for chain_id in chain_dir.keys():
            # go projection-by-projection: fixed gene, loop over chains
            block_starts = []
            block_sizes = []
            ranges = ranges_chain[chain_id]
            name = f"{gene}.{chain_id}"  # projection name for bed file

            if len(ranges) == 0:  # this projection is completely missing
                skipped.append(f"{name}\tall exons are deleted.")
                continue
            direct = chain_dir[chain_id]
            exon_nums = sorted(ranges.keys()) if direct else sorted(ranges.keys(), reverse=True)

            # get basic coordinates
            chrom = ranges[exon_nums[0]]["chrom"]
            chrom_start = ranges[exon_nums[0]]["start"] if direct else ranges[exon_nums[0]]["end"]
            chrom_end = ranges[exon_nums[-1]]["end"] if direct else ranges[exon_nums[-1]]["start"]

            # we do not predict UTRs: thickStart/End = chrom_start/End
            thickStart = chrom_start
            thick_end = chrom_end
            strand = "+" if direct else "-"
            block_count = len(exon_nums)

            # need to convert to "block starts" \ "block sizes" format
            for exon_num in exon_nums:
                ex_range = ranges[exon_num]
                block_sizes.append(abs(ex_range["end"] - ex_range["start"]))
                blockStart = ex_range["start"] - chrom_start if direct else ex_range["end"] - chrom_start
                block_starts.append(blockStart)

            # need this as strings to save it in a text file
            block_starts_str = ",".join(map(str, block_starts)) + ","
            block_sizes_str = ",".join(map(str, block_sizes)) + ","

            # join in a bed line
            bed_list = map(str, [chrom, chrom_start, chrom_end, name,
                                 DEFAULT_SCORE, strand, thickStart, thick_end, BLACK,
                                 block_count, block_sizes_str, block_starts_str])
            bed_line = "\t".join(bed_list)
            bed_lines.append(bed_line)

    # arrange fasta content
    fasta_lines_lst = []
    for gene, chain_exon_seq in pred_seq_chain.items():
        # write target gene info
        t_gene_seq_dct = t_exon_seqs.get(gene)
        if t_gene_seq_dct is None:
            # no sequence data for this transcript?
            eprint(f"Warning! Missing data for {gene}")
            skipped.append(f"{gene}\tmissing data after cesar stage")
            continue
        # We have sequence fragments split between different exons
        t_exon_nums = sorted(t_gene_seq_dct.keys())
        t_header = ">ref_{0}\n".format(gene)
        t_seq = "".join([t_gene_seq_dct[num] for num in t_exon_nums]) + "\n"
        # append data to fasta strings
        fasta_lines_lst.append(t_header)
        fasta_lines_lst.append(t_seq)

        # and query info
        for chain_id, exon_seq in chain_exon_seq.items():
            track_header = ">{0}.{1}\n".format(gene, chain_id)
            exon_nums = sorted(exon_seq.keys())
            # also need to assemble different exon sequences
            seq = "".join([exon_seq[num] for num in exon_nums]) + "\n"
            fasta_lines_lst.append(track_header)
            fasta_lines_lst.append(seq)

    # save corrupted exons as bed-6 track
    # to make it possible to save them and visualize in the browser
    trash_exons = []
    for elem in wrong_exons:
        elem_fields = [s.replace(" ", "") for s in elem.split("|")]
        # need to fill the following:
        # chrom, start, end, name, score, strand
        gene_name = elem_fields[0]
        exon_num = elem_fields[1]
        chain_id = elem_fields[2]
        label = ".".join([gene_name, exon_num, chain_id])
        grange = elem_fields[3].split(":")
        chrom, (start, end) = grange[0], grange[1].split("-")
        strand = "+"
        score = str(int(float(elem_fields[4]) * 10))
        bed_6 = "\t".join([chrom, start, end, label, score, strand]) + "\n"
        trash_exons.append(bed_6)

    # join output strings
    meta_str = "\n".join(all_meta_data) + "\n"
    skipped_str = "\n".join(skipped) + "\n"
    prot_fasta = "".join(prot_data)
    codon_fasta = "".join(codon_data)
    fasta_lines = "".join(fasta_lines_lst)
    return bed_lines, trash_exons, fasta_lines, meta_str, prot_fasta, codon_fasta, skipped_str


def merge_cesar_output(input_dir, output_bed, output_fasta,
                       meta_data_arg, skipped_arg, prot_arg,
                       codon_arg, output_trash):
    """Merge multiple CESAR output files."""
    # check that input dir is correct
    die(f"Error! {input_dir} is not a dir!") \
        if not os.path.isdir(input_dir) else None
    # get list of bdb files (output of CESAR part)
    bdbs = [x for x in os.listdir(input_dir) if x.endswith(".txt")]

    # initiate lists for different types of output:
    bed_summary = []
    fasta_summary = []
    trash_summary = []
    meta_summary = []
    prot_summary = []
    codon_summary = []
    skipped = []
    all_ok = True

    task_size = len(bdbs)
    # extract data for all the files
    for num, bdb_file in enumerate(bdbs):
        # parse bdb files one by one
        bdb_path = os.path.join(input_dir, bdb_file)
        try:  # try to parse data
            parsed_data = parse_cesar_bdb(bdb_path)
        except AssertionError:
            # if this happened: some assertion was violated
            # probably CESAR output data is corrupted
            sys.exit(f"Error! Failed reading file {bdb_file}")

        # unpack parsed data tuple:
        bed_lines = parsed_data[0]
        trash_exons = parsed_data[1]
        fasta_lines = parsed_data[2]
        meta_data = parsed_data[3]
        prot_fasta = parsed_data[4]
        codon_fasta = parsed_data[5]
        skip = parsed_data[6]

        if len(bed_lines) == 0:
            # actually should not happen, but can
            eprint(f"Warning! {bdb_file} is empty")
            all_ok = False
            continue  # it is empty
        
        # append data to lists
        bed_summary.append("\n".join(bed_lines) + "\n")
        fasta_summary.append(fasta_lines)
        trash_summary.append("".join(trash_exons))
        meta_summary.append(meta_data)
        skipped.append(skip)
        prot_summary.append(prot_fasta)
        codon_summary.append(codon_fasta)
        eprint(f"Reading file {num + 1}/{task_size}", end="\r")

    # save output
    eprint("Saving the output")
    if len(bed_summary) == 0:
        # if so, no need to continue
        eprint("! merge_cesar_output.py:")
        die("No projections found! Abort.")

    # save bed, fasta and the rest
    with open(output_bed, "w") as f:
        f.write("".join(bed_summary))
    with open(output_fasta, "w") as f:
        f.write("".join(fasta_summary))
    with open(meta_data_arg, "w") as f:
        f.write("\n".join(meta_summary))
    with open(skipped_arg, "w") as f:
        f.write("\n".join(skipped))
    with open(prot_arg, "w") as f:
        f.write("\n".join(prot_summary))
    with open(codon_arg, "w") as f:
        f.write("\n".join(codon_summary))

    if output_trash:
        # if requested: provide trash annotation
        f = open(output_trash, "w")
        f.write("".join(trash_summary))
        f.close()
    return all_ok


def main():
    """Entry point."""
    args = parse_args()
    merge_cesar_output(args.input_dir,
                       args.output_bed,
                       args.output_fasta,
                       args.meta_data,
                       args.skipped,
                       args.prot_fasta,
                       args.codon_fasta,
                       args.output_trash)


if __name__ == "__main__":
    main()
