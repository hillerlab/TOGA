#!/usr/bin/env python3
"""Extract data from CESAR output."""
import argparse
from collections import defaultdict
import sys
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
DEFAULT_SCORE = 100

# header for exons meta data file
META_HEADER = "\t".join("gene exon_num chain_id act_region exp_region"
                        " in_exp pid blosum gap class paralog q_mark".split())


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("input", help="Bdb CESAR output or a fasta.")
    app.add_argument("--bed", help="Save bed file.", default=None)
    app.add_argument("--trash_exons", help="Save exons declined by a "
                                           "low pID into a separated bed",
                     default=None)
    app.add_argument("--fasta", help="Save fasta file if required.", default=None)
    app.add_argument("--prot_fasta", help="Save protein fasta", default=None)
    app.add_argument("--exons_metadata", help="Save exon data", default=None)
    app.add_argument("--skipped", help="Save skipped genes", default=None)
    app.add_argument("--exon_class", help="Produce exon/class track", default=None)
    app.add_argument("--verbose", action="store_true", dest="verbose", help="Show verbose messages")
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
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


def read_paralogs_log(paralog_log_file):
    """Get list of paralogous gene-chains."""
    f = open(paralog_log_file)
    gene_chains = [x for x in f.read().split("\n") if x != ""]
    f.close()
    return set(gene_chains)


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
        exon_lens, chain_raw_scores = {}, defaultdict(list)

        # split fasta headers in different classes
        # query, ref and prot sequence headers are explicitly marked
        query_headers = [h for h in order if h.endswith("query_exon")]
        ref_headers = [h for h in order if h.endswith("reference_exon")]
        prot_ids = [h for h in order if "PROT" in h]

        # parse reference exons, quite simple
        for header in ref_headers:
            # one header for one exon
            # fields look like this:
            # FIELD_1 | FIELD_2 | FIELD_3\n
            header_fields = [s.replace(" ", "") for s in header.split("|")]
            exon_num = int(header_fields[1])  # 0-based!
            exon_seq = sequences[header].replace("-", "")  # header is also a key for seq dict
            seq_len = len(exon_seq)
            t_exon_seqs[gene][exon_num] = exon_seq
            exon_lens[exon_num] = seq_len
        
        # save protein data
        for prot_id in prot_ids:
            prot_seq = sequences[prot_id]
            prot_line = f">{prot_id}\n{prot_seq}\n"
            prot_data.append(prot_line)

        # will need this for score computation (deprecated)
        gene_len = sum(exon_lens.values())
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

            try:  # time to get exon score (normalize to exon length)
                exon_score = int(pid / 100 * exon_lens[exon_num])
            except KeyError:
                pid, exon_score = 0, 0
            except ZeroDivisionError:
                # should never happen, because exon length cannot be 0
                # but just in case
                pid, exon_score = 0, 0
            
            if exon_decision is False:
                # exon is deleted/missing
                wrong_exons.append(header)  # save this data
                gene_chain_exon_status[stat_key][exon_num] = False
                chain_raw_scores[chain_id].append(0)
            else:  # exon is not deleted
                chain_raw_scores[chain_id].append(exon_score)
                # get/write necessary info
                gene_chain_exon_status[stat_key][exon_num] = True
                directed = exon_region["end"] > exon_region["start"]
                chain_dir[chain_id] = directed
                ranges_chain[chain_id][exon_num] = exon_region
                pred_seq_chain[gene][chain_id][exon_num] = sequences[header]
            # collect exon meta-data
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

        # make bed regions
        for chain_id in chain_dir.keys():
            # go projection-by-projection
            # in this chunk we operate with the same projection
            # change only chains
            chain_r_score = sum(chain_raw_scores[chain_id])
            chain_score = int(MAX_SCORE * (chain_r_score / gene_len))

            # assemble bed data
            blockStarts, blockSizes = [], []
            ranges = ranges_chain[chain_id]
            name = f"{gene}.{chain_id}"  # projection name for bed file
            # stat_key = (gene, chain_id)  # projection identifier for dicts

            if len(ranges) == 0:
                # this gene is completely missing
                skipped.append(f"{name}\tall exons are deleted.")
                continue
            direct = chain_dir[chain_id]
            exon_nums = sorted(ranges.keys()) if direct \
                else sorted(ranges.keys(), reverse=True)

            # get basic coordinates
            chrom = ranges[exon_nums[0]]["chrom"]
            chromStart = ranges[exon_nums[0]]["start"] if direct \
                else ranges[exon_nums[0]]["end"]
            chromEnd = ranges[exon_nums[-1]]["end"] if direct \
                else ranges[exon_nums[-1]]["start"]

            score = chain_score  # TODO: check whether we need that
            # we do not predict UTRs: thickStart/End = chromStart/End
            thickStart = chromStart
            thickEnd = chromEnd
            itemRgb = BLACK
            strand = "+" if direct else "-"
            blockCount = len(exon_nums)

            # need to convert to "block starts" \ "block sizes" format
            for exon_num in exon_nums:
                ex_range = ranges[exon_num]
                blockSizes.append(abs(ex_range["end"] - ex_range["start"]))
                blockStart = ex_range["start"] - chromStart if direct \
                    else ex_range["end"] - chromStart
                blockStarts.append(blockStart)

            # need this as strings to save it in a text file
            blockStarts_str = ",".join(map(str, blockStarts)) + ","
            blockSizes_str = ",".join(map(str, blockSizes)) + ","

            # join in a bed line
            bed_list = map(str, [chrom, chromStart, chromEnd,
                                 name, score, strand,
                                 thickStart, thickEnd, itemRgb,
                                 blockCount, blockSizes_str, blockStarts_str])
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
    fasta_lines = "".join(fasta_lines_lst)
    return bed_lines, trash_exons, fasta_lines, meta_str, prot_fasta, skipped_str


def main():
    """Entry point as stand-alone tool."""
    args = parse_args()
    # basically call the function that extracts everything
    parsed_data = parse_cesar_bdb(args.input, v=args.verbose)

    bed_lines = parsed_data[0]
    trash_exons = parsed_data[1]
    fasta_lines = parsed_data[2]
    exons_meta = parsed_data[3]
    prot_fasta = parsed_data[4]
    skipped = parsed_data[5]

    # and save this (as a standalone script)
    if args.bed:
        f = open(args.bed, "w") if args.bed != "stdout" else sys.stdout
        f.write("\n".join(bed_lines) + "\n")
        f.close() if args.bed != "stdout" else None

    if args.trash_exons:  # save bed
        f = open(args.trash_exons, "w") if args.trash_exons != "stdout" else sys.stdout
        f.write("\n".join(trash_exons) + "\n")
        f.close() if args.trash_exons != "stdout" else None

    if args.fasta:  # save fasta
        f = open(args.fasta, "w") if args.fasta != "stdout" else sys.stdout
        f.write(fasta_lines)
        f.close() if args.fasta != "stdout" else None

    if args.exons_metadata:  # save metadata
        f = open(args.exons_metadata, "w") if args.exons_metadata != "stdout" else sys.stdout
        f.write(exons_meta)
        f.close() if args.exons_metadata != "stdout" else None

    if args.prot_fasta:  # save protein fasta
        f = open(args.prot_fasta, "w") if args.prot_fasta != "stdout" else sys.stdout
        f.write(prot_fasta)
        f.close() if args.prot_fasta != "stdout" else None

    if args.skipped:  # save data of skipped projections
        f = open(args.skipped, "w") if args.skipped != "stdout" else sys.stdout
        f.write(skipped)
        f.close()
    sys.exit(0)


if __name__ == "__main__":
    main()
