#!/usr/bin/env python3
"""Extract data from CESAR output."""
import argparse
from collections import defaultdict
import sys
try:
    from modules.parse_cesar_output import classify_exon
except ImportError:
    from parse_cesar_output import classify_exon

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]

# TODO: revise this
PID_ANN_THR = 45
BLOSUM_ANN_THR = 25
MAX_SCORE = 1000
MAX_COLOR = 255
PID_HQ_THR = 65
BLOSUM_HQ_THR = 35

Q_HEADER_FIELDS_NUM = 12

# thresholds for exon classification
HQ_PID = 75
HQ_BLOSUM = 65

AB_INCL_PID = 25
AB_INCL_BLOSUM = 25

C_INCL_PID = 65
C_INCL_BLOSUM = 50

META_HEADER = "\t".join("gene exon_num chain_id act_region exp_region"
                        " in_exp pid blosum gap class paralog q_mark".split())


def eprint(msg, end="\n"):
    """Like print but for stderr."""
    sys.stderr.write(str(msg) + end)


def die(msg, rc=0):
    """Write msg to stderr and abort program."""
    eprint(msg)
    sys.exit(rc)


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
    app.add_argument("--exons_left", help="Exons that are left, a tsv", default=None)
    app.add_argument("--skipped", help="Save skipped genes", default=None)
    app.add_argument("--exon_class", help="Produce exon/class track", default=None)
    app.add_argument("--verbose", action="store_true", dest="verbose", help="Show verbose messages")
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def read_fasta(fasta_line, v=False, show_headers=False, rm=""):
    """Read fasta, return dict and type."""
    fasta_data = fasta_line.split(">")
    eprint(f"fasta_data[0] is:\n{fasta_data[0]}") if v else None
    eprint(f"fasta_data[1] is:\n{fasta_data[1]}") if v else None
    if fasta_data[0] != "":
        # this is a bug
        eprint("ERROR! Cesar output is corrupted")
        eprint(f"Issue detected in the folling string:\n{fasta_line}")
        die("Abort")
    # assert fasta_data[0] == ""  # if a file starts with > it should be empty
    del fasta_data[0]  # remove it "" we don't need that
    sequences = {}  # accumulate data here
    order = []  # to have ordered list

    for elem in fasta_data:
        raw_lines = elem.split("\n")
        # it must be first ['capHir1', 'ATGCCGCGCCAATTCCCCAAGCTGA... ]
        header = raw_lines[0]
        # separate nucleotide-containing lines
        lines = [x for x in raw_lines[1:] if x != "" and not x.startswith("!")]
        if len(lines) == 0:  # it is a mistake - empty sequene --> get rid of
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


def get_itemRgb(score, seq_type=None):
    """Return color for a gene according the score computed."""
    # should be RED,GREEN,BLUE each 0:255
    # bigger score --> reddish one
    score_coeff = score / MAX_SCORE
    if seq_type == "LOST":
        red = int(MAX_COLOR * score_coeff)
        green = int(red / 4)
        blue = int(red / 4)
    elif seq_type == "GREY":
        red = int(MAX_COLOR * score_coeff)
        green = red
        blue = red
    else:  # normal case, usual intact sequence
        blue = int(MAX_COLOR * score_coeff)
        green = int(blue / 4)
        red = int(blue / 4)

    # normalize just in case
    red = red if red <= MAX_COLOR else MAX_COLOR
    blue = blue if blue <= MAX_COLOR else MAX_COLOR
    green = green if green <= MAX_COLOR else MAX_COLOR
    itemRgb = f"{red},{green},{blue}"
    return itemRgb


def read_paralogs_log(paralog_log_file):
    """Get list of paralogous gene-chains."""
    f = open(paralog_log_file)
    gene_chains = [x for x in f.read().split("\n") if x != ""]
    f.close()
    return set(gene_chains)


def parse_cesar_bdb(arg_input, v=False):
    """Entry point."""
    in_ = open(arg_input, "r")
    # two \n\n divide each unit of information
    content = [x for x in in_.read().split("\n\n") if x]
    in_.close()

    # initiate collectors
    bed_lines = []
    skipped = []
    pred_seq_chain = {}
    t_exon_seqs = defaultdict(dict)
    exons_left = defaultdict(list)
    wrong_exons = []
    all_meta_data = [META_HEADER]
    prot_data = []
    exon_class_data = []

    for elem in content:
        # one elem - one CESAR call (one gene and >=1 chains)
        # now loop gene-by-genew
        gene = elem.split("\n")[0][1:]
        eprint(f"Reading gene {gene}") if v else None
        # error = False
        cesar_out = "\n".join(elem.split("\n")[1:])
        sequences, order = read_fasta(cesar_out, v=v)
        ranges_chain, chain_dir = defaultdict(dict), {}
        pred_seq_chain[gene] = defaultdict(dict)
        exon_lens, chain_raw_scores = {}, defaultdict(list)
        chain_pid_scores, chain_blosum_scores = defaultdict(list), defaultdict(list)
        chain_classes = defaultdict(list)

        query_headers = [h for h in order if h.endswith("query_exon")]
        ref_headers = [h for h in order if h.endswith("reference_exon")]
        prot_ids = [h for h in order if "PROT" in h]

        # parse reference exons
        for header in ref_headers:
            # one header for one exon
            header_fields = [s.replace(" ", "") for s in header.split("|")]
            exon_num = int(header_fields[1])
            exon_seq = sequences[header]
            seq_len = len(exon_seq)
            t_exon_seqs[gene][exon_num] = exon_seq
            exon_lens[exon_num] = seq_len
        
        # save protein data
        for prot_id in prot_ids:
            prot_seq = sequences[prot_id]
            prot_line = f">{prot_id}\n{prot_seq}\n"
            prot_data.append(prot_line)

        gene_len = sum(exon_lens.values())

        # get gene: exons dict to trace deleted exons
        gene_chain_exon_status = defaultdict(dict)

        for num, header in enumerate(query_headers):
            header_fields = [s.replace(" ", "") for s in header.split("|")]
            if len(header_fields) != Q_HEADER_FIELDS_NUM:
                continue  # ref exon

            # extract metadata
            trans = header_fields[0]
            exon_num = int(header_fields[1])
            chain_id = int(header_fields[2])
            exon_region = read_region(header_fields[3])
            pid = float(header_fields[4])
            blosum = float(header_fields[5])

            is_gap = header_fields[6]
            exon_class = header_fields[7]
            exp_region_str = header_fields[8]
            in_exp = header_fields[9]
            in_exp_b = True if in_exp == "INC" else False

            chain_pid_scores[chain_id].append(pid)
            chain_blosum_scores[chain_id].append(blosum)
            chain_classes[chain_id].append(exon_class)
            para_annot = True if header_fields[10] == "True" else False
            stat_key = (trans, chain_id)
            exon_decision, q_mark = classify_exon(exon_class, in_exp_b, pid, blosum)
            exon_class_track = (trans, str(chain_id), str(exon_num), header_fields[3], q_mark)
            exon_class_data.append(exon_class_track)

            # time to classify exons
            try:
                exon_score = int(pid / 100 * exon_lens[exon_num])
            except KeyError:
                # error = True
                pid, exon_score = 0, 0
            
            if exon_decision is False:
                # too bad exon
                wrong_exons.append(header)
                gene_chain_exon_status[stat_key][exon_num] = False
                chain_raw_scores[chain_id].append(0)
            else:
                chain_raw_scores[chain_id].append(exon_score)
                # get/write necessary info
                gene_chain_exon_status[stat_key][exon_num] = True
                directed = exon_region["end"] > exon_region["start"]
                chain_dir[chain_id] = directed
                ranges_chain[chain_id][exon_num] = exon_region
                pred_seq_chain[gene][chain_id][exon_num] = sequences[header]
                exons_left[(gene, chain_id)].append(str(exon_num))

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
            # no exons left for the particular name
            # need to log it
            name_ = f"{name[0]}.{name[1]}"
            skipped.append(f"{name_}\tall exons are deleted.")

        # make bed regions
        for chain_id in chain_dir.keys():
            # go gene.chain-by-gene.chain
            # get chain score
            chain_r_score = sum(chain_raw_scores[chain_id])
            chain_score = int(MAX_SCORE * (chain_r_score / gene_len))
            chain_pids = chain_pid_scores[chain_id]
            chain_blosums = chain_blosum_scores[chain_id]
            chain_cls = chain_classes[chain_id]

            # check if gene is high confident\quality
            pid_HQ = all(x >= PID_HQ_THR for x in chain_pids)
            blosum_HQ = all(x >= BLOSUM_HQ_THR for x in chain_blosums)
            class_HQ = all(x == "A" for x in chain_cls)
            HQ = pid_HQ and blosum_HQ and class_HQ

            blockStarts, blockSizes = [], []
            ranges = ranges_chain[chain_id]
            name = f"{gene}.{chain_id}"
            stat_key = (gene, chain_id)
            if HQ:  # add mark that it's high quality gene
                # name = f"HQ_{name}"
                # TODO: implement this better
                pass

            if len(ranges) == 0:
                # this gene is completely missed
                skipped.append(f"{name}\tall exons are deleted.")
                continue
            direct = chain_dir[chain_id]
            exon_nums = sorted(ranges.keys()) if direct \
                else sorted(ranges.keys(), reverse=True)

            # naming convention as on UCSC browser site
            chrom = ranges[exon_nums[0]]["chrom"]
            chromStart = ranges[exon_nums[0]]["start"] if direct \
                else ranges[exon_nums[0]]["end"]
            chromEnd = ranges[exon_nums[-1]]["end"] if direct \
                else ranges[exon_nums[-1]]["start"]

            score = chain_score
            thickStart = chromStart  # for now
            thickEnd = chromEnd
            itemRgb = get_itemRgb(score)
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
    fasta_lines = ""
    for gene, chain_exon_seq in pred_seq_chain.items():
        # write target gene info
        t_gene_seq_dct = t_exon_seqs.get(gene)
        if t_gene_seq_dct is None:
            eprint(f"Warning! Missed data for {gene}")
            skipped.append(f"{gene}\tmissed data after cesar stage")
            continue
        t_exon_nums = sorted(t_gene_seq_dct.keys())
        t_header = ">ref_{0}\n".format(gene)
        t_seq = "".join([t_gene_seq_dct[num] for num in t_exon_nums]) + "\n"
        fasta_lines += t_header
        fasta_lines += t_seq

        # and query info
        for chain_id, exon_seq in chain_exon_seq.items():
            track_header = ">{0}.{1}\n".format(gene, chain_id)
            exon_nums = sorted(exon_seq.keys())
            seq = "".join([exon_seq[num] for num in exon_nums]) + "\n"
            fasta_lines += track_header
            fasta_lines += seq

    # re-format wrong exons to bed-6 lines
    # to make it possible to save them and visualize
    # in the browser
    trash_exons = []
    for elem in wrong_exons:
        elem_fields = [s.replace(" ", "") for s in elem.split("|")]
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

    meta_str = "\n".join(all_meta_data) + "\n"
    skipped_str = "\n".join(skipped) + "\n"
    prot_fasta = "".join(prot_data)
    return bed_lines, trash_exons, fasta_lines, meta_str, prot_fasta, skipped_str


def main():
    """Enrty point as stand-alone tool."""
    args = parse_args()
    # basically call the function that extracts everything
    parsed_data = parse_cesar_bdb(args.input, v=args.verbose)

    bed_lines = parsed_data[0]
    trash_exons = parsed_data[1]
    fasta_lines = parsed_data[2]
    exons_meta = parsed_data[3]
    prot_fasta = parsed_data[4]
    skipped = parsed_data[5]

    if args.bed:  # save bed
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

    if args.prot_fasta:
        f = open(args.prot_fasta, "w") if args.prot_fasta != "stdout" else sys.stdout
        f.write(prot_fasta)
        f.close() if args.prot_fasta != "stdout" else None

    if args.skipped:
        f = open(args.skipped, "w") if args.skipped != "stdout" else sys.stdout
        f.write(skipped)
        f.close()
    sys.exit(0)


if __name__ == "__main__":
    main()
