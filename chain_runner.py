#!/usr/bin/env python3
"""Script to run chain classification job.

Extract features from each chain to gene intersection.
Loads a list of chain: genes tasks and calls
modules.processor.unit for each chain: genes task.
"""
import argparse
import sys
import os
from datetime import datetime as dt
from modules.overlap_select import overlap_select
from modules.common import bed_extract_id
from modules.common import make_cds_track
from modules.common import die
from modules.common import load_chain_dict
from modules.common import setup_logger
from modules.common import to_log
from version import __version__

__author__ = "Bogdan M. Kirilenko"

FLANK_SIZE = 10000  # gene flank size -> for flank ali feature
COMBINED_BED_ID = "COMBINED"  # placeholder gene name for intermediate tracks
ALL_EXONS_COMBINED = "ALL_EXONS_COMBINED"


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument(
        "input_file",
        type=str,
        help="File containing chain to genes lines."
        'Also you can use "chain [genes]" as a single argument.',
    )
    app.add_argument(
        "bed_file", type=str, help="BDB file containing annotation tracks."
    )
    app.add_argument("chain_file", type=str, help="Chain file.")
    app.add_argument("--log_file", type=str, help="Log file")
    app.add_argument(
        "--verbose", "-v", action="store_true", dest="verbose", help="Verbose messages."
    )
    app.add_argument(
        "--extended",
        "-e",
        action="store_true",
        dest="extended",
        help="Write the output in extended (human readable) format. "
        "Is not recommended for genome-wide scale.",
    )
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def intersect(range_1, range_2):
    """Return intersection size."""
    return min(range_1[1], range_2[1]) - max(range_1[0], range_2[0])


def marge_ranges(range_1, range_2):
    """Return merged range."""
    return min(range_1[0], range_2[0]), max(range_1[1], range_2[1])


def extract_chain(chain_file, chain_dict, chain):
    """Extract chain string.

    We have: chain file, chain_id, start byte and offset.
    """
    f = open(chain_file, "rb")
    start, offset = chain_dict.get(int(chain))
    f.seek(start)  # jump to start_byte_position
    chain = f.read(offset).decode("utf-8")  # read OFFSET bytes
    f.close()
    return chain


def check_args(
    chain_id, genes, chain_file, chain_dict, bed_file, verbose_level, work_data, result
):
    """Check if arguments are correct, extract initial data if so."""
    work_data["chain_id"] = chain_id
    # check genes
    raw_genes = [x for x in genes.split(",") if x != ""]
    # bed_lines = bedExtractSqlite(raw_genes, bed_index, bed_file)
    bed_lines = bed_extract_id(bed_file, raw_genes)
    work_data["bed"] = bed_lines  # save it
    work_data["genes"] = [x.split("\t")[3] for x in bed_lines.split("\n")[:-1]]

    # check if numbers of genes are equal
    if len(raw_genes) != len(bed_lines.split("\n")[:-1]):
        to_log("Warning. Not all the transcripts were found!\n")
        need_ = len(raw_genes)
        extracted_ = len(bed_lines.split("\n")[:-1])
        to_log(f"Expected {need_} transcripts, extracted {extracted_}")
        missing_genes = ",".join([x for x in raw_genes if x not in work_data["genes"]])
        to_log(f"Missing transcripts:\n{missing_genes}")

    # extract chain body from the file
    work_data["chain"] = extract_chain(chain_file, chain_dict, chain_id)

    # parse chain header
    chain_header = work_data["chain"].split("\n")[0].split()
    q_start = int(chain_header[10])
    q_end = int(chain_header[11])
    q_len = abs(q_end - q_start)
    work_data["chain_QLen"] = q_len
    work_data["chain_Tstarts"] = int(chain_header[5])
    work_data["chain_Tends"] = int(chain_header[6])
    result["chain_global_score"] = int(chain_header[1])
    result["chain_len"] = work_data["chain_Tends"] - work_data["chain_Tstarts"]


def read_input(input_file):
    """Read input."""
    # it must be chain TAB genes line
    if os.path.isfile(input_file):
        tasks = {}
        f = open(input_file)
        for line in f:
            line_info = line.rstrip().split("\t")
            chain = line_info[0]
            genes = line_info[1]
            tasks[chain] = genes
        f.close()
        return tasks
    elif len(input_file.split()) == 2:
        # it is not a file but chain<space>[,-sep list of genes]
        chain = input_file.split()[0]
        genes = input_file.split()[1]
        return {chain: genes}
    else:
        err_msg = (
            "Error! Wrong input. Please provide either a file containing chain to genes\n"
            'list or a "chain<space>[comma-separated list of genes]" formatted-file'
        )
        die(err_msg)
        return


def bed12_to_ranges(bed):
    """Convert bed-12 file to set of sorted ranges."""
    ranges_unsort, chrom = [], None
    for line in bed.split("\n")[:-1]:
        # parse line and extract blocks
        line_info = line.split("\t")
        chrom = line_info[0]
        glob_start = int(line_info[1])
        blocks_num = int(line_info[9])
        block_sizes = [int(x) for x in line_info[10].split(",") if x != ""]
        block_starts = [
            glob_start + int(x) for x in line_info[11].split(",") if x != ""
        ]
        block_ends = [block_starts[i] + block_sizes[i] for i in range(blocks_num)]
        for i in range(blocks_num):  # save the range for each exon
            ranges_unsort.append((block_starts[i], block_ends[i]))
    # return sorted ranges
    die("(bed12_to_ranges) error, cannot read bed properly") if not chrom else None
    return chrom, sorted(ranges_unsort, key=lambda x: x[0])


def bedcov_ranges(ranges, chrom):
    """Return a set of exons without overlaps.

    Python re-implementation of bedCov (kent) functionality.
    """
    ranges_filtered, pointer = [ranges[0]], 0  # initial values for filter
    gene = COMBINED_BED_ID  # if there is a mixture of genes - no ID anyway
    nested = False  # default value
    for i in range(1, len(ranges)):  # ranges are sorted so we can
        # compare each with only the next one
        if intersect(ranges[i], ranges_filtered[pointer]) <= 0:
            pointer += 1  # we have no intersection
            ranges_filtered.append(ranges[i])
        else:  # intersect - add merged range to the pointer
            # pointer = pointer - don't move it
            # replace the last range with merged last + new one
            nested = True  # at least one pair intersected
            ranges_filtered[pointer] = marge_ranges(ranges_filtered[pointer], ranges[i])

    # chr | start | end | gene | - bed4 structure
    # now make bed4 file
    exons, template = [], "{0}\t{1}\t{2}\t{3}\n"
    for grange in ranges_filtered:
        exons.append(template.format(chrom, grange[0], grange[1], gene))
    return exons, nested


def check_nest(work_data, cds_bed):
    """Return True if genes are nested."""
    chrom, ranges = bed12_to_ranges(cds_bed)
    exons, nested = bedcov_ranges(ranges, chrom)
    work_data["exons"] = exons
    return nested


def get_tot_exons_track(work_data):
    """Get all exons including UTR and collapse them."""
    chrom, ranges = bed12_to_ranges(work_data["bed"])
    exons, _ = bedcov_ranges(ranges, chrom)
    bed_template = "{0}\t{1}\t{2}\t{6}\t1000\t+\t{1}\t{2}\t0,0,0\t{3}\t{4}\t{5}"
    gene = ALL_EXONS_COMBINED
    blocks_uns = [(int(x.split("\t")[1]), int(x.split("\t")[2])) for x in exons]
    blocks = sorted(
        blocks_uns, key=lambda x: x[0]
    )  # no guarantee that it is sorted initially
    bed_12_start = min([x[0] for x in blocks])
    bed_12_end = max(x[1] for x in blocks)
    block_starts = ",".join([str(x[0] - bed_12_start) for x in blocks]) + ","
    block_sizes = ",".join([str(x[1] - x[0]) for x in blocks]) + ","
    bed_12 = bed_template.format(
        chrom, bed_12_start, bed_12_end, len(exons), block_sizes, block_starts, gene
    )
    return bed_12


def collapse_exons(work_data):
    """Compensate nested genes."""
    # how bed12 looks like:
    # chr15	19964665	19988117	O	1000	+	19964665	19988117
    # 0,0,0	3	307,46,313,	0,390,22990,
    # I need to fill it with chrom, start and end and blocks info
    bed_template = "{0}\t{1}\t{2}\t{6}\t1000\t+\t{1}\t{2}\t0,0,0\t{3}\t{4}\t{5}"
    chrom, _, _, gene = work_data["exons"][0][:-1].split("\t")
    blocks_uns = [
        (int(x.split("\t")[1]), int(x.split("\t")[2])) for x in work_data["exons"]
    ]
    # TODO: fix duplicated code fragment
    blocks = sorted(
        blocks_uns, key=lambda x: x[0]
    )  # no guarantee that it is sorted initially
    bed_12_start = min([x[0] for x in blocks])
    bed_12_end = max(x[1] for x in blocks)
    block_starts = ",".join([str(x[0] - bed_12_start) for x in blocks]) + ","
    block_sizes = ",".join([str(x[1] - x[0]) for x in blocks]) + ","
    bed_12 = bed_template.format(
        chrom,
        bed_12_start,
        bed_12_end,
        len(work_data["exons"]),
        block_sizes,
        block_starts,
        gene,
    )
    work_data["nested"] = bed_12


def extend_bed_lines(bed_lines):
    """Create bed tracks for overlapSelect."""
    bed_lines_extended = ""  # init the variable to store the extended bed lines
    for line in bed_lines.split("\n")[:-1]:
        # verbose(f"Extending line:\n{line}")
        bed_lines_extended += line + "\n"  # first, I add the original bed line
        grange_track = line.split("\t")  # tab-separated file
        # create the second track for the genomic region of the same gene
        # also known as "gene body"
        grange_track[3] = (
            grange_track[3] + "_grange"
        )  # I add _grange for the gene name, mark it
        grange_track[11] = "0"  # one block --> one start, starts from 0
        # size of block == size of the gene
        grange_track[10] = str(int(grange_track[2]) - int(grange_track[1]))
        grange_track[9] = "1"  # it means that it will be only one block
        bed_lines_extended += "\t".join(grange_track) + "\n"

        # create a separate track for FLANKS
        flanks_track = line.split("\t")
        flanks_track[3] = flanks_track[3] + "_flanks"
        flanks_track[11] = "0"
        flanks_track[9] = "1"
        # need to avoid negative value here!
        # used 1 just for robustness
        flank_start = int(flanks_track[1]) - FLANK_SIZE
        flank_start = 1 if flank_start < 1 else flank_start
        flanks_track[1] = str(flank_start)
        flanks_track[6] = flanks_track[1]
        # no need to care about flank_end > chrSize
        # if not exists -> will be no blocks
        flanks_track[2] = str(int(flanks_track[2]) + FLANK_SIZE)
        flanks_track[7] = flanks_track[2]
        flanks_track[10] = str(int(flanks_track[2]) - int(flanks_track[1]))
        bed_lines_extended += "\t".join(flanks_track) + "\n"

        # add CDS track
        cds_track = make_cds_track(line)
        bed_lines_extended += cds_track + "\n"
    return bed_lines_extended


# def cound_cds_exons(bed_lines_extended):
#     """Count CDS exons in each gene/transcript."""
#     bed_lines = [x.rstrip().split("\t") for x in bed_lines_extended.split("\n") if x != ""]
#     cds_lines = [x for x in bed_lines if x[3].endswith("_CDS")]
#     ret = {x[3][:-4]: int(x[9]) for x in cds_lines}
#     return ret


def get_features(work_data, result, bed_lines_extended, nested=False):
    """Compute local exon overlap score.

    For every line in the bed file (X - exonic base, I - intronic base)
    the new line will be created, represents the genomic region (called gene_grange).
    After the overlapSelect have been called, it returns the number of bases in chain blocks
    overlapping the gene exons:
    ----XXXXIIIIXXXXIIIIXXXX----- gene A - 5 overlapped bases / in exons and blocks
    ----XXXXXXXXXXXXXXXXXXXX----- gene A_grange - 9 overlapped bases / in all genomic region and blocks
    --bbbb----bbb-----bbbb------- chain N
    Here we consider the entire gene as a single exon
    In this toy example local_fractionExonOverlap score would be 5/9

    The raw overlapSelectOutput looks like:
    #inId   selectId        inOverlap       selectOverlap   overBases       similarity      inBases selectBases
    ENSG00000167232 chr17   0.112   1       399     0.201   3576    399
    ENSG00000167232_grange    chr17   0.0111  1       399     0.022   35952   399
    """
    # call overlap select
    chain_glob_bases, local_exo_dict, bed_cov_times = overlap_select(
        bed_lines_extended, work_data["chain"]
    )
    nums_of_cds_exons_covered = [
        len(v) for k, v in bed_cov_times.items() if k.endswith("_CDS")
    ]
    max_num_of_cds_exons_covered = (
        max(nums_of_cds_exons_covered) if len(nums_of_cds_exons_covered) > 0 else 0
    )
    # compute for each gene finally
    chain_cds_bases = 0  # summarize global set here

    for gene in work_data["genes"]:
        # pick the data from overlap select table
        blocks_v_exons = local_exo_dict[gene]
        blocks_v_cds = local_exo_dict[gene + "_CDS"]
        blocks_v_gene = local_exo_dict[gene + "_grange"]
        blocks_v_flanks_and_gene = local_exo_dict[gene + "_flanks"]
        # cds_exons_num = gene_to_cds_exons[gene]

        # all exons - CDS exons -> UTR exons
        blocks_v_utr_exons = blocks_v_exons - blocks_v_cds
        # gene blocks - UTR exons -> gene without UTR exons
        blocks_v_no_utr_exons = blocks_v_gene - blocks_v_utr_exons

        # if something like this happened -> there is a bug
        assert blocks_v_exons >= blocks_v_utr_exons
        assert blocks_v_exons >= blocks_v_cds

        # blocks gene + flanks - blocks gene -> blocks X flanks
        blocks_v_flanks = blocks_v_flanks_and_gene - blocks_v_gene
        # blocks gene - blocks exons -> blocks introns
        blocks_v_introns = blocks_v_gene - blocks_v_exons
        assert blocks_v_introns >= 0

        flank_feature = blocks_v_flanks / (FLANK_SIZE * 2)

        # global counters
        # CDS bases increase with blocks V cds in the gene
        chain_cds_bases += blocks_v_cds
        # increase number of UTR exons
        # chain_utr_exon_bases += blocks_v_utr_exons

        # get local results
        result["gene_coverage"] += f"{gene}={blocks_v_cds},"
        result["gene_introns"] += f"{gene}={blocks_v_introns},"
        result["flanks_cov"] += f"{gene}={flank_feature},"
        local_exo = (
            blocks_v_cds / blocks_v_no_utr_exons
            if blocks_v_no_utr_exons != 0.0
            else 0.0
        )
        assert local_exo >= 0
        assert local_exo <= 1
        result["local_exons"] += "{0}={1},".format(gene, local_exo)
        # increase synteny if > 0 CDS bases covered
        if blocks_v_cds > 0:
            result["chain_synteny"] += 1
            ov_block = f"{gene}={work_data['chain_id']}"
            result["gene_overlaps"].append(ov_block)
        else:
            # verbose(f"Chain don't overlap any exons in {gene}")
            # TO CHECK - it was like this here:
            result["gene_overlaps"].append(f"{gene}=None")
    # do not forget about global feature
    # chain_glob_bases -= chain_utr_exon_bases  # ignore UTR exons!
    chain_v_all_exons = local_exo_dict[ALL_EXONS_COMBINED]
    chain_cds_bases = local_exo_dict[COMBINED_BED_ID] if nested else chain_cds_bases
    chain_v_utr_exons = chain_v_all_exons - chain_cds_bases
    q_len_corrected = work_data["chain_QLen"] - chain_v_utr_exons
    assert (
        q_len_corrected >= 0
    )  # chain length in query - blocks cover UTR cannot be a negative number
    result["global_exo"] = (
        chain_cds_bases / chain_glob_bases if chain_glob_bases != 0 else 0
    )
    # here we consider this transcript separately
    # nested genes do not affect this feature
    # chain_exon_bases = local_exo_dict[COMBINED_BED_ID] if nested else chain_exon_bases
    if max_num_of_cds_exons_covered > 1:
        result["Exlen_to_Qlen"] = (
            chain_cds_bases / q_len_corrected if q_len_corrected != 0 else 0
        )
    else:  # if chain covers at most 1 CDS exon -> this feature is not applicable
        result["Exlen_to_Qlen"] = 0
    assert result["Exlen_to_Qlen"] <= 1  # blocklen / qlen cannot be > 1


def extract_cds_lines(all_bed_lines):
    """Extract bed lines with names end with _CDS."""
    selected = []
    for line in all_bed_lines.split("\n"):
        if line == "":
            continue
        if line.split("\t")[3].endswith("_CDS"):
            selected.append(line)
    return "\n".join(selected) + "\n"


def make_output(work_data, result, t0):
    """Arrange the output."""
    # verbose("Making the output...")
    chain_fields = [
        "chain",
        work_data["chain_id"],
        result["chain_synteny"],
        result["chain_global_score"],
        result["global_exo"],
        result["Exlen_to_Qlen"],
        result["local_exons"],
        result["gene_coverage"],
        result["gene_introns"],
        result["flanks_cov"],
        result["chain_len"],
    ]
    chain_output = "\t".join([str(x) for x in chain_fields]) + "\n"
    genes_output = "genes\t{0}\n".format("\t".join(result["gene_overlaps"]))
    time_output = f"#estimated time: {dt.now() - t0}\n"
    return chain_output, genes_output, time_output


def extended_output(result, t0):
    """Make human-readable output for small tests."""
    chain_output = "Chain-related features:\n"
    for key, value in result.items():
        if key == "gene_overlaps":
            continue
        chain_output += f'"{key}": {value}\n'
    genes_output = "These genes are overlapped by these chains:\n{0}".format(
        "\t".join(result["gene_overlaps"])
    )
    time_output = f"#estimated time: {dt.now() - t0}\n"
    return chain_output, genes_output, time_output


def chain_feat_extractor(
    chain_id, transcripts, chain_file, bed_file, chain_dict, verbose_arg=None, extended=False
):
    """Chain features extractor entry point."""
    # global vars
    t0 = dt.now()
    to_log(f"processing chain_id: {chain_id} transcripts: {transcripts}")
    # global work_data: bed, chain, etc
    work_data = {
        "bed": "",
        "chain": "",
        "nested": None,
        "chain_id": "",
        "chain_Qlen": 0,
        "genes": [],
        "chain_len": 0,
        "chain_Tstarts": 0,
        "chain_Tends": 0,
        "chain_global_score": 0,
    }
    # structure to collect the results
    result = {
        "global_exo": 0.0,
        "flanks_cov": "",
        "gene_coverage": "",
        "gene_introns": "",
        "chain_synteny": 0,
        "local_exons": "",
        "gene_overlaps": [],
        "Exlen_to_Qlen": 0,
    }
    # check if all the files, dependencies etc are correct
    check_args(
        chain_id,
        transcripts,
        chain_file,
        chain_dict,
        bed_file,
        verbose_arg,
        work_data,
        result,
    )

    # the main part, computations
    bed_lines_extended = extend_bed_lines(work_data["bed"])
    cds_bed_lines = extract_cds_lines(bed_lines_extended)
    tot_track = get_tot_exons_track(work_data)
    bed_lines_extended += f"{tot_track}\n"
    nested = check_nest(work_data, cds_bed_lines)  # check if the genes are nested

    if not nested:
        # 99% cases go here
        # there are no nested genes
        get_features(work_data, result, bed_lines_extended)
    else:
        # another case, firstly need to make bed track with no intersections
        # and only after that call this function with flag NESTED for updated bed file
        collapse_exons(work_data)
        bed_lines_extended += f"{work_data['nested']}\n"
        get_features(work_data, result, bed_lines_extended, nested=True)
    # make a tuple with chain, genes and time output
    if not extended:
        # provide short version of output
        output = make_output(work_data, result, t0)
    else:
        # provide extended output
        # human-readable version
        output = extended_output(result, t0)
    return output


def main():
    """Entry point."""
    t0 = dt.now()
    args = parse_args()
    setup_logger(args.log_file, write_to_console=False)
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    # read input: meaning chain ids and gene names
    # there are 2 ways how they could be provided:
    # 1) Just a file, each line contains chain id and genes
    # 2) an argument: "chain ,-sep list of genes"
    batch = read_input(args.input_file)
    task_size = len(batch)
    to_log(f"processing {task_size} chains")
    # TODO: check whether I don't need .bst
    # load chains dict; it would be much faster to load chain_ID: (start_byte, offset)
    # python dict once than ask HDF5 database each time TOGA needs another chain
    index_file = args.chain_file.replace(".chain", ".chain_ID_position")
    chain_dict = load_chain_dict(index_file)

    # call main processing tool
    # TODO: rename genes to transcripts where appropropriate
    for job_num, (chain, transcripts) in enumerate(batch.items(), 1):
        # one unit: one chain + intersected genes
        # call routine that extracts chain feature
        unit_output = chain_feat_extractor(
            chain,
            transcripts,
            args.chain_file,
            args.bed_file,
            chain_dict,
            verbose_arg=args.verbose,
            extended=args.extended,
        )
        chain_output, genes_output, time_output = unit_output
        # stdout is used by subsequent TOGA commands
        sys.stdout.write(chain_output)
        sys.stdout.write(genes_output)
        sys.stdout.write(time_output)
        # sys.stderr.write(f"Job {job_num}/{task_size} done\r") if args.verbose else None
    to_log(f"Total job time: {dt.now() - t0}")


if __name__ == "__main__":
    main()
