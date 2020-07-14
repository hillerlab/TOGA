"""Common functions."""
# in progress
import sys
import sqlite3
import bsddb3

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]


def parts(lst, n=3):
    """Split an iterable into parts with size n."""
    return [lst[i:i + n] for i in iter(range(0, len(lst), n))]


def bedExtractID(index_file, gene_ids):
    """Extract chain from BDB file."""
    db = bsddb3.btopen(index_file, "r")
    if type(gene_ids) != str:
        keys = [str(gene_id).encode() for gene_id in gene_ids]
    else:
        keys = [str(gene_ids).encode()]
    bed_lines = "".join([db.get(key).decode("utf-8") for key in keys if db.get(key)])
    db.close()
    if len(bed_lines) == 0:
        sys.stderr.write("Error! Bed tracks not found! (bedExtractID)")
        sys.exit(1)
    return bed_lines


def bedExtractIDText(bed_file, gene_ids_param):
    """Extract bed-12 tracks according the gene names."""
    if type(gene_ids_param) != str:  # so it's list
        gene_ids = set(gene_ids_param)
    else:
        gene_ids = set([gene_ids_param])
    output = []
    f = open(bed_file, "r")
    for line in f:
        gene_id = line.split("\t")[3]
        if gene_id not in gene_ids:
            continue
        output.append(line)
    f.close()
    return "".join(output)


def make_cds_track(line):
    """Trim UTRs from bed track."""
    line_data = line.rstrip().split("\t")
    if len(line_data) != 12:
        sys.exit(f"Error! Bed line:\n{line}\nis a not bed-12 formatted line!")
    chrom = line_data[0]
    chromStart = int(line_data[1])
    chromEnd = int(line_data[2])
    name = line_data[3]  # gene_name usually
    name += "_CDS"
    bed_score = int(line_data[4])  # never used
    strand = line_data[5]  # otherwise:
    # strand = True if line_data[5] == '+' else False
    thickStart = int(line_data[6])
    thickEnd = int(line_data[7])
    itemRgb = line_data[8]  # never used
    blockCount = int(line_data[9])
    blockSizes = [int(x) for x in line_data[10].split(',') if x != '']
    blockStarts = [int(x) for x in line_data[11].split(',') if x != '']
    blockEnds = [blockStarts[i] + blockSizes[i] for i in range(blockCount)]
    blockAbsStarts = [blockStarts[i] + chromStart for i in range(blockCount)]
    blockAbsEnds = [blockEnds[i] + chromStart for i in range(blockCount)]
    blockNewStarts, blockNewEnds = [], []

    for block_num in range(blockCount):
        blockStart = blockAbsStarts[block_num]
        blockEnd = blockAbsEnds[block_num]

        # skip the block if it is entirely UTR
        if blockEnd <= thickStart:
            continue
        elif blockStart >= thickEnd:
            continue

        # remove UTRs
        blockNewStart = blockStart if blockStart >= thickStart else thickStart
        blockNewEnd = blockEnd if blockEnd <= thickEnd else thickEnd
        blockNewStarts.append(blockNewStart - thickStart)
        blockNewEnds.append(blockNewEnd - thickStart)

    blockNewCount = len(blockNewStarts)
    blockNewSizes = [blockNewEnds[i] - blockNewStarts[i]
                        for i in range(blockNewCount)]

    new_track = [chrom, thickStart, thickEnd, name, bed_score,
                    strand, thickStart, thickEnd, itemRgb, blockNewCount,
                    ",".join([str(x) for x in blockNewSizes]) + ",",
                    ",".join([str(x) for x in blockNewStarts]) + ","]
    new_line = "\t".join([str(x) for x in new_track])
    return new_line


def chainExtractID(index_file, chain_id):
    """Extract chain from BDB file."""
    db = bsddb3.btopen(index_file, "r")
    key = str(chain_id).encode()
    chain_data_enc = db.get(key)
    if not chain_data_enc:
        sys.stderr.write("Error! Chain {} not found in the bdb. Abort\n".format(chain_id))
        sys.exit(1)
    chain_data = chain_data_enc.decode("utf-8")
    db.close()
    return chain_data


def flatten(lst):
    """Flat list out of list of lists."""
    return [item for sublist in lst for item in sublist]


def split_proj_name(proj_name):
    """Split projection name."""
    proj_name_split = proj_name.split(".")
    q_num_str = proj_name_split[-1]
    trans_name = ".".join(proj_name_split[:-1])
    return trans_name, q_num_str
