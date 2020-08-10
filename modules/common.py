"""TOGA common functions."""
import os
import sys
import numpy as np
import ctypes
import h5py

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]

SLIB_NAME = "chain_bst_lib.so"


def parts(lst, n=3):
    """Split an iterable into parts with size n."""
    return [lst[i:i + n] for i in iter(range(0, len(lst), n))]


def bedExtractID(index_file, gene_ids):
    """Extract a bed track from a BDB file."""
    h = h5py.File(index_file, "r")
    # accept both a list of gene_ids (type=list) and a single gene_id (type=str)
    if type(gene_ids) != str:
        keys = [str(gene_id) for gene_id in gene_ids]
    else:
        keys = [str(gene_ids), ]
    bed_lines = []
    for key in keys:
        try:  # catch keyerror
            b_bed_line = h[key][()]
            u_type = f"U{len(b_bed_line)}"
            bed_line = b_bed_line.astype(u_type)
            bed_lines.append(bed_line)
        except KeyError:
            # a single item not found
            # crush only if ALL tracks not found
            continue
    h.close()
    if len(bed_lines) == 0:
        # if nothing found -> there must be an error
        sys.stderr.write("Error! Bed tracks not found! (bedExtractID)")
        sys.exit(1)
    return "".join(bed_lines)


def bedExtractIDText(bed_file, gene_ids_param):
    """Extract bed-12 tracks from a text bed file."""
    # the function accepts both a list of gene_ids (type=list)
    # and a single gene_id (type=string)
    if type(gene_ids_param) != str:  # so it's list
        gene_ids = set(gene_ids_param)
    else:
        gene_ids = set([gene_ids_param])
    # gene_ids is a set for consistency; even if a single gene id provided
    output = []
    f = open(bed_file, "r")
    for line in f:
        # read line-by-line and catch lines with gene_id in the gene_ids set
        gene_id = line.split("\t")[3]
        if gene_id not in gene_ids:
            continue
        output.append(line)
    f.close()
    return "".join(output)


def make_cds_track(line):
    """Trim UTRs from a bed track."""
    line_data = line.rstrip().split("\t")
    if len(line_data) != 12:
        sys.exit(f"Error! Bed line:\n{line}\nis a not bed-12 formatted line!")
    # parse bed12 line according to the specification
    chrom = line_data[0]
    chromStart = int(line_data[1])
    chromEnd = int(line_data[2])
    name = line_data[3]  # gene_name usually
    name += "_CDS"  # mark that UTRs are trimmed
    bed_score = int(line_data[4])  # never used
    strand = line_data[5]
    thickStart = int(line_data[6])
    thickEnd = int(line_data[7])
    itemRgb = line_data[8]  # never used
    blockCount = int(line_data[9])
    # chrom start and end define the entire transcript location
    # this includes both UTRs and CDS
    # thick start and end limit the CDS only
    blockSizes = [int(x) for x in line_data[10].split(',') if x != '']
    blockStarts = [int(x) for x in line_data[11].split(',') if x != '']
    blockEnds = [blockStarts[i] + blockSizes[i] for i in range(blockCount)]
    # block starts are given in the relative coordinates -> need to convert them
    # into absolute coordinates using chromstart
    blockAbsStarts = [blockStarts[i] + chromStart for i in range(blockCount)]
    blockAbsEnds = [blockEnds[i] + chromStart for i in range(blockCount)]
    # arrays for blocks with trimmed UTRs
    blockNewStarts, blockNewEnds = [], []

    for block_num in range(blockCount):
        # go block-by-block
        blockStart = blockAbsStarts[block_num]
        blockEnd = blockAbsEnds[block_num]

        # skip the block if it is entirely UTR
        if blockEnd <= thickStart:
            continue
        elif blockStart >= thickEnd:
            continue
        
        # if we are here: this is not an entirely UTR exon
        # it migth intersect the CDS border or to be in the CDS entirely
        # remove UTRs: block start must be >= CDS_start (thickStart)
        # block end must be <= CDS_end (thickEnd)
        blockNewStart = blockStart if blockStart >= thickStart else thickStart
        blockNewEnd = blockEnd if blockEnd <= thickEnd else thickEnd
        # save blocks with updated coordinates
        # also convert them back to relative coordinates with - thickStart
        # after the update thickStart/End are equal to chromStart/End
        blockNewStarts.append(blockNewStart - thickStart)
        blockNewEnds.append(blockNewEnd - thickStart)

    # blockCount could change due to entirely UTR exons
    blockNewCount = len(blockNewStarts)
    # this is also a subject to change
    blockNewSizes = [blockNewEnds[i] - blockNewStarts[i]
                     for i in range(blockNewCount)]

    # save the updated bed line with trimmed UTRs
    new_track = [chrom, thickStart, thickEnd, name, bed_score,
                 strand, thickStart, thickEnd, itemRgb, blockNewCount,
                 ",".join([str(x) for x in blockNewSizes]) + ",",
                 ",".join([str(x) for x in blockNewStarts]) + ","]
    new_line = "\t".join([str(x) for x in new_track])
    return new_line


def chainExtractID(index_file, chain_id, chain_file=None):
    """Extract chain text using index file."""
    # within TOGA should be fine:
    chain_file = chain_file if chain_file else index_file.replace(".bst", ".chain")
    if not os.path.isfile(chain_file):
        # need this check anyways
        sys.exit(f"chainExtractID error: cannot find {chain_file} file")
    # connect shared library
    script_location = os.path.dirname(__file__)
    slib_location = os.path.join(script_location, SLIB_NAME)
    sh_lib = ctypes.CDLL(slib_location)
    sh_lib.get_s_byte.argtypes = [ctypes.c_char_p,
                                ctypes.c_uint64,
                                ctypes.POINTER(ctypes.c_uint64),
                                ctypes.POINTER(ctypes.c_uint64)]
    sh_lib.get_s_byte.restype = ctypes.c_uint64

    # call library: find chain start byte and offset
    c_index_path = ctypes.c_char_p(index_file.encode())
    c_chain_id = ctypes.c_uint64(int(chain_id))
    c_sb = ctypes.c_uint64(0)
    c_of = ctypes.c_uint64(0)

    _ = sh_lib.get_s_byte(c_index_path,
                          c_chain_id,
                          ctypes.byref(c_sb),
                          ctypes.byref(c_of))
    
    if c_sb.value == c_of.value == 0:
        # if they are 0: nothing found then
        sys.stderr.write(f"Error, chain {chain_id} ")
        sys.stderr.write("not found\n")
        sys.exit(1)

    # we got start byte and offset
    f = open(chain_file, "rb")
    f.seek(c_sb.value)  # jump to start_byte_position
    chain = f.read(c_of.value).decode("utf-8")  # read OFFSET bytes
    f.close()
    return chain


def flatten(lst):
    """Flat list out of list of lists."""
    return [item for sublist in lst for item in sublist]


def split_proj_name(proj_name):
    """Split projection name.
    
    Projections named as follows: ${transcript_ID}.{$chain_id}.
    This function splits projection back into transcript and chain ids.
    We cannot just use split("."), because there migth be dots
    in the original transcript ID.
    """
    proj_name_split = proj_name.split(".")
    q_num_str = proj_name_split[-1]
    trans_name = ".".join(proj_name_split[:-1])
    return trans_name, q_num_str
