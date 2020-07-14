#!/usr/bin/env python3
"""Find all genic intersections for chain X bed.

Writes to stdout the following table:
chain_id<tab>comma-separated list of ovrelapped genes.
"""
import sys
from collections import defaultdict
import subprocess

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]


def flatten(lst):
    """Flat list out of list of lists."""
    return [item for sublist in lst for item in sublist]


def parse_chain(chain):
    """Return chrom: ranges from chain file."""
    # I need the headers only
    chrom_range = defaultdict(list)
    cmd = "cat {0} | grep chain".format(chain)
    headers = subprocess.check_output(cmd, shell=True).decode("utf-8")
    for header in headers.split("\n"):
        if len(header) == 0:
            continue
        header_info = header.split()
        chrom = header_info[2]
        start = int(header_info[5])
        end = int(header_info[6])
        chain_id = header_info[12]
        chrom_range[chrom].append((chain_id, start, end))
    return chrom_range


def parse_bed(bed):
    """Return chrom: ranges from bed file."""
    chrom_range = defaultdict(list)
    f = open(bed, "r")
    for line in f:
        line_info = line.split("\t")
        chrom = line_info[0]
        # start = int(line_info[1])
        # end = int(line_info[2])
        # count only CDS!
        start = int(line_info[6])
        end = int(line_info[7])
        gene = line_info[3]
        chrom_range[chrom].append((gene, start, end))
    f.close()
    return chrom_range


def intersect(range_1, range_2):
    """Return intersection size."""
    return min(range_1[2], range_2[2]) - max(range_1[1], range_2[1])


def find_first(chains, beds):
    """Find indexes for the first intersection."""
    if len(chains) == 0 or len(beds) == 0:
        return 0, 0
    if intersect(chains[0], beds[0]) > 0:  # no need to search
        return 0, 0
    first_bed_start, _ = beds[0][1], beds[0][2]
    _, first_chain_end = chains[0][1], chains[0][2]
    if first_chain_end < first_bed_start:
        # we have lots of chains in the beginning not intersecting beds
        for i in range(len(chains)):
            if intersect(chains[i], beds[0]) > 0:
                return i, 0
            elif chains[i][1] > beds[0][2]:
                return i, 1
            else:  # no intersection
                return 0, 0
    else:  # lots of beds in the beginning not covered by chains
        for i in range(len(beds)):
            if intersect(chains[i], beds[0]) > 0:
                return 0, i
            elif beds[i][1] > chains[0][2]:
                return 1, i
            else:  # no ibtersection
                return 0, 0


def overlap(chains, beds):
    """Return intersections for chain: bed."""
    # init state, find FIRST bed intersecting the FIRST chain
    chain_beds = defaultdict(list)
    chain_init, start_with = find_first(chains, beds)
    chains_num = len(chains)
    bed_len = len(beds)
    bed_num = None

    if chains_num == 0 or bed_len == 0:
        # check if there are chains and beds in this chrom
        return {}
    for i in range(chain_init, chains_num):
        FLAG = False  # was intersection or not?
        FIRST = True
        chain = chains[i]
        while True:
            if FIRST:  # start with previous start, first iteration
                bed_num = start_with
                FIRST = False  # guarantee that this condition works ONCE per loop
            else:  # just increase the pointer
                bed_num += 1  # to avoid inf loop

            if bed_num >= bed_len:
                break  # beds are over
            # pick the bed range
            bed = beds[bed_num]

            if chain[2] < bed[1]:  # too late
                break  # means that bed is "righter" than chain

            if intersect(chain, bed) > 0:
                if not FLAG:  # the FIRST intersection of this chain
                    start_with = bed_num  # guarantee that I will assign to starts with
                    # only the FIRST intersection (if it took place)
                    FLAG = True  # otherwise starts with will be preserved
                # save the intersection
                chain_beds[chain[0]].append(bed[0])

            else:  # we recorded all the region with intersections
                if chain[1] > bed[2]:  # too early
                    # in case like:
                    # gene A EEEEE----------------------------------------EEEEEE #
                    # gene B               EEEEEEEEEE                            #
                    # gene C                               EEEEEEEEE             #
                    # chain                                    ccccc             #
                    # at gene A I will get FLAG = True and NO intersection with gene B
                    # --> I will miss gene C in this case without this condition.
                    continue

                elif FLAG:
                    # this is not a nested gene
                    # and all intersections are saved
                    # --> proceed to the next chain
                    break

    return chain_beds


def chain_bed_intersect(chain, bed):
    """Entry point."""
    # get list of chrom: ranges for both
    skipped = []
    chain_data = parse_chain(chain)
    bed_data = parse_bed(bed)

    chroms = list(set(bed_data.keys()).intersection(chain_data.keys()))
    # to save beds that are skipped at this stage:
    only_bed_chroms = list(set(bed_data.keys()).difference(chain_data.keys()))
    genes_rejected = [x[0] for x in flatten([bed_data[chr] for chr in only_bed_chroms])]
    for gene in genes_rejected:
        skipped.append((gene, "chromosome is not aligned"))

    if len(chroms) == 0:
        sys.exit("Error! No common chromosomes between the bed and chain files found!")
    chain_bed_dict = {}  # out answer

    # main loop
    for chrom in chroms:
        # sort the ranges
        bed_ranges = sorted(bed_data[chrom], key=lambda x: x[1])
        gene_at_chrom = set(x[0] for x in bed_ranges)
        chain_ranges = sorted(chain_data[chrom], key=lambda x: x[1])
        chrom_chain_beds = overlap(chain_ranges, bed_ranges)
        # get genes that are not intersected with any chain
        genes_in = set(flatten(v for v in chrom_chain_beds.values()))
        genes_out = gene_at_chrom.difference(genes_in)
        for gene in genes_out:
            skipped.append((gene, "no intersecting chains"))
        chain_bed_dict.update(chrom_chain_beds)
    return chain_bed_dict, skipped


def save(dct, output="stdout"):
    """Save output in the file given."""
    f = open(output, "w") if output != "stdout" else sys.stdout
    for k, v in dct.items():
        f.write("{0}\t{1}\n".format(k, ",".join(v) + ","))
    f.close()


if __name__ == "__main__":
    try:  # read args
        chain_file = sys.argv[1]
        bed_file = sys.argv[2]
    except IndexError:
        sys.stderr.write("Usage: {} [chain_file] [bed file]\n".format(sys.argv[0]))
        sys.stderr.write("Output goes to stdout.\n")
        sys.exit(0)
    chain_bed_dict, skipped = chain_bed_intersect(chain_file, bed_file)
    save(chain_bed_dict)
