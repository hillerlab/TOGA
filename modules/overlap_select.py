"""Python replacement for overlapSelect.

For a chain and a set of genes returns the following:
gene: how many bases this chain overlap in exons.
"""

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]



def make_bed_ranges(bed_line):
    """Convert a bed line to a set of exon ranges."""
    line_info = bed_line.split("\t")
    # parse bed-12 line according to specification
    chrom_start = int(line_info[1])
    gene_name = line_info[3]
    # basically get exon absolute coordinates
    block_starts = [chrom_start + int(x) for x in line_info[11].split(",") if x != ""]
    block_sizes = [int(x) for x in line_info[10].split(",") if x != ""]
    blocks_num = int(line_info[9])
    block_ends = [block_starts[i] + block_sizes[i] for i in range(blocks_num)]
    # create and return (exon start, exon end, gene name) tuples
    genes = [gene_name for _ in range(blocks_num)]
    return list(zip(block_starts, block_ends, genes))


def parse_bed(bed_lines):
    """Return sorted genic regions."""
    ranges = []
    for bed_line in bed_lines.split("\n")[:-1]:
        gene_ranges = make_bed_ranges(bed_line)
        ranges.extend(gene_ranges)
    return list(sorted(ranges, key=lambda x: x[0]))


def chain_reader(chain):
    """Yield chain blocks one by one."""
    chain_data = chain.split("\n")
    chain_head = chain_data[0]
    del chain_data[0]  # keep blocks only in the chain_data list
    # define starting point
    progress = int(chain_head.split()[5])
    blocks_num = len(chain_data)
    for i in range(blocks_num):
        # read block-by-block
        block_info = chain_data[i].split()
        if len(block_info) == 1:
            # only the last block has length == 1
            block_size = int(block_info[0])
            block_start = progress
            block_end = block_start + block_size
            yield block_start, block_end
            break  # end the loop
        
        # get intermediate block data
        block_size = int(block_info[0])
        dt = int(block_info[1])
        block_start = progress
        block_end = block_start + block_size
        progress = block_end + dt
        yield block_start, block_end


def intersect(ch_block, be_block):
    """Return intersection size."""
    return min(ch_block[1], be_block[1]) - max(ch_block[0], be_block[0])


def overlap_select(bed, chain):
    """Python implementation of some overlapSelect (kent) functionality."""
    ranges = parse_bed(bed)  # list of exon ranges
    genes = [x[2] for x in ranges]

    # bed overlaps: our results, count overlapped bases per gene
    bed_overlaps = {gene: 0 for gene in genes}
    chain_len = 0  # sum of chain blocks
    start_with = 0  # bed tracks counter

    # init blocks generator
    for block in chain_reader(chain):
        # add to len
        chain_len += block[1] - block[0]  # block[1] - block[0] is block length
        FLAG = False  # was there an intersection or not?
        FIRST = True

        while True:
            if FIRST:  # start with previous start, first iteration
                bed_num = start_with
                FIRST = False  # guarantee that this condition works ONCE per loop
            else:  # just increase the pointer
                bed_num += 1  # to avoid inf loop

            if bed_num >= len(ranges):
                # we intersected all beds, they are over
                break

            exon = ranges[bed_num]  # pick the exon
            if block[1] < exon[0]:  # too late
                break  # means that bed is "righter" than chain block
            
            # get intersection between the chain block and exon
            block_vs_exon = intersect(block, (exon[0], exon[1]))
            if block_vs_exon > 0:  # they intersect
                if not FLAG:  # the FIRST intersection of this chain block
                    # shift the start: all exons leftside will never be reached
                    start_with = bed_num  # guarantee that I will assign to starts with
                                          # only the FIRST intersection (if it took place)
                    FLAG = True  # otherwise starts_with will be preserved
                # add the intersection size
                bed_overlaps[exon[2]] += block_vs_exon

            else:  # we recorded all the region with intersections
                if block[0] > exon[1]:  # too early
                    # in case like:
                    # gene A EEEEE----------------------------------------EEEEEE #
                    # gene B               EEEEEEEEEE                            #
                    # gene C                               EEEEEEEEE             #
                    # chain                                    ccccc             #
                    # at gene A I will get FLAG = True and NO intersection with gene B
                    # --> I would miss gene C in this case without this condition.
                    continue

                elif FLAG:  # this is not a nested gene
                    break  # and all intersections are saved --> proceed to the next chain

    # return the required values
    return chain_len, bed_overlaps
