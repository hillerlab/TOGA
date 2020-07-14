#!/usr/bin/env python3
"""Join projection in genes for query annotation."""
import argparse
import sys
import os
from collections import defaultdict
try:
    from modules.common import flatten
except ImportError:
    from common import flatten


def parse_args():
    """Read CMD args."""
    app = argparse.ArgumentParser()
    app.add_argument("query_bed", help="Query annotation bed file.")
    app.add_argument("output", help="Output containing query isoforms data")
    app.add_argument("--genes_track", help="Save gene borders track")
    if len(sys.argv) < 3:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def intersect_ranges(region_1, region_2):
    """Return TRUE if ranges intersect."""
    inter_size = min(region_1[1], region_2[1]) - max(region_1[0], region_2[0])
    return True if inter_size > 0 else False


def read_query_bed(bed_file):
    """Extract everything what is possible from a query bed."""
    chrom_dir_enrties = defaultdict(list)
    gene_exons_coord = defaultdict(dict)

    f = open(bed_file, "r")
    for line in f:
        bed_info = line[:-1].split("\t")
        chrom = bed_info[0]
        chromStart = int(bed_info[1])
        name = bed_info[3]
        strand = bed_info[5]
        thickStart = int(bed_info[6])
        thickEnd = int(bed_info[7])
        blockCount = int(bed_info[9])
        blockSizes = [int(x) for x in bed_info[10].split(',') if x != '']
        blockStarts = [int(x) for x in bed_info[11].split(',') if x != '']
        chrom_dir_enrties[(chrom, strand)].append((name, thickStart, thickEnd))
        # add exon ranges
        for i in range(blockCount):
            exon_i_start = blockStarts[i] + chromStart
            exon_i_end = exon_i_start + blockSizes[i]
            gene_exons_coord[name][i] = (exon_i_start, exon_i_end)
    f.close()
    return chrom_dir_enrties, gene_exons_coord


def intersect_genes(gene_1, gene_2):
    """Return true if at least one exon intersects in both genes."""
    exons_1 = gene_1.values()
    exons_2 = gene_2.values()
    for exon_1 in exons_1:
        for exon_2 in exons_2:
            # compare exons one-to-one
            intersect = intersect_ranges(exon_1, exon_2)
            if intersect:
                return True
    # no intersection detected
    return False


def merge_regions(chrom_dir_enrties, gene_exons_coord):
    """Split query genes in regions.

    Return query_gene: Region ID dict
    Also region_ID: coordinates."""
    region_id = 0
    region_to_projections = defaultdict(set)
    # not the best solution ever, but it works:
    region_range_buffer = {}

    # split genes in region IDs
    for _, genes in chrom_dir_enrties.items():
        # chrom dir -> chromosome + strain!
        region_id += 1  # if we look to a new chrom --> region must be different
        if len(genes) == 1:  # nothing to intersect with
            gene = genes[0]  # anyway there is only one gene
            region_to_projections[region_id].add(gene[0])
            continue
        # there is something to intersect!
        genes_sorted = sorted(genes, key=lambda x: x[1])
        # add the first element
        # first_here = genes[0]
        first_here = genes_sorted[0]
        region_to_projections[region_id].add(first_here[0])

        for i in range(1, len(genes_sorted)):
            # get current and previous transcript
            # beginning of the gene-oriented
            current = genes_sorted[i]
            current_range = (current[1], current[2])
            prev = genes_sorted[i - 1]

            if not region_range_buffer.get(region_id):
                # nothing in this region before -> save prev gene region to buffer
                prev_range = (prev[1], prev[2])
                region_range_buffer[region_id] = prev_range
            else:
                # add previous to buffer, merge these regions
                buffered = region_range_buffer[region_id]
                prev_range = (min(buffered[0], prev[1]), max(buffered[1], prev[2]))
                region_range_buffer[region_id] = prev_range

            ranges_intersect = intersect_ranges(current_range, prev_range)

            # even ranges do not intersect --> goto next region ID
            if not ranges_intersect:
                region_id += 1
                region_to_projections[region_id].add(current[0])
                continue

            # we are here - let's check if genes intersects with their exons
            current_exons = gene_exons_coord[current[0]]
            prev_exons = gene_exons_coord[prev[0]]
            genes_intersect = intersect_genes(current_exons, prev_exons)

            # genes dont intersect --> goto the next region then
            if not genes_intersect:
                region_id += 1
                region_to_projections[region_id].add(current[0])
                continue
            # they do intersect actually
            region_to_projections[region_id].add(current[0])

    return region_to_projections


def save_regions_data(reg_to_proj, out_file):
    """Save regions to projections data."""
    f = open(out_file, "w") if out_file != "stdout" else sys.stdout
    f.write("Region_ID\tProjection_ID\n")
    for reg_num, projections in reg_to_proj.items():
        for proj_id in projections:
            f.write(f"reg_{reg_num}\t{proj_id}\n")
    f.close() if out_file != "stdout" else None


def get_regions_ranges(region_to_proj, chrom_to_ranges):
    """We like to save region ID (gene) to range data."""
    proj_to_grange = {}
    proj_to_chrom = {}
    for chr_reg, p_ranges in chrom_to_ranges.items():
        chrom = chr_reg[0]
        for p_range_named in p_ranges:
            p_name = p_range_named[0]
            p_range = (p_range_named[1], p_range_named[2])
            proj_to_grange[p_name] = p_range
            proj_to_chrom[p_name] = chrom
    region_id_to_range = {}
    for region_id, projections in region_to_proj.items():
        chrom = proj_to_chrom[iter(projections).__next__()]
        granges = [proj_to_grange[x] for x in projections]
        granges_vals = flatten(granges)
        reg_start = min(granges_vals)
        reg_end = max(granges_vals)
        region_id_to_range[region_id] = (chrom, reg_start, reg_end)
    return region_id_to_range


def save_bed4(region_id_to_range, bed_out):
    """Save bed file containing predicted genes ranges."""
    f = open(bed_out, "w") if bed_out != "stdout" else sys.stdout
    for region_id, grange in region_id_to_range.items():
        reg_str = f"reg_{region_id}"
        chrom = grange[0]
        start = grange[1]
        end = grange[2]
        f.write(f"{chrom}\t{start}\t{end}\t{reg_str}\n")
    f.close() if bed_out != "stdout" else None

    
def get_query_isoforms_data(query_bed, query_isoforms, save_genes_track=None):
    """Create isoforms track for query."""
    chrom_dir_enrties, gene_exons_coord = read_query_bed(query_bed)
    region_to_proj = merge_regions(chrom_dir_enrties, gene_exons_coord)
    save_regions_data(region_to_proj, query_isoforms)
    if save_genes_track:
        region_to_range = get_regions_ranges(region_to_proj, chrom_dir_enrties)
        save_bed4(region_to_range, save_genes_track)

if __name__ == "__main__":
    args = parse_args()
    get_query_isoforms_data(args.query_bed, args.output, save_genes_track=args.genes_track)
