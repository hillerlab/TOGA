#!/usr/bin/env python3
"""Perform all operations to create BigBed file for UCSC browser."""
import os
from re import sub
import sys
import argparse
import subprocess
from collections import defaultdict
import getpass


GENERATE_TAB_FILES = "generate_tab_files.py"
GENERATE_PLOTS = "make_togaPlot.py"
LOCATION = os.path.dirname(__file__)
SCHEMA_LOCATION = os.path.join(LOCATION, "bb_schema.as")

IFEAT_PLACEHOLDER = ["0.0" for _ in range(6)]
PROT_PLACEHOLDER = ["<TT>NO_DATA</TT>"]
UNDEF = "UNDEFINED"
SVG_PLACEHOLDER = """<svg>
<rect fill="#aaa" stroke="#000" x="0" y="0" width="400" height="100"/>
<line x1="0" y1="0" x2="400" y2="100" stroke="red" stroke-width="4" />
<line x1="0" y1="100" x2="400" y2="0" stroke="red" stroke-width="4" />
</svg>
"""
PLOT_PLACEHOLDER = [SVG_PLACEHOLDER, ]


def parse_args():
    """Parse args."""
    app = argparse.ArgumentParser()
    app.add_argument("project_dir", help="Directory containing TOGA output")
    app.add_argument("--do_not_cleanup",
                     action="store_true",
                     dest="do_not_cleanup",
                     help="Do not clean tabs dir up.")
    app.add_argument("--bb_version",
                     default="v1",
                     help="Version tag")
    app.add_argument("--no_plots",
                     "--np",
                     dest="no_plots",
                     action="store_true",
                     help="If inactivating mutation plots are already generated, do not recreate them")
    # TODO: chroms sizes if not standard or not HL
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def call_gen_tab_files(project_dir):
    """Call Generate tab files subprocess."""
    print(f"Calling {GENERATE_TAB_FILES}")
    out_dir = os.path.join(project_dir, "tabs")
    os.mkdir(out_dir) if not os.path.isdir(out_dir) else None
    exe_ = os.path.join(LOCATION, GENERATE_TAB_FILES)
    cmd = f"{exe_} {project_dir} {out_dir}"
    rc = subprocess.call(cmd, shell=True)
    if rc != 0:
        raise ValueError(f"Command {cmd} died - cannot generate tab files")
    print("5 of 6 tab files generated")


def call_gen_plot_files(project_dir, no_plots=False):
    """Generate svg plots."""
    print(f"Calling {GENERATE_PLOTS}")
    out_file = os.path.join(project_dir, "tabs", "togaPlot.tab")
    if no_plots is True and os.path.isfile(out_file):
        print(f"--no_plots flag is on, found {out_file} --> skip making plots")
        return
    elif no_plots is True and not os.path.isfile(out_file):
        print(f"Error! --no_plots flag is on, but {out_file} not found\nAbort")
        sys.exit(1)
    exe_ = os.path.join(LOCATION, GENERATE_PLOTS)
    cmd = f"{exe_} {project_dir} {out_file}"
    rc = subprocess.call(cmd, shell=True)
    if rc != 0:
        raise ValueError(f"Command {cmd} died - cannot generate plots")
    print("All tab files generated")


def read_tab_file(tab_file):
    """Make trans_id: data track."""
    trans_to_dat = {}
    f = open(tab_file, "r")
    for line in f:
        line_data = line.rstrip().split("\t")
        trans = line_data[0]
        data = line_data[1:]
        trans_to_dat[trans] = data
    f.close()
    return trans_to_dat



def __gen_info_placeholder(p_name):
    """Generate info placeholder."""
    t_name = ".".join(p_name.split(".")[:-1])
    placeholder = [t_name, UNDEF, UNDEF, "0.0", "0", "0.0", "0.0", "0.0", "0.0", "0.0", UNDEF]
    return placeholder


def make_toga_data(project_dir):
    """Merge togaPlot, togaProt, togaInfo and togaInactFeat tables into togaData."""
    tabs_dir = os.path.join(project_dir, "tabs")
    toga_info_file = os.path.join(tabs_dir, "togaInfo.tab")
    toga_inact_feat_file = os.path.join(tabs_dir, "togaInactFeat.tab")
    toga_plot_file = os.path.join(tabs_dir, "togaPlot.tab")
    toga_prot_file = os.path.join(tabs_dir, "togaProt.tab")
    out_file = os.path.join(tabs_dir, "togaData.tab")

    print("Reading tab files to merge")
    trans_to_info = read_tab_file(toga_info_file)
    trans_to_ifeat = read_tab_file(toga_inact_feat_file)
    trans_to_plot = read_tab_file(toga_plot_file)
    trans_to_prot = read_tab_file(toga_prot_file)

    transcripts = trans_to_info.keys()
    print(f"Got data for {len(transcripts)} transcripts")

    f = open(out_file, "w")
    for t in transcripts:
        info = trans_to_info.get(t)
        if info is None:
            info = __gen_info_placeholder(t)
        ifeat = trans_to_ifeat.get(t, IFEAT_PLACEHOLDER.copy())
        prot = trans_to_prot.get(t, PROT_PLACEHOLDER.copy())
        plot = trans_to_plot.get(t, PLOT_PLACEHOLDER.copy())
        data_lst = [t] + info + ifeat + prot + plot
        tab_str = "\t".join(data_lst)
        f.write(tab_str)
        f.write("\n")
    f.close()
    # do cleanup
    os.remove(toga_info_file)
    os.remove(toga_inact_feat_file)
    os.remove(toga_plot_file)
    os.remove(toga_prot_file)
    print("Done")


def _check_seq_of_intervals_intersect(intervals):
    """Check whether any interval in the seq intersects"""
    intervals_num = len(intervals)
    for i in range(intervals_num - 1):
        # (start, end)
        curr_one = intervals[i]
        next_one = intervals[i + 1]
        # sorted by beginning
        # if start of the next < end of the curr
        # -> they intersect
        if next_one[0] < curr_one[1]:
            return (curr_one[0], curr_one[1])
    return None  # nothing suspicious found


def _filter_bed(in_bed, out_bed):
    """Filter out bed entries with:
     - None chromosome
     - Intersecting blocks."""
    _in = open(in_bed, "r")
    _out = open(out_bed, "w")
    for line in _in:
        ld = line.rstrip().split("\t")
        chrom = ld[0]
        if chrom == "None":
            continue
        block_sizes = [int(x) for x in ld[10].split(",") if x != ""]
        block_starts = [int(x) for x in ld[11].split(",") if x != ""]
        block_count = len(block_starts)
        if block_count == 1:
            # no need to check whether exons intersect if there is only 1 exon
            _out.write(line)
            continue
        block_ends = [block_starts[i] + block_sizes[i] for i in range(block_count)]
        exon_ranges = sorted(zip(block_starts, block_ends), key=lambda x: x[0])
        ranges_intersect = _check_seq_of_intervals_intersect(exon_ranges)
        if ranges_intersect:
            continue
        _out.write(line)
    _in.close()
    _out.close()


def make_sorted_bed(project_dir):
    """For bigBed, we need properly sorted bed file."""
    non_sorted_bed = os.path.join(project_dir, "query_annotation.bed")
    tabs_dir = os.path.join(project_dir, "tabs")
    bed_no_none = os.path.join(project_dir, "query_annotation_intact_chroms.bed")
    _filter_bed(non_sorted_bed, bed_no_none)
    sorted_bed = os.path.join(tabs_dir, "query_annot_sorted.bed")
    cmd = f"sort -k1,1 -k2,2n {bed_no_none} > {sorted_bed}"
    subprocess.call(cmd, shell=True)
    # sorted bed file, good


def __get_ts_from_bed(bed):
    """Get transcript IDs from bed file."""
    f = open(bed, 'r')
    ts = []
    ts_to_bed_line = {}
    for line in f:
        ld = line.rstrip().split("\t")
        t = ld[3]
        ts.append(t)
        ts_to_bed_line[t] = ld
    f.close()
    return ts, ts_to_bed_line


def __get_ts_to_toga_data(toga_data_file):
    """Get transID: togaData fields."""
    ret = {}
    f = open(toga_data_file, "r")
    for line in f:
        ld = line.rstrip().split("\t")
        t = ld[0]
        fields = ld[1:]
        ret[t] = fields
    f.close()
    return ret


def __get_ts_to_toga_inact(toga_inact_file):
    """Read inact mut data."""
    ts_to_exon_to_data = defaultdict(dict)
    f = open(toga_inact_file, 'r')
    for line in f:
        line_data = line.rstrip().split("\t")
        ts = line_data[0]
        exon_num = int(line_data[1])
        rest = line_data[2:]
        ts_to_exon_to_data[ts][exon_num] = rest
    f.close()
    # now: TS to HTML_formatted line
    trans_to_inact_line = {}
    for t in ts_to_exon_to_data.keys():
        exon_nums_unsort = ts_to_exon_to_data[t].keys()
        exon_nums = sorted(exon_nums_unsort)
        t_lines = []
        for exon_num in exon_nums:
            inact_mut_data = ts_to_exon_to_data[t][exon_num]
            """
            transcript varchar(50) not null,  -- unique projection ID: ${ref_transcript_ID}.${chain_ID}
            exon_num int unsigned not null,  -- exon number
            position int unsigned not null,  -- possition where mutation happened
            mut_class varchar(15) not null,  -- mutation class such as FS deletion
            mutation varchar(20) not null,  -- what exactly happened
            is_inact tinyint unsigned not null,  -- is this mutation inactivating, yes 1 or not 0
            mut_id varchar(10) not null  -- mut identifier
            """
            # template:
            """
            struct togaInactMut *info = NULL;
            info = togaInactMutLoad(row);
            printf("<tr>\n");
            printf("<td>%s</td>\n", info->exon_num);
            printf("<td>%s</td>\n", info->position);
            printf("<td>%s</td>\n", info->mut_class);
            printf("<td>%s</td>\n", info->mutation);
            if (sameWord(info->is_inact, ONE_)){
                printf("<td>%s</td>\n", YES_);
            } else {
                printf("<td>%s</td>\n", NO_);
            }
            printf("<td>%s</td>\n", info->mut_id);
            printf("</tr>\n");
            togaInactMutFree(&info);
            """
            position = inact_mut_data[0]
            mut_class  = inact_mut_data[1]
            mutation = inact_mut_data[2]
            is_inact = "YES" if inact_mut_data[3] == "1" else "NO"
            mut_id = inact_mut_data[4]
            lines = []
            lines.append("<tr>")
            lines.append(f"<td>{exon_num}</td>")
            lines.append(f"<td>{position}</td>")
            lines.append(f"<td>{mut_class}</td>")
            lines.append(f"<td>{mutation}</td>")
            lines.append(f"<td>{is_inact}</td>")
            lines.append(f"<td>{mut_id}</td>")
            lines.append("</tr>")
            inact_mut_line = "".join(lines)
            t_lines.append(inact_mut_line)
        t_line = "".join(t_lines)
        trans_to_inact_line[t] = t_line
    return trans_to_inact_line


def __get_ts_to_toga_exons(toga_exons_file):
    trans_to_exon_line = {}
    """Table schema
    
    transcript varchar(50) not null,  -- unique projection ID: ${ref_transcript_ID}.${chain_ID}
    exon_num int unsigned not null,  -- exon number
    exon_region varchar(100) not null,  -- region where exon was detected
    pid float not null,  -- nucleotide %id
    blosum float not null,  -- normalized blosum score
    gaps tinyint unsigned not null,  -- are there any asm gaps near? 1 - yes 0 - no
    ali_class varchar(4) not null,  -- alignemnt class: A, B, C, A+
    exp_region varchar(50) not null,  -- where exon was expected
    in_exp_region tinyint unsigned not null,  -- detected in expected region or not 1 yes 0 no
    is_del_or_no char -- D - deleted, M - missing, I - intact
    alignment longblob not null  -- exon sequence in query
    """
    ts_to_exon_to_data = defaultdict(dict)
    f = open(toga_exons_file, "r")
    for line in f:
        ld = line.rstrip().split("\t")
        t = ld[0]
        ex_num = int(ld[1])
        datum = ld[2:]
        ts_to_exon_to_data[t][ex_num] = datum
    f.close()

    
    for transcript in ts_to_exon_to_data.keys():
        exons_unsorted = ts_to_exon_to_data[transcript].keys()
        exon_nums = sorted(exons_unsorted)
        transcript_lines = []
        for exon_num in exon_nums:
            exon_lines = []
            datum = ts_to_exon_to_data[transcript][exon_num]
            """To convert this:
            struct togaNucl *info = NULL;
            info = togaNuclLoad(row);
            printf("<h5>Exon number: %s</h5><BR>\n", info->exon_num);
            printf("<B>Exon region:</B> %s<BR>\n", info->exon_region);
            printf("<B>Nucleotide percent identity:</B> %s | <B>BLOSUM:</B> %s <BR>\n", info->pid, info->blosum);
            if (sameWord(info->gaps, ONE_)){
                printf("<B>Intersects assembly gaps:</B> %s<BR>\n", YES_);
            } else {
                printf("<B>Intersects assembly gaps:</B> %s<BR>\n", NO_);
            }
            printf("<B>Exon alignment class:</B> %s<BR>\n", info->ali_class);
            if (sameWord(info->in_exp_region, ONE_)){
                printf("<B>Detected within expected region (%s):</B> %s<BR>\n", info->exp_region, YES_);
            } else {
                printf("<B>Detected within expected region (%s):</B> %s<BR>\n", info->exp_region, NO_);
            }
            // printf("<B>Expected region:</B> %s<BR>\n", info->exp_region);
            printf("<BR>\n");
            printf("<B>Sequence alignment between reference and query exon:</B><BR>\n");
            printf("%s<BR>\n", info->alignment);
            togaNuclFree(&info);
            """
            exon_region = datum[0]
            pid = datum[1]
            blosum = datum[2]
            gaps = "YES" if datum[3] == "1" else "NO"
            ali_class = datum[4]
            exp_reg = datum[5]
            in_exp_region = "YES" if datum[6] == "1" else "NO"
            is_del_or_no = datum[7]
            ali = datum[8]

            if is_del_or_no == "D":
                exon_lines.append(f"<h5>Exon number: {exon_num} - Deleted</h5><BR>")
            elif is_del_or_no == "M":
                exon_lines.append(f"<h5>Exon number: {exon_num} - Missing</h5><BR>")
            else:
                exon_lines.append(f"<h5>Exon number: {exon_num}</h5><BR>")

            exon_lines.append(f"<B>Exon region:</B> {exon_region}</B><BR>")
            exon_lines.append(f"<B>Nucleotide percent identity:</B> {pid} | <B>BLOSUM:</B> {blosum} <BR>")
            exon_lines.append(f"<B>Intersects assembly gaps:</B> {gaps}<BR>")
            exon_lines.append(f"<B>Exon alignment class:</B> {ali_class}<BR>")
            exon_lines.append(f"<B>Detected within expected region ({exp_reg}):</B> {in_exp_region}<BR>")
            exon_lines.append(f"<BR>")
            exon_lines.append(f"<B>Sequence alignment between reference and query exon:</B><BR>")
            exon_lines.append(f"{ali}<BR>")
            exon_line = "".join(exon_lines)
            transcript_lines.append(exon_line)
        transcript_line = "".join(transcript_lines)
        trans_to_exon_line[transcript] = transcript_line
    return trans_to_exon_line


def __get_ref_transcript(q_trans):
    return ".".join(q_trans.split(".")[:-1])


def get_ref_ts_to_link(p_dir):
    """Get ref_transcript: link dict."""
    ret = {}
    ref_2bit_link = os.path.join(p_dir, "t2bit.link")
    ref_2bit_path = os.readlink(ref_2bit_link)
    ref_2bit_basename = os.path.basename(ref_2bit_path)
    ref_genome_gname = ".".join(ref_2bit_basename.split(".")[:-1])
    ref_toga_dir = f"/projects/hillerlab/genome/gbdb-HL/{ref_genome_gname}/TOGA/"
    ref_ts_to_link_path = os.path.join(ref_toga_dir, "toga.transcript2Link.txt")
    # dirname = os.path.dirname(ref_2bit_path)
    # chrom_sizes_path = os.path.join(dirname, "chrom.sizes")
    if not os.path.isfile(ref_ts_to_link_path):
        return ret, ref_genome_gname
    f = open(ref_ts_to_link_path, "r")
    for line in f:
        ld = line.rstrip().split("\t")
        ret[ld[0]] = ld[1]
    f.close()
    return ret, ref_genome_gname


def get_que_chrom_sizes(p_dir):
    que_2bit_link = os.path.join(p_dir, "q2bit.link")
    que_2bit_path = os.readlink(que_2bit_link)
    dirname = os.path.dirname(que_2bit_path)
    chrom_sizes_path = os.path.join(dirname, "chrom.sizes")
    que_2bit_basename = os.path.basename(que_2bit_path)
    que_genome_gname = ".".join(que_2bit_basename.split(".")[:-1])
    return chrom_sizes_path, que_genome_gname


def merge_all_tables(project_dir, r_ts_to_link):
    """Merge all tabs in ready-to-merge bed file."""
    print("Joining tsv for bigbed")
    tabs_dir = os.path.join(project_dir, "tabs")
    sorted_bed = os.path.join(tabs_dir, "query_annot_sorted.bed")
    toga_data_file = os.path.join(tabs_dir, "togaData.tab")
    toga_inact_file = os.path.join(tabs_dir, "togaInactMut.tab")
    toga_nucl_file = os.path.join(tabs_dir, "togaNucl.tab")
    transcripts_ordered, t_to_bed = __get_ts_from_bed(sorted_bed)
    ts_to_data_track = __get_ts_to_toga_data(toga_data_file)
    ts_to_inact_data = __get_ts_to_toga_inact(toga_inact_file)
    ts_to_exon_data = __get_ts_to_toga_exons(toga_nucl_file)
    

    bigbed_material = os.path.join(tabs_dir, "query_annot_for_big.tsv")
    f = open(bigbed_material, "w")
    for transcript in transcripts_ordered:
        ref_trans = __get_ref_transcript(transcript)
        bed_fields = t_to_bed[transcript]
        data_track = ts_to_data_track[transcript]
        inact_line = ts_to_inact_data.get(transcript, "<BR>")
        exon_line = ts_to_exon_data[transcript]
        link = r_ts_to_link.get(ref_trans, ref_trans)  # just reference transcript if no link
        combined_datum_lst = bed_fields + data_track + [link] + [inact_line] + [exon_line]
        combined_datum = "\t".join(combined_datum_lst)
        f.write(f"{combined_datum}\n")
    f.close()


def sort_and_make_bb(project_dir, chrom_sizes):
    """bed12+22 as the result
    sort -k1,1 -k2,2n query_annot_for_big.tsv > query_annot_for_big_sorta.tsv
    bedToBigBed -type=bed12+22 query_annot_for_big_sorta.tsv
    /projects/hillerlab/genome/gbdb-HL/mm10/chrom.sizes query_annot.bb -tab
    """
    print("Making bigbed...")
    tabs_dir = os.path.join(project_dir, "tabs")
    bigbed_material_not_sorted = os.path.join(tabs_dir, "query_annot_for_big.tsv")
    bigbed_material = os.path.join(tabs_dir, "query_annot_for_big__sorted.tsv")
    sort_cmd = f"sort -k1,1 -k2,2n {bigbed_material_not_sorted} > {bigbed_material}"
    rc = subprocess.call(sort_cmd, shell=True)
    if rc != 0:
        sys.exit(f"Error! Command:\n{sort_cmd}\ncrashed")
    output = os.path.join(tabs_dir, "query_annotation.bb")
    # now transform into bigbed
    print("Calling:")
    bb_cmd = f"bedToBigBed -type=bed12+22 {bigbed_material} {chrom_sizes} {output} -tab -extraIndex=name -as={SCHEMA_LOCATION}"
    print(bb_cmd)
    rc = subprocess.call(bb_cmd, shell=True)
    if rc != 0:
        sys.exit(f"Error! Command:\n{bb_cmd}\ncrashed")

    # make ix.txt
    query_annot = os.path.join(project_dir, "query_annotation.bed")
    make_ix_script = os.path.join(LOCATION, "get_names_from_bed.py")
    ix_txt = os.path.join(tabs_dir, "query_annotation.ix.txt")
    ix_cmd = f"{make_ix_script} {query_annot} | sort -u > {ix_txt}"
    rc = subprocess.call(ix_cmd, shell=True)
    if rc != 0:
        sys.exit(f"Error! Command:\n{ix_cmd}\ncrashed")
    
    # make .ix and .ixx indexes
    ix_path = os.path.join(tabs_dir, "query_annotation.bb.ix")
    ixx_path = os.path.join(tabs_dir, "query_annotation.bb.ixx")
    ixixx_cmd = f"ixIxx {ix_txt} {ix_path} {ixx_path}"
    rc = subprocess.call(ixixx_cmd, shell=True)
    if rc != 0:
        sys.exit(f"Error! Command:\n{ixixx_cmd}\ncrashed")


def cleanup(project_dir, dont):
    """Clean project dir up."""
    if dont:  # do not cleanup
        return
    tabs_dir = os.path.join(project_dir, "tabs")
    bigbed_material_not_sorted = os.path.join(tabs_dir, "query_annot_for_big.tsv")
    bigbed_material = os.path.join(tabs_dir, "query_annot_for_big__sorted.tsv")
    toga_info_file = os.path.join(tabs_dir, "togaInfo.tab")
    toga_inact_feat_file = os.path.join(tabs_dir, "togaInactFeat.tab")
    toga_plot_file = os.path.join(tabs_dir, "togaPlot.tab")
    toga_prot_file = os.path.join(tabs_dir, "togaProt.tab")
    sorted_bed = os.path.join(tabs_dir, "query_annot_sorted.bed")
    toga_data_file = os.path.join(tabs_dir, "togaData.tab")
    toga_inact_file = os.path.join(tabs_dir, "togaInactMut.tab")
    toga_nucl_file = os.path.join(tabs_dir, "togaNucl.tab")
    to_rm = [bigbed_material_not_sorted, bigbed_material, toga_info_file, 
             toga_inact_feat_file, toga_plot_file, toga_prot_file,
             sorted_bed, toga_data_file, toga_inact_file, toga_nucl_file]
    for path in to_rm:
        os.remove(path) if os.path.isfile(path) else None


def __ssh_mkdir(uname, dirname):
    cmd = f"ssh {uname}@genome mkdir -p {dirname}"
    print(f"Calling {cmd}")
    subprocess.call(cmd, shell=True)


def __ssh_ln(uname, src, dest):
    """Call ln on delta"""
    cmd = f"ssh {uname}@genome ln -sf {src} {dest}"
    print(f"Calling {cmd}")
    subprocess.call(cmd, shell=True)


def __get_ccase_refname(refname):
    if refname == "hg38":
        return "Hg38"
    elif refname == "mm10":
        return "Mm10"
    elif refname == "mm39":
        return "Mm39"
    else:
        return refname


def load_bb_track(project_dir, dont, ref, que, bb_ver):
    if dont:  # do not load anything
        return
    tabs_dir = os.path.join(project_dir, "tabs")
    bb_file = os.path.join(tabs_dir, "query_annotation.bb")
    uname = getpass.getuser()
    # scp query_annotation.bb genome:/genome/gbdb-HL/$ref/TOGA/vs_$query/HLTOGAannotVs$refv1.bb
    ref_dir__genome = f"/genome/gbdb-HL/{ref}"
    toga_ref_dir__genome = f"{ref_dir__genome}/TOGA"
    # make if needed
    __ssh_mkdir(uname, toga_ref_dir__genome)
    vs_que_dir__genome = f"{toga_ref_dir__genome}/vs_{que}"
    __ssh_mkdir(uname, vs_que_dir__genome)
    # transfer data
    bb_file = os.path.join(tabs_dir, "query_annotation.bb")
    ix_file = os.path.join(tabs_dir, "query_annotation.bb.ix")
    ixx_file = os.path.join(tabs_dir, "query_annotation.bb.ixx")

    ref_camel_case = __get_ccase_refname(ref)
    bb_filename__genome = f"HLTOGAannotVs{ref_camel_case}{bb_ver}.bb"
    ix_filename__genome = f"HLTOGAannotVs{ref_camel_case}{bb_ver}.ix"
    ixx_filename__genome = f"HLTOGAannotVs{ref_camel_case}{bb_ver}.ixx"

    bb_path_genome = f"{vs_que_dir__genome}/{bb_filename__genome}"
    ix_path_genome = f"{vs_que_dir__genome}/{ix_filename__genome}"
    ixx_path_genome = f"{vs_que_dir__genome}/{ixx_filename__genome}"

    rsync_cmd = f"rsync -av {bb_file} {uname}@genome:{bb_path_genome}"
    subprocess.call(rsync_cmd, shell=True)

    rsync_cmd = f"rsync -av {ix_file} {uname}@genome:{ix_path_genome}"
    subprocess.call(rsync_cmd, shell=True)

    rsync_cmd = f"rsync -av {ixx_file} {uname}@genome:{ixx_path_genome}"
    subprocess.call(rsync_cmd, shell=True)

    # make link in /var/www
    www_data__genome = "/var/www/data"
    www_que__genome = f"{www_data__genome}/{que}"
    __ssh_mkdir(uname, www_que__genome)

    www_bb_link_dest = f"{www_que__genome}/{bb_filename__genome}"
    www_ix_link_dest = f"{www_que__genome}/{ix_filename__genome}"
    www_ixx_link_dest = f"{www_que__genome}/{ixx_filename__genome}"

    __ssh_ln(uname, bb_path_genome, www_bb_link_dest)
    __ssh_ln(uname, ix_path_genome, www_ix_link_dest)
    __ssh_ln(uname, ixx_path_genome, www_ixx_link_dest)


def main():
    """Entry point."""
    args = parse_args()
    r_ts_to_link, ref_name = get_ref_ts_to_link(args.project_dir)
    chrom_sizes, que_name = get_que_chrom_sizes(args.project_dir)
    # call generate_tab_files.py and make_togaPlot.py
    call_gen_tab_files(args.project_dir)
    call_gen_plot_files(args.project_dir, no_plots=args.no_plots)
    # merge 4 transcriptID-related tables into a single one
    make_toga_data(args.project_dir)

    # OK, now making bigBed
    # togaNucl.tab and togaInactMut.tab
    # -> merge into HTML_formatted strings
    # -> add to togaData
    # save to as.tab for bedToBigBed
    make_sorted_bed(args.project_dir)

    merge_all_tables(args.project_dir, r_ts_to_link)
    sort_and_make_bb(args.project_dir, chrom_sizes)

    cleanup(args.project_dir, args.do_not_cleanup)


if __name__ == "__main__":
    main()
