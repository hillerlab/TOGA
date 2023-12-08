#!/usr/bin/env python3
"""Plot TOGA-detected mutations.

Original script provided by Björn Langer, 2014.
Adapted for TOGA by Bogdan Kirilenko, 2020.
"""
import argparse
import sys
import os
import re
from copy import copy
from collections import defaultdict
import h5py

# Mutation classes
EX_DEL = "Deleted exon"
EX_MIS = "Missing exon"
START_MIS = "START_MISSED"
FS_DEL = "FS_DEL"
FS_INS = "FS_INS"
BIG_DEL = "BIG_DEL"
BIG_INS = "BIG_INS"
INDELS = {FS_DEL, FS_INS, BIG_DEL, BIG_INS}
INS = {FS_INS, BIG_INS}
DELS = {FS_DEL, BIG_DEL}
STOP = "STOP"
SSM = "SSM"
# (ag)acceptor-EXON-donor(gt)
SSM_D = "SSMD"  # Donor, right, GT,GC
SSM_A = "SSMA"  # Acceptor, left, AG 
SSM_TYPES = {SSM_D, SSM_A}

COMP = "COMPENSATION"

TEMPLATE_PATH_1 = "svg_template.txt"
TEMPLATE_PATH_2 = "supply/svg_template.txt"

OPACITY = 100
HORIZONTAL = "horizontal"
MUT_LINE_FIELDS = 8
MUT_LINE_FIELDS_SP = MUT_LINE_FIELDS

BLACK = "#121212"  # almost black
MISS_SEQ_COLOR = "#878787"  # Grey
MASKED_MUT_COLOR = "#878787"  # Grey
INACT_MUT_COLOR = "#cf232b"  # Dark Red
FP_EXON_DEL_COLOR = "#0c7bdc"  # blue-ish
BACKGROUND_COLOR = "#87bcbc"  # cyan-ish
STOP_CODON_COLOR = "#121212"  # almost black

MAX_VERTICAL_SPACE = 100  # pixels
HALF_EXON_HEIGHT = 15  # pixels
HALF_UTR_HEIGHT = 5  # pixels
UTR_WIDTH = 0  # pixels

EXON_BASE_SIZE = 0.8  # pixels per base
INTRON_BASE_SIZE = 0.02  # pixels per base
MAX_INTRON_SIZE = 2000  # bases
MIN_INTRON_SIZE = 2000  # bases

GAP_WIDTH = 0  # pixels
HALF_GAP_HEIGHT = 0  # pixels

MIN_ARROW_SIZE = 10  # pixels
MAX_ARROW_SIZE = 20  # pixels
ARROW_SIZE = 1  # pixels per base

PRINT_EXON_SCHEME = False
ROUND_EDGES = False
INTRON_STYLE = 'style="stroke:#999; stroke-width:3;" '
EXON_ANC_STYLE = "stroke-width:3;"
EXON_NON_ANC_STYLE = "stroke: black; stroke-width:3; stroke-dasharray: 5,5;"
INSERT_STYLE = 'style="fill:{0}; stroke-opacity:1; fill-opacity:1"'
DEL_STYLE = 'style="stroke:{0}; stroke-width:{1}; stroke-opacity:1"'
MISSSEQ_STYLE = 'style="stroke:{0}; stroke-width:{1}; stroke-opacity:1"'
STOP_CODON_STYLE = 'style="stroke:{0};stroke-width:3;"'
COMP_INDEL_STYLE = 'style="fill:none;stroke-width:1;stroke:green;"'
TEXT_HIGHLIGHT_STYLE = f'style="fill:{INACT_MUT_COLOR};"'
FONTFAMILY = "Courier New"
STOP_LABEL_FONTSIZE = 18
SS_LABEL_FONTSIZE = 18
MO_FONTSIZE = 18
STACKING_THRESHOLD = 35  # pixels
J_STYLE = False
picture_width = 0


SVG_PLACEHOLDER = """
<svg>
<rect fill="#aaa" stroke="#000" x="0" y="0" width="400" height="100"/>
<line x1="0" y1="0" x2="400" y2="100" stroke="red" stroke-width="4" />
<line x1="0" y1="100" x2="400" y2="0" stroke="red" stroke-width="4" />
</svg>
"""


DEFAULT_MODE = "DEFAULT"
PUB_MODE_I = "PUB_MODE_I"
DRAW_MODES = {
    DEFAULT_MODE,
    PUB_MODE_I,
}

# drawing modes where we don't show compensations
# for now only one, can be extended
DRAW_MODES_NO_COMPENSATIONS_SHOWN = {PUB_MODE_I}


class MouseOver:
    id_counter = 0

    def __init__(self, pos, texts):
        """Probably doesn't really work with TOGA output!"""
        MouseOver.id_counter += 1
        self.id_counter = MouseOver.id_counter
        lens_ = [
            len(x) - 4 * x.count("<cc>") for x in texts if not x.startswith("http://")
        ]
        self.width = max(lens_) * 0.65 * MO_FONTSIZE
        self.height = MO_FONTSIZE * sum(
            [(1, 0)[t.startswith("http://")] for t in texts]
        )
        # align text horizontally and vertically to image boarder
        self.pos = (
            (pos[0], max(0, picture_width - self.width))[
                pos[0] + self.width > picture_width
            ],
            (pos[1], pos[1] + self.height)[pos[1] - self.height < 0],
        )
        self.texts = texts

    def __str__(self):
        if not J_STYLE:
            elements, i, linecounter = "", 0, 0
            while i < len(self.texts):
                y_ = self.pos[1] - linecounter * MO_FONTSIZE
                elem_list = [
                    f'<tspan x="{self.pos[0]}" y="{y_}" ',
                ]
                if self.texts[i].startswith("http://"):
                    # create link
                    elem_list.append(">")
                    xlink_ = self.texts[i].replace("&", "&amp;")
                    elem_list.append(f'<a xlink:href="{xlink_}" target="_blank">')
                    i += 1
                    elem_list.append("<tspan>")
                    tparts = self.texts[i].split("<cc>")
                    elem_list.append(tparts[0])
                    for j in range(1, len(tparts)):
                        s_1 = ("", TEXT_HIGHLIGHT_STYLE)[j % 2]  # B: wtf is that?
                        s_2 = tparts[j]
                        elem_list.append(f"</tspan><tspan {s_1}>{s_2}")
                    elem_list.append("</tspan></a>")
                else:
                    tparts = self.texts[i].split("<cc>")
                    start_index = tparts[0] == ""
                    s_1 = ("", TEXT_HIGHLIGHT_STYLE)[start_index % 2]
                    elem_list.append(f"{s_1}>")
                    elem_list.append(tparts[start_index])
                    for j in range(start_index + 1, len(tparts)):
                        s_1 = ("", TEXT_HIGHLIGHT_STYLE)[j % 2]
                        s_2 = tparts[j]
                        elem_list.append(f"</tspan><tspan {s_1}>{s_2}")
                elem_list.append("</tspan>")
                elements = "".join(elem_list)
                i += 1
                linecounter += 1
            out = [
                f'  <g id="Mouseover{self.id_counter}" visibility="hidden">\n',
            ]
            out.append(
                f'    <text style="font-size:{MO_FONTSIZE}px">{elements}</text>\n'
            )
            out.append("  </g>\n")
            out_str = "".join(out)
            return out_str
        else:
            return ""


class Textstack:
    def __init__(
        self,
        pos,
        direction,
        label,
        mouseover=None,
        fontsize=STOP_LABEL_FONTSIZE,
        style="",
        anchor="",
        color=None,
    ):
        self.pos = pos
        self.direction = direction
        self.size = fontsize
        self.labels = [(label, mouseover)]
        self.style = style
        self.color = [color]
        anchor_cond = anchor == "start" or anchor == "middle" or anchor == "end"
        self.anchor = f'text-anchor="{anchor}"' if anchor_cond else ""

    def add_label(self, label, y, mouseover=None, color=None):
        if "up" in self.direction:
            self.pos = (self.pos[0], min(self.pos[1], y))
        elif "down" in self.direction:
            self.pos = (self.pos[0], max(self.pos[1], y))
        self.labels.append((label, mouseover))
        if color is not None and self.color is not None:
            self.color.append(color)
        elif color is None and self.color is None:
            pass
        else:
            print("Error: either all or no color for a stack must be provide")

    def shift_vertical(self, distance):
        self.pos = (self.pos[0], self.pos[1] + distance)

    def toString(self):
        global picture_width  # TODO: needs to be refactored
        out_lines = []
        mouseovers = []
        if self.direction == "upwards" or self.direction == "downwards":
            print("deprecated text direction")
            # some hacks because 'writing-mode="tb"' and 'glyph-orientation-vertical:360'
            # are not supported by firefox yet
            self.labels = [(" ".join(label[0]), label[1]) for label in self.labels]
            length = len(self.labels) + sum([len(label[0]) for label in self.labels])
            offset = HALF_EXON_HEIGHT + 0.5 * self.size

            if self.direction == "upwards":
                self.labels = [
                    (label[0][::-1], label[1]) for label in self.labels
                ]  # reverse each label
                rotation = 90
                pos_of_ = self.pos[1] - offset
                pos_of_len_ = self.pos[1] - offset - length * self.size
                path = f"M{self.pos[0]} {pos_of_} V{pos_of_len_}"
            else:
                rotation = 270
                pos_of_ = self.pos[1] + offset
                pos_of_len_ = self.pos[1] + offset + length * self.size
                path = f"M{self.pos[0]} {pos_of_} V{pos_of_len_}"

            out_lines.append(
                f'  <path id="textpath{self.pos[0]}{self.direction}" d="{path}" fill="none"/>\n'
            )
            out_lines.append(f'  <text style="{self.style}" rotate="{rotation}">\n')
            out_lines.append(
                f'    <textPath xlink:href="#textpath{self.pos[0]}{self.direction}">\n'
            )

            for i in range(len(self.labels)):
                fill_ = (f"fill:{BLACK};", f"fill:{BLACK};")[i % 2]
                out_lines.append(
                    f'<tspan style="font-size:{self.size};{fill_}">{self.labels[i][0]} </tspan>'
                )
            out_lines.append("    </textPath>\n")
            height = 0.5 * self.size + length * self.size
            out_lines.append("  </text>\n")
            out_str = "".join(out_lines)
            return out_str, height
        else:
            out_lines.append(f'  <text style="{self.style}" {self.anchor}>\n')
            y = self.pos[1]
            c = 1
            if self.direction == "down":
                y += HALF_EXON_HEIGHT + 1.5 * self.size
            elif self.direction == "up":
                y -= HALF_EXON_HEIGHT + 0.5 * self.size
                c = -1
            for i in range(len(self.labels)):
                if self.labels[i][1] is not None:
                    mouseovers.append(
                        MouseOver((self.pos[0], y - self.size), self.labels[i][1])
                    )
                    picture_width = max(picture_width, mouseovers[-1].width)
                    mo_ref = f"onmouseover=\"mouseover(evt, 'Mouseover{MouseOver.id_counter}')\""
                else:
                    mo_ref = ""
                if self.color == [None]:
                    switch_ = (self.style, f"fill:{INACT_MUT_COLOR};")[
                        i % 2
                    ]  # B: don't understand this
                    label_ = self.labels[i][0]
                    out_lines.append(
                        f'<tspan x="{self.pos[0]}" y="{y}" style="font-size:{self.size}px;'
                    )
                    out_lines.append(f'{switch_}" {mo_ref}>{label_}</tspan>')
                else:
                    label_ = self.labels[i][0]
                    out_lines.append(
                        f'<tspan x="{self.pos[0]}" y="{y}" style="font-size:{self.size}px;'
                    )
                    out_lines.append(
                        f'fill:{self.color[i]};" {mo_ref}>{label_}</tspan>'
                    )
                y += c * self.size
            height = y - self.pos[1] - (0, self.size)[self.direction == "up"]

            mouseovers_join = "".join([str(mo) for mo in mouseovers])
            out_lines_merge = "".join(out_lines)
            outstr = f"{mouseovers_join}{out_lines_merge}  </text>\n"
            return outstr, -1 * height


class Exon:
    def __init__(
        self,
        pos,
        width,
        is_exon=True,
        opacity=OPACITY,
        color=None,
        ancestral=True,
        annotation=None,
        draw_mode=DEFAULT_MODE,
    ):
        self.pos = [pos[0], pos[1]]
        self.width = width
        self.half_height = HALF_EXON_HEIGHT if is_exon else HALF_UTR_HEIGHT
        self.final_height = 2 * self.half_height
        self.opacity = opacity
        self.data = []
        self.labels = (
            {}
        )  # keys: '(direction, base_number)'   values: '(tag, y-value, mouseover)'
        self.mouseover = ""
        self.annotation = annotation
        self.ancestral = ancestral
        self.draw_mode = draw_mode

        if color is None and ancestral:
            self.color = BACKGROUND_COLOR
        elif color is None:
            self.color = "none"
        else:
            self.color = color

    def __copy__(self):
        """Copy method implementation."""
        pos_copy = [self.pos[0], self.pos[1]]
        new_one = Exon(
            pos_copy,
            self.width,
            opacity=self.opacity,
            color=self.color,
            ancestral=self.ancestral,
            annotation=self.annotation,
            draw_mode=self.draw_mode,
        )
        return new_one

    def shift_vertical(self, distance):
        self.pos[1] += distance

    def cluster_labels(self):
        keys = self.labels.keys()  # keys have the form '(direction, base_number)'
        height = self.half_height

        if keys:
            keys = sorted(keys)
            first_key, last_key, cluster = keys[0], keys[0], []
            for k in keys:
                if k[1] > first_key[1] + STACKING_THRESHOLD or k[0] != first_key[0]:
                    cluster_center = (
                        self.pos[0] + (first_key[1] + last_key[1]) / 2 * EXON_BASE_SIZE
                    )
                    t = Textstack(
                        (cluster_center, self.labels[cluster[0]][1]),
                        first_key[0],
                        self.labels[cluster[0]][0],
                        mouseover=self.labels[cluster[0]][2],
                        anchor="middle",
                        color=self.labels[cluster[0]][3],
                    )
                    for i in range(1, len(cluster)):
                        t.add_label(
                            self.labels[cluster[i]][0],
                            self.labels[cluster[i]][1],
                            self.labels[cluster[i]][2],
                            color=self.labels[cluster[i]][3],
                        )

                    s, h = t.toString()
                    self.data.append(s)
                    # print self.labels[cluster[0]][1], self.pos[1], h
                    height = max(height, self.pos[1] - self.labels[cluster[0]][1] + h)
                    first_key, last_key, cluster = k, k, [k]
                else:
                    cluster.append(k)
                    last_key = k
            cluster_center = (
                self.pos[0] + (first_key[1] + last_key[1]) / 2 * EXON_BASE_SIZE
            )

            t = Textstack(
                (cluster_center, self.labels[cluster[0]][1]),
                first_key[0],
                self.labels[cluster[0]][0],
                mouseover=self.labels[cluster[0]][2],
                anchor="middle",
                color=self.labels[cluster[0]][3],
            )
            for i in range(1, len(cluster)):
                t.add_label(
                    self.labels[cluster[i]][0],
                    self.labels[cluster[i]][1],
                    self.labels[cluster[i]][2],
                    color=self.labels[cluster[i]][3],
                )
            s, h = t.toString()
            self.data.append(s)
            height = max(height, self.pos[1] - self.labels[cluster[0]][1] + h)
        self.height = height

    def add_mouseover(self, labels):
        self.mouseover = MouseOver(
            (self.pos[0] + self.width / 2, self.pos[1] - 1.1 * self.half_height), labels
        )

    def add_stop_codon(self, start_base_pos, tag, mouseover, is_masked):
        x = self.pos[0] + (start_base_pos + 1) * EXON_BASE_SIZE
        y = self.pos[1]
        color_ = MASKED_MUT_COLOR if is_masked else STOP_CODON_COLOR
        style_ = STOP_CODON_STYLE.format(color_)
        self.data.append(
            draw_line((x, y - HALF_EXON_HEIGHT), (x, y + HALF_EXON_HEIGHT), style_)
        )
        self.labels[("up", start_base_pos)] = (tag, y, mouseover, color_)

    def add_deletion(self, start_base_pos, length, mouseover, is_masked):
        if self.draw_mode == PUB_MODE_I:
            # all mutations in black
            color = BLACK
        else:
            color = MASKED_MUT_COLOR if is_masked else INACT_MUT_COLOR
        style = DEL_STYLE.format(color, length * EXON_BASE_SIZE)
        self.data.append(
            draw_line(
                (
                    self.pos[0] + (start_base_pos + length / 2) * EXON_BASE_SIZE,
                    self.pos[1] - HALF_EXON_HEIGHT,
                ),
                (
                    self.pos[0] + (start_base_pos + length / 2) * EXON_BASE_SIZE,
                    self.pos[1] + HALF_EXON_HEIGHT,
                ),
                style,
            )
        )
        self.labels[("up", start_base_pos)] = (
            "-" + str(length),
            self.pos[1],
            mouseover,
            color,
        )

    def add_missing_sequence(self, start_base_pos, length, mouseover):
        style = MISSSEQ_STYLE.format(MISS_SEQ_COLOR, length * EXON_BASE_SIZE)
        self.data.append(
            draw_line(
                (
                    self.pos[0] + (start_base_pos + length / 2) * EXON_BASE_SIZE,
                    self.pos[1] - HALF_EXON_HEIGHT,
                ),
                (
                    self.pos[0] + (start_base_pos + length / 2) * EXON_BASE_SIZE,
                    self.pos[1] + HALF_EXON_HEIGHT,
                ),
                style,
            )
        )
        self.labels[("up", start_base_pos)] = (
            str(length),
            self.pos[1],
            mouseover,
            MISS_SEQ_COLOR,
        )

    def add_insertion(self, base_pos, length, mouseover, is_masked):
        if self.draw_mode == PUB_MODE_I:
            color = BLACK
        else:
            color = MASKED_MUT_COLOR if is_masked else INACT_MUT_COLOR
        style = INSERT_STYLE.format(color)
        x = self.pos[0] + base_pos * EXON_BASE_SIZE
        y = self.pos[1] - HALF_EXON_HEIGHT
        height = min(max(length * ARROW_SIZE, MIN_ARROW_SIZE), MAX_ARROW_SIZE)
        # B: get points coordinates, for point 1 it's pretty simple: x,y
        width = height / 2
        p2x_ = x + width / 2
        p2y_ = y - height
        p3x_ = x - width / 2
        p3y_ = y - height
        self.data.append(
            f'  <polygon points="{x},{y} {p2x_},{p2y_} {p3x_},{p3y_}" {style}/>\n'
        )
        plus_len_ = f"+{str(length)}"
        self.labels[("up", base_pos)] = (
            plus_len_,
            self.pos[1] - height,
            mouseover,
            color,
        )

    def add_altSS(self, offset, width, style):
        rect_x_ = self.pos[0] + offset
        rect_y_ = self.pos[1] - self.half_height
        height_ = 2 * self.half_height  # B: I don't get it
        to_insert = f'  <rect x="{rect_x_}" y="{rect_y_}" width ="{width}" height="{height_}" style="{style}"/>\n'
        self.data.insert(0, to_insert)

    def toString(self):
        self.cluster_labels()
        if self.mouseover == "":
            attr_lines = []
        else:
            # mo_id_str_ = str(self.mouseover.id_counter)
            mo_id_str_ = ""
            attr_lines = [
                f"onmouseover=\"mouseover(evt, 'Mouseover{mo_id_str_}')\"",
            ]

        extra_attr = f' rx="{HALF_EXON_HEIGHT / 3}" ry ="{HALF_EXON_HEIGHT / 3}" '  # B: don't know what it stands for
        attr_lines.append(("", extra_attr)[ROUND_EDGES])
        if self.annotation is not None:
            mouseover = MouseOver(
                (self.pos[0] + self.width / 2, self.pos[1] + 1.1 * self.half_height),
                self.annotation,
            )
            mo_str_ = str(mouseover)
            dline_str = f'style="stroke:white; stroke-width:5;" onmouseover="mouseover(evt, \'Mouseover{mo_str_}\')"'
            dline_ = draw_line(
                [self.pos[0], self.pos[1] + self.half_height],
                [self.pos[0] + self.width, self.pos[1] + self.half_height],
                dline_str,
            )
            annotation = str(mouseover) + dline_
        else:
            annotation = ""

        # time to assemble the result
        self_mo_str_ = str(self.mouseover)
        rect_class_ = ("non_", "")[self.ancestral]
        x_ = self.pos[0]
        y_ = self.pos[1] - self.half_height
        height_ = 2 * self.half_height
        fill_opac_ = self.opacity / 100.0
        attributes = "".join(attr_lines)
        data_str = "".join(self.data)
        rect_line = (
            f'  <rect class="{rect_class_}anc_exon" x="{x_}" y="{y_}" width="{self.width}" '
            f'height="{height_}" style="fill:{self.color};fill-opacity:{fill_opac_:1.2f}" {attributes}/>\n'
        )
        result_line = f"{annotation}{self_mo_str_}{rect_line}{data_str}"
        return result_line, self.height


class Intron:
    def __init__(self, pos, width, seperated=False):
        self.pos = [pos[0], pos[1]]
        self.width = width
        self.seperated = seperated

    def shift_vertical(self, distance):
        self.pos[1] += distance

    def toString(self):
        if not self.seperated:
            line = draw_line(
                self.pos, (self.pos[0] + self.width, self.pos[1]), INTRON_STYLE
            )
        else:
            gap_l = self.pos[0] + self.width // 2 - GAP_WIDTH // 2
            gap_r = gap_l + GAP_WIDTH
            line_pieces = []
            line_pieces.append(draw_line(self.pos, (gap_l, self.pos[1]), INTRON_STYLE))
            line_pieces.append(
                draw_line(
                    (gap_l, self.pos[1] - HALF_GAP_HEIGHT),
                    (gap_l, self.pos[1] + HALF_GAP_HEIGHT),
                    INTRON_STYLE,
                )
            )
            line_pieces.append(
                draw_line(
                    (gap_r, self.pos[1] - HALF_GAP_HEIGHT),
                    (gap_r, self.pos[1] + HALF_GAP_HEIGHT),
                    INTRON_STYLE,
                )
            )
            line_pieces.append(
                draw_line(
                    (gap_r, self.pos[1]),
                    (self.pos[0] + self.width, self.pos[1]),
                    INTRON_STYLE,
                )
            )
            line = "".join(line_pieces)
        return line, 2 * HALF_GAP_HEIGHT

    def __copy__(self):
        upd_pos = [self.pos[0], self.pos[1]]
        return Intron(upd_pos, self.width, self.seperated)


def draw_line(pos1, pos2, otherAttributes):
    return f'  <line x1="{pos1[0]}" y1="{pos1[1]}" x2="{pos2[0]}" y2="{pos2[1]}" {otherAttributes}/>\n'


def draw_comp_indel(pos1, pos2):
    center = " ".join(
        2 * [str((pos1[0] + pos2[0]) // 2), str(pos1[1] - HALF_UTR_HEIGHT)]
    )
    d = f"M{pos1[0]} {pos1[1] - HALF_EXON_HEIGHT} C {center} {pos2[0]} {pos2[1] - HALF_EXON_HEIGHT}"
    return f'  <path d="{d}" {COMP_INDEL_STYLE} />\n'


def parse_args():
    """Read arguments."""
    app = argparse.ArgumentParser()
    app.add_argument("bed_file", help="Reference annotation, bed12 file")
    app.add_argument("toga_mut_file", help="TOGA mutations")
    app.add_argument(
        "transcripts", help="Transcripts of interest (comma-sep list or a single one)"
    )
    app.add_argument("output", help="Path to output image")
    app.add_argument("--chain", default=None, help="Limit to a particular projection")
    app.add_argument(
        "--isoforms_file",
        "-i",
        default=None,
        help="If you provide geneID instead of transcript, please also provide "
        "an isoforms file.",
    )
    app.add_argument(
        "--multi_species",
        "--ms",
        action="store_true",
        dest="multi_species",
        help="Please see documentation which is not yet written",
    )
    app.add_argument("--alt_template", "-a", help="Alternative template path")
    app.add_argument(
        "--publication_mode_heni",
        "--pmh",
        help="Publication mode; Indrischek et al. 2021 version",
        dest="publication_mode_heni",
        action="store_true",
    )
    app.add_argument("--verbose", "-v", action="store_true", dest="verbose")
    if len(sys.argv) < 4:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    if not os.path.isfile(args.bed_file):
        sys.exit(f"Error! Bed file {args.bed_file} not found")
    elif not os.path.isfile(args.toga_mut_file):
        sys.exit(f"Error! Mut file {args.toga_mut_file} not found")
    return args


def get_transcript(bed_file, trans):
    """Find the requested transcript."""
    if bed_file.endswith(".hdf5"):
        # TODO: smarter check for hdf5 file
        h = h5py.File(bed_file, "r")
        try:
            b_bed_file = h[str(trans)][()]
            u_type = f"U{len(b_bed_file)}"
            line_req = b_bed_file.astype(u_type)
        except KeyError:
            line_req = None
        h.close()
        return line_req
    else:
        f = open(bed_file, "r")
        line_req = None
        for line in f:
            if line.split("\t")[3] == trans:
                line_req = line.rstrip()
                break
        f.close()
        return line_req


def get_mut_list(mut_file, trans):
    """Get a list of mutations."""
    if mut_file.endswith(".hdf5"):
        os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"  # otherwise hdf5
        # TODO: a smarter check for hdf5 file
        h = h5py.File(mut_file, "r")
        try:
            b_mut_lines = h[trans][()]
            u_type = f"U{len(b_mut_lines)}"
            trans_lines_str = b_mut_lines.astype(u_type)
            trans_related_lines = trans_lines_str.split("\n")
        except KeyError:
            trans_related_lines = [
                "",
            ]
        h.close()
        return trans_related_lines
    else:  # ordinary text file
        f = open(mut_file, "r")
        trans_related_lines = []
        for line in f:
            if not line.startswith("#"):
                continue
            line_data = line.rstrip().split("\t")
            trans_line = line_data[0][2:]
            if trans_line == trans:
                trans_related_lines.append(line)
        f.close()
        return trans_related_lines


def get_rel_starts(block_sizes):
    """Like accumulated sum."""
    ans = [0]
    for block in block_sizes[:-1]:
        val = ans[-1] + block
        ans.append(val)
    return ans


def parse_trans_data(bed_line):
    """Parse trans data."""
    line_data = bed_line.rstrip().split("\t")
    # get rid of UTR's
    chrom = line_data[0]
    chromStart = int(line_data[1])
    _ = int(line_data[2])
    name = line_data[3]  # gene_name usually
    _ = int(line_data[4])  # never used
    strand = line_data[5]
    thickStart = int(line_data[6])
    thickEnd = int(line_data[7])
    _ = line_data[8]  # never used
    blockCount = int(line_data[9])
    blockSizes = [int(x) for x in line_data[10].split(",") if x != ""]
    blockStarts = [int(x) for x in line_data[11].split(",") if x != ""]
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

    new_block_sizes = [x[1] - x[0] for x in zip(blockNewStarts, blockNewEnds)]
    if strand == "-":
        new_block_sizes = new_block_sizes[::-1]
    rel_starts = get_rel_starts(new_block_sizes)

    exons = list(zip(blockNewStarts, blockNewEnds))

    bed_data = {
        "chrom": chrom,
        "name": name,
        "strand": strand,
        "exons": exons,
        "rel_starts": rel_starts,
        "exon_num": len(exons),
    }
    return bed_data


def comp_width(trans_data):
    """Compute picture width."""
    exon_coords = trans_data["exons"]
    gene_width = (
        UTR_WIDTH
        + int((exon_coords[0][1] - exon_coords[0][0]) * EXON_BASE_SIZE)
        + max(UTR_WIDTH, 2 * SS_LABEL_FONTSIZE)
    )
    for i in range(1, trans_data["exon_num"]):
        intron_length = exon_coords[i][0] - exon_coords[i - 1][1]
        intron_length = min(max(intron_length, MIN_INTRON_SIZE), MAX_INTRON_SIZE)
        gene_width += int(intron_length * INTRON_BASE_SIZE)
        gene_width += int((exon_coords[i][1] - exon_coords[i][0]) * EXON_BASE_SIZE)
    return gene_width


def get_rel_pos(abs_pos, ex_num, rel_starts):
    """Convert absolute position to relativ."""
    if len(rel_starts) == 1:
        return abs_pos
    ex_num_0 = int(ex_num) - 1
    rel_pos_uncorr = abs_pos - rel_starts[ex_num_0]
    rel_pos = rel_pos_uncorr if rel_pos_uncorr >= 0 else 0
    return rel_pos


def get_comp_pairs(mut_ids):
    """Get pairs of compensated IDs."""
    ans = []
    # cannot be < 2
    if len(mut_ids) == 2:
        pair = (mut_ids[0], mut_ids[1])
        ans.append(pair)
        return ans
    # if 4: need 3 pairs:
    # (1-2) (2-3) (3-4)
    for i in range(len(mut_ids) - 1):
        curr_id = mut_ids[i]
        next_id = mut_ids[i+1]
        pair = (curr_id, next_id)
        ans.append(pair)
    return ans


def prepare_mut_data(mut_lines, sp, rel_starts, chain_lim=None):
    """Adjust mutations data.

    Absolute codon positions -> relative nucl.
    Split compensations in pairs.
    """
    proj_lines = defaultdict(list)
    for line in mut_lines:
        if not line.startswith("#"):
            continue
        line_data = line.rstrip().split("\t")
        upd_line_data = line_data.copy()
        # update trans, it was "# ENST0000..."
        trans = line_data[0].replace("# ", "")
        chain = line_data[1]
        if chain_lim and chain != chain_lim:
            # limit to a projection
            continue
        projection = f"{trans}.{chain}"
        upd_line_data[0] = trans
        if len(line.split("\t")) != MUT_LINE_FIELDS:
            to_merge = [str(x) for x in upd_line_data]
            upd_line = "\t".join(to_merge)
            proj_lines[projection].append(upd_line)
            continue
        exon_num = int(line_data[2])
        abs_codon_pos = int(line_data[3])
        mclass = line_data[4]
        if mclass in INDELS or mclass == STOP:
            # correct positioning
            abs_nucl_pos = (abs_codon_pos * 3) - 2
            rel_nucl_pos = str(get_rel_pos(abs_nucl_pos, exon_num, rel_starts))
            upd_line_data[3] = rel_nucl_pos
            to_merge = [str(x) for x in upd_line_data]
            upd_line = "\t".join(to_merge)
            proj_lines[projection].append(upd_line)
        elif mclass == COMP:
            # split compensations if required
            comp_field = line_data[5]
            try:
                # comp_ids = comp_field.split("_")[1].split(",")
                comp_ids_range_str = comp_field.split("_")[1].split("-")
                # fmt: FS_{start}-{end}
                _comp_start = int(comp_ids_range_str[0])
                _comp_end = int(comp_ids_range_str[1])
            except ValueError:
                # old-style formatting
                comp_ids_range_str = comp_field.split("_")[1].split(",")
                _comp_start = int(comp_ids_range_str[0])
                _comp_end = int(comp_ids_range_str[-1])
            comp_ids = list(range(_comp_start, _comp_end + 1))
            comp_pairs = get_comp_pairs(comp_ids)
            for pair in comp_pairs:
                upd_mark = f"FS_{pair[0]},{pair[1]}"
                pair_line = upd_line_data.copy()
                pair_line[5] = upd_mark
                upd_line = "\t".join([str(x) for x in pair_line])
                proj_lines[projection].append(upd_line)
        else:
            proj_lines[projection].append("\t".join([str(x) for x in upd_line_data]))

    # compensations must be in the end of each projection mutations list
    upd_lines = []
    for projection, lines in proj_lines.items():
        comp_lines, no_comp_lines = [], []
        for line in lines:
            line_split = line.split("\t")
            if len(line_split) < MUT_LINE_FIELDS:
                no_comp_lines.append(line)
            elif line.split("\t")[4] == COMP:
                comp_lines.append(line)
            else:
                no_comp_lines.append(line)
        proj_list = no_comp_lines + comp_lines
        upd_lines.extend(proj_list)

    to_ret = []
    for line in upd_lines:
        # add sp identifier
        # line_ret = f"{line}\t{sp}"
        line_ret = (line, sp)
        to_ret.append(line_ret)
    return to_ret


def write_pic(filebuffer, picture_width, picture_height, vbox_h, mcv, alt_template):
    """Write picture data."""
    this_file_loc = os.path.dirname(__file__)
    templ_3 = os.path.join(this_file_loc, "svg_template.txt")
    if os.path.isfile(TEMPLATE_PATH_1):
        template_file = TEMPLATE_PATH_1
    elif os.path.isfile(TEMPLATE_PATH_2):
        template_file = TEMPLATE_PATH_2
    elif alt_template:
        template_file = alt_template
    elif templ_3:
        template_file = templ_3
    else:
        sys.exit(f"Error! Cannot locate svg template")
    with open(template_file, "r") as f:
        template = f.read()

    # replace all placeholders with our values
    vbox_val = min(0, MAX_VERTICAL_SPACE / 2 - vbox_h)
    rep = {
        "FILE_BUFFER": str(filebuffer),
        "PICTURE_WIDTH": str(picture_width),
        "PICTURE_HEIGHT": str(picture_height),
        "VBOX": str(vbox_val),
        "EXON_ANC_STYLE": EXON_ANC_STYLE,
        "EXON_NON_ANC_STYLE": EXON_NON_ANC_STYLE,
        "STOP_LABEL_FONTSIZE": str(STOP_LABEL_FONTSIZE),
        "FONT_FAMILY": FONTFAMILY,
        "MOUSEOVERCOUNTER": str(mcv),
    }
    rep_d = {re.escape(k): v for k, v in rep.items()}
    pattern = re.compile("|".join(rep_d.keys()))
    text = pattern.sub(lambda m: rep[re.escape(m.group(0))], template)
    return text


def trans_data_direction(trans_data):
    """Adjust trans data direction depending on strand."""
    if trans_data["strand"] == "+":
        return
    max_coord = trans_data["exons"][-1][1]
    for i in range(trans_data["exon_num"]):
        trans_data["exons"][i] = (
            max_coord - trans_data["exons"][i][1],
            max_coord - trans_data["exons"][i][0],
        )
    trans_data["exons"].reverse()


def init_human_exons(x, y, ancestral_exons, trans_data, drawing_mode=DEFAULT_MODE):
    """Just init human exons list."""
    human_exons = [
        Exon((0, y), x[0], is_exon=False, draw_mode=drawing_mode),
    ]
    # 1 is always in ancestrals!
    # one_in_anc_exons = "1" in ancestral_exons
    human_exons.append(
        Exon(
            (x[0], y),
            x[1] - x[0],
            opacity=OPACITY,
            ancestral=True,
            annotation=None,
            draw_mode=drawing_mode,
        )
    )

    for i in range(1, trans_data["exon_num"]):
        intron_length = trans_data["exons"][i][0] - trans_data["exons"][i - 1][1]
        seperated = intron_length > MAX_INTRON_SIZE
        intron_length = min(max(intron_length, MIN_INTRON_SIZE), MAX_INTRON_SIZE)

        # adding element depending on the current last element
        to_add_ = int(x[-1] + intron_length * INTRON_BASE_SIZE)
        x.append(to_add_)
        # doing this twice
        to_add_ = int(
            x[-1]
            + (trans_data["exons"][i][1] - trans_data["exons"][i][0]) * EXON_BASE_SIZE
        )
        x.append(to_add_)

        human_exons.append(Intron((x[-3], y), x[-2] - x[-3], seperated))  # ith Intron
        is_ancestral = str(i + 1) in ancestral_exons
        human_exons.append(
            Exon(
                (x[-2], y),
                x[-1] - x[-2],
                opacity=OPACITY,
                ancestral=is_ancestral,
                annotation=None,
                draw_mode=drawing_mode,
            )
        )  # (i+1)th Exon
    human_exons.append(
        Exon((x[-1], y), UTR_WIDTH, is_exon=False, draw_mode=drawing_mode)
    )  # UTR
    return human_exons


def generate_filebuffer(mut_lines, human_exons_dct, x, y, drawing_mode=DEFAULT_MODE, verbose=False):
    """Generate the body of our SVG file."""
    print("Generating SVG line") if verbose else None
    curr_projection = None  # name of the current handled species
    curr_sp = None
    speciesdata = {}  # mutation data: {mutationnumber:(exon, start, stop)}
    comp_substrings = []
    filebuffer_substrings = []
    vert_offset = 0
    uniq_projections = set()

    for defect_line, sp_num_tup in mut_lines:
        defect_line_dat = defect_line.rstrip().split("\t")
        trans = defect_line_dat[0]
        chain = defect_line_dat[1]
        # sp = defect_line_dat[-1]
        sp_name = sp_num_tup[0]
        print(f"{defect_line_dat} {sp_num_tup} curr sp: {curr_sp}") if verbose else None
        human_exons = human_exons_dct[trans]
        projection = f"{trans}.{chain}"
        uniq_projections.add(projection)
        is_defect = len(defect_line_dat) == MUT_LINE_FIELDS_SP

        # check if we are at the same projection or now
        if projection != curr_projection or sp_num_tup != curr_sp:
            # initiate new plot!
            if PRINT_EXON_SCHEME or curr_projection is not None:
                height = 0
                for e in exons:  # draw genome of last projection
                    s, h = e.toString()
                    filebuffer_substrings.append(s)
                    height = max(height, h)
                for label in labels:
                    # draw defects and labels of last species
                    filebuffer_substrings.append(label.toString()[0])
                vert_offset += MAX_VERTICAL_SPACE
            # initiate new drawing!
            compensations = "".join(comp_substrings)
            filebuffer_substrings.append(compensations)
            exons = [copy(h) for h in human_exons]
            for e in exons:
                e.shift_vertical(vert_offset)

            label = f"{sp_name} {trans}.{chain}"
            labels = [
                (
                    Textstack(
                        (x[0], y + HALF_EXON_HEIGHT + MO_FONTSIZE),
                        HORIZONTAL,
                        label,
                        fontsize=MO_FONTSIZE,
                    )
                )
            ]
            for label in labels:
                label.shift_vertical(vert_offset)
            curr_projection = projection
            curr_sp = sp_num_tup
            speciesdata = {}
            comp_substrings = []

        # okay, now deal with what we have
        # 0 - trans 1 - chain 2 - exon_num
        # 3 - position  4 - mut class 5 - mut itself
        # 6 - trash 7 - mut ID
        if is_defect is False:
            continue
        defect_class = defect_line_dat[4]
        is_masked = True if defect_line_dat[6] == "masked" else False

        if defect_class != COMP:
            # define where we put the mutation
            exonnumber = int(defect_line_dat[2])
            exon = exons[2 * exonnumber - 1]
            pos = int(defect_line_dat[3])
            mut = defect_line_dat[5]
            mut_id = defect_line_dat[7]
        else:  # compensation track is a bit special
            if drawing_mode in DRAW_MODES_NO_COMPENSATIONS_SHOWN:
                # show it only if drawing mode allows so
                # otherwise simply skip it
                continue
            comp_ids_num = defect_line_dat[5].split("_")[1].split(",")
            fs_1_id = f"FS_{comp_ids_num[0]}"
            fs_2_id = f"FS_{comp_ids_num[1]}"
            try:
                data1, data2 = speciesdata[fs_1_id], speciesdata[fs_2_id]
                start = (
                    exons[2 * int(data1[0]) - 1].pos[0] + data1[2] * EXON_BASE_SIZE,
                    exons[2 * int(data1[0]) - 1].pos[1],
                )
                end = (
                    exons[2 * int(data2[0]) - 1].pos[0] + data2[1] * EXON_BASE_SIZE,
                    exons[2 * int(data2[0]) - 1].pos[1],
                )
                comp_substrings.append(draw_comp_indel(start, end))
            except KeyError:
                # TODO: check why FS2 could be deleted from list
                # if it belongs to a deleted exon: why compensation is
                # still there?
                pass
            continue

        # deal with splice site mutations
        if defect_class in SSM_TYPES:
            mouseover = None
            to_what = mut.split("->")[1]

            # color depends on the
            if drawing_mode == PUB_MODE_I:
                color_ = BLACK
            else:  # mode: default or non-affecting this mutation type
                color_ = MASKED_MUT_COLOR if is_masked else INACT_MUT_COLOR

            # if pos == 1:  # donor
            if defect_class == SSM_D:
                labels.append(
                    Textstack(
                        (exon.pos[0] + exon.width, exon.pos[1]),
                        HORIZONTAL,
                        to_what,
                        mouseover,
                        SS_LABEL_FONTSIZE,
                        f"fill:{color_};",
                    )
                )
            else:  # acceptor
                pos_tup_ = (
                    exon.pos[0] - len(to_what) * 0.64 * SS_LABEL_FONTSIZE,
                    exon.pos[1] + SS_LABEL_FONTSIZE,
                )
                labels.append(
                    Textstack(
                        pos_tup_,
                        HORIZONTAL,
                        to_what,
                        mouseover,
                        SS_LABEL_FONTSIZE,
                        f"fill:{color_};",
                    )
                )
        # deal with insertions
        elif defect_class in INS:
            start = pos
            length = int(mut.replace("+", "").replace("-", ""))
            mouseover = None
            exon.add_insertion(start, length, mouseover, is_masked)
            speciesdata[mut_id] = (exonnumber, start, start + length - 1)
        # the same for deletions
        elif defect_class in DELS:
            start = pos
            length = int(mut.replace("+", "").replace("-", ""))
            mouseover = None
            exon.add_deletion(start, length, mouseover, is_masked)
            speciesdata[mut_id] = (exonnumber, start, start + length - 1)
        # add deleted exon
        elif defect_class == EX_DEL:

            if drawing_mode == PUB_MODE_I:
                # in this mode: all deleted exons in red, no blue ones anymore
                exon.color = INACT_MUT_COLOR
            else:  # default mode or not affecting deleted exons color
                exon.color = FP_EXON_DEL_COLOR if is_masked else INACT_MUT_COLOR

            speciesdata[mut_id] = (exonnumber, 0, 0)
        # and also with missing exons
        elif defect_class == EX_MIS:
            exon.color = MISS_SEQ_COLOR
            speciesdata[mut_id] = (exonnumber, 0, 0)
        # draw stop
        elif defect_class == STOP:
            if "->" in mut:
                codon = mut.split("->")[1]
            else:
                # split stop codon
                codon = mut
            mouseover = None
            exon.add_stop_codon(pos, codon, mouseover, is_masked)
        else:
            pass

    if curr_projection is not None:
        vert_offset += MAX_VERTICAL_SPACE
        vbox_h = 0
        for e in exons:
            s, h = e.toString()
            filebuffer_substrings.append(s)
            vbox_h = max(vbox_h, h)
        compensations = "".join(comp_substrings)
        filebuffer_substrings.append(compensations)
        for label in labels:
            filebuffer_substrings.append(label.toString()[0])

    filebuffer = "".join(filebuffer_substrings)
    return filebuffer, len(uniq_projections), vbox_h


def compute_pic_and_vbox_h(up_num, vbox_h):
    """Compute vbox and picture h."""
    # B: don't really understand what it stands for.
    picture_height = MAX_VERTICAL_SPACE if PRINT_EXON_SCHEME else 0
    picture_height += up_num * MAX_VERTICAL_SPACE
    if up_num == 1:
        picture_height = max(picture_height, MAX_VERTICAL_SPACE / 2 + vbox_h)
    else:
        vbox_h = MAX_VERTICAL_SPACE / 2
    return picture_height, vbox_h


def get_transcripts(gene_id, isoforms_file):
    """Get transcript for gene."""
    transcripts = []
    f = open(isoforms_file, "r")
    for line in f:
        line_data = line.rstrip().split("\t")
        gene = line_data[0]
        trans = line_data[1]
        if gene == gene_id:
            transcripts.append(trans)
    f.close()
    return transcripts


def get_sp_to_mbd_data(mut_files_list):
    """Read species: mut file data."""
    sp_to_file = {}
    f = open(mut_files_list, "r")
    for num, line in enumerate(f, 1):
        line_data = line.rstrip().split()
        if len(line_data) < 2:
            err_msg = (
                f"Error! {mut_files_list} file is incorrect. The following format "
                f"applies:\n[Species name] <space> [path to corresponding inact mut "
                f"file]. Species name can also be space-separated."
            )
            print(err_msg)
            sys.exit(1)
        sp = " ".join(line_data[:-1])
        rel_path = line_data[-1]
        dir_name = os.path.dirname(mut_files_list)
        abs_path = os.path.join(dir_name, rel_path)
        sp_to_file[(sp, num)] = abs_path
    f.close()
    return sp_to_file


def make_plot(
    bed_file,
    toga_mut_file,
    transcripts,
    chain,
    isoforms_file,
    multi_species,
    alt_template,
    inact_lines=None,
    drawing_mode=DEFAULT_MODE,
    verbose=False
):
    """Enrty point."""
    if isoforms_file is None:
        transcripts = [t for t in transcripts.split(",") if t != ""]
    else:
        transcripts = get_transcripts(transcripts, isoforms_file)

    if multi_species is True:
        sp_to_mdb = get_sp_to_mbd_data(toga_mut_file)
    else:
        sp_to_mdb = {("SP", 0): toga_mut_file}

    if verbose:
        print(f"SP to Mut data: {sp_to_mdb}\n\n")

    all_mut_lines = []
    trans_to_h_exons = {}
    picture_width = 0

    for transcript in transcripts:
        if verbose:
            print(f"Parsing data for transcript {transcript}")

        trans_bed = get_transcript(bed_file, transcript)
        if trans_bed is None:
            # not found
            sys.stderr.write(
                f"Warning! Cannot find {transcript} in the bed {bed_file}\n"
            )
            continue

        trans_data = parse_trans_data(trans_bed)
        trans_data_direction(trans_data)

        for sp, mdb in sp_to_mdb.items():
            print(f"# Extracting data for {sp} {mdb}") if verbose else None
            if inact_lines is None:
                # ~100% cases
                mut_lines_raw = get_mut_list(mdb, transcript)
            else:
                # called by external program: inact data already extracted
                mut_lines_raw = inact_lines
            print(f"Extracted {len(mut_lines_raw)} raw mut lines") if verbose else None

            mut_lines = prepare_mut_data(
                mut_lines_raw, sp, trans_data["rel_starts"], chain_lim=chain
            )
            print(f"After filter: {len(mut_lines_raw)} lines") if verbose else None

            all_mut_lines.extend(mut_lines)

        t_picture_width = comp_width(trans_data)
        picture_width = (
            t_picture_width if t_picture_width > picture_width else picture_width
        )

        print(f"Picture width: {picture_width}") if verbose else None

        # something we need
        x = [
            UTR_WIDTH,
            UTR_WIDTH
            + (trans_data["exons"][0][1] - trans_data["exons"][0][0]) * EXON_BASE_SIZE,
        ]
        y = MAX_VERTICAL_SPACE / 2

        ancestral_exons = [str(x) for x in range(1, trans_data["exon_num"] + 1)]
        human_exons = init_human_exons(
            x, y, ancestral_exons, trans_data, drawing_mode=drawing_mode
        )
        trans_to_h_exons[transcript] = human_exons

    if len(trans_to_h_exons.keys()) == 0:
        sys.exit("Not found any data for your transcripts")
    if len(all_mut_lines) == 0 and inact_lines is None:
        sys.exit(f"{str(transcripts)}: Please check your filters -> nothing to plot")
    elif len(all_mut_lines) == 0:
        return SVG_PLACEHOLDER
    filebuffer, up_num, vbox_h = generate_filebuffer(
        all_mut_lines, trans_to_h_exons, x, y, drawing_mode=drawing_mode, verbose=verbose
    )

    # compute height
    picture_height, vbox_h = compute_pic_and_vbox_h(up_num, vbox_h)

    # save figure
    mouse_counter_val = MouseOver.id_counter + 1
    svg_line = write_pic(
        filebuffer,
        picture_width,
        picture_height,
        vbox_h,
        mouse_counter_val,
        alt_template,
    )
    return svg_line


def __infer_draw_mode(args):
    """Parse args to infer drawing mode.

    If nothing specified: use DEFAULT.
    """
    # TODO: if several modes available, exit if > 1 modes assigned
    if args.publication_mode_heni:
        return PUB_MODE_I
    return DEFAULT_MODE


if __name__ == "__main__":
    args = parse_args()
    _d_mode = __infer_draw_mode(args)
    svg_line = make_plot(
        args.bed_file,
        args.toga_mut_file,
        args.transcripts,
        args.chain,
        args.isoforms_file,
        args.multi_species,
        args.alt_template,
        drawing_mode=_d_mode,
        verbose=args.verbose
    )
    f = open(args.output, "w") if args.output != "stdout" else sys.stdout
    f.write(svg_line)
    f.close() if args.output != "stdout" else None
