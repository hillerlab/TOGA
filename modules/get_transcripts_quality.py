#!/usr/bin/env python3
"""Using metadata given, create transcript confidence level."""
import argparse
import sys
from collections import defaultdict

META_DATA_FIELDS_NUM = 12


def parse_args():
    """Parse and check args."""
    app = argparse.ArgumentParser()
    app.add_argument("meta_data", help="Exons meta data")
    app.add_argument("orthology_scores", help="Orthology scores file")
    app.add_argument("--hq_threshold", "--hqt", type=float, default=0.95,
                     help="Orthology score for high-confidence transcripts")
    app.add_argument("output", help="Save output to")
    if len(sys.argv) < 3:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    if args.hq_threshold > 1.0:
        sys.exit("Error! HQ threshold must be <= 1.0!")
    return args


def get_exon_marks(meta_data_file):
    """Collect exon marks for each projection."""
    transcript_exon_marks = defaultdict(list)
    f = open(meta_data_file, "r")
    for num, line in enumerate(f, 1):
        line_data = line.rstrip().split("\t")
        if line_data[0] == "":
            # a gap
            continue
        elif len(line_data) != META_DATA_FIELDS_NUM:
            # if not true -> file is corrupted
            sys.stderr.write(f"Line caused errors:\n")
            sys.stderr.write(f"{line}\n")
            sys.stderr.write(f"Line num: {num}\n")
            raise ValueError("Input file is corrupted")
        elif line_data[0] == "gene":
            # header
            continue
        transcript = line_data[0]
        chain = line_data[2]
        transcript_key = f"{transcript}.{chain}"
        exon_mark = line_data[11]
        transcript_exon_marks[transcript_key].append(exon_mark)
    f.close()
    return transcript_exon_marks



def get_transcript_score(ort_scores_file):
    """Read orthology scores."""
    result = {}
    f = open(ort_scores_file, "r")
    _ = f.__next__()  # skip header
    for line in f:
        line_data = line.rstrip().split("\t")
        transcript = line_data[0]
        chain = line_data[1]
        score = float(line_data[2])
        projection = f"{transcript}.{chain}"
        result[projection] = score
    f.close()
    return result


def classify_transcripts(meta_data_file, ort_score_file, hq_threshold, output):
    """Classify transcripts."""
    # first, get exon marks
    transcript_exon_marks = get_exon_marks(meta_data_file)
    transcripts_scores = get_transcript_score(ort_score_file)
    # ther get transcript classes
    result = {}
    for trans, marks in transcript_exon_marks.items():
        trans_score = transcripts_scores[trans]
        if "LQ" in marks:
            t_mark = "low_confidence"
        elif "NA" in marks:
            t_mark = "partial"
        elif all(x == "HQ" for x in marks) and trans_score >= hq_threshold:
            t_mark = "high_confidence"
        else:
            t_mark = "average_confidence"
        result[trans] = t_mark  
    save(result, output)

def save(transcript_class, output_file):
    """Just save the table."""
    f = open(output_file, "w")
    f.write("Projection_ID\tconfidence_level\n")
    for k, v in transcript_class.items():
        f.write(f"{k}\t{v}\n")
    f.close()


def main():
    """Entry point."""
    args = parse_args()
    classify_transcripts(args.meta_data,
                         args.orthology_scores,
                         args.hq_threshold,
                         args.output)

if __name__ == "__main__":
    main()
