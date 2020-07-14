#!/usr/bin/env python3
"""ML version of classify chains script."""
import argparse
import sys
from collections import defaultdict
from numpy import log10
import pandas as pd
import joblib


__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]

SE_MODEL = "models/se_model.dat"
ME_MODEL = "models/me_model.dat"


se_drop_in_X = ["gene",
                "gene_overs",
                "chain",
                "synt",
                "gl_score", 
                "chain_len",
                "exon_cover",
                "intr_cover",
                "gene_len",
                "ex_num",
                "ex_fract",
                "intr_fract",
                "chain_len_log",
                "single_exon",
                "intr_perc",
                "loc_exo",
                "cds_qlen"]
me_drop_in_X = ["gene",
                "gene_overs",
                "chain",
                "synt",
                "gl_score", 
                "chain_len",
                "exon_cover",
                "intr_cover",
                "gene_len",
                "ex_num",
                "ex_fract",
                "intr_fract",
                "chain_len_log",
                "exon_perc",
                "single_exon",
                "cds_qlen"]


def eprint(*lines):
    """Like print but for stderr."""
    for line in lines:
        sys.stderr.write(line + "\n")


def die(msg):
    """Write msg to stderr and abort program."""
    eprint(msg)
    sys.exit(1)


def verbose(msg):
    """Eprint for verbose messages."""
    eprint(msg + "\n")


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("table", type=str,
                     help="Table containing merged results.")
    app.add_argument("output", type=str,
                     help="Write the table here.")
    app.add_argument("se_model", help="A trained XGBoost model for single-exon genes.")
    app.add_argument("me_model", help="A trained XGBoost model for multi-exon genes.")
    app.add_argument("--raw_model_out", "--ro", default=None,
                     help="Save gene: chain xgboost output")
    # app.add_argument("--verbose", "-v", action="store_true", dest="verbose",
    #                  help="Verbose messages.")
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def classify_chains(table, output, se_model_path, me_model_path,
                    raw_out=None, rejected=None,
                    annot_threshold=0.5):
    """Entry point."""
    # global VERBOSE  # set verbosity level
    # VERBOSE = True if args.verbose else False
    verbose("Loading dataset...")
    df = pd.read_csv(table, header=0, sep="\t")
    init_genes_set = set(df["gene"])

    # remove processed pseudogenes
    # p_pseudogenes = df[(df["synt"] == 1) &
    #                    (df["intr_cover"] == 0) &
    #                    (df["flank_cov"] == 0) &
    #                    (df["ex_num"] > 1)].index
    # remove processed pseudogenes
    p_pseudogenes = df[(df["synt"] == 1) &
                       (df["cds_qlen"] > 0.95) &
                       (df["ex_num"] > 1)].index
    verbose(f"Found {len(p_pseudogenes)} processed pseudogenes")
    gene_to_pp_genes = defaultdict(list)
    for ind in p_pseudogenes:
        p_line = df.iloc[ind]
        gene = p_line["gene"]
        chain = p_line["chain"]
        gene_to_pp_genes[gene].append(chain)
    df.drop(p_pseudogenes, inplace=True)
    # remove trans chains
    trans_lines = df[df["exon_cover"] == 0]
    # df_no_trans = df[df["exon_cover"] > 0]
    trans_lines = trans_lines[trans_lines["synt"] > 1]
    gene_trans = defaultdict(list)
    df = df[df["exon_cover"] > 0]  # no need if no CDS
    df = df[df["synt"] > 0]  # also skip this

    verbose("Extracting TRANS chains...")
    for row in trans_lines.itertuples():
        gene_trans[row.gene].append(row.chain)

    df_final = df.copy()
    df_final["exon_perc"] = df_final["exon_cover"] / df_final["ex_fract"]
    df_final["chain_len_log"] = log10(df_final["chain_len"])
    df_final["synt_log"] = log10(df_final["synt"])
    df_final["intr_perc"] = df_final["intr_cover"] / df_final["intr_fract"]
    df_final.loc[df_final["ex_num"] == 1, "single_exon"] = 1
    df_final.loc[df_final["ex_num"] != 1, "single_exon"] = 0
    df_final = df_final.fillna(0.0)  # fill NA values

    # need > 0, because if no df_no_trans, then will be Value Error
    if len(df_final) > 0:
        df_final.loc[df_final["ex_num"] == 1, "single_exon"] = 1
        df_final.loc[df_final["ex_num"] != 1, "single_exon"] = 0
    else:
        # this df is empty anyway, so any value fits
        df_final["single_exon"] = 0

    df_se = df_final[df_final["single_exon"] == 1]
    df_me = df_final[df_final["single_exon"] == 0]
    X_se = df_se.copy()
    X_me = df_me.copy()
    X_se = df_se.drop(se_drop_in_X, axis=1)
    X_me = df_me.drop(me_drop_in_X, axis=1)

    # load model and apply
    verbose("Load and apply model")
    se_model = joblib.load(se_model_path)
    me_model = joblib.load(me_model_path)

    se_pred = se_model.predict_proba(X_se)[:,1]
    me_pred = me_model.predict_proba(X_me)[:,1]

    df_se_result = df_se.copy()
    df_me_result = df_me.copy()
    trans_result = trans_lines.copy()
    df_se_result["pred"] = se_pred
    df_me_result["pred"] = me_pred
    trans_result["pred"] = -1

    df_se_result = df_se_result.loc[:, ["gene", "chain", "pred"]]
    df_me_result = df_me_result.loc[:, ["gene", "chain", "pred"]]
    trans_result = trans_result.loc[:, ["gene", "chain", "pred"]]
    overall_result = pd.concat([df_se_result, df_me_result, trans_result])

    if raw_out:  # save raw scores
        overall_result.to_csv(raw_out, sep="\t", index=False)

    gene_class_chains = {}
    genes = set(overall_result["gene"])
    for gene in genes:
        gene_class_chains[gene] = {"ORTH": [], "PARA": [], "TRANS": [], "P_PGENES": []}

    for num, data in enumerate(overall_result.itertuples()):
        gene = data.gene
        chain = data.chain
        pred = data.pred

        if pred == -1:  # spanning chain
            gene_class_chains[gene]["TRANS"].append(chain)
        elif pred < annot_threshold:
            gene_class_chains[gene]["PARA"].append(chain)
        else:  # > annot_threshold
            gene_class_chains[gene]["ORTH"].append(chain)
        # verbose
        if gene_to_pp_genes[gene]:
            gene_class_chains[gene]["P_PGENES"] = gene_to_pp_genes[gene]
        if num % 10000 == 0:
            print(num)

    # save orthologs output
    f = open(output, "w") if output != "stdout" else sys.stdout
    f.write("GENE\tORTH\tPARA\tTRANS\tP_PGENES\n")
    i = 0
    for k, v in gene_class_chains.items():
        i += 1
        orth = v["ORTH"]
        para = v["PARA"]
        trans = v["TRANS"]
        pp = v["P_PGENES"]
        orth_f = ",".join(str(x) for x in orth) if orth else "0"
        para_f = ",".join(str(x) for x in para) if para else "0"
        trans_f = ",".join(str(x) for x in trans) if trans else "0"
        p_pgenes_f = ",".join(str(x) for x in pp) if pp else "0"
        f.write("\t".join([k, orth_f, para_f, trans_f, p_pgenes_f]) + "\n")
    f.close() if output != "stdout" else None
    genes_missed = list(init_genes_set.difference(genes))

    # save rejected genes
    if rejected:
        f = open(rejected, "w")
        for gene in genes_missed:
            f.write(f"{gene}\tNo classifiable chains\n")
        f.close()



def main():
    args = parse_args()
    classify_chains(args.table,
                    args.output,
                    args.se_model,
                    args.me_model,
                    raw_out=args.raw_model_out)

if __name__ == "__main__":
    main()
