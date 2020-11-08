#!/usr/bin/env python3
"""Classify chain-gene pairs as orthologous, paralogous etc.

Each chain-gene pair has a set of features.
We have a XGBoost pre-trained model that can classify them.
"""
import argparse
import sys
from collections import defaultdict
import functools
from numpy import log10
import pandas as pd
import joblib
import xgboost as xgb

try:  # for robustness
    from modules.common import eprint
    from modules.common import die
except ImportError:
    from common import eprint
    from common import die


__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]

# paths to single (SE) and multi-exon (ME) models
SE_MODEL = "models/se_model.dat"
ME_MODEL = "models/me_model.dat"

ORTH = "ORTH"
PARA = "PARA"
TRANS = "TRANS"
P_PGENES = "P_PGENES"

# we actually extract a redundant amount of features
# lists of features required by single and multi exon models
SE_MODEL_FEATURES = ['gl_exo', 'flank_cov', 'exon_perc', 'synt_log']
ME_MODEL_FEATURES = ['gl_exo', 'loc_exo', 'flank_cov', 'synt_log', 'intr_perc']
print = functools.partial(print, flush=True)


def verbose(msg):
    """Eprint for verbose messages."""
    print(msg + "\n")


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
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def classify_chains(table, output, se_model_path, me_model_path,
                    raw_out=None, rejected=None,
                    annot_threshold=0.5):
    """Core chain classifier function."""
    verbose("Loading dataset...")
    # read dataframe
    df = pd.read_csv(table, header=0, sep="\t")
    init_genes_set = set(df["gene"])

    # get indexes of processed pseudogene chains
    # if a chain has synteny = 1, introns are deleted and exon_num > 1
    # -> then this is a proc pseudogene

    # move trans chains to a different dataframe
    # trans chain -> a syntenic chain that passes throw the gene body
    #                but has no aligning bases in the CDS
    trans_lines = df[(df["exon_cover"] == 0) & (df["synt"] > 1)]
    # remove from dataframe: (this includes trans chains)
    # 1) chains that don't cover CDS
    # 2) remove chains that have synteny == 0
    df = df[(df["exon_cover"] > 0) & (df["synt"] > 0)]
    df_final = df.copy()  # filtered dataframe: what we will classify
    # compute some necessary features
    df_final["exon_perc"] = df_final["exon_cover"] / df_final["ex_fract"]
    df_final["chain_len_log"] = log10(df_final["chain_len"])
    df_final["synt_log"] = log10(df_final["synt"])
    df_final["intr_perc"] = df_final["intr_cover"] / df_final["intr_fract"]
    df_final = df_final.fillna(0.0)  # fill NA values with 0.0

    if len(df_final) > 0:
        # add "is single exon" column -> to separate it for different models
        # need to do this only if the df is not empty
        df_final.loc[df_final["ex_num"] == 1, "single_exon"] = 1
        df_final.loc[df_final["ex_num"] != 1, "single_exon"] = 0
    else:  # this df is empty anyway, so any value fits
        df_final["single_exon"] = 0

    # split df into two: for single and multi exon models
    df_se = df_final[df_final["single_exon"] == 1]
    df_me = df_final[df_final["single_exon"] == 0]
    # create X dataframes: skip unnecessary columns
    X_se = df_se.copy()
    X_me = df_me.copy()
    X_se = X_se[SE_MODEL_FEATURES]
    X_me = X_me[ME_MODEL_FEATURES]

    # load models
    verbose("Load and apply model")
    try:  # there are 2 potential problems:
        # no files at all and
        # cannot load files
        se_model = joblib.load(se_model_path)
        me_model = joblib.load(me_model_path)
    except FileNotFoundError:
        err_msg = f"Cannot find models {se_model_path} and {me_model_path}\n" \
                  f"Please call train_model.py to create them.\nAbort."
        raise FileNotFoundError(err_msg)
    except (xgb.core.XGBoostError, AttributeError):
        xgboost_version = xgb.__version__
        err_msg = f"Cannot load models {se_model_path} and {me_model_path} " \
                  f"Probably, models were trained with a different version of " \
                  f"XGBoost. You used XBGoost version: {xgboost_version}; " \
                  f"Please make sure you called train_model.py with the same version."
        raise ValueError(err_msg)
    # and apply them
    se_pred = se_model.predict_proba(X_se)[:, 1]
    me_pred = me_model.predict_proba(X_me)[:, 1]

    # add predictions to the dataframe
    # prediction is basically a single-column
    df_se_result = df_se.copy()
    df_me_result = df_me.copy()
    trans_result = trans_lines.copy()
    df_se_result["pred"] = se_pred
    df_me_result["pred"] = me_pred
    # model prediction is a float from 0 to 1, -1 -> for trans chains
    trans_result["pred"] = -1
    # identify processed pseudogenes, they satisfy the following criteria:
    # 1) multi-exon (single-exon ones are out of score of the method)
    # 2) synteny == 1
    # 3) cds_to_qlen > 0.95
    # set them score -2
    df_me_result.loc[(df_me_result["synt"] == 1)
                     & (df_me_result["exon_qlen"] > 0.95)
                     & (df_me_result["pred"] < 0.5), "pred"] = -2

    # we need gene -> chain -> prediction from each row
    df_se_result = df_se_result.loc[:, ["gene", "chain", "pred"]]
    df_me_result = df_me_result.loc[:, ["gene", "chain", "pred"]]
    trans_result = trans_result.loc[:, ["gene", "chain", "pred"]]
    # concatenate the results
    overall_result = pd.concat([df_se_result, df_me_result, trans_result])

    if raw_out:  # save raw scores if required
        overall_result.to_csv(raw_out, sep="\t", index=False)

    # create a different TSV
    # transcript: lists of different chain classes
    # such as transcript A: [orthologous chains] [paralogous chains] etc
    # 0 -> placeholder, means "the class if empty"
    gene_class_chains = {}
    genes = set(overall_result["gene"])
    for gene in genes:
        gene_class_chains[gene] = {ORTH: [], PARA: [], TRANS: [], P_PGENES: []}

    for num, data in enumerate(overall_result.itertuples()):
        gene = data.gene
        chain = data.chain
        pred = data.pred

        if pred == -1:  # spanning (trans) chain
            gene_class_chains[gene][TRANS].append(chain)
        elif pred == -2:  # processed pseudogene
            gene_class_chains[gene][P_PGENES].append(chain)
        elif pred < annot_threshold:
            gene_class_chains[gene][PARA].append(chain)
        else:  # > annot_threshold
            gene_class_chains[gene][ORTH].append(chain)
        # verbose
        if num % 10000 == 0:
            print(num)

    # save orthologs output
    f = open(output, "w") if output != "stdout" else sys.stdout
    f.write(f"GENE\t{ORTH}\t{PARA}\t{TRANS}\t{P_PGENES}\n")
    i = 0
    for k, v in gene_class_chains.items():
        i += 1
        orth = v[ORTH]
        para = v[PARA]
        trans = v[TRANS]
        pp = v[P_PGENES]
        orth_f = ",".join(str(x) for x in orth) if orth else "0"
        para_f = ",".join(str(x) for x in para) if para else "0"
        trans_f = ",".join(str(x) for x in trans) if trans else "0"
        p_pgenes_f = ",".join(str(x) for x in pp) if pp else "0"
        f.write("\t".join([k, orth_f, para_f, trans_f, p_pgenes_f]) + "\n")
    f.close() if output != "stdout" else None
    # for some genes there are no classifiable chains, save them
    genes_missing = list(init_genes_set.difference(genes))

    if rejected:  # we requested to save transcripts without classifiable chains
        f = open(rejected, "w")
        for gene in genes_missing:
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
