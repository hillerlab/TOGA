#!/usr/bin/env python3
"""Classify chain-gene pairs as orthologous, paralogous etc.

Each chain-gene pair has a set of features.
We have a XGBoost pre-trained model that can classify them.
"""
import argparse
import sys
import numpy as np
import pandas as pd
import joblib
import xgboost as xgb
from version import __version__

try:  # TODO: check whether it's needed
    from modules.common import setup_logger
    from modules.common import to_log
    from modules.common import die
except ImportError:
    from modules.common import setup_logger
    from modules.common import to_log
    from common import die

__author__ = "Bogdan M. Kirilenko"

# paths to single (SE) and multi-exon (ME) models
SE_MODEL = "models/se_model.dat"
ME_MODEL = "models/me_model.dat"
LD_MODEL = "models/long_dist_model.dat"

ORTH = "ORTH"
PARA = "PARA"
SPAN = "SPAN"
P_PGENES = "P_PGENES"

SPANNING_SCORE = -1.0
PPGENE_SCORE = -2.0

# we actually extract a redundant amount of features
# lists of features required by single and multi exon models
SE_MODEL_FEATURES = ["gl_exo", "flank_cov", "exon_perc", "synt_log"]
ME_MODEL_FEATURES = ["gl_exo", "loc_exo", "flank_cov", "synt_log", "intr_perc"]
LD_MODEL_FEATURES = ["gl_exo", "flank_cov", "exon_perc", "synt_log", "loc_exo",
                     "intr_perc", "score", "single_exon"]


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser()
    app.add_argument("table", type=str, help="Table containing merged results.")
    app.add_argument("output", type=str, help="Write the table here.")
    app.add_argument("se_model", help="A trained XGBoost model for single-exon transcripts.")
    app.add_argument("me_model", help="A trained XGBoost model for multi-exon transcripts.")
    app.add_argument(
        "--raw_model_out", "--ro", default=None, help="Save gene: chain xgboost output"
    )
    app.add_argument("--ld_model", action="store_true", dest="ld_model", help="Apply LD model")
    app.add_argument("--log_file", default=None, help="Path to the log file")
    # print help if there are no args
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def classify_chains(
    table,
    output,
    se_model_path,
    me_model_path,
    raw_out=None,
    rejected=None,
    annot_threshold=0.5,
    ld_model=None
):
    """Core chain classifier function."""
    # read dataframe
    df = pd.read_csv(table, header=0, sep="\t")
    init_transcripts_set = set(df["gene"])
    to_log(f"classify_chains: loaded dataframe of size {len(df)}")
    to_log(f"classify_chains: total number of transcripts: {len(df)}")

    # get indexes of processed pseudogene chains
    # if a chain has synteny = 1, introns are deleted and exon_num > 1
    # -> then this is a proc pseudogene
    # move spanning chains to a different dataframe
    # TODO: rename trans -> spanning
    # trans chain -> a syntenic chain that passes throw the gene body
    #                but has no aligning bases in the CDS
    spanning_ts_lines = df[(df["exon_cover"] == 0) & (df["synt"] > 1)]
    to_log(f"classify_chains: {len(spanning_ts_lines)} rows with spanning chains")
    # remove from dataframe: (this includes trans chains)
    # 1) chains that don't cover CDS
    # 2) remove chains that have synteny == 0
    df = df[(df["exon_cover"] > 0) & (df["synt"] > 0)]
    df_final = df.copy()  # filtered dataframe: what we will classify
    to_log(f"classify_chains: filtered dataset contains {len(df_final)} records")
    # compute some necessary features
    to_log("classify_chains: omputing additional features...")
    df_final["exon_perc"] = df_final["exon_cover"] / df_final["ex_fract"]
    df_final["chain_len_log"] = np.log10(df_final["chain_len"])
    df_final["synt_log"] = np.log10(df_final["synt"])
    df_final["intr_perc"] = df_final["intr_cover"] / df_final["intr_fract"]
    df_final = df_final.fillna(0.0)  # fill NA values with 0.0

    if len(df_final) > 0:
        # add "is single exon" column -> to separate it for different models
        # need to do this only if the df is not empty
        df_final.loc[df_final["ex_num"] == 1, "single_exon"] = 1
        df_final.loc[df_final["ex_num"] != 1, "single_exon"] = 0
    else:  # this df is empty anyway, so any value fits
        to_log("classify_chains: WARNING! The final df for classification is empty")
        df_final["single_exon"] = 0

    # split df into two: for single and multi exon models
    df_se = df_final[df_final["single_exon"] == 1]
    df_me = df_final[df_final["single_exon"] == 0]
    to_log(f"classify_chains: df for single-exon model contains {len(df_se)} records")
    to_log(f"classify_chains: df for multi-exon model contains {len(df_me)} records")

    # create X dataframes: skip unnecessary columns
    X_se = df_se.copy()
    X_me = df_me.copy()
    X_se = X_se[SE_MODEL_FEATURES]
    X_me = X_me[ME_MODEL_FEATURES]

    # load models
    try:  # there are 2 potential problems:
        # no files at all and
        # cannot load files
        to_log(f"classify_chains: loading models at {se_model_path} (SE) and {me_model_path} (ME)")
        se_model = joblib.load(se_model_path)
        me_model = joblib.load(me_model_path)
    except FileNotFoundError:
        err_msg = (
            f"Cannot find models {se_model_path} and {me_model_path}\n"
            f"Please call train_model.py to create them.\nAbort."
        )
        to_log(f"classify_chains: ERROR! {err_msg}")
        raise FileNotFoundError(err_msg)
    except (xgb.core.XGBoostError, AttributeError):
        xgboost_version = xgb.__version__
        err_msg = (
            f"Cannot load models {se_model_path} and {me_model_path} "
            f"Probably, models were trained with a different version of "
            f"XGBoost. You used XBGoost version: {xgboost_version}; "
            f"Please make sure you called train_model.py with the same version."
        )
        to_log(f"classify_chains: ERROR! {err_msg}")
        raise ValueError(err_msg)

    # and apply them
    to_log(f"classify_chains: applying models to SE and ME datasets...")
    me_pred = me_model.predict_proba(X_me)[:, 1] if len(X_me) > 0 else np.array([])
    se_pred = se_model.predict_proba(X_se)[:, 1] if len(X_se) > 0 else np.array([])

    # add predictions to the dataframe
    # prediction is basically a single-column
    df_se_result = df_se.copy()
    df_me_result = df_me.copy()
    spanning_chains_result = spanning_ts_lines.copy()

    df_se_result["pred"] = se_pred
    df_me_result["pred"] = me_pred

    if ld_model:
        to_log(f"classify_chains: applying model for higher molecular distances")
        to_log(f"classify_chains: WARNING! This is an experimental feature")
        to_log(f"classify_chains: it is not recommended to use for research purposes yet")
        # apply LD model in addition
        # score from previous prediction is a feature for the LD model
        ld_model = joblib.load(ld_model)
        df_se_result = df_se_result.rename({"pred": "score"}, axis="columns")
        df_me_result = df_me_result.rename({"pred": "score"}, axis="columns")
        X_LD_me = df_me_result[LD_MODEL_FEATURES].copy()
        X_LD_se = df_se_result[LD_MODEL_FEATURES].copy()
        me_ld_pred = ld_model.predict_proba(X_LD_me)[:, 1] if len(X_LD_me) > 0 else np.array([])
        se_ld_pred = ld_model.predict_proba(X_LD_se)[:, 1] if len(X_LD_se) > 0 else np.array([])
        df_se_result["pred"] = se_ld_pred
        df_me_result["pred"] = me_ld_pred

    # model prediction is a float from 0 to 1, -1 -> for trans chains
    to_log(f"classify_chains: applying -1.0 score to the spanning chains")
    spanning_chains_result["pred"] = SPANNING_SCORE
    # identify processed pseudogenes, they satisfy the following criteria:
    # 1) multi-exon (single-exon ones are out of score of the method)
    # 2) synteny == 1
    # 3) cds_to_qlen > 0.95
    # set them score -2
    to_log(f"classify_chains: applying -2.0 score to the processed pseudogene alignments")
    df_me_result.loc[
        (df_me_result["synt"] == 1)
        & (df_me_result["exon_qlen"] > 0.95)
        & (df_me_result["pred"] < annot_threshold)
        & (df_me_result["exon_perc"] > 0.65),
        "pred",
    ] = PPGENE_SCORE

    pp_gene_count = df_me_result[df_me_result["pred"] == PPGENE_SCORE].shape[0]
    to_log(f"classify_chains: number of processed pseudogene alignments: {pp_gene_count}")

    # we need gene -> chain -> prediction from each row
    to_log(f"classify_chains: arranging the final output")
    df_se_result = df_se_result.loc[:, ["gene", "chain", "pred"]]
    df_me_result = df_me_result.loc[:, ["gene", "chain", "pred"]]
    spanning_chains_result = spanning_chains_result.loc[:, ["gene", "chain", "pred"]]
    # concatenate the results
    overall_result = pd.concat([df_se_result, df_me_result, spanning_chains_result])
    # some stats
    paralogs_count = overall_result[overall_result["pred"] < annot_threshold].shape[0]
    orthologs_count = overall_result[overall_result["pred"] >= annot_threshold].shape[0]
    spanning_chains_count = overall_result[overall_result["pred"] == 1.0].shape[0]
    pp_genes_count = overall_result[overall_result["pred"] == 2.0].shape[0]
    stats_msg = (
        f"* orthologs: {orthologs_count}\n"
        f"* paralogs: {paralogs_count}\n"
        f"* spanning chains: {spanning_chains_count}\n"
        f"* processed pseudogenes: {pp_genes_count}"
    )
    to_log(f"classify_chains: classification result stats:\n{stats_msg}")
    to_log(f"classify_chains: using {annot_threshold} as a threshold to separate orthologs from paralogs")

    if raw_out:  # save raw scores if required
        overall_result.to_csv(raw_out, sep="\t", index=False)

    # create a different TSV
    # transcript: lists of different chain classes
    # such as transcript A: [orthologous chains] [paralogous chains] etc
    # 0 -> placeholder, means "the class if empty"
    gene_class_chains = {}
    transcripts = set(overall_result["gene"])
    for transcript in transcripts:
        gene_class_chains[transcript] = {ORTH: [], PARA: [], SPAN: [], P_PGENES: []}

    to_log(f"classify_chains: combining results for {len(transcripts)} individual transcripts")
    for data in overall_result.itertuples():
        gene = data.gene
        chain = data.chain
        pred = data.pred

        if pred == SPANNING_SCORE:  # spanning (trans) chain
            gene_class_chains[gene][SPAN].append(chain)
        elif pred == PPGENE_SCORE:  # processed pseudogene
            gene_class_chains[gene][P_PGENES].append(chain)
        elif pred < annot_threshold:
            gene_class_chains[gene][PARA].append(chain)
        else:  # > annot_threshold
            gene_class_chains[gene][ORTH].append(chain)

    # save orthologs output
    to_log(f"classify_chains: saving the classification to {output}")
    f = open(output, "w") if output != "stdout" else sys.stdout
    f.write(f"GENE\t{ORTH}\t{PARA}\t{SPAN}\t{P_PGENES}\n")
    i = 0
    for k, v in gene_class_chains.items():
        i += 1
        orth = v[ORTH]
        para = v[PARA]
        trans = v[SPAN]
        pp = v[P_PGENES]
        orth_f = ",".join(str(x) for x in orth) if orth else "0"
        para_f = ",".join(str(x) for x in para) if para else "0"
        trans_f = ",".join(str(x) for x in trans) if trans else "0"
        p_pgenes_f = ",".join(str(x) for x in pp) if pp else "0"
        f.write("\t".join([k, orth_f, para_f, trans_f, p_pgenes_f]) + "\n")
    f.close() if output != "stdout" else None
    # for some transcripts there are no classifiable chains, save them
    transcripts_missing = list(init_transcripts_set.difference(transcripts))
    to_log(f"classify_chains: found no classifiable chains for {len(transcripts_missing)} transcripts")

    if rejected:  # we requested to save transcripts without classifiable chains
        to_log(f"classify_chains: saving these transcripts to: {rejected}")
        f = open(rejected, "w")
        for transcript in transcripts_missing:
            f.write(f"{transcript}\tNo classifiable chains\n")
        f.close()


def main():
    args = parse_args()
    setup_logger(args.log_file)
    classify_chains(
        args.table,
        args.output,
        args.se_model,
        args.me_model,
        raw_out=args.raw_model_out,
        ld_model=args.ld_model,
    )


if __name__ == "__main__":
    main()
