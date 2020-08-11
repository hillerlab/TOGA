# Chain classification train data

This directory contains:

- default train dataset (human vs rat)
- ipython notebook to create train dataset

## Training dataset structure

Train.tsv is a tab-separated file containing the following fields:

- gene: gene or transcript identifier, string
- gene_overs: number of chains that intersect this gene, positive integer
- chain: chain identifier, positive integer
- synt: chain synteny, positive integer
- gl_score: chain alignment score, positive integer
- gl_exo: global CDS fraction, float
- chain_len: chain length in the reference genome, positive integer
- loc_exo: local CDS fraction, float
- exon_cover: number of exon bases intersecting chain aligning blocks, non-negative integer
- intr_cover: number of intron bases intersecting chain aligning blocks, non-negative integer
- gene_len: gene or transcript length, positive integer
- ex_num: number of exons, positive integer
- ex_fract: number of reference bases in exons, positive integer
- intr_fract: number of reference bases in introns, non-negative integer
- flank_cov: flank fraction, float
- exon_perc: basically exon_cover / ex_fract * 100, float
- chain_len_log: chain legnth in reference log10, float
- synt_log: chain synteny log10, float
- intr_perc: intr_cover / intr_fract * 100, float
- single_exon: 1 (is single exon) or 0 (multi-exon)
- y: target, 1 (ortholog) or 0 (paralog)

One line corresponds to one gene/chain intersection.
Please note that not all of these features are necessary for training the model!

TOGA uses two models: one for single-exon (SE) genes, another for multi-exon (ME) ones.
single_exon feature is mandatory to split the training dataset for SE and ME rows.

Features required to train ME model:

- gl_exo
- loc_exo
- flank_cov
- synt_log
- intr_perc

Features required to train SE model:

- gl_exo
- flank_cov
- exon_perc
- synt_log

## Creating custom training datasets

You can use "create_train.ipynb" notebook as a reference
to create custom training datasets.
Fields in the train.tsv look a bit reduntant, also not all of these features
are used in the model training or prediction.
Please feel free to use any subset of features you like and/or create your own.

To call our default example you will need:

- Reference annotation in bed-12 format
- Dataset containing chain features (to obtain this call toga.py with --sac flag).
this dataset will be saved into $project_dir/chain_results_df.tsv file.
- Biomart output containing orthologs data. An example for rat would contain
the following fields:
  - Gene stable ID
  - Transcript stable ID
  - Rat gene stable ID
  - Rat gene name
  - Rat homology type
  - Rat protein or transcript stable ID

In general, TOGA classifies neither genes nor chains but gene/chain intersections.
So basically the training dataset should be: gene_ID + chain_ID + features + class
For classification most interesting classes are "ortholog" and "paralog".
