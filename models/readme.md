This directory contains default XGBoost model, trained on human vs mouse dataset.
Model is saved in the "model.dat" file.

However, it is possible to create custom models.
Jupyter notebook "create_train.ipynb" contains all necessary commands for that.
To make your own custom model you need the following:

- merged_data.tsv - intermediate TOGA file containing chain features.
Call toga with --stop_at_chain_class, then it will produce this file and stop.
- reference genome annotation track in bed12 format
- Ensembl orthologs data downloaded from Biomart

