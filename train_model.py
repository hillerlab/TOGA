#!/usr/bin/env python3
"""Just create a trained model."""
import xgboost as xgb
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.metrics import mean_squared_error
import joblib
from datetime import datetime as dt
import os

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]


def train_on(X, y, save_to, name=None):
    """Train model on the X and y given."""
    # models parametets
    n_trees = 50
    max_depth = 3
    learning_rate = 0.1
    # create and cross-validate model
    model = xgb.XGBClassifier(n_estimators=n_trees,
                              max_depth=max_depth,
                              learning_rate=learning_rate)
    kfold = StratifiedKFold(n_splits=5, random_state=777, shuffle=True)
    results = cross_val_score(model, X, y, cv=kfold)
    model.fit(X, y)
    if name:
        print(f"{name} model: ")
        print(f"Training on {len(X)} samples")
        print(f"Using features: {X.columns}")
    print("Accuracy: {0:.3f} {1:.3f}".format(results.mean() * 100, results.std() * 100))
    joblib.dump(model, save_to)
    print(f"Model saved to: {save_to}")


# load dataset, defile where is what
t0 = dt.now()
file_location = os.path.dirname(__file__)
models_dir = "models"
train_tsv = "train.tsv"
se_model_dat = "se_model.dat"
me_model_dat = "me_model.dat"

train_set = os.path.join(file_location, models_dir, train_tsv)
se_model_path = os.path.join(file_location, models_dir, se_model_dat)
me_model_path = os.path.join(file_location, models_dir, me_model_dat)
df = pd.read_csv(train_set, header=0, sep="\t")
print(f"Training dataset size: {len(df)}")

# define which columns to drop in SE and ME models
se_drop_in_X = ["gene", "gene_overs", "chain", "synt", "gl_score", 
                "chain_len", "exon_cover", "intr_cover", "gene_len",
                "ex_num", "ex_fract", "intr_fract", "chain_len_log",
                "y", "single_exon", "intr_perc", "loc_exo"]

me_drop_in_X = ["gene", "gene_overs", "chain", "synt", "gl_score", 
                "chain_len", "exon_cover", "intr_cover", "gene_len",
                "ex_num", "ex_fract", "intr_fract", "chain_len_log",
                "y", "exon_perc", "single_exon"]

# get X, y for SE and ME models
df_se = df[df["single_exon"] == 1]
df_me = df[df["single_exon"] == 0]

X_se = df_se.copy()
X_se = X_se.drop(se_drop_in_X, axis=1)
y_se = df_se["y"]

X_me = df_me.copy()
X_me = X_me.drop(me_drop_in_X, axis=1)
y_me = df_me["y"]

print(f"Single exon train length: {len(X_se)}; multi exon: {len(X_me)}")

# train and save models
train_on(X_se, y_se, se_model_path, name="Single exon")
train_on(X_me, y_me, me_model_path, name="Multi exon")

print(f"Done in {dt.now() - t0}")
