#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 00:25:17 2021

@author: aahmed
"""

import sys
import pandas as pd
import tqdm
import matplotlib.pyplot as plt
import matplotlib
import argparse


# set these parameters to 42 to make the text editable
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
new_rc_params = {'text.usetex': False, "svg.fonttype": 'none'}
matplotlib.rcParams.update(new_rc_params)


####get arguments

HELP=f"""USAGE: ONE of the following

    To keep the best class found for each gene or transcript
            {sys.argv[0]} <assemblies_file> -m merge  (optional: -ances <ancestral_gene_file>, -pre <TOGAclass1#TOGAclass2#TOGAclass3...> )

    To make the summary statistics of TOGA classification
            {sys.argv[0]} <assemblies_file> -m stats -ances <ancestral_gene_file>  (optional: -aN <assembly_names_file>, -d)

            -pre is I#PI#UL#L#M#PM#PG#abs by default, it's to change the order in which classes are considered
            -ances is a file where each line is a gene you want to keep
            -aN is to specify the names of the assemblies, otherwise TOGA directory names will be used
            -d is to display the statistics for each class TOGA recognizes

          """

parser = argparse.ArgumentParser(usage=HELP)
parser.add_argument('file')
parser.add_argument('-m', '--mode')
parser.add_argument('-ances', '--ancestral',default=False)
parser.add_argument('-pre', '--precedence',default="I#PI#UL#L#M#PM#PG#abs")
parser.add_argument('-aN', '--assemblyNames',default=False)
parser.add_argument('-d', '--detailed',action='store_true')  # on/off flag
ARGS = parser.parse_args()



####get file name for saving
if len(ARGS.file.split("."))<2:
    filename=ARGS.file
else:
    filename=".".join(ARGS.file.split(".")[0:-1])

####create list of species then dataframe of classes
ASSEMBLIES=[]
with open(ARGS.file,"r") as f:
    for line in f:
        if line.startswith("vs_"):
            ASSEMBLIES.append(line.rstrip())
        else:
            ASSEMBLIES.append("vs_"+line.rstrip())

NAMES={x:x for x in ASSEMBLIES}
if ARGS.assemblyNames!=False:
    with open(ARGS.assemblyNames,"r") as f:
        for line in f:
            line=line.rstrip().split("\t")
            NAMES["vs_"+line[0]]=line[1]


labels=set()
def get_classes(X,typ):
    new={}

    for x in tqdm.tqdm(X):
        temp={}
        try:
            df=pd.read_csv(x+"/loss_summ_data.tsv",sep="\t",header=0,index_col=0).loc[typ.upper()]
            for l in range(len(df.iloc[:,0])):
                temp[df.iloc[l,0]]=df.iloc[l,1]
                labels.add(df.iloc[l,1])
            new[x]=temp
        except:
            print("no "+typ+"s in "+x)

    new=pd.DataFrame(new)
    new=new.fillna("abs")
    for x in new.index:
        for m in new.columns:
            labels.add(new.loc[x,m])

    return new

print("Getting TOGA classes\n")
CLASSES=get_classes(ASSEMBLIES,"gene")
if ARGS.mode=="merge":
    CLASSES_trans=get_classes(ASSEMBLIES,"transcript")

ANCESTRAL=[]
if ARGS.ancestral!=False:
    with open(ARGS.ancestral,"r") as f:
        for line in f:
            if line.rstrip() in CLASSES.index:
                ANCESTRAL.append(line.rstrip())
    CLASSES=CLASSES.loc[ANCESTRAL]

####process dataframe
def merge(df,typ):
    l=ARGS.precedence.split("#")
    l={x:l[x] for x in range(len(l))}
    l[len(l)]="abs"

    order={l[x]:x for x in range(len(l))}

    inner={x:{"status":len(order)+1} for x in df.index}
    for x in tqdm.tqdm(df.index):
        for m in df.columns:
            if order[df.loc[x,m]]<inner[x]["status"]:
                inner[x]["status"]=order[df.loc[x,m]]
    inner=pd.DataFrame(inner).replace(to_replace=l).T

    #inner.to_csv(filename+"_merge.tsv",sep="\t")

    stats={x:0 for x in set(inner["status"])}
    for x in inner.index:
        stats[inner.loc[x,"status"]]+=1
    print("Merge finished, the numbers per class for "+typ+" are:")
    for x in stats:
        print(x,stats[x])
    return inner


def stats(df):
    stats={x:{"L":0,"UL":0,"M":0,"PM":0,"I":0,"PI":0,"abs":0,"PG":0} for x in df.columns}
    for x in df.index:
        for m in df.columns:
            stats[m][df.loc[x,m]]+=1
    #stats={NAMES[x]:stats[x] for x in stats}
    stats=pd.DataFrame(stats).T
    stats=stats.loc[ASSEMBLIES]

    if ARGS.detailed:
        stats.rename(index=NAMES).to_csv(filename+"_stats.tsv",sep="\t")
        stats.loc[list(stats.index)[::-1],["I","PI","UL","L","M","PM","PG","abs"]].rename(index=NAMES).plot.barh(stacked=True,color={"I":"#0247FE","PI":"#7ABFFF","UL":"#FF823B","M":"#B2BEB5","L":"#CE1E16","PM":"#A2A2D0","PG":"#000000","abs":"#FEFEFA"},rot=0).get_figure().savefig(filename+"_statsplot.pdf", bbox_inches="tight")
    else:
        stats["genes with missing sequence"]=stats.loc[:,["PI","M","PM","PG","abs"]].sum(axis=1)
        stats["intact genes"]=stats["I"]
        stats["genes with inactivating mutations"]=stats.loc[:,["L","UL"]].sum(axis=1)

        stats=stats.loc[:,["intact genes","genes with inactivating mutations","genes with missing sequence"]]
        stats.rename(index=NAMES).to_csv(filename+"_stats.tsv",sep="\t")

        stats.loc[list(stats.index)[::-1]].rename(index=NAMES).plot.barh(stacked=True,color={"intact genes":"#0247FE","genes with inactivating mutations":"#FE7B15","genes with missing sequence":"#BFBFBF"},rot=0).get_figure().savefig(filename+"_statsplot.pdf", bbox_inches="tight")


if ARGS.mode=="stats":
    print("Calculating statistics")
    stats(CLASSES)
elif ARGS.mode=="merge":
    print("Merging in progress")
    genes=merge(CLASSES, "genes")
    transcripts=merge(CLASSES_trans, "transcripts")
    with open(filename+"_merge.tsv","w") as f:
        for n in genes.index:
            f.write("\t".join(["GENE",n,genes.loc[n,"status"]])+"\n")
        for n in transcripts.index:
            f.write("\t".join(["TRANSCRIPT",n,transcripts.loc[n,"status"]])+"\n")

        

