# TOGA - Tool to infer Orthologs from Genome Alignments

TOGA is a new method that integrates gene annotation, inferring orthologs and classifying genes as intact or lost.

TOGA implements a novel machine learning based paradigm to infer orthologous genes between related species and to accurately distinguish orthologs from paralogs or processed pseudogenes.

This tutorial explains how to get started using TOGA.
It shows how to install and execute TOGA, and how to handle possible issues that may occur.

## Installation and preparing

TOGA supports both Linux and MacOS systems.
The package was properly tested on Python version 3.6.5. and 3.7.3.

It is highly recommended to have access to computational cluster, but
for small or partial genomes with short genes a desktop PC will be enough.

To get TOGA do the following:

```shell
# clone the repository
git clone https://github.com/hillerlab/TOGA.git
cd TOGA
# install CESAR, compile C code, install necessary libraries
./configure.sh
# run a test, it will take about a minute
./run_test.sh micro
```

If you see something like this at the very end, TOGA is almost ready to go:

```txt
Orthology class sizes:
one2one: 3
Done! Estimated time: 0:01:02.800084
Program finished with exit code 0
JH567521 299723 336583 ENST00000618101.1169 879 + 299723 336583 0,0,200 7 28,923,130,173,200,179,248, 0,1256,6085,6677,19146,21311,36612,
JH567521 463144 506100 ENST00000262455.1169 711 - 463144 506100 0,200,255 8 102,103,142,112,117,58,116,185, 0,1982,30295,31351,36911,38566,41322,42771,
JH567521 395878 449234 ENST00000259400.1169 942 + 395878 449234 0,0,200 7 123,66,226,116,51,87,240, 0,11871,38544,45802,45994,52305,53116,
Success!
```

If something is wrong and you get something different at the end please see the "troubleshooting" section.

If you are planning to run TOGA on cluster and your cluster is managed by SLURM or LSF then do the following:

```shell
git clone https://github.com/hillerlab/Parasol_LSF_Slurm.git
# for a SLURM cluster
ln -s Parasol_LSF_Slurm/paraSlurm.perl para
# for a LSF cluster:
ln -s Parasol_LSF_Slurm/paraLSF.perl para
```

See "Para on cluster" section for details.

This repository also contains sample data to perform a wide-scale test.
To do so, please download genome sequences for human (GRCh38/hg38) and mouse (GRCm38) in the 2bit format.
Then call the following:

```shell
./run_test.sh normal  ${path_to_human_2bit} ${path_to_mouse_2bit}
```

If something doesn't work, most likely you:

- don't have installed berkeley DB properly
- have berkeley DB installation, but python pip cannot link it to the library
- xgboost was installed for GPU but cannot use it
- TOGA cannot find "para"

### Troubleshooting

TOGA's configure script automatically tries to install all dependencies.
If you encounter error messages related to these three dependencies, please see below for help.

1) BerkeleyDB (database itself and python package)
2) XGBoost
3) Para

#### BerkeleyDB

If the configure script gives an error like

```txt
Complete output from command python setup.py egg_info:
Can't find a local Berkeley DB installation.
(suggestion: try the --berkeley-db=/path/to/bsddb option)
Command "python setup.py egg_info" failed with error code 1
```

please have a look at these links for a likely solution:

- for Linux:
  - <https://github.com/DOsinga/deep_learning_cookbook/issues/59>
  - <https://88plug.com/linux/install-berkeley-4-8-db-libs-on-ubuntu-16-04/>
- for MacOS:
  - <https://stackoverflow.com/questions/16003224/installing-bsddb-package-python>
  - <https://github.com/scrapy-plugins/scrapy-deltafetch/issues/33>

#### XGboost

Sometimes xgboost installation with pip doesn't work and shows a message like:

```txt
Command "/usr/bin/python3 -u -c "import setuptools, tokenize;__file__='/genome/scratch/tmp/pip-install-g6qbjl5j/xgboost/setup.py';f=getattr(tokenize, 'open', open)(__file__);code=f.read().replace('\r\n', '\n');f.close();exec(compile(code, __file__, 'exec'))" install --record /genome/scratch/tmp/pip-record-4dhjvr_9/install-record.txt --single-version-externally-managed --compile" failed with error code 1 in /genome/scratch/tmp/pip-install-g6qbjl5j/xgboost/
```

One of solutions to install XGBoost properly is to compile it from sources, as explained here:

<https://xgboost.readthedocs.io/en/latest/build.html>

#### Para on cluster

Two steps of the TOGA pipeline could be run on cluster:

1) Chain features extraction.
2) Realigning of reference exons in the query genome with CESAR2.0.

To execute these jobs in parallel TOGA uses Para,
a wrapper around LSF or Slurm for managing batches of jobs.
You can find para here:

<https://github.com/hillerlab/Parasol_LSF_Slurm>

If your cluster workload manager is SLURM, just create a symlink to the paraSlurm.perl script.
Like this:

```shell
git clone https://github.com/hillerlab/Parasol_LSF_Slurm.git
ln -s Parasol_LSF_Slurm/paraSlurm.perl para
```

In turn, for LSF cluster use paraLSF.perl

If something doesn't work and you like to configure managing cluster jobs yourself then
please have a look at the following functions in the toga.py script:

1) __chain_genes_run: this function pushes cluster jobs to extract chain features
2) __run_cesar_jobs: this one is responsible for calling exon realign jobs

To execute a batch of jobs in parallel TOGA creates a temporary text file
containing commands that might be executed independently, it looks like this:

```txt
./script.py input/part_1.txt output/part_1.txt
./script.py input/part_2.txt output/part_2.txt
./script.py input/part_3.txt output/part_3.txt
./script.py input/part_4.txt output/part_4.txt
...
```

Let's say we call this file "toga_jobs.txt".
Then TOGA calls a command that looks like this:

```shell
para make joblist_name toga_jobs.txt
```

Then TOGA waits until these jobs are done (para process returns 0) and merges the output files.

Please note that joblist_name should be unique for each batch of jobs.

If you have created a symlink to para but TOGA still cannot find it then
add directory with "para" to $PATH.
TOGA uses "para", not "./para" or "paraSlurm.perl".

## Usage

This section explains TOGA usage, especially toga.py arguments and input files format.

### Input files

TOGA is a reference-based genome annotation tool, which means that needs the following data as input:

1) Gene annotation of the reference genome
2) Genome alignment between the reference and query genome(s)
3) Reference and query genome sequences

#### Gene annotation of reference genome

##### Bed-12 file

TOGA accepts a bed12-formatted file as a reference genome annotation.
This file is mandatory for running TOGA.

Please find bed12 format specification under:
<https://genome.ucsc.edu/FAQ/FAQformat.html#format1>

Example for human gene MAP1S which has two transcripts

```txt
user@user$ grep ENST00000544059 supply/hg38.wgEncodeGencodeCompV34.bed
chr19 17720155 17734490 ENST00000544059 0 + 17720390 17734428 0 7 275,102,83,141,2344,236,218, 0,780,3970,4893,5673,13037,14117,
user@user$ grep ENST00000324096 supply/hg38.wgEncodeGencodeCompV34.bed
chr19 17719479 17734513 ENST00000324096 0 + 17719502 17734428 0 7 141,102,83,141,2344,236,241, 0,1456,4646,5569,6349,13713,14793,
```

This repository contains examples of bed12 file for human and mouse:

1) Human genome annotation: supply/hg38.wgEncodeGencodeCompV34.bed
2) Mouse genome annotation: supply/mm10.wgEncodeGencodeCompVM25.bed

Some advice about your reference annotation:

- Please make sure that the length of the CDS of your annotations is divisible by 3.
TOGA will skip transcripts that do not satisfy this criteria.
- This is highly recommended that CDS of your transcripts start with ATG and end with a canonical stop codon.
- Your transcripts are coding, meaning that thickStart and thickEnd are not equal.
TOGA would skip non-coding transcripts.
- Avoid any pseudogenes in the reference annotations.
- Also, try to avoid merged and incomplete transcripts.

##### Optional but highly recommended: Isoform data

One gene can have multiple isoforms.
TOGA can handle more than one isoform per gene, meaning it is not necessary to reduce the
transcript data to the isoform with the longest CDS.
Isoform data is optional, but if available increases annotation completeness and gene loss determination accuracy.
If you do not provide isoforms data, TOGA will treat each transcript in the bed12 file as a separate gene.

Isoforms can be provided to TOGA in a single two-column tab-separated file in the following format:
GeneIdentifier {tab} TranscriptIdentifier
The first line can be a header.

Example: The human gene MAP1S (ENSG00000130479) has 2 isoforms: ENST00000544059 and ENST00000324096.

```txt
user@user$ grep ENSG00000130479 supply/hg38.wgEncodeGencodeCompV34.isoforms.txt
ENSG00000130479 ENST00000544059
ENSG00000130479 ENST00000324096
```

Importantly, all transcripts listed in the bed12 file have to occur in this isoform file, otherwise TOGA throw an error.
Examples for human and mouse Gencode annotations:

1) For human: hg38.wgEncodeGencodeCompV34.isoforms.txt
2) For mouse: mm10.wgEncodeGencodeCompVM25.isoforms.txt

The simplest way to obtain isoforms file is:

- Visit <https://www.ensembl.org/biomart/martview>
- Choose Ensembl Genes N dataset and then [species of interest] genes
- Go to Filters tab, select "gene type" - protein coding
- Go to Attributes tab, select:
  - Gene stable ID
  - Transcript stable ID
  - Uncheck all other marks!
- Download the results as a tsv file

##### U12 introns data

You also can provide data of U12 exons in the reference genome,
it would facilitate gene loss detection process.
However, this is not mandatory.

There are examples for U12 data files:

1) For human: supply/hg38.U12sites.tsv
2) And for mouse: supply/mm10.U12sites.tsv

There are tab-separated files containing 3 columns:

1) First column contains transcript identifier
2) Second one: exon number
3) Third column contains A/D letter, which means "acceptor" or "donor".

If gene loss pipeline detects mutations of the splice sites listed in this file,
it would not consider them as inactivating.

#### Genome alignment

This section explains what "please provide a genome alignment" actually means.

##### Chain file

Chain file is a mandatory (more precisely, essential) for running TOGA.

According to its name, chain file is a text file that describes chains.
Chains represent co-linear local alignments that occur in the same order
on a reference and a query chromosome.
Thus, a collection of chains could represent a whole-genome pairwise alignment.

Here is a more detailed explanation of chains:

<http://genomewiki.ucsc.edu/index.php/Chains_Nets>

You can find chain file format specification here:

<https://genome.ucsc.edu/goldenPath/help/chain.html>

You can also provide gzipped chain file, then please make sure that the filename
ends with ".chain.gz".

Please make sure that each chain you provide has a unique identifier!

#### Reference and query genome sequences

TOGA accepts reference and query genomes sequences in the 2bit format.

Here is the specification for 2bit file format:

<http://genome.ucsc.edu/FAQ/FAQformat.html#format7>.

Please make sure that chromosome/scaffold names in the bed12, chain and 2bit files are consistent.
For example, if the first chromosome 1 is called "chr1" in the 2bit file, it should also be "chr1" in
the bed12 file, not "chromosome_1" or just "1".

### Arguments

To call TOGA use toga.py script.
It accepts the following arguments:

#### Positional, mandatory arguments

1) Chain file containing the genome alignment.
2) Bed file containing reference genome annotation.
3) Path to reference genome 2bit file.
4) Path to query genome 2bit file.

#### Optional arguments

##### -h, --help

Show help message and exit.
Calling ./toga.py without arguments does the same.

##### --project_folder PROJECT_FOLDER

Project directory.
TOGA will save all intermediate and output files exactly in this directory.
If not specified, use CURRENT_DIR/PROJECT_NAME as default (see below).

##### --project_name PROJECT_NAME

If you don't like to provide a full path to the project directory with --project_folder
you can use this parameter.
In this case TOGA will create project directory in the current directory as
"CURRENT_DIR/PROJECT_NAME".
If not provided, TOGA will try to extract the project name from chain filename,
which is not recommended.

##### --min_score MIN_SCORE

Do not consider chains that have a score lower than this threshold.
Default value is 15000.

##### --no_chain_filter, --ncf

A flag.
Do not filter the chain file (be sure you specified a .chain but not .gz file in this case)

###### --isoforms, -i ISOFORMS_FILE

Path to isoforms data file, details are above.

##### --keep_temp, --kt

A flag.
With this flag TOGA will not remove temporary and intermediate files.
Highly recommended for the first times: it will be easier to trace issues.

##### --no_para

A flag.
If you have access to a cluster but don't like to use it for some reason,
with this flag TOGA would use it as an ordinary PC.

##### --limit_to_chrom

Limit analysis to a single reference chromosome.
If you provide whole genome alignment for human and mouse but
would like to perform the analysis on the chr11 only, then add
"--limit_to_chrom chr11" parameter.

##### --chain_jobs_num, --chn NUMBER_OF_JOBS

Number of cluster jobs for extracting chain features.
Recommended from 10 to 50 jobs.

##### --cesar_jobs_num, --cjn NUMBER_OF_JOBS

Number of CESAR cluster jobs.
Recommended from 300 to 1500, depending on your patience and cluster size.

##### --mask_stops, --ms

A flag.
CESAR cannot process coding sequences containing in-frame stop codons.
However, sometimes we need to project these genes and they are intact,
like selenocysteine-coding ones.
With this parameter TOGA will mask stop codons and CESAR could process them.
Without this parameter TOGA will crash if meet any in-frame stop codon in the
reference.
Use this flag only if you are sure that your reference sequences might
contain inframe-stop codons.

##### --cesar_mem_limit CESAR_MEM_LIMIT

CESAR2.0 might be very memory-consuming.
For most of the genes it requires less than 5Gb of RAM, but some exceptionally long
genes might require more than 100Gb.

TOGA will skip CESAR jobs that would require more that CESAR_MEM_LIMIT GB of memory.
Of course, you will find such genes in the log.

##### --cesar_buckets, --cb CESAR_BUCKETS

Where CESAR_BUCKETS is a comma-separated list of integers, such as "5,15,50".
This is the evolution of previous parameter.
You can split CESAR jobs into different buckets according to their memory requirements.
Imagine you need to call CESAR for 1000 genes.
900 of the cluster jobs will require less than 10Gb of RAM and the rest
require from 10 to 100Gb.
For sure you can split these genes into 100 cluster jobs and require 100Gb
of RAM for each of them.
But then 90% of the time your jobs will consume just a tiny bit of the
requested memory.
But if you call TOGA with the following parameter:

"--cesar_buckets 10,100"

then TOGA will create and push separately two joblists:
one for cluster jobs that require less then 10Gb,
and the second one for jobs that require from 10 to 100Gb.

It will automatically set --cesar_mem_limit, in case you provided
"5,15,50" as --cesar_buckets, the --cesar_mem_limit would be 50Gb.

###### --cesar_chain_limit CESAR_CHAIN_LIMIT

Skip genes that have more that CESAR_CHAIN_LIMIT orthologous chains.
Recommended values are a 50-100.

##### --u12 U12

Path to U12 introns data.

##### --stop_at_chain_class, --sac

A flag.
If set, TOGA halts after chain classification step.

##### --o2o_only, --o2o

Process only the genes that have a single orthologous chain.
Please not that many-2-one orthologs could not be filtered out
at this stage!

##### --no_fpi

A flag.
Consider long frame-preserving indels as inactivating mutations.

## Output reading

This section describes TOGA output.

### Transcript naming convention

To predict a transcript in the query genome TOGA projects a reference transcript
via an orthologous chain.
Since each reference transcript might have more than one orthologous chains,
TOGA uses the following convention for naming predicted transcripts:
"reference_transcript_name"."orthologous_chain_id".

For example, if TOGA identified two orthologous chains (1 and 2) for a transcript A,
you will find two transcripts in the query genome: A.1 and A.2.
If you see two identical annotations in the query genome named as A.1 and B.2, it means
that predicted transcript is an ortholog to reference transcripts A and B.

In the code predicted transcripts are also called "projections".

### query_annotation.bed

Bed12 formatted file containing annotation tracks for the query genome.
Could be loaded to UCSC genome browser.

Predictions have different colors according to gene loss pipeline classes
(please see "loss_summ_data.tsv" section for details).
The color code is:

1) Black - N, no data
2) Brown - PG, paralogous projection
3) Grey - M and PM, missing and partially missing
4) Red - L, clearly lost.
5) Salmon - "grey", neither intact nor lost.
6) Light blue - PI, partially intact.
7) Blue - I, intact.

### query_isoforms.tsv

Isoforms file for the query.

### Fasta files

TOGA produces 2 fasta files: prot.fasta and nucleotide.fasta.
It saves both the reference and predicted query sequences.

### orthology_classification.tsv

File containing orthology data.
It contains 5 columns:

1) t_gene - gene name in the reference
2) t_transcript - transcript identifier in the reference
3) q_gene - gene identifier in the query
4) q_transcript - transcript identifier in the query
5) orthology_class - class of orthology relationships, such as:
one2one, one2many, many2one, many2many and one2zero.

### loss_summ_data.tsv

This file contains gene loss pipeline classification for each
projection, transcript and gene.
TOGA identifies the following classes:

1) N - no data due to technical reasons (such as CESAR memory requirements)
2) PG - no orthologous chains identified, TOGA projected transcripts via paralogous
chains and cannot make any conclusion.
3) PM - partially missing. Most of the projection lies outside scaffold borders.
4) L - clearly lost.
5) M - missing, assembly gaps mask >50% of the prediction CDS.
6) G - "grey", there are inactivating mutations but not enough evidence
for being list.
7) PI - partially intact: some fraction of CDS is missing, but most likely
this is intact.
8) I - clearly intact.

## Inactivating mutations visualization

There is a possibility to visualise inactivating mutations detected
in the projected transcript.
There are 3 levels available:

1) Projection level: visualize one projection only, which is named as
"transcript_id"."chain_id".
2) Visualize all projections of a particular transcript. If a transcript has only
one orthologous chain, that it's the same with level 1. But if TOGA identified
several orthologous chains you can visualize them all at once.
3) Visualize the entire gene. If a gene has several transcripts, then it makes sense.

To make visualisations you need:

Merge all inactivating mutations data into a single file.

```shell
cat ${PROJECT_DIR}/inact_mut_data/* > ${PROJECT_DIR}/inact_mut_data.txt
```

Since the merged file could be huge you can index it.
It could also be helpful if you like to plot numerous genes.
You can do it using "mut_index.py" script in the "supply" directory:

```shell
./supply/mut_index.py ${PROJECT_DIR}/inact_mut_data.txt ${PROJECT_DIR}/inact_mut_data.bdb
```

Then use "./supply/plot_mutations.py" script to create a visualization.
This script requires the following:

1) Reference bed file.
2) List of mutations: either text file of an indexed one.
3) Transcript identifier to plot (or a gene identifier, see below)
4) Path to output file: script creates a svg figure

For example:

```shell
/supply/plot_mutations.py ${REFERENCE_BED_FILE} ${PROJECT_DIR}/inact_mut_data.bdb ENST0000011111 test.svg
```

This will create a plot of all inactivating mutations detected for all projections of the ENST0000011111
transcript.

If you like to visualize all projections of a gene then:

1) Provide gene name instead of transcript ID
2) Also provide isoforms file using --isoforms_file parameter.

For example:

```shell
/supply/plot_mutations.py ${REFERENCE_BED_FILE} ${PROJECT_DIR}/inact_mut_data.bdb ENSG0000011111 test.svg -i ${ISOFORMS FILE}
```

The script will look for all transcripts of the ENSG0000011111 gene in the isoforms file.
Then it will make a plot for each of these transcripts.

If you are interested only in a particular projection then
provide the chain of interest with --chain parameter:

```shell
/supply/plot_mutations.py ${REFERENCE_BED_FILE} ${PROJECT_DIR}/inact_mut_data.bdb ENST0000011111 test.svg --chain 222
```

## Citation

Kirilenko et al, 2020.
In preparation.
