# TOGA

<img src="https://github.com/hillerlab/TOGA/blob/master/supply/logo.png" width="500">

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![version](https://img.shields.io/badge/version-1.1.7.dev-blue)
[![DOI](https://zenodo.org/badge/277817661.svg)](https://zenodo.org/badge/latestdoi/277817661)
[![License](https://img.shields.io/github/license/hillerlab/TOGA.svg)](https://github.com/hillerlab/TOGA/blob/master/LICENSE)
[![made-with-Nextflow](https://img.shields.io/badge/Made%20with-Nextflow-23aa62.svg)](https://www.nextflow.io/)
[![Published in Science](https://img.shields.io/badge/Published%20in-Science-blue.svg)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10193443/)

TOGA is a new method that integrates gene annotation, inferring orthologs and classifying
genes as intact or lost.

TOGA implements a novel machine learning based paradigm to infer orthologous genes between
related species and to accurately distinguish orthologs from paralogs or processed pseudogenes.

This tutorial explains how to get started using TOGA.
It shows how to install and execute TOGA, and how to handle possible issues that may occur.

For more details, please check out the [TOGA wiki](https://github.com/hillerlab/TOGA/wiki).


## GitHub discussions section

Interested in contributing, have questions, or want to discuss the science behind TOGA?
Head over to our [Discussions section](https://github.com/hillerlab/TOGA/discussions).
It's a new, experimental space (authors did not have a chance to try this GitHub function yet) where we can talk
about anything that doesn't quite fit into the Issues framework.

## Changelog

For a detailed history of changes made to the TOGA project, please refer to the [Changelog](https://github.com/hillerlab/TOGA/blob/master/VersionHistory.md).
This document provides version-specific updates, including new features, bug fixes, and other modifications.

## Installation

TOGA is compatible with Linux and MacOS, including M1-based systems.
It is recommended to use Python version 3.11.

It is highly recommended to have access to computational cluster, but
for small or partial genomes with short genes a desktop PC will be enough.

TOGA requires [Nextflow](https://www.nextflow.io), which in turn requires java >=8.
Check your version of java and install nextflow using one of the following commands:

```shell
curl -fsSL https://get.nextflow.io | bash
# OR
conda install -c bioconda nextflow
```

If you've downloaded nextflow using curl, move the nextflow executable to a
directory accessible by your $PATH variable.

To get TOGA do the following:

```shell
# clone the repository
git clone https://github.com/hillerlab/TOGA.git
cd TOGA
# install necessary python packages:
python3 -m pip install -r requirements.txt --user
# call configure to:
# 1) train xgboost models
# 2) download CESAR2.0
# 3) compile C code
./configure.sh
# run a test, it will take a couple of minutes
./run_test.sh micro
```

If you see something like this at the very end, then TOGA is almost ready to go:

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

If you experience any problems installing TOGA, please visit the [troubleshooting](https://github.com/hillerlab/TOGA/blob/master/TroubleShooting.md) section.

### Configuring TOGA for cluster

TOGA uses nextflow to run cluster-dependent steps.
To run a pipeline on cluster nextflow requires a configuration file defining "executors" component.
This repository contains configuration files for slurm cluster,
please find them in the nextflow_config_files directory.

To create configuration files for non-slurm cluster do the following:

1) Find [here](https://www.nextflow.io/docs/latest/executor.html) what parameters are available for your cluster.
Most likely, you can use slurm configuration files as a reference.
2) Create a separate directory for configuration files, or re-use nextflow_config_files dir.
3) Create "extract_chain_features_config.nf" file.
This file contains configuration for chain features extraction step.
These jobs are expected to be short and not memory consuming, so 1 hour of runtime limit
and 10Gb of memory would be enough.
4) Create "call_cesar_config_template.nf" file.
This configuration file is for CESAR jobs.
These jobs usually take much longer that chain feature extraction, it's recommended to request 24 hours for them.
You don't have to provide an exact amount of memory for these jobs, TOGA will compute this itself.
Please write a placeholder instead, as follows: process.memory = "${\_MEMORY\_}G".

### Final test

This repository also contains sample data to perform a wide-scale test.
To do so, please download genome sequences for human (GRCh38/hg38) and mouse (GRCm38) in the 2bit format.
You can download these 2bit files using the following links:

Human 2bit: wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit

Mouse 2bit: wget https://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.2bit

Then call the following:

```shell
./toga.py test_input/hg38.mm10.chr11.chain test_input/hg38.genCode27.chr11.bed ${path_to_human_2bit} ${path_to_mouse_2bit} --kt --pn test -i supply/hg38.wgEncodeGencodeCompV34.isoforms.txt --nc ${path_to_nextflow_config_dir} --cb 3,5 --cjn 500 --u12 supply/hg38.U12sites.tsv --ms
```

This will take about 20 minutes on 500 cores cluster.

### Troubleshooting

Please see [here](https://github.com/hillerlab/TOGA/blob/master/TroubleShooting.md).

## Usage

This section explains TOGA usage, especially toga.py arguments and input files format.

### Input files

TOGA is a reference-based genome annotation tool, which means that it needs the following data as input:

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
- Make sure that transcript identifiers are unique, e.g. avoid cases where two or more
transcripts have the same identifier.

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
Chain-formatted genome alignments suitable for running TOGA can
be generated using our [make_lastz_chain pipeline](https://github.com/hillerlab/make_lastz_chains).

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

##### --project_dir PROJECT_DIR

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

Chain score threshold.
Exclude chains that have a lower score from the analysis.
Default value is 15000.

##### --no_chain_filter, --ncf

A flag.
Do not filter the chain file (make sure you specified a .chain but not .gz file in this case)

###### --isoforms, -i ISOFORMS_FILE

Path to isoforms data file, details are above.

##### --keep_temp, --kt

A flag.
With this flag TOGA will not remove temporary and intermediate files.
Highly recommended for the first times: it will be easier to trace issues.

##### --limit_to_ref_chrom

Limit analysis to a single reference chromosome.
If you provide whole genome alignment for human and mouse but
would like to perform the analysis on the human chr11 only, then add
"--limit_to_chrom chr11" parameter.

##### --chain_jobs_num, --chn NUMBER_OF_JOBS

Number of cluster jobs for extracting chain features.
Recommended from 20 to 50 jobs.

##### --cesar_jobs_num, --cjn NUMBER_OF_JOBS

Number of CESAR cluster jobs.
Recommended from 300 to 1500, depending on your patience and a number
of available cluster cores.

##### --mask_stops, --ms

A flag.
CESAR cannot process coding sequences containing in-frame stop codons.
However, sometimes we need to project these genes and they are
actually intact, like selenocysteine-coding ones.
With this parameter TOGA will mask stop codons and CESAR could process them.
Without this parameter TOGA will crash if meet any in-frame stop codon in the
reference.
Use this flag only if you are sure that your reference sequences might
contain in-frame stop codons.

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
requested memory, which is likely a unsustainable use of cluster resources.
But if you call TOGA with the following parameter:

```txt
--cesar_buckets 10,100
```

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
So it doesn't annotate query genome, just produces a list
of orthologous chains for each reference gene.

##### --o2o_only, --o2o

Process only the genes that have a single orthologous chain.
Please note that many-2-one orthologs could not be filtered out
at this stage!

##### --no_fpi

A flag.
Consider long frame-preserving indels as inactivating mutations.

##### --parallelization_strategy STRATEGY, --ps STRATEGY

This option lets you choose the strategy for running jobs in parallel.

📅 As of September 9, 2023, the available strategies are:

* `nextflow` (default)
* `para` (internal Hillerlab script to manage sbatch)
* `custom`
* 
🛠 Custom Strategy:
If you're feeling adventurous, you can roll your own strategy.
Just dive into `parallel_jobs_manager.py` and implement the required methods.

🔮 Coming Soon:
We're planning to add `snakemake` support in future releases.

##### --nextflow_dir NEXTFLOW_DIR, --nd NEXTFLOW_DIR

Nextflow working directory: from this directory
nextflow is executed, also there all nextflow log
files are kept

##### --nextflow_config_dir NEXTFLOW_CONFIG_DIR, --nc NEXTFLOW_CONFIG_DIR

Directory containing nextflow configuration files for
cluster, pls see nextflow_config_files/readme.txt for
details.

##### --quiet, -q

Run without printing messages to console. 
This does not affect the writing of log _files_: they will be written as usual. 

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

TOGA produces 3 fasta files: prot.fasta, codon.fasta nucleotide.fasta.
It saves both the reference and predicted query sequences.
- prot.fasta contains protein sequences of reference genes and predicted transcripts.
- codon.fasta contains codon alignments, corrected for frameshifting insertions and deletions
- nucleotide.fasta contains exon nucleotide alignments

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
3) PM - partial & missing. Most of the projection lies outside scaffold borders.
4) M - missing, assembly gaps mask >50% of the prediction CDS.
5) L - clearly lost.
6) UL - "uncertain loss", there are inactivating mutations but not enough evidence
for "clearly lost" class. In other words: neither lost nor intact.
7) PI - partially intact: some fraction of CDS is missing, but most likely
this is intact.
8) I - clearly intact.

### proc_pseudogenes.bed

Bed-formatted file containing annotation of processed pseudogenes in the query.

### genes_rejection_reason.tsv

If TOGA skips any gene, transcript or projection, it writes about that in this file.
Also this file shows a reason, why this happened.

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
It would also be helpful if you plan to plot numerous genes.
You can do it using "mut_index.py" script in the "supply" directory:

```shell
./supply/mut_index.py ${PROJECT_DIR}/inact_mut_data.txt ${PROJECT_DIR}/inact_mut_data.hdf5
```

Then use "./supply/plot_mutations.py" script to create a visualization.
This script requires the following:

1) Reference bed file.
2) List of mutations: either text file of an indexed one.
3) Transcript identifier to plot (or a gene identifier, see below)
4) Path to output file: script creates a svg figure

For example:

```shell
/supply/plot_mutations.py ${REFERENCE_BED_FILE} ${PROJECT_DIR}/inact_mut_data.hdf5 ENST0000011111 test.svg
```

This will create a plot of all inactivating mutations detected for all projections of the ENST0000011111
transcript.

If you like to visualize all projections of a gene then:

1) Provide gene name instead of transcript ID
2) Also provide isoforms file using --isoforms_file parameter.

For example:

```shell
/supply/plot_mutations.py ${REFERENCE_BED_FILE} ${PROJECT_DIR}/inact_mut_data.hdf5 ENSG0000011111 test.svg -i ${ISOFORMS FILE}
```

The script will look for all transcripts of the ENSG0000011111 gene in the isoforms file.
Then it will make a plot for each of these transcripts.

If you are interested only in a particular projection then
provide the chain of interest with --chain parameter:

```shell
/supply/plot_mutations.py ${REFERENCE_BED_FILE} ${PROJECT_DIR}/inact_mut_data.hdf5 ENST0000011111 test.svg --chain 222
```

## Getting assembly quality statistics

TOGA also provides a powerful way to benchmark assembly completeness and quality. TOGA’s gene classification explicitly distinguishes between genes with missing sequences (indicative of assembly incompleteness) and genes with inactivating mutations (an excess of genes with inactivating mutations indicates a higher base error rate). We used a set of 18430 ancestral placental mammal genes (file is [here](https://github.com/hillerlab/TOGA/blob/master/TOGAInput/human_hg38/Ancestral_placental.txt)) to compute assembly quality statistics, as illustrated previously for the [vampire bat](https://www.science.org/doi/10.1126/sciadv.abm6494) or [Rhinolophid bat](https://europepmc.org/article/ppr/ppr616538) genomes.  

To use TOGA to benchmark assembly quality, use the "supply/TOGA_assemblyStats.py" script. This script produces a TSV file with a summary of the number of genes that are intact (classified as I), having missing sequence (TOGA status PI, M, PM, PG, or absent) or inactivating mutations (L and UL). The script will also generate a PDF image of a stacked plot of the statistics (used in Figure 1 in our previous studies). The two output files will be named respectively ${TOGA_DIRS_FILE}\_stats.tsv and ${TOGA_DIRS_FILE}\_statsplot.pdf

The input is a text file listing the names of the TOGA directories. Each directory must contain the loss_summ_data.tsv output file. 

Use this command to consider all genes.
```shell
./TOGA_assemblyStats.py ${TOGA_DIRS_FILE} -m stats
```

You can restrict the statistics to a subset of genes (for example, the Placental ancestral gene set in "./TOGAInput/human_hg38/Ancestral_placental.txt") by providing the -ances/--ancestral parameter with the path to the file listing them like this:
```shell
./TOGA_assemblyStats.py ${TOGA_DIRS_FILE} -m stats -ances ./TOGAInput/human_hg38/Ancestral_placental.txt
```

If you want the script to provide the details of all the classes, without grouping for example L and UL under "genes with inactivating mutations", simply add the -d/--detailed flag.

If you want to change the names of the TOGA directories to something more meaningful (e.g. instead of "vs_speNam" listing the species common or latin name), you can provide a TSV file mapping the names listed in the ${TOGA_DIRS_FILE} to new ones with the -aN/--assemblyNames parameter like this:
```shell
./TOGA_assemblyStats.py ${TOGA_DIRS_FILE} -m stats -ances ./TOGAInput/human_hg38/Ancestral_placental.txt -aN ${NAME_MAP_FILE}
```

#### Assembly statistics for haplotype resolved assemblies

If you have multiple assemblies of the same species, or you have haplotype-resolved assemblies, where e.g. sex chromosomes are only contained in one assembly, 
you may want to combine the results to get a more accurate view. For example, X-linked genes will be deemed missing (M) by TOGA in the Haplotype with Y-chromosome.
Similarly, in case of a polymorphic frameshift or a base error, a gene may be classified as 'uncertain loss' (UL) in haplotype 1 but Intact in haplotype 2.

TOGA_assemblyStats.py, can merge different assemblies, using a precedence map: I>PI>UL>L>M>PM>PG>abs, to keep the "best" class for a given gene out of a range of TOGA runs:
```shell
./TOGA_assemblyStats.py ${TOGA_DIRS_FILE} -m merge
```

This will output a file named ${TOGA_DIRS_FILE}\_merge.tsv with the same structure as a loss_summ_data.tsv file for genes and transcripts
(but not projections of course, which are run-specific).

If for some reason you want to change the precedence map (for example, to focus on how many genes are ever part of an assembly gap), 
you can set the parameter -pre/--precedence like so:

```shell
./TOGA_assemblyStats.py ${TOGA_DIRS_FILE} -m merge -pre M#PM#PG#abs#I#PI#UL#L
```

Note that the outputs, with the suffix \_merge.tsv, can be renamed loss_summ_data.tsv and put in a directory, to act as the output of a fictional TOGA run.
This can be useful when using the same script as outlined in the previous section to get summary statistics.

## Contributing

TOGA uses the python built-in [doctest](https://docs.python.org/3/library/doctest.html) framework.
In this framework, tests are embedded within each function's docstring.
To run all the tests in a file, you could run:

```shell
python3 -m doctest toga.py
# Or, for verbose output: python3 -m doctest -v toga.py
```

To run all the tests in all the python files in this project, try:

```shell
python3 -m doctest $(find . -name "*py")
# Or, for verbose output: python3 -m doctest -v $(find . -name "*py") 
```

## Citation

Kirilenko BM, Munegowda C, Osipova E, Jebb D, Sharma V, Blumer M, Morales AE, Ahmed AW, Kontopoulos DG, Hilgers L, Lindblad-Toh K, Karlsson EK, Zoonomia Consortium, Hiller M. [Integrating gene annotation with orthology inference at scale.](https://www.science.org/stoken/author-tokens/ST-1161/full) Science, 380(6643), eabn3107, 2023
