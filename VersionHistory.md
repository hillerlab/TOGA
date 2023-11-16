# TOGA 1.0 #
This was used to produce all data in the original manuscript: 

Kirilenko BM, Munegowda C, Osipova E, Jebb D, Sharma V, Blumer M, Morales AE, Ahmed AW, Kontopoulos DG, Hilgers L, Lindblad-Toh K, Karlsson EK, Zoonomia Consortium, Hiller M.
[Integrating gene annotation with orthology inference at scale](https://www.biorxiv.org/content/10.1101/2022.09.08.507143v1). _Science_, in press, 2022

# TOGA 1.1.0 #


### Improved transcript classification ###
* Improved handling of compensating frameshifts. If two or more frameshifts (such as +1, +1, +1) compensate each other and restore the reading frame without having a stop codon in the frame-shifted region, TOGA does not consider these frameshifts as inactivating. However, the region between these frameshifts is translated into a protein sequence that is likely entirely different from the reference proteins. We observed extreme cases where compensating frameshifts cover essentially the entire transcript, such that considering the transcript as intact is not correct. 
In TOGA 1.1, we now consider the region between compensating frameshifts as deleted. This improves the transcript classification step. For example, the case where compensating frameshifts cover the entire transcript is now classified as lost. 

* Improved classification of N-terminal mutations. In TOGA 1.0, we masked inactivating mutations if they are located in the first 10% of the reading frame. The reason was that often there is a downstream ATG in proximity, but this is not always the case. 
In TOGA 1.1, we now determine where the next inframe ATG codon is located that is downstream of a potential inactivating mutation. If this ATG is located in the first 10% of the reading frame, the mutation is not considered as inactivating (as in 1.0). However, if this ATG is located outside the first 10% of the reading frame, a truncated protein would lack at least 10% of the N-terminus and TOGA 1.1 now considers this mutation as inactivating in the transcript classification step. 

* Improved classification of borderline cases where a projection has considerable amount of both deleted and missing regions.

* Fixed intron retention module bug. Now, TOGA should correctly identify them.

### Codon alignments ###
* The codonAlignment.fasta file previously contained the raw CESAR output, including exons that were classified as deleted or missing in later processing steps. TOGA 1.1 now generates a codonAlignment.fasta where such exons are excluded. 


### bigBed browser tracks ###
* TOGA 1.0 generated bed and tab files that could be loaded as a total of 4 SQL tables to display the query annotation and provide rich information about each transcript via the TOGA handler function upon clicking on a transcript in the UCSC genome browser. TOGA 1.1 utilizes UCSC's bigBed functionality and produces a single bigBed file that contains all that information, enabling a SQL-free TOGA display. The TOGA handler function was updated to work with bigBed. The bigBed TOGA track often displays detailed information faster than the SQL version. 


### Minor changes ###
* Longer lists of several compensating frameshifts (e.g. FS_1,2,3,4,5) could result in warnings when loading SQL tables if the string gets too long. To avoid this, we now write a condensed format FS_1-5. 
* Avoiding overlapping exon annotations in fragmented query assemblies. If a gene is split into several query chains, there were rare cases where several query chains overlap in the query genome, resulting in exon annotations that overlap each other. This creates invalid bed or genePred records. In TOGA 1.1, we now classify such genes as missing and do not output a bed entry for this gene.
* Sometimes a gene had several orthologous chains and one of them is a 'runaway' chain that spans megabases in the query locus, which is not computable in the CESAR step. TOGA 1.1 now ignores the runaway chain but runs CESAR for the other (normal) orthologous chains. This adds a few orthologs to the annotation.
* Added TOGA_assemblyStats.py
* Bug fixes in merge_cesar_output.py and CESAR output parser.
* Minor changes in plot_mutations.py.


### Nomenculature updates ###
* `MIDDLE_80%_INTACT` projection's feature renamed to `MIDDLE_IS_PRESENT`
* `SSM` mutation type (splice site mutation) is split into more precise `SSMD` (donor) and `SSMA` (acceptor)


# TOGA 1.1.1 #

### Critical bugfixes ###

* Fixed codon alignment generation procedure - all codons are 3bp long. Previous version could mistakenly include dinucleoties into the codon alignment. 
* Fixed parsing the protein and codon alignment output. Previous version could add a sequence to protein alignment if its ID contained "PROT".

# TOGA 1.1.2 #

Documentation improvements.

# TOGA 1.1.3 #

* Added support for Nextflow DSL2.
* Enabled generation of the processed pseudogenes annotation track.
* Removed `-f/--fragmented_genome` flag from toga.py CLI - the "assemble transcript from fragments" feature is enabled by default. Please use `--disable_fragments_joining/--dfj` option to disable this feature.
* Fixed package versions in the `requirements.txt`

# TOGA 1.1.4 #

* CESAR_wrapper.py does not create temporary input files for CESAR. Instead, passes the input directly to CESARs stdin.
* Fixed GCC flags for ARM architecture.
* A little improvement in the `run_test.sh` script.
* Better versions management - added version.py

# TOGA 1.1.5 #

* extract_codon_alignment.py: "!" characters inserted my MACSE2.0 to compensate frameshifts are replaced with "N". 
* CESAR_wrapper.py -> does not use `/dev/shm/` or `/tmp/` partitions anymore.
* Significantly improved logging - all logs are automatically stored in the output directory.
* Nomenclature update: identifiers of TOGA-annotated genes now start with "TOGA_" instead of potentially confusing "reg_"
* Nomenclature update: TRANS chain class (considered confusing) renamed to SPAN (better fits spanning chains)
* Added `bed2gtf` submodule to facilitate results post-processing.
* Improved parallelization strategy: using abstract class to handle different ways to parallelize computations.
* (started): better code organisation - constants class
* Replaced CESAR2.0 submodule with Kirilenko's lightweight fork (without Kent, etc.)

# TOGA 1.1.6 (current release) #

* Removed obsolete optimised CESAR path.
* Fixed a minor bug with re-running CESAR jobs memory parameter if no buckets were specified.
* (by shjenkins94) added gene_prefix argument
* (by shjenkins94) fixed U12 argument check bug
* minor code style improvements

# TOGA 1.1.7 (planned) #

* (planned) Step management - now user can rerun TOGA from any desired step.
* (planned) Using CESAR in the single exon mode for most of the exons to solve the RAM bottleneck (maybe TOGA 1.2.x?)
* Added logging for filtered out reference transcripts
* Updated library versions - now compatible with python3.11
* Updated configure.sh -> fixed inability to install CESAR if TOGA is installed from a zip archive + cleanup scenario
* Fixed bug with parsing intron retention in the CESAR output -> in the prev versions, query chars located vs > were skipped
