# TOGA plugin for UCSC genome browser

Please note that this part of TOGA framework is at an early development stage.

## Purpose

This plugin handles the click event on a TOGA annotation track.
It creates a custom HTML pages showing the following data for each TOGA projection:

1) General projection info: gene loss classification, chain orthology probability, etc.
2) Detected inactivating mutations.
3) Information per each exon, such as nucleotide %id, BLOSUM score, alignment to reference.
4) Protein sequence alignment.

Indeed, this radically extends the default bed-annotation functionality.

What's planned: switch to genePred-formatter tracks.

## Usage

This guide is built on the premise that you have your own UCSC browser mirror
and are familiar with Kent source code.

### Modify kent source code

Please copy togaClick.h and togaClick.c to the following directory:
${KENTHOME}/src/hg/hgc
Also you need to add the "togaClick.o" target to the following file:
${KENTHOME}/src/hg/hgc/makefile

Makefile O variable should contain "togaClick.o" as shown here:

```makefile
O = $A.o bamClick.o barChartClick.o bigBedClick.o bigDbSnpClick.o \
    ccdsClick.o cgapSageClick.o dbRIP.o \
    encodeClick.o expClick.o \
    gencodeClick.o gtexEqtlClusterClick.o gtexClick.o gvfClick.o hgdpClick.o \
    interactClick.o lowelab.o lrgClick.o \
    mafClick.o makeItemsClick.o mgcClick.o \
    parClick.o peakClusters.o regMotif.o retroClick.o rmskJoinedClick.o \
    rnaFoldClick.o togaClick.o \
    transMapClick.o txCdsInfo.o pubs.o vcfClick.o wiggleClick.o \
    wikiTrack.o variomeClick.o numtsClick.o geneReviewsClick.o
```

Then open the following file:
${KENTHOME}/src/hg/hgc/hgc.c

There are 2 places to change:

First: add #include "togaClick.h" directive, source code should look like this:

```c
#include "aveStats.h"
#include "trix.h"
#include "bPlusTree.h"
#include "customFactory.h"
#include "iupac.h"
#include "togaClick.h"

static char *rootDir = "hgcData";
```

Second: find a conditional expression encapsuled between the following comments:

```c
/* Start of 1000+ line dispatch on table involving 100+ if/elses. */

/* End of 1000+ line dispatch on table involving 100+ if/elses. */
```

There will be quite a lot of repetitive conditionals.
They check whether the track name is equal or starts with some string and then
call the corresponding function to handle this track.
We will add another condition to handle TOGA tracks.

Put another condition as shown here:

```c
else if (startsWith("HLTOGA", table) && hTableExists(database, "TOGAData"))
{
    doHillerLabTOGAGene(tdb, item);
}
```

Then re-build your browser.

### Loading TOGA tables

#### Generate tab files

${project_dir}: directory containing TOGA results and intermediate data.

```shell
./ucsc_browser_visualisation/make_sql_data.py ${project_dir}```
# Wait, it will take a few minutes, most likely less than an hour
```

#### Load tab files to browser database

${tab_files_dir} = ${project_dir}/tabs
${query} - annotated genome identifier.
Use *.sql files located in the ucsc_browser_visualisation directory.
Call the following commands:

```shell
hgLoadSqlTab ${query} TOGAData togaData.sql ${tab_files_dir}/togaInfo.tab
hgLoadSqlTab ${query} TOGANucl togaNucl.sql ${tab_files_dir}/togaNucl.tab
hgLoadSqlTab ${query} TOGAInactMut togaInactMut.sql ${tab_files_dir}/togaInactMut.tab
```

Also load bed file using hgLoadBed command.
Create a track starting with HLTOGA, for instance HLTOGAannotation.

```shell script
hgLoadBed ${query} HLTOGAannotation ${project_dir}/query_annotation.bed
```
