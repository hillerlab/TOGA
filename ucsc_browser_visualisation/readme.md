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
else if (startsWith("HLTOGAannot", trackHubSkipHubName(table)))
    {
    doHillerLabTOGAGene(database, tdb, item, table);
    }
```


[//]: # (```c)
[//]: # (else if (startsWith("HLTOGA", table) && hTableExists(database, "TOGAData")))
[//]: # ({)
[//]: # (    doHillerLabTOGAGene(database, tdb, item, table);)
[//]: # (}
[//]: # (```)

Then re-build your browser.

### Loading TOGA tables

#### Generate tab files

${project_dir}: directory containing TOGA results and intermediate data.

```shell
./ucsc_browser_visualisation/make_bigbed_data_public.py ${project_dir}```
# Wait, it will take a few minutes, most likely less than an hour
```

#### Load tab files to browser database

Transfer created *.bb, *.ix, and *.ixx files to the machine hosting your instance of UCSC genome browser, if need be.
Create a table schema for bigBed tracks, using bigDataUrl field to specify the bigBed URL, and searchTrix field to specify the *.ix file URL (no need to specify *.ixx file separately, it just should be located in the same directory with the *.ix file.)

Please also specify the following fields:
```
type bigBed 12 +
labelFields name
searchIndex name
```
