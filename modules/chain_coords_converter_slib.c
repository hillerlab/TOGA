/*
To be compiled as a shared library
chain_coords_converter function
For a chain and a list of regions in the reference genome
make a list of corresponding regions in the query genome
Briefly: project regions through a chain

Input:
chain -> string containing chain (header + blocks)
shift -> integer > 0: number of flanking blocks; the bigger is shift
         the bigger are corresponding regions in the query
regions num -> integer > 0, number of reference genome regions
granges -> list of strings, each string: reference genome region

Author: Bogdan Kirilenko, 2020;
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include "chain.h"

#define MAXCHAR 255
#define DEFAULT_SHIFT 2
#define NUM_BASE 10


char ** chain_coords_converter(char *chain, int shift, int regions_num, char **granges)
{
    // parse region from argv[2]; looks like chromN:start-end
    struct Regions t_regions[regions_num];

    for (int i = 0; i < regions_num; i++)
    {
        // parse regions; their format should be chrom:start-end
        int reg_index = i;
        char *current_region = granges[i];
        // format requires 2 splits: by : and then by -
        char *reg_split = strtok(current_region, ":");
        strcpy(t_regions[reg_index].chrom, reg_split);
        reg_split = strtok(NULL, ":");

        char *split_range = strtok(reg_split, "-");
        t_regions[reg_index].start = strtol(split_range, NULL, NUM_BASE);
        split_range = strtok(NULL, "-");
        t_regions[reg_index].end = strtol(split_range, NULL, NUM_BASE);
    }  // end regions parsing

    for (int i = 0; i<regions_num; i++)
    {
        // check that start < end in each region
        // if that not true -> switch them
        if (t_regions[i].start > t_regions[i].end)
        {
            // this is an extraordinary situation, so let user know about this
            fprintf(stderr, "Warning! End < Start in the region num %d! Switching...\n", i);
            int temp = t_regions[i].start;
            t_regions[i].start = t_regions[i].end;
            t_regions[i].end = temp;
        }
    }  // end regions checking

    // global search iterables
    // about most of them: they are arrays, one per region
    bool head_caught = false;  // flag; true if header was parsed
    bool starts_found[regions_num];
    bool ends_found[regions_num];
    bool finished[regions_num];

    // need to do memset here
    memset(starts_found, false, regions_num * sizeof(bool));
    memset(ends_found, false, regions_num * sizeof(bool));
    memset(finished, false, regions_num * sizeof(bool));

    struct Chain_info head;  // parse head information there

    // variables to collect the results I need
    int starts_in_q[regions_num];
    int ends_in_q[regions_num];
    memset(starts_in_q, 0, regions_num * sizeof(int));
    memset(ends_in_q, 0, regions_num * sizeof(int));

    // hold previous SHIFT + 1 values and NEXT SHIFT + 1 values
    int *starts_buff = (int*)calloc((shift + 2), sizeof(int));
    int *ks = (int*)calloc(regions_num, sizeof(int));

    // variables describing each block
    // intermediate loop values
    // variables called as pointers BUT
    // they are not pointers!
    int t_start_pointer = 0;
    int q_start_pointer = 0;
    int blockTStarts = 0;
    int blockTEnds = 0;
    int blockQStarts = 0;
    int blockQEnds = 0;
    int prevBlockQEnds = 0;

    // this is how we can read a string
    // line-by-line (lines separated by \n)
    char * curLine = chain;

    while(curLine)
    {
        char * nextLine = strchr(curLine, '\n');
        if (nextLine)
        {
            *nextLine = '\0';
        }

        // parse head if it wasn't
        // this should happen at the first iteration
        if (!head_caught)
        {
            // chain head should be caught only once
            // because we operate with a single chain
            // which is also produced by TOGA
            head_caught = true;
            head = parse_head(curLine);
            t_start_pointer = head.tStart;
            q_start_pointer = head.qStart;

            prevBlockQEnds = head.qStart;  // init previous block end with the chain's start

            for (int i = 0; i < regions_num; i++)
            // go region-by-region
            {
                if (strcmp(t_regions[i].chrom, head.tName) != 0)
                {
                    // if the range starts with chrN, but chain is aligned to the chrM,
                    // it is wrong, so those regions must be skipped
                    // likely an error, it also might happen if some gene identifier
                    // belongs to some different chromosomes
                    fprintf(stderr, "Error! Requested range num %d lies in the chromosome %s "
                                    "meanwhile the chain covers chrom %s in target genome\n",
                            i, t_regions[i].chrom, head.tName);
                    starts_found[i] = true;
                    ends_found[i] = true;
                    starts_in_q[i] = 0;
                    ends_in_q[i] = 0;

                // it the region lies outside the chain, it also should be noticed
                } else if (t_regions[i].end < head.tStart) {
                    fprintf(stderr, "Warning! Query %d lies outside the chain borders (to the left)\n", i);
                    starts_found[i] = true;
                    ends_found[i] = true;
                    starts_in_q[i] = (!head.qStrand) ? head.qSize - head.qEnd : head.qStart;
                    ends_in_q[i] = starts_in_q[i] + 1;
                    // ends_in_q[i] = (!head.qStrand) ? head.qSize - head.qStart : head.qEnd;
                
                // the same situation here
                } else if (t_regions[i].start > head.tEnd) {
                    fprintf(stderr, "Warning! Query %d lies outside the chain borders (to the right)\n", i);
                    starts_found[i] = true;
                    ends_found[i] = true;
                    // start_in_q = (!head.qStrand) ? head.qSize - head.qEnd : head.qStart;
                    ends_in_q[i] = (!head.qStrand) ? head.qSize - head.qStart : head.qEnd;
                    starts_in_q[i] = ends_in_q[i] - 1;
                }

            }  // end region checking
        }  // end chain header parsing

        // each chain ends with double /n
        else if (strcmp(curLine, "\n") == 0) 
        {
            // if the chain is over -> break
            // because the input consists on a single chain
            break;
        }

        // now we're parsing each block
        struct Block block = parse_block(curLine);
        // convert relative block coordinates to
        // absolute coords
        blockTStarts = t_start_pointer;
        blockQStarts = q_start_pointer;
        blockTEnds = blockTStarts + block.size;
        blockQEnds = blockQStarts + block.size;
        t_start_pointer = blockTEnds + block.dt;
        q_start_pointer = blockQEnds + block.dq;

        // then restore newline-char, just to be tidy
        if (nextLine) *nextLine = '\n';
        curLine = nextLine ? (nextLine + 1) : NULL;

        // need to make >= 1 not > 0 because unsigned
        // otherwise unsigned int behaves weird
        for (int k = shift + 1; k >= 1; k--)
        {
            starts_buff[k] = starts_buff[k - 1];
        }

        starts_buff[0] = blockQStarts;
        // all_finished == true is default value
        // it could be changed further
        bool all_finished = true;

        // then intersect this chain block with each region
        for (int i = 0; i < regions_num; i++)
        {
            // if the region was mapped -> skip it
            // if not -> all_finished flag to false
            if (finished[i]) {continue;} else {all_finished = false;}

            if ((!starts_found[i]) && (blockTStarts > t_regions[i].start))
            {  // start in between of blocks
                // block-------T.start----------------block
                starts_found[i] = true;
                starts_in_q[i] = (shift == 0) ? blockQStarts : starts_buff[shift + 1];
                // starts_in_q[i] = blockQStarts;
                // block------block-----block--TSTART--block
                // _____ - get this block
            }
            else if ((!starts_found[i]) && (blockTEnds > t_regions[i].start))
            {  // start intersects the block
                starts_found[i] = true;
                int delta = t_regions[i].start - blockTStarts;
                starts_in_q[i] = (shift == 0) ? blockQStarts + delta : starts_buff[shift];

            }  // do not search the end before we find the start

            // parse ends
            // if shift = N is > 0 then we need to add +/- N blocks
            // up and downstream
            if (shift > 0)
            {
                if (starts_found[i] && !ends_found[i] && blockTEnds > t_regions[i].end)
                {
                    ends_found[i] = true;
                    ks[i] = 0;  // I will use this pointer to get this position + SHIFT
                }

                if (ends_found[i])
                {
                    if (ks[i] >= shift)
                    {
                        ends_in_q[i] = blockQEnds;
                        finished[i] = true;
                    }
                    else {ks[i] += 1;}
                }  // I need index of end found + SHIFT
            }
            else
            {
                // a bit more complicated in case if no shift required
                if (starts_found[i] && !ends_found[i] && blockTStarts > t_regions[i].end)
                {  // reference region ends between the chain blocks
                    ends_found[i] = true;
                    ends_in_q[i] = prevBlockQEnds;
                    finished[i] = true;
                // get the end of the previous block
                }
                else if (starts_found[i] && !ends_found[i] && blockTEnds >= t_regions[i].end)
                {
                    // reference region ends inside a block
                    ends_found[i] = true;
                    int delta = t_regions[i].end - blockTStarts;
                    ends_in_q[i] = blockQStarts + delta;
                    finished[i] = true;
                }
            }
        }  // end mapping this block with regions
        // current block will be the previous in the next iteration
        prevBlockQEnds = blockQEnds;

        // if all blocks are finished -> nothing could switch this flag to false
        if (all_finished) {break;}
    }

    // post-process
    for (int i = 0; i < regions_num; i++)
    {
        if (ends_in_q[i] == 0)
        { 
            // we didn't reach +SHIFT block
            // then assign it to the chain end
            ends_in_q[i] = head.qEnd;
        }
        if (starts_in_q[i] == 0)
        {
            // out of borders too
            // assign block to the chain start
            starts_in_q[i] = head.qStart;
        }
        if (!head.qStrand)
        {   
            // if strand is negative we need to reverse coordinates
            int temp = starts_in_q[i];
            // not the most elegant way, yes
            starts_in_q[i] = head.qSize - ends_in_q[i];
            ends_in_q[i] = head.qSize - temp;
        }
    }

    // convert 1\0 (strand) to +\- (human-readable)
    char *tStrand = (head.tStrand == 1) ? "+" : "-";
    char *qStrand = (head.qStrand == 1) ? "+" : "-";
    
    // allocate memory for answer -> array of strings
    char ** answer = malloc (sizeof (char *) * (regions_num + 1));
    for (int i = 0; i < regions_num + 1; i++)
    {
        answer[i] = malloc(MAXCHAR * sizeof(char));
    }

    // save general chain info in the first line of answer
    sprintf(answer[0], "chain %s %s %d %d %d %s %s %d %d %d\n",
            head.tName, tStrand, head.tSize, head.tStart, head.tEnd,
            head.qName, qStrand, head.qSize, head.qStart, head.qEnd);

    for (int i = 0; i < regions_num; i++)
    {
        // write regions strings to answer array
        // each line contains:
        // region number (ordered as in the granges array)
        // region in the reference
        // corresponding region in the query
        // fields are tab-separated
        sprintf(answer[i + 1],
                "%d\t%s:%d-%d\t%s:%d-%d\n",
                i, head.qName, starts_in_q[i], ends_in_q[i],
                t_regions[i].chrom, t_regions[i].start,
                t_regions[i].end);
    }
    return answer;
}

int main()
{
    // should be compiled with -fPIC and -shared flags!
    printf("Warning! This code is not intended to be compiled as a standalone tool.\n");
    printf("Please compile this file with -fPIC and -shared flags.\n");
    printf("and create the /modules/chain_coords_converter_slib.so shared library file.\n");
    return 0;
}
