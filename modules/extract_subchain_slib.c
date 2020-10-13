/*
To be compiled as a shared library.
Function to be called: extract_subchain

Extract chain blocks that intersect the requested region.
You can request subchain for the target or the query genome.
Usage: extract_subchain [chain] [q/t] [target/query genome range]

Arguments:
chain -> string, a chain itself
mode_cd -> char, "q" or "t" -> stands or "query" or "target" genome
search_region -> string formatted as chrom:start-end
                 Region in the query or reference (depending on the mode param)
                 Function extracts chain blocks that intersect this region

Output structure:
list of strings
each string except the last one contains the following:
    chain block start and end in the reference, chain block start and end in the query
    values are space-separated
    each of those chain blocks intersect the region of interest (which is specified)
The last string is "END"

Author: Bogdan Kirilenko, 2020;
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <stdint.h>
#include "chain.h"
#define MAXCHAR 255
#define ALLOC_START 500
#define ALLOC_EXT 75
#define MEM_CHUNK 5


char ** extract_subchain(char *chain, char *mode_ch, char *search_region)
{
    // read the chain file OR stdin
    // read genomic range mode
    bool mode;
    // the region we are interested in might lie in the target
    // as well as in the query genome
    if (strcmp(mode_ch, "t") == 0) {
        mode = true;  // meaning it is target genome
    } else if (strcmp(mode_ch, "q") == 0) {
        mode = false;  // so query genome
    } else {
        // wrong parameter, need either t or q
        fprintf(stderr, "Error! Either q or t parameter requested!\n");
        fprintf(stderr, "Assigned to T automatically\n");
        mode = true;
    }

    // parse the range | chrN:X-Y
    struct Regions request_region;
    char *reg_split = strtok(search_region, ":");
    strcpy(request_region.chrom, reg_split);
    reg_split = strtok(NULL, ":");
    char *split_range = strtok(reg_split, "-");
    request_region.start = strtol(split_range, NULL, 10);
    split_range = strtok(NULL, "-");
    request_region.end = strtol(split_range, NULL, 10);

    // state of the program
    // if reading -> we're going throw the region of interest
    bool head_caught = false;
    bool reading = false;

    // variables describing the block
    // variables called as pointers BUT
    // they are not pointers! (in C-meaning)
    int t_start_pointer = 0;
    int q_start_pointer = 0;
    int blockTStarts = 0; 
    int blockTEnds = 0;
    int blockQStarts = 0;
    int blockQEnds = 0;
    int req_block_starts = 0;
    int req_block_ends = 0;

    struct Chain_info in_chain_head;

    unsigned long long int rows_added = 0;
    char * curLine = chain;
    // answer is an array of strings (array)
    // we don't know the size of this array now
    // so will dynamically reallocate memory for this
    unsigned long long int arr_size = ALLOC_START;
    char ** answer = (char**)malloc(sizeof(char*) * arr_size);

    while (curLine)
    {
        // deal with reading line-by-line business
        // we are not reading a file, it's just a string
        char * nextLine = strchr(curLine, '\n');
        if (nextLine) 
        {
            *nextLine = '\0';
        }

        // parse head if it wasn't
        // it must happen to the first line
        if (!head_caught)
        {
            head_caught = true;  // should happen only once
            in_chain_head = parse_head(curLine);

            // mark starting position
            t_start_pointer = in_chain_head.tStart;
            q_start_pointer = in_chain_head.qStart;

            // invert coordinates if needed
            // if strand is negative, we need to do this operation
            // read chain file docs for more information
            if ((!in_chain_head.tStrand && mode) || (!in_chain_head.qStrand && !mode))
            {
                int temp = request_region.start;
                int req_chrom_size = mode ? in_chain_head.tSize : in_chain_head.qSize;
                request_region.start = req_chrom_size - request_region.end;
                request_region.end = req_chrom_size - temp;
            }

            // check chrom, it the requested region has a different chrom,
            // it is definitely wrong!
            // if we are interested in the reference region at chrom 1
            // and provided chain aligned to chrom 2 -> something is definitely wrong!
            if (mode && strcmp(request_region.chrom, in_chain_head.tName) != 0)
            {
                fprintf(stderr, "Error! Target genome chrom differs in chain and requested region!\n");
                break;
            }
            // check this also for the query
            else if (!mode && strcmp(request_region.chrom, in_chain_head.qName) != 0)
            {
                fprintf(stderr, "Error! Query genome chrom differs in chain and requested region!\n");
                break;
            }
        }  // end header parsing

        // we are reading a block, so update variables
        struct Block block = parse_block(curLine);
        // convert to absolute coordinates
        blockTStarts = t_start_pointer;
        blockQStarts = q_start_pointer;
        blockTEnds = blockTStarts + block.size;
        blockQEnds = blockQStarts + block.size;
        t_start_pointer = blockTEnds + block.dt;
        q_start_pointer = blockQEnds + block.dq;

        // then restore newline-char, just to be tidy   
        if (nextLine) *nextLine = '\n'; 
        curLine = nextLine ? (nextLine + 1) : NULL;

        // depending on mode we are interested in we look at
        // either the target or the query genome
        req_block_starts = mode ? blockTStarts : blockQStarts;
        req_block_ends = mode ? blockTEnds : blockQEnds;

        // these conditions are here to check whether anything of this is true:
        // 1) we haven't read the chain yet but just reached the region of interest
        // 2) we finished reading already and are outside the region of interest
        // conceptually, these conditions are mutually exclusive

        // if we're not reading yet but we just reached the region of interest
        // then start reading
        if (!reading && req_block_ends >= request_region.start)
        {
            // ---chainchainchainchainchain------
            // ----regionregion------------------
            //     ^ we are here
            reading = true;
            // for each row we allocate a string in the answer
            answer[rows_added] = (char*)malloc(MAXCHAR * sizeof(char));
            // and write this block to answer
            sprintf(answer[rows_added], "%d %d %d %d\n", blockTStarts, blockTEnds, blockQStarts, blockQEnds);
            rows_added++;
            // it is the first row --> a answer has place a priori
            continue;
        } 
        else if (req_block_starts >= request_region.end) 
        {
            // we have already read the region of interest
            // ---chainchainchainchainchain------
            // ----regionregion------------------
            //                     ^ we are here
            // we can simply stop reading chain then
            break;
        }

        if (reading)
        {
            // we are reading the chain now
            // need to add a block in the result
            // ---chainchainchainchainchain------
            // ----regionregion------------------
            //          ^ we are here
            answer[rows_added] = (char*)malloc(MAXCHAR * sizeof(char));
            sprintf(answer[rows_added], "%d %d %d %d\n", blockTStarts, blockTEnds, blockQStarts, blockQEnds);
            rows_added++;

            if (rows_added + MEM_CHUNK > arr_size)
            // make sure that we have allocated enough memory
            {
                // need to realloc answer array
                arr_size += ALLOC_EXT;
                answer = (char**)realloc(answer, sizeof(char*) * arr_size);
            }
        }

        // if not reading: just go to the next block
        // no need to write this condition
    }
    // finish answer with END
    // will be simpler to parse the output in the python func
    answer[rows_added] = (char*)malloc(MAXCHAR * sizeof(char));
    sprintf(answer[rows_added], "END");
    return answer;
}


int main()
{
    // should be compiled with -fPIC and -shared flags!
    printf("Warning! This code is not intended to be compiled as a standalone tool.\n");
    printf("Please compile this file with -fPIC and -shared flags.\n");
    printf("and create the /modules/extract_subchain_slib.so shared library file.\n");
    return 0;
}
