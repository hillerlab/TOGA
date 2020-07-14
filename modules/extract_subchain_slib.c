/*
For a chain and target range return blocks that
cover only the requested range; some sort of subchain in other words.
You can request subchain for the target or the query genome.
Usage: extract_subchain [chain] [q/t] [target/query genome range]

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
#define MEM_CHUNCK 5


char ** extract_subchain(char *chain, char *mode_ch, char *search_region)
{
    // read the chain file OR stdin
    // read genic range mode
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
    bool head_catched = false;
    bool reading = false;

    // variables describing the block
    // varialbes called as pointers BUT
    // they are not pointers!
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
        // it must happen at the first line
        if (!head_catched)
        {
            head_catched = true;  // should happen only once
            in_chain_head = parse_head(curLine);

            // mark starting position
            t_start_pointer = in_chain_head.tStart;
            q_start_pointer = in_chain_head.qStart;

            // invert coordinates if needed
            // if strand is negative, we need to do this operation
            // read chain file docs for more information
            if ((!in_chain_head.tStrand && mode) || (!in_chain_head.qStrand && !mode))
            {
                int temp = request_region.start;  // yes I know the way is stupid
                int req_chrom_size = mode ? in_chain_head.tSize : in_chain_head.qSize;
                request_region.start = req_chrom_size - request_region.end;
                request_region.end = req_chrom_size - temp;
            }

            // check chrom, it the requested region has a different chrom,
            // it is defenitely wrong!
            if (mode && strcmp(request_region.chrom, in_chain_head.tName) != 0)
            {
                fprintf(stderr, "Error! Target genome chrom differs in chain and requested region!\n");
                break;
            }
            // check for both target and query genome ranges
            else if (!mode && strcmp(request_region.chrom, in_chain_head.qName) != 0)
            {
                fprintf(stderr, "Error! Query genome chrom differs in chain and requested region!\n");
                break;
            }
        }  // end header parsing

        // we are reading a block, so update variables
        struct Block block = parse_block(curLine);
        blockTStarts = t_start_pointer;
        blockQStarts = q_start_pointer;
        blockTEnds = blockTStarts + block.size;
        blockQEnds = blockQStarts + block.size;
        t_start_pointer = blockTEnds + block.dt;
        q_start_pointer = blockQEnds + block.dq;

        // then restore newline-char, just to be tidy   
        if (nextLine) *nextLine = '\n'; 
        curLine = nextLine ? (nextLine + 1) : NULL;

        // depending on mode we are interested in eather the target
        // or the query genome
        req_block_starts = mode ? blockTStarts : blockQStarts;
        req_block_ends = mode ? blockTEnds : blockQEnds;

        // if we're not reading but we just reached the region of interest, then
        // start reading
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
            //                 ^ we are here
            break;
        }

        if (reading)
        {
            // just continue reading the chain, now we are in the overlap
            // we have already read the region of interest
            // ---chainchainchainchainchain------
            // ----regionregion------------------
            //          ^ we are here
            answer[rows_added] = (char*)malloc(MAXCHAR * sizeof(char));
            sprintf(answer[rows_added], "%d %d %d %d\n", blockTStarts, blockTEnds, blockQStarts, blockQEnds);
            rows_added++;

            if (rows_added + MEM_CHUNCK > arr_size)
            // make sure that we have allocated enough memory
            {
                // need to realloc answer array
                arr_size += ALLOC_EXT;
                answer = (char**)realloc(answer, sizeof(char*) * arr_size);
            }
        }
    }
    // finish answer with END
    // will be simpler to parse the output
    answer[rows_added] = (char*)malloc(MAXCHAR * sizeof(char));
    sprintf(answer[rows_added], "END");
    return answer;
}


int main()
{
    // should be compiled with -fPIC and -shared flags!
    printf("Warining! This code is not intended to be compiled as a standalone tool.\n");
    printf("Please compile this file with -fPIC and -shared flags.\n");
    printf("and create the /modules/extract_subchain_slib.so shared library file.\n");
    return 0;
}
