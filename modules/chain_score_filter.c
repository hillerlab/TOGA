/*
Filter chain file: extract only those that have a score > min_score
Write output to stdout
Usage: chain_score_filter [chain_file] [minimal_score]

Author: Bogdan Kirilenko, 2020;
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <stdbool.h>

#define MAXCHAR 255
#define CHAIN_WORD_LEN 5
#define ARGS_NUM 3
#define NUM_BASE 10


uint64_t get_chain_score(char *header_string)
{
    // parse chain header, we need a field num 1.
    // fields are separated with space
    char *split_head = strtok(header_string, " ");
    // initiate vars with default values
    uint64_t score = 0;
    uint8_t field = 0;

    while(split_head)
    {
        switch (field)
        {
            case 0: break;
            case 1:  // this is the score
                score = strtoul(split_head, NULL, NUM_BASE);
                break;
            default: break;
        }
        split_head = strtok(NULL, " ");
        field += 1;
    }
    return score;
}


int main(int argc, char **argv)
{
    // what if the number of arguments is wrong?
    // then show usage and quit
    if (argc != ARGS_NUM)
    {
        fprintf(stderr, "Usage: %s [chain_file] [min_score]\n", argv[0]);
        return 0;
    }

    // check that the second arg is numeric
    if (strspn(argv[2], "0123456789") != strlen(argv[2]))
    {
        fprintf(stderr, "Error! Second arg must be numeric, got %s\n", argv[2]);
        return 1;
    }

    uint32_t min_score = strtoul(argv[2], NULL, 10);
    // read the chain file OR stdin
    FILE *chain;

    if (strcmp(argv[1], "stdin") == 0) {
        chain = stdin;  // if stdin --> ok
    } else if (access(argv[1], F_OK) != -1) {
        // not stdin; check that file exists
        chain = fopen(argv[1], "r");
    } else {
        // not stderr, non-existent file -> kill the process
        fprintf(stderr, "Error! File %s doesn't exist!\n", argv[1]);
        return 1;
    }

    // we will read file line-by-line
    // some lines are headers, the rest -> are blocks
    // we read a chain header and extract the score
    // if score is > min_score, we print the header and all blocks
    // beneath this header, until we get to the next one
    // then we repeat the first step again

    char buff[MAXCHAR];  // file reading buffer
    bool reading = false;  // if true: we are reading a chain with
                           // score > required minimal
    uint64_t score = 0;  // current chain score, initiate with 0
    uint64_t passed = 0;  // count chains that passed the filter

    while (fgets(buff, MAXCHAR, chain) != NULL)
    {
        if (strncmp("chain", buff, CHAIN_WORD_LEN) == 0)
        {
            // it is a header, starts with "chain"
            // need to extract chain score
            char header[MAXCHAR];
            strcpy(header, buff);
            score = get_chain_score(header);
            
            // compare score to the min score
            if (score > min_score){
                // if this happened -> print the header
                // as a part of the chain
                printf("%s", buff);
                // reading true -> will print all blocks
                // belonging to this chain
                reading = true;
                // also count chains that passed the filter
                ++passed;
            // if score < required then
            // reading flag is false
            // skip all lines until we reach the next header
            } else {reading = false;}

        } else if (reading) {
            // print if reading -> score of
            // this chain was > required
            printf("%s", buff);
        }
        // else: no reading -> do nothing
    }  // stop reading file, we can close it now
    fclose(chain);

    if (passed == 0)
    {
        // if we counted 0 chains that passed the filter, the user should be warned at least
        fprintf(stderr, "Error! No chains passed the filter, probably the input is empty!\n");
        return 1;
    }
    return 0;
}
