/*
For a chain file and minimal score writes to stdout
only those chains that have score more than required.
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


uint64_t chain_score(char *header_string)
{
    // parse chain header, we need a field num 1.
    // fields are separated with space
    char *split_head = strtok(header_string, " ");
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
    if (argc != ARGS_NUM)
    {
        fprintf(stderr, "Usage: %s [chain_file] [min_score]\n", argv[0]);
        return 0;
    }

    // check that second arg is a number
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

    char buff[MAXCHAR];  // file reading buffer
    bool reading = false;  // if true: we are reading a chain with
                           // score > required minimal
    uint64_t score = 0;
    uint64_t passed = 0;

    while (fgets(buff, MAXCHAR, chain) != NULL)
    {
        if (strncmp("chain", buff, CHAIN_WORD_LEN) == 0)
        {
            // it is a header
            char header[MAXCHAR];
            strcpy(header, buff);
            score = chain_score(header);
            
            // if this happened -> print the header
            // as a part of the chain
            if (score > min_score){
                printf("%s", buff);
                // reading true -> will print all blocks
                // belonging to this chain
                reading = true;
                // also count chains that passed the filter
                ++passed;
            // if score < required then
            // reading flag is false
            } else {reading = false;}

        } else if (reading) {
            // print if reading -> score of
            // this chain was > required
            printf("%s", buff);
        }
        // else: no reading -> do nothing
    }
    fclose(chain);
    if (passed == 0)
    {
        // if we counted 0 chains that passed the filter, the user should be warned at least
        fprintf(stderr, "Warning! No chains passed the filter, probably the input is empty!\n");
    }

    return 0;
}
