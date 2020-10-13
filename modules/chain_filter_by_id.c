/*
Extract a single chain of ID requested from a chain file.
Based on assumption, that each chain has an unique ID.
Also chain_ids should start with 1.

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


// parse chain header to get an ID value
uint64_t get_chain_id(char *header_string)
{
    // header is space-separated
    char *split_head = strtok(header_string, " ");
    // initiate values such as field and chain id
    uint64_t chain_id = 0;
    uint8_t field = 0;

    while(split_head)
    {
        switch (field)
        {
            case 0: break;
            case 12:  // this field contains the score, num 12
                chain_id = strtoul(split_head, NULL, NUM_BASE); break;
            default: break;  // other fields are out of our interest
        }
        split_head = strtok(NULL, " ");
        field += 1;
    }
    return chain_id;
}


int main(int argc, char **argv)
{
    if (argc != ARGS_NUM)
    {
        // just show usage
        fprintf(stderr, "Usage: %s [chain_file] [chain_id]\n", argv[0]);
        return 0;
    }

    // check that second arg is a number
    // maybe not the most elegant way
    // negative chain scores are not allowed
    if (strspn(argv[2], "0123456789") != strlen(argv[2]))
    {
        fprintf(stderr, "Error! Second arg must be numeric, got %s\n", argv[2]);
        return 1;  // clearly an error
    }

    // and read the second arg then
    uint64_t chain_id_req = strtol(argv[2], NULL, NUM_BASE);

    // read the chain file OR stdin
    FILE *chain;
    if (strcmp(argv[1], "stdin") == 0) {
        chain = stdin;  // if stdin --> ok
    } else if (access(argv[1], F_OK) != -1) {
        // not stdin; check that file exists
        chain = fopen(argv[1], "r");
    } else {
        // OR there is no file, raise an error then
        fprintf(stderr, "Error! File %s doesn't exist!\n", argv[1]);
        return 1;
    }

    // file reading buffer, 255 characters must be enough
    char buff[MAXCHAR];
    uint64_t chain_id = 0;
    // if reading == TRUE -> print what in the buffer
    // it must happen only if we reading a chain
    // of the requested ID
    // init with false: we haven't read any chain yet
    bool reading = false;

    // read chain file line-by-line
    while (fgets(buff, MAXCHAR, chain) != NULL)
    {
        if (strncmp("\n", buff, 1) == 0 && reading)
        {
            // a single \n line means that chain is over
            // if reading -> that we're reading the chain we need
            // in sum, we should stop here, because the chain 
            // we are interested in is over
            printf("\n");
            break;
        }
        else if (strncmp("chain", buff, CHAIN_WORD_LEN) == 0)
        {
            // a new chain header, need to extract the ID and compare with requested
            char header[MAXCHAR];
            strcpy(header, buff);
            chain_id = get_chain_id(header);
            if (chain_id_req == chain_id)
            {
                // if our chain -> need to print the header
                // and turn reading (printing the chain) on
                printf("%s", buff);
                reading = true;
            // else/ if it's not our chain -> just do nothing
            } else {reading = false;}

        }
        else if (reading)
        {
            // if reading -> this is a part of our chain
            printf("%s", buff);
        }

        // if not reading: skip the line, goto the next one
        // no need to write this condition
    }  // stop reading file, we can close it now
    fclose(chain);

    if (!reading)
    {
        // if we are here -> we gone through the entire file, but didn't catch the chain we need
        // meaning there was no chain with the requested ID
        // most likely this is an error: return code 1 just in case
        fprintf(stderr, "Warning! No chains passed the filter, probably the input is empty!\n");
        return 1;
    }
    // result is in the stdout now:
    // just return 0
    return 0;
}
