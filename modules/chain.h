/* Header file containing chain-related structures and functions
*/
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#define MAXCHAR 255

// all the necessary data from chain header
struct Chain_info{
    char tName[MAXCHAR];
    int tSize;
    bool tStrand;
    int tStart;
    int tEnd;
    char qName[MAXCHAR];
    int qSize;
    bool qStrand;
    int qStart;
    int qEnd;
};

// also chain block data struct
struct Block{
    int size;
    int dt;
    int dq;
};

// struct describing genomic region CHROM START END
struct Regions
{
    char chrom[MAXCHAR];
    int start;
    int end;
};

// parse chain header and return all values from it
// read chain file format explanation for extra info
struct Chain_info parse_head(char *header_string)
{
    struct Chain_info head;  // init empty struct
    // space-separated fields:
    char *split_head = strtok(header_string, " ");
    short int field = 0;
    // field 0 - just word "chain"
    // field 1 - chain score -> ignore this
    while(split_head) {
        switch (field)
        {
            case 2:  // extract tName
                strcpy(head.tName, split_head);
                break;
            case 3:
                head.tSize = strtol(split_head, NULL, 10);
                break;
            case 4:
                if (strcmp(split_head, "+") == 0) {head.tStrand = true;}
                else {head.tStrand = false;}
                break;
            case 5:
                head.tStart = strtol(split_head, NULL, 10);
                break;
            case 6:
                head.tEnd = strtol(split_head, NULL, 10);
                break;
            case 7:
                strcpy(head.qName, split_head);
                break;
            case 8:
                head.qSize = strtol(split_head, NULL, 10);
                break;
            case 9:
                if (strcmp(split_head, "+") == 0) {head.qStrand = true;}
                else {head.qStrand = false;}
                break;
            case 10:
                head.qStart = strtol(split_head, NULL, 10);
                break;
            case 11:
                head.qEnd = strtol(split_head, NULL, 10); 
                break;
        }
        split_head = strtok(NULL, " ");
        field += 1;
    }
    return head;
}

// parse chain block
// pls read chain docs for more information
struct Block parse_block(char *block_string)
{
    struct Block block;
    // default values
    block.size = 0;
    block.dq = 0;
    block.dt = 0;

    char *split_block = strtok(block_string, "\t");
    short int field = 0;
    while (split_block) {
        switch (field) {
            case 0: block.size = strtol(split_block, NULL, 10); break;
            case 1: block.dt = strtol(split_block, NULL, 10); break;
            case 2: block.dq = strtol(split_block, NULL, 10); break;
        }
        split_block = strtok(NULL, "\t");
        field += 1;
    }
    return block;
}
