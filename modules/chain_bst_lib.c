/*
To be compiled as a shared library
make_index and  get_s_byte functions

Store binary search tree (index) for a chain_file.
In detail BST contains start byte and offset for each
chain_id. Which means that we can rapidly extract
any chain from a given chain file.

Author: Bogdan Kirilenko, 2020;
*/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define MIN_TREE_SIZE 8

typedef struct{
// this is a BST node
// key: chain_id
// start byte and offset: our value
    uint64_t chain_id;
    uint64_t start_byte;
    uint64_t offset;
    bool is_null;
    bool is_terminal;
} Point;

// initiate empty node
// Point NOTHING = {};


int compare_structs
// compare chain_ids of two points
(const void *a, const void *b)
{
    Point *ia = (Point*)a;
    Point *ib = (Point*)b;
    if (ia->chain_id > ib->chain_id){
        return 1;
    } else {
        return -1;
    }
}   


// void print_point
// for debug: print point data
// (Point *point)
// {
//     printf("Chain ID %llu;", point->chain_id);
//     printf("Start byte %llu;", point->start_byte);
//     printf("Offset %llu\n", point->offset);
// }


void copy_point
// copy point structure (except is_null and is_terminal)
(Point *from, Point *to)
{
    to->chain_id = from->chain_id;
    to->start_byte = from->start_byte;
    to->offset = from->offset;
    to->is_null = false;
    to->is_terminal = false;
}


void add_node
// add node to BSD (recursive function!)
(Point *points, Point *bst, uint64_t bst_size, 
uint64_t low, uint64_t high, uint64_t i)
{
    // printf("Adding node %llu\n", i);
    if (i >= bst_size)
    // we are out of borders: skip
    {
        // printf("Node out of borders: return\n");
        return;
    }

    if (low > high)
    // this branch doesn't actually exist (nothing is here)
    {
        // printf("Non-existent branch\n");
        bst[i].is_null = true;
        bst[i].is_terminal = true;
        return;
    }

    uint64_t mid = low + (high - low) / 2;
    // printf("Lo: %llu; hi: %llu; mid: %llu\n", low, high, mid);
    copy_point(&points[mid], &bst[i]);

    if (2 * i > bst_size)
    // this point has no children
    {
        // printf("Terminal branch (second half of array)\n");
        bst[i].is_terminal = true;
        return;
    }
    // add left branch
    // to avoid unsigned overflow
    uint64_t left_mid;
    if (mid > 0)
    // we can decrement safely
    {
        left_mid = mid - 1;
    } else
    // we cannot increment -> unsigned overflow
    {
        left_mid = 0;
    }
    add_node(points, bst, bst_size, low, left_mid, 2 * i);
    // add right branch
    add_node(points, bst, bst_size, mid + 1, high, (2 * i) + 1);
    return;
}


Point *make_bst
// function to make bst array itself
(Point *array, uint64_t arr_size, uint64_t *bst_size)
{
    *bst_size = 0;
    uint64_t i = 0;

    // define the size of BST array
    // must be a number == 2 ** N and > initial array size
    while (*bst_size < arr_size)
    {
        *bst_size = pow(2, i);
        ++i;
    }
    *bst_size = *bst_size > MIN_TREE_SIZE ? *bst_size : MIN_TREE_SIZE;
    *bst_size += 1;

    // printf("BST size is %llu\n", *bst_size);
    Point *bst = (Point*)malloc(sizeof(Point) * *bst_size);
    // fill bst with zeros
    for (uint64_t i = 0; i < *bst_size; ++i)
    {
        bst[i].chain_id = 0;
        bst[i].start_byte = 0;
        bst[i].offset = 0;
        bst[i].is_null = 0;
        bst[i].is_terminal = 0;
    }
    uint64_t init_i = 1;

    // recursively fill the BST array:
    add_node(array, bst, *bst_size, 0, arr_size - 1, init_i);
    return bst;
}


int make_index
// shared library entry point 1
// loads all data into an index file
// python provides 3 arrays:
// chain identifiers, starting positions and offsets
// + length of these arrays (they must be same)
// this function packs them into an array of Point structures
// and saves into the BST table
(uint64_t *chain_ids, uint64_t *start_bytes, uint64_t *offsets,
uint64_t arr_size, char *table_path)
{
    // struct of arrays to array of structs
    Point *array = (Point *)malloc(sizeof(Point) * arr_size);

    // load arrays data to structs
    for (uint64_t i = 0; i < arr_size; ++i)
    {
        array[i].chain_id = chain_ids[i];
        array[i].start_byte = start_bytes[i];
        array[i].offset = offsets[i];
        array[i].is_null = false;
        array[i].is_terminal = false;
    }

    // we need chains to be sorted (search is binary)
    qsort(array, arr_size, sizeof(Point), compare_structs);

    // organize these points into a BST
    uint64_t bst_size = 0;
    Point *bst = make_bst(array, arr_size, &bst_size);

    // initiate root point
    bst[0].chain_id = 0;
    bst[0].start_byte = 0;
    bst[0].offset = 0;
    bst[0].is_terminal = false;
    bst[0].is_null = false;

    // and now we can save this array
    // for (uint64_t i = 0; i < bst_size; ++i)
    // {
    //     print_point(&bst[i]);
    // }
    FILE *fp = NULL;
    fp = fopen(table_path, "wb");
    fwrite(bst, sizeof(Point), bst_size, fp);
    fclose(fp);
    // free(bst);
    // free(array);
    return 0;
}


uint64_t get_s_byte
// find value for the chain_id given
// value is: start_byte + offset
// they are provided by reference
(char *index_path, uint64_t chain_id,
uint64_t *start_byte, uint64_t *offset)
{
    FILE *fp = NULL;  // initiate file reading
    fp = fopen(index_path, "rb");
    Point reader;

    uint64_t cur = 0;  // position of the next structure
    uint64_t n = 0;  // number of structure in the file reader

    while (fread(&reader, sizeof(Point), 1, fp))
    // file is basically a sequence of Point structures
    // read them one-by-one
    {
        if (cur != n){
            // cur -> where are we going
            // n -> where are we now (number of struct in the file we read)
            // if cur != n: skip this chunk, goto the next one
            ++n;
            continue;
        }

        // basically a binary search procedure
        if (chain_id > reader.chain_id)
        // requested chain_id is greater than current
        // goto the right branch
        {
            cur = (2 * n + 1);
        } else if (chain_id < reader.chain_id)
        // go to the left branch
        {
            cur = (2 * n);
        } else
        // chain_id == current_chain_id: got it
        // write values and break
        {
            *start_byte = reader.start_byte;
            *offset = reader.offset;
            break;
        }
        ++n;  // increment current struct number
    }
    fclose(fp);
    // if we found nothing: offset and start byte are not changed
    // upstream function will catch this
    return 0;
}
