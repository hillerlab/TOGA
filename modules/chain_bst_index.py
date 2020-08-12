#!/usr/bin/env python3
"""Save chain index file.

Chain index is just a BST saved into a separate file.
Using this index file for any chain ID we can extract:
1) Start byte of this chain in the chain file
2) Length of this chain in the file.
And then simply extract it.
"""
import sys
import os
import ctypes

__author__ = "Bogdan Kirilenko, 2020."
__version__ = "1.0"
__email__ = "kirilenk@mpi-cbg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]

SLIB_NAME = "chain_bst_lib.so"


def chain_bst_index(chain_file, index_file, txt_index=None):
    """Create index file for chain."""
    # assume that shared lib is in the same dir
    script_location = os.path.dirname(__file__)
    slib_location = os.path.join(script_location, SLIB_NAME)
    if not os.path.isfile(slib_location):
        sys.exit("chain_bst_lib.so lib not found, call ./configure to compile")
    # connect shared lib
    sh_lib = ctypes.CDLL(slib_location)
    sh_lib.make_index.argtypes = [ctypes.POINTER(ctypes.c_uint64),
                                  ctypes.POINTER(ctypes.c_uint64),
                                  ctypes.POINTER(ctypes.c_uint64),
                                  ctypes.c_uint64,
                                  ctypes.c_char_p]
    sh_lib.make_index.restype = ctypes.c_int

    # read chain file, get necessary data: start bytes and offsets
    print("Reading chain...")
    chain_ids, start_bytes, offsets = [0, ], [0, ], [0, ]
    byte_num, offset = 0, 0

    f = open(chain_file, "rb")
    for line in f:
        if not line.startswith(b"chain"):
            # just count these bytes
            byte_num += len(line)
            offset += len(line)
            continue
        # if we're here -> this is a chain header
        offsets.append(offset)
        header = line.decode("utf-8").rstrip()
        chain_id = header.split()[-1]
        chain_ids.append(int(chain_id))
        start_bytes.append(int(byte_num))
        # byte num is absolute and offset (chain size) is relative
        byte_num += len(line)  # si just continue incrementing byte num
        offset = len(line)  # and reset (to line length) offset
    f.close()

    offsets.append(offset)
    arr_size = len(chain_ids)
    del offsets[0]

    if arr_size == 0:
        sys.exit("Error! No chains found. Abort")
    print(f"There are {arr_size} chains to index")

    if txt_index:
        # save text (non-binary) dict for chain ids and positions in the file
        # in some cases this is more efficient that extracting data from BST
        # via shared library
        f = open(txt_index, "w")
        for x in zip(chain_ids, start_bytes, offsets):
            f.write(f"{x[0]}\t{x[1]}\t{x[2]}\n")
        f.close()

    # call shared lib
    print("Saved chain index data...")
    c_chain_ids = (ctypes.c_uint64 * (arr_size + 1))()
    c_chain_ids[:-1] = chain_ids
    c_s_bytes = (ctypes.c_uint64 * (arr_size + 1))()
    c_s_bytes[:-1] = start_bytes
    c_offsets = (ctypes.c_uint64 * (arr_size + 1))()
    c_offsets[:-1] = offsets

    c_arr_size = ctypes.c_uint64(arr_size)
    c_table_path = ctypes.c_char_p(str(index_file).encode())
    _ = sh_lib.make_index(c_chain_ids,
                          c_s_bytes,
                          c_offsets,
                          c_arr_size,
                          c_table_path)
    print(f"Saved chain {chain_file} index to {index_file}")


if __name__ == "__main__":
    try:  # read args
        in_chain_arg = sys.argv[1]
        index_file_arg = sys.argv[2]
    except IndexError:  # if args are wrong: show usage and quit
        sys.stderr.write("Recommended file extension for index file is bst.\n")
        sys.exit(f"Usage: {sys.argv[0]} [in_chain] [index_file]\n")
    chain_bst_index(in_chain_arg, index_file_arg)
