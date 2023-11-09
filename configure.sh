#!/usr/bin/env bash
# exit_status=0
my_dir="${0%/*}"


OVERRIDE=false
BUILD_BED2GTF=false
CLEANUP=false

### Parsing arguments ###
for arg in "$@"
do
    case $arg in
        --override)
            OVERRIDE=true
            echo "Overriding existing CESAR installation and models"
            ;;
        --bed2gtf)
            BUILD_BED2GTF=true
            echo "Building bed2gtf"
            ;;
        --cleanup)
            CLEANUP=true
            echo "Cleanup activated"
            ;;
    esac
done

### Cleanup scenario ###
### Cleanup ###
if $CLEANUP; then
    printf "Cleaning up...\n"
    rm -f "${my_dir}"/modules/*.so
    rm -f "${my_dir}"/modules/chain_score_filter
    rm -f "${my_dir}"/modules/chain_filter_by_id
    rm -rf ./CESAR2.0
    rm -f ./models/*.dat
    echo "Cleanup completed"
fi

### C modules ###
# Check machine architecture and set appropriate flags
printf "Compiling C code...\n"

if [ "$(uname -m)" = "arm64" ]; then
    CFLAGS="-Wall -Wextra -O2 -g -std=c99 -arch arm64" # adjust flags for M1 if necessary
else
    CFLAGS="-Wall -Wextra -O2 -g -std=c99" # original flags for x86
fi

gcc $CFLAGS -o "${my_dir}"/modules/chain_score_filter "${my_dir}"/modules/chain_score_filter.c
gcc $CFLAGS -o "${my_dir}"/modules/chain_filter_by_id "${my_dir}"/modules/chain_filter_by_id.c
gcc $CFLAGS -fPIC -shared -o "${my_dir}"/modules/chain_coords_converter_slib.so "${my_dir}"/modules/chain_coords_converter_slib.c
gcc $CFLAGS -fPIC -shared -o "${my_dir}"/modules/extract_subchain_slib.so "${my_dir}"/modules/extract_subchain_slib.c
gcc $CFLAGS -fPIC -shared -o "${my_dir}"/modules/chain_bst_lib.so "${my_dir}"/modules/chain_bst_lib.c


### XGBoost models ###
if ! $OVERRIDE && { [[ -f "./models/se_model.dat" ]] || [[ -f "./models/me_model.dat" ]]; }
then
    printf "Model found\n";
else
    printf "XGBoost model not found\nTraining...\n"
    eval "python3 train_model.py"
    printf "Model created\n"
fi

### CESAR2.0 ###
if ! $OVERRIDE && [[ -f "./CESAR2.0/cesar" ]]
then
    printf "CESAR installation found\n"
else
    if [ -d ".git" ]
    then
        printf "Git repo detected, cloning CESAR using git submodule...\n"
        git submodule init CESAR2.0
        git submodule update CESAR2.0
        cd CESAR2.0 && make
    else
        printf "No git repo detected, downloading CESAR using wget...\n"
        wget -q https://github.com/kirilenkobm/CESAR2.0/archive/refs/tags/toga_submodule_1.zip
        unzip -q toga_submodule_1.zip -d .
        mv CESAR2.0-toga_submodule_1 CESAR2.0
        cd CESAR2.0 && make
    fi
    printf "Don't worry about '*** are the same file' message if you see it\n"
fi

### Bed2Gtf ###
if $BUILD_BED2GTF; then
    if ! $OVERRIDE && [[ -f "./bed2gtf/some_check_file" ]]
    then
        printf "bed2gtf installation found\n"
    else
        printf "bed2gtf installation not found, cloning and building\n"
        git submodule init bed2gtf
        git submodule update bed2gtf
        # ACTUAL BUILD COMMAND HERE
    fi
fi
