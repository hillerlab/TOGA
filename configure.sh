#!/usr/bin/env bash
exit_status=0
mydir="${0%/*}"

printf "Compiling C code...\n"
CFLAGS="-Wall -Wextra -O2 -g -std=c99"

gcc $CFLAGS -o ${mydir}/modules/chain_score_filter ${mydir}/modules/chain_score_filter.c 
gcc $CFLAGS -o ${mydir}/modules/chain_filter_by_id ${mydir}/modules/chain_filter_by_id.c 
gcc $CFLAGS -fPIC -shared -o ${mydir}/modules/chain_coords_converter_slib.so ${mydir}/modules/chain_coords_converter_slib.c
gcc $CFLAGS -fPIC -shared -o ${mydir}/modules/extract_subchain_slib.so ${mydir}/modules/extract_subchain_slib.c
gcc $CFLAGS -fPIC -shared -o ${mydir}/modules/chain_bst_lib.so ${mydir}/modules/chain_bst_lib.c
exit ${exit_status}

if [[ -f "./cesar" ]]
then
    printf "CESAR installation found\n"
else
    printf "CESAR installation not found, cloning\n"
    git clone https://github.com/hillerlab/CESAR2.0/
    cd CESAR2.0/
    make
    cd ..
    echo $'#!/usr/bin/env bash' > cesar
    echo $'mydir="${0%/*}"' >> cesar
    echo $'exe="$mydir/CESAR2.0/cesar"' >> cesar
    echo $'$exe "$@"' >> cesar
    chmod +x ./cesar
    printf "CESAR installed\n"
fi


if [[ -f "./models/se_model.dat" ]] || [[ -f "./models/me_model.dat" ]]
then
    printf "Model found\n";
else
    printf "XGBoost model not found\nTraining...\n"
    eval "python3 train_model.py"
    printf "Model created\n"
fi
