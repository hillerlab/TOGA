#!/usr/bin/env bash
set -e

if [ "$#" -ne 1 ];
then
    echo "Usage: $0 [MODE: micro/normal]"
    exit 0
fi

eval ./configure.sh
if [[ $? -ne 0 ]]
then
    echo "Configure + build failed"
    exit 1
else
    echo "Configure successful!"
fi

test_out="test_out"
eval rm -rf ${test_out}
eval mkdir -p ${test_out}

if [ $1 = "micro" ];
then
    echo "Running test in the micro mode"
    eval "rm -rf micro_test_out"
    cmd="./toga.py test_input/hg38.mm10.chrM.chain test_input/hg38.chrM.bed \
         test_input/hg38.chrM.2bit test_input/mm10.chrM.2bit \
         --pn micro_test_out  --kt \
         --cjn 1 --chn 1 --ms \
         --no_para
         "
    echo "Calling ${cmd}"
    eval $cmd
    eval "cat micro_test_out/query_annotation.bed"

elif [ $1 = "normal" ];
then
    echo "Running widescale test"
    echo "Run on cluster only!"
    test_bed="test_input/hg38.genCode27.chr11.bed"
    test_chain="test_input/hg38.mm10.chr11.chain"

    cmd="./toga.py  ${test_chain} ${test_bed} \
         hg38 mm10 \
         --chn 10 --project_name ${test_out} \
         --kt -i supply/hg38.wgEncodeGencodeCompV34.isoforms.txt \
         --cjn 200 --ms --cb 2,4 \
         "
    echo "Calling ${cmd}"
    eval $cmd
else
    echo "Please choose between 'micro' and 'normal' modes"
    echo "Usage: ./run_test.sh [micro/normal]"
    exit 0
fi

if [[ $? -ne 0 ]]
then
    echo "TOGA failed"
    exit 1
else
    echo "Success!"
fi
