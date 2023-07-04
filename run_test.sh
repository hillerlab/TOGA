#!/usr/bin/env bash
set -e

if [ "$#" -lt 1 ] || ( [ "$1" != "micro" ] && [ "$1" != "normal" ] ); then
    echo "Usage: $0 [MODE: micro/normal] [if normal: paths to human and mouse 2 bit files]"
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

if [ $1 = "micro" ];
then
    echo "Running test in the micro mode"
    eval "rm -rf micro_test_out"
    cmd="./toga.py test_input/align_micro_sample.chain test_input/annot_micro_sample.bed \
         test_input/hg38.micro_sample.2bit test_input/q2bit_micro_sample.2bit \
         --pn micro_test_out  --kt \
         --cjn 1 --chn 1 --ms \
         "
    echo "Calling ${cmd}"
    eval $cmd
    eval "cat micro_test_out/query_annotation.bed"

elif [ $1 = "normal" ];
then
    echo "Running widescale test"
    eval "rm -rf test_out"

    hg38_2bit=$2
    mm10_2bit=$3
    if [ "$#" -lt 3 ]; then
        echo "For widescale test please define paths to hg38 and mm10 2bit files"
        echo "Such as: $0 normal hg38.2bit mm10.2bit"
        exit 0
    fi
    test_bed="test_input/hg38.genCode27.chr11.bed"
    test_chain="test_input/hg38.mm10.chr11.chain"

    cmd="./toga.py  ${test_chain} ${test_bed} \
         ${hg38_2bit} ${mm10_2bit} \
         --chn 10 --project_name test_out \
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
