#!/usr/bin/env bash
# orthology_type_map.py wrapper
# calls the script on a toga directory
set -e

if [ "$#" -ne 1 ];
then
    echo "Usage: $0 [TOGA OUTPUT DIRECTORY]"
    exit 0
fi

ref_bed="toga_filt_ref_annot.bed"
que_bed="query_annotation.bed"
out="orthology_classification.tsv"
riforms="isoforms.tsv"
qiforms="query_isoforms.tsv"
paralogs="paralogs.txt"
loss_data="loss_summ_data.tsv"
o_scores="orthology_scores.tsv"

cur_dir=$(dirname "$(readlink -f "$0")")
omap_script="${cur_dir}/orthology_type_map.py"
cmd="${omap_script} ${1}/${ref_bed} ${1}/${que_bed} ${1}/${out} --ri ${1}/${riforms} --qi ${1}/${qiforms} \
    -p ${1}/${paralogs} -l ${1}/${loss_data} -o ${1}/${o_scores}"

echo "Calling ${cmd}"
eval "${cmd}"

echo "DONE"
