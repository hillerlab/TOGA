#!/usr/bin/env bash
# Call a sequence of scripts to continue TOGA process
# after the CESAR step
set -e

if [ "$#" -ne 1 ];
then
    echo "Usage: $0 [OUTPUT DIRECTORY]"
    exit 0
fi

cesar_res_dir="cesar_results/"
intermediate_bed="intermediate.bed"
nucl_fasta="nucleotide.fasta"
exons_meta="exons_meta_data.tsv"
prot_fasta="prot.fasta"
codon_fasta="codon.fasta"
inact_data_dir="inact_mut_data/"
ref_bed="toga_filt_ref_annot.bed"
query_annot="query_annotation.bed"
loss_summ="loss_summ_data.tsv"
isoforms="isoforms.tsv"
paralogs="paralogs.txt"
query_isoforms="query_isoforms.tsv"
orth_class="orthology_classification.tsv"


cmd="./modules/merge_cesar_output.py ${1}/${cesar_res_dir} ${1}/${intermediate_bed} \
     ${1}/${nucl_fasta} ${1}/${exons_meta} ${1}/${prot_fasta} ${1}/${codon_fasta} \
     ${1}/rejected/cesar_step_rejected.txt"
echo "Merging cesar output"
echo "Calling ${cmd}"
eval ${cmd}

cmd="./modules/gene_losses_summary.py ${1}/${inact_data_dir} ${1}/${ref_bed} \
     ${1}/${intermediate_bed} ${1}/${query_annot} ${1}/${loss_summ} \
     -i ${1}/${isoforms} --paral_projections ${1}/${paralogs}"
echo "Gene loss part..."
echo "Calling ${cmd}"
eval ${cmd}

echo "Making query isoforms data"
cmd="./modules/make_query_isoforms.py ${1}/${query_annot} ${1}/${query_isoforms} \
    --genes_track ${1}/query_gene_spans.bed"
echo "Calling ${cmd}"
eval ${cmd}

echo "Orthology classification"
cmd="./modules/orthology_type_map.py ${1}/${ref_bed} ${1}/${query_annot} ${1}/${orth_class} \
    --ri ${1}/${isoforms} --qi ${1}/${query_isoforms} -p ${1}/${paralogs}  \
    -l ${1}/${loss_summ} -s ${1}/rejected/ref_orphan_transcripts.txt"
echo "Calling ${cmd}"
eval ${cmd}

cat ${1}/${cesar_res_dir}/* > ${1}/cesar_results.txt
cat ${1}/${inact_data_dir}/* > ${1}/inact_mut_data.txt

echo "DONE"
