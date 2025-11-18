#!/bin/bash

source activate cellbender-test

#mutant or wildtype
genotype=$1

# e.g. p1 or a0
stage=$2

# s1 or s2
sample=$3

#raw for unique counts, raw_multi_floor for unique+multimappers, default raw_multi_floor
#don't use raw_multi, CellBender doesnt't like non-integer counts
quant_dir=${4:-"raw_multi_floor"}

seq_dir="/mnt/z/Caelen/snRNAseq_v2"

seq_run="concat_all_multimapper_em_1st_pass"

ID="${genotype}_${stage}_${sample}"

out_dir="cellbender/${seq_run}_auto_tuned/${ID}"
mkdir -p $out_dir

cellbender remove-background \
--cuda \
--input "${seq_dir}/${seq_run}/${ID}/Solo.out/GeneFull_Ex50pAS/${quant_dir}" \
--output "${out_dir}/cellbender.h5" \
--checkpoint-mins 1000000 \
--estimator-multiple-cpu \
--total-droplets-included 25000 \

rm ckpt.tar.gz