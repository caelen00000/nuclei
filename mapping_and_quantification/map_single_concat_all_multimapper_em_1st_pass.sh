#!/bin/bash

#mutant or wildtype
genotype=$1

# e.g. p1 or a0
stage=$2

# s1 or s2
sample=$3

out_prefix="concat_all_multimapper_em_1st_pass"

# Increase maximum number of open file descriptors. Avoids errors with --outSAMtype BAM SortedByCoordinate
# May have to increase further for parallel runs, or do stuff from https://github.com/alexdobin/STAR/issues/1292
ulimit -n 100000

#not enough room in /tmp/, rather than trying to resize it, I used the C drive for speed (theoretically).
#temp_dir="/tmp/temp_star/${stage}_${genotype}_${sample}"
temp_dir="/mnt/c/Users/camplain/Documents/STARtemp/${stage}_${genotype}_${sample}"

# Have to use a nonstandard tempdir due to running under WSL
rm -rf $temp_dir
mkdir -p $temp_dir
chmod 777 $temp_dir

mkdir $out_prefix
out_dir="${out_prefix}/${genotype}_${stage}_${sample}/"


if [[ "${stage}" == "wp" ]]; then
	in_dir="./fastqs/wp_lanes_concat"
else
	in_dir="./fastqs/small_medium_large_concat"
fi

STAR \
--runThreadN 36 \
--genomeDir genome_dir_large \
--outTmpDir $temp_dir/temp \
--genomeChrSetMitochondrial mitochondrion_genome \
--clipAdapterType CellRanger4 \
--soloUMIlen 12 \
--soloCBwhitelist barcode_whitelists/3Pv3.txt \
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
--soloType CB_UMI_Simple \
--soloFeatures Gene GeneFull GeneFull_Ex50pAS SJ Velocyto \
--soloCellFilter EmptyDrops_CR 10000 0.99 10 45000 90000 500 0.01 20000 0.01 10000 \
--soloMultiMappers EM \
--soloUMIdedup 1MM_CR \
--soloUMIfiltering MultiGeneUMI_CR \
--outReadsUnmapped Fastx \
--outSAMtype None \
--winAnchorMultimapNmax 100 \
--outFilterMultimapNmax 100 \
--outMultimapperOrder Random \
--limitBAMsortRAM 70000000000 \
--outFileNamePrefix $out_dir \
--readFilesIn "${in_dir}/${genotype}_${stage}_${sample}_r2.fastq" "${in_dir}/${genotype}_${stage}_${sample}_r1.fastq"

rm -rf $temp_dir
