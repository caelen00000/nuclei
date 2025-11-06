#!/bin/bash

source activate cellbender-test

#mutant, wildtype, or both
genotype=$1

#auto or expect_10k_tuned, default auto
run_type=${2:-"auto"}

#raw for unique counts, raw_multi_floor for unique+multimappers, default raw_multi_floor
#don't use raw_multi, CellBender doesnt't like non-integer counts
quant_dir=${3:-"raw_multi_floor"}

if [[ "${genotype}" == "wildtype" || "${genotype}" == "both" ]]; then
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh wildtype a0 s1 ${quant_dir}
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh wildtype a0 s2 ${quant_dir}
	
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh wildtype p1 s1 ${quant_dir}
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh wildtype p1 s2 ${quant_dir}
	
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh wildtype p2 s1 ${quant_dir}
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh wildtype p2 s2 ${quant_dir}
	
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh wildtype p3 s1 ${quant_dir}
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh wildtype p3 s2 ${quant_dir}
	
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh wildtype p4 s1 ${quant_dir}
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh wildtype p4 s2 ${quant_dir}
	
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh wildtype wp s1 ${quant_dir}
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh wildtype wp s2 ${quant_dir}
fi

if [[ "${genotype}" == "mutant" || "${genotype}" == "both" ]]; then
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh mutant a0 s1 ${quant_dir}
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh mutant a0 s2 ${quant_dir}
	
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh mutant p1 s1 ${quant_dir}
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh mutant p1 s2 ${quant_dir}
	
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh mutant p2 s1 ${quant_dir}
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh mutant p2 s2 ${quant_dir}
	
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh mutant p3 s1 ${quant_dir}
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh mutant p3 s2 ${quant_dir}
	
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh mutant p4 s1 ${quant_dir}
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh mutant p4 s2 ${quant_dir}
	
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh mutant wp s1 ${quant_dir}
	../empty_droplet_calling_and_background_removal/run_cellbender_${run_type}.sh mutant wp s2 ${quant_dir}
fi

