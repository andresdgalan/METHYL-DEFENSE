#!/bin/bash
set -e  # stop if any script fails

# Master script to run VARIMETH pipeline in order
# To run this pipeline we need a REFERENCE genome and FILTERED, DEMULTIPLEXED SAMPLES
# In DEMULTIQUAL pipeline we obtained SAMPLES and a MultiQC quality reports which we need to evaluate before starting
# In DENOVO pipeline we obtained a de novo constructed REFERENCE, but we can also use available references in ncbi

# Default start is VARIMETH1, if you want to start at any other point run something like this:
# ./VARIMETH0_masterscript.sh VARIMETH3

start_at=${1:-"VARIMETH1"}

run_step() {
    step_name=$1
    shift
    if [[ "$start_at" == "$step_name" || "$started" == "yes" ]]; then
        started="yes"
        echo "Running $step_name..."
        "$@"
        echo "Finished $step_name"
    else
        echo "Skipping $step_name"
    fi
}

run_step VARIMETH1 ./VARIMETH1_bismark.sh
run_step VARIMETH2 ./VARIMETH2_MDtag.sh
run_step VARIMETH3 ./VARIMETH3_masking.sh
run_step VARIMETH4 ./VARIMETH4_freebayes.sh
run_step VARIMETH5 ./VARIMETH5_filtering_SNPs.sh

echo "All VARIMETH scripts completed successfully!"
