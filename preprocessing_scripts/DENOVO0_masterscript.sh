#!/bin/bash
set -e  # stop if any script fails

# Master script to run DENOVO pipeline in order
# To run this script this pipeline we need input from the previous pipeline:
# DEMULTIQUAL
# In DEMULTIQUAL pipeline we obtained MultiQC quality reports which we need to evaluate before starting

# IMPORTANT ---
# Bad quality samples will be removed in DENOVO1_trunclen.sh and DENOVO5_PL_Map.txt

# Default start is DENOVO1, if you want to start at any other point run something like this:
# ./DENOVO0_masterscript.sh DENOVO3

start_at=${1:-"DENOVO1"}

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

run_step DENOVO1 ./DENOVO1_trunclen.sh
run_step DENOVO2 ./DENOVO2_intersection.sh
run_step DENOVO3 ./DENOVO3_join.sh
run_step DENOVO4 ./DENOVO4_ustacks.sh
run_step DENOVO5 ./DENOVO5_cstacks.sh
run_step DENOVO6 ./DENOVO6_genome_construction.sh
run_step DENOVO7 ./DENOVO7_genome_preparation.sh
run_step DENOVO8 ./DENOVO8_genome_position_map.sh

echo "All DENOVO scripts completed successfully!"

