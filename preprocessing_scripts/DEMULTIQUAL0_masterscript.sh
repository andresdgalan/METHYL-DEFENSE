#!/bin/bash
set -e  # stop if any script fails

# Master script to run DEMULTIQUAL pipeline in order
# To run this pipeline we ONLY NEED RAW SEQUENCING DATA
# In DEMULTIQUAL pipeline we will obtain demultiplexed SAMPLES and a MultiQC quality reports which we need to evaluate before starting

# Default start is DEMULTIQUAL1, if you want to start at any other point run something like this:
# ./DEMULTIQUAL0_masterscript.sh DEMULTIQUAL3

start_at=${1:-"DEMULTIQUAL1"}

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

run_step DEMULTIQUAL1 ./DEMULTIQUAL1_extract_UMI_to_header.sh
run_step DEMULTIQUAL2 ./DEMULTIQUAL2_process_radtags.sh
run_step DEMULTIQUAL3 ./DEMULTIQUAL3_clone_filter.sh
run_step DEMULTIQUAL4 ./DEMULTIQUAL4_quality_adapter_trim.sh
run_step DEMULTIQUAL5 ./DEMULTIQUAL5_multiqc.sh

echo "All DEMULTIQUAL scripts completed successfully!"
