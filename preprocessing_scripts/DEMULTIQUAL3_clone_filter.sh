# In this script we remove PCR clones using clone_filter from the stacks software
#!/bin/bash
exec > >(tee -i DEMULTIQUAL3_clone_filter.log) 2>&1

IN_DIR="./DEMULTIQUAL2_process_radtags"
OUT_DIR="./DEMULTIQUAL3_clone_filter"

mkdir -p "$OUT_DIR"

echo "==== Starting clone_filter ===="

# Loop through *.1.fq files that don't contain "rem"
for R1 in "$IN_DIR"/*.1.fq; do

    # Skip if "rem"
    if [[ "$R1" == *rem* ]]; then
        continue
    fi

    filename=$(basename "$R1")     # example: PL11_2.1.fq
    sample=${filename%.1.fq}       # removes ".1.fq" → PL11_2

    R2="${IN_DIR}/${sample}.2.fq"  # R2 pair construction
    

    # Make sure R2 exists
    if [[ ! -f "$R2" ]]; then
        echo "R2 not found for $R1."
        continue
    fi

    echo "Processing: $R1 y $R2"

    clone_filter \
        -1 "$R1" \
        -2 "$R2" \
        -P \
        --index_index \
        --oligo-len-1 3 \
        --oligo-len-2 3 \
        -o "$OUT_DIR" \
        -i fastq \
        -y fastq

done

echo "==== clone_filter finalized ===="

