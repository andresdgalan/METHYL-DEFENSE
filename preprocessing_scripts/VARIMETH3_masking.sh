#!/bin/bash
# Processing bam file with a double-masking procedure.
# Ensuring that variant calling is not affected by the bisulfite conversion, as explained in Nunn et al. (2021)

# Activar conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate variant_env

exec > >(tee -i VARIMETH3_masking.log) 2>&1

set -e
set -o pipefail

INPUT_DIR="VARIMETH2_MDtag"
SCRIPT="./VARIMETH3_revelio.py"
TEMP_DIR="${INPUT_DIR}/tmp_revelio"
OUTPUT_DIR="VARIMETH3_masking"

# Create output folder
mkdir -p "$OUTPUT_DIR"

# Create a temporal folder
mkdir -p "$TEMP_DIR"

echo "Running ReVeLio on $INPUT_DIR/merged.calmd.bam..."

# Run revelio.py
    python3 "$SCRIPT" \
        "$INPUT_DIR"/merged.calmd.bam \
        "$OUTPUT_DIR"/merged_masked.bam \
        -t "$TEMP_DIR" \
        -T 15  # You may adjust threads

    echo "Output written to $OUTPUT_DIR/merged_masked.bam"

# Cleaning
# rm -rf "$TEMP_DIR" # he comentado esto para mirar qué hay, porque estoy teniendo problemas
echo "All files processed successfully!"

