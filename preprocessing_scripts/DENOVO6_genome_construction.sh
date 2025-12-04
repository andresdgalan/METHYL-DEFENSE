#!/bin/bash
# Script to create pseudogenome.fa from catalog.tags.tsv
exec > >(tee -i DENOVO6_genome_construction.log) 2>&1
# Define output folder
OUTPUT_DIR="./DENOVO67_genome_construction"

# Create the folder if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Define input and output files
INPUT_FILE="DENOVO45_stacks/catalog.tags.tsv"
OUTPUT_FILE="$OUTPUT_DIR/pseudogenome.fa"

# Generate the FASTA
awk '$3=="consensus"{print ">Plantago_"$2"\n"$6}' "$INPUT_FILE" > "$OUTPUT_FILE"

# Optional message
echo "pseudogenome.fa has been created in $OUTPUT_DIR"

