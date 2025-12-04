#!/bin/bash
# After indexing the merged bam file, we add the MD tag: a string encoding mismatched and deleted reference bases
# MD tag is added using SamToolS calmd
exec > >(tee -i VARIMETH2_MDtag.log) 2>&1

set -e
set -o pipefail

# Input directory with BAMs
INPUT_DIR="VARIMETH1_bismark"
# Original genome
GENOME_ORIG="DENOVO67_genome_construction/pseudogenome.fa"
# Fixed genome output
GENOME_FIXED="DENOVO67_genome_construction/pseudogenome_fixed.fa"
# Output directory with tagged BAMs
OUTPUT_DIR="VARIMETH2_MDtag"
mkdir -p "$OUTPUT_DIR"

# Fix FASTA format
echo "Reformatting FASTA for samtools..."
{
    # Write the header
    grep '^>' "$GENOME_ORIG"
    # Concatenate all sequence lines and wrap to 60 bases per line
    grep -v '^>' "$GENOME_ORIG" | tr -d '\n' | fold -w 60
} > "$GENOME_FIXED"

# Index genome if index does not exist
    echo "Indexing genome..."
    samtools faidx "$GENOME_FIXED"

# Add MD tag and index output
    samtools calmd -b "$INPUT_DIR"/merged_all_samples.bam "$GENOME_FIXED" > "$OUTPUT_DIR"/merged.calmd.bam
    samtools index "$OUTPUT_DIR"/merged.calmd.bam


    echo "Output written to $OUTPUT_FILE"

echo "All files processed successfully!"

