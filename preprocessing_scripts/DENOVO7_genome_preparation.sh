#!/bin/bash
# Script to collapse pseudogenome FASTA to 1 chromosome and run Bismark genome preparation
exec > >(tee -i DENOVO7_genome_preparation.log) 2>&1

# Set variables
GENOME_DIR="./DENOVO67_genome_construction"
INPUT_FASTA="$GENOME_DIR/pseudogenome.fa"
MULTI_CHR_BACKUP="./pseudogenome_contigs.fa"  # outside genome folder
FIXED_FASTA="$GENOME_DIR/pseudogenome_single_chr.fa"

# Make a backup of the original multi-chromosome fasta
cp "$INPUT_FASTA" "$MULTI_CHR_BACKUP"
echo "Backup of multi-chromosome FASTA saved as $MULTI_CHR_BACKUP"

# Collapse to single chromosome
grep -v "^>" "$INPUT_FASTA" > tmp.fa
echo ">Plantago_genome" > "$FIXED_FASTA"
cat tmp.fa >> "$FIXED_FASTA"
rm tmp.fa
echo "Collapsed FASTA created: $FIXED_FASTA"

# replace original FASTA in genome folder with single chromosome version
mv "$FIXED_FASTA" "$INPUT_FASTA"
echo "Single chromosome FASTA ready: $INPUT_FASTA"

# 4 Activate bismark environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bismark_env

# 5 Run Bismark genome preparation
bismark_genome_preparation "$GENOME_DIR"
echo "Bismark genome preparation completed for $GENOME_DIR"

