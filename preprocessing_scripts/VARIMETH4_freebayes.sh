#!/bin/bash
# Variant calling using freebayes with epiGBS2 default parameters
exec > >(tee -i VARIMETH4_freebayes.log) 2>&1

# Activar conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bismark_env


set -e
set -o pipefail

INPUT_DIR="VARIMETH3_masking"
GENOME="DENOVO67_genome_construction/pseudogenome_fixed.fa"  # Use the same FASTA as with samtools calmd
FREEBAYES_OUT_DIR="VARIMETH4_freebayes"

# Create output directory
mkdir -p "$FREEBAYES_OUT_DIR"

    echo "Running FreeBayes on $INPUT_DIR/merged_masked.bam..."

    freebayes -f "$GENOME" \
        --no-partial-observations \
        --report-genotype-likelihood-max \
        --genotype-qualities \
        --min-coverage 0 \
        --min-base-quality 1 \
        --min-mapping-quality 10 \
        --no-population-priors \
        "$INPUT_DIR"/merged_masked.bam > "$FREEBAYES_OUT_DIR"/variants.vcf

    echo "Output written to $OUTPUT_VCF"

echo "All FreeBayes runs completed successfully!"

