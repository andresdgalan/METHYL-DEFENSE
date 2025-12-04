#!/bin/bash
# In this script we align our samples with the de novo constructed reference using bismark.
# All per-sample alignment files are indexed and merged to a single combined alignment file using SamToolS.
# The indexed bams will be used for methylation calling with methylKit in R.
# The merged bam will be used for variant calling in the next steps of this pipeline.
#exec > >(tee -i VARIMETH1_bismark.log) 2>&1
set -euo pipefail

FASTQ_DIR=./DEMULTIQUAL4_quality_adapter_trim
OUTPUT_DIR=./VARIMETH1_bismark
GENOME=./DENOVO67_genome_construction
THREADS=15

mkdir -p "$OUTPUT_DIR"
: << 'FIN'
# echo "[INFO] Starting Bismark mapping..."

# -----------------------------
# 1. Run Bismark for each sample
# -----------------------------

#  Activate bismark environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bismark_env

#for R1 in "$FASTQ_DIR"/*_R1_trimmed.fq; do
#
    # Skip non-treated replicates
#    if [[ "$R1" == *"_2_"* ]]; then
#        continue
#    fi

 #   R2="${R1/_R1_trimmed.fq/_R2_trimmed.fq}"
#    if [[ ! -f "$R2" ]]; then
#        echo "[WARNING] R2 not found for $R1, skipping..."
#        continue
#    fi

#    echo "[INFO] Mapping: $R1 and $R2"
#    bismark --non_directional \
#            --genome "$GENOME" \
#            --output_dir "$OUTPUT_DIR" \
#            -p "$THREADS" \
#            -1 "$R1" \
#            -2 "$R2"
#done

# echo "[INFO] Bismark mapping finished."
echo "[INFO] Starting sorting and indexing..."

# -----------------------------
# 2. Sort and index each BAM produced by Bismark
# -----------------------------
cd "$OUTPUT_DIR"

for bam in *_bismark_bt2_pe.bam; do
    base=$(basename "$bam" .bam)
    echo "[INFO] Sorting $bam"
    samtools sort -o "${base}.sorted.bam" "$bam"

    echo "[INFO] Indexing ${base}.sorted.bam"
    samtools index "${base}.sorted.bam"
done

echo "[INFO] Sorting and indexing completed."
FIN
# -----------------------------
# 3. Add Read Groups and merge sorted BAMs
# -----------------------------
echo "[INFO] Adding Read Groups to sorted BAM files (excluding those containing '_2')..."
cd "$OUTPUT_DIR"
RG_DIR="bam_with_rg"
mkdir -p "$RG_DIR"

for bam in *.sorted.bam; do
    if [[ "$bam" == *_2* ]]; then
        echo "[INFO] Skipping $bam because it contains '_2'"
        continue
    fi

    # Remove the final .bam
    sample="${bam%.bam}"

    echo "[INFO] Adding RG to $bam (sample=$sample)"

    samtools addreplacerg \
        -r "ID:$sample" \
        -r "SM:$sample" \
        -r "PL:ILLUMINA" \
        -o "$RG_DIR/${sample}.RG.bam" \
        "$bam"

    samtools index -c "$RG_DIR/${sample}.RG.bam"
done

echo "[INFO] Read Groups added to all selected BAM files."

# -----------------------------
# Merge de BAMs con RG
# -----------------------------
echo "[INFO] Merging BAMs with RG information..."
samtools merge -@ "$THREADS" -f merged_all_samples.bam "$RG_DIR"/*.RG.bam

echo "[INFO] Indexing merged BAM..."
samtools index -c merged_all_samples.bam

echo "[INFO] DONE. Combined BAM generated with proper Read Groups:"
echo "       $OUTPUT_DIR/merged_all_samples.bam"
