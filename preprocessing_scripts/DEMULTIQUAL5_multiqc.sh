# In this script we create a MultiQC report to evaluate quality filtering
# and decide if we need to remove  any bad quality samples before continuing with the pipeline
#!/bin/bash
exec > >(tee -i DEMULTIQUAL5_multiqc.log) 2>&1

INPUT="./DEMULTIQUAL4_quality_adapter_trim"
OUTPUT="./DEMULTIQUAL5_multiqc"
mkdir -p "$OUTPUT"

# FastQC
echo "FastQC post-trimming..."
fastqc "$INPUT"/*_R1_trimmed.fq "$INPUT"/*_R2_trimmed.fq -o "$OUTPUT" -t 8

# MultiQC
echo "MultiQC post-trimming..."
multiqc "$OUTPUT" -o "$OUTPUT"

echo "FastQC and MultiQC completed."
