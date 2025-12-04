#!/bin/bash
# Script to join R1 and R2 FASTQ files in current directory using vsearch
exec > >(tee -i DENOVO3_join.log) 2>&1

mkdir -p DENOVO3_join

for r1 in ./DENOVO2_intersection/*_2_intersection_R1.fq; do
    r2="${r1/_R1/_R2}"   # cambia _R1_ por _R2_
    base=$(basename "$r1" _intersection_R1.fq)  # quita la ruta y el sufijo
    if [[ -f "$r2" ]]; then
        echo "Joining $r1 and $r2 ..."
        vsearch --fastq_join "$r1" \
                --reverse "$r2" \
                --join_padgap ATATATAT \
                --join_padgapq IIIIIIII \
                --fastqout "DENOVO3_join/${base}_joined.fq"
    else
        echo "Warning: Matching R2 file for $r1 not found."
    fi
done

