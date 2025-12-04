# In this script I keep only the reads in the intersection of R1 and R2
# This step is needed after quality filtering with fastp
# Our R2 reads have worse quality, and the quality trimming is assymetric
# This would results in an error when trying to join R1 and R2 in the next step (more R1 than R2, so it doesn't work properly)
exec > >(tee -i DENOVO2_intersection.log) 2>&1

#!/bin/bash
set -euo pipefail

mkdir -p DENOVO2_intersection

echo "Searching R1/R2 pairs in directory DENOVO1_trunclen/..."

# Loop through all R1 files in the input directory
for R1 in DENOVO1_trunclen/*_R1_trunc.fq; do
    # Get sample name (ej. PL22_2)
    prefix=$(basename "$R1" _R1_trunc.fq)
    R2="DENOVO1_trunclen/${prefix}_R2_trunc.fq"

    # Make sure R2 exist for each R1
    if [[ ! -f "$R2" ]]; then
        echo "o s2 pfound for $prefix"
        continue
    fi

    echo "Processing pair: $prefix"

    # Intermediate and output files
    R1_nosuf="DENOVO2_intersection/${prefix}_R1_nosuf.fq"
    R2_nosuf="DENOVO2_intersection/${prefix}_R2_nosuf.fq"
    ids_R1="DENOVO2_intersection/${prefix}_ids_R1.txt"
    ids_R2="DENOVO2_intersection/${prefix}_ids_R2.txt"
    filter="DENOVO2_intersection/${prefix}_intersection.txt"
    R1_out="DENOVO2_intersection/${prefix}_intersection_R1.fq"
    R2_out="DENOVO2_intersection/${prefix}_intersection_R2.fq"

    # Remove suffix that contains UMIs making impossible to recognise pairs as pairs
    sed 's/ .*//' "$R1" > "$R1_nosuf"
    sed 's/ .*//' "$R2" > "$R2_nosuf"


    # Extract iDs
    grep '^@' "$R1_nosuf" | sort > "$ids_R1"
    grep '^@' "$R2_nosuf" | sort > "$ids_R2"

    # IDs intersection
    comm -12 "$ids_R1" "$ids_R2" > "$filter"
    echo "   - IDs comunes: $(wc -l < "$filter")"

    # Filter R1 (and restore /1)
    awk -v ids="$filter" '
    BEGIN{
      while((getline line < ids) > 0) seen[line]=1
      close(ids)
    }
    NR%4==1 {
      hdr=$0; getline seq; getline plus; getline qual
      split(hdr,a," "); id=a[1]
      if(id in seen){
        print id"/1"; print seq; print plus; print qual
      }
    }' "$R1_nosuf" > "$R1_out"

    # Filter R2 (and restore /2)
    awk -v ids="$filter" '
    BEGIN{
      while((getline line < ids) > 0) seen[line]=1
      close(ids)
    }
    NR%4==1 {
      hdr=$0; getline seq; getline plus; getline qual
      split(hdr,a," "); id=a[1]
      if(id in seen){
        print id"/2"; print seq; print plus; print qual
      }
    }' "$R2_nosuf" > "$R2_out"

    echo "Intersection filter completed:"
    echo " $R1_out"
    echo " $R2_out"
    echo
done

echo "All pairs processed. Results in 'DENOVO2_intersection/'."

