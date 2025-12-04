# De novo reference construction with the stacks pipeline requires homogeneity in read lenght
# We have first checked that almost all the R1 reads have 129 bases and the R2 reads have 135 bases
# In this script we truncate all R1 reads to 129 bases and R2 reads to 135 bases
exec > >(tee -i DENOVO1_trunclen.log) 2>&1

mkdir -p DENOVO1_trunclen

for f in trimmed/*_2_R[12]_trimmed.fq; do
    base=$(basename "$f" _trimmed.fq)

# exclude bad quality samples
    if [[ "$base" == PL49_2_R1 || "$base" == PL49_2_R2 || \
          "$base" == PL6_2_R1 || "$base" == PL6_2_R2 || \
          "$base" == PL25_2_R1 || "$base" == PL25_2_R2 ]]; then
        echo "Excluyendo $f por baja calidad"
        continue
    fi

    if [[ "$f" == *_2_R1_trimmed.fq ]]; then
        len=129
    elif [[ "$f" == *_2_R2_trimmed.fq ]]; then
        len=135
    else
        echo "Skipping $f (not R1 or R2)"
        continue
    fi

    echo "Truncating $f to $len bp..."
    vsearch --fastq_filter "$f" \
        --fastq_trunclen "$len" \
        --fastqout "DENOVO1_trunclen/${base}_trunc.fq"
done

