# In this script we use fastp to trim adapters, polyGs derived from Ilumina errors and bad quality regions
# --trim_poly_g only removes polyGs in the tails so we used --low_complexity_filter to remove polyGs out of the tails
#!/bin/bash
exec > >(tee -i DEMULTIQUAL4_quality_adapter_trim.log) 2>&1

input_folder="./DEMULTIQUAL3_clone_filter"       
output_folder="./DEMULTIQUAL4_quality_adapter_trim"
mkdir -p "$output_folder"

# For some reason, during clone_filter, an extra .1 and .2 suffix was added to files that already had .1 and .2
# Removing redundant suffixes from file names
for f in "$input_folder"/*.fq; do
    if [[ "$f" =~ \.1\.1\.fq$ ]]; then
        mv "$f" "${f/.1.1.fq/.1.fq}"
    elif [[ "$f" =~ \.2\.2\.fq$ ]]; then
        mv "$f" "${f/.2.2.fq/.2.fq}"
    fi
done

# Loop through all R1 files (except those with ".rem.")
for r1 in "$input_folder"/*.1.fq; do
    # Saltar si el archivo contiene ".rem."
    if [[ "$r1" == *".rem."* ]]; then
        continue
    fi

    # Sample name (without .1.fq)
    base=$(basename "$r1" .1.fq)

    # Define corresponding R2 file
    r2="$input_folder/${base}.2.fq"

    # Make sure R2 exists
    if [[ ! -f "$r2" ]]; then
        echo "R2 file for $base not found, skipping..."
        continue
    fi

    # output files
    out_r1="$output_folder/${base}_R1_trimmed.fq"
    out_r2="$output_folder/${base}_R2_trimmed.fq"
    report_html="$output_folder/${base}_fastp_report.html"
    report_json="$output_folder/${base}_fastp_report.json"

    echo "Processing sample $base"
    fastp \
      -i "$r1" -I "$r2" \
      -o "$out_r1" -O "$out_r2" \
      --trim_poly_g \
      --trim_front1 5 \
      --trim_front2 9 \
      --low_complexity_filter \
      --complexity_threshold 30 \
      --disable_adapter_trimming \
      --dont_eval_duplication \
      --cut_front \
      --cut_tail \
      --cut_window_size 4 \
      --cut_mean_quality 20 \
      --qualified_quality_phred 15 \
      --length_required 30 \
      --thread 8 \
      --html "$report_html" \
      --json "$report_json"
done

