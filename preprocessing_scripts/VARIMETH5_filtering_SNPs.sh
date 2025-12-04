#!/bin/bash
exec > >(tee -i VARIMETH5_filtering_SNPs.log) 2>&1
# Activar conda dentro de script con el entorno donde está bcftools
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bcft_env

set -euo pipefail

VCF="VARIMETH4_freebayes/variants.vcf"
OUTPUT="filtered_variants_cov5-1000_biallelic_MAF0.01_80pct.vcf.gz"
TREATMENT_FILE="VARIMETH5_plantago_treatment.csv"

echo "Filtrado de SNPs: $VCF"

# 1. Coverage filter (DP 5 - 1000)

vcf_step1=$(mktemp --suffix=.vcf.gz)
echo "=== Coverage Filter (DP 5–1000) ==="
bcftools view -i 'MIN(FMT/DP)>=5 && MAX(FMT/DP)<=1000' -Oz -o "$vcf_step1" "$VCF"
sites_step1=$(bcftools view -H "$vcf_step1" | wc -l)
echo "Sites after coverage filtering (DP 5–1000): $sites_step1"

# 2. Keep only biallelic SNPs

vcf_step2=$(mktemp --suffix=.vcf.gz)
echo "=== Keeping biallelic SNPs ==="
bcftools view -m2 -M2 -v snps -Oz -o "$vcf_step2" "$vcf_step1"
sites_step2=$(bcftools view -H "$vcf_step2" | wc -l)
echo "Sites after keeping only biallelic SNPs: $sites_step2"

# 3. Minimum Allele Frequency MAF = 0.01

vcf_step3=$(mktemp --suffix=.vcf.gz)
echo "=== MAF >= 0.01 ==="
bcftools view -i 'MAF>=0.01' -Oz -o "$vcf_step3" "$vcf_step2"
sites_step3=$(bcftools view -H "$vcf_step3" | wc -l)
echo "Sites after MAF >= 0.01: $sites_step3"

# 4. Keep SNPs present in at least 80% of the samples per treatment

echo "=== Filter out SNPs with missingness in >80% of the samples per treatment ==="

declare -a treatments=("graz_clip" "graz_unclip" "ungraz_clip" "ungraz_unclip")

# Normalizar saltos de línea del CSV

dos2unix "$TREATMENT_FILE" 2>/dev/null || true

# Create sample lists for each treatment

for t in "${treatments[@]}"; do
maternal="${t%%_*}"
simulation="${t##*_}"
awk -F';' -v m="$maternal" -v s="$simulation" '
BEGIN {OFS=""}
NR>1 {
gsub(/\r/,"");
gsub(/^[ \t]+|[ \t]+$/,"",$1);
gsub(/^[ \t]+|[ \t]+$/,"",$2);
gsub(/^[ \t]+|[ \t]+$/,"",$3);
if ($2==m && $3==s) print $1, "_R1_trimmed_bismark_bt2_pe.sorted";
}' "$TREATMENT_FILE" > "$t.txt"
echo "  → $(wc -l < "$t.txt") samples for $t"
done

# Create sub-VCFs for each treatment and filter each one

for t in "${treatments[@]}"; do
echo "Filtrando tratamiento $t..."
bcftools view -S "$t.txt" -Oz -o "$t.vcf.gz" "$vcf_step3"
bcftools index "$t.vcf.gz"
bcftools view -i 'COUNT(GT!="mis")/N_SAMPLES>=0.8' -Oz -o "$t.filtered.vcf.gz" "$t.vcf.gz"
bcftools query -f '%CHROM\t%POS\n' "$t.filtered.vcf.gz" > "$t.pos"
done

# Intersection of SNPs present in all 4 treatments

sort graz_clip.pos graz_unclip.pos ungraz_clip.pos ungraz_unclip.pos | uniq -c | awk '$1==4 {print $2 "\t" $3}' > shared.pos

# Keep only those SNPs from the MAF-filtered vcf

vcf_step4=$(mktemp --suffix=.vcf.gz)
bcftools view -T shared.pos -Oz -o "$vcf_step4" "$vcf_step3"
sites_step4=$(bcftools view -H "$vcf_step4" | wc -l)
echo "Sites after per-treatment ≥80% genotyping: $sites_step4"

# Create and index final output

mv "$vcf_step4" "$OUTPUT"
bcftools index -f "$OUTPUT"

# Move temporal files to an intermediate folder

# Create folder
mkdir -p intermediate_filtering_steps

# Move temporal files
mv -f "$vcf_step1" "$vcf_step2" "$vcf_step3" \
   graz_clip.* graz_unclip.* ungraz_clip.* ungraz_unclip.* shared.pos \
   intermediate_filtering_steps/ 2>/dev/null


echo "Filtring completed: $OUTPUT ready."

# ================================================================
# Calculate missingness per sample in the final VCF
# ================================================================

echo
echo "===Summary of per sample missingness==="

samples=($(bcftools query -l "$OUTPUT"))
n_sites=$(bcftools view -H "$OUTPUT" | wc -l)

printf "%-20s %-10s %-10s %s\n" "Muestra" "Missing" "Total" "Porcentaje"

for s in "${samples[@]}"; do
    missing=$(bcftools query -s "$s" -f '[%GT\n]' "$OUTPUT" | grep -Ewc '\./\.|\.')
    pct=$(echo "scale=2; 100 * $missing / $n_sites" | bc)
    printf "%-20s %-10s %-10s %s%%\n" "$s" "$missing" "$n_sites" "$pct"
done | sort -k4 -n

echo "=============================================================="
