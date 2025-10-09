#!/bin/bash
set -euo pipefail

# -----------------------------
# Variables
# -----------------------------
WD=~/nuevo_datos/Plantago_lib
RAW_FASTQ="/home/andres/nuevo_datos/Plantago_lib/samples/PL26_2.1.fq"
SUBSAMPLED_FASTQ="$WD/PL26_2_1.subsampled.fq"
SUBSAMPLED_FASTA="$WD/PL26_2_1.subsampled.fa"
BLAST_OUT="$WD/PL26_2_1_blast.tsv"
ACCESSIONS="$WD/accessions_324_2.1.txt"
ACCESSION_TAXID="$WD/accession_324_2.1_taxid.tsv"
TAXID2NAME="$WD/taxid2name.tsv"
ACCESSION_TAXID_NAME="$WD/accession_taxid_name_324_2.1.tsv"
SPECIES_COUNTS="$WD/counts_by_species_324_2.1.tsv"

# -----------------------------
# 1. Subsample reads and convert to FASTA
# -----------------------------
echo "[INFO] Subsampling 1000 reads..."
seqtk sample -s100 "$RAW_FASTQ" 1000 > "$SUBSAMPLED_FASTQ"

echo "[INFO] Converting to FASTA..."
seqtk seq -a "$SUBSAMPLED_FASTQ" > "$SUBSAMPLED_FASTA"

# -----------------------------
# 2. Run BLASTn remotely
# -----------------------------
echo "[INFO] Running BLASTn..."
blastn -query "$SUBSAMPLED_FASTA" -db nt -remote \
       -max_target_seqs 5 \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids ssciname" \
       -out "$BLAST_OUT"

# -----------------------------
# 3. Download accession→taxid mapping
# -----------------------------
echo "[INFO] Downloading NCBI accession2taxid table..."
wget -O "$WD/nucl_gb.accession2taxid.gz" ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gunzip -f "$WD/nucl_gb.accession2taxid.gz"

# -----------------------------
# 4. Extract accessions from BLAST output
# -----------------------------
echo "[INFO] Extracting accessions from BLAST..."
cut -f2 "$BLAST_OUT" \
  | sed 's/gi|[0-9]*|[a-z]*|//' \
  | sed 's/|.*//' \
  > "$ACCESSIONS"

# -----------------------------
# 5. Map accessions → taxid
# -----------------------------
echo "[INFO] Mapping accessions to taxids..."
join -t $'\t' \
  <(sort -k1,1 "$ACCESSIONS") \
  <(cut -f2,3 "$WD/nucl_gb.accession2taxid" | sort -k1,1) \
  > "$ACCESSION_TAXID"

# -----------------------------
# 6. Download NCBI taxonomy database
# -----------------------------
echo "[INFO] Downloading NCBI taxdump..."
wget -O "$WD/new_taxdump.tar.gz" https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar -xvzf "$WD/new_taxdump.tar.gz" -C "$WD" names.dmp nodes.dmp

# -----------------------------
# 7. Extract taxid → scientific name
# -----------------------------
echo "[INFO] Extracting taxid → scientific name..."
awk -F'\t|\t' '$7=="scientific name" {print $1 "\t" $3}' "$WD/names.dmp" \
  > "$TAXID2NAME"

# -----------------------------
# 8. Clean & sort for join
# -----------------------------
echo "[INFO] Cleaning and sorting files for join..."
awk '{$1=$1; $2=$2; print}' OFS='\t' "$ACCESSION_TAXID" \
  | sort -k2,2 > "$WD/accession_taxid_clean_sorted.tsv"

awk '{$1=$1; $2=$2; print}' OFS='\t' "$TAXID2NAME" \
  | sort -k1,1 > "$WD/taxid_name_clean_sorted.tsv"

awk -F $'\t' 'NF==2' "$ACCESSION_TAXID" > "$WD/accession_taxid_clean.tsv"

# -----------------------------
# 9. Join accession→taxid with taxid→name
# -----------------------------
echo "[INFO] Joining accession→taxid with taxid→name..."
join -t $'\t' -1 2 -2 1 \
  <(cat "$WD/accession_taxid_clean_sorted.tsv") \
  "$WD/taxid_name_clean_sorted.tsv" \
  > "$ACCESSION_TAXID_NAME"

echo "[INFO] Example lines from joined file:"
head -n 10 "$ACCESSION_TAXID_NAME"
wc -l "$ACCESSION_TAXID_NAME"

# -----------------------------
# 10. Count hits per species
# -----------------------------
echo "[INFO] Counting hits per species..."
cut -f3 "$ACCESSION_TAXID_NAME" \
  | sort \
  | uniq -c \
  | sort -nr \
  > "$SPECIES_COUNTS"

echo "[INFO] Top species counts:"
head -n 20 "$SPECIES_COUNTS"

echo "[INFO] All steps completed. Output files:"
echo "  - $ACCESSION_TAXID_NAME"
echo "  - $SPECIES_COUNTS"

