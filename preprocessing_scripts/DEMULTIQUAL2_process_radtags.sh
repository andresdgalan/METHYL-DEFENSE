#!/bin/bash
# Script to Demultiplex samples according to adapter barcodes using process_radtags from stacks software
exec > >(tee -i DEMULTIQUAL2_process_radtags.log) 2>&1

# Directorio de entrada
INPUT_DIR="./DEMULTIQUAL1_extract_UMI_to_header"

# Archivo de códigos de barras
BARCODES_FILE="./DEMULTIQUAL2_barcodes.txt"

# Directorio de salida (relativo al directorio actual)
OUTPUT_DIR="./DEMULTIQUAL2_process_radtags"

# Crear directorio de salida si no existe
mkdir -p "$OUTPUT_DIR"

echo "==== Iniciando process_radtags ===="

process_radtags \
    -P \
    -p "$INPUT_DIR" \
    -o "$OUTPUT_DIR" \
    -b "$BARCODES_FILE" \
    -e bamHI \
    -c -q -r \
    --inline-null \
    --disable-rad-check \
    -y fastq \
    --retain-header

echo "==== process_radtags finalizado ===="
