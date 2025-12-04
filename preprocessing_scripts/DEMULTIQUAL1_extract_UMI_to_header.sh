#!/bin/bash
set -euo pipefail
exec > >(tee -i DEMULTIQUAL1_extract_UMI_to_header.log) 2>&1

RAW_DIR="/home/andres/datos/INTERTROPHIC/EPIGENÉTICA/NOVOGENE_DOWNLOAD/X204SC24100669-Z01-F001_02/01.RawData/Plantago_Lib/raw"
OUT_DIR="./DEMULTIQUAL1_extract_UMI_to_header"

mkdir -p "$OUT_DIR"

# Lanes a procesar
LANES=("L005" "L007")

for L in "${LANES[@]}"; do

    R1_IN="${RAW_DIR}/Plantago_Lib_${L}_R1_${L##L}.fastq.gz"
    R2_IN="${RAW_DIR}/Plantago_Lib_${L}_R2_${L##L}.fastq.gz"

    R1_OUT="${OUT_DIR}/Plantago_Lib_${L}_R1_${L##L}.umi.fastq.gz"
    R2_OUT="${OUT_DIR}/Plantago_Lib_${L}_R2_${L##L}.umi.fastq.gz"

    echo "Procesando lane ${L} → extrayendo UMI (3 bases) inline → header…"

    #############################################
    # R1
    #############################################
    awk '
        NR % 4 == 1 { hdr = $0 }
        NR % 4 == 2 { umi = substr($0,1,3); seq = substr($0,4) }
        NR % 4 == 3 { plus = $0 }
        NR % 4 == 0 {
            qual = substr($0,4)
            print hdr "+" umi
            print seq
            print plus
            print qual
        }
    ' <(zcat "$R1_IN") | gzip > "${R1_OUT}.tmp"

    #############################################
    # R2
    #############################################
    awk '
        NR % 4 == 1 { hdr = $0 }
        NR % 4 == 2 { umi = substr($0,1,3); seq = substr($0,4) }
        NR % 4 == 3 { plus = $0 }
        NR % 4 == 0 {
            qual = substr($0,4)
            print hdr "+" umi
            print seq
            print plus
            print qual
        }
    ' <(zcat "$R2_IN") | gzip > "${R2_OUT}.tmp"

    mv -f "${R1_OUT}.tmp" "$R1_OUT"
    mv -f "${R2_OUT}.tmp" "$R2_OUT"

    echo "Listo:"
    echo "  $R1_OUT"
    echo "  $R2_OUT"
    echo

done

echo ">>> Todos los lanes procesados correctamente."

