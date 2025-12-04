#!/bin/bash
# Script para procesar archivos truncados con ustacks
exec > >(tee -i DENOVO4_ustacks.log) 2>&1

# Carpetas
input_folder="./DENOVO3_join"
output_folder="./DENOVO45_stacks"

# Crear carpeta de salida
mkdir -p "$output_folder"

# Número de muestras, que va incrementando por cada archivo
counter=1

# Bucle para recorrer todos los archivos .fq truncados
for file in "$input_folder"/*_joined.fq; do
    echo "Procesando $file con ustacks..."
    ustacks -f "$file" --disable-gapped -p 15 -M 3 -i $counter -o "$output_folder"
    counter=$((counter + 1))
done

