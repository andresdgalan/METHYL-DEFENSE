# En este sscript voy a hacer pruebas para el GO enrichment de Avena
# web para utilizar el paquete: https://yulab-smu.top/biomedical-knowledge-mining-book/021-go.html


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")

library(here) # reproducibility
library(methylKit)
library(clusterProfiler) # enrichment (function enricher)
library(dplyr) # tidy
library(tidyr) # tidy
library(purrr) # tidy
library(genomation) # anotación de promotores y visualización de features
library(rtracklayer) # importar gff


#load data ----

avena_CpG_GUvsUU <- readRDS(here("avena_ncbi","DMPs_CpG_GUvsUU.rds")) # background
avena_CpG_UUvsUC <- readRDS(here("avena_ncbi","DMPs_CpG_UUvsUC.rds")) # background
avena_CpG_GUvsGC <- readRDS(here("avena_ncbi","DMPs_CpG_GUvsGC.rds")) # background

avena_CHG_GUvsUU <- readRDS(here("avena_ncbi","DMPs_CHG_GUvsUU.rds")) # background
avena_CHG_UUvsUC <- readRDS(here("avena_ncbi","DMPs_CHG_UUvsUC.rds")) # background
avena_CHG_GUvsGC <- readRDS(here("avena_ncbi","DMPs_CHG_GUvsGC.rds")) # background

avena_CHH_GUvsUU <- readRDS(here("avena_ncbi","DMPs_CHH_GUvsUU.rds")) # background
avena_CHH_UUvsUC <- readRDS(here("avena_ncbi","DMPs_CHH_UUvsUC.rds")) # background
avena_CHH_GUvsGC <- readRDS(here("avena_ncbi","DMPs_CHH_GUvsGC.rds")) # background

avena_all_GUvsUU <- read.table(here("avena_ncbi/DMPs","DMPs_all_GUvsUU.txt"), header = T) # sig
avena_all_UUvsUC <- read.table(here("avena_ncbi/DMPs","DMPs_all_UUvsUC.txt"), header = T) # sig
avena_all_GUvsGC <- read.table(here("avena_ncbi/DMPs","DMPs_all_GUvsGC.txt"), header = T) # sig

GOfeature <- read.table(here("GO_enrichment/Avena/AVBAR.PI388828.10000a.pgsb.r1.Feb2023.feature_table.csv"), header = T, sep = ",")

# ENRICHMENT BACKGROUND ----

# 1. Crear una lista con todos los sitios evaluados (background)
lista_objetos <- list(
  avena_CpG_GUvsGC$all, avena_CpG_GUvsUU$all, avena_CpG_UUvsUC$all,
  avena_CHG_GUvsGC$all, avena_CHG_GUvsUU$all, avena_CHG_UUvsUC$all,
  avena_CHH_GUvsGC$all, avena_CHH_GUvsUU$all, avena_CHH_UUvsUC$all
)

# 2. Convertir a data frame, seleccionar 'chr' y 'start', unir y eliminar filas redundantes completas
background_Cs <- lista_objetos %>%
  map_dfr(~ getData(.x) %>% dplyr::select(chr, start, end, strand)) %>%
  distinct(chr, start, .keep_all = TRUE)



# CAMBIO DE NOMENCLATURA ----
# El objeto gff tiene una nomenclatura distinta a la de nuestros DMPs para los cromosomas
# Vamos a cambiar los objetos de DMPs para que tengan el nombre real del cromosoma,
# y no el de genbank

# Function to convert GenBank sequence names to chromosome/contig names
convert_genbank_to_chr <- function(df, mapping_file) {
  
  # Rename the original chromosome column
  df$genbank <- df$chr
  df$chr <- NULL
  
  # Read GenBank-to-chr mapping
  seqmap <- read.delim(
    mapping_file,
    header = FALSE,
    sep = "\t",
    col.names = c("genbank", "chr"),
    stringsAsFactors = FALSE
  )
  
  # Check that GenBank IDs are unique in the mapping
  stopifnot(!anyDuplicated(seqmap$genbank))
  
  # Check that all GenBank IDs in the data have a mapping
  missing_genbank <- setdiff(
    unique(df$genbank),
    seqmap$genbank
  )
  
  if (length(missing_genbank) > 0) {
    stop(
      "These GenBank IDs are missing from the mapping: ",
      paste(missing_genbank, collapse = ", ")
    )
  }
  
  # Create named lookup vector
  genbank_to_chr <- setNames(seqmap$chr, seqmap$genbank)
  
  # Create chr from genbank
  df$chr <- unname(
    genbank_to_chr[df$genbank]
  )
  
  # Check that no chromosome names are missing
  stopifnot(!anyNA(df$chr))
  
  # Move chr to the first column and genbank to the second
  df <- df[, c(
    "chr",
    "genbank",
    setdiff(names(df), c("chr", "genbank"))
  )]
  
  return(df)
}

# MAPPING FILE - correspondencia de genbank con el nombre de cada contig
mapping_file <- here(
  "GO_enrichment",
  "Avena",
  "Abarbata_genbank_to_chr.txt"
)

# APLICAMOS LA FUNCIÓN PARA CADA ARCHIVO
avena_all_GUvsUU <- convert_genbank_to_chr(
  avena_all_GUvsUU,
  mapping_file
)

avena_all_UUvsUC <- convert_genbank_to_chr(
  avena_all_UUvsUC,
  mapping_file
)

avena_all_GUvsGC <- convert_genbank_to_chr(
  avena_all_GUvsGC,
  mapping_file
)

background_Cs <- convert_genbank_to_chr(
  background_Cs,
  mapping_file
)


# ANNOTATION ----


# Como lo hace McNew
# Incluir promotores 2kb upstream del TSS
gene.obj <- readTranscriptFeatures(here("GO_enrichment/Avena/longest.bed"),remove.unusual=FALSE,
                                   up.flank=2000,down.flank=0,unique.prom=TRUE)

GUvsUU_annotation <- annotateWithGeneParts(as(avena_all_GUvsUU, "GRanges"), gene.obj)
UUvsUC_annotation <- annotateWithGeneParts(as(avena_all_UUvsUC, "GRanges"), gene.obj)
GUvsGC_annotation <- annotateWithGeneParts(as(avena_all_GUvsGC, "GRanges"), gene.obj)
background_annotation <- annotateWithGeneParts(as(background_Cs, "GRanges"), gene.obj)

# Corregimos feature.name para que sea el nombre del gen
GUvsUU_annotation@dist.to.TSS$feature.name <- sub(
  "^AVBAR\\.10000a\\.r1\\..*G([0-9]+)\\.\\d+$",
  "AVBAR.10000a.r1.1AG\\1",
  GUvsUU_annotation@dist.to.TSS$feature.name
)
UUvsUC_annotation@dist.to.TSS$feature.name <- sub(
  "^AVBAR\\.10000a\\.r1\\..*G([0-9]+)\\.\\d+$",
  "AVBAR.10000a.r1.1AG\\1",
  UUvsUC_annotation@dist.to.TSS$feature.name
)
GUvsGC_annotation@dist.to.TSS$feature.name <- sub(
  "^AVBAR\\.10000a\\.r1\\..*G([0-9]+)\\.\\d+$",
  "AVBAR.10000a.r1.1AG\\1",
  GUvsGC_annotation@dist.to.TSS$feature.name
)
background_annotation@dist.to.TSS$feature.name <- sub(
  "^AVBAR\\.10000a\\.r1\\..*G([0-9]+)\\.\\d+$",
  "AVBAR.10000a.r1.1AG\\1",
  background_annotation@dist.to.TSS$feature.name
)

# Function to create annotation table and calculate feature proportions
summarise_annotation <- function(annotation, name) {
  
  anno <- annotation@members %>%
    as.data.frame() %>%
    mutate(target.row = seq_len(nrow(annotation@members))) %>%
    merge(annotation@dist.to.TSS, by = "target.row")
  
  # Save annotation table
  write.csv(
    anno,
    here("GO_enrichment", "Avena", paste0(name, "_annotation.csv")),
    row.names = FALSE
  )
  
  # Number and percentage of cytosines in each feature
  summary <- data.frame(
    Feature = c("Promoter", "Exon", "Intron"),
    N = c(
      sum(anno$prom == 1),
      sum(anno$exon == 1),
      sum(anno$intron == 1)
    ),
    Percentage = c(
      mean(anno$prom == 1) * 100,
      mean(anno$exon == 1) * 100,
      mean(anno$intron == 1) * 100
    )
  )
  
  cat("\n", name, "\n", sep = "")
  print(summary)
  
  invisible(anno)
}

# Generate tables and summaries
GUvsUU_anno <- summarise_annotation(GUvsUU_annotation, "GUvsUU")
UUvsUC_anno <- summarise_annotation(UUvsUC_annotation, "UUvsUC")
GUvsGC_anno <- summarise_annotation(GUvsGC_annotation, "GUvsGC")
background_anno <- summarise_annotation(background_annotation, "background")



# ENRICHMENT ----

# Creamos TERM2GENE seleccionando las columnas gID y GO
library(dplyr)
library(tidyr)

TERM2GENE <- GOfeature %>%
  select(gID, GO) %>%
  filter(!is.na(GO), GO != "", GO != "nan") %>%
  separate_rows(GO, sep = ",") %>%
  filter(grepl("^GO:\\d+$", GO)) %>%
  select(term = GO, gene = gID)

# Overrepresented biological processes, molecular functions, and cellular 
# components were identified with an FDR-adjusted threshold of alpha-value 
# < 0.05 (Benjamini & Hochberg, 1995).

run_go_enrichment <- function(gene_set) {
  
  enricher(
    gene = unique(gene_set),
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe = background_anno$feature.name,
    TERM2GENE = TERM2GENE
  )
}


GUvsUU_GO <- run_go_enrichment(
  GUvsUU_anno$feature.name)

UUvsUC_GO <- run_go_enrichment(
  UUvsUC_anno$feature.name)

GUvsGC_GO <- run_go_enrichment(
  GUvsGC_anno$feature.name)


