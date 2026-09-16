# En este sscript voy a hacer pruebas para el GO enrichment de Avena

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")

library(clusterProfiler)

# web para utilizar el paquete: https://yulab-smu.top/biomedical-knowledge-mining-book/021-go.html

#load data
library(here)

avena_CpG_GUvsUU <- readRDS(here("avena_ncbi","DMPs_CpG_GUvsUU.rds"))
avena_CpG_UUvsUC <- readRDS(here("avena_ncbi","DMPs_CpG_UUvsUC.rds"))
avena_CpG_GUvsGC <- readRDS(here("avena_ncbi","DMPs_CpG_GUvsGC.rds"))

avena_CHG_GUvsUU <- readRDS(here("avena_ncbi","DMPs_CHG_GUvsUU.rds"))
avena_CHG_UUvsUC <- readRDS(here("avena_ncbi","DMPs_CHG_UUvsUC.rds"))
avena_CHG_GUvsGC <- readRDS(here("avena_ncbi","DMPs_CHG_GUvsGC.rds"))

avena_CHH_GUvsUU <- readRDS(here("avena_ncbi","DMPs_CHH_GUvsUU.rds"))
avena_CHH_UUvsUC <- readRDS(here("avena_ncbi","DMPs_CHH_UUvsUC.rds"))
avena_CHH_GUvsGC <- readRDS(here("avena_ncbi","DMPs_CHH_GUvsGC.rds"))

avena_all_GUvsUU <- read.table(here("avena_ncbi/DMPs","DMPs_all_GUvsUU.txt"), header = T)
avena_all_UUvsUC <- read.table(here("avena_ncbi/DMPs","DMPs_all_UUvsUC.txt"), header = T)
avena_all_GUvsGC <- read.table(here("avena_ncbi/DMPs","DMPs_all_GUvsGC.txt"), header = T)



library(dplyr)
library(purrr)

# 1. Crear una lista con todos los sitios evaluados (background)
lista_objetos <- list(
  avena_CpG_GUvsGC$all, avena_CpG_GUvsUU$all, avena_CpG_UUvsUC$all,
  avena_CHG_GUvsGC$all, avena_CHG_GUvsUU$all, avena_CHG_UUvsUC$all,
  avena_CHH_GUvsGC$all, avena_CHH_GUvsUU$all, avena_CHH_UUvsUC$all
)

# 2. Convertir a data frame, seleccionar 'chr' y 'start', unir y eliminar filas redundantes completas
resultado_unico <- lista_objetos %>%
  map_dfr(~ getData(.x) %>% select(chr, start)) %>%
  distinct(chr, start, .keep_all = TRUE)

# 3. Mirar si todavía quedan filas repetidas analizando SOLO la columna 'start'
repetidos_solo_start <- resultado_unico %>%
  filter(duplicated(start) | duplicated(start, fromLast = TRUE)) %>%
  arrange(start)

background_Cs <- resultado_unico$start
sig_MCs_GUvsUU <- avena_all_GUvsUU$start
sig_MCs_UUvsUC <- avena_all_UUvsUC$start 
sig_MCs_GUvsGC <- avena_all_GUvsGC$start


# Overrepresented biological processes, molecular functions, and cellular 
# components were identified with an FDR-adjusted threshold of alpha-value 
# < 0.05 (Benjamini & Hochberg, 1995).

enricher( gene, # a vector of gene id
          pvalueCutoff = 0.05, 
          pAdjustMethod = "BH", # "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
          universe = NULL, # background genes. If missing, the all genes listed in the database
          minGSSize = 10, # minimal size of genes annotated for testing
          maxGSSize = 500, # maximal size of genes annotated for testing
          qvalueCutoff = 0.2, # qvalue cutoff on enrichment tests to report as significant.
          gson = NULL, # annotation data
          TERM2GENE, # Only used when gson is NULL
          TERM2NAME = NA ) # Only used when gson is NULL



