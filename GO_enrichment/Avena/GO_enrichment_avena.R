# En este sscript voy a hacer pruebas para el GO enrichment de Avena
# web para utilizar el paquete: https://yulab-smu.top/biomedical-knowledge-mining-book/021-go.html


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")
# BiocManager::install("genomation")

library(here) # reproducibility
library(methylKit)
library(clusterProfiler) # enrichment (function enricher)
library(dplyr) # tidy
library(tidyr) # tidy
library(purrr) # tidy
library(genomation) # anotación de promotores y visualización de features
library(rtracklayer) # importar gff
library(ggplot2)


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
<<<<<<< HEAD

# ESTOS PASOS LLEVAN BASTANTE TIEMPO POR ESO LOS HE DEJADO EN COMENTARIO Y HE CREADO OBJETOS
# CON LOS QUE TRABAJAMOS DIRECTAMENTE

gene.obj <- readTranscriptFeatures(here("GO_enrichment/Avena/gff3.bed"),remove.unusual=FALSE, # También probé con longest.bed
                                   up.flank=2000,down.flank=0,unique.prom=TRUE)
=======
>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910

# ESTOS PASOS LLEVAN BASTANTE TIEMPO POR ESO LOS HE DEJADO EN COMENTARIO Y HE CREADO OBJETOS
# CON LOS QUE TRABAJAMOS DIRECTAMENTE

<<<<<<< HEAD
saveRDS(GUvsUU_annotation, file = "GO_enrichment/Avena/GUvsUU_annotation.rds")
saveRDS(UUvsUC_annotation, file = "GO_enrichment/Avena/UUvsUC_annotation.rds")
saveRDS(GUvsGC_annotation, file = "GO_enrichment/Avena/GUvsGC_annotation.rds")
saveRDS(background_annotation, file = "GO_enrichment/Avena/background_annotation.rds")
=======
# gene.obj <- readTranscriptFeatures(here("GO_enrichment/Avena/gff3.bed"),remove.unusual=FALSE, # También probé con longest.bed
#                                    up.flank=2000,down.flank=0,unique.prom=TRUE)
# 
# GUvsUU_annotation <- annotateWithGeneParts(as(avena_all_GUvsUU, "GRanges"), gene.obj)
# UUvsUC_annotation <- annotateWithGeneParts(as(avena_all_UUvsUC, "GRanges"), gene.obj)
# GUvsGC_annotation <- annotateWithGeneParts(as(avena_all_GUvsGC, "GRanges"), gene.obj)
# background_annotation <- annotateWithGeneParts(as(background_Cs, "GRanges"), gene.obj)
# 
# # Corregimos feature.name para que sea el nombre del TRANSCRITO (TAMBIÉN LO HICE CON EL DEL GEN)
# # GUvsUU_annotation@dist.to.TSS$feature.name <- sub(
# #   "^AVBAR\\.10000a\\.r1\\..*G([0-9]+)\\.\\d+$",         # eSTA VERSIÓN SERÍA PARA QUITAR EL NÚMERO DE TRANSCRITO
# #   "AVBAR.10000a.r1.1AG\\1",
# #   GUvsUU_annotation@dist.to.TSS$feature.name
# # )
# GUvsUU_annotation@dist.to.TSS$feature.name <- sub(
#   "^AVBAR\\.10000a\\.r1\\..*G([0-9]+)(\\.\\d+)$",
#   "AVBAR.10000a.r1.1AG\\1\\2",
#   GUvsUU_annotation@dist.to.TSS$feature.name
# )
# UUvsUC_annotation@dist.to.TSS$feature.name <- sub(
#   "^AVBAR\\.10000a\\.r1\\..*G([0-9]+)(\\.\\d+)$",
#   "AVBAR.10000a.r1.1AG\\1\\2",
#   UUvsUC_annotation@dist.to.TSS$feature.name
# )
# GUvsGC_annotation@dist.to.TSS$feature.name <- sub(
#   "^AVBAR\\.10000a\\.r1\\..*G([0-9]+)(\\.\\d+)$",
#   "AVBAR.10000a.r1.1AG\\1\\2",
#   GUvsGC_annotation@dist.to.TSS$feature.name
# )
# background_annotation@dist.to.TSS$feature.name <- sub(
#   "^AVBAR\\.10000a\\.r1\\..*G([0-9]+)(\\.\\d+)$",
#   "AVBAR.10000a.r1.1AG\\1\\2",
#   background_annotation@dist.to.TSS$feature.name
# )
# 
# saveRDS(GUvsUU_annotation, file = "GO_enrichment/Avena/GUvsUU_annotation.rds")
# saveRDS(UUvsUC_annotation, file = "GO_enrichment/Avena/UUvsUC_annotation.rds")
# saveRDS(GUvsGC_annotation, file = "GO_enrichment/Avena/GUvsGC_annotation.rds")
# saveRDS(background_annotation, file = "GO_enrichment/Avena/background_annotation.rds")
>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910

GUvsUU_annotation <- readRDS(here("GO_enrichment/Avena", "GUvsUU_annotation.rds"))
UUvsUC_annotation <- readRDS(here("GO_enrichment/Avena", "UUvsUC_annotation.rds"))
GUvsGC_annotation <- readRDS(here("GO_enrichment/Avena", "GUvsGC_annotation.rds"))
background_annotation <- readRDS(here("GO_enrichment/Avena", "background_annotation.rds"))



# Function to create annotation table and calculate feature proportions
summarise_annotation <- function(annotation, name) {
  
  anno <- annotation@members %>%
    as.data.frame() %>%
    mutate(target.row = seq_len(nrow(annotation@members))) %>%
    merge(annotation@dist.to.TSS, by = "target.row")
  
  # Save annotation table
  write.csv(
    anno,
    here("GO_enrichment", "Avena", paste0(name, "_annotation_TRANSCRIPT.csv")),
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


# FIGURA DE ANOTACIÓN ----



# 1. Construir el data frame con los porcentajes de cada categoría
#    Intergenic = 100 - (Promoter + Exon + Intron)
prop_df <- data.frame(
  Condition = rep(c("GUvsUU", "UUvsUC", "GUvsGC", "background"), each = 4),
  Feature   = rep(c("Promoter", "Exon", "Intron", "Intergenic"), times = 4),
  Percentage = c(
    # GUvsUU
    3.101628, 49.406071, 4.377475, 100 - (3.101628 + 49.406071 + 4.377475),
    # UUvsUC
    2.803738, 49.532710, 1.869159, 100 - (2.803738 + 49.532710 + 1.869159),
    # GUvsGC
    3.370787, 40.449438, 1.872659, 100 - (3.370787 + 40.449438 + 1.872659),
    # background
    2.498747, 53.416979, 3.803564, 100 - (2.498747 + 53.416979 + 3.803564)
  )
)

# 2. Fijar orden de condiciones y features
prop_df$Condition <- factor(prop_df$Condition,
                            levels = c("GUvsUU", "UUvsUC", "GUvsGC", "background"))
prop_df$Feature <- factor(prop_df$Feature,
                          levels = c("Promoter", "Exon", "Intron", "Intergenic"))

# 3. Paleta de colores bonitos (paleta tipo Nature / pastel elegante)
nice_colors <- c(
  "Promoter"   = "#E64B35",  # rojo coral
  "Exon"       = "#4DBBD5",  # azul cielo
  "Intron"     = "#00A087",  # verde azulado
  "Intergenic" = "#F39B7F"   # salmón suave
)

# Alternativa: paleta pastel
# nice_colors <- c("Promoter"="#F28E8E", "Exon"="#7EB6E0",
#                  "Intron"="#7DC9A9", "Intergenic"="#D3D3D3")

# 4. Figura de barras apiladas al 100%
p <- ggplot(prop_df, aes(x = Condition, y = Percentage, fill = Feature)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)),
            position = position_stack(vjust = 0.5),
            size = 3.2, color = "white", fontface = "bold") +
  scale_fill_manual(values = nice_colors) +
  scale_y_continuous(expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  labs(
    x = NULL,
    y = "Percentage of cytosines (%)",
    fill = "Feature",
    title = "Distribution of cytosines across genomic features"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(face = "bold", color = "grey20"),
    legend.position    = "right",
    plot.title         = element_text(face = "bold", hjust = 0.5)
  )

print(p)

# 5. Guardar
ggsave(
  here::here("GO_enrichment", "Avena", "feature_proportions_barplot.pdf"),
  plot = p, width = 8, height = 5, dpi = 300
)


# ENRICHMENT ----

# Como lleva tiempo, hemos guardado los resultados en csv y los cargamos directamente
# Para la figura

# Creamos TERM2GENE seleccionando las columnas gID y GO

TERM2GENE <- GOfeature %>%
  select(tID, GO) %>%                                 # sería gID para genes
  filter(!is.na(GO), GO != "", GO != "nan") %>%
  separate_rows(GO, sep = ",") %>%
  filter(grepl("^GO:\\d+$", GO)) %>%
  select(term = GO, gene = tID)                     # sería gID para genes

# Overrepresented biological processes, molecular functions, and cellular
# components were identified with an FDR-adjusted threshold of alpha-value
# < 0.05 (Benjamini & Hochberg, 1995).

run_go_enrichment <- function(gene_set) {

  enricher(
    gene = unique(gene_set),
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    universe = background_anno$feature.name,
    qvalueCutoff = 1,
    TERM2GENE = TERM2GENE
  )
}


GUvsUU_GO <- run_go_enrichment(
  GUvsUU_anno$feature.name)

UUvsUC_GO <- run_go_enrichment(
  UUvsUC_anno$feature.name)

GUvsGC_GO <- run_go_enrichment(
  GUvsGC_anno$feature.name)

as.data.frame(GUvsUU_GO)
as.data.frame(UUvsUC_GO)
as.data.frame(GUvsGC_GO)

write.table(GUvsUU_GO, file = "GO_enrichment/Avena/GUvsUU_GO.csv", sep = ";", row.names = T)
write.table(UUvsUC_GO, file = "GO_enrichment/Avena/UUvsUC_GO.csv", sep = ";", row.names = T)
write.table(GUvsGC_GO, file = "GO_enrichment/Avena/GUvsGC_GO.csv", sep = ";", row.names = T)
<<<<<<< HEAD


# FIGURAS ----

GUvsUU_GO <- read.table("GO_enrichment/Avena/GUvsUU_GO.csv", sep = ";", header = T)
UUvsUC_GO <- read.table("GO_enrichment/Avena/UUvsUC_GO.csv", sep = ";", header = T)
GUvsGC_GO <- read.table("GO_enrichment/Avena/GUvsGC_GO.csv", sep = ";", header = T)


# Asociar cada GOid con una descripción

library(GO.db)
library(AnnotationDbi)

# GO ID -> GO term and ontology
GO_terms <- AnnotationDbi::select(
  GO.db,
  keys = unique(c(
    UUvsUC_GO$ID,
    GUvsGC_GO$ID,
    GUvsUU_GO$ID
  )),
  keytype = "GOID",
  columns = c("GOID", "TERM", "ONTOLOGY")
) %>%
  distinct(GOID, .keep_all = TRUE)


# Function to prepare GO enrichment results
prepare_GO_fig <- function(df, label) {
  
  df %>%
    separate_wider_delim(
      GeneRatio,
      delim = "/",
      names = c("Nr Test", "Total Test")
    ) %>%
    separate_wider_delim(
      BgRatio,
      delim = "/",
      names = c("Nr Reference", "Total Reference")
    ) %>%
    mutate(
      `Nr Test` = as.numeric(`Nr Test`),
      `Total Test` = as.numeric(`Total Test`),
      `Nr Reference` = as.numeric(`Nr Reference`),
      `Total Reference` = as.numeric(`Total Reference`),
      
      `Not Annot Test` = `Total Test` - `Nr Test`,
      `Not Annot Ref` = `Total Reference` - `Nr Reference`,
      
      Tag = "OVER",
      `GO Term` = ID,
      `Adj. P-value` = p.adjust,
      `P-value` = pvalue
    ) %>%
    left_join(
      GO_terms,
      by = c("GO Term" = "GOID")
    ) %>%
    rename(
      `GO Name` = TERM,
      `GO Category` = ONTOLOGY
    ) %>%
    mutate(
      logAdjP = -log10(`Adj. P-value`),
      File = label
    ) %>%
    dplyr::select(
      Tag,
      `GO Term`,
      `GO Name`,
      `GO Category`,
      `Adj. P-value`,
      `P-value`,
      `Nr Test`,
      `Nr Reference`,
      `Not Annot Test`,
      `Not Annot Ref`,
      logAdjP,
      File
    )
}


# Prepare the three comparisons
UUvsUC_fig <- prepare_GO_fig(UUvsUC_GO, "UUvsUC")
GUvsGC_fig <- prepare_GO_fig(GUvsGC_GO, "GUvsGC")
GUvsUU_fig <- prepare_GO_fig(GUvsUU_GO, "GUvsUU")
=======
>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910


# FIGURAS ----

<<<<<<< HEAD
=======
GUvsUU_GO <- read.table("GO_enrichment/Avena/GUvsUU_GO.csv", sep = ";", header = T)
UUvsUC_GO <- read.table("GO_enrichment/Avena/UUvsUC_GO.csv", sep = ";", header = T)
GUvsGC_GO <- read.table("GO_enrichment/Avena/GUvsGC_GO.csv", sep = ";", header = T)


# Asociar cada GOid con una descripción

library(GO.db)
library(AnnotationDbi)

# GO ID -> GO term and ontology
GO_terms <- AnnotationDbi::select(
  GO.db,
  keys = unique(c(
    UUvsUC_GO$ID,
    GUvsGC_GO$ID,
    GUvsUU_GO$ID
  )),
  keytype = "GOID",
  columns = c("GOID", "TERM", "ONTOLOGY")
) %>%
  distinct(GOID, .keep_all = TRUE)


# Function to prepare GO enrichment results
prepare_GO_fig <- function(df, label) {
  
  df %>%
    separate_wider_delim(
      GeneRatio,
      delim = "/",
      names = c("Nr Test", "Total Test")
    ) %>%
    separate_wider_delim(
      BgRatio,
      delim = "/",
      names = c("Nr Reference", "Total Reference")
    ) %>%
    mutate(
      `Nr Test` = as.numeric(`Nr Test`),
      `Total Test` = as.numeric(`Total Test`),
      `Nr Reference` = as.numeric(`Nr Reference`),
      `Total Reference` = as.numeric(`Total Reference`),
      
      `Not Annot Test` = `Total Test` - `Nr Test`,
      `Not Annot Ref` = `Total Reference` - `Nr Reference`,
      
      Tag = "OVER",
      `GO Term` = ID,
      `Adj. P-value` = p.adjust,
      `P-value` = pvalue
    ) %>%
    left_join(
      GO_terms,
      by = c("GO Term" = "GOID")
    ) %>%
    rename(
      `GO Name` = TERM,
      `GO Category` = ONTOLOGY
    ) %>%
    mutate(
      logAdjP = -log10(`Adj. P-value`),
      File = label
    ) %>%
    dplyr::select(
      Tag,
      `GO Term`,
      `GO Name`,
      `GO Category`,
      `Adj. P-value`,
      `P-value`,
      `Nr Test`,
      `Nr Reference`,
      `Not Annot Test`,
      `Not Annot Ref`,
      logAdjP,
      File
    )
}


# Prepare the three comparisons
UUvsUC_fig <- prepare_GO_fig(UUvsUC_GO, "UUvsUC")
GUvsGC_fig <- prepare_GO_fig(GUvsGC_GO, "GUvsGC")
GUvsUU_fig <- prepare_GO_fig(GUvsUU_GO, "GUvsUU")


>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910
# -----------------------------
# Prepare combined data
# -----------------------------
go_data <- list(
  GUvsUU = GUvsUU_fig,
  UUvsUC = UUvsUC_fig,
  GUvsGC = GUvsGC_fig
)


# -----------------------------
# Common X-axis limits
# -----------------------------
<<<<<<< HEAD
=======
bubble_x_max <- max(
  c(
    GUvsUU_fig$logAdjP,
    UUvsUC_fig$logAdjP,
    GUvsGC_fig$logAdjP
  ),
  na.rm = TRUE
)

bubble_x_max <- ceiling(bubble_x_max)


>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910
bar_x_max <- max(
  c(
    100 * GUvsUU_fig$`Nr Test` /
      (GUvsUU_fig$`Nr Test` + GUvsUU_fig$`Not Annot Test`),
    100 * UUvsUC_fig$`Nr Test` /
      (UUvsUC_fig$`Nr Test` + UUvsUC_fig$`Not Annot Test`),
    100 * GUvsGC_fig$`Nr Test` /
      (GUvsGC_fig$`Nr Test` + GUvsGC_fig$`Not Annot Test`),
    100 * GUvsUU_fig$`Nr Reference` /
      (GUvsUU_fig$`Nr Reference` + GUvsUU_fig$`Not Annot Ref`),
    100 * UUvsUC_fig$`Nr Reference` /
      (UUvsUC_fig$`Nr Reference` + UUvsUC_fig$`Not Annot Ref`),
    100 * GUvsGC_fig$`Nr Reference` /
      (GUvsGC_fig$`Nr Reference` + GUvsGC_fig$`Not Annot Ref`)
  ),
  na.rm = TRUE
)

bar_x_max <- ceiling(bar_x_max / 5) * 5


# -----------------------------
<<<<<<< HEAD
=======
# Common size scale for Nr Test
# Same scale as the original Plantago figure
# -----------------------------
size_limits <- c(0, 60)


# -----------------------------
# Bubble plot
# -----------------------------
plot_bubble <- function(df, top_n_terms = 14) {
  
  top_terms <- df %>%
    arrange(desc(logAdjP)) %>%
    slice_head(n = top_n_terms)
  
  top_terms$`GO Name` <- factor(
    top_terms$`GO Name`,
    levels = rev(unique(top_terms$`GO Name`))
  )
  
  ggplot(
    top_terms,
    aes(
      x = logAdjP,
      y = `GO Name`,
      size = `Nr Test`
    )
  ) +
    geom_point(
      alpha = 0.8,
      color = "#CD69C9"
    ) +
    scale_size(
      limits = size_limits,
      range = c(2, 10),
      breaks = c(0, 10, 20, 30, 40, 50, 60),
      name = "Nr Test"
    ) +
    scale_x_continuous(
      limits = c(0, bubble_x_max)
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    xlab(expression(-log[10](p)))
}


# -----------------------------
>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910
# Two-bar percentage plot
# -----------------------------
plot_two_bar <- function(df, top_n_terms = 14) {
  
<<<<<<< HEAD
  # Keep only GO terms with valid names
  # and select the most significant terms
  top_terms <- df %>%
    filter(
      !is.na(`GO Name`),
      `GO Name` != ""
    ) %>%
=======
  top_terms <- df %>%
>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910
    arrange(desc(logAdjP)) %>%
    slice_head(n = top_n_terms) %>%
    mutate(
      pct_Test = 100 * (
        `Nr Test` /
          (`Nr Test` + `Not Annot Test`)
      ),
      pct_Ref = 100 * (
        `Nr Reference` /
          (`Nr Reference` + `Not Annot Ref`)
<<<<<<< HEAD
      ),
      p_adjust = 10^(-logAdjP),
      significance = case_when(
        p_adjust < 0.001 ~ "***",
        p_adjust < 0.01  ~ "**",
        p_adjust < 0.05  ~ "*",
        p_adjust < 0.1   ~ "·",
        TRUE ~ "ns"
      )
    )
  
  
  # Convert to long format for the two bars
=======
      )
    )
  
>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910
  long <- top_terms %>%
    dplyr::select(
      `GO Name`,
      pct_Test,
      pct_Ref
    ) %>%
    pivot_longer(
      cols = c(pct_Ref, pct_Test),
      names_to = "Series",
      values_to = "% of sequences"
    ) %>%
    mutate(
      Series = factor(
        Series,
        levels = c("pct_Ref", "pct_Test")
      )
    )
  
<<<<<<< HEAD
  
  # Most significant term at the top
  long$`GO Name` <- factor(
    long$`GO Name`,
    levels = rev(top_terms$`GO Name`)
  )
  
  
  # Data frame for significance symbols
  significance_df <- top_terms %>%
    mutate(
      x = pmax(pct_Test, pct_Ref) + 2
    )
  
  significance_df$`GO Name` <- factor(
    significance_df$`GO Name`,
    levels = rev(top_terms$`GO Name`)
  )
  
  
  # Plot
=======
  long$`GO Name` <- factor(
    long$`GO Name`,
    levels = rev(unique(top_terms$`GO Name`))
  )
  
>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910
  ggplot(
    long,
    aes(
      x = `% of sequences`,
      y = `GO Name`,
      fill = Series
    )
  ) +
<<<<<<< HEAD
    geom_col(
      position = "dodge"
    ) +
    
    geom_text(
      data = significance_df,
      aes(
        x = x,
        y = `GO Name`,
        label = significance
      ),
      inherit.aes = FALSE,
      hjust = 0,
      size = 5
    ) +
    
=======
    geom_col(position = "dodge") +
>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910
    scale_fill_manual(
      values = c(
        "pct_Ref" = "#008B45",
        "pct_Test" = "#CD69C9"
      ),
      labels = c(
        "Reference",
        "Test"
      )
    ) +
<<<<<<< HEAD
    
    scale_x_continuous(
      limits = c(0, bar_x_max + 8)
    ) +
    
    theme_bw(
      base_size = 14
    ) +
    
=======
    scale_x_continuous(
      limits = c(0, bar_x_max)
    ) +
    theme_bw(base_size = 14) +
>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 12)
    ) +
<<<<<<< HEAD
    
=======
>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910
    xlab("% of sequences")
}


<<<<<<< HEAD

# -----------------------------
# Generate plots
# -----------------------------
plots_GUvsUU <- plot_two_bar(
  GUvsUU_fig
)

plots_UUvsUC <- plot_two_bar(
  UUvsUC_fig
)

plots_GUvsGC <- plot_two_bar(
  GUvsGC_fig
=======
# -----------------------------
# Generate plots
# -----------------------------
plots_GUvsUU <- list(
  bubble = plot_bubble(GUvsUU_fig),
  twobar = plot_two_bar(GUvsUU_fig)
)

plots_UUvsUC <- list(
  bubble = plot_bubble(UUvsUC_fig),
  twobar = plot_two_bar(UUvsUC_fig)
)

plots_GUvsGC <- list(
  bubble = plot_bubble(GUvsGC_fig),
  twobar = plot_two_bar(GUvsGC_fig)
)


# -----------------------------
# Prepare plots for combination
# -----------------------------
prep_for_combination <- function(bar_plot, bubble_plot) {
  
  bubble_clean <- bubble_plot +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  list(
    bars = bar_plot,
    bubble = bubble_clean
  )
}


GUvsUU_clean <- prep_for_combination(
  plots_GUvsUU$twobar,
  plots_GUvsUU$bubble
)

UUvsUC_clean <- prep_for_combination(
  plots_UUvsUC$twobar,
  plots_UUvsUC$bubble
)

GUvsGC_clean <- prep_for_combination(
  plots_GUvsGC$twobar,
  plots_GUvsGC$bubble
>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910
)


# -----------------------------
# Side labels
# -----------------------------
row_title_right <- function(text) {
  
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = text,
      angle = 270,
      size = 6
    ) +
    theme_void()
}


GUvsUU_label <- row_title_right("GUvsUU")
UUvsUC_label <- row_title_right("UUvsUC")
GUvsGC_label <- row_title_right("GUvsGC")


# -----------------------------
# Build rows
# -----------------------------
<<<<<<< HEAD
row1 <- plots_GUvsUU +
  GUvsUU_label +
  plot_layout(
    widths = c(1, 0.15)
=======
row1 <- GUvsUU_clean$bars +
  GUvsUU_clean$bubble +
  GUvsUU_label +
  plot_layout(
    widths = c(1, 1, 0.15)
>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910
  ) &
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


<<<<<<< HEAD
row2 <- plots_UUvsUC +
  UUvsUC_label +
  plot_layout(
    widths = c(1, 0.15)
=======
row2 <- UUvsUC_clean$bars +
  UUvsUC_clean$bubble +
  UUvsUC_label +
  plot_layout(
    widths = c(1, 1, 0.15)
>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910
  ) &
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


<<<<<<< HEAD
row3 <- plots_GUvsGC +
  GUvsGC_label +
  plot_layout(
    widths = c(1, 0.15)
=======
row3 <- GUvsGC_clean$bars +
  GUvsGC_clean$bubble +
  GUvsGC_label +
  plot_layout(
    widths = c(1, 1, 0.15)
>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910
  )


# -----------------------------
# Combine rows
# -----------------------------
final_fig <- row1 / row2 / row3 +
  plot_layout(
    heights = c(1, 1, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )


# -----------------------------
# Show X-axis only in bottom row
# -----------------------------
final_fig[[3]][[1]] <- final_fig[[3]][[1]] +
  theme(
    axis.title.x = element_text(),
    axis.text.x = element_text(),
    axis.ticks.x = element_line()
  )

<<<<<<< HEAD
=======
final_fig[[3]][[2]] <- final_fig[[3]][[2]] +
  theme(
    axis.title.x = element_text(),
    axis.text.x = element_text(),
    axis.ticks.x = element_line()
  )

>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910

# -----------------------------
# Display figure
# -----------------------------
<<<<<<< HEAD
final_fig
=======
final_fig
>>>>>>> 24c812fa5f44b6ab2ab1aca934400a5cb5014910
