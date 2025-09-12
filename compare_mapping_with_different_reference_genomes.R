# Imports multiple Bismark single-end summary reports, extracts key alignment and methylation metrics,
# combines them into a single dataframe with mean ± standard error, creates a metadata table explaining each variable,
# and generates summary plots showing distribution of values by species and reference genome.


# Librerías
library(tidyverse)
library(here)

# Directorio de archivos
input_dir <- here("SE_reports")

# Listar archivos
files <- list.files(input_dir, pattern = "bismark_summary_.*\\.txt", full.names = TRUE)

# Función para procesar un archivo
read_summary_alt <- function(file) {
  # Leer todas las líneas
  lines <- readLines(file)
  
  # Filtrar solo líneas que contienen ":"
  lines_colon <- lines[str_detect(lines, ":")]
  
  # Leer como tabla de dos columnas, sin interpretar comillas
  dat <- read.table(
    text = lines_colon,
    sep = ":",
    fill = TRUE,
    stringsAsFactors = FALSE,
    strip.white = TRUE,
    quote = ""   # <--- Desactiva el parsing de comillas
  )
  colnames(dat) <- c("variable","value")
  
  # Normalizar espacios Unicode y extraer media ± SE
  dat <- dat %>%
    mutate(
      value = str_replace_all(value, "\u00A0", " "),
      value = str_trim(value),
      mean = as.numeric(str_replace(str_extract(value, "^[0-9]+[,.][0-9]+"), ",", ".")),
      se   = as.numeric(str_replace(str_extract(value, "(?<=± )[0-9]+[,.][0-9]+"), ",", "."))
    )
  
  # Extraer especie y referencia del nombre del archivo
  fname <- basename(file)
  dat$sp <- str_extract(fname, "avena|hordeum|plantago")
  dat$ref <- str_extract(fname, "ncbi|andres")
  
  return(dat)
}


# Leer todos los archivos y unir
all_data_long <- map_dfr(files, read_summary_alt)

# Pivotar a formato ancho: mean_XXX y se_XXX
all_data <- all_data_long %>%
  select(sp, ref, variable, mean, se) %>%
  pivot_wider(names_from = variable, values_from = c(mean,se))

# Guardar dataframe
write_csv(all_data, file.path(input_dir, "summary_dataframe.csv"))

# Metadata de variables
metadata <- all_data_long %>%
  select(variable) %>%
  distinct() %>%
  mutate(
    description = c(
      "Total sequences analysed",
      "Unique best alignments",
      "Mapping efficiency",
      "Sequences without alignments",
      "Non-unique mappings",
      "Sequences discarded (genomic sequence could not be extracted)",
      "Number of sequences with unique best (first) alignment from Bowtie output",
      "Converted top strand (CT/CT)",
      "Converted bottom strand (CT/GA)",
      "Complementary to converted top strand (GA/CT)",
      "Complementary to converted bottom strand (GA/GA)",
      "Total cytosines analysed",
      "Methylated C in CpG",
      "Methylated C in CHG",
      "Methylated C in CHH",
      "Methylated C in unknown context",
      "Unmethylated C in CpG",
      "Unmethylated C in CHG",
      "Unmethylated C in CHH",
      "Unmethylated C in unknown context",
      "C methylation in CpG",
      "C methylation in CHG",
      "C methylation in CHH",
      "C methylation in unknown context"
    )[1:nrow(.)],
    units = c(
      "count",  # Sequences analysed in total
      "count",  # Unique best alignments
      "%",      # Mapping efficiency
      "count",  # Sequences without alignments
      "count",  # Non-unique mappings
      "count",  # Sequences discarded
      "count",  # Unique best first alignment (Bowtie)
      "count",  # CT/CT
      "count",  # CT/GA
      "count",  # GA/CT
      "count",  # GA/GA
      "count",  # Total cytosines analysed
      "count",  # Methylated C in CpG
      "count",  # Methylated C in CHG
      "count",  # Methylated C in CHH
      "count",  # Methylated C in unknown
      "count",  # Unmethylated C in CpG
      "count",  # Unmethylated C in CHG
      "count",  # Unmethylated C in CHH
      "count",  # Unmethylated C in unknown
      "%",      # C methylation in CpG
      "%",      # C methylation in CHG
      "%",      # C methylation in CHH
      "%"       # C methylation in unknown
    )[1:nrow(.)],
    method_SE = "Standard error calculated across reads or alignments per Bismark summary"
  )



write_csv(metadata, file.path(input_dir, "metadata_variables.csv"))

# ---- GRÁFICOS ----
plot_dir <- file.path(input_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE)

# Graficar cada variable (media ± SE) por especie y referencia
for (v in unique(all_data_long$variable)) {
  mean_col <- paste0("mean_", v)
  se_col <- paste0("se_", v)
  if (!(mean_col %in% names(all_data))) next
  
  df_plot <- all_data %>%
    select(sp, ref, !!mean_col, !!se_col) %>%
    rename(mean = !!mean_col, se = !!se_col)
  
  p <- ggplot(df_plot, aes(x = sp, y = mean, fill = ref)) +
    geom_col(position = position_dodge(width = 0.8)) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(width = 0.8), width = 0.3) +
    labs(title = v,
         y = paste0(v," (", metadata$units[metadata$variable==v],")"),
         x = "Species") +
    theme_minimal(base_size = 14) +
    scale_fill_brewer(palette = "Set2")
  
  ggsave(filename = file.path(plot_dir, paste0(v,".png")),
         plot = p, width = 7, height = 5)
}
