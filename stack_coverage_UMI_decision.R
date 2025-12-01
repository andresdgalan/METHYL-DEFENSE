library(ggplot2)
library(dplyr)
library(scales)
library(readr)

# Cargar tabla
avena <- read_tsv("coverage_summary_avena.tsv")
plantago <- read_tsv("coverage_summary_plantago.tsv")
hordeum <- read_tsv("coverage_summary_hordeum.tsv")

species_list <- list(avena, plantago, hordeum)
species_names <- c("Avena", "Plantago", "Hordeum")

# Layout 3 filas x 2 columnas
par(mfrow = c(3, 2))

for (j in seq_along(species_list)) {
  
  df <- species_list[[j]]
  samples <- unique(df$sample)
  
  ## 1️⃣ Densidad log10
  plot(NULL, xlim = c(1, max(df$coverage_primary, na.rm = TRUE)),
       ylim = c(0, max(density(df$coverage_primary)$y)),
       log = "x",
       xlab = "Coverage",
       ylab = "Density",
       main = paste(species_names[j], "- Densidad log10"))
  
  cols <- scales::hue_pal()(length(samples))
  for (i in seq_along(samples)) {
    vals <- df$coverage_primary[df$sample == samples[i]]
    dens <- density(vals)
    lines(dens$x, dens$y, col = cols[i], lwd = 2)
  }
  legend("topright", legend = samples, col = cols, lwd = 2, cex = 0.8)
  
  ## 2️⃣ Curvas suavizadas de frecuencia
  breaks_seq <- seq(0, 70, by = 1)
  # calcular máximo
  max_count <- 0
  for (s in samples) {
    vals <- df$coverage_primary[
      df$sample == s &
        df$coverage_primary >= 0 &
        df$coverage_primary <= 70
    ]
    h <- hist(vals, breaks = breaks_seq, plot = FALSE)
    max_count <- max(max_count, max(h$counts, na.rm = TRUE))
  }
  y_max <- max_count + 10
  
  plot(NULL,
       xlim = c(0, 70),
       ylim = c(0, y_max),
       xlab = "Coverage",
       ylab = "Frecuencia",
       main = paste(species_names[j], "- Frecuencia suavizada")
  )
  
  for (i in seq_along(samples)) {
    vals <- df$coverage_primary[
      df$sample == samples[i] &
        df$coverage_primary >= 0 &
        df$coverage_primary <= 70
    ]
    h <- hist(vals, breaks = breaks_seq, plot = FALSE)
    smooth_line <- spline(h$mids, h$counts)
    lines(smooth_line$x, smooth_line$y,
          col = cols[i],
          lwd = 2)
  }
  
  legend("topright", legend = samples, col = cols, lwd = 2, cex = 0.8)
}
``


library(dplyr)

# Lista de especies y nombres
species_list <- list(avena = avena, plantago = plantago, hordeum = hordeum)

# Crear tabla vacía
summary_table <- data.frame()

# Loop sobre especies
for (sp in names(species_list)) {
  
  df <- species_list[[sp]]
  
  tmp <- df %>%
    group_by(sample) %>%
    summarise(
      n_stacks = n(),
      max_coverage = max(coverage_primary, na.rm = TRUE)
    ) %>%
    mutate(species = sp)  # añadir columna de especie
  
  summary_table <- bind_rows(summary_table, tmp)
}

# Reordenar columnas
summary_table <- summary_table %>%
  select(species, sample, n_stacks, max_coverage)

summary_table
