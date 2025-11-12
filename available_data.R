library(tidyverse)
library(here)
library(ggtext)

# -----------------------------
# Cargar datos
# -----------------------------
df <- read_delim(
  here("MultiQC_reports", "read_counts.csv"),
  delim = ";"
)

# -----------------------------
# Filtrar y limpiar columnas
# -----------------------------
df_clean <- df %>%
  filter(!str_ends(sample, ".rem")) %>%      # eliminar *.rem
  select(-R2_reads, -`difference `) %>%         # eliminar columnas
  mutate(
    group = case_when(
      str_starts(sample, "HO") ~ "HO",
      str_starts(sample, "AV") ~ "AV",
      str_starts(sample, "PL") ~ "PL",
      TRUE ~ "Other"
    ),
    is_2 = if_else(str_ends(sample, "_2"), TRUE, FALSE)
  )

# -----------------------------
# Ordenar samples por grupo
# -----------------------------
df_clean <- df_clean %>%
  arrange(group, is_2, sample) %>%
  mutate(sample = factor(sample, levels = sample))

# -----------------------------
# Etiquetas del eje X con _2 en negrita
# -----------------------------
df_clean <- df_clean %>%
  mutate(
    sample_label = if_else(
      is_2,
      paste0("<b>", sample, "</b>"),
      sample
    )
  )

# -----------------------------
# Gráfico
# -----------------------------
ggplot(df_clean, aes(x = sample_label, y = R1_reads, fill = group)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(
    x = "sample",
    y = "nº reads before filtering",
    fill = "Grupo"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_markdown(size = 8),  # permitir negrita selectiva
    axis.title.x = element_text(margin = margin(t = 10))
  )
