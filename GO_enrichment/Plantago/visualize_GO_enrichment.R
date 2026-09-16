# -----------------------------
# Load libraries
# -----------------------------
library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(forcats)
library(tidyr)
library(patchwork)


# -----------------------------
# Prepare GO data
# -----------------------------
prepare_go_data <- function(file_path, label) {
  data <- read_tsv(file_path)
  
  data$Tag <- ifelse(data$Tag == "", "NOT_SIG", data$Tag)
  
  data <- data %>%
    mutate(logAdjP = -log10(`Adj. P-value`),
           File = label) %>%
    filter(Tag %in% c("OVER", "UNDER"))
  
  return(data)
}

# -----------------------------
# File paths
# -----------------------------
file_GUvsUU <- here("GO_enrichment", "results_all_GUvsUU.txt")
file_UUvsUC <- here("GO_enrichment", "results_all_UUvsUC.txt")

data_GUvsUU <- prepare_go_data(file_GUvsUU, "GUvsUU")
data_UUvsUC <- prepare_go_data(file_UUvsUC, "UUvsUC")






# -----------------------------
# Bubble plot with unified size scale
# -----------------------------
plot_bubble <- function(df, top_n_terms = 14, size_range = c(0, 60)) {
  
  top_terms <- df %>%
    arrange(desc(logAdjP)) %>%
    slice(1:min(top_n_terms, n()))
  
  top_terms$`GO Name` <- factor(top_terms$`GO Name`, levels = rev(unique(top_terms$`GO Name`)))
  
  ggplot(top_terms, aes(
    x = logAdjP,
    y = `GO Name`,
    size = `Nr Test`,
    color = Tag
  )) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("OVER" = "#CD69C9", "UNDER" = "#008B45")) +
    scale_size(limits = size_range) +
    theme_bw(base_size = 14) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    xlab(expression(-log[10](p))) +
    xlim(0, 36)  # fija escala X para todas las comparaciones
}

# -----------------------------
# Two-bar percentage plot with unified X scale
# -----------------------------
plot_two_bar <- function(df, top_n_terms = 14) {
  
  top_terms <- df %>%
    arrange(desc(logAdjP)) %>%
    slice(1:min(top_n_terms, n())) %>%
    mutate(
      pct_Test = 100 * (`Nr Test` / (`Nr Test` + `Not Annot Test`)),
      pct_Ref  = 100 * (`Nr Reference` / (`Nr Reference` + `Not Annot Ref`))
    )
  
  long <- top_terms %>%
    select(`GO Name`, pct_Test, pct_Ref) %>%
    pivot_longer(cols = c(pct_Test, pct_Ref), names_to = "Series", values_to = "% of sequences") %>%
    mutate(
      Series = factor(Series, levels = c("pct_Ref","pct_Test"))  # Reference primero, Test segundo
    )
  
  long$`GO Name` <- factor(long$`GO Name`, levels = rev(unique(top_terms$`GO Name`)))
  
  ggplot(long, aes(
    x = `% of sequences`,
    y = `GO Name`,
    fill = Series
  )) +
    geom_col(position = "dodge") +
    scale_fill_manual(
      values = c("pct_Ref" = "#008B45", "pct_Test" = "#CD69C9"),
      labels = c("Reference","Test")
    ) +
    xlim(0, 30) +  # fija escala para todas las barras
    theme_bw(base_size = 14) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 12)
    )
  
}


# -----------------------------
# Generate plots for each comparison
# -----------------------------
plots_GUvsUU <- list(
  bubble = plot_bubble(data_GUvsUU),
  twobar = plot_two_bar(data_GUvsUU)
)

plots_UUvsUC <- list(
  bubble = plot_bubble(data_UUvsUC),
  twobar = plot_two_bar(data_UUvsUC)
)



# -----------------------------
# Combine plots
# -----------------------------
prep_for_combination <- function(bar_plot, bubble_plot) {
  
  # Bubble plot limpio
  bubble_clean <- bubble_plot +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  list(bars = bar_plot, bubble = bubble_clean)
}

GU_clean <- prep_for_combination(plots_GUvsUU$twobar, plots_GUvsUU$bubble)
UU_clean <- prep_for_combination(plots_UUvsUC$twobar, plots_UUvsUC$bubble)

# -----------------------------
# Side labels (comparisons)
# -----------------------------
row_title_right <- function(text) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = text, angle = 270, size = 6) +
    theme_void()
}

row1_label <- row_title_right("Maternal Environment")
row2_label <- row_title_right("Real-time Herbivory Simulation")

# -----------------------------
# Ocultar ejes X solo en fila superior
# -----------------------------
row1 <- GU_clean$bars + GU_clean$bubble + row1_label +
  plot_layout(widths = c(1,1,0.15)) &
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

row2 <- UU_clean$bars + UU_clean$bubble + row2_label +
  plot_layout(widths = c(1,1,0.15))  # aquí no tocamos los ejes X

# -----------------------------
# Combinar filas en final_fig
# -----------------------------
final_fig <- row1 / row2 +
  plot_layout(heights = c(1,1), guides = "collect") &
  theme(
    legend.position = "right"
  )

# Mostrar etiquetas X solo en fila inferior
final_fig[[2]][[1]] <- final_fig[[2]][[1]] + theme(axis.title.x = element_text(), axis.text.x = element_text())
final_fig[[2]][[2]] <- final_fig[[2]][[2]] + theme(axis.title.x = element_text(), axis.text.x = element_text())

# -----------------------------
# Mostrar figura
# -----------------------------
final_fig
