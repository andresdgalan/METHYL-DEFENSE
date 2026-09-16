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
  
  data <- read_tsv(
    file_path,
    show_col_types = FALSE
  )
  
  data <- data %>%
    mutate(
      Tag = ifelse(Tag == "", "NOT_SIG", Tag),
      logAdjP = -log10(`Adj. P-value`),
      File = label
    ) %>%
    filter(Tag %in% c("OVER", "UNDER"))
  
  return(data)
}


# -----------------------------
# File paths
# -----------------------------
go_dir <- here("GO_enrichment/Avena")


# -----------------------------
# Select comparison and context
# -----------------------------
comparison_GUvsUU <- "GUvsUU"
comparison_UUvsUC <- "UUvsUC"

context <- "COMBINED"


file_GUvsUU <- file.path(
  go_dir,
  paste0(
    "results_",
    context,
    "_",
    comparison_GUvsUU,
    "_BP_ARGOT.txt"
  )
)

file_UUvsUC <- file.path(
  go_dir,
  paste0(
    "results_",
    context,
    "_",
    comparison_UUvsUC,
    "_BP_ARGOT.txt"
  )
)


data_GUvsUU <- prepare_go_data(
  file_GUvsUU,
  comparison_GUvsUU
)

data_UUvsUC <- prepare_go_data(
  file_UUvsUC,
  comparison_UUvsUC
)


# -----------------------------
# Bubble plot with unified size scale
# -----------------------------
plot_bubble <- function(
    df,
    top_n_terms = 14,
    size_range = c(0, 60)
) {
  
  top_terms <- df %>%
    filter(`Nr Test` > 0) %>%
    arrange(desc(logAdjP)) %>%
    slice(1:min(top_n_terms, n()))
  
  top_terms$`GO Name` <- factor(
    top_terms$`GO Name`,
    levels = rev(unique(top_terms$`GO Name`))
  )
  
  ggplot(
    top_terms,
    aes(
      x = logAdjP,
      y = `GO Name`,
      size = `Nr Test`,
      color = Tag
    )
  ) +
    
    geom_point(alpha = 0.8) +
    
    scale_color_manual(
      values = c(
        "OVER" = "#CD69C9",
        "UNDER" = "#008B45"
      )
    ) +
    
    scale_size(
      limits = size_range
    ) +
    
    theme_bw(
      base_size = 14
    ) +
    
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    
    xlab(
      expression(-log[10](q))
    ) +
    
    xlim(
      0,
      36
    )
}


# -----------------------------
# Two-bar percentage plot
# -----------------------------
plot_two_bar <- function(df, top_n_terms = 14) {
  
  top_terms <- df %>%
    filter(`Nr Test` > 0) %>%
    arrange(desc(logAdjP)) %>%
    slice(1:min(top_n_terms, n())) %>%
    mutate(
      pct_Test = 100 * (`Nr Test` / (`Nr Test` + `Not Annot Test`)),
      pct_Ref  = 100 * (`Nr Reference` / (`Nr Reference` + `Not Annot Ref`))
    )
  
  long <- top_terms %>%
    select(`GO Name`, pct_Test, pct_Ref) %>%
    pivot_longer(
      cols = c(pct_Test, pct_Ref),
      names_to = "Series",
      values_to = "% of sequences"
    ) %>%
    mutate(
      Series = factor(Series, levels = c("pct_Ref", "pct_Test"))
    )
  
  long$`GO Name` <- factor(
    long$`GO Name`,
    levels = rev(unique(top_terms$`GO Name`))
  )
  
  ggplot(long, aes(
    x = `% of sequences`,
    y = `GO Name`,
    fill = Series
  )) +
    geom_col(position = "dodge") +
    scale_fill_manual(
      values = c("pct_Ref" = "#008B45", "pct_Test" = "#CD69C9"),
      labels = c("Reference", "Test")
    ) +
    xlim(0, 30) +
    theme_bw(base_size = 14) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 12)
    )
}


# -----------------------------
# Generate plots
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
# Prepare plots for combination
# -----------------------------
prep_for_combination <- function(
    bar_plot,
    bubble_plot
) {
  
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


GU_clean <- prep_for_combination(
  plots_GUvsUU$twobar,
  plots_GUvsUU$bubble
)

UU_clean <- prep_for_combination(
  plots_UUvsUC$twobar,
  plots_UUvsUC$bubble
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


row1_label <- row_title_right(
  "Maternal Environment"
)

row2_label <- row_title_right(
  "Real-time Herbivory Simulation"
)


# -----------------------------
# Hide X axes in upper row
# -----------------------------
row1 <- GU_clean$bars +
  GU_clean$bubble +
  row1_label +
  
  plot_layout(
    widths = c(1, 1, 0.15)
  ) &
  
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )


row2 <- UU_clean$bars +
  UU_clean$bubble +
  row2_label +
  
  plot_layout(
    widths = c(1, 1, 0.15)
  )


# -----------------------------
# Combine final figure
# -----------------------------
final_fig <- row1 / row2 +
  
  plot_layout(
    heights = c(1, 1),
    guides = "collect"
  ) &
  
  theme(
    legend.position = "right"
  )


# -----------------------------
# Show X-axis labels only
# in the bottom row
# -----------------------------
final_fig[[2]][[1]] <-
  final_fig[[2]][[1]] +
  theme(
    axis.title.x = element_text(),
    axis.text.x = element_text()
  )

final_fig[[2]][[2]] <-
  final_fig[[2]][[2]] +
  theme(
    axis.title.x = element_text(),
    axis.text.x = element_text()
  )


# -----------------------------
# Display figure
# -----------------------------
final_fig
