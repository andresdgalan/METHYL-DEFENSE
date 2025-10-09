# CARGAR AVENA CpG primero

# extract genomic positions
coords <- getData(CpG_meth)$start

# assign them as rownames to meth_matrix
rownames(CpG_meth_matrix) <- coords

grp <- as.factor(treatments)


# Run DAPC keeping 14 principal components (and all possible discriminant functions). I used 14 because it's a safe maximum: 42/3
dapc_res <- dapc(t(CpG_meth_matrix), grp, n.pca = 14, n.da = length(levels(grp)) - 1)

# Optimize the number of PCs using the alpha-score method (balances discrimination vs overfitting).
optim_a <- optim.a.score(dapc_res)

# Print the optimal number of PCs suggested by optim.a.score
cat("Optimal number of PCs (alpha-score):", optim_a$best, "\n")

# Re-run DAPC using the optimal number of PCs suggested automatically
dapc_res <- dapc(t(CpG_meth_matrix), grp, 
                 n.pca = optim_a$best, 
                 n.da = length(levels(grp)) - 1)

# Extract full loadings matrix (all CpGs × LDs)
loadings <- dapc_res$var.contr  

# Convert to tidy data frame: CpG, LD, Loading
loadings_df <- melt(loadings,
                    varnames = c("CpG", "LD"),
                    value.name = "Loading")

# cREATING A DATAFRAME WITH DAPC, TREATMENT AND SAMPLE ID
dapc_scores <- as.data.frame(dapc_res$ind.coord)
dapc_scores$id <- sample_ids
dapc_scores$combined <- grp

# Correct factor levels and labels (matching your desired order and names)
dapc_scores$combined <- factor(dapc_scores$combined,
                               levels = c("0", "2", "1", "3"),
                               labels = c("Grazed + Clipped",
                                          "Grazed + Unclipped",
                                          "Ungrazed + Clipped",
                                          "Ungrazed + Unclipped"))

# Compute group centroids
centroids <- dapc_scores %>%
  group_by(combined) %>%
  summarise(LD1 = mean(LD1), LD2 = mean(LD2))

# Plot
ggplot(dapc_scores, aes(x = LD1, y = LD2, color = combined)) +
  # Vectors from centroid to each point
  geom_segment(data = dapc_scores %>%
                 left_join(centroids, by = "combined"),
               aes(x = LD1.y, y = LD2.y, xend = LD1.x, yend = LD2.x, color = combined),
               alpha = 0.4, size = 0.6, inherit.aes = FALSE) +
  # Points
  geom_point(size = 1.8) +
  # Open inertia/confidence ellipses
  stat_ellipse(type = "norm", level = 0.67, size = 0.6) +
  # Colors for each group
  scale_color_manual(values = c(
    "Grazed + Clipped" = "#FF69B4",
    "Grazed + Unclipped" = "#8B3A62",
    "Ungrazed + Clipped" = "#00FF7F",
    "Ungrazed + Unclipped" = "#008B45"
  )) +
  # Legend title
  labs(title = "DAPC of CpG Methylation in A.barbata",
       x = "LD1", y = "LD2",
       color = "Herbivory treatment") +
  # Theme adjustments
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),   # remove grid
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1), # box around plot
    aspect.ratio = 1  # make the plot square
  )





# Filtrar solo los dos grupos de interés
dapc_scores_sub <- dapc_scores %>%
  filter(combined %in% c("Grazed + Unclipped", "Ungrazed + Unclipped"))

centroids_sub <- centroids %>%
  filter(combined %in% c("Grazed + Unclipped", "Ungrazed + Unclipped"))

# Plot solo dos grupos
ggplot(dapc_scores_sub, aes(x = LD1, y = LD2, color = combined)) +
  geom_segment(data = dapc_scores_sub %>%
                 left_join(centroids_sub, by = "combined"),
               aes(x = LD1.y, y = LD2.y, xend = LD1.x, yend = LD2.x, color = combined),
               alpha = 0.4, size = 0.6, inherit.aes = FALSE) +
  geom_point(size = 1.8) +
  stat_ellipse(type = "norm", level = 0.67, size = 0.6) +
  scale_color_manual(values = c(
    "Grazed + Unclipped" = "#8B3A62",
    "Ungrazed + Unclipped" = "#008B45"
  )) +
  labs(title = "DAPC of CpG Methylation in A.barbata (Unclipped only)",
       x = "LD1", y = "LD2",
       color = "Herbivory treatment") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    aspect.ratio = 1
  )






# CARGAR CHH PLANTAGO PRIMERO

# extract genomic positions
coords <- getData(CHH_meth)$start

# assign them as rownames to meth_matrix
rownames(CHH_meth_matrix) <- coords

grp <- as.factor(treatments)

# Run DAPC with initial max PCs
dapc_res <- dapc(t(CHH_meth_matrix), grp, n.pca = 14, n.da = length(levels(grp)) - 1)

# Optimize PCs
optim_a <- optim.a.score(dapc_res)
cat("Optimal number of PCs (alpha-score):", optim_a$best, "\n")

# Re-run with optimal number of PCs
dapc_res <- dapc(t(CHH_meth_matrix), grp, 
                 n.pca = optim_a$best,
                 n.da = length(levels(grp)) - 1)

# Extract loadings
loadings <- dapc_res$var.contr  
loadings_df <- melt(loadings,
                    varnames = c("CHH", "LD"),
                    value.name = "Loading")

# Scores dataframe
dapc_scores <- as.data.frame(dapc_res$ind.coord)
dapc_scores$id <- sample_ids
dapc_scores$combined <- grp

# Factor levels with labels
dapc_scores$combined <- factor(dapc_scores$combined,
                               levels = c("0", "2", "1", "3"),
                               labels = c("Grazed + Clipped",
                                          "Grazed + Unclipped",
                                          "Ungrazed + Clipped",
                                          "Ungrazed + Unclipped"))

# Centroids
centroids <- dapc_scores %>%
  group_by(combined) %>%
  summarise(LD1 = mean(LD1), LD2 = mean(LD2), .groups = "drop")

# --- FILTRAR SOLO LOS DOS GRUPOS DE INTERÉS ---
dapc_scores_sub <- dapc_scores %>%
  filter(combined %in% c("Grazed + Clipped", "Grazed + Unclipped"))

centroids_sub <- centroids %>%
  filter(combined %in% c("Grazed + Clipped", "Grazed + Unclipped"))

# Plot
ggplot(dapc_scores_sub, aes(x = LD1, y = LD2, color = combined)) +
  geom_segment(data = dapc_scores_sub %>%
                 left_join(centroids_sub, by = "combined"),
               aes(x = LD1.y, y = LD2.y, xend = LD1.x, yend = LD2.x, color = combined),
               alpha = 0.4, size = 0.6, inherit.aes = FALSE) +
  geom_point(size = 1.8) +
  stat_ellipse(type = "norm", level = 0.67, size = 0.6) +
  scale_color_manual(values = c(
    "Grazed + Clipped" = "#FF69B4",
    "Grazed + Unclipped" = "#8B3A62"
  )) +
  labs(title = "DAPC of CHH Methylation in P.lagopus",
       x = "LD1", y = "LD2",
       color = "Herbivory treatment") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    aspect.ratio = 1
  )







# Volcano plot function
plot_volcano <- function(diff_obj, min_diff = 10, q_cutoff = 0.05, title = "") {
  df <- as.data.frame(diff_obj)
  
  # Calculate variables
  df$logq <- -log10(df$qvalue)
  df$status <- ifelse(abs(df$meth.diff) >= min_diff & df$qvalue <= q_cutoff, "sig", "ns")
  
  # Plot
  p <- ggplot(df, aes(x = meth.diff, y = logq)) +
    geom_point(aes(color = status), alpha = 0.6, size = 1.5, show.legend = FALSE) +
    scale_color_manual(values = c("sig" = "red", "ns" = "black")) +
    geom_vline(xintercept = c(-min_diff, min_diff), linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -log10(q_cutoff), linetype = "dashed", color = "blue") +
    labs(
      x = "% Methylation Difference",
      y = expression(-log[10]("q-value")),
      title = title
    ) +
    theme_bw() +
    theme(
      axis.title = element_blank(),  
      axis.text = element_text(size = 9),
      plot.title = element_text(size = 11, hjust = 0.5)
    ) +
    # Fixed axes
    scale_x_continuous(limits = c(-100, 100)) +
    scale_y_continuous(limits = c(0, 55), expand = expansion(mult = c(0.02, 0.02)))
  
  return(p)
}

# Generate volcano plots
p1 <- plot_volcano(res_GUvsUU$all, title = "Ungrazed")
p2 <- plot_volcano(res_GUvsGC$all, title = "Clipping (Grazed)")
p3 <- plot_volcano(res_UUvsUC$all, title = "Clipping (Ungrazed)")


grid.arrange(
  p1, p2, p3, ncol = 3,
  left   = textGrob(expression(-log[10](q-value)), rot = 90, gp = gpar(fontsize = 12)),
  bottom = textGrob("% Methylation Difference", gp = gpar(fontsize = 12))
)





# I will create a function to then compare all against all

# Function for pairwise comparison
pairwise_diff <- function(meth_obj, g1, g2, min_diff = 10, q_cutoff = 0.05) {
  
  # Subset the object to only keep the two groups of interest
  sub_obj <- reorganize(
    meth_obj,
    sample.ids = getSampleID(meth_obj)[meth_obj@treatment %in% c(g1, g2)],
    treatment = meth_obj@treatment[meth_obj@treatment %in% c(g1, g2)]
  )
  
  # Run differential methylation
  diff <- calculateDiffMeth(
    sub_obj,
    overdispersion = "MN",
    test = "Chisq"
  )
  
  # Extract significant CHHs
  sig <- getMethylDiff(diff, difference = min_diff, qvalue = q_cutoff)
  
  return(list(all = diff, significant = sig))
}


# MAKE COMPARISONS

# FIRST 3 ARE USING GU AS CONTROL AND COMPARED TO THE 2 TREATMENTS AND THE COMBINATION OF BOTH TREATMENTS
# GRAZED VS. UNGRAZED (WITHIN UNCLIPPED)
res_GUvsUU <- pairwise_diff(CHH_meth, g1 = 2, g2 = 3)
# UNCLIPPED VS. CLIPPED (WITHIN GRAZED)
res_GUvsGC <- pairwise_diff(CHH_meth, g1 = 2, g2 = 0)
# GRAZED+UNCLIPPED VS. UNGRAZED+CLIPPED 
res_GUvsUC <- pairwise_diff(CHH_meth, g1 = 2, g2 = 1)
