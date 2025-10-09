
library(ggplot2)
library(ggrepel)
library(ggbreak)
library(dplyr)

# Semi-transparent point colors
point_colors <- c(
  "graz-unclip"   = alpha("#8B3A62", 0.5),
  "graz-clip"     = alpha("#FF69B4", 0.5),
  "ungraz-unclip" = alpha("#008B45", 0.5),
  "ungraz-clip"   = alpha("#00FF7F", 0.5)
)

# Outliers
CpGs_per_sample <- CpGs_per_sample %>%
  group_by(combined) %>%
  mutate(outlier = find_outliers(CpG)) %>%
  ungroup()

# Plot
ggplot(CpGs_per_sample, aes(x = combined, y = CpG)) +
  geom_boxplot(aes(fill = combined), outlier.shape = NA) + # boxplots opaque
  geom_point(aes(fill = combined), shape = 21, color = "black",
             size = 3, alpha = 0.5) + # semi-transparent points
  geom_text_repel(
    data = CpGs_per_sample %>% filter(outlier),
    aes(label = id),
    size = 3,
    max.overlaps = 10
  ) +
  scale_fill_manual(values = c(
    "graz-unclip"   = "#8B3A62",
    "graz-clip"     = "#FF69B4",
    "ungraz-unclip" = "#008B45",
    "ungraz-clip"   = "#00FF7F"
  )) +
  # Two breaks
  scale_y_break(c(6500, 750000), scales = 0.5) +
  scale_y_break(c(870000, 1300000), scales = 0.5) +
  # Optional: manual zig-zag at break positions
  geom_segment(aes(x = 0.6, xend = 1.4, y = 750000, yend = 750000),
               color = "black", size = 0.8, linetype = "22") +
  geom_segment(aes(x = 0.6, xend = 1.4, y = 1300000, yend = 1300000),
               color = "black", size = 0.8, linetype = "22") +
  theme_minimal() +
  labs(title = "CpG counts per treatment group (broken axis)",
       x = "Treatment", y = "CpG count")



# Your raw list of filenames (copy-paste or read from file)
files <- c(
  "HO301.1_bismark_bt2_SE_report.txt",
  "HO302.1_bismark_bt2_SE_report.txt",
  "HO305.1_bismark_bt2_SE_report.txt",
  "HO307.1_bismark_bt2_SE_report.txt",
  "HO308.1_bismark_bt2_SE_report.txt",
  "HO315.1_bismark_bt2_SE_report.txt",
  "HO316B.1_bismark_bt2_SE_report.txt",
  "HO317.1_bismark_bt2_SE_report.txt",
  "HO319.1_bismark_bt2_SE_report.txt",
  "HO322.1_bismark_bt2_SE_report.txt",
  "HO323.1_bismark_bt2_SE_report.txt",
  "HO326.1_bismark_bt2_SE_report.txt",
  "HO334.1_bismark_bt2_SE_report.txt",
  "HO337.1_bismark_bt2_SE_report.txt",
  "HO340A.1_bismark_bt2_SE_report.txt",
  "HO340B.1_bismark_bt2_SE_report.txt",
  "HO342.1_bismark_bt2_SE_report.txt",
  "HO346.1_bismark_bt2_SE_report.txt",
  "HO351.1_bismark_bt2_SE_report.txt",
  "HO353.1_bismark_bt2_SE_report.txt",
  "HO359.1_bismark_bt2_SE_report.txt"
)

# Remove the suffix to keep only the sample ID
ids <- sub("\\.1_bismark_bt2_SE_report\\.txt$", "", files)

# Subset your dataframe
subset_df <- CpGs_per_sample[CpGs_per_sample$id %in% ids, ]




# Example: your Bash output as a character vector
mapping_output <- c(
  "HO301.1_bismark_bt2_SE_report.txt:Mapping 0.4%",
  "HO302.1_bismark_bt2_SE_report.txt:Mapping 1.3%",
  "HO305.1_bismark_bt2_SE_report.txt:Mapping 0.1%",
  "HO307.1_bismark_bt2_SE_report.txt:Mapping 1.0%",
  "HO308.1_bismark_bt2_SE_report.txt:Mapping 0.4%",
  "HO315.1_bismark_bt2_SE_report.txt:Mapping 0.0%",
  "HO316B.1_bismark_bt2_SE_report.txt:Mapping 0.1%",
  "HO317.1_bismark_bt2_SE_report.txt:Mapping 0.7%",
  "HO319.1_bismark_bt2_SE_report.txt:Mapping 0.1%",
  "HO322.1_bismark_bt2_SE_report.txt:Mapping 0.0%",
  "HO323.1_bismark_bt2_SE_report.txt:Mapping 1.2%",
  "HO326.1_bismark_bt2_SE_report.txt:Mapping 0.0%",
  "HO334.1_bismark_bt2_SE_report.txt:Mapping 0.3%",
  "HO337.1_bismark_bt2_SE_report.txt:Mapping 1.0%",
  "HO340A.1_bismark_bt2_SE_report.txt:Mapping 0.3%",
  "HO340B.1_bismark_bt2_SE_report.txt:Mapping 0.1%",
  "HO342.1_bismark_bt2_SE_report.txt:Mapping 0.0%",
  "HO346.1_bismark_bt2_SE_report.txt:Mapping 0.2%",
  "HO351.1_bismark_bt2_SE_report.txt:Mapping 1.0%",
  "HO353.1_bismark_bt2_SE_report.txt:Mapping 0.0%",
  "HO359.1_bismark_bt2_SE_report.txt:Mapping 0.0%"
)

# Split into two parts: filename and percentage
parts <- do.call(rbind, strsplit(mapping_output, ":Mapping "))
filenames <- parts[,1]
percentages <- as.numeric(sub("%", "", parts[,2]))

# Keep only those with < 1%
low_map <- filenames[percentages < 1]

# Remove suffix to get sample IDs
ids <- sub("\\.1_bismark_bt2_SE_report\\.txt$", "", low_map)

# Subset your dataframe
subset_df <- CpGs_per_sample[CpGs_per_sample$id %in% ids, ]






# Step 1: CpG < 500 subset
low_cpg <- CpGs_per_sample[CpGs_per_sample$CpG < 500, ]

# Step 2: check if all IDs in low_cpg are also in subset_df
all_in <- all(low_cpg$id %in% subset_df$id)

# Step 3 (optional): list the ones that are NOT in subset_df
not_in <- setdiff(low_cpg$id, subset_df$id)

all_in
not_in





rest_cpg <- CpGs_per_sample[!(CpGs_per_sample$id %in% low_cpg$id), ]

