# I wrote this script to:
# 1) combine CpG/CHG/CHH files into DMPs_all_<comparison>.txt (rbind + add CTXT)
# 2) for EVERY DMPs_*.txt file in the folder (the 9 context files + the 3 combined),
#    create:
#      - <basename>_ids.txt       (one ID per line, unique, no header)
#      - <basename>_up_ids.txt    (meth.diff > 0)
#      - <basename>_down_ids.txt  (meth.diff < 0)
#      - <basename>_ranked.txt    (contig_id <TAB> agg_meth.diff) -> suitable for GSEA

# Load packages
library(here)
library(tidyverse)  # for bind_rows, readr helpers

setwd(here("plantago_pseudo"))

# ----- Part A: combine contexts into DMPs_all_<comparison>.txt (as before) -----
comparisons <- c("GUvsGC", "GUvsUU", "UUvsUC")
contexts <- c("CpG", "CHG", "CHH")

# I will loop comparisons and read the three context files, add CTXT column and bind_rows
for (comp in comparisons) {
  dfs <- list()
  for (ctx in contexts) {
    file_name <- paste0("DMPs_", ctx, "_", comp, ".txt")
    if (!file.exists(file_name)) {
      warning("Expected file not found: ", file_name, " — skipping this context for ", comp)
      next
    }
    # I read with read.table because the files are tab-delimited and may not be huge.
    df <- read.table(file_name, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # Add CTXT column so we keep track of the context origin
    df$CTXT <- ctx
    dfs[[ctx]] <- df
    message("Read ", nrow(df), " rows from ", file_name)
  }
  
  # If no context files read, skip
  if (length(dfs) == 0) {
    warning("No context files found for comparison: ", comp, " — no combined file created.")
    next
  }
  
  # Combine
  df_all <- bind_rows(dfs)
  
  # Write combined file
  out_name <- paste0("DMPs_all_", comp, ".txt")
  write.table(df_all, file = out_name, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote combined file: ", out_name, " (", nrow(df_all), " rows)")
}

# ----- Part B: create ID lists (and up/down/ranked) for every DMPs_*.txt file -----
# I will detect all files that start with "DMPs_" and end with ".txt"
all_dmp_files <- list.files(pattern = "^DMPs_.*\\.txt$")

if (length(all_dmp_files) == 0) {
  stop("No DMPs_*.txt files found in the current directory. Aborting.")
}

message("Found ", length(all_dmp_files), " DMPs_*.txt files. Will create ID lists for each.")

# Helper function: safe column access (in case of slightly different headers)
get_col_safe <- function(df, candidates) {
  for (c in candidates) {
    if (c %in% names(df)) return(df[[c]])
  }
  return(NULL)
}

# Process each file
for (f in all_dmp_files) {
  message("\nProcessing file: ", f)
  df <- tryCatch(
    read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
    error = function(e) {
      warning("Failed to read ", f, " : ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(df)) next
  
  # Ensure contig_id exists (some files might call it 'contig_id' exactly)
  contig_col <- get_col_safe(df, c("contig_id", "contig", "seq_id", "ID"))
  if (is.null(contig_col)) {
    warning("No contig_id-like column in ", f, ". Columns found: ", paste(names(df), collapse=", "), " — skipping.")
    next
  }
  # Normalize contig IDs as character and trim whitespace
  contig_vec <- trimws(as.character(contig_col))
  # Unique IDs for the basic list
  ids <- unique(contig_vec)
  ids_out <- sub("\\.txt$", "_ids.txt", f)
  writeLines(ids, con = ids_out)
  message("  Wrote IDs: ", ids_out, " (", length(ids), " unique IDs)")
  
  # Check for meth.diff column (could be named exactly 'meth.diff' or similar)
  methdiff_col <- get_col_safe(df, c("meth.diff", "meth_diff", "methdiff"))
  if (is.null(methdiff_col)) {
    # If no meth.diff, still create empty up/down files (informative)
    message("  No meth.diff column found in ", f, " — creating empty up/down lists.")
    writeLines(character(0), con = sub("\\.txt$", "_up_ids.txt", f))
    writeLines(character(0), con = sub("\\.txt$", "_down_ids.txt", f))
    # create an empty ranked file
    write.table(data.frame(), file = sub("\\.txt$", "_ranked.txt", f), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    next
  }
  
  # Build a data.frame pairing contig_id and meth.diff (coerce to numeric)
  md <- as.numeric(methdiff_col)
  pair_df <- data.frame(contig_id = contig_vec, meth.diff = md, stringsAsFactors = FALSE)
  
  # Remove NA contig_id rows
  pair_df <- pair_df[!is.na(pair_df$contig_id) & pair_df$contig_id != "", , drop = FALSE]
  
  # UP (meth.diff > 0) and DOWN (meth.diff < 0)
  up_ids <- unique(pair_df$contig_id[!is.na(pair_df$meth.diff) & pair_df$meth.diff > 0])
  down_ids <- unique(pair_df$contig_id[!is.na(pair_df$meth.diff) & pair_df$meth.diff < 0])
  
  writeLines(up_ids, con = sub("\\.txt$", "_up_ids.txt", f))
  writeLines(down_ids, con = sub("\\.txt$", "_down_ids.txt", f))
  message("  Wrote up IDs:   ", sub("\\.txt$", "_up_ids.txt", f), " (", length(up_ids), ")")
  message("  Wrote down IDs: ", sub("\\.txt$", "_down_ids.txt", f), " (", length(down_ids), ")")
  
  # Create ranked file: aggregate by contig_id (mean meth.diff) so we have 1 value per contig
  ranked_agg <- pair_df %>%
    filter(!is.na(meth.diff)) %>%
    group_by(contig_id) %>%
    summarise(rank_value = mean(meth.diff, na.rm = TRUE)) %>%
    ungroup()
  
  ranked_out <- sub("\\.txt$", "_ranked.txt", f)
  # For GSEA we typically want "ID TAB value" with no header
  if (nrow(ranked_agg) > 0) {
    write.table(ranked_agg, file = ranked_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    message("  Wrote ranked file: ", ranked_out, " (", nrow(ranked_agg), " IDs)")
  } else {
    # write an empty file to keep consistent outputs
    write.table(data.frame(), file = ranked_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    message("  No valid ranked entries; wrote empty file: ", ranked_out)
  }
}

message("\nAll done. You should now have for each DMPs_*.txt file:")
message("  - *_ids.txt      (unique contig IDs, one per line)")
message("  - *_up_ids.txt   (contigs with meth.diff > 0)")
message("  - *_down_ids.txt (contigs with meth.diff < 0)")
message("  - *_ranked.txt   (contig_id <TAB> mean_meth.diff)")

# End of script
