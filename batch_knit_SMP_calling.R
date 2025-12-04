# batch_render_parallel.R

library(rmarkdown)
library(future)

# Use multisession for parallel on Windows
plan(multisession)

# Project root
library(here)

# Directories to search
dirs <- c("avena_ncbi", "plantago_pseudo")

# Find all *_results.Rmd in those directories
files <- unlist(lapply(dirs, function(d) {
  list.files(
    path = here(d),
    pattern = "_results\\.Rmd$",
    recursive = TRUE,
    full.names = TRUE
  )
}))

# Print files to check
cat("Found the following Rmd files:\n")
print(files)

# Knit all files in parallel
futures <- lapply(files, function(f) {
  future({
    cat("Rendering:", f, "\n")
    rmarkdown::render(f, quiet = TRUE)
  })
})

# Wait for all to finish
lapply(futures, value)

cat("All reports rendered!\n")
