# Extraer matriz interna del methylBase
meth <- getData(CpG_meth_filtered_mef)

# Índices de columnas de coverage
cov_idx <- CpG_meth_filtered_mef@coverage.index

# Vector con tratamientos (uno por muestra)
treatment <- CpG_meth_filtered_mef@treatment

# Crear lista de columnas de cobertura por tratamiento
group_idx <- split(cov_idx, treatment)

# Generar tabla por tratamiento
for (grp in names(group_idx)) {
  idx <- group_idx[[grp]]
  n_present <- rowSums(!is.na(meth[, idx]))
  total_sites <- length(n_present)
  
  res <- data.frame(
    Muestras = paste0("≥", 1:10),
    Sitios = sapply(1:10, function(n) sum(n_present >= n))
  )
  res$Porcentaje <- round((res$Sitios / total_sites) * 100, 2)
  
  cat("\n=== Tratamiento", grp, "===\n")
  print(res, row.names = FALSE)
}

