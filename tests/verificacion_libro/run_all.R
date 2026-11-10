# =============================================================================
# tests/verificacion_libro/run_all.R
# Corre TODOS los análisis de la app sobre los ejemplos del libro y compara cada
# cifra con la impresa en el libro. Uso (desde la carpeta de la app):
#   Rscript tests/verificacion_libro/run_all.R
# Salida: tests/verificacion_libro/resultados.csv  (una fila por cifra comparada)
#   match = TRUE  → la app reproduce la cifra del libro
#   match = FALSE → difiere, o la app no la muestra (ver columna note)
# La columna `line` es la línea del texto del libro (versión Final101126) donde
# aparece la cifra.
# =============================================================================
here <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), winslash = "/")
Sys.setenv(VERIFY_DIR = here)
rs <- file.path(R.home("bin"), "Rscript")
scripts <- c(t = "verify_t.R", anova = "verify_anova.R", rm = "verify_rm.R", cat = "verify_cat.R", reg = "verify_reg.R")
all <- list()
for (k in names(scripts)) {
  cat("Corriendo", scripts[[k]], "...\n")
  st <- system2(rs, shQuote(file.path(here, scripts[[k]])), stdout = FALSE, stderr = FALSE)
  f <- file.path(here, paste0("results_", k, ".csv"))
  if (st != 0 || !file.exists(f)) { cat("  ERROR al correr", scripts[[k]], "\n"); next }
  r <- read.csv(f, fileEncoding = "UTF-8", stringsAsFactors = FALSE); r$familia <- k; all[[k]] <- r
}
R <- do.call(rbind, all)
R$estado <- ifelse(R$match, "coincide", ifelse(grepl("no lo muestra", R$note) | grepl("^NA", R$app), "la app no lo muestra", "no coincide"))
write.csv(R, file.path(here, "resultados.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("\n== RESUMEN POR ANÁLISIS ==\n"); print(table(R$analysis, R$estado))
cat("\n== TOTAL ==\n"); print(table(R$estado))
cat("\n== NO COINCIDEN ==\n"); print(R[R$estado == "no coincide", c("analysis", "item", "book", "app", "line", "note")], row.names = FALSE)
