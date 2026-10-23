# =============================================================================
# 00_config.R — Configuración global de la app ADO
# (Análisis de Datos en Odontología. Un enfoque pragmático con R)
# Estilo: R base, como en el libro.
# =============================================================================

ALPHA <- 0.05

# Paleta y tema del libro (ver cap. Resumen de datos y familia de medias)
ADO_COL  <- c("#2E86C1", "#E74C3C", "#27AE60", "#F39C12", "#8E44AD", "#555555")
ADO_FILL <- c("#AED6F1", "#FADBD8", "#A9DFBF", "#FDEBD0", "#E8DAEF", "#D5D8DC")

ado_theme <- function(base_size = 12) ggplot2::theme_bw(base_size = base_size)

# pandoc para rmarkdown cuando la app corre fuera de RStudio (reporte Word)
if (!rmarkdown::pandoc_available()) {
  cand <- c("C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
            "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/x86_64",
            "C:/Program Files/Quarto/bin/tools")
  for (p in cand) if (dir.exists(p)) { Sys.setenv(RSTUDIO_PANDOC = p); break }
}

# ---- Formateadores para el texto "Cómo reportar" --------------------------
fmt_num <- function(x, d = 2) {
  if (is.null(x) || length(x) == 0) return("…")
  vapply(as.numeric(x), function(xi) if (is.na(xi)) "…" else formatC(xi, digits = d, format = "f"), character(1))
}
# Devuelve "= 0.034" o "< 0.001" para escribir "p {p}"
fmt_p <- function(p, d = 3) {
  if (is.null(p) || length(p) == 0) return("= …")
  vapply(as.numeric(p), function(pi) if (is.na(pi)) "= …" else if (pi < 0.001) "< 0.001" else paste0("= ", formatC(pi, digits = d, format = "f")), character(1))
}
# Interpretación del tamaño de efecto r (Cohen 1988 / libro)
efecto_r <- function(r) {
  if (is.na(r)) return("…")
  r <- abs(r)
  if (r < 0.1) "trivial" else if (r < 0.3) "pequeño" else if (r < 0.5) "mediano" else "grande"
}
# eta² / omega² (Cohen): 0.01 pequeño, 0.06 mediano, 0.14 grande
efecto_eta <- function(e) {
  if (is.na(e)) return("…")
  if (e < 0.01) "trivial" else if (e < 0.06) "pequeño" else if (e < 0.14) "mediano" else "grande"
}
# W de Kendall (libro cap. Friedman)
efecto_W <- function(w) {
  if (is.na(w)) return("…")
  if (w < 0.1) "trivial" else if (w < 0.3) "pequeño" else if (w < 0.5) "mediano" else "grande"
}

# Placeholder que el usuario debe completar con el contexto de sus datos
CTX <- function(txt) paste0("[", txt, "]")
