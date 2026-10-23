# =============================================================================
# 04_descriptives.R — Estadística descriptiva con table1 (como en el libro)
# Patrón del libro: table1(~ respuesta | factor, data = dat)
# =============================================================================

desc_table1 <- function(d, response, by = NULL) {
  response <- response[nzchar(response)]
  rhs <- paste(bq(response), collapse = " + ")
  f <- if (is.null(by) || !length(by) || !nzchar(by[1])) as.formula(paste("~", rhs))
       else as.formula(paste("~", rhs, "|", paste(bq(by), collapse = " + ")))
  dd <- d
  if (!is.null(by)) for (b in by) if (nzchar(b)) dd[[b]] <- fac(dd, b)
  tryCatch(table1::table1(f, data = dd), error = function(e) NULL)
}

# HTML para Shiny
table1_html <- function(t1) {
  if (is.null(t1)) return(htmltools::HTML("<em>No fue posible generar la tabla descriptiva.</em>"))
  htmltools::HTML(paste(unclass(t1), collapse = ""))
}

# Resumen numérico por grupo (n, media, DE, EE, mediana, IQR) — R base
desc_numeric <- function(y, g = NULL) {
  y <- as.numeric(y)
  if (is.null(g)) g <- factor(rep("Total", length(y)))
  ok <- !is.na(y) & !is.na(g)
  y <- y[ok]; g <- droplevels(factor(g[ok]))
  f <- function(v) c(n = length(v), Media = mean(v), DE = sd(v), EE = sd(v) / sqrt(length(v)),
                     Mediana = median(v), IQR = IQR(v), Min = min(v), Max = max(v))
  m <- t(sapply(split(y, g), f))
  data.frame(Grupo = rownames(m), round(m, 3), row.names = NULL, check.names = FALSE)
}
