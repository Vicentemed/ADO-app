# =============================================================================
# analyses/descriptiva.R - Tabla 1 de estadística descriptiva
# Sigue la lógica de table1::table1() tal como se usa en el libro:
#   table1(~ explicativa1 + explicativa2 + ... | respuesta, data = dat)
# Las variables de la izquierda son las filas de la tabla y la variable a la
# derecha de "|" define las columnas (una por categoría, más "Overall").
# Si no se indica variable de columnas: table1(~ x1 + x2, data = dat).
# Con dos variables de columnas: table1(~ x | g1*g2, data = dat) (anidadas).
# =============================================================================

.t1_stats <- c("Mean (SD) y Median [Min, Max] (table1 por defecto)" = "default",
               "Mean (SD)" = "mean_sd",
               "Median [Q1, Q3]" = "median_iqr",
               "Mean (SD) y Median [Q1, Q3]" = "mean_median_iqr")

# Estadísticos de las variables continuas (render.continuous de table1)
.t1_render_cont <- function(kind) {
  if (is.null(kind) || identical(kind, "default")) return(NULL)
  function(x) {
    s <- table1::stats.apply.rounding(table1::stats.default(x), digits = 3)
    out <- c("")
    if (kind %in% c("mean_sd", "mean_median_iqr")) out <- c(out, "Mean (SD)" = sprintf("%s (%s)", s$MEAN, s$SD))
    if (kind %in% c("median_iqr", "mean_median_iqr")) out <- c(out, "Median [Q1, Q3]" = sprintf("%s [%s, %s]", s$MEDIAN, s$Q1, s$Q3))
    out
  }
}

# Prepara los datos y la fórmula de table1
.t1_build <- function(d, v) {
  rows <- unique(c(v$explanatory, v$factors)); rows <- rows[nzchar(rows)]
  cols <- v$response; cols <- cols[nzchar(cols)]
  if (!length(rows)) stop("Seleccione al menos una variable explicativa (filas de la tabla).")
  if (length(cols) > 2) stop("Seleccione como máximo dos variables de columnas.")
  if (length(intersect(rows, cols))) stop("Una misma variable no puede estar en las filas y en las columnas: ", paste(intersect(rows, cols), collapse = ", "))
  dd <- d[, c(rows, cols), drop = FALSE]
  for (x in rows) {
    k <- var_kind(dd[[x]])
    if (x %in% v$factors || !k %in% c("numeric", "numeric_small")) dd[[x]] <- factor(dd[[x]]) else dd[[x]] <- as.numeric(dd[[x]])
  }
  for (g in cols) dd[[g]] <- fac(dd, g)
  if (length(cols)) dd <- dd[complete.cases(dd[, cols, drop = FALSE]), , drop = FALSE]
  f <- as.formula(paste("~", paste(bq(rows), collapse = " + "), if (length(cols)) paste("|", paste(bq(cols), collapse = " * ")) else ""))
  list(data = dd, formula = f, rows = rows, cols = cols,
       call = paste0("table1(~ ", paste(rows, collapse = " + "), if (length(cols)) paste0(" | ", paste(cols, collapse = " * ")) else "", ", data = dat)"))
}
.t1_make <- function(d, v, opts = list()) {
  b <- .t1_build(d, v)
  overall <- if (isFALSE(opts$overall) && length(b$cols)) FALSE else "Overall"   # sin variable de columnas, «Overall» es la única columna
  rc <- .t1_render_cont(opts$stats)
  args <- list(b$formula, data = b$data, overall = overall)
  if (!is.null(rc)) args$render.continuous <- rc
  list(t1 = do.call(table1::table1, args), b = b)
}

ado_register(list(
  id = "tabla1", family = "descriptiva", order = 1,
  label = "Tabla 1 (estadística descriptiva)",
  description = "Tabla descriptiva con la lógica de table1(~ explicativas | respuesta): las variables explicativas son las filas y las categorías de la variable de respuesta (o de grupo) son las columnas.",
  chapter = "Resumen de datos (función table1)",
  roles = list(
    explanatory = role("Variables explicativas (filas de la tabla)", "any", multiple = TRUE, min = 1,
                       help = "Las variables numéricas se resumen con media (DE) y mediana; las categóricas, con n (%)."),
    response    = role("Variable de respuesta o de grupo (columnas; opcional, hasta 2)", "factor", multiple = TRUE, required = FALSE,
                       help = "Cada categoría forma una columna, como en table1(~ edad + sexo | grupo). Vacío: una sola columna con todos los datos. Para describir una respuesta numérica por grupo, ponga la respuesta en las filas y el grupo aquí."),
    factors     = role("Variables explicativas categóricas codificadas con números (opcional)", "factor", multiple = TRUE, required = FALSE,
                       help = "Se añaden a las filas y se muestran como categorías (n, %) en lugar de media y mediana.")),
  example = list(file = "datAsthma.txt", vars = list(explanatory = c("Age", "Sex", "OrParaF"), response = "Asthma"),
                 note = "Edad, sexo y parafunción oral según el diagnóstico de asma (cap. 40)."),
  alternative = NULL,
  descriptives = function(d, v) tryCatch(.t1_make(d, v)$t1, error = function(e) NULL),
  assumptions = function(d, v) {
    b <- tryCatch(.t1_build(d, v), error = function(e) e)
    if (inherits(b, "error")) return(list(checks = list(check_obj("tabla1", "Selección de variables", FALSE, conditionMessage(b))),
                                          decision = list(parametric = TRUE, all_met = FALSE, violated = character(0), reason = conditionMessage(b))))
    nmiss <- sapply(b$data[b$rows], function(x) sum(is.na(x)))
    txt <- paste0("Llamada equivalente: ", b$call, ". Filas analizadas: ", nrow(b$data),
                  if (length(b$cols) && nrow(b$data) < nrow(d)) paste0(" (se excluyen ", nrow(d) - nrow(b$data), " sin dato en la variable de columnas)") else "", ". ",
                  if (any(nmiss > 0)) paste0("Valores faltantes: ", paste0(names(nmiss)[nmiss > 0], " = ", nmiss[nmiss > 0], collapse = "; "), " (table1 los muestra en la fila «Missing»).") else "Sin valores faltantes en las variables de las filas.")
    list(checks = list(check_obj("tabla1", "Estructura de la tabla", NA, txt)),
         decision = list(parametric = TRUE, all_met = TRUE, violated = character(0), reason = "La Tabla 1 es descriptiva: no requiere supuestos ni pruebas de hipótesis."))
  },
  run = function(d, v, opts) {
    m <- .t1_make(d, v, opts); b <- m$b
    num <- b$rows[sapply(b$data[b$rows], is.numeric)]; cat_ <- setdiff(b$rows, num)
    list(main = m$t1, t1 = m$t1,
         stats = list(n = nrow(b$data), rows = b$rows, cols = b$cols, num = num, cat = cat_, call = b$call,
                      stats = if (is.null(opts$stats)) "default" else opts$stats, overall = !isFALSE(opts$overall) && length(b$cols) > 0,
                      levels = if (length(b$cols)) lapply(b$data[b$cols], levels) else list()),
         tables = list(),
         notes = paste0("Llamada equivalente: ", b$call, ". La tabla se incluye en el reporte Word en la sección de estadística descriptiva."))
  },
  posthoc = NULL,
  plots = list(),
  opts_ui = list(
    stats   = list(label = "Estadísticos de las variables numéricas", type = "select", choices = .t1_stats, value = "default"),
    overall = list(label = "Incluir la columna «Overall» (total)", type = "checkbox", value = TRUE)),
  report = function(res, d, v, checks) {
    s <- res$stats
    st <- switch(s$stats, mean_sd = "media (desviación estándar)", median_iqr = "mediana [cuartil 1, cuartil 3]",
                 mean_median_iqr = "media (desviación estándar) y mediana [cuartil 1, cuartil 3]", "media (desviación estándar) y mediana [mínimo, máximo]")
    paste0("La Tabla 1 presenta la estadística descriptiva de ", CTX(paste(s$rows, collapse = ", ")),
           if (length(s$cols)) paste0(" según ", CTX(paste(s$cols, collapse = " y ")), " (", paste(unlist(lapply(s$levels, paste, collapse = ", ")), collapse = "; "), ")") else "",
           " en ", s$n, " ", CTX("sujetos"), ". ",
           if (length(s$num)) paste0("Las variables numéricas se resumen como ", st, ". ") else "",
           if (length(s$cat)) "Las variables categóricas se resumen como frecuencia (porcentaje). " else "",
           CTX("Describir las características relevantes de la muestra."))
  }
))
