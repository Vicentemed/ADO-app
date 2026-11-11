# =============================================================================
# analyses/descriptiva.R - Tabla 1 de estadística descriptiva
# Sigue la lógica de table1::table1() tal como se usa en el libro:
#   table1(~ respuesta1 + respuesta2 + ... | explicativa1 * explicativa2, data = dat)
# Las variables de respuesta (izquierda de "|") son las filas de la tabla y las
# variables explicativas o de grupo (derecha de "|") definen las columnas.
#   - Sin explicativas:  table1(~ y1 + y2, data = dat)   → una sola columna.
#   - Una explicativa:   table1(~ VH | Group)            → una columna por nivel + Overall.
#   - Dos explicativas:  table1(~ VH | Group*Cure.time)  → columnas anidadas (ANOVA de dos vías).
#   - Tres o más: table1 solo admite un nivel de anidamiento, así que la primera queda
#     como nivel superior y las demás se combinan en una sola variable (p. ej. "Dual · F").
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
  rows <- unique(c(v$response, v$factors)); rows <- rows[nzchar(rows)]
  cols <- unique(v$explanatory); cols <- cols[nzchar(cols)]
  if (!length(rows)) stop("Seleccione al menos una variable de respuesta (filas de la tabla).")
  if (length(intersect(rows, cols))) stop("Una misma variable no puede estar a ambos lados de «|»: ", paste(intersect(rows, cols), collapse = ", "))
  dd <- d[, c(rows, cols), drop = FALSE]
  for (x in rows) {
    k <- var_kind(dd[[x]])
    if (x %in% v$factors || !k %in% c("numeric", "numeric_small")) dd[[x]] <- factor(dd[[x]]) else dd[[x]] <- as.numeric(dd[[x]])
  }
  for (g in cols) dd[[g]] <- fac(dd, g)
  if (length(cols)) dd <- dd[complete.cases(dd[, cols, drop = FALSE]), , drop = FALSE]
  strata <- cols; comb <- NULL
  if (length(cols) > 2) {   # table1: «Only 1 level of nesting is supported»
    comb <- paste(cols[-1], collapse = " · ")
    dd[[comb]] <- droplevels(interaction(dd[cols[-1]], sep = " · ", lex.order = TRUE))
    strata <- c(cols[1], comb)
  }
  f <- as.formula(paste("~", paste(bq(rows), collapse = " + "), if (length(strata)) paste("|", paste(bq(strata), collapse = " * ")) else ""))
  list(data = dd, formula = f, rows = rows, cols = cols, strata = strata, comb = comb,
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
  description = "Tabla descriptiva con la lógica de table1(~ respuestas | explicativas): las variables de respuesta son las filas y las variables explicativas o de grupo forman las columnas. Admite una o más variables a cada lado de «|».",
  chapter = "Resumen de datos (función table1)",
  roles = list(
    response    = role("Variables de respuesta (izquierda de «|»; filas de la tabla)", "any", multiple = TRUE, min = 1,
                       help = "Una o más. Las numéricas se resumen con media (DE) y mediana; las categóricas, con n (%)."),
    explanatory = role("Variables explicativas o de grupo (derecha de «|»; columnas; opcional)", "factor", multiple = TRUE, required = FALSE,
                       help = "Una: una columna por nivel, como en table1(~ VH | Group). Dos: columnas anidadas, como en table1(~ VH | Group*Cure.time) (ANOVA de dos vías, o grupo y sexo). Tres o más: la primera es el nivel superior y las demás se combinan. Vacío: una sola columna con todos los datos."),
    factors     = role("Variables de respuesta categóricas codificadas con números (opcional)", "factor", multiple = TRUE, required = FALSE,
                       help = "Se añaden a las filas y se muestran como categorías (n, %) en lugar de media y mediana.")),
  example = list(file = "VickersF.txt", vars = list(response = "VH", explanatory = c("Group", "Cure.time")),
                 note = "Dureza Vickers por cemento y protocolo de polimerización: table1(~ VH | Group*Cure.time), como en el ANOVA de dos vías (cap. 24)."),
  alternative = NULL,
  descriptives = function(d, v) tryCatch(.t1_make(d, v)$t1, error = function(e) NULL),
  assumptions = function(d, v) {
    b <- tryCatch(.t1_build(d, v), error = function(e) e)
    if (inherits(b, "error")) return(list(checks = list(check_obj("tabla1", "Selección de variables", FALSE, conditionMessage(b))),
                                          decision = list(parametric = TRUE, all_met = FALSE, violated = character(0), reason = conditionMessage(b))))
    nmiss <- sapply(b$data[b$rows], function(x) sum(is.na(x)))
    txt <- paste0("Llamada equivalente: ", b$call, ". Filas analizadas: ", nrow(b$data),
                  if (length(b$cols) && nrow(b$data) < nrow(d)) paste0(" (se excluyen ", nrow(d) - nrow(b$data), " sin dato en las variables explicativas)") else "", ". ",
                  if (!is.null(b$comb)) paste0("table1 admite un solo nivel de anidamiento: ", b$cols[1], " queda como nivel superior y ", paste(b$cols[-1], collapse = ", "), " se combinan en una variable. ") else "",
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
                      levels = if (length(b$cols)) lapply(b$data[b$cols], levels) else list(), comb = b$comb),
         tables = list(),
         notes = paste0("Llamada equivalente: ", b$call, ". ", if (!is.null(b$comb)) paste0("Con más de dos variables explicativas, ", paste(b$cols[-1], collapse = " y "), " se combinan en una sola. ") else "", "La tabla se incluye en el reporte Word en la sección de estadística descriptiva."))
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
           if (length(s$cols)) paste0(" según ", CTX(paste(s$cols, collapse = ", ")), " (", paste(unlist(lapply(s$levels, paste, collapse = ", ")), collapse = "; "), ")") else "",
           " en ", s$n, " ", CTX("sujetos"), ". ",
           if (length(s$num)) paste0("Las variables numéricas se resumen como ", st, ". ") else "",
           if (length(s$cat)) "Las variables categóricas se resumen como frecuencia (porcentaje). " else "",
           CTX("Describir las características relevantes de la muestra."))
  }
))
