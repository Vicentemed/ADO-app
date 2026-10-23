# =============================================================================
# 01_registry.R — Registro genérico de análisis
# Cada análisis se DECLARA con cdo_register(); la interfaz es genérica.
#
# Esquema de una entrada (spec):
#   id, family, label, order, description, chapter
#   roles      : lista nombrada de role() → qué variables pide el análisis
#   example    : list(file, prep = function(d) d, vars = list(rol = "columna"))
#   descriptives(d, v)          → objeto table1 (o data.frame)
#   assumptions(d, v)           → list(checks = list(...), decision = list(...))
#   alternative                 → id de la alternativa no paramétrica (o NULL)
#   run(d, v, opts)             → res: list(main=, stats=list(), tables=list(), notes=)
#   posthoc(res, d, v, opts)    → list(tables = list(), notes =) o NULL
#   plots      : lista nombrada de list(label, fn = function(d, v, res, cos))
#   report(res, d, v, checks)   → texto markdown "Cómo reportar" con valores
#   opts_ui    : lista de opciones (selectInput) que la UI ofrece → opts
# =============================================================================

# El registro vive en options(): es único para todo el proceso, sin importar
# en qué entorno haya sido cargado cada archivo (Shiny carga app.R y los
# archivos anidados en un entorno local, no en el global).
.cdo_env <- function() {
  e <- getOption("cdo.registry")
  if (is.null(e)) { e <- new.env(parent = emptyenv()); options(cdo.registry = e) }
  e
}
.CDO_REGISTRY <- .cdo_env()

FAMILIES <- c(
  medias        = "Comparación de medias",
  correlacion   = "Correlación",
  categoricos   = "Datos categóricos",
  concordancia  = "Concordancia",
  regresion     = "Regresión",
  supervivencia = "Supervivencia"
)

role <- function(label, type = c("numeric", "factor", "binary", "id", "time", "event", "any"),
                 multiple = FALSE, required = TRUE, nlevels = NULL, min = NULL, help = NULL) {
  type <- match.arg(type)
  list(label = label, type = type, multiple = multiple, required = required,
       nlevels = nlevels, min = min, help = help)
}

cdo_register <- function(spec) {
  stopifnot(is.list(spec), !is.null(spec$id), !is.null(spec$family))
  if (is.null(spec$order)) spec$order <- 99
  if (is.null(spec$plots)) spec$plots <- list()
  if (is.null(spec$opts_ui)) spec$opts_ui <- list()
  assign(spec$id, spec, envir = .cdo_env())
  invisible(spec)
}
cdo_get  <- function(id) get(id, envir = .cdo_env(), inherits = FALSE)
cdo_has  <- function(id) exists(id, envir = .cdo_env(), inherits = FALSE)
cdo_list <- function(family = NULL) {
  ids <- ls(.cdo_env())
  specs <- lapply(ids, cdo_get)
  if (!is.null(family)) specs <- Filter(function(s) identical(s$family, family), specs)
  ord <- order(sapply(specs, function(s) s$family), sapply(specs, function(s) s$order))
  specs[ord]
}
cdo_choices <- function(family) {
  sp <- cdo_list(family)
  setNames(sapply(sp, `[[`, "id"), sapply(sp, `[[`, "label"))
}

# Rellena {clave} en una plantilla
fill_tpl <- function(tpl, vals) {
  for (k in names(vals)) {
    tpl <- gsub(paste0("{", k, "}"), as.character(vals[[k]]), tpl, fixed = TRUE)
  }
  tpl
}

# Objeto de verificación de supuesto (estructura común)
check_obj <- function(kind, label, met, text, ...) {
  c(list(kind = kind, label = label, met = met, text = text), list(...))
}

# Decisión paramétrico / no paramétrico según el libro:
# si algún supuesto se viola y existe alternativa → no paramétrico
decide_path <- function(checks, has_alternative) {
  mets <- sapply(checks, function(cch) cch$met)
  mets <- mets[!is.na(mets)]
  all_met <- length(mets) == 0 || all(mets)
  violated <- sapply(checks, function(cch) isFALSE(cch$met))
  list(
    parametric = all_met || !has_alternative,
    all_met    = all_met,
    violated   = sapply(checks, `[[`, "label")[violated],
    reason     = if (all_met) "Se cumplen los supuestos: se aplica la prueba paramétrica."
                 else if (has_alternative) paste0("Se viola: ", paste(sapply(checks, `[[`, "label")[violated], collapse = ", "),
                                                  ". Según el libro, se opta por la alternativa no paramétrica.")
                 else paste0("Se viola: ", paste(sapply(checks, `[[`, "label")[violated], collapse = ", "),
                             ". Este análisis no tiene alternativa no paramétrica en el libro; interpretar con cautela.")
  )
}
