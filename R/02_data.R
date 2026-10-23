# =============================================================================
# 02_data.R — Carga de datos (datasets del libro + archivo propio) y variables
# =============================================================================

# Llamada canónica del libro para leer los .txt
read_book_data <- function(path) {
  read.table(path, header = TRUE, sep = "\t", na.strings = "NA", dec = ".",
             strip.white = TRUE, check.names = TRUE)
}

# Lee .txt/.tsv (tab), .csv (, o ;) o Excel
load_any <- function(path, name = basename(path)) {
  ext <- tolower(tools::file_ext(name))
  d <- switch(ext,
    txt = , tsv = read_book_data(path),
    csv = {
      d1 <- tryCatch(read.csv(path, strip.white = TRUE, na.strings = c("NA", "")), error = function(e) NULL)
      if (is.null(d1) || ncol(d1) < 2) d1 <- read.csv2(path, strip.white = TRUE, na.strings = c("NA", ""))
      d1
    },
    xlsx = , xls = as.data.frame(readxl::read_excel(path)),
    stop("Formato no soportado: .", ext, " (use .txt, .csv o .xlsx)")
  )
  # limpiar espacios en texto
  for (nm in names(d)) if (is.character(d[[nm]])) d[[nm]] <- trimws(d[[nm]])
  sanitize_names(d)
}

# Nombres de columna válidos para R: quita símbolos (%, /, -, espacios…),
# evita empezar por dígito y asegura unicidad. Deja en attr "renamed" el mapa viejo → nuevo.
sanitize_names <- function(d) {
  old <- names(d)
  new <- gsub("[^A-Za-z0-9_.]+", "_", trimws(old))      # símbolos y espacios → "_"
  new <- gsub("^_+|_+$", "", new)                        # sin guiones bajos en los extremos
  new <- ifelse(grepl("^[0-9.]", new) | !nzchar(new), paste0("V_", new), new)
  new <- make.unique(new, sep = "_")
  chg <- old != new
  names(d) <- new
  attr(d, "renamed") <- if (any(chg)) setNames(new[chg], old[chg]) else NULL
  d
}

bundled_datasets <- function(dir = "data") {
  f <- list.files(dir, pattern = "\\.txt$")
  setNames(file.path(dir, f), sub("\\.txt$", "", f))
}

# Tipo "práctico" de una columna
var_kind <- function(x) {
  x <- x[!is.na(x)]
  nu <- length(unique(x))
  if (nu == 0) return("empty")
  if (is.numeric(x)) {
    if (nu == 2) "binary"
    else if (nu <= 8 && all(x == round(x))) "numeric_small"   # entero con pocos valores: sirve como factor o numérica
    else "numeric"
  } else {
    if (nu == 2) "binary" else if (nu <= 30) "factor" else "text"
  }
}

vars_for_role <- function(d, r) {
  kinds <- sapply(d, var_kind)
  ok <- switch(r$type,
    numeric = kinds %in% c("numeric", "numeric_small", "binary"),
    factor  = kinds %in% c("factor", "binary", "numeric_small"),
    binary  = kinds == "binary",
    id      = rep(TRUE, length(kinds)),
    time    = kinds %in% c("numeric", "numeric_small"),
    event   = kinds %in% c("binary", "numeric_small"),
    any     = kinds != "empty"
  )
  names(d)[ok]
}

# Protege nombres no sintácticos (p. ej. "%VIAB", "Cure time") para usarlos en fórmulas
bq <- function(x) ifelse(make.names(x) == x, x, paste0("`", gsub("`", "", x), "`"))

# Accesores con coerción segura
num <- function(d, x) as.numeric(d[[x]])
fac <- function(d, x) droplevels(if (is.factor(d[[x]])) d[[x]] else factor(d[[x]]))

# Evento 0/1 a partir de binaria (numérica o texto)
as_event01 <- function(x, event_level = NULL) {
  if (is.numeric(x)) {
    u <- sort(unique(na.omit(x)))
    if (all(u %in% c(0, 1))) return(as.integer(x))
    return(as.integer(x == max(u)))
  }
  f <- factor(x)
  if (is.null(event_level)) event_level <- levels(f)[length(levels(f))]
  as.integer(f == event_level)
}

# Formato ancho (varias columnas de medición) → largo (Px, Times, Response)
wide_to_long <- function(d, measures, id = NULL) {
  Px <- if (is.null(id) || !nzchar(id)) seq_len(nrow(d)) else d[[id]]
  long <- data.frame(
    Px       = factor(rep(Px, times = length(measures))),
    Times    = factor(rep(measures, each = nrow(d)), levels = measures),
    Response = as.numeric(unlist(d[measures], use.names = FALSE))
  )
  long[!is.na(long$Response), ]
}

# Valida la selección de variables contra los roles del spec
validate_roles <- function(d, spec, v) {
  probs <- character(0)
  for (rn in names(spec$roles)) {
    r <- spec$roles[[rn]]
    sel <- v[[rn]]
    sel <- sel[nzchar(sel)]
    if (length(sel) == 0) {
      if (isTRUE(r$required)) probs <- c(probs, paste0("Falta: ", r$label))
      next
    }
    if (!is.null(r$min) && length(sel) < r$min)
      probs <- c(probs, paste0(r$label, ": se requieren al menos ", r$min, " variables"))
    if (!is.null(r$nlevels)) {
      nl <- length(unique(na.omit(d[[sel[1]]])))
      if (!is.null(r$nlevels) && nl != r$nlevels)
        probs <- c(probs, paste0(r$label, ": debe tener exactamente ", r$nlevels, " niveles (tiene ", nl, ")"))
    }
  }
  # variables repetidas en roles distintos
  allv <- unlist(v); allv <- allv[nzchar(allv)]
  dup <- unique(allv[duplicated(allv)])
  if (length(dup)) probs <- c(probs, paste0("Variable usada en más de un rol: ", paste(dup, collapse = ", ")))
  probs
}
