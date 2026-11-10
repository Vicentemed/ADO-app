# =============================================================================
# verify_lib.R — Harness común: corre la app (sus funciones registradas) sobre los
# ejemplos del libro y compara cada cifra con la impresa en el libro.
# Uso en cada verify_<familia>.R:
#   source(file.path(Sys.getenv("VERIFY_DIR"), "verify_lib.R"))
#   x <- run_app("kruskal", data = book_data("ReMat.txt"), vars = list(...), opts = list(...))
#   chk("kruskal", "H", "51.86", x$res$stats$H, line = 6890)
# =============================================================================
APP_DIR  <- normalizePath(file.path(Sys.getenv("VERIFY_DIR"), "..", ".."), winslash = "/")   # carpeta de la app
BOOK_TXT <- file.path(Sys.getenv("VERIFY_DIR"), "libro_Final101126.txt")   # texto del libro extraído del .docx
BOOK_DATA_DIRS <- c(file.path(APP_DIR, "data"), Sys.getenv("ADO_BOOK_DATA", unset = file.path(APP_DIR, "data")))   # ADO_BOOK_DATA: carpeta opcional con los sets de datos del libro

if (!exists(".verify_loaded")) {
  owd <- setwd(APP_DIR)
  suppressPackageStartupMessages({ library(ggplot2) })
  for (f in c("00_config.R", "01_registry.R", "02_data.R", "03_assumptions.R", "04_descriptives.R",
              "05_plots.R", "06_report.R")) source(file.path("R", f), encoding = "UTF-8")
  for (f in list.files(file.path("R", "analyses"), full.names = TRUE)) source(f, encoding = "UTF-8")
  setwd(owd)
  .verify_loaded <- TRUE
  RESULTS <- data.frame(analysis = character(), item = character(), book = character(), app = character(),
                        match = logical(), line = integer(), note = character(), stringsAsFactors = FALSE)
}

# Lee un set de datos del libro (busca en data/ de la app y en la carpeta del libro)
book_data <- function(file) {
  for (dd in BOOK_DATA_DIRS) { p <- file.path(dd, file); if (file.exists(p)) return(read_book_data(p)) }
  stop("No encontré ", file)
}

# Corre una entrada registrada de la app exactamente como lo hace la app:
# supuestos → run → post-hoc → reporte. Devuelve todo para inspección.
run_app <- function(id, data, vars, opts = NULL) {
  s <- ado_get(id)
  if (is.null(opts)) opts <- lapply(s$opts_ui, `[[`, "value")
  else { def <- lapply(s$opts_ui, `[[`, "value"); for (n in names(opts)) def[[n]] <- opts[[n]]; opts <- def }
  pr <- validate_roles(data, s, vars); if (length(pr)) stop(id, ": roles inválidos: ", paste(pr, collapse = "; "))
  a <- s$assumptions(data, vars)
  r <- s$run(data, vars, opts)
  ph <- if (!is.null(s$posthoc)) tryCatch(s$posthoc(r, data, vars, opts), error = function(e) list(error = conditionMessage(e))) else NULL
  rep <- tryCatch(s$report(r, data, vars, a$checks), error = function(e) paste("ERROR reporte:", conditionMessage(e)))
  list(spec = s, assumptions = a, res = r, posthoc = ph, report = rep, opts = opts)
}

# Compara un valor del libro (tal como está impreso, en texto) con el de la app.
# Coincide si la app, redondeada a los decimales impresos en el libro, da lo mismo.
# Para p impresos como "< 0.001" / "< 2e-16" pase book = "<0.001" y la app debe ser menor.
# Texto: pase app como character y book como character → igualdad exacta (sin espacios extremos).
chk <- function(analysis, item, book, app, line = NA_integer_, note = "", tol = NULL) {
  book <- trimws(as.character(book))
  if (is.null(app) || length(app) == 0 || (length(app) == 1 && is.na(app))) {
    ok <- FALSE; app_s <- "NA / no encontrado"
  } else if (is.character(app) && !grepl("^<", book) && is.na(suppressWarnings(as.numeric(book)))) {
    ok <- identical(trimws(app[1]), book); app_s <- app[1]
  } else if (grepl("^<", book)) {
    lim <- as.numeric(sub("^<\\s*", "", book)); a <- as.numeric(app[1]); ok <- a < lim; app_s <- format(signif(a, 6))
  } else {
    b <- as.numeric(book); a <- as.numeric(app[1])
    dec <- if (grepl("e", book, ignore.case = TRUE)) NA else { s <- sub("^-", "", book); if (grepl("\\.", s)) nchar(sub("^[^.]*\\.", "", s)) else 0 }
    if (!is.null(tol)) ok <- abs(a - b) <= tol
    else if (is.na(dec)) ok <- abs(a - b) <= abs(b) * 1e-3
    else ok <- isTRUE(abs(round(a, dec) - b) < 10^(-dec) / 2 + 1e-12) || isTRUE(abs(a - b) <= 10^(-dec) / 2 + 1e-12)
    app_s <- format(signif(a, 7))
  }
  RESULTS[nrow(RESULTS) + 1, ] <<- list(analysis, item, book, app_s, isTRUE(ok), as.integer(line), note)
  invisible(ok)
}

# Utilidad: busca un valor en una tabla de la app por nombre de tabla (regex) y columna/fila
tab_val <- function(x, table_regex, col, row = 1) {
  tbs <- c(x$res$tables, if (!is.null(x$posthoc)) x$posthoc$tables)
  nm <- grep(table_regex, names(tbs), value = TRUE)
  if (!length(nm)) return(NA)
  t <- as.data.frame(tbs[[nm[1]]], check.names = FALSE)
  if (is.character(row)) { k <- which(apply(t, 1, function(r) any(trimws(as.character(r)) == row))); if (!length(k)) return(NA); row <- k[1] }
  if (!col %in% names(t)) return(NA)
  t[[col]][row]
}

write_results <- function(file) { utils::write.csv(RESULTS, file, row.names = FALSE, fileEncoding = "UTF-8"); invisible(RESULTS) }
