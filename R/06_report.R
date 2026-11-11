# =============================================================================
# 06_report.R - Generación del reporte Word (rmarkdown → docx)
# ctx: list(title, analysis, family, chapter, contexto, dataset, vars,
#           t1 (table1), checks (list), decision (list), results_md, tables (list df),
#           posthoc_tables (list df), plots (list ggplot/función), report_md, notes)
# =============================================================================

build_report <- function(ctx, out_docx, template = file.path("report", "report_template.Rmd")) {
  stopifnot(file.exists(template))
  # copia temporal del template para no escribir junto al original
  tmpdir <- tempfile("ado_rep_"); dir.create(tmpdir)
  tpl <- file.path(tmpdir, "report.Rmd"); file.copy(template, tpl)
  rmarkdown::render(tpl, output_file = out_docx, output_dir = dirname(out_docx),
                    params = list(ctx = ctx), envir = new.env(parent = globalenv()), quiet = TRUE)
  out_docx
}

# Texto estándar que abre el reporte: el usuario debe describir el contexto
contexto_placeholder <- function() {
  paste0("**Contexto del estudio (completar por el usuario):** ",
         CTX("Describa aquí el objetivo del análisis, la población, el diseño, el origen de los datos y el significado de cada variable"),
         ". El texto de resultados que sigue fue prellenado automáticamente con los valores calculados; ",
         "sustituya los corchetes por el contexto de sus datos y revise la interpretación clínica.")
}

df_to_md <- function(df, digits = 3) {
  if (is.null(df) || !nrow(df)) return("")
  num <- sapply(df, is.numeric)
  df[num] <- lapply(df[num], function(x) round(x, digits))
  paste(knitr::kable(df, format = "markdown"), collapse = "\n")
}

# ---- Texto del reporte sin referencias al libro -------------------------------------
# El reporte Word solo menciona el libro en el subtítulo y en la nota final; en el cuerpo
# se eliminan las remisiones a capítulos y las frases del tipo «según el libro».
sin_libro <- function(x) {
  if (is.null(x)) return(x)
  one <- function(t) {
    t <- as.character(t); if (is.na(t) || !nzchar(t)) return(t)
    rp <- list(
      c("\\s*\\((Tabla [0-9.]+ del libro|regla del libro|formato del libro|libro|como (en )?el libro|como en el reporte del libro)\\)", ""),
      c("\\s*\\([Cc]aps?\\. ?[0-9][^()]*\\)", ""),                                   # (cap. 15.1.6) · (caps. 19.8 y 27.11)
      c("\\s*\\([Cc]ap\\. ?[0-9.]+ lo muestra a manera de ejemplo\\)", ""),
      c(",? ?(donde se denomina[^)]*)", ""),
      c("[,;]? ?[Cc]aps?\\. ?[0-9][0-9.]*( y [0-9][0-9.]*)?", ""),                    # …, cap. 18.2
      c("Criterio\\s*:", "Criterio:"), c("Nota\\s*:", "Nota:"),
      c("[,;] como (en )?el libro( \\([^)]*\\))?", ""), c("[,;] como en el reporte del libro", ""), c("Como en el libro, ", ""),
      c(": según el libro, ", ": "), c("Según el libro, ([a-záéíóúñ])", "\\U\\1"), c("[Ss]egún el libro, ", ""), c(" según el libro", ""),
      c("[Ee]l libro (la )?recomienda", "se \\1recomienda"), c("[Ee]l libro sugiere", "se sugiere"), c("[Ee]l libro exige", "se exige"),
      c("[Ee]l libro recurre entonces al", "se recurre entonces al"), c("[Ee]l libro usa", "Se usa"), c("[Ee]l libro centra", "Se centra"),
      c("[Ee]l libro especifica los factores al inicio \\(as\\.factor\\)", "Los factores se especifican al inicio (as.factor)"),
      c("[Ee]l libro convierte las covariables categóricas a factor \\(factor\\(\\)\\)", "Las covariables categóricas se convierten a factor"),
      c("la regla del libro pide", "la regla práctica pide"), c("(Flujo|Procedimiento|Alternativa) del libro:", "\\1:"),
      c(" para la lectura visual del libro", " para la lectura visual"), c(" de las plantillas de reporte del libro", ""),
      c(" no tiene alternativa no paramétrica en el libro", " no tiene alternativa no paramétrica"),
      c("\\s*\\(regla del libro\\)", ""), c(" del libro", ""), c(" en el libro", ""),
      c("\\s\\(\\s*\\)", ""), c("\\(\\s*[,;]\\s*", "("), c("\\s+([.,;:])", "\\1"), c("  +", " "))   # paréntesis vacíos sueltos; no llamadas como lme()
    for (r in rp) t <- gsub(r[1], r[2], t, perl = TRUE)
    t <- gsub("(^|\\. )se (la )?(recomienda|sugiere|exige|recurre)", "\\1Se \\2\\3", t, perl = TRUE)   # mayúscula solo al inicio de oración
    trimws(t)
  }
  if (is.list(x)) return(lapply(x, sin_libro))
  nm <- names(x); out <- vapply(x, one, NA_character_, USE.NAMES = FALSE); names(out) <- nm; out
}

# ---- Paquetes de R que usa un análisis (a partir de su propio código) ------------------
analysis_packages <- function(specs) {
  seen <- character(0); pk <- character(0)
  scan <- function(txt, depth) {
    pk <<- union(pk, unique(regmatches(txt, gregexpr("[A-Za-z][A-Za-z0-9.]*(?=::)", txt, perl = TRUE))[[1]]))
    if (depth >= 4) return(invisible())
    sy <- setdiff(unique(regmatches(txt, gregexpr("[.A-Za-z_][.A-Za-z0-9_]*(?=\\()", txt, perl = TRUE))[[1]]), seen)
    for (f in sy) if (exists(f, envir = globalenv(), inherits = FALSE) && is.function(get(f, envir = globalenv()))) {
      seen <<- c(seen, f); scan(paste(deparse(get(f, envir = globalenv())), collapse = "\n"), depth + 1) }
  }
  for (s in specs) scan(paste(unlist(lapply(s[c("descriptives", "assumptions", "run", "posthoc", "report", "plots", "qq")], function(z) deparse(z))), collapse = "\n"), 0)
  pk <- union(pk, c("table1", "car"))                                                # tabla descriptiva y gráfico Q-Q
  pk <- setdiff(pk, c("stats", "utils", "graphics", "grDevices", "tools", "methods", "base", "shiny", "bslib", "htmltools", "DT"))
  pk <- sort(pk[vapply(pk, requireNamespace, TRUE, quietly = TRUE)])
  setNames(vapply(pk, function(p) as.character(utils::packageVersion(p)), ""), pk)
}
report_session <- function(specs) {
  inf <- c("rmarkdown", "knitr", "flextable")
  list(r = R.version.string, platform = R.version$platform, pkgs = analysis_packages(specs),
       infra = setNames(vapply(inf, function(p) as.character(utils::packageVersion(p)), ""), inf))
}
