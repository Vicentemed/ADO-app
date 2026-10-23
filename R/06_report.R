# =============================================================================
# 06_report.R — Generación del reporte Word (rmarkdown → docx)
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
