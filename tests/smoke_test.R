# =============================================================================
# tests/smoke_test.R — Ejecuta CADA análisis del registro sobre su ejemplo del
# libro: descriptiva → supuestos → run → post-hoc → gráficos → texto → (1 docx)
# Uso:  Rscript tests/smoke_test.R   (desde la carpeta ADO-app)
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2) })
for (f in c("00_config.R", "01_registry.R", "02_data.R", "03_assumptions.R", "04_descriptives.R",
            "05_plots.R", "06_report.R")) source(file.path("R", f), encoding = "UTF-8")
for (f in list.files(file.path("R", "analyses"), full.names = TRUE)) source(f, encoding = "UTF-8")

specs <- ado_list()
cat("Análisis registrados:", length(specs), "\n\n")
results <- list(); first_ctx <- NULL

for (s in specs) {
  t0 <- Sys.time(); status <- "OK"; msg <- ""; nplots <- 0
  tryCatch({
    d <- read_book_data(file.path("data", s$example$file))
    if (is.function(s$example$prep)) d <- s$example$prep(d)
    v <- s$example$vars
    pr <- validate_roles(d, s, v); if (length(pr)) stop("roles: ", paste(pr, collapse = "; "))
    t1 <- s$descriptives(d, v)
    a  <- s$assumptions(d, v)
    stopifnot(is.list(a$checks), is.list(a$decision))
    # ruta según la decisión del libro
    s_run <- if (!is.null(s$alternative) && !isTRUE(a$decision$parametric)) ado_get(s$alternative) else s
    opts <- lapply(s_run$opts_ui, `[[`, "value")
    r <- s_run$run(d, v, opts)
    stopifnot(is.list(r$stats), is.list(r$tables))
    ph <- if (!is.null(s_run$posthoc)) s_run$posthoc(r, d, v, opts) else NULL
    # gráficos → png temporal
    png(tempfile(fileext = ".png"), width = 900, height = 600)
    for (pl in s_run$plots) { p <- pl$fn(d, v, r, list()); draw_plot(p); nplots <- nplots + 1 }
    dev.off()
    txt <- s_run$report(r, d, v, a$checks)
    stopifnot(is.character(txt), nchar(txt) > 50)
    if (is.null(first_ctx)) first_ctx <- list(
      title = paste("Reporte:", s_run$label), analysis = s_run$label, family = FAMILIES[[s_run$family]], chapter = s_run$chapter,
      contexto = contexto_placeholder(), dataset = s$example$file, vars = v, t1 = t1, checks = a$checks, decision = a$decision,
      tables = r$tables, posthoc_tables = if (!is.null(ph)) ph$tables else list(),
      plots = lapply(s_run$plots, function(pl) pl$fn(d, v, r, list())), plot_labels = sapply(s_run$plots, `[[`, "label"),
      report_md = txt, notes = r$notes)
    msg <- paste0(if (!identical(s_run$id, s$id)) paste0("→ alternativa: ", s_run$id, " | ") else "", "gráficos: ", nplots, " | tablas: ", length(r$tables))
  }, error = function(e) { status <<- "FALLA"; msg <<- conditionMessage(e); try(dev.off(), silent = TRUE) })
  results[[s$id]] <- data.frame(id = s$id, familia = s$family, estado = status, detalle = substr(msg, 1, 120), seg = round(as.numeric(Sys.time() - t0, units = "secs"), 1))
  cat(sprintf("%-14s %-14s %-6s %s\n", s$id, s$family, status, msg))
}

res <- do.call(rbind, results)
cat("\n== RESUMEN ==\n"); print(table(res$estado))
if (any(res$estado == "FALLA")) { cat("\nFALLAS:\n"); print(res[res$estado == "FALLA", c("id", "detalle")], row.names = FALSE) }

# Reporte Word con el primer análisis exitoso
cat("\n== Reporte Word de prueba ==\n")
tryCatch({ out <- file.path(tempdir(), "ado_smoke.docx"); build_report(first_ctx, out); cat("docx generado:", out, "|", round(file.size(out) / 1024), "KB\n") },
         error = function(e) cat("FALLA docx:", conditionMessage(e), "\n"))
