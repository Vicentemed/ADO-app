# Reporte Word de prueba con historial de preprocesamiento (t apareada sobre datos filtrados)
suppressPackageStartupMessages(library(ggplot2))
for (f in c("00_config.R", "01_registry.R", "02_data.R", "03_assumptions.R", "04_descriptives.R",
            "05_plots.R", "06_report.R")) source(file.path("R", f), encoding = "UTF-8")
for (f in list.files(file.path("R", "analyses"), full.names = TRUE)) source(f, encoding = "UTF-8")
s <- ado_get("wilcox_sr"); raw <- s$example$prep(read_book_data(file.path("data", s$example$file)))
d <- raw[raw$Ed > 20, ]; d <- d[setdiff(names(d), c("Tot8", "Tremed"))]
lg <- c("Filtro (expresión): Ed > 20", "Removidas: Tot8, Tremed")
v <- s$example$vars; a <- s$assumptions(d, v); r <- s$run(d, v, list())
ctx <- list(title = "Reporte: Wilcoxon rangos con signo", analysis = s$label, family = FAMILIES[[s$family]], chapter = s$chapter,
            contexto = contexto_placeholder(), dataset = s$example$file, vars = v, preproc_log = lg,
            n_raw = nrow(raw), n_proc = nrow(d), k_raw = ncol(raw), k_proc = ncol(d),
            t1 = s$descriptives(d, v), checks = a$checks, decision = a$decision, tables = r$tables, posthoc_tables = list(),
            plots = lapply(s$plots, function(pl) pl$fn(d, v, r, list())), plot_labels = sapply(s$plots, `[[`, "label"),
            report_md = s$report(r, d, v, a$checks), notes = r$notes)
out <- file.path(normalizePath("report"), "Ejemplo_reporte_preproc.docx"); build_report(ctx, out)
cat("OK", out, round(file.size(out) / 1024), "KB\n")
