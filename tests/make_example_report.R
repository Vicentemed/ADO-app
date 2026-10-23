# Genera un reporte Word de ejemplo (ANOVA de una vía, datos del libro) en report/
# Uso: Rscript tests/make_example_report.R   (desde la carpeta CDO-app)
suppressPackageStartupMessages(library(ggplot2))
for (f in c("00_config.R", "01_registry.R", "02_data.R", "03_assumptions.R", "04_descriptives.R",
            "05_plots.R", "06_report.R")) source(file.path("R", f), encoding = "UTF-8")
for (f in list.files(file.path("R", "analyses"), full.names = TRUE)) source(f, encoding = "UTF-8")

s  <- cdo_get("anova1")
d  <- s$example$prep(read_book_data(file.path("data", s$example$file)))
v  <- s$example$vars
t1 <- s$descriptives(d, v)
a  <- s$assumptions(d, v)
r  <- s$run(d, v, list(posthoc = "tukey"))
ph <- s$posthoc(r, d, v, list(posthoc = "tukey"))
ctx <- list(title = "Reporte: ANOVA de una vía", analysis = s$label, family = FAMILIES[[s$family]], chapter = s$chapter,
            contexto = contexto_placeholder(), dataset = s$example$file, vars = v, t1 = t1,
            checks = a$checks, decision = a$decision, tables = r$tables, posthoc_tables = ph$tables,
            plots = lapply(s$plots, function(pl) pl$fn(d, v, r, list())), plot_labels = sapply(s$plots, `[[`, "label"),
            report_md = s$report(r, d, v, a$checks), notes = r$notes)
out <- file.path(normalizePath("report"), "Ejemplo_reporte_ANOVA.docx")
build_report(ctx, out)
cat("OK", out, round(file.size(out) / 1024), "KB\n")
