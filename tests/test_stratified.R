# Análisis estratificado: t independiente de MAT dentro de cada PROP (dat24h.xlsx) + reporte Word por estratos
suppressPackageStartupMessages(library(ggplot2))
for (f in c("00_config.R", "01_registry.R", "02_data.R", "03_assumptions.R", "04_descriptives.R",
            "05_plots.R", "06_report.R")) source(file.path("R", f), encoding = "UTF-8")
for (f in list.files(file.path("R", "analyses"), full.names = TRUE)) source(f, encoding = "UTF-8")
ok <- function(cond, msg) cat(if (isTRUE(cond)) "OK   " else "FALLA", msg, "\n")

d <- load_any("C:/Escritorio 2026/CARPETAS/Protocolos/Dr. Raul Valdivia/Tesis Raul/dat24h.xlsx", "dat24h.xlsx")
resp <- grep("VIAB", names(d), value = TRUE, ignore.case = TRUE)[1]
s <- ado_get("t_indep"); v <- list(response = resp, group = "MAT"); st <- "PROP"
L <- split(d, fac(d, st)); names(L) <- paste0(st, " = ", names(L))
ok(length(L) == 4, "4 estratos de PROP")

spec_for <- function(a) if (isTRUE(a$decision$parametric)) s else ado_get(s$alternative)
blocks <- lapply(names(L), function(nm) { dd <- L[[nm]]
  a <- s$assumptions(dd, v); sr <- spec_for(a); r <- sr$run(dd, v, list()); ph <- if (!is.null(sr$posthoc)) sr$posthoc(r, dd, v, list()) else NULL
  cat("   ", nm, "→", sr$id, "| p =", fmt_num(r$stats$p, 3), "\n")
  list(name = nm, n = nrow(dd), analysis = sr$label, t1 = s$descriptives(dd, v), checks = a$checks, decision = a$decision,
       tables = r$tables, posthoc_tables = if (!is.null(ph)) ph$tables else list(),
       plots = lapply(sr$plots, function(pl) pl$fn(dd, v, r, list())), plot_labels = sapply(sr$plots, `[[`, "label"),
       report_md = sr$report(r, dd, v, a$checks), notes = r$notes) })
ok(all(sapply(blocks, function(b) !is.null(b$t1))), "descriptiva en cada estrato")
ok(all(sapply(blocks, function(b) nchar(b$report_md) > 50)), "texto 'cómo reportar' en cada estrato")

ctx <- list(title = paste("Reporte:", s$label, "— estratificado por PROP"), analysis = s$label, family = FAMILIES[[s$family]], chapter = s$chapter,
            contexto = contexto_placeholder(), dataset = "dat24h.xlsx", vars = v, strata = st, preproc_log = character(0),
            n_raw = nrow(d), n_proc = nrow(d), k_raw = ncol(d), k_proc = ncol(d), blocks = blocks)
out <- file.path(normalizePath("report"), "Ejemplo_reporte_estratificado.docx"); build_report(ctx, out)
ok(file.exists(out) && file.size(out) > 20000, paste("docx estratificado generado:", round(file.size(out) / 1024), "KB"))
