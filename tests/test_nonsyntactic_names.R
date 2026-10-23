# Prueba con nombres de columna no sintácticos ("VIAB") — archivo real del usuario
suppressPackageStartupMessages(library(ggplot2))
for (f in c("00_config.R", "01_registry.R", "02_data.R", "03_assumptions.R", "04_descriptives.R",
            "05_plots.R", "06_report.R")) source(file.path("R", f), encoding = "UTF-8")
for (f in list.files(file.path("R", "analyses"), full.names = TRUE)) source(f, encoding = "UTF-8")
ok <- function(cond, msg) cat(if (isTRUE(cond)) "OK   " else "FALLA", msg, "\n")

d <- load_any("C:/Escritorio 2026/CARPETAS/Protocolos/Dr. Raul Valdivia/Tesis Raul/dat24h.xlsx", "dat24h.xlsx")
cat("   nombres:", paste(names(d), collapse = ", "), "\n   renombradas:", paste(paste0(names(attr(d, "renamed")), "→", attr(d, "renamed")), collapse = "; "), "\n")
ok("VIAB" %in% names(d) || "X_VIAB" %in% names(d) || any(grepl("VIAB", names(d))), "columna con símbolo fue saneada")
v <- list(response = grep("VIAB", names(d), value = TRUE)[1], group = "GRUPO")
s <- cdo_get("anova1")
ok(length(validate_roles(d, s, v)) == 0, "roles válidos")
t1 <- s$descriptives(d, v)
ok(!is.null(t1), "table1 con %VIAB | GRUPO se genera")
if (!is.null(t1)) cat("   columnas table1:", length(t1) , "| primer encabezado:", substr(gsub("<[^>]+>", " ", unclass(t1)[[1]]), 1, 60), "\n")
a <- s$assumptions(d, v); ok(length(a$checks) == 2, "supuestos calculados")
s_run <- if (!isTRUE(a$decision$parametric)) cdo_get(s$alternative) else s
cat("   ruta:", s_run$id, "|", a$decision$reason, "\n")
r <- s_run$run(d, v, list(posthoc = "tukey")); ok(is.list(r$stats), paste("run", s_run$id))
ph <- s_run$posthoc(r, d, v, list(posthoc = "tukey")); ok(length(ph$tables) > 0, "post-hoc")
png(tempfile(fileext = ".png")); for (pl in s_run$plots) draw_plot(pl$fn(d, v, r, list())); dev.off(); ok(TRUE, "gráficos")
txt <- s_run$report(r, d, v, a$checks); ok(nchar(txt) > 50, "texto 'cómo reportar'")

# Subanálisis: MAT dentro de cada PROP (t independiente) — misma variable con símbolo
for (p in unique(d$PROP)) { dd <- d[d$PROP == p, ]; s2 <- cdo_get("t_indep"); v2 <- list(response = "VIAB", group = "MAT")
  r2 <- tryCatch(s2$run(dd, v2, list()), error = function(e) NULL); ok(!is.null(r2) && !is.null(s2$descriptives(dd, v2)), paste("subanálisis MAT | PROP =", p)) }

# Regresión lineal con respuesta %VIAB (fórmula con backticks)
s3 <- cdo_get("lineal"); v3 <- list(response = "VIAB", predictors = c("MAT", "PROP"))
r3 <- tryCatch(s3$run(d, v3, list()), error = function(e) { cat("  ", conditionMessage(e), "\n"); NULL }); ok(!is.null(r3), "regresión lineal con %VIAB")
