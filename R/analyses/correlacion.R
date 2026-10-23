# =============================================================================
# analyses/correlacion.R — Pearson (alt: Spearman/Kendall), parcial, PPS
# Reglas del libro (cap. Correlaciones): Pearson requiere normalidad de ambas
# variables (Q-Q/Shapiro); si no → Spearman; con muchos empates/ordinal → Kendall.
# =============================================================================

efecto_cor <- function(r) { r <- abs(r); if (is.na(r)) "…" else if (r < 0.1) "trivial" else if (r < 0.3) "pequeño" else if (r < 0.5) "mediano" else "grande" }

.cor_run <- function(d, v, method) {
  dd <- data.frame(x = num(d, v$x), y = num(d, v$y)); dd <- dd[complete.cases(dd), ]
  ct <- suppressWarnings(cor.test(dd$x, dd$y, alternative = "two.sided", method = method))
  est <- unname(ct$estimate); ci <- if (!is.null(ct$conf.int)) ct$conf.int else c(NA, NA)
  sym <- switch(method, pearson = "r", spearman = "ρ", kendall = "τ")
  list(main = ct,
       stats = list(est = est, p = ct$p.value, lo = ci[1], hi = ci[2], n = nrow(dd), sym = sym, r2 = if (method == "pearson") est^2 else NA,
                    stat = unname(ct$statistic), df = if (!is.null(ct$parameter)) unname(ct$parameter) else NA),
       tables = list("Resultado de cor.test()" = data.frame(Coeficiente = sym, Estimación = est, `IC 95% inf` = ci[1], `IC 95% sup` = ci[2], p = ct$p.value, n = nrow(dd), check.names = FALSE)),
       notes = if (method == "pearson") "R² = r² indica la proporción de varianza compartida." else "Spearman/Kendall no proporcionan IC directamente (el libro sugiere bootstrap).")
}
.cor_report <- function(res, d, v, checks, method) {
  s <- res$stats; nm <- c(pearson = "Pearson", spearman = "Spearman", kendall = "Kendall")[method]
  swx <- checks[[1]]; swy <- checks[[2]]
  paste0("La normalidad de ", CTX(v$x), " (Shapiro-Wilk, ", .sw_txt(swx), ") y de ", CTX(v$y), " (", .sw_txt(swy), ") se evaluó junto con el gráfico Q-Q",
         if (isTRUE(swx$met) && isTRUE(swy$met)) ", sin evidencia de desviación de la normalidad" else ", con evidencia de desviación de la normalidad en al menos una variable",
         "; por ello se aplicó la correlación de ", nm, ". La relación entre ", CTX(v$x), " y ", CTX(v$y), " fue ",
         if (s$est > 0) "positiva" else "negativa", " y ", if (s$p < ALPHA) "estadísticamente significativa" else "no significativa",
         " (", s$sym, " = ", fmt_num(s$est, 3), if (!is.na(s$lo)) paste0(", IC 95 %: ", fmt_num(s$lo, 3), " a ", fmt_num(s$hi, 3)) else "", ", p ", fmt_p(s$p), ", n = ", s$n, "), un efecto ", efecto_cor(s$est),
         if (!is.na(s$r2)) paste0("; R² = ", fmt_num(s$r2, 3), " (", round(100 * s$r2), " % de varianza compartida)") else "", " [Cohen, 1988; Field, 2012].")
}
.cor_roles <- list(x = role("Variable X (numérica)", "numeric"), y = role("Variable Y (numérica)", "numeric"))
.cor_assump <- function(d, v) list(check_normality(num(d, v$x), v$x), check_normality(num(d, v$y), v$y))
.cor_plots <- function(jit) list(
  scatter = list(label = "Dispersión + recta de regresión", fn = function(d, v, res, cos) plot_scatter(d, v$x, v$y, cos, r_label = paste0(res$stats$sym, " = ", fmt_num(res$stats$est, 2)), jitter = jit)),
  qqx = list(label = "Q-Q de X", fn = function(d, v, res, cos) plot_qq(num(d, v$x), paste("Q-Q:", v$x), cos)),
  qqy = list(label = "Q-Q de Y", fn = function(d, v, res, cos) plot_qq(num(d, v$y), paste("Q-Q:", v$y), cos)))

cdo_register(list(
  id = "pearson", family = "correlacion", order = 1, label = "Correlación de Pearson",
  description = "Relación lineal entre dos variables numéricas normales. Alternativa: Spearman.",
  chapter = "Correlaciones", roles = .cor_roles,
  example = list(file = "DXKTw.txt", vars = list(x = "Imc", y = "Tqx"), note = "IMC y tiempo quirúrgico (min)."),
  alternative = "spearman",
  descriptives = function(d, v) desc_table1(d, c(v$x, v$y)),
  assumptions = function(d, v) { ch <- .cor_assump(d, v); list(checks = ch, decision = decide_path(ch, TRUE)) },
  run = function(d, v, opts) .cor_run(d, v, "pearson"), posthoc = NULL, plots = .cor_plots(FALSE),
  report = function(res, d, v, checks) .cor_report(res, d, v, checks, "pearson")
))
cdo_register(list(
  id = "spearman", family = "correlacion", order = 2, label = "Correlación de Spearman (ρ)",
  description = "Alternativa no paramétrica de Pearson: variables no normales u ordinales.",
  chapter = "Correlaciones", roles = .cor_roles,
  example = list(file = "datcron.txt", vars = list(x = "Dental_Age", y = "Cron_Age"), note = "Edad dental y cronológica (meses)."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, c(v$x, v$y)),
  assumptions = function(d, v) { ch <- .cor_assump(d, v); list(checks = ch, decision = list(parametric = FALSE, all_met = FALSE, violated = character(0), reason = "No paramétrica.")) },
  run = function(d, v, opts) .cor_run(d, v, "spearman"), posthoc = NULL, plots = .cor_plots(TRUE),
  report = function(res, d, v, checks) .cor_report(res, d, v, checks, "spearman")
))
cdo_register(list(
  id = "kendall", family = "correlacion", order = 3, label = "Tau de Kendall (τ)",
  description = "Recomendada con datos ordinales de pocos niveles o muchos empates.",
  chapter = "Correlaciones", roles = .cor_roles,
  example = list(file = "datcron.txt", vars = list(x = "Cron_Age", y = "Bone_Maturation"), note = "Edad cronológica y maduración ósea (ordinal 1-5)."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, c(v$x, v$y)),
  assumptions = function(d, v) { ch <- .cor_assump(d, v); list(checks = ch, decision = list(parametric = FALSE, all_met = FALSE, violated = character(0), reason = "No paramétrica.")) },
  run = function(d, v, opts) .cor_run(d, v, "kendall"), posthoc = NULL, plots = .cor_plots(TRUE),
  report = function(res, d, v, checks) .cor_report(res, d, v, checks, "kendall")
))

cdo_register(list(
  id = "parcial", family = "correlacion", order = 4, label = "Correlación parcial (controlando Z)",
  description = "Correlación entre X e Y eliminando la influencia de una tercera variable Z (ggm::pcor).",
  chapter = "Correlaciones",
  roles = list(x = role("Variable X", "numeric"), y = role("Variable Y", "numeric"), z = role("Variable de control Z", "numeric")),
  example = list(file = "DXKTw.txt", vars = list(x = "Tqx", y = "Imc", z = "Ed"), note = "Tiempo quirúrgico e IMC controlando por edad."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, c(v$x, v$y, v$z)),
  assumptions = function(d, v) { ch <- .cor_assump(d, v); list(checks = ch, decision = decide_path(ch, FALSE)) },
  run = function(d, v, opts) {
    dd <- data.frame(x = num(d, v$x), y = num(d, v$y), z = num(d, v$z)); dd <- dd[complete.cases(dd), ]
    S <- var(dd); rp <- ggm::pcor(c("x", "y", "z"), S)
    biv <- cor.test(dd$x, dd$y); c1 <- cor.test(dd$x, dd$z); c2 <- cor.test(dd$y, dd$z)
    n <- nrow(dd); tval <- rp * sqrt((n - 3) / (1 - rp^2)); p_partial <- 2 * pt(-abs(tval), n - 3)
    list(main = rp,
         stats = list(r_partial = rp, p_partial = p_partial, r_biv = unname(biv$estimate), p_biv = biv$p.value, r_xz = unname(c1$estimate), r_yz = unname(c2$estimate), n = n),
         tables = list("Correlaciones" = data.frame(Relación = c(paste(v$x, "~", v$y, "(bivariada)"), paste(v$x, "~", v$y, "| ", v$z, "(parcial)"), paste(v$x, "~", v$z), paste(v$y, "~", v$z)),
                                                    r = c(unname(biv$estimate), rp, unname(c1$estimate), unname(c2$estimate)), p = c(biv$p.value, p_partial, c1$p.value, c2$p.value))),
         notes = "pcor() no proporciona IC; el p de la parcial se obtiene con t = r·√((n−3)/(1−r²)).")
  },
  posthoc = NULL,
  plots = list(scatter = list(label = "Dispersión X-Y", fn = function(d, v, res, cos) plot_scatter(d, v$x, v$y, cos, r_label = paste0("r parcial = ", fmt_num(res$stats$r_partial, 2))))),
  report = function(res, d, v, checks) {
    s <- res$stats
    paste0("La correlación bivariada entre ", CTX(v$x), " y ", CTX(v$y), " fue r = ", fmt_num(s$r_biv, 3), " (p ", fmt_p(s$p_biv), "). Al controlar por ", CTX(v$z),
           " —posible variable de confusión (", v$x, " ~ ", v$z, ": r = ", fmt_num(s$r_xz, 2), "; ", v$y, " ~ ", v$z, ": r = ", fmt_num(s$r_yz, 2), ")—, la correlación parcial fue r = ",
           fmt_num(s$r_partial, 3), " (p ", fmt_p(s$p_partial), ", n = ", s$n, "), ", if (abs(s$r_partial - s$r_biv) > 0.1) "lo que modifica sustancialmente la relación observada." else "similar a la bivariada, por lo que la tercera variable no explica la relación.",
           " Siempre debe considerarse la posible influencia de variables de confusión al calcular correlaciones [Field, 2012].")
  }
))

cdo_register(list(
  id = "pps", family = "correlacion", order = 5, label = "Predictive Power Score (PPS)",
  description = "Capacidad predictiva (asimétrica, no lineal) de varias variables sobre una respuesta (ppsr).",
  chapter = "Correlaciones",
  roles = list(target = role("Variable objetivo (respuesta)", "any"), predictors = role("Predictores candidatos", "any", multiple = TRUE, min = 1)),
  example = list(file = "DXKTw.txt", vars = list(target = "Tqx", predictors = c("Classcx", "Imc", "Ed", "Sexo", "Grupo")), note = "¿Qué variables predicen el tiempo quirúrgico?"),
  alternative = NULL,
  descriptives = function(d, v) NULL,
  assumptions = function(d, v) list(checks = list(check_obj("info", "PPS", NA, "El PPS no exige supuestos distribucionales; se estima por validación cruzada (árbol de decisión), por lo que puede variar ligeramente entre ejecuciones.")), decision = list(parametric = TRUE, all_met = TRUE, violated = character(0), reason = "Sin supuestos.")),
  run = function(d, v, opts) {
    dd <- d[c(v$target, v$predictors)]; for (nm in names(dd)) if (is.character(dd[[nm]])) dd[[nm]] <- factor(dd[[nm]])
    dd <- dd[complete.cases(dd), ]
    sc <- do.call(rbind, lapply(v$predictors, function(x) { r <- ppsr::score(dd, x = x, y = v$target); data.frame(x = x, y = v$target, pps = r$pps) }))
    sc <- sc[order(-sc$pps), ]
    list(main = sc, stats = list(top = sc$x[1], top_pps = sc$pps[1]),
         tables = list("PPS de cada predictor sobre la respuesta" = sc), notes = "PPS 0 = sin poder predictivo; 1 = predicción perfecta. Es asimétrico (x→y ≠ y→x).")
  },
  posthoc = NULL,
  plots = list(pps = list(label = "PPS por predictor", fn = function(d, v, res, cos) {
    cos <- cos_merge(cos); sc <- res$main; sc$x <- factor(sc$x, levels = rev(sc$x))
    p <- ggplot(sc, aes(x, pps)) + geom_col(fill = cos$colour[1], alpha = 0.85) + coord_flip() + ylim(0, 1)
    .finish(p, cos, "Predictor", paste("PPS →", v$target)) })),
  report = function(res, d, v, checks) {
    s <- res$stats; sc <- res$main
    paste0("Se calculó el Predictive Power Score (PPS) de ", length(v$predictors), " variables sobre ", CTX(v$target), ". El mayor poder predictivo correspondió a ", CTX(s$top), " (PPS = ", fmt_num(s$top_pps), "); ",
           paste(paste0(sc$x, " = ", fmt_num(sc$pps)), collapse = ", "), ". El PPS complementa a la correlación: informa de la capacidad predictiva —incluso no lineal— sin indicar dirección del efecto; ",
           "para un IC 95 % puede emplearse remuestreo bootstrap [Wetschoreck et al., 2021].")
  }
))
