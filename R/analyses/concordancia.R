# =============================================================================
# analyses/concordancia.R — Kappa (nominal/ordinal), ICC, Bland-Altman (alt: CCC)
# Regla del libro: categórica → kappa (ordinal → ponderado); continua con
# varios evaluadores → ICC; dos métodos continuos → Bland-Altman si las
# diferencias son normales, CCC de Lin si no (o ambos).
# =============================================================================

landis <- function(k) { if (is.na(k)) return("…"); if (k < 0) "sin acuerdo" else if (k <= 0.20) "ligero" else if (k <= 0.40) "razonable" else if (k <= 0.60) "moderado" else if (k <= 0.80) "sustancial" else "casi perfecto" }
koo_li <- function(i) { if (is.na(i)) return("…"); if (i < 0.5) "pobre" else if (i < 0.75) "moderada" else if (i < 0.9) "buena" else "excelente" }
ccc_band <- function(r) { if (is.na(r)) return("…"); if (r < 0.90) "pobre" else if (r < 0.95) "moderada" else if (r < 0.99) "sustancial" else "casi perfecta" }

ado_register(list(
  id = "kappa", family = "concordancia", order = 1, label = "Kappa de Cohen (2 evaluadores; ponderado si ordinal)",
  description = "Acuerdo entre dos evaluadores en una escala categórica. Con > 2 evaluadores calcula además Fleiss.",
  chapter = "Kappa en R",
  roles = list(raters = role("Columnas de los evaluadores (≥ 2)", "factor", multiple = TRUE, min = 2)),
  example = list(file = "PGMax.txt", vars = list(raters = c("Obs2", "Obs3")), note = "Clasificación ordinal (AI < AII < BII < BIII) de 20 radiografías por 2 observadores."),
  alternative = NULL,
  descriptives = function(d, v) NULL,
  assumptions = function(d, v) {
    lv <- unique(unlist(lapply(v$raters, function(r) as.character(d[[r]]))))
    list(checks = list(check_obj("escala", "Tipo de escala", NA, paste0("Categorías observadas: ", paste(sort(lv), collapse = ", "), ". Si la escala es ORDINAL, el libro recomienda el kappa ponderado (lineal o cuadrático); si es nominal, el no ponderado. Recuerde la paradoja de prevalencia: con prevalencias extremas el kappa puede ser artificialmente bajo."))),
         decision = list(parametric = TRUE, all_met = TRUE, violated = character(0), reason = ""))
  },
  run = function(d, v, opts) {
    lv <- sort(unique(unlist(lapply(v$raters, function(r) as.character(d[[r]])))))
    if (!is.null(opts$levels) && nzchar(opts$levels)) { l2 <- trimws(strsplit(opts$levels, ",")[[1]]); if (all(lv %in% l2)) lv <- l2 }
    M <- as.data.frame(lapply(d[v$raters], function(x) factor(as.character(x), levels = lv))); M <- M[complete.cases(M), ]
    k_unw <- irr::kappa2(M[, 1:2], "unweighted"); k_lin <- irr::kappa2(M[, 1:2], "equal"); k_sq <- irr::kappa2(M[, 1:2], "squared")
    fl <- if (ncol(M) > 2) irr::kappam.fleiss(M) else NULL
    ci <- tryCatch({ ck <- psych::cohen.kappa(as.matrix(M[, 1:2])); ck$confid }, error = function(e) NULL)
    tab <- table(M[, 1], M[, 2]); names(dimnames(tab)) <- v$raters[1:2]
    st <- list(k_unw = k_unw$value, z_unw = k_unw$statistic, p_unw = k_unw$p.value, k_lin = k_lin$value, z_lin = k_lin$statistic, p_lin = k_lin$p.value,
               k_sq = k_sq$value, z_sq = k_sq$statistic, p_sq = k_sq$p.value, n = nrow(M), r1 = v$raters[1], r2 = v$raters[2], k_fleiss = if (!is.null(fl)) fl$value else NA, nr = ncol(M))
    tabs <- list("Tabla de acuerdo (evaluador 1 × evaluador 2)" = .tab_df(tab),
                 "Kappa de Cohen" = data.frame(Versión = c("No ponderado", "Ponderado lineal", "Ponderado cuadrático"), κ = c(st$k_unw, st$k_lin, st$k_sq), z = c(st$z_unw, st$z_lin, st$z_sq), p = c(st$p_unw, st$p_lin, st$p_sq), Interpretación = sapply(c(st$k_unw, st$k_lin, st$k_sq), landis)))
    if (!is.null(ci)) tabs[["IC 95 % (psych::cohen.kappa)"]] <- data.frame(Tipo = rownames(ci), ci, check.names = FALSE, row.names = NULL)
    if (!is.null(fl)) tabs[["Kappa de Fleiss (todos los evaluadores)"]] <- data.frame(κ = fl$value, z = fl$statistic, p = fl$p.value, Evaluadores = ncol(M), Interpretación = landis(fl$value))
    list(main = k_lin, stats = st, tables = tabs, notes = "Escala de Landis y Koch. El kappa de Fleiss es no ponderado (nominal): subestima el acuerdo en escalas ordinales.")
  },
  posthoc = NULL,
  plots = list(matriz = list(label = "Matriz de acuerdo (mapa de calor)", fn = function(d, v, res, cos) {
    cos <- cos_merge(cos); M <- d[v$raters[1:2]]; M <- M[complete.cases(M), ]
    df <- as.data.frame(table(factor(as.character(M[, 1])), factor(as.character(M[, 2])))); names(df) <- c("r1", "r2", "n")
    p <- ggplot(df, aes(r2, r1, fill = n)) + geom_tile(colour = "white") + geom_text(aes(label = n), size = 4) + scale_fill_gradient(low = "white", high = cos$colour[1])
    cos$legend <- TRUE; .finish(p, cos, v$raters[2], v$raters[1]) })),
  opts_ui = list(levels = list(label = "Orden de categorías (ordinal), separadas por coma", type = "text", value = "")),
  report = function(res, d, v, checks) {
    s <- res$stats
    paste0("La concordancia entre los evaluadores ", CTX(s$r1), " y ", CTX(s$r2), " para ", CTX("la clasificación"), " se evaluó con el Kappa de Cohen (n = ", s$n, "). ",
           "El kappa no ponderado fue κ = ", fmt_num(s$k_unw), " (z = ", fmt_num(s$z_unw), "; p ", fmt_p(s$p_unw), "), acuerdo ", landis(s$k_unw), "; ",
           "si la escala es ordinal procede la versión ponderada: lineal κ = ", fmt_num(s$k_lin), " (z = ", fmt_num(s$z_lin), "; p ", fmt_p(s$p_lin), ") y cuadrático κ = ", fmt_num(s$k_sq), " (z = ", fmt_num(s$z_sq), "; p ", fmt_p(s$p_sq), "), acuerdo ", landis(s$k_sq), " (escala de Landis y Koch)",
           if (!is.na(s$k_fleiss)) paste0(". Con los ", s$nr, " evaluadores, el kappa de Fleiss fue κ = ", fmt_num(s$k_fleiss)) else "", " [Koo y Li, 2016; Heo, 2008; Landis y Koch, 1977].")
  }
))

ado_register(list(
  id = "icc", family = "concordancia", order = 2, label = "Coeficiente de correlación intraclase (ICC)",
  description = "Acuerdo entre ≥ 2 evaluadores en una variable continua (irr::icc, modelos de Shrout-Fleiss).",
  chapter = "ICC en R",
  roles = list(raters = role("Columnas de los evaluadores (≥ 2, numéricas)", "numeric", multiple = TRUE, min = 2)),
  example = list(file = "Concor.txt", vars = list(raters = c("A", "B", "C")), note = "Tres evaluadores, escala continua."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, v$raters),
  assumptions = function(d, v) {
    M <- d[v$raters]; M <- M[complete.cases(M), ]; vs <- sapply(M, var)
    list(checks = list(check_obj("varianzas", "Varianzas por evaluador", NA, paste0("Varianzas: ", paste(paste0(names(vs), " = ", fmt_num(vs)), collapse = "; "), ". Defina de antemano el modelo: consistencia ICC(3,1) para calibrar a evaluadores fijos; acuerdo absoluto ICC(2,1) para generalizar a otros evaluadores."))),
         decision = list(parametric = TRUE, all_met = TRUE, violated = character(0), reason = ""))
  },
  run = function(d, v, opts) {
    M <- as.data.frame(lapply(d[v$raters], as.numeric)); M <- M[complete.cases(M), ]
    unit <- if (is.null(opts$unit)) "single" else opts$unit
    ic <- irr::icc(M, model = "twoway", type = "consistency", unit = unit); ia <- irr::icc(M, model = "twoway", type = "agreement", unit = unit)
    i1 <- irr::icc(M, model = "oneway", unit = unit)
    st <- list(icc_c = ic$value, lo_c = ic$lbound, hi_c = ic$ubound, p_c = ic$p.value, icc_a = ia$value, lo_a = ia$lbound, hi_a = ia$ubound, p_a = ia$p.value, icc_1 = i1$value, n = nrow(M), k = ncol(M), unit = unit)
    list(main = ic, stats = st,
         tables = list("ICC (irr::icc, twoway)" = data.frame(Modelo = c("ICC(1) un factor", "ICC consistencia (3,1)", "ICC acuerdo absoluto (2,1)"), ICC = c(i1$value, ic$value, ia$value), `IC inf` = c(i1$lbound, ic$lbound, ia$lbound), `IC sup` = c(i1$ubound, ic$ubound, ia$ubound), p = c(i1$p.value, ic$p.value, ia$p.value), Interpretación = sapply(c(i1$value, ic$value, ia$value), koo_li), check.names = FALSE)),
         notes = "Interpretación de Koo y Li: < 0.50 pobre, 0.50–0.75 moderada, 0.75–0.90 buena, > 0.90 excelente.")
  },
  posthoc = NULL,
  plots = list(pares = list(label = "Dispersión por pares con línea de identidad", fn = function(d, v, res, cos) {
    cos <- cos_merge(cos); M <- as.data.frame(lapply(d[v$raters], as.numeric)); M <- M[complete.cases(M), ]; prs <- combn(names(M), 2, simplify = FALSE)
    function() { op <- par(mfrow = c(1, length(prs))); on.exit(par(op))
      for (pr in prs) { plot(M[[pr[1]]], M[[pr[2]]], pch = 16, col = cos$colour[1], xlab = pr[1], ylab = pr[2], main = paste(pr, collapse = " vs "))
        abline(0, 1, col = cos$colour[2], lty = 2, lwd = 2); abline(lm(M[[pr[2]]] ~ M[[pr[1]]]), col = cos$colour[1], lwd = 1.5) } } })),
  opts_ui = list(unit = list(label = "Unidad", type = "select", choices = c("Medición única (single)" = "single", "Promedio de k mediciones (average)" = "average"), value = "single")),
  report = function(res, d, v, checks) {
    s <- res$stats
    paste0("La concordancia entre los ", s$k, " evaluadores para ", CTX("la variable medida"), " se evaluó con el ICC (n = ", s$n, "). ", CTX("Como el objetivo es calibrar a los evaluadores específicos del estudio (fijos)"),
           ", se reporta el ICC de consistencia —ICC(3,1)— que fue de ", fmt_num(s$icc_c, 3), " (IC 95 %: ", fmt_num(s$lo_c), " a ", fmt_num(s$hi_c), "; p ", fmt_p(s$p_c), "), una confiabilidad ", koo_li(s$icc_c), " según Koo y Li. ",
           "El ICC de acuerdo absoluto —ICC(2,1)— fue ", fmt_num(s$icc_a, 3), if (abs(s$icc_a - s$icc_c) < 0.05) ", muy similar, lo que indica un sesgo sistemático pequeño entre evaluadores" else ", menor, lo que sugiere un sesgo sistemático entre evaluadores", " [Koo y Li, 2016].")
  }
))

.ba_stats <- function(m, r) { ok <- !is.na(m) & !is.na(r); m <- m[ok]; r <- r[ok]; dif <- m - r; prom <- (m + r) / 2
  bias <- mean(dif); sd_d <- sd(dif); n <- length(dif); se <- sd_d / sqrt(n); tc <- qt(0.975, n - 1)
  list(bias = bias, sd = sd_d, loa_l = bias - 1.96 * sd_d, loa_u = bias + 1.96 * sd_d, ci_lo = bias - tc * se, ci_hi = bias + tc * se, ref = mean(prom), n = n, dif = dif) }

ado_register(list(
  id = "bland_altman", family = "concordancia", order = 3, label = "Bland-Altman (dos métodos continuos)",
  description = "Sesgo y límites de acuerdo entre un método y una referencia. Alternativa: CCC de Lin si las diferencias no son normales.",
  chapter = "Bland-Altman en R",
  roles = list(method = role("Método a evaluar (numérica)", "numeric"), reference = role("Referencia / gold standard (numérica)", "numeric")),
  example = list(file = "datBF.txt", vars = list(method = "Inter", reference = "Basal"), note = "Dos mediciones continuas del mismo sujeto (ejemplo ilustrativo)."),
  alternative = "ccc",
  descriptives = function(d, v) desc_table1(d, c(v$method, v$reference)),
  assumptions = function(d, v) { ch <- list(check_normality(num(d, v$method) - num(d, v$reference), "las diferencias método − referencia")); list(checks = ch, decision = decide_path(ch, TRUE)) },
  run = function(d, v, opts) {
    b <- .ba_stats(num(d, v$method), num(d, v$reference)); ccc <- tryCatch(epiR::epi.ccc(num(d, v$method), num(d, v$reference))$rho.c, error = function(e) NULL)
    st <- list(bias = b$bias, sd = b$sd, loa_l = b$loa_l, loa_u = b$loa_u, ci_lo = b$ci_lo, ci_hi = b$ci_hi, ref = b$ref, n = b$n, pct_l = 100 * b$loa_l / b$ref, pct_u = 100 * b$loa_u / b$ref, ccc = if (!is.null(ccc)) ccc$est else NA)
    list(main = b, stats = st,
         tables = list("Bland-Altman" = data.frame(Sesgo = b$bias, `IC 95% sesgo inf` = b$ci_lo, `IC 95% sesgo sup` = b$ci_hi, DE = b$sd, `LoA inferior` = b$loa_l, `LoA superior` = b$loa_u, `Nivel medio` = b$ref, n = b$n, `CCC de Lin` = st$ccc, check.names = FALSE)),
         notes = "LoA = sesgo ± 1.96·DE; requieren diferencias normales. Cada línea también se expresa como % del nivel medio de medición.")
  },
  posthoc = NULL,
  plots = list(ba = list(label = "Gráfico de Bland-Altman", fn = function(d, v, res, cos) plot_ba(num(d, v$method), num(d, v$reference), v$method, cos)),
               qq = list(label = "Q-Q de las diferencias", fn = function(d, v, res, cos) plot_qq(num(d, v$method) - num(d, v$reference), "Q-Q: diferencias", cos))),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- checks[[1]]
    paste0("La concordancia entre ", CTX(v$method), " y ", CTX(v$reference), " se evaluó mediante el análisis de Bland-Altman, ",
           if (isTRUE(sw$met)) "tras verificar que las diferencias no muestran evidencia de una distribución diferente de la normal" else "aunque las diferencias se apartan de la normalidad, por lo que los límites de acuerdo deben interpretarse con cautela y se recomienda el CCC de Lin",
           " (Shapiro-Wilk, p ", fmt_p(sw$p), "). El sesgo (diferencia media) fue de ", fmt_num(s$bias), " ", CTX("unidades"), " (IC 95 %: ", fmt_num(s$ci_lo), " a ", fmt_num(s$ci_hi), ")",
           if (s$ci_lo > 0 || s$ci_hi < 0) ", lo que indica una pequeña diferencia sistemática" else ", sin diferencia sistemática apreciable",
           ". Los límites de acuerdo al 95 % fueron de ", fmt_num(s$loa_l, 1), " a ", fmt_num(s$loa_u, 1), " (", fmt_num(s$pct_l, 1), " % a ", fmt_num(s$pct_u, 1), " % del nivel medio), el rango dentro del cual se espera que caiga el 95 % de las diferencias individuales; ",
           "si esa magnitud es clínicamente aceptable, ambos métodos pueden considerarse intercambiables [Bland y Altman, 1986; Giavarina, 2015].")
  }
))

ado_register(list(
  id = "ccc", family = "concordancia", order = 4, label = "Coeficiente de concordancia de Lin (CCC)",
  description = "Índice único de concordancia (precisión × exactitud) entre dos métodos continuos; robusto si las diferencias no son normales.",
  chapter = "CCC de Lin en R",
  roles = list(method = role("Método a evaluar (numérica)", "numeric"), reference = role("Referencia / gold standard (numérica)", "numeric")),
  example = list(file = "datBF.txt", vars = list(method = "Inter", reference = "Basal"), note = "Dos mediciones continuas del mismo sujeto (ejemplo ilustrativo)."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, c(v$method, v$reference)),
  assumptions = function(d, v) { ch <- list(check_normality(num(d, v$method) - num(d, v$reference), "las diferencias método − referencia")); list(checks = ch, decision = list(parametric = FALSE, all_met = FALSE, violated = character(0), reason = "El CCC no requiere normalidad de las diferencias.")) },
  run = function(d, v, opts) {
    m <- num(d, v$method); r <- num(d, v$reference); ok <- !is.na(m) & !is.na(r)
    cc <- epiR::epi.ccc(m[ok], r[ok], conf.level = 0.95, rep.measure = FALSE); pr <- cor(m[ok], r[ok])
    st <- list(ccc = cc$rho.c$est, lo = cc$rho.c$lower, hi = cc$rho.c$upper, Cb = cc$C.b, r = pr, bias = cc$sblalt$est, loa_l = cc$sblalt$lower, loa_u = cc$sblalt$upper, n = sum(ok))
    list(main = cc, stats = st,
         tables = list("CCC de Lin (epiR::epi.ccc)" = data.frame(`ρc` = st$ccc, `IC 95% inf` = st$lo, `IC 95% sup` = st$hi, `Cb (corrección de sesgo)` = st$Cb, `r de Pearson` = st$r, Sesgo = st$bias, `LoA inf` = st$loa_l, `LoA sup` = st$loa_u, Interpretación = ccc_band(st$ccc), check.names = FALSE)),
         notes = "ρc = r · Cb. Interpretación (McBride): < 0.90 pobre, 0.90–0.95 moderada, 0.95–0.99 sustancial, > 0.99 casi perfecta.")
  },
  posthoc = NULL,
  plots = list(identidad = list(label = "Gráfico de identidad (línea 1:1)", fn = function(d, v, res, cos) plot_identity(num(d, v$method), num(d, v$reference), v$method, v$reference, cos)),
               ba = list(label = "Bland-Altman (complementario)", fn = function(d, v, res, cos) plot_ba(num(d, v$method), num(d, v$reference), v$method, cos))),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- checks[[1]]
    paste0("La concordancia entre ", CTX(v$method), " y ", CTX(v$reference), " se evaluó mediante el Coeficiente de Concordancia de Lin",
           if (isFALSE(sw$met)) ", apropiado porque las diferencias muestran una distribución diferente de la normal (Shapiro-Wilk, p " else " (Shapiro-Wilk de las diferencias, p ", fmt_p(sw$p), if (isFALSE(sw$met)) ") y los límites de Bland-Altman pierden validez" else ")",
           ". El CCC fue ρc = ", fmt_num(s$ccc, 3), " (IC 95 %: ", fmt_num(s$lo, 3), " a ", fmt_num(s$hi, 3), "), una concordancia ", ccc_band(s$ccc), ". El coeficiente de correlación de Pearson fue r = ", fmt_num(s$r, 3),
           " y el factor de corrección de sesgo Cb = ", fmt_num(s$Cb, 3), if (s$Cb > 0.98) ", lo que indica una desviación prácticamente nula respecto de la línea de identidad (sin sesgo sistemático apreciable)" else ", lo que revela desviación respecto de la línea de identidad (sesgo sistemático)",
           " [Lin, 1989; Giavarina, 2015].")
  }
))
