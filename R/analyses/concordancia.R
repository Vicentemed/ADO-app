# =============================================================================
# analyses/concordancia.R - Kappa (nominal/ordinal), ICC, Bland-Altman (alt: CCC)
# Regla del libro (caps. 43–47): categórica → kappa (ordinal → ponderado;
# > 2 evaluadores → Fleiss y AC1/AC2 de Gwet); continua con varios
# evaluadores → ICC; dos métodos continuos → Bland-Altman si las diferencias
# son normales, CCC de Lin si no (se recomienda usar ambos).
# =============================================================================

landis <- function(k) { if (is.na(k)) return("…"); if (k < 0) "sin acuerdo" else if (k <= 0.20) "ligero" else if (k <= 0.40) "razonable" else if (k <= 0.60) "moderado" else if (k <= 0.80) "sustancial" else "casi perfecto" }
koo_li <- function(i) { if (is.na(i)) return("…"); if (i < 0.5) "pobre" else if (i < 0.75) "moderada" else if (i < 0.9) "buena" else "excelente" }
koo_li_rango <- function(i) { if (is.na(i)) return(""); if (i < 0.5) "< 0.50" else if (i < 0.75) "0.50–0.75" else if (i < 0.9) "0.75–0.90" else "> 0.90" }
ccc_band <- function(r) { if (is.na(r)) return("…"); if (r < 0.90) "pobre" else if (r < 0.95) "moderada" else if (r < 0.99) "sustancial" else "casi perfecta" }
.conc_tab_df <- function(tab) { m <- addmargins(tab); data.frame(Fila = rownames(m), as.data.frame.matrix(m), check.names = FALSE, row.names = NULL) }

# Orden de categorías: el indicado por el usuario (si cubre todas las observadas)
# o el alfabético, que es el que usa irr::kappa2() sobre texto (como en el libro).
.kap_levels <- function(d, raters, txt) {
  obs <- unique(unlist(lapply(raters, function(r) as.character(d[[r]])))); obs <- obs[!is.na(obs) & nzchar(obs)]
  lv <- sort(obs); note <- NULL
  if (!is.null(txt) && nzchar(trimws(txt))) {
    l2 <- trimws(strsplit(txt, ",")[[1]]); l2 <- unique(l2[nzchar(l2)])
    if (all(obs %in% l2)) lv <- l2
    else note <- paste0("El orden indicado no incluye todas las categorías observadas (faltan: ", paste(setdiff(obs, l2), collapse = ", "), "); se usó el orden alfabético.")
  }
  list(lv = lv, note = note)
}
.gwet_ci <- function(s) { x <- suppressWarnings(as.numeric(strsplit(gsub("[()]", "", s), ",")[[1]])); if (length(x) == 2) x else c(NA, NA) }

ado_register(list(
  id = "kappa", family = "concordancia", order = 1, label = "Kappa de Cohen (2 evaluadores; ponderado si ordinal)",
  description = "Acuerdo entre dos evaluadores en una escala categórica. Con > 2 evaluadores calcula además Fleiss y AC1/AC2 de Gwet.",
  chapter = "Kappa en R",
  roles = list(raters = role("Columnas de los evaluadores (≥ 2; Cohen usa las dos primeras)", "factor", multiple = TRUE, min = 2)),
  example = list(file = "PGMax.txt", vars = list(raters = c("Obs2", "Obs3")), note = "Clasificación ordinal (AI < AII < BII < BIII) de 20 radiografías por 2 observadores."),
  alternative = NULL,
  descriptives = function(d, v) NULL,
  assumptions = function(d, v) {
    lv <- unique(unlist(lapply(v$raters, function(r) as.character(d[[r]])))); lv <- lv[!is.na(lv)]
    list(checks = list(check_obj("escala", "Tipo de escala", NA, paste0("Categorías observadas: ", paste(sort(lv), collapse = ", "), ". Si la escala es ORDINAL, el libro recomienda el kappa ponderado (lineal o cuadrático) e indicar el orden de las categorías; si es nominal, el no ponderado. Paradoja de prevalencia: con prevalencias extremas el kappa puede ser artificialmente bajo aunque el acuerdo observado sea alto (en tablas 2×2 conviene reportar PA+ y PA−)."))),
         decision = list(parametric = TRUE, all_met = TRUE, violated = character(0), reason = ""))
  },
  run = function(d, v, opts) {
    ordinal <- !identical(opts$escala, "nominal")
    L <- .kap_levels(d, v$raters, opts$levels); lv <- L$lv
    # Cohen: casos completos SOLO en los dos evaluadores comparados
    P <- data.frame(a = as.character(d[[v$raters[1]]]), b = as.character(d[[v$raters[2]]]), stringsAsFactors = FALSE); P <- P[complete.cases(P), ]
    lvp <- lv[lv %in% c(P$a, P$b)]                                   # categorías observadas en el par (como kappa2)
    C <- data.frame(a = sprintf("%03d", match(P$a, lvp)), b = sprintf("%03d", match(P$b, lvp)))   # códigos que respetan el orden indicado
    k_unw <- irr::kappa2(C, "unweighted"); k_lin <- irr::kappa2(C, "equal"); k_sq <- irr::kappa2(C, "squared")
    ci <- tryCatch(psych::cohen.kappa(cbind(match(P$a, lvp), match(P$b, lvp)))$confid, error = function(e) NULL)
    tab <- table(factor(P$a, levels = lvp), factor(P$b, levels = lvp)); names(dimnames(tab)) <- v$raters[1:2]
    st <- list(k_unw = k_unw$value, z_unw = k_unw$statistic, p_unw = k_unw$p.value, k_lin = k_lin$value, z_lin = k_lin$statistic, p_lin = k_lin$p.value,
               k_sq = k_sq$value, z_sq = k_sq$statistic, p_sq = k_sq$p.value, n = nrow(P), r1 = v$raters[1], r2 = v$raters[2], lv = lv, lvp = lvp, ordinal = ordinal,
               ci_unw = if (!is.null(ci)) ci[1, c(1, 3)] else c(NA, NA), ci_w = if (!is.null(ci)) ci[2, c(1, 3)] else c(NA, NA), nr = length(v$raters),
               k_fleiss = NA, ac1 = NA, ac2 = NA)
    tabs <- list("Tabla de acuerdo (evaluador 1 × evaluador 2)" = .conc_tab_df(tab),
                 "Kappa de Cohen (irr::kappa2)" = data.frame(Versión = c("No ponderado", "Ponderado lineal (equal)", "Ponderado cuadrático (squared)"), κ = c(st$k_unw, st$k_lin, st$k_sq), z = c(st$z_unw, st$z_lin, st$z_sq), p = c(st$p_unw, st$p_lin, st$p_sq), `Landis y Koch` = sapply(c(st$k_unw, st$k_lin, st$k_sq), landis), check.names = FALSE))
    if (!is.null(ci)) tabs[["IC 95 % (psych::cohen.kappa; ponderado = cuadrático)"]] <- data.frame(Tipo = c("No ponderado", "Ponderado"), ci, check.names = FALSE, row.names = NULL)
    # Paradoja de prevalencia (cap. 43.4): PA+ y PA− en tablas 2×2
    if (length(lvp) == 2) {
      a <- tab[1, 1]; b <- tab[1, 2]; cc <- tab[2, 1]; dd <- tab[2, 2]
      st$pa_pos <- 2 * a / (2 * a + b + cc); st$pa_neg <- 2 * dd / (2 * dd + b + cc); st$pos <- lvp[1]
      tabs[["Acuerdo positivo y negativo (Heo, 2008)"]] <- data.frame(`Categoría positiva` = lvp[1], `PA+ = 2a/(2a+b+c)` = st$pa_pos, `PA− = 2d/(2d+b+c)` = st$pa_neg, check.names = FALSE)
    }
    # > 2 evaluadores: Fleiss (no ponderado) y AC1/AC2 de Gwet sobre los evaluadores seleccionados
    if (length(v$raters) > 2) {
      R <- as.data.frame(lapply(d[v$raters], as.character), stringsAsFactors = FALSE)
      Rc <- R[complete.cases(R), ]
      fl <- irr::kappam.fleiss(Rc); st$k_fleiss <- fl$value; st$z_fleiss <- fl$statistic; st$p_fleiss <- fl$p.value; st$n_fleiss <- nrow(Rc)
      g1 <- tryCatch(irrCAC::gwet.ac1.raw(R, weights = "unweighted", categ.labels = lv)$est, error = function(e) NULL)
      g2 <- if (ordinal) tryCatch(irrCAC::gwet.ac1.raw(R, weights = "quadratic", categ.labels = lv)$est, error = function(e) NULL) else NULL
      mult <- data.frame(Coeficiente = "Kappa de Fleiss (no ponderado)", Valor = fl$value, `IC 95% inf` = NA, `IC 95% sup` = NA, z = fl$statistic, p = fl$p.value, `Landis y Koch` = landis(fl$value), check.names = FALSE)
      for (g in list(g1, g2)) if (!is.null(g)) { ci_g <- .gwet_ci(g$conf.int)
        mult <- rbind(mult, data.frame(Coeficiente = if (g$coeff.name == "AC1") "AC1 de Gwet (no ponderado)" else "AC2 de Gwet (pesos cuadráticos)", Valor = g$coeff.val, `IC 95% inf` = ci_g[1], `IC 95% sup` = ci_g[2], z = NA, p = g$p.value, `Landis y Koch` = landis(g$coeff.val), check.names = FALSE)) }
      if (!is.null(g1)) st$ac1 <- g1$coeff.val
      if (!is.null(g2)) { st$ac2 <- g2$coeff.val; st$ac2_ci <- .gwet_ci(g2$conf.int); st$ac2_p <- g2$p.value }
      tabs[[paste0("Múltiples evaluadores (", length(v$raters), "): Fleiss y Gwet (irrCAC)")]] <- mult
    }
    list(main = k_lin, stats = st, tables = tabs,
         notes = paste0("Orden de categorías usado: ", paste(lv, collapse = " < "), ". Escala de Landis y Koch. Cohen usa los casos completos de los dos primeros evaluadores (n = ", nrow(P), "). ",
                        if (length(v$raters) > 2) "El kappa de Fleiss es no ponderado (nominal): en escalas ordinales subestima el acuerdo; el AC2 de Gwet admite pesos ordinales y es robusto a la prevalencia. " else "",
                        if (!is.null(L$note)) L$note else ""))
  },
  posthoc = NULL,
  plots = list(matriz = list(label = "Matriz de acuerdo (mapa de calor)", fn = function(d, v, res, cos) {
    cos <- cos_merge(cos); M <- d[v$raters[1:2]]; M <- M[complete.cases(M), ]
    lv <- if (!is.null(res$stats$lvp)) res$stats$lvp else sort(unique(c(as.character(M[, 1]), as.character(M[, 2]))))
    df <- as.data.frame(table(factor(as.character(M[, 1]), levels = lv), factor(as.character(M[, 2]), levels = lv))); names(df) <- c("r1", "r2", "n")
    p <- ggplot(df, aes(r2, r1, fill = n)) + geom_tile(colour = "white") + geom_text(aes(label = n), size = 4) + scale_fill_gradient(low = "white", high = cos$colour[1], name = "Frecuencia")
    cos$legend <- TRUE; .finish(p, cos, v$raters[2], v$raters[1]) })),
  opts_ui = list(escala = list(label = "Tipo de escala", type = "select", choices = c("Ordinal (kappa ponderado)" = "ordinal", "Nominal (kappa no ponderado)" = "nominal"), value = "ordinal"),
                 levels = list(label = "Orden de categorías (ordinal), separadas por coma; vacío = alfabético", type = "text", value = "")),
  report = function(res, d, v, checks) {
    s <- res$stats
    lk <- function(k) { l <- landis(k); if (l == "sin acuerdo") "ausencia de acuerdo" else paste("un acuerdo", l) }
    ci_txt <- function(ci) if (any(is.na(ci))) "" else paste0("; IC 95 %: ", fmt_num(ci[1]), " a ", fmt_num(ci[2]))
    base <- paste0("La concordancia entre los evaluadores ", CTX(s$r1), " y ", CTX(s$r2), " para ", CTX("la clasificación"), " se evaluó con el Kappa de Cohen (n = ", s$n, "). ")
    body <- if (s$ordinal && length(s$lvp) > 2) {
      paste0("Como la clasificación es ordinal (", paste(s$lvp, collapse = " < "), "), procede la versión ponderada: el kappa no ponderado fue κ = ", fmt_num(s$k_unw), " (z = ", fmt_num(s$z_unw), "; p ", fmt_p(s$p_unw), "; ", lk(s$k_unw),
             "), mientras que el ponderado lineal (κ = ", fmt_num(s$k_lin), "; z = ", fmt_num(s$z_lin), "; p ", fmt_p(s$p_lin), ") y el cuadrático (κ = ", fmt_num(s$k_sq), "; z = ", fmt_num(s$z_sq), "; p ", fmt_p(s$p_sq), ci_txt(s$ci_w), ") mostraron ",
             if (landis(s$k_lin) == landis(s$k_sq)) paste0("ambos ", lk(s$k_sq)) else paste0(lk(s$k_lin), " y ", lk(s$k_sq), ", respectivamente"), " (escala de Landis y Koch).")
    } else {
      paste0(if (length(s$lvp) == 2) "Con dos categorías los kappas ponderados coinciden con el no ponderado, por lo que se reporta este último: " else "Como la clasificación es nominal, se reporta el kappa no ponderado: ",
             "κ = ", fmt_num(s$k_unw), " (z = ", fmt_num(s$z_unw), "; p ", fmt_p(s$p_unw), ci_txt(s$ci_unw), "), ", lk(s$k_unw), " (escala de Landis y Koch).")
    }
    pa <- if (!is.null(s$pa_pos)) paste0(" Ante la posible paradoja de prevalencia se reportan además la proporción de acuerdo positivo (PA+ = ", fmt_num(s$pa_pos), ", categoría ", s$pos, ") y negativo (PA− = ", fmt_num(s$pa_neg), ").") else ""
    mult <- if (!is.na(s$k_fleiss)) paste0(" Con los ", s$nr, " evaluadores, el kappa de Fleiss (no ponderado) fue κ = ", fmt_num(s$k_fleiss), " (z = ", fmt_num(s$z_fleiss), "; p ", fmt_p(s$p_fleiss), ")",
                                            if (!is.na(s$ac1)) paste0("; el AC1 de Gwet, robusto a la paradoja de prevalencia, fue ", fmt_num(s$ac1)) else "",
                                            if (!is.na(s$ac2)) paste0(", y el AC2 de Gwet con pesos cuadráticos, ", fmt_num(s$ac2), " (IC 95 %: ", fmt_num(s$ac2_ci[1]), "–", fmt_num(s$ac2_ci[2]), "; p ", fmt_p(s$ac2_p), "), ", lk(s$ac2)) else "", ".") else ""
    txt <- sub("[.[:space:]]+$", "", paste0(base, body, pa, mult))
    paste0(txt, " (Heo, 2008; Koo & Li, 2016; Landis & Koch, 1977", if (!is.na(s$ac1)) "; Gwet, 2008" else "", ").")
  }
))

ado_register(list(
  id = "icc", family = "concordancia", order = 2, label = "Coeficiente de correlación intraclase (ICC)",
  description = "Acuerdo entre ≥ 2 evaluadores en una variable continua (irr::icc y psych::ICC, modelos de Shrout-Fleiss).",
  chapter = "ICC en R",
  roles = list(raters = role("Columnas de los evaluadores (≥ 2, numéricas)", "numeric", multiple = TRUE, min = 2)),
  example = list(file = "Concor.txt", vars = list(raters = c("A", "B", "C")), note = "Tres evaluadores, escala continua."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, v$raters),
  assumptions = function(d, v) {
    M <- d[v$raters]; M <- M[complete.cases(M), ]; vs <- sapply(M, var)
    list(checks = list(check_obj("varianzas", "Varianzas por evaluador", NA, paste0("Varianzas: ", paste(paste0(names(vs), " = ", fmt_num(vs, 4)), collapse = "; "), ". Defina de antemano el modelo: consistencia ICC(3,1) para calibrar a evaluadores fijos; acuerdo absoluto ICC(2,1) para generalizar a otros evaluadores."))),
         decision = list(parametric = TRUE, all_met = TRUE, violated = character(0), reason = ""))
  },
  run = function(d, v, opts) {
    M <- as.data.frame(lapply(d[v$raters], as.numeric)); M <- M[complete.cases(M), ]
    unit <- if (is.null(opts$unit)) "single" else opts$unit; u <- if (unit == "single") "1" else "k"
    ic <- irr::icc(M, model = "twoway", type = "consistency", unit = unit); ia <- irr::icc(M, model = "twoway", type = "agreement", unit = unit)
    ps <- tryCatch(suppressMessages(suppressWarnings(psych::ICC(M)))$results, error = function(e) NULL)
    st <- list(icc_c = ic$value, lo_c = ic$lbound, hi_c = ic$ubound, p_c = ic$p.value, icc_a = ia$value, lo_a = ia$lbound, hi_a = ia$ubound, p_a = ia$p.value,
               n = nrow(M), k = ncol(M), unit = unit, u = u)
    tabs <- list("ICC (irr::icc, twoway)" = data.frame(
      Modelo = c(paste0("ICC(3,", u, ") consistencia - dos factores, mixto"), paste0("ICC(2,", u, ") acuerdo absoluto - dos factores, aleatorio")),
      ICC = c(ic$value, ia$value), F = c(ic$Fvalue, ia$Fvalue), gl1 = c(ic$df1, ia$df1), gl2 = c(ic$df2, ia$df2), p = c(ic$p.value, ia$p.value),
      `IC 95% inf` = c(ic$lbound, ia$lbound), `IC 95% sup` = c(ic$ubound, ia$ubound), `Koo y Li` = sapply(c(ic$value, ia$value), koo_li), check.names = FALSE))
    if (!is.null(ps)) tabs[["ICC de Shrout y Fleiss (psych::ICC)"]] <- data.frame(`Shrout-Fleiss` = c("ICC(1,1)", "ICC(2,1)", "ICC(3,1)", "ICC(1,k)", "ICC(2,k)", "ICC(3,k)"),
      psych = ps$type, ICC = ps$ICC, F = ps$F, gl1 = ps$df1, gl2 = ps$df2, p = ps$p, `IC 95% inf` = ps[["lower bound"]], `IC 95% sup` = ps[["upper bound"]], check.names = FALSE)
    list(main = ic, stats = st, tables = tabs,
         notes = paste0("Unidad: ", if (unit == "single") "medición única (ICC(·,1))" else "promedio de k mediciones (ICC(·,k))", ". Interpretación de Koo y Li: < 0.50 pobre, 0.50–0.75 moderada, 0.75–0.90 buena, > 0.90 excelente."))
  },
  posthoc = NULL,
  plots = list(pares = list(label = "Dispersión por pares con línea de identidad", fn = function(d, v, res, cos) {
    cos <- cos_merge(cos); M <- as.data.frame(lapply(d[v$raters], as.numeric)); M <- M[complete.cases(M), ]; prs <- combn(names(M), 2, simplify = FALSE)
    function() { op <- par(mfrow = c(1, length(prs))); on.exit(par(op))
      for (pr in prs) { plot(M[[pr[1]]], M[[pr[2]]], pch = 16, col = cos$colour[1], xlab = pr[1], ylab = pr[2], main = paste(pr[1], "vs.", pr[2]))
        abline(0, 1, col = cos$colour[2], lty = 2, lwd = 2); abline(lm(M[[pr[2]]] ~ M[[pr[1]]]), col = cos$colour[1], lwd = 1.5) } } })),
  opts_ui = list(unit = list(label = "Unidad", type = "select", choices = c("Medición única (single)" = "single", "Promedio de k mediciones (average)" = "average"), value = "single")),
  report = function(res, d, v, checks) {
    s <- res$stats
    paste0("La concordancia entre los ", s$k, " evaluadores para ", CTX("la variable medida"), " se evaluó con el Coeficiente de Correlación Intraclase (ICC) (n = ", s$n, if (s$unit == "average") "; unidad: promedio de las k mediciones" else "", "). ",
           CTX("Como el objetivo es calibrar a los evaluadores específicos que participarán en el estudio (evaluadores fijos)"),
           ", se reporta el ICC de consistencia, ICC(3,", s$u, "), que fue de ", fmt_num(s$icc_c, 3), " (IC 95 %: ", fmt_num(s$lo_c), " a ", fmt_num(s$hi_c), "; p ", fmt_p(s$p_c), "), una confiabilidad ", koo_li(s$icc_c), " (", koo_li_rango(s$icc_c), ") según los criterios de Koo y Li. ",
           "El ICC de acuerdo absoluto, ICC(2,", s$u, "), fue de ", fmt_num(s$icc_a, 3), " (IC 95 %: ", fmt_num(s$lo_a), " a ", fmt_num(s$hi_a), "); a diferencia del de consistencia, penaliza las diferencias sistemáticas (sesgo) entre evaluadores, por lo que la distancia entre ambos refleja la magnitud de ese sesgo (Koo & Li, 2016).")
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
    m <- num(d, v$method); r <- num(d, v$reference); ok <- !is.na(m) & !is.na(r)
    b <- .ba_stats(m, r); ccc <- tryCatch(epiR::epi.ccc(m[ok], r[ok])$rho.c, error = function(e) NULL)
    st <- list(bias = b$bias, sd = b$sd, loa_l = b$loa_l, loa_u = b$loa_u, ci_lo = b$ci_lo, ci_hi = b$ci_hi, ref = b$ref, n = b$n,
               pct_b = 100 * b$bias / b$ref, pct_l = 100 * b$loa_l / b$ref, pct_u = 100 * b$loa_u / b$ref,
               ccc = if (!is.null(ccc)) ccc$est else NA, ccc_lo = if (!is.null(ccc)) ccc$lower else NA, ccc_hi = if (!is.null(ccc)) ccc$upper else NA)
    list(main = b, stats = st,
         tables = list("Bland-Altman" = data.frame(Sesgo = b$bias, `Sesgo (% nivel medio)` = st$pct_b, `IC 95% sesgo inf` = b$ci_lo, `IC 95% sesgo sup` = b$ci_hi, DE = b$sd,
                                                   `LoA inferior` = b$loa_l, `LoA inf (%)` = st$pct_l, `LoA superior` = b$loa_u, `LoA sup (%)` = st$pct_u, `Nivel medio` = b$ref, n = b$n, check.names = FALSE),
                       "CCC de Lin (epiR::epi.ccc, complementario)" = data.frame(`ρc` = st$ccc, `IC 95% inf` = st$ccc_lo, `IC 95% sup` = st$ccc_hi, check.names = FALSE)),
         notes = "LoA = sesgo ± 1.96·DE; requieren diferencias normales. IC del sesgo con t(n − 1). Cada valor también se expresa como % del nivel medio de medición.")
  },
  posthoc = NULL,
  plots = list(ba = list(label = "Gráfico de Bland-Altman", fn = function(d, v, res, cos) { m <- num(d, v$method); r <- num(d, v$reference); ok <- !is.na(m) & !is.na(r); plot_ba(m[ok], r[ok], v$method, cos) }),
               qq = list(label = "Q-Q de las diferencias", fn = function(d, v, res, cos) plot_qq(num(d, v$method) - num(d, v$reference), "Q-Q: diferencias", cos))),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- checks[[1]]; sist <- s$ci_lo > 0 || s$ci_hi < 0
    paste0("La concordancia entre ", CTX(v$method), " y ", CTX(v$reference), " se evaluó mediante el análisis de Bland-Altman, ",
           if (isTRUE(sw$met)) "tras verificar que las diferencias no muestran evidencia de una distribución diferente de la normal" else "aunque las diferencias se apartan de la normalidad, por lo que sus límites de acuerdo no son confiables y la concordancia se resume mejor con el CCC de Lin",
           " (Shapiro-Wilk, p ", fmt_p(sw$p), "). El sesgo (diferencia media) fue de ", fmt_num(s$bias), " ", CTX("unidades"), " (IC 95 %: ", fmt_num(s$ci_lo), " a ", fmt_num(s$ci_hi), ")",
           if (sist) paste0(", lo que indica una diferencia sistemática (", v$method, " mide ", if (s$bias > 0) "por encima" else "por debajo", " de ", v$reference, ")") else "; el intervalo incluye el 0, por lo que no se detecta una diferencia sistemática",
           ". Los límites de acuerdo al 95 % fueron de ", fmt_num(s$loa_l, 1), " a ", fmt_num(s$loa_u, 1), " ", CTX("unidades"), " (de ", fmt_num(s$pct_l, 1), " % a ", fmt_num(s$pct_u, 1), " % del nivel medio de medición), el rango dentro del cual se espera que caiga el 95 % de las diferencias individuales; ",
           "si esa magnitud es clínicamente aceptable, ambos métodos pueden considerarse intercambiables (Bland & Altman, 1986; Giavarina, 2015).")
  }
))

ado_register(list(
  id = "ccc", family = "concordancia", order = 4, label = "Coeficiente de concordancia de Lin (CCC)",
  description = "Índice único de concordancia (precisión × exactitud) entre dos métodos continuos; no exige normalidad de las diferencias.",
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
         notes = "ρc = r · Cb (Cb = 1: la recta de regresión coincide con la línea de identidad). Interpretación (Tabla 47.1 del libro): < 0.90 pobre, 0.90–0.95 moderada, 0.95–0.99 sustancial, > 0.99 casi perfecta.")
  },
  posthoc = NULL,
  plots = list(identidad = list(label = "Gráfico de identidad (línea 1:1)", fn = function(d, v, res, cos) plot_identity(num(d, v$method), num(d, v$reference), v$method, v$reference, cos)),
               ba = list(label = "Bland-Altman (complementario)", fn = function(d, v, res, cos) { m <- num(d, v$method); r <- num(d, v$reference); ok <- !is.na(m) & !is.na(r); plot_ba(m[ok], r[ok], v$method, cos) })),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- checks[[1]]
    paste0("La concordancia entre ", CTX(v$method), " y ", CTX(v$reference), " se evaluó mediante el Coeficiente de Concordancia de Lin",
           if (isFALSE(sw$met)) paste0(", apropiado en este caso porque las diferencias muestran una distribución diferente de la normal (Shapiro-Wilk, p ", fmt_p(sw$p), ") y los límites de Bland-Altman pierden validez")
           else paste0(" (Shapiro-Wilk de las diferencias, p ", fmt_p(sw$p), "; con diferencias normales se recomienda acompañarlo del análisis de Bland-Altman)"),
           ". El CCC fue ρc = ", fmt_num(s$ccc, 3), " (IC 95 %: ", fmt_num(s$lo, 3), " a ", fmt_num(s$hi, 3), "), una concordancia ", ccc_band(s$ccc), ". El coeficiente de correlación de Pearson fue r = ", fmt_num(s$r, 3),
           " y el factor de corrección de sesgo, Cb = ", fmt_num(s$Cb, 3), " (Cb = 1 indica ausencia de desviación respecto de la línea de identidad; cuanto menor, mayor el sesgo sistemático)",
           " (Giavarina, 2015; Lin, 1989).")
  }
))
