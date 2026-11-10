# =============================================================================
# R/analyses/medias_rm.R - Diseños de medidas repetidas (Caps. 21, 22, 35, 36)
#   7. ANOVA de medidas repetidas: ezANOVA (F + Mauchly + GG/HF) o multinivel
#      lme() comparado con el modelo nulo por razón de verosimilitud (L.Ratio)
#   8. ANOVA de Friedman (alternativa para respuestas ordinales o no normales, Cap. 35)
#   8.5 Curvas de crecimiento (drc, LL.2) comparadas por razón de verosimilitud
# =============================================================================

# ---- utilidades comunes ----------------------------------------------------------
.rm_long <- function(d, v) { long <- wide_to_long(d, v$measures, v$subject); long <- long[complete.cases(long), ]
  # sujetos completos
  full <- names(which(table(long$Px) == length(v$measures))); long <- long[long$Px %in% full, ]; long$Px <- droplevels(long$Px)
  long$Times <- factor(long$Times, levels = v$measures); long }

.rm_ez <- function(long) suppressWarnings(ez::ezANOVA(data = long, dv = Response, wid = Px, within = Times, detailed = TRUE, type = 3))

# Esfericidad según el libro (Cap. 21.3): Mauchly p > 0.05 → F sin corrección;
# p < 0.05 → corrección: Greenhouse-Geisser si ε < 0.75, Huynh-Feldt si ε ≥ 0.75.
.rm_sph <- function(ez) {
  m <- ez[["Mauchly's Test for Sphericity"]]; cr <- ez[["Sphericity Corrections"]]
  if (is.null(m)) return(list(ok = FALSE))
  W <- as.numeric(m$W[1]); p <- as.numeric(m$p[1]); GGe <- as.numeric(cr$GGe[1]); HFe <- min(1, as.numeric(cr$HFe[1]))
  corr <- if (p < ALPHA) (if (GGe < 0.75) "GG" else "HF") else "none"
  list(ok = TRUE, W = W, p = p, GGe = GGe, HFe = HFe, pGG = as.numeric(cr[["p[GG]"]][1]), pHF = as.numeric(cr[["p[HF]"]][1]), corr = corr,
       eps = switch(corr, GG = GGe, HF = HFe, 1))
}
.rm_sph_check <- function(s) {
  if (!isTRUE(s$ok)) return(check_obj("esfericidad", "Esfericidad (Mauchly)", NA, "La prueba de Mauchly no aplica (dos niveles) o no pudo calcularse."))
  met <- s$p > ALPHA
  txt <- paste0("Prueba de Mauchly: W = ", fmt_num(s$W, 3), ", p ", fmt_p(s$p), ". ",
                if (met) "La esfericidad no se viola: se usa el estadístico F sin corrección."
                else paste0("La esfericidad se viola (ε de Greenhouse-Geisser = ", fmt_num(s$GGe, 3), "; ε de Huynh-Feldt = ", fmt_num(s$HFe, 3),
                            "). Si se mantiene el ANOVA de medidas repetidas se aplica la corrección de ",
                            if (s$corr == "GG") "Greenhouse-Geisser (ε < 0.75)" else "Huynh-Feldt (ε ≥ 0.75)", " (Cap. 21)."))
  check_obj("esfericidad", "Esfericidad (Mauchly)", met, txt, W = s$W, p = s$p, GGe = s$GGe, HFe = s$HFe)
}

# Modelo multinivel (Cap. 22.7): lme(Response ~ Times, random = ~1 | Px/Times, ML) vs. modelo nulo
.rm_lme <- function(long) {
  m1 <- nlme::lme(Response ~ Times, random = ~1 | Px/Times, data = long, method = "ML")
  m0 <- nlme::lme(Response ~ 1,     random = ~1 | Px/Times, data = long, method = "ML")
  a <- as.data.frame(anova(m0, m1)); a <- a[, intersect(c("df", "AIC", "BIC", "logLik", "L.Ratio", "p-value"), names(a))]
  tab <- data.frame(Modelo = c("Nulo (Response ~ 1)", "Tiempo (Response ~ Times)"), a, check.names = FALSE, row.names = NULL)
  fx <- as.data.frame(summary(m1)$tTable); fx <- data.frame(Término = rownames(fx), fx, check.names = FALSE, row.names = NULL)
  list(model = m1, null = m0, table = tab, fixed = fx, L = as.numeric(a$L.Ratio[2]), df = as.numeric(a$df[2] - a$df[1]), p = as.numeric(a[["p-value"]][2]))
}

# Media ± DE (mediana) por tiempo, como en las plantillas de reporte del libro
.rm_desc_txt <- function(long, digits = 2) {
  s <- split(long$Response, long$Times)
  paste(sapply(names(s), function(k) paste0(k, ": media = ", fmt_num(mean(s[[k]]), digits), " ± ", fmt_num(sd(s[[k]]), digits),
                                            " (mediana ", fmt_num(median(s[[k]]), digits), ")")), collapse = "; ")
}

# =============================================================================
# 7. ANOVA de medidas repetidas  (alt: Friedman)
# =============================================================================
ado_register(list(
  id = "rm_anova", family = "medias", order = 7,
  label = "ANOVA de medidas repetidas (≥ 3 tiempos)",
  description = "Compara ≥ 3 mediciones del mismo sujeto: ezANOVA (F, Mauchly, corrección GG/HF) o modelo multinivel lme() con razón de verosimilitud. Alternativa: ANOVA de Friedman.",
  chapter = "ANOVA de medidas repetidas en R",
  roles = list(measures = role("Columnas de medición (≥ 3 tiempos/condiciones)", "numeric", multiple = TRUE, min = 3),
               subject  = role("Identificador del sujeto (opcional)", "id", required = FALSE)),
  example = list(file = "DXKTw.txt", prep = function(d) subset(d, Grupo == 1), vars = list(measures = c("X6", "X7", "X8"), subject = "Px"),
                 note = "Dolor (NRS) a las 6, 7 y 8 h en el grupo Dexketoprofeno."),
  alternative = "friedman",
  descriptives = function(d, v) desc_table1(.rm_long(d, v), "Response", "Times"),
  assumptions = function(d, v) {
    long <- .rm_long(d, v)
    ez <- tryCatch(.rm_ez(long), error = function(e) NULL)
    s <- if (is.null(ez)) list(ok = FALSE) else .rm_sph(ez)
    checks <- list(.rm_sph_check(s))
    # Caps. 21-22 y 35.3: si Mauchly es significativo, el F sin corrección no es apropiado, pero el ANOVA de
    # medidas repetidas se mantiene con la corrección GG/HF (o el modelo multinivel). Friedman queda como
    # alternativa para respuestas ordinales o marcadamente no normales, a elección del usuario.
    viol <- isTRUE(s$ok) && s$p < ALPHA
    dec <- list(parametric = TRUE, all_met = !viol, violated = if (viol) "Esfericidad (Mauchly)" else character(0),
                reason = if (!viol) "Se cumplen los supuestos: se aplica el ANOVA de medidas repetidas con el estadístico F sin corrección."
                         else paste0("La prueba de Mauchly indica violación de la esfericidad: el estadístico F sin corrección no es apropiado. ",
                                     "Se mantiene el ANOVA de medidas repetidas con la corrección de ", if (identical(s$corr, "GG")) "Greenhouse-Geisser" else "Huynh-Feldt",
                                     " (Cap. 21); también puede usarse el enfoque multinivel lme(), que no requiere esfericidad (Cap. 22.7). ",
                                     "Si la respuesta es ordinal o su distribución es marcadamente asimétrica, elija la ruta del ANOVA de Friedman (Cap. 35)."))
    list(checks = checks, decision = dec)
  },
  run = function(d, v, opts) {
    long <- .rm_long(d, v); enfoque <- if (identical(opts$enfoque, "lme")) "lme" else "ez"
    ez <- .rm_ez(long); A <- ez$ANOVA; row <- which(A$Effect == "Times"); s <- .rm_sph(ez)
    F <- A$F[row]; DFn <- A$DFn[row]; DFd <- A$DFd[row]; p <- A$p[row]; ges <- A$ges[row]
    corr <- if (isTRUE(s$ok)) s$corr else "none"; eps <- if (isTRUE(s$ok)) s$eps else 1
    p_use <- switch(corr, GG = s$pGG, HF = s$pHF, p)
    lm_ <- tryCatch(.rm_lme(long), error = function(e) NULL)
    tabs <- list("Descriptiva por tiempo" = desc_numeric(long$Response, long$Times),
                 "ANOVA (ezANOVA, tipo III)" = A)
    if (isTRUE(s$ok)) {
      tabs[["Prueba de Mauchly"]] <- ez[["Mauchly's Test for Sphericity"]]
      tabs[["Correcciones de esfericidad"]] <- ez[["Sphericity Corrections"]]
    }
    if (corr != "none")
      tabs[[paste0("ANOVA con corrección de ", if (corr == "GG") "Greenhouse-Geisser" else "Huynh-Feldt")]] <-
        data.frame(ε = eps, F = F, `gl1 corregidos` = DFn * eps, `gl2 corregidos` = DFd * eps, `p corregida` = p_use, check.names = FALSE)
    if (!is.null(lm_)) {
      tabs[["Multinivel: modelo nulo vs. modelo con tiempo (lme, razón de verosimilitud)"]] <- lm_$table
      tabs[["Multinivel: efectos fijos (lme)"]] <- lm_$fixed
    }
    list(main = if (enfoque == "lme" && !is.null(lm_)) lm_$model else ez,
         stats = list(enfoque = if (enfoque == "lme" && !is.null(lm_)) "lme" else "ez",
                      F = F, df1 = DFn, df2 = DFd, p_raw = p, ges = ges, corr = corr, eps = eps, p_corr = p_use,
                      p = if (enfoque == "lme" && !is.null(lm_)) lm_$p else p_use,
                      L = if (!is.null(lm_)) lm_$L else NA, L_df = if (!is.null(lm_)) lm_$df else NA, L_p = if (!is.null(lm_)) lm_$p else NA,
                      W_m = if (isTRUE(s$ok)) s$W else NA, p_m = if (isTRUE(s$ok)) s$p else NA, k = nlevels(long$Times),
                      desc_txt = .rm_desc_txt(long), posthoc = if (identical(opts$posthoc, "bh")) "Benjamini-Hochberg" else "Bonferroni"),
         tables = tabs,
         notes = paste0(if (corr == "none") "Mauchly no significativo: F sin corrección. "
                        else paste0("Mauchly significativo: se reporta F con gl y p corregidos (", if (corr == "GG") "Greenhouse-Geisser, ε < 0.75" else "Huynh-Feldt, ε ≥ 0.75", "). "),
                        "El modelo multinivel lme(Response ~ Times, random = ~1 | Px/Times) no requiere esfericidad; se compara con el modelo nulo mediante la razón de verosimilitud (L.Ratio)."))
  },
  posthoc = function(res, d, v, opts) {
    long <- .rm_long(d, v); meth <- if (identical(opts$posthoc, "bh")) "bh" else "bonferroni"
    pw <- pairwise.t.test(long$Response, long$Times, p.adjust.method = if (meth == "bh") "BH" else "bonferroni", paired = TRUE)
    df <- as.data.frame(as.table(pw$p.value)); df <- df[!is.na(df$Freq), ]; names(df) <- c("Tiempo A", "Tiempo B", "p ajustada")
    list(tables = setNames(list(df), paste0("pairwise.t.test apareado (", if (meth == "bh") "Benjamini-Hochberg" else "Bonferroni", ")")),
         notes = "Bonferroni es el recomendado cuando la esfericidad se viola (Cap. 22.6).")
  },
  plots = list(
    perfiles = list(label = "Perfiles por sujeto + boxplot (medias)", fn = function(d, v, res, cos) plot_paired(.rm_long(d, v), cos, "mean")),
    nulo     = list(label = "Modelo nulo vs. modelo del experimento (multinivel)", fn = function(d, v, res, cos) {
      cos <- cos_merge(cos); long <- .rm_long(d, v)
      mt <- data.frame(Times = factor(levels(long$Times), levels = levels(long$Times)), media = as.numeric(tapply(long$Response, long$Times, mean)))
      gm <- data.frame(y = mean(long$Response))
      p <- ggplot(long, aes(Times, Response)) +
        geom_boxplot(aes(fill = Times), alpha = 0.55, outlier.shape = NA, width = 0.45) +
        geom_jitter(width = 0.1, alpha = 0.25, size = 1.5, colour = "gray40") +
        geom_line(data = mt, aes(Times, media, group = 1, colour = "Modelo del experimento"), linewidth = 1.6) +
        geom_point(data = mt, aes(Times, media, colour = "Modelo del experimento"), size = 4, shape = 18) +
        geom_hline(data = gm, aes(yintercept = y, colour = "Modelo nulo"), linewidth = 1.4, linetype = "dashed") +
        scale_fill_manual(values = .fillN(cos, nlevels(long$Times)), guide = "none") +
        scale_colour_manual(name = "", values = c("Modelo del experimento" = "#E74C3C", "Modelo nulo" = "#2E86C1"))
      .finish(p, cos, "Tiempo", "Respuesta") + theme(legend.position = "bottom") })
  ),
  opts_ui = list(
    enfoque = list(label = "Enfoque", type = "select",
                   choices = c("ANOVA de medidas repetidas (ezANOVA: F, Mauchly, GG/HF)" = "ez", "Modelo multinivel lme() (razón de verosimilitud, sin esfericidad)" = "lme"), value = "ez"),
    posthoc = list(label = "Método post-hoc", type = "select", choices = c("Bonferroni" = "bonferroni", "Benjamini-Hochberg" = "bh"), value = "bonferroni")),
  report = function(res, d, v, checks) {
    s <- res$stats; ph <- if (is.null(s$posthoc)) "Bonferroni" else s$posthoc
    sig <- s$p < ALPHA
    if (identical(s$enfoque, "lme")) {
      return(paste0(CTX("La variable"), " fue evaluada en ", s$k, " tiempos ", CTX("(describir)"), " (", s$desc_txt, "). ",
             "Las diferencias entre tiempos se analizaron con un modelo multinivel (modelo lineal mixto; lme) con el sujeto como efecto aleatorio, que no requiere el supuesto de esfericidad, ",
             "comparado con un modelo nulo mediante la prueba de razón de verosimilitud. El modelo con el efecto del tiempo ",
             if (sig) "explicó significativamente mejor" else "no explicó significativamente mejor", " la respuesta que el modelo nulo (L.Ratio = ", fmt_num(s$L),
             ", gl = ", s$L_df, ", p ", fmt_p(s$L_p), "). Los análisis post-hoc con corrección de ", ph, " ", CTX("confirmaron/revelaron …"),
             " [Field et al., 2012; Vittinghoff et al., 2012]."))
    }
    gl <- if (s$corr == "none") paste0(fmt_num(s$df1, 0), ", ", fmt_num(s$df2, 0)) else paste0(fmt_num(s$df1 * s$eps, 2), ", ", fmt_num(s$df2 * s$eps, 2))
    paste0(CTX("La variable"), " fue evaluada en ", s$k, " tiempos ", CTX("(describir)"), ". El supuesto de esfericidad fue evaluado con la prueba de Mauchly (W = ", fmt_num(s$W_m, 3), ", p ", fmt_p(s$p_m), ") y ",
           if (s$corr == "none") "no fue violado" else paste0("fue violado, por lo que se aplicó la corrección de ", if (s$corr == "GG") "Greenhouse-Geisser" else "Huynh-Feldt", " (ε = ", fmt_num(s$eps, 3), ")"),
           ". El ANOVA de medidas repetidas ", if (sig) "mostró" else "no mostró", " diferencia significativa entre los tiempos de evaluación (F[", gl, "] = ", fmt_num(s$F), ", p ", fmt_p(s$p), ", η²g = ", fmt_num(s$ges, 3), "). ",
           "Los análisis post-hoc con corrección de ", ph, " ", CTX("confirmaron/revelaron …"), " [Field et al., 2012].")
  }
))

# =============================================================================
# 8. ANOVA de Friedman
# =============================================================================
.fr_M <- function(d, v) { M <- as.matrix(d[v$measures]); storage.mode(M) <- "double"; M[complete.cases(M), , drop = FALSE] }

ado_register(list(
  id = "friedman", family = "medias", order = 8,
  label = "ANOVA de Friedman (≥ 3 tiempos)",
  description = "Alternativa no paramétrica al ANOVA de medidas repetidas para respuestas ordinales o marcadamente no normales; no requiere esfericidad.",
  chapter = "ANOVA de Friedman en R",
  roles = list(measures = role("Columnas de medición (≥ 3 tiempos/condiciones)", "numeric", multiple = TRUE, min = 3),
               subject  = role("Identificador del sujeto (opcional)", "id", required = FALSE)),
  example = list(file = "DXKTw.txt", prep = function(d) subset(d, Grupo == 1), vars = list(measures = c("X0", "X1", "X2"), subject = "Px"),
                 note = "Dolor (NRS) a las 0, 1 y 2 h en el grupo Dexketoprofeno."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(.rm_long(d, v), "Response", "Times"),
  assumptions = function(d, v) {
    M <- .fr_M(d, v)
    prs <- combn(colnames(M), 2, simplify = FALSE)
    dv <- sapply(prs, function(ij) var(M[, ij[1]] - M[, ij[2]])); names(dv) <- sapply(prs, paste, collapse = "−")
    ez <- tryCatch(.rm_ez(.rm_long(d, v)), error = function(e) NULL); s <- if (is.null(ez)) list(ok = FALSE) else .rm_sph(ez)
    txt <- paste0("Varianzas de las diferencias: ", paste0(names(dv), " = ", fmt_num(dv, 2), collapse = "; "), ". ",
                  if (isTRUE(s$ok)) paste0("Prueba de Mauchly: W = ", fmt_num(s$W, 3), ", p ", fmt_p(s$p), if (s$p < ALPHA) " → la esfericidad se viola (justifica Friedman)." else " → no hay evidencia de violación de la esfericidad.")
                  else "La prueba de Mauchly no pudo calcularse.")
    checks <- list(check_obj("esfericidad", "Esfericidad (varianzas de las diferencias y Mauchly)", if (isTRUE(s$ok)) s$p > ALPHA else NA, txt,
                             min = min(dv), max = max(dv), W = if (isTRUE(s$ok)) s$W else NA, p = if (isTRUE(s$ok)) s$p else NA))
    list(checks = checks, decision = list(parametric = FALSE, all_met = FALSE, violated = character(0),
                                          reason = "Prueba no paramétrica: no requiere el supuesto de esfericidad (Cap. 35.2)."))
  },
  run = function(d, v, opts) {
    M <- .fr_M(d, v); Fr <- friedman.test(M); N <- nrow(M); k <- ncol(M); W <- unname(Fr$statistic) / (N * (k - 1))
    R <- colSums(t(apply(M, 1, rank)))
    long <- data.frame(Response = as.vector(M), Times = factor(rep(colnames(M), each = N), levels = colnames(M)))
    mc <- tryCatch(pgirmess::friedmanmc(M), error = function(e) NULL)
    prs <- combn(colnames(M), 2, simplify = FALSE)
    deltas <- paste(sapply(prs, function(ij) paste0("Δ ", ij[1], "−", ij[2], " = ", fmt_num(mean(M[, ij[1]]) - mean(M[, ij[2]]), 2))), collapse = "; ")
    list(main = Fr,
         stats = list(chi2 = unname(Fr$statistic), df = unname(Fr$parameter), p = Fr$p.value, W = W, N = N, k = k,
                      desc_txt = .rm_desc_txt(long), deltas = deltas,
                      mc_sig = if (!is.null(mc)) rownames(mc$dif.com)[mc$dif.com$stat.signif] else NA),
         tables = list("Resultado de friedman.test()" = data.frame(`χ²F` = unname(Fr$statistic), gl = unname(Fr$parameter), p = Fr$p.value, `W de Kendall` = W, check.names = FALSE),
                       "Suma de rangos por condición (rangos dentro de cada sujeto)" = data.frame(Condición = colnames(M), `Suma de rangos` = unname(R), check.names = FALSE),
                       "Descriptiva por tiempo" = desc_numeric(long$Response, long$Times)),
         notes = "W de Kendall = χ²F / [N(k − 1)]: 0–0.1 trivial, 0.1–0.3 pequeño, 0.3–0.5 mediano, > 0.5 grande. friedman.test() corrige por empates.")
  },
  posthoc = function(res, d, v, opts) {
    M <- .fr_M(d, v)
    mc <- tryCatch(pgirmess::friedmanmc(M), error = function(e) NULL)
    if (is.null(mc)) return(list(tables = list(), notes = "pgirmess::friedmanmc no disponible."))
    list(tables = list("friedmanmc (pgirmess)" = data.frame(Comparación = rownames(mc$dif.com), mc$dif.com, row.names = NULL, check.names = FALSE)),
         notes = "Post-hoc solo si Friedman global es significativo (Cap. 36.4 lo muestra a manera de ejemplo).")
  },
  plots = list(perfiles = list(label = "Perfiles por sujeto + boxplot (mediana)", fn = function(d, v, res, cos) plot_paired(.rm_long(d, v), cos, "median"))),
  report = function(res, d, v, checks) {
    s <- res$stats; ch <- checks[[1]]
    viol <- isTRUE(!is.na(ch$p) && ch$p < ALPHA)
    intro <- if (viol) paste0("Dado que las varianzas de las diferencias entre los tiempos fueron heterogéneas (rango: ", fmt_num(ch$min), " a ", fmt_num(ch$max),
                              ") y la prueba de Mauchly confirmó la violación del supuesto de esfericidad (W = ", fmt_num(ch$W, 2), ", p ", fmt_p(ch$p), "), ")
             else paste0("Las varianzas de las diferencias entre los tiempos oscilaron entre ", fmt_num(ch$min), " y ", fmt_num(ch$max),
                         if (!is.na(ch$p)) paste0(" (prueba de Mauchly: W = ", fmt_num(ch$W, 2), ", p ", fmt_p(ch$p), ")"), "; ",
                         CTX("justificar el uso de la prueba no paramétrica"), ", ")
    ph <- if (length(s$mc_sig) == 1 && is.na(s$mc_sig)) CTX("no reveló/reveló diferencias entre …")
          else if (length(s$mc_sig) == 0) "no reveló diferencias significativas entre ningún par de tiempos"
          else paste0("reveló diferencias significativas entre ", paste(s$mc_sig, collapse = ", "))
    paste0(intro, "la diferencia en ", CTX("la variable"), " entre las ", s$k, " evaluaciones se evaluó con el ANOVA de Friedman. ",
           CTX("La variable"), " por tiempo fue: ", s$desc_txt, ", con diferencias de medias entre pares (", s$deltas, "). La diferencia ",
           if (s$p < ALPHA) "fue estadísticamente significativa" else "no fue estadísticamente significativa",
           " (χ²F[", s$df, "] = ", fmt_num(s$chi2, 1), ", p ", fmt_p(s$p), ", W = ", fmt_num(s$W, 2), ", efecto ", efecto_W(s$W), "). ",
           "El análisis post-hoc (friedmanmc) ", ph, " [Field et al., 2012; Friedman, 1937].")
  }
))

# =============================================================================
# 8.5 Curvas de crecimiento (drc): grupos a lo largo del tiempo (Cap. 22.9)
# =============================================================================
.drc_df <- function(d, v) {
  dd <- data.frame(y = num(d, v$response), t1 = num(d, v$time_lower), t2 = num(d, v$time), g = fac(d, v$group))
  dd <- dd[complete.cases(dd), ]; dd$g <- droplevels(dd$g); dd
}
.drc_fit <- function(dd) {
  full <- drc::drm(y ~ t1 + t2, fct = drc::LL.2(), data = dd, type = "event", curveid = g)
  null <- drc::drm(y ~ t1 + t2, fct = drc::LL.2(), data = dd, type = "event")
  list(full = full, null = null)
}

ado_register(list(
  id = "crecimiento", family = "medias", order = 8.5,
  label = "Curvas de crecimiento por grupo (drc, razón de verosimilitud)",
  description = "Modela la trayectoria de la respuesta en el tiempo (log-logística LL.2, datos de intervalo) y compara curvas por grupo contra una curva común con la prueba de razón de verosimilitud.",
  chapter = "ANOVA de medidas repetidas en R (22.9: modelos de crecimiento)",
  roles = list(response   = role("Respuesta (acotada entre 0 y 1, p. ej. absorbancia)", "numeric"),
               time       = role("Tiempo: límite superior del intervalo", "numeric"),
               time_lower = role("Tiempo: límite inferior del intervalo", "numeric"),
               group      = role("Variable de grupo", "factor")),
  example = list(file = "viabSBF.txt", vars = list(response = "Absorbance", time = "Time", time_lower = "Timeb", group = "Group"),
                 note = "Viabilidad celular (absorbancia) durante 21 días con tres sueros (Preciado-Gutiérrez et al., 2023)."),
  alternative = NULL,
  descriptives = function(d, v) {
    dd <- .drc_df(d, v)
    a <- aggregate(y ~ g + t2, data = dd, FUN = function(x) c(n = length(x), Media = mean(x), DE = sd(x)))
    data.frame(Grupo = a$g, Tiempo = a$t2, n = a$y[, "n"], Media = round(a$y[, "Media"], 4), DE = round(a$y[, "DE"], 4))
  },
  assumptions = function(d, v) {
    dd <- .drc_df(d, v)
    checks <- list(check_obj("rango", "Respuesta acotada entre 0 y 1 (LL.2)", all(dd$y >= 0 & dd$y <= 1),
                             paste0("Rango de la respuesta: ", fmt_num(min(dd$y), 4), " a ", fmt_num(max(dd$y), 4),
                                    ". LL.2() fija los límites inferior y superior en 0 y 1.")))
    list(checks = checks, decision = list(parametric = TRUE, all_met = TRUE, violated = character(0),
                                          reason = "Modelo de crecimiento: se comparan curvas por grupo contra una curva común (razón de verosimilitud)."))
  },
  run = function(d, v, opts) {
    dd <- .drc_df(d, v); f <- .drc_fit(dd)
    pe <- summary(f$full)$coefficients
    invisible(utils::capture.output(lr <- suppressMessages(anova(f$full, f$null, test = "Chisq")))); lr <- as.data.frame(unclass(lr))
    LR <- as.numeric(lr[2, 4]); df <- as.numeric(lr[2, 3]); p <- as.numeric(lr[2, 5])
    list(main = f$full,
         stats = list(LR = LR, df = df, p = p, groups = levels(dd$g), tmax = max(dd$t2)),
         tables = list("Parámetros por grupo (b = pendiente; e = tiempo al 50 % del máximo)" = data.frame(Parámetro = rownames(pe), pe, check.names = FALSE, row.names = NULL),
                       "Razón de verosimilitud: curvas por grupo vs. curva común" = data.frame(Modelo = c("Curvas por grupo", "Curva común"), lr, check.names = FALSE, row.names = NULL)),
         notes = "p > 0.05: permitir curvas distintas por grupo no mejora el ajuste (los grupos comparten la trayectoria). p < 0.05: al menos un grupo sigue una curva diferente.")
  },
  posthoc = NULL,
  plots = list(curvas = list(label = "Curvas log-logísticas ajustadas por grupo", fn = function(d, v, res, cos) {
    cos <- cos_merge(cos); dd <- .drc_df(d, v); f <- .drc_fit(dd)
    function() plot(f$full, log = "", xlab = .lab(cos, "xlab", "Tiempo"), ylab = .lab(cos, "ylab", v$response),
                    col = rep_len(cos$colour, nlevels(dd$g)), main = if (nzchar(cos$title)) cos$title else "")
  })),
  report = function(res, d, v, checks) {
    s <- res$stats
    paste0(CTX("La variable"), " se modeló a lo largo de ", fmt_num(s$tmax, 0), " ", CTX("días"), " mediante una función de crecimiento log-logística de dos parámetros (LL.2) con el paquete drc (Ritz et al., 2015). ",
           "Se comparó un modelo con curvas específicas por grupo (", paste(s$groups, collapse = ", "), ") frente a un modelo con una curva común mediante una prueba de razón de verosimilitud. ",
           if (s$p < ALPHA) "Se encontró evidencia de que las curvas difirieron entre grupos" else "No se encontró evidencia de que las curvas difirieran entre grupos",
           " (LR = ", fmt_num(s$LR), ", gl = ", fmt_num(s$df, 0), ", p ", if (s$p >= 0.9995) paste0("= ", fmt_num(s$p, 4)) else fmt_p(s$p), ")", if (s$p < ALPHA) ". " else ", por lo que los grupos comparten una misma trayectoria de crecimiento. ",
           CTX("Interpretación clínica."))
  }
))
