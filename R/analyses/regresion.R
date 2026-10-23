# =============================================================================
# analyses/regresion.R — Regresión lineal múltiple y regresión logística
# Reglas del libro: VIF < 5 (5–10 moderada, > 10 grave); Cook > 1 revisar;
# D-W ≈ 2; ≥ 10–20 sujetos por gl (lineal) / EPV ≥ 10 (logística);
# logística: Hosmer-Lemeshow p > 0.05 = calibración adecuada; AUC > 0.70 aceptable.
# =============================================================================

.prep_pred <- function(d, resp, preds) {
  dd <- d[c(resp, preds)]
  for (p in preds) if (!is.numeric(dd[[p]]) || var_kind(dd[[p]]) %in% c("binary") && !is.numeric(dd[[p]])) dd[[p]] <- factor(dd[[p]])
  for (p in preds) if (is.character(dd[[p]])) dd[[p]] <- factor(dd[[p]])
  dd[complete.cases(dd), ]
}
.coef_table <- function(m, expo = FALSE) {
  s <- summary(m)$coefficients; ci <- suppressMessages(confint(m))
  df <- data.frame(Término = rownames(s), Estimación = s[, 1], EE = s[, 2], `Estadístico` = s[, 3], p = s[, 4], `IC 95% inf` = ci[, 1], `IC 95% sup` = ci[, 2], check.names = FALSE, row.names = NULL)
  if (expo) { df$OR <- exp(df$Estimación); df$`OR IC inf` <- exp(df$`IC 95% inf`); df$`OR IC sup` <- exp(df$`IC 95% sup`) }
  df
}

ado_register(list(
  id = "lineal", family = "regresion", order = 1, label = "Regresión lineal múltiple",
  description = "Predice una respuesta continua a partir de varios predictores (lm). Diagnósticos: VIF, residuos, Cook, Durbin-Watson.",
  chapter = "Regresión lineal múltiple en R",
  roles = list(response = role("Variable de respuesta (numérica)", "numeric"), predictors = role("Predictores (numéricos o factores)", "any", multiple = TRUE, min = 1)),
  example = list(file = "datBF.txt", vars = list(response = "Final", predictors = c("Basal", "Group", "Age", "Sex")), note = "Fluorescencia final en función de la basal, grupo, edad y sexo."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, c(v$response, v$predictors)),
  assumptions = function(d, v) {
    dd <- .prep_pred(d, v$response, v$predictors); f <- as.formula(paste(bq(v$response), "~", paste(bq(v$predictors), collapse = " + ")))
    m <- lm(f, data = dd); res <- residuals(m); ck <- cooks.distance(m)
    n <- nrow(dd); k <- length(coef(m)) - 1
    dw <- tryCatch(car::durbinWatsonTest(m), error = function(e) NULL)
    checks <- list(
      check_normality(res, "los residuos"),
      if (k > 1) check_vif(m) else check_obj("vif", "Multicolinealidad (VIF)", NA, "Un solo predictor: VIF no aplica."),
      check_obj("cook", "Casos influyentes (Cook > 1)", max(ck) <= 1, paste0("Distancia de Cook máxima = ", fmt_num(max(ck), 3), if (max(ck) > 1) " — hay casos con influencia excesiva que merecen revisión." else " — sin casos con influencia excesiva."), max_cook = max(ck)),
      check_obj("n", "Tamaño muestral (10–20 sujetos por gl)", n / max(k, 1) >= 10, paste0(n, " observaciones para ", k, " parámetros (", fmt_num(n / max(k, 1), 1), " por gl); la regla del libro pide 10 a 20."), ratio = n / max(k, 1)),
      if (!is.null(dw)) check_obj("dw", "Independencia (Durbin-Watson ≈ 2)", NA, paste0("D-W = ", fmt_num(dw$dw), ", p ", fmt_p(dw$p), ". Solo es informativa si los datos tienen orden secuencial; en datos transversales refleja el orden de las filas."), dw = dw$dw, p = dw$p) else NULL)
    checks <- Filter(Negate(is.null), checks)
    list(checks = checks, decision = decide_path(checks, FALSE))
  },
  run = function(d, v, opts) {
    dd <- .prep_pred(d, v$response, v$predictors); f <- as.formula(paste(bq(v$response), "~", paste(bq(v$predictors), collapse = " + ")))
    m <- lm(f, data = dd); s <- summary(m); fs <- s$fstatistic
    eta <- tryCatch({ a <- car::Anova(m, type = "II"); ss <- a[, "Sum Sq"]; data.frame(Predictor = rownames(a), `η² parcial` = ss / (ss + a["Residuals", "Sum Sq"]), check.names = FALSE)[rownames(a) != "Residuals", ] }, error = function(e) NULL)
    list(main = m,
         stats = list(F = unname(fs[1]), df1 = unname(fs[2]), df2 = unname(fs[3]), p = pf(fs[1], fs[2], fs[3], lower.tail = FALSE), R2 = s$r.squared, aR2 = s$adj.r.squared, n = nrow(dd), k = length(coef(m)) - 1, AIC = AIC(m), BIC = BIC(m), formula = deparse(f)),
         tables = c(list("Coeficientes (lm)" = .coef_table(m), "Ajuste global" = data.frame(F = unname(fs[1]), gl1 = unname(fs[2]), gl2 = unname(fs[3]), p = pf(fs[1], fs[2], fs[3], lower.tail = FALSE), `R²` = s$r.squared, `R² ajustado` = s$adj.r.squared, AIC = AIC(m), BIC = BIC(m), check.names = FALSE)), if (!is.null(eta)) list("Tamaño de efecto por predictor (η² parcial, SS tipo II)" = eta)),
         notes = "Para simplificar el modelo, el libro elimina un término a la vez (el de mayor p) con car::Anova tipo II, comparando con anova()/AIC.")
  },
  posthoc = function(res, d, v, opts) {
    m <- res$main; dd <- m$model
    a2 <- car::Anova(m, type = "II"); tab <- data.frame(Término = rownames(a2), a2, check.names = FALSE, row.names = NULL)
    st <- tryCatch(step(m, trace = 0), error = function(e) NULL)
    tabs <- list("Anova tipo II (contribución de cada término)" = tab)
    if (!is.null(st)) tabs[["Modelo simplificado por AIC (step)"]] <- data.frame(Fórmula = deparse(formula(st)), AIC = AIC(st), `R² ajustado` = summary(st)$adj.r.squared, check.names = FALSE)
    list(tables = tabs, notes = "La selección automática debe complementarse con criterio clínico y parsimonia.")
  },
  plots = list(diag = list(label = "Diagnósticos (residuos, Q-Q, escala, influencia)", fn = function(d, v, res, cos) plot_diag(res$main)),
               scatter = list(label = "Observado vs. ajustado", fn = function(d, v, res, cos) { cos <- cos_merge(cos); df <- data.frame(ajustado = fitted(res$main), observado = res$main$model[[1]])
                 p <- ggplot(df, aes(ajustado, observado)) + geom_point(colour = cos$colour[1], alpha = 0.7, size = cos$point_size) + geom_abline(intercept = 0, slope = 1, linetype = 2, colour = cos$colour[2]); .finish(p, cos, "Valor ajustado", "Valor observado") })),
  report = function(res, d, v, checks) {
    s <- res$stats; vif <- .find(checks, "vif"); sw <- .find(checks, "normalidad"); ck <- .find(checks, "cook")
    paste0("Se construyó un modelo de regresión lineal múltiple para predecir ", CTX(v$response), " a partir de ", CTX(paste(v$predictors, collapse = ", ")), " sobre ", s$n, " observaciones ",
           "(", s$k, " parámetros; la regla de 10 a 20 sujetos por grado de libertad ", if (s$n / max(s$k, 1) >= 10) "se cumple" else "no se cumple", "). ",
           if (!is.na(vif$met)) paste0(if (isTRUE(vif$met)) "Se descartó la multicolinealidad (todos los VIF < 5). " else paste0("Se detectó multicolinealidad (VIF máximo = ", fmt_num(vif$max), "). ")) else "",
           "El modelo ", if (s$p < ALPHA) "resultó significativo" else "no resultó significativo", " (F(", fmt_num(s$df1, 0), ", ", fmt_num(s$df2, 0), ") = ", fmt_num(s$F), ", p ", fmt_p(s$p), "), con R² = ", fmt_num(s$R2, 3), " y R² ajustado = ", fmt_num(s$aR2, 3),
           ": explica alrededor del ", round(100 * s$R2), " % de la variabilidad de la respuesta. Los residuos ", if (isTRUE(sw$met)) "no mostraron desviación de la normalidad" else "mostraron desviación de la normalidad", " (Shapiro-Wilk, p ", fmt_p(sw$p), ")",
           if (isTRUE(ck$met)) " y no hubo casos con influencia excesiva (Cook ≤ 1)" else " y se identificaron casos influyentes (Cook > 1)", ". ", CTX("Describir los predictores significativos con β (IC 95 %) y su interpretación clínica"), " [Vittinghoff et al., 2012; Field, 2012].")
  }
))

ado_register(list(
  id = "logistica", family = "regresion", order = 2, label = "Regresión logística",
  description = "Predice una respuesta binaria (glm binomial). OR con IC, razón de verosimilitud, Hosmer-Lemeshow, ROC/AUC.",
  chapter = "Regresión logística en R",
  roles = list(response = role("Variable de respuesta (binaria)", "binary"), predictors = role("Predictores (numéricos o factores)", "any", multiple = TRUE, min = 1)),
  example = list(file = "datAsthma.txt", vars = list(response = "Bruxism", predictors = c("Asthma", "Age", "Sex", "OrParaF")), note = "Bruxismo en función de asma, edad, sexo y parafunción oral."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, v$predictors, v$response),
  assumptions = function(d, v) {
    dd <- .prep_pred(d, v$response, v$predictors); y <- as_event01(dd[[v$response]]); dd$.y <- y
    f <- as.formula(paste(".y ~", paste(bq(v$predictors), collapse = " + "))); m <- glm(f, data = dd, family = binomial(link = "logit"))
    k <- length(coef(m)) - 1; ev <- min(sum(y == 1), sum(y == 0)); epv <- ev / max(k, 1)
    checks <- list(
      check_obj("binaria", "Respuesta binaria", TRUE, paste0("Respuesta codificada 0/1 (evento = ", if (is.numeric(dd[[v$response]])) "valor máximo" else levels(factor(dd[[v$response]]))[2], "); ", sum(y == 1), " eventos y ", sum(y == 0), " no eventos.")),
      check_obj("epv", "Eventos por variable (EPV ≥ 10)", epv >= 10, paste0("EPV = ", fmt_num(epv, 1), " (", ev, " eventos del grupo menor / ", k, " parámetros). ", if (epv >= 10) "Tamaño adecuado." else "Modelo posiblemente sobreparametrizado: estimaciones inestables e IC amplios."), epv = epv),
      if (k > 1) check_vif(m) else check_obj("vif", "Multicolinealidad (VIF)", NA, "Un solo predictor: VIF no aplica."))
    list(checks = checks, decision = decide_path(checks, FALSE))
  },
  run = function(d, v, opts) {
    dd <- .prep_pred(d, v$response, v$predictors); y <- as_event01(dd[[v$response]]); dd$.y <- y
    f <- as.formula(paste(".y ~", paste(bq(v$predictors), collapse = " + "))); m <- glm(f, data = dd, family = binomial(link = "logit"))
    chi <- m$null.deviance - m$deviance; chdf <- m$df.null - m$df.residual; pchi <- pchisq(chi, chdf, lower.tail = FALSE)
    ll0 <- as.numeric(logLik(update(m, . ~ 1))); ll1 <- as.numeric(logLik(m)); n <- nrow(dd)
    r2_cs <- 1 - exp(2 * (ll0 - ll1) / n); r2_n <- r2_cs / (1 - exp(2 * ll0 / n)); r2_hl <- chi / (-2 * ll0)
    prob <- fitted(m); roc <- pROC::roc(y, prob, levels = c(0, 1), direction = "<", quiet = TRUE); auc <- as.numeric(pROC::auc(roc))
    hl <- tryCatch(ResourceSelection::hoslem.test(y, prob, g = 10), error = function(e) NULL)
    cut <- tryCatch(pROC::coords(roc, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity")), error = function(e) NULL)
    st <- list(chi = chi, chdf = chdf, pchi = pchi, r2_cs = r2_cs, r2_n = r2_n, r2_hl = r2_hl, auc = auc, hl_chi = if (!is.null(hl)) unname(hl$statistic) else NA, hl_p = if (!is.null(hl)) hl$p.value else NA, n = n, k = length(coef(m)) - 1, AIC = AIC(m),
               thr = if (!is.null(cut)) cut$threshold[1] else NA, sens = if (!is.null(cut)) cut$sensitivity[1] else NA, spec = if (!is.null(cut)) cut$specificity[1] else NA)
    list(main = m, roc = roc, stats = st,
         tables = list("Coeficientes y odds ratios (glm binomial)" = .coef_table(m, expo = TRUE),
                       "Ajuste global" = data.frame(`χ² razón de verosimilitud` = chi, gl = chdf, p = pchi, `R² Hosmer-Lemeshow` = r2_hl, `R² Cox-Snell` = r2_cs, `R² Nagelkerke` = r2_n, AIC = AIC(m), check.names = FALSE),
                       "Discriminación y calibración" = data.frame(AUC = auc, `Umbral óptimo (Youden)` = st$thr, Sensibilidad = st$sens, Especificidad = st$spec, `Hosmer-Lemeshow χ²` = st$hl_chi, `Hosmer-Lemeshow p` = st$hl_p, check.names = FALSE)),
         notes = "OR = exp(β). AUC: 0.5 azar, > 0.70 aceptable, 1 perfecta. Hosmer-Lemeshow p > 0.05 = calibración adecuada.")
  },
  posthoc = function(res, d, v, opts) {
    m <- res$main; a2 <- car::Anova(m, type = "II"); tab <- data.frame(Término = rownames(a2), a2, check.names = FALSE, row.names = NULL)
    st <- tryCatch(step(m, trace = 0), error = function(e) NULL); tabs <- list("Anova tipo II (razón de verosimilitud por término)" = tab)
    if (!is.null(st)) tabs[["Modelo simplificado por AIC (step)"]] <- data.frame(Fórmula = deparse(formula(st)), AIC = AIC(st), check.names = FALSE)
    list(tables = tabs, notes = "Eliminación hacia atrás un término a la vez; complementar con criterio clínico.")
  },
  plots = list(roc = list(label = "Curva ROC", fn = function(d, v, res, cos) plot_roc(res$roc, cos)),
               or = list(label = "Odds ratios con IC 95 %", fn = function(d, v, res, cos) { cos <- cos_merge(cos); tb <- res$tables[[1]]; tb <- tb[tb$Término != "(Intercept)", ]
                 p <- ggplot(tb, aes(x = OR, y = Término)) + geom_vline(xintercept = 1, linetype = 2, colour = cos$colour[2]) + geom_pointrange(aes(xmin = `OR IC inf`, xmax = `OR IC sup`), colour = cos$colour[1]) + scale_x_log10(); .finish(p, cos, "Odds ratio (escala log)", "") })),
  report = function(res, d, v, checks) {
    s <- res$stats; vif <- .find(checks, "vif"); epv <- .find(checks, "epv")
    paste0("Se ajustó un modelo de regresión logística para predecir ", CTX(v$response), " a partir de ", CTX(paste(v$predictors, collapse = ", ")), " (n = ", s$n, "; EPV = ", fmt_num(epv$epv, 1), if (isTRUE(epv$met)) ", adecuado" else ", insuficiente: interpretar con cautela", "). ",
           if (!is.na(vif$met)) paste0(if (isTRUE(vif$met)) "No hubo multicolinealidad (VIF < 5). " else "Se detectó multicolinealidad. ") else "",
           "El modelo ", if (s$pchi < ALPHA) "resultó significativo" else "no resultó significativo", " en la prueba de razón de verosimilitud (χ² = ", fmt_num(s$chi), ", gl = ", s$chdf, ", p ", fmt_p(s$pchi), "). Los pseudo-R² oscilaron entre ", fmt_num(min(s$r2_hl, s$r2_cs, s$r2_n)), " y ", fmt_num(max(s$r2_hl, s$r2_cs, s$r2_n)), ". ",
           "El modelo alcanzó un AUC de ", fmt_num(s$auc, 3), " (", if (s$auc > 0.7) "discriminación aceptable" else "discriminación limitada", ")",
           if (!is.na(s$hl_p)) paste0(" y una calibración ", if (s$hl_p > ALPHA) "adecuada" else "deficiente", " según la prueba de Hosmer-Lemeshow (p ", fmt_p(s$hl_p), ")") else "", ". ",
           CTX("Describir los factores asociados significativamente con sus razones de momios (OR, IC 95 %) — recuerde que el OR es una razón de momios, no de probabilidades"), " [Hosmer et al., 2013; Harrell, 2015].")
  }
))
