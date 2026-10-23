# =============================================================================
# analyses/supervivencia.R — Kaplan-Meier + log-rank y regresión de Cox
# Reglas del libro: status entero 0/1; log-rank H0: S1(t) = ... = Sk(t);
# Cox: ≥ 10 eventos por covariable, VIF < 5, Schoenfeld p > 0.05 (PH);
# si PH se viola → estratificar.
# =============================================================================

.surv_df <- function(d, v) {
  dd <- data.frame(time = num(d, v$time), status = as_event01(d[[v$event]]))
  if (!is.null(v$group) && nzchar(v$group)) dd$group <- fac(d, v$group)
  if (!is.null(v$predictors)) for (p in v$predictors) { x <- d[[p]]; dd[[p]] <- if (is.numeric(x) && var_kind(x) %in% c("numeric", "numeric_small")) as.numeric(x) else factor(x) }
  dd[complete.cases(dd), ]
}
.surv_at <- function(fit, t) { s <- tryCatch(summary(fit, times = t, extend = TRUE), error = function(e) NULL)
  if (is.null(s)) return(NULL); data.frame(Estrato = if (is.null(s$strata)) "Global" else as.character(s$strata), Tiempo = s$time, `En riesgo` = s$n.risk, Eventos = s$n.event, `S(t)` = s$surv, `IC inf` = s$lower, `IC sup` = s$upper, check.names = FALSE) }

cdo_register(list(
  id = "km", family = "supervivencia", order = 1, label = "Kaplan-Meier y prueba log-rank",
  description = "Curvas de supervivencia (global y por grupo) y comparación con log-rank.",
  chapter = "Kaplan-Meier y prueba log-rank en R",
  roles = list(time = role("Tiempo hasta el evento", "time"), event = role("Estatus (1 = evento, 0 = censura)", "event"), group = role("Variable de grupo (opcional)", "factor", required = FALSE)),
  example = list(file = "CRIT1.txt", vars = list(time = "time", event = "status", group = "brux"), note = "Supervivencia de coronas de acero según bruxismo (402 pacientes, 11 eventos)."),
  alternative = NULL,
  descriptives = function(d, v) { dd <- d; dd$.estatus <- factor(as_event01(d[[v$event]]), labels = c("Sin evento", "Evento")); vars <- c(v$time, if (!is.null(v$group) && nzchar(v$group)) v$group); desc_table1(dd, vars, ".estatus") },
  assumptions = function(d, v) {
    dd <- .surv_df(d, v); ne <- sum(dd$status)
    checks <- list(check_obj("codificacion", "Codificación del estatus", TRUE, paste0("Estatus convertido a entero 0/1 (as.integer): ", ne, " eventos y ", sum(dd$status == 0), " censuras (", round(100 * mean(dd$status == 0)), " %). El KM no ajusta por covariables: las diferencias entre grupos pueden estar afectadas por confusión.")),
                   check_obj("eventos", "Número de eventos", ne >= 10, paste0(ne, " eventos. ", if (ne < 10) "Baja tasa de eventos: interpretar con cautela." else "Suficientes para estimaciones estables.")))
    list(checks = checks, decision = list(parametric = TRUE, all_met = TRUE, violated = character(0), reason = "Método no paramétrico (descriptivo)."))
  },
  run = function(d, v, opts) {
    dd <- .surv_df(d, v); has_g <- !is.null(dd$group)
    km0 <- survival::survfit(survival::Surv(time, status) ~ 1, data = dd)
    tt <- if (!is.null(opts$t_fixed) && !is.na(as.numeric(opts$t_fixed))) as.numeric(opts$t_fixed) else stats::median(dd$time)
    st <- list(n = nrow(dd), events = sum(dd$status), t_fixed = tt, med = unname(summary(km0)$table["median"]))
    tabs <- list("Supervivencia global en tiempo fijo" = .surv_at(km0, tt))
    fitg <- NULL
    if (has_g) {
      fitg <- survival::survfit(survival::Surv(time, status) ~ group, data = dd)
      lr <- survival::survdiff(survival::Surv(time, status) ~ group, data = dd); p_lr <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
      st$chisq <- lr$chisq; st$df_lr <- length(lr$n) - 1; st$p_lr <- p_lr; st$k <- length(lr$n); st$gname <- v$group
      tabs[["Supervivencia por grupo en tiempo fijo"]] <- .surv_at(fitg, tt)
      tabs[["Prueba log-rank"]] <- data.frame(`χ²` = lr$chisq, gl = st$df_lr, p = p_lr, check.names = FALSE)
      tabs[["Eventos observados vs. esperados por grupo"]] <- data.frame(Grupo = levels(dd$group), N = as.numeric(lr$n), Observados = as.numeric(lr$obs), Esperados = round(as.numeric(lr$exp), 2))
    }
    list(main = km0, fitg = fitg, dd = dd, stats = st, tables = tabs, notes = "La supervivencia se estima hasta el último tiempo con eventos; el tiempo fijo se elige en las opciones.")
  },
  posthoc = function(res, d, v, opts) {
    dd <- res$dd; if (is.null(dd$group) || nlevels(dd$group) < 3) return(NULL)
    prs <- combn(levels(dd$group), 2, simplify = FALSE)
    out <- do.call(rbind, lapply(prs, function(pr) { s <- dd[dd$group %in% pr, ]; s$group <- droplevels(s$group); lr <- survival::survdiff(survival::Surv(time, status) ~ group, data = s)
      data.frame(Comparación = paste(pr, collapse = " vs "), `χ²` = lr$chisq, p = pchisq(lr$chisq, 1, lower.tail = FALSE), check.names = FALSE) }))
    out$`p ajustada (Bonferroni)` <- pmin(1, out$p * nrow(out))
    list(tables = list("Log-rank por pares" = out), notes = "Comparaciones por pares con ajuste de Bonferroni.")
  },
  plots = list(km = list(label = "Curva(s) de Kaplan-Meier", fn = function(d, v, res, cos) if (!is.null(res$fitg)) plot_km(res$fitg, res$dd, cos, pval = TRUE, risk_table = FALSE, legend_labs = levels(res$dd$group)) else plot_km(res$main, res$dd, cos)),
               km_risk = list(label = "KM con tabla de riesgo", fn = function(d, v, res, cos) plot_km(if (!is.null(res$fitg)) res$fitg else res$main, res$dd, cos, pval = !is.null(res$fitg), risk_table = TRUE, legend_labs = if (!is.null(res$fitg)) levels(res$dd$group) else NULL))),
  opts_ui = list(t_fixed = list(label = "Tiempo fijo para S(t) (p. ej., 365.25 = 1 año)", type = "numeric", value = NA)),
  report = function(res, d, v, checks) {
    s <- res$stats; sf <- res$tables[[1]]
    paste0("La supervivencia de ", CTX("el evento estudiado"), " se estimó mediante el método de Kaplan-Meier en ", s$n, " sujetos, de los cuales ", s$events, " presentaron el evento. ",
           "La probabilidad de supervivencia global en t = ", fmt_num(s$t_fixed, 1), " fue de ", fmt_num(sf$`S(t)`[1], 3), " (IC 95 %: ", fmt_num(sf$`IC inf`[1], 3), "–", fmt_num(sf$`IC sup`[1], 3), ")",
           if (!is.na(s$med)) paste0("; la mediana de supervivencia fue ", fmt_num(s$med, 1)) else "; la mediana de supervivencia no se alcanzó",
           if (!is.null(s$p_lr)) paste0(". La prueba log-rank ", if (s$p_lr < ALPHA) "mostró diferencias significativas" else "no mostró diferencias significativas", " en la supervivencia según ", CTX(s$gname), " (χ² = ", fmt_num(s$chisq, 1), ", gl = ", s$df_lr, ", p ", fmt_p(s$p_lr), ")") else "",
           if (s$events < 10) ". Debido a la baja tasa de eventos, este resultado debe interpretarse con cautela" else "", " [Kaplan y Meier, 1958; Mantel, 1966].")
  }
))

cdo_register(list(
  id = "cox", family = "supervivencia", order = 2, label = "Regresión de Cox (riesgos proporcionales)",
  description = "Hazard ratios ajustados por covariables; verifica VIF, eventos por covariable y riesgos proporcionales (Schoenfeld).",
  chapter = "Regresión de Cox en R",
  roles = list(time = role("Tiempo hasta el evento", "time"), event = role("Estatus (1 = evento, 0 = censura)", "event"), predictors = role("Covariables", "any", multiple = TRUE, min = 1)),
  example = list(file = "CRIT1.txt", vars = list(time = "time", event = "status", predictors = c("asa", "age", "sex", "brux")), note = "Falla de coronas según ASA, edad, sexo y bruxismo."),
  alternative = NULL,
  descriptives = function(d, v) { dd <- d; dd$.estatus <- factor(as_event01(d[[v$event]]), labels = c("Sin evento", "Evento")); desc_table1(dd, c(v$time, v$predictors), ".estatus") },
  assumptions = function(d, v) {
    dd <- .surv_df(d, v); f <- as.formula(paste("survival::Surv(time, status) ~", paste(bq(v$predictors), collapse = " + ")))
    m <- survival::coxph(f, data = dd); k <- length(coef(m)); ne <- sum(dd$status)
    checks <- list(check_obj("epc", "Eventos por covariable (≥ 10)", ne / k >= 10, paste0(ne, " eventos para ", k, " coeficientes (", fmt_num(ne / k, 1), " por covariable). ", if (ne / k >= 10) "Adecuado." else "Modelo sobreparametrizado: estimaciones inestables e IC amplios (simplificar)."), epc = ne / k),
                   if (k > 1) check_vif(m) else check_obj("vif", "Multicolinealidad (VIF)", NA, "Una sola covariable: VIF no aplica."),
                   check_ph(m))
    list(checks = checks, decision = decide_path(checks[3], FALSE))
  },
  run = function(d, v, opts) {
    dd <- .surv_df(d, v); f <- as.formula(paste("survival::Surv(time, status) ~", paste(bq(v$predictors), collapse = " + ")))
    m <- survival::coxph(f, data = dd); s <- summary(m); co <- s$coefficients; ci <- s$conf.int
    tab <- data.frame(Término = rownames(co), β = co[, "coef"], HR = co[, "exp(coef)"], `HR IC inf` = ci[, "lower .95"], `HR IC sup` = ci[, "upper .95"], z = co[, "z"], p = co[, "Pr(>|z|)"], check.names = FALSE, row.names = NULL)
    st <- list(n = nrow(dd), events = sum(dd$status), lr_chi = s$logtest[["test"]], lr_df = s$logtest[["df"]], lr_p = s$logtest[["pvalue"]], conc = s$concordance[1], AIC = AIC(m),
               top = tab$Término[which.min(tab$p)], top_HR = tab$HR[which.min(tab$p)], top_lo = tab$`HR IC inf`[which.min(tab$p)], top_hi = tab$`HR IC sup`[which.min(tab$p)], top_p = min(tab$p))
    list(main = m, dd = dd, stats = st,
         tables = list("Hazard ratios (coxph)" = tab, "Ajuste global" = data.frame(`χ² razón de verosimilitud` = st$lr_chi, gl = st$lr_df, p = st$lr_p, Concordancia = st$conc, AIC = AIC(m), check.names = FALSE)),
         notes = "HR = 1 sin asociación; > 1 mayor riesgo; < 1 factor protector. Para muestras pequeñas se prefiere la prueba de razón de verosimilitud.")
  },
  posthoc = function(res, d, v, opts) {
    m <- res$main; a <- car::Anova(m); tab <- data.frame(Término = rownames(a), a, check.names = FALSE, row.names = NULL)
    st <- tryCatch(step(m, trace = 0), error = function(e) NULL); tabs <- list("Anova (razón de verosimilitud por término)" = tab)
    if (!is.null(st)) tabs[["Modelo simplificado por AIC (step)"]] <- data.frame(Fórmula = deparse(formula(st)), AIC = AIC(st), check.names = FALSE)
    list(tables = tabs, notes = "Simplificación hacia atrás: se prefiere el modelo más parsimonioso si anova() no es significativo.")
  },
  plots = list(schoenfeld = list(label = "Residuos de Schoenfeld (PH)", fn = function(d, v, res, cos) { z <- survival::cox.zph(res$main); function() { k <- nrow(z$table) - 1; op <- par(mfrow = c(1, max(1, k))); on.exit(par(op)); for (i in seq_len(k)) plot(z[i], main = rownames(z$table)[i]) } }),
               hr = list(label = "Hazard ratios con IC 95 %", fn = function(d, v, res, cos) { cos <- cos_merge(cos); tb <- res$tables[[1]]
                 p <- ggplot(tb, aes(x = HR, y = Término)) + geom_vline(xintercept = 1, linetype = 2, colour = cos$colour[2]) + geom_pointrange(aes(xmin = `HR IC inf`, xmax = `HR IC sup`), colour = cos$colour[1]) + scale_x_log10(); .finish(p, cos, "Hazard ratio (escala log)", "") }),
               km = list(label = "KM por la covariable más significativa (si es factor)", fn = function(d, v, res, cos) { dd <- res$dd; tb <- res$tables[[1]]; top <- v$predictors[sapply(v$predictors, function(p) any(startsWith(tb$Término[which.min(tb$p)], p)))][1]
                 if (is.na(top) || !is.factor(dd[[top]])) return(plot_km(survival::survfit(survival::Surv(time, status) ~ 1, data = dd), dd, cos))
                 dd$group <- dd[[top]]; plot_km(survival::survfit(survival::Surv(time, status) ~ group, data = dd), dd, cos, pval = TRUE, legend_labs = levels(dd$group)) })),
  report = function(res, d, v, checks) {
    s <- res$stats; epc <- checks[[1]]; vif <- .find(checks, "vif"); ph <- .find(checks, "ph")
    paste0("Se ajustó un modelo de riesgos proporcionales de Cox para identificar los factores asociados con ", CTX("el evento"), " en ", s$n, " sujetos, de los cuales ", s$events, " presentaron el evento. ",
           if (!is.na(vif$met)) paste0("Tras evaluar la multicolinealidad (VIF ", if (isTRUE(vif$met)) "< 5" else paste0("máximo = ", fmt_num(vif$max)), " en las covariables), ") else "",
           "el modelo ", if (s$lr_p < ALPHA) "resultó significativo" else "no resultó significativo", " (razón de verosimilitud χ² = ", fmt_num(s$lr_chi), ", gl = ", s$lr_df, ", p ", fmt_p(s$lr_p), "; concordancia = ", fmt_num(s$conc, 3), "). ",
           "El factor con mayor evidencia de asociación fue ", CTX(s$top), " (HR = ", fmt_num(s$top_HR), "; IC 95 %: ", fmt_num(s$top_lo), "–", fmt_num(s$top_hi), "; p ", fmt_p(s$top_p), ")",
           if (s$top_HR > 1) ", lo que indica un mayor riesgo" else ", lo que indica un menor riesgo (factor protector)", ". El supuesto de riesgos proporcionales se verificó con la prueba de Schoenfeld (p global ", fmt_p(ph$p), "), que ",
           if (isTRUE(ph$met)) "se cumplió" else "no se cumplió (considerar estratificación)", if (!isTRUE(epc$met)) ". Dado el reducido número de eventos, las estimaciones deben interpretarse con cautela" else "", " [Cox, 1972; Harrell, 2015].")
  }
))
