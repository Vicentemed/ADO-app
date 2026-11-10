# =============================================================================
# analyses/supervivencia.R - Kaplan-Meier + log-rank (cap. 52–53) y regresión
# de Cox (cap. 54–55), siguiendo el flujo de trabajo del libro.
# Reglas del libro: status entero 0/1 (as.integer); S(t) solo hasta el último tiempo
# observado (sin extrapolar); log-rank H0: S1(t) = ... = Sk(t);
# Cox: ≈ 10 eventos por covariable, VIF < 5, simplificación con Anova() y anova() entre
# modelos anidados, cph()/validate() (rms), cox.zph() del modelo final p > 0.05
# (por variable y global); si PH se viola → estratificar.
# =============================================================================

.sv_or <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

.surv_df <- function(d, v) {
  dd <- data.frame(time = num(d, v$time), status = as.integer(as_event01(d[[v$event]])))
  if (!is.null(v$group) && length(v$group) && nzchar(v$group)) dd$group <- fac(d, v$group)
  facs <- .sv_or(v$factors, character(0))
  for (p in c(v$predictors, facs)) { x <- d[[p]]
    dd[[p]] <- if (!(p %in% facs) && is.numeric(x) && var_kind(x) %in% c("numeric", "numeric_small")) as.numeric(x) else factor(x) }
  dd <- dd[complete.cases(dd), , drop = FALSE]
  for (p in names(dd)) if (is.factor(dd[[p]])) dd[[p]] <- droplevels(dd[[p]])
  dd
}
.sv_tab <- function(s, strata = TRUE) {
  if (is.null(s) || !length(s$time)) return(NULL)
  data.frame(Estrato = if (is.null(s$strata)) "Global" else as.character(s$strata), Tiempo = s$time, `En riesgo` = s$n.risk, Eventos = s$n.event,
             `S(t)` = s$surv, EE = s$std.err, `IC inf` = s$lower, `IC sup` = s$upper, check.names = FALSE)
}
# S(t) en tiempos fijos, sin extrapolar: cada estrato solo hasta su último tiempo observado
.surv_at <- function(fit, t, dd) {
  if (!length(t)) return(NULL)
  grp <- if (!is.null(fit$strata)) rep(names(fit$strata), fit$strata) else rep("Global", length(fit$time))
  tmax <- tapply(fit$time, grp, max)
  out <- lapply(names(tmax), function(g) { tt <- t[t <= tmax[[g]]]; if (!length(tt)) return(NULL)
    f1 <- if (g == "Global") fit else fit[g]
    tab <- .sv_tab(summary(f1, times = tt, extend = FALSE)); if (!is.null(tab)) tab$Estrato <- g; tab })
  do.call(rbind, out)
}
.sv_t <- function(x) format(round(x, 2), trim = TRUE, drop0trailing = TRUE)
.sv_parse_times <- function(txt) { txt <- .sv_or(txt, ""); x <- suppressWarnings(as.numeric(trimws(unlist(strsplit(as.character(txt), "[;,[:space:]]+"))))); sort(unique(x[!is.na(x) & x >= 0])) }

# =============================================================================
# KAPLAN-MEIER Y LOG-RANK
# =============================================================================
ado_register(list(
  id = "km", family = "supervivencia", order = 1, label = "Kaplan-Meier y prueba log-rank",
  description = "Tabla y curvas de supervivencia (global y por grupo), S(t) en tiempos clínicamente relevantes y comparación con log-rank.",
  chapter = "Kaplan-Meier y prueba log-rank en R",
  roles = list(time = role("Tiempo hasta el evento", "time"), event = role("Estatus (1 = evento, 0 = censura)", "event"), group = role("Variable de grupo (opcional)", "factor", required = FALSE)),
  example = list(file = "CRIT1.txt", vars = list(time = "time", event = "status", group = "brux"), note = "Supervivencia de coronas de acero según bruxismo (402 pacientes, 11 eventos)."),
  alternative = NULL,
  descriptives = function(d, v) { dd <- d; dd$.estatus <- factor(as_event01(d[[v$event]]), levels = 0:1, labels = c("Sin evento", "Evento")); vars <- c(v$time, if (!is.null(v$group) && length(v$group) && nzchar(v$group)) v$group); desc_table1(dd, vars, ".estatus") },
  assumptions = function(d, v) {
    dd <- .surv_df(d, v); ne <- sum(dd$status)
    checks <- list(check_obj("codificacion", "Codificación del estatus", TRUE, paste0("Estatus convertido a entero 0/1 (as.integer): ", ne, " eventos y ", sum(dd$status == 0), " censuras (", round(100 * mean(dd$status == 0)), " %). El KM es descriptivo y no ajusta por covariables: las diferencias entre grupos pueden estar afectadas por confusión.")),
                   check_obj("eventos", "Número de eventos", ne >= 10, paste0(ne, " eventos. ", if (ne < 10) "Baja tasa de eventos: interpretar con cautela." else "Suficientes para estimaciones estables.")))
    list(checks = checks, decision = list(parametric = TRUE, all_met = TRUE, violated = character(0), reason = "Método no paramétrico (descriptivo)."))
  },
  run = function(d, v, opts) {
    dd <- .surv_df(d, v); has_g <- !is.null(dd$group)
    km0 <- survival::survfit(survival::Surv(time, status) ~ 1, data = dd)
    tev <- dd$time[dd$status == 1]
    tt <- .sv_parse_times(opts$tiempos); user_t <- length(tt) > 0
    if (!user_t && length(tev)) tt <- unique(c(min(tev), max(tev)))       # el libro: menor y mayor tiempo con evento
    tmax <- max(dd$time); out_t <- tt[tt > tmax]
    st <- list(n = nrow(dd), events = sum(dd$status), times = tt[tt <= tmax], user_t = user_t, out_t = out_t, tmax = tmax,
               last_ev = if (length(tev)) max(tev) else NA)
    tabs <- list("Tabla de supervivencia global (tiempos con evento)" = .sv_tab(summary(km0)),
                 "S(t) global en los tiempos seleccionados" = .surv_at(km0, st$times, dd))
    fitg <- NULL
    if (has_g) {
      fitg <- survival::survfit(survival::Surv(time, status) ~ group, data = dd)
      lr <- survival::survdiff(survival::Surv(time, status) ~ group, data = dd); p_lr <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
      st$chisq <- lr$chisq; st$df_lr <- length(lr$n) - 1; st$p_lr <- p_lr; st$k <- length(lr$n); st$gname <- v$group
      tabs[["Prueba log-rank"]] <- data.frame(`χ²` = lr$chisq, gl = st$df_lr, p = p_lr, check.names = FALSE)
      tabs[["Eventos observados vs. esperados por grupo"]] <- data.frame(Grupo = levels(dd$group), N = as.numeric(lr$n), Observados = as.numeric(lr$obs), Esperados = as.numeric(lr$exp),
                                                                         `(O-E)^2/E` = (lr$obs - lr$exp)^2 / lr$exp, `(O-E)^2/V` = (lr$obs - lr$exp)^2 / diag(lr$var), check.names = FALSE)
      tabs[["Tabla de supervivencia por grupo (tiempos con evento)"]] <- .sv_tab(summary(fitg))
      tabs[["S(t) por grupo en los tiempos seleccionados"]] <- .surv_at(fitg, st$times, dd)
    }
    if (has_g) for (nm in names(tabs)) if (!is.null(tabs[[nm]]$Estrato)) tabs[[nm]]$Estrato <- sub("^group=", paste0(v$group, "="), tabs[[nm]]$Estrato)
    if (has_g) st$fixed_g <- tabs[["S(t) por grupo en los tiempos seleccionados"]]
    st$fixed <- tabs[["S(t) global en los tiempos seleccionados"]]
    list(main = km0, fitg = fitg, dd = dd, stats = st, tables = Filter(Negate(is.null), tabs),
         notes = paste0(if (user_t) "S(t) en los tiempos indicados" else "Sin tiempos indicados: S(t) en el menor y el mayor tiempo con evento (como el libro); indique tiempos clínicamente relevantes (p. ej., 365.25 = 1 año)",
                        ". La probabilidad solo se estima hasta el último tiempo observado de cada grupo (sin extrapolar)",
                        if (length(out_t)) paste0("; se omitieron por superar el seguimiento máximo (", fmt_num(tmax, 1), "): ", paste(out_t, collapse = ", ")) else "", "."))
  },
  posthoc = NULL,
  plots = list(km = list(label = "Curva(s) de Kaplan-Meier", fn = function(d, v, res, cos) if (!is.null(res$fitg)) plot_km(res$fitg, res$dd, cos, pval = TRUE, risk_table = FALSE, legend_labs = levels(res$dd$group)) else plot_km(res$main, res$dd, cos)),
               km_risk = list(label = "KM con tabla de sujetos en riesgo", fn = function(d, v, res, cos) plot_km(if (!is.null(res$fitg)) res$fitg else res$main, res$dd, cos, pval = !is.null(res$fitg), risk_table = TRUE, legend_labs = if (!is.null(res$fitg)) levels(res$dd$group) else NULL)),
               loglog = list(label = "Gráfico log(-log(S(t))) por grupo (riesgos proporcionales)", fn = function(d, v, res, cos) { cos <- cos_merge(cos)
                 function() { if (is.null(res$fitg)) { plot.new(); text(0.5, 0.5, "Requiere una variable de grupo."); return(invisible()) }
                   k <- length(res$fitg$strata); cl <- rep_len(cos$colour, k)
                   plot(res$fitg, fun = "cloglog", col = cl, lwd = 2, xlab = "log(tiempo)", ylab = "log(-log(S(t)))", main = if (nzchar(cos$title)) cos$title else "Gráfico log(-log) - Verificación de riesgos proporcionales")
                   legend("topleft", legend = levels(res$dd$group), col = cl, lwd = 2, bty = "n") } })),
  opts_ui = list(tiempos = list(label = "Tiempos clínicamente relevantes para S(t), separados por comas (p. ej., 365.25, 730.5). Vacío = menor y mayor tiempo con evento", type = "text", value = "")),
  report = function(res, d, v, checks) {
    s <- res$stats; sf <- s$fixed
    st_txt <- if (!is.null(sf) && nrow(sf)) paste0("La probabilidad de supervivencia global ", paste0("a t = ", .sv_t(sf$Tiempo), " fue de ", fmt_num(sf$`S(t)`, 3), " (IC 95 %: ", fmt_num(sf$`IC inf`, 3), "–", fmt_num(sf$`IC sup`, 3), ")", collapse = "; "),
                                                   if (!s$user_t) paste0(" ", CTX("sustituir por tiempos clínicamente relevantes, p. ej., a un año")) else "", ". ") else ""
    gtxt <- if (!is.null(s$p_lr)) {
      fg <- s$fixed_g
      paste0("La prueba log-rank ", if (s$p_lr < ALPHA) "mostró diferencias significativas" else "no mostró diferencias significativas", " en la supervivencia según ", CTX(s$gname), " (χ² = ", fmt_num(s$chisq, 1), ", gl = ", s$df_lr, ", p ", fmt_p(s$p_lr), ")",
             if (!is.null(fg) && nrow(fg)) paste0("; a t = ", .sv_t(max(fg$Tiempo)), " la supervivencia fue de ", paste0(fmt_num(fg$`S(t)`[fg$Tiempo == max(fg$Tiempo)], 3), " (IC 95 %: ", fmt_num(fg$`IC inf`[fg$Tiempo == max(fg$Tiempo)], 3), "–", fmt_num(fg$`IC sup`[fg$Tiempo == max(fg$Tiempo)], 3), ") en ", fg$Estrato[fg$Tiempo == max(fg$Tiempo)], collapse = " y ")) else "", ". ")
    } else ""
    paste0("La supervivencia de ", CTX("el tratamiento / evento estudiado"), " se estimó mediante el método de Kaplan-Meier en ", s$n, " sujetos, de los cuales ", s$events, " presentaron el evento ", CTX("definir el evento"), ". ",
           st_txt, gtxt, if (s$events < 10) "Debido a la baja tasa de eventos, estos resultados deben interpretarse con cautela. " else "",
           "(Harrington & Fleming, 1982; Kaplan & Meier, 1958; Mantel, 1966).")
  }
))

# =============================================================================
# REGRESIÓN DE COX
# =============================================================================
.sv_formula <- function(terms, dd) as.formula(paste("survival::Surv(time, status) ~", if (length(terms)) paste(bq(terms), collapse = " + ") else "1"), env = list2env(list(dd = dd), parent = globalenv()))   # survfit()/Anova() reevalúan data = dd
# Simplificación del libro: se elimina la covariable no significativa de mayor p (Anova, razón de
# verosimilitud), un término a la vez, y se reajusta; se conserva al menos una covariable.
.sv_backward <- function(dd, terms, alpha = ALPHA) {
  steps <- list()
  repeat {
    m <- survival::coxph(.sv_formula(terms, dd), data = dd); a <- car::Anova(m)
    # con una sola covariable car::Anova() devuelve la tabla secuencial (columnas Chisq / Pr(>|Chi|))
    pcol <- grep("^Pr", colnames(a), value = TRUE)[1]; ccol <- grep("Chisq", colnames(a), value = TRUE)[1]
    p <- setNames(a[bq(terms), pcol], terms)
    if (length(terms) <= 1 || max(p, na.rm = TRUE) <= alpha) break
    i <- which.max(p)
    steps[[length(steps) + 1]] <- data.frame(Paso = length(steps) + 1, `Covariable eliminada` = terms[i], `LR χ²` = a[bq(terms[i]), ccol], p = unname(p[i]), check.names = FALSE)
    terms <- terms[-i]
  }
  list(terms = terms, model = m, anova = a, steps = if (length(steps)) do.call(rbind, steps) else NULL)
}
.sv_hr <- function(m) { s <- summary(m); co <- s$coefficients; ci <- s$conf.int
  data.frame(Término = rownames(co), β = co[, "coef"], HR = co[, "exp(coef)"], `HR IC inf` = ci[, "lower .95"], `HR IC sup` = ci[, "upper .95"], EE = co[, "se(coef)"], z = co[, "z"], p = co[, "Pr(>|z|)"], check.names = FALSE, row.names = NULL) }
.sv_global <- function(m, label) { s <- summary(m)
  data.frame(Modelo = label, `LR χ²` = s$logtest[["test"]], `gl` = s$logtest[["df"]], `p (LR)` = s$logtest[["pvalue"]], `Wald` = s$waldtest[["test"]], `p (Wald)` = s$waldtest[["pvalue"]],
             `Score (log-rank)` = s$sctest[["test"]], `p (Score)` = s$sctest[["pvalue"]], Concordancia = s$concordance[1], AIC = AIC(m), check.names = FALSE) }
# cox.zph() por covariable y global, sobre el modelo final (cap. 55.7)
.sv_check_ph <- function(m) {
  z <- tryCatch(survival::cox.zph(m), error = function(e) NULL)
  if (is.null(z)) return(check_obj("ph", "Riesgos proporcionales (Schoenfeld, modelo final)", NA, "No fue posible evaluar cox.zph()."))
  tb <- z$table; pg <- tb["GLOBAL", "p"]; pv <- tb[rownames(tb) != "GLOBAL", "p"]; names(pv) <- rownames(tb)[rownames(tb) != "GLOBAL"]
  met <- all(c(pv, pg) > ALPHA); bad <- names(pv)[pv <= ALPHA]
  txt <- paste0("cox.zph() del modelo final: ", paste0(names(pv), " p ", fmt_p(pv), collapse = "; "), "; global p ", fmt_p(pg), ". ",
                if (met) "Sin relación entre los residuos de Schoenfeld y el tiempo: el supuesto de riesgos proporcionales se cumple."
                else paste0("El supuesto se viola", if (length(bad)) paste0(" (", paste(bad, collapse = ", "), ")") else "", ". Alternativa del libro: estratificar por la covariable problemática (strata())."))
  check_obj("ph", "Riesgos proporcionales (Schoenfeld, modelo final)", met, txt, p = pg, zph = z,
            table = data.frame(Término = rownames(tb), `χ²` = tb[, "chisq"], gl = tb[, "df"], p = tb[, "p"], check.names = FALSE, row.names = NULL))
}
.sv_cph_stats <- function(f, dd, label) { c1 <- rms::cph(f, x = TRUE, y = TRUE, data = dd); s <- c1$stats
  list(fit = c1, table = data.frame(Modelo = label, Obs = s[["Obs"]], Eventos = s[["Events"]], `LR χ²` = s[["Model L.R."]], gl = s[["d.f."]], p = s[["P"]], `Score χ²` = s[["Score"]], `p (Score)` = s[["Score P"]], R2 = s[["R2"]], Dxy = s[["Dxy"]], check.names = FALSE)) }
.sv_newdata <- function(txt, dd, vars) {
  txt <- .sv_or(txt, ""); if (!nzchar(trimws(txt))) return(NULL)
  recs <- trimws(strsplit(txt, ";")[[1]]); recs <- recs[nzchar(recs)]
  rows <- lapply(recs, function(r) { kv <- lapply(strsplit(strsplit(r, ",")[[1]], "="), trimws); setNames(lapply(kv, `[`, 2), sapply(kv, `[`, 1)) })
  nd <- data.frame(row.names = seq_along(rows))
  for (vn in vars) {
    vals <- sapply(rows, function(r) .sv_or(r[[vn]], NA_character_))
    if (any(is.na(vals))) stop("Pacientes nuevos: falta el valor de «", vn, "».")
    if (is.factor(dd[[vn]])) { nd[[vn]] <- factor(vals, levels = levels(dd[[vn]])); if (anyNA(nd[[vn]])) stop("Pacientes nuevos: nivel no válido para «", vn, "» (niveles: ", paste(levels(dd[[vn]]), collapse = ", "), ").") }
    else nd[[vn]] <- as.numeric(vals)
  }
  nd
}
.sv_preds <- function(v) { p <- c(v$predictors, .sv_or(v$factors, character(0))); p[nzchar(p)] }

ado_register(list(
  id = "cox", family = "supervivencia", order = 2, label = "Regresión de Cox (riesgos proporcionales)",
  description = "Hazard ratios ajustados por covariables. Flujo del libro: VIF, eventos por covariable, Anova, simplificación, anova() entre modelos anidados, cph()/validate() (rms), cox.zph() del modelo final y predicción para pacientes nuevos.",
  chapter = "Regresión de Cox en R",
  roles = list(time = role("Tiempo hasta el evento", "time"), event = role("Estatus (1 = evento, 0 = censura)", "event"),
               predictors = role("Covariables (continuas o ya categóricas en texto)", "any", multiple = TRUE, min = 1),
               factors = role("Covariables categóricas codificadas con números (se convierten a factor)", "factor", multiple = TRUE, required = FALSE,
                              help = "El libro convierte las covariables categóricas a factor (factor()) antes de ajustar el modelo.")),
  example = list(file = "CRIT1.txt", vars = list(time = "time", event = "status", predictors = c("asa", "age", "sex", "brux")), note = "Falla de coronas según ASA, edad, sexo y bruxismo."),
  alternative = NULL,
  descriptives = function(d, v) { dd <- d; dd$.estatus <- factor(as_event01(d[[v$event]]), levels = 0:1, labels = c("Sin evento", "Evento")); desc_table1(dd, c(v$time, .sv_preds(v)), ".estatus") },
  assumptions = function(d, v) {
    dd <- .surv_df(d, v); pr <- .sv_preds(v)
    m <- survival::coxph(.sv_formula(pr, dd), data = dd); k <- length(coef(m)); ne <- sum(dd$status)
    bw <- .sv_backward(dd, pr)
    checks <- list(check_obj("epc", "Eventos por covariable (≈ 10)", ne / k >= 10, paste0(ne, " eventos para ", k, " coeficientes del modelo saturado (", fmt_num(ne / k, 1), " por covariable). ",
                                     if (ne / k >= 10) "Adecuado." else "Modelo saturado sobreparametrizado: estimaciones inestables e IC amplios; interpretar con cautela y priorizar el modelo más parsimonioso."), epc = ne / k),
                   if (k > 1) check_vif(m) else check_obj("vif", "Multicolinealidad (VIF)", NA, "Una sola covariable: VIF no aplica."),
                   .sv_check_ph(bw$model))
    list(checks = checks, decision = decide_path(checks[3], FALSE))
  },
  run = function(d, v, opts) {
    dd <- .surv_df(d, v); pr <- .sv_preds(v)
    modS <- survival::coxph(.sv_formula(pr, dd), data = dd)
    bw <- .sv_backward(dd, pr); modF <- bw$model
    cmp <- if (length(bw$terms) < length(pr)) anova(modF, modS) else NULL
    cS <- tryCatch(.sv_cph_stats(.sv_formula(pr, dd), dd, "Saturado"), error = function(e) NULL)
    cF <- tryCatch(.sv_cph_stats(.sv_formula(bw$terms, dd), dd, "Final"), error = function(e) NULL)
    val <- NULL
    if (isTRUE(opts$boot) && !is.null(cS)) val <- tryCatch({ set.seed(1); vb <- rms::validate(cS$fit, B = 200)
      data.frame(Índice = rownames(vb), unclass(vb)[, c("index.orig", "training", "test", "optimism", "index.corrected", "n")], check.names = FALSE, row.names = NULL) }, error = function(e) NULL)
    ph <- .sv_check_ph(modF)
    nd <- tryCatch(.sv_newdata(opts$nuevos, dd, pr), error = function(e) e)
    pnew <- NULL; sf_new <- NULL
    if (inherits(nd, "error")) pnew <- data.frame(Aviso = conditionMessage(nd))
    else if (!is.null(nd)) { pp <- predict(modS, newdata = nd, type = "risk", se.fit = TRUE)
      pnew <- data.frame(Paciente = seq_len(nrow(nd)), nd, Riesgo = pp$fit, EE = pp$se.fit, `IC inf` = pmax(0, pp$fit - 1.96 * pp$se.fit), `IC sup` = pp$fit + 1.96 * pp$se.fit, check.names = FALSE)
      sf_new <- survival::survfit(modS, newdata = nd) }
    hrF <- .sv_hr(modF); sF <- summary(modF)
    top <- which.min(hrF$p)
    st <- list(n = nrow(dd), events = sum(dd$status), final_terms = bw$terms, formula = paste("Surv(time, status) ~", paste(bw$terms, collapse = " + ")),
               lr_chi = sF$logtest[["test"]], lr_df = sF$logtest[["df"]], lr_p = sF$logtest[["pvalue"]], conc = sF$concordance[1],
               sig = hrF[hrF$p < ALPHA, , drop = FALSE], top = hrF[top, ], ph_p = ph$p, ph_met = ph$met,
               cmp_p = if (!is.null(cmp)) cmp[2, "Pr(>|Chi|)"] else NA, cmp_chi = if (!is.null(cmp)) cmp[2, "Chisq"] else NA, cmp_df = if (!is.null(cmp)) cmp[2, "Df"] else NA)
    tabs <- list(
      "Modelo saturado: Anova (razón de verosimilitud)" = { aS <- car::Anova(modS); data.frame(Término = rownames(aS), as.data.frame(aS), check.names = FALSE, row.names = NULL) },
      "Modelo saturado: coeficientes y hazard ratios" = .sv_hr(modS),
      "Modelo final: coeficientes y hazard ratios" = hrF,
      "Pruebas globales y concordancia" = rbind(.sv_global(modS, paste("Saturado:", paste(pr, collapse = " + "))), .sv_global(modF, paste("Final:", paste(bw$terms, collapse = " + ")))),
      "Índices de ajuste (rms::cph)" = if (!is.null(cS) && !is.null(cF)) rbind(cS$table, cF$table) else NULL,
      "Validación interna por bootstrap del modelo saturado (rms::validate, B = 200, set.seed(1))" = val,
      "Riesgos proporcionales: cox.zph() del modelo final" = ph$table,
      "Predicción para pacientes nuevos (risk score del modelo saturado, IC 95 % por aproximación normal)" = pnew)
    list(main = modF, full = modS, bw = bw, cmp = cmp, zph = ph$zph, dd = dd, sf_new = sf_new, stats = st, tables = Filter(Negate(is.null), tabs),
         notes = "HR = 1 sin asociación; > 1 mayor riesgo; < 1 factor protector. Para muestras pequeñas se prefiere la prueba de razón de verosimilitud. Como en el libro, la validación bootstrap y la predicción usan el modelo saturado.")
  },
  posthoc = function(res, d, v, opts) {
    tabs <- list("Simplificación: covariables eliminadas (una a la vez, mayor p)" = res$bw$steps,
                 "Modelo final: Anova (razón de verosimilitud)" = data.frame(Término = rownames(res$bw$anova), as.data.frame(res$bw$anova), check.names = FALSE, row.names = NULL),
                 "Modelo saturado vs. modelo final (anova de modelos anidados)" = if (!is.null(res$cmp)) data.frame(Modelo = c("Final", "Saturado"), loglik = res$cmp$loglik, `χ²` = res$cmp$Chisq, gl = res$cmp$Df, p = res$cmp$`Pr(>|Chi|)`, check.names = FALSE) else NULL)
    list(tables = Filter(Negate(is.null), tabs), notes = "Si la comparación no es significativa (p > 0.05), el modelo reducido, más parsimonioso, es preferible.")
  },
  plots = list(
    schoenfeld = list(label = "Residuos de Schoenfeld escalados (modelo final)", fn = function(d, v, res, cos) { cos <- cos_merge(cos)
      if (is.null(res$zph)) return(function() { plot.new(); text(0.5, 0.5, "No fue posible evaluar cox.zph().") })
      survminer::ggcoxzph(res$zph, ggtheme = ado_theme(cos$base_size)) }),
    loglog = list(label = "Gráfico log(-log(S(t))) por covariable categórica del modelo final", fn = function(d, v, res, cos) { cos <- cos_merge(cos); dd <- res$dd
      fv <- res$stats$final_terms[sapply(res$stats$final_terms, function(p) is.factor(dd[[p]]))]
      function() { if (!length(fv)) { plot.new(); text(0.5, 0.5, "El modelo final no contiene covariables categóricas."); return(invisible()) }
        op <- par(mfrow = c(1, length(fv))); on.exit(par(op))
        for (p in fv) { dd$.g <- dd[[p]]; f <- survival::survfit(survival::Surv(time, status) ~ .g, data = dd); cl <- rep_len(cos$colour, nlevels(dd$.g))
          plot(f, fun = "cloglog", col = cl, lwd = 2, xlab = "log(tiempo)", ylab = "log(-log(S(t)))", main = paste("log(-log):", p))
          legend("topleft", legend = levels(dd$.g), col = cl, lwd = 2, bty = "n") } } }),
    pred = list(label = "Supervivencia predicha para pacientes nuevos", fn = function(d, v, res, cos) { cos <- cos_merge(cos); sf <- res$sf_new
      function() { if (is.null(sf)) { plot.new(); text(0.5, 0.5, "Indique pacientes nuevos en las opciones."); return(invisible()) }
        k <- ncol(as.matrix(sf$surv)); cl <- rep_len(cos$colour, k)
        plot(sf, col = cl, lwd = 2, lty = seq_len(k), ylim = c(0, 1), xlab = .lab(cos, "xlab", "Tiempo"), ylab = .lab(cos, "ylab", "Probabilidad de supervivencia S(t)"))
        legend("bottomleft", bty = "n", legend = paste("Paciente", seq_len(k)), col = cl, lwd = 2, lty = seq_len(k)) } })),
  opts_ui = list(
    boot = list(label = "Validación interna por bootstrap (rms::validate, B = 200; como el libro)", type = "checkbox", value = TRUE),
    nuevos = list(label = "Pacientes nuevos para predecir (todas las covariables: var=valor separados por comas; pacientes separados por «;»)", type = "text", value = "")),
  report = function(res, d, v, checks) {
    s <- res$stats; epc <- .find(checks, "epc"); vif <- .find(checks, "vif"); sig <- s$sig
    fac_txt <- if (nrow(sig)) paste0(if (nrow(sig) == 1) "el único factor asociado de forma significativa fue " else "los factores asociados de forma significativa fueron ",
                                     paste0(CTX(sig$Término), " (HR = ", fmt_num(sig$HR), "; IC 95 %: ", fmt_num(sig$`HR IC inf`), "–", fmt_num(sig$`HR IC sup`), "; p ", fmt_p(sig$p), ")", collapse = ", "),
                                     ", lo que indica ", paste0(ifelse(sig$HR > 1, "un mayor riesgo", "un menor riesgo (factor protector)"), collapse = " / "), " ", CTX("describir el evento"))
               else paste0("ninguna covariable se asoció de forma significativa con el evento (la de mayor evidencia fue ", CTX(s$top$Término), ": HR = ", fmt_num(s$top$HR), "; IC 95 %: ", fmt_num(s$top$`HR IC inf`), "–", fmt_num(s$top$`HR IC sup`), "; p ", fmt_p(s$top$p), ")")
    paste0("Se ajustó un modelo de riesgos proporcionales de Cox para identificar los factores asociados con ", CTX("el evento"), " en ", s$n, " sujetos, de los cuales ", s$events, " presentaron el evento. ",
           "Tras evaluar la multicolinealidad", if (!is.na(vif$met)) paste0(" (", if (isTRUE(vif$met)) "VIF < 5 en todas las covariables" else paste0("VIF máximo = ", fmt_num(vif$max)), ")") else "", " y simplificar el modelo",
           if (!is.na(s$cmp_p)) paste0(" (comparación con el modelo saturado: χ² = ", fmt_num(s$cmp_chi), ", gl = ", s$cmp_df, ", p ", fmt_p(s$cmp_p), ")") else "", ", ", fac_txt, ". ",
           "El supuesto de riesgos proporcionales se verificó con la prueba de Schoenfeld (p global ", fmt_p(s$ph_p), "), que ", if (isTRUE(s$ph_met)) "se cumplió" else "no se cumplió (considerar estratificación)", ". ",
           if (!isTRUE(epc$met)) "Dado el reducido número de eventos, las estimaciones deben interpretarse con cautela " else "", "(Cox, 1972; Harrell, 2015).")
  }
))
