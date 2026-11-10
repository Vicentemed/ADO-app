# =============================================================================
# analyses/regresion.R - Regresión lineal múltiple (cap. 48–49) y regresión
# logística (cap. 50–51), siguiendo el flujo de trabajo del libro:
#  1) especificar los factores (as.factor) → modelo inicial con todos los predictores;
#  2) supuestos/diagnóstico del modelo inicial (VIF < 5; Cook > 1; DFFITS > 0.5; D-W ≈ 2);
#  3) valores únicos y GAM (s()) para evaluar tendencias curvilíneas → términos I() centrados;
#  4) simplificación hacia atrás: un término a la vez (mayor p, car::Anova tipo II),
#     respetando la jerarquía polinómica; comparación de modelos anidados con anova()/AIC/BIC;
#     método de todos los subconjuntos (leaps::regsubsets, nbest = 5, nvmax = 10);
#  5) modelo final + validación (partición 80/20, validación cruzada; logística: bootstrap rms)
#     + tamaño del efecto / predicción.
# =============================================================================

.rg_or <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# ---- Preparación: factores, polinomios centrados y términos ------------------------
# opts$poly: "age=2; lwt=4" (variable = grado, máx. 4). El libro centra la variable (x - mean(x))
# antes de añadir los términos I(x^k) (cap. 49.5).
.rg_parse_poly <- function(txt, ok_vars) {
  txt <- .rg_or(txt, ""); if (!nzchar(trimws(txt))) return(integer(0))
  it <- trimws(unlist(strsplit(txt, "[;,]"))); it <- it[nzchar(it)]
  out <- integer(0)
  for (x in it) { kv <- trimws(strsplit(x, "[=:]")[[1]]); if (length(kv) != 2) next
    g <- suppressWarnings(as.integer(kv[2])); if (kv[1] %in% ok_vars && !is.na(g) && g >= 2) out[kv[1]] <- min(g, 4L) }
  out
}

.rg_build <- function(d, v, opts, binary = FALSE) {
  resp <- v$response; preds <- .rg_or(v$predictors, character(0)); facs <- .rg_or(v$factors, character(0))
  preds <- preds[nzchar(preds)]; facs <- facs[nzchar(facs)]; allp <- c(preds, facs)
  if (!length(allp)) stop("Seleccione al menos un predictor.")
  dd <- d[c(resp, allp)]
  for (p in facs) dd[[p]] <- factor(dd[[p]])                                  # el libro: within(dat, x <- as.factor(x))
  for (p in preds) if (!is.numeric(dd[[p]])) dd[[p]] <- factor(dd[[p]])
  dd <- dd[complete.cases(dd), , drop = FALSE]
  for (p in allp) if (is.factor(dd[[p]])) dd[[p]] <- droplevels(dd[[p]])
  y <- resp
  if (binary) { dd$.y <- as_event01(dd[[resp]]); y <- ".y" }
  numv <- allp[sapply(allp, function(p) is.numeric(dd[[p]]))]
  poly <- .rg_parse_poly(opts$poly, numv)
  rows <- list()
  for (p in allp) {
    if (p %in% names(poly)) {
      pc <- make.names(paste0(p, "_c")); dd[[pc]] <- dd[[p]] - mean(dd[[p]])
      for (g in seq_len(poly[[p]])) rows[[length(rows) + 1]] <- data.frame(c = if (g == 1) bq(pc) else paste0("I(", bq(pc), "^", g, ")"),
                                                                            o = if (g == 1) bq(p) else paste0("I(", bq(p), "^", g, ")"), var = p, deg = g)
    } else rows[[length(rows) + 1]] <- data.frame(c = bq(p), o = bq(p), var = p, deg = 1)
  }
  tt <- do.call(rbind, rows); tt$lab <- attr(terms(.rg_formula(y, tt)), "term.labels")
  # Modelo inicial del libro (mod_full / mod_log): todos los predictores, sin términos polinómicos
  f_full <- as.formula(paste(bq(y), "~", paste(bq(allp), collapse = " + ")), env = globalenv())
  list(dd = dd, y = y, resp = resp, tt = tt, poly = poly, preds = preds, facs = facs, allp = allp, numv = numv, f_full = f_full)
}
.rg_formula <- function(y, tt, scale = "c") as.formula(paste(bq(y), "~", if (NROW(tt)) paste(tt[[scale]], collapse = " + ") else "1"), env = globalenv())
.rg_ftxt <- function(y, tt, scale = "o") paste(y, "~", if (NROW(tt)) paste(tt[[scale]], collapse = " + ") else "1")

# ---- Simplificación hacia atrás (cap. 49.6 / 51.7) -----------------------------------
# Un término a la vez: se retira el de mayor p (Anova tipo II, car); si es un término polinómico,
# se retiran también los de mayor grado. Un término de grado inferior no se retira mientras un
# término de mayor grado de la misma variable sea significativo (queda «protegido»).
.rg_backward <- function(fit, y, tt, alpha = ALPHA) {
  steps <- list(); prot <- character(0)
  repeat {
    m <- fit(.rg_formula(y, tt)); a <- car::Anova(m, type = "II")
    pcol <- grep("^Pr", colnames(a), value = TRUE)[1]; p <- setNames(a[tt$lab, pcol], tt$lab)
    prot <- character(0); removed <- FALSE
    if (nrow(tt) > 1) for (i in order(p, decreasing = TRUE)) {
      if (is.na(p[i]) || p[i] <= alpha) break
      hi <- which(tt$var == tt$var[i] & tt$deg > tt$deg[i])
      if (length(hi) && any(p[hi] <= alpha, na.rm = TRUE)) { prot <- c(prot, tt$lab[i]); next }
      drop <- c(i, hi)
      steps[[length(steps) + 1]] <- data.frame(Paso = length(steps) + 1, `Término(s) eliminado(s)` = paste(tt$lab[drop], collapse = ", "), `p (Anova tipo II)` = unname(p[i]), check.names = FALSE)
      tt <- tt[-drop, , drop = FALSE]; removed <- TRUE; break
    }
    if (!removed) break
  }
  list(tt = tt, model = m, anova = a, protected = prot, steps = if (length(steps)) do.call(rbind, steps) else NULL)
}
# Variables «protegidas» no significativas: se compara el modelo con y sin la variable (todos sus
# grados) mediante anova() y AIC/BIC; si p > 0.05 se prefiere el modelo más parsimonioso (cap. 49.6.1 y 49.7).
.rg_blocks <- function(fit, y, bw, test = NULL, alpha = ALPHA) {
  tt <- bw$tt; m <- bw$model; rows <- list(); vars <- unique(tt$var[tt$lab %in% bw$protected])
  for (vr in vars) {
    tr <- tt[tt$var != vr, , drop = FALSE]; mr <- fit(.rg_formula(y, tr))
    an <- if (is.null(test)) anova(mr, m) else anova(mr, m, test = test)
    pc <- grep("^Pr", colnames(an))[1]; p <- an[2, pc]
    rows[[length(rows) + 1]] <- data.frame(Variable = vr, `Modelo sin la variable` = .rg_ftxt(y, tr, "c"), gl = abs(an[2, which(colnames(an) == "Df")[1]]), p = p,
                                           `AIC sin` = AIC(mr), `AIC con` = AIC(m), `BIC sin` = BIC(mr), `BIC con` = BIC(m),
                                           Decisión = if (p > alpha) "Se elimina (modelo más parsimonioso)" else "Se conserva", check.names = FALSE)
    if (p > alpha) { tt <- tr; m <- mr }
  }
  list(tt = tt, model = m, table = if (length(rows)) do.call(rbind, rows) else NULL)
}

# ---- Todos los subconjuntos (leaps::regsubsets, nbest = 5, nvmax = 10) -------------------
.rg_leaps <- function(y, tt, dd, yvec = NULL) {
  f <- .rg_formula(y, tt); d2 <- dd; if (!is.null(yvec)) d2[[y]] <- yvec
  X <- model.matrix(lm(f, data = d2)); asg <- attr(X, "assign"); cn <- colnames(X)
  rs <- leaps::regsubsets(f, data = d2, nbest = 5, nvmax = 10); s <- summary(rs)
  best <- c(`R² ajustado (máximo)` = which.max(s$adjr2), `Cp de Mallows (mínimo)` = which.min(s$cp), `BIC (mínimo)` = which.min(s$bic))
  tl <- lapply(best, function(i) { cols <- colnames(s$which)[s$which[i, ]]; cols <- setdiff(cols, "(Intercept)"); unique(tt$lab[asg[match(cols, cn)]]) })
  tab <- data.frame(Criterio = names(best), `Modelo (índice)` = paste0(rownames(s$which)[best], " [", best, "]"),
                    Términos = sapply(tl, paste, collapse = " + "), check.names = FALSE, row.names = NULL)
  list(table = tab, terms = tl)
}

# ---- Tablas ---------------------------------------------------------------------------------
.coef_table <- function(m, expo = FALSE) {
  s <- summary(m)$coefficients; ci <- suppressMessages(confint(m))
  df <- data.frame(Término = rownames(s), Estimación = s[, 1], EE = s[, 2], `Estadístico` = s[, 3], p = s[, 4], `IC 95% inf` = ci[, 1], `IC 95% sup` = ci[, 2], check.names = FALSE, row.names = NULL)
  if (expo) { df$OR <- exp(df$Estimación); df$`OR IC inf` <- exp(df$`IC 95% inf`); df$`OR IC sup` <- exp(df$`IC 95% sup`) }
  df
}
.rg_anova_df <- function(a) data.frame(Término = rownames(a), as.data.frame(a), check.names = FALSE, row.names = NULL)
.rg_unique <- function(b) if (length(b$numv)) data.frame(Predictor = b$numv, `Valores únicos` = sapply(b$numv, function(p) length(unique(b$dd[[p]]))),
  `Suavizado s() en el GAM` = ifelse(sapply(b$numv, function(p) length(unique(b$dd[[p]]))) > 8, "Sí", "No (pocos valores únicos)"), check.names = FALSE, row.names = NULL) else NULL
# GAM temporal (paquete gam): s() solo en predictores continuos con suficientes valores únicos (cap. 49.4)
.rg_gam <- function(b, family = gaussian()) {
  sv <- b$numv[sapply(b$numv, function(p) length(unique(b$dd[[p]])) > 8)]
  if (!length(sv)) return(NULL)
  rhs <- sapply(b$allp, function(p) if (p %in% sv) paste0("s(", bq(p), ")") else bq(p))
  e <- list2env(list(s = gam::s, .gdat = b$dd, .gfam = family), parent = globalenv())   # plot.Gam reevalúa la llamada en este entorno
  e$.gf <- as.formula(paste(bq(b$y), "~", paste(rhs, collapse = " + ")), env = e)
  g <- tryCatch(eval(quote(gam::gam(.gf, data = .gdat, family = .gfam)), e), error = function(e) NULL); if (is.null(g)) return(NULL)
  an <- summary(g)$anova; an <- as.data.frame(an[grepl("^s\\(", rownames(an)), , drop = FALSE])
  list(model = g, sterms = paste0("s(", bq(sv), ")"), table = data.frame(Término = rownames(an), an, check.names = FALSE, row.names = NULL))
}
# Sujetos nuevos: "x=1, g=A; x=2, g=B" (punto decimal)
.rg_newdata <- function(txt, dd, vars) {
  txt <- .rg_or(txt, ""); if (!nzchar(trimws(txt))) return(NULL)
  recs <- trimws(strsplit(txt, ";")[[1]]); recs <- recs[nzchar(recs)]
  rows <- lapply(recs, function(r) { kv <- lapply(strsplit(strsplit(r, ",")[[1]], "="), trimws); setNames(lapply(kv, `[`, 2), sapply(kv, `[`, 1)) })
  nd <- data.frame(row.names = seq_along(rows))
  for (vn in vars) {
    vals <- sapply(rows, function(r) .rg_or(r[[vn]], NA_character_))
    if (any(is.na(vals))) stop("Sujetos nuevos: falta el valor de «", vn, "».")
    if (is.factor(dd[[vn]])) { nd[[vn]] <- factor(vals, levels = levels(dd[[vn]])); if (anyNA(nd[[vn]])) stop("Sujetos nuevos: nivel no válido para «", vn, "» (niveles: ", paste(levels(dd[[vn]]), collapse = ", "), ").") }
    else nd[[vn]] <- as.numeric(vals)
  }
  nd
}
.rg_poly_plot <- function(b, cos) {
  cos <- cos_merge(cos); sv <- b$numv[sapply(b$numv, function(p) length(unique(b$dd[[p]])) > 8)]
  if (!length(sv)) return(function() { plot.new(); text(0.5, 0.5, "No hay predictores continuos con suficientes valores únicos.") })
  yv <- b$dd[[b$y]]
  function() { op <- par(mfrow = c(1, length(sv))); on.exit(par(op))
    for (p in sv) { x <- b$dd[[p]]; xs <- seq(min(x), max(x), length.out = 200)
      plot(x, yv, pch = 16, col = "#88888860", xlab = p, ylab = b$resp, main = p)
      r2 <- numeric(0)
      for (g in 1:4) { m <- lm(yv ~ poly(x, g, raw = TRUE)); r2[g] <- summary(m)$r.squared
        lines(xs, predict(m, newdata = data.frame(x = xs)), col = rep_len(cos$colour, 4)[g], lwd = 2, lty = g) }
      legend("topright", bty = "n", cex = 0.8, col = rep_len(cos$colour, 4), lty = 1:4, lwd = 2,
             legend = paste0(c("Lineal", "Cuadrática", "Cúbica", "Cuártica"), " (R² = ", round(r2, 3), ")")) } }
}
.rg_vif_txt <- function(vif) if (is.na(vif$met)) "" else if (isTRUE(vif$met)) "se descartó la multicolinealidad (todos los VIF < 5)" else paste0("se detectó multicolinealidad (VIF máximo = ", fmt_num(vif$max), ")")
.rg_poly_txt <- function(poly) if (length(poly)) paste0("se evaluaron términos polinómicos centrados (", paste0(names(poly), " hasta grado ", poly, collapse = "; "), ")") else "se evaluó la linealidad de los predictores continuos mediante un GAM"

.rg_opts_common <- list(
  poly = list(label = "Términos polinómicos (variable=grado; p. ej., age=2; lwt=4). Vacío = solo términos lineales; decídalo con el GAM y los valores únicos", type = "text", value = "")
)

# =============================================================================
# REGRESIÓN LINEAL MÚLTIPLE
# =============================================================================
ado_register(list(
  id = "lineal", family = "regresion", order = 1, label = "Regresión lineal múltiple",
  description = "Predice una respuesta continua (lm). Flujo del libro: modelo inicial, diagnóstico (VIF, Cook, DFFITS, D-W), GAM, simplificación hacia atrás con Anova tipo II, regsubsets, validación (80/20 y validación cruzada) y η².",
  chapter = "Regresión lineal múltiple en R",
  roles = list(response = role("Variable de respuesta (numérica)", "numeric"),
               predictors = role("Predictores continuos (o ya categóricos en texto)", "any", multiple = TRUE, min = 1),
               factors = role("Predictores categóricos codificados con números (se convierten a factor)", "factor", multiple = TRUE, required = FALSE,
                              help = "El libro especifica los factores al inicio (as.factor), p. ej., race 1/2/3. Las variables numéricas que no se indiquen aquí se tratan como continuas.")),
  example = list(file = "datBF.txt", vars = list(response = "Final", predictors = c("Basal", "Group", "Age", "Sex")), note = "Fluorescencia final en función de la basal, grupo, edad y sexo."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, c(v$response, v$predictors, v$factors)),
  assumptions = function(d, v) {
    b <- .rg_build(d, v, list()); m <- lm(b$f_full, data = b$dd)
    ck <- cooks.distance(m); dff <- dffits(m); n <- nrow(b$dd); k <- length(coef(m)) - 1
    set.seed(123); dw <- tryCatch(car::durbinWatsonTest(m), error = function(e) NULL)
    nd <- sum(dff > 0.5, na.rm = TRUE)
    checks <- list(
      check_normality(residuals(m), "los residuos del modelo inicial"),
      if (k > 1) check_vif(m) else check_obj("vif", "Multicolinealidad (VIF)", NA, "Un solo predictor: VIF no aplica."),
      check_obj("cook", "Casos influyentes (Cook > 1)", max(ck) <= 1, paste0("Distancia de Cook máxima = ", fmt_num(max(ck), 3), " (observaciones con mayor valor: ", paste(names(sort(ck, decreasing = TRUE))[1:min(3, n)], collapse = ", "), ")",
        if (max(ck) > 1) ": hay casos con influencia excesiva que merecen revisión." else ": ningún valor supera 1."), max_cook = max(ck)),
      check_obj("dffits", "DFFITS > 0.5", NA, paste0(nd, " observaciones con DFFITS > 0.5", if (nd) paste0(" (", paste(head(names(dff)[dff > 0.5], 10), collapse = ", "), if (nd > 10) ", …" else "", ")") else "",
        ". Si el DFFIT es cercano a 0, el caso probablemente no influye de forma importante en los coeficientes.")),
      check_obj("n", "Tamaño muestral (10–20 sujetos por gl)", n / max(k, 1) >= 10, paste0(n, " observaciones para ", k, " gl de los predictores (", fmt_num(n / max(k, 1), 1), " por gl); la regla del libro pide 10 a 20 sujetos por grado de libertad."), ratio = n / max(k, 1)),
      if (!is.null(dw)) check_obj("dw", "Independencia (Durbin-Watson ≈ 2)", NA, paste0("D-W = ", fmt_num(dw$dw), ", p ", fmt_p(dw$p), " (set.seed(123)). Un valor cercano a 2 indica ausencia de autocorrelación; la prueba solo es informativa si los datos tienen orden secuencial (en datos transversales refleja el orden de las filas)."), dw = dw$dw, p = dw$p) else NULL)
    checks <- Filter(Negate(is.null), checks)
    list(checks = checks, decision = decide_path(checks, FALSE))
  },
  run = function(d, v, opts) {
    b <- .rg_build(d, v, opts); dd <- b$dd; fit <- function(f) lm(f, data = dd)
    full <- lm(b$f_full, data = dd); sf <- summary(full); fsf <- sf$fstatistic
    g <- .rg_gam(b)
    # Simplificación (modelo saturado centrado → modS → modelo final)
    sat <- fit(.rg_formula(b$y, b$tt)); bw <- .rg_backward(fit, b$y, b$tt); bl <- .rg_blocks(fit, b$y, bw)
    modS <- lm(.rg_formula(b$y, bw$tt, "o"), data = dd); fin <- lm(.rg_formula(b$y, bl$tt, "o"), data = dd)
    modB <- lm(as.formula(paste(bq(b$y), "~ 1")), data = dd)
    lp <- tryCatch(.rg_leaps(b$y, b$tt, dd), error = function(e) NULL)
    cmp <- list("Eliminación hacia atrás" = bw$tt, "Modelo final" = bl$tt)
    if (!is.null(lp)) cmp <- c(cmp, setNames(lapply(lp$terms, function(l) b$tt[b$tt$lab %in% l, , drop = FALSE]), paste("regsubsets:", names(lp$terms))))
    cmp_tab <- do.call(rbind, lapply(names(cmp), function(nm) { m <- fit(.rg_formula(b$y, cmp[[nm]], "o")); data.frame(Modelo = nm, Fórmula = .rg_ftxt(b$y, cmp[[nm]], "o"), `R² ajustado` = summary(m)$adj.r.squared, AIC = AIC(m), check.names = FALSE) }))
    s <- summary(fin); fs <- s$fstatistic; pF <- pf(fs[1], fs[2], fs[3], lower.tail = FALSE)
    set.seed(123); dw <- tryCatch(car::durbinWatsonTest(fin), error = function(e) NULL)
    sw <- check_normality(residuals(fin), "los residuos del modelo final")
    vif_f <- if (length(coef(fin)) > 2) check_vif(fin) else NULL
    eta <- tryCatch({ e <- lsr::etaSquared(fin); data.frame(Término = rownames(e), `η²` = e[, "eta.sq"], `η² parcial` = e[, "eta.sq.part"], check.names = FALSE, row.names = NULL) }, error = function(e) NULL)
    eta1 <- tryCatch(sum(lsr::etaSquared(fin, type = 1)[, "eta.sq"]), error = function(e) NA)
    # Validación: partición 80/20 (set.seed(123)), CV 10-fold y CV repetida 10 × 20
    yv <- dd[[b$y]]; same <- identical(bw$tt$o, bl$tt$o)
    set.seed(123); idx <- caret::createDataPartition(yv, p = 0.8, list = FALSE); tr <- dd[idx, ]; te <- dd[-idx, ]
    met <- function(tt) { m <- lm(.rg_formula(b$y, tt, "o"), data = tr); pr <- predict(m, te); ob <- te[[b$y]]
      data.frame(R2 = caret::R2(pr, ob), RMSE = caret::RMSE(pr, ob), MAE = caret::MAE(pr, ob), PER = caret::RMSE(pr, ob) / mean(ob)) }
    vsp <- tryCatch(rbind(if (!same) cbind(Modelo = paste("Hacia atrás:", .rg_ftxt(b$y, bw$tt)), met(bw$tt)), cbind(Modelo = paste("Final:", .rg_ftxt(b$y, bl$tt)), met(bl$tt))), error = function(e) NULL)
    set.seed(123); cv <- tryCatch(caret::train(.rg_formula(b$y, bl$tt, "o"), data = dd, method = "lm", trControl = caret::trainControl(method = "cv", number = 10))$results, error = function(e) NULL)
    cvr <- NULL
    if (isTRUE(opts$cv_rep)) {
      cvf <- function(tt) { set.seed(123); caret::train(.rg_formula(b$y, tt, "o"), data = dd, method = "lm", trControl = caret::trainControl(method = "repeatedcv", number = 10, repeats = 20))$results[, c("RMSE", "Rsquared", "MAE")] }
      cvr <- tryCatch(rbind(if (!same) cbind(Modelo = paste("Hacia atrás:", .rg_ftxt(b$y, bw$tt)), cvf(bw$tt)), cbind(Modelo = paste("Final:", .rg_ftxt(b$y, bl$tt)), cvf(bl$tt))), error = function(e) NULL)
    }
    pred <- tryCatch({ pi <- predict(fin, newdata = head(te, 6), interval = "confidence"); data.frame(Fila = rownames(head(te, 6)), head(te[, unique(bl$tt$var), drop = FALSE], 6), Predicho = pi[, "fit"], `IC inf` = pi[, "lwr"], `IC sup` = pi[, "upr"], check.names = FALSE, row.names = NULL) }, error = function(e) NULL)
    dff <- dffits(full)
    st <- list(n = nrow(dd), k_full = length(coef(full)) - 1, preds_txt = paste(b$allp, collapse = ", "), poly = b$poly, formula = .rg_ftxt(b$resp, bl$tt),
               F = unname(fs[1]), df1 = unname(fs[2]), df2 = unname(fs[3]), p = pF, R2 = s$r.squared, aR2 = s$adj.r.squared, AIC = AIC(fin), BIC = BIC(fin),
               vif_max = if (!is.null(vif_f)) vif_f$max else NA, dw = if (!is.null(dw)) dw$dw else NA, dw_p = if (!is.null(dw)) dw$p else NA, sw = sw, eta_sum = eta1,
               cv_RMSE = if (!is.null(cv)) cv$RMSE else NA, cv_R2 = if (!is.null(cv)) cv$Rsquared else NA, cv_MAE = if (!is.null(cv)) cv$MAE else NA,
               sp_R2 = if (!is.null(vsp)) tail(vsp$R2, 1) else NA, sp_MAE = if (!is.null(vsp)) tail(vsp$MAE, 1) else NA)
    tabs <- list(
      "Modelo inicial (todos los predictores): coeficientes e IC 95 %" = .coef_table(full),
      "Modelo inicial: ajuste global" = data.frame(F = unname(fsf[1]), gl1 = unname(fsf[2]), gl2 = unname(fsf[3]), p = pf(fsf[1], fsf[2], fsf[3], lower.tail = FALSE), `R²` = sf$r.squared, `R² ajustado` = sf$adj.r.squared, check.names = FALSE),
      "Observaciones con DFFITS > 0.5 (modelo inicial)" = if (any(dff > 0.5)) data.frame(Fila = names(dff)[dff > 0.5], DFFITS = dff[dff > 0.5], dd[dff > 0.5, c(b$resp, b$allp), drop = FALSE], check.names = FALSE, row.names = NULL) else NULL,
      "Valores únicos de los predictores numéricos" = .rg_unique(b),
      "GAM: efectos no paramétricos (tendencia curvilínea si p < 0.05)" = if (!is.null(g)) g$table else NULL,
      "Modelo final: coeficientes e IC 95 %" = .coef_table(fin),
      "Modelo final: ajuste global" = data.frame(Fórmula = st$formula, F = st$F, gl1 = st$df1, gl2 = st$df2, p = pF, `R²` = st$R2, `R² ajustado` = st$aR2, AIC = st$AIC, BIC = st$BIC, check.names = FALSE),
      "Modelo final: Durbin-Watson (set.seed(123))" = if (!is.null(dw)) data.frame(Autocorrelación = dw$r, `D-W` = dw$dw, p = dw$p, check.names = FALSE) else NULL,
      "Tamaño del efecto (etaSquared, lsr)" = eta,
      "Validación: partición 80/20 (set.seed(123); métricas sobre el conjunto de prueba)" = vsp,
      "Validación cruzada 10-fold del modelo final (set.seed(123))" = if (!is.null(cv)) data.frame(cv[, c("RMSE", "Rsquared", "MAE")], `Error porcentual (RMSE/media, %)` = 100 * cv$RMSE / mean(yv), check.names = FALSE) else NULL,
      "Validación cruzada repetida (10 particiones × 20 repeticiones)" = cvr,
      "Predicción con IC 95 % (primeros sujetos del conjunto de prueba)" = pred)
    list(main = fin, full = full, gam = g, b = b, stats = st, tables = Filter(Negate(is.null), tabs),
         simp = list(sat = sat, bw = bw, bl = bl, modS = modS, modB = modB, lp = lp, cmp = cmp_tab),
         notes = paste0("Factores: ", if (length(b$facs)) paste(b$facs, collapse = ", ") else "ninguno indicado", ". Validación cruzada y partición con set.seed(123), como el libro. Si el GAM sugiere curvatura, indique los términos polinómicos en las opciones."))
  },
  posthoc = function(res, d, v, opts) {
    sp <- res$simp; b <- res$b
    a_bas <- anova(sp$modB, sp$modS)
    tabs <- list(
      "Modelo saturado: Anova tipo II" = .rg_anova_df(car::Anova(sp$sat, type = "II")),
      "Eliminación hacia atrás (un término a la vez, mayor p)" = sp$bw$steps,
      "Modelo simplificado: Anova tipo II" = .rg_anova_df(sp$bw$anova),
      "Modelo simplificado vs. modelo basal (anova, AIC, BIC)" = data.frame(Modelo = c("Basal (~ 1)", .rg_ftxt(b$y, sp$bw$tt)), `gl residuales` = a_bas$Res.Df, RSS = a_bas$RSS, F = a_bas$F, p = a_bas$`Pr(>F)`,
                                                                          AIC = c(AIC(sp$modB), AIC(sp$modS)), BIC = c(BIC(sp$modB), BIC(sp$modS)), check.names = FALSE),
      "Efecto de las variables no significativas que quedaron por jerarquía (anova de modelos anidados)" = sp$bl$table,
      "Todos los subconjuntos (regsubsets): mejor modelo por criterio" = if (!is.null(sp$lp)) sp$lp$table else NULL,
      "Comparación de modelos candidatos (R² ajustado y AIC)" = sp$cmp)
    list(tables = Filter(Negate(is.null), tabs), notes = "Procedimiento del libro: Anova tipo II, eliminación del término de mayor p (con la jerarquía polinómica), anova()/AIC/BIC entre modelos anidados y regsubsets. La elección final debe combinar criterio clínico y parsimonia.")
  },
  plots = list(
    diag = list(label = "Diagnósticos del modelo inicial (residuos, Q-Q, escala-localización, leverage)", fn = function(d, v, res, cos) plot_diag(res$full, c(1, 2, 3, 5))),
    cook0 = list(label = "Distancia de Cook (modelo inicial)", fn = function(d, v, res, cos) function() plot(res$full, 4)),
    gam = list(label = "GAM: efectos parciales suavizados s()", fn = function(d, v, res, cos) { g <- res$gam
      if (is.null(g)) return(function() { plot.new(); text(0.5, 0.5, "Sin predictores continuos para suavizar.") })
      function() { op <- par(mfrow = c(1, length(g$sterms))); on.exit(par(op)); plot(g$model, se = TRUE, terms = g$sterms) } }),
    poly = list(label = "Tendencias polinómicas bivariadas (R²)", fn = function(d, v, res, cos) .rg_poly_plot(res$b, cos)),
    cook = list(label = "Distancia de Cook (modelo final)", fn = function(d, v, res, cos) function() plot(res$main, 4, id.n = 5))),
  opts_ui = c(.rg_opts_common, list(cv_rep = list(label = "Validación cruzada repetida 10 × 20 (como el libro; tarda unos segundos)", type = "checkbox", value = TRUE))),
  report = function(res, d, v, checks) {
    s <- res$stats; vif <- .find(checks, "vif"); nn <- .find(checks, "n")
    sw <- s$sw; norm_txt <- if (is.na(sw$met)) "se evaluaron con el gráfico Q-Q" else if (isTRUE(sw$met)) "no mostraron desviación de la normalidad" else "mostraron desviación de la normalidad"
    norm_how <- if (identical(sw$rule, "skew_kurt")) paste0(" (n > 50: sesgo y curtosis, |skew.2SE| = ", fmt_num(abs(sw$skew2SE)), ", |kurt.2SE| = ", fmt_num(abs(sw$kurt2SE)), "; véase el gráfico Q-Q)") else if (!is.na(sw$p)) paste0(" (Shapiro-Wilk, p ", fmt_p(sw$p), "; véase el gráfico Q-Q)") else ""
    paste0("Se construyó un modelo de regresión lineal múltiple para predecir ", CTX(res$b$resp), " a partir de ", CTX("factores estudiados"), ", sobre un set de datos de ", s$n, " observaciones. ",
           "Partiendo de un modelo con todos los predictores (", s$preds_txt, "), cuyo número de parámetros (", s$k_full, " gl) ", if (isTRUE(nn$met)) "es adecuado" else "no es adecuado", " para el tamaño de muestra, según la regla de 10 a 20 sujetos por grado de libertad, ",
           .rg_vif_txt(vif), if (!is.na(vif$met)) " y " else "", .rg_poly_txt(s$poly), ". Tras la simplificación por eliminación hacia atrás, el modelo final fue:\n\n`", s$formula, "`\n\n",
           "El modelo ", if (s$p < ALPHA) "resultó significativo" else "no resultó significativo", " (F(", fmt_num(s$df1, 0), ", ", fmt_num(s$df2, 0), ") = ", fmt_num(s$F), ", p ", fmt_p(s$p), "), con R² = ", fmt_num(s$R2, 3), " y R² ajustado = ", fmt_num(s$aR2, 3), ". ",
           if (!is.na(s$vif_max)) paste0("El VIF máximo fue ", fmt_num(s$vif_max), if (s$vif_max < 5) " (< 5), lo que descarta problemas de multicolinealidad. " else ", lo que indica multicolinealidad. ") else "",
           "Los residuos ", norm_txt, norm_how,
           if (!is.na(s$dw)) paste0("; la prueba de Durbin-Watson (D-W = ", fmt_num(s$dw), ", p ", fmt_p(s$dw_p), ") ", if (s$dw_p < ALPHA) paste0("señaló autocorrelación ", CTX("si los datos son transversales, atribuible al orden de las filas")) else "no señaló autocorrelación") else "", ". ",
           if (!is.na(s$eta_sum)) paste0("Los tamaños de efecto (η²) de los predictores ", CTX("describir su magnitud"), " y, en conjunto, equivalen al R² del modelo (", fmt_num(s$eta_sum), "). ") else "",
           "El modelo se validó por dos métodos: (a) validación cruzada de 10 particiones (RMSE = ", fmt_num(s$cv_RMSE, 1), ", R² = ", fmt_num(s$cv_R2, 3), " y error absoluto medio (MAE) = ", fmt_num(s$cv_MAE, 1), "); y (b) partición 80/20 (entrenamiento/prueba): sobre el conjunto de prueba, el R² (correlación al cuadrado entre valores predichos y observados) fue de ", fmt_num(s$sp_R2, 3), ", con un MAE de ", fmt_num(s$sp_MAE, 1), ". ",
           "En conjunto, ", CTX(res$b$resp), " se explica ", CTX("de forma modesta/sustancial"), " por los predictores del modelo final: el modelo capta alrededor del ", round(100 * s$R2), " % de su variabilidad, con un error de predicción típico (MAE) cercano a ", fmt_num(s$cv_MAE, 1), " (Field et al., 2012; Vittinghoff et al., 2012).")
  }
))

# =============================================================================
# REGRESIÓN LOGÍSTICA
# =============================================================================
.lg_pseudo <- function(m) { dev <- m$deviance; nd <- m$null.deviance; n <- length(m$fitted.values)   # LogisticPseudoR2s() del libro
  rcs <- 1 - exp(-(nd - dev) / n); c(hl = 1 - dev / nd, cs = rcs, n = rcs / (1 - exp(-nd / n))) }

ado_register(list(
  id = "logistica", family = "regresion", order = 2, label = "Regresión logística",
  description = "Predice una respuesta binaria (glm binomial). Flujo del libro: modelo inicial, EPV, VIF, simplificación hacia atrás, pseudo-R², bondad de ajuste (le Cessie), sobredispersión, función de enlace, ROC/AUC, punto de corte, Hosmer-Lemeshow, validación, nomograma.",
  chapter = "Regresión logística en R",
  roles = list(response = role("Variable de respuesta (binaria)", "binary"),
               predictors = role("Predictores continuos (o ya categóricos en texto)", "any", multiple = TRUE, min = 1),
               factors = role("Predictores categóricos codificados con números (se convierten a factor)", "factor", multiple = TRUE, required = FALSE,
                              help = "El libro especifica los factores al inicio (as.factor). Las variables numéricas que no se indiquen aquí se tratan como continuas.")),
  example = list(file = "datAsthma.txt", vars = list(response = "Bruxism", predictors = c("Asthma", "Age", "Sex", "OrParaF")), note = "Bruxismo en función de asma, edad, sexo y parafunción oral."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, c(v$predictors, v$factors), v$response),
  assumptions = function(d, v) {
    b <- .rg_build(d, v, list(), binary = TRUE); y <- b$dd$.y
    m <- glm(b$f_full, data = b$dd, family = binomial(link = "logit"))
    df <- length(coef(m)) - 1; ev <- min(sum(y == 1), sum(y == 0)); P <- ev / length(y); epv <- ev / max(df, 1); nreq <- 10 * df / P; n_ok <- length(y) >= nreq
    checks <- list(
      check_obj("binaria", "Respuesta binaria", TRUE, paste0("Respuesta codificada 0/1 (evento = ", if (is.numeric(d[[v$response]])) "valor máximo" else levels(factor(b$dd[[v$response]]))[2], "; la otra categoría es la referencia); ", sum(y == 1), " eventos y ", sum(y == 0), " no eventos.")),
      check_obj("epv", "Eventos por variable (EPV ≥ 10; Peduzzi)", n_ok,
                paste0("Modelo inicial con ", df, " gl (parámetros). Evento menos frecuente: ", ev, " (P = ", fmt_num(P, 3), "); EPV = ", fmt_num(epv, 1), ". Según Peduzzi (n = 10·gl/P) se requieren ≈ ", round(nreq), " sujetos y hay ", length(y), ": ",
                       if (n_ok) "tamaño adecuado." else "el modelo inicial está sobreparametrizado (priorizar el modelo más parsimonioso)."), epv = epv, nreq = nreq, df = df),
      if (df > 1) check_vif(m) else check_obj("vif", "Multicolinealidad (VIF)", NA, "Un solo predictor: VIF no aplica."))
    list(checks = checks, decision = decide_path(checks, FALSE))
  },
  run = function(d, v, opts) {
    b <- .rg_build(d, v, opts, binary = TRUE); dd <- b$dd; y <- dd$.y; n <- nrow(dd)
    fit <- function(f) glm(f, data = dd, family = binomial(link = "logit"))
    full <- glm(b$f_full, data = dd, family = binomial(link = "logit"))
    chi0 <- full$null.deviance - full$deviance; df0 <- full$df.null - full$df.residual
    g <- .rg_gam(b, binomial())
    sat <- fit(.rg_formula(b$y, b$tt)); bw <- .rg_backward(fit, b$y, b$tt); bl <- .rg_blocks(fit, b$y, bw, test = "Chisq")
    fin <- glm(.rg_formula(b$y, bl$tt, "o"), data = dd, family = binomial(link = "logit")); modB <- glm(.y ~ 1, data = dd, family = binomial(link = "logit"))
    lp <- tryCatch(.rg_leaps(b$y, b$tt, dd), error = function(e) NULL)
    chi <- fin$null.deviance - fin$deviance; chdf <- fin$df.null - fin$df.residual; pchi <- pchisq(chi, chdf, lower.tail = FALSE)
    pr2 <- .lg_pseudo(fin); ffo <- .rg_formula(b$y, bl$tt, "o")
    gof <- tryCatch(resid(rms::lrm(ffo, x = TRUE, y = TRUE, data = dd), "gof"), error = function(e) NULL)
    od <- tryCatch({ invisible(capture.output(o <- qcc::qcc.overdispersion.test(c(0, 1), as.vector(table(y))))); o }, error = function(e) NULL)
    links <- lapply(c(logit = "logit", probit = "probit", cloglog = "cloglog"), function(l) tryCatch(glm(ffo, data = dd, family = binomial(link = l)), error = function(e) NULL))
    links <- Filter(Negate(is.null), links)
    prob <- fitted(fin); roc <- pROC::roc(y, prob, levels = c(0, 1), direction = "<", quiet = TRUE); auc <- as.numeric(pROC::auc(roc))
    cut <- tryCatch(pROC::coords(roc, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity")), error = function(e) NULL)
    thr <- if (!is.null(cut)) cut$threshold[1] else 0.5
    clas <- function(t) { cl <- as.integer(prob >= t); c(Sensibilidad = sum(cl == 1 & y == 1) / sum(y == 1), Especificidad = sum(cl == 0 & y == 0) / sum(y == 0), Error = mean(cl != y)) }
    cm <- table(Predicho = factor(as.integer(prob >= thr), levels = 0:1), Observado = factor(y, levels = 0:1))
    hl <- tryCatch(ResourceSelection::hoslem.test(y, prob, g = 10), error = function(e) NULL)
    # Validación por división de la muestra (set.seed(123)): modelo final vs. modelo inicial completo
    set.seed(123); idx <- caret::createDataPartition(factor(y), p = 0.8, list = FALSE); tr <- dd[idx, ]; te <- dd[-idx, ]
    auc_te <- function(f) { m <- glm(f, data = tr, family = binomial(link = "logit")); pr <- predict(m, newdata = te, type = "response"); pROC::roc(te$.y, pr, levels = c(0, 1), direction = "<", quiet = TRUE) }
    roc_te <- tryCatch(list(final = auc_te(ffo), inicial = auc_te(b$f_full)), error = function(e) NULL)
    val <- NULL
    if (isTRUE(opts$boot)) val <- tryCatch({ set.seed(42); vb <- rms::validate(rms::lrm(ffo, x = TRUE, y = TRUE, data = dd), B = 200)
      data.frame(Índice = rownames(vb), unclass(vb)[, c("index.orig", "training", "test", "optimism", "index.corrected", "n")], check.names = FALSE, row.names = NULL) }, error = function(e) NULL)
    ci_lr <- tryCatch(car::Confint(fin, level = 0.95, type = "LR"), error = function(e) NULL)
    nd <- tryCatch(.rg_newdata(opts$nuevos, dd, unique(bl$tt$var)), error = function(e) e)
    pnew <- if (inherits(nd, "error")) data.frame(Aviso = conditionMessage(nd)) else if (!is.null(nd)) data.frame(nd, `Probabilidad predicha` = predict(fin, newdata = nd, type = "response"), check.names = FALSE) else NULL
    dfull <- length(coef(full)) - 1
    st <- list(n = n, nvars = ncol(d), chi = chi, chdf = chdf, pchi = pchi, r2 = pr2, auc = auc, thr = thr, formula = .rg_ftxt(b$resp, bl$tt),
               gof_z = if (!is.null(gof)) unname(gof["Z"]) else NA, gof_p = if (!is.null(gof)) unname(gof["P"]) else NA,
               od_p = if (!is.null(od)) od[1, "p-value"] else NA, hl_p = if (!is.null(hl)) hl$p.value else NA, hl_chi = if (!is.null(hl)) unname(hl$statistic) else NA,
               auc_te = if (!is.null(roc_te)) as.numeric(pROC::auc(roc_te$final)) else NA, terms = unique(bl$tt$var))
    tabs <- list(
      "Modelo inicial (todos los predictores): coeficientes y OR (IC 95 %)" = .coef_table(full, expo = TRUE),
      "Modelo inicial: razón de verosimilitud" = data.frame(`χ²` = chi0, gl = df0, p = pchisq(chi0, df0, lower.tail = FALSE), check.names = FALSE),
      "Valores únicos de los predictores numéricos" = .rg_unique(b),
      "GAM (binomial): efectos no paramétricos (tendencia curvilínea si p < 0.05)" = if (!is.null(g)) g$table else NULL,
      "Modelo final: coeficientes (escala logOR)" = .coef_table(fin),
      "Modelo final: logOR y OR con IC 95 % por razón de verosimilitud (car::Confint, type = \"LR\")" = if (!is.null(ci_lr)) data.frame(Término = rownames(ci_lr), logOR = ci_lr[, 1], `IC inf` = ci_lr[, 2], `IC sup` = ci_lr[, 3], OR = exp(ci_lr[, 1]), `OR IC inf` = exp(ci_lr[, 2]), `OR IC sup` = exp(ci_lr[, 3]), check.names = FALSE, row.names = NULL) else NULL,
      "Modelo final: razón de verosimilitud y pseudo-R²" = data.frame(Fórmula = st$formula, `χ²` = chi, gl = chdf, p = pchi, `R² Hosmer-Lemeshow` = pr2[["hl"]], `R² Cox-Snell` = pr2[["cs"]], `R² Nagelkerke` = pr2[["n"]], AIC = AIC(fin), BIC = BIC(fin), check.names = FALSE),
      "Bondad de ajuste (le Cessie-van Houwelingen; rms::resid(lrm, \"gof\"))" = if (!is.null(gof)) data.frame(`Suma de cuadrados del error` = gof[1], `Valor esperado|H0` = gof[2], DE = gof[3], Z = gof[4], P = gof[5], check.names = FALSE, row.names = NULL) else NULL,
      "Sobredispersión (qcc.overdispersion.test)" = if (!is.null(od)) data.frame(`Obs.Var/Theor.Var` = od[1, 1], Estadístico = od[1, 2], p = od[1, 3], check.names = FALSE) else NULL,
      "Función de enlace: comparación (devianza residual, AIC, BIC)" = data.frame(Enlace = names(links), `Devianza residual` = sapply(links, deviance), AIC = sapply(links, AIC), BIC = sapply(links, BIC), check.names = FALSE, row.names = NULL),
      "Discriminación: AUC y punto de corte óptimo (Youden)" = data.frame(AUC = auc, `Umbral óptimo` = thr, Sensibilidad = if (!is.null(cut)) cut$sensitivity[1] else NA, Especificidad = if (!is.null(cut)) cut$specificity[1] else NA, check.names = FALSE),
      "Matriz de confusión (punto de corte óptimo)" = data.frame(Predicho = c("0", "1"), `Observado 0` = cm[, "0"], `Observado 1` = cm[, "1"], check.names = FALSE, row.names = NULL),
      "Sensibilidad, especificidad y error según el punto de corte" = data.frame(Umbral = c("0.5", "0.2", paste0("Óptimo (", round(thr, 3), ")")), rbind(clas(0.5), clas(0.2), clas(thr)), check.names = FALSE),
      "Calibración: prueba de Hosmer-Lemeshow (g = 10)" = if (!is.null(hl)) data.frame(`χ²` = st$hl_chi, gl = unname(hl$parameter), p = st$hl_p, check.names = FALSE) else NULL,
      "Validación por división 80/20 (AUC en el conjunto de prueba)" = if (!is.null(roc_te)) data.frame(Modelo = c(paste("Final:", st$formula), paste("Inicial:", b$resp, "~", paste(b$allp, collapse = " + "))), `AUC (prueba)` = c(as.numeric(pROC::auc(roc_te$final)), as.numeric(pROC::auc(roc_te$inicial))), check.names = FALSE) else NULL,
      "Validación interna por bootstrap (rms::validate, B = 200, set.seed(42))" = val,
      "Predicción para sujetos nuevos" = pnew)
    list(main = fin, full = full, gam = g, b = b, roc = roc, roc_te = roc_te, stats = st, tables = Filter(Negate(is.null), tabs),
         simp = list(sat = sat, bw = bw, bl = bl, modB = modB, lp = lp),
         notes = paste0("Factores: ", if (length(b$facs)) paste(b$facs, collapse = ", ") else "ninguno indicado", ". OR = exp(β); recuerde que es una razón de momios, no de probabilidades. AUC corregido por optimismo = 0.5 + Dxy/2 (fila Dxy del bootstrap). La prueba de sobredispersión es solo complementaria con respuesta binaria individual."))
  },
  posthoc = function(res, d, v, opts) {
    sp <- res$simp; b <- res$b; fin <- res$main
    an <- anova(sp$modB, fin, test = "Chisq")
    tabs <- list(
      "Modelo saturado: Anova tipo II (razón de verosimilitud)" = .rg_anova_df(car::Anova(sp$sat, type = "II")),
      "Eliminación hacia atrás (un término a la vez, mayor p)" = sp$bw$steps,
      "Modelo simplificado: Anova tipo II" = .rg_anova_df(sp$bw$anova),
      "Efecto de las variables no significativas que quedaron por jerarquía (anova de modelos anidados)" = sp$bl$table,
      "Modelo final vs. modelo basal (anova, AIC, BIC)" = data.frame(Modelo = c("Basal (~ 1)", res$stats$formula), `gl residuales` = an$`Resid. Df`, `Devianza residual` = an$`Resid. Dev`, Devianza = an$Deviance, p = an$`Pr(>Chi)`,
                                                                    AIC = c(AIC(sp$modB), AIC(fin)), BIC = c(BIC(sp$modB), BIC(fin)), check.names = FALSE),
      "Todos los subconjuntos (regsubsets): mejor modelo por criterio" = if (!is.null(sp$lp)) sp$lp$table else NULL)
    list(tables = Filter(Negate(is.null), tabs), notes = "Los criterios automáticos (regsubsets) pueden no coincidir entre sí ni con la simplificación: complementar con criterio clínico y parsimonia.")
  },
  plots = list(
    roc = list(label = "Curva ROC del modelo final", fn = function(d, v, res, cos) plot_roc(res$roc, cos)),
    roc_val = list(label = "Validación: curvas ROC en datos de prueba (80/20)", fn = function(d, v, res, cos) { cos <- cos_merge(cos); r <- res$roc_te
      function() { if (is.null(r)) { plot.new(); text(0.5, 0.5, "No fue posible la partición 80/20."); return(invisible()) }
        plot(r$final, col = cos$colour[1], lwd = 2, main = "Validación: curvas ROC en datos de prueba")
        plot(r$inicial, col = cos$colour[2], lwd = 2, lty = 2, add = TRUE)
        legend("bottomright", bty = "n", col = cos$colour[1:2], lwd = 2, lty = 1:2,
               legend = c(paste0("Modelo final (AUC = ", round(pROC::auc(r$final), 3), ")"), paste0("Modelo inicial (AUC = ", round(pROC::auc(r$inicial), 3), ")"))) } }),
    nomo = list(label = "Nomograma (rms)", fn = function(d, v, res, cos) { dd <- res$b$dd
      assign(".ado_rg_ddist", rms::datadist(dd[setdiff(names(dd), ".y")]), envir = globalenv()); op <- options(datadist = ".ado_rg_ddist"); on.exit(options(op))
      nom <- rms::nomogram(rms::lrm(.rg_formula(res$b$y, res$simp$bl$tt, "o"), data = dd), fun = plogis, funlabel = "Probabilidad del evento")
      function() plot(nom) }),
    dca = list(label = "Análisis de curva de decisión (DCA)", fn = function(d, v, res, cos) { dd <- data.frame(y = res$main$y, pred = fitted(res$main))
      plot(dcurves::dca(y ~ pred, data = dd, thresholds = seq(0, 0.5, 0.01)), smooth = TRUE) })),
  opts_ui = c(.rg_opts_common, list(
    boot = list(label = "Validación interna por bootstrap (rms::validate, B = 200; como el libro)", type = "checkbox", value = TRUE),
    nuevos = list(label = "Sujetos nuevos para predecir (var=valor separados por comas; sujetos separados por «;»)", type = "text", value = ""))),
  report = function(res, d, v, checks) {
    s <- res$stats; vif <- .find(checks, "vif"); r2 <- s$r2
    paste0("El set de datos empleado incluye ", s$nvars, " variables y ", s$n, " observaciones, y el objetivo del modelo es predecir ", CTX(res$b$resp), " a partir de ", CTX("distintos factores clínicos"), ". ",
           "Tras evaluar la multicolinealidad", if (!is.na(vif$met)) paste0(" (", if (isTRUE(vif$met)) "VIF < 5" else paste0("VIF máximo = ", fmt_num(vif$max)), ")") else "", " y simplificar el modelo mediante eliminación hacia atrás (backward), el modelo final fue `", s$formula, "`. ",
           "El modelo ", if (s$pchi < ALPHA) "resultó significativo" else "no resultó significativo", " en la prueba de razón de verosimilitud (χ² = ", fmt_num(s$chi), ", gl = ", s$chdf, ", p ", fmt_p(s$pchi), "). Los pseudo-R² oscilaron entre ", fmt_num(min(r2)), " y ", fmt_num(max(r2)), ". ",
           if (!is.na(s$gof_p)) paste0("La bondad de ajuste, evaluada con la prueba de le Cessie-Hosmer, ", if (s$gof_p > ALPHA) "no mostró evidencia de falta de ajuste" else "mostró evidencia de falta de ajuste", " (z = ", fmt_num(s$gof_z), ", p ", fmt_p(s$gof_p), ")") else "",
           if (!is.na(s$od_p)) paste0(", y ", if (s$od_p > ALPHA) "tampoco se detectó sobredispersión" else "se detectó sobredispersión", " (p ", fmt_p(s$od_p), ")") else "", ". ",
           "El modelo alcanzó un AUC de ", fmt_num(s$auc, 3), if (!is.na(s$hl_p)) paste0(" y una calibración ", if (s$hl_p > ALPHA) "adecuada" else "deficiente", " según la prueba de Hosmer-Lemeshow (p ", fmt_p(s$hl_p), ")") else "", ". ",
           "Los factores asociados de forma significativa con ", CTX(res$b$resp), " fueron ", paste(s$terms, collapse = ", "), ", con sus respectivas razones de momios e intervalos de confianza ", CTX("OR (IC 95 %) de cada factor"), " (Harrell, 2015; Hosmer et al., 2013).")
  }
))
