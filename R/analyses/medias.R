# =============================================================================
# analyses/medias.R - Familia "Comparación de medias" (caps. 15, 16, 18, 19,
# 23, 24, 25, 26, 27, 28–34): t independiente, t apareada, ANOVA de una vía,
# ANOVA factorial, ANCOVA y sus alternativas no paramétricas (Wilcoxon,
# Kruskal-Wallis). ANOVA de medidas repetidas y Friedman → medias_rm.R.
# Fiel a las llamadas, tablas, gráficos y plantillas de reporte del libro.
# =============================================================================

r_from_t <- function(t, df) sqrt(t^2 / (t^2 + df))
r_from_z <- function(p, N) abs(qnorm(p / 2)) / sqrt(N)

# Frases de supuestos para el reporte (a partir de checks) - usadas también por otras familias
.sw_txt  <- function(ch) {
  if (is.null(ch)) return("")
  s <- paste0("W = ", fmt_num(ch$W, 2), ", p ", fmt_p(ch$p))
  if (identical(ch$rule, "skew_kurt")) s <- paste0(s, "; con n > 50: skew.2SE = ", fmt_num(ch$skew2SE), ", kurt.2SE = ", fmt_num(ch$kurt2SE))
  s
}
.lev_txt <- function(ch) if (is.null(ch)) "" else paste0("F = ", fmt_num(ch$F, 3), ", p ", fmt_p(ch$p))
.find <- function(checks, kind) { for (ch in checks) if (identical(ch$kind, kind)) return(ch); NULL }
.msd <- function(v) paste0("Media = ", fmt_num(mean(v, na.rm = TRUE)), ", DE = ", fmt_num(sd(v, na.rm = TRUE)))
.mmd <- function(v) paste0("mediana = ", fmt_num(median(v, na.rm = TRUE)), "; media = ", fmt_num(mean(v, na.rm = TRUE)), " ± ", fmt_num(sd(v, na.rm = TRUE)))

# ---- Utilidades internas de esta familia (prefijo .md_) ------------------------------
# Veredicto de normalidad según la regla aplicada (Cap. 9): Shapiro-Wilk (n ≤ 50) o sesgo/curtosis (n > 50)
.md_norm_how <- function(ch) if (identical(ch$rule, "skew_kurt")) "la prueba de Shapiro-Wilk, el sesgo y la curtosis estandarizados (n > 50)" else "la prueba de Shapiro-Wilk"
.md_norm_verdict <- function(ch, yes = ", sin evidencia de distribución diferente a la normal",
                             no = ", con evidencia de desviación de la normalidad", na = "") {
  if (is.null(ch) || is.na(ch$met)) na else if (isTRUE(ch$met)) yes else no
}
.md_join <- function(x, last = " y ") { x <- as.character(x); n <- length(x)
  if (n <= 1) paste(x, collapse = "") else paste0(paste(x[-n], collapse = ", "), last, x[n]) }
.md_dd <- function(d, v) { dd <- data.frame(y = num(d, v$response), g = fac(d, v$group)); dd <- dd[complete.cases(dd), ]; dd$g <- droplevels(dd$g); dd }
.md_lev <- function(y, g) tryCatch(car::leveneTest(y ~ g, center = median)[1, "Pr(>F)"], error = function(e) NA)

# Descriptiva por grupo como en las tablas dplyr::summarise del libro
.md_desc <- function(y, g, cols = c("n", "Media", "DE", "Mediana", "IQR"), lab = "Grupo") {
  ok <- !is.na(y) & !is.na(g); y <- y[ok]; g <- droplevels(factor(g[ok]))
  f <- function(v) c(n = length(v), Media = mean(v), DE = sd(v), EE = sd(v) / sqrt(length(v)),
                     Mediana = median(v), IQR = IQR(v), Min = min(v), Max = max(v))
  m <- t(sapply(split(y, g), f))[, cols, drop = FALSE]
  out <- data.frame(levels(g), round(m, 3), row.names = NULL, check.names = FALSE); names(out)[1] <- lab; out
}

# Grupo control elegido por el usuario (texto); vacío o inexistente → primer nivel del factor
.md_ctrl <- function(g, ctrl) { lv <- levels(g); ctrl <- trimws(if (is.null(ctrl)) "" else as.character(ctrl)[1])
  if (nzchar(ctrl) && ctrl %in% lv) ctrl else lv[1] }

# Contrastes ortogonales del libro (caps. 18–19, 24, 27): C1 = control (−(k−1)) vs. resto (+1);
# los siguientes comparan los grupos restantes entre sí (Helmert), con 0 para el control.
# Con 3 niveles y control primero reproduce cbind(c(-2,1,1), c(0,-1,1)); con 2 niveles, c(-1, 1).
.md_orth <- function(lv, ctrl) {
  k <- length(lv); oth <- setdiff(lv, ctrl)
  M <- matrix(0, k, k - 1, dimnames = list(lv, NULL))
  M[, 1] <- ifelse(lv == ctrl, -(k - 1), 1)
  nm <- paste0(ctrl, "_vs_", paste(oth, collapse = "+"))
  if (k > 2) { H <- contr.helmert(k - 1)
    for (j in seq_len(k - 2)) { M[oth, j + 1] <- H[, j]; nm <- c(nm, paste0(paste(oth[1:j], collapse = "+"), "_vs_", oth[j + 1])) } }
  colnames(M) <- nm; M
}

# Coeficientes personalizados: "-2,1,1; 0,-1,1" (un contraste por bloque, en el orden de los niveles)
.md_parse_contr <- function(txt, lv) {
  if (is.null(txt) || !nzchar(trimws(txt))) stop("Escriba los coeficientes de los contrastes personalizados (p. ej. -2,1,1; 0,-1,1).")
  blocks <- strsplit(strsplit(txt, ";")[[1]], ",")
  M <- do.call(cbind, lapply(blocks, function(b) suppressWarnings(as.numeric(trimws(b)))))
  k <- length(lv)
  if (nrow(M) != k || anyNA(M)) stop("Cada contraste debe tener ", k, " coeficientes numéricos (niveles: ", paste(lv, collapse = ", "), ").")
  if (ncol(M) > k - 1) stop("Se admiten como máximo ", k - 1, " contrastes (k − 1).")
  if (any(abs(colSums(M)) > 1e-8)) stop("Los coeficientes de cada contraste deben sumar 0.")
  cp <- crossprod(M); orth <- all(abs(cp[upper.tri(cp)]) < 1e-8)
  dimnames(M) <- list(lv, paste0("Contraste ", seq_len(ncol(M)))); attr(M, "orthogonal") <- orth; M
}

# Tabla de contrastes a partir de summary.lm() (filas idx) con r = √(t²/(t²+gl)) (caps. 19.8 y 27.11)
.md_contr_tab <- function(m, idx, names, with_r = TRUE) {
  cf <- summary.lm(m)$coefficients[idx, , drop = FALSE]; dfR <- df.residual(m)
  out <- data.frame(Contraste = names, Estimación = cf[, 1], EE = cf[, 2], t = cf[, 3], gl = dfR, p = cf[, 4], check.names = FALSE, row.names = NULL)
  if (with_r) out$r <- r_from_t(cf[, 3], dfR)
  out
}

# Comparaciones por pares en formato común: g1, g2, est (= media g2 − media g1), lwr, upr, p
.md_pairs_glht <- function(g, lv) {
  s <- summary(g); ci <- confint(g)$confint; rn <- rownames(ci)
  pr <- t(sapply(rn, function(r) { for (a in lv) for (b in lv) if (a != b && identical(r, paste(b, "-", a))) return(c(a, b)); c(NA, NA) }))
  data.frame(g1 = pr[, 1], g2 = pr[, 2], est = ci[, "Estimate"], lwr = ci[, "lwr"], upr = ci[, "upr"],
             EE = as.numeric(s$test$sigma), t = as.numeric(s$test$tstat), p = as.numeric(s$test$pvalues), row.names = NULL)
}
.md_pairs_tab <- function(pr, padj_lab = "p ajustada") {
  out <- data.frame(Comparación = paste(pr$g2, "-", pr$g1), Estimación = pr$est, check.names = FALSE)
  for (cc in c("EE", "t", "lwr", "upr")) if (!is.null(pr[[cc]]) && !all(is.na(pr[[cc]]))) out[[switch(cc, lwr = "IC 95% inf", upr = "IC 95% sup", cc)]] <- pr[[cc]]
  out[[padj_lab]] <- pr$p; out
}
.md_glht <- function(model, gname, type) {
  lf <- do.call(multcomp::mcp, setNames(list(type), gname))
  multcomp::glht(model, linfct = lf)
}
.glht_mcp <- .md_glht

# Games-Howell (rstatix::games_howell_test, cap. 18.6.1) - implementación nativa equivalente
.md_games_howell <- function(y, g) {
  lv <- levels(g); k <- length(lv); m <- tapply(y, g, mean); v <- tapply(y, g, var); n <- tapply(y, g, length)
  pr <- utils::combn(k, 2)
  do.call(rbind, lapply(seq_len(ncol(pr)), function(c) { i <- pr[1, c]; j <- pr[2, c]
    d <- unname(m[j] - m[i]); se2 <- unname(v[i] / n[i] + v[j] / n[j])
    df <- se2^2 / (unname((v[i] / n[i])^2 / (n[i] - 1) + (v[j] / n[j])^2 / (n[j] - 1)))
    q <- abs(d) / sqrt(se2 / 2); hw <- qtukey(1 - ALPHA, k, df) * sqrt(se2 / 2)
    data.frame(g1 = lv[i], g2 = lv[j], est = d, lwr = d - hw, upr = d + hw, gl = df, p = ptukey(q, k, df, lower.tail = FALSE)) }))
}

# pairwise.t.test → formato común (Δ = media g2 − media g1)
.md_pairs_pw <- function(y, g, adj) {
  pw <- pairwise.t.test(y, g, p.adjust.method = adj); pm <- pw$p.value; m <- tapply(y, g, mean)
  out <- list(); for (cn in colnames(pm)) for (rn in rownames(pm)) if (!is.na(pm[rn, cn]))
    out[[length(out) + 1]] <- data.frame(g1 = cn, g2 = rn, est = unname(m[rn] - m[cn]), p = pm[rn, cn])
  do.call(rbind, out)
}

# Frase post-hoc: "reveló diferencia significativa únicamente entre A y B (Δ = …, p = …)…"
.md_pairs_sentence <- function(pr, ptxt, list_ns = TRUE, unit = "") {
  item <- function(i) { hi <- if (pr$est[i] >= 0) pr$g2[i] else pr$g1[i]; lo <- if (pr$est[i] >= 0) pr$g1[i] else pr$g2[i]
    list(ab = paste0(hi, " y ", lo), par = paste0("(Δ = ", fmt_num(abs(pr$est[i])), unit, ", ", ptxt[i], ")")) }
  s <- which(pr$sig); ns <- which(!pr$sig)
  if (!length(s)) return(paste0("no reveló diferencias significativas entre ningún par de grupos (",
                                paste(sapply(ns, function(i) { it <- item(i); paste(it$ab, it$par) }), collapse = "; "), ")"))
  out <- paste0("reveló ", if (length(s) > 1) "diferencias significativas" else "diferencia significativa",
                if (length(s) == 1 && nrow(pr) > 1 && !list_ns) " únicamente" else "", " ",
                .md_join(sapply(s, function(i) { it <- item(i); paste("entre", it$ab, it$par) })))
  if (list_ns && length(ns)) out <- paste0(out, ", mientras que ",
                .md_join(sapply(ns, function(i) { it <- item(i); paste(it$ab, "no difirieron significativamente", it$par) }), "; "))
  out
}
.md_ptxt <- function(p) paste0("p ", fmt_p(p))
# Cierra el párrafo con la cita del libro: "… texto (cita)."
.md_cite <- function(txt, cite) paste0(sub("[.[:space:]]+$", "", txt), " (", cite, ").")

# =============================================================================
# 1. t de Student independiente  (alt: Wilcoxon suma de rangos) - caps. 14–15
# =============================================================================
ado_register(list(
  id = "t_indep", family = "medias", order = 1,
  label = "Prueba t independiente (dos grupos)",
  description = "Compara la media de dos grupos independientes. Alternativa: Wilcoxon suma de rangos.",
  chapter = "Prueba t independiente en R",
  roles = list(response = role("Variable de respuesta (numérica)", "numeric"),
               group    = role("Variable de grupo (2 niveles)", "factor", nlevels = 2)),
  example = list(file = "VickersF.txt", prep = function(d) subset(d, Cure.time == "Dual"),
                 vars = list(response = "VH", group = "Group"),
                 note = "Dureza Vickers (VH) de dos cementos (Group) en protocolo Dual."),
  alternative = "wilcox_rs",
  descriptives = function(d, v) desc_table1(d, v$response, v$group),
  assumptions = function(d, v) {
    y <- num(d, v$response); g <- fac(d, v$group)
    checks <- list(check_normality(y, v$response), check_homogeneity(y, g))
    dec <- decide_path(checks[1], TRUE)   # decide la normalidad; la heterocedasticidad la cubre Welch (cap. 15.1.5)
    if (isFALSE(checks[[2]]$met)) dec$reason <- paste(dec$reason, "Varianzas heterogéneas (Levene): la corrección de Welch (var.equal = FALSE) las ajusta.")
    list(checks = checks, decision = dec)
  },
  run = function(d, v, opts) {
    dd <- .md_dd(d, v); veq <- isTRUE(opts$var_equal)
    itt <- t.test(y ~ g, data = dd, alternative = "two.sided", var.equal = veq)
    t <- unname(itt$statistic); df <- unname(itt$parameter)
    lv <- levels(dd$g); m <- tapply(dd$y, dd$g, mean); s <- tapply(dd$y, dd$g, sd)
    # Libro 15.1.6: Δ = diff(itt$estimate) (segundo − primer nivel) e IC en la misma dirección
    delta <- unname(diff(itt$estimate)); ci <- c(-itt$conf.int[2], -itt$conf.int[1])
    lab <- paste(lv[2], "−", lv[1])
    list(main = itt,
         stats = list(t = t, df = df, p = itt$p.value, r = r_from_t(t, df), delta = delta, lo = ci[1], hi = ci[2], dlab = lab, var_equal = veq,
                      g1 = lv[1], g2 = lv[2], m1 = m[[1]], m2 = m[[2]], s1 = s[[1]], s2 = s[[2]], n = nrow(dd)),
         tables = list("Estadística descriptiva por grupo" = .md_desc(dd$y, dd$g),
                       "Resultado de t.test()" = data.frame(Prueba = if (veq) "t de Student (var.equal = TRUE)" else "t de Welch (var.equal = FALSE)",
                                                            t = t, gl = df, p = itt$p.value, check.names = FALSE),
                       "Tamaño del efecto y diferencia de medias" = data.frame(r = r_from_t(t, df), Diferencia = lab, `Δ medias` = delta,
                                                            `IC 95% inf` = ci[1], `IC 95% sup` = ci[2], check.names = FALSE)),
         notes = paste0(if (veq) "t de Student clásica (var.equal = TRUE)." else "Corrección de Welch (var.equal = FALSE): opción por defecto y más robusta según el libro.",
                        " t se calcula como ", lv[1], " − ", lv[2], " (salida de t.test); Δ e IC se presentan como ", lab, " (cap. 15.1.6). r = √(t²/(t²+gl))."))
  },
  posthoc = NULL,
  plots = list(
    boxplot = list(label = "Boxplot por grupo", fn = function(d, v, res, cos) plot_box(d, v$response, v$group, cos)),
    qq      = list(label = "Q-Q de la respuesta", fn = function(d, v, res, cos) plot_qq(num(d, v$response), paste("Q-Q:", v$response), cos)),
    effect  = list(label = "Diferencia de medias (IC 95 %)", fn = function(d, v, res, cos) plot_effect(res$stats$delta, res$stats$lo, res$stats$hi, res$stats$dlab, cos))
  ),
  opts_ui = list(var_equal = list(label = "Varianzas iguales (t de Student clásica, var.equal = TRUE)", type = "checkbox", value = FALSE)),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- .find(checks, "normalidad"); lev <- .find(checks, "homogeneidad")
    paste0("La normalidad se evaluó mediante ", .md_norm_how(sw), " (", .sw_txt(sw), ") y la evaluación gráfica (Q-Q plot)", .md_norm_verdict(sw),
           ". La homogeneidad de varianzas se evaluó con la prueba de Levene (", .lev_txt(lev), "). ",
           "La ", CTX("variable de respuesta"), " entre ", CTX(s$g1), " (Media = ", fmt_num(s$m1), ", DE = ", fmt_num(s$s1), ") y ",
           CTX(s$g2), " (Media = ", fmt_num(s$m2), ", DE = ", fmt_num(s$s2), ") ",
           if (s$p < ALPHA) "mostró una diferencia estadísticamente significativa" else "no mostró una diferencia estadísticamente significativa",
           " (", if (s$var_equal) "t de Student" else "t de Welch", ": t = ", fmt_num(s$t), ", df = ", fmt_num(s$df, 1), ", p ", fmt_p(s$p),
           "), con una diferencia de medias (", s$dlab, ") de Δ = ", fmt_num(s$delta),
           " (IC 95 %: ", fmt_num(s$lo), ", ", fmt_num(s$hi), "), que representa un tamaño de efecto ", efecto_r(s$r), " (r = ", fmt_num(s$r), ") (Cohen, 1992).")
  }
))

# =============================================================================
# 2. Wilcoxon suma de rangos (Mann-Whitney) - caps. 29–30
# =============================================================================
ado_register(list(
  id = "wilcox_rs", family = "medias", order = 2,
  label = "Wilcoxon suma de rangos (Mann-Whitney)",
  description = "Alternativa no paramétrica a la prueba t independiente.",
  chapter = "Wilcoxon suma de rangos en R",
  roles = list(response = role("Variable de respuesta (numérica)", "numeric"),
               group    = role("Variable de grupo (2 niveles)", "factor", nlevels = 2)),
  example = list(file = "AUCFil.txt", vars = list(response = "AUC", group = "Group"),
                 note = "AUC de microfiltración en dos composites (Group)."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, v$response, v$group),
  assumptions = function(d, v) {
    y <- num(d, v$response); g <- fac(d, v$group)
    checks <- list(check_normality(y, v$response), check_homogeneity(y, g))
    list(checks = checks, decision = list(parametric = FALSE, all_met = FALSE, violated = character(0), reason = "Prueba no paramétrica: los supuestos se reportan como justificación."))
  },
  run = function(d, v, opts) {
    dd <- .md_dd(d, v); exacto <- identical(opts$metodo, "exacto")
    wmsg <- character(0)
    W <- withCallingHandlers(if (exacto) wilcox.test(y ~ g, data = dd, alternative = "two.sided", exact = TRUE)
                             else wilcox.test(y ~ g, data = dd, alternative = "two.sided", exact = FALSE, correct = TRUE),
                             warning = function(w) { wmsg <<- c(wmsg, conditionMessage(w)); invokeRestart("muffleWarning") })
    N <- nrow(dd); r <- r_from_z(W$p.value, N); lv <- levels(dd$g)
    # Ranking manual (cap. 29.3.2) y r por rangos medios (cap. 30.6)
    rk <- rank(dd$y); Rs <- tapply(rk, dd$g, sum); nk <- tapply(dd$y, dd$g, length); Wk <- Rs - nk * (nk + 1) / 2; Rm <- tapply(rk, dd$g, mean)
    r_rk <- 2 * abs(diff(Rm)) / sum(nk)
    med <- tapply(dd$y, dd$g, median); m <- tapply(dd$y, dd$g, mean); s <- tapply(dd$y, dd$g, sd)
    hay_emp <- length(wmsg) > 0
    list(main = W,
         stats = list(W = unname(W$statistic), p = W$p.value, r = r, r_rk = unname(r_rk), N = N, g1 = lv[1], g2 = lv[2],
                      Rm1 = Rm[[1]], Rm2 = Rm[[2]], med1 = med[[1]], med2 = med[[2]], m1 = m[[1]], m2 = m[[2]], s1 = s[[1]], s2 = s[[2]],
                      dmed = abs(unname(med[2] - med[1])), dmean = abs(unname(m[2] - m[1])), exacto = exacto),
         tables = list("Estadística descriptiva (mediana e IQR; media y DE)" = .md_desc(dd$y, dd$g, c("n", "Mediana", "IQR", "Min", "Max", "Media", "DE")),
                       "Rangos por grupo" = data.frame(Grupo = lv, n = as.numeric(nk), `Suma de rangos` = as.numeric(Rs), `Rango medio` = as.numeric(Rm),
                                                       `W del grupo` = as.numeric(Wk), check.names = FALSE),
                       "Resultado de wilcox.test()" = data.frame(`W (salida de R)` = unname(W$statistic), `W menor (estadístico de prueba)` = min(Wk), p = W$p.value,
                                                                 `r (a partir de Z)` = r, `r (rangos medios)` = unname(r_rk),
                                                                 `|Δ medianas|` = abs(unname(med[2] - med[1])), `|Δ medias|` = abs(unname(m[2] - m[1])), check.names = FALSE)),
         notes = paste0(if (exacto) paste0("Valor p exacto (exact = TRUE)", if (hay_emp) paste0("; aviso de R: ", paste(unique(wmsg), collapse = " / "), ".") else paste0(": ", W$method, "."))
                        else "Aproximación normal con corrección de continuidad (exact = FALSE, correct = TRUE), como en el libro.",
                        " W del grupo = suma de rangos − n(n+1)/2; el menor se usa como estadístico (cap. 29.3.2). r = |Z|/√N con Z = qnorm(p/2) (cap. 30.6)."))
  },
  posthoc = NULL,
  plots = list(boxplot = list(label = "Boxplot por grupo (mediana)", fn = function(d, v, res, cos) plot_box(d, v$response, v$group, cos, center = "median"))),
  opts_ui = list(metodo = list(label = "Cálculo del valor p (cap. 29.3.3)", type = "select",
                               choices = c("Aproximación normal con corrección de continuidad (libro)" = "normal", "Exacto (muestras pequeñas sin empates)" = "exacto"), value = "normal")),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- .find(checks, "normalidad"); lev <- .find(checks, "homogeneidad")
    hi1 <- s$Rm1 >= s$Rm2
    A <- if (hi1) list(g = s$g1, m = s$m1, s = s$s1, md = s$med1) else list(g = s$g2, m = s$m2, s = s$s2, md = s$med2)
    B <- if (hi1) list(g = s$g2, m = s$m2, s = s$s2, md = s$med2) else list(g = s$g1, m = s$m1, s = s$s1, md = s$med1)
    paste0("La normalidad de los datos fue evaluada con ", .md_norm_how(sw), " (", .sw_txt(sw), ")",
           .md_norm_verdict(sw, ", sin evidencia clara de desviación de la normalidad", ", lo que evidencia una distribución no normal"),
           ". La homogeneidad de varianzas fue evaluada con la prueba de Levene (", .lev_txt(lev), "). ",
           "La ", CTX("variable de respuesta"), " fue ", if (s$Rm1 == s$Rm2) "similar en " else "mayor en ", CTX(A$g), " (media = ", fmt_num(A$m), " ± ", fmt_num(A$s), "; mediana = ", fmt_num(A$md), ")",
           if (s$Rm1 == s$Rm2) " y en " else " que en ", CTX(B$g), " (media = ", fmt_num(B$m), " ± ", fmt_num(B$s), "; mediana = ", fmt_num(B$md), "); ",
           "la diferencia entre los grupos (Δ de medias = ", fmt_num(s$dmean), ") ",
           if (s$p < ALPHA) "resultó estadísticamente significativa" else "no resultó estadísticamente significativa",
           " (Wilcoxon: W = ", fmt_num(s$W, 1), ", p ", fmt_p(s$p), "), con un tamaño de efecto ", efecto_r(s$r), " (r = ", fmt_num(s$r), ") (Field et al., 2012; Wilcoxon, 1945).")
  }
))

# =============================================================================
# 3. t apareada (dos mediciones del mismo sujeto)  (alt: Wilcoxon rangos con signo) - cap. 16
# =============================================================================
.paired_diff <- function(d, v) { a <- num(d, v$measures[1]); b <- num(d, v$measures[2]); ok <- !is.na(a) & !is.na(b); list(a = a[ok], b = b[ok], dif = a[ok] - b[ok]) }
.md_paired_desc <- function(pd, v, cols) { y <- c(pd$a, pd$b); g <- factor(rep(v$measures[1:2], each = length(pd$a)), levels = v$measures[1:2]); .md_desc(y, g, cols, "Condición") }

ado_register(list(
  id = "t_paired", family = "medias", order = 3,
  label = "Prueba t apareada (antes/después)",
  description = "Compara dos mediciones del mismo sujeto. Alternativa: Wilcoxon de rangos con signo.",
  chapter = "Prueba t dependiente (apareada) en R",
  roles = list(measures = role("Dos columnas de medición (antes y después)", "numeric", multiple = TRUE, min = 2),
               subject  = role("Identificador del sujeto (opcional)", "id", required = FALSE)),
  example = list(file = "Carbuse.txt",
                 prep = function(d) { w <- reshape(d, direction = "wide", idvar = "ID", timevar = "Groups"); names(w) <- sub("^Ra\\.", "Ra_", names(w)); w },
                 vars = list(measures = c("Ra_0_uses", "Ra_1_use"), subject = "ID"),
                 note = "Rugosidad (Ra) de fresas antes (0 usos) y después (1 uso)."),
  alternative = "wilcox_sr",
  descriptives = function(d, v) { long <- wide_to_long(d, v$measures[1:2], v$subject); desc_table1(long, "Response", "Times") },
  assumptions = function(d, v) {
    pd <- .paired_diff(d, v)
    checks <- list(check_normality(pd$dif, "las diferencias intra-sujeto"))
    list(checks = checks, decision = decide_path(checks, TRUE))
  },
  run = function(d, v, opts) {
    pd <- .paired_diff(d, v)
    tt <- t.test(pd$a, pd$b, paired = TRUE, alternative = "two.sided")
    t <- unname(tt$statistic); df <- unname(tt$parameter); ci <- tt$conf.int
    list(main = tt,
         stats = list(t = t, df = df, p = tt$p.value, r = r_from_t(t, df), delta = unname(tt$estimate), lo = ci[1], hi = ci[2],
                      m1 = mean(pd$a), s1 = sd(pd$a), m2 = mean(pd$b), s2 = sd(pd$b), n = length(pd$a), c1 = v$measures[1], c2 = v$measures[2]),
         tables = list("Estadística descriptiva por condición" = .md_paired_desc(pd, v, c("n", "Media", "DE", "Mediana", "IQR")),
                       "Resultado de t.test(paired = TRUE)" = data.frame(t = t, gl = df, p = tt$p.value, Diferencia = paste(v$measures[1], "−", v$measures[2]),
                                                                         `Δ medias` = unname(tt$estimate), `IC 95% inf` = ci[1], `IC 95% sup` = ci[2], r = r_from_t(t, df), check.names = FALSE)),
         notes = "En diseño apareado no se evalúa homogeneidad de varianzas; la normalidad se evalúa sobre las diferencias (cap. 16.2.4). r = √(t²/(t²+gl)).")
  },
  posthoc = NULL,
  plots = list(
    pares = list(label = "Líneas por sujeto + boxplot", fn = function(d, v, res, cos) plot_paired(wide_to_long(d, v$measures[1:2], v$subject), cos, "mean")),
    qq    = list(label = "Q-Q de las diferencias", fn = function(d, v, res, cos) plot_qq(.paired_diff(d, v)$dif, "Q-Q: diferencias", cos))
  ),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- .find(checks, "normalidad")
    paste0("La normalidad de las diferencias fue evaluada mediante ", .md_norm_how(sw), " (", .sw_txt(sw), ") y evaluación gráfica (Q-Q plot)", .md_norm_verdict(sw),
           ". La diferencia en ", CTX("la variable medida"), " entre la condición ", CTX(s$c1), " (Media = ", fmt_num(s$m1, 3), ", DE = ", fmt_num(s$s1, 3),
           ") y ", CTX(s$c2), " (Media = ", fmt_num(s$m2, 3), ", DE = ", fmt_num(s$s2, 3), ") ",
           if (s$p < ALPHA) "fue estadísticamente significativa" else "no fue estadísticamente significativa",
           " (t = ", fmt_num(s$t), ", df = ", s$df, ", p ", fmt_p(s$p), "), con una diferencia media de Δ = ", fmt_num(s$delta),
           " (IC 95 %: ", fmt_num(s$lo), ", ", fmt_num(s$hi), ") y un tamaño de efecto ", efecto_r(s$r), " (r = ", fmt_num(s$r), ") (Cohen, 1988; Tomczak y Tomczak, 2014).")
  }
))

# =============================================================================
# 4. Wilcoxon rangos con signo - caps. 31–32
# =============================================================================
ado_register(list(
  id = "wilcox_sr", family = "medias", order = 4,
  label = "Wilcoxon de rangos con signo (apareado)",
  description = "Alternativa no paramétrica a la prueba t apareada.",
  chapter = "Wilcoxon rangos con signo en R",
  roles = list(measures = role("Dos columnas de medición (antes y después)", "numeric", multiple = TRUE, min = 2),
               subject  = role("Identificador del sujeto (opcional)", "id", required = FALSE)),
  example = list(file = "DXKTw.txt", prep = function(d) subset(d, Grupo == 2), vars = list(measures = c("X0", "X2"), subject = "Px"),
                 note = "Dolor (NRS) inmediato (X0) y a las 2 horas (X2), grupo control."),
  alternative = NULL,
  descriptives = function(d, v) { long <- wide_to_long(d, v$measures[1:2], v$subject); desc_table1(long, "Response", "Times") },
  assumptions = function(d, v) {
    pd <- .paired_diff(d, v); checks <- list(check_normality(pd$dif, "las diferencias intra-sujeto"))
    list(checks = checks, decision = list(parametric = FALSE, all_met = FALSE, violated = character(0), reason = "Prueba no paramétrica."))
  },
  run = function(d, v, opts) {
    pd <- .paired_diff(d, v)
    W <- wilcox.test(pd$a, pd$b, paired = TRUE, alternative = "two.sided", exact = FALSE, correct = TRUE)
    N <- length(pd$a); r <- r_from_z(W$p.value, N)
    # Ranking con signo manual (cap. 31.3.2–31.3.3): se excluyen las diferencias = 0
    dz <- pd$dif[pd$dif != 0]; n0 <- length(dz); rk <- rank(abs(dz)) * sign(dz)
    Wp <- sum(rk[rk > 0]); Wn <- abs(sum(rk[rk < 0])); Tst <- min(Wp, Wn)
    muT <- n0 * (n0 + 1) / 4; seT <- sqrt(n0 * (n0 + 1) * (2 * n0 + 1) / 24); zT <- (Tst - muT) / seT
    list(main = W,
         stats = list(V = unname(W$statistic), p = W$p.value, r = r, N = N, c1 = v$measures[1], c2 = v$measures[2], Wp = Wp, Wn = Wn,
                      med1 = median(pd$a), med2 = median(pd$b), m1 = mean(pd$a), s1 = sd(pd$a), m2 = mean(pd$b), s2 = sd(pd$b),
                      dmed = abs(median(pd$a) - median(pd$b)), dmean = abs(mean(pd$a) - mean(pd$b))),
         tables = list("Estadística descriptiva por condición" = .md_paired_desc(pd, v, c("n", "Mediana", "IQR", "Media", "DE")),
                       "Rangos con signo" = data.frame(`Pares con diferencia ≠ 0` = n0, `W+ (suma de rangos positivos)` = Wp, `W− (suma de rangos negativos)` = Wn,
                                                       `T (el menor)` = Tst, `Media de T` = muT, `EE de T` = seT, z = zT, check.names = FALSE),
                       "Resultado de wilcox.test(paired = TRUE)" = data.frame(V = unname(W$statistic), p = W$p.value, `r (a partir de Z)` = r,
                                                                              `|Δ medianas|` = abs(median(pd$a) - median(pd$b)), `|Δ medias|` = abs(mean(pd$a) - mean(pd$b)), check.names = FALSE)),
         notes = paste0("Diferencias = ", v$measures[1], " − ", v$measures[2], ". Aproximación normal con corrección de continuidad (exact = FALSE, correct = TRUE). ",
                        "z de la tabla de rangos: (T − n(n+1)/4)/√(n(n+1)(2n+1)/24), sin corrección por empates (cap. 31.3.3); el valor p proviene de wilcox.test(). r = |Z|/√N con N = número de pares."))
  },
  posthoc = NULL,
  plots = list(pares = list(label = "Líneas por sujeto + boxplot (mediana)", fn = function(d, v, res, cos) plot_paired(wide_to_long(d, v$measures[1:2], v$subject), cos, "median"))),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- .find(checks, "normalidad")
    # Dirección: W− > W+ ⇒ las diferencias (c1 − c2) son mayormente negativas ⇒ c2 mayor
    hi2 <- s$Wn > s$Wp
    A <- if (hi2) list(c = s$c2, md = s$med2, m = s$m2, s = s$s2) else list(c = s$c1, md = s$med1, m = s$m1, s = s$s1)
    B <- if (hi2) list(c = s$c1, md = s$med1, m = s$m1, s = s$s1) else list(c = s$c2, md = s$med2, m = s$m2, s = s$s2)
    paste0("La normalidad de la diferencia entre ", CTX(s$c1), " y ", CTX(s$c2), " fue evaluada con ", .md_norm_how(sw), " (", .sw_txt(sw), ")",
           .md_norm_verdict(sw, ", sin evidencia clara de desviación de la normalidad", ", lo que evidencia una distribución no normal"),
           ". La diferencia en ", CTX("la variable medida"), " entre ", CTX(s$c1), " y ", CTX(s$c2), " fue evaluada con la prueba de Wilcoxon de rangos con signo. ",
           CTX("La variable"), if (s$Wn == s$Wp) " fue similar en " else " fue mayor en ", CTX(A$c), " (mediana = ", fmt_num(A$md), "; media = ", fmt_num(A$m), " ± ", fmt_num(A$s), ")",
           if (s$Wn == s$Wp) " y en " else " que en ", CTX(B$c), " (mediana = ", fmt_num(B$md), "; media = ", fmt_num(B$m), " ± ", fmt_num(B$s), "), ",
           "con una diferencia de medianas de ", fmt_num(s$dmed), " y una diferencia de medias de ", fmt_num(s$dmean), ". Esta diferencia ",
           if (s$p < ALPHA) "fue estadísticamente significativa" else "no fue estadísticamente significativa",
           " (V = ", fmt_num(s$V, 1), ", p ", fmt_p(s$p), "), con un tamaño de efecto ", efecto_r(s$r), " (r = ", fmt_num(s$r), ") (Field et al., 2012; Wilcoxon, 1945).")
  }
))

# =============================================================================
# 5. ANOVA de una vía  (alt: Kruskal-Wallis) - caps. 18–19
# =============================================================================
.md_a1_methods <- c(tukey = "Tukey HSD", dunnett = "Dunnett", bonferroni = "Bonferroni", holm = "Holm", bh = "Benjamini-Hochberg", games_howell = "Games-Howell")

# Contrastes planeados / tendencia (summary.lm de aov con contrastes; caps. 18.5, 19.5, 19.6, 19.8)
.md_a1_contrasts <- function(dd, opts) {
  tipo <- if (is.null(opts$contrastes)) "ninguno" else opts$contrastes
  if (identical(tipo, "ninguno")) return(NULL)
  lv <- levels(dd$g); k <- length(lv); ctrl <- .md_ctrl(dd$g, opts$control); nt <- NULL
  if (tipo == "tendencia") { M <- contr.poly(k); nm <- c("Lineal", "Cuadrática", "Cúbica", if (k > 4) paste0("Grado ", 4:(k - 1)))[seq_len(k - 1)]
    lab <- "Análisis de tendencia (contr.poly)"; nt <- "Contrastes polinomiales: niveles tratados como ordenados y equiespaciados, en el orden del factor (cap. 19.6)."
  } else if (tipo == "personalizado") { M <- .md_parse_contr(opts$coef, lv); nm <- colnames(M); lab <- "Contrastes planeados personalizados"
    nt <- paste0("Coeficientes en el orden de niveles: ", paste(lv, collapse = ", "), ".",
                 if (!isTRUE(attr(M, "orthogonal"))) " ATENCIÓN: los contrastes no son ortogonales (el libro exige ortogonalidad para interpretarlos de forma independiente)." else " Contrastes ortogonales.")
  } else { M <- .md_orth(lv, ctrl); nm <- colnames(M); lab <- paste0("Contrastes planeados ortogonales (control = ", ctrl, ")")
    nt <- "C1 compara el control contra el promedio de los demás grupos; los siguientes comparan los grupos restantes entre sí (cap. 18.5.3)." }
  dc <- dd; contrasts(dc$g) <- M; mc <- aov(y ~ g, data = dc)
  tab <- .md_contr_tab(mc, 1 + seq_along(nm), nm)
  coefs <- data.frame(Grupo = lv, round(unclass(M)[, seq_along(nm), drop = FALSE], 4), check.names = FALSE, row.names = NULL); names(coefs)[-1] <- nm
  list(tipo = tipo, label = lab, table = tab, coefs = coefs, note = nt)
}

# Post-hoc elegido (caps. 18.6–18.7, 19.7)
.md_a1_posthoc <- function(dd, opts) {
  meth <- if (is.null(opts$posthoc) || !opts$posthoc %in% names(.md_a1_methods)) "tukey" else opts$posthoc
  lv <- levels(dd$g); ctrl <- .md_ctrl(dd$g, opts$control); gl <- NULL; nt <- NULL
  if (meth == "tukey") {
    gl <- .md_glht(aov(y ~ g, data = dd), "g", "Tukey"); pr <- .md_pairs_glht(gl, lv); lab <- "Post-hoc de Tukey HSD (multcomp::glht + confint)"
  } else if (meth == "dunnett") {
    d2 <- dd; d2$g <- relevel(d2$g, ref = ctrl); gl <- .md_glht(aov(y ~ g, data = d2), "g", "Dunnett"); pr <- .md_pairs_glht(gl, levels(d2$g))
    lab <- paste0("Post-hoc de Dunnett: cada grupo vs. control (", ctrl, ")")
  } else if (meth == "games_howell") {
    pr <- .md_games_howell(dd$y, dd$g); lab <- "Post-hoc de Games-Howell (no asume varianzas iguales)"
  } else {
    adj <- switch(meth, bonferroni = "bonferroni", holm = "holm", bh = "BH"); pr <- .md_pairs_pw(dd$y, dd$g, adj)
    lab <- paste0("pairwise.t.test (", .md_a1_methods[[meth]], ", DE combinada)")
  }
  pr$sig <- pr$p < ALPHA
  tab <- .md_pairs_tab(pr); if (meth == "games_howell") tab$gl <- pr$gl
  p_lev <- .md_lev(dd$y, dd$g)
  if (!is.na(p_lev)) nt <- if (p_lev <= ALPHA && meth != "games_howell") "Levene indica varianzas heterogéneas: el libro recomienda Games-Howell (cap. 18.7)."
                            else if (p_lev > ALPHA && meth == "games_howell") "Levene indica varianzas homogéneas: el libro sugiere Tukey HSD u Holm en ese escenario (cap. 18.7)." else NULL
  list(method = meth, name = .md_a1_methods[[meth]], label = lab, table = tab, pairs = pr, glht = gl, ctrl = ctrl, note = nt)
}

ado_register(list(
  id = "anova1", family = "medias", order = 5,
  label = "ANOVA de una vía (≥ 3 grupos)",
  description = "Compara la media de tres o más grupos independientes. Alternativa: Kruskal-Wallis.",
  chapter = "ANOVA de una vía en R",
  roles = list(response = role("Variable de respuesta (numérica)", "numeric"),
               group    = role("Variable de grupo (≥ 3 niveles)", "factor")),
  example = list(file = "VickersF.txt", prep = function(d) subset(d, Group == "Multilink"),
                 vars = list(response = "VH", group = "Cure.time"),
                 note = "Dureza Vickers del cemento Multilink en tres protocolos de polimerización."),
  alternative = "kruskal",
  descriptives = function(d, v) desc_table1(d, v$response, v$group),
  assumptions = function(d, v) {
    y <- num(d, v$response); g <- fac(d, v$group)
    checks <- list(check_normality(y, v$response), check_homogeneity(y, g))
    dec <- decide_path(checks, TRUE)
    nm <- isFALSE(checks[[1]]$met); hm <- isFALSE(checks[[2]]$met)
    if (nm && !hm) dec$reason <- paste(dec$reason, "Nota (cap. 18.4 y 33.3.1): el ANOVA es robusto ante violaciones moderadas de la normalidad, pero más sensible a la heterocedasticidad; con varianzas homogéneas puede mantenerse la ruta paramétrica.")
    if (hm && !nm) dec$reason <- paste(dec$reason, "Si se mantiene el ANOVA, el libro recomienda Games-Howell como post-hoc (cap. 18.7).")
    list(checks = checks, decision = dec)
  },
  run = function(d, v, opts) {
    dd <- .md_dd(d, v)
    m <- aov(y ~ g, data = dd); sm <- summary(m)[[1]]
    SS_M <- sm[1, "Sum Sq"]; SS_R <- sm[2, "Sum Sq"]; SS_T <- SS_M + SS_R; df_M <- sm[1, "Df"]; MS_R <- sm[2, "Mean Sq"]
    eta2 <- SS_M / SS_T; omega2 <- (SS_M - df_M * MS_R) / (SS_T + MS_R)
    ct <- .md_a1_contrasts(dd, opts); ph <- .md_a1_posthoc(dd, opts)
    tabs <- list("Estadística descriptiva por grupo" = .md_desc(dd$y, dd$g, c("n", "Media", "DE", "EE")),
                 "Tabla ANOVA" = data.frame(Fuente = c(v$group, "Residuales"), sm, check.names = FALSE, row.names = NULL),
                 "Tamaño de efecto" = data.frame(`η²` = eta2, `ω²` = omega2, `r (√η²)` = sqrt(max(eta2, 0)), check.names = FALSE))
    if (!is.null(ct)) { tabs[[paste0(ct$label, ": coeficientes")]] <- ct$coefs; tabs[[paste0(ct$label, " (summary.lm)")]] <- ct$table }
    list(main = m, ph = ph, ct = ct,
         stats = list(F = sm[1, "F value"], df1 = df_M, df2 = sm[2, "Df"], p = sm[1, "Pr(>F)"], eta2 = eta2, omega2 = omega2, r = sqrt(max(eta2, 0)), k = nlevels(dd$g), gname = v$group),
         tables = tabs,
         notes = paste(c("Modelo aov(respuesta ~ grupo). η² = SS_M/SS_T; ω² = (SS_M − df_M·MS_R)/(SS_T + MS_R) (cap. 19.8).",
                         if (!is.null(ct)) paste0(ct$note, " r por contraste = √(t²/(t²+gl)).")), collapse = " "))
  },
  posthoc = function(res, d, v, opts) {
    ph <- if (!is.null(res$ph)) res$ph else .md_a1_posthoc(.md_dd(d, v), opts)
    list(tables = setNames(list(ph$table), ph$label), glht = ph$glht,
         notes = paste(c("Interpretar las comparaciones por pares solo si el ANOVA global es significativo (prueba ómnibus, cap. 18.2).", ph$note), collapse = " "))
  },
  plots = list(
    boxplot = list(label = "Boxplot por grupo", fn = function(d, v, res, cos) plot_box(d, v$response, v$group, cos)),
    qq      = list(label = "Q-Q de la respuesta", fn = function(d, v, res, cos) plot_qq(num(d, v$response), paste("Q-Q:", v$response), cos)),
    tukey   = list(label = "IC 95 % del post-hoc (Tukey o Dunnett)", fn = function(d, v, res, cos) {
      gl <- res$ph$glht
      if (is.null(gl)) gl <- .md_glht(aov(y ~ g, data = .md_dd(d, v)), "g", "Tukey")
      plot_tukey(gl, cos) })
  ),
  opts_ui = list(
    posthoc = list(label = "Método post-hoc", type = "select",
                   choices = c("Tukey HSD (glht)" = "tukey", "Dunnett (vs. grupo control)" = "dunnett", "Bonferroni (pairwise.t.test)" = "bonferroni",
                               "Holm (pairwise.t.test)" = "holm", "Benjamini-Hochberg (pairwise.t.test)" = "bh", "Games-Howell (varianzas desiguales)" = "games_howell"),
                   value = "tukey"),
    contrastes = list(label = "Contrastes planeados", type = "select",
                      choices = c("Ninguno" = "ninguno", "Ortogonales: control vs. resto, luego entre los demás" = "ortogonal",
                                  "Personalizados (coeficientes abajo)" = "personalizado", "Tendencia polinómica (contr.poly)" = "tendencia"),
                      value = "ninguno"),
    control = list(label = "Grupo control (Dunnett y contrastes ortogonales; vacío = primer nivel)", type = "text", value = ""),
    coef = list(label = "Coeficientes personalizados: un contraste por bloque separado por ';', en el orden de los niveles (p. ej. -2,1,1; 0,-1,1)", type = "text", value = "")
  ),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- .find(checks, "normalidad"); lev <- .find(checks, "homogeneidad"); ph <- res$ph; ct <- res$ct
    ct_txt <- if (!is.null(ct)) {
      it <- sapply(seq_len(nrow(ct$table)), function(i) { r <- ct$table[i, ]
        paste0(if (ct$tipo == "tendencia") paste0("la tendencia ", tolower(r$Contraste)) else paste0("el contraste ", r$Contraste),
               if (r$p < ALPHA) " fue significativ" else " no fue significativ", if (ct$tipo == "tendencia") "a" else "o",
               " (t = ", fmt_num(r$t), ", df = ", r$gl, ", p ", fmt_p(r$p), ", r = ", fmt_num(r$r), ")") })
      it[1] <- paste0(toupper(substr(it[1], 1, 1)), substring(it[1], 2))
      paste0(" ", paste(it, collapse = "; "), ".")
    } else ""
    ph_txt <- if (!is.null(ph)) paste0(" El análisis post-hoc de ", ph$name, if (ph$method == "dunnett") paste0(" (vs. ", ph$ctrl, ")") else "", " ",
                                       .md_pairs_sentence(ph$pairs, .md_ptxt(ph$pairs$p), list_ns = FALSE), ".") else ""
    .md_cite(paste0("La normalidad de ", CTX(v$response), " se evaluó con ", .md_norm_how(sw), " (", .sw_txt(sw), ") y Q-Q plot", .md_norm_verdict(sw),
           ". La homogeneidad de varianzas ", if (isTRUE(lev$met)) "fue confirmada" else if (isFALSE(lev$met)) "no se cumplió" else "fue evaluada", " por la prueba de Levene (", .lev_txt(lev), "). ",
           "La ", CTX("variable de respuesta"), " entre los ", s$k, " grupos de ", CTX(v$group), " ",
           if (s$p < ALPHA) "fue significativamente diferente" else "no mostró diferencias significativas",
           " (F(", s$df1, ", ", s$df2, ") = ", fmt_num(s$F), ", p ", fmt_p(s$p), ", ω² = ", fmt_num(s$omega2), ", efecto ", efecto_eta(s$omega2), ").",
           ct_txt, if (s$p >= ALPHA && nzchar(ph_txt)) " Dado que el ANOVA global no fue significativo, las comparaciones por pares se presentan solo de forma descriptiva." else ph_txt),
           "Cohen, 1988; Tomczak y Tomczak, 2014")
  }
))

# =============================================================================
# 6. Kruskal-Wallis - caps. 33–34
# =============================================================================
# Método post-hoc elegido (valores desconocidos, p. ej. "tukey" heredado del ANOVA → Siegel y Castellan)
.kw_meth <- function(opts) { m <- opts$posthoc; if (is.null(m) || !m %in% c("sc", "wilcox_bonf")) "sc" else m }

# Comparaciones múltiples de Siegel y Castellan (1988, p. 213), igual que pgirmess::kruskalmc()
.siegel_castellan <- function(y, g, alpha = 0.05) {
  g <- droplevels(as.factor(g)); r <- rank(y); N <- length(y); k <- nlevels(g)
  Rm <- tapply(r, g, mean); n <- tapply(r, g, length)
  z_crit <- qnorm(alpha / (k * (k - 1)), lower.tail = FALSE); m <- k * (k - 1) / 2
  pr <- utils::combn(levels(g), 2)
  do.call(rbind, lapply(seq_len(ncol(pr)), function(j) { a <- pr[1, j]; b <- pr[2, j]
    se <- sqrt(N * (N + 1) / 12 * (1 / n[[a]] + 1 / n[[b]])); dif <- abs(Rm[[a]] - Rm[[b]]); z <- dif / se
    data.frame(Comparación = paste(a, "-", b), `Rango medio A` = Rm[[a]], `Rango medio B` = Rm[[b]],
               `Dif. observada` = dif, `Dif. crítica` = z_crit * se, z = z,
               `p ajustada (Bonferroni)` = min(1, 2 * pnorm(z, lower.tail = FALSE) * m),
               Significativa = if (dif > z_crit * se) "Sí" else "No", check.names = FALSE) }))
}

.md_kw_posthoc <- function(dd, opts) {
  meth <- .kw_meth(opts); m <- tapply(dd$y, dd$g, mean); pc <- utils::combn(levels(dd$g), 2)
  pr <- data.frame(g1 = pc[1, ], g2 = pc[2, ], est = unname(m[pc[2, ]] - m[pc[1, ]]))
  if (meth == "sc") {
    sc <- .siegel_castellan(dd$y, dd$g, ALPHA); pr$sig <- sc$Significativa == "Sí"
    ptxt <- ifelse(pr$sig, paste0("p < ", ALPHA), paste0("p > ", ALPHA))
    tab <- sc[, c("Comparación", "Rango medio A", "Rango medio B", "Dif. observada", "Dif. crítica", "Significativa")]
    lab <- "Siegel y Castellan (kruskalmc): diferencia de rangos medios vs. diferencia crítica"
    nt <- paste0("Siegel y Castellan (1988): |R̄i − R̄j| > z[α/(k(k−1))]·√(N(N+1)/12·(1/ni + 1/nj)), α = ", ALPHA, ". Equivale a pgirmess::kruskalmc(probs = ", ALPHA, ").")
    name <- "Siegel y Castellan"
  } else {
    pw <- pairwise.wilcox.test(dd$y, dd$g, p.adjust.method = "bonferroni", exact = FALSE)
    pr$p <- mapply(function(a, b) { x <- pw$p.value; if (b %in% rownames(x) && a %in% colnames(x) && !is.na(x[b, a])) x[b, a] else x[a, b] }, pr$g1, pr$g2)
    pr$sig <- pr$p < ALPHA; ptxt <- .md_ptxt(pr$p)
    tab <- data.frame(Comparación = paste(pr$g2, "-", pr$g1), `Δ medias` = pr$est, `p ajustada (Bonferroni)` = pr$p, check.names = FALSE)
    lab <- "pairwise.wilcox.test (Bonferroni, exact = FALSE)"; nt <- NULL; name <- "Wilcoxon por pares con corrección de Bonferroni"
  }
  list(method = meth, name = name, label = lab, table = tab, pairs = pr, ptxt = ptxt, note = nt)
}

ado_register(list(
  id = "kruskal", family = "medias", order = 6,
  label = "Kruskal-Wallis (≥ 3 grupos)",
  description = "Alternativa no paramétrica al ANOVA de una vía.",
  chapter = "Kruskal-Wallis en R",
  roles = list(response = role("Variable de respuesta (numérica)", "numeric"),
               group    = role("Variable de grupo (≥ 3 niveles)", "factor")),
  example = list(file = "ReMat.txt", vars = list(response = "Response", group = "Groups"),
                 note = "Resistencia mecánica (MPa) de tres materiales."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, v$response, v$group),
  assumptions = function(d, v) {
    y <- num(d, v$response); g <- fac(d, v$group)
    checks <- list(check_normality(y, v$response), check_homogeneity(y, g))
    list(checks = checks, decision = list(parametric = FALSE, all_met = FALSE, violated = character(0), reason = "Prueba no paramétrica."))
  },
  run = function(d, v, opts) {
    dd <- .md_dd(d, v)
    KW <- kruskal.test(y ~ g, data = dd); N <- nrow(dd); K <- nlevels(dd$g)
    eta2 <- (unname(KW$statistic) - K + 1) / (N - K)
    rk <- rank(dd$y); ph <- .md_kw_posthoc(dd, opts)
    list(main = KW, ph = ph,
         stats = list(H = unname(KW$statistic), df = unname(KW$parameter), p = KW$p.value, eta2 = eta2, N = N, K = K, posthoc = ph$method,
                      means = tapply(dd$y, dd$g, mean), sds = tapply(dd$y, dd$g, sd)),
         tables = list("Estadística descriptiva (media, DE, mediana e IQR)" = .md_desc(dd$y, dd$g),
                       "Rangos por grupo" = data.frame(Grupo = levels(dd$g), n = as.numeric(tapply(rk, dd$g, length)), `Suma de rangos (R_k)` = as.numeric(tapply(rk, dd$g, sum)),
                                                       `Rango medio` = as.numeric(tapply(rk, dd$g, mean)), check.names = FALSE),
                       "Resultado de kruskal.test()" = data.frame(H = unname(KW$statistic), gl = unname(KW$parameter), p = KW$p.value, `η²_H` = eta2, check.names = FALSE)),
         notes = "η²_H = (H − k + 1)/(n − k) (cap. 34.7). El libro recomienda reportar media ± DE junto con mediana e IQR (cap. 28.5).")
  },
  posthoc = function(res, d, v, opts) {
    ph <- if (!is.null(res$ph)) res$ph else .md_kw_posthoc(.md_dd(d, v), opts)
    list(tables = setNames(list(ph$table), ph$label), notes = paste(c("Post-hoc solo si Kruskal-Wallis global es significativo (prueba ómnibus, cap. 33.3.2).", ph$note), collapse = " "))
  },
  opts_ui = list(posthoc = list(label = "Método post-hoc", type = "select",
                                choices = c("Siegel y Castellan (kruskalmc)" = "sc", "Wilcoxon por pares (Bonferroni)" = "wilcox_bonf"),
                                value = "sc")),
  plots = list(boxplot = list(label = "Boxplot por grupo (mediana)", fn = function(d, v, res, cos) plot_box(d, v$response, v$group, cos, "median"))),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- .find(checks, "normalidad"); lev <- .find(checks, "homogeneidad"); ph <- res$ph
    viol <- c(if (isFALSE(sw$met)) "normalidad", if (isFALSE(lev$met)) "homogeneidad de varianzas")
    o <- order(s$means, decreasing = TRUE); gl <- names(s$means)[o]
    grp <- sapply(seq_along(o), function(i) paste0(CTX(gl[i]), " (", if (i == 1) "media = " else "", fmt_num(s$means[[gl[i]]]), " ± ", fmt_num(s$sds[[gl[i]]]), ")"))
    .md_cite(paste0("La normalidad de ", CTX(v$response), " fue evaluada con ", .md_norm_how(sw), " (", .sw_txt(sw), ") y la homogeneidad de varianzas con la prueba de Levene (", .lev_txt(lev), "); ",
           if (length(viol) == 2) "ambos supuestos paramétricos se incumplen, por lo que " else if (length(viol) == 1) paste0("el supuesto de ", viol, " se incumple, por lo que ") else "",
           "la comparación entre los ", s$K, " grupos se realizó con la prueba de Kruskal-Wallis. La diferencia de ", CTX("la variable de respuesta"), " ",
           if (s$p < ALPHA) "fue significativa" else "no fue significativa", " entre los grupos ", .md_join(grp),
           " (H[", s$df, "] = ", fmt_num(s$H), ", p ", fmt_p(s$p), ", η² = ", fmt_num(s$eta2), ", efecto ", efecto_eta(s$eta2), "). ",
           if (s$p < ALPHA) paste0("El análisis post-hoc (", ph$name, ") ", .md_pairs_sentence(ph$pairs, ph$ptxt, list_ns = TRUE), ". ")
           else "Dado que la prueba global no fue significativa, no se interpretan las comparaciones post-hoc. "),
           paste0("Field et al., 2012; Kruskal y Wallis, 1952", if (ph$method == "sc") "; Siegel y Castellan, 1988" else ""))
  }
))

# =============================================================================
# 9. ANOVA factorial (dos vías) - caps. 23, 24, 26 - sin alternativa NP en el libro
# =============================================================================
.md_fact_dd <- function(d, v) { dd <- data.frame(y = num(d, v$response), A = fac(d, v$factorA), B = fac(d, v$factorB)); dd <- dd[complete.cases(dd), ]
  dd$A <- droplevels(dd$A); dd$B <- droplevels(dd$B); dd }
.md_term_name <- function(rn, vars) sapply(strsplit(rn, ":", fixed = TRUE), function(p) paste(substring(p, 2), collapse = " × "))

# Efectos simples (cap. 24.5): factor combinado A×B con contrastes ortogonales:
# contrastes de B (promediando A) + contrastes de A dentro de cada nivel de B
.md_fact_simple <- function(dd, MA, MB) {
  la <- levels(dd$A); lb <- levels(dd$B)
  cells <- as.vector(outer(la, lb, paste, sep = " · "))
  dd$celda <- factor(paste(dd$A, dd$B, sep = " · "), levels = cells)
  if (any(table(dd$celda) == 0)) stop("Hay celdas A×B vacías: no es posible el análisis de efectos simples.")
  ca <- rep(la, times = length(lb)); cb <- rep(lb, each = length(la))
  M <- cbind(MB[cb, , drop = FALSE]); nm <- colnames(MB)
  for (l in lb) { W <- MA[ca, , drop = FALSE] * (cb == l); M <- cbind(M, W); nm <- c(nm, paste0(colnames(MA), " en ", l)) }
  dimnames(M) <- list(cells, nm)
  dc <- dd; contrasts(dc$celda) <- M; mS <- aov(y ~ celda, data = dc)
  tab <- .md_contr_tab(mS, 1 + seq_along(nm), nm, with_r = FALSE)
  tab$Tipo <- c(rep(paste("Contraste de", "B"), ncol(MB)), rep("Efecto simple de A dentro de B", length(lb) * ncol(MA)))
  tab$nivelB <- c(rep(NA, ncol(MB)), rep(lb, each = ncol(MA)))
  list(table = tab, dd = dd)
}

ado_register(list(
  id = "factorial", family = "medias", order = 9,
  label = "ANOVA factorial (dos factores)",
  description = "Efectos principales e interacción de dos factores con contrastes ortogonales y suma de cuadrados tipo III (car::Anova).",
  chapter = "ANOVA de dos vías en R",
  roles = list(response = role("Variable de respuesta (numérica)", "numeric"),
               factorA  = role("Factor A", "factor"), factorB = role("Factor B", "factor")),
  example = list(file = "VickersF.txt", vars = list(response = "VH", factorA = "Group", factorB = "Cure.time"),
                 note = "Dureza Vickers según cemento (Group) y protocolo de polimerización (Cure.time)."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, v$response, c(v$factorA, v$factorB)),
  assumptions = function(d, v) {
    y <- num(d, v$response); A <- fac(d, v$factorA); B <- fac(d, v$factorB)
    checks <- list(check_normality(y, v$response), check_homogeneity(y, interaction(B, A), "las celdas A×B"))
    list(checks = checks, decision = decide_path(checks, FALSE))
  },
  run = function(d, v, opts) {
    dd <- .md_fact_dd(d, v)
    cA <- .md_ctrl(dd$A, opts$ctrlA); cB <- .md_ctrl(dd$B, opts$ctrlB)
    MA <- .md_orth(levels(dd$A), cA); MB <- .md_orth(levels(dd$B), cB)
    contrasts(dd$A) <- MA; contrasts(dd$B) <- MB
    m <- aov(y ~ A * B, data = dd); a3 <- car::Anova(m, type = "III")
    a3 <- a3[rownames(a3) != "(Intercept)", ]
    # ω² como en el libro (cap. 24.6): SS de summary(aov) (secuencial, tipo I) y SS_T = suma total
    sm <- summary(m)[[1]]; rownames(sm) <- trimws(rownames(sm))
    MS_R <- sm["Residuals", "Mean Sq"]; SS_T <- sum(sm[, "Sum Sq"])
    eff <- c("A", "B", "A:B"); omega <- setNames((sm[eff, "Sum Sq"] - sm[eff, "Df"] * MS_R) / (SS_T + MS_R), eff)
    fn <- c(v$factorA, v$factorB, paste(v$factorA, "×", v$factorB))
    tab <- data.frame(Fuente = c(fn, "Residuales"), a3[c(eff, "Residuals"), ], check.names = FALSE, row.names = NULL)
    cf <- summary.lm(m)$coefficients; rn <- rownames(cf)[-1]
    ctab <- .md_contr_tab(m, 1 + seq_along(rn), .md_term_name(rn), with_r = FALSE)
    celdas <- .md_desc(dd$y, interaction(dd$A, dd$B, sep = " · ", lex.order = TRUE), c("n", "Media", "DE"), paste(v$factorA, "·", v$factorB))
    meth <- if (is.null(opts$posthoc) || !opts$posthoc %in% c("ambos", "contrastes", "tukey")) "ambos" else opts$posthoc
    se <- if (meth %in% c("ambos", "contrastes")) tryCatch(.md_fact_simple(dd, MA, MB), error = function(e) NULL) else NULL
    tk <- if (meth %in% c("ambos", "tukey")) {
      dt <- dd; dt$celda <- factor(paste(dt$A, dt$B, sep = " · ")); tkh <- TukeyHSD(aov(y ~ celda, data = dt))$celda
      lv <- levels(dt$celda); prs <- t(sapply(rownames(tkh), function(r) { for (a in lv) for (b in lv) if (a != b && identical(r, paste0(b, "-", a))) return(c(a, b)); c(NA, NA) }))
      list(table = data.frame(Comparación = rownames(tkh), Diferencia = tkh[, "diff"], `IC 95% inf` = tkh[, "lwr"], `IC 95% sup` = tkh[, "upr"], `p ajustada` = tkh[, "p adj"], check.names = FALSE, row.names = NULL),
           pairs = data.frame(g1 = prs[, 1], g2 = prs[, 2], est = tkh[, "diff"], p = tkh[, "p adj"], row.names = NULL))
    } else NULL
    list(main = m, simple = se, tukey = tk,
         stats = list(FA = a3["A", "F value"], pA = a3["A", "Pr(>F)"], dfA = a3["A", "Df"],
                      FB = a3["B", "F value"], pB = a3["B", "Pr(>F)"], dfB = a3["B", "Df"],
                      FAB = a3["A:B", "F value"], pAB = a3["A:B", "Pr(>F)"], dfAB = a3["A:B", "Df"], dfR = a3["Residuals", "Df"],
                      wA = omega[["A"]], wB = omega[["B"]], wAB = omega[["A:B"]], ctrlA = cA, ctrlB = cB, lA = levels(dd$A), lB = levels(dd$B),
                      balanced = length(unique(as.vector(table(dd$A, dd$B)))) == 1),
         tables = list("Estadística descriptiva por celda" = celdas,
                       "ANOVA tipo III (car::Anova) con contrastes ortogonales" = tab,
                       "Contrastes planeados (summary.lm)" = ctab,
                       "Tamaño de efecto ω²" = data.frame(Fuente = fn, SS = sm[eff, "Sum Sq"], gl = sm[eff, "Df"], `ω²` = as.numeric(omega), check.names = FALSE, row.names = NULL)),
         notes = paste0("Contrastes ortogonales por factor (cap. 24.3): control vs. resto y luego entre los demás niveles (control ", v$factorA, " = ", cA, "; control ", v$factorB, " = ", cB,
                        "). ω² = (SS − gl·MS_R)/(SS_T + MS_R) con las SS de summary(aov) (secuenciales), como en el libro (cap. 24.6, donde se denomina ω² parcial)",
                        if (length(unique(as.vector(table(dd$A, dd$B)))) == 1) "; en este diseño balanceado coinciden con las SS tipo III." else "; ATENCIÓN: el diseño está desbalanceado y las SS secuenciales dependen del orden de los factores (cap. 26)."))
  },
  posthoc = function(res, d, v, opts) {
    tabs <- list()
    if (!is.null(res$simple)) { t <- res$simple$table; t$Tipo <- sub("B", v$factorB, sub("A", v$factorA, t$Tipo)); tabs[["Efectos simples: contrastes ortogonales sobre el factor combinado A×B (summary.lm)"]] <- t[, setdiff(names(t), "nivelB")] }
    if (!is.null(res$tukey)) tabs[["Efectos simples: Tukey HSD sobre el factor combinado A×B"]] <- res$tukey$table
    list(tables = tabs, notes = "Con muchas celdas la interpretación de los contrastes se vuelve confusa; el libro recurre entonces al post-hoc de Tukey sobre el factor combinado (cap. 24.5).")
  },
  plots = list(
    interaccion = list(label = "Gráfico de interacción (medias ± IC 95 % bootstrap)", fn = function(d, v, res, cos) plot_interaction(d, v$response, v$factorB, v$factorA, cos)),
    facetas = list(label = "Boxplot por factor B, en paneles por factor A", fn = function(d, v, res, cos) {
      cos <- cos_merge(cos); dd <- .md_fact_dd(d, v)
      p <- ggplot(dd, aes(B, y, fill = B)) + geom_boxplot(alpha = 0.7, outlier.shape = if (cos$jitter) NA else 16) +
        facet_wrap(~ A) + scale_fill_manual(values = .fillN(cos, nlevels(dd$B)))
      if (cos$jitter) p <- p + geom_jitter(width = 0.1, alpha = 0.4, size = 1.5)
      .finish(p, cos, v$factorB, v$response) }),
    qq   = list(label = "Q-Q de la respuesta", fn = function(d, v, res, cos) plot_qq(num(d, v$response), paste("Q-Q:", v$response), cos))
  ),
  opts_ui = list(
    posthoc = list(label = "Análisis de efectos simples", type = "select",
                   choices = c("Contrastes ortogonales y Tukey sobre celdas" = "ambos", "Solo contrastes ortogonales (summary.lm)" = "contrastes", "Solo Tukey HSD sobre celdas" = "tukey"), value = "ambos"),
    ctrlA = list(label = "Nivel control del factor A (vacío = primer nivel)", type = "text", value = ""),
    ctrlB = list(label = "Nivel control del factor B (vacío = primer nivel)", type = "text", value = "")
  ),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- .find(checks, "normalidad"); lev <- .find(checks, "homogeneidad")
    se_txt <- ""
    if (s$pAB < ALPHA && length(s$lA) == 2) {
      # El reporte del libro (cap. 24.7) se apoya en el Tukey sobre el factor combinado; si no se calculó, en los contrastes
      if (!is.null(res$tukey)) { pr <- res$tukey$pairs; spl <- function(x) strsplit(x, " · ", fixed = TRUE)[[1]]
        keep <- mapply(function(a, b) spl(a)[2] == spl(b)[2], pr$g1, pr$g2); pr <- pr[keep, ]
        sigB <- sapply(strsplit(pr$g1[pr$p < ALPHA], " · ", fixed = TRUE), `[`, 2); how <- "Los análisis de efectos simples (Tukey HSD sobre las celdas)" }
      else { t <- res$simple$table; t <- t[!is.na(t$nivelB), ]; sigB <- t$nivelB[t$p < ALPHA]; how <- "Los análisis de efectos simples (contrastes ortogonales)" }
      se_txt <- paste0(" ", how, " revelaron que la diferencia entre ", CTX(s$lA[1]), " y ", CTX(s$lA[2]),
                       if (!length(sigB)) " no fue significativa en ningún nivel de " else paste0(" fue significativa", if (length(sigB) < length(s$lB)) " únicamente" else "",
                                                                                                 if (length(sigB) > 1) " en los niveles " else " en el nivel ", .md_join(sigB), " de "),
                       CTX(v$factorB), ".")
    } else if (s$pAB < ALPHA) se_txt <- paste0(" Los análisis de efectos simples ", CTX("revelaron … (ver tablas de efectos simples)"), ".")
    paste0("La normalidad de ", CTX(v$response), if (isTRUE(sw$met)) " fue confirmada" else if (isFALSE(sw$met)) " no se cumplió" else " no pudo determinarse", " (evaluada con ", .md_norm_how(sw), ": ", .sw_txt(sw),
           ") y la homogeneidad de varianzas fue evaluada con la prueba de Levene (", .lev_txt(lev), "). ",
           "El ANOVA factorial reveló ", if (s$pA < ALPHA) "un efecto principal significativo" else "ausencia de efecto principal significativo", " de ", CTX(v$factorA),
           " (F[", s$dfA, ",", s$dfR, "] = ", fmt_num(s$FA, 1), ", p ", fmt_p(s$pA), ", ω² = ", fmt_num(s$wA), ") y ",
           if (s$pB < ALPHA) "un efecto principal significativo" else "ausencia de efecto principal significativo", " de ", CTX(v$factorB),
           " (F[", s$dfB, ",", s$dfR, "] = ", fmt_num(s$FB, 1), ", p ", fmt_p(s$pB), ", ω² = ", fmt_num(s$wB), "). ",
           "La interacción ", v$factorA, " × ", v$factorB, if (s$pAB < ALPHA) " fue significativa" else " no fue significativa",
           " (F[", s$dfAB, ",", s$dfR, "] = ", fmt_num(s$FAB, 1), ", p ", fmt_p(s$pAB), ", ω² = ", fmt_num(s$wAB), ")",
           if (s$pAB < ALPHA) paste0(", indicando que ", CTX("el efecto de un factor depende del nivel del otro (describir)"), ".") else ".",
           se_txt) |> .md_cite("Cohen, 1988; Field et al., 2012")
  }
))

# =============================================================================
# 10. ANCOVA - caps. 25–27 - sin alternativa NP en el libro
# =============================================================================
.md_anc_dd <- function(d, v) { dd <- data.frame(y = num(d, v$response), g = fac(d, v$group), x = num(d, v$covariate)); dd <- dd[complete.cases(dd), ]; dd$g <- droplevels(dd$g); dd }

ado_register(list(
  id = "ancova", family = "medias", order = 10,
  label = "ANCOVA (grupo + covariable)",
  description = "Compara grupos ajustando por una covariable. Verifica independencia de la covariable y homogeneidad de pendientes.",
  chapter = "ANCOVA en R",
  roles = list(response  = role("Variable de respuesta (numérica)", "numeric"),
               group     = role("Variable de grupo", "factor"),
               covariate = role("Covariable (numérica)", "numeric")),
  example = list(file = "datBF.txt", prep = function(d) { d$Group <- factor(d$Group, levels = c("G2", "G1", "G3")); d },
                 vars = list(response = "Final", group = "Group", covariate = "Basal"),
                 note = "Fluorescencia final por grupo, ajustando por la basal (control = G2)."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, c(v$covariate, v$response), v$group),
  assumptions = function(d, v) {
    dd <- .md_anc_dd(d, v)
    ind <- summary(aov(x ~ g, data = dd))[[1]]; p_ind <- ind[1, "Pr(>F)"]
    rs <- car::Anova(aov(y ~ g + x + g:x, data = dd), type = "III"); p_slope <- rs["g:x", "Pr(>F)"]
    checks <- list(
      check_normality(dd$y, v$response), check_normality(dd$x, v$covariate),
      check_homogeneity(dd$y, dd$g, paste(v$response, "entre grupos")), check_homogeneity(dd$x, dd$g, paste(v$covariate, "(covariable) entre grupos")),
      check_obj("independencia", "Independencia de la covariable", p_ind > ALPHA,
                paste0("ANOVA de ", v$covariate, " ~ ", v$group, ": F = ", fmt_num(ind[1, "F value"], 3), ", p ", fmt_p(p_ind), if (p_ind > ALPHA) ": la covariable no difiere entre grupos (independiente del factor)." else ": la covariable difiere entre grupos: condición necesaria del ANCOVA comprometida."), F = ind[1, "F value"], p = p_ind),
      check_obj("pendientes", "Homogeneidad de pendientes", p_slope > ALPHA,
                paste0("Interacción ", v$group, " × ", v$covariate, " (Anova tipo III): F = ", fmt_num(rs["g:x", "F value"], 3), ", p ", fmt_p(p_slope), if (p_slope > ALPHA) ": el supuesto de homogeneidad de pendientes se cumple." else ": las pendientes difieren entre grupos; el ANCOVA clásico no es apropiado."), F = rs["g:x", "F value"], p = p_slope))
    dec <- decide_path(checks[5:6], FALSE)
    if (isFALSE(checks[[2]]$met) && !isFALSE(checks[[4]]$met)) dec$reason <- paste(dec$reason, "La covariable no es normal, pero sus varianzas son homogéneas; el estadístico F es robusto ante desviaciones de la normalidad (cap. 27.5).")
    list(checks = checks, decision = dec)
  },
  run = function(d, v, opts) {
    dd <- .md_anc_dd(d, v); ctrl <- .md_ctrl(dd$g, opts$control); M <- .md_orth(levels(dd$g), ctrl)
    contrasts(dd$g) <- M
    m <- aov(y ~ g + x, data = dd); a3 <- car::Anova(m, type = "III"); a3 <- a3[rownames(a3) != "(Intercept)", ]
    SS_R <- summary(m)[[1]]["Residuals", "Sum Sq"]
    eta_g <- a3["g", "Sum Sq"] / (a3["g", "Sum Sq"] + SS_R); eta_x <- a3["x", "Sum Sq"] / (a3["x", "Sum Sq"] + SS_R)
    ctab <- .md_contr_tab(m, 1 + seq_len(ncol(M)), colnames(M))
    # ANOVA vs. ANCOVA (cap. 27.7)
    m0 <- aov(y ~ g, data = dd); s0 <- summary(m0)[[1]]; s1 <- summary(m)[[1]]; ssr0 <- s0[2, "Sum Sq"]; ssr1 <- s1["Residuals", "Sum Sq"]
    comp <- data.frame(Modelo = c(paste(v$response, "~", v$group, "(ANOVA)"), paste(v$response, "~", v$group, "+", v$covariate, "(ANCOVA)")),
                       `SS grupo` = c(s0[1, "Sum Sq"], s1[1, "Sum Sq"]), `F grupo` = c(s0[1, "F value"], s1[1, "F value"]), `p grupo` = c(s0[1, "Pr(>F)"], s1[1, "Pr(>F)"]),
                       `SS residual` = c(ssr0, ssr1), `gl residual` = c(s0[2, "Df"], s1["Residuals", "Df"]), `Reducción SSR (%)` = c(NA, 100 * (1 - ssr1 / ssr0)), check.names = FALSE)
    adj <- tryCatch({ dd2 <- dd; names(dd2) <- c("y", v$group, "x"); m2 <- aov(as.formula(paste("y ~", bq(v$group), "+ x")), data = dd2)
      e <- effects::effect(v$group, m2, se = TRUE); s <- summary(e)
      data.frame(Grupo = names(s$effect), `Media ajustada` = as.numeric(s$effect), `IC inf` = as.numeric(s$lower), `IC sup` = as.numeric(s$upper), check.names = FALSE) }, error = function(e) NULL)
    gl <- .md_glht(m, "g", "Tukey"); pr <- .md_pairs_glht(gl, levels(dd$g)); pr$sig <- pr$p < ALPHA
    d_desc <- cbind(.md_desc(dd$x, dd$g, c("n", "Media", "DE")), .md_desc(dd$y, dd$g, c("Media", "DE"))[, -1])
    names(d_desc) <- c("Grupo", "n", paste(v$covariate, c("Media", "DE")), paste(v$response, c("Media", "DE")))
    list(main = m, glht = gl, pairs = pr, ctab = ctab,
         stats = list(Fg = a3["g", "F value"], dfg = a3["g", "Df"], pg = a3["g", "Pr(>F)"], Fx = a3["x", "F value"], dfx = a3["x", "Df"], px = a3["x", "Pr(>F)"],
                      dfR = a3["Residuals", "Df"], eta_g = eta_g, eta_x = eta_x, red = 100 * (1 - ssr1 / ssr0), ctrl = ctrl),
         tables = c(list("Estadística descriptiva de covariable y respuesta por grupo" = d_desc,
                         "ANCOVA tipo III (car::Anova) con contrastes ortogonales" = data.frame(Fuente = c(v$group, v$covariate, "Residuales"), a3, check.names = FALSE, row.names = NULL),
                         "Contrastes planeados (summary.lm) y r por contraste" = ctab,
                         "η² parcial" = data.frame(Efecto = c(v$group, v$covariate), `η² parcial` = c(eta_g, eta_x), check.names = FALSE),
                         "Comparación ANOVA vs. ANCOVA (reducción de la SS residual)" = comp),
                    if (!is.null(adj)) list("Medias ajustadas (effects)" = adj)),
         notes = paste0("Contrastes ortogonales (cap. 27.6): control = ", ctrl, " vs. resto, luego entre los grupos restantes. η² parcial = SS_efecto (tipo III)/(SS_efecto + SS_R); r = √(t²/(t²+gl)) por contraste (cap. 27.11)."))
  },
  posthoc = function(res, d, v, opts) {
    gl <- res$glht; pr <- res$pairs
    if (is.null(gl)) { dd <- .md_anc_dd(d, v); gl <- .md_glht(aov(y ~ g + x, data = dd), "g", "Tukey"); pr <- .md_pairs_glht(gl, levels(dd$g)) }
    list(tables = list("Tukey sobre medias ajustadas (multcomp::glht + confint)" = .md_pairs_tab(pr)), glht = gl, notes = "Comparaciones de Tukey sobre las medias ajustadas por la covariable (cap. 27.9).")
  },
  plots = list(
    boxplot = list(label = "Boxplot de la respuesta por grupo", fn = function(d, v, res, cos) plot_box(d, v$response, v$group, cos)),
    pendientes = list(label = "Dispersión covariable–respuesta con rectas por grupo (homogeneidad de pendientes)", fn = function(d, v, res, cos) {
      cos <- cos_merge(cos); dd <- .md_anc_dd(d, v)
      p <- ggplot(dd, aes(x, y, colour = g)) + geom_point(alpha = 0.6, size = 2) +
        geom_smooth(method = "lm", se = cos$show_ci, alpha = 0.10, formula = y ~ x) +
        scale_colour_manual(values = .colN(cos, nlevels(dd$g)), name = v$group)
      cos$legend <- TRUE; .finish(p, cos, v$covariate, v$response) }),
    ajustadas = list(label = "Medias ajustadas (effects)", fn = function(d, v, res, cos) {
      dd <- .md_anc_dd(d, v); names(dd)[2] <- v$group; m <- aov(as.formula(paste("y ~", bq(v$group), "+ x")), data = dd)
      plot_adjusted(effects::effect(v$group, m, se = TRUE), cos) }),
    qq  = list(label = "Q-Q de la respuesta", fn = function(d, v, res, cos) plot_qq(num(d, v$response), paste("Q-Q:", v$response), cos)),
    qqx = list(label = "Q-Q de la covariable", fn = function(d, v, res, cos) plot_qq(num(d, v$covariate), paste("Q-Q:", v$covariate), cos))
  ),
  opts_ui = list(control = list(label = "Grupo control para los contrastes ortogonales (vacío = primer nivel)", type = "text", value = "")),
  report = function(res, d, v, checks) {
    s <- res$stats; swy <- checks[[1]]; swx <- checks[[2]]; levy <- checks[[3]]; levx <- checks[[4]]; ind <- checks[[5]]; pen <- checks[[6]]
    lev_ok <- isTRUE(levy$met) && isTRUE(levx$met)
    ct <- res$ctab; ct_txt <- paste(sapply(seq_len(nrow(ct)), function(i) paste0("el contraste ", ct$Contraste[i], if (ct$p[i] < ALPHA) " fue significativo" else " no fue significativo",
                                                                           " (t = ", fmt_num(ct$t[i]), ", p ", fmt_p(ct$p[i]), ", r = ", fmt_num(ct$r[i]), ", efecto ", efecto_r(ct$r[i]), ")")), collapse = "; ")
    paste0("La normalidad fue evaluada con ", .md_norm_how(swy), " para ", CTX(v$response), " (", .sw_txt(swy), ") y ", CTX(v$covariate), " (", .sw_txt(swx), "); ",
           "la homogeneidad de varianzas ", if (lev_ok) "se confirmó" else "se evaluó", " con la prueba de Levene (", v$response, ": ", .lev_txt(levy), "; ", v$covariate, ": ", .lev_txt(levx), "). ",
           "La covariable (", CTX(v$covariate), ") ", if (isTRUE(ind$met)) "no difirió" else "difirió", " entre grupos (ANOVA: F = ", fmt_num(ind$F, 3), ", p ", fmt_p(ind$p), ")",
           if (isTRUE(ind$met)) ", confirmando su independencia del factor experimental" else "", ". El modelo ANCOVA reveló ",
           if (s$pg < ALPHA) "un efecto significativo" else "ausencia de efecto significativo", " del grupo (F[", s$dfg, ",", s$dfR, "] = ", fmt_num(s$Fg), ", p ", fmt_p(s$pg), ", η²parcial = ", fmt_num(s$eta_g), ") y ",
           if (s$px < ALPHA) "un efecto significativo" else "ausencia de efecto significativo", " de la covariable (F[", s$dfx, ",", s$dfR, "] = ", fmt_num(s$Fx), ", p ", fmt_p(s$px), ", η²parcial = ", fmt_num(s$eta_x), "); ",
           "la covariable redujo la suma de cuadrados residual en un ", fmt_num(s$red, 1), " % respecto al ANOVA sin covariable. ",
           "El supuesto de homogeneidad de pendientes ", if (isTRUE(pen$met)) "fue confirmado" else "no se cumplió", " (interacción ", v$group, " × ", v$covariate, ": p ", fmt_p(pen$p), "). ",
           "En los contrastes planeados, ", ct_txt, ". El análisis de Tukey sobre las medias ajustadas ", .md_pairs_sentence(res$pairs, .md_ptxt(res$pairs$p), list_ns = TRUE),
           " (Field et al., 2012; Miller y Chapman, 2001).")
  }
))
