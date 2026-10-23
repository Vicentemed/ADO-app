# =============================================================================
# analyses/medias.R — Familia "Comparación de medias" (caps. t, ANOVA, ANCOVA
# y sus alternativas no paramétricas: Wilcoxon, Kruskal-Wallis, Friedman)
# Fiel a las llamadas y reglas del libro. R base.
# =============================================================================

r_from_t <- function(t, df) sqrt(t^2 / (t^2 + df))
r_from_z <- function(p, N) abs(qnorm(p / 2)) / sqrt(N)

# Frases de supuestos para el reporte (a partir de checks)
.sw_txt  <- function(ch) if (is.null(ch)) "" else paste0("W = ", fmt_num(ch$W, 2), ", p ", fmt_p(ch$p))
.lev_txt <- function(ch) if (is.null(ch)) "" else paste0("F = ", fmt_num(ch$F, 3), ", p ", fmt_p(ch$p))
.find <- function(checks, kind) { for (ch in checks) if (identical(ch$kind, kind)) return(ch); NULL }
.msd <- function(v) paste0("Media = ", fmt_num(mean(v, na.rm = TRUE)), ", DE = ", fmt_num(sd(v, na.rm = TRUE)))
.mmd <- function(v) paste0("mediana = ", fmt_num(median(v, na.rm = TRUE)), "; media = ", fmt_num(mean(v, na.rm = TRUE)), " ± ", fmt_num(sd(v, na.rm = TRUE)))

# =============================================================================
# 1. t de Student independiente  (alt: Wilcoxon suma de rangos)
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
    list(checks = checks, decision = decide_path(checks[1], TRUE))   # decide por normalidad; Welch cubre varianzas
  },
  run = function(d, v, opts) {
    dd <- data.frame(y = num(d, v$response), g = fac(d, v$group)); dd <- dd[complete.cases(dd), ]
    itt <- t.test(y ~ g, data = dd, alternative = "two.sided", var.equal = isTRUE(opts$var_equal))
    t <- unname(itt$statistic); df <- unname(itt$parameter)
    lv <- levels(dd$g); m <- tapply(dd$y, dd$g, mean); s <- tapply(dd$y, dd$g, sd)
    delta <- unname(m[1] - m[2]); ci <- itt$conf.int
    list(main = itt,
         stats = list(t = t, df = df, p = itt$p.value, r = r_from_t(t, df), delta = delta, lo = ci[1], hi = ci[2],
                      g1 = lv[1], g2 = lv[2], m1 = m[1], m2 = m[2], s1 = s[1], s2 = s[2], n = nrow(dd)),
         tables = list("Resultado de t.test()" = data.frame(t = t, gl = df, p = itt$p.value, `Δ medias` = delta, `IC 95% inf` = ci[1], `IC 95% sup` = ci[2], r = r_from_t(t, df), check.names = FALSE)),
         notes = if (isTRUE(opts$var_equal)) "t de Student clásica (var.equal = TRUE)." else "Corrección de Welch (var.equal = FALSE), opción por defecto y más robusta según el libro.")
  },
  posthoc = NULL,
  plots = list(
    boxplot = list(label = "Boxplot por grupo", fn = function(d, v, res, cos) plot_box(d, v$response, v$group, cos)),
    qq      = list(label = "Q-Q de la respuesta", fn = function(d, v, res, cos) plot_qq(num(d, v$response), paste("Q-Q:", v$response), cos)),
    effect  = list(label = "Diferencia de medias (IC 95 %)", fn = function(d, v, res, cos) plot_effect(res$stats$delta, res$stats$lo, res$stats$hi, "Δ de medias", cos))
  ),
  opts_ui = list(var_equal = list(label = "Varianzas iguales (t clásica)", type = "checkbox", value = FALSE)),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- .find(checks, "normalidad"); lev <- .find(checks, "homogeneidad")
    paste0("La normalidad se evaluó mediante la prueba de Shapiro-Wilk (", .sw_txt(sw), ") y la evaluación gráfica (Q-Q plot)",
           if (isTRUE(sw$met)) ", sin evidencia de distribución diferente a la normal" else ", con evidencia de desviación de la normalidad",
           ". La homogeneidad de varianzas se evaluó con la prueba de Levene (", .lev_txt(lev), "). ",
           "La ", CTX("variable de respuesta"), " entre ", CTX(s$g1), " (", .msd(num(d, v$response)[fac(d, v$group) == s$g1]), ") y ",
           CTX(s$g2), " (", .msd(num(d, v$response)[fac(d, v$group) == s$g2]), ") ",
           if (s$p < ALPHA) "mostró una diferencia estadísticamente significativa" else "no mostró una diferencia estadísticamente significativa",
           " (t = ", fmt_num(s$t), ", df = ", fmt_num(s$df, 1), ", p ", fmt_p(s$p), "), con una diferencia de medias de Δ = ", fmt_num(s$delta),
           " (IC 95 %: ", fmt_num(s$lo), ", ", fmt_num(s$hi), "), que representa un tamaño de efecto ", efecto_r(s$r), " (r = ", fmt_num(s$r), ") [Cohen, 1992].")
  }
))

# =============================================================================
# 2. Wilcoxon suma de rangos (Mann-Whitney)
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
    dd <- data.frame(y = num(d, v$response), g = fac(d, v$group)); dd <- dd[complete.cases(dd), ]
    W <- wilcox.test(y ~ g, data = dd, alternative = "two.sided", exact = FALSE, correct = TRUE)
    N <- nrow(dd); r <- r_from_z(W$p.value, N); lv <- levels(dd$g)
    med <- tapply(dd$y, dd$g, median); m <- tapply(dd$y, dd$g, mean)
    list(main = W,
         stats = list(W = unname(W$statistic), p = W$p.value, r = r, N = N, g1 = lv[1], g2 = lv[2],
                      med1 = med[1], med2 = med[2], dmed = unname(med[1] - med[2]), dmean = unname(m[1] - m[2])),
         tables = list("Resultado de wilcox.test()" = data.frame(W = unname(W$statistic), p = W$p.value, r = r, `Δ medianas` = unname(med[1] - med[2]), `Δ medias` = unname(m[1] - m[2]), check.names = FALSE),
                       "Descriptiva por grupo" = desc_numeric(dd$y, dd$g)),
         notes = "Aproximación normal con corrección de continuidad (exact = FALSE, correct = TRUE), como en el libro.")
  },
  posthoc = NULL,
  plots = list(boxplot = list(label = "Boxplot por grupo (mediana)", fn = function(d, v, res, cos) plot_box(d, v$response, v$group, cos, center = "median"))),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- .find(checks, "normalidad"); lev <- .find(checks, "homogeneidad")
    y <- num(d, v$response); g <- fac(d, v$group)
    paste0("La normalidad de los datos fue evaluada con la prueba de Shapiro-Wilk (", .sw_txt(sw), "), ",
           if (isFALSE(sw$met)) "lo que evidencia una distribución no normal" else "sin evidencia clara de desviación",
           ". La homogeneidad de varianzas fue evaluada con la prueba de Levene (", .lev_txt(lev), "). ",
           "La ", CTX("variable de respuesta"), " en ", CTX(s$g1), " (", .mmd(y[g == s$g1]), ") y en ", CTX(s$g2), " (", .mmd(y[g == s$g2]), ") ",
           "presentó una diferencia (Δ de medias = ", fmt_num(s$dmean), "; Δ de medianas = ", fmt_num(s$dmed), ") que ",
           if (s$p < ALPHA) "resultó estadísticamente significativa" else "no resultó estadísticamente significativa",
           " (Wilcoxon: W = ", fmt_num(s$W, 1), ", p ", fmt_p(s$p), "), con un tamaño de efecto ", efecto_r(s$r), " (r = ", fmt_num(s$r), ") [Wilcoxon, 1945; Field, 2012].")
  }
))

# =============================================================================
# 3. t apareada (dos mediciones del mismo sujeto)  (alt: Wilcoxon rangos con signo)
# =============================================================================
.paired_diff <- function(d, v) { a <- num(d, v$measures[1]); b <- num(d, v$measures[2]); ok <- !is.na(a) & !is.na(b); list(a = a[ok], b = b[ok], dif = a[ok] - b[ok]) }

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
         tables = list("Resultado de t.test(paired = TRUE)" = data.frame(t = t, gl = df, p = tt$p.value, `Δ medias` = unname(tt$estimate), `IC 95% inf` = ci[1], `IC 95% sup` = ci[2], r = r_from_t(t, df), check.names = FALSE)),
         notes = "En diseño apareado no se evalúa homogeneidad de varianzas; la normalidad se evalúa sobre las diferencias (libro).")
  },
  posthoc = NULL,
  plots = list(
    pares = list(label = "Líneas por sujeto + boxplot", fn = function(d, v, res, cos) plot_paired(wide_to_long(d, v$measures[1:2], v$subject), cos, "mean")),
    qq    = list(label = "Q-Q de las diferencias", fn = function(d, v, res, cos) plot_qq(.paired_diff(d, v)$dif, "Q-Q: diferencias", cos))
  ),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- .find(checks, "normalidad")
    paste0("La normalidad de las diferencias fue evaluada mediante la prueba de Shapiro-Wilk (", .sw_txt(sw), ") y evaluación gráfica (Q-Q plot)",
           if (isTRUE(sw$met)) ", sin evidencia de distribución diferente a la normal" else ", con evidencia de desviación de la normalidad",
           ". La diferencia en ", CTX("la variable medida"), " entre la condición ", CTX(s$c1), " (Media = ", fmt_num(s$m1), ", DE = ", fmt_num(s$s1),
           ") y ", CTX(s$c2), " (Media = ", fmt_num(s$m2), ", DE = ", fmt_num(s$s2), ") ",
           if (s$p < ALPHA) "fue estadísticamente significativa" else "no fue estadísticamente significativa",
           " (t = ", fmt_num(s$t), ", df = ", s$df, ", p ", fmt_p(s$p), "), con una diferencia media de Δ = ", fmt_num(s$delta),
           " (IC 95 %: ", fmt_num(s$lo), ", ", fmt_num(s$hi), ") y un tamaño de efecto ", efecto_r(s$r), " (r = ", fmt_num(s$r), ") [Cohen, 1988; Tomczak, 2014].")
  }
))

# =============================================================================
# 4. Wilcoxon rangos con signo
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
    list(main = W,
         stats = list(V = unname(W$statistic), p = W$p.value, r = r, N = N, c1 = v$measures[1], c2 = v$measures[2],
                      med1 = median(pd$a), med2 = median(pd$b), m1 = mean(pd$a), s1 = sd(pd$a), m2 = mean(pd$b), s2 = sd(pd$b),
                      dmed = median(pd$a) - median(pd$b), dmean = mean(pd$a) - mean(pd$b)),
         tables = list("Resultado de wilcox.test(paired = TRUE)" = data.frame(V = unname(W$statistic), p = W$p.value, r = r, `Δ medianas` = median(pd$a) - median(pd$b), `Δ medias` = mean(pd$a) - mean(pd$b), check.names = FALSE)),
         notes = "Aproximación normal con corrección de continuidad (exact = FALSE, correct = TRUE).")
  },
  posthoc = NULL,
  plots = list(pares = list(label = "Líneas por sujeto + boxplot (mediana)", fn = function(d, v, res, cos) plot_paired(wide_to_long(d, v$measures[1:2], v$subject), cos, "median"))),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- .find(checks, "normalidad")
    paste0("La normalidad de la diferencia entre ", CTX(s$c1), " y ", CTX(s$c2), " fue evaluada con la prueba de Shapiro-Wilk (", .sw_txt(sw), "), ",
           if (isFALSE(sw$met)) "lo que evidencia una distribución no normal" else "sin evidencia clara de desviación",
           ". La diferencia en ", CTX("la variable medida"), " entre ", CTX(s$c1), " (mediana = ", fmt_num(s$med1), "; media = ", fmt_num(s$m1), " ± ", fmt_num(s$s1), ") y ",
           CTX(s$c2), " (mediana = ", fmt_num(s$med2), "; media = ", fmt_num(s$m2), " ± ", fmt_num(s$s2), ") fue evaluada con la prueba de Wilcoxon de rangos con signo, ",
           "con una diferencia de medianas de ", fmt_num(s$dmed), " y una diferencia de medias de ", fmt_num(s$dmean), ". Esta diferencia ",
           if (s$p < ALPHA) "fue estadísticamente significativa" else "no fue estadísticamente significativa",
           " (V = ", fmt_num(s$V, 1), ", p ", fmt_p(s$p), "), con un tamaño de efecto ", efecto_r(s$r), " (r = ", fmt_num(s$r), ") [Wilcoxon, 1945; Field, 2012].")
  }
))

# =============================================================================
# 5. ANOVA de una vía  (alt: Kruskal-Wallis)
# =============================================================================
.glht_mcp <- function(model, gname, type) {
  lf <- do.call(multcomp::mcp, setNames(list(type), gname))
  multcomp::glht(model, linfct = lf)
}
.glht_table <- function(g) {
  s <- summary(g); ci <- confint(g)$confint
  data.frame(Comparación = rownames(ci), Estimación = ci[, "Estimate"], `IC inf` = ci[, "lwr"], `IC sup` = ci[, "upr"],
             p = as.numeric(s$test$pvalues), check.names = FALSE, row.names = NULL)
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
    list(checks = checks, decision = decide_path(checks, TRUE))
  },
  run = function(d, v, opts) {
    dd <- data.frame(y = num(d, v$response), g = fac(d, v$group)); dd <- dd[complete.cases(dd), ]
    m <- aov(y ~ g, data = dd); sm <- summary(m)[[1]]
    SS_M <- sm[1, "Sum Sq"]; SS_R <- sm[2, "Sum Sq"]; SS_T <- SS_M + SS_R; df_M <- sm[1, "Df"]; MS_R <- sm[2, "Mean Sq"]
    eta2 <- SS_M / SS_T; omega2 <- (SS_M - df_M * MS_R) / (SS_T + MS_R)
    list(main = m,
         stats = list(F = sm[1, "F value"], df1 = df_M, df2 = sm[2, "Df"], p = sm[1, "Pr(>F)"], eta2 = eta2, omega2 = omega2, r = sqrt(max(eta2, 0)), k = nlevels(dd$g), gname = v$group),
         tables = list("Tabla ANOVA" = data.frame(Fuente = trimws(rownames(sm)), sm, check.names = FALSE, row.names = NULL),
                       "Tamaño de efecto" = data.frame(`η²` = eta2, `ω²` = omega2, r = sqrt(max(eta2, 0)), check.names = FALSE)),
         notes = "Modelo aov(respuesta ~ grupo). η² = SS_M/SS_T; ω² = (SS_M − df_M·MS_R)/(SS_T + MS_R).")
  },
  posthoc = function(res, d, v, opts) {
    dd <- data.frame(y = num(d, v$response), g = fac(d, v$group)); dd <- dd[complete.cases(dd), ]
    m <- aov(y ~ g, data = dd); meth <- if (is.null(opts$posthoc)) "tukey" else opts$posthoc
    tabs <- list()
    if (meth %in% c("tukey", "dunnett")) {
      g <- .glht_mcp(m, "g", if (meth == "tukey") "Tukey" else "Dunnett")
      tabs[[if (meth == "tukey") "Post-hoc de Tukey (multcomp::glht)" else "Dunnett vs. primer nivel (control)"]] <- .glht_table(g)
      attr(tabs, "glht") <- g
    } else {
      pw <- pairwise.t.test(dd$y, dd$g, p.adjust.method = if (meth == "bh") "BH" else "bonferroni")
      pm <- pw$p.value; df <- as.data.frame(as.table(pm)); df <- df[!is.na(df$Freq), ]; names(df) <- c("Grupo A", "Grupo B", "p ajustada")
      tabs[[paste0("pairwise.t.test (", if (meth == "bh") "Benjamini-Hochberg" else "Bonferroni", ")")]] <- df
    }
    list(tables = tabs, glht = attr(tabs, "glht"), notes = "Interpretar el post-hoc solo si el ANOVA global es significativo.")
  },
  plots = list(
    boxplot = list(label = "Boxplot por grupo", fn = function(d, v, res, cos) plot_box(d, v$response, v$group, cos)),
    qq      = list(label = "Q-Q de la respuesta", fn = function(d, v, res, cos) plot_qq(num(d, v$response), paste("Q-Q:", v$response), cos)),
    tukey   = list(label = "IC 95 % de Tukey", fn = function(d, v, res, cos) {
      dd <- data.frame(y = num(d, v$response), g = fac(d, v$group)); dd <- dd[complete.cases(dd), ]
      plot_tukey(.glht_mcp(aov(y ~ g, data = dd), "g", "Tukey"), cos) })
  ),
  opts_ui = list(posthoc = list(label = "Método post-hoc", type = "select", choices = c("Tukey (glht)" = "tukey", "Dunnett (vs. control)" = "dunnett", "Bonferroni" = "bonferroni", "Benjamini-Hochberg" = "bh"), value = "tukey")),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- .find(checks, "normalidad"); lev <- .find(checks, "homogeneidad")
    paste0("La normalidad de ", CTX(v$response), " se evaluó con la prueba de Shapiro-Wilk (", .sw_txt(sw), ") y Q-Q plot",
           if (isTRUE(sw$met)) ", sin evidencia de distribución diferente a la normal" else ", con evidencia de desviación",
           ". La homogeneidad de varianzas fue evaluada por la prueba de Levene (", .lev_txt(lev), "). ",
           "La ", CTX("variable de respuesta"), " entre los ", s$k, " grupos de ", CTX(v$group), " ",
           if (s$p < ALPHA) "fue significativamente diferente" else "no mostró diferencias significativas",
           " (F(", s$df1, ", ", s$df2, ") = ", fmt_num(s$F), ", p ", fmt_p(s$p), ", ω² = ", fmt_num(s$omega2, 3), ", efecto ", efecto_eta(s$omega2), "). ",
           "El análisis post-hoc ", CTX("(Tukey/Bonferroni) reveló diferencias entre …"), " [Cohen, 1988; Tomczak, 2014].")
  }
))

# =============================================================================
# 6. Kruskal-Wallis
# =============================================================================
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
    dd <- data.frame(y = num(d, v$response), g = fac(d, v$group)); dd <- dd[complete.cases(dd), ]
    KW <- kruskal.test(y ~ g, data = dd); N <- nrow(dd); K <- nlevels(dd$g)
    eta2 <- (unname(KW$statistic) - K + 1) / (N - K)
    list(main = KW,
         stats = list(H = unname(KW$statistic), df = unname(KW$parameter), p = KW$p.value, eta2 = eta2, N = N, K = K),
         tables = list("Resultado de kruskal.test()" = data.frame(H = unname(KW$statistic), gl = unname(KW$parameter), p = KW$p.value, `η²_H` = eta2, check.names = FALSE),
                       "Descriptiva por grupo" = desc_numeric(dd$y, dd$g)),
         notes = "η²_H = (H − k + 1)/(n − k). El libro recomienda reportar media ± DE junto con mediana e IQR.")
  },
  posthoc = function(res, d, v, opts) {
    dd <- data.frame(y = num(d, v$response), g = fac(d, v$group)); dd <- dd[complete.cases(dd), ]
    tabs <- list()
    mc <- tryCatch(pgirmess::kruskalmc(y ~ g, data = dd, probs = 0.05), error = function(e) NULL)
    if (!is.null(mc)) tabs[["Siegel y Castellan (pgirmess::kruskalmc)"]] <- data.frame(Comparación = rownames(mc$dif.com), mc$dif.com, row.names = NULL, check.names = FALSE)
    pw <- pairwise.wilcox.test(dd$y, dd$g, p.adjust.method = "bonferroni", exact = FALSE)
    df <- as.data.frame(as.table(pw$p.value)); df <- df[!is.na(df$Freq), ]; names(df) <- c("Grupo A", "Grupo B", "p ajustada (Bonferroni)")
    tabs[["pairwise.wilcox.test (Bonferroni)"]] <- df
    list(tables = tabs, notes = "Post-hoc solo si Kruskal-Wallis global es significativo.")
  },
  plots = list(boxplot = list(label = "Boxplot por grupo (mediana)", fn = function(d, v, res, cos) plot_box(d, v$response, v$group, cos, "median"))),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- .find(checks, "normalidad"); lev <- .find(checks, "homogeneidad")
    paste0("La normalidad de ", CTX(v$response), " fue evaluada con la prueba de Shapiro-Wilk (", .sw_txt(sw), ") y la homogeneidad de varianzas con la prueba de Levene (", .lev_txt(lev), "); ",
           if (isFALSE(sw$met) || isFALSE(lev$met)) "los supuestos paramétricos se incumplen, por lo que " else "",
           "la comparación entre los ", s$K, " grupos se realizó con la prueba de Kruskal-Wallis. La diferencia ",
           if (s$p < ALPHA) "fue significativa" else "no fue significativa", " entre los grupos ", CTX("(media ± DE de cada grupo)"),
           " (H[", s$df, "] = ", fmt_num(s$H), ", p ", fmt_p(s$p), ", η² = ", fmt_num(s$eta2), ", efecto ", efecto_eta(s$eta2), "). ",
           "El análisis post-hoc (Siegel y Castellan) ", CTX("reveló diferencias significativas entre …"), " [Kruskal y Wallis, 1952; Field, 2012; Siegel y Castellan, 1988].")
  }
))

# =============================================================================
# 7. ANOVA de medidas repetidas (ezANOVA)  (alt: Friedman)
# =============================================================================
.rm_long <- function(d, v) { long <- wide_to_long(d, v$measures, v$subject); long <- long[complete.cases(long), ]
  # sujetos completos
  full <- names(which(table(long$Px) == length(v$measures))); long <- long[long$Px %in% full, ]; long$Px <- droplevels(long$Px); long }

ado_register(list(
  id = "rm_anova", family = "medias", order = 7,
  label = "ANOVA de medidas repetidas (≥ 3 tiempos)",
  description = "Compara ≥ 3 mediciones del mismo sujeto. Alternativa: ANOVA de Friedman.",
  chapter = "ANOVA de medidas repetidas en R",
  roles = list(measures = role("Columnas de medición (≥ 3 tiempos/condiciones)", "numeric", multiple = TRUE, min = 3),
               subject  = role("Identificador del sujeto (opcional)", "id", required = FALSE)),
  example = list(file = "DXKTw.txt", prep = function(d) subset(d, Grupo == 1), vars = list(measures = c("X6", "X7", "X8"), subject = "Px"),
                 note = "Dolor (NRS) a las 6, 7 y 8 h en el grupo Dexketoprofeno."),
  alternative = "friedman",
  descriptives = function(d, v) desc_table1(.rm_long(d, v), "Response", "Times"),
  assumptions = function(d, v) {
    long <- .rm_long(d, v)
    ez <- tryCatch(suppressWarnings(ez::ezANOVA(data = long, dv = Response, wid = Px, within = Times, detailed = TRUE, type = 3)), error = function(e) NULL)
    checks <- list(check_normality(long$Response, "la respuesta"), if (!is.null(ez)) check_sphericity(ez) else check_obj("esfericidad", "Esfericidad (Mauchly)", NA, "No se pudo calcular ezANOVA."))
    # regla del libro: esfericidad violada → corrección GG (sigue paramétrico); normalidad violada → Friedman
    list(checks = checks, decision = decide_path(checks[1], TRUE))
  },
  run = function(d, v, opts) {
    long <- .rm_long(d, v)
    ez <- suppressWarnings(ez::ezANOVA(data = long, dv = Response, wid = Px, within = Times, detailed = TRUE, type = 3))
    A <- ez$ANOVA; row <- which(A$Effect == "Times")
    F <- A$F[row]; DFn <- A$DFn[row]; DFd <- A$DFd[row]; p <- A$p[row]; ges <- A$ges[row]
    sph <- ez[["Mauchly's Test for Sphericity"]]; corr <- ez[["Sphericity Corrections"]]
    p_gg <- if (!is.null(corr)) as.numeric(corr[["p[GG]"]][1]) else NA
    use_gg <- !is.null(sph) && as.numeric(sph$p[1]) < ALPHA
    list(main = ez,
         stats = list(F = F, df1 = DFn, df2 = DFd, p = if (use_gg) p_gg else p, p_raw = p, p_gg = p_gg, ges = ges, use_gg = use_gg,
                      W_m = if (!is.null(sph)) as.numeric(sph$W[1]) else NA, p_m = if (!is.null(sph)) as.numeric(sph$p[1]) else NA, k = length(v$measures)),
         tables = c(list("ANOVA (ezANOVA, tipo III)" = A), if (!is.null(sph)) list("Prueba de Mauchly" = sph), if (!is.null(corr)) list("Correcciones de esfericidad" = corr)),
         notes = if (use_gg) "Mauchly significativo: se reporta p corregida de Greenhouse-Geisser." else "Esfericidad no violada: F sin corrección.")
  },
  posthoc = function(res, d, v, opts) {
    long <- .rm_long(d, v); meth <- if (is.null(opts$posthoc)) "bonferroni" else opts$posthoc
    pw <- pairwise.t.test(long$Response, long$Times, p.adjust.method = if (meth == "bh") "BH" else "bonferroni", paired = TRUE)
    df <- as.data.frame(as.table(pw$p.value)); df <- df[!is.na(df$Freq), ]; names(df) <- c("Tiempo A", "Tiempo B", "p ajustada")
    list(tables = setNames(list(df), paste0("pairwise.t.test apareado (", if (meth == "bh") "Benjamini-Hochberg" else "Bonferroni", ")")),
         notes = "Bonferroni recomendado cuando la esfericidad se viola (libro).")
  },
  plots = list(perfiles = list(label = "Perfiles por sujeto + boxplot", fn = function(d, v, res, cos) plot_paired(.rm_long(d, v), cos, "mean"))),
  opts_ui = list(posthoc = list(label = "Método post-hoc", type = "select", choices = c("Bonferroni" = "bonferroni", "Benjamini-Hochberg" = "bh"), value = "bonferroni")),
  report = function(res, d, v, checks) {
    s <- res$stats
    paste0(CTX("La variable"), " fue evaluada en ", s$k, " tiempos ", CTX("(describir)"), ". El supuesto de esfericidad fue evaluado con la prueba de Mauchly (W = ", fmt_num(s$W_m, 3), ", p ", fmt_p(s$p_m), ") y ",
           if (isTRUE(s$use_gg)) "fue violado, por lo que se reporta la corrección de Greenhouse-Geisser" else "no fue violado",
           ". El ANOVA de medidas repetidas ", if (s$p < ALPHA) "mostró" else "no mostró", " diferencia significativa entre los tiempos (F[", fmt_num(s$df1, 0), ", ", fmt_num(s$df2, 0), "] = ", fmt_num(s$F), ", p ", fmt_p(s$p), ", η²g = ", fmt_num(s$ges, 3), "). ",
           "Los análisis post-hoc con corrección de Bonferroni ", CTX("confirmaron/revelaron …"), " [Field, 2012].")
  }
))

# =============================================================================
# 8. ANOVA de Friedman
# =============================================================================
ado_register(list(
  id = "friedman", family = "medias", order = 8,
  label = "ANOVA de Friedman (≥ 3 tiempos)",
  description = "Alternativa no paramétrica al ANOVA de medidas repetidas.",
  chapter = "ANOVA de Friedman en R",
  roles = list(measures = role("Columnas de medición (≥ 3 tiempos/condiciones)", "numeric", multiple = TRUE, min = 3),
               subject  = role("Identificador del sujeto (opcional)", "id", required = FALSE)),
  example = list(file = "DXKTw.txt", prep = function(d) subset(d, Grupo == 1), vars = list(measures = c("X0", "X1", "X2"), subject = "Px"),
                 note = "Dolor (NRS) a las 0, 1 y 2 h en el grupo Dexketoprofeno."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(.rm_long(d, v), "Response", "Times"),
  assumptions = function(d, v) {
    M <- as.matrix(d[v$measures]); M <- M[complete.cases(M), , drop = FALSE]
    dif_vars <- combn(ncol(M), 2, function(ij) var(M[, ij[1]] - M[, ij[2]]))
    txt <- paste0("Las varianzas de las diferencias entre pares de tiempos oscilan entre ", fmt_num(min(dif_vars)), " y ", fmt_num(max(dif_vars)),
                  if (max(dif_vars) / max(min(dif_vars), 1e-9) > 2) ": heterogéneas, lo que compromete la esfericidad (justifica Friedman)." else ": relativamente homogéneas.")
    checks <- list(check_obj("esfericidad", "Varianzas de las diferencias", NA, txt, min = min(dif_vars), max = max(dif_vars)))
    list(checks = checks, decision = list(parametric = FALSE, all_met = FALSE, violated = character(0), reason = "Prueba no paramétrica."))
  },
  run = function(d, v, opts) {
    M <- as.matrix(d[v$measures]); M <- M[complete.cases(M), , drop = FALSE]
    Fr <- friedman.test(M); N <- nrow(M); k <- ncol(M); W <- unname(Fr$statistic) / (N * (k - 1))
    list(main = Fr,
         stats = list(chi2 = unname(Fr$statistic), df = unname(Fr$parameter), p = Fr$p.value, W = W, N = N, k = k),
         tables = list("Resultado de friedman.test()" = data.frame(`χ²F` = unname(Fr$statistic), gl = unname(Fr$parameter), p = Fr$p.value, `W de Kendall` = W, check.names = FALSE),
                       "Descriptiva por tiempo" = desc_numeric(as.vector(M), rep(colnames(M), each = N))),
         notes = "W de Kendall = χ²F / [N(k − 1)]: 0–0.1 trivial, 0.1–0.3 pequeño, 0.3–0.5 mediano, > 0.5 grande.")
  },
  posthoc = function(res, d, v, opts) {
    M <- as.matrix(d[v$measures]); M <- M[complete.cases(M), , drop = FALSE]
    mc <- tryCatch(pgirmess::friedmanmc(M), error = function(e) NULL)
    if (is.null(mc)) return(list(tables = list(), notes = "pgirmess::friedmanmc no disponible."))
    list(tables = list("friedmanmc (pgirmess)" = data.frame(Comparación = rownames(mc$dif.com), mc$dif.com, row.names = NULL, check.names = FALSE)),
         notes = "Post-hoc solo si Friedman global es significativo (el libro lo muestra a modo de ejemplo).")
  },
  plots = list(perfiles = list(label = "Perfiles por sujeto + boxplot (mediana)", fn = function(d, v, res, cos) plot_paired(.rm_long(d, v), cos, "median"))),
  report = function(res, d, v, checks) {
    s <- res$stats; ch <- checks[[1]]
    paste0("Dado que las varianzas de las diferencias entre los tiempos fueron heterogéneas (rango: ", fmt_num(ch$min), " a ", fmt_num(ch$max), "), lo que compromete el supuesto de esfericidad, ",
           "la diferencia en ", CTX("la variable"), " entre las ", s$k, " evaluaciones se evaluó con el ANOVA de Friedman. ", CTX("Describir media ± DE (mediana) por tiempo"), ". La diferencia ",
           if (s$p < ALPHA) "fue estadísticamente significativa" else "no fue estadísticamente significativa",
           " (χ²F[", s$df, "] = ", fmt_num(s$chi2), ", p ", fmt_p(s$p), ", W = ", fmt_num(s$W, 3), ", efecto ", efecto_W(s$W), "). ",
           "El análisis post-hoc (friedmanmc) ", CTX("no reveló/reveló diferencias entre …"), " [Friedman, 1937; Field, 2012].")
  }
))

# =============================================================================
# 9. ANOVA factorial (dos vías) — sin alternativa NP en el libro
# =============================================================================
ado_register(list(
  id = "factorial", family = "medias", order = 9,
  label = "ANOVA factorial (dos factores)",
  description = "Efectos principales e interacción de dos factores. Suma de cuadrados tipo III (car::Anova).",
  chapter = "ANOVA de dos vías en R",
  roles = list(response = role("Variable de respuesta (numérica)", "numeric"),
               factorA  = role("Factor A", "factor"), factorB = role("Factor B", "factor")),
  example = list(file = "VickersF.txt", vars = list(response = "VH", factorA = "Group", factorB = "Cure.time"),
                 note = "Dureza Vickers según cemento (Group) y protocolo de polimerización (Cure.time)."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, v$response, c(v$factorA, v$factorB)),
  assumptions = function(d, v) {
    y <- num(d, v$response); A <- fac(d, v$factorA); B <- fac(d, v$factorB)
    checks <- list(check_normality(y, v$response), check_homogeneity(y, interaction(A, B), "las celdas A×B"))
    list(checks = checks, decision = decide_path(checks, FALSE))
  },
  run = function(d, v, opts) {
    dd <- data.frame(y = num(d, v$response), A = fac(d, v$factorA), B = fac(d, v$factorB)); dd <- dd[complete.cases(dd), ]
    contrasts(dd$A) <- contr.sum(nlevels(dd$A)); contrasts(dd$B) <- contr.sum(nlevels(dd$B))
    m <- aov(y ~ A * B, data = dd); a3 <- car::Anova(m, type = "III")
    a3 <- a3[rownames(a3) != "(Intercept)", ]
    SS_R <- a3["Residuals", "Sum Sq"]; MS_R <- SS_R / a3["Residuals", "Df"]; SS_T <- sum(a3[, "Sum Sq"])
    eff <- rownames(a3)[rownames(a3) != "Residuals"]
    omega <- sapply(eff, function(e) (a3[e, "Sum Sq"] - a3[e, "Df"] * MS_R) / (SS_T + MS_R))
    tab <- data.frame(Fuente = c(v$factorA, v$factorB, paste(v$factorA, "×", v$factorB), "Residuales")[seq_len(nrow(a3))], a3, check.names = FALSE, row.names = NULL)
    list(main = m,
         stats = list(FA = a3["A", "F value"], pA = a3["A", "Pr(>F)"], dfA = a3["A", "Df"],
                      FB = a3["B", "F value"], pB = a3["B", "Pr(>F)"], dfB = a3["B", "Df"],
                      FAB = a3["A:B", "F value"], pAB = a3["A:B", "Pr(>F)"], dfAB = a3["A:B", "Df"], dfR = a3["Residuals", "Df"],
                      wA = omega["A"], wB = omega["B"], wAB = omega["A:B"]),
         tables = list("ANOVA tipo III (car::Anova)" = tab, "ω² parcial" = data.frame(Efecto = c(v$factorA, v$factorB, "Interacción"), `ω² parcial` = as.numeric(omega), check.names = FALSE)),
         notes = "Contrastes suma (contr.sum) para SS tipo III. ω² parcial = (SS − df·MS_R)/(SS_T + MS_R).")
  },
  posthoc = function(res, d, v, opts) {
    dd <- data.frame(y = num(d, v$response), A = fac(d, v$factorA), B = fac(d, v$factorB)); dd <- dd[complete.cases(dd), ]
    dd$celda <- interaction(dd$A, dd$B, sep = " × ")
    tk <- TukeyHSD(aov(y ~ celda, data = dd))$celda
    list(tables = list("Efectos simples: Tukey HSD sobre celdas A×B" = data.frame(Comparación = rownames(tk), tk, row.names = NULL, check.names = FALSE)),
         notes = "Con muchas celdas, el libro prefiere el post-hoc sobre el factor combinado a los contrastes planeados.")
  },
  plots = list(
    interaccion = list(label = "Gráfico de interacción (medias ± IC)", fn = function(d, v, res, cos) plot_interaction(d, v$response, v$factorA, v$factorB, cos)),
    boxA = list(label = "Boxplot por factor A", fn = function(d, v, res, cos) plot_box(d, v$response, v$factorA, cos)),
    boxB = list(label = "Boxplot por factor B", fn = function(d, v, res, cos) plot_box(d, v$response, v$factorB, cos)),
    qq   = list(label = "Q-Q de la respuesta", fn = function(d, v, res, cos) plot_qq(num(d, v$response), paste("Q-Q:", v$response), cos))
  ),
  report = function(res, d, v, checks) {
    s <- res$stats; sw <- .find(checks, "normalidad"); lev <- .find(checks, "homogeneidad")
    paste0("La normalidad de ", CTX(v$response), " fue evaluada (Shapiro-Wilk, ", .sw_txt(sw), ") y la homogeneidad de varianzas con la prueba de Levene (", .lev_txt(lev), "). ",
           "El ANOVA factorial reveló ", if (s$pA < ALPHA) "un efecto principal significativo" else "ausencia de efecto principal", " de ", CTX(v$factorA), " (F[", s$dfA, ",", s$dfR, "] = ", fmt_num(s$FA), ", p ", fmt_p(s$pA), ", ω² = ", fmt_num(s$wA, 3), ") y ",
           if (s$pB < ALPHA) "un efecto significativo" else "ausencia de efecto", " de ", CTX(v$factorB), " (F[", s$dfB, ",", s$dfR, "] = ", fmt_num(s$FB), ", p ", fmt_p(s$pB), ", ω² = ", fmt_num(s$wB, 3), "). ",
           "La interacción ", v$factorA, " × ", v$factorB, if (s$pAB < ALPHA) " fue significativa" else " no fue significativa", " (F[", s$dfAB, ",", s$dfR, "] = ", fmt_num(s$FAB), ", p ", fmt_p(s$pAB), ", ω² = ", fmt_num(s$wAB, 3), ")",
           if (s$pAB < ALPHA) paste0(", indicando que ", CTX("el efecto de un factor depende del nivel del otro (describir)"), ". Los análisis de efectos simples ", CTX("revelaron …")) else ".",
           " [Field, 2012; Cohen, 1988].")
  }
))

# =============================================================================
# 10. ANCOVA — sin alternativa NP en el libro
# =============================================================================
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
    y <- num(d, v$response); g <- fac(d, v$group); x <- num(d, v$covariate)
    dd <- data.frame(y, g, x); dd <- dd[complete.cases(dd), ]
    ind <- summary(aov(x ~ g, data = dd))[[1]]; p_ind <- ind[1, "Pr(>F)"]
    rs <- car::Anova(aov(y ~ g + x + g:x, data = dd), type = "III"); p_slope <- rs["g:x", "Pr(>F)"]
    checks <- list(
      check_normality(y, v$response), check_normality(x, v$covariate), check_homogeneity(y, g),
      check_obj("independencia", "Independencia de la covariable", p_ind > ALPHA,
                paste0("ANOVA de ", v$covariate, " ~ ", v$group, ": F = ", fmt_num(ind[1, "F value"], 3), ", p ", fmt_p(p_ind), if (p_ind > ALPHA) ": la covariable no difiere entre grupos (independiente del factor)." else ": la covariable difiere entre grupos — condición necesaria del ANCOVA comprometida."), F = ind[1, "F value"], p = p_ind),
      check_obj("pendientes", "Homogeneidad de pendientes", p_slope > ALPHA,
                paste0("Interacción ", v$group, " × ", v$covariate, ": p ", fmt_p(p_slope), if (p_slope > ALPHA) ": el supuesto de homogeneidad de pendientes se cumple." else ": las pendientes difieren entre grupos; el ANCOVA clásico no es apropiado."), p = p_slope))
    list(checks = checks, decision = decide_path(checks[4:5], FALSE))
  },
  run = function(d, v, opts) {
    dd <- data.frame(y = num(d, v$response), g = fac(d, v$group), x = num(d, v$covariate)); dd <- dd[complete.cases(dd), ]
    contrasts(dd$g) <- contr.sum(nlevels(dd$g))
    m <- aov(y ~ g + x, data = dd); a3 <- car::Anova(m, type = "III"); a3 <- a3[rownames(a3) != "(Intercept)", ]
    SS_R <- a3["Residuals", "Sum Sq"]
    eta_g <- a3["g", "Sum Sq"] / (a3["g", "Sum Sq"] + SS_R); eta_x <- a3["x", "Sum Sq"] / (a3["x", "Sum Sq"] + SS_R)
    adj <- tryCatch({ dd2 <- dd; names(dd2) <- c("y", v$group, "x"); m2 <- aov(as.formula(paste("y ~", bq(v$group), "+ x")), data = dd2)
      e <- effects::effect(v$group, m2, se = TRUE); s <- summary(e)
      data.frame(Grupo = names(s$effect), `Media ajustada` = as.numeric(s$effect), `IC inf` = as.numeric(s$lower), `IC sup` = as.numeric(s$upper), check.names = FALSE) }, error = function(e) NULL)
    list(main = m,
         stats = list(Fg = a3["g", "F value"], dfg = a3["g", "Df"], pg = a3["g", "Pr(>F)"], Fx = a3["x", "F value"], dfx = a3["x", "Df"], px = a3["x", "Pr(>F)"], dfR = a3["Residuals", "Df"], eta_g = eta_g, eta_x = eta_x),
         tables = c(list("ANCOVA tipo III (car::Anova)" = data.frame(Fuente = c(v$group, v$covariate, "Residuales")[seq_len(nrow(a3))], a3, check.names = FALSE, row.names = NULL),
                         "η² parcial" = data.frame(Efecto = c(v$group, v$covariate), `η² parcial` = c(eta_g, eta_x), check.names = FALSE)),
                    if (!is.null(adj)) list("Medias ajustadas (effects)" = adj)),
         notes = "η² parcial = SS_efecto/(SS_efecto + SS_R).")
  },
  posthoc = function(res, d, v, opts) {
    dd <- data.frame(y = num(d, v$response), g = fac(d, v$group), x = num(d, v$covariate)); dd <- dd[complete.cases(dd), ]
    m <- aov(y ~ g + x, data = dd); g <- .glht_mcp(m, "g", "Tukey")
    list(tables = list("Tukey sobre medias ajustadas (multcomp::glht)" = .glht_table(g)), glht = g, notes = "")
  },
  plots = list(
    boxplot = list(label = "Boxplot de la respuesta por grupo", fn = function(d, v, res, cos) plot_box(d, v$response, v$group, cos)),
    ajustadas = list(label = "Medias ajustadas (effects)", fn = function(d, v, res, cos) {
      dd <- data.frame(y = num(d, v$response), g = fac(d, v$group), x = num(d, v$covariate)); dd <- dd[complete.cases(dd), ]
      names(dd)[2] <- v$group; m <- aov(as.formula(paste("y ~", bq(v$group), "+ x")), data = dd)
      plot_adjusted(effects::effect(v$group, m, se = TRUE), cos) }),
    qq = list(label = "Q-Q respuesta", fn = function(d, v, res, cos) plot_qq(num(d, v$response), paste("Q-Q:", v$response), cos))
  ),
  report = function(res, d, v, checks) {
    s <- res$stats; swy <- checks[[1]]; swx <- checks[[2]]; lev <- checks[[3]]; ind <- checks[[4]]; pen <- checks[[5]]
    paste0("La normalidad fue evaluada con la prueba de Shapiro-Wilk para ", CTX(v$response), " (", .sw_txt(swy), ") y ", CTX(v$covariate), " (", .sw_txt(swx), "); la homogeneidad de varianzas se evaluó con la prueba de Levene (", .lev_txt(lev), "). ",
           "La covariable (", CTX(v$covariate), ") ", if (isTRUE(ind$met)) "no difirió" else "difirió", " entre grupos (ANOVA: F = ", fmt_num(ind$F, 3), ", p ", fmt_p(ind$p), ")",
           if (isTRUE(ind$met)) ", confirmando su independencia del factor experimental" else "", ". El modelo ANCOVA reveló ",
           if (s$pg < ALPHA) "un efecto significativo" else "ausencia de efecto", " del grupo (F[", s$dfg, ",", s$dfR, "] = ", fmt_num(s$Fg), ", p ", fmt_p(s$pg), ", η²parcial = ", fmt_num(s$eta_g, 3), ") y ",
           if (s$px < ALPHA) "un efecto significativo" else "ausencia de efecto", " de la covariable (F[", s$dfx, ",", s$dfR, "] = ", fmt_num(s$Fx), ", p ", fmt_p(s$px), ", η²parcial = ", fmt_num(s$eta_x, 3), "). ",
           "El supuesto de homogeneidad de pendientes ", if (isTRUE(pen$met)) "fue confirmado" else "no se cumplió", " (interacción grupo × covariable: p ", fmt_p(pen$p), "). ",
           "Las comparaciones de Tukey sobre las medias ajustadas ", CTX("mostraron …"), " [Field, 2012; Miller y Chapman, 2001].")
  }
))
