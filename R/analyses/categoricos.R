# =============================================================================
# analyses/categoricos.R — Chi-cuadrado (alt: Fisher), residuos (r×c), McNemar
# Reglas del libro: esperadas ≥ 5 en todas las celdas y N ≥ 40 → chi²;
# si no → Fisher exacta. Yates solo con n < 40. ASR |>1.96| = celda que difiere.
# =============================================================================

.ctab <- function(d, v) {
  dd <- data.frame(a = fac(d, v$rowvar), b = fac(d, v$colvar)); dd <- dd[complete.cases(dd), ]
  tab <- table(dd$a, dd$b); names(dimnames(tab)) <- c(v$rowvar, v$colvar); tab
}
.cramer <- function(tab) { ct <- suppressWarnings(chisq.test(tab, correct = FALSE)); N <- sum(tab); k <- min(dim(tab)); sqrt(as.numeric(ct$statistic) / (N * (k - 1))) }
.tab_df <- function(tab) { m <- addmargins(tab); data.frame(Fila = rownames(m), as.data.frame.matrix(m), check.names = FALSE, row.names = NULL) }

ado_register(list(
  id = "chi2", family = "categoricos", order = 1, label = "Chi-cuadrado de Pearson (asociación)",
  description = "Asociación entre dos variables categóricas (2×2 o r×c). Alternativa: Fisher exacta.",
  chapter = "Chi-cuadrado con Prueba exacta de Fisher en R",
  roles = list(rowvar = role("Variable explicativa (filas)", "factor"), colvar = role("Variable de respuesta (columnas)", "factor")),
  example = list(file = "datAsthma.txt", prep = function(d) { d$Asthma <- factor(d$Asthma, levels = c("Yes", "No")); d$OrParaF <- factor(d$OrParaF, levels = c("Yes", "No")); d },
                 vars = list(rowvar = "Asthma", colvar = "OrParaF"), note = "Asma (explicativa) y parafunción oral (respuesta), 228 niños."),
  alternative = "fisher",
  descriptives = function(d, v) desc_table1(d, v$colvar, v$rowvar),
  assumptions = function(d, v) { ch <- list(check_expected_counts(.ctab(d, v))); list(checks = ch, decision = decide_path(ch, TRUE)) },
  run = function(d, v, opts) {
    tab <- .ctab(d, v); N <- sum(tab); is22 <- all(dim(tab) == 2)
    yates <- isTRUE(opts$yates)
    ct <- suppressWarnings(chisq.test(tab, correct = yates))
    V <- .cramer(tab); phi <- if (is22) sqrt(as.numeric(ct$statistic) / N) else NA
    OR <- if (is22) fisher.test(tab) else NULL
    st <- list(chi2 = unname(ct$statistic), df = unname(ct$parameter), p = ct$p.value, N = N, V = V, phi = phi, min_exp = min(ct$expected),
               OR = if (!is.null(OR)) unname(OR$estimate) else NA, OR_lo = if (!is.null(OR)) OR$conf.int[1] else NA, OR_hi = if (!is.null(OR)) OR$conf.int[2] else NA, is22 = is22)
    pw <- if (is22 && !is.na(phi)) tryCatch(pwr::pwr.chisq.test(w = phi, N = N, df = 1, sig.level = 0.05)$power, error = function(e) NA) else NA
    st$power <- pw
    tabs <- list("Tabla de contingencia (con márgenes)" = .tab_df(tab),
                 "Frecuencias esperadas" = data.frame(Fila = rownames(ct$expected), round(as.data.frame.matrix(ct$expected), 2), check.names = FALSE, row.names = NULL),
                 "Chi-cuadrado" = data.frame(`χ²` = st$chi2, gl = st$df, p = st$p, `V de Cramér` = V, φ = phi, check.names = FALSE))
    if (is22) tabs[["Odds ratio (fisher.test)"]] <- data.frame(OR = st$OR, `IC 95% inf` = st$OR_lo, `IC 95% sup` = st$OR_hi, `Poder post hoc` = pw, check.names = FALSE)
    list(main = ct, stats = st, tables = tabs, notes = if (yates) "Con corrección de Yates (recomendada solo si n < 40)." else "Sin corrección de Yates (correct = FALSE), como en el libro.")
  },
  posthoc = function(res, d, v, opts) {
    tab <- .ctab(d, v); ct <- suppressWarnings(chisq.test(tab)); asr <- ct$stdres
    df <- as.data.frame(as.table(round(asr, 3))); names(df) <- c(v$rowvar, v$colvar, "ASR"); df$Significativo <- ifelse(abs(df$ASR) > 1.96, "sí (|ASR| > 1.96)", "no")
    list(tables = list("Residuos estandarizados ajustados (chisq.test()$stdres)" = df), notes = "Los ASR son post hoc: interpretar solo si la prueba global es significativa. Signo + = más de lo esperado.")
  },
  plots = list(barras = list(label = "Barras apiladas de proporción", fn = function(d, v, res, cos) plot_bar_prop(d, v$rowvar, v$colvar, cos)),
               asr = list(label = "Mapa de calor de residuos (ASR)", fn = function(d, v, res, cos) { tab <- .ctab(d, v); plot_heat_asr(tab, suppressWarnings(chisq.test(tab))$stdres, cos) })),
  opts_ui = list(yates = list(label = "Corrección de Yates (n < 40)", type = "checkbox", value = FALSE)),
  report = function(res, d, v, checks) {
    s <- res$stats
    paste0("La relación entre ", CTX(v$rowvar), " y ", CTX(v$colvar), " se evaluó mediante la prueba de chi-cuadrado de Pearson. La frecuencia mínima esperada fue de ", fmt_num(s$min_exp),
           if (s$min_exp >= 5) ", cumpliendo el supuesto de frecuencias esperadas ≥ 5" else " (< 5; se recomienda la prueba exacta de Fisher)", ". La asociación ",
           if (s$p < ALPHA) "resultó estadísticamente significativa" else "no resultó estadísticamente significativa", " (χ²[", s$df, "] = ", fmt_num(s$chi2), ", p ", fmt_p(s$p),
           if (s$is22) paste0(", OR = ", fmt_num(s$OR), ", IC 95 % [", fmt_num(s$OR_lo), ", ", fmt_num(s$OR_hi), "], φ = ", fmt_num(s$phi)) else paste0(", V de Cramér = ", fmt_num(s$V)), ")",
           if (s$is22 && !is.na(s$power)) paste0(". El poder estadístico para detectar un efecto de esta magnitud fue ", if (s$power < 0.8) "bajo" else "adecuado", " (≈ ", round(100 * s$power), " %)") else "",
           " [Field, 2012; Ranganathan, 2015].")
  }
))

ado_register(list(
  id = "fisher", family = "categoricos", order = 2, label = "Prueba exacta de Fisher",
  description = "Para frecuencias esperadas < 5 o N < 40 (tablas r×c: Freeman-Halton).",
  chapter = "Prueba exacta de Fisher",
  roles = list(rowvar = role("Variable explicativa (filas)", "factor"), colvar = role("Variable de respuesta (columnas)", "factor")),
  example = list(file = "datAsthma.txt", prep = function(d) { d$Asthma <- factor(d$Asthma, levels = c("Yes", "No")); d$Bruxism <- factor(d$Bruxism); d },
                 vars = list(rowvar = "Asthma", colvar = "Bruxism"), note = "Asma y bruxismo."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, v$colvar, v$rowvar),
  assumptions = function(d, v) { ch <- list(check_expected_counts(.ctab(d, v))); list(checks = ch, decision = list(parametric = FALSE, all_met = FALSE, violated = character(0), reason = "Prueba exacta.")) },
  run = function(d, v, opts) {
    tab <- .ctab(d, v); is22 <- all(dim(tab) == 2)
    ft <- fisher.test(tab, simulate.p.value = !is22 && prod(dim(tab)) > 10)
    st <- list(p = ft$p.value, N = sum(tab), is22 = is22, OR = if (is22) unname(ft$estimate) else NA, OR_lo = if (is22) ft$conf.int[1] else NA, OR_hi = if (is22) ft$conf.int[2] else NA, V = .cramer(tab))
    list(main = ft, stats = st,
         tables = list("Tabla de contingencia" = .tab_df(tab), "Fisher" = data.frame(p = st$p, OR = st$OR, `IC 95% inf` = st$OR_lo, `IC 95% sup` = st$OR_hi, `V de Cramér` = st$V, check.names = FALSE)),
         notes = "El OR de fisher.test() es el estimador condicional de máxima verosimilitud.")
  },
  posthoc = NULL,
  plots = list(barras = list(label = "Barras apiladas de proporción", fn = function(d, v, res, cos) plot_bar_prop(d, v$rowvar, v$colvar, cos))),
  report = function(res, d, v, checks) {
    s <- res$stats; ch <- checks[[1]]
    paste0("Dado que ", if (ch$min_expected < 5) "alguna frecuencia esperada fue < 5" else if (s$N < 40) "N < 40" else "se optó por una prueba exacta", " (mínima esperada = ", fmt_num(ch$min_expected), ", N = ", s$N, "), ",
           "la asociación entre ", CTX(v$rowvar), " y ", CTX(v$colvar), " se evaluó con la prueba exacta de Fisher, que ", if (s$p < ALPHA) "resultó significativa" else "no resultó significativa", " (p ", fmt_p(s$p),
           if (s$is22) paste0("; OR = ", fmt_num(s$OR), ", IC 95 % [", fmt_num(s$OR_lo), ", ", fmt_num(s$OR_hi), "]") else paste0("; V de Cramér = ", fmt_num(s$V)), ") [Field, 2012].")
  }
))

ado_register(list(
  id = "mcnemar", family = "categoricos", order = 3, label = "Prueba de McNemar (apareada)",
  description = "Cambio en una respuesta binaria del mismo sujeto entre dos condiciones (tabla apareada).",
  chapter = "Prueba de McNemar en R",
  roles = list(measures = role("Dos columnas binarias (condición 1 y condición 2)", "binary", multiple = TRUE, min = 2)),
  example = list(file = "datMusicT.txt",
                 prep = function(d) { m <- d[d$Group == "Music", c("Px", "Behavior")]; cc <- d[d$Group == "Control", c("Px", "Behavior")]
                   names(m)[2] <- "Music"; names(cc)[2] <- "Control"; w <- merge(m, cc, by = "Px")
                   w$Music <- factor(w$Music, levels = c(1, 2), labels = c("Colaborador", "No colaborador")); w$Control <- factor(w$Control, levels = c(1, 2), labels = c("Colaborador", "No colaborador")); w },
                 vars = list(measures = c("Music", "Control")), note = "Comportamiento (colaborador/no) de 59 niños con y sin musicoterapia."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, v$measures[1:2]),
  assumptions = function(d, v) list(checks = list(check_obj("apareado", "Diseño apareado", TRUE, "Cada sujeto aporta una fila con su respuesta en ambas condiciones; las observaciones no son independientes, por lo que corresponde McNemar y no chi-cuadrado.")), decision = list(parametric = TRUE, all_met = TRUE, violated = character(0), reason = "")),
  run = function(d, v, opts) {
    a <- fac(d, v$measures[1]); b <- fac(d, v$measures[2]); ok <- !is.na(a) & !is.na(b); a <- a[ok]; b <- b[ok]
    lv <- union(levels(a), levels(b)); a <- factor(a, levels = lv); b <- factor(b, levels = lv)
    tab <- table(a, b); names(dimnames(tab)) <- v$measures[1:2]
    mc <- mcnemar.test(tab, correct = FALSE)
    bb <- tab[1, 2]; cc <- tab[2, 1]; OR <- if (cc > 0) bb / cc else Inf; OR_h <- (bb + 0.5) / (cc + 0.5)
    p1 <- sum(tab[1, ]) / sum(tab); p2 <- sum(tab[, 1]) / sum(tab)
    list(main = mc,
         stats = list(chi2 = unname(mc$statistic), df = unname(mc$parameter), p = mc$p.value, b = bb, c = cc, OR = OR, OR_h = OR_h, N = sum(tab), p1 = p1, p2 = p2, lvl = lv[1], c1 = v$measures[1], c2 = v$measures[2]),
         tables = list("Tabla apareada" = .tab_df(tab), "McNemar" = data.frame(`χ²` = unname(mc$statistic), gl = unname(mc$parameter), p = mc$p.value, `Discordantes b` = bb, `Discordantes c` = cc, `OR condicional` = OR, `OR (Haldane)` = OR_h, check.names = FALSE)),
         notes = "OR condicional = b/c (pares discordantes). Sin corrección de continuidad, como en el libro.")
  },
  posthoc = NULL,
  plots = list(barras = list(label = "Proporción por condición", fn = function(d, v, res, cos) {
    long <- data.frame(cond = factor(rep(v$measures[1:2], each = nrow(d)), levels = v$measures[1:2]), resp = factor(c(as.character(d[[v$measures[1]]]), as.character(d[[v$measures[2]]]))))
    long <- long[complete.cases(long), ]; names(long) <- c("Condición", "Respuesta")
    plot_bar_prop(long, "Condición", "Respuesta", cos) })),
  report = function(res, d, v, checks) {
    s <- res$stats
    paste0("Un mismo conjunto de ", s$N, " sujetos participó en ambas condiciones (", CTX(s$c1), " y ", CTX(s$c2), "), por lo que la respuesta se analizó con la prueba de McNemar sobre la tabla apareada. ",
           "La proporción de ", CTX(s$lvl), " fue de ", round(100 * s$p1, 1), " % bajo ", s$c1, " y de ", round(100 * s$p2, 1), " % bajo ", s$c2, "; la diferencia ",
           if (s$p < ALPHA) "resultó estadísticamente significativa" else "no resultó estadísticamente significativa", " (χ²[", s$df, "] = ", fmt_num(s$chi2, 1), ", p ", fmt_p(s$p), "). ",
           s$b, " sujetos cambiaron en una dirección y ", s$c, " en la contraria; el odds ratio condicional fue ", if (is.finite(s$OR)) fmt_num(s$OR) else "indefinido (infinito)", " (con corrección de Haldane ≈ ", fmt_num(s$OR_h, 1), ") [McNemar, 1947].")
  }
))
