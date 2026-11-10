# =============================================================================
# analyses/categoricos.R - Chi-cuadrado (alt: Fisher), residuos (r×c), McNemar
# Reglas del libro (caps. 37–42):
#  - Esperadas ≥ 5 en todas las celdas → chi-cuadrado de Pearson (cap. 38.4).
#  - Alguna esperada < 5 → Fisher exacta. El N total no es criterio
#    (indicación del autor, 2026-11-10; versión Final101126 del libro).
#  - Se reporta el chi-cuadrado sin Yates; φ y el poder siempre sin corrección
#    (cap. 40.5.1). Yates solo con n < 40 (cap. 38.6).
#  - 2×2: residuos estandarizados de CrossTable (sresid) en la tabla;
#    r×c: residuos estandarizados ajustados (ASR, |ASR| > 1.96) como post-hoc
#    (cap. 41), interpretables solo si la prueba global es significativa.
# =============================================================================

.ctab <- function(d, v) {
  dd <- data.frame(a = fac(d, v$rowvar), b = fac(d, v$colvar)); dd <- dd[complete.cases(dd), ]
  tab <- table(dd$a, dd$b); names(dimnames(tab)) <- c(v$rowvar, v$colvar); tab
}
# V de Cramér (cap. 41.3): sqrt(chi² / (N (k − 1))), chi² sin corrección
.cramer <- function(tab) { ct <- suppressWarnings(chisq.test(tab, correct = FALSE)); N <- sum(tab); k <- min(dim(tab)); sqrt(as.numeric(ct$statistic) / (N * (k - 1))) }
.tab_df <- function(tab) { m <- addmargins(tab); data.frame(Fila = rownames(m), as.data.frame.matrix(m), check.names = FALSE, row.names = NULL) }
.cat_mat_df <- function(m, digits = 3) data.frame(Fila = rownames(m), round(as.data.frame.matrix(m), digits), check.names = FALSE, row.names = NULL)

# Regla chi² vs. Fisher del libro (sustituye a check_expected_counts compartido)
.cat_check <- function(tab) {
  ex <- suppressWarnings(chisq.test(tab, correct = FALSE))$expected
  mn <- min(ex); N <- sum(tab); is22 <- all(dim(tab) == 2)
  met <- mn >= 5
  motivos <- if (!met) "alguna frecuencia esperada fue < 5" else character(0)
  txt <- paste0("La frecuencia esperada mínima es ", fmt_num(mn), " (N = ", N, "; tabla ", nrow(tab), "×", ncol(tab), "). ",
                if (met) "Todas las frecuencias esperadas son ≥ 5: procede el chi-cuadrado de Pearson."
                else "Alguna frecuencia esperada es < 5: se usa la prueba exacta de Fisher.",
                " Criterio (cap. 39): frecuencias esperadas ≥ 5 en todas las celdas → chi-cuadrado; alguna < 5 → Fisher exacta.")
  check_obj("esperadas", "Frecuencias esperadas ≥ 5 (chi² vs. Fisher)", met, txt,
            min_expected = mn, N = N, is22 = is22, motivos = motivos)
}

# fisher.test() exacto, como en el libro; si la tabla r×c agota la memoria de
# trabajo se reintenta con más workspace y, solo si aun así falla, se simula.
.cat_fisher <- function(tab) {
  tryCatch(fisher.test(tab), error = function(e)
    tryCatch(fisher.test(tab, workspace = 2e7), error = function(e2) {
      r <- fisher.test(tab, simulate.p.value = TRUE, B = 1e5); r$simulado <- TRUE; r }))
}

# Gráfico de mosaico (cap. 40, Figura 40.1): ancho ∝ total de fila, altura ∝
# frecuencia condicional; cada combinación con su color, % del total y n.
.cat_mosaic <- function(tab, rowvar, colvar, cos = list()) {
  cos <- cos_merge(cos); N <- sum(tab); rs <- rowSums(tab); xr <- c(0, cumsum(rs / N))
  nr <- nrow(tab); nc <- ncol(tab); recs <- NULL; pal <- character(0)
  for (i in seq_len(nr)) {
    yr <- c(0, cumsum(tab[i, ] / max(rs[i], 1)))
    sh <- grDevices::colorRampPalette(c(.colN(cos, nr)[i], "white"))(nc + 1)[seq_len(nc)]
    for (j in seq_len(nc)) {
      lab <- paste0(rowvar, " ", rownames(tab)[i], " · ", colvar, " ", colnames(tab)[j])
      recs <- rbind(recs, data.frame(xmin = xr[i], xmax = xr[i + 1], ymin = yr[j], ymax = yr[j + 1], combo = lab,
                                     pct = round(tab[i, j] / N * 100, 1), n = as.numeric(tab[i, j]), stringsAsFactors = FALSE))
      pal[lab] <- sh[j]
    }
  }
  recs$combo <- factor(recs$combo, levels = names(pal))
  centros <- (xr[-length(xr)] + xr[-1]) / 2
  p <- ggplot(recs) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = combo), colour = "white", linewidth = 1.2) +
    geom_text(aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = paste0(pct, "%\n(n=", n, ")")), size = 3.3, fontface = "bold") +
    scale_fill_manual(values = pal, name = NULL, guide = guide_legend(nrow = nr, byrow = TRUE)) +
    scale_x_continuous(breaks = centros, labels = paste0(rowvar, ": ", rownames(tab), "\n(n=", rs, ")"), expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0))
  cos$legend <- TRUE
  .finish(p, cos, NULL, paste0("Proporción de ", colvar, " (dentro de cada grupo)")) + theme(panel.grid = element_blank(), axis.text.x = element_text(face = "bold"))
}

# Cierra el texto con la cita entre paréntesis antes del punto final (formato del libro)
.cat_cite <- function(txt, cita) paste0(sub("[.[:space:]]+$", "", txt), " (", cita, ").")

.cat_msg_plot <- function(msg) ggplot() + annotate("text", x = 0, y = 0, label = msg, size = 4.5) + theme_void()

ado_register(list(
  id = "chi2", family = "categoricos", order = 1, label = "Chi-cuadrado de Pearson (asociación)",
  description = "Asociación entre dos variables categóricas (2×2 o r×c). Alternativa: Fisher exacta si alguna frecuencia esperada < 5.",
  chapter = "Chi-cuadrado con Prueba exacta de Fisher en R",
  roles = list(rowvar = role("Variable explicativa (filas)", "factor"), colvar = role("Variable de respuesta (columnas)", "factor")),
  example = list(file = "datAsthma.txt", prep = function(d) { d$Asthma <- factor(d$Asthma, levels = c("Yes", "No")); d$OrParaF <- factor(d$OrParaF, levels = c("Yes", "No")); d },
                 vars = list(rowvar = "Asthma", colvar = "OrParaF"), note = "Asma (explicativa) y parafunción oral (respuesta), 228 niños."),
  alternative = "fisher",
  descriptives = function(d, v) desc_table1(d, v$colvar, v$rowvar),
  assumptions = function(d, v) { ch <- list(.cat_check(.ctab(d, v))); list(checks = ch, decision = decide_path(ch, TRUE)) },
  run = function(d, v, opts) {
    tab <- .ctab(d, v); N <- sum(tab); is22 <- all(dim(tab) == 2); k <- min(dim(tab))
    ct  <- suppressWarnings(chisq.test(tab, correct = FALSE))                    # Pearson (el que reporta el libro)
    cty <- if (is22) suppressWarnings(chisq.test(tab, correct = TRUE)) else NULL # Yates (CrossTable lo muestra en 2×2)
    yates <- is22 && isTRUE(opts$yates)
    mn <- if (yates) cty else ct
    ft <- .cat_fisher(tab)                                                       # Fisher exacta junto al chi² (CrossTable fisher = TRUE)
    chi_u <- as.numeric(ct$statistic)
    V   <- sqrt(chi_u / (N * (k - 1)))
    phi <- if (is22) sqrt(chi_u / N) else NA                                     # φ siempre con correct = FALSE (cap. 40.5.1)
    pw  <- if (is22) tryCatch(pwr::pwr.chisq.test(w = phi, N = N, df = 1, sig.level = 0.05)$power, error = function(e) NA) else NA
    n80 <- if (is22) tryCatch(pwr::pwr.chisq.test(w = phi, power = 0.80, df = 1, sig.level = 0.05)$N, error = function(e) NA) else NA
    asr <- ct$stdres
    sig_cells <- which(abs(asr) > 1.96, arr.ind = TRUE)
    sig_txt <- if (nrow(sig_cells)) apply(sig_cells, 1, function(ix) paste0(rownames(asr)[ix[1]], " × ", colnames(asr)[ix[2]], " (ASR = ", fmt_num(asr[ix[1], ix[2]]), ")")) else character(0)
    st <- list(chi2 = unname(mn$statistic), df = unname(mn$parameter), p = mn$p.value, yates = yates,
               chi2_u = chi_u, p_u = ct$p.value, chi2_y = if (is22) unname(cty$statistic) else NA, p_y = if (is22) cty$p.value else NA,
               p_fisher = ft$p.value, fisher_sim = isTRUE(ft$simulado),
               N = N, V = V, k = k, phi = phi, min_exp = min(ct$expected), power = pw, n80 = n80, is22 = is22,
               OR = if (is22) unname(ft$estimate) else NA, OR_lo = if (is22) ft$conf.int[1] else NA, OR_hi = if (is22) ft$conf.int[2] else NA,
               sig_cells = sig_txt)
    pruebas <- data.frame(Prueba = c("Chi-cuadrado de Pearson", if (is22) "Chi-cuadrado con corrección de Yates", "Prueba exacta de Fisher"),
                          `Estadístico` = c(chi_u, if (is22) st$chi2_y, NA), gl = c(unname(ct$parameter), if (is22) unname(cty$parameter), NA),
                          p = c(ct$p.value, if (is22) cty$p.value, ft$p.value), check.names = FALSE)
    tabs <- list("Tabla de contingencia (con márgenes)" = .tab_df(tab),
                 "Frecuencias esperadas" = .cat_mat_df(ct$expected, 3),
                 "Porcentaje por fila" = .cat_mat_df(100 * prop.table(tab, 1), 1))
    if (is22) tabs[["Residuos estandarizados (CrossTable, sresid)"]] <- .cat_mat_df(ct$residuals, 3)
    tabs[["Pruebas (como CrossTable: chi², Yates y Fisher)"]] <- pruebas
    if (is22) {
      tabs[["Odds ratio (fisher.test) y φ"]] <- data.frame(OR = st$OR, `IC 95% inf` = st$OR_lo, `IC 95% sup` = st$OR_hi, `φ` = phi, check.names = FALSE)
      tabs[["Poder (pwr.chisq.test, w = φ, gl = 1)"]] <- data.frame(w = phi, N = N, `Poder post hoc` = pw, `N para 80 % de poder` = ceiling(n80), check.names = FALSE)
    } else {
      tabs[["V de Cramér"]] <- data.frame(`V de Cramér` = V, N = N, k = k, check.names = FALSE)
    }
    notes <- c(if (yates) "Estadístico principal con corrección de Yates (el libro la recomienda solo con n < 40)." else "Chi-cuadrado sin corrección de Yates (correct = FALSE), como en el reporte del libro.",
               if (is22) "φ y el poder se calculan siempre con el chi-cuadrado sin corrección. El OR es el de fisher.test() (máxima verosimilitud condicional). El poder post hoc solo contextualiza un resultado no significativo." ,
               if (isTRUE(ft$simulado)) "La tabla es demasiado grande para el cálculo exacto de Fisher: su p se obtuvo por simulación (B = 100 000).")
    list(main = mn, stats = st, tables = tabs, notes = paste(notes, collapse = " "))
  },
  posthoc = function(res, d, v, opts) {
    tab <- .ctab(d, v)
    if (all(dim(tab) == 2)) return(NULL)   # el libro aplica los ASR a tablas mayores que 2×2 (cap. 41)
    asr <- suppressWarnings(chisq.test(tab))$stdres
    df <- as.data.frame(as.table(round(asr, 3))); names(df) <- c(v$rowvar, v$colvar, "ASR"); df$Significativo <- ifelse(abs(df$ASR) > 1.96, "sí (|ASR| > 1.96)", "no")
    list(tables = list("Residuos estandarizados ajustados (chisq.test()$stdres)" = df),
         notes = "Los ASR son post hoc: interpretar solo si la prueba global es significativa. |ASR| > 1.96 → la celda difiere de lo esperado (α = 0.05); signo + = más observaciones de lo esperado, − = menos.")
  },
  plots = list(barras = list(label = "Barras apiladas de proporción", fn = function(d, v, res, cos) plot_bar_prop(d, v$rowvar, v$colvar, cos)),
               mosaico = list(label = "Gráfico de mosaico", fn = function(d, v, res, cos) .cat_mosaic(.ctab(d, v), v$rowvar, v$colvar, cos)),
               asr = list(label = "Mapa de calor de residuos (ASR, tablas r×c)", fn = function(d, v, res, cos) { tab <- .ctab(d, v)
                 if (all(dim(tab) == 2)) return(.cat_msg_plot("El libro usa el mapa de ASR solo en tablas mayores que 2×2."))
                 plot_heat_asr(tab, suppressWarnings(chisq.test(tab))$stdres, cos) })),
  opts_ui = list(yates = list(label = "Reportar con corrección de Yates (solo 2×2; el libro la recomienda si n < 40)", type = "checkbox", value = FALSE)),
  report = function(res, d, v, checks) {
    s <- res$stats; sig <- s$p < ALPHA
    intro <- paste0("La relación entre ", CTX(v$rowvar), " y ", CTX(v$colvar), " se evaluó mediante la prueba de chi-cuadrado de Pearson", if (s$yates) " con corrección de Yates" else "", ". ")
    supu <- if (s$min_exp >= 5) paste0("La frecuencia mínima esperada fue de ", fmt_num(s$min_exp), ", cumpliendo el supuesto de frecuencias esperadas ≥ 5. ")
            else paste0("La frecuencia mínima esperada fue de ", fmt_num(s$min_exp), " (< 5), por lo que el libro recomienda la prueba exacta de Fisher (p ", fmt_p(s$p_fisher), "). ")
    if (s$is22) {
      .cat_cite(paste0(intro, supu, "La asociación entre ambas variables ", if (sig) "resultó estadísticamente significativa" else "no resultó estadísticamente significativa",
             " (χ²[", s$df, "] = ", fmt_num(s$chi2), ", p ", fmt_p(s$p), ", OR = ", fmt_num(s$OR), ", IC 95 % [", fmt_num(s$OR_lo), ", ", fmt_num(s$OR_hi), "], φ = ", fmt_num(s$phi),
             "), lo que indica un tamaño del efecto ", efecto_r(s$phi), " (referencia de Cohen: φ ≈ 0.1 pequeño, 0.3 mediano, 0.5 grande)", if (!sig) "; asociación no concluyente" else "", ".",
             if (!sig && !is.na(s$power)) paste0(" El poder estadístico para detectar un efecto de esta magnitud fue ", if (s$power < 0.8) "bajo" else "adecuado", " (≈ ", round(100 * s$power), " %)",
                                                 if (s$power < 0.8) paste0(", por lo que la ausencia de significancia no descarta una asociación real; detectar el efecto observado con un 80 % de poder habría requerido N ≈ ", ceiling(s$n80)) else "", ".") else "",
             ""), "Field et al., 2012; Ranganathan et al., 2015")
    } else {
      .cat_cite(paste0(intro, "En esta muestra (N = ", s$N, ") ", if (s$min_exp >= 5) paste0("todas las frecuencias esperadas son ≥ 5 (la menor es ", fmt_num(s$min_exp), "), por lo que la prueba apropiada es el chi-cuadrado de Pearson. ") else paste0("la menor frecuencia esperada es ", fmt_num(s$min_exp), " (< 5), por lo que el libro recomienda la prueba exacta de Fisher (p ", fmt_p(s$p_fisher), "). "),
             "La relación entre ", CTX(v$rowvar), " y ", CTX(v$colvar), if (sig) " resulta estadísticamente significativa" else " no resulta estadísticamente significativa",
             " (χ² = ", fmt_num(s$chi2), ", gl = ", s$df, ", p ", fmt_p(s$p), "; V de Cramér = ", fmt_num(s$V), ")",
             if (s$k == 2) paste0(", un tamaño del efecto ", efecto_r(s$V), " (referencia: V ≈ 0.1 pequeño, 0.3 mediano, 0.5 grande)") else "", ". ",
             if (sig) { if (length(s$sig_cells)) paste0("Los residuos estandarizados ajustados (|ASR| > 1.96) identifican las celdas que contribuyen a la asociación: ", paste(s$sig_cells, collapse = "; "), " (signo positivo = más observaciones de lo esperado; negativo = menos). ")
                        else "Ninguna celda presenta un residuo estandarizado ajustado significativo (|ASR| > 1.96). " }
             else "Los residuos estandarizados ajustados son un análisis post hoc, por lo que solo se interpretan cuando la prueba global resulta significativa. ",
             ""), "Field et al., 2012")
    }
  }
))

ado_register(list(
  id = "fisher", family = "categoricos", order = 2, label = "Prueba exacta de Fisher",
  description = "Si alguna frecuencia esperada < 5 (tablas r×c: Freeman-Halton).",
  chapter = "Prueba exacta de Fisher",
  roles = list(rowvar = role("Variable explicativa (filas)", "factor"), colvar = role("Variable de respuesta (columnas)", "factor")),
  example = list(file = "datAsthma.txt", prep = function(d) { d$Asthma <- factor(d$Asthma, levels = c("Yes", "No")); d$Bruxism <- factor(d$Bruxism); d },
                 vars = list(rowvar = "Asthma", colvar = "Bruxism"), note = "Asma y bruxismo."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, v$colvar, v$rowvar),
  assumptions = function(d, v) { ch <- list(.cat_check(.ctab(d, v))); list(checks = ch, decision = list(parametric = FALSE, all_met = FALSE, violated = character(0), reason = "Prueba exacta.")) },
  run = function(d, v, opts) {
    tab <- .ctab(d, v); is22 <- all(dim(tab) == 2)
    ft <- .cat_fisher(tab); ex <- suppressWarnings(chisq.test(tab, correct = FALSE))$expected
    st <- list(p = ft$p.value, N = sum(tab), is22 = is22, sim = isTRUE(ft$simulado), min_exp = min(ex),
               OR = if (is22) unname(ft$estimate) else NA, OR_lo = if (is22) ft$conf.int[1] else NA, OR_hi = if (is22) ft$conf.int[2] else NA)
    tabs <- list("Tabla de contingencia" = .tab_df(tab), "Frecuencias esperadas" = .cat_mat_df(ex, 2),
                 "Fisher" = if (is22) data.frame(p = st$p, OR = st$OR, `IC 95% inf` = st$OR_lo, `IC 95% sup` = st$OR_hi, check.names = FALSE)
                            else data.frame(`Método` = "Freeman-Halton (r × c)", p = st$p, check.names = FALSE))
    list(main = ft, stats = st, tables = tabs,
         notes = paste0(if (is22) "El OR de fisher.test() es el estimador de máxima verosimilitud condicional (difiere ligeramente del OR muestral)." else "fisher.test() aplica el método de Freeman-Halton a tablas r × c.",
                        if (st$sim) " La tabla es demasiado grande para el cálculo exacto: p obtenido por simulación (B = 100 000)." else ""))
  },
  posthoc = NULL,
  plots = list(barras = list(label = "Barras apiladas de proporción", fn = function(d, v, res, cos) plot_bar_prop(d, v$rowvar, v$colvar, cos)),
               mosaico = list(label = "Gráfico de mosaico", fn = function(d, v, res, cos) .cat_mosaic(.ctab(d, v), v$rowvar, v$colvar, cos))),
  report = function(res, d, v, checks) {
    s <- res$stats; ch <- checks[[1]]; mot <- ch$motivos
    paste0("Dado que ", if (length(mot)) paste(mot, collapse = " y ") else "se optó por una prueba exacta", " (frecuencia mínima esperada = ", fmt_num(ch$min_expected), ", N = ", s$N, "), ",
           "la asociación entre ", CTX(v$rowvar), " y ", CTX(v$colvar), " se evaluó con la prueba exacta de Fisher", if (!s$is22) " (método de Freeman-Halton para tablas r × c)" else "",
           ", que ", if (s$p < ALPHA) "resultó estadísticamente significativa" else "no resultó estadísticamente significativa", " (p ", fmt_p(s$p),
           if (s$is22) paste0("; OR = ", fmt_num(s$OR), ", IC 95 % [", fmt_num(s$OR_lo), ", ", fmt_num(s$OR_hi), "]") else "", ") (Field et al., 2012; West & Hankin, 2008).")
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
    lv <- union(levels(a), levels(b)); if (length(lv) != 2) stop("McNemar requiere una respuesta binaria (2 categorías) en ambas columnas; se encontraron: ", paste(lv, collapse = ", "))
    a <- factor(a, levels = lv); b <- factor(b, levels = lv)
    tab <- table(a, b); names(dimnames(tab)) <- v$measures[1:2]
    mc <- mcnemar.test(tab, correct = FALSE)
    bb <- tab[1, 2]; cc <- tab[2, 1]; OR <- if (cc > 0) bb / cc else Inf
    cero <- bb == 0 || cc == 0; OR_h <- (bb + 0.5) / (cc + 0.5)
    n1 <- sum(tab[1, ]); n2 <- sum(tab[, 1]); N <- sum(tab)
    mc_tab <- data.frame(`χ²` = unname(mc$statistic), gl = unname(mc$parameter), p = mc$p.value, `Discordantes b` = bb, `Discordantes c` = cc, `OR condicional (b/c)` = OR, check.names = FALSE)
    if (cero) mc_tab[["OR con corrección de continuidad (+0.5)"]] <- OR_h
    list(main = mc,
         stats = list(chi2 = unname(mc$statistic), df = unname(mc$parameter), p = mc$p.value, b = bb, c = cc, OR = OR, OR_h = OR_h, cero = cero, N = N,
                      n1 = n1, n2 = n2, p1 = n1 / N, p2 = n2 / N, lvl = lv[1], c1 = v$measures[1], c2 = v$measures[2]),
         tables = list("Tabla apareada" = .tab_df(tab), "McNemar" = mc_tab),
         notes = paste0("OR condicional = b/c (solo pares discordantes). Sin corrección de continuidad en la prueba (correct = FALSE), como en el libro.",
                        if (cero) " Como una celda discordante es 0, se añade el OR con corrección de continuidad (sumar 0.5 a cada celda)." else ""))
  },
  posthoc = NULL,
  plots = list(barras = list(label = "Proporción por condición", fn = function(d, v, res, cos) {
    long <- data.frame(cond = factor(rep(v$measures[1:2], each = nrow(d)), levels = v$measures[1:2]), resp = factor(c(as.character(d[[v$measures[1]]]), as.character(d[[v$measures[2]]]))))
    long <- long[complete.cases(long), ]; names(long) <- c("Condición", "Respuesta")
    plot_bar_prop(long, "Condición", "Respuesta", cos) })),
  report = function(res, d, v, checks) {
    s <- res$stats
    paste0("Un mismo conjunto de ", s$N, " sujetos participó en ambas condiciones (", CTX(s$c1), " y ", CTX(s$c2), "), por lo que la respuesta se analizó con la prueba de McNemar sobre la tabla apareada. ",
           "La proporción de ", CTX(s$lvl), " fue de ", s$n1, "/", s$N, " (", fmt_num(100 * s$p1, 1), " %) bajo ", s$c1, " y de ", s$n2, "/", s$N, " (", fmt_num(100 * s$p2, 1), " %) bajo ", s$c2, "; la diferencia ",
           if (s$p < ALPHA) "resultó estadísticamente significativa" else "no resultó estadísticamente significativa", " (χ²[", s$df, "] = ", fmt_num(s$chi2, 1), ", p ", fmt_p(s$p), "). ",
           "De los ", s$N, " sujetos, ", s$b, " cambiaron solo en una dirección (", s$lvl, " únicamente bajo ", s$c1, ") y ", s$c, " en la contraria",
           if (s$cero && (s$b + s$c) > 0) ", un patrón de cambio completamente unidireccional" else "", ". ",
           if (is.finite(s$OR) && s$b > 0) paste0("El odds ratio condicional fue ", fmt_num(s$OR))
           else paste0("En consecuencia, el odds ratio condicional es ", if (s$b == 0 && s$c == 0) "indeterminado" else if (s$c == 0) "indefinido (infinito)" else "0"),
           if (s$cero) paste0("; con corrección de continuidad se aproxima a ", fmt_num(s$OR_h, 1)) else "", " (McNemar, 1947).")
  }
))
