# =============================================================================
# verify_cat.R — Verificación de correlaciones (cap. 10), datos categóricos
# (caps. 37–42) y concordancia (caps. 43–47) contra las cifras impresas en el libro.
# Solo se usan salidas de la app (stats, tablas, supuestos, reporte); nunca se recalcula.
# =============================================================================
source(file.path(Sys.getenv("VERIFY_DIR"), "verify_lib.R"))
NA_APP <- "la app no lo muestra"

# ---- utilidades ---------------------------------------------------------------
# Primer número capturado por un patrón (grupo 1) en el texto del reporte
rx <- function(txt, pattern, group = 1) {
  m <- regmatches(txt, regexec(pattern, txt, perl = TRUE))[[1]]
  if (length(m) < group + 1) return(NA)
  as.numeric(gsub("−", "-", m[group + 1]))
}
NUM <- "(-?[0-9]+(?:\\.[0-9]+)?)"
# Celdas de una fila de un table1 (html) por etiqueta; k = ocurrencia
t1_row <- function(html, label, k = 1) {
  html <- paste(as.character(html), collapse = "\n")
  rows <- regmatches(html, gregexpr("(?s)<tr>.*?</tr>", html, perl = TRUE))[[1]]
  lab <- sub("(?s)^<tr>\\s*<t[dh][^>]*>(.*?)</t[dh]>.*$", "\\1", rows, perl = TRUE)
  idx <- which(trimws(gsub("<[^>]+>", "", lab)) == label)
  if (length(idx) < k) return(NA)
  cells <- regmatches(rows[idx[k]], gregexpr("(?s)<td[^>]*>.*?</td>", rows[idx[k]], perl = TRUE))[[1]]
  trimws(gsub("<[^>]+>", "", cells))[-1]
}
# "68 (89.5%)" → c(68, 89.5) ; "35.6 (12.8)" → c(35.6, 12.8)
t1_nums <- function(cell) as.numeric(regmatches(cell, gregexpr("-?[0-9]+(\\.[0-9]+)?", cell))[[1]])
# Valor de una tabla de matriz (Fila + columnas) por fila y columna
mat_val <- function(x, table_regex, row, col) tab_val(x, table_regex, col, row)
# Matriz → data.frame de casos (factores en el orden de la matriz)
mk <- function(m) {
  rn <- rownames(m); cn <- colnames(m); g <- expand.grid(i = seq_along(rn), j = seq_along(cn))
  n <- m[cbind(g$i, g$j)]
  data.frame(F = factor(rep(rn[g$i], n), levels = rn), C = factor(rep(cn[g$j], n), levels = cn))
}
# Dos evaluadores a partir de a, b, c, d (tabla 2×2 de acuerdo)
mk_k <- function(a, b, cc, dd) data.frame(r1 = rep(c("Si", "Si", "No", "No"), c(a, b, cc, dd)),
                                          r2 = rep(c("Si", "No", "Si", "No"), c(a, b, cc, dd)), stringsAsFactors = FALSE)
na_chk <- function(an, item, book, line, note = NA_APP) chk(an, item, book, NA, line = line, note = note)

# =============================================================================
# CAPÍTULO 10 — CORRELACIONES
# =============================================================================
# ---- 10.1–10.2 Edad vs. AINEs (5 pacientes) ----------------------------------
pat <- data.frame(Edad = c(22, 24, 35, 48, 49), Consumo_AINEs = c(2, 3, 5, 11, 13))
x <- run_app("pearson", pat, list(x = "Edad", y = "Consumo_AINEs"))
s <- x$res$stats; rep <- x$report
ds <- tryCatch(x$spec$descriptives(pat, list(x = "Edad", y = "Consumo_AINEs")), error = function(e) NULL)
ms1 <- t1_nums(t1_row(ds, "Mean (SD)", 1)); ms2 <- t1_nums(t1_row(ds, "Mean (SD)", 2))
chk("pearson", "Edad/AINEs: media de Edad (descriptivos)", "35.6", ms1[1], line = 1911)
chk("pearson", "Edad/AINEs: media de AINEs (descriptivos)", "6.8", ms2[1], line = 1913)
na_chk("pearson", "Edad/AINEs: desviaciones y productos cruzados (tabla)", "65.28", line = 1926)
na_chk("pearson", "Edad/AINEs: covarianza (manual)", "61.4", line = 1934)
na_chk("pearson", "Edad/AINEs: covarianza (cov())", "61.4", line = 1937)
chk("pearson", "Edad/AINEs: DE Edad (descriptivos, table1 a 3 cifras)", "12.779", ms1[2], line = 1948, note = "table1 muestra 3 cifras significativas (12.8)")
chk("pearson", "Edad/AINEs: DE AINEs (descriptivos, table1 a 3 cifras)", "4.919", ms2[2], line = 1950, note = "table1 muestra 3 cifras significativas (4.92)")
chk("pearson", "Edad/AINEs: r (manual)", "0.9767", s$est, line = 1952)
chk("pearson", "Edad/AINEs: r (cor())", "0.9767", s$est, line = 1955)
chk("pearson", "Edad/AINEs: r en leyenda Fig. 10.1", "0.98", s$est, line = 1981)
chk("pearson", "Edad/AINEs: t", "7.8852", tab_val(x, "cor.test", "t"), line = 1990)
chk("pearson", "Edad/AINEs: df", "3", tab_val(x, "cor.test", "gl"), line = 1990)
chk("pearson", "Edad/AINEs: p", "0.004251", tab_val(x, "cor.test", "p"), line = 1990)
chk("pearson", "Edad/AINEs: IC 95% inf", "0.6830314", tab_val(x, "cor.test", "IC 95% inf"), line = 1993)
chk("pearson", "Edad/AINEs: IC 95% sup", "0.9985274", tab_val(x, "cor.test", "IC 95% sup"), line = 1993)
chk("pearson", "Edad/AINEs: r (cor.test)", "0.9767144", tab_val(x, "cor.test", "r"), line = 1996)
chk("pearson", "Edad/AINEs: r en texto/reporte", "0.98", rx(rep, paste0("\\(r = ", NUM)), line = 1997)
chk("pearson", "Edad/AINEs: p en texto/reporte (libro: p < 0.001)", "<0.001", s$p, line = 1997,
    note = "reporte de la app: 'p = 0.004'; el propio libro imprime p-value = 0.004251 (línea 1990)")
chk("pearson", "Edad/AINEs: R² (tabla)", "0.96", tab_val(x, "cor.test", "R²"), line = 1998, note = "libro: R² = 0.98² = 0.96 (eleva r ya redondeado)")
chk("pearson", "Edad/AINEs: R² (reporte)", "0.96", rx(rep, paste0("R² = ", NUM)), line = 1998)
chk("pearson", "Edad/AINEs: % de variación explicada (reporte)", "96", rx(rep, paste0("el ", NUM, " % de la variación")), line = 1998)

# ---- 10.2.1 IMC vs. tiempo quirúrgico ----------------------------------------
dx <- book_data("DXKTw.txt")
x <- run_app("pearson", dx, list(x = "Imc", y = "Tqx")); s <- x$res$stats; rep <- x$report
chk("pearson", "IMC/Tqx: decisión (Q-Q → Pearson)", "TRUE", as.character(isTRUE(x$assumptions$decision$all_met)), line = 2044)
chk("pearson", "IMC/Tqx: t", "0.2107", tab_val(x, "cor.test", "t"), line = 2052)
chk("pearson", "IMC/Tqx: df", "58", tab_val(x, "cor.test", "gl"), line = 2052)
chk("pearson", "IMC/Tqx: p", "0.8339", tab_val(x, "cor.test", "p"), line = 2052)
chk("pearson", "IMC/Tqx: IC 95% inf", "-0.2278691", tab_val(x, "cor.test", "IC 95% inf"), line = 2055)
chk("pearson", "IMC/Tqx: IC 95% sup", "0.2796169", tab_val(x, "cor.test", "IC 95% sup"), line = 2055)
chk("pearson", "IMC/Tqx: r", "0.02765578", tab_val(x, "cor.test", "r"), line = 2058)
chk("pearson", "IMC/Tqx: r en leyenda Fig. 10.3 (≈)", "0.03", s$est, line = 2065)
chk("pearson", "IMC/Tqx: r (reporte)", "0.028", rx(rep, paste0("\\(r = ", NUM)), line = 2066)
chk("pearson", "IMC/Tqx: p (reporte)", "0.83", rx(rep, paste0("p = ", NUM)), line = 2066)

# ---- 10.3 Spearman: edad dental vs. cronológica ------------------------------
cr <- book_data("datcron.txt")
x <- run_app("pearson", cr, list(x = "Dental_Age", y = "Cron_Age"))
chk("spearman", "Dental/Cron: decisión Q-Q (Pearson → Spearman)", "FALSE", as.character(isTRUE(x$assumptions$decision$all_met)), line = 2092)
x <- run_app("spearman", cr, list(x = "Dental_Age", y = "Cron_Age")); s <- x$res$stats; rep <- x$report
chk("spearman", "Dental/Cron: S", "3304546", tab_val(x, "spearman", "S"), line = 2100)
chk("spearman", "Dental/Cron: p", "<2.2e-16", tab_val(x, "spearman", "p"), line = 2100)
chk("spearman", "Dental/Cron: rho", "0.8556836", tab_val(x, "spearman", "rho"), line = 2104)
chk("spearman", "Dental/Cron: ρ en leyenda Fig. 10.5 (≈)", "0.86", s$est, line = 2111)
chk("spearman", "Dental/Cron: ρ (reporte)", "0.856", rx(rep, paste0("ρ = ", NUM)), line = 2112)
chk("spearman", "Dental/Cron: p < 0.0001 (texto)", "<0.0001", s$p, line = 2112, note = "el reporte de la app escribe 'p < 0.001'")

# ---- 10.4 Kendall: edad cronológica vs. maduración ósea ----------------------
x <- run_app("kendall", cr, list(x = "Cron_Age", y = "Bone_Maturation")); s <- x$res$stats; rep <- x$report
chk("kendall", "τ de Kendall", "0.5933", tab_val(x, "Kendall", "Estimación", "Kendall τ"), line = 2131)
chk("kendall", "p de Kendall (round 4)", "0", tab_val(x, "Kendall", "p", "Kendall τ"), line = 2131)
chk("kendall", "ρ de Spearman (comparativa)", "0.7202", tab_val(x, "Kendall", "Estimación", "Spearman ρ"), line = 2134)
chk("kendall", "p de Spearman (round 4)", "0", tab_val(x, "Kendall", "p", "Spearman ρ"), line = 2134)
chk("kendall", "N pacientes (leyenda Fig. 10.6)", "516", s$n, line = 2149)
chk("kendall", "τ en leyenda Fig. 10.6", "0.59", s$tau, line = 2149)
chk("kendall", "ρ en leyenda Fig. 10.6", "0.72", s$rho, line = 2149)
chk("kendall", "niveles de maduración ósea", "5", s$k_y, line = 2116)
chk("kendall", "τ (reporte)", "0.593", rx(rep, paste0("τ = ", NUM)), line = 2131)
chk("kendall", "ρ (reporte)", "0.720", rx(rep, paste0("ρ = ", NUM)), line = 2134)

# ---- 10.5 Correlación parcial: Tqx ~ Imc | Ed -------------------------------
x <- run_app("parcial", dx, list(x = "Tqx", y = "Imc", z = "Ed")); s <- x$res$stats; rep <- x$report
xb <- run_app("pearson", dx, list(x = "Tqx", y = "Imc"))
nb <- "la app de correlación parcial no lo muestra; verificado con la entrada Pearson (Tqx, Imc)"
chk("parcial", "bivariada Tqx~Imc: t", "0.2107", tab_val(xb, "cor.test", "t"), line = 2165, note = nb)
chk("parcial", "bivariada Tqx~Imc: df", "58", tab_val(xb, "cor.test", "gl"), line = 2165, note = nb)
chk("parcial", "bivariada Tqx~Imc: p", "0.8339", tab_val(x, "Correlaciones", "p", 1), line = 2165)
chk("parcial", "bivariada Tqx~Imc: IC 95% inf", "-0.2278691", tab_val(xb, "cor.test", "IC 95% inf"), line = 2168, note = nb)
chk("parcial", "bivariada Tqx~Imc: IC 95% sup", "0.2796169", tab_val(xb, "cor.test", "IC 95% sup"), line = 2168, note = nb)
chk("parcial", "bivariada Tqx~Imc: r", "0.02765578", tab_val(x, "Correlaciones", "r", 1), line = 2171)
chk("parcial", "r parcial (Tqx ~ Imc | Ed)", "-0.0452", tab_val(x, "Correlaciones", "r", 2), line = 2181)
chk("parcial", "Tqx ~ Ed: r", "0.329", tab_val(x, "Correlaciones", "r", 3), line = 2186)
chk("parcial", "Tqx ~ Ed: p", "0.01", tab_val(x, "Correlaciones", "p", 3), line = 2186)
chk("parcial", "Imc ~ Ed: r", "0.211", tab_val(x, "Correlaciones", "r", 4), line = 2188)
chk("parcial", "Imc ~ Ed: p", "0.105", tab_val(x, "Correlaciones", "p", 4), line = 2188)
chk("parcial", "texto: r bivariada ≈ 0.03 (reporte)", "0.03", rx(rep, paste0("fue r = ", NUM)), line = 2189)
chk("parcial", "texto: r parcial ≈ −0.05", "-0.05", s$r_partial, line = 2189, note = "reporte de la app: r = -0.045")
chk("parcial", "texto: r parcial (reporte, 3 dec.)", "-0.05", rx(rep, paste0("parcial fue r = ", NUM)), line = 2189)
chk("parcial", "texto: Ed ~ Tqx r ≈ 0.33 (reporte)", "0.33", rx(rep, paste0("\\[Tqx\\] \\(r = ", NUM)), line = 2189)
chk("parcial", "texto: Ed ~ Imc r ≈ 0.21 (reporte)", "0.21", rx(rep, paste0("\\[Imc\\] \\(r = ", NUM)), line = 2189)

# ---- 10.6 PPS -----------------------------------------------------------------
x <- run_app("pps", dx, list(target = "Tqx", predictors = c("Classcx", "Imc", "Ed", "Sexo", "Grupo"))); rep <- x$report
ntol <- "PPS por validación cruzada sin semilla en el libro: se compara con tolerancia ±0.02"
chk("pps", "PPS Classcx → Tqx", "0.2788096", tab_val(x, "PPS", "PPS x → y", "Classcx"), line = 2200, tol = 0.02, note = ntol)
chk("pps", "PPS Classcx → Tqx (texto ≈)", "0.28", tab_val(x, "PPS", "PPS x → y", "Classcx"), line = 2201, tol = 0.02, note = ntol)
chk("pps", "PPS Classcx → Tqx (texto asimetría ≈)", "0.28", rx(rep, paste0("Classcx = ", NUM)), line = 2205, tol = 0.02, note = ntol)
chk("pps", "PPS Tqx → Classcx (inverso ≈)", "0.49", tab_val(x, "PPS", "PPS y → x", "Classcx"), line = 2205, tol = 0.02, note = ntol)
chk("pps", "PPS Tqx → Classcx (reporte ≈)", "0.49", rx(rep, paste0("Tqx → Classcx = ", NUM)), line = 2205, tol = 0.02, note = ntol)

# =============================================================================
# CAPÍTULOS 37–39 — tabla 2×2 de 200 pacientes (28, 10, 48, 114)
# =============================================================================
obs <- matrix(c(28, 10, 48, 114), 2, byrow = TRUE, dimnames = list(c("Experimental", "Control"), c("Evento: Sí", "Evento: No")))
x <- run_app("chi2", mk(obs), list(rowvar = "F", colvar = "C")); s <- x$res$stats
TC <- "Tabla de contingencia"
chk("chi2", "37.1 marginal fila Experimental", "38", mat_val(x, TC, "Experimental", "Sum"), line = 7367)
chk("chi2", "37.1 marginal fila Control", "162", mat_val(x, TC, "Control", "Sum"), line = 7371)
chk("chi2", "37.1 marginal col Sí", "76", mat_val(x, TC, "Sum", "Evento: Sí"), line = 7373)
chk("chi2", "37.1 marginal col No", "124", mat_val(x, TC, "Sum", "Evento: No"), line = 7374)
chk("chi2", "37.1 N", "200", mat_val(x, TC, "Sum", "Sum"), line = 7375)
FE <- "Frecuencias esperadas"
chk("chi2", "38.2 esperada Exp/Sí", "14.44", mat_val(x, FE, "Experimental", "Evento: Sí"), line = 7511)
chk("chi2", "38.2 esperada Exp/No", "23.56", mat_val(x, FE, "Experimental", "Evento: No"), line = 7511)
chk("chi2", "38.2 esperada Con/Sí", "61.56", mat_val(x, FE, "Control", "Evento: Sí"), line = 7512)
chk("chi2", "38.2 esperada Con/No", "100.44", mat_val(x, FE, "Control", "Evento: No"), line = 7512)
PR <- "Pruebas"
chk("chi2", "38.2 χ² manual", "25.356", tab_val(x, PR, "Estadístico", 1), line = 7520)
chk("chi2", "38.2 df manual", "1", tab_val(x, PR, "gl", 1), line = 7520)
chk("chi2", "38.2 p manual (round 4)", "0", tab_val(x, PR, "p", 1), line = 7520)
chk("chi2", "38.2 chisq.test χ²", "25.356", s$chi2, line = 7530)
chk("chi2", "38.2 chisq.test df", "1", s$df, line = 7530)
chk("chi2", "38.2 chisq.test p", "4.767e-07", s$p, line = 7530)
na_chk("chi2", "38.3 OR muestral (A/B)/(C/D)", "6.65", line = 7540, note = "la app no lo muestra: solo da el OR condicional de fisher.test() (6.58)")
na_chk("chi2", "38.3 OR muestral en frase", "6.65", line = 7543, note = "la app no lo muestra: solo da el OR condicional de fisher.test() (6.58)")
chk("chi2", "38.6 χ² sin corregir (texto)", "25.36", s$chi2_u, line = 7571)
chk("chi2", "38.6 χ² Yates (texto)", "23.52", s$chi2_y, line = 7571)
chk("chi2", "38.6 Yates χ² (tabla Pruebas)", "23.52", tab_val(x, PR, "Estadístico", 2), line = 7578)
chk("chi2", "38.6 Yates df", "1", tab_val(x, PR, "gl", 2), line = 7578)
chk("chi2", "38.6 Yates p", "1.236e-06", tab_val(x, PR, "p", 2), line = 7578)
xy <- run_app("chi2", mk(obs), list(rowvar = "F", colvar = "C"), list(yates = TRUE))
chk("chi2", "38.6 opción Yates: χ² principal", "23.52", xy$res$stats$chi2, line = 7578)
chk("chi2", "38.6 opción Yates: p principal", "1.236e-06", xy$res$stats$p, line = 7578)
chk("chi2", "38.6 opción Yates: χ² en reporte", "23.52", rx(xy$report, paste0("χ²\\[1\\] = ", NUM)), line = 7578)
na_chk("chi2", "38.7 L χ² (razón de verosimilitud)", "24.932", line = 7589)
na_chk("chi2", "38.7 p de L χ² (round 4)", "0", line = 7589)
na_chk("chi2", "38.7.1 devianza modelo nulo", "265.63", line = 7610)
na_chk("chi2", "38.7.1 devianza modelo completo", "240.69", line = 7611)
na_chk("chi2", "38.7.1 LRT deviance", "24.932", line = 7611)
na_chk("chi2", "38.7.1 LRT p", "5.94e-07", line = 7611)
na_chk("chi2", "38.7.1 diferencia de devianza (texto)", "24.93", line = 7614)

# ---- 39.2 Fisher sobre la misma tabla ----------------------------------------
x <- run_app("fisher", mk(obs), list(rowvar = "F", colvar = "C")); s <- x$res$stats
chk("fisher", "39.2 p", "1.312e-06", tab_val(x, "^Fisher$", "p"), line = 7655)
chk("fisher", "39.2 IC 95% inf", "2.837773", tab_val(x, "^Fisher$", "IC 95% inf"), line = 7658)
chk("fisher", "39.2 IC 95% sup", "16.429686", tab_val(x, "^Fisher$", "IC 95% sup"), line = 7658)
chk("fisher", "39.2 OR (fisher.test)", "6.579265", tab_val(x, "^Fisher$", "OR"), line = 7661)
chk("fisher", "39.2 OR condicional (texto ≈)", "6.58", rx(x$report, paste0("OR = ", NUM)), line = 7662)
na_chk("fisher", "39.2 OR muestral (texto)", "6.65", line = 7662, note = "la app no lo muestra: solo da el OR condicional de fisher.test()")
chk("fisher", "39.2 todas las esperadas ≥ 5 (texto)", "TRUE", as.character(isTRUE(x$assumptions$checks[[1]]$met)), line = 7662)

# ---- 39.4 tabla con frecuencia esperada baja ---------------------------------
ts <- matrix(c(4, 2, 1, 10), 2, byrow = TRUE, dimnames = list(c("Tratamiento A", "Tratamiento B"), c("Éxito", "Fracaso")))
x <- run_app("fisher", mk(ts), list(rowvar = "F", colvar = "C"))
chk("fisher", "39.4 esperada A/Éxito", "1.76", mat_val(x, FE, "Tratamiento A", "Éxito"), line = 7695)
chk("fisher", "39.4 esperada A/Fracaso", "4.24", mat_val(x, FE, "Tratamiento A", "Fracaso"), line = 7695)
chk("fisher", "39.4 esperada B/Éxito", "3.24", mat_val(x, FE, "Tratamiento B", "Éxito"), line = 7696)
chk("fisher", "39.4 esperada B/Fracaso", "7.76", mat_val(x, FE, "Tratamiento B", "Fracaso"), line = 7696)
chk("fisher", "39.4 ¿alguna esperada < 5?", "TRUE", as.character(isFALSE(x$assumptions$checks[[1]]$met)), line = 7699)
chk("fisher", "39.4 p", "0.02763", tab_val(x, "^Fisher$", "p"), line = 7709)
chk("fisher", "39.4 IC 95% inf", "0.9571449", tab_val(x, "^Fisher$", "IC 95% inf"), line = 7712)
chk("fisher", "39.4 IC 95% sup", "1068.8016204", tab_val(x, "^Fisher$", "IC 95% sup"), line = 7712)
chk("fisher", "39.4 OR", "15.51469", tab_val(x, "^Fisher$", "OR"), line = 7715)
xc <- run_app("chi2", mk(ts), list(rowvar = "F", colvar = "C"))
chk("chi2", "39.4 chi² → Fisher (alternativa por esperada < 5)", "FALSE", as.character(isTRUE(xc$assumptions$decision$all_met)), line = 7699)
chk("chi2", "39.4 p chisq.test(correct = TRUE) (comparación)", "0.05325428", tab_val(xc, PR, "p", 2), line = 7702)
chk("chi2", "39.4 p Fisher mostrado por chi²", "0.02763", xc$res$stats$p_fisher, line = 7709)

# ---- 39.5 tabla 3×2 (Freeman-Halton) -----------------------------------------
t32 <- matrix(c(3, 1, 5, 2, 1, 8), 3, byrow = TRUE, dimnames = list(c("G1", "G2", "G3"), c("Sí", "No")))
x <- run_app("fisher", mk(t32), list(rowvar = "F", colvar = "C"))
chk("fisher", "39.5 p Freeman-Halton (3×2)", "0.02202", tab_val(x, "^Fisher$", "p"), line = 7729)

# =============================================================================
# CAPÍTULO 40 — Asma vs. parafunción oral
# =============================================================================
da <- book_data("datAsthma.txt"); da$Asthma <- factor(da$Asthma, levels = c("Yes", "No")); da$OrParaF <- factor(da$OrParaF, levels = c("Yes", "No"))
x <- run_app("chi2", da, list(rowvar = "Asthma", colvar = "OrParaF")); s <- x$res$stats; rep <- x$report
cells <- list(c("Yes", "Yes", 7777, "68"), c("Yes", "No", 7778, "8"), c("Yes", "Sum", 7779, "76"), c("No", "Yes", 7781, "123"),
              c("No", "No", 7782, "29"), c("No", "Sum", 7783, "152"), c("Sum", "Yes", 7785, "191"), c("Sum", "No", 7786, "37"), c("Sum", "Sum", 7787, "228"))
for (cc in cells) chk("chi2", paste0("40.2 tabla ", cc[1], "/", cc[2]), cc[4], mat_val(x, TC, cc[1], cc[2]), line = as.integer(cc[3]))
dsa <- x$spec$descriptives(da, list(rowvar = "Asthma", colvar = "OrParaF"))
for (it in list(c("Age Mean (SD) Yes", "7.61 (1.98)", 7802), c("Age Mean (SD) No", "7.61 (1.89)", 7803), c("Age Mean (SD) Overall", "7.61 (1.92)", 7804),
                c("Sex F Yes", "30 (39.5 %)", 7811), c("Sex M Yes", "46 (60.5 %)", 7815), c("Asthma Yes Overall", "76 (33.3 %)", 7822)))
  na_chk("chi2", paste0("40.3 table1 ", it[1]), it[2], line = as.integer(it[3]), note = "la app no lo muestra: su table1 solo incluye la variable de respuesta por grupo")
t1y <- t1_row(dsa, "Yes"); t1n <- t1_row(dsa, "No")
chk("chi2", "40.3 table1 OrParaF Yes | Asma Yes (%)", "89.5", t1_nums(t1y[1])[2], line = 7829)
chk("chi2", "40.3 table1 OrParaF Yes | Asma No (%)", "80.9", t1_nums(t1y[2])[2], line = 7830)
chk("chi2", "40.3 table1 OrParaF Yes | Overall n", "191", t1_nums(t1y[3])[1], line = 7831)
chk("chi2", "40.3 table1 OrParaF Yes | Overall (%)", "83.8", t1_nums(t1y[3])[2], line = 7831)
chk("chi2", "40.3 table1 OrParaF No | Asma Yes (%)", "10.5", t1_nums(t1n[1])[2], line = 7833)
chk("chi2", "40.3 table1 OrParaF No | Asma No (%)", "19.1", t1_nums(t1n[2])[2], line = 7834)
chk("chi2", "40.3 table1 OrParaF No | Overall (%)", "16.2", t1_nums(t1n[3])[2], line = 7835)
chk("chi2", "40.3 texto 89.5 %", "89.5", tab_val(x, "Porcentaje por fila", "Yes", "Yes"), line = 7836)
chk("chi2", "40.3 texto 80.9 %", "80.9", tab_val(x, "Porcentaje por fila", "Yes", "No"), line = 7836)
# CrossTable
chk("chi2", "40.4 esperada Yes/Yes", "63.667", mat_val(x, FE, "Yes", "Yes"), line = 7866)
chk("chi2", "40.4 esperada Yes/No", "12.333", mat_val(x, FE, "Yes", "No"), line = 7866)
chk("chi2", "40.4 esperada No/Yes", "127.333", mat_val(x, FE, "No", "Yes"), line = 7871)
chk("chi2", "40.4 esperada No/No", "24.667", mat_val(x, FE, "No", "No"), line = 7871)
PF <- "Porcentaje por fila"; npf <- "la app muestra el % por fila con 1 decimal"
chk("chi2", "40.4 % fila Yes/Yes", "89.474", mat_val(x, PF, "Yes", "Yes"), line = 7867, note = npf)
chk("chi2", "40.4 % fila Yes/No", "10.526", mat_val(x, PF, "Yes", "No"), line = 7867, note = npf)
na_chk("chi2", "40.4 % total fila Yes", "33.333", line = 7867)
chk("chi2", "40.4 % fila No/Yes", "80.921", mat_val(x, PF, "No", "Yes"), line = 7872, note = npf)
chk("chi2", "40.4 % fila No/No", "19.079", mat_val(x, PF, "No", "No"), line = 7872, note = npf)
na_chk("chi2", "40.4 % total fila No", "66.667", line = 7872)
SR <- "sresid"
chk("chi2", "40.4 sresid Yes/Yes", "0.543", mat_val(x, SR, "Yes", "Yes"), line = 7868)
chk("chi2", "40.4 sresid Yes/No", "-1.234", mat_val(x, SR, "Yes", "No"), line = 7868)
chk("chi2", "40.4 sresid No/Yes", "-0.384", mat_val(x, SR, "No", "Yes"), line = 7873)
chk("chi2", "40.4 sresid No/No", "0.873", mat_val(x, SR, "No", "No"), line = 7873)
chk("chi2", "40.4 χ² Pearson", "2.726192", tab_val(x, PR, "Estadístico", 1), line = 7884)
chk("chi2", "40.4 df Pearson", "1", tab_val(x, PR, "gl", 1), line = 7884)
chk("chi2", "40.4 p Pearson", "0.09871439", tab_val(x, PR, "p", 1), line = 7884)
chk("chi2", "40.4 χ² Yates", "2.133366", tab_val(x, PR, "Estadístico", 2), line = 7888)
chk("chi2", "40.4 df Yates", "1", tab_val(x, PR, "gl", 2), line = 7888)
chk("chi2", "40.4 p Yates", "0.1441239", tab_val(x, PR, "p", 2), line = 7888)
OT <- "Odds ratio"
chk("chi2", "40.4 Fisher OR", "1.998463", tab_val(x, OT, "OR"), line = 7893)
chk("chi2", "40.4 Fisher p bilateral", "0.1273169", tab_val(x, PR, "p", 3), line = 7896)
chk("chi2", "40.4 Fisher IC inf", "0.8321903", tab_val(x, OT, "IC 95% inf"), line = 7897)
chk("chi2", "40.4 Fisher IC sup", "5.349523", tab_val(x, OT, "IC 95% sup"), line = 7897)
na_chk("chi2", "40.4 Fisher p unilateral (less)", "0.9705404", line = 7900)
na_chk("chi2", "40.4 Fisher IC sup (less)", "4.577022", line = 7901)
na_chk("chi2", "40.4 Fisher p unilateral (greater)", "0.06923538", line = 7904)
na_chk("chi2", "40.4 Fisher IC inf (greater)", "0.9389371", line = 7905)
chk("chi2", "40.4 mínima esperada", "12.33333", s$min_exp, line = 7909)
chk("chi2", "40.5 mínima esperada (texto)", "12.33", rx(x$assumptions$checks[[1]]$text, paste0("mínima es ", NUM)), line = 7912)
chk("chi2", "40.5 decisión: chi² de Pearson", "TRUE", as.character(isTRUE(x$assumptions$decision$all_met)), line = 7912)
chk("chi2", "40.5 OR (round 4)", "1.9985", s$OR, line = 7970)
chk("chi2", "40.5 IC inf (round 3)", "0.832", s$OR_lo, line = 7972)
chk("chi2", "40.5 IC sup (round 3)", "5.35", s$OR_hi, line = 7972)
chk("chi2", "40.5 p Fisher (round 4)", "0.1273", s$p_fisher, line = 7974)
chk("chi2", "40.5.1 φ", "0.109", tab_val(x, OT, "φ"), line = 7982)
chk("chi2", "40.5.1 φ (texto ≈)", "0.11", s$phi, line = 7983)
PW <- "Poder"
chk("chi2", "40.6 w", "0.109348", tab_val(x, PW, "w"), line = 7993)
chk("chi2", "40.6 N", "228", tab_val(x, PW, "N"), line = 7994)
chk("chi2", "40.6 df", "1", s$df, line = 7995)
chk("chi2", "40.6 poder", "0.378872", tab_val(x, PW, "Poder post hoc"), line = 7997)
chk("chi2", "40.6 poder (texto)", "0.38", s$power, line = 8000)
chk("chi2", "40.6 poder % (texto)", "38", rx(rep, paste0("\\(≈ ", NUM, " %\\)")), line = 8000)
chk("chi2", "40.6 N para 80 % de poder (tabla)", "656.4248", tab_val(x, PW, "N para 80 % de poder"), line = 8008,
    note = "la app muestra ceiling(N) = 657; el valor sin redondear (s$n80) es 656.4248")
chk("chi2", "40.6 N para 80 % de poder (reporte)", "656.4248", rx(rep, paste0("N ≈ ", NUM)), line = 8008)
# Cómo reportar (40.7)
chk("chi2", "40.7 reporte: mínima esperada", "12.33", rx(rep, paste0("esperada fue de ", NUM)), line = 8019)
chk("chi2", "40.7 reporte: gl", "1", rx(rep, "χ²\\[([0-9]+)\\]"), line = 8019)
chk("chi2", "40.7 reporte: χ²", "2.73", rx(rep, paste0("χ²\\[1\\] = ", NUM)), line = 8019)
chk("chi2", "40.7 reporte: p", "0.099", rx(rep, paste0(", p = ", NUM)), line = 8019)
chk("chi2", "40.7 reporte: OR", "2.00", rx(rep, paste0("OR = ", NUM)), line = 8019)
chk("chi2", "40.7 reporte: IC inf", "0.83", rx(rep, paste0("IC 95 % \\[", NUM)), line = 8019)
chk("chi2", "40.7 reporte: IC sup", "5.35", rx(rep, paste0("IC 95 % \\[[^,]+, ", NUM)), line = 8019)
chk("chi2", "40.7 reporte: φ", "0.11", rx(rep, paste0("φ = ", NUM)), line = 8019)
chk("chi2", "40.7 reporte: poder %", "38", rx(rep, paste0("\\(≈ ", NUM, " %\\)")), line = 8019)

# =============================================================================
# CAPÍTULO 41 — Sexo vs. diagnóstico periodontal (r×c, ASR)
# =============================================================================
dp <- book_data("datPerio.txt"); dp$Perio_Dx <- factor(dp$Perio_Dx, levels = c("Ging", "Low_P", "Mod_P", "Adv_P"))
x <- run_app("chi2", dp, list(rowvar = "Sex", colvar = "Perio_Dx")); s <- x$res$stats; rep <- x$report
lv <- c("Ging", "Low_P", "Mod_P", "Adv_P", "Sum")
cnt <- list(F = c(18, 7, 10, 5, 40), M = c(26, 16, 10, 8, 60), Sum = c(44, 23, 20, 13, 100)); cl <- c(F = 8076, M = 8082, Sum = 8088)
for (r in names(cnt)) for (j in 1:5) chk("chi2", paste0("41.2 tabla ", r, "/", lv[j]), as.character(cnt[[r]][j]), mat_val(x, TC, r, lv[j]), line = cl[[r]] + j - 1)
dsp <- x$spec$descriptives(dp, list(rowvar = "Sex", colvar = "Perio_Dx"))
pct1 <- list(Ging = c(45.0, 43.3, 44.0, 8108), Low_P = c(17.5, 26.7, 23.0, 8112), Mod_P = c(25.0, 16.7, 20.0, 8116), Adv_P = c(12.5, 13.3, 13.0, 8120))
for (l in names(pct1)) { rr <- t1_row(dsp, l); for (j in 1:3) chk("chi2", paste0("41.2 table1 ", l, " col ", c("F", "M", "Overall")[j], " (%)"), sprintf("%.1f", pct1[[l]][j]), t1_nums(rr[j])[2], line = pct1[[l]][4] + j - 1) }
ex <- list(F = c("17.600", "9.200", "8.000", "5.200"), M = c("26.400", "13.800", "12.000", "7.800")); el <- c(F = 8147, M = 8155)
for (r in names(ex)) for (j in 1:4) chk("chi2", paste0("41.2 esperada ", r, "/", lv[j]), ex[[r]][j], mat_val(x, FE, r, lv[j]), line = el[[r]])
na_chk("chi2", "41.2 contribución al χ² por celda (F/Ging…)", "0.009", line = 8148)
pf <- list(F = c("45.000", "17.500", "25.000", "12.500"), M = c("43.333", "26.667", "16.667", "13.333")); pl <- c(F = 8149, M = 8157)
for (r in names(pf)) for (j in 1:4) chk("chi2", paste0("41.2 % fila ", r, "/", lv[j]), pf[[r]][j], mat_val(x, PF, r, lv[j]), line = pl[[r]], note = npf)
na_chk("chi2", "41.2 % por columna (F/Ging…)", "40.909", line = 8150)
na_chk("chi2", "41.2 % del total (F/Ging…)", "18.000", line = 8151)
asr_df <- x$posthoc$tables[[1]]
asr <- function(sx, dxv) asr_df$ASR[asr_df$Sex == sx & asr_df$Perio_Dx == dxv]
ab <- list(F = c("0.164", "-1.067", "1.021", "-0.121"), M = c("-0.164", "1.067", "-1.021", "0.121")); al <- c(F = 8152, M = 8160); al2 <- c(F = 8205, M = 8206)
for (r in names(ab)) for (j in 1:4) {
  chk("chi2", paste0("41.2 ASR (CrossTable) ", r, "/", lv[j]), ab[[r]][j], asr(r, lv[j]), line = al[[r]])
  chk("chi2", paste0("41.4 ASR (stdres) ", r, "/", lv[j]), ab[[r]][j], asr(r, lv[j]), line = al2[[r]])
}
chk("chi2", "41.2 χ²", "1.738117", tab_val(x, PR, "Estadístico", 1), line = 8172)
chk("chi2", "41.2 df", "3", tab_val(x, PR, "gl", 1), line = 8172)
chk("chi2", "41.2 p", "0.628493", tab_val(x, PR, "p", 1), line = 8172)
chk("chi2", "41.2 Fisher p", "0.6456615", tab_val(x, PR, "p", 2), line = 8179)
chk("chi2", "41.2 mínima esperada", "5.2", s$min_exp, line = 8182)
chk("chi2", "41.3 V de Cramér", "0.132", tab_val(x, "Cram", "V de Cramér"), line = 8192)
chk("chi2", "41.3 V (texto ≈)", "0.13", rx(rep, paste0("V de Cramér = ", NUM)), line = 8193)
chk("chi2", "41.3 p global (texto)", "0.63", s$p, line = 8193)
chk("chi2", "41.4 nº celdas con |ASR| > 1.96", "0", sum(asr_df$Significativo != "no"), line = 8216)
chk("chi2", "41.4 texto ASR F/Ging", "0.16", asr("F", "Ging"), line = 8240)
chk("chi2", "41.4 texto ASR F/Low_P", "-1.07", asr("F", "Low_P"), line = 8240)
chk("chi2", "41.4 texto ASR F/Mod_P", "1.02", asr("F", "Mod_P"), line = 8240)
chk("chi2", "41.4 texto ASR F/Adv_P", "-0.12", asr("F", "Adv_P"), line = 8240)
chk("chi2", "41.4 reporte: N", "100", rx(rep, paste0("\\(N = ", NUM)), line = 8241)
chk("chi2", "41.4 reporte: menor esperada", "5.2", rx(rep, paste0("la menor es ", NUM)), line = 8241)
chk("chi2", "41.4 reporte: χ²", "1.74", rx(rep, paste0("χ² = ", NUM)), line = 8241)
chk("chi2", "41.4 reporte: gl", "3", rx(rep, paste0("gl = ", NUM)), line = 8241)
chk("chi2", "41.4 reporte: p", "0.63", rx(rep, paste0("gl = [0-9]+, p = ", NUM)), line = 8241)

# =============================================================================
# CAPÍTULO 42 — McNemar (musicoterapia)
# =============================================================================
dm <- ado_get("mcnemar")$example$prep(book_data("datMusicT.txt"))   # mismo reacomodo que el libro (merge por Px)
x <- run_app("mcnemar", dm, list(measures = c("Music", "Control"))); s <- x$res$stats; rep <- x$report
TP <- "Tabla apareada"
mt <- list(c("Colaborador", "Colaborador", 8339, "23"), c("Colaborador", "No colaborador", 8340, "31"), c("Colaborador", "Sum", 8341, "54"),
           c("No colaborador", "Colaborador", 8343, "0"), c("No colaborador", "No colaborador", 8344, "5"), c("No colaborador", "Sum", 8345, "5"),
           c("Sum", "Colaborador", 8347, "23"), c("Sum", "No colaborador", 8348, "36"), c("Sum", "Sum", 8349, "59"))
for (cc in mt) chk("mcnemar", paste0("42.2.2 tabla ", cc[1], "/", cc[2]), cc[4], mat_val(x, TP, cc[1], cc[2]), line = as.integer(cc[3]))
chk("mcnemar", "42.2.2 texto b = 31", "31", s$b, line = 8350)
chk("mcnemar", "42.2.2 texto c = 0", "0", s$c, line = 8350)
chk("mcnemar", "42.2.3 χ² McNemar", "31", tab_val(x, "^McNemar$", "χ²"), line = 8360)
chk("mcnemar", "42.2.3 df", "1", tab_val(x, "^McNemar$", "gl"), line = 8360)
chk("mcnemar", "42.2.3 p", "2.58e-08", tab_val(x, "^McNemar$", "p"), line = 8360)
chk("mcnemar", "42.2.4 b", "31", tab_val(x, "^McNemar$", "Discordantes b"), line = 8367)
chk("mcnemar", "42.2.4 c", "0", tab_val(x, "^McNemar$", "Discordantes c"), line = 8367)
chk("mcnemar", "42.2.4 OR b/c", "infinito", if (is.infinite(tab_val(x, "^McNemar$", "OR condicional (b/c)"))) "infinito" else as.character(s$OR), line = 8369)
chk("mcnemar", "42.2.4 OR con corrección +0.5", "63", tab_val(x, "^McNemar$", "OR con corrección de continuidad (+0.5)"), line = 8372)
chk("mcnemar", "42.4 reporte: N", "59", rx(rep, paste0("conjunto de ", NUM)), line = 8389)
chk("mcnemar", "42.4 reporte: n colaborador bajo música", "54", rx(rep, paste0("fue de ", NUM, "/")), line = 8389)
chk("mcnemar", "42.4 reporte: % bajo música", "91.5", rx(rep, paste0("fue de [0-9]+/[0-9]+ \\(", NUM)), line = 8389)
chk("mcnemar", "42.4 reporte: n colaborador bajo control", "23", rx(rep, paste0("y de ", NUM, "/")), line = 8389)
chk("mcnemar", "42.4 reporte: % bajo control", "39.0", rx(rep, paste0("y de [0-9]+/[0-9]+ \\(", NUM)), line = 8389)
chk("mcnemar", "42.4 reporte: χ²", "31.0", rx(rep, paste0("χ²\\[1\\] = ", NUM)), line = 8389)
chk("mcnemar", "42.4 reporte: p < 0.001", "<0.001", s$p, line = 8389)
chk("mcnemar", "42.4 reporte: p < 0.001 (texto)", "p < 0.001", paste0("p ", sub("^.*χ²\\[1\\] = [0-9.]+, p (< 0\\.001|= [0-9.]+).*$", "\\1", rep)), line = 8389)
chk("mcnemar", "42.4 reporte: 31 solo con musicoterapia", "31", rx(rep, paste0("sujetos, ", NUM, " cambiaron")), line = 8389)
chk("mcnemar", "42.4 reporte: OR corregido ~63", "63", rx(rep, paste0("se aproxima a ", NUM)), line = 8389)

# =============================================================================
# CAPÍTULO 43 — ejemplos de la paradoja de prevalencia (κ de tablas 2×2)
# =============================================================================
for (e in list(list(a = 1, b = 4, c = 6, d = 89, k = "0.115", p0 = "0.900", pe = "0.887", l = 8460),
               list(a = 89, b = 2, c = 8, d = 1, k = "0.127", p0 = "0.900", pe = "0.885", l = 8465),
               list(a = 55, b = 9, c = 1, d = 35, k = "0.793", p0 = "0.900", pe = "0.517", l = 8470))) {
  xk <- run_app("kappa", mk_k(e$a, e$b, e$c, e$d), list(raters = c("r1", "r2")))
  tag <- sprintf("43.3 a=%d b=%d c=%d d=%d", e$a, e$b, e$c, e$d)
  na_chk("kappa", paste(tag, "P0"), e$p0, line = e$l)
  na_chk("kappa", paste(tag, "Pe"), e$pe, line = e$l)
  chk("kappa", paste(tag, "κ"), e$k, xk$res$stats$k_unw, line = e$l)
}

# =============================================================================
# CAPÍTULO 44 — Kappa (PGMax: clasificación de terceros molares)
# =============================================================================
pg <- book_data("PGMax.txt")
x <- run_app("kappa", pg, list(raters = c("Obs2", "Obs3"))); s <- x$res$stats; rep <- x$report
KT <- "Kappa de Cohen"
chk("kappa", "44.2 n sujetos", "20", s$n, line = 8566)
chk("kappa", "44.2 κ no ponderado", "0.0291", tab_val(x, KT, "κ", 1), line = 8568)
chk("kappa", "44.2 z no ponderado", "0.277", tab_val(x, KT, "z", 1), line = 8570)
chk("kappa", "44.2 p no ponderado", "0.781", tab_val(x, KT, "p", 1), line = 8571)
chk("kappa", "44.3 κ lineal", "0.251", tab_val(x, KT, "κ", 2), line = 8586)
chk("kappa", "44.3 z lineal", "1.57", tab_val(x, KT, "z", 2), line = 8588)
chk("kappa", "44.3 p lineal", "0.115", tab_val(x, KT, "p", 2), line = 8589)
chk("kappa", "44.3 κ cuadrático", "0.38", tab_val(x, KT, "κ", 3), line = 8597)
chk("kappa", "44.3 z cuadrático", "1.75", tab_val(x, KT, "z", 3), line = 8599)
chk("kappa", "44.3 p cuadrático", "0.0793", tab_val(x, KT, "p", 3), line = 8600)
xn <- run_app("kappa", pg, list(raters = c("Obs2", "Obs3")), list(escala = "nominal"))
chk("kappa", "44.2 opción nominal: κ en reporte", "0.03", rx(xn$report, paste0("κ = ", NUM)), line = 8568)
# Cómo reportar (44.6)
chk("kappa", "44.6 reporte: κ no ponderado", "0.03", rx(rep, paste0("no ponderado fue κ = ", NUM)), line = 8657)
chk("kappa", "44.6 reporte: z no ponderado", "0.28", rx(rep, paste0("no ponderado fue κ = [-0-9.]+ \\(z = ", NUM)), line = 8657)
chk("kappa", "44.6 reporte: p no ponderado", "0.781", rx(rep, paste0("no ponderado fue κ = [-0-9.]+ \\(z = [-0-9.]+; p = ", NUM)), line = 8657)
chk("kappa", "44.6 reporte: κ lineal", "0.25", rx(rep, paste0("lineal \\(κ = ", NUM)), line = 8657)
chk("kappa", "44.6 reporte: z lineal", "1.57", rx(rep, paste0("lineal \\(κ = [-0-9.]+; z = ", NUM)), line = 8657)
chk("kappa", "44.6 reporte: p lineal", "0.115", rx(rep, paste0("lineal \\(κ = [-0-9.]+; z = [-0-9.]+; p = ", NUM)), line = 8657)
chk("kappa", "44.6 reporte: κ cuadrático", "0.38", rx(rep, paste0("cuadrático \\(κ = ", NUM)), line = 8657)
chk("kappa", "44.6 reporte: z cuadrático", "1.75", rx(rep, paste0("cuadrático \\(κ = [-0-9.]+; z = ", NUM)), line = 8657)
chk("kappa", "44.6 reporte: p cuadrático", "0.079", rx(rep, paste0("cuadrático \\(κ = [-0-9.]+; z = [-0-9.]+; p = ", NUM)), line = 8657)
chk("kappa", "44.6 reporte: interpretación ponderados", "razonable", landis(s$k_sq), line = 8657, note = "libro: 'acuerdo aceptable (escala de Landis y Koch)'; la app usa la etiqueta de Landis y Koch")

# ---- 44.4 múltiples evaluadores (10) -----------------------------------------
x <- run_app("kappa", pg, list(raters = names(pg))); s <- x$res$stats; rep <- x$report
MT <- "Múltiples evaluadores"
chk("kappa", "44.4 Fleiss: evaluadores", "10", s$nr, line = 8609)
chk("kappa", "44.4 Fleiss: sujetos", "20", s$n_fleiss, line = 8608)
chk("kappa", "44.4 Fleiss κ", "0.0486", tab_val(x, MT, "Valor", 1), line = 8610)
chk("kappa", "44.4 Fleiss z", "3.13", tab_val(x, MT, "z", 1), line = 8612)
chk("kappa", "44.4 Fleiss p", "0.00173", tab_val(x, MT, "p", 1), line = 8613)
chk("kappa", "44.4 Fleiss κ (texto ≈)", "0.05", rx(rep, paste0("Fleiss \\(no ponderado\\) fue κ = ", NUM)), line = 8615)
ng <- "el libro pasa categ.labels = c('AI','AII','BII','BIII') (4 categorías); la app usa todas las categorías observadas (9) y no permite restringirlas"
na_chk("kappa", "44.4 AC1 pa", "0.3468651", line = 8623)
na_chk("kappa", "44.4 AC1 pe", "0.2374243", line = 8623)
chk("kappa", "44.4 AC1 de Gwet", "0.14351", tab_val(x, MT, "Valor", 2), line = 8623, note = ng)
na_chk("kappa", "44.4 AC1 EE", "0.03421", line = 8623)
chk("kappa", "44.4 AC1 IC inf", "0.072", tab_val(x, MT, "IC 95% inf", 2), line = 8623, note = ng)
chk("kappa", "44.4 AC1 IC sup", "0.215", tab_val(x, MT, "IC 95% sup", 2), line = 8623, note = ng)
chk("kappa", "44.4 AC1 p", "0.0002453858", tab_val(x, MT, "p", 2), line = 8623, note = ng)
na_chk("kappa", "44.4 AC1 tot.obs", "200", line = 8624)
na_chk("kappa", "44.4 AC2 pa", "0.8461199", line = 8630)
na_chk("kappa", "44.4 AC2 pe", "0.6858924", line = 8630)
chk("kappa", "44.4 AC2 de Gwet (cuadrático)", "0.5101", tab_val(x, MT, "Valor", 3), line = 8630, note = ng)
na_chk("kappa", "44.4 AC2 EE", "0.08978", line = 8630)
chk("kappa", "44.4 AC2 IC inf", "0.322", tab_val(x, MT, "IC 95% inf", 3), line = 8630, note = ng)
chk("kappa", "44.4 AC2 IC sup", "0.698", tab_val(x, MT, "IC 95% sup", 3), line = 8630, note = ng)
chk("kappa", "44.4 AC2 p", "8.875079e-06", tab_val(x, MT, "p", 3), line = 8630, note = ng)
chk("kappa", "44.4 texto: Fleiss ≈ 0.05", "0.05", s$k_fleiss, line = 8634)
chk("kappa", "44.4 texto/reporte: AC1 ≈ 0.14", "0.14", rx(rep, paste0("prevalencia, fue ", NUM)), line = 8635, note = ng)
chk("kappa", "44.4 texto/reporte: AC2 ≈ 0.51", "0.51", rx(rep, paste0("cuadráticos, ", NUM)), line = 8636, note = ng)
chk("kappa", "44.4 texto/reporte: AC2 IC inf", "0.32", rx(rep, paste0("cuadráticos, [-0-9.]+ \\(IC 95 %: ", NUM)), line = 8636, note = ng)
chk("kappa", "44.4 texto/reporte: AC2 IC sup", "0.70", rx(rep, paste0("cuadráticos, [-0-9.]+ \\(IC 95 %: [-0-9.]+–", NUM)), line = 8636, note = ng)
chk("kappa", "44.4 texto: AC2 p < 0.001", "<0.001", s$ac2_p, line = 8636)
chk("kappa", "44.4 texto: AC2 'acuerdo moderado'", "moderado", landis(s$ac2), line = 8636)
xl <- run_app("kappa", pg, list(raters = names(pg)), list(levels = "AI, AII, BII, BIII"))
chk("kappa", "44.4 AC2 con orden 'AI, AII, BII, BIII' indicado en la app", "0.5101", tab_val(xl, MT, "Valor", 3), line = 8630,
    note = "la app descarta el orden indicado porque no cubre todas las categorías observadas (concordancia.R:22-23)")

# =============================================================================
# CAPÍTULO 45 — ICC (irr::anxiety)
# =============================================================================
data("anxiety", package = "irr", envir = environment())
x <- run_app("icc", anxiety, list(raters = c("rater1", "rater2", "rater3"))); s <- x$res$stats; rep <- x$report
vt <- x$assumptions$checks[[1]]$text
chk("icc", "45.2 varianza rater1", "2.1342", rx(vt, paste0("rater1 = ", NUM)), line = 8699)
chk("icc", "45.2 varianza rater2", "1.6079", rx(vt, paste0("rater2 = ", NUM)), line = 8701)
chk("icc", "45.2 varianza rater3", "1.8", rx(vt, paste0("rater3 = ", NUM)), line = 8703)
IT <- "irr::icc"
chk("icc", "45.3 sujetos", "20", s$n, line = 8719)
chk("icc", "45.3 evaluadores", "3", s$k, line = 8720)
chk("icc", "45.3 ICC(C,1)", "0.216", tab_val(x, IT, "ICC", 1), line = 8721)
chk("icc", "45.3 ICC(C,1) gl1", "19", tab_val(x, IT, "gl1", 1), line = 8724)
chk("icc", "45.3 ICC(C,1) gl2", "38", tab_val(x, IT, "gl2", 1), line = 8724)
chk("icc", "45.3 ICC(C,1) F", "1.83", tab_val(x, IT, "F", 1), line = 8724)
chk("icc", "45.3 ICC(C,1) p", "0.0562", tab_val(x, IT, "p", 1), line = 8724)
chk("icc", "45.3 ICC(C,1) IC inf", "-0.046", tab_val(x, IT, "IC 95% inf", 1), line = 8727)
chk("icc", "45.3 ICC(C,1) IC sup", "0.522", tab_val(x, IT, "IC 95% sup", 1), line = 8727)
chk("icc", "45.3 ICC(A,1)", "0.198", tab_val(x, IT, "ICC", 2), line = 8738)
chk("icc", "45.3 ICC(A,1) gl1", "19", tab_val(x, IT, "gl1", 2), line = 8741)
nv <- "irr 0.85 instalado da F(19,38), p = 0.0562 para el ICC(A,1); el libro se generó con otra versión de irr (gl2 de McGraw-Wong = 39.7)"
chk("icc", "45.3 ICC(A,1) gl2", "39.7", tab_val(x, IT, "gl2", 2), line = 8741, note = nv)
chk("icc", "45.3 ICC(A,1) F", "1.83", tab_val(x, IT, "F", 2), line = 8741)
chk("icc", "45.3 ICC(A,1) p", "0.0543", tab_val(x, IT, "p", 2), line = 8741, note = nv)
chk("icc", "45.3 ICC(A,1) IC inf", "-0.039", tab_val(x, IT, "IC 95% inf", 2), line = 8744)
chk("icc", "45.3 ICC(A,1) IC sup", "0.494", tab_val(x, IT, "IC 95% sup", 2), line = 8744)
PS <- "psych::ICC"
psb <- list(list("ICC1", "0.18", "1.6", "19", "40", "0.094", "-0.077", "0.48", 8751), list("ICC2", "0.20", "1.8", "19", "38", "0.056", "-0.039", "0.49", 8752),
            list("ICC3", "0.22", "1.8", "19", "38", "0.056", "-0.046", "0.52", 8753), list("ICC1k", "0.39", "1.6", "19", "40", "0.094", "-0.275", "0.74", 8754),
            list("ICC2k", "0.43", "1.8", "19", "38", "0.056", "-0.127", "0.75", 8755), list("ICC3k", "0.45", "1.8", "19", "38", "0.056", "-0.153", "0.77", 8756))
cols <- c("ICC", "F", "gl1", "gl2", "p", "IC 95% inf", "IC 95% sup")
for (b in psb) for (j in seq_along(cols)) chk("icc", paste0("45.3 psych ", b[[1]], " ", cols[j]), b[[j + 1]], tab_val(x, PS, cols[j], b[[1]]), line = b[[9]])
chk("icc", "45.3 texto: ICC(3,1) = 0.216", "0.216", s$icc_c, line = 8778)
chk("icc", "45.3 texto: ICC(2,1) = 0.198", "0.198", s$icc_a, line = 8778)
# opción unidad = promedio (equivalente a ICC2k / ICC3k de psych)
xa <- run_app("icc", anxiety, list(raters = c("rater1", "rater2", "rater3")), list(unit = "average"))
na_ <- "opción unit = 'average' (irr::icc) frente a ICC3k/ICC2k de psych impresos en el libro"
chk("icc", "45.3 average: ICC(3,k)", "0.45", tab_val(xa, IT, "ICC", 1), line = 8756, note = na_)
chk("icc", "45.3 average: ICC(3,k) IC inf", "-0.153", tab_val(xa, IT, "IC 95% inf", 1), line = 8756, note = na_)
chk("icc", "45.3 average: ICC(3,k) IC sup", "0.77", tab_val(xa, IT, "IC 95% sup", 1), line = 8756, note = na_)
chk("icc", "45.3 average: ICC(2,k)", "0.43", tab_val(xa, IT, "ICC", 2), line = 8755, note = na_)
chk("icc", "45.3 average: ICC(2,k) IC inf", "-0.127", tab_val(xa, IT, "IC 95% inf", 2), line = 8755, note = na_)
chk("icc", "45.3 average: ICC(2,k) IC sup", "0.75", tab_val(xa, IT, "IC 95% sup", 2), line = 8755, note = na_)
# Cómo reportar (45.5)
chk("icc", "45.5 reporte: ICC(3,1)", "0.216", rx(rep, paste0("que fue de ", NUM)), line = 8801)
chk("icc", "45.5 reporte: IC inf", "-0.05", rx(rep, paste0("que fue de [-0-9.]+ \\(IC 95 %: ", NUM)), line = 8801)
chk("icc", "45.5 reporte: IC sup", "0.52", rx(rep, paste0("que fue de [-0-9.]+ \\(IC 95 %: [-0-9.]+ a ", NUM)), line = 8801)
chk("icc", "45.5 reporte: p", "0.056", rx(rep, paste0("que fue de [-0-9.]+ \\(IC 95 %: [-0-9.]+ a [-0-9.]+; p = ", NUM)), line = 8801)
chk("icc", "45.5 reporte: umbral pobre (< 0.50)", "pobre", koo_li(s$icc_c), line = 8801)
chk("icc", "45.5 reporte: ICC(2,1)", "0.198", rx(rep, paste0("ICC\\(2,1\\)— fue de ", NUM)), line = 8801)

# =============================================================================
# CAPÍTULOS 46–47 — radiómetros simulados (set.seed(123), n = 40)
# =============================================================================
set.seed(123); n <- 40
Gold <- round(rnorm(n, mean = 1000, sd = 150))
A    <- round(Gold + rnorm(n, mean = 15, sd = 30))
B    <- round(Gold + (rexp(n, rate = 1 / 40) - 40))
lamp <- data.frame(Gold, A, B)
x <- run_app("bland_altman", lamp, list(method = "A", reference = "Gold")); s <- x$res$stats; rep <- x$report; sw <- x$assumptions$checks[[1]]
BT <- "^Bland-Altman$"
chk("bland_altman", "46.3 Shapiro A−Gold: W", "0.98726", sw$W, line = 8858)
chk("bland_altman", "46.3 Shapiro A−Gold: p", "0.9256", sw$p, line = 8858)
chk("bland_altman", "46.3 texto: p = 0.93", "0.93", sw$p, line = 8866)
chk("bland_altman", "46.4.1 A: sesgo", "14.75", tab_val(x, BT, "Sesgo"), line = 8911)
chk("bland_altman", "46.4.1 A: sesgo %", "1.5", tab_val(x, BT, "Sesgo (% nivel medio)"), line = 8911)
chk("bland_altman", "46.4.1 A: LoA inf", "-41.7", tab_val(x, BT, "LoA inferior"), line = 8911)
chk("bland_altman", "46.4.1 A: LoA inf %", "-4.1", tab_val(x, BT, "LoA inf (%)"), line = 8911)
chk("bland_altman", "46.4.1 A: LoA sup", "71.2", tab_val(x, BT, "LoA superior"), line = 8911)
chk("bland_altman", "46.4.1 A: LoA sup %", "7.0", tab_val(x, BT, "LoA sup (%)"), line = 8911)
chk("bland_altman", "46.4.1 texto: sesgo 14.75", "14.75", s$bias, line = 8913)
chk("bland_altman", "46.4.1 texto: sesgo 1.5 %", "1.5", s$pct_b, line = 8913)
chk("bland_altman", "46.4.1 texto: LoA −41.7", "-41.7", s$loa_l, line = 8913)
chk("bland_altman", "46.4.1 texto: LoA 71.2", "71.2", s$loa_u, line = 8913)
chk("bland_altman", "46.4.1 texto: −4.1 %", "-4.1", s$pct_l, line = 8913)
chk("bland_altman", "46.4.1 texto: 7.0 %", "7.0", s$pct_u, line = 8913)
chk("bland_altman", "46.5 IC 95% sesgo inf", "5.54", tab_val(x, BT, "IC 95% sesgo inf"), line = 8934)
chk("bland_altman", "46.5 IC 95% sesgo sup", "23.96", tab_val(x, BT, "IC 95% sesgo sup"), line = 8934)
chk("bland_altman", "46.5 texto IC inf", "5.54", s$ci_lo, line = 8935)
chk("bland_altman", "46.5 texto IC sup", "23.96", s$ci_hi, line = 8935)
CC <- "CCC de Lin"
chk("bland_altman", "46.6 CCC A vs Gold", "0.9721226", tab_val(x, CC, "ρc"), line = 8943)
chk("bland_altman", "46.6 CCC A IC inf", "0.9487943", tab_val(x, CC, "IC 95% inf"), line = 8943)
chk("bland_altman", "46.6 CCC A IC sup", "0.9849053", tab_val(x, CC, "IC 95% sup"), line = 8943)
chk("bland_altman", "46.6 texto: CCC A ≈ 0.97", "0.97", s$ccc, line = 8949)
# Cómo reportar (46.7)
chk("bland_altman", "46.7 reporte: Shapiro p", "0.93", rx(rep, paste0("Shapiro-Wilk, p = ", NUM)), line = 8953)
chk("bland_altman", "46.7 reporte: sesgo", "14.75", rx(rep, paste0("fue de ", NUM)), line = 8953)
chk("bland_altman", "46.7 reporte: IC inf", "5.54", rx(rep, paste0("IC 95 %: ", NUM)), line = 8953)
chk("bland_altman", "46.7 reporte: IC sup", "23.96", rx(rep, paste0("IC 95 %: [-0-9.]+ a ", NUM)), line = 8953)
chk("bland_altman", "46.7 reporte: LoA inf", "-41.7", rx(rep, paste0("fueron de ", NUM)), line = 8953)
chk("bland_altman", "46.7 reporte: LoA sup", "71.2", rx(rep, paste0("fueron de [-0-9.]+ a ", NUM)), line = 8953)
chk("bland_altman", "46.7 reporte: LoA inf %", "-4.1", rx(rep, paste0("\\(de ", NUM, " %")), line = 8953)
chk("bland_altman", "46.7 reporte: LoA sup %", "7.0", rx(rep, paste0("\\(de [-0-9.]+ % a ", NUM, " %")), line = 8953)
xB <- run_app("bland_altman", lamp, list(method = "B", reference = "Gold")); sB <- xB$res$stats; swB <- xB$assumptions$checks[[1]]
chk("bland_altman", "46.3 Shapiro B−Gold: W", "0.88752", swB$W, line = 8865)
chk("bland_altman", "46.3 Shapiro B−Gold: p", "0.0008452", swB$p, line = 8865)
chk("bland_altman", "46.3 texto: p < 0.001 (B)", "<0.001", swB$p, line = 8866)
chk("bland_altman", "46.3 B: decisión → CCC (diferencias no normales)", "FALSE", as.character(isTRUE(xB$assumptions$decision$all_met)), line = 8866)
chk("bland_altman", "46.4.2 B: sesgo", "-6.80", tab_val(xB, BT, "Sesgo"), line = 8917)
chk("bland_altman", "46.4.2 B: sesgo %", "-0.7", tab_val(xB, BT, "Sesgo (% nivel medio)"), line = 8917)
chk("bland_altman", "46.4.2 B: LoA inf", "-61.4", tab_val(xB, BT, "LoA inferior"), line = 8917)
chk("bland_altman", "46.4.2 B: LoA inf %", "-6.1", tab_val(xB, BT, "LoA inf (%)"), line = 8917)
chk("bland_altman", "46.4.2 B: LoA sup", "47.8", tab_val(xB, BT, "LoA superior"), line = 8917)
chk("bland_altman", "46.4.2 B: LoA sup %", "4.8", tab_val(xB, BT, "LoA sup (%)"), line = 8917)
chk("bland_altman", "46.4.2 texto: LoA −61.4", "-61.4", sB$loa_l, line = 8919)
chk("bland_altman", "46.4.2 texto: LoA 47.8", "47.8", sB$loa_u, line = 8919)
chk("bland_altman", "46.4.2 texto: −6.1 %", "-6.1", sB$pct_l, line = 8919)
chk("bland_altman", "46.4.2 texto: 4.8 %", "4.8", sB$pct_u, line = 8919)
chk("bland_altman", "46.6 CCC B vs Gold", "0.9775555", tab_val(xB, CC, "ρc"), line = 8948)
chk("bland_altman", "46.6 CCC B IC inf", "0.9581052", tab_val(xB, CC, "IC 95% inf"), line = 8948)
chk("bland_altman", "46.6 CCC B IC sup", "0.9880309", tab_val(xB, CC, "IC 95% sup"), line = 8948)
chk("bland_altman", "46.6 texto: CCC B ≈ 0.98", "0.98", sB$ccc, line = 8949)
chk("bland_altman", "46.7 reporte B: Shapiro p < 0.001", "p < 0.001", paste0("p ", sub("^.*Shapiro-Wilk, p (< 0\\.001|= [0-9.]+).*$", "\\1", xB$report)), line = 8953)

# ---- 47 CCC de Lin: B vs. Gold ------------------------------------------------
x <- run_app("ccc", lamp, list(method = "B", reference = "Gold")); s <- x$res$stats; rep <- x$report; sw <- x$assumptions$checks[[1]]
LT <- "epi.ccc"
chk("ccc", "47.3 Shapiro B−Gold: W", "0.88752", sw$W, line = 9012)
chk("ccc", "47.3 Shapiro B−Gold: p", "0.0008452", sw$p, line = 9012)
chk("ccc", "47.5 ρc", "0.9775555", tab_val(x, LT, "ρc"), line = 9033)
chk("ccc", "47.5 ρc IC inf", "0.9581052", tab_val(x, LT, "IC 95% inf"), line = 9033)
chk("ccc", "47.5 ρc IC sup", "0.9880309", tab_val(x, LT, "IC 95% sup"), line = 9033)
na_chk("ccc", "47.5 s.shift", "0.9896866", line = 9036)
na_chk("ccc", "47.5 l.shift", "0.05086971", line = 9039)
chk("ccc", "47.5 C.b", "0.9986542", tab_val(x, LT, "Cb (corrección de sesgo)"), line = 9042)
na_chk("ccc", "47.5 $blalt (40 pares media/diferencia; p. ej. fila 1: 896.0, −40)", "896.0", line = 9046, note = "la app no lo muestra como tabla (solo en el gráfico de Bland-Altman)")
chk("ccc", "47.5 sblalt est", "-6.8", tab_val(x, LT, "Sesgo"), line = 9089)
na_chk("ccc", "47.5 sblalt delta.sd", "27.8634", line = 9089, note = "la app no lo muestra en el CCC (sí aparece como DE en la tabla de Bland-Altman B vs Gold)")
chk("ccc", "47.5 sblalt delta.sd (vía Bland-Altman B)", "27.8634", tab_val(xB, BT, "DE"), line = 9089, note = "verificado en la entrada Bland-Altman")
chk("ccc", "47.5 sblalt lower", "-61.41127", tab_val(x, LT, "LoA inf"), line = 9089)
chk("ccc", "47.5 sblalt upper", "47.81127", tab_val(x, LT, "LoA sup"), line = 9089)
na_chk("ccc", "47.5 nmissing", "0", line = 9092)
chk("ccc", "47.5 CCC (round 4)", "0.9776", s$ccc, line = 9096)
chk("ccc", "47.5 IC inf (round 3)", "0.958", s$lo, line = 9099)
chk("ccc", "47.5 IC sup (round 3)", "0.988", s$hi, line = 9099)
chk("ccc", "47.5 r de Pearson", "0.9789", tab_val(x, LT, "r de Pearson"), line = 9101)
chk("ccc", "47.5 Cb (round 4)", "0.9987", s$Cb, line = 9103)
chk("ccc", "47.5 sesgo medio", "-6.8", s$bias, line = 9105)
chk("ccc", "47.5 LoA inferior (round 3)", "-61.411", s$loa_l, line = 9107)
chk("ccc", "47.5 LoA superior (round 3)", "47.811", s$loa_u, line = 9109)
chk("ccc", "47.6 interpretación (Tabla 47.1)", "sustancial", tab_val(x, LT, "Interpretación"), line = 9128)
# Cómo reportar (47.7)
chk("ccc", "47.7 reporte: ρc", "0.978", rx(rep, paste0("ρc = ", NUM)), line = 9128)
chk("ccc", "47.7 reporte: IC inf", "0.958", rx(rep, paste0("ρc = [-0-9.]+ \\(IC 95 %: ", NUM)), line = 9128)
chk("ccc", "47.7 reporte: IC sup", "0.988", rx(rep, paste0("ρc = [-0-9.]+ \\(IC 95 %: [-0-9.]+ a ", NUM)), line = 9128)
chk("ccc", "47.7 reporte: r", "0.979", rx(rep, paste0("Pearson fue r = ", NUM)), line = 9128)
chk("ccc", "47.7 reporte: Cb", "0.999", rx(rep, paste0("Cb = ", NUM)), line = 9128)
xA <- run_app("ccc", lamp, list(method = "A", reference = "Gold"))
chk("ccc", "46.6 CCC A vs Gold (entrada CCC)", "0.9721226", tab_val(xA, LT, "ρc"), line = 8943)
chk("ccc", "46.6 CCC A IC inf (entrada CCC)", "0.9487943", tab_val(xA, LT, "IC 95% inf"), line = 8943)
chk("ccc", "46.6 CCC A IC sup (entrada CCC)", "0.9849053", tab_val(xA, LT, "IC 95% sup"), line = 8943)

# =============================================================================
out <- write_results(file.path(Sys.getenv("VERIFY_DIR"), "results_cat.csv"))
noapp <- grepl("^la app no lo muestra", out$note) & out$app == "NA / no encontrado"
cat("\nTotal:", nrow(out), "| coinciden:", sum(out$match), "| no coinciden:", sum(!out$match & !noapp), "| la app no lo muestra:", sum(noapp), "\n\n")
mm <- out[!out$match & !noapp, ]
if (nrow(mm)) for (i in seq_len(nrow(mm))) cat(sprintf("[%s] %s | libro %s (l. %s) | app %s%s\n", mm$analysis[i], mm$item[i], mm$book[i], mm$line[i], mm$app[i],
                                                        if (nzchar(mm$note[i])) paste0(" | ", mm$note[i]) else ""))
