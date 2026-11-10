# =============================================================================
# verify_t.R — Verifica que la app reproduce las cifras del libro en:
#   cap. 9 (supuestos, DXKTw X8 ~ Grupo), caps. 14–15 (t independiente, VickersF Dual),
#   cap. 16 (t apareada, Carbuse), caps. 29–30 (Wilcoxon suma de rangos, AUCFil),
#   caps. 31–32 (Wilcoxon rangos con signo, DXKTw grupo 2, X0 vs X2).
#   (Cap. 28 solo contiene un ejemplo didáctico de rangos con un vector inventado: no
#    corresponde a ninguna salida de la app, por lo que no se verifica.)
# Solo se usan salidas de la app (stats, tablas, supuestos, reporte y table1 de descriptivos).
# =============================================================================
source(file.path(Sys.getenv("VERIFY_DIR"), "verify_lib.R"))

BL <- readLines(BOOK_TXT, encoding = "UTF-8", warn = FALSE)
# Línea del libro: primera coincidencia de `pat` a partir de `from` (n-ésima si nth > 1)
L <- function(pat, from, nth = 1, fixed = TRUE) {
  k <- which(grepl(pat, BL[from:length(BL)], fixed = fixed))
  if (length(k) < nth) { warning("Línea no encontrada: ", pat); return(NA_integer_) }
  as.integer(k[nth] + from - 1)
}
# chk con el texto del libro normalizado (signo menos tipográfico U+2212 → "-")
ck <- function(analysis, item, book, app, line, note = "", tol = NULL)
  chk(analysis, item, gsub("−", "-", book), app, line = line, note = note, tol = tol)
NOAPP <- "la app no lo muestra"
# Extrae del reporte el grupo `k` de la regex (numérico) o como texto (rxs)
rx  <- function(txt, pat, k = 1) { m <- regmatches(txt, regexec(pat, txt, perl = TRUE))[[1]]; if (length(m) < k + 1) NA else as.numeric(m[k + 1]) }
rxs <- function(txt, pat, k = 1) { m <- regmatches(txt, regexec(pat, txt, perl = TRUE))[[1]]; if (length(m) < k + 1) NA else m[k + 1] }
NUM <- "(-?[0-9]+(?:\\.[0-9]+)?)"
# table1 que muestra la app (pestaña Descriptivos): celda → números
T1 <- function(spec, d, v) { t <- as.data.frame(spec$descriptives(d, v)); names(t)[1] <- "fila"; t$fila <- trimws(t$fila); t }
t1n <- function(t, fila, col, k) { s <- t[[col]][t$fila == fila][1]; as.numeric(regmatches(s, gregexpr("-?[0-9]+(\\.[0-9]+)?", s))[[1]][k]) }
TOL_T1 <- "table1 de la app (3 cifras significativas); se compara dentro de su precisión"
sw_chk <- function(x) Find(function(c) identical(c$kind, "normalidad"), x$assumptions$checks)
lv_chk <- function(x) Find(function(c) identical(c$kind, "homogeneidad"), x$assumptions$checks)

# =============================================================================
# Cap. 9 — Supuestos: DXKTw, X8 ~ Grupo (la app: t independiente con esas variables)
# =============================================================================
A <- "t_indep (cap. 9, DXKTw X8)"
d9 <- book_data("DXKTw.txt"); d9$Grupo <- factor(d9$Grupo); levels(d9$Grupo) <- c("Dexketoprofeno", "Control")
v9 <- list(response = "X8", group = "Grupo")
x9 <- run_app("t_indep", d9, v9); t9 <- T1(x9$spec, d9, v9)
TD <- "Estad.stica descriptiva por grupo"
for (g in c("Dexketoprofeno", "Control")) {
  l0 <- L(paste0("^", g, "$"), 1700, fixed = FALSE)
  vals <- BL[l0 + 1:6]
  ck(A, paste(g, "n"),       vals[1], tab_val(x9, TD, "n", g), l0 + 1)
  ck(A, paste(g, "Media"),   vals[2], tab_val(x9, TD, "Media", g), l0 + 2)
  ck(A, paste(g, "DE"),      vals[3], tab_val(x9, TD, "DE", g), l0 + 3)
  ck(A, paste(g, "EE"),      vals[4], tab_val(x9, TD, "EE", g), l0 + 4, NOAPP)
  ck(A, paste(g, "Mediana"), vals[5], tab_val(x9, TD, "Mediana", g), l0 + 5)
  ck(A, paste(g, "IQR"),     vals[6], tab_val(x9, TD, "IQR", g), l0 + 6)
}
s9 <- sw_chk(x9); l9 <- lv_chk(x9)
ls <- L("     30.0000      35.7667", 1700); lk <- ls + 2; lp <- ls + 4
ck(A, "stat.desc median",       "30.0000",  t1n(t9, "Median [Min, Max]", "Overall", 1), ls, TOL_T1, tol = 0.05)
ck(A, "stat.desc mean",         "35.7667",  t1n(t9, "Mean (SD)", "Overall", 1), ls, TOL_T1, tol = 0.05)
ck(A, "stat.desc SE.mean",      "3.4418",   NA, ls, NOAPP)
ck(A, "stat.desc CI.mean.0.95", "6.8870",   NA, ls, NOAPP)
ck(A, "stat.desc var",          "710.7582", NA, ls, NOAPP)
ck(A, "stat.desc std.dev",      "26.6600",  t1n(t9, "Mean (SD)", "Overall", 2), ls, TOL_T1, tol = 0.05)
ck(A, "stat.desc coef.var",     "0.7454",   NA, lk, NOAPP)
ck(A, "stat.desc skewness",     "0.2331",   NA, lk, NOAPP)
ck(A, "stat.desc skew.2SE",     "0.3776",   s9$skew2SE, lk)
ck(A, "stat.desc kurtosis",     "-1.3064",  NA, lk, NOAPP)
ck(A, "stat.desc kurt.2SE",     "-1.0735",  s9$kurt2SE, lk)
ck(A, "stat.desc normtest.W",   "0.9263",   s9$W, lk)
ck(A, "stat.desc normtest.p",   "0.0014",   s9$p, lp)
l <- L("W = 0.9263, p-value", 1800)
ck(A, "Shapiro W", "0.9263", s9$W, l); ck(A, "Shapiro p", "0.001388", s9$p, l)
l <- L("En este ejemplo, la variable X8", 1800)
ck(A, "Texto: Shapiro p (reporte)", "0.001", rx(x9$report, paste0("\\(W = ", NUM, ", p = ", NUM), 2), l)
l <- L("group  1  0.1285 0.7213", 1800)
ck(A, "Levene gl1", "1", l9$df1, l); ck(A, "Levene F", "0.1285", l9$F, l); ck(A, "Levene p", "0.7213", l9$p, l)
ck(A, "Levene gl2", "58", l9$df2, l + 1)
ck(A, "Shapiro log(X8) W (cap. 9.4.1)", "0.88346", NA, L("W = 0.88346", 1800), "la app no transforma la respuesta")
ck(A, "Shapiro log(X8) p (cap. 9.4.1)", "8.005e-05", NA, L("W = 0.88346", 1800), "la app no transforma la respuesta")

# =============================================================================
# Caps. 14–15 — t independiente: VickersF, protocolo Dual, VH ~ Group
# =============================================================================
s <- ado_get("t_indep"); dV <- s$example$prep(book_data("VickersF.txt")); dV$Group <- factor(dV$Group)
vV <- list(response = "VH", group = "Group")
x <- run_app("t_indep", dV, vV); xS <- run_app("t_indep", dV, vV, list(var_equal = TRUE)); tV <- T1(x$spec, dV, vV)
st <- x$res$stats; sw <- sw_chk(x); lv <- lv_chk(x)

for (cap in c(14, 15)) {
  A <- paste0("t_indep (cap. ", cap, ")")
  base <- if (cap == 14) 2700 else 2860
  # Tabla dplyr
  for (g in c("Multilink", "Ultimate")) {
    l0 <- L(paste0("^", g, "$"), base, fixed = FALSE)
    if (cap == 14) { cols <- c("n", "Media", "DE", "EE") } else { cols <- c("n", "Media", "DE", "Mediana", "IQR") }
    for (j in seq_along(cols)) ck(A, paste(g, cols[j]), BL[l0 + j], tab_val(x, TD, cols[j], g), l0 + j, if (cols[j] == "EE") NOAPP else "")
  }
  # table1
  lm1 <- L("41.9 (2.52)", base)
  ck(A, "table1 Multilink Mean", "41.9", t1n(tV, "Mean (SD)", "Multilink", 1), lm1); ck(A, "table1 Multilink SD", "2.52", t1n(tV, "Mean (SD)", "Multilink", 2), lm1)
  ck(A, "table1 Ultimate Mean", "42.5", t1n(tV, "Mean (SD)", "Ultimate", 1), lm1 + 1); ck(A, "table1 Ultimate SD", "2.21", t1n(tV, "Mean (SD)", "Ultimate", 2), lm1 + 1)
  ck(A, "table1 Overall Mean", "42.2", t1n(tV, "Mean (SD)", "Overall", 1), lm1 + 2); ck(A, "table1 Overall SD", "2.35", t1n(tV, "Mean (SD)", "Overall", 2), lm1 + 2)
  lmd <- L("41.9 [37.3, 46.0]", base)
  bk <- list(Multilink = c("41.9", "37.3", "46.0"), Ultimate = c("41.3", "39.8", "46.7"), Overall = c("41.7", "37.3", "46.7"))
  for (j in 1:3) for (k in 1:3) ck(A, paste("table1", names(bk)[j], c("Median", "Min", "Max")[k]), bk[[j]][k], t1n(tV, "Median [Min, Max]", names(bk)[j], k), lmd + j - 1)
  # Supuestos
  l <- L("W = 0.96541, p-value = 0.4223", base)
  ck(A, "Shapiro W", "0.96541", sw$W, l); ck(A, "Shapiro p", "0.4223", sw$p, l)
  l <- L("group  1   0.077 0.7834", base)
  ck(A, "Levene gl1", "1", lv$df1, l); ck(A, "Levene F", "0.077", lv$F, l); ck(A, "Levene p", "0.7834", lv$p, l); ck(A, "Levene gl2", "28", lv$df2, l + 1)
}

# Cap. 14: modelo lineal lm(VH ~ Group) — la app no ajusta lm; b0 = media Multilink y b1 = Δ sí se muestran
A <- "t_indep (cap. 14)"
l <- L("(Intercept)    41.8879", 2740)
ck(A, "lm b0 Estimate", "41.8879", st$m1, l, "b0 = media de Multilink (stats m1)")
ck(A, "lm b0 Std. Error", "0.6124", NA, l, NOAPP); ck(A, "lm b0 t value", "68.402", NA, l, NOAPP)
ck(A, "lm b1 Estimate", "0.5795", st$delta, l + 1, "b1 = Δ Ultimate − Multilink")
ck(A, "lm b1 Std. Error", "0.8660", NA, l + 1, NOAPP)
ck(A, "lm b1 |t value|", "0.669", abs(xS$res$stats$t), l + 1, "t de Student (opción var.equal = TRUE); la app muestra t con signo Multilink − Ultimate")
ck(A, "lm b1 Pr(>|t|)", "0.509", xS$res$stats$p, l + 1, "p de Student (var.equal = TRUE)")
lr <- L("    Min      1Q  Median      3Q     Max", 2740) + 1
for (q in c("-4.6025", "-1.7335", "-0.1388", "1.5779", "4.2240")) ck(A, paste("lm residuo", q), q, NA, lr, NOAPP)
l <- L("Residual standard error: 2.372", 2740)
ck(A, "lm Residual SE", "2.372", NA, l, NOAPP); ck(A, "lm gl residuales", "28", xS$res$stats$df, l, "gl de la t de Student (var.equal = TRUE)")
ck(A, "lm R2", "0.01574", NA, l + 1, NOAPP); ck(A, "lm R2 ajustado", "-0.01941", NA, l + 1, NOAPP)
ck(A, "lm F", "0.4477", NA, l + 2, NOAPP); ck(A, "lm F gl1", "1", NA, l + 2, NOAPP); ck(A, "lm F gl2", "28", xS$res$stats$df, l + 2, "gl de la t de Student")
ck(A, "lm F p-value", "0.5089", xS$res$stats$p, l + 2, "p de Student (var.equal = TRUE)")
l <- L("b0 (media Multilink): 41.8879", 2740)
ck(A, "b0 (media Multilink)", "41.8879", st$m1, l); ck(A, "b1 (Ultimate - Multilink)", "0.5795", st$delta, l + 1); ck(A, "Media Ultimate", "42.4674", st$m2, l + 2)
l <- L("El coeficiente  = 41.888", 2740)
ck(A, "Texto b0", "41.888", st$m1, l); ck(A, "Texto b1", "0.579", st$delta, l)

# Cap. 15: t.test, tamaño de efecto y reporte
A <- "t_indep (cap. 15)"
l <- L("t = -0.66914, df = 27.514", 2900)
ck(A, "t", "-0.66914", tab_val(x, "Resultado de t.test", "t"), l); ck(A, "df", "27.514", tab_val(x, "Resultado de t.test", "gl"), l)
ck(A, "p", "0.509", tab_val(x, "Resultado de t.test", "p"), l)
l <- L(" -2.354895  1.195907", 2900)
ck(A, "IC 95% inf (Multilink − Ultimate)", "-2.354895", -st$hi, l, "la app muestra el IC en dirección Ultimate − Multilink (signo invertido)")
ck(A, "IC 95% sup (Multilink − Ultimate)", "1.195907", -st$lo, l, "la app muestra el IC en dirección Ultimate − Multilink (signo invertido)")
l <- L("               41.88791                42.46741", 2900)
ck(A, "mean in group Multilink", "41.88791", st$m1, l); ck(A, "mean in group Ultimate", "42.46741", st$m2, l)
TE <- "Tama.o del efecto"
l <- L("t = -0.6691 ", 2950); ck(A, "t (cat)", "-0.6691", st$t, l)
l <- L("df = 27.51 ", 2950); ck(A, "df (cat)", "27.51", st$df, l)
l <- L("r (tamaño de efecto) = 0.127", 2950); ck(A, "r", "0.127", tab_val(x, TE, "r"), l)
l <- L("Diferencia de medias (Ultimate − Multilink) = 0.579", 2950); ck(A, "Δ Ultimate − Multilink", "0.579", tab_val(x, TE, "Δ medias"), l)
l <- L("IC 95%: [ -1.196 , 2.355 ]", 2950)
ck(A, "IC 95% inf (Ultimate − Multilink)", "-1.196", tab_val(x, TE, "IC 95% inf"), l); ck(A, "IC 95% sup (Ultimate − Multilink)", "2.355", tab_val(x, TE, "IC 95% sup"), l)
l <- L("Interpretación: La prueba de Shapiro-Wilk muestra W = 0.97", 2900)
ck(A, "Interpretación W", "0.97", sw$W, l); ck(A, "Interpretación p Shapiro", "0.42", sw$p, l)
ck(A, "Interpretación F", "0.077", lv$F, l); ck(A, "Interpretación p Levene", "0.78", lv$p, l)
# Reporte
R <- x$report; l <- L("La normalidad se evaluó mediante la prueba de Shapiro-Wilk (W = 0.97", 2950)
ck(A, "Reporte W", "0.97", rx(R, paste0("\\(W = ", NUM)), l)
ck(A, "Reporte p Shapiro", "0.42", rx(R, paste0("\\(W = ", NUM, ", p = ", NUM), 2), l)
ck(A, "Reporte F", "0.077", rx(R, paste0("\\(F = ", NUM)), l)
ck(A, "Reporte p Levene", "0.78", rx(R, paste0("\\(F = ", NUM, ", p = ", NUM), 2), l)
ck(A, "Reporte Media Multilink", "41.89", rx(R, paste0("\\[Multilink\\] \\(Media = ", NUM)), l)
ck(A, "Reporte DE Multilink", "2.52", rx(R, paste0("\\[Multilink\\] \\(Media = ", NUM, ", DE = ", NUM), 2), l)
ck(A, "Reporte Media Ultimate", "42.47", rx(R, paste0("\\[Ultimate\\] \\(Media = ", NUM)), l)
ck(A, "Reporte DE Ultimate", "2.21", rx(R, paste0("\\[Ultimate\\] \\(Media = ", NUM, ", DE = ", NUM), 2), l)
ck(A, "Reporte t", "-0.67", rx(R, paste0(" t = ", NUM)), l)
ck(A, "Reporte df", "27.5", rx(R, paste0("df = ", NUM)), l)
ck(A, "Reporte p", "0.51", rx(R, paste0("df = ", NUM, ", p = ", NUM), 2), l)
ck(A, "Reporte Δ", "0.58", rx(R, paste0("Δ = ", NUM)), l)
ck(A, "Reporte IC inf", "-1.2", rx(R, paste0("IC 95 %: ", NUM, ", ", NUM)), l)
ck(A, "Reporte IC sup", "2.35", rx(R, paste0("IC 95 %: ", NUM, ", ", NUM), 2), l)
ck(A, "Reporte r", "0.13", rx(R, paste0("\\(r = ", NUM, "\\)")), l)

# =============================================================================
# Cap. 16 — t apareada: Carbuse (formato wide, Ra 0_uses vs 1_use)
# =============================================================================
A <- "t_paired (cap. 16)"
s <- ado_get("t_paired"); dP <- s$example$prep(book_data("Carbuse.txt")); vP <- s$example$vars
x <- run_app("t_paired", dP, vP); tP <- T1(x$spec, dP, vP); st <- x$res$stats; sw <- sw_chk(x)
TDc <- "Estad.stica descriptiva por condici.n"
for (g in c("0_uses", "1_use")) {
  l0 <- L(paste0("^", g, "$"), 3060, fixed = FALSE); gg <- paste0("Ra_", g)
  cols <- c("n", "Media", "DE", "Mediana", "IQR")
  for (j in seq_along(cols)) ck(A, paste(g, cols[j]), BL[l0 + j], tab_val(x, TDc, cols[j], gg), l0 + j)
}
lm1 <- L("6.01 (1.61)", 3060)
cn <- c("Ra_0_uses", "Ra_1_use", "Overall"); bm <- list(c("6.01", "1.61"), c("5.74", "1.21"), c("5.87", "1.41"))
for (j in 1:3) { ck(A, paste("table1", cn[j], "Mean"), bm[[j]][1], t1n(tP, "Mean (SD)", cn[j], 1), lm1 + j - 1); ck(A, paste("table1", cn[j], "SD"), bm[[j]][2], t1n(tP, "Mean (SD)", cn[j], 2), lm1 + j - 1) }
lmd <- L("6.15 [3.60, 8.40]", 3060); bk <- list(c("6.15", "3.60", "8.40"), c("5.75", "3.90", "8.20"), c("5.95", "3.60", "8.40"))
for (j in 1:3) for (k in 1:3) ck(A, paste("table1", cn[j], c("Median", "Min", "Max")[k]), bk[[j]][k], t1n(tP, "Median [Min, Max]", cn[j], k), lmd + j - 1)
l <- L("W = 0.93632, p-value = 0.2041", 3100); ck(A, "Shapiro W (diferencias)", "0.93632", sw$W, l); ck(A, "Shapiro p (diferencias)", "0.2041", sw$p, l)
l <- L("muestra W = 0.936, p = 0.204", 3100); ck(A, "Texto W", "0.936", sw$W, l); ck(A, "Texto p", "0.204", sw$p, l)
TP <- "Resultado de t.test\\(paired"
l <- L("t = 0.66609, df = 19, p-value = 0.5134", 3100)
ck(A, "t", "0.66609", tab_val(x, TP, "t"), l); ck(A, "df", "19", tab_val(x, TP, "gl"), l); ck(A, "p", "0.5134", tab_val(x, TP, "p"), l)
l <- L(" -0.5676912  1.0976912", 3100); ck(A, "IC 95% inf", "-0.5676912", tab_val(x, TP, "IC 95% inf"), l); ck(A, "IC 95% sup", "1.0976912", tab_val(x, TP, "IC 95% sup"), l)
l <- L("          0.265 ", 3100); ck(A, "mean difference", "0.265", tab_val(x, TP, "Δ medias"), l)
l <- L("t = 0.6661 ", 3190); ck(A, "t (cat)", "0.6661", st$t, l)
l <- L("df = 19 ", 3190); ck(A, "df (cat)", "19", st$df, l)
l <- L("Delta (diferencia de medias) = 0.265", 3190); ck(A, "Δ (cat)", "0.265", st$delta, l)
l <- L("IC 95%: [ -0.568 , 1.098 ]", 3190); ck(A, "IC inf (cat)", "-0.568", st$lo, l); ck(A, "IC sup (cat)", "1.098", st$hi, l)
l <- L("r (tamaño de efecto) = 0.151", 3190); ck(A, "r", "0.151", tab_val(x, TP, "r"), l)
R <- x$report; l <- L("La normalidad de las diferencias fue evaluada mediante", 3190)
ck(A, "Reporte W", "0.94", rx(R, paste0("\\(W = ", NUM)), l)
ck(A, "Reporte p Shapiro", "0.20", rx(R, paste0("\\(W = ", NUM, ", p = ", NUM), 2), l)
ck(A, "Reporte Media antes", "6.005", rx(R, paste0("\\[Ra_0_uses\\] \\(Media = ", NUM)), l)
ck(A, "Reporte DE antes", "1.607", rx(R, paste0("\\[Ra_0_uses\\] \\(Media = ", NUM, ", DE = ", NUM), 2), l)
ck(A, "Reporte Media después", "5.740", rx(R, paste0("\\[Ra_1_use\\] \\(Media = ", NUM)), l)
ck(A, "Reporte DE después", "1.215", rx(R, paste0("\\[Ra_1_use\\] \\(Media = ", NUM, ", DE = ", NUM), 2), l)
ck(A, "Reporte t", "0.67", rx(R, paste0("\\(t = ", NUM)), l)
ck(A, "Reporte df", "19", rx(R, paste0("df = ", NUM)), l)
ck(A, "Reporte p", "0.51", rx(R, paste0("df = ", NUM, ", p = ", NUM), 2), l)
ck(A, "Reporte Δ", "0.26", rx(R, paste0("Δ = ", NUM)), l)
ck(A, "Reporte IC inf", "−0.57", rx(R, paste0("IC 95 %: ", NUM, ", ", NUM)), l)
ck(A, "Reporte IC sup", "1.10", rx(R, paste0("IC 95 %: ", NUM, ", ", NUM), 2), l)
ck(A, "Reporte r", "0.15", rx(R, paste0("\\(r = ", NUM, "\\)")), l)

# =============================================================================
# Caps. 29–30 — Wilcoxon suma de rangos: AUCFil, AUC ~ Group
# =============================================================================
s <- ado_get("wilcox_rs"); dW <- book_data("AUCFil.txt"); dW$Group <- factor(dW$Group); vW <- s$example$vars
x <- run_app("wilcox_rs", dW, vW); tW <- T1(x$spec, dW, vW); st <- x$res$stats; sw <- sw_chk(x); lv <- lv_chk(x)
RG <- "Rangos por grupo"; RW <- "Resultado de wilcox.test"; TDw <- "Estad.stica descriptiva \\(mediana"
for (cap in c(29, 30)) {
  A <- paste0("wilcox_rs (cap. ", cap, ")"); base <- if (cap == 29) 6100 else 6300
  l <- L("W = 0.45348, p-value = 8.363e-10", base); ck(A, "Shapiro W", "0.45348", sw$W, l); ck(A, "Shapiro p", "8.363e-10", sw$p, l)
  l <- L("group  1  4.0222  0.054", base)
  ck(A, "Levene gl1", "1", lv$df1, l); ck(A, "Levene F", "4.0222", lv$F, l); ck(A, "Levene p", "0.054", lv$p, l); ck(A, "Levene gl2", "30", lv$df2, l + 1)
  l <- L(if (cap == 29) "(Shapiro-Wilk p < 0.001)" else "resulta significativa (p < 0.001)", base); ck(A, "Texto: Shapiro p < 0.001", "<0.001", sw$p, l)
}
A <- "wilcox_rs (cap. 29)"
ck(A, "Rangos por observación (head(auc_ord, 12))", "10.5", NA, L("27  KERR 667.6", 6150), "la app no muestra el rango de cada observación (solo sumas y rangos medios)")
l <- L("319.5 208.5", 6150); ck(A, "Suma de rangos 3M", "319.5", tab_val(x, RG, "Suma de rangos", "3M"), l); ck(A, "Suma de rangos KERR", "208.5", tab_val(x, RG, "Suma de rangos", "KERR"), l)
l <- L("183.5  72.5", 6150); ck(A, "W 3M", "183.5", tab_val(x, RG, "W del grupo", "3M"), l); ck(A, "W KERR", "72.5", tab_val(x, RG, "W del grupo", "KERR"), l)
l <- L("W (el menor se usa como estadístico de prueba): 72.5", 6150); ck(A, "W menor", "72.5", tab_val(x, RW, "W menor (estadístico de prueba)"), l)
l <- L("La suma de rangos es de 319.5", 6150); ck(A, "Texto suma 3M", "319.5", tab_val(x, RG, "Suma de rangos", "3M"), l); ck(A, "Texto suma KERR", "208.5", tab_val(x, RG, "Suma de rangos", "KERR"), l)

A <- "wilcox_rs (cap. 30)"
for (g in c("3M", "KERR")) {
  l0 <- L(paste0("^", g, "$"), 6280, fixed = FALSE)
  cols <- c("n", "Mediana", "IQR", "Min", "Max")
  for (j in seq_along(cols)) ck(A, paste(g, cols[j]), BL[l0 + j], tab_val(x, TDw, cols[j], g), l0 + j)
}
lm1 <- L("686 (34.8)", 6300); cn <- c("3M", "KERR", "Overall"); bm <- list(c("686", "34.8"), c("668", "2.62"), c("677", "25.8"))
for (j in 1:3) { ck(A, paste("table1", cn[j], "Mean"), bm[[j]][1], t1n(tW, "Mean (SD)", cn[j], 1), lm1 + j - 1); ck(A, paste("table1", cn[j], "SD"), bm[[j]][2], t1n(tW, "Mean (SD)", cn[j], 2), lm1 + j - 1) }
lmd <- L("674 [664, 805]", 6300); bk <- list(c("674", "664", "805"), c("668", "663", "674"), c("669", "663", "805"))
for (j in 1:3) for (k in 1:3) ck(A, paste("table1", cn[j], c("Median", "Min", "Max")[k]), bk[[j]][k], t1n(tW, "Median [Min, Max]", cn[j], k), lmd + j - 1)
l <- L("W = 183.5, p-value = 0.03811", 6300); ck(A, "W (wilcox.test)", "183.5", tab_val(x, RW, "W (salida de R)"), l)
ck(A, "p", "0.03811", tab_val(x, RW, "p"), l)
l <- L("Efecto r = 0.367", 6300); ck(A, "r (a partir de Z)", "0.367", tab_val(x, RW, "r (a partir de Z)"), l)
l <- L("Efecto r (por rangos medios) = 0.434", 6300); ck(A, "r (rangos medios)", "0.434", tab_val(x, RW, "r (rangos medios)"), l)
l <- L("Media 685.69 668.15", 6300)
ck(A, "Media 3M", "685.69", tab_val(x, TDw, "Media", "3M"), l); ck(A, "Media KERR", "668.15", tab_val(x, TDw, "Media", "KERR"), l)
ck(A, "DE 3M", "34.77", tab_val(x, TDw, "DE", "3M"), l + 1); ck(A, "DE KERR", "2.62", tab_val(x, TDw, "DE", "KERR"), l + 1)
l <- L("Diferencia de medianas de AUC = 5.6", 6300); ck(A, "|Δ medianas|", "5.6", tab_val(x, RW, "|Δ medianas|"), l)
l <- L("Diferencia de medias de AUC   = 17.54", 6300); ck(A, "|Δ medias|", "17.54", tab_val(x, RW, "|Δ medias|"), l)
R <- x$report; l <- L("La normalidad de los datos fue evaluada con la prueba de Shapiro-Wilk (W = 0.45", 6300)
ck(A, "Reporte W Shapiro", "0.45", rx(R, paste0("Shapiro-Wilk \\(W = ", NUM)), l)
ck(A, "Reporte p Shapiro", "p < 0.001", rxs(R, paste0("Shapiro-Wilk \\(W = ", NUM, ", (p [<=] [0-9.]+)"), 2), l)
ck(A, "Reporte F", "4.02", rx(R, paste0("\\(F = ", NUM)), l)
ck(A, "Reporte p Levene", "0.054", rx(R, paste0("\\(F = ", NUM, ", p = ", NUM), 2), l)
ck(A, "Reporte media 3M", "685.69", rx(R, paste0("\\[3M\\] \\(media = ", NUM)), l)
ck(A, "Reporte DE 3M", "34.77", rx(R, paste0("\\[3M\\] \\(media = ", NUM, " ± ", NUM), 2), l)
ck(A, "Reporte mediana 3M", "673.85", rx(R, paste0("\\[3M\\] \\(media = ", NUM, " ± ", NUM, "; mediana = ", NUM), 3), l)
ck(A, "Reporte media KERR", "668.15", rx(R, paste0("\\[KERR\\] \\(media = ", NUM)), l)
ck(A, "Reporte DE KERR", "2.62", rx(R, paste0("\\[KERR\\] \\(media = ", NUM, " ± ", NUM), 2), l)
ck(A, "Reporte mediana KERR", "668.25", rx(R, paste0("\\[KERR\\] \\(media = ", NUM, " ± ", NUM, "; mediana = ", NUM), 3), l)
ck(A, "Reporte Δ medias", "17.54", rx(R, paste0("Δ de medias = ", NUM)), l)
ck(A, "Reporte W", "183.5", rx(R, paste0("Wilcoxon: W = ", NUM)), l)
ck(A, "Reporte p", "0.038", rx(R, paste0("Wilcoxon: W = ", NUM, ", p = ", NUM), 2), l)
ck(A, "Reporte r", "0.37", rx(R, paste0("\\(r = ", NUM, "\\)")), l)

# =============================================================================
# Caps. 31–32 — Wilcoxon rangos con signo: DXKTw grupo 2 (control), X0 vs X2
# =============================================================================
s <- ado_get("wilcox_sr"); dS <- s$example$prep(book_data("DXKTw.txt")); vS <- s$example$vars
x <- run_app("wilcox_sr", dS, vS); tS <- T1(x$spec, dS, vS); st <- x$res$stats; sw <- sw_chk(x)
RS <- "Rangos con signo"; RWs <- "Resultado de wilcox.test\\(paired"; TDs <- "Estad.stica descriptiva por condici.n"
for (cap in c(31, 32)) {
  A <- paste0("wilcox_sr (cap. ", cap, ")"); base <- if (cap == 31) 6440 else 6610
  l <- L("W = 0.92713, p-value = 0.04122", base); ck(A, "Shapiro W (diferencias)", "0.92713", sw$W, l); ck(A, "Shapiro p (diferencias)", "0.04122", sw$p, l)
}
A <- "wilcox_sr (cap. 31)"
ck(A, "Rangos con signo por observación (head(sr, 15))", "11.5", NA, L("40 P42          6", 6440), "la app no muestra el rango de cada par (solo W+, W−, T)")
l <- L("W+ (suma de rangos positivos): 122.5", 6480); ck(A, "W+", "122.5", tab_val(x, RS, "W+ (suma de rangos positivos)"), l)
l <- L("W- (suma de rangos negativos): 202.5", 6480); ck(A, "W−", "202.5", tab_val(x, RS, "W− (suma de rangos negativos)"), l)
l <- L("T (estadístico de prueba, el menor): 122.5", 6480); ck(A, "T", "122.5", tab_val(x, RS, "T (el menor)"), l)
l <- L("aquí T = 122.5", 6480); ck(A, "Texto T", "122.5", tab_val(x, RS, "T (el menor)"), l)
l <- L("donde n = 25 es el número de pares", 6480)
ck(A, "n pares con diferencia ≠ 0", "25", tab_val(x, RS, "Pares con diferencia ≠ 0"), l)
ck(A, "Diferencias iguales a 0 omitidas", "5", st$N - tab_val(x, RS, "Pares con diferencia ≠ 0"), l, "derivado: N (pares) − pares con diferencia ≠ 0, ambos mostrados por la app")

A <- "wilcox_sr (cap. 32)"
l <- L("la mediana pasa de 5 a 10", 6545); ck(A, "Figura: mediana X0", "5", st$med1, l); ck(A, "Figura: mediana X2", "10", st$med2, l)
for (g in c("Inmediato (X0)", "2 horas (X2)")) {
  l0 <- L(paste0("^", gsub("([()])", "\\\\\\1", g), "$"), 6560, fixed = FALSE); gg <- if (grepl("X0", g)) "X0" else "X2"
  cols <- c("n", "Mediana", "IQR")
  for (j in seq_along(cols)) ck(A, paste(gg, cols[j]), BL[l0 + j], tab_val(x, TDs, cols[j], gg), l0 + j)
}
lm1 <- L("10.6 (12.6)", 6560); cn <- c("X0", "X2", "Overall"); bm <- list(c("10.6", "12.6"), c("13.3", "11.6"), c("12.0", "12.1"))
for (j in 1:3) { ck(A, paste("table1", cn[j], "Mean"), bm[[j]][1], t1n(tS, "Mean (SD)", cn[j], 1), lm1 + j - 1); ck(A, paste("table1", cn[j], "SD"), bm[[j]][2], t1n(tS, "Mean (SD)", cn[j], 2), lm1 + j - 1) }
lmd <- L("5.00 [0, 40.0]", 6560); bk <- list(c("5.00", "0", "40.0"), c("10.0", "0", "40.0"), c("9.50", "0", "40.0"))
for (j in 1:3) for (k in 1:3) ck(A, paste("table1", cn[j], c("Median", "Min", "Max")[k]), bk[[j]][k], t1n(tS, "Median [Min, Max]", cn[j], k), lmd + j - 1)
l <- L("V = 122.5, p-value = 0.2871", 6560); ck(A, "V", "122.5", tab_val(x, RWs, "V"), l); ck(A, "p", "0.2871", tab_val(x, RWs, "p"), l)
l <- L("r (tamaño de efecto) = 0.194", 6560); ck(A, "r (a partir de Z)", "0.194", tab_val(x, RWs, "r (a partir de Z)"), l)
l <- L("Mediana           5.00        10.00", 6560)
ck(A, "Mediana X0", "5.00", tab_val(x, TDs, "Mediana", "X0"), l); ck(A, "Mediana X2", "10.00", tab_val(x, TDs, "Mediana", "X2"), l)
ck(A, "Media X0", "10.63", tab_val(x, TDs, "Media", "X0"), l + 1); ck(A, "Media X2", "13.33", tab_val(x, TDs, "Media", "X2"), l + 1)
ck(A, "DE X0", "12.57", tab_val(x, TDs, "DE", "X0"), l + 2); ck(A, "DE X2", "11.62", tab_val(x, TDs, "DE", "X2"), l + 2)
l <- L("Diferencia de medianas = 5", 6560); ck(A, "|Δ medianas|", "5", tab_val(x, RWs, "|Δ medianas|"), l)
l <- L("Diferencia de medias   = 2.7", 6560); ck(A, "|Δ medias|", "2.7", tab_val(x, RWs, "|Δ medias|"), l)
R <- x$report; l <- L("La normalidad de la diferencia entre X0 y X2", 6560)
ck(A, "Reporte W", "0.93", rx(R, paste0("\\(W = ", NUM)), l)
ck(A, "Reporte p Shapiro", "0.041", rx(R, paste0("\\(W = ", NUM, ", p = ", NUM), 2), l)
ck(A, "Reporte mediana X2", "10", rx(R, paste0("\\[X2\\] \\(mediana = ", NUM)), l)
ck(A, "Reporte media X2", "13.33", rx(R, paste0("\\[X2\\] \\(mediana = ", NUM, "; media = ", NUM), 2), l)
ck(A, "Reporte DE X2", "11.62", rx(R, paste0("\\[X2\\] \\(mediana = ", NUM, "; media = ", NUM, " ± ", NUM), 3), l)
ck(A, "Reporte mediana X0", "5", rx(R, paste0("\\[X0\\] \\(mediana = ", NUM)), l)
ck(A, "Reporte media X0", "10.63", rx(R, paste0("\\[X0\\] \\(mediana = ", NUM, "; media = ", NUM), 2), l)
ck(A, "Reporte DE X0", "12.57", rx(R, paste0("\\[X0\\] \\(mediana = ", NUM, "; media = ", NUM, " ± ", NUM), 3), l)
ck(A, "Reporte dif. medianas", "5", rx(R, paste0("diferencia de medianas de ", NUM)), l)
ck(A, "Reporte dif. medias", "2.7", rx(R, paste0("diferencia de medias de ", NUM)), l)
ck(A, "Reporte V", "122.5", rx(R, paste0("\\(V = ", NUM)), l)
ck(A, "Reporte p", "0.29", rx(R, paste0("\\(V = ", NUM, ", p = ", NUM), 2), l)
ck(A, "Reporte r", "0.19", rx(R, paste0("\\(r = ", NUM, "\\)")), l)

# =============================================================================
res <- write_results(file.path(Sys.getenv("VERIFY_DIR"), "results_t.csv"))
na_app <- res$app == "NA / no encontrado"
cat("\nTotal:", nrow(res), "| coinciden:", sum(res$match), "| no coinciden (app con valor):", sum(!res$match & !na_app),
    "| app NA:", sum(na_app), "\n")
cat("\n--- NO COINCIDEN ---\n"); print(res[!res$match & !na_app, c("analysis", "item", "book", "app", "line", "note")], right = FALSE)
cat("\n--- APP NA ---\n"); print(res[na_app, c("analysis", "item", "book", "line", "note")], right = FALSE)
if (any(is.na(res$line))) { cat("\n--- SIN LÍNEA ---\n"); print(res[is.na(res$line), 1:3]) }
