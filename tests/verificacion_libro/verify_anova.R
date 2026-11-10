# =============================================================================
# verify_anova.R — Verifica que la app reproduzca las cifras impresas en el libro
# para ANOVA de una vía (caps. 18–19), ANOVA factorial (23–24, 26), ANCOVA (25, 27)
# y Kruskal-Wallis (33–34). Solo usa salidas de la app (tablas, stats, checks, reporte).
# =============================================================================
source(file.path(Sys.getenv("VERIFY_DIR"), "verify_lib.R"))
options(width = 200)

NO   <- "la app no lo muestra"
GLHT <- "glht (multcomp) usa cuantiles por simulación: p/IC varían en el 4.º decimal entre corridas; tol = 0.001"

# ---- utilidades ---------------------------------------------------------------
# Valores impresos en notación científica (p. ej. "4.53e-05"): se compara con media unidad
# de la última cifra significativa impresa (el harness usa tolerancia relativa 1e-3, más
# estricta que la precisión con la que R imprime 2–3 cifras significativas).
.chk_h <- chk
chk <- function(analysis, item, book, app, line = NA_integer_, note = "", tol = NULL) {
  b <- trimws(as.character(book))
  if (is.null(tol) && grepl("^-?[0-9.]+e[-+]?[0-9]+$", b, ignore.case = TRUE)) {
    m <- sub("e.*$", "", b, ignore.case = TRUE); dg <- nchar(gsub("[^0-9]", "", m))
    ex <- as.integer(sub("^.*e", "", b, ignore.case = TRUE)); tol <- 0.5 * 10^(ex - (dg - 1)) * (1 + 1e-9)
  }
  .chk_h(analysis, item, book, app, line, note, tol)
}
# NA explícito: cifra impresa por el libro que la app no muestra
nas <- function(an, item, book, line, note = NO) chk(an, item, book, NA, line, note)
# Captura el grupo g de la n-ésima coincidencia de 'pat' en el reporte
rx <- function(rep, pat, g = 1, n = 1) {
  m <- regmatches(rep, gregexpr(pat, rep, perl = TRUE))[[1]]
  if (length(m) < n) return(NA_character_)
  regmatches(m[n], regexec(pat, m[n], perl = TRUE))[[1]][g + 1]
}
rnum <- function(rep, pat, g = 1, n = 1) as.numeric(rx(rep, pat, g, n))
# p del reporte: app escribe "= 0.003" o "< 0.001"; libro "0.003", "< 0.001" o "> 0.05"
chk_rp <- function(an, item, book, ptxt, line, note = "") {
  if (is.na(ptxt)) return(chk(an, item, book, NA, line, note))
  if (grepl("^[<>]", book)) {
    if (grepl("^>", book) && grepl("^=", ptxt)) {   # libro "p > lim": la app da el valor exacto
      lim <- as.numeric(sub("^>\\s*", "", book)); v <- as.numeric(sub("^=\\s*", "", ptxt))
      return(chk(an, item, paste("p", book), if (v > lim) paste("p", book) else paste("p", ptxt), line, note))
    }
    return(chk(an, item, paste("p", book), paste("p", ptxt), line, note))
  }
  if (grepl("^=", ptxt)) chk(an, item, book, as.numeric(sub("^=\\s*", "", ptxt)), line, note)
  else chk(an, item, paste("p =", book), paste("p", ptxt), line, note)
}
# Libro "p > 0.05" contra un p numérico de la app
chk_gt <- function(an, item, book, p, line, note = "") {
  lim <- as.numeric(sub("^.*>\\s*", "", book))
  chk(an, item, book, if (isTRUE(p > lim)) book else paste0("p = ", signif(p, 4)), line, note)
}
# Celda de table1 (descriptives de la app) como texto
t1 <- function(x, data, vars) as.data.frame(x$spec$descriptives(data, vars), stringsAsFactors = FALSE)
t1cell <- function(t, rowpat, col, nth = 1) { r <- grep(rowpat, t[[1]], fixed = TRUE); trimws(t[r[nth], col]) }
# Fila de una tabla de la app contra varios valores del libro
crow <- function(an, x, tabpat, row, cols, book, line, lab, note = "", tol = NULL) {
  for (i in seq_along(cols)) {
    if (is.na(book[i])) next
    tl <- if (is.list(tol)) tol[[i]] else tol
    nt <- if (length(note) > 1) note[i] else note
    chk(an, paste0(lab, " — ", cols[i]), book[i], tab_val(x, tabpat, cols[i], row), line, nt, tol = tl)
  }
}

# =============================================================================
# 1. ANOVA de una vía — caps. 18–19 (VickersF.txt, Multilink)
# =============================================================================
vf <- book_data("VickersF.txt"); ML <- subset(vf, Group == "Multilink")
v1 <- list(response = "VH", group = "Cure.time")
a_or <- run_app("anova1", ML, v1, list(posthoc = "tukey", contrastes = "ortogonal"))
a_pe <- run_app("anova1", ML, v1, list(posthoc = "tukey", contrastes = "personalizado", coef = "-2,1,1; 0,-1,1"))
a_tr <- run_app("anova1", ML, v1, list(contrastes = "tendencia"))
a_du <- run_app("anova1", ML, v1, list(posthoc = "dunnett", control = "Control"))
a_bo <- run_app("anova1", ML, v1, list(posthoc = "bonferroni"))
a_ho <- run_app("anova1", ML, v1, list(posthoc = "holm"))
a_bh <- run_app("anova1", ML, v1, list(posthoc = "bh"))
a_gh <- run_app("anova1", ML, v1, list(posthoc = "games_howell"))
AN <- "anova1"
grp <- c("Control", "Dual", "Light(20s)")
sw1 <- a_or$assumptions$checks[[1]]; lv1 <- a_or$assumptions$checks[[2]]

# -- 18.3 cálculo manual
nas(AN, "Media general", "41.90756", 3401)
nas(AN, "Varianza general", "5.18186", 3403)
nas(AN, "SST", "228.0019", 3405); nas(AN, "gl total", "44", 3405)
for (i in 1:3) {
  ln <- c(3431, 3435, 3439)[i]
  chk(AN, paste("n", grp[i]), "15", tab_val(a_or, "descriptiva", "n", grp[i]), ln)
  chk(AN, paste("media (5 dec.)", grp[i]), c("43.27447", "41.88791", "40.56029")[i], tab_val(a_or, "descriptiva", "Media", grp[i]), ln + 1,
      "la tabla descriptiva de la app redondea a 3 decimales")
  nas(AN, paste("varianza", grp[i]), c("2.14758", "6.37273", "3.81845")[i], ln + 2)
}
nas(AN, "Media general (texto)", "41.91", 3415)
for (i in 1:3) chk(AN, paste("media texto", grp[i]), c("43.27", "41.89", "40.56")[i], tab_val(a_or, "descriptiva", "Media", grp[i]), 3415)
chk(AN, "SSM (texto)", "55.26", tab_val(a_or, "Tabla ANOVA", "Sum Sq", 1), 3416)
chk(AN, "SSM", "55.25928", tab_val(a_or, "Tabla ANOVA", "Sum Sq", 1), 3445)
chk(AN, "df_M", "2", tab_val(a_or, "Tabla ANOVA", "Df", 1), 3445)
chk(AN, "SSR", "172.7426", tab_val(a_or, "Tabla ANOVA", "Sum Sq", 2), 3454)
chk(AN, "df_R", "42", tab_val(a_or, "Tabla ANOVA", "Df", 2), 3454)
chk(AN, "SSR (verificación)", "172.7426", tab_val(a_or, "Tabla ANOVA", "Sum Sq", 2), 3458)
chk(AN, "F (texto, 'unas 6.7 veces')", "6.7", a_or$res$stats$F, 3466)
# Tabla ANOVA manual
crow(AN, a_or, "Tabla ANOVA", 1, c("Sum Sq", "Df", "Mean Sq", "F value", "Pr(>F)"), c("55.2593", "2", "27.6296", "6.7178", "0.0029"), 3490, "Tabla ANOVA 18.3 Modelo")
crow(AN, a_or, "Tabla ANOVA", 2, c("Sum Sq", "Df", "Mean Sq"), c("172.7426", "42", "4.1129"), 3496, "Tabla ANOVA 18.3 Residuos")
nas(AN, "Tabla ANOVA 18.3 Total — SS", "228.0019", 3502); nas(AN, "Tabla ANOVA 18.3 Total — df", "44", 3503)

# -- summary.lm con contrastes ortogonales (cap. 18.5.3 con opción 'ortogonal'; cap. 19.5 con 'personalizado')
slm_block <- function(x, tabpat, l0, lab, rows, book_rows, intercept, footer, R2_is_eta = TRUE) {
  nas(AN, paste(lab, "— Intercepto (Estimate, SE, t, p)"), intercept, l0)
  for (i in seq_along(rows))
    crow(AN, x, tabpat, i, c("Estimación", "EE", "t", "p"), book_rows[[i]], l0 + i, paste(lab, "—", rows[i]))
  nas(AN, paste(lab, "— Residual standard error"), footer[1], l0 + 6)
  chk(AN, paste(lab, "— Multiple R-squared (= η²)"), footer[2], x$res$stats$eta2, l0 + 7, "la app lo muestra como η² en 'Tamaño de efecto'")
  nas(AN, paste(lab, "— Adjusted R-squared"), footer[3], l0 + 7)
  chk(AN, paste(lab, "— F-statistic"), footer[4], x$res$stats$F, l0 + 8)
  chk(AN, paste(lab, "— F gl1"), "2", x$res$stats$df1, l0 + 8)
  chk(AN, paste(lab, "— F gl2"), "42", x$res$stats$df2, l0 + 8)
  chk(AN, paste(lab, "— p global"), footer[5], x$res$stats$p, l0 + 8)
}
ft <- c("2.028", "0.2424", "0.2063", "6.718", "0.002942")
slm_block(a_or, "summary.lm", 3572, "summary.lm contrastes (cap. 18, opción ortogonal)", c("Control_vs_Exp", "Dual_vs_Light"),
          list(c("-0.6835", "0.2138", "-3.197", "0.00264"), c("-0.6638", "0.3703", "-1.793", "0.08021")), "41.9076 / 0.3023 / 138.619 / <2e-16", ft)
chk(AN, "p contraste 1 (texto, ≈)", "0.003", tab_val(a_or, "summary.lm", "p", 1), 3581)
chk(AN, "p contraste 2 (texto, ≈)", "0.080", tab_val(a_or, "summary.lm", "p", 2), 3581)

# -- TukeyHSD (cap. 18.6) → la app usa glht + confint (mismas estimaciones; p/IC por simulación)
tk_rows <- c("Dual - Control", "Light(20s) - Control", "Light(20s) - Dual")
tk_b <- list(c("-1.386552", "-3.185672", "0.4125684", "0.1593010"), c("-2.714175", "-4.513295", "-0.9150546", "0.0019510"),
             c("-1.327623", "-3.126743", "0.4714974", "0.1843476"))
for (i in 1:3) crow(AN, a_or, "Tukey HSD", tk_rows[i], c("Estimación", "IC 95% inf", "IC 95% sup", "p ajustada"), tk_b[[i]], 3635 + i,
                    paste("TukeyHSD", tk_rows[i]), note = c("", rep(paste("libro: TukeyHSD (ptukey exacto); app:", GLHT), 3)),
                    tol = list(NULL, 0.001, 0.001, 0.001))

# -- Tabla 18.2 (p por método)
chk(AN, "Tabla 18.2 caption: p Contraste 1", "0.003", tab_val(a_or, "summary.lm", "p", 1), 3670)
t182 <- list(Bonferroni = list(a_bo, "Bonferroni", c("0.204", "0.002", "0.241"), 3679),
             Tukey = list(a_or, "Tukey HSD", c("0.159", "0.002", "0.184"), 3680),
             Holm = list(a_ho, "Holm", c("0.136", "0.002", "0.136"), 3681),
             GamesHowell = list(a_gh, "Games-Howell", c("0.180", "0.001", "0.259"), 3682))
for (m in names(t182)) { z <- t182[[m]]
  for (i in 1:3) chk(AN, paste("Tabla 18.2", m, tk_rows[i]), z[[3]][i], tab_val(z[[1]], z[[2]], "p ajustada", tk_rows[i]), z[[4]] + 6 * (i - 1),
                     if (m == "Tukey") GLHT else "", tol = if (m == "Tukey") 0.001 else NULL) }
chk(AN, "Tabla 18.2 Contraste planeado Dual vs Light", "0.080", tab_val(a_or, "summary.lm", "p", 2), 3690)
chk(AN, "Texto 18.6.1: p contraste Dual vs Light", "0.080", tab_val(a_or, "summary.lm", "p", 2), 3695)
chk(AN, "Texto 18.6.1: Holm Dual vs Light", "0.136", tab_val(a_ho, "Holm", "p ajustada", "Light(20s) - Dual"), 3695)
chk(AN, "Texto 18.6.1: Bonferroni Dual vs Light", "0.241", tab_val(a_bo, "Bonferroni", "p ajustada", "Light(20s) - Dual"), 3695)

# -- 19.3 descriptiva
for (i in 1:3) crow(AN, a_or, "descriptiva", grp[i], c("n", "Media", "DE", "EE"),
                    list(c("15", "43.274", "1.465", "0.378"), c("15", "41.888", "2.524", "0.652"), c("15", "40.560", "1.954", "0.505"))[[i]],
                    3773 + 5 * (i - 1), paste("Descriptiva 19.3", grp[i]))
T1 <- t1(a_or, ML, v1)
for (i in 1:4) {
  chk(AN, paste("table1 Mean (SD)", c(grp, "Overall")[i]), c("43.3 (1.47)", "41.9 (2.52)", "40.6 (1.95)", "41.9 (2.28)")[i], t1cell(T1, "Mean (SD)", i + 1), 3799 + i)
  chk(AN, paste("table1 Median [Min, Max]", c(grp, "Overall")[i]), c("43.2 [40.4, 45.2]", "41.9 [37.3, 46.0]", "40.4 [36.5, 43.2]", "41.9 [36.5, 46.0]")[i],
      t1cell(T1, "Median [Min, Max]", i + 1), 3804 + i)
}
# -- 19.4 supuestos
chk(AN, "Fig 19.2: Shapiro W", "0.975", sw1$W, 3816); chk(AN, "Fig 19.2: Shapiro p", "0.44", sw1$p, 3816); chk(AN, "Fig 19.2: Levene p", "0.21", lv1$p, 3816)
chk(AN, "Shapiro W", "0.97519", sw1$W, 3823); chk(AN, "Shapiro p", "0.4398", sw1$p, 3823)
chk(AN, "Levene gl1", "2", lv1$df1, 3828); chk(AN, "Levene F", "1.6241", lv1$F, 3828); chk(AN, "Levene p", "0.2092", lv1$p, 3828); chk(AN, "Levene gl2", "42", lv1$df2, 3829)
chk(AN, "Texto 19.4: W", "0.975", sw1$W, 3830); chk(AN, "Texto 19.4: p SW", "0.44", sw1$p, 3830)
chk(AN, "Texto 19.4: Levene F", "1.624", lv1$F, 3830); chk(AN, "Texto 19.4: Levene p", "0.21", lv1$p, 3830)
# -- 19.5 contrastes (opción personalizada con los coeficientes del libro)
slm_block(a_pe, "summary.lm", 3851, "summary.lm contrastes (cap. 19.5, opción personalizado -2,1,1; 0,-1,1)", c("Control_vs_Exp", "Dual_vs_Light"),
          list(c("-0.6835", "0.2138", "-3.197", "0.00264"), c("-0.6638", "0.3703", "-1.793", "0.08021")), "41.9076 / 0.3023 / 138.619 / <2e-16", ft)
# -- 19.6 tendencia
slm_block(a_tr, "tendencia.*summary.lm", 3876, "summary.lm tendencia (contr.poly)", c("Lineal (.L)", "Cuadrática (.Q)"),
          list(c("-1.91921", "0.52364", "-3.665", "0.000688"), c("0.02406", "0.52364", "0.046", "0.963573")), "41.90756 / 0.30232 / 138.619 / <2e-16", ft)
# -- 19.7 post-hoc
for (i in 1:3) chk(AN, paste("Bonferroni", tk_rows[i]), c("0.2044", "0.0021", "0.2406")[i], tab_val(a_bo, "Bonferroni", "p ajustada", tk_rows[i]), c(3900, 3901, 3901)[i])
for (i in 1:3) chk(AN, paste("BH", tk_rows[i]), c("0.0802", "0.0021", "0.0802")[i], tab_val(a_bh, "Benjamini", "p ajustada", tk_rows[i]), c(3912, 3913, 3913)[i])
gl_b <- list(c("-1.3866", "0.7405", "-1.872", "0.15939"), c("-2.7142", "0.7405", "-3.665", "0.00194"), c("-1.3276", "0.7405", "-1.793", "0.18436"))
for (i in 1:3) crow(AN, a_or, "Tukey HSD", tk_rows[i], c("Estimación", "EE", "t", "p ajustada"), gl_b[[i]], 3931 + i, paste("glht Tukey", tk_rows[i]),
                    note = c("", "", "", GLHT), tol = list(NULL, NULL, NULL, 0.001))
nas(AN, "confint Tukey: Quantile", "2.4293", 3947)
ci_b <- list(c("-1.3866", "-3.1855", "0.4124"), c("-2.7142", "-4.5132", "-0.9152"), c("-1.3276", "-3.1266", "0.4714"))
for (i in 1:3) crow(AN, a_or, "Tukey HSD", tk_rows[i], c("Estimación", "IC 95% inf", "IC 95% sup"), ci_b[[i]], 3952 + i, paste("confint Tukey", tk_rows[i]),
                    note = c("", GLHT, GLHT), tol = list(NULL, 0.001, 0.001))
chk(AN, "Fig 19.3: p Control vs Light", "0.002", tab_val(a_or, "Tukey HSD", "p ajustada", "Light(20s) - Control"), 3959, GLHT, tol = 0.001)
du_b <- list(c("-1.3866", "0.7405", "-1.872", "0.12116"), c("-2.7142", "0.7405", "-3.665", "0.00134"))
for (i in 1:2) crow(AN, a_du, "Dunnett", tk_rows[i], c("Estimación", "EE", "t", "p ajustada"), du_b[[i]], 3972 + i, paste("Dunnett", tk_rows[i]),
                    note = c("", "", "", GLHT), tol = list(NULL, NULL, NULL, 0.001))
# -- 19.8 tamaño del efecto
chk(AN, "η²", "0.2424", tab_val(a_or, "Tamaño de efecto", "η²"), 3994)
chk(AN, "ω²", "0.2026", tab_val(a_or, "Tamaño de efecto", "ω²"), 3996)
chk(AN, "r (√η²)", "0.4923", tab_val(a_or, "Tamaño de efecto", "r (√η²)"), 3998)
chk(AN, "r Contraste 1", "0.442", tab_val(a_or, "summary.lm", "r", 1), 4011)
chk(AN, "r Contraste 2", "0.267", tab_val(a_or, "summary.lm", "r", 2), 4013)
# -- 19.9 Cómo reportar
R <- a_or$report; L <- 4017; RA <- "Reporte ANOVA 1 vía"
chk(AN, paste(RA, "W"), "0.975", rnum(R, "Shapiro-Wilk \\(W = ([0-9.]+)"), L, "la app escribe W con 2 decimales (fmt_num(W, 2))")
chk_rp(AN, paste(RA, "p SW"), "0.44", rx(R, "Shapiro-Wilk \\(W = [0-9.]+, p ([<=] [0-9.]+)"), L)
chk(AN, paste(RA, "Levene F"), "1.624", rnum(R, "Levene \\(F = ([0-9.]+)"), L)
chk_rp(AN, paste(RA, "Levene p"), "0.21", rx(R, "Levene \\(F = [0-9.]+, p ([<=] [0-9.]+)"), L)
chk(AN, paste(RA, "F gl1"), "2", rnum(R, "F\\((\\d+), (\\d+)\\)", 1), L); chk(AN, paste(RA, "F gl2"), "42", rnum(R, "F\\((\\d+), (\\d+)\\)", 2), L)
chk(AN, paste(RA, "F"), "6.72", rnum(R, "F\\(\\d+, \\d+\\) = ([0-9.]+)"), L)
chk_rp(AN, paste(RA, "p F"), "0.003", rx(R, "F\\(\\d+, \\d+\\) = [0-9.]+, p ([<=] [0-9.]+)"), L)
chk(AN, paste(RA, "ω²"), "0.20", rnum(R, "ω² = ([0-9.]+)"), L)
cpat <- "t = (-?[0-9.]+), df = (\\d+), p ([<=] [0-9.]+), r = ([0-9.]+)"
for (n in 1:2) {
  cl <- c("Control vs Exp", "Dual vs Light")[n]
  chk(AN, paste(RA, "t", cl), c("-3.20", "-1.79")[n], rnum(R, cpat, 1, n), L)
  chk(AN, paste(RA, "df", cl), "42", rnum(R, cpat, 2, n), L)
  chk_rp(AN, paste(RA, "p", cl), c("0.003", "0.08")[n], rx(R, cpat, 3, n), L)
  chk(AN, paste(RA, "r", cl), c("0.44", "0.27")[n], rnum(R, cpat, 4, n), L)
}
chk(AN, paste(RA, "Tukey Δ Control-Light"), "2.71", rnum(R, "Δ = ([0-9.]+)"), L)
chk_rp(AN, paste(RA, "Tukey p Control-Light"), "0.002", rx(R, "Δ = [0-9.]+, p ([<=] [0-9.]+)"), L)

# =============================================================================
# 2. ANOVA factorial — caps. 23, 24 y 26 (VickersF.txt completo)
# =============================================================================
FA <- "factorial"
vF <- list(response = "VH", factorA = "Group", factorB = "Cure.time")
f_am <- run_app("factorial", vf, vF, list(posthoc = "ambos"))
f_co <- run_app("factorial", vf, vF, list(posthoc = "contrastes"))
f_tk <- run_app("factorial", vf, vF, list(posthoc = "tukey"))
swF <- f_am$assumptions$checks[[1]]; lvF <- f_am$assumptions$checks[[2]]
A3 <- "ANOVA tipo III"; OM <- "Tamaño de efecto ω²"; SL <- "^Contrastes planeados"
src <- c("Group", "Cure.time", "Group × Cure.time")

# -- 23.4 tabla manual
mt <- list(c("109.555", "1", "109.555", "18.519"), c("306.499", "2", "153.250", "25.904"), c("106.466", "2", "53.233", "8.998"))
for (i in 1:3) {
  l0 <- 4804 + 5 * (i - 1)
  chk(FA, paste("Tabla 23.4", src[i], "— SS"), mt[[i]][1], tab_val(f_am, A3, "Sum Sq", i), l0, "app: SS tipo III (= secuencial en diseño balanceado)")
  chk(FA, paste("Tabla 23.4", src[i], "— df"), mt[[i]][2], tab_val(f_am, A3, "Df", i), l0 + 1)
  nas(FA, paste("Tabla 23.4", src[i], "— MS"), mt[[i]][3], l0 + 2)
  chk(FA, paste("Tabla 23.4", src[i], "— F"), mt[[i]][4], tab_val(f_am, A3, "F value", i), l0 + 3)
}
chk(FA, "Tabla 23.4 Residuos — SS", "496.939", tab_val(f_am, A3, "Sum Sq", 4), 4819)
chk(FA, "Tabla 23.4 Residuos — df", "84", tab_val(f_am, A3, "Df", 4), 4820)
nas(FA, "Tabla 23.4 Residuos — MS", "5.916", 4821)
nas(FA, "Tabla 23.4 Total — SS", "1019.459", 4824); nas(FA, "Tabla 23.4 Total — df", "89", 4825)
cellnm <- function(g, c) paste(g, "·", c)
cm <- list(Control = c("43.2745", "40.8044"), Dual = c("41.8879", "42.4674"), `Light(20s)` = c("40.5603", "35.8311"))
for (j in seq_along(cm)) for (g in 1:2)
  chk(FA, paste("Media de celda (4 dec.)", c("Multilink", "Ultimate")[g], names(cm)[j]), cm[[j]][g],
      tab_val(f_am, "descriptiva por celda", "Media", cellnm(c("Multilink", "Ultimate")[g], names(cm)[j])), 4840 + 3 * (j - 1) + g,
      "la tabla descriptiva de la app redondea a 3 decimales")
# -- 24.1 Levene
chk(FA, "Levene gl1", "5", lvF$df1, 4903); chk(FA, "Levene F", "1.8166", lvF$F, 4903); chk(FA, "Levene p", "0.1183", lvF$p, 4903); chk(FA, "Levene gl2", "84", lvF$df2, 4904)
# -- 24.2 descriptiva
dsc <- list(c("15", "43.274", "1.465"), c("15", "41.888", "2.524"), c("15", "40.560", "1.954"), c("15", "40.804", "2.317"), c("15", "42.467", "2.208"), c("15", "35.831", "3.593"))
cl6 <- c(cellnm("Multilink", grp), cellnm("Ultimate", grp))
for (i in 1:6) crow(FA, f_am, "descriptiva por celda", cl6[i], c("n", "Media", "DE"), dsc[[i]], 4920 + 5 * (i - 1), paste("Descriptiva 24.2", cl6[i]))
TF <- t1(f_am, vf, vF)
msd <- c("43.3 (1.47)", "41.9 (2.52)", "40.6 (1.95)", "40.8 (2.32)", "42.5 (2.21)", "35.8 (3.59)", "42.0 (2.28)", "42.2 (2.35)", "38.2 (3.72)")
mmm <- c("43.2 [40.4, 45.2]", "41.9 [37.3, 46.0]", "40.4 [36.5, 43.2]", "40.7 [37.3, 44.5]", "41.3 [39.8, 46.7]", "36.2 [30.5, 43.5]",
         "42.5 [37.3, 45.2]", "41.7 [37.3, 46.7]", "38.9 [30.5, 43.5]")
c9 <- c(cellnm("Multilink", grp), cellnm("Ultimate", grp), paste("Overall ·", grp))
for (i in 1:9) {
  chk(FA, paste("table1 Mean (SD)", c9[i]), msd[i], t1cell(TF, "Mean (SD)", i + 1), 4973 + i)
  chk(FA, paste("table1 Median [Min, Max]", c9[i]), mmm[i], t1cell(TF, "Median [Min, Max]", i + 1), 4983 + i)
}
# -- 24.3 Anova tipo III
nas(FA, "Anova III (Intercept) — SS/F/p", "149849 / 25329.6591 / <2.2e-16", 5014)
a3b <- list(c("110", "1", "18.5185", "4.525e-05"), c("306", "2", "25.9045", "1.725e-09"), c("106", "2", "8.9982", "0.0002878"))
for (i in 1:3) crow(FA, f_am, A3, i, c("Sum Sq", "Df", "F value", "Pr(>F)"), a3b[[i]], 5014 + i, paste("Anova III", src[i]))
crow(FA, f_am, A3, 4, c("Sum Sq", "Df"), c("497", "84"), 5018, "Anova III Residuals")
# summary.lm del modelo factorial
nas(FA, "summary.lm (Intercept) — Estimate/SE/t/p", "40.80426 / 0.25638 / 159.153 / <2e-16", 5032)
slf <- list(c("-1.10330", "0.25638", "-4.303", "4.53e-05"), c("-0.61759", "0.18129", "-3.407", "0.00101"), c("-1.99100", "0.31400", "-6.341", "1.10e-08"),
            c("0.06586", "0.18129", "0.363", "0.71729"), c("-1.32718", "0.31400", "-4.227", "6.00e-05"))
sln <- c("Ulti_vs_Multi", "Control_vs_Exp", "Dual_vs_Light", "Ulti_vs_Multi:Control_vs_Exp", "Ulti_vs_Multi:Dual_vs_Light")
for (i in 1:5) crow(FA, f_am, SL, i, c("Estimación", "EE", "t", "p"), slf[[i]], 5032 + i, paste("summary.lm", sln[i]))
nas(FA, "summary.lm Residual standard error", "2.432", 5048)
nas(FA, "summary.lm Multiple R-squared", "0.5125", 5049); nas(FA, "summary.lm Adjusted R-squared", "0.4835", 5049)
nas(FA, "summary.lm F-statistic (modelo global)", "17.66", 5050); nas(FA, "summary.lm p global", "6.345e-12", 5050)
# Texto 24.3 (≈)
chk(FA, "Texto: F Group ≈", "18.5", tab_val(f_am, A3, "F value", 1), 5052); chk(FA, "Texto: p Group", "<0.001", tab_val(f_am, A3, "Pr(>F)", 1), 5052)
chk(FA, "Texto: F Cure.time ≈", "25.9", tab_val(f_am, A3, "F value", 2), 5052); chk(FA, "Texto: p Cure.time", "<0.001", tab_val(f_am, A3, "Pr(>F)", 2), 5052)
chk(FA, "Texto: F interacción ≈", "9.0", tab_val(f_am, A3, "F value", 3), 5052); chk(FA, "Texto: p interacción", "<0.001", tab_val(f_am, A3, "Pr(>F)", 3), 5052)
chk(FA, "Texto: t Ulti_vs_Multi ≈", "-4.3", tab_val(f_am, SL, "t", 1), 5054); chk(FA, "Texto: p Ulti_vs_Multi", "<0.001", tab_val(f_am, SL, "p", 1), 5054)
chk(FA, "Texto: t Control_vs_Exp ≈", "-3.4", tab_val(f_am, SL, "t", 2), 5055); chk(FA, "Texto: p Control_vs_Exp", "<0.05", tab_val(f_am, SL, "p", 2), 5055)
chk(FA, "Texto: t Dual_vs_Light ≈", "-6.3", tab_val(f_am, SL, "t", 3), 5056); chk(FA, "Texto: p Dual_vs_Light", "<0.001", tab_val(f_am, SL, "p", 3), 5056)
chk(FA, "Texto: t Ulti×Dual_vs_Light ≈", "-4.2", tab_val(f_am, SL, "t", 5), 5057); chk(FA, "Texto: p Ulti×Dual_vs_Light", "<0.001", tab_val(f_am, SL, "p", 5), 5057)
# 24.3.1 interpretación detallada
nas(FA, "Media marginal Multilink", "41.908", 5061, "la app no muestra medias marginales por factor A (table1 solo da marginales del factor B)")
nas(FA, "Media marginal Ultimate", "39.701", 5061, "la app no muestra medias marginales por factor A")
chk(FA, "24.3.1 t Ulti_vs_Multi", "-4.303", tab_val(f_am, SL, "t", 1), 5061); chk(FA, "24.3.1 p Ulti_vs_Multi", "<0.001", tab_val(f_am, SL, "p", 1), 5061)
nas(FA, "Media marginal Control", "42.039", 5062, "table1 la muestra con 1 decimal (42.0)")
chk(FA, "24.3.1 t Control_vs_Exp", "-3.407", tab_val(f_am, SL, "t", 2), 5062); chk(FA, "24.3.1 p Control_vs_Exp", "0.001", tab_val(f_am, SL, "p", 2), 5062)
nas(FA, "Media marginal Dual", "42.178", 5063, "table1 la muestra con 1 decimal (42.2)"); nas(FA, "Media marginal Light(20s)", "38.196", 5063, "table1 la muestra con 1 decimal (38.2)")
chk(FA, "24.3.1 t Dual_vs_Light", "-6.341", tab_val(f_am, SL, "t", 3), 5063); chk(FA, "24.3.1 p Dual_vs_Light", "<0.001", tab_val(f_am, SL, "p", 3), 5063)
nas(FA, "Contraste Control vs Exp en Multilink", "-2.050", 5067); nas(FA, "Contraste Control vs Exp en Ultimate", "-1.655", 5067)
chk(FA, "24.3.1 t Ulti×Control_vs_Exp", "0.363", tab_val(f_am, SL, "t", 4), 5067); chk(FA, "24.3.1 p Ulti×Control_vs_Exp", "0.717", tab_val(f_am, SL, "p", 4), 5067)
nas(FA, "Contraste Dual vs Light en Multilink", "-1.33", 5071); nas(FA, "Contraste Dual vs Light en Ultimate", "-6.64", 5071)
chk(FA, "24.3.1 t Ulti×Dual_vs_Light", "-4.227", tab_val(f_am, SL, "t", 5), 5071); chk(FA, "24.3.1 p Ulti×Dual_vs_Light", "<0.001", tab_val(f_am, SL, "p", 5), 5071)
chk(FA, "Fig 24.2 t Ulti×Control_vs_Exp", "0.363", tab_val(f_am, SL, "t", 4), 5119); chk(FA, "Fig 24.2 p Ulti×Control_vs_Exp", "0.717", tab_val(f_am, SL, "p", 4), 5119)
chk(FA, "Fig 24.2 t Ulti×Dual_vs_Light", "-4.227", tab_val(f_am, SL, "t", 5), 5119); chk(FA, "Fig 24.2 p Ulti×Dual_vs_Light", "<0.001", tab_val(f_am, SL, "p", 5), 5119)
# -- 24.5 efectos simples (opción 'contrastes'): se emparejan por nombre
SE <- "Efectos simples: contrastes"
nas(FA, "Efectos simples (Intercept) — Estimate/SE/t/p", "40.8043 / 0.2564 / 159.153 / <2e-16", 5207)
seb <- list(list("Control_vs_Dual+Light(20s)", "Con_vs_Exp", c("-0.6176", "0.1813", "-3.407", "0.00101"), 5208),
            list("Dual_vs_Light(20s)", "Dual_vs_Light", c("1.9910", "0.3140", "6.341", "1.10e-08"), 5209),
            list("Multilink_vs_Ultimate en Dual", "Dual_Multi_vs_Dual_Ulti", c("0.2897", "0.4441", "0.652", "0.51587"), 5210),
            list("Multilink_vs_Ultimate en Control", "Con_Multi_vs_Con_Ulti", c("-1.2350", "0.4441", "-2.781", "0.00669"), 5211),
            list("Multilink_vs_Ultimate en Light(20s)", "Light_Multi_vs_Light_Ulti", c("-2.3646", "0.4441", "-5.325", "8.29e-07"), 5212))
for (z in seb) crow(FA, f_co, SE, z[[1]], c("Estimación", "EE", "t", "p"), z[[3]], z[[4]], paste("Efectos simples", z[[2]]),
                    note = if (z[[2]] == "Dual_vs_Light") "el libro codifica Dual_vs_Light como c(-1,1,0,...) (Light −1, Dual +1); la app usa Dual −1, Light +1 (como en el modelo factorial del libro)" else "")
nas(FA, "Efectos simples Residual SE / R² / R² aj.", "2.432 / 0.5125 / 0.4835", 5216)
nas(FA, "Efectos simples F global / p", "17.66 / 6.345e-12", 5218)
# Tukey sobre el factor combinado (opción 'tukey')
lab6 <- c(Mult_L20 = "Multilink · Light(20s)", Mult_Dual = "Multilink · Dual", Mult_Con = "Multilink · Control",
          Ulti_L20 = "Ultimate · Light(20s)", Ulti_Dual = "Ultimate · Dual", Ulti_Con = "Ultimate · Control")
tkb <- read.table(text = "
Mult_Dual Mult_L20   1.3276230  -1.2626767  3.9179228 0.6684770
Mult_Con Mult_L20    2.7141750   0.1238753  5.3044748 0.0344588
Ulti_L20 Mult_L20   -4.7292426  -7.3195424 -2.1389429 0.0000121
Ulti_Dual Mult_L20   1.9071169  -0.6831829  4.4974166 0.2737021
Ulti_Con Mult_L20    0.2441194  -2.3461803  2.8344192 0.9997804
Mult_Con Mult_Dual   1.3865520  -1.2037477  3.9768518 0.6260396
Ulti_L20 Mult_Dual  -6.0568656  -8.6471654 -3.4665659 0.0000000
Ulti_Dual Mult_Dual  0.5794938  -2.0108059  3.1697936 0.9864693
Ulti_Con Mult_Dual  -1.0835036  -3.6738034  1.5067961 0.8259324
Ulti_L20 Mult_Con   -7.4434177 -10.0337174 -4.8531179 0.0000000
Ulti_Dual Mult_Con  -0.8070582  -3.3973579  1.7832416 0.9432626
Ulti_Con Mult_Con   -2.4700556  -5.0603554  0.1202441 0.0705239
Ulti_Dual Ulti_L20   6.6363595   4.0460597  9.2266592 0.0000000
Ulti_Con Ulti_L20    4.9733620   2.3830623  7.5636618 0.0000039
Ulti_Con Ulti_Dual  -1.6629975  -4.2532972  0.9273023 0.4261581", colClasses = "character")
TK <- f_tk$posthoc$tables[["Efectos simples: Tukey HSD sobre el factor combinado A×B"]]
for (i in seq_len(nrow(tkb))) {
  b <- lab6[[tkb$V1[i]]]; a <- lab6[[tkb$V2[i]]]; nm <- paste0(tkb$V1[i], "-", tkb$V2[i]); ln <- 5228 + i
  k <- match(paste0(b, "-", a), TK$Comparación); flip <- FALSE
  if (is.na(k)) { k <- match(paste0(a, "-", b), TK$Comparación); flip <- TRUE }
  nt <- if (flip) "la app presenta la comparación en orden inverso (orden alfabético de celdas); se invierte el signo de diff/IC para comparar" else ""
  d <- TK$Diferencia[k]; lw <- TK$`IC 95% inf`[k]; up <- TK$`IC 95% sup`[k]
  if (flip) { d <- -d; tmp <- -up; up <- -lw; lw <- tmp }
  chk(FA, paste("Tukey celdas", nm, "— diff"), tkb$V3[i], d, ln, nt)
  chk(FA, paste("Tukey celdas", nm, "— lwr"), tkb$V4[i], lw, ln, nt)
  chk(FA, paste("Tukey celdas", nm, "— upr"), tkb$V5[i], up, ln, nt)
  chk(FA, paste("Tukey celdas", nm, "— p adj"), tkb$V6[i], TK$`p ajustada`[k], ln)
}
# -- 24.6 ω²
omb <- list(c("109.5545", "1", "0.1011"), c("306.4992", "2", "0.2874"), c("106.4659", "2", "0.0923"))
for (i in 1:3) crow(FA, f_am, OM, i, c("SS", "gl", "ω²"), omb[[i]], 5266 + 4 * (i - 1), paste("ω² parcial", src[i]))
# -- 24.7 Cómo reportar
R <- f_am$report; L <- 5280; RF <- "Reporte factorial"
chk_rp(FA, paste(RF, "Shapiro p"), "> 0.05", rx(R, "W = [0-9.]+, p ([<=] [0-9.]+)"), L)
chk_rp(FA, paste(RF, "Levene p"), "> 0.05", rx(R, "Levene \\(F = [0-9.]+, p ([<=] [0-9.]+)"), L)
fpat <- "F\\[(\\d+),(\\d+)\\] = ([0-9.]+), p ([<=] [0-9.]+), ω² = ([0-9.]+)"
fb <- list(c("1", "84", "18.5", "< 0.001", "0.10"), c("2", "84", "25.9", "< 0.001", "0.29"), c("2", "84", "9.0", "< 0.001", "0.09"))
for (n in 1:3) {
  chk(FA, paste(RF, src[n], "gl1"), fb[[n]][1], rnum(R, fpat, 1, n), L); chk(FA, paste(RF, src[n], "gl2"), fb[[n]][2], rnum(R, fpat, 2, n), L)
  chk(FA, paste(RF, src[n], "F"), fb[[n]][3], rnum(R, fpat, 3, n), L); chk_rp(FA, paste(RF, src[n], "p"), fb[[n]][4], rx(R, fpat, 4, n), L)
  chk(FA, paste(RF, src[n], "ω²"), fb[[n]][5], rnum(R, fpat, 5, n), L)
}
# -- Cap. 26: tipos de SS
nas(FA, "dat_ex anova(Y ~ A + B): SS/MS/F/p A", "1.600 / 1.6000 / 0.1807 / 0.6835", 5387, "modelo aditivo sin interacción: la app (factorial) siempre ajusta A*B")
nas(FA, "dat_ex anova(Y ~ A + B): SS/MS/F/p B", "0.017 / 0.0167 / 0.0019 / 0.9666", 5388, "modelo aditivo sin interacción: la app (factorial) siempre ajusta A*B")
nas(FA, "dat_ex Residuals SS/MS", "61.983 / 8.8548", 5389, "modelo aditivo sin interacción: la app (factorial) siempre ajusta A*B")
# Balanceado: Tipo I (SS secuenciales = tabla ω² de la app), Tipo II (no), Tipo III
for (i in 1:3) chk(FA, paste("Tipo I balanceado SS", src[i]), c("109.55", "306.50", "106.47")[i], tab_val(f_am, OM, "SS", i), 5433 + i, "app: SS secuenciales en la tabla ω²")
nas(FA, "Tipo I balanceado MS/F/p (tabla secuencial)", "109.555 / 18.5185 / 4.525e-05 ...", 5434, "la app no muestra la tabla tipo I (solo sus SS en la tabla ω²)")
chk(FA, "Tipo I balanceado Residuals SS", "496.94", tab_val(f_am, A3, "Sum Sq", 4), 5437)
nas(FA, "Tipo II balanceado (tabla completa)", "109.55 / 306.50 / 106.47 ...", 5449, "la app no ofrece SS tipo II")
for (i in 1:3) crow(FA, f_am, A3, i, c("Sum Sq", "Df", "F value", "Pr(>F)"), a3b[[i]], 5466 + i, paste("Tipo III balanceado (contr.sum)", src[i]))
# Desbalanceado (se quitan 5 + 8 observaciones como en el libro)
q <- c(which(vf$Group == "Multilink" & vf$Cure.time == "Control")[1:5], which(vf$Group == "Ultimate" & vf$Cure.time == "Light(20s)")[1:8])
db <- vf[-q, ]
u_ab <- run_app("factorial", db, vF, list(posthoc = "contrastes"))
u_ba <- run_app("factorial", db, list(response = "VH", factorA = "Cure.time", factorB = "Group"), list(posthoc = "contrastes"))
for (i in 1:3) chk(FA, paste("Desbalanceado Tipo I (Group primero) SS", src[i]), c("42.60", "217.66", "115.11")[i], tab_val(u_ab, OM, "SS", i), 5500 + i,
                   "app: SS secuenciales en la tabla ω² (factorA = Group)")
nas(FA, "Desbalanceado Tipo I (Group primero) MS/F/p", "42.597 / 8.1623 / 0.005607 ...", 5501, "la app no muestra la tabla tipo I")
chk(FA, "Desbalanceado Residuals SS", "370.53", tab_val(u_ab, A3, "Sum Sq", 4), 5504); chk(FA, "Desbalanceado Residuals df", "71", tab_val(u_ab, A3, "Df", 4), 5504)
for (i in 1:3) chk(FA, paste("Desbalanceado Tipo I (Cure.time primero) SS", c("Cure.time", "Group", "Cure.time:Group")[i]), c("171.88", "88.38", "115.11")[i],
                   tab_val(u_ba, OM, "SS", i), 5514 + i, "app: SS secuenciales en la tabla ω² (factorA = Cure.time)")
nas(FA, "Desbalanceado Tipo I (Cure.time primero) MS/F/p", "85.941 / 16.468 / 1.332e-06 ...", 5515, "la app no muestra la tabla tipo I")
nas(FA, "Desbalanceado Tipo II (tabla completa)", "88.38 / 217.66 / 115.11 ...", 5531, "la app no ofrece SS tipo II")
nas(FA, "Desbalanceado Tipo III (Intercept)", "117388 / 22493.537", 5548)
a3u <- list(c("125", "1", "23.919", "6.047e-06"), c("279", "2", "26.737", "2.210e-09"), c("115", "2", "11.028", "6.749e-05"))
for (i in 1:3) crow(FA, u_ab, A3, i, c("Sum Sq", "Df", "F value", "Pr(>F)"), a3u[[i]], 5548 + i, paste("Desbalanceado Tipo III", src[i]),
                    note = "app: contrastes ortogonales del libro; libro: contr.sum (equivalentes para tipo III)")
crow(FA, u_ab, A3, 4, c("Sum Sq", "Df"), c("371", "71"), 5552, "Desbalanceado Tipo III Residuals")
chk(FA, "Texto 26.3.1: SS Group Tipo I (primero)", "42.60", tab_val(u_ab, OM, "SS", 1), 5555)
nas(FA, "Texto 26.3.1: p Group Tipo I (primero)", "0.006", 5555, "la app no muestra el p de la tabla tipo I")
chk(FA, "Texto 26.3.1: SS Group Tipo I (después)", "88.38", tab_val(u_ba, OM, "SS", 2), 5555)
nas(FA, "Texto 26.3.1: p Group Tipo I (después)", "<0.001", 5555, "la app no muestra el p de la tabla tipo I")
nas(FA, "Texto 26.3.1: SS Group Tipo II", "88.38", 5555, "la app no ofrece SS tipo II")
chk(FA, "Texto 26.3.1: SS Group Tipo III", "124.83", tab_val(u_ab, A3, "Sum Sq", 1), 5555)

# =============================================================================
# 3. ANCOVA — caps. 25 y 27 (datBF.txt, control = G2)
# =============================================================================
AC <- "ancova"
dB <- book_data("datBF.txt"); dB$Group <- factor(dB$Group, levels = c("G2", "G1", "G3"))
vA <- list(response = "Final", group = "Group", covariate = "Basal")
an <- run_app("ancova", dB, vA, list(control = "G2"))
ck <- an$assumptions$checks  # 1 SW Final, 2 SW Basal, 3 Lev Final, 4 Lev Basal, 5 independencia, 6 pendientes
g3 <- c("G2", "G1", "G3")
# -- 27.3 descriptiva
adb <- list(c("45", "17.104", "1.902", "14.944", "2.079"), c("39", "17.137", "2.052", "12.773", "2.341"), c("39", "17.125", "2.226", "12.329", "2.452"))
for (i in 1:3) crow(AC, an, "descriptiva de covariable", g3[i], c("n", "Basal Media", "Basal DE", "Final Media", "Final DE"), adb[[i]], 5664 + 6 * (i - 1),
                    paste("Descriptiva 27.3", g3[i]))
TA <- t1(an, dB, vA)
cA <- c(g3, "Overall")
b_ms <- c("17.1 (1.90)", "17.1 (2.05)", "17.1 (2.23)", "17.1 (2.04)"); b_md <- c("17.2 [12.3, 20.0]", "17.3 [13.7, 20.3]", "16.7 [13.2, 20.0]", "17.2 [12.3, 20.3]")
f_ms <- c("14.9 (2.08)", "12.8 (2.34)", "12.3 (2.45)", "13.4 (2.55)"); f_md <- c("15.3 [10.3, 18.3]", "12.3 [7.33, 18.3]", "13.2 [7.16, 17.3]", "13.7 [7.16, 18.3]")
for (i in 1:4) {
  chk(AC, paste("table1 Basal Mean (SD)", cA[i]), b_ms[i], t1cell(TA, "Mean (SD)", i + 1, 1), 5693 + i)
  chk(AC, paste("table1 Basal Median [Min, Max]", cA[i]), b_md[i], t1cell(TA, "Median [Min, Max]", i + 1, 1), 5698 + i)
  chk(AC, paste("table1 Final Mean (SD)", cA[i]), f_ms[i], t1cell(TA, "Mean (SD)", i + 1, 2), 5704 + i)
  chk(AC, paste("table1 Final Median [Min, Max]", cA[i]), f_md[i], t1cell(TA, "Median [Min, Max]", i + 1, 2), 5709 + i)
}
# -- 27.4 supuestos
chk(AC, "Shapiro Basal W", "0.94035", ck[[2]]$W, 5730); chk(AC, "Shapiro Basal p", "3.636e-05", ck[[2]]$p, 5730)
chk(AC, "Shapiro Final W", "0.98461", ck[[1]]$W, 5736); chk(AC, "Shapiro Final p", "0.1768", ck[[1]]$p, 5736)
chk(AC, "Levene Basal gl1", "2", ck[[4]]$df1, 5740); chk(AC, "Levene Basal F", "1.2202", ck[[4]]$F, 5740); chk(AC, "Levene Basal p", "0.2988", ck[[4]]$p, 5740)
chk(AC, "Levene Basal gl2", "120", ck[[4]]$df2, 5741)
chk(AC, "Levene Final gl1", "2", ck[[3]]$df1, 5745); chk(AC, "Levene Final F", "0.7124", ck[[3]]$F, 5745); chk(AC, "Levene Final p", "0.4925", ck[[3]]$p, 5745)
chk(AC, "Levene Final gl2", "120", ck[[3]]$df2, 5746)
# -- 27.5 independencia de la covariable
nas(AC, "ANOVA Basal ~ Group: Df / SS / MS", "2 / 0.0 / 0.012", 5754, "la app solo muestra F y p de la independencia")
chk(AC, "ANOVA Basal ~ Group: F", "0.003", ck[[5]]$F, 5754); chk(AC, "ANOVA Basal ~ Group: p", "0.997", ck[[5]]$p, 5754)
nas(AC, "ANOVA Basal ~ Group: Residuals Df / SS / MS", "120 / 507.5 / 4.229", 5755, "la app solo muestra F y p de la independencia")
chk(AC, "Texto 27.5: p independencia", "0.997", ck[[5]]$p, 5756); chk(AC, "Texto 27.5: W Basal", "0.94", ck[[2]]$W, 5756)
chk(AC, "Texto 27.5: p SW Basal", "<0.001", ck[[2]]$p, 5756); chk(AC, "Texto 27.5: Levene F Basal", "1.22", ck[[4]]$F, 5756)
chk(AC, "Texto 27.5: Levene p Basal", "0.299", ck[[4]]$p, 5756)
# -- 27.6 ANCOVA tipo III y summary.lm
A3A <- "ANCOVA tipo III"; CT <- "Contrastes planeados"
nas(AC, "Anova III (Intercept) SS/Df/F/p", "0.86 / 1 / 0.2927 / 0.5895", 5786)
crow(AC, an, A3A, 1, c("Sum Sq", "Df", "F value", "Pr(>F)"), c("169.78", "2", "28.8218", "6.205e-11"), 5787, "Anova III Group")
crow(AC, an, A3A, 2, c("Sum Sq", "Df", "F value", "Pr(>F)"), c("276.50", "1", "93.8780", "<2.2e-16"), 5788, "Anova III Basal")
crow(AC, an, A3A, 3, c("Sum Sq", "Df"), c("350.49", "119"), 5789, "Anova III Residuals")
nas(AC, "summary.lm (Intercept) Estimate/SE/t/p", "0.71066 / 1.31352 / 0.541 / 0.589", 5803)
crow(AC, an, CT, 1, c("Estimación", "EE", "t", "p"), c("-0.80417", "0.10709", "-7.509", "1.2e-11"), 5804, "summary.lm Control_vs_Exp")
crow(AC, an, CT, 2, c("Estimación", "EE", "t", "p"), c("-0.21766", "0.19432", "-1.120", "0.265"), 5805, "summary.lm G1_vs_G3")
nas(AC, "summary.lm Basal (pendiente) Estimate/SE/t/p", "0.73811 / 0.07618 / 9.689 / <2e-16", 5806)
nas(AC, "summary.lm Residual SE / R² / R² aj.", "1.716 / 0.5587 / 0.5476", 5810)
nas(AC, "summary.lm F global / p", "50.22 / <2.2e-16", 5812)
CMP <- "Comparación ANOVA vs. ANCOVA"
nas(AC, "Fig 27.3: SST", "794.20", 5813)
chk(AC, "Fig 27.3: SSM (grupo, ANOVA)", "167.21", tab_val(an, CMP, "SS grupo", 1), 5813)
chk(AC, "Fig 27.3: SSR ANOVA", "626.99", tab_val(an, CMP, "SS residual", 1), 5813)
chk(AC, "Fig 27.3: SS Basal", "276.50", tab_val(an, A3A, "Sum Sq", 2), 5813)
chk(AC, "Fig 27.3: SSR ANCOVA", "350.49", tab_val(an, CMP, "SS residual", 2), 5813)
chk(AC, "Fig 27.3: reducción ≈ %", "44", tab_val(an, CMP, "Reducción SSR (%)", 2), 5813)
# -- 27.7 ANOVA vs ANCOVA
crow(AC, an, CMP, 1, c("SS grupo", "F grupo", "p grupo"), c("167.2", "16", "6.92e-07"), 5823, "ANOVA sin covariable Group")
nas(AC, "ANOVA sin covariable Group Df / MS", "2 / 83.60", 5823)
crow(AC, an, CMP, 1, c("SS residual", "gl residual"), c("627.0", "120"), 5824, "ANOVA sin covariable Residuals")
nas(AC, "ANOVA sin covariable Residuals MS", "5.22", 5824)
crow(AC, an, CMP, 2, c("SS grupo", "F grupo", "p grupo"), c("167.2", "28.39", "8.33e-11"), 5832, "ANCOVA (secuencial) Group")
nas(AC, "ANCOVA (secuencial) Group Df / MS", "2 / 83.60", 5832)
chk(AC, "ANCOVA (secuencial) Basal SS", "276.5", tab_val(an, A3A, "Sum Sq", 2), 5833, "app: tabla tipo III (igual a la secuencial para el último término)")
chk(AC, "ANCOVA (secuencial) Basal F", "93.88", tab_val(an, A3A, "F value", 2), 5833, "app: tabla tipo III (igual a la secuencial para el último término)")
chk(AC, "ANCOVA (secuencial) Basal p", "<2e-16", tab_val(an, A3A, "Pr(>F)", 2), 5833)
nas(AC, "ANCOVA (secuencial) Basal MS", "276.50", 5833)
crow(AC, an, CMP, 2, c("SS residual", "gl residual"), c("350.5", "119"), 5834, "ANCOVA (secuencial) Residuals")
nas(AC, "ANCOVA (secuencial) Residuals MS", "2.95", 5834)
chk(AC, "SSR ANOVA", "626.99", tab_val(an, CMP, "SS residual", 1), 5844)
chk(AC, "SSR ANCOVA", "350.49", tab_val(an, CMP, "SS residual", 2), 5845)
chk(AC, "Reducción SSR %", "44.1", tab_val(an, CMP, "Reducción SSR (%)", 2), 5846)
chk(AC, "Texto: reducción ~%", "44", tab_val(an, CMP, "Reducción SSR (%)", 2), 5847)
# -- 27.8 medias ajustadas
MA <- "Medias ajustadas"
for (i in 1:3) {
  chk(AC, paste("Media ajustada", g3[i]), c("14.95615", "12.76130", "12.32597")[i], tab_val(an, MA, "Media ajustada", g3[i]), 5858)
  chk(AC, paste("IC inf media ajustada", g3[i]), c("14.44956", "12.21714", "11.78182")[i], tab_val(an, MA, "IC inf", g3[i]), 5863)
  chk(AC, paste("IC sup media ajustada", g3[i]), c("15.46273", "13.30546", "12.87012")[i], tab_val(an, MA, "IC sup", g3[i]), 5868)
}
# -- 27.9 Tukey sobre medias ajustadas
TKA <- "Tukey sobre medias ajustadas"; ar <- c("G1 - G2", "G3 - G2", "G3 - G1")
tka <- list(c("-2.1948", "0.3755", "-5.846", "<1e-04"), c("-2.6302", "0.3755", "-7.005", "<1e-04"), c("-0.4353", "0.3886", "-1.120", "0.504"))
for (i in 1:3) crow(AC, an, TKA, ar[i], c("Estimación", "EE", "t", "p ajustada"), tka[[i]], 5892 + i, paste("glht Tukey", ar[i]),
                    note = c("", "", "", GLHT), tol = list(NULL, NULL, NULL, if (i == 3) 0.001 else NULL))
nas(AC, "confint Tukey: Quantile", "2.3736", 5908)
tkc <- list(c("-2.1948", "-3.0861", "-1.3036"), c("-2.6302", "-3.5214", "-1.7390"), c("-0.4353", "-1.3578", "0.4871"))
for (i in 1:3) crow(AC, an, TKA, ar[i], c("Estimación", "IC 95% inf", "IC 95% sup"), tkc[[i]], 5913 + i, paste("confint Tukey", ar[i]),
                    note = c("", GLHT, GLHT), tol = list(NULL, 0.001, 0.001))
# -- 27.10 homogeneidad de pendientes
nas(AC, "Pendientes: (Intercept) / Group / Basal (SS, F, p)", "0.57 / 1.66 0.2808 0.7557 / 280.66 95.2120 <2e-16", 5927, "la app solo muestra F y p de la interacción Group:Basal")
nas(AC, "Pendientes: Group:Basal SS / Df", "5.61 / 2", 5930, "la app solo muestra F y p de la interacción")
chk(AC, "Pendientes: Group:Basal F", "0.9515", ck[[6]]$F, 5930); chk(AC, "Pendientes: Group:Basal p", "0.3891", ck[[6]]$p, 5930)
nas(AC, "Pendientes: Residuals SS / Df", "344.88 / 117", 5931)
# -- 27.11 tamaño del efecto
chk(AC, "η² parcial Group", "0.3263", tab_val(an, "η² parcial", "η² parcial", 1), 5946)
chk(AC, "η² parcial Basal", "0.441", tab_val(an, "η² parcial", "η² parcial", 2), 5948)
chk(AC, "r Control vs Exp", "0.567", tab_val(an, CT, "r", 1), 5958)
chk(AC, "r G1 vs G3", "0.102", tab_val(an, CT, "r", 2), 5960)
# -- 27.12 Cómo reportar
R <- an$report; L <- 5964; RC <- "Reporte ANCOVA"
chk(AC, paste(RC, "W Final"), "0.98", rnum(R, "Final\\]? \\(W = ([0-9.]+)"), L)
chk_rp(AC, paste(RC, "p SW Final"), "0.18", rx(R, "Final\\]? \\(W = [0-9.]+, p ([<=] [0-9.]+)"), L)
chk(AC, paste(RC, "W Basal"), "0.94", rnum(R, "Basal\\]? \\(W = ([0-9.]+)"), L)
chk_rp(AC, paste(RC, "p SW Basal"), "< 0.001", rx(R, "Basal\\]? \\(W = [0-9.]+, p ([<=] [0-9.]+)"), L)
chk_rp(AC, paste(RC, "Levene p Final"), "> 0.05", rx(R, "Final: F = [0-9.]+, p ([<=] [0-9.]+)"), L)
chk_rp(AC, paste(RC, "Levene p Basal"), "> 0.05", rx(R, "Basal: F = [0-9.]+, p ([<=] [0-9.]+)"), L)
chk(AC, paste(RC, "F independencia"), "0.003", rnum(R, "ANOVA: F = ([0-9.]+)"), L)
chk_rp(AC, paste(RC, "p independencia"), "0.997", rx(R, "ANOVA: F = [0-9.]+, p ([<=] [0-9.]+)"), L)
epat <- "F\\[(\\d+),(\\d+)\\] = ([0-9.]+), p ([<=] [0-9.]+), η²parcial = ([0-9.]+)"
eb <- list(c("2", "119", "28.82", "< 0.001", "0.33"), c("1", "119", "93.88", "< 0.001", "0.44")); en <- c("grupo", "covariable")
for (n in 1:2) {
  chk(AC, paste(RC, en[n], "gl1"), eb[[n]][1], rnum(R, epat, 1, n), L); chk(AC, paste(RC, en[n], "gl2"), eb[[n]][2], rnum(R, epat, 2, n), L)
  chk(AC, paste(RC, en[n], "F"), eb[[n]][3], rnum(R, epat, 3, n), L); chk_rp(AC, paste(RC, en[n], "p"), eb[[n]][4], rx(R, epat, 4, n), L)
  chk(AC, paste(RC, en[n], "η² parcial"), eb[[n]][5], rnum(R, epat, 5, n), L)
}
chk_rp(AC, paste(RC, "p pendientes"), "0.39", rx(R, "Basal: p ([<=] [0-9.]+)"), L)
chk_rp(AC, paste(RC, "Tukey p G2 vs G1"), "< 0.001", rx(R, "entre G2 y G1 \\(Δ = [0-9.]+, p ([<=] [0-9.]+)"), L)
chk_rp(AC, paste(RC, "Tukey p G2 vs G3"), "< 0.001", rx(R, "entre G2 y G3 \\(Δ = [0-9.]+, p ([<=] [0-9.]+)"), L)
chk(AC, paste(RC, "r Control vs Exp"), "0.57", rnum(R, "r = ([0-9.]+)", 1, 1), L)
chk(AC, paste(RC, "r G1 vs G3"), "0.10", rnum(R, "r = ([0-9.]+)", 1, 2), L)

# =============================================================================
# 4. Kruskal-Wallis — caps. 33–34 (ReMat.txt)
# =============================================================================
KW <- "kruskal"
rm <- book_data("ReMat.txt"); vK <- list(response = "Response", group = "Groups")
k_sc <- run_app("kruskal", rm, vK, list(posthoc = "sc"))
k_wb <- run_app("kruskal", rm, vK, list(posthoc = "wilcox_bonf"))
swK <- k_sc$assumptions$checks[[1]]; lvK <- k_sc$assumptions$checks[[2]]
gk <- c("CenN", "EQForteF", "RivaLC")
for (ln in c(6740, 6883)) { chk(KW, "Shapiro W", "0.93895", swK$W, ln); chk(KW, "Shapiro p", "0.0003702", swK$p, ln) }
for (ln in c(6744, 6887)) { chk(KW, "Levene gl1", "2", lvK$df1, ln); chk(KW, "Levene F", "4.1702", lvK$F, ln); chk(KW, "Levene p", "0.01864", lvK$p, ln)
  chk(KW, "Levene gl2", "87", lvK$df2, ln + 1) }
RG <- "Rangos por grupo"
for (i in 1:3) {
  chk(KW, paste("Suma de rangos", gk[i]), c("2166.0", "741.5", "1187.5")[i], tab_val(k_sc, RG, "Suma de rangos (R_k)", gk[i]), 6763)
  chk(KW, paste("n", gk[i]), "30", tab_val(k_sc, RG, "n", gk[i]), 6767)
  chk(KW, paste("Rango medio", gk[i]), c("72.20000", "24.71667", "39.58333")[i], tab_val(k_sc, RG, "Rango medio", gk[i]), 6779)
}
chk(KW, "N", "90", k_sc$res$stats$N, 6769)
chk(KW, "H (cálculo manual)", "51.8613", k_sc$res$stats$H, 6771, "el libro calcula H sin corrección por empates; la app muestra kruskal.test() (con corrección)")
chk(KW, "p (cálculo manual)", "5.4761e-12", k_sc$res$stats$p, 6773, "p del H manual sin corrección por empates")
# 34.2 descriptiva
KD <- "Estadística descriptiva"
kdb <- list(c("30", "166.17", "38.35", "171.71", "49.05"), c("30", "78.04", "23.08", "74.09", "26.36"), c("30", "96.36", "22.25", "100.22", "26.48"))
for (i in 1:3) crow(KW, k_sc, KD, gk[i], c("n", "Media", "DE", "Mediana", "IQR"), kdb[[i]], 6837 + 6 * (i - 1), paste("Descriptiva 34.2", gk[i]))
TK1 <- t1(k_sc, rm, vK); ck4 <- c(gk, "Overall")
for (i in 1:4) {
  chk(KW, paste("table1 Mean (SD)", ck4[i]), c("166 (38.4)", "78.0 (23.1)", "96.4 (22.2)", "114 (47.7)")[i], t1cell(TK1, "Mean (SD)", i + 1), 6866 + i)
  chk(KW, paste("table1 Median [Min, Max]", ck4[i]), c("172 [79.3, 250]", "74.1 [38.8, 131]", "100 [44.5, 137]", "102 [38.8, 250]")[i], t1cell(TK1, "Median [Min, Max]", i + 1), 6871 + i)
}
# 34.4 kruskal.test
chk(KW, "H (kruskal.test)", "51.863", tab_val(k_sc, "kruskal.test", "H"), 6899)
chk(KW, "gl", "2", tab_val(k_sc, "kruskal.test", "gl"), 6899)
chk(KW, "p (kruskal.test)", "5.47e-12", tab_val(k_sc, "kruskal.test", "p"), 6899)
# 34.5 Siegel y Castellan
SC <- "Siegel"
nas(KW, "z crítico [α/(k(k−1))]", "2.39", 6903)
chk(KW, "Texto: rango medio CenN", "72.20", tab_val(k_sc, RG, "Rango medio", "CenN"), 6904)
chk(KW, "Texto: rango medio RivaLC", "39.58", tab_val(k_sc, RG, "Rango medio", "RivaLC"), 6904)
chk(KW, "Texto: rango medio EQForteF", "24.72", tab_val(k_sc, RG, "Rango medio", "EQForteF"), 6904)
chk(KW, "Texto: dif. CenN − EQForteF", "47.48", tab_val(k_sc, SC, "Dif. observada", "CenN - EQForteF"), 6904)
chk(KW, "Texto: dif. CenN − RivaLC", "32.62", tab_val(k_sc, SC, "Dif. observada", "CenN - RivaLC"), 6904)
chk(KW, "Texto: dif. RivaLC − EQForteF", "14.86", tab_val(k_sc, SC, "Dif. observada", "EQForteF - RivaLC"), 6904)
chk(KW, "Texto: dif. crítica kruskalmc", "16.15", tab_val(k_sc, SC, "Dif. crítica", 1), 6905)
scr <- c("CenN - EQForteF", "CenN - RivaLC", "EQForteF - RivaLC")
for (i in 1:3) {
  chk(KW, paste("kruskalmc obs.dif", scr[i]), c("47.48333", "32.61667", "14.86667")[i], tab_val(k_sc, SC, "Dif. observada", scr[i]), 6913 + i)
  chk(KW, paste("kruskalmc critical.dif", scr[i]), "16.14828", tab_val(k_sc, SC, "Dif. crítica", scr[i]), 6913 + i)
  chk(KW, paste("kruskalmc stat.signif", scr[i]), c("Sí", "Sí", "No")[i], tab_val(k_sc, SC, "Significativa", scr[i]), 6913 + i, "libro TRUE/FALSE; app Sí/No")
}
# pairwise Wilcoxon (Bonferroni)
WB <- "pairwise.wilcox"
chk(KW, "Wilcoxon Bonf. EQForteF vs CenN", "1.8e-09", tab_val(k_wb, WB, "p ajustada (Bonferroni)", "EQForteF - CenN"), 6927)
chk(KW, "Wilcoxon Bonf. RivaLC vs CenN", "5.1e-08", tab_val(k_wb, WB, "p ajustada (Bonferroni)", "RivaLC - CenN"), 6928, "R imprime 2 cifras significativas; tolerancia relativa del harness 1e-3")
chk(KW, "Wilcoxon Bonf. RivaLC vs EQForteF", "0.0077", tab_val(k_wb, WB, "p ajustada (Bonferroni)", "RivaLC - EQForteF"), 6928)
# 34.7 η²_H
chk(KW, "η²_H (texto, ≈ 57 %)", "0.57", tab_val(k_sc, "kruskal.test", "η²_H"), 6945)
chk(KW, "η²_H", "0.5731", tab_val(k_sc, "kruskal.test", "η²_H"), 6951)
# 34.8 Cómo reportar
R <- k_sc$report; L <- 6955; RK <- "Reporte Kruskal-Wallis"
chk(KW, paste(RK, "W"), "0.94", rnum(R, "W = ([0-9.]+)"), L)
chk_rp(KW, paste(RK, "p SW"), "< 0.001", rx(R, "W = [0-9.]+, p ([<=] [0-9.]+)"), L)
chk(KW, paste(RK, "Levene F"), "4.17", rnum(R, "Levene \\(F = ([0-9.]+)"), L)
chk_rp(KW, paste(RK, "Levene p"), "0.019", rx(R, "Levene \\(F = [0-9.]+, p ([<=] [0-9.]+)"), L)
mpat <- "([0-9.]+) ± ([0-9.]+)"; mo <- c("CenN", "RivaLC", "EQForteF")
mb <- list(c("166.17", "38.35"), c("96.36", "22.25"), c("78.04", "23.08"))
for (n in 1:3) { chk(KW, paste(RK, "media", mo[n]), mb[[n]][1], rnum(R, mpat, 1, n), L); chk(KW, paste(RK, "DE", mo[n]), mb[[n]][2], rnum(R, mpat, 2, n), L) }
hpat <- "H\\[(\\d+)\\] = ([0-9.]+), p ([<=] [0-9.]+), η² = ([0-9.]+)"
chk(KW, paste(RK, "gl H"), "2", rnum(R, hpat, 1), L); chk(KW, paste(RK, "H"), "51.86", rnum(R, hpat, 2), L)
chk_rp(KW, paste(RK, "p H"), "< 0.001", rx(R, hpat, 3), L); chk(KW, paste(RK, "η²"), "0.57", rnum(R, hpat, 4), L)
dpat <- "Δ = ([0-9.]+), p ([<>] [0-9.]+)"; dn <- c("CenN vs EQForteF", "CenN vs RivaLC", "RivaLC vs EQForteF")
db_ <- list(c("88.12", "< 0.05"), c("69.81", "< 0.05"), c("18.31", "> 0.05"))
for (n in 1:3) {
  chk(KW, paste(RK, "Δ", dn[n]), db_[[n]][1], rnum(R, dpat, 1, n), L, "el libro añade la unidad (MPa); la app deja la unidad al usuario")
  chk(KW, paste(RK, "p", dn[n]), paste("p", db_[[n]][2]), paste("p", rx(R, dpat, 2, n)), L)
}

# =============================================================================
out <- file.path(Sys.getenv("VERIFY_DIR"), "results_anova.csv")
res <- write_results(out)
isna <- grepl("^NA", res$app)
cat("\nTotal:", nrow(res), "| coinciden:", sum(res$match), "| no coinciden:", sum(!res$match & !isna), "| la app no lo muestra:", sum(isna), "\n")
cat("\n=== NO COINCIDEN ===\n")
print(res[!res$match & !isna, c("analysis", "item", "book", "app", "line", "note")], right = FALSE, row.names = FALSE)
