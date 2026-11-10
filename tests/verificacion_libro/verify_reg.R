# =============================================================================
# verify_reg.R — Verifica que la app reproduce las cifras del libro en:
#   cap. 11 (regresión lineal simple), 48–49 (lineal múltiple), 50–51 (logística),
#   52–53 (Kaplan-Meier / log-rank) y 54–55 (Cox).
# Solo se usan salidas de la app (stats, tablas, checks, reporte); no se recalcula nada.
# Validación cruzada / partición / bootstrap: la app usa las mismas semillas que el libro
# (lineal: set.seed(123); logística: partición set.seed(123) y bootstrap set.seed(42);
#  Cox: bootstrap set.seed(1)).
# =============================================================================
source(file.path(Sys.getenv("VERIFY_DIR"), "verify_lib.R"))
suppressPackageStartupMessages(library(MASS))

# --- utilidades ----------------------------------------------------------------
# k-ésimo grupo de captura de la primera coincidencia de `pattern` en `txt` (numérico)
rx <- function(txt, pattern, k = 1) {
  m <- regmatches(txt, regexec(pattern, txt, perl = TRUE))[[1]]
  if (length(m) < k + 1) return(NA_real_)
  s <- m[k + 1]; if (grepl("^<", trimws(s))) return(0)   # "< 0.001" → valor por debajo del límite
  suppressWarnings(as.numeric(sub("\\.$", "", s)))       # quita el punto final de la oración
}
# p impresos en notación científica ("6e-04", "1.14e-04"): tolerancia = media unidad de la última
# cifra impresa (el harness usa 0.1 % relativo, demasiado estricto para 1 cifra significativa)
chk0 <- chk
chk <- function(analysis, item, book, app, line = NA_integer_, note = "", tol = NULL) {
  b <- trimws(as.character(book))
  if (is.null(tol) && grepl("^-?[0-9]+(\\.[0-9]+)?e[-+]?[0-9]+$", b)) {
    mant <- sub("e.*", "", b); ex <- as.integer(sub(".*e", "", b))
    dec <- if (grepl("\\.", mant)) nchar(sub(".*\\.", "", mant)) else 0
    tol <- 0.5 * 10^(ex - dec) * (1 + 1e-9)
  }
  chk0(analysis, item, book, app, line, note, tol)
}
# Valor p escrito en el reporte tras `prefix` ("p = 0.054" o "p < 0.001")
rp <- function(txt, prefix) {
  m <- regmatches(txt, regexec(paste0(prefix, "p (=|<) ([0-9.]+)"), txt, perl = TRUE))[[1]]
  if (length(m) < 3) return(NA_real_)
  if (m[2] == "<") 0 else as.numeric(m[3])
}
NOAPP <- "la app no lo muestra"
na <- function(A, item, book, line, note = NOAPP) chk(A, item, book, NA, line, note)
tv <- function(x, re, col, row = 1) tab_val(x, re, col, row)
# Coeficientes: tabla de la app con columnas Término/Estimación/EE/Estadístico/p
coefs <- function(A, x, re, rows, line0, pref = "") {
  for (i in seq_len(nrow(rows))) { r <- rows[i, ]; ln <- line0 + i - 1
    chk(A, paste0(pref, r$term, " b"),  r$b,  tv(x, re, "Estimación", r$term), ln)
    chk(A, paste0(pref, r$term, " EE"), r$se, tv(x, re, "EE", r$term), ln)
    chk(A, paste0(pref, r$term, " t/z"), r$t, tv(x, re, "Estadístico", r$term), ln)
    chk(A, paste0(pref, r$term, " p"),  r$p,  tv(x, re, "p", r$term), ln) }
}
cdf <- function(...) { m <- rbind(...); data.frame(term = m[, 1], b = m[, 2], se = m[, 3], t = m[, 4], p = m[, 5], stringsAsFactors = FALSE) }
# Captura las leyendas de un gráfico base de la app (la app las dibuja con legend())
capture_legends <- function(f) {
  .LEG <- list()
  assign("legend", function(...) { a <- list(...); .LEG[[length(.LEG) + 1]] <<- a$legend; invisible(NULL) }, envir = globalenv())
  on.exit(rm("legend", envir = globalenv()))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  f(); .LEG
}
normterms <- function(s) paste(sort(trimws(strsplit(gsub("_c", "", gsub("`", "", s)), "\\+")[[1]])), collapse = " + ")

# =============================================================================
# CAP. 11 — Regresión lineal simple (curva de calibración)
# =============================================================================
A <- "lineal (cap. 11)"
Conc <- data.frame(Absorbancia = c(0.1472, 0.3949, 0.5635, 0.5861, 1.6347), Concentracion = c(1, 3, 5, 7, 15))
x <- run_app("lineal", Conc, list(response = "Concentracion", predictors = "Absorbancia"), list(cv_rep = FALSE))
coefs(A, x, "Modelo inicial \\(todos", cdf(c("(Intercept)", "-0.02724", "0.74553", "-0.037", "0.97315"),
                                         c("Absorbancia", "9.36033", "0.88969", "10.521", "0.00183")), 2295)
for (i in 1:5) na(A, paste("Residuo", i), c("-0.3506", "-0.6692", "-0.2473", "1.5412", "-0.2741")[i], 2291)
na(A, "Error estándar residual", "1.014", 2300)
chk(A, "R²", "0.9736", tv(x, "Modelo inicial: ajuste", "R²"), 2301)
chk(A, "R² ajustado", "0.9648", tv(x, "Modelo inicial: ajuste", "R² ajustado"), 2301)
chk(A, "F", "110.7", tv(x, "Modelo inicial: ajuste", "F"), 2302)
chk(A, "gl1", "1", tv(x, "Modelo inicial: ajuste", "gl1"), 2302)
chk(A, "gl2", "3", tv(x, "Modelo inicial: ajuste", "gl2"), 2302)
chk(A, "p (F)", "0.001834", tv(x, "Modelo inicial: ajuste", "p"), 2302)
na(A, "SS total", "116.8", 2318); na(A, "SS modelo", "113.7179", 2320); na(A, "SS residual", "3.0821", 2322)
chk(A, "R² (manual)", "0.9736", x$res$stats$R2, 2324)
chk(A, "F calculado", "110.688", x$res$stats$F, 2339)
chk(A, "p-value F", "0.00183", x$res$stats$p, 2341)
chk(A, "Texto: R² = 0.974", "0.974", rx(x$report, "R² = ([0-9.]+)"), 2327, "del reporte de la app")
chk(A, "Texto: F = 110.7", "110.7", rx(x$report, "F\\(1, 3\\) = ([0-9.]+)"), 2342, "del reporte de la app (imprime 110.69)")
chk(A, "Texto: p = 0.002", "0.002", rp(x$report, "110\\.69, "), 2342, "del reporte de la app")
na(A, "Predicción Abs = 0.5635", "5.2473", 2386, "la app no lo muestra: solo predice el conjunto de prueba 80/20, que con n = 5 queda vacío")
na(A, "Predicción IC inf", "3.776", 2388, "la app no lo muestra (sin predicción para valores nuevos en lineal)")
na(A, "Predicción IC sup", "6.718", 2388, "la app no lo muestra (sin predicción para valores nuevos en lineal)")

# =============================================================================
# CAP. 48 — Regresión múltiple: coeficientes simple vs. múltiple
# =============================================================================
A <- "lineal (cap. 48)"
bw <- MASS::birthwt
x <- run_app("lineal", bw, list(response = "bwt", predictors = "lwt"), list(cv_rep = FALSE))
chk(A, "b(lwt) simple", "4.429", tv(x, "Modelo inicial \\(todos", "Estimación", "lwt"), 9191)
x <- run_app("lineal", bw, list(response = "bwt", predictors = c("lwt", "age"), factors = "smoke"), list(cv_rep = FALSE))
chk(A, "b(lwt) múltiple", "4.016", tv(x, "Modelo inicial \\(todos", "Estimación", "lwt"), 9195)
chk(A, "b(age) múltiple", "7.154", tv(x, "Modelo inicial \\(todos", "Estimación", "age"), 9197)
chk(A, "b(smoke1) múltiple", "-269.257", tv(x, "Modelo inicial \\(todos", "Estimación", "smoke1"), 9199)

# =============================================================================
# CAP. 48–49 — Regresión lineal múltiple: flujo completo (birthwt)
# Opciones del libro: factores ht, race, smoke, ui; polinomios age=2, lwt=4; CV repetida.
# =============================================================================
A <- "lineal"
vL <- list(response = "bwt", predictors = c("age", "ftv", "lwt", "ptl"), factors = c("ht", "race", "smoke", "ui"))
x <- run_app("lineal", bw, vL, list(poly = "age=2; lwt=4", cv_rep = TRUE))
ck <- setNames(x$assumptions$checks, sapply(x$assumptions$checks, `[[`, "kind"))
st <- x$res$stats

# VIF del modelo completo (cap. 48.4 y 49.3): GVIF y GVIF^(1/(2*Df))
gv <- c(age = "1.155069", ftv = "1.077053", ht = "1.087725", lwt = "1.252135", ptl = "1.125009", race = "1.335314", smoke = "1.207031", ui = "1.087875")
gv2 <- c(age = "1.074741", ftv = "1.037811", ht = "1.042941", lwt = "1.118988", ptl = "1.060665", race = "1.074969", smoke = "1.098650", ui = "1.043012")
for (k in seq_along(gv)) { nm <- names(gv)[k]
  chk(A, paste("VIF (GVIF)", nm), gv[[k]], ck$vif$vif[[nm]], 9245 + k)
  chk(A, paste("VIF (GVIF)", nm, "— cap. 49"), gv[[k]], ck$vif$vif[[nm]], 9689 + k)
  na(A, paste("GVIF^(1/(2Df))", nm), gv2[[k]], 9245 + k, "la app no lo muestra (solo GVIF)") }

# AIC/BIC de modelos (cap. 48.9.1): sin age = modS1 (bloque de age), con age+age² = modS
bl <- function(col) tv(x, "jerarquía", col)
chk(A, "AIC sin age (mod1)", "2991.15319687154", bl("AIC sin"), 9301)
chk(A, "BIC sin age (mod1)", "3017.08717299202", bl("BIC sin"), 9302)
na(A, "AIC con age lineal (mod2)", "2992.89188238441", 9304, "la app no lo muestra (no ajusta el modelo con age solo lineal)")
na(A, "BIC con age lineal (mod2)", "3022.06760551995", 9305, "la app no lo muestra (no ajusta el modelo con age solo lineal)")
chk(A, "AIC con age+age² (mod3)", "2989.01581985068", bl("AIC con"), 9307)
chk(A, "BIC con age+age² (mod3)", "3021.43329000127", bl("BIC con"), 9308)

# Modelo completo (cap. 49.2)
RE0 <- "Modelo inicial \\(todos"
coefs(A, x, RE0, cdf(c("(Intercept)", "2927.962", "312.904", "9.357", "<2e-16"),
                    c("age", "-3.570", "9.620", "-0.371", "0.711012"),
                    c("ftv", "-14.058", "46.468", "-0.303", "0.762598"),
                    c("ht1", "-592.827", "202.321", "-2.930", "0.003830"),
                    c("lwt", "4.354", "1.736", "2.509", "0.013007"),
                    c("ptl", "-48.402", "101.972", "-0.475", "0.635607"),
                    c("race2", "-488.428", "149.985", "-3.257", "0.001349"),
                    c("race3", "-355.077", "114.753", "-3.094", "0.002290"),
                    c("smoke1", "-352.045", "106.476", "-3.306", "0.001142"),
                    c("ui1", "-516.081", "138.885", "-3.716", "0.000271")), 9626, "inicial ")
for (q in 1:5) na(A, paste("Inicial: residuo", c("Min", "1Q", "Mediana", "3Q", "Max")[q]), c("-1825.26", "-435.21", "55.91", "473.46", "1701.20")[q], 9622)
na(A, "Inicial: error estándar residual", "650.3", 9639)
chk(A, "Inicial: gl residual", "179", tv(x, "Modelo inicial: ajuste", "gl2"), 9639)
chk(A, "Inicial: R²", "0.2427", tv(x, "Modelo inicial: ajuste", "R²"), 9640)
chk(A, "Inicial: R² ajustado", "0.2047", tv(x, "Modelo inicial: ajuste", "R² ajustado"), 9640)
chk(A, "Inicial: F", "6.376", tv(x, "Modelo inicial: ajuste", "F"), 9641)
chk(A, "Inicial: gl1", "9", tv(x, "Modelo inicial: ajuste", "gl1"), 9641)
chk(A, "Inicial: p (F)", "7.891e-08", tv(x, "Modelo inicial: ajuste", "p"), 9641)
ci0 <- rbind(c("(Intercept)", "2310.5062595", "3545.417614"), c("age", "-22.5535895", "15.413721"), c("ftv", "-105.7536827", "77.637574"),
             c("ht1", "-992.0688887", "-193.586000"), c("lwt", "0.9291722", "7.778853"), c("ptl", "-249.6231426", "152.819074"),
             c("race2", "-784.3928361", "-192.462241"), c("race3", "-581.5204622", "-128.633751"), c("smoke1", "-562.1550309", "-141.934036"), c("ui1", "-790.1441980", "-242.017757"))
for (i in seq_len(nrow(ci0))) { chk(A, paste("Inicial IC inf", ci0[i, 1]), ci0[i, 2], tv(x, RE0, "IC 95% inf", ci0[i, 1]), 9647 + i)
  chk(A, paste("Inicial IC sup", ci0[i, 1]), ci0[i, 3], tv(x, RE0, "IC 95% sup", ci0[i, 1]), 9647 + i) }

# Diagnóstico del modelo completo (cap. 49.3)
cook_top <- sub("\\).*", "", sub(".*observaciones con mayor valor: ", "", ck$cook$text))
chk(A, "Cook: 3 observaciones con mayor valor", "11, 188, 226", paste(sort(as.integer(strsplit(cook_top, ", ")[[1]])), collapse = ", "), 9681, "conjunto (orden ascendente)")
chk(A, "Cook máxima < 1 (sin casos problemáticos)", "1", as.numeric(ck$cook$max_cook < 1), 9680, "regla del libro: Cook > 1 merece revisión")
na(A, "D-W modelo inicial: autocorrelación", "0.6973801", 9685, "la app no lo muestra (el check del modelo inicial solo da D-W y p)")
chk(A, "D-W modelo inicial: estadístico", "0.6049808", ck$dw$dw, 9685)
chk(A, "D-W modelo inicial: p", "0", ck$dw$p, 9685)
dfr <- c("188", "197", "202", "210", "226")
chk(A, "DFFITS > 0.5: observaciones", paste(dfr, collapse = ", "), paste(tv(x, "DFFITS", "Fila", 1:5), collapse = ", "), 9701)
dfv <- rbind(c(25, 95, 1, 1, 3, 0, 1, 0, 3637), c(19, 184, 1, 1, 0, 1, 0, 0, 3756), c(25, 241, 2, 0, 0, 1, 0, 0, 3790),
             c(33, 117, 1, 0, 0, 0, 1, 1, 3912), c(45, 123, 1, 0, 0, 0, 0, 1, 4990))
cols <- c("age", "lwt", "race", "smoke", "ptl", "ht", "ui", "ftv", "bwt")
for (i in 1:5) for (j in seq_along(cols)) chk(A, paste0("DFFITS fila ", dfr[i], ": ", cols[j]), as.character(dfv[i, j]), as.numeric(as.character(tv(x, "DFFITS", cols[j], i))), 9701 + i)

# Valores únicos y GAM (cap. 49.4)
uq <- c(age = 24, ftv = 6, lwt = 75, ptl = 4); uql <- c(age = 9715, ftv = 9717, lwt = 9719, ptl = 9721)
for (nm in names(uq)) chk(A, paste("Valores únicos", nm), as.character(uq[[nm]]), tv(x, "Valores únicos", "Valores únicos", nm), uql[[nm]])
chk(A, "GAM s(age) Npar F", "3.2043", tv(x, "GAM", "Npar F", "s(age)"), 9759)
chk(A, "GAM s(age) Npar Df", "3", tv(x, "GAM", "Npar Df", "s(age)"), 9759)
chk(A, "GAM s(age) p", "0.02462", tv(x, "GAM", "Pr(F)", "s(age)"), 9759)
chk(A, "GAM s(lwt) Npar F", "2.5131", tv(x, "GAM", "Npar F", "s(lwt)"), 9762)
chk(A, "GAM s(lwt) Npar Df", "3", tv(x, "GAM", "Npar Df", "s(lwt)"), 9762)
chk(A, "GAM s(lwt) p", "0.06016", tv(x, "GAM", "Pr(F)", "s(lwt)"), 9762)
na(A, "GAM: dispersión", "399126", 9734); na(A, "GAM: devianza nula", "99969656", 9736)
na(A, "GAM: devianza residual", "69048781", 9737); na(A, "GAM: AIC", "2991.18", 9738)
gp <- rbind(c("s(age)", "2.2864", "0.1323343"), c("ftv", "0.1594", "0.6902018"), c("ht", "5.4272", "0.0209803"), c("s(lwt)", "11.0282", "0.0010946"),
            c("ptl", "3.0478", "0.0826199"), c("race", "5.5425", "0.0046439"), c("smoke", "12.5349", "0.0005133"), c("ui", "14.8865", "0.0001609"))
for (i in seq_len(nrow(gp))) { na(A, paste("GAM efectos paramétricos F", gp[i, 1]), gp[i, 2], 9743 + i); na(A, paste("GAM efectos paramétricos p", gp[i, 1]), gp[i, 3], 9743 + i) }
# R² de las tendencias polinómicas bivariadas (leyenda del gráfico de la app)
leg <- tryCatch(capture_legends(x$spec$plots$poly$fn(bw, vL, x$res, list())), error = function(e) list())
lg <- function(var, k) { i <- match(var, c("age", "lwt")); if (length(leg) < i) return(NA); rx(leg[[i]][k], "R² = ([0-9.]+)") }
chk(A, "R² bivariado age lineal", "0.008", lg("age", 1), 9791, "leyenda del gráfico «Tendencias polinómicas»")
chk(A, "R² bivariado age cuadrática", "0.04", lg("age", 2), 9791, "leyenda del gráfico «Tendencias polinómicas»")
chk(A, "R² bivariado lwt lineal", "0.034", lg("lwt", 1), 9810, "leyenda del gráfico «Tendencias polinómicas»")
chk(A, "R² bivariado lwt cuártica", "0.069", lg("lwt", 4), 9810, "leyenda del gráfico «Tendencias polinómicas»")

# Modelo centrado / modS (cap. 49.5–49.6.1): coeficientes no mostrados por la app
mc <- cdf(c("(Intercept)", "3305.622", "96.340", "34.312", "<2e-16"), c("age_c", "-15.133", "10.210", "-1.482", "0.140024"),
          c("I(age_c^2)", "2.788", "1.169", "2.384", "0.018159"), c("ht1", "-563.371", "198.063", "-2.844", "0.004965"),
          c("lwt_c", "4.470", "1.685", "2.652", "0.008718"), c("race2", "-506.767", "147.450", "-3.437", "0.000731"),
          c("race3", "-338.877", "112.657", "-3.008", "0.003006"), c("smoke1", "-348.152", "102.842", "-3.385", "0.000873"),
          c("ui1", "-519.553", "133.426", "-3.894", "0.000139"))
for (i in seq_len(nrow(mc))) for (j in 2:5) na(A, paste("modS (centrado)", mc$term[i], c("", "b", "EE", "t", "p")[j]), mc[i, j], 9830 + i,
                                               "la app no lo muestra (solo Anova tipo II del modelo simplificado; sus coeficientes no se tabulan)")
mF <- cdf(c("(Intercept)", "4582.643", "756.723", "6.056", "7.95e-09"), c("age", "-144.713", "59.457", "-2.434", "0.015914"))
for (i in 1:2) for (j in 2:5) na(A, paste("modF (sin centrar)", mF$term[i], c("", "b", "EE", "t", "p")[j]), mF[i, j], 9993 + i,
                                 "la app no lo muestra (coeficientes del modelo con age+age² en escala original)")
na(A, "modS: error estándar residual", "639.1", 9843); na(A, "modS: R²", "0.2647", 9844)
chk(A, "modS: R² ajustado", "0.232", tv(x, "candidatos", "R² ajustado", "Eliminación hacia atrás"), 9844)
chk(A, "modS: F", "8.098", tv(x, "vs. modelo basal", "F", 2), 9845)
chk(A, "modS: p (F)", "2.542e-09", tv(x, "vs. modelo basal", "p", 2), 9845)

# Anova tipo II del modelo saturado (cap. 49.6)
sa <- rbind(c("age_c", "678382", "1", "1.6636", "0.1988144"), c("I(age_c^2)", "2526455", "1", "6.1958", "0.0137395"),
            c("ftv", "240654", "1", "0.5902", "0.4433888"), c("ht", "3496561", "1", "8.5748", "0.0038620"),
            c("lwt_c", "713204", "1", "1.7490", "0.1877239"), c("I(lwt_c^2)", "1215296", "1", "2.9803", "0.0860459"),
            c("I(lwt_c^3)", "37520", "1", "0.0920", "0.7619953"), c("I(lwt_c^4)", "394608", "1", "0.9677", "0.3266069"),
            c("ptl", "28948", "1", "0.0710", "0.7902123"), c("race", "6781627", "2", "8.3155", "0.0003549"),
            c("smoke", "4529203", "1", "11.1073", "0.0010493"), c("ui", "5054204", "1", "12.3947", "0.0005490"))
for (i in seq_len(nrow(sa))) { ln <- 9874 + i
  chk(A, paste("Anova saturado SS", sa[i, 1]), sa[i, 2], tv(x, "Modelo saturado: Anova", "Sum Sq", sa[i, 1]), ln)
  chk(A, paste("Anova saturado Df", sa[i, 1]), sa[i, 3], tv(x, "Modelo saturado: Anova", "Df", sa[i, 1]), ln)
  chk(A, paste("Anova saturado F", sa[i, 1]), sa[i, 4], tv(x, "Modelo saturado: Anova", "F value", sa[i, 1]), ln)
  chk(A, paste("Anova saturado p", sa[i, 1]), sa[i, 5], tv(x, "Modelo saturado: Anova", "Pr(>F)", sa[i, 1]), ln) }
chk(A, "Anova saturado SS residual", "71359733", tv(x, "Modelo saturado: Anova", "Sum Sq", "Residuals"), 9887)
chk(A, "Anova saturado gl residual", "175", tv(x, "Modelo saturado: Anova", "Df", "Residuals"), 9887)
# Pasos de eliminación (cap. 49.6, lista)
stp <- c("ptl", "I(lwt_c^3), I(lwt_c^4)", "ftv", "I(lwt_c^2)")
for (i in 1:4) chk(A, paste("Paso", i, "de eliminación"), stp[i], as.character(tv(x, "Eliminación hacia atrás", "Término(s) eliminado(s)", i)), 9890 + i)
chk(A, "Término no significativo restante (age_c protegido)", "age", as.character(tv(x, "jerarquía", "Variable")), 9895)
# Anova tipo II de modS
ss <- rbind(c("age_c", "897266", "1", "2.1971", "0.1400236"), c("I(age_c^2)", "2321383", "1", "5.6842", "0.0181588"),
            c("ht", "3304159", "1", "8.0906", "0.0049650"), c("lwt_c", "2872143", "1", "7.0327", "0.0087177"),
            c("race", "6614302", "2", "8.0979", "0.0004290"), c("smoke", "4680356", "1", "11.4603", "0.0008728"),
            c("ui", "6192428", "1", "15.1628", "0.0001388"))
for (i in seq_len(nrow(ss))) { ln <- 9901 + i
  chk(A, paste("Anova modS SS", ss[i, 1]), ss[i, 2], tv(x, "Modelo simplificado: Anova", "Sum Sq", ss[i, 1]), ln)
  chk(A, paste("Anova modS Df", ss[i, 1]), ss[i, 3], tv(x, "Modelo simplificado: Anova", "Df", ss[i, 1]), ln)
  chk(A, paste("Anova modS F", ss[i, 1]), ss[i, 4], tv(x, "Modelo simplificado: Anova", "F value", ss[i, 1]), ln)
  chk(A, paste("Anova modS p", ss[i, 1]), ss[i, 5], tv(x, "Modelo simplificado: Anova", "Pr(>F)", ss[i, 1]), ln) }
chk(A, "Anova modS SS residual", "73511202", tv(x, "Modelo simplificado: Anova", "Sum Sq", "Residuals"), 9909)
chk(A, "Anova modS gl residual", "180", tv(x, "Modelo simplificado: Anova", "Df", "Residuals"), 9909)

# modB vs modS (cap. 49.6.1)
VB <- "vs. modelo basal"
chk(A, "modB: gl residual", "188", tv(x, VB, "gl residuales", 1), 9922); chk(A, "modB: RSS", "99969656", tv(x, VB, "RSS", 1), 9922)
chk(A, "modS: gl residual", "180", tv(x, VB, "gl residuales", 2), 9923); chk(A, "modS: RSS", "73511202", tv(x, VB, "RSS", 2), 9923)
na(A, "modB vs modS: Df", "8", 9923); na(A, "modB vs modS: Sum of Sq", "26458454", 9923)
chk(A, "modB vs modS: F", "8.0983", tv(x, VB, "F", 2), 9923); chk(A, "modB vs modS: p", "2.542e-09", tv(x, VB, "p", 2), 9923)
chk(A, "AIC modB", "3031.120", tv(x, VB, "AIC", 1), 9928); chk(A, "AIC modS", "2989.016", tv(x, VB, "AIC", 2), 9929)
chk(A, "BIC modB", "3037.603", tv(x, VB, "BIC", 1), 9932); chk(A, "BIC modS", "3021.433", tv(x, VB, "BIC", 2), 9933)
# modS1 vs modS (efecto de age)
na(A, "modS1: gl residual", "182", 9943); na(A, "modS1: RSS", "75937505", 9943)
chk(A, "modS1 vs modS: Df", "2", bl("gl"), 9944, "columna «gl» de la tabla de bloques")
na(A, "modS1 vs modS: Sum of Sq", "2426303", 9944); na(A, "modS1 vs modS: F", "2.9705", 9944)
chk(A, "modS1 vs modS: p", "0.0538", bl("p"), 9944)
chk(A, "AIC modS1", "2991.153", bl("AIC sin"), 9949); chk(A, "AIC modS (vs modS1)", "2989.016", bl("AIC con"), 9950)
chk(A, "Texto: p de age = 0.054", "0.054", bl("p"), 9951)

# regsubsets (cap. 49.6.2)
RS <- "regsubsets"
chk(A, "regsubsets: índice máx R² ajustado", "46", rx(tv(x, RS, "Modelo (índice)", 1), "\\[(\\d+)\\]"), 10150)
chk(A, "regsubsets: índice mín Cp", "36", rx(tv(x, RS, "Modelo (índice)", 2), "\\[(\\d+)\\]"), 10150)
chk(A, "regsubsets: índice mín BIC", "26", rx(tv(x, RS, "Modelo (índice)", 3), "\\[(\\d+)\\]"), 10150)
chk(A, "regsubsets: modelo máx R² ajustado", normterms("age + I(age^2) + ht + lwt + I(lwt^2) + I(lwt^4) + race + smoke + ui"), normterms(tv(x, RS, "Términos", 1)), 10152, "términos (sin el sufijo _c, orden alfabético)")
chk(A, "regsubsets: modelo mín Cp", normterms("I(age^2) + ht + I(lwt^2) + I(lwt^3) + race + smoke + ui"), normterms(tv(x, RS, "Términos", 2)), 10154, "términos (sin el sufijo _c, orden alfabético)")
chk(A, "regsubsets: modelo mín BIC", normterms("ht + I(lwt^3) + race + smoke + ui"), normterms(tv(x, RS, "Términos", 3)), 10156, "términos (sin el sufijo _c, orden alfabético)")
# Comparación de candidatos (escala original)
CM <- "candidatos"
chk(A, "Candidatos: modF R² ajustado", "0.2320", tv(x, CM, "R² ajustado", "Eliminación hacia atrás"), 10175)
chk(A, "Candidatos: modF AIC", "2989.016", tv(x, CM, "AIC", "Eliminación hacia atrás"), 10176)
na(A, "Candidatos: modR (+lwt²) R² ajustado", "0.2281", 10178, "la app no lo muestra (no ajusta modF + lwt²)"); na(A, "Candidatos: modR AIC", "2990.922", 10179, "la app no lo muestra (no ajusta modF + lwt²)")
na(A, "Candidatos: modC (reducido) R² ajustado", "0.1922", 10181, "la app no lo muestra (modelo ht + race + smoke + ui)"); na(A, "Candidatos: modC AIC", "2995.695", 10182, "la app no lo muestra (modelo ht + race + smoke + ui)")
chk(A, "Candidatos: modBI (menor BIC) R² ajustado", "0.2163", tv(x, CM, "R² ajustado", "regsubsets: BIC (mínimo)"), 10184)
chk(A, "Candidatos: modBI AIC", "2990.936", tv(x, CM, "AIC", "regsubsets: BIC (mínimo)"), 10185)

# anova(mod1, mod2, mod3) (cap. 49.7)
na(A, "anova 3 modelos: mod2 (+lwt²) RSS", "73474705", 10204); na(A, "anova 3 modelos: mod2 F", "0.0889", 10204); na(A, "anova 3 modelos: mod2 p", "0.7659", 10204)
na(A, "anova 3 modelos: mod3 Sum of Sq", "-2462800", 10205); na(A, "anova 3 modelos: mod3 F", "2.0000", 10205); na(A, "anova 3 modelos: mod3 p", "0.1157", 10205)
chk(A, "anova 3 modelos: mod1 RSS", "73511202", tv(x, VB, "RSS", 2), 10203)

# Modelo final (cap. 49.7)
chk(A, "Fórmula del modelo final", normterms("ht + lwt + race + smoke + ui"), normterms(sub(".*~", "", st$formula)), 10209, "términos (orden alfabético); libro: bwt ~ ht + lwt + race + smoke + ui")
REF <- "Modelo final: coeficientes"
coefs(A, x, REF, cdf(c("(Intercept)", "2837.264", "243.676", "11.644", "<2e-16"), c("ht1", "-585.193", "199.644", "-2.931", "0.003810"),
                    c("lwt", "4.242", "1.675", "2.532", "0.012198"), c("race2", "-475.058", "145.603", "-3.263", "0.001318"),
                    c("race3", "-348.150", "112.361", "-3.099", "0.002254"), c("smoke1", "-356.321", "103.444", "-3.445", "0.000710"),
                    c("ui1", "-525.524", "134.675", "-3.902", "0.000134")), 10221, "final ")
for (q in 1:5) na(A, paste("Final: residuo", c("Min", "1Q", "Mediana", "3Q", "Max")[q]), c("-1842.14", "-433.19", "67.09", "459.21", "1631.03")[q], 10217)
na(A, "Final: error estándar residual", "645.9", 10231)
chk(A, "Final: R²", "0.2404", st$R2, 10232); chk(A, "Final: R² ajustado", "0.2154", st$aR2, 10232)
chk(A, "Final: F", "9.6", st$F, 10233); chk(A, "Final: gl1", "6", st$df1, 10233); chk(A, "Final: gl2", "182", st$df2, 10233)
chk(A, "Final: p (F)", "3.601e-09", st$p, 10233)
ciF <- rbind(c("(Intercept)", "2356.4706569", "3318.057183"), c("ht1", "-979.1080814", "-191.278160"), c("lwt", "0.9358509", "7.547249"),
             c("race2", "-762.3440159", "-187.771193"), c("race3", "-569.8476393", "-126.453123"), c("smoke1", "-560.4237850", "-152.218115"), c("ui1", "-791.2496587", "-259.798136"))
for (i in seq_len(nrow(ciF))) { chk(A, paste("Final IC inf", ciF[i, 1]), ciF[i, 2], tv(x, REF, "IC 95% inf", ciF[i, 1]), 10235 + i)
  chk(A, paste("Final IC sup", ciF[i, 1]), ciF[i, 3], tv(x, REF, "IC 95% sup", ciF[i, 1]), 10235 + i) }
chk(A, "D-W final: autocorrelación", "0.7007329", tv(x, "Durbin", "Autocorrelación"), 10249)
chk(A, "D-W final: estadístico", "0.5984074", tv(x, "Durbin", "D-W"), 10249)
chk(A, "D-W final: p", "0", tv(x, "Durbin", "p"), 10249)

# Validación (cap. 49.9) — set.seed(123) como el libro
SP <- "partición 80/20"
vv <- rbind(c("0.3769", "707.8", "546.9", "0.2414"), c("0.3326", "727.3", "552.8", "0.2480")); mm <- c("R2", "RMSE", "MAE", "PER")
for (i in 1:2) for (j in 1:4) chk(A, paste0("Partición 80/20 ", c("modF (con age)", "modNF1 (final)")[i], ": ", mm[j]), vv[i, j], tv(x, SP, mm[j], i), c(10340, 10345)[i] + j - 1, "set.seed(123), como el libro")
CR <- "repetida"; cr <- rbind(c("657.294", "0.218", "538.983"), c("652.822", "0.226", "532.284")); mc2 <- c("RMSE", "Rsquared", "MAE")
for (i in 1:2) for (j in 1:3) chk(A, paste0("CV repetida 10×20 ", c("con age + age²", "sin age")[i], ": ", mc2[j]), cr[i, j], tv(x, CR, mc2[j], i), c(10373, 10377)[i] + j - 1, "set.seed(123), como el libro")
chk(A, "CV 10-fold: RMSE", "649.5529", tv(x, "10-fold", "RMSE"), 10398, "set.seed(123), como el libro")
chk(A, "CV 10-fold: R²", "0.2192741", tv(x, "10-fold", "Rsquared"), 10398, "set.seed(123), como el libro")
chk(A, "CV 10-fold: MAE", "533.5506", tv(x, "10-fold", "MAE"), 10398, "set.seed(123), como el libro")
chk(A, "CV 10-fold: error porcentual", "22.06", tv(x, "10-fold", "Error porcentual (RMSE/media, %)"), 10401, "set.seed(123), como el libro")

# Tamaño del efecto (cap. 49.10)
et <- rbind(c("ht", "0.03585927", "0.04507964"), c("lwt", "0.02675041", "0.03401820"), c("race", "0.06632135", "0.08029930"),
            c("smoke", "0.04952136", "0.06120345"), c("ui", "0.06355146", "0.07720454"))
for (i in 1:5) { chk(A, paste("η²", et[i, 1]), et[i, 2], tv(x, "etaSquared", "η²", et[i, 1]), 10412 + i)
  chk(A, paste("η² parcial", et[i, 1]), et[i, 3], tv(x, "etaSquared", "η² parcial", et[i, 1]), 10412 + i) }
chk(A, "Suma de η² (tipo 1)", "0.2404", st$eta_sum, 10420); chk(A, "R² del modelo (η²)", "0.2404", st$R2, 10422)

# Predicción (cap. 49.11)
pr <- rbind(c("85", "2608.645", "2225.947", "2991.342"), c("95", "2960.238", "2784.837", "3135.639"), c("99", "2417.435", "2145.991", "2688.880"),
            c("109", "2998.100", "2830.096", "3166.103"), c("117", "2841.501", "2551.503", "3131.500"), c("121", "2892.400", "2615.574", "3169.226"))
xr <- rbind(c(182, 0, 2, 0, 1), c(113, 0, 1, 1, 0), c(107, 0, 3, 0, 1), c(120, 0, 3, 0, 0), c(113, 0, 2, 0, 0), c(125, 0, 2, 0, 0))
for (i in 1:6) { ln <- 10440 + i
  chk(A, paste("Predicción: fila", i), pr[i, 1], as.character(tv(x, "Predicción con IC", "Fila", i)), ln)
  for (j in 1:5) chk(A, paste0("Predicción fila ", pr[i, 1], ": ", c("lwt", "ht", "race", "smoke", "ui")[j]), as.character(xr[i, j]),
                     as.numeric(as.character(tv(x, "Predicción con IC", c("lwt", "ht", "race", "smoke", "ui")[j], i))), 10430 + i)
  chk(A, paste("Predicción fila", pr[i, 1], "fit"), pr[i, 2], tv(x, "Predicción con IC", "Predicho", i), ln)
  chk(A, paste("Predicción fila", pr[i, 1], "IC inf"), pr[i, 3], tv(x, "Predicción con IC", "IC inf", i), ln)
  chk(A, paste("Predicción fila", pr[i, 1], "IC sup"), pr[i, 4], tv(x, "Predicción con IC", "IC sup", i), ln) }
chk(A, "Texto: sujeto 85 predicho (g)", "2609", tv(x, "Predicción con IC", "Predicho", 1), 10447)
chk(A, "Texto: sujeto 85 IC inf", "2226", tv(x, "Predicción con IC", "IC inf", 1), 10447)
chk(A, "Texto: sujeto 85 IC sup", "2991", tv(x, "Predicción con IC", "IC sup", 1), 10447)

# Cómo reportar (cap. 49.12)
r <- x$report
chk(A, "Reporte: n observaciones", "189", rx(r, "set de datos de (\\d+) observaciones"), 10451)
chk(A, "Reporte: gl del modelo inicial", "9", rx(r, "parámetros \\((\\d+) gl\\)"), 10451)
chk(A, "Reporte: VIF < 5 (texto)", "se descartó la multicolinealidad (todos los VIF < 5)", regmatches(r, regexpr("se descartó la multicolinealidad \\(todos los VIF < 5\\)", r)), 10451)
chk(A, "Reporte: fórmula final", normterms("ht + lwt + race + smoke + ui"), normterms(sub(".*~", "", sub(".*`(.*)`.*", "\\1", gsub("\n", " ", r)))), 10452, "términos (orden alfabético); libro: bwt ~ ht + lwt + race + smoke + ui")
chk(A, "Reporte: F gl1", "6", rx(r, "F\\((\\d+), (\\d+)\\)", 1), 10453); chk(A, "Reporte: F gl2", "182", rx(r, "F\\((\\d+), (\\d+)\\)", 2), 10453)
chk(A, "Reporte: F", "9.6", rx(r, "F\\(\\d+, \\d+\\) = ([0-9.]+)"), 10453)
chk(A, "Reporte: p del modelo", "<0.001", rp(r, "F\\(\\d+, \\d+\\) = [0-9.]+, "), 10453)
chk(A, "Reporte: R²", "0.24", rx(r, "con R² = ([0-9.]+)"), 10453)
chk(A, "Reporte: R² ajustado", "0.215", rx(r, "R² ajustado = ([0-9.]+)"), 10453)
chk(A, "Reporte: VIF máximo", "1.26", rx(r, "VIF máximo fue ([0-9.]+)"), 10453)
chk(A, "Reporte: D-W", "0.6", rx(r, "D-W = ([0-9.]+)"), 10453)
chk(A, "Reporte: p D-W", "<0.001", rp(r, "D-W = [0-9.]+, "), 10453)
chk(A, "Reporte: Σ η² = R²", "0.24", rx(r, "equivalen al R² del modelo \\(([0-9.]+)\\)"), 10453)
chk(A, "Reporte: CV RMSE", "649.6", rx(r, "RMSE = ([0-9.]+)"), 10454)
chk(A, "Reporte: CV R²", "0.219", rx(r, "RMSE = [0-9.]+, R² = ([0-9.]+)"), 10454)
chk(A, "Reporte: CV MAE", "533.6", rx(r, "\\(MAE\\) = ([0-9.]+)"), 10454)
chk(A, "Reporte: 80/20 R²", "0.333", rx(r, "observados\\) fue de ([0-9.]+)"), 10454)
chk(A, "Reporte: 80/20 MAE", "552.8", rx(r, "con un MAE de ([0-9.]+)"), 10454)
chk(A, "Reporte: % variabilidad", "24", rx(r, "alrededor del (\\d+) %"), 10455)
chk(A, "Reporte: MAE típico", "533.6", rx(r, "cercano a ([0-9.]+) \\(Field"), 10455)

# =============================================================================
# CAP. 50 — Regresión logística: fundamentos
# =============================================================================
A <- "logistica (cap. 50)"
x <- run_app("logistica", bw, list(response = "low", predictors = "smoke"), list(boot = FALSE))
chk(A, "β intercepto (low ~ smoke)", "-1.0870515", tv(x, "Modelo inicial \\(todos", "Estimación", "(Intercept)"), 10559)
chk(A, "β smoke1", "0.7040592", tv(x, "Modelo inicial \\(todos", "Estimación", "smoke"), 10559, "smoke 0/1 como predictor numérico (mismo coeficiente que el factor smoke1)")
chk(A, "OR intercepto", "0.3372093", tv(x, "Modelo inicial \\(todos", "OR", "(Intercept)"), 10565)
chk(A, "OR smoke1", "2.0219436", tv(x, "Modelo inicial \\(todos", "OR", "smoke"), 10565)

# =============================================================================
# CAP. 51 — Regresión logística: flujo completo (birthwt, low)
# Opciones del libro: factores ht, race, smoke, ui; polinomio age=2; bootstrap; 2 sujetos nuevos.
# =============================================================================
A <- "logistica"
vG <- list(response = "low", predictors = c("age", "lwt", "ftv", "ptl"), factors = c("ht", "race", "smoke", "ui"))
x <- run_app("logistica", bw, vG, list(poly = "age=2", boot = TRUE,
             nuevos = "lwt=120, ht=0, race=1, smoke=1, ui=0; lwt=155, ht=0, race=3, smoke=0, ui=0"))
ck <- setNames(x$assumptions$checks, sapply(x$assumptions$checks, `[[`, "kind")); st <- x$res$stats

# EPV (cap. 50.10 y 51.5)
chk(A, "Eventos (bajo peso)", "59", rx(ck$binaria$text, "(\\d+) eventos"), 10584)
na(A, "Predictores (EPV cap. 50)", "7", 10586, "la app no lo muestra: usa los gl del modelo inicial (9), no 7 predictores")
chk(A, "EPV (cap. 50)", "8.4", ck$epv$epv, 10588, "el libro usa 7 predictores; la app divide entre los 9 gl del modelo inicial")
chk(A, "Conteo low = 0", "130", rx(ck$binaria$text, "(\\d+) no eventos"), 11089)
chk(A, "Conteo low = 1", "59", rx(ck$binaria$text, "(\\d+) eventos"), 11089)
chk(A, "gl del modelo inicial (Peduzzi)", "9", ck$epv$df, 10863)
chk(A, "n requerido (Peduzzi, 10 EPV)", "288", ck$epv$nreq, 10867)

# Modelo inicial (cap. 51.3.1)
coefs(A, x, RE0, cdf(c("(Intercept)", "0.480623", "1.196888", "0.402", "0.68801"), c("age", "-0.029549", "0.037031", "-0.798", "0.42489"),
                    c("lwt", "-0.015424", "0.006919", "-2.229", "0.02580"), c("ftv", "0.065302", "0.172394", "0.379", "0.70484"),
                    c("ht1", "1.863303", "0.697533", "2.671", "0.00756"), c("ptl", "0.543337", "0.345403", "1.573", "0.11571"),
                    c("race2", "1.272260", "0.527357", "2.413", "0.01584"), c("race3", "0.880496", "0.440778", "1.998", "0.04576"),
                    c("smoke1", "0.938846", "0.402147", "2.335", "0.01957"), c("ui1", "0.767648", "0.459318", "1.671", "0.09467")), 10745, "inicial ")
na(A, "Inicial: devianza nula", "234.67", 10760); na(A, "Inicial: devianza residual", "201.28", 10761); na(A, "Inicial: AIC", "221.28", 10762)
chk(A, "Texto: β lwt", "-0.0154", tv(x, RE0, "Estimación", "lwt"), 10765)
chk(A, "Texto: OR lwt", "0.985", tv(x, RE0, "OR", "lwt"), 10789)
chk(A, "Texto: EE lwt", "0.0069", tv(x, RE0, "EE", "lwt"), 10790)
orT <- rbind(c("(Intercept)", "1.617", "0.159", "17.769"), c("age", "0.971", "0.901", "1.043"), c("lwt", "0.985", "0.971", "0.998"),
             c("ftv", "1.067", "0.753", "1.490"), c("ht1", "6.445", "1.703", "27.694"), c("ptl", "1.722", "0.884", "3.477"),
             c("race2", "3.569", "1.273", "10.238"), c("race3", "2.412", "1.027", "5.842"), c("smoke1", "2.557", "1.175", "5.743"), c("ui1", "2.155", "0.866", "5.317"))
for (i in seq_len(nrow(orT))) { ln <- 10802 + 4 * (i - 1)
  chk(A, paste("Inicial OR", orT[i, 1]), orT[i, 2], tv(x, RE0, "OR", orT[i, 1]), ln + 1)
  chk(A, paste("Inicial OR IC inf", orT[i, 1]), orT[i, 3], tv(x, RE0, "OR IC inf", orT[i, 1]), ln + 2)
  chk(A, paste("Inicial OR IC sup", orT[i, 1]), orT[i, 4], tv(x, RE0, "OR IC sup", orT[i, 1]), ln + 3) }
# Valor p del modelo inicial (cap. 51.4)
chk(A, "Inicial: χ² razón de verosimilitud", "33.3872", tv(x, "Modelo inicial: razón", "χ²"), 10848)
chk(A, "Inicial: gl", "9", tv(x, "Modelo inicial: razón", "gl"), 10852)
chk(A, "Inicial: p", "0.0001143272", tv(x, "Modelo inicial: razón", "p"), 10856)
chk(A, "Texto: χ² = 33.39", "33.39", tv(x, "Modelo inicial: razón", "χ²"), 10857)
chk(A, "Texto: p = 1.14e-4", "1.14e-04", tv(x, "Modelo inicial: razón", "p"), 10857)

# VIF (cap. 51.6)
gvl <- c(age = "1.100003", lwt = "1.303489", ftv = "1.087144", ht = "1.168419", ptl = "1.087847", race = "1.510253", smoke = "1.348243", ui = "1.063061")
gvl2 <- c(age = "1.048810", lwt = "1.141704", ftv = "1.042662", ht = "1.080934", ptl = "1.042999", race = "1.108568", smoke = "1.161139", ui = "1.031048")
for (k in seq_along(gvl)) { nm <- names(gvl)[k]; chk(A, paste("VIF (GVIF)", nm), gvl[[k]], ck$vif$vif[[nm]], 10874 + k)
  na(A, paste("GVIF^(1/(2Df))", nm), gvl2[[k]], 10874 + k, "la app no lo muestra (solo GVIF)") }

# Anova tipo II saturado (cap. 51.7)
SAL <- "Modelo saturado: Anova"
la <- rbind(c("age_c", "0.2978", "1", "0.585275"), c("I(age_c^2)", "0.8185", "1", "0.365633"), c("ftv", "0.1802", "1", "0.671192"),
            c("ht", "7.2828", "1", "0.006962"), c("lwt", "5.6761", "1", "0.017198"), c("ptl", "2.2703", "1", "0.131878"),
            c("race", "7.6983", "2", "0.021298"), c("smoke", "5.7392", "1", "0.016590"), c("ui", "2.6651", "1", "0.102573"))
for (i in seq_len(nrow(la))) { ln <- 10897 + i
  chk(A, paste("Anova saturado LR χ²", la[i, 1]), la[i, 2], tv(x, SAL, "LR Chisq", la[i, 1]), ln)
  chk(A, paste("Anova saturado Df", la[i, 1]), la[i, 3], tv(x, SAL, "Df", la[i, 1]), ln)
  chk(A, paste("Anova saturado p", la[i, 1]), la[i, 4], tv(x, SAL, "Pr(>Chisq)", la[i, 1]), ln) }
chk(A, "Primer término eliminado (mayor p)", "ftv", as.character(tv(x, "Eliminación hacia atrás", "Término(s) eliminado(s)", 1)), 10909)
fa <- rbind(c("ht", "7.6980", "0.005528"), c("lwt", "6.9493", "0.008385"), c("race", "8.6091", "0.013507"), c("smoke", "7.2931", "0.006922"), c("ui", "4.0308", "0.044677"))
for (i in 1:5) { chk(A, paste("Anova final LR χ²", fa[i, 1]), fa[i, 2], tv(x, "Modelo simplificado: Anova", "LR Chisq", fa[i, 1]), 10915 + i)
  chk(A, paste("Anova final p", fa[i, 1]), fa[i, 3], tv(x, "Modelo simplificado: Anova", "Pr(>Chisq)", fa[i, 1]), 10915 + i) }
chk(A, "Anova final race Df", "2", tv(x, "Modelo simplificado: Anova", "Df", "race"), 10918)
# modB vs LogmF
VB <- "Modelo final vs. modelo basal"
chk(A, "modB: gl residual", "188", tv(x, VB, "gl residuales", 1), 10932); chk(A, "modB: devianza residual", "234.67", tv(x, VB, "Devianza residual", 1), 10932)
chk(A, "LogmF: gl residual", "182", tv(x, VB, "gl residuales", 2), 10933); chk(A, "LogmF: devianza residual", "204.22", tv(x, VB, "Devianza residual", 2), 10933)
chk(A, "modB vs LogmF: Df", "6", st$chdf, 10933); chk(A, "modB vs LogmF: Deviance", "30.455", tv(x, VB, "Devianza", 2), 10933)
chk(A, "modB vs LogmF: p", "3.22e-05", tv(x, VB, "p", 2), 10933)
chk(A, "AIC modB", "236.6720", tv(x, VB, "AIC", 1), 10942); chk(A, "AIC LogmF", "218.2166", tv(x, VB, "AIC", 2), 10943)
chk(A, "BIC modB", "239.9137", tv(x, VB, "BIC", 1), 10946); chk(A, "BIC LogmF", "240.9088", tv(x, VB, "BIC", 2), 10947)
# regsubsets (cap. 51.7.1)
chk(A, "regsubsets: índice máx R² ajustado", "31", rx(tv(x, RS, "Modelo (índice)", 1), "\\[(\\d+)\\]"), 11027)
chk(A, "regsubsets: índice mín Cp", "31", rx(tv(x, RS, "Modelo (índice)", 2), "\\[(\\d+)\\]"), 11027)
chk(A, "regsubsets: índice mín BIC", "11", rx(tv(x, RS, "Modelo (índice)", 3), "\\[(\\d+)\\]"), 11027)
chk(A, "regsubsets: modelo máx R² ajustado", normterms("ht + lwt + ptl + race + smoke + ui"), normterms(tv(x, RS, "Términos", 1)), 11029, "términos (orden alfabético)")
chk(A, "regsubsets: modelo mín Cp", normterms("ht + lwt + ptl + race + smoke + ui"), normterms(tv(x, RS, "Términos", 2)), 11031, "términos (orden alfabético)")
chk(A, "regsubsets: modelo mín BIC", normterms("ht + lwt + ptl"), normterms(tv(x, RS, "Términos", 3)), 11033, "términos (orden alfabético)")
# Pseudo-R² (cap. 51.7.2)
PS <- "pseudo-R²"
chk(A, "R² Hosmer-Lemeshow", "0.13", tv(x, PS, "R² Hosmer-Lemeshow"), 11053)
chk(A, "R² Cox-Snell", "0.149", tv(x, PS, "R² Cox-Snell"), 11054)
chk(A, "R² Nagelkerke", "0.209", tv(x, PS, "R² Nagelkerke"), 11055)
# Bondad de ajuste le Cessie (cap. 51.8.1)
GF <- "le Cessie"
chk(A, "GOF: suma de cuadrados del error", "34.4618771", tv(x, GF, "Suma de cuadrados del error"), 11075)
chk(A, "GOF: valor esperado|H0", "34.2609949", tv(x, GF, "Valor esperado|H0"), 11075)
chk(A, "GOF: DE", "0.3299107", tv(x, GF, "DE"), 11075)
chk(A, "GOF: Z", "0.6088986", tv(x, GF, "Z"), 11077); chk(A, "GOF: P", "0.5425916", tv(x, GF, "P"), 11077)
chk(A, "Texto: z = 0.61", "0.61", st$gof_z, 11078); chk(A, "Texto: p = 0.543", "0.543", st$gof_p, 11078)
# Sobredispersión (cap. 51.8.2)
chk(A, "Sobredispersión: p", "0.316", tv(x, "Sobredispersión", "p"), 11094)
# Funciones de enlace (cap. 51.8.3)
LK <- "Función de enlace"
lk <- rbind(c("logit", "218.2166", "240.9088", "204.22"), c("probit", "217.9177", "240.6099", "203.92"), c("cloglog", "218.6597", "241.3519", "204.66"))
for (i in 1:3) { chk(A, paste("Enlace", lk[i, 1], "AIC"), lk[i, 2], tv(x, LK, "AIC", lk[i, 1]), c(11131, 11137, 11143)[i])
  chk(A, paste("Enlace", lk[i, 1], "BIC"), lk[i, 3], tv(x, LK, "BIC", lk[i, 1]), c(11131, 11137, 11143)[i])
  chk(A, paste("Enlace", lk[i, 1], "devianza residual"), lk[i, 4], tv(x, LK, "Devianza residual", lk[i, 1]), 11151 + i) }
na(A, "anova enlaces: Deviance probit vs logit", "0.29891", 11153); na(A, "anova enlaces: Deviance cloglog vs probit", "-0.74198", 11154)
# ROC, Youden, confusión (cap. 51.9–51.10)
chk(A, "AUC", "0.735", st$auc, 11173)
chk(A, "Umbral óptimo (Youden)", "0.3148913", tv(x, "Discriminación", "Umbral óptimo"), 11181)
chk(A, "Sensibilidad óptima", "0.6949153", tv(x, "Discriminación", "Sensibilidad"), 11181)
chk(A, "Especificidad óptima", "0.7", tv(x, "Discriminación", "Especificidad"), 11181)
chk(A, "Texto: umbral 0.315", "0.315", st$thr, 11182)
chk(A, "Confusión Pred 0 / Obs 0", "91", tv(x, "Matriz de confusión", "Observado 0", 1), 11191)
chk(A, "Confusión Pred 0 / Obs 1", "18", tv(x, "Matriz de confusión", "Observado 1", 1), 11192)
chk(A, "Confusión Pred 1 / Obs 0", "39", tv(x, "Matriz de confusión", "Observado 0", 2), 11194)
chk(A, "Confusión Pred 1 / Obs 1", "41", tv(x, "Matriz de confusión", "Observado 1", 2), 11195)
um <- rbind(c("0.356", "0.900", "0.270"), c("0.864", "0.354", "0.487"), c("0.695", "0.700", "0.302")); um_n <- c("Sensibilidad", "Especificidad", "Error")
for (i in 1:3) for (j in 1:3) chk(A, paste0(c("Umbral 0.5", "Umbral 0.2", "Umbral óptimo")[i], ": ", um_n[j]), um[i, j], tv(x, "Sensibilidad, especificidad y error", um_n[j], i), c(11215, 11219, 11223)[i] + j - 1)
# Hosmer-Lemeshow (cap. 51.11)
chk(A, "Hosmer-Lemeshow χ²", "11.715", tv(x, "Hosmer-Lemeshow", "χ²"), 11238)
chk(A, "Hosmer-Lemeshow gl", "8", tv(x, "Hosmer-Lemeshow", "gl"), 11238)
chk(A, "Hosmer-Lemeshow p", "0.1644", tv(x, "Hosmer-Lemeshow", "p"), 11238)
# Validación 80/20 (cap. 51.12.1) — set.seed(123)
chk(A, "AUC prueba modelo final (80/20)", "0.636", tv(x, "Validación por división", "AUC (prueba)", 1), 11271, "set.seed(123), como el libro")
na(A, "AUC prueba modelo + age (80/20)", "0.64", 11271, "la app no lo muestra: compara con el modelo inicial completo, no con modelo final + age")
# Bootstrap (cap. 51.12.2) — set.seed(42)
vb <- list(Dxy = c("0.4701", "0.5018", "0.4411", "0.0607", "0.4094", "0.2684", "0.5511"), R2 = c("0.2093", "0.2413", "0.1793", "0.0620", "0.1473", "0.0024", "0.2655"),
           Intercept = c("0.0000", "0.0000", "-0.1164", "0.1164", "-0.1164", "-0.5777", "0.4577"), Slope = c("1.0000", "1.0000", "0.8296", "0.1704", "0.8296", "0.4201", "1.2925"),
           Emax = c("0.0000", "0.0000", "0.0924", "-0.0924", "0.0924", "-0.0103", "0.3329"), D = c("0.1558", "0.1839", "0.1312", "0.0527", "0.1031", "-0.0242", "0.2045"),
           U = c("-0.0106", "-0.0106", "0.0135", "-0.0241", "0.0135", "-0.0261", "0.3785"), Q = c("0.1664", "0.1945", "0.1177", "0.0768", "0.0897", "-0.2985", "0.2156"),
           B = c("0.1823", "0.1755", "0.1903", "-0.0149", "0.1972", "0.1678", "0.2326"), g = c("1.1100", "1.2332", "0.9892", "0.2441", "0.8660", "0.0423", "1.3183"),
           gp = c("0.2043", "0.2155", "0.1873", "0.0282", "0.1761", "0.0926", "0.2409"))
vcols <- c("index.orig", "training", "test", "optimism", "index.corrected")
for (k in seq_along(vb)) { ln <- 11281 + k; nm <- names(vb)[k]
  for (j in 1:5) chk(A, paste("validate", nm, vcols[j]), vb[[nm]][j], tv(x, "bootstrap", vcols[j], nm), ln, "set.seed(42), B = 200, como el libro")
  na(A, paste("validate", nm, "Lower"), vb[[nm]][6], ln, "la app no lo muestra (columnas Lower/Upper de rms::validate)")
  na(A, paste("validate", nm, "Upper"), vb[[nm]][7], ln, "la app no lo muestra (columnas Lower/Upper de rms::validate)")
  chk(A, paste("validate", nm, "n"), "200", tv(x, "bootstrap", "n", nm), 11293 + k) }
# Coeficientes y OR por razón de verosimilitud (cap. 51.13)
LR <- "car::Confint"
ci <- rbind(c("(Intercept)", "0.0563", "-1.7339", "1.9597", "1.0579", "0.1766", "7.0971"), c("ht1", "1.8714", "0.5513", "3.3137", "6.4975", "1.7355", "27.4860"),
            c("lwt", "-0.0167", "-0.0309", "-0.0041", "0.9834", "0.9695", "0.9959"), c("race2", "1.3246", "0.3048", "2.3663", "3.7605", "1.3563", "10.6576"),
            c("race3", "0.9262", "0.0965", "1.7939", "2.5249", "1.1013", "6.0127"), c("smoke1", "1.0358", "0.2808", "1.8290", "2.8174", "1.3242", "6.2278"),
            c("ui1", "0.9050", "0.0218", "1.7891", "2.4719", "1.0220", "5.9840"))
lc <- c("logOR", "IC inf", "IC sup", "OR", "OR IC inf", "OR IC sup")
for (i in seq_len(nrow(ci))) for (j in 1:6) chk(A, paste("Confint LR", ci[i, 1], lc[j]), ci[i, j + 1], tv(x, LR, lc[j], ci[i, 1]), (if (j <= 3) 11314 else 11324) + i)
chk(A, "Texto: OR ht ≈ 6.5", "6.5", tv(x, LR, "OR", "ht1"), 11332); chk(A, "Texto: OR lwt 0.983", "0.983", tv(x, LR, "OR", "lwt"), 11332)
# Nomograma y predicciones (cap. 51.14–51.15)
na(A, "Nomograma: predictor lineal (lwt=155, race=3)", "-1.6111", 11341, "la app no lo muestra (solo la probabilidad)")
chk(A, "Nomograma: probabilidad (lwt=155, race=3)", "0.1664", tv(x, "sujetos nuevos", "Probabilidad predicha", 2), 11343, "sujeto nuevo 2 de las opciones")
chk(A, "Texto: probabilidad 0.166", "0.166", tv(x, "sujetos nuevos", "Probabilidad predicha", 2), 11344)
na(A, "Nomograma: puntos lwt = 155", "58.3", 11377, "solo en el gráfico del nomograma"); na(A, "Nomograma: puntos race = 3", "30.8", 11379, "solo en el gráfico del nomograma")
na(A, "Nomograma: puntos totales", "89.1", 11388, "solo en el gráfico del nomograma")
chk(A, "Nuevo paciente (lwt=120, race=1, smoke=1)", "0.2858", tv(x, "sujetos nuevos", "Probabilidad predicha", 1), 11403)
chk(A, "Nuevo paciente (%)", "28.6", 100 * tv(x, "sujetos nuevos", "Probabilidad predicha", 1), 11403)
# Cómo reportar (cap. 51.17)
r <- x$report
chk(A, "Reporte: n variables", "10", rx(r, "incluye (\\d+) variables"), 11418)
chk(A, "Reporte: n observaciones", "189", rx(r, "variables y (\\d+) observaciones"), 11418)
chk(A, "Reporte: fórmula final", normterms("ht + lwt + race + smoke + ui"), normterms(sub(".*~", "", sub(".*modelo final fue `([^`]*)`.*", "\\1", r))), 11418, "términos (orden alfabético); libro: low ~ ht + lwt + race + smoke + ui")
chk(A, "Fórmula del modelo final", normterms("ht + lwt + race + smoke + ui"), normterms(sub(".*~", "", st$formula)), 10925, "términos (orden alfabético); libro: LogmF <- glm(low ~ ht + lwt + race + smoke + ui)")
chk(A, "Reporte: χ²", "30.46", rx(r, "χ² = ([0-9.]+)"), 11418); chk(A, "Reporte: gl", "6", rx(r, "gl = (\\d+)"), 11418)
chk(A, "Reporte: p del modelo", "<0.001", rp(r, "gl = \\d+, "), 11418)
chk(A, "Reporte: pseudo-R² mínimo", "0.13", rx(r, "oscilaron entre ([0-9.]+)"), 11418)
chk(A, "Reporte: pseudo-R² máximo", "0.21", rx(r, "oscilaron entre [0-9.]+ y ([0-9.]+)"), 11418)
chk(A, "Reporte: z le Cessie", "0.61", rx(r, "z = ([0-9.]+)"), 11418); chk(A, "Reporte: p le Cessie", "0.54", rp(r, "z = [0-9.]+, "), 11418)
chk(A, "Reporte: p sobredispersión", "0.32", rp(r, "sobredispersión \\("), 11418)
chk(A, "Reporte: AUC", "0.735", rx(r, "AUC de ([0-9.]+)"), 11418)
chk(A, "Reporte: p Hosmer-Lemeshow", "0.16", rp(r, "Hosmer-Lemeshow \\("), 11418)

# =============================================================================
# CAP. 53 — Kaplan-Meier y log-rank (CRIT1)
# =============================================================================
A <- "km"
crit <- book_data("CRIT1.txt")
x0 <- run_app("km", crit, list(time = "time", event = "status"), list(tiempos = ""))   # sin tiempos: menor y mayor tiempo con evento
gt <- rbind(c(288, 269, 0.996, "0.00371", "0.989", "1.000"), c(409, 235, 0.992, "0.00562", "0.981", "1.000"), c(414, 233, 0.988, "0.00702", "0.974", "1.000"),
            c(476, 212, 0.983, "0.00839", "0.967", "1.000"), c(754, 153, 0.977, "0.01052", "0.956", "0.998"), c(862, 132, 0.969, "0.01278", "0.945", "0.995"),
            c(887, 127, 0.962, "0.01478", "0.933", "0.991"), c(966, 113, 0.953, "0.01692", "0.921", "0.987"), c(999, 107, 0.944, "0.01897", "0.908", "0.982"),
            c(1141, 78, 0.932, "0.02225", "0.890", "0.977"), c(2078, 16, 0.874, "0.06014", "0.764", "1.000"))
kc <- c("Tiempo", "En riesgo", "S(t)", "EE", "IC inf", "IC sup")
for (i in seq_len(nrow(gt))) { chk(A, paste("Global t =", gt[i, 1], "eventos"), "1", tv(x0, "global \\(tiempos con evento", "Eventos", i), 11633 + i)
  for (j in 1:6) chk(A, paste("Global t =", gt[i, 1], kc[j]), as.character(gt[i, j]), tv(x0, "global \\(tiempos con evento", kc[j], i), 11633 + i) }
chk(A, "S(288) (texto)", "0.996", tv(x0, "global en los tiempos", "S(t)", 1), 11646); chk(A, "S(288) IC inf", "0.989", tv(x0, "global en los tiempos", "IC inf", 1), 11646)
chk(A, "S(288) IC sup", "1", tv(x0, "global en los tiempos", "IC sup", 1), 11646)
chk(A, "S(2078) (texto)", "0.874", tv(x0, "global en los tiempos", "S(t)", 2), 11646); chk(A, "S(2078) IC inf", "0.764", tv(x0, "global en los tiempos", "IC inf", 2), 11646)
chk(A, "S(2078) IC sup", "1", tv(x0, "global en los tiempos", "IC sup", 2), 11646)
# S(t) a 1 año (365.25)
x1 <- run_app("km", crit, list(time = "time", event = "status"), list(tiempos = "365.25"))
s1 <- c("En riesgo" = "245", Eventos = "1", "S(t)" = "0.996", EE = "0.00371", "IC inf" = "0.989", "IC sup" = "1")
for (nm in names(s1)) chk(A, paste("S(365.25) global", nm), s1[[nm]], tv(x1, "global en los tiempos", nm, 1), 11662)
chk(A, "S(365) texto", "0.996", tv(x1, "global en los tiempos", "S(t)", 1), 11663)
# log-rank por sexo
xs <- run_app("km", crit, list(time = "time", event = "status", group = "sex"), list(tiempos = "365.25"))
OE <- "observados vs. esperados"
lr <- rbind(c("F", "188", "2", "4.87", "1.69", "3.04"), c("M", "214", "9", "6.13", "1.34", "3.04")); lrc <- c("N", "Observados", "Esperados", "(O-E)^2/E", "(O-E)^2/V")
for (i in 1:2) for (j in 1:5) chk(A, paste0("log-rank sexo ", lr[i, 1], ": ", lrc[j]), lr[i, j + 1], tv(xs, OE, lrc[j], i), 11672 + i)
chk(A, "log-rank sexo χ²", "3", tv(xs, "Prueba log-rank", "χ²"), 11676); chk(A, "log-rank sexo gl", "1", tv(xs, "Prueba log-rank", "gl"), 11676)
chk(A, "log-rank sexo p", "0.08", tv(xs, "Prueba log-rank", "p"), 11676)
# log-rank por bruxismo
xb <- run_app("km", crit, list(time = "time", event = "status", group = "brux"), list(tiempos = "730.5"))
lr <- rbind(c("NO", "317", "4", "8.98", "2.76", "15.1"), c("SI", "85", "7", "2.02", "12.27", "15.1"))
for (i in 1:2) for (j in 1:5) chk(A, paste0("log-rank brux ", lr[i, 1], ": ", lrc[j]), lr[i, j + 1], tv(xb, OE, lrc[j], i), 11695 + i)
chk(A, "log-rank brux χ²", "15.1", tv(xb, "Prueba log-rank", "χ²"), 11699); chk(A, "log-rank brux gl", "1", tv(xb, "Prueba log-rank", "gl"), 11699)
chk(A, "log-rank brux p", "1e-04", tv(xb, "Prueba log-rank", "p"), 11699)
chk(A, "Texto: brux p = 0.0001", "0.0001", tv(xb, "Prueba log-rank", "p"), 11700)
# Tabla por grupo (brux)
bt <- rbind(c("brux=NO", 288, 210, 0.995, "0.00475", "0.986", "1"), c("brux=NO", 476, 169, 0.989, "0.00754", "0.975", "1"),
            c("brux=NO", 966, 96, 0.979, "0.01268", "0.955", "1"), c("brux=NO", 999, 92, 0.968, "0.01641", "0.937", "1"),
            c("brux=SI", 409, 49, "0.980", "0.0202", "0.941", "1.000"), c("brux=SI", 414, 48, 0.959, "0.0283", "0.905", "1.000"),
            c("brux=SI", 754, 29, 0.926, "0.0424", "0.847", "1.000"), c("brux=SI", 862, 24, 0.888, "0.0555", "0.785", "1.000"),
            c("brux=SI", 887, 23, 0.849, "0.0651", "0.730", "0.987"), c("brux=SI", 1141, 12, 0.778, "0.0903", "0.620", "0.977"),
            c("brux=SI", 2078, 3, 0.519, "0.2202", "0.226", "1.000"))
bl_ln <- c(11722:11725, 11729:11735)
for (i in seq_len(nrow(bt))) { chk(A, paste(bt[i, 1], "t =", bt[i, 2], "estrato"), bt[i, 1], as.character(tv(xb, "por grupo \\(tiempos con evento", "Estrato", i)), bl_ln[i])
  for (j in 1:6) chk(A, paste(bt[i, 1], "t =", bt[i, 2], kc[j]), as.character(bt[i, j + 1]), tv(xb, "por grupo \\(tiempos con evento", kc[j], i), bl_ln[i]) }
# S(t) a 730.5 y 2889 por grupo
SG <- "por grupo en los tiempos"
b7 <- rbind(c("126", "2", "0.989", "7.54e-03", "9.75e-01", "1"), c("31", "2", "0.9592", "0.0283", "0.9054", "1.0000")); bc <- c("En riesgo", "Eventos", "S(t)", "EE", "IC inf", "IC sup")
for (i in 1:2) for (j in 1:6) chk(A, paste0("S(730.5) ", c("brux=NO", "brux=SI")[i], ": ", bc[j]), b7[i, j], tv(xb, SG, bc[j], i), c(11743, 11749)[i])
xb2 <- run_app("km", crit, list(time = "time", event = "status", group = "brux"), list(tiempos = "2889"))
b2 <- rbind(c("1", "4", "0.968", "1.64e-02", "9.37e-01", "1"), c("1", "7", "0.519", "0.220", "0.226", "1.000"))
for (i in 1:2) for (j in 1:6) chk(A, paste0("S(2889) ", c("brux=NO", "brux=SI")[i], ": ", bc[j]), b2[i, j], tv(xb2, SG, bc[j], i), c(11758, 11764)[i])
chk(A, "Texto: S(2889) brux=SI", "0.519", tv(xb2, SG, "S(t)", 2), 11767); chk(A, "Texto: S(2889) brux=SI IC inf", "0.226", tv(xb2, SG, "IC inf", 2), 11767)
chk(A, "Texto: S(2889) brux=NO", "0.968", tv(xb2, SG, "S(t)", 1), 11767); chk(A, "Texto: S(2889) brux=NO IC inf", "0.937", tv(xb2, SG, "IC inf", 1), 11767)
chk(A, "Último tiempo observado (censura)", "2920", x0$res$stats$tmax, 11752)
chk(A, "Último evento", "2078", x0$res$stats$last_ev, 11752)
# Cómo reportar (cap. 53.7)
r <- x1$report
chk(A, "Reporte: n pacientes", "402", rx(r, "Kaplan-Meier en (\\d+) sujetos"), 11771)
chk(A, "Reporte: eventos", "11", rx(r, "de los cuales (\\d+) presentaron"), 11771)
chk(A, "Reporte: S a un año", "0.996", rx(r, "fue de ([0-9.]+) \\(IC"), 11771)
chk(A, "Reporte: S a un año IC inf", "0.989", rx(r, "IC 95 %: ([0-9.]+)–"), 11771)
chk(A, "Reporte: S a un año IC sup", "1", rx(r, "IC 95 %: [0-9.]+–([0-9.]+)\\)"), 11771)
chk(A, "Reporte: log-rank sexo χ²", "3", rx(xs$report, "χ² = ([0-9.]+)"), 11771)
chk(A, "Reporte: log-rank sexo p", "0.08", rp(xs$report, "gl = \\d+, "), 11771)
chk(A, "Reporte: log-rank brux χ²", "15.1", rx(xb$report, "χ² = ([0-9.]+)"), 11771)
chk(A, "Reporte: log-rank brux p (texto)", "p = 0.0001", { m <- regmatches(xb$report, regexpr("gl = \\d+, p [<=] [0-9.]+", xb$report)); sub("\\.$", "", sub(".*(p [<=] [0-9.]+)", "\\1", m)) }, 11771,
    "comparación de texto: el valor de la app es p = 0.000103 (véase «log-rank brux p»)")

# =============================================================================
# CAP. 55 — Regresión de Cox (CRIT1)
# =============================================================================
A <- "cox"
vC <- list(time = "time", event = "status", predictors = c("asa", "age", "sex", "brux"))
x <- run_app("cox", crit, vC, list(boot = TRUE, nuevos = "asa=I, age=4, sex=F, brux=NO; asa=II, age=6, sex=F, brux=NO; asa=II, age=5, sex=F, brux=SI"))
ck <- setNames(x$assumptions$checks, sapply(x$assumptions$checks, `[[`, "kind")); st <- x$res$stats
vf <- c(asaII = "1.012333", age = "1.015239", sexM = "1.019302", bruxSI = "1.010003")
for (k in seq_along(vf)) { nm <- names(vf)[k]; av <- ck$vif$vif; chk(A, paste("VIF", nm), vf[[k]], if (nm %in% names(av)) av[[nm]] else av[k], 11904) }
an <- rbind(c("asa", "1.8720", "0.171248"), c("age", "1.6798", "0.194945"), c("sex", "3.1423", "0.076284"), c("brux", "8.7521", "0.003092"))
for (i in 1:4) { chk(A, paste("Anova saturado LR χ²", an[i, 1]), an[i, 2], tv(x, "Modelo saturado: Anova", "LR Chisq", an[i, 1]), 11911 + i)
  chk(A, paste("Anova saturado p", an[i, 1]), an[i, 3], tv(x, "Modelo saturado: Anova", "Pr(>Chisq)", an[i, 1]), 11911 + i) }
HS <- "Modelo saturado: coeficientes"
hs <- rbind(c("asaII", "1.24029", "3.45662", "1.06023", "1.170", "0.24207", "0.2893", "0.4327", "27.614"),
            c("age", "0.11760", "1.12479", "0.08725", "1.348", "0.17771", "0.8891", "0.9480", "1.335"),
            c("sexM", "1.25626", "3.51225", "0.79023", "1.590", "0.11189", "0.2847", "0.7463", "16.528"),
            c("bruxSI", "1.84075", "6.30126", "0.63276", "2.909", "0.00362", "0.1587", "1.8231", "21.779"))
hc <- c("β", "HR", "EE", "z", "p")
for (i in 1:4) { for (j in 1:5) chk(A, paste("Saturado", hs[i, 1], hc[j]), hs[i, j + 1], tv(x, HS, hc[j], hs[i, 1]), 11926 + i)
  na(A, paste("Saturado", hs[i, 1], "exp(-coef)"), hs[i, 7], 11934 + i)
  chk(A, paste("Saturado", hs[i, 1], "HR IC inf"), hs[i, 8], tv(x, HS, "HR IC inf", hs[i, 1]), 11934 + i)
  chk(A, paste("Saturado", hs[i, 1], "HR IC sup"), hs[i, 9], tv(x, HS, "HR IC sup", hs[i, 1]), 11934 + i) }
GL <- "Pruebas globales"
chk(A, "Saturado: n", "402", st$n, 11924); chk(A, "Saturado: eventos", "11", st$events, 11924)
chk(A, "Saturado: concordancia", "0.796", tv(x, GL, "Concordancia", 1), 11940); na(A, "Saturado: EE concordancia", "0.062", 11940)
chk(A, "Saturado: LR", "16.97", tv(x, GL, "LR χ²", 1), 11941); chk(A, "Saturado: LR gl", "4", tv(x, GL, "gl", 1), 11941); chk(A, "Saturado: LR p", "0.002", tv(x, GL, "p (LR)", 1), 11941)
chk(A, "Saturado: Wald", "14.65", tv(x, GL, "Wald", 1), 11942); chk(A, "Saturado: Wald p", "0.005", tv(x, GL, "p (Wald)", 1), 11942)
chk(A, "Saturado: Score", "19.75", tv(x, GL, "Score (log-rank)", 1), 11943); chk(A, "Saturado: Score p", "6e-04", tv(x, GL, "p (Score)", 1), 11943)
chk(A, "Texto: HR brux 6.3", "6.3", tv(x, HS, "HR", "bruxSI"), 11947); chk(A, "Texto: IC 1.82", "1.82", tv(x, HS, "HR IC inf", "bruxSI"), 11947)
chk(A, "Texto: IC 21.78", "21.78", tv(x, HS, "HR IC sup", "bruxSI"), 11947)
# Simplificación (cap. 55.4)
SS <- "Simplificación: covariables"
chk(A, "Paso 1 eliminado", "age", as.character(tv(x, SS, "Covariable eliminada", 1)), 11953)
chk(A, "Paso 2 eliminado", "asa", as.character(tv(x, SS, "Covariable eliminada", 2)), 11953)
chk(A, "modS2: LR χ² sex", "2.2042", tv(x, SS, "LR χ²", 3), 11958); chk(A, "modS2: p sex", "0.137636", tv(x, SS, "p", 3), 11958)
na(A, "modS2: LR χ² brux", "9.8299", 11959, "la app no lo muestra (los pasos solo listan la covariable eliminada)"); na(A, "modS2: p brux", "0.001717", 11959, "la app no lo muestra (los pasos solo listan la covariable eliminada)")
chk(A, "Paso 3 eliminado", "sex", as.character(tv(x, SS, "Covariable eliminada", 3)), 11962)
HF <- "Modelo final: coeficientes"
chk(A, "Final bruxSI β", "2.0637", tv(x, HF, "β", "bruxSI"), 11971); chk(A, "Final bruxSI HR", "7.8749", tv(x, HF, "HR", "bruxSI"), 11971)
chk(A, "Final bruxSI EE", "0.6289", tv(x, HF, "EE", "bruxSI"), 11971); chk(A, "Final bruxSI z", "3.281", tv(x, HF, "z", "bruxSI"), 11971)
chk(A, "Final bruxSI p", "0.00103", tv(x, HF, "p", "bruxSI"), 11971); na(A, "Final bruxSI exp(-coef)", "0.127", 11976)
chk(A, "Final HR IC inf", "2.296", tv(x, HF, "HR IC inf", "bruxSI"), 11976); chk(A, "Final HR IC sup", "27.01", tv(x, HF, "HR IC sup", "bruxSI"), 11976)
chk(A, "Final: concordancia", "0.696", tv(x, GL, "Concordancia", 2), 11978); na(A, "Final: EE concordancia", "0.084", 11978)
chk(A, "Final: LR", "10.99", tv(x, GL, "LR χ²", 2), 11979); chk(A, "Final: LR p", "9e-04", tv(x, GL, "p (LR)", 2), 11979)
chk(A, "Final: Wald", "10.77", tv(x, GL, "Wald", 2), 11980); chk(A, "Final: Wald p", "0.001", tv(x, GL, "p (Wald)", 2), 11980)
chk(A, "Final: Score", "15.09", tv(x, GL, "Score (log-rank)", 2), 11981); chk(A, "Final: Score p", "1e-04", tv(x, GL, "p (Score)", 2), 11981)
chk(A, "Texto: HR final 7.87", "7.87", tv(x, HF, "HR", "bruxSI"), 11982); chk(A, "Texto: IC 2.3", "2.3", tv(x, HF, "HR IC inf", "bruxSI"), 11982)
CMP <- "Modelo saturado vs. modelo final"
chk(A, "anova: loglik saturado", "-44.666", tv(x, CMP, "loglik", "Saturado"), 11989); chk(A, "anova: loglik final", "-47.656", tv(x, CMP, "loglik", "Final"), 11990)
chk(A, "anova: χ²", "5.9793", tv(x, CMP, "χ²", "Saturado"), 11990); chk(A, "anova: gl", "3", tv(x, CMP, "gl", "Saturado"), 11990)
chk(A, "anova: p", "0.1126", tv(x, CMP, "p", "Saturado"), 11990)
# cph (cap. 55.5)
IA <- "Índices de ajuste"
cp <- rbind(c("402", "11", "16.97", "4", "0.0020", "19.75", "0.0006", "0.178", "0.592"), c("402", "11", "10.99", "1", "0.0009", "15.09", "0.0001", "0.116", "0.392"))
cpc <- c("Obs", "Eventos", "LR χ²", "gl", "p", "Score χ²", "p (Score)", "R2", "Dxy")
for (i in 1:2) for (j in seq_along(cpc)) chk(A, paste0("cph ", c("saturado", "final")[i], ": ", cpc[j]), cp[i, j], tv(x, IA, cpc[j], i), c(12007, 12028)[i] + c(0, 1, 0, 1, 2, 3, 4, 0, 3)[j])
na(A, "cph saturado: Center", "2.6988", 12009); na(A, "cph saturado: R2(4,402)", "0.032", 12008); na(A, "cph saturado: R2(4,11)", "0.692", 12009)
na(A, "cph final: Center", "0.4363", 12030); na(A, "cph final: R2(1,402)", "0.025", 12029); na(A, "cph final: R2(1,11)", "0.597", 12030)
cpS <- rbind(c("asaII", "1.2403", "1.0602", "1.17", "0.2421"), c("age", "0.1176", "0.0872", "1.35", "0.1777"), c("sexM", "1.2563", "0.7902", "1.59", "0.1119"), c("bruxSI", "1.8407", "0.6328", "2.91", "0.0036"))
for (i in 1:4) for (j in 1:4) chk(A, paste("cph saturado", cpS[i, 1], c("Coef", "S.E.", "Wald Z", "p")[j]), cpS[i, j + 1], tv(x, HS, c("β", "EE", "z", "p")[j], cpS[i, 1]), 12013 + i)
chk(A, "cph final bruxSI Coef", "2.0637", tv(x, HF, "β", "bruxSI"), 12035); chk(A, "cph final bruxSI S.E.", "0.6289", tv(x, HF, "EE", "bruxSI"), 12035)
chk(A, "cph final bruxSI Wald Z", "3.28", tv(x, HF, "z", "bruxSI"), 12035); chk(A, "cph final bruxSI p", "0.0010", tv(x, HF, "p", "bruxSI"), 12035)
chk(A, "Texto: R² saturado", "0.178", tv(x, IA, "R2", 1), 12018); chk(A, "Texto: Dxy saturado", "0.592", tv(x, IA, "Dxy", 1), 12018)
chk(A, "Texto: R² final", "0.116", tv(x, IA, "R2", 2), 12036); chk(A, "Texto: Dxy final", "0.392", tv(x, IA, "Dxy", 2), 12036)
# validate (cap. 55.6) — set.seed(1)
vc <- list(Dxy = c("0.5919", "0.6531", "0.5356", "0.1175", "0.4745", "0.0811", "0.7668"), R2 = c("0.1779", "0.2273", "0.1280", "0.0993", "0.0785", "-0.1310", "0.2738"),
           Slope = c("1.0000", "1.0000", "0.6476", "0.3524", "0.6476", "-0.1042", "1.5066"), D = c("0.1502", "0.1962", "0.1050", "0.0913", "0.0589", "-0.1460", "0.2370"),
           U = c("-0.0188", "-0.0215", "0.0922", "-0.1137", "0.0949", "-0.0909", "1.1151"), Q = c("0.1690", "0.2177", "0.0128", "0.2050", "-0.0360", "-1.0321", "0.3071"),
           g = c("1.3536", "2.9260", "1.3770", "1.5490", "-0.1955", "-4.4841", "2.5660"))
for (k in seq_along(vc)) { ln <- 12042 + k; nm <- names(vc)[k]
  for (j in 1:5) chk(A, paste("validate", nm, vcols[j]), vc[[nm]][j], tv(x, "bootstrap", vcols[j], nm), ln, "set.seed(1), B = 200, como el libro")
  na(A, paste("validate", nm, "Lower"), vc[[nm]][6], ln, "la app no lo muestra (columnas Lower/Upper de rms::validate)")
  na(A, paste("validate", nm, "Upper"), vc[[nm]][7], ln, "la app no lo muestra (columnas Lower/Upper de rms::validate)")
  chk(A, paste("validate", nm, "n"), "199", tv(x, "bootstrap", "n", nm), ln) }
# cox.zph (cap. 55.7)
ZP <- "cox.zph"
chk(A, "zph brux χ²", "0.338", tv(x, ZP, "χ²", "brux"), 12057); chk(A, "zph brux gl", "1", tv(x, ZP, "gl", "brux"), 12057); chk(A, "zph brux p", "0.56", tv(x, ZP, "p", "brux"), 12057)
chk(A, "zph GLOBAL χ²", "0.338", tv(x, ZP, "χ²", "GLOBAL"), 12058); chk(A, "zph GLOBAL gl", "1", tv(x, ZP, "gl", "GLOBAL"), 12058); chk(A, "zph GLOBAL p", "0.56", tv(x, ZP, "p", "GLOBAL"), 12058)
chk(A, "Texto: p global 0.561", "0.561", tv(x, ZP, "p", "GLOBAL"), 12061)
chk(A, "zph p global (check de supuestos)", "0.561", ck$ph$p, 12061)
# Pacientes nuevos (cap. 55.8)
PN <- "pacientes nuevos"
rk <- rbind(c("0.741", "0.192", "0.365", "1.116"), c("3.238", "1.915", "0.000", "6.992"), c("18.142", "5.163", "8.022", "28.262")); rc <- c("Riesgo", "EE", "IC inf", "IC sup")
for (i in 1:3) for (j in 1:4) chk(A, paste0("Paciente ", i, ": ", rc[j]), rk[i, j], tv(x, PN, rc[j], i), 12092 + 4 * (i - 1) + j - 1)
chk(A, "Texto: riesgo paciente 1", "0.74", tv(x, PN, "Riesgo", 1), 12104); chk(A, "Texto: riesgo paciente 3", "18.14", tv(x, PN, "Riesgo", 3), 12104)
# Cómo reportar (cap. 55.9)
r <- x$report
chk(A, "Reporte: n pacientes", "402", rx(r, "en (\\d+) sujetos"), 12118); chk(A, "Reporte: eventos", "11", rx(r, "de los cuales (\\d+) presentaron"), 12118)
chk(A, "Reporte: VIF < 5 (texto)", "VIF < 5 en todas las covariables", regmatches(r, regexpr("VIF < 5 en todas las covariables", r)), 12118)
chk(A, "Reporte: HR", "7.87", rx(r, "HR = ([0-9.]+)"), 12118); chk(A, "Reporte: IC inf", "2.3", rx(r, "IC 95 %: ([0-9.]+)–"), 12118)
chk(A, "Reporte: IC sup", "27.01", rx(r, "IC 95 %: [0-9.]+–([0-9.]+)"), 12118); chk(A, "Reporte: p HR", "0.001", rp(r, "–[0-9.]+; "), 12118)
chk(A, "Reporte: p global Schoenfeld", "0.561", rx(r, "p global = ([0-9.]+)"), 12118)

# =============================================================================
write_results(file.path(Sys.getenv("VERIFY_DIR"), "results_reg.csv"))
cat("\nTotal:", nrow(RESULTS), " coinciden:", sum(RESULTS$match),
    " no coinciden (con valor de la app):", sum(!RESULTS$match & RESULTS$note != NOAPP & !grepl("^la app no lo muestra|^solo en", RESULTS$note)),
    " la app no lo muestra:", sum(grepl("^la app no lo muestra|^solo en", RESULTS$note)), "\n")
options(width = 250)
mm <- RESULTS[!RESULTS$match & !grepl("^la app no lo muestra|^solo en", RESULTS$note), c("analysis", "item", "book", "app", "line", "note")]
print(mm, right = FALSE)
