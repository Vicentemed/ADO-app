# =============================================================================
# 03_assumptions.R — Motor de supuestos (cap. "Supuestos del análisis paramétrico")
# Reglas del libro:
#  - Shapiro-Wilk: p > 0.05 → sin evidencia de desviación de la normalidad
#    (no equivale a demostrarla); p < 0.05 → supuesto violado.
#  - n > 50: Shapiro tiende a ser significativo; apoyarse en Q-Q y sesgo/curtosis.
#  - |skew.2SE| o |kurt.2SE| > 1.96 → desviación significativa (α = 0.05).
#  - Levene (center = median): p > 0.05 → varianzas homogéneas.
# =============================================================================

check_normality <- function(x, label = "normalidad") {
  x <- as.numeric(x); x <- x[!is.na(x)]
  n <- length(x)
  sw <- if (n >= 3 && n <= 5000) tryCatch(shapiro.test(x), error = function(e) NULL) else NULL
  sk <- tryCatch(pastecs::stat.desc(x, basic = FALSE, norm = TRUE), error = function(e) NULL)
  W <- if (!is.null(sw)) unname(sw$statistic) else NA
  p <- if (!is.null(sw)) sw$p.value else NA
  skew2 <- if (!is.null(sk)) unname(sk["skew.2SE"]) else NA
  kurt2 <- if (!is.null(sk)) unname(sk["kurt.2SE"]) else NA
  met <- if (is.na(p)) NA else p > ALPHA
  txt <- if (is.na(p)) {
    paste0("No fue posible aplicar Shapiro-Wilk a ", label, " (n = ", n, ").")
  } else if (met) {
    paste0("La prueba de Shapiro-Wilk sobre ", label, " (W = ", fmt_num(W, 3), ", p ", fmt_p(p),
           ") no muestra evidencia de una distribución diferente a la normal (lo cual no equivale a demostrarla).")
  } else {
    paste0("La prueba de Shapiro-Wilk sobre ", label, " (W = ", fmt_num(W, 3), ", p ", fmt_p(p),
           ") indica una desviación significativa de la normalidad: el supuesto se considera violado.")
  }
  if (n > 50 && !is.na(p))
    txt <- paste0(txt, " Con n > 50 la prueba tiende a ser significativa ante desviaciones pequeñas; ",
                  "conviene apoyarse en el gráfico Q-Q y en los estadísticos de sesgo (skew.2SE = ",
                  fmt_num(skew2), ") y curtosis (kurt.2SE = ", fmt_num(kurt2), "): |valor| > 1.96 indica desviación significativa.")
  check_obj("normalidad", paste0("Normalidad (", label, ")"), met, txt,
            W = W, p = p, n = n, skew2SE = skew2, kurt2SE = kurt2, large_n = n > 50)
}

check_homogeneity <- function(y, g, label = "varianzas entre grupos") {
  ok <- tryCatch({
    lt <- car::leveneTest(y ~ g, center = median)
    list(F = lt[1, "F value"], df1 = lt[1, "Df"], df2 = lt[2, "Df"], p = lt[1, "Pr(>F)"])
  }, error = function(e) NULL)
  if (is.null(ok)) return(check_obj("homogeneidad", "Homogeneidad de varianzas", NA,
                                    "No fue posible aplicar la prueba de Levene."))
  met <- ok$p > ALPHA
  txt <- paste0("La prueba de Levene (centro = mediana) para ", label, " (F(", ok$df1, ", ", ok$df2, ") = ",
                fmt_num(ok$F, 3), ", p ", fmt_p(ok$p), ") ",
                if (met) "indica varianzas homogéneas: el supuesto se cumple." else "indica varianzas significativamente diferentes: el supuesto se viola.")
  check_obj("homogeneidad", "Homogeneidad de varianzas", met, txt, F = ok$F, df1 = ok$df1, df2 = ok$df2, p = ok$p)
}

# Frecuencias esperadas para tablas de contingencia (regla del libro: ≥ 5; N < 40 → Fisher)
check_expected_counts <- function(tab) {
  ct <- suppressWarnings(chisq.test(tab, correct = FALSE))
  mn <- min(ct$expected); N <- sum(tab)
  met <- mn >= 5 && N >= 40
  txt <- paste0("La frecuencia esperada mínima es ", fmt_num(mn), " (N = ", N, "). ",
                if (mn >= 5) "Se cumple el supuesto de frecuencias esperadas ≥ 5" else "Alguna celda tiene frecuencia esperada < 5",
                if (N < 40) "; además N < 40, por lo que se recomienda la prueba exacta de Fisher." else ".")
  check_obj("esperadas", "Frecuencias esperadas (≥ 5, N ≥ 40)", met, txt, min_expected = mn, N = N)
}

# Esfericidad (Mauchly) a partir de la salida de ezANOVA
check_sphericity <- function(ez) {
  m <- ez[["Mauchly's Test for Sphericity"]]
  if (is.null(m)) return(check_obj("esfericidad", "Esfericidad (Mauchly)", NA,
                                   "La prueba de Mauchly no aplica (dos niveles) o no pudo calcularse."))
  W <- as.numeric(m$W[1]); p <- as.numeric(m$p[1]); met <- p > ALPHA
  txt <- paste0("La prueba de Mauchly (W = ", fmt_num(W, 3), ", p ", fmt_p(p), ") ",
                if (met) "no indica violación de la esfericidad: se usa el estadístico F sin corrección."
                else "indica violación de la esfericidad: se reportan los valores p corregidos (Greenhouse-Geisser).")
  check_obj("esfericidad", "Esfericidad (Mauchly)", met, txt, W = W, p = p)
}

# Colinealidad (VIF < 5 sin problema; 5–10 moderada; > 10 grave)
check_vif <- function(model) {
  v <- tryCatch(suppressWarnings(car::vif(model)), error = function(e) NULL)
  if (is.null(v)) return(check_obj("vif", "Multicolinealidad (VIF)", NA, "VIF no aplicable (un solo predictor)."))
  vv <- if (is.matrix(v)) v[, 1] else v
  mx <- max(vv)
  met <- mx < 5
  txt <- paste0("VIF máximo = ", fmt_num(mx), ": ",
                if (mx < 5) "sin problemas de multicolinealidad." else if (mx <= 10) "multicolinealidad moderada — precaución." else "multicolinealidad grave — considerar eliminar o combinar predictores.")
  check_obj("vif", "Multicolinealidad (VIF)", met, txt, vif = vv, max = mx)
}

# Riesgos proporcionales (Schoenfeld): p > 0.05 global y por variable → se cumple
check_ph <- function(cox) {
  z <- tryCatch(survival::cox.zph(cox), error = function(e) NULL)
  if (is.null(z)) return(check_obj("ph", "Riesgos proporcionales (Schoenfeld)", NA, "No fue posible evaluar cox.zph()."))
  pg <- z$table["GLOBAL", "p"]; met <- pg > ALPHA
  txt <- paste0("La prueba de Schoenfeld (p global ", fmt_p(pg), ") ",
                if (met) "no muestra relación entre los residuos y el tiempo: el supuesto de riesgos proporcionales se cumple."
                else "indica relación con el tiempo: el supuesto se viola. Alternativa del libro: estratificar por la covariable problemática.")
  check_obj("ph", "Riesgos proporcionales (Schoenfeld)", met, txt, p = pg, zph = z)
}

# Resumen en prosa de todos los chequeos (para el reporte)
assumption_sentences <- function(checks) {
  if (length(checks) == 0) return("")
  paste(sapply(checks, `[[`, "text"), collapse = " ")
}
