# =============================================================================
# 03_assumptions.R - Motor de supuestos (cap. "Supuestos del análisis paramétrico")
# Reglas del libro:
#  - Shapiro-Wilk: p > 0.05 → sin evidencia de desviación de la normalidad
#    (no equivale a demostrarla); p < 0.05 → supuesto violado.
#  - n > 50: Shapiro tiende a ser significativo; la decisión se toma con
#    sesgo/curtosis estandarizados: |skew.2SE| o |kurt.2SE| > 1 → violado (Field et al., 2012).
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
  # Cap. 9 (9.2.5–9.2.6): con n > 50 Shapiro-Wilk tiende a ser significativo ante desviaciones pequeñas;
  # la decisión se apoya en sesgo y curtosis estandarizados (sesgo y curtosis divididos entre 2·EE: |skew.2SE| o |kurt.2SE| > 1 → desviación, p < 0.05).
  large <- n > 50 && !is.na(skew2) && !is.na(kurt2)
  sw_txt <- if (is.na(p)) paste0("No fue posible aplicar Shapiro-Wilk a ", label, " (n = ", n, ").")
            else paste0("Shapiro-Wilk sobre ", label, ": W = ", fmt_num(W, 3), ", p ", fmt_p(p), ".")
  if (large) {
    met <- abs(skew2) <= 1 && abs(kurt2) <= 1
    txt <- paste0(sw_txt, " Con n > 50 (n = ", n, ") la prueba de Shapiro-Wilk tiende a ser significativa ante desviaciones pequeñas, por lo que la decisión se apoya en el gráfico Q-Q ",
                  "y en el sesgo y la curtosis estandarizados (skew.2SE = ", fmt_num(skew2), "; kurt.2SE = ", fmt_num(kurt2), "): ",
                  if (met) "|valores| ≤ 1, sin evidencia de una desviación relevante de la normalidad."
                  else "|valor| > 1, desviación significativa de la normalidad: el supuesto se considera violado.")
  } else {
    met <- if (is.na(p)) NA else p > ALPHA
    txt <- if (is.na(p)) sw_txt
           else if (met) paste0(sw_txt, " No hay evidencia de una desviación de la normalidad (lo cual no equivale a demostrarla).")
           else paste0(sw_txt, " Indica una desviación significativa de la normalidad: el supuesto se considera violado.")
  }
  check_obj("normalidad", paste0("Normalidad (", label, ")"), met, txt,
            W = W, p = p, n = n, skew2SE = skew2, kurt2SE = kurt2, large_n = n > 50, rule = if (large) "skew_kurt" else "shapiro")
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

# Frecuencias esperadas para tablas de contingencia (regla: alguna esperada < 5 → Fisher)
check_expected_counts <- function(tab) {
  ct <- suppressWarnings(chisq.test(tab, correct = FALSE))
  mn <- min(ct$expected); N <- sum(tab)
  # Criterio (cap. 39): solo las frecuencias esperadas; alguna < 5 → Fisher
  met <- mn >= 5
  txt <- paste0("La frecuencia esperada mínima es ", fmt_num(mn), " (N = ", N, "). ",
                if (met) "Se cumple el supuesto de frecuencias esperadas ≥ 5: chi-cuadrado de Pearson."
                else "Alguna celda tiene frecuencia esperada < 5: se usa la prueba exacta de Fisher.")
  check_obj("esperadas", "Frecuencias esperadas (≥ 5)", met, txt, min_expected = mn, N = N)
}

# Esfericidad (Mauchly) a partir de la salida de ezANOVA
check_sphericity <- function(ez) {
  m <- ez[["Mauchly's Test for Sphericity"]]
  if (is.null(m)) return(check_obj("esfericidad", "Esfericidad (Mauchly)", NA,
                                   "La prueba de Mauchly no aplica (dos niveles) o no pudo calcularse."))
  W <- as.numeric(m$W[1]); p <- as.numeric(m$p[1]); met <- p > ALPHA
  txt <- paste0("La prueba de Mauchly (W = ", fmt_num(W, 3), ", p ", fmt_p(p), ") ",
                if (met) "no indica violación de la esfericidad: se usa el estadístico F sin corrección."
                else "indica violación de la esfericidad: se aplica la corrección de Greenhouse-Geisser (ε < 0.75) o Huynh-Feldt (ε ≥ 0.75).")
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
                if (mx < 5) "sin problemas de multicolinealidad." else if (mx <= 10) "multicolinealidad moderada: precaución." else "multicolinealidad grave: considerar eliminar o combinar predictores.")
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
