# =============================================================================
# analyses/correlacion.R - Pearson (alt: Spearman), Kendall, parcial, PPS
# Reglas del libro (cap. 10 Correlaciones):
#  - Pearson requiere datos continuos, normalidad de ambas variables y relación
#    lineal (10.2.1). La normalidad se evalúa con el gráfico Q-Q de car::qqPlot():
#    puntos dentro de la banda de confianza al 95 % → se aproximan a la normal
#    (Fig. 10.2); puntos fuera de la banda → no normal → Spearman (Fig. 10.4).
#  - Spearman: alternativa no paramétrica (variables no normales u ordinales);
#    no da IC directamente (sería necesario bootstrap) (10.3).
#  - Kendall: datos ordinales con pocos valores posibles / muchos empates; el
#    libro la presenta junto a Spearman como comparación (10.4).
#  - Parcial: ggm::pcor() sobre var(); se muestran la bivariada y la relación de
#    Z con X e Y (10.5).
#  - PPS: ppsr::score() / visualize_pps() (10.6).
# =============================================================================

# ---- Banda Q-Q (idéntica a car::qqPlot: línea por cuartiles, banda al 95 %) ----
# El libro decide visualmente. Para automatizar esa lectura la app cuenta los
# puntos fuera de la banda; con una banda puntual al 95 % se espera que ~5 % de
# los puntos caigan fuera por azar, de modo que «≤ 5 % fuera» se toma como
# «los puntos permanecen dentro de la banda». Con este criterio se reproducen
# las dos decisiones del libro (IMC/Tqx → normal; edad dental/cronológica → no).
.qq_band_check <- function(x, label) {
  x <- sort(as.numeric(x)[!is.na(x)]); n <- length(x)
  if (n < 3) return(check_obj("normalidad", paste0("Normalidad - Q-Q (", label, ")"), NA,
                              paste0("No es posible construir el gráfico Q-Q de ", label, " (n = ", n, ").")))
  P <- ppoints(n); z <- qnorm(P)
  Qx <- quantile(x, c(0.25, 0.75), names = FALSE); Qz <- qnorm(c(0.25, 0.75))
  b <- (Qx[2] - Qx[1]) / (Qz[2] - Qz[1]); a <- Qx[1] - b * Qz[1]
  SE <- (b / dnorm(z)) * sqrt(P * (1 - P) / n); fit <- a + b * z
  out <- sum(x > fit + qnorm(0.975) * SE | x < fit - qnorm(0.975) * SE)
  pct <- 100 * out / n; met <- if (b <= 0) FALSE else pct <= 5
  txt <- paste0("Gráfico Q-Q de ", label, " (car::qqPlot, banda de confianza al 95 %): ", out, " de ", n, " puntos (",
                fmt_num(pct, 1), " %) quedan fuera de la banda. ",
                if (met) "Los puntos permanecen dentro de la banda: la variable se aproxima a la distribución normal."
                else "Hay puntos fuera de la banda de confianza: la variable se aparta de la distribución normal.",
                " (Criterio operativo de la app para la lectura visual del libro: ≤ 5 % de puntos fuera de la banda; confírmelo en el gráfico Q-Q.)")
  check_obj("normalidad", paste0("Normalidad - Q-Q (", label, ")"), met, txt, out = out, n = n, pct_out = pct)
}

# Magnitud según los criterios de Cohen del libro (0.10 pequeño, 0.30 mediano, ≥ 0.50 grande);
# por debajo de 0.10 el libro escribe «sumamente pequeña» (r = 0.028).
.mag_cor <- function(r) { r <- abs(r); if (is.na(r)) "…" else if (r < 0.1) "sumamente pequeña" else if (r < 0.3) "pequeña" else if (r < 0.5) "mediana" else "grande" }
.dir_txt <- function(est, x, y) if (est > 0) paste0("a mayor ", x, ", mayor ", y) else paste0("a mayor ", x, ", menor ", y)

.cor_dd <- function(d, v) { dd <- data.frame(x = num(d, v$x), y = num(d, v$y)); dd[complete.cases(dd), ] }

.cor_run <- function(d, v, method) {
  dd <- .cor_dd(d, v)
  ct <- suppressWarnings(cor.test(dd$x, dd$y, alternative = "two.sided", method = method))
  est <- unname(ct$estimate); n <- nrow(dd)
  if (method == "pearson") {
    ci <- ct$conf.int
    list(main = ct,
         stats = list(est = est, p = ct$p.value, lo = ci[1], hi = ci[2], n = n, sym = "r", r2 = est^2,
                      stat = unname(ct$statistic), df = unname(ct$parameter)),
         tables = list("Resultado de cor.test(method = \"pearson\")" = data.frame(
           t = unname(ct$statistic), gl = unname(ct$parameter), p = ct$p.value, `IC 95% inf` = ci[1], `IC 95% sup` = ci[2],
           r = est, `R²` = est^2, n = n, check.names = FALSE)),
         notes = "R² = r² indica la proporción de varianza de una variable explicada por la otra.")
  } else {  # spearman
    list(main = ct,
         stats = list(est = est, p = ct$p.value, lo = NA, hi = NA, n = n, sym = "ρ", r2 = NA, stat = unname(ct$statistic), df = NA),
         tables = list("Resultado de cor.test(method = \"spearman\")" = data.frame(S = unname(ct$statistic), p = ct$p.value, rho = est, n = n, check.names = FALSE)),
         notes = "En Spearman no se obtienen intervalos de confianza directamente; para calcularlos sería necesario usar métodos bootstrap.")
  }
}

# Frase de supuestos a partir de los chequeos Q-Q (los dos primeros: X e Y)
.qq_sentence <- function(checks, v) {
  qx <- checks[[1]]; qy <- checks[[2]]
  both <- isTRUE(qx$met) && isTRUE(qy$met)
  paste0("La normalidad de ", CTX(v$x), " y de ", CTX(v$y), " se evaluó mediante gráficos Q-Q con banda de confianza al 95 % (",
         v$x, ": ", qx$out %||% "…", " de ", qx$n %||% "…", " puntos fuera de la banda; ", v$y, ": ", qy$out %||% "…", " de ", qy$n %||% "…", "); ",
         if (both) "ambas variables se aproximan a la distribución normal"
         else "al menos una variable se aparta de la normalidad (puntos fuera de la banda de confianza)")
}

.cor_report <- function(res, d, v, checks, method) {
  s <- res$stats; sig <- s$p < ALPHA
  pos <- if (s$est > 0) "positiva" else "negativa"
  both <- length(checks) >= 2 && isTRUE(checks[[1]]$met) && isTRUE(checks[[2]]$met)
  intro <- if (length(checks) >= 2 && !is.null(checks[[1]]$out)) paste0(.qq_sentence(checks, v), "; ") else ""
  if (method == "pearson") {
    paste0(intro, "por ello se aplicó la correlación de Pearson. ",
           if (sig) paste0("La salida muestra una relación ", .mag_cor(s$est), ", ", pos, " y significativa (r = ", fmt_num(s$est, 2), ", p ", fmt_p(s$p),
                           ", n = ", s$n, ") entre ", CTX(v$x), " y ", CTX(v$y), " (", .dir_txt(s$est, v$x, v$y), "), con un intervalo de confianza al 95 % (",
                           fmt_num(s$lo, 3), " a ", fmt_num(s$hi, 3), ") que ", if (s$lo > 0 || s$hi < 0) "no incluye" else "incluye", " el cero.")
           else paste0("La salida indica una correlación ", .mag_cor(s$est), " y no significativa (r = ", fmt_num(s$est, 3), ", p ", fmt_p(s$p),
                       ", n = ", s$n, "; IC 95 %: ", fmt_num(s$lo, 3), " a ", fmt_num(s$hi, 3), ") entre ", CTX(v$x), " y ", CTX(v$y),
                       ". Dado que p > 0.05, no se puede rechazar la hipótesis nula de que no existe relación entre ambas variables."),
           " El coeficiente de determinación es R² = ", fmt_num(s$r2, 3), ": el ", fmt_num(100 * s$r2, 1), " % de la variación de ", CTX(v$y),
           " es explicada por ", CTX(v$x), " [Cohen, 1988; Field et al., 2012].")
  } else {
    paste0(intro, if (both) "se aplicó la correlación de Spearman, alternativa no paramétrica de Pearson también indicada cuando las variables son ordinales. "
                  else "por ello se aplicó la correlación de Spearman, alternativa no paramétrica de Pearson. ",
           if (sig) paste0("La salida indica una relación ", pos, " ", .mag_cor(s$est), " y significativa (ρ = ", fmt_num(s$est, 3), ", p ", fmt_p(s$p),
                           ", n = ", s$n, ") entre ", CTX(v$x), " y ", CTX(v$y), ": ", .dir_txt(s$est, v$x, v$y), ".")
           else paste0("La salida indica una relación ", .mag_cor(s$est), " y no significativa (ρ = ", fmt_num(s$est, 3), ", p ", fmt_p(s$p),
                       ", n = ", s$n, ") entre ", CTX(v$x), " y ", CTX(v$y), ". Dado que p > 0.05, no se puede rechazar la hipótesis nula de que no existe relación entre ambas variables."),
           " En Spearman no se obtienen intervalos de confianza directamente; para calcularlos sería necesario usar métodos bootstrap [Field et al., 2012].")
  }
}

.cor_roles <- list(x = role("Variable X (numérica)", "numeric"), y = role("Variable Y (numérica)", "numeric"))
# Supuestos de Pearson (10.2.1): continuos, normalidad (Q-Q) y relación lineal
.cor_assump <- function(d, v) list(
  .qq_band_check(num(d, v$x), v$x), .qq_band_check(num(d, v$y), v$y),
  check_obj("info", "Datos continuos y relación lineal", NA,
            "Pearson requiere además datos continuos (escala de intervalo o razón) en ambas variables y una relación lineal entre ellas: verifíquelo en el diagrama de dispersión."))

.scatter_plot <- function(d, v, res, cos) {
  plot_scatter(d, v$x, v$y, cos, r_label = paste0(res$stats$sym, " = ", fmt_num(res$stats$est, 2)), jitter = isTRUE(cos_merge(cos)$jitter))
}
.cor_plots <- list(
  scatter = list(label = "Dispersión (jitter) + recta de regresión", fn = .scatter_plot),
  qqx = list(label = "Q-Q de X", fn = function(d, v, res, cos) plot_qq(num(d, v$x), paste("Q-Q:", v$x), cos)),
  qqy = list(label = "Q-Q de Y", fn = function(d, v, res, cos) plot_qq(num(d, v$y), paste("Q-Q:", v$y), cos)))

ado_register(list(
  id = "pearson", family = "correlacion", order = 1, label = "Correlación de Pearson",
  description = "Relación lineal entre dos variables numéricas continuas y normales (Q-Q). Alternativa: Spearman.",
  chapter = "Correlaciones", roles = .cor_roles,
  example = list(file = "DXKTw.txt", vars = list(x = "Imc", y = "Tqx"), note = "IMC y tiempo quirúrgico (min)."),
  alternative = "spearman",
  descriptives = function(d, v) desc_table1(d, c(v$x, v$y)),
  assumptions = function(d, v) {
    ch <- .cor_assump(d, v); dec <- decide_path(ch, TRUE)
    dec$reason <- if (dec$all_met) "Los gráficos Q-Q muestran que ambas variables se aproximan a la distribución normal (puntos dentro de la banda de confianza al 95 %): se aplica la correlación de Pearson."
                  else "En el gráfico Q-Q, al menos una variable se aparta de la normalidad (puntos fuera de la banda de confianza al 95 %): según el libro, se aplica la correlación de Spearman."
    list(checks = ch, decision = dec) },
  run = function(d, v, opts) .cor_run(d, v, "pearson"), posthoc = NULL, plots = .cor_plots,
  report = function(res, d, v, checks) .cor_report(res, d, v, checks, "pearson")
))
ado_register(list(
  id = "spearman", family = "correlacion", order = 2, label = "Correlación de Spearman (ρ)",
  description = "Alternativa no paramétrica de Pearson: una o ambas variables no normales, o variables ordinales.",
  chapter = "Correlaciones", roles = .cor_roles,
  example = list(file = "datcron.txt", vars = list(x = "Dental_Age", y = "Cron_Age"), note = "Edad dental y cronológica (meses)."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, c(v$x, v$y)),
  assumptions = function(d, v) { ch <- .cor_assump(d, v)[1:2]
    list(checks = ch, decision = list(parametric = FALSE, all_met = FALSE, violated = character(0),
                                      reason = "Spearman: alternativa no paramétrica de Pearson cuando una o ambas variables no cumplen la normalidad (Q-Q) o son ordinales.")) },
  run = function(d, v, opts) .cor_run(d, v, "spearman"), posthoc = NULL, plots = .cor_plots,
  report = function(res, d, v, checks) .cor_report(res, d, v, checks, "spearman")
))

# ---- Kendall: τ junto a ρ de Spearman (comparativa), como en el libro ----------
.ties_check <- function(x, label) {
  x <- as.numeric(x); x <- x[!is.na(x)]; n <- length(x); k <- length(unique(x))
  check_obj("info", paste0("Valores distintos y empates (", label, ")"), NA,
            paste0(label, " toma ", k, " valores distintos en ", n, " observaciones (", n - k, " observaciones repiten un valor ya presente, ",
                   fmt_num(100 * (n - k) / max(n, 1), 1), " %)."), k = k, n = n)
}
ado_register(list(
  id = "kendall", family = "correlacion", order = 3, label = "Tau de Kendall (τ)",
  description = "Recomendada con rangos de muy pocos valores posibles (ordinales con pocas categorías, muchos empates). Se muestra junto a Spearman.",
  chapter = "Correlaciones", roles = .cor_roles,
  example = list(file = "datcron.txt", vars = list(x = "Cron_Age", y = "Bone_Maturation"), note = "Edad cronológica y maduración ósea (ordinal 1-5)."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, c(v$x, v$y)),
  assumptions = function(d, v) list(checks = list(.ties_check(num(d, v$x), v$x), .ties_check(num(d, v$y), v$y)),
                                    decision = list(parametric = FALSE, all_met = FALSE, violated = character(0),
                                                    reason = "Según el libro, Kendall se recomienda cuando los datos son rangos con muy pocos valores posibles (muchos empates), p. ej. ordinales con pocas categorías; si hay pocos empates, use Spearman.")),
  run = function(d, v, opts) {
    dd <- .cor_dd(d, v)
    kt <- suppressWarnings(cor.test(dd$x, dd$y, alternative = "two.sided", method = "kendall"))
    sp <- suppressWarnings(cor.test(dd$x, dd$y, alternative = "two.sided", method = "spearman"))
    list(main = kt,
         stats = list(tau = unname(kt$estimate), p = kt$p.value, rho = unname(sp$estimate), p_rho = sp$p.value, n = nrow(dd),
                      k_x = length(unique(dd$x)), k_y = length(unique(dd$y))),
         tables = list("Kendall y Spearman (comparativa)" = data.frame(
           Coeficiente = c("Kendall τ", "Spearman ρ"), Estimación = c(unname(kt$estimate), unname(sp$estimate)),
           p = c(kt$p.value, sp$p.value), n = nrow(dd), check.names = FALSE)),
         notes = "Como en el libro, la rho de Spearman se calcula como comparativa de la tau de Kendall.")
  },
  posthoc = NULL,
  plots = list(scatter = list(label = "Dispersión (jitter) con τ y ρ", fn = function(d, v, res, cos) {
    cos <- cos_merge(cos); dd <- .cor_dd(d, v); s <- res$stats
    p <- ggplot(dd, aes(x, y)) +
      geom_jitter(height = 0.15, alpha = 0.3, colour = cos$colour[1], size = 1.5) +
      geom_smooth(method = "lm", se = cos$show_ci, colour = cos$colour[2], alpha = 0.15, formula = y ~ x) +
      annotate("text", x = min(dd$x) + 0.05 * diff(range(dd$x)), y = max(dd$y) - 0.1 * diff(range(dd$y)),
               label = paste0("Kendall τ = ", round(s$tau, 2), "\nSpearman ρ = ", round(s$rho, 2)), hjust = 0, size = 3.5, colour = "#333333")
    if (s$k_y <= 10 && all(dd$y == round(dd$y))) p <- p + scale_y_continuous(breaks = sort(unique(dd$y)))
    .finish(p, cos, v$x, v$y) })),
  report = function(res, d, v, checks) {
    s <- res$stats; sig <- s$p < ALPHA
    few <- if (s$k_y <= s$k_x) list(var = v$y, k = s$k_y) else list(var = v$x, k = s$k_x)
    paste0("Dado que ", CTX(few$var), " toma solo ", few$k, " valores distintos y, por tanto, presenta numerosos empates, se aplicó la tau de Kendall, más apropiada que la rho de Spearman en esta situación. ",
           "La relación entre ", CTX(v$x), " y ", CTX(v$y), " fue ", if (s$tau > 0) "positiva" else "negativa", " y ", if (sig) "significativa" else "no significativa",
           " (τ = ", fmt_num(s$tau, 3), ", p ", fmt_p(s$p), ", n = ", s$n, "); a modo comparativo, la rho de Spearman fue ρ = ", fmt_num(s$rho, 3), " (p ", fmt_p(s$p_rho), "). ",
           if (abs(s$tau) < abs(s$rho)) "La tau de Kendall es más conservadora que la rho de Spearman, lo cual es esperable: Kendall es más robusta ante empates y generalmente produce estimaciones más conservadoras" else "Kendall es más robusta ante empates",
           " [Field et al., 2012].")
  }
))

ado_register(list(
  id = "parcial", family = "correlacion", order = 4, label = "Correlación parcial (controlando Z)",
  description = "Correlación entre X e Y controlando el efecto de una tercera variable Z, posible variable de confusión (ggm::pcor).",
  chapter = "Correlaciones",
  roles = list(x = role("Variable X", "numeric"), y = role("Variable Y", "numeric"), z = role("Variable de control Z", "numeric")),
  example = list(file = "DXKTw.txt", vars = list(x = "Tqx", y = "Imc", z = "Ed"), note = "Tiempo quirúrgico e IMC controlando por edad."),
  alternative = NULL,
  descriptives = function(d, v) desc_table1(d, c(v$x, v$y, v$z)),
  assumptions = function(d, v) {
    ch <- list(check_obj("info", "Base de la correlación parcial", NA,
                         "pcor() (ggm) calcula la correlación parcial a partir de la matriz de covarianzas var(); se apoya en correlaciones de Pearson (véanse sus supuestos)."))
    list(checks = ch, decision = list(parametric = TRUE, all_met = NA, violated = character(0),
                                      reason = "Se calcula la correlación bivariada sin ajuste, la correlación parcial controlando Z y la relación de Z con X y con Y.")) },
  run = function(d, v, opts) {
    dd <- data.frame(x = num(d, v$x), y = num(d, v$y), z = num(d, v$z)); dd <- dd[complete.cases(dd), ]
    S <- var(dd); rp <- ggm::pcor(c("x", "y", "z"), S)
    biv <- cor.test(dd$x, dd$y, alternative = "two.sided", method = "pearson")
    c1 <- cor.test(dd$x, dd$z, method = "pearson"); c2 <- cor.test(dd$y, dd$z, method = "pearson")
    list(main = rp,
         stats = list(r_partial = rp, r_biv = unname(biv$estimate), p_biv = biv$p.value, r_xz = unname(c1$estimate), p_xz = c1$p.value,
                      r_yz = unname(c2$estimate), p_yz = c2$p.value, n = nrow(dd)),
         tables = list("Correlaciones" = data.frame(
           Relación = c(paste(v$x, "~", v$y, "(bivariada, sin ajuste)"), paste0(v$x, " ~ ", v$y, " | ", v$z, " (parcial)"), paste(v$x, "~", v$z), paste(v$y, "~", v$z)),
           r = c(unname(biv$estimate), rp, unname(c1$estimate), unname(c2$estimate)),
           p = c(biv$p.value, NA, c1$p.value, c2$p.value), check.names = FALSE)),
         notes = "Como en el libro, pcor() se reporta sin valor p ni IC; los valores p corresponden a cor.test() de Pearson.")
  },
  posthoc = NULL,
  plots = list(),
  report = function(res, d, v, checks) {
    s <- res$stats
    paste0("La correlación bivariada entre ", CTX(v$x), " y ", CTX(v$y), " fue r = ", fmt_num(s$r_biv, 3), " (p ", fmt_p(s$p_biv), "). Al controlar por ", CTX(v$z),
           ", la correlación parcial fue r = ", fmt_num(s$r_partial, 3), " (n = ", s$n, "); por su parte, ", CTX(v$z), " se asocia con ", CTX(v$x),
           " (r = ", fmt_num(s$r_xz, 2), ", p ", fmt_p(s$p_xz), ") y con ", CTX(v$y), " (r = ", fmt_num(s$r_yz, 2), ", p ", fmt_p(s$p_yz), "). ",
           CTX(paste0("Interprete la comparación entre la correlación bivariada y la parcial: si la asociación entre ", v$x, " y ", v$y,
                      " se mantiene una vez ajustado el efecto de ", v$z, ", esta no induce una asociación espuria; si cambia, ", v$z, " actúa como variable de confusión.")),
           " Siempre se debe considerar la posible influencia de variables de confusión al calcular correlaciones [Austin, 2011; Field et al., 2012].")
  }
))

ado_register(list(
  id = "pps", family = "correlacion", order = 5, label = "Predictive Power Score (PPS)",
  description = "Poder predictivo (0 a 1; asimétrico, no lineal, admite variables categóricas) de varias variables sobre una respuesta (ppsr).",
  chapter = "Correlaciones",
  roles = list(target = role("Variable objetivo (respuesta)", "any"), predictors = role("Predictores candidatos", "any", multiple = TRUE, min = 1)),
  example = list(file = "DXKTw.txt", vars = list(target = "Tqx", predictors = c("Classcx", "Imc", "Ed", "Sexo", "Grupo")), note = "¿Qué variables predicen el tiempo quirúrgico?"),
  alternative = NULL,
  descriptives = function(d, v) NULL,
  assumptions = function(d, v) list(checks = list(check_obj("info", "PPS", NA, "El PPS se calcula con árboles de clasificación y regresión (CART) comparados con un modelo base (mediana o clase más frecuente) y se estima mediante validación cruzada, por lo que su valor puede variar ligeramente entre ejecuciones.")),
                                    decision = list(parametric = TRUE, all_met = TRUE, violated = character(0), reason = "Se calcula el PPS de cada predictor sobre la respuesta (ppsr::score).")),
  run = function(d, v, opts) {
    # Como en el libro: se pasan las columnas tal cual; score() usa los casos completos de cada par x-y
    dd <- d[c(v$target, v$predictors)]
    sc <- do.call(rbind, lapply(v$predictors, function(x) {
      data.frame(x = x, y = v$target, pps = ppsr::score(dd, x = x, y = v$target, verbose = FALSE)$pps, pps_inverso = ppsr::score(dd, x = v$target, y = x, verbose = FALSE)$pps) }))
    sc <- sc[order(-sc$pps), ]; rownames(sc) <- NULL
    tab <- setNames(sc, c("x (predictor)", "y (respuesta)", "PPS x → y", "PPS y → x"))
    list(main = sc, stats = list(top = sc$x[1], top_pps = sc$pps[1]),
         tables = list("PPS de cada predictor sobre la respuesta (y su inverso)" = tab),
         notes = "PPS 0 = ausencia de poder predictivo; 1 = predicción perfecta. Es asimétrico (x → y ≠ y → x) y no indica la dirección del efecto.")
  },
  posthoc = NULL,
  plots = list(pps = list(label = "PPS por predictor (visualize_pps)", fn = function(d, v, res, cos) {
    cos <- cos_merge(cos)
    p <- ppsr::visualize_pps(d[c(v$target, v$predictors)], y = v$target, verbose = FALSE) + theme_bw(base_size = cos$base_size)
    if (nzchar(cos$title)) p <- p + labs(title = cos$title)
    p })),
  report = function(res, d, v, checks) {
    s <- res$stats; sc <- res$main
    paste0("Se calculó el Predictive Power Score (PPS) de ", length(v$predictors), " variables para predecir ", CTX(v$target), ". El mayor poder predictivo correspondió a ",
           CTX(s$top), " (PPS ", s$top, " → ", v$target, " = ", fmt_num(s$top_pps), "); ", paste(paste0(sc$x, " = ", fmt_num(sc$pps)), collapse = ", "), ". ",
           "La naturaleza asimétrica del PPS se aprecia al invertir los papeles: ", paste(paste0(v$target, " → ", sc$x, " = ", fmt_num(sc$pps_inverso)), collapse = ", "), ". ",
           "El PPS no sustituye a la correlación, sino que la complementa: informa sobre la capacidad predictiva (incluso no lineal) de una variable sobre otra, sin indicar la dirección del efecto. ",
           "Dado que se estima mediante validación cruzada, su valor puede variar ligeramente entre ejecuciones; para reportar un intervalo de confianza al 95 % puede emplearse remuestreo bootstrap [Laken, 2021].")
  }
))
