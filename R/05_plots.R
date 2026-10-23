# =============================================================================
# 05_plots.R — Gráficos del libro con parámetros cosméticos (cos)
# Convención: cada builder devuelve un objeto ggplot, o una FUNCIÓN sin
# argumentos que dibuja (para gráficos base como Q-Q, TukeyHSD, ROC).
# =============================================================================
library(ggplot2)

cos_default <- function() list(
  fill      = ADO_FILL, colour = ADO_COL,
  base_size = 12, title = "", xlab = NULL, ylab = NULL,
  jitter = TRUE, show_center = TRUE, show_ci = TRUE, legend = FALSE,
  point_alpha = 0.55, point_size = 2.5
)
cos_merge <- function(cos) { d <- cos_default(); for (n in names(cos)) d[[n]] <- cos[[n]]; d }
.lab <- function(cos, what, default) if (is.null(cos[[what]]) || !nzchar(cos[[what]])) default else cos[[what]]
.finish <- function(p, cos, xlab, ylab) {
  p + labs(title = if (nzchar(cos$title)) cos$title else NULL, x = .lab(cos, "xlab", xlab), y = .lab(cos, "ylab", ylab)) +
    ado_theme(cos$base_size) + theme(legend.position = if (isTRUE(cos$legend)) "bottom" else "none")
}
.fillN <- function(cos, k) rep_len(cos$fill, k)
.colN  <- function(cos, k) rep_len(cos$colour, k)

# ---- Boxplot + jitter (patrón de todos los capítulos de medias) --------------
plot_box <- function(d, y, g, cos = list(), center = "mean") {
  cos <- cos_merge(cos)
  dd <- data.frame(y = num(d, y), g = fac(d, g)); dd <- dd[complete.cases(dd), ]
  p <- ggplot(dd, aes(g, y, fill = g)) +
    geom_boxplot(alpha = 0.6, outlier.shape = if (cos$jitter) NA else 16) +
    scale_fill_manual(values = .fillN(cos, nlevels(dd$g)))
  if (cos$jitter) p <- p + geom_jitter(width = 0.12, alpha = cos$point_alpha, size = cos$point_size, colour = "gray20")
  if (cos$show_center) p <- p + stat_summary(fun = if (center == "mean") mean else median, geom = "point", size = 3.5, colour = "black", shape = 18)
  .finish(p, cos, g, y)
}

# ---- Q-Q (car::qqPlot, estilo del libro) — función base -----------------------
plot_qq <- function(x, main = "Q-Q", cos = list()) {
  cos <- cos_merge(cos)
  x <- as.numeric(x); x <- x[!is.na(x)]
  function() {
    car::qqPlot(x, main = if (nzchar(cos$title)) cos$title else main,
                xlab = "Cuantiles teóricos (Normal)", ylab = "Cuantiles observados",
                pch = 16, col = cos$colour[1], col.lines = cos$colour[2], id = FALSE)
  }
}

# ---- Diferencia de medias con IC (pointrange, cap. t independiente) ------------
plot_effect <- function(est, lo, hi, label = "Diferencia de medias", cos = list()) {
  cos <- cos_merge(cos)
  dd <- data.frame(x = label, est = est, lo = lo, hi = hi)
  p <- ggplot(dd, aes(x, est)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = cos$colour[2]) +
    geom_pointrange(aes(ymin = lo, ymax = hi), colour = cos$colour[1], size = 0.9) + coord_flip()
  .finish(p, cos, "", "Δ (IC 95 %)")
}

# ---- Diseños apareados / medidas repetidas: líneas por sujeto + boxplot --------
plot_paired <- function(long, cos = list(), center = "mean") {
  cos <- cos_merge(cos)
  p <- ggplot(long, aes(Times, Response)) +
    geom_line(aes(group = Px), alpha = 0.3, colour = "gray50") +
    geom_boxplot(aes(fill = Times), width = 0.4, alpha = 0.6, outlier.shape = NA) +
    scale_fill_manual(values = .fillN(cos, nlevels(long$Times)))
  if (cos$show_center) p <- p +
    stat_summary(fun = if (center == "mean") mean else median, geom = "point", size = 3.5, colour = "black", aes(group = 1)) +
    stat_summary(fun = if (center == "mean") mean else median, geom = "line", linewidth = 1, colour = "black", aes(group = 1))
  .finish(p, cos, "Condición / tiempo", "Respuesta")
}

# ---- Interacción (ANOVA factorial): medias ± IC bootstrap ---------------------
plot_interaction <- function(d, y, A, B, cos = list()) {
  cos <- cos_merge(cos)
  dd <- data.frame(y = num(d, y), A = fac(d, A), B = fac(d, B)); dd <- dd[complete.cases(dd), ]
  p <- ggplot(dd, aes(A, y, colour = B, group = B)) +
    stat_summary(fun = mean, geom = "point", size = 3.5) +
    stat_summary(fun = mean, geom = "line", linewidth = 1)
  if (cos$show_ci) p <- p + stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.18)
  p <- p + scale_colour_manual(values = .colN(cos, nlevels(dd$B)), name = B)
  cos$legend <- TRUE
  .finish(p, cos, A, paste("Media de", y))
}

# ---- Dispersión + recta (correlación) ------------------------------------------
plot_scatter <- function(d, x, y, cos = list(), r_label = NULL, jitter = FALSE) {
  cos <- cos_merge(cos)
  dd <- data.frame(x = num(d, x), y = num(d, y)); dd <- dd[complete.cases(dd), ]
  p <- ggplot(dd, aes(x, y))
  p <- if (jitter) p + geom_jitter(width = 0.2, height = 0.2, alpha = cos$point_alpha, colour = cos$colour[1], size = cos$point_size)
       else p + geom_point(size = cos$point_size, colour = cos$colour[1], alpha = 0.85)
  p <- p + geom_smooth(method = "lm", se = cos$show_ci, alpha = 0.15, colour = cos$colour[2], formula = y ~ x)
  if (!is.null(r_label)) p <- p + annotate("label", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, label = r_label, size = 4)
  .finish(p, cos, x, y)
}

# ---- Barras apiladas de proporción (categóricos) -------------------------------
plot_bar_prop <- function(d, rowvar, colvar, cos = list()) {
  cos <- cos_merge(cos)
  dd <- data.frame(a = fac(d, rowvar), b = fac(d, colvar)); dd <- dd[complete.cases(dd), ]
  p <- ggplot(dd, aes(a, fill = b)) + geom_bar(position = "fill", alpha = 0.85) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = .colN(cos, nlevels(dd$b)), name = colvar)
  cos$legend <- TRUE
  .finish(p, cos, rowvar, "Proporción")
}

# ---- Mapa de calor de residuos estandarizados ajustados (ASR) -------------------
plot_heat_asr <- function(tab, asr, cos = list()) {
  cos <- cos_merge(cos)
  df <- as.data.frame(as.table(asr)); names(df) <- c("fila", "col", "ASR")
  ob <- as.data.frame(as.table(tab)); names(ob) <- c("fila", "col", "n")
  df <- merge(df, ob)
  df$lab <- paste0("n=", df$n, "\nASR=", round(df$ASR, 2))
  p <- ggplot(df, aes(col, fila, fill = ASR)) + geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = lab), size = 3.3, lineheight = 0.95) +
    scale_fill_gradient2(low = cos$colour[2], mid = "white", high = cos$colour[1], midpoint = 0)
  cos$legend <- TRUE
  .finish(p, cos, "", "")
}

# ---- Bland-Altman (función del libro) --------------------------------------------
plot_ba <- function(metodo, referencia, nombre = "Método", cos = list(), unidad = "") {
  cos <- cos_merge(cos)
  prom <- (metodo + referencia) / 2; dif <- metodo - referencia
  bias <- mean(dif, na.rm = TRUE); sd_d <- sd(dif, na.rm = TRUE)
  loa_u <- bias + 1.96 * sd_d; loa_l <- bias - 1.96 * sd_d; ref <- mean(prom, na.rm = TRUE)
  pct <- function(v) sprintf("%.1f %%", 100 * v / ref)
  dd <- data.frame(prom, dif)
  xr <- max(prom, na.rm = TRUE)
  p <- ggplot(dd, aes(prom, dif)) +
    geom_point(size = cos$point_size, colour = cos$colour[1], alpha = 0.8) +
    geom_hline(yintercept = bias, colour = cos$colour[2], linewidth = 1.1) +
    geom_hline(yintercept = c(loa_u, loa_l), colour = cos$colour[2], linetype = 2, linewidth = 0.8) +
    annotate("text", x = xr, y = bias,  label = paste0("Sesgo = ", round(bias, 2), " (", pct(bias), ")"), hjust = 1, vjust = -0.5, size = 3.4) +
    annotate("text", x = xr, y = loa_u, label = paste0("+1.96 SD = ", round(loa_u, 2), " (", pct(loa_u), ")"), hjust = 1, vjust = -0.5, size = 3.4) +
    annotate("text", x = xr, y = loa_l, label = paste0("-1.96 SD = ", round(loa_l, 2), " (", pct(loa_l), ")"), hjust = 1, vjust = 1.5, size = 3.4) +
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.15)))
  .finish(p, cos, paste0("Promedio de las dos mediciones", if (nzchar(unidad)) paste0(" (", unidad, ")") else ""),
          paste0("Diferencia (", nombre, " − referencia)"))
}

# ---- Gráfico de identidad (CCC de Lin) --------------------------------------------
plot_identity <- function(x, y, xlab = "Método", ylab = "Referencia", cos = list()) {
  cos <- cos_merge(cos)
  dd <- data.frame(x, y); dd <- dd[complete.cases(dd), ]
  p <- ggplot(dd, aes(x, y)) + geom_point(size = cos$point_size, colour = cos$colour[1], alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, colour = cos$colour[2], linetype = 2, linewidth = 1.2) +
    geom_smooth(method = "lm", se = cos$show_ci, alpha = 0.15, colour = cos$colour[1], fill = cos$fill[1], formula = y ~ x) +
    coord_equal()
  .finish(p, cos, xlab, ylab)
}

# ---- Kaplan-Meier (survminer::ggsurvplot) -----------------------------------------
plot_km <- function(fit, d, cos = list(), pval = FALSE, risk_table = FALSE, legend_labs = NULL) {
  cos <- cos_merge(cos)
  k <- max(1, length(fit$strata))
  args <- list(fit, data = d, conf.int = cos$show_ci, pval = pval, risk.table = risk_table,
               palette = .colN(cos, k), xlab = .lab(cos, "xlab", "Tiempo"),
               ylab = .lab(cos, "ylab", "Probabilidad de supervivencia S(t)"),
               ggtheme = ado_theme(cos$base_size), title = if (nzchar(cos$title)) cos$title else NULL)
  if (!is.null(legend_labs)) args$legend.labs <- legend_labs
  if (risk_table) { args$tables.height <- 0.25; args$tables.theme <- survminer::theme_cleantable() }
  do.call(survminer::ggsurvplot, args)
}

# ---- Curva ROC (pROC, base) --------------------------------------------------------
plot_roc <- function(roc_obj, cos = list()) {
  cos <- cos_merge(cos)
  function() {
    plot(roc_obj, col = cos$colour[1], lwd = 2, main = if (nzchar(cos$title)) cos$title else paste0("Curva ROC — AUC = ", round(as.numeric(pROC::auc(roc_obj)), 3)))
    abline(a = 0, b = 1, lty = 2, col = "gray50")
  }
}

# ---- Diagnósticos de regresión (base plot(mod, k)) ----------------------------------
plot_diag <- function(model, which = c(1, 2, 3, 5)) {
  function() { op <- par(mfrow = c(2, 2)); on.exit(par(op)); for (k in which) plot(model, k) }
}

# ---- IC de Tukey (multcomp::confint(glht)) ------------------------------------------
plot_tukey <- function(glht_obj, cos = list()) {
  cos <- cos_merge(cos)
  function() { op <- par(mar = c(4, 12, 3, 2)); on.exit(par(op))
    plot(confint(glht_obj), col = cos$colour[1], main = if (nzchar(cos$title)) cos$title else "Tukey HSD — IC 95 %") }
}

# ---- Medias ajustadas (effects) ----------------------------------------------------
plot_adjusted <- function(eff, cos = list()) {
  cos <- cos_merge(cos)
  function() plot(eff, main = if (nzchar(cos$title)) cos$title else "Medias ajustadas — ANCOVA",
                  lines = list(col = cos$colour[1]), symbols = list(col = cos$colour[1], pch = 16), confint = list(col = cos$colour[1]))
}

# Dibuja cualquier objeto devuelto por los builders
draw_plot <- function(p) { if (is.function(p)) p() else print(p) }
