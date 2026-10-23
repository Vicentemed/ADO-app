# =============================================================================
# 07_ui_helpers.R — Constructores de UI genéricos sobre el registro
# =============================================================================
library(shiny)

# ---- Selectores de variables por rol ----------------------------------------
role_input_ui <- function(prefix, spec, d, preset = list()) {
  if (is.null(d)) return(tags$em("Cargue un set de datos primero."))
  tagList(lapply(names(spec$roles), function(rn) {
    r <- spec$roles[[rn]]; ch <- vars_for_role(d, r); sel <- preset[[rn]]
    sel <- sel[sel %in% ch]
    if (!isTRUE(r$required) && !isTRUE(r$multiple)) ch <- c("(ninguna)" = "", ch)
    selectInput(paste0(prefix, "_", rn), label = tags$span(r$label, if (isTRUE(r$required)) tags$span(" *", style = "color:#E74C3C")),
                choices = ch, selected = if (length(sel)) sel else if (isTRUE(r$multiple)) NULL else ch[1], multiple = isTRUE(r$multiple), width = "100%")
  }))
}
read_roles <- function(input, prefix, spec) {
  v <- lapply(names(spec$roles), function(rn) { x <- input[[paste0(prefix, "_", rn)]]; if (is.null(x)) character(0) else x[nzchar(x)] })
  setNames(v, names(spec$roles))
}

# ---- Opciones específicas del análisis ------------------------------------------
opts_input_ui <- function(prefix, spec) {
  if (!length(spec$opts_ui)) return(NULL)
  tagList(lapply(names(spec$opts_ui), function(on) {
    o <- spec$opts_ui[[on]]; id <- paste0(prefix, "_opt_", on)
    switch(o$type,
      checkbox = checkboxInput(id, o$label, value = isTRUE(o$value)),
      select   = selectInput(id, o$label, choices = o$choices, selected = o$value, width = "100%"),
      text     = textInput(id, o$label, value = o$value, width = "100%"),
      numeric  = numericInput(id, o$label, value = o$value, width = "100%"))
  }))
}
read_opts <- function(input, prefix, spec) {
  if (!length(spec$opts_ui)) return(list())
  setNames(lapply(names(spec$opts_ui), function(on) input[[paste0(prefix, "_opt_", on)]]), names(spec$opts_ui))
}

# ---- Parámetros cosméticos de los gráficos ----------------------------------------
PALETTES <- list(
  "Libro (azul/rojo/verde)" = list(colour = CDO_COL, fill = CDO_FILL),
  "Escala de grises"        = list(colour = c("#222222", "#666666", "#999999", "#BBBBBB", "#DDDDDD", "#444444"), fill = c("#D9D9D9", "#BFBFBF", "#A6A6A6", "#8C8C8C", "#737373", "#595959")),
  "Viridis"                 = list(colour = c("#440154", "#31688E", "#35B779", "#FDE725", "#21918C", "#90D743"), fill = c("#B5A3CC", "#9DBBD6", "#A9E1BD", "#F8F2A3", "#A6D8D3", "#D4EDB5")),
  "Pastel"                  = list(colour = c("#8E9BFF", "#FF9E9E", "#8ED8A0", "#FFD08E", "#C99EFF", "#9EDCE8"), fill = c("#D3D8FF", "#FFD3D3", "#D3F0DA", "#FFEBD3", "#E7D3FF", "#D3F1F6"))
)
cosmetic_ui <- function(prefix = "cos") tagList(
  selectInput(paste0(prefix, "_palette"), "Paleta de colores", choices = names(PALETTES), selected = names(PALETTES)[1]),
  sliderInput(paste0(prefix, "_base"), "Tamaño de fuente", min = 8, max = 20, value = 12, step = 1),
  textInput(paste0(prefix, "_title"), "Título", ""),
  textInput(paste0(prefix, "_xlab"), "Etiqueta eje X (vacío = automática)", ""),
  textInput(paste0(prefix, "_ylab"), "Etiqueta eje Y (vacío = automática)", ""),
  checkboxInput(paste0(prefix, "_jitter"), "Mostrar puntos individuales (jitter)", TRUE),
  checkboxInput(paste0(prefix, "_center"), "Marcar media/mediana", TRUE),
  checkboxInput(paste0(prefix, "_ci"), "Mostrar IC / banda de confianza", TRUE),
  checkboxInput(paste0(prefix, "_legend"), "Mostrar leyenda", FALSE),
  sliderInput(paste0(prefix, "_alpha"), "Transparencia de puntos", min = 0.1, max = 1, value = 0.55, step = 0.05),
  sliderInput(paste0(prefix, "_psize"), "Tamaño de puntos", min = 1, max = 5, value = 2.5, step = 0.5)
)
read_cos <- function(input, prefix = "cos") {
  pal <- PALETTES[[input[[paste0(prefix, "_palette")]] %||% names(PALETTES)[1]]]
  list(colour = pal$colour, fill = pal$fill, base_size = input[[paste0(prefix, "_base")]] %||% 12,
       title = input[[paste0(prefix, "_title")]] %||% "", xlab = input[[paste0(prefix, "_xlab")]] %||% "", ylab = input[[paste0(prefix, "_ylab")]] %||% "",
       jitter = isTRUE(input[[paste0(prefix, "_jitter")]]), show_center = isTRUE(input[[paste0(prefix, "_center")]]), show_ci = isTRUE(input[[paste0(prefix, "_ci")]]),
       legend = isTRUE(input[[paste0(prefix, "_legend")]]), point_alpha = input[[paste0(prefix, "_alpha")]] %||% 0.55, point_size = input[[paste0(prefix, "_psize")]] %||% 2.5)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Presentación de chequeos de supuestos ------------------------------------------
checks_ui <- function(checks) {
  if (!length(checks)) return(tags$p("Este análisis no requiere verificación de supuestos."))
  tags$ul(class = "cdo-checks", lapply(checks, function(ch) {
    ico <- if (is.na(ch$met)) tags$span("ℹ", class = "cdo-ico cdo-info") else if (isTRUE(ch$met)) tags$span("✔", class = "cdo-ico cdo-ok") else tags$span("✘", class = "cdo-ico cdo-bad")
    tags$li(ico, tags$b(ch$label), tags$br(), tags$span(ch$text))
  }))
}

# ---- Pestañas por estrato ------------------------------------------------------------------
# names: etiquetas de los estratos; content_fn(i) devuelve la UI del estrato i.
# Con un solo estrato devuelve el contenido directamente (sin pestañas).
strata_tabs <- function(names, content_fn) {
  if (length(names) <= 1) return(content_fn(1))
  do.call(bslib::navset_pill, c(list(id = NULL),
    lapply(seq_along(names), function(i) bslib::nav_panel(names[i], tags$div(class = "pt-2", content_fn(i))))))
}

# ---- Tablas ---------------------------------------------------------------------------
df_table <- function(df, digits = 3) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(NULL)
  num <- sapply(df, is.numeric); df[num] <- lapply(df[num], function(x) signif(x, digits + 1))
  DT::datatable(df, rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE), class = "compact stripe")
}
tables_ui <- function(tabs, ns_id) {
  if (!length(tabs)) return(tags$em("Sin tablas."))
  tagList(lapply(seq_along(tabs), function(i) tagList(tags$h6(names(tabs)[i]), DT::renderDT(df_table(tabs[[i]]), server = FALSE))))
}
