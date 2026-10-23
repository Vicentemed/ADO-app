# =============================================================================
# 08_preproc.R — Preprocesamiento de datos (paso 2 del flujo)
#   · Evaluar / eliminar valores faltantes (NA)
#   · Remover variables
#   · Subset de filas (filtro por variable o expresión R)
#   · Transformación ancho → largo (wide → long) y largo → ancho (long → wide)
#   · Cambiar tipo de variable (factor / numérica)
# Estilo R base, como en el libro (cap. "Primeros pasos con un set de datos").
# =============================================================================

# ---- Utilidades ----------------------------------------------------------------
na_summary <- function(d) {
  if (is.null(d) || !ncol(d)) return(data.frame())
  n_na <- sapply(d, function(x) sum(is.na(x)))
  data.frame(Variable = names(d), Tipo = sapply(d, function(x) class(x)[1]),
             `n NA` = as.integer(n_na), `% NA` = round(100 * n_na / nrow(d), 1),
             `Valores únicos` = sapply(d, function(x) length(unique(x[!is.na(x)]))),
             check.names = FALSE, row.names = NULL)
}

# ancho → largo (R base). Conserva las columnas id; si no hay id crea Px.
pp_wide_to_long <- function(d, measures, ids = character(0), key = "Times", value = "Response") {
  measures <- measures[measures %in% names(d)]; ids <- ids[ids %in% names(d)]
  if (length(measures) < 2) stop("Seleccione al menos dos columnas de medición.")
  n <- nrow(d); k <- length(measures)
  base <- if (length(ids)) d[rep(seq_len(n), times = k), ids, drop = FALSE] else data.frame(Px = factor(rep(seq_len(n), times = k)))
  long <- cbind(base, setNames(data.frame(factor(rep(measures, each = n), levels = measures)), key),
                setNames(data.frame(unlist(lapply(measures, function(m) d[[m]]), use.names = FALSE)), value))
  rownames(long) <- NULL
  long
}

# largo → ancho (R base, reshape). Nombres resultantes: <valor>_<nivel> como en el libro.
pp_long_to_wide <- function(d, id, key, value) {
  stopifnot(all(c(id, key, value) %in% names(d)))
  dd <- d[c(id, key, value)]
  w <- reshape(dd, direction = "wide", idvar = id, timevar = key, v.names = value)
  names(w) <- sub(paste0("^", value, "\\."), paste0(value, "_"), names(w))
  rownames(w) <- NULL
  attr(w, "reshapeWide") <- NULL
  w
}

# ---- UI --------------------------------------------------------------------------
preproc_ui <- function() {
  layout_columns(col_widths = c(5, 7),
    card(card_header("Operaciones"),
      tags$p(class = "small text-muted", "Las operaciones se aplican en orden sobre los datos actuales. Los pasos siguientes (variables, análisis, reporte) usan el resultado."),
      accordion(id = "pp_acc", open = FALSE, multiple = FALSE,
        accordion_panel("1 · Valores faltantes (NA)", uiOutput("pp_na_ui")),
        accordion_panel("2 · Remover variables", uiOutput("pp_drop_ui")),
        accordion_panel("3 · Subset de filas (filtrar)", uiOutput("pp_subset_ui")),
        accordion_panel("4 · Ancho → largo (wide → long)", uiOutput("pp_w2l_ui")),
        accordion_panel("5 · Largo → ancho (long → wide)", uiOutput("pp_l2w_ui")),
        accordion_panel("6 · Tipo de variable", uiOutput("pp_type_ui"))),
      tags$hr(),
      div(class = "d-flex gap-2 flex-wrap",
          actionButton("pp_undo", "↶ Deshacer última", class = "btn-outline-secondary btn-sm"),
          actionButton("pp_reset", "Restablecer datos originales", class = "btn-outline-danger btn-sm"))),
    card(card_header("Estado de los datos"),
      uiOutput("pp_status"),
      tags$h6("Valores faltantes por variable"), DTOutput("pp_na_table"),
      tags$h6("Historial de operaciones"), uiOutput("pp_log"),
      tags$h6("Vista previa (datos actuales)"), DTOutput("pp_preview"))
  )
}

# ---- Servidor --------------------------------------------------------------------
# Devuelve list(data = reactive con los datos procesados, log = reactive con el historial)
preproc_server <- function(input, output, session, d_raw) {
  wd    <- reactiveVal(NULL)      # datos de trabajo
  hist  <- reactiveVal(list())    # pila de estados previos (para deshacer)
  lg    <- reactiveVal(character(0))

  observeEvent(d_raw(), { wd(d_raw()); hist(list()); lg(character(0)) }, ignoreNULL = FALSE)

  apply_op <- function(fn, msg) {
    d <- wd(); if (is.null(d)) return(invisible())
    res <- tryCatch(fn(d), error = function(e) { showNotification(paste("No se aplicó:", conditionMessage(e)), type = "error", duration = 8); NULL })
    if (is.null(res)) return(invisible())
    if (!nrow(res)) { showNotification("La operación dejaría 0 filas; no se aplicó.", type = "warning"); return(invisible()) }
    hist(c(hist(), list(d))); wd(res); lg(c(lg(), msg))
    showNotification(msg, type = "message", duration = 4)
  }
  observeEvent(input$pp_undo, { h <- hist(); if (!length(h)) return(); wd(h[[length(h)]]); hist(h[-length(h)]); l <- lg(); lg(l[-length(l)]) })
  observeEvent(input$pp_reset, { wd(d_raw()); hist(list()); lg(character(0)); showNotification("Datos restablecidos.", type = "message") })

  vars <- reactive({ d <- wd(); if (is.null(d)) character(0) else names(d) })
  nvars <- reactive({ d <- wd(); if (is.null(d)) character(0) else names(d)[sapply(d, is.numeric)] })

  # -- 1. NA ------------------------------------------------------------------------
  output$pp_na_ui <- renderUI({ req(wd()); tagList(
    radioButtons("pp_na_scope", "Eliminar filas con NA en:", choices = c("cualquier variable" = "all", "las variables seleccionadas" = "sel")),
    conditionalPanel("input.pp_na_scope == 'sel'", selectInput("pp_na_vars", NULL, choices = vars(), multiple = TRUE, width = "100%")),
    actionButton("pp_na_rows", "Eliminar filas con NA", class = "btn-primary btn-sm"),
    tags$hr(),
    numericInput("pp_na_pct", "Eliminar variables con más de este % de NA:", value = 50, min = 0, max = 100, step = 5, width = "60%"),
    actionButton("pp_na_cols", "Eliminar variables", class = "btn-primary btn-sm")) })
  observeEvent(input$pp_na_rows, {
    sc <- input$pp_na_scope; vs <- if (identical(sc, "sel")) input$pp_na_vars else NULL
    apply_op(function(d) { keep <- if (is.null(vs) || !length(vs)) complete.cases(d) else complete.cases(d[vs]); d[keep, , drop = FALSE] },
             paste0("Eliminadas filas con NA", if (!is.null(vs) && length(vs)) paste0(" en: ", paste(vs, collapse = ", ")) else " en cualquier variable",
                    " (", sum(!(if (is.null(vs) || !length(vs)) complete.cases(wd()) else complete.cases(wd()[vs]))), " filas)")) })
  observeEvent(input$pp_na_cols, {
    th <- input$pp_na_pct
    apply_op(function(d) { pct <- 100 * sapply(d, function(x) mean(is.na(x))); drop <- names(d)[pct > th]
      if (!length(drop)) stop("Ninguna variable supera ese porcentaje de NA."); d[setdiff(names(d), drop)] },
             paste0("Eliminadas variables con > ", th, " % NA: ", paste(names(wd())[100 * sapply(wd(), function(x) mean(is.na(x))) > th], collapse = ", "))) })

  # -- 2. Remover variables --------------------------------------------------------------
  output$pp_drop_ui <- renderUI({ req(wd()); tagList(
    selectInput("pp_drop_vars", "Variables a remover", choices = vars(), multiple = TRUE, width = "100%"),
    actionButton("pp_drop", "Remover", class = "btn-primary btn-sm")) })
  observeEvent(input$pp_drop, { vs <- input$pp_drop_vars; req(length(vs))
    apply_op(function(d) { if (length(vs) >= ncol(d)) stop("No se pueden remover todas las variables."); d[setdiff(names(d), vs)] }, paste0("Removidas: ", paste(vs, collapse = ", "))) })

  # -- 3. Subset ---------------------------------------------------------------------------
  output$pp_subset_ui <- renderUI({ req(wd()); tagList(
    selectInput("pp_sub_var", "Filtrar por variable", choices = vars(), width = "100%"),
    uiOutput("pp_sub_cond"),
    actionButton("pp_sub_apply", "Aplicar filtro", class = "btn-primary btn-sm"),
    tags$hr(),
    textInput("pp_sub_expr", "O bien, expresión R (como en subset()):", placeholder = 'p. ej.  Age > 30 & Sex == "F"', width = "100%"),
    actionButton("pp_sub_expr_apply", "Aplicar expresión", class = "btn-outline-primary btn-sm")) })
  output$pp_sub_cond <- renderUI({ d <- wd(); v <- input$pp_sub_var; req(d, v, v %in% names(d)); x <- d[[v]]
    if (is.numeric(x) && length(unique(x[!is.na(x)])) > 8) {
      r <- range(x, na.rm = TRUE); sliderInput("pp_sub_range", "Conservar valores en el rango", min = floor(r[1]), max = ceiling(r[2]), value = c(floor(r[1]), ceiling(r[2])), width = "100%")
    } else {
      lv <- sort(unique(as.character(x[!is.na(x)])))
      checkboxGroupInput("pp_sub_levels", "Conservar categorías", choices = lv, selected = lv)
    } })
  observeEvent(input$pp_sub_apply, { d <- wd(); v <- input$pp_sub_var; req(d, v); x <- d[[v]]
    if (is.numeric(x) && length(unique(x[!is.na(x)])) > 8) { rg <- input$pp_sub_range; req(rg)
      apply_op(function(d) d[!is.na(d[[v]]) & d[[v]] >= rg[1] & d[[v]] <= rg[2], , drop = FALSE], paste0("Filtro: ", v, " en [", rg[1], ", ", rg[2], "]"))
    } else { lv <- input$pp_sub_levels; req(length(lv))
      apply_op(function(d) { r <- d[!is.na(d[[v]]) & as.character(d[[v]]) %in% lv, , drop = FALSE]; if (is.factor(r[[v]])) r[[v]] <- droplevels(r[[v]]); r }, paste0("Filtro: ", v, " ∈ {", paste(lv, collapse = ", "), "}")) } })
  observeEvent(input$pp_sub_expr_apply, { ex <- trimws(input$pp_sub_expr); req(nzchar(ex))
    apply_op(function(d) { e <- parse(text = ex)[[1]]; keep <- eval(e, d, parent.frame()); if (!is.logical(keep)) stop("La expresión debe producir TRUE/FALSE."); d[keep & !is.na(keep), , drop = FALSE] }, paste0("Filtro (expresión): ", ex)) })

  # -- 4. Ancho → largo ----------------------------------------------------------------------
  output$pp_w2l_ui <- renderUI({ req(wd()); tagList(
    selectInput("pp_w2l_meas", "Columnas de medición (se apilarán)", choices = nvars(), multiple = TRUE, width = "100%"),
    selectInput("pp_w2l_ids", "Columnas identificadoras a conservar (opcional; si no hay, se crea Px)", choices = vars(), multiple = TRUE, width = "100%"),
    layout_columns(textInput("pp_w2l_key", "Nombre de la columna de tiempo/condición", "Times"), textInput("pp_w2l_val", "Nombre de la columna de respuesta", "Response")),
    actionButton("pp_w2l", "Convertir a formato largo", class = "btn-primary btn-sm")) })
  observeEvent(input$pp_w2l, { m <- input$pp_w2l_meas; ids <- input$pp_w2l_ids; k <- trimws(input$pp_w2l_key); vl <- trimws(input$pp_w2l_val); req(length(m) >= 2)
    apply_op(function(d) pp_wide_to_long(d, m, ids, if (nzchar(k)) k else "Times", if (nzchar(vl)) vl else "Response"),
             paste0("Ancho → largo: ", paste(m, collapse = ", "), " → ", if (nzchar(k)) k else "Times", "/", if (nzchar(vl)) vl else "Response")) })

  # -- 5. Largo → ancho --------------------------------------------------------------------------
  output$pp_l2w_ui <- renderUI({ req(wd()); tagList(
    selectInput("pp_l2w_id", "Identificador del sujeto", choices = vars(), width = "100%"),
    selectInput("pp_l2w_key", "Variable de tiempo/condición", choices = vars(), width = "100%"),
    selectInput("pp_l2w_val", "Variable de respuesta (valores)", choices = nvars(), width = "100%"),
    actionButton("pp_l2w", "Convertir a formato ancho", class = "btn-primary btn-sm")) })
  observeEvent(input$pp_l2w, { id <- input$pp_l2w_id; k <- input$pp_l2w_key; vl <- input$pp_l2w_val; req(id, k, vl)
    if (length(unique(c(id, k, vl))) < 3) { showNotification("Las tres variables deben ser distintas.", type = "warning"); return() }
    apply_op(function(d) pp_long_to_wide(d, id, k, vl), paste0("Largo → ancho: ", vl, " por ", k, " (id = ", id, ")")) })

  # -- 6. Tipo de variable -------------------------------------------------------------------------
  output$pp_type_ui <- renderUI({ req(wd()); tagList(
    selectInput("pp_type_var", "Variable", choices = vars(), width = "100%"),
    radioButtons("pp_type_to", "Convertir a", choices = c("factor (categórica)" = "factor", "numérica" = "numeric"), inline = TRUE),
    actionButton("pp_type", "Convertir", class = "btn-primary btn-sm")) })
  observeEvent(input$pp_type, { v <- input$pp_type_var; to <- input$pp_type_to; req(v)
    apply_op(function(d) { d[[v]] <- if (to == "factor") factor(d[[v]]) else { x <- suppressWarnings(as.numeric(as.character(d[[v]]))); if (all(is.na(x))) stop("No se puede convertir a numérica."); x }; d },
             paste0("Tipo: ", v, " → ", to)) })

  # -- Estado / salidas ------------------------------------------------------------------------------
  output$pp_status <- renderUI({ d <- wd(); r <- d_raw(); if (is.null(d)) return(tags$em("Cargue un set de datos."))
    mod <- length(lg()) > 0
    tags$div(class = if (mod) "alert alert-info py-2 small" else "alert alert-light py-2 small",
             tags$b(nrow(d)), " filas × ", tags$b(ncol(d)), " columnas",
             if (mod) tags$span(" (original: ", nrow(r), " × ", ncol(r), ") — ", tags$b(length(lg())), " operación(es) aplicada(s)") else " — sin modificaciones.",
             tags$br(), "Total de NA: ", tags$b(sum(is.na(d)))) })
  output$pp_na_table <- renderDT({ d <- wd(); req(d); datatable(na_summary(d), rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollY = "220px", scrollX = TRUE), class = "compact stripe") })
  output$pp_log <- renderUI({ l <- lg(); if (!length(l)) return(tags$p(class = "text-muted small", "Ninguna operación aplicada.")); tags$ol(class = "small", lapply(l, tags$li)) })
  output$pp_preview <- renderDT({ d <- wd(); req(d); datatable(d, options = list(pageLength = 6, scrollX = TRUE, lengthMenu = c(6, 15, 50)), class = "compact stripe", rownames = FALSE) })

  list(data = reactive(wd()), log = reactive(lg()))
}
