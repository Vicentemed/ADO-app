# =============================================================================
# app.R — CDO: asistente de análisis guiado (Ciencia de Datos en Odontología)
# Flujo: análisis → datos → preprocesar → variables (+ estrato opcional) →
#        descriptiva → supuestos/decisión → análisis → post-hoc → gráficos → reporte Word
# Si se elige una variable de estratificación, TODA la secuencia se ejecuta
# por separado dentro de cada nivel (subanálisis) y se muestra en pestañas.
# =============================================================================
suppressPackageStartupMessages({ library(shiny); library(bslib); library(DT); library(ggplot2) })
for (f in c("00_config.R", "01_registry.R", "02_data.R", "03_assumptions.R", "04_descriptives.R",
            "05_plots.R", "06_report.R", "07_ui_helpers.R", "08_preproc.R")) source(file.path("R", f), encoding = "UTF-8")
for (f in list.files(file.path("R", "analyses"), full.names = TRUE)) source(f, encoding = "UTF-8")

STEPS <- c("1. Datos", "2. Preprocesar", "3. Variables", "4. Descriptiva", "5. Supuestos", "6. Análisis", "7. Post-hoc", "8. Gráficos", "9. Reporte")

ui <- page_sidebar(
  title = tags$span(tags$b("CDO"), " · Asistente de análisis — ", tags$i("Ciencia de Datos en Odontología")),
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2E86C1"),
  tags$head(tags$link(rel = "stylesheet", href = "style.css")),
  sidebar = sidebar(width = 330, open = TRUE,
    tags$h6("Tipo de análisis"),
    selectInput("family", "Familia", choices = setNames(names(FAMILIES), FAMILIES), selectize = FALSE),
    uiOutput("analysis_ui"),
    uiOutput("analysis_desc"),
    tags$hr(),
    tags$h6("Fuente de datos"),
    radioButtons("src", NULL, choices = c("Ejemplo del libro" = "example", "Subir mi archivo" = "upload"), inline = TRUE),
    conditionalPanel("input.src == 'upload'", fileInput("file", "Archivo (.txt/.tsv tab, .csv, .xlsx)", accept = c(".txt", ".tsv", ".csv", ".xlsx", ".xls"))),
    conditionalPanel("input.src == 'example'", uiOutput("example_note")),
    uiOutput("data_info"),
    uiOutput("renamed_ui"),
    uiOutput("strata_info"),
    tags$hr(),
    div(class = "d-flex gap-2", actionButton("prev", "◀ Anterior", class = "btn-outline-secondary btn-sm"), actionButton("next", "Siguiente ▶", class = "btn-primary btn-sm"))
  ),
  navset_card_pill(id = "steps",
    nav_panel(STEPS[1], value = "s1",
      card(card_header("Vista de los datos"), DTOutput("data_preview"))),
    nav_panel(STEPS[2], value = "s2", preproc_ui()),
    nav_panel(STEPS[3], value = "s3",
      layout_columns(col_widths = c(5, 7),
        card(card_header("Selección de variables"), uiOutput("roles_ui"), uiOutput("roles_msg"),
             tags$hr(), tags$h6("Análisis estratificado (opcional)"),
             uiOutput("strata_ui"),
             tags$p(class = "small text-muted", "Si eliges una variable, el análisis completo se repite por separado dentro de cada uno de sus niveles (subanálisis), con su propia descriptiva, supuestos, decisión, post-hoc, gráficos y reporte.")),
        card(card_header("Qué pide este análisis"), uiOutput("roles_help")))),
    nav_panel(STEPS[4], value = "s4",
      card(card_header("Estadística descriptiva (table1)"), uiOutput("desc_ui"))),
    nav_panel(STEPS[5], value = "s5",
      card(card_header("Verificación de supuestos y decisión (regla del libro)"), uiOutput("path_ui"), uiOutput("checks_ui"))),
    nav_panel(STEPS[6], value = "s6",
      layout_columns(col_widths = c(3, 9),
        card(card_header("Opciones"), uiOutput("opts_ui"), actionButton("run", "Ejecutar análisis", class = "btn-primary")),
        card(card_header(uiOutput("run_title")), uiOutput("run_ui")))),
    nav_panel(STEPS[7], value = "s7",
      card(card_header("Contrastes / post-hoc"), uiOutput("posthoc_ui"))),
    nav_panel(STEPS[8], value = "s8",
      layout_columns(col_widths = c(3, 9),
        card(card_header("Gráfico y parámetros cosméticos"), uiOutput("plot_stratum_ui"), uiOutput("plot_select"), cosmetic_ui("cos")),
        card(card_header("Vista"), plotOutput("plot", height = "520px")))),
    nav_panel(STEPS[9], value = "s9",
      layout_columns(col_widths = c(5, 7),
        card(card_header("Título y contexto del reporte"),
             textInput("rep_title", "Título del análisis (encabezado del Word)", value = "", width = "100%",
                       placeholder = "p. ej. Viabilidad celular a 24 h según material y proporción"),
             tags$label(class = "control-label", "Contexto de los datos (obligatorio para el reporte)"),
             textAreaInput("contexto", NULL, rows = 8, placeholder = "Describa el objetivo, la población, el diseño, el origen de los datos y el significado de cada variable…", width = "100%"),
             checkboxGroupInput("rep_parts", "Incluir en el reporte", choices = c("Preprocesamiento" = "preproc", "Descriptiva" = "desc", "Supuestos" = "checks", "Tablas del análisis" = "tables", "Post-hoc" = "posthoc", "Gráficos" = "plots"), selected = c("preproc", "desc", "checks", "tables", "posthoc", "plots")),
             downloadButton("download", "Generar reporte Word (.docx)", class = "btn-success")),
        card(card_header("Cómo reportar (prellenado con sus resultados)"), tags$p(class = "text-muted small", "Los corchetes [ ] indican contexto que usted debe especificar."), uiOutput("report_ui"))))
  )
)

server <- function(input, output, session) {
  # ---- análisis --------------------------------------------------------------
  output$analysis_ui <- renderUI({ req(input$family)
    selectInput("analysis", "Análisis", choices = cdo_choices(input$family), selectize = FALSE, width = "100%") })
  spec <- reactive({ req(input$analysis); req(cdo_has(input$analysis)); cdo_get(input$analysis) })
  output$analysis_desc <- renderUI({ s <- spec(); tags$p(class = "small text-muted", s$description, tags$br(), tags$i("Capítulo: ", s$chapter)) })
  output$example_note <- renderUI({ s <- spec(); tags$p(class = "small", tags$b(s$example$file), " — ", s$example$note) })

  # ---- datos -------------------------------------------------------------------
  d_raw <- reactive({
    s <- spec()
    if (input$src == "example") {
      d <- read_book_data(file.path("data", s$example$file))
      if (is.function(s$example$prep)) d <- s$example$prep(d)
      d
    } else { req(input$file); tryCatch(load_any(input$file$datapath, input$file$name), error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }) }
  })
  # Aviso de columnas renombradas (nombres con símbolos, espacios o que empiezan por dígito)
  renamed <- reactive({ d <- d_raw(); if (is.null(d)) NULL else attr(d, "renamed") })
  observeEvent(renamed(), { rn <- renamed(); req(length(rn))
    showNotification(tags$div(tags$b("Se renombraron columnas con nombres no válidos en R:"), tags$br(),
                              paste(paste0("«", names(rn), "» → ", rn), collapse = "; ")), type = "warning", duration = 15) })
  output$renamed_ui <- renderUI({ rn <- renamed(); if (!length(rn)) return(NULL)
    tags$div(class = "alert alert-warning py-2 small mt-2", tags$b("Columnas renombradas: "),
             tags$ul(class = "mb-0", lapply(names(rn), function(o) tags$li("«", o, "» → ", tags$code(rn[[o]])))),
             "Se quitaron símbolos/espacios para que R pueda usarlas en las fórmulas.") })
  # ---- preprocesamiento (paso 2): d() son los datos que usan todos los pasos siguientes ----
  pp <- preproc_server(input, output, session, d_raw)
  d  <- pp$data
  output$data_info <- renderUI({ dd <- d(); if (is.null(dd)) return(NULL); n_ops <- length(pp$log())
    tags$p(class = "small", tags$b(nrow(dd)), " filas × ", tags$b(ncol(dd)), " columnas", if (n_ops) tags$span(class = "text-info", tags$br(), "(", n_ops, " op. de preprocesamiento)")) })
  output$data_preview <- renderDT({ dd <- d(); req(dd); datatable(dd, options = list(pageLength = 10, scrollX = TRUE), class = "compact stripe", rownames = FALSE) })

  # ---- variables + estrato --------------------------------------------------------
  output$roles_ui <- renderUI({ s <- spec(); dd <- d(); role_input_ui("var", s, dd, preset = if (input$src == "example") s$example$vars else list()) })
  v <- reactive({ read_roles(input, "var", spec()) })
  output$strata_ui <- renderUI({ dd <- d(); req(dd); used <- unlist(v())
    cand <- setdiff(names(dd)[sapply(dd, function(x) { k <- var_kind(x); k %in% c("factor", "binary", "numeric_small") })], used)
    selectInput("strata", "Variable de estratificación", choices = c("(ninguna — análisis global)" = "", cand), selected = isolate(input$strata) %||% "", selectize = FALSE, width = "100%") })
  strata <- reactive({ x <- input$strata; if (is.null(x) || !nzchar(x)) NULL else x })
  problems <- reactive({ dd <- d(); req(dd); p <- validate_roles(dd, spec(), v()); st <- strata()
    if (!is.null(st) && st %in% unlist(v())) p <- c(p, "La variable de estratificación no puede ser también una variable del análisis."); p })
  output$roles_msg <- renderUI({ p <- problems(); if (length(p)) tags$div(class = "alert alert-warning small", tags$ul(lapply(p, tags$li))) else tags$div(class = "alert alert-success small", "Selección válida.") })
  output$roles_help <- renderUI({ s <- spec(); tags$ul(lapply(names(s$roles), function(rn) { r <- s$roles[[rn]]; tags$li(tags$b(r$label), " — tipo: ", r$type, if (isTRUE(r$multiple)) " (varias)", if (!isTRUE(r$required)) " (opcional)") })) })
  ok <- reactive({ length(problems()) == 0 })

  # Lista de subconjuntos: uno por estrato (o uno solo con todos los datos)
  d_list <- reactive({ dd <- d(); req(dd, ok()); st <- strata()
    if (is.null(st)) return(list(`Todos los datos` = dd))
    g <- fac(dd, st); L <- split(dd, g); L <- L[sapply(L, nrow) > 0]
    setNames(L, paste0(st, " = ", names(L))) })
  is_strat <- reactive({ length(d_list()) > 1 })
  output$strata_info <- renderUI({ if (!isTRUE(tryCatch(is_strat(), error = function(e) FALSE))) return(NULL)
    tags$p(class = "small text-info", "Estratificado por ", tags$b(strata()), ": ", length(d_list()), " subanálisis") })

  # ---- descriptiva (por estrato) -------------------------------------------------------
  output$desc_ui <- renderUI({ L <- d_list(); s <- spec()
    strata_tabs(names(L), function(i) { t1 <- tryCatch(s$descriptives(L[[i]], v()), error = function(e) NULL)
      tagList(tags$p(class = "small text-muted", "n = ", nrow(L[[i]])), if (is.null(t1)) tags$em("Este análisis no genera tabla descriptiva.") else table1_html(t1)) }) })

  # ---- supuestos y decisión (por estrato) ------------------------------------------------
  assump_list <- reactive({ L <- d_list(); s <- spec()
    lapply(L, function(dd) tryCatch(s$assumptions(dd, v()), error = function(e) list(checks = list(check_obj("error", "Error", NA, conditionMessage(e))), decision = list(parametric = TRUE, all_met = NA, violated = character(0), reason = paste("No se pudieron evaluar los supuestos:", conditionMessage(e)))))) })
  output$path_ui <- renderUI({ s <- spec(); if (is.null(s$alternative)) return(NULL)
    radioButtons("path", "Ruta a seguir", inline = TRUE,
                 choices = setNames(c("auto", "primary", "alt"), c("Automática según los supuestos (regla del libro)", s$label, cdo_get(s$alternative)$label)), selected = "auto") })
  spec_for <- function(a) { s <- spec(); if (is.null(s$alternative)) return(s); p <- input$path %||% "auto"
    if (identical(p, "alt")) cdo_get(s$alternative) else if (identical(p, "primary")) s else if (isTRUE(a$decision$parametric)) s else cdo_get(s$alternative) }
  output$checks_ui <- renderUI({ A <- assump_list(); s <- spec()
    strata_tabs(names(A), function(i) { a <- A[[i]]; sr <- spec_for(a)
      tagList(tags$div(class = if (isTRUE(a$decision$parametric)) "alert alert-success py-2" else "alert alert-warning py-2", tags$b(a$decision$reason), tags$br(), "Prueba que se ejecutará: ", tags$b(sr$label)),
              checks_ui(a$checks)) }) })

  # ---- análisis (por estrato) -----------------------------------------------------------------
  spec_opts <- reactive({ A <- assump_list(); spec_for(A[[1]]) })   # opciones según la prueba del primer estrato
  output$opts_ui <- renderUI({ opts_input_ui("opt", spec_opts()) })
  output$run_title <- renderUI({ tags$span("Resultados: ", tags$b(spec()$label), if (is_strat()) tags$span(" · estratificado por ", tags$b(strata()))) })
  res_list <- eventReactive(input$run, { req(ok()); L <- d_list(); A <- assump_list()
    lapply(seq_along(L), function(i) { sr <- spec_for(A[[i]]); opts <- read_opts(input, "opt", sr)
      r <- tryCatch(sr$run(L[[i]], v(), opts), error = function(e) { showNotification(paste0(names(L)[i], ": ", conditionMessage(e)), type = "error", duration = 10); NULL })
      list(spec = sr, res = r, opts = opts, name = names(L)[i]) }) })
  output$run_ui <- renderUI({ R <- res_list(); if (is.null(R)) return(tags$em("Pulse «Ejecutar análisis»."))
    strata_tabs(sapply(R, `[[`, "name"), function(i) { x <- R[[i]]; if (is.null(x$res)) return(tags$em("No se pudo ejecutar en este estrato."))
      tagList(tags$p(class = "small", tags$b(x$spec$label), if (!is.null(x$res$notes)) tags$span(class = "text-muted", " — ", x$res$notes)), tables_ui(x$res$tables)) }) })

  # ---- post-hoc (por estrato) --------------------------------------------------------------
  ph_list <- reactive({ R <- res_list(); req(R); L <- d_list()
    lapply(seq_along(R), function(i) { x <- R[[i]]; if (is.null(x$res) || is.null(x$spec$posthoc)) return(NULL)
      tryCatch(x$spec$posthoc(x$res, L[[i]], v(), x$opts), error = function(e) { showNotification(paste("Post-hoc:", conditionMessage(e)), type = "warning"); NULL }) }) })
  output$posthoc_ui <- renderUI({ R <- res_list(); if (is.null(R)) return(tags$em("Ejecute primero el análisis.")); P <- ph_list()
    strata_tabs(sapply(R, `[[`, "name"), function(i) { p <- P[[i]]; if (is.null(p)) return(tags$p("Este análisis no requiere post-hoc (dos grupos o modelo global)."))
      tagList(if (!is.null(p$notes)) tags$p(class = "small text-muted", p$notes), tables_ui(p$tables)) }) })

  # ---- gráficos ----------------------------------------------------------------------
  output$plot_stratum_ui <- renderUI({ R <- res_list(); req(R); if (length(R) < 2) return(NULL)
    selectInput("plot_stratum", "Estrato", choices = setNames(seq_along(R), sapply(R, `[[`, "name")), selectize = FALSE) })
  cur_idx <- reactive({ R <- res_list(); req(R); i <- if (length(R) < 2) 1 else as.integer(input$plot_stratum %||% 1); min(max(i, 1), length(R)) })
  output$plot_select <- renderUI({ R <- res_list(); req(R); pl <- R[[cur_idx()]]$spec$plots; if (!length(pl)) return(tags$em("Sin gráficos."))
    selectInput("plot_key", "Gráfico", choices = setNames(names(pl), sapply(pl, `[[`, "label"))) })
  output$plot <- renderPlot({ R <- res_list(); req(R, input$plot_key); i <- cur_idx(); x <- R[[i]]; req(x$res); fn <- x$spec$plots[[input$plot_key]]$fn; req(fn)
    p <- tryCatch(fn(d_list()[[i]], v(), x$res, read_cos(input, "cos")), error = function(e) { plot.new(); text(0.5, 0.5, conditionMessage(e)); NULL }); if (!is.null(p)) draw_plot(p) }, res = 96)

  # ---- reporte -------------------------------------------------------------------------
  report_list <- reactive({ R <- res_list(); req(R); A <- assump_list(); L <- d_list()
    lapply(seq_along(R), function(i) { x <- R[[i]]; if (is.null(x$res)) return("No se pudo ejecutar el análisis en este estrato.")
      tryCatch(x$spec$report(x$res, L[[i]], v(), A[[i]]$checks), error = function(e) paste("No se pudo generar el texto:", conditionMessage(e))) }) })
  output$report_ui <- renderUI({ R <- res_list(); req(R); T <- report_list()
    strata_tabs(sapply(R, `[[`, "name"), function(i) tags$div(class = "cdo-report", markdown(T[[i]]))) })

  output$download <- downloadHandler(
    filename = function() { tt <- trimws(input$rep_title %||% "")
      base <- if (nzchar(tt)) gsub("[^A-Za-z0-9_.-]+", "_", iconv(tt, to = "ASCII//TRANSLIT")) else paste0("Reporte_", spec()$id, if (is_strat()) paste0("_por_", strata()) else "")
      paste0(base, "_", format(Sys.time(), "%Y%m%d_%H%M"), ".docx") },
    content = function(file) {
      R <- res_list(); A <- assump_list(); L <- d_list(); P <- ph_list(); T <- report_list(); cos <- read_cos(input, "cos"); parts <- input$rep_parts
      withProgress(message = "Generando reporte Word…", value = 0.2, {
        blocks <- lapply(seq_along(R), function(i) { x <- R[[i]]; if (is.null(x$res)) return(NULL)
          list(name = x$name, n = nrow(L[[i]]), analysis = x$spec$label,
               t1 = if ("desc" %in% parts) tryCatch(x$spec$descriptives(L[[i]], v()), error = function(e) NULL) else NULL,
               checks = if ("checks" %in% parts) A[[i]]$checks else list(), decision = A[[i]]$decision,
               tables = if ("tables" %in% parts) x$res$tables else list(), posthoc_tables = if ("posthoc" %in% parts && !is.null(P[[i]])) P[[i]]$tables else list(),
               plots = if ("plots" %in% parts) Filter(Negate(is.null), lapply(x$spec$plots, function(pl) tryCatch(pl$fn(L[[i]], v(), x$res, cos), error = function(e) NULL))) else list(),
               plot_labels = sapply(x$spec$plots, `[[`, "label"), report_md = T[[i]], notes = x$res$notes) })
        blocks <- Filter(Negate(is.null), blocks); s <- spec()
        tt <- trimws(input$rep_title %||% "")
        ctx <- list(title = if (nzchar(tt)) tt else paste("Reporte:", s$label, if (is_strat()) paste("— estratificado por", strata()) else ""),
                    subtitle = if (nzchar(tt)) paste(s$label, if (is_strat()) paste("— estratificado por", strata()) else "") else NULL,
                    analysis = s$label, family = FAMILIES[[s$family]], chapter = s$chapter,
                    contexto = if (nzchar(trimws(input$contexto))) input$contexto else contexto_placeholder(),
                    dataset = if (input$src == "example") s$example$file else input$file$name, vars = v(), strata = strata(),
                    preproc_log = if ("preproc" %in% parts) pp$log() else character(0),
                    n_raw = nrow(d_raw()), n_proc = nrow(d()), k_raw = ncol(d_raw()), k_proc = ncol(d()), blocks = blocks)
        incProgress(0.4); build_report(ctx, file); incProgress(0.4) })
    })

  # ---- navegación -------------------------------------------------------------------------
  step_ids <- paste0("s", seq_along(STEPS))
  observeEvent(input[["next"]], { i <- match(input$steps, step_ids); if (!is.na(i) && i < length(step_ids)) nav_select("steps", step_ids[i + 1]) })
  observeEvent(input[["prev"]], { i <- match(input$steps, step_ids); if (!is.na(i) && i > 1) nav_select("steps", step_ids[i - 1]) })
}

shinyApp(ui, server)
