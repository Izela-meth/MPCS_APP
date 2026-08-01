# ==============================================================================
# MPCS_APP — Aplicacion Shiny del Modelo Predictivo de Cambio Conductual por Sistemas
# ==============================================================================
# Repositorio: https://github.com/Izela-meth/MPCS_APP
# ==============================================================================

# --- Cargar librerias ---
library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(readxl)
library(haven)
library(ggplot2)
library(igraph)
library(tidyr)
library(RColorBrewer)
library(patchwork)

# --- Cargar funciones modulares ---
source("functions/mpcs_functions.R", local = TRUE)

# ==============================================================================
# GENERAR DATOS DE DEMOSTRACION INICIALES (SOLO ENDES)
# ==============================================================================
if (!file.exists("data/demo_data.csv")) {
  dir.create("data", showWarnings = FALSE)
  datos_demo <- generar_demo_data(n = 1000, seed = 123)
  write.csv(datos_demo, "data/demo_data.csv", row.names = FALSE)
}

# ==============================================================================
# UI — INTERFAZ DE USUARIO
# ==============================================================================

ui <- page_navbar(
  title = tags$div(
    tags$img(src = "https://img.icons8.com/color/48/000000/artificial-intelligence.png", 
             height = "30px", style = "margin-right: 10px;"),
    "MPCS Calculator"
  ),
  theme = bs_theme(bootswatch = "flatly", version = 5),
  
  # ============================================================================
  # Pestaña 1: Carga de Datos
  # ============================================================================
  nav_panel(
    "1. Data Upload",
    fluidRow(
      column(
        width = 4,
        wellPanel(
          h4("Upload file"),
          fileInput("file", "Select a file",
                    accept = c(".csv", ".xlsx", ".xls", ".dta"),
                    buttonLabel = "Browse",
                    placeholder = "No file selected"),
          tags$small("Supported formats: CSV, Excel (.xlsx, .xls), Stata (.dta)"),
          hr(),
          h4("Or use demo data"),
          p("Load simulated health data (ENDES) to test the application."),
          actionButton("load_demo", "Load ENDES Data (Health)", 
                       class = "btn-primary w-100",
                       icon = icon("heart"))
        )
      ),
      column(
        width = 8,
        wellPanel(
          h5("Preview (first 10 rows)"),
          DTOutput("data_preview"),
          hr(),
          h5("Basic statistics"),
          verbatimTextOutput("data_stats")
        )
      )
    )
  ),
  
  # ============================================================================
  # Pestaña 2: Configuracion
  # ============================================================================
  nav_panel(
    "2. Configuration",
    fluidRow(
      column(
        width = 6,
        wellPanel(
          h4("Graph variables"),
          helpText("Select at least 5 numeric variables to build the behavioral graph."),
          uiOutput("graph_vars_ui"),
          hr(),
          h4("Grouping variable (optional)"),
          helpText("If selected, the analysis will be performed separately for each group."),
          uiOutput("group_var_ui")
        )
      ),
      column(
        width = 6,
        wellPanel(
          h4("Markov configuration"),
          helpText("Select the column containing the behavioral states."),
          uiOutput("markov_var_ui"),
          hr(),
          h4("Adjustable parameters"),
          sliderInput("threshold", "Correlation threshold",
                      min = 0.05, max = 0.30, value = 0.10, step = 0.01,
                      post = tags$span("  (|r| > value)")),
          
          # Campo para tasas de avance de Markov
          hr(),
          h4("Markov transition rates (optional)"),
          helpText("Enter the forward transition probabilities between consecutive states."),
          helpText("Example for 5 states (4 transitions): 0.20, 0.48, 0.89, 0.62")
textInput("tasas_avance_input", "Forward transition rates (comma separated):",
          placeholder = "e.g., 0.20, 0.48, 0.89, 0.62")
          helpText("If left empty, the generic heuristic matrix will be used."),
          
          # Modulo de Bootstrap para selección de umbral
          hr(),
          h4("Bootstrap threshold selection (optional)"),
          helpText("Evaluates the stability of the optimal node across different correlation thresholds."),
          helpText("Note: 20 replicates take ~30-60 seconds on Render."),
          fluidRow(
            column(6,
              numericInput("boot_n_sim", "Number of bootstrap replicates:",
                           value = 20,
                           min = 10,
                           max = 50,
                           step = 5)
            ),
            column(6,
              actionButton("run_bootstrap", "Run Bootstrap",
                           class = "btn-info w-100",
                           icon = icon("random"))
            )
          ),
          conditionalPanel(
            condition = "input.run_bootstrap > 0",
            div(style = "padding: 10px; background-color: #f8f9fa; border-radius: 5px; margin-top: 10px;",
                h5("Bootstrap Results"),
                uiOutput("bootstrap_results_ui"),
                plotOutput("bootstrap_plot", height = "250px")
            )
          ),
          div(
            class = "alert alert-info mt-2",
            style = "font-size: 0.85em;",
            icon("circle-info"),
            HTML(
              "<b>Tip:</b> Bootstrap is optional. If it takes too long, 
              just use the default threshold of 0.10 (validated in the article)."
            )
          ),
          
          hr(),
          h5("MPCS Index Weights"),
          fluidRow(
            column(4, numericInput("w1", "Graph (w1)", 
                                   value = 0.35, min = 0, max = 1, step = 0.05)),
            column(4, numericInput("w2", "Markov (w2)", 
                                   value = 0.40, min = 0, max = 1, step = 0.05)),
            column(4, numericInput("w3", "Games (w3)", 
                                   value = 0.25, min = 0, max = 1, step = 0.05))
          ),
          tags$small("Weights must sum to 1. Current: ", 
                     textOutput("suma_ponderadores", inline = TRUE)),
          hr(),
          sliderInput("R_factor", "Resource factor (R)",
                      min = 0, max = 1, value = 0.65, step = 0.05,
                      post = tags$span("  (higher = more resources available)")),
          hr(),
          h4("Contextual indicator (alpha)"),
          helpText("alpha represents the favorability of the context for behavior adoption."),
          fluidRow(
            column(6,
              radioButtons("alpha_mode", "Method for alpha:",
                           choices = c("Select variables" = "vars",
                                      "Enter manual value" = "manual",
                                      "Automatic" = "auto"),
                           selected = "auto")
            ),
            column(6,
              conditionalPanel(
                condition = "input.alpha_mode == 'vars'",
                uiOutput("alpha_vars_ui")
              ),
              conditionalPanel(
                condition = "input.alpha_mode == 'manual'",
                numericInput("alpha_manual", "Alpha value (0-1):",
                             value = 0.60, min = 0, max = 1, step = 0.01)
              ),
              conditionalPanel(
                condition = "input.alpha_mode == 'auto'",
                helpText("The app will use the average of access/resource related variables.")
              )
            )
          ),
          div(
            class = "alert alert-info mt-2",
            style = "font-size: 0.85em;",
            icon("circle-info"),
            HTML(
              "<b>Methodological note.</b><br><br>
              <b>I_Markov:</b> Calculated as the accumulated probability in advanced states 
              at a fixed horizon of 10 periods.<br><br>
              <b>alpha (context):</b> If you select variables, they are normalized to 0-1 and averaged. 
              If you enter manually, that value is used. In automatic mode, it looks for common 
              variables like 'Acceso_salud' or uses the first numeric variable."
            )
          )
        )
      )
    ),
    fluidRow(
      column(
        width = 12,
        wellPanel(
          actionButton("run_mpcs", "Run MPCS", 
                       class = "btn-success btn-lg w-100",
                       icon = icon("play")),
          uiOutput("validation_msg")
        )
      )
    )
  ),
  
  # ============================================================================
  # Pestaña 3: Resultados
  # ============================================================================
  nav_panel(
    "3. Results",
    fluidRow(
      column(
        width = 12,
        wellPanel(
          h4("Results table"),
          DTOutput("results_table")
        )
      )
    ),
    fluidRow(
      column(
        width = 6,
        wellPanel(
          h5("Behavioral graph"),
          plotOutput("plot_graph", height = "500px")
        )
      ),
      column(
        width = 6,
        wellPanel(
          h5("State distribution by group"),
          plotOutput("plot_states", height = "500px")
        )
      )
    ),
    fluidRow(
      column(
        width = 6,
        wellPanel(
          h5("Markov transitions"),
          plotOutput("plot_arbol_markov", height = "450px")
        )
      ),
      column(
        width = 6,
        wellPanel(
          h5("Game theory dynamics"),
          plotOutput("plot_juego_evolutivo", height = "450px")
        )
      )
    ),
    fluidRow(
      column(
        width = 6,
        wellPanel(
          h5("Markov trajectories"),
          plotOutput("plot_markov", height = "400px")
        )
      ),
      column(
        width = 6,
        wellPanel(
          h5("I_MPCS ranking"),
          plotOutput("plot_ranking", height = "400px")
        )
      )
    ),
    # ==========================================================================
    # Analisis de Sensibilidad
    # ==========================================================================
    fluidRow(
      column(
        width = 12,
        wellPanel(
          h4("Sensitivity Analysis"),
          p("Evaluates the robustness of the recommendation to changes in weights w1, w2, w3."),
          hr(),
          fluidRow(
            column(
              width = 3,
              numericInput("sens_n_sim", 
                           "Number of simulations:", 
                           value = 10000, 
                           min = 100, 
                           max = 50000,
                           step = 100),
              br(),
              actionButton("run_sensitivity", 
                           "Run analysis",
                           icon = icon("chart-line"),
                           class = "btn-primary",
                           style = "width: 100%;")
            ),
            column(
              width = 9,
              conditionalPanel(
                condition = "input.run_sensitivity > 0",
                div(style = "padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
                    h5("MPCS Index Distribution"),
                    plotOutput("sens_plot", height = "300px"),
                    verbatimTextOutput("sens_summary")
                )
              )
            )
          )
        )
      )
    ),
    # ==========================================================================
    # Interpretacion
    # ==========================================================================
    fluidRow(
      column(
        width = 12,
        wellPanel(
          h4("Automatic interpretation"),
          uiOutput("interpretation_text")
        )
      )
    )
  ),
  
  # ============================================================================
  # Pestaña 4: Reporte y Descarga
  # ============================================================================
  nav_panel(
    "4. Report & Download",
    fluidRow(
      column(
        width = 6,
        wellPanel(
          h4("Download data"),
          p("Download the results table as CSV."),
          downloadButton("download_csv", "Download CSV", 
                         class = "btn-primary w-100",
                         icon = icon("file-csv"))
        )
      ),
      column(
        width = 6,
        wellPanel(
          h4("Generate report"),
          p("Download a complete report in HTML format."),
          downloadButton("download_report", "Download HTML Report", 
                         class = "btn-danger w-100",
                         icon = icon("file-code"))
        )
      )
    )
  )
)

# ==============================================================================
# SERVER — LOGICA DE LA APLICACION
# ==============================================================================

server <- function(input, output, session) {
  
  rv <- reactiveValues(
    data = NULL,
    results = NULL,
    results_df = NULL,
    plots = NULL,
    sens_resultado = NULL,
    bootstrap_resultado = NULL,
    dominio = "salud"
  )
  
  # ==========================================================================
  # CARGA DE DATOS DE ARCHIVO
  # ==========================================================================
  observeEvent(input$file, {
    req(input$file)
    ext <- tools::file_ext(input$file$datapath)
    datos <- tryCatch({
      switch(ext,
             csv = read.csv(input$file$datapath),
             xlsx = readxl::read_excel(input$file$datapath),
             xls = readxl::read_excel(input$file$datapath),
             dta = haven::read_dta(input$file$datapath),
             NULL
      )
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      return(NULL)
    })
    
    if (is.null(datos)) {
      showNotification("Unsupported format.", type = "error")
      return()
    }
    
    rv$data <- datos
    rv$dominio <- "salud"
    showNotification(paste("Data loaded:", nrow(datos), "rows,", ncol(datos), "columns"), type = "message")
    
    vars <- names(datos)
    vars_num <- names(datos)[sapply(datos, is.numeric)]
    updateSelectInput(session, "graph_vars", choices = vars_num, selected = vars_num[1:min(5, length(vars_num))])
    updateSelectInput(session, "markov_var", choices = vars, selected = if ("Estado_Markov" %in% vars) "Estado_Markov" else vars[length(vars)])
    updateSelectInput(session, "group_var", choices = c("None (Global)", vars), selected = if ("Region" %in% vars) "Region" else "None (Global)")
    updateSelectInput(session, "alpha_vars", choices = vars_num, selected = vars_num[1:min(2, length(vars_num))])
  })
  
  # ==========================================================================
  # CARGA DE DATOS DE DEMOSTRACION (ENDES - Salud)
  # ==========================================================================
  observeEvent(input$load_demo, {
    showNotification("Loading demo data (ENDES - Health)...", type = "message")
    
    generar_y_cargar_demo <- function() {
      datos_demo <- generar_demo_data(1000, 123)
      rv$data <- datos_demo
      rv$dominio <- "salud"
      showNotification("Demo data generated successfully (n=1000).", type = "message")
    }
    
    if (file.exists("data/demo_data.csv")) {
      tryCatch({
        rv$data <- read.csv("data/demo_data.csv")
        rv$dominio <- "salud"
        showNotification("Demo data loaded from file.", type = "message")
      }, error = function(e) {
        showNotification(paste("Error reading file:", e$message), type = "error")
        generar_y_cargar_demo()
      })
    } else {
      showNotification("File not found. Generating demo data...", type = "message")
      generar_y_cargar_demo()
    }
    
    vars <- names(rv$data)
    vars_num <- names(rv$data)[sapply(rv$data, is.numeric)]
    updateSelectInput(session, "graph_vars", choices = vars_num, selected = vars_num[1:min(5, length(vars_num))])
    updateSelectInput(session, "markov_var", choices = vars, selected = "Estado_Markov")
    updateSelectInput(session, "group_var", choices = c("None (Global)", vars), selected = "Region")
    updateSelectInput(session, "alpha_vars", choices = vars_num, selected = c("Acceso_salud", "Tiene_seguro"))
  })
  
  # ==========================================================================
  # SALIDAS DE VISTA PREVIA
  # ==========================================================================
  output$data_preview <- renderDT({
    req(rv$data)
    datatable(
      head(rv$data, 10), 
      options = list(
        scrollX = TRUE, 
        dom = 't', 
        pageLength = 10,
        columnDefs = list(
          list(className = 'dt-center', targets = '_all')
        )
      ), 
      rownames = FALSE,
      class = 'display compact stripe hover'
    )
  })
  
  output$data_stats <- renderPrint({
    req(rv$data)
    cat("=== DATA STATISTICS ===\n\n")
    cat("Observations (n):", nrow(rv$data), "\n")
    cat("Variables:", ncol(rv$data), "\n")
    cat("Missing values:", sum(is.na(rv$data)), "\n\n")
    cat("--- Variable types ---\n")
    print(table(sapply(rv$data, class)))
    cat("\n--- First variables ---\n")
    print(names(rv$data)[1:min(10, ncol(rv$data))])
  })
  
  # ==========================================================================
  # UI DINAMICA
  # ==========================================================================
  output$graph_vars_ui <- renderUI({
    req(rv$data)
    vars_num <- names(rv$data)[sapply(rv$data, is.numeric)]
    if (length(vars_num) == 0) {
      return(div(class = "alert alert-warning", "No numeric variables found."))
    }
    checkboxGroupInput("graph_vars", "System variables:", choices = vars_num, 
                       selected = vars_num[1:min(5, length(vars_num))])
  })
  
  output$alpha_vars_ui <- renderUI({
    req(rv$data)
    vars_num <- names(rv$data)[sapply(rv$data, is.numeric)]
    if (length(vars_num) == 0) {
      return(helpText("No numeric variables available for alpha."))
    }
    selectInput("alpha_vars", "Variables for alpha (select 2-3):",
                choices = vars_num, multiple = TRUE,
                selected = vars_num[1:min(2, length(vars_num))])
  })
  
  output$markov_var_ui <- renderUI({
    req(rv$data)
    vars <- names(rv$data)
    default <- if ("Estado_Markov" %in% vars) {
      "Estado_Markov"
    } else if ("Estado_Participacion" %in% vars) {
      "Estado_Participacion"
    } else {
      vars[1]
    }
    selectInput("markov_var", "Markov state variable:", choices = vars, selected = default)
  })
  
  output$group_var_ui <- renderUI({
    req(rv$data)
    vars <- c("None (Global)", names(rv$data))
    default <- if ("Region" %in% names(rv$data)) {
      "Region"
    } else if ("Grupo" %in% names(rv$data)) {
      "Grupo"
    } else {
      "None (Global)"
    }
    selectInput("group_var", "Grouping variable:", choices = vars, selected = default)
  })
  
  output$suma_ponderadores <- renderText({
    total <- sum(input$w1, input$w2, input$w3)
    if (abs(total - 1) < 0.01) paste0(total, " OK") else paste0(total, " (must sum to 1)")
  })
  
  # ==========================================================================
  # FUNCION PARA DETECTAR ORDEN LOGICO DE ESTADOS
  # ==========================================================================
  detectar_orden_estados <- function(estados_unicos) {
    
    if (all(grepl("^E[0-9]+$", estados_unicos))) {
      return(estados_unicos[order(as.numeric(gsub("E", "", estados_unicos)))])
    }
    
    orden_aula <- c("Nunca", "Rara_vez", "A_veces", "Casi_siempre")
    if (any(orden_aula %in% estados_unicos)) {
      return(orden_aula[orden_aula %in% estados_unicos])
    }
    
    orden_bma <- c("Bajo", "Medio", "Alto")
    if (any(orden_bma %in% estados_unicos)) {
      return(orden_bma[orden_bma %in% estados_unicos])
    }
    
    if (all(grepl("^[0-9]+$", estados_unicos))) {
      return(as.character(sort(as.numeric(estados_unicos))))
    }
    
    if (all(grepl("^Nivel [0-9]+$", estados_unicos))) {
      return(estados_unicos[order(as.numeric(gsub("Nivel ", "", estados_unicos)))])
    }
    
    warning("Could not determine logical order of states. Using alphabetical order.")
    return(sort(estados_unicos))
  }
  
  # ==========================================================================
  # FUNCION PARA CALCULAR ALPHA
  # ==========================================================================
  calcular_alpha <- function(sub) {
    alpha <- 0.60
    
    if (input$alpha_mode == "manual") {
      alpha <- input$alpha_manual
      if (is.na(alpha) || alpha < 0 || alpha > 1) alpha <- 0.60
      
    } else if (input$alpha_mode == "vars" && !is.null(input$alpha_vars) && length(input$alpha_vars) > 0) {
      vars_alpha <- intersect(input$alpha_vars, names(sub))
      if (length(vars_alpha) > 0) {
        alpha_values <- sapply(vars_alpha, function(var) {
          vals <- sub[[var]]
          if (is.numeric(vals)) {
            rng <- range(vals, na.rm = TRUE)
            if (rng[2] > rng[1]) {
              (vals - rng[1]) / (rng[2] - rng[1])
            } else {
              rep(0.5, length(vals))
            }
          } else {
            rep(0.5, length(vals))
          }
        })
        if (is.matrix(alpha_values)) {
          alpha <- mean(rowMeans(alpha_values, na.rm = TRUE), na.rm = TRUE)
        } else {
          alpha <- mean(alpha_values, na.rm = TRUE)
        }
        if (is.na(alpha) || is.nan(alpha)) alpha <- 0.60
      }
      
    } else {
      if ("Acceso_salud" %in% names(sub)) {
        alpha <- mean(sub$Acceso_salud, na.rm = TRUE)
      } else if (all(c("Peer_Influence", "Teacher_Encouragement") %in% names(sub))) {
        peer_norm <- (sub$Peer_Influence - 1) / 3
        teacher_norm <- sub$Teacher_Encouragement / 10
        alpha <- mean(rowMeans(cbind(peer_norm, teacher_norm), na.rm = TRUE), na.rm = TRUE)
      } else {
        vars_num <- names(sub)[sapply(sub, is.numeric)]
        if (length(vars_num) > 0) {
          vals <- sub[[vars_num[1]]]
          rng <- range(vals, na.rm = TRUE)
          if (rng[2] > rng[1]) {
            alpha <- mean((vals - rng[1]) / (rng[2] - rng[1]), na.rm = TRUE)
          }
        }
      }
      if (is.na(alpha) || is.nan(alpha)) alpha <- 0.60
    }
    
    return(max(0, min(1, alpha)))
  }
  
  # ==========================================================================
  # FUNCION PARA PROCESAR TASAS DE AVANCE
  # ==========================================================================
  procesar_tasas_avance <- function(orden_estados) {
    tasas_avance <- NULL
    if (!is.null(input$tasas_avance_input) && input$tasas_avance_input != "") {
      tasas_raw <- strsplit(trimws(input$tasas_avance_input), ",")[[1]]
      tasas_raw <- trimws(tasas_raw)
      tasas_num <- suppressWarnings(as.numeric(tasas_raw))
      if (!any(is.na(tasas_num)) && length(tasas_num) > 0) {
        m <- length(orden_estados)
        if (length(tasas_num) == (m - 1)) {
          tasas_avance <- tasas_num
        } else {
          showNotification(
            paste("Expected", m - 1, "rates for", m, "states. Using heuristic."),
            type = "warning"
          )
        }
      }
    }
    return(tasas_avance)
  }
  
  # ==========================================================================
  # ANALISIS DE BOOTSTRAP PARA RENDER
  # ==========================================================================
  observeEvent(input$run_bootstrap, {
    req(rv$data)
    
    if (is.null(input$graph_vars) || length(input$graph_vars) < 5) {
      showNotification("Select at least 5 variables for the graph.", type = "error")
      return()
    }
    
    showNotification("Running bootstrap analysis. This may take 30-60 seconds...", 
                     type = "message", duration = 5)
    
    withProgress(message = 'Running bootstrap...', value = 0.1, {
      
      datos <- rv$data
      variables <- input$graph_vars
      n_boot <- input$boot_n_sim
      
      umbrales <- c(0.05, 0.07, 0.10, 0.12, 0.15, 0.20)
      
      incProgress(0.2, detail = "Processing thresholds...")
      
      resultado_boot <- seleccionar_umbral_bootstrap(
        datos = datos,
        variables = variables,
        umbrales = umbrales,
        n_boot = n_boot,
        seed = 123,
        criterio_jaccard = 0.70
      )
      
      incProgress(0.9, detail = "Generating results...")
      
      rv$bootstrap_resultado <- resultado_boot
      
      incProgress(1.0, detail = "Completed!")
    })
    
    showNotification("Bootstrap completed successfully!", type = "message")
  })
  
  # ==========================================================================
  # SALIDA DE BOOTSTRAP - MANEJO DE ERRORES
  # ==========================================================================
  output$bootstrap_results_ui <- renderUI({
    req(rv$bootstrap_resultado)
    
    res <- rv$bootstrap_resultado
    
    if (is.null(res$umbral_optimo) || is.na(res$umbral_optimo)) {
      return(div(
        class = "alert alert-warning",
        icon("exclamation-triangle"),
        " Could not determine an optimal threshold. Try with more replicates or a different range."
      ))
    }
    
    umbral_optimo <- round(res$umbral_optimo, 2)
    
    nodo_original <- res$nodo_original
    if (is.null(nodo_original) || is.na(nodo_original)) nodo_original <- "No disponible"
    
    nodo_mas_frecuente <- res$nodo_mas_frecuente
    if (is.null(nodo_mas_frecuente) || is.na(nodo_mas_frecuente)) nodo_mas_frecuente <- "No disponible"
    
    pct_estable <- res$pct_nodo_mas_frecuente
    if (is.null(pct_estable) || is.na(pct_estable) || !is.numeric(pct_estable)) {
      pct_estable <- 0
    }
    pct_estable <- round(pct_estable, 1)
    
    jaccard_optimo <- NA
    if (!is.null(res$resultados) && nrow(res$resultados) > 0) {
      idx <- which(res$resultados$Umbral == res$umbral_optimo)
      if (length(idx) > 0 && !is.na(res$resultados$Jaccard_Promedio[idx[1]])) {
        jaccard_optimo <- round(res$resultados$Jaccard_Promedio[idx[1]], 3)
      }
    }
    
    html_parts <- c(
      "<div style='background-color: #f0f8ff; padding: 12px; border-radius: 5px; border-left: 4px solid #3498DB;'>",
      "<b>Optimal threshold:</b> <span style='color:#C0392B;font-size:20px;font-weight:bold;'>", 
      umbral_optimo, "</span><br>",
      "<b>Original node:</b> ", nodo_original, "<br>",
      "<b>Most frequent node (bootstrap):</b> ", nodo_mas_frecuente, 
      " (", pct_estable, "% stability)<br>"
    )
    
    if (!is.na(jaccard_optimo)) {
      html_parts <- c(html_parts, "<b>Jaccard similarity:</b> ", jaccard_optimo)
    } else {
      html_parts <- c(html_parts, "<b>Jaccard similarity:</b> Not available")
    }
    
    html_parts <- c(html_parts, "</div>")
    
    HTML(paste(html_parts, collapse = ""))
  })
  
  # ==========================================================================
  # GRAFICO DE BOOTSTRAP
  # ==========================================================================
  output$bootstrap_plot <- renderPlot({
    req(rv$bootstrap_resultado)
    
    res <- rv$bootstrap_resultado
    df <- res$resultados
    
    if (is.null(df) || nrow(df) == 0) {
      return(ggplot() + theme_void() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No hay datos para el grafico",
                        size = 5, color = "#7F8C8D"))
    }
    
    if (!"Jaccard_Promedio" %in% names(df)) {
      return(ggplot() + theme_void() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "Datos incompletos para el grafico",
                        size = 5, color = "#7F8C8D"))
    }
    
    ggplot(df, aes(x = Umbral)) +
      geom_line(aes(y = Jaccard_Promedio, color = "Jaccard"), linewidth = 1.2) +
      geom_point(aes(y = Jaccard_Promedio, color = "Jaccard"), size = 2) +
      geom_vline(xintercept = res$umbral_optimo, linetype = "dashed", color = "#C0392B", linewidth = 1) +
      annotate("text", x = res$umbral_optimo + 0.008, y = max(df$Jaccard_Promedio, na.rm = TRUE) * 0.9,
               label = paste0("Optimal = ", round(res$umbral_optimo, 2)),
               color = "#C0392B", fontface = "bold", hjust = 0) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      labs(title = "Graph stability by threshold (Jaccard)",
           x = "Correlation threshold", y = "Jaccard similarity") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")
  })
  
  # ==========================================================================
  # EJECUTAR MPCS
  # ==========================================================================
  observeEvent(input$run_mpcs, {
    req(rv$data)
    
    if (is.null(input$graph_vars) || length(input$graph_vars) < 5) {
      showNotification("Select at least 5 variables.", type = "error")
      output$validation_msg <- renderUI({ 
        tags$div(class = "alert alert-danger mt-2", 
                 icon("exclamation-triangle"), 
                 " Select at least 5 variables.")
      })
      return()
    }
    output$validation_msg <- renderUI({ NULL })
    
    if (is.null(input$markov_var) || input$markov_var == "") {
      showNotification("Select a Markov state variable.", type = "error")
      output$validation_msg <- renderUI({ 
        tags$div(class = "alert alert-danger mt-2", 
                 icon("exclamation-triangle"), 
                 " Select a Markov state variable.")
      })
      return()
    }
    
    estados <- rv$data[[input$markov_var]]
    if (length(unique(na.omit(estados))) < 3) {
      showNotification("Markov variable must have at least 3 states.", type = "error")
      output$validation_msg <- renderUI({ 
        tags$div(class = "alert alert-danger mt-2", 
                 icon("exclamation-triangle"), 
                 " Markov variable must have at least 3 states.")
      })
      return()
    }
    output$validation_msg <- renderUI({ NULL })
    
    withProgress(message = 'Running MPCS...', value = 0, {
      data_analysis <- rv$data
      
      if (input$group_var == "None (Global)") {
        data_analysis$Group <- "Global"
        grupos <- "Global"
      } else {
        data_analysis$Group <- data_analysis[[input$group_var]]
        grupos <- unique(data_analysis$Group)
        grupos <- grupos[!is.na(grupos) & grupos != ""]
      }
      
      results_list <- list()
      
      for (i in seq_along(grupos)) {
        g <- grupos[i]
        incProgress(1 / length(grupos), detail = paste("Processing group:", g))
        
        sub <- data_analysis[data_analysis$Group == g, ]
        
        if (nrow(sub) < 30) {
          showNotification(paste("Group", g, "has less than 30 observations. Skipping."), type = "message")
          next
        }
        
        # --- Grafo ---
        graph_data <- sub[, input$graph_vars, drop = FALSE]
        graph_res <- calcular_grafo(graph_data, input$graph_vars, input$threshold)
        
        if (is.null(graph_res$graph) || vcount(graph_res$graph) < 2) {
          showNotification(paste("Group", g, "has insufficient connections in the graph. Using default value."), type = "message")
          graph_res$score <- 0.5
          graph_res$optimal_node <- "No node"
          graph_res$graph <- NULL
        }
        
        alpha_grupo <- calcular_alpha(sub)
        
        # --- Detectar orden de estados ---
        estados_unicos <- unique(sub[[input$markov_var]])
        estados_unicos <- estados_unicos[!is.na(estados_unicos) & estados_unicos != ""]
        orden_estados <- detectar_orden_estados(estados_unicos)
        
        # --- Procesar tasas de avance ---
        tasas_avance <- procesar_tasas_avance(orden_estados)
        
        # --- Markov ---
        markov_res <- calcular_markov(
          sub[[input$markov_var]], 
          orden_estados = orden_estados,
          tasas_avance = tasas_avance,
          umbral_objetivo = 0.50, 
          horizonte = 10
        )
        
        # --- Juegos ---
        games_res <- calcular_juegos(markov_res$mat, input$R_factor, alpha = alpha_grupo)
        
        # --- Indice integrado ---
        index_res <- calcular_indice(
          I_grafo = graph_res$score,
          I_markov = markov_res$score,
          I_juegos = games_res$score,
          w1 = input$w1,
          w2 = input$w2,
          w3 = input$w3,
          R_factor = input$R_factor
        )
        
        results_list[[as.character(g)]] <- list(
          grupo = g,
          n = nrow(sub),
          nodo_optimo = if (is.null(graph_res$optimal_node)) "No node" else graph_res$optimal_node,
          I_grafo = graph_res$score,
          I_markov = markov_res$score,
          I_juegos = games_res$score,
          alpha = alpha_grupo,
          I_MPCS = index_res$I_MPCS,
          k = index_res$k,
          tipo = index_res$nudge_type,
          graph = graph_res$graph,
          graph_data = graph_data,
          markov_mat = markov_res$mat,
          sim_base = markov_res$sim_base,
          dist_actual = markov_res$dist_actual,
          T_base = markov_res$T_base,
          estados = colnames(markov_res$mat)
        )
      }
      
      if (length(results_list) == 0) {
        showNotification("Could not process any group. Verify that at least one group has >30 observations.", type = "error")
        return()
      }
      
      # --- Calcular sim_nudge para cada grupo ---
      for (g in names(results_list)) {
        r <- results_list[[g]]
        if (!is.null(r$markov_mat) && !is.null(r$sim_base)) {
          P_n <- r$markov_mat
          
          k_nudge <- r$k
          if (is.na(k_nudge) || k_nudge < 0 || k_nudge > 1) k_nudge <- 0.4
          
          for (i in 1:(nrow(P_n)-1)) {
            av <- P_n[i, i] * k_nudge
            P_n[i, i] <- P_n[i, i] - av
            P_n[i, i+1] <- P_n[i, i+1] + av
            P_n[i, ] <- P_n[i, ] / sum(P_n[i, ])
          }
          
          sim_nudge <- r$sim_base
          for (j in 2:nrow(sim_nudge)) {
            sim_nudge[j, ] <- sim_nudge[j-1, ] %*% P_n
          }
          results_list[[g]]$sim_nudge <- sim_nudge
        }
      }
      
      rv$results <- results_list
      
      results_df <- do.call(rbind, lapply(results_list, function(r) {
        data.frame(
          Grupo = r$grupo,
          n = r$n,
          Alpha = round(r$alpha, 3),
          I_Markov_H10 = round(r$I_markov, 4),
          I_Games = round(r$I_juegos, 4),
          I_MPCS = round(r$I_MPCS, 4),
          Nodo_Optimo = r$nodo_optimo,
          k = round(r$k, 4),
          Tipo_Nudge = r$tipo,
          stringsAsFactors = FALSE
        )
      }))
      
      rv$results_df <- results_df
      
      rv$plots <- generate_plots(
        data = data_analysis,
        graph_vars = input$graph_vars,
        markov_var = input$markov_var,
        results = results_df,
        results_list = results_list,
        threshold = input$threshold,
        dominio = rv$dominio
      )
      
      showNotification(paste("MPCS executed successfully for", nrow(results_df), "groups."), type = "message")
    })
  })
  
  # ==========================================================================
  # FUNCIONES PARA GRAFICOS
  # ==========================================================================
  generate_plots <- function(data, graph_vars, markov_var, results, results_list, threshold, dominio = "salud") {
    p_graph <- NULL
    p_states <- NULL
    p_markov <- NULL
    p_rank <- NULL
    p_arbol <- NULL
    p_juego <- NULL
    
    if (is.null(data) || nrow(data) == 0 || is.null(results) || nrow(results) == 0) {
      return(list(graph = NULL, states = NULL, markov = NULL, rank = NULL, 
                  arbol = NULL, juego = NULL))
    }
    
    first_group <- results$Grupo[1]
    r <- results_list[[as.character(first_group)]]
    
    # --- Grafo ---
    if (!is.null(r$graph) && vcount(r$graph) > 0) {
      V(r$graph)$color <- ifelse(V(r$graph)$name == r$nodo_optimo, "#C0392B", "#F0DFC0")
      V(r$graph)$size <- ifelse(V(r$graph)$name == r$nodo_optimo, 25, 15)
      p_graph <- function() {
        plot(r$graph, layout = layout_with_fr(r$graph), vertex.label.cex = 0.8,
             vertex.label.color = "black", vertex.label.dist = 2, edge.color = "gray60",
             edge.width = 1.5, main = paste("Behavioral Graph -", first_group), cex.main = 0.9)
        legend("topright", legend = c("Optimal node", "Other nodes"), 
               fill = c("#C0392B", "#F0DFC0"), cex = 0.8, bty = "n")
      }
    }
    
    # --- Distribucion de estados ---
    if (!is.null(markov_var) && markov_var %in% names(data)) {
      p_states <- function() {
        if (!"Group" %in% names(data)) data$Group <- "Global"
        plot_data <- data[!is.na(data$Group) & !is.na(data[[markov_var]]), ]
        if (nrow(plot_data) == 0) {
          return(ggplot() + theme_void() + 
                   annotate("text", x = 0.5, y = 0.5, label = "No data"))
        }
        ggplot(plot_data, aes(x = .data[["Group"]], fill = .data[[markov_var]])) + 
          geom_bar(position = "fill") +
          scale_fill_brewer(palette = "Set2") + theme_minimal() +
          labs(x = "Group", y = "Proportion", fill = "State") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
    }
    
    # --- Ranking ---
    if (!is.null(results) && nrow(results) > 0) {
      p_rank <- function() {
        ggplot(results, aes(x = reorder(Grupo, I_MPCS), y = I_MPCS, fill = Tipo_Nudge)) +
          geom_col(width = 0.7) + coord_flip() + theme_minimal(base_size = 12) +
          labs(x = "Group", y = "MPCS Index", fill = "Nudge Type") +
          scale_fill_manual(
            values = c("Informational" = "#74B3CE", "Structural" = "#2E86AB",
                       "Normative" = "#E84855", "Systemic multi-nudge" = "#1A3A5C")
          ) +
          geom_text(aes(label = round(I_MPCS, 3)), hjust = -0.2, size = 3.5) +
          theme(legend.position = "bottom")
      }
    }
    
    # --- Trayectorias de Markov ---
    p_markov <- function() {
      if (!is.null(r$sim_base) && nrow(r$sim_base) > 0) {
        y_label <- if (dominio == "educacion" || any(grepl("Participacion", colnames(r$sim_base)))) {
          "Participation probability"
        } else {
          "Occupation probability"
        }
        return(graficar_trayectorias_markov(
          sim_base = r$sim_base, 
          sim_nudge = r$sim_nudge,
          estados = colnames(r$sim_base),
          titulo = paste("Markov Trajectories -", first_group),
          y_label = y_label
        ))
      }
      return(ggplot() + theme_void() + 
               annotate("text", x = 0.5, y = 0.5, label = "No data for trajectories"))
    }
    
    # --- Arbol de Markov ---
    p_arbol <- function() {
      if (!is.null(r$markov_mat)) {
        estados_arbol <- if (!is.null(r$estados) && length(r$estados) == ncol(r$markov_mat)) {
          r$estados
        } else if (!is.null(colnames(r$markov_mat))) {
          colnames(r$markov_mat)
        } else {
          paste0("E", 1:ncol(r$markov_mat))
        }
        
        return(graficar_arbol_markov(
          P = r$markov_mat,
          estados = estados_arbol,
          umbral_prob = 0.05,
          titulo = paste("Markov Transitions -", first_group)
        ))
      }
      return(ggplot() + theme_void() + 
               annotate("text", x = 0.5, y = 0.5, label = "No data for Markov tree"))
    }
    
    # --- Juego evolutivo ---
    p_juego <- function() {
      if (!is.null(r$alpha)) {
        return(graficar_juego_evolutivo(
          alpha = r$alpha,
          p_star = r$I_juegos,
          titulo = paste("Replicator Dynamics -", first_group)
        ))
      }
      return(ggplot() + theme_void() + 
               annotate("text", x = 0.5, y = 0.5, label = "No data for game theory"))
    }
    
    return(list(graph = p_graph, states = p_states, markov = p_markov, rank = p_rank,
                arbol = p_arbol, juego = p_juego))
  }
  
  # ==========================================================================
  # SALIDAS DE RESULTADOS
  # ==========================================================================
  output$results_table <- renderDT({
    req(rv$results_df)
    datatable(rv$results_df, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE) %>%
      formatRound(columns = c("I_Markov_H10", "I_Games", "I_MPCS", "k"), digits = 4)
  })
  
  output$plot_graph <- renderPlot({
    req(rv$plots)
    if (!is.null(rv$plots$graph)) rv$plots$graph() else {
      plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(0, 0, "Could not generate graph.")
    }
  })
  
  output$plot_states <- renderPlot({
    req(rv$plots)
    if (!is.null(rv$plots$states)) rv$plots$states() else {
      ggplot() + theme_void() + 
        annotate("text", x = 0.5, y = 0.5, label = "Could not generate chart.")
    }
  })
  
  output$plot_markov <- renderPlot({
    req(rv$plots)
    if (!is.null(rv$plots$markov)) rv$plots$markov() else {
      ggplot() + theme_void() + 
        annotate("text", x = 0.5, y = 0.5, label = "Could not generate chart.")
    }
  })
  
  output$plot_ranking <- renderPlot({
    req(rv$plots)
    if (!is.null(rv$plots$rank)) rv$plots$rank() else {
      ggplot() + theme_void() + 
        annotate("text", x = 0.5, y = 0.5, label = "Could not generate chart.")
    }
  })
  
  output$plot_arbol_markov <- renderPlot({
    req(rv$plots)
    if (!is.null(rv$plots$arbol)) {
      rv$plots$arbol()
    } else {
      ggplot() + theme_void() + 
        annotate("text", x = 0.5, y = 0.5, label = "Could not generate Markov tree.")
    }
  })
  
  output$plot_juego_evolutivo <- renderPlot({
    req(rv$plots)
    if (!is.null(rv$plots$juego)) rv$plots$juego() else {
      ggplot() + theme_void() + 
        annotate("text", x = 0.5, y = 0.5, label = "Could not generate game theory chart.")
    }
  })
  
  # ==========================================================================
  # ANALISIS DE SENSIBILIDAD
  # ==========================================================================
  observeEvent(input$run_sensitivity, {
    
    if (is.null(rv$results_df) || nrow(rv$results_df) == 0) {
      showNotification("First run MPCS to get results.", type = "error")
      return()
    }
    
    top_group_idx <- which.max(rv$results_df$I_MPCS)
    top_group <- rv$results_df$Grupo[top_group_idx]
    
    r <- rv$results[[as.character(top_group)]]
    
    if (is.null(r)) {
      showNotification("No results found for selected group.", type = "error")
      return()
    }
    
    I_Grafo <- r$I_grafo
    I_Markov <- r$I_markov
    I_Juegos <- r$I_juegos
    
    if (is.na(I_Grafo)) I_Grafo <- 0.5
    if (is.na(I_Markov)) I_Markov <- 0.5
    if (is.na(I_Juegos)) I_Juegos <- 0.5
    
    showNotification(
      paste("Running sensitivity analysis for group:", top_group),
      type = "message",
      duration = 2
    )
    
    withProgress(message = 'Analyzing sensitivity...', value = 0, {
      
      resultado_sens <- analisis_sensibilidad(
        I_Grafo = I_Grafo,
        I_Markov = I_Markov,
        I_Juegos = I_Juegos,
        n_sim = input$sens_n_sim,
        R = input$R_factor,
        escala = 1.5,
        seed = 123
      )
      
      incProgress(0.8, detail = "Generating charts...")
      
      rv$sens_resultado <- resultado_sens
      
      incProgress(1, detail = "Completed")
    })
    
    showNotification("Sensitivity analysis completed.", type = "message")
  })
  
  output$sens_plot <- renderPlot({
    req(rv$sens_resultado)
    
    df <- data.frame(
      I_MPCS = rv$sens_resultado$I_MPCS_sim,
      Tipo = rv$sens_resultado$tipo_sim
    )
    
    colores <- c(
      "Informational" = "#74B3CE",
      "Structural" = "#2E86AB",
      "Normative" = "#E84855",
      "Systemic multi-nudge" = "#1A3A5C"
    )
    
    tipos_presentes <- unique(df$Tipo)
    colores_filtrados <- colores[names(colores) %in% tipos_presentes]
    
    ggplot(df, aes(x = I_MPCS, fill = Tipo)) +
      geom_histogram(bins = 50, color = "white", alpha = 0.85) +
      geom_vline(xintercept = rv$sens_resultado$I_MPCS_default,
                 color = "#C0392B", linetype = "dashed", linewidth = 1.2) +
      geom_vline(xintercept = rv$sens_resultado$ic_inf,
                 color = "#2C3E50", linetype = "dotted", linewidth = 0.8, alpha = 0.5) +
      geom_vline(xintercept = rv$sens_resultado$ic_sup,
                 color = "#2C3E50", linetype = "dotted", linewidth = 0.8, alpha = 0.5) +
      annotate("text", 
               x = rv$sens_resultado$I_MPCS_default + 0.008, 
               y = max(table(cut(df$I_MPCS, breaks = 50))) * 0.9,
               label = paste0("Default weights\nI_MPCS = ", 
                             round(rv$sens_resultado$I_MPCS_default, 4)),
               hjust = 0, size = 3.5, color = "#C0392B") +
      scale_fill_manual(values = colores_filtrados) +
      labs(
        title = "Sensitivity Analysis of MPCS Index",
        subtitle = paste0(rv$sens_resultado$n_sim, " random combinations of weights w1, w2, w3"),
        x = "MPCS Index",
        y = "Frequency",
        fill = "Nudge type"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
  })
  
  output$sens_summary <- renderPrint({
    req(rv$sens_resultado)
    
    s <- rv$sens_resultado
    
    cat("=== SENSITIVITY ANALYSIS RESULTS ===\n\n")
    cat(sprintf("Number of simulations: %d\n", s$n_sim))
    cat(sprintf("Mean I_MPCS: %.4f\n", s$stats$media))
    cat(sprintf("Standard deviation: %.4f\n", s$stats$sd))
    cat(sprintf("Coefficient of variation: %.1f%%\n", s$stats$cv))
    cat(sprintf("95%% Confidence interval: [%.4f, %.4f]\n", 
                s$ic_inf, s$ic_sup))
    cat(sprintf("I_MPCS with default weights: %.4f\n", s$I_MPCS_default))
    cat(sprintf("Default nudge type: %s\n", s$tipo_default))
    cat(sprintf("\nMost frequent nudge type in simulations: %s (%.1f%%)\n",
                s$tipo_mas_frecuente, s$pct_mas_frecuente))
    cat(sprintf("Agreement with default type: %.1f%%\n", 
                s$pct_coincidencia))
    
    cat("\n=== INTERPRETATION ===\n")
    if (s$pct_coincidencia >= 70) {
      cat("The nudge recommendation is HIGHLY ROBUST.\n")
      cat("  The nudge type remains in more than 70% of simulations.\n")
    } else if (s$pct_coincidencia >= 50) {
      cat("The nudge recommendation is ROBUST.\n")
      cat("  The nudge type remains in more than 50% of simulations.\n")
    } else {
      cat("The nudge recommendation is SENSITIVE to weights.\n")
      cat("  Consider reviewing the weighting or conducting additional analysis.\n")
    }
  })
  
  # ==========================================================================
  # INTERPRETACION
  # ==========================================================================
  output$interpretation_text <- renderUI({
    req(rv$results_df)
    
    top_group <- rv$results_df[which.max(rv$results_df$I_MPCS), ]
    bottom_group <- rv$results_df[which.min(rv$results_df$I_MPCS), ]
    
    HTML(paste0(
      "<div class='well'>",
      "<p><b>Results Summary</b></p>",
      "<p>Analyzed <b>", nrow(rv$results_df), " groups</b> with a total of <b>", 
      sum(rv$results_df$n), " observations</b>.</p>",
      "<p><b>Domain:</b> ",
      ifelse(any(grepl("Participacion", names(rv$data))), "Education - Classroom participation", "Health - Medication adherence"),
      "</p>",
      "<hr>",
      "<p><b>Highest priority group:</b> <span style='color:#C0392B;font-weight:bold;'>", 
      top_group$Grupo, "</span> (I_MPCS = ", round(top_group$I_MPCS, 4), ")</p>",
      "<p>The optimal intervention node is <b>", top_group$Nodo_Optimo, "</b>.</p>",
      "<p>Recommended nudge type: <b>", top_group$Tipo_Nudge, 
      "</b> with intensity k = ", round(top_group$k, 4), ".</p>",
      "<p><b>Estimated critical mass:</b> ", round(top_group$I_Games * 100, 1), 
      "% of adopters needed for self-sustaining change.</p>",
      "<hr>",
      "<p><b>Lowest priority group:</b> <span style='color:#2ECC71;font-weight:bold;'>", 
      bottom_group$Grupo, "</span> (I_MPCS = ", round(bottom_group$I_MPCS, 4), ")</p>",
      "<p><b>Recommendation:</b> Focus efforts on the highest priority group (", 
      top_group$Grupo, ") and apply a <b>", top_group$Tipo_Nudge, 
      "</b> nudge with intensity ", round(top_group$k * 100, 1), "%.</p>",
      "</div>"
    ))
  })
  
  # ==========================================================================
  # DESCARGAS
  # ==========================================================================
  output$download_csv <- downloadHandler(
    filename = function() {
      paste0("MPCS_Results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(rv$results_df)
      write.csv(rv$results_df, file, row.names = FALSE)
    }
  )
  
  output$download_report <- downloadHandler(
    filename = function() {
      paste0("MPCS_Report_", Sys.Date(), ".html")
    },
    content = function(file) {
      if (is.null(rv$results_df) || nrow(rv$results_df) == 0) {
        showNotification("No results to generate report.", type = "error")
        return()
      }
      
      top_group <- rv$results_df[which.max(rv$results_df$I_MPCS), ]
      
      html_lines <- c(
        "<!DOCTYPE html>",
        "<html>",
        "<head>",
        "<meta charset='UTF-8'>",
        "<title>MPCS Report</title>",
        "<style>",
        "body { font-family: Arial, sans-serif; margin: 40px; max-width: 1000px; margin-left: auto; margin-right: auto; }",
        "h1 { color: #2C3E50; border-bottom: 3px solid #3498DB; padding-bottom: 10px; }",
        "h2 { color: #34495E; border-bottom: 2px solid #3498DB; padding-bottom: 5px; margin-top: 30px; }",
        "table { border-collapse: collapse; width: 100%; margin: 20px 0; }",
        "th { background-color: #3498DB; color: white; padding: 12px; text-align: left; }",
        "td { padding: 10px; border-bottom: 1px solid #ddd; }",
        "tr:nth-child(even) { background-color: #f9f9f9; }",
        "tr:hover { background-color: #f5f5f5; }",
        ".highlight { background-color: #FFEAA7; padding: 15px; border-radius: 5px; border-left: 5px solid #F39C12; }",
        ".footer { margin-top: 50px; font-size: 12px; color: #7F8C8D; border-top: 2px solid #ddd; padding-top: 15px; text-align: center; }",
        "</style>",
        "</head>",
        "<body>",
        "<h1>MPCS Report</h1>",
        "<p><b>Generation date:</b> ", format(Sys.Date(), "%B %d, %Y"), 
        " at ", format(Sys.time(), "%H:%M"), "</p>",
        "<p><b>Domain:</b> ",
        ifelse(any(grepl("Participacion", names(rv$data))), "Education - Classroom participation", "Health - Medication adherence"),
        "</p>",
        "<hr>",
        "<h2>1. Results Summary</h2>",
        "<p><b>Groups analyzed:</b> ", nrow(rv$results_df), "</p>",
        "<p><b>Total observations:</b> ", sum(rv$results_df$n), "</p>",
        "<br>",
        "<h2>2. Results Table</h2>",
        "<table>",
        "<thead><tr><th>Group</th><th>n</th><th>Alpha</th><th>I_Markov_H10</th><th>I_Games</th><th>I_MPCS</th><th>Optimal Node</th><th>k</th><th>Nudge Type</th></tr></thead>",
        "<tbody>"
      )
      
      for (i in 1:nrow(rv$results_df)) {
        row <- rv$results_df[i, ]
        is_top <- row$I_MPCS == max(rv$results_df$I_MPCS)
        bg_color <- if (is_top) " style='background-color: #FFEAA7;'" else ""
        icono <- if (is_top) " [PRIORITY]" else ""
        html_lines <- c(
          html_lines,
          paste0("<tr", bg_color, ">"),
          paste0("<td><b>", row$Grupo, icono, "</b></td>"),
          paste0("<td>", row$n, "</td>"),
          paste0("<td>", round(row$Alpha, 3), "</td>"),
          paste0("<td>", round(row$I_Markov_H10, 4), "</td>"),
          paste0("<td>", round(row$I_Games, 4), "</td>"),
          paste0("<td><b>", round(row$I_MPCS, 4), "</b></td>"),
          paste0("<td>", row$Nodo_Optimo, "</td>"),
          paste0("<td>", round(row$k, 4), "</td>"),
          paste0("<td>", row$Tipo_Nudge, "</td>"),
          "</tr>"
        )
      }
      
      html_lines <- c(
        html_lines,
        "</tbody></table>",
        "<h2>3. Interpretation</h2>",
        "<div class='highlight'>",
        "<p><b>Highest priority group for intervention:</b></p>",
        "<ul>",
        "<li><b>Group:</b> ", top_group$Grupo, "</li>",
        "<li><b>MPCS Index:</b> ", round(top_group$I_MPCS, 4), "</li>",
        "<li><b>Optimal intervention node:</b> ", top_group$Nodo_Optimo, "</li>",
        "<li><b>Recommended intensity (k):</b> ", round(top_group$k, 4), "</li>",
        "<li><b>Recommended nudge type:</b> <b>", top_group$Tipo_Nudge, "</b></li>",
        "<li><b>Critical mass:</b> ", round(top_group$I_Games * 100, 1), "% of adopters needed</li>",
        "</ul>",
        "</div>",
        "<h2>4. Citation</h2>",
        "<p>MPCS: A method for the systematic design of behavioral nudges... (Zela Llanque, 2026).</p>",
        "<div class='footer'>",
        "<p>Report automatically generated by <b>MPCS Calculator</b></p>",
        "<p>Repository: <a href='https://github.com/Izela-meth/MPCS_APP' target='_blank'>",
        "https://github.com/Izela-meth/MPCS_APP</a></p>",
        "</div>",
        "</body>",
        "</html>"
      )
      
      writeLines(html_lines, file)
      showNotification("HTML report generated successfully.", type = "message")
    }
  )
}

# ==============================================================================
# EJECUTAR LA APLICACION CON FOOTER
# ==============================================================================

ui <- tagList(
  ui,
  tags$style(HTML("
    .navbar + .container-fluid {
      min-height: calc(100vh - 200px);
    }
    .mpcs-footer {
      background-color: #f8f9fa;
      padding: 15px;
      text-align: center;
      border-top: 1px solid #dee2e6;
      margin-top: 30px;
      font-size: 14px;
      color: #6c757d;
      width: 100%;
    }
    .mpcs-footer a {
      color: #2C3E50;
      text-decoration: none;
    }
    .mpcs-footer a:hover {
      text-decoration: underline;
    }
    .well-sens {
      background-color: #f8f9fa;
      border-radius: 5px;
      padding: 15px;
    }
  ")),
  tags$footer(
    class = "mpcs-footer",
    "MPCS Calculator v1.0 - Developed by the MPCS Team | ",
    tags$a(href = "https://github.com/Izela-meth/MPCS_APP", target = "_blank", "GitHub"),
    " | ",
    tags$a(href = "mailto:izela@unsa.edu.pe", "Contact")
  )
)

# ==============================================================================
# EJECUTAR LA APLICACION
# ==============================================================================
shinyApp(
  ui = ui,
  server = server
)
