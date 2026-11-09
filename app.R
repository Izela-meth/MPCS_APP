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

# --- Cargar funciones modulares (Asegurate de que mpcs_functions.R esté actualizado) ---
source("functions/mpcs_functions.R", local = TRUE)

# ==============================================================================
# GENERAR DATOS DE DEMOSTRACION INICIALES
# ==============================================================================
if (!file.exists("data/demo_data3.csv")) {
  dir.create("data", showWarnings = FALSE)
  # La nueva función generar_demo_data ya crea datos con estructura real (no ruido)
  datos_demo <- generar_demo_data(n = 1000, seed = 123)
  write.csv(datos_demo, "data/demo_data3.csv", row.names = FALSE)
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
          p("Load simulated health data to test the application."),
          actionButton("load_demo", "Load Demo Data", 
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
          uiOutput("group_var_ui")
        )
      ),
      column(
        width = 6,
        wellPanel(
          h4("Markov configuration"),
          uiOutput("markov_var_ui"),
          hr(),
          h4("Adjustable parameters"),
          sliderInput("threshold", "Correlation threshold",
                      min = 0.05, max = 0.30, value = 0.10, step = 0.01),
          
          hr(),
          h4("Markov Transition Rates"),
          radioButtons("markov_rates_mode", "Method to define transition probabilities:",
             choices = c(
               "Estimate from data (Conversion Ratio)" = "estimate",
               "Use calibrated values (Manuscript)" = "calibrated",
               "Enter manually (4 rates)" = "manual"
             ),
             selected = "calibrated"),
          
          conditionalPanel(
            condition = "input.markov_rates_mode == 'estimate'",
            div(class = "alert alert-warning", style = "font-size: 0.85em;",
                icon("exclamation-triangle"),
                HTML("<b>Warning:</b> Monte Carlo simulations show that the conversion ratio rule 
                     introduces a massive error (approx. -0.22) in populations dominated by 
                     undiagnosed individuals (like ENDES 2024). Use calibrated mode for reliable results.")
            )
          ),
          
          conditionalPanel(
            condition = "input.markov_rates_mode == 'calibrated'",
            div(class = "alert alert-info", style = "font-size: 0.85em;",
                icon("info-circle"),
                HTML("<b>Calibrated mode:</b> Uses t1=0.20, t2=0.48, t3=0.89, t4=0.62. 
                     These are design assumptions. You should report a sensitivity analysis.")
            )
          ),
          
          conditionalPanel(
            condition = "input.markov_rates_mode == 'manual'",
            textInput("tasas_avance_input", "Enter 4 rates (comma separated):",
                      placeholder = "e.g., 0.20, 0.48, 0.89, 0.62"),
            helpText("Provide rates for E1→E2, E2→E3, E3→E4, E4→E5.")
          ),
          
          hr(),
          h4("Bootstrap threshold selection (optional)"),
          fluidRow(
            column(6, numericInput("boot_n_sim", "Replicates:", value = 20, min = 10, max = 50, step = 5)),
            column(6, actionButton("run_bootstrap", "Run", class = "btn-info w-100"))
          ),
          conditionalPanel(
            condition = "input.run_bootstrap > 0",
            div(style = "padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
                uiOutput("bootstrap_results_ui"),
                plotOutput("bootstrap_plot", height = "250px"))
          ),
          
          hr(),
          h5("MPCS Index Weights"),
          fluidRow(
            column(4, numericInput("w1", "Graph (w1)", value = 0.35, min = 0, max = 1, step = 0.05)),
            column(4, numericInput("w2", "Markov (w2)", value = 0.40, min = 0, max = 1, step = 0.05)),
            column(4, numericInput("w3", "Games (w3)", value = 0.25, min = 0, max = 1, step = 0.05))
          ),
          tags$small("Weights must sum to 1. Current: ", textOutput("suma_ponderadores", inline = TRUE)),
          hr(),
          sliderInput("R_factor", "Resource factor (R)", min = 0, max = 1, value = 0.65, step = 0.05),
          hr(),
          h4("Contextual indicator (alpha)"),
          fluidRow(
            column(6, radioButtons("alpha_mode", "Method:", choices = c("Select vars" = "vars", "Manual" = "manual", "Auto" = "auto"), selected = "auto")),
            column(6, 
              conditionalPanel(condition = "input.alpha_mode == 'vars'", uiOutput("alpha_vars_ui")),
              conditionalPanel(condition = "input.alpha_mode == 'manual'", numericInput("alpha_manual", "Value (0-1):", 0.60, 0, 1, 0.01))
            )
          )
        )
      )
    ),
    fluidRow(
      column(12, wellPanel(actionButton("run_mpcs", "Run MPCS", class = "btn-success btn-lg w-100", icon = icon("play")), uiOutput("validation_msg")))
    )
  ),
  
  # ============================================================================
  # Pestaña 3: Resultados
  # ============================================================================
  nav_panel(
    "3. Results",
    fluidRow(column(width = 12, wellPanel(h4("Results table"), DTOutput("results_table")))),
    fluidRow(
      column(6, wellPanel(h5("Behavioral graph"), plotOutput("plot_graph", height = "400px"))),
      column(6, wellPanel(h5("State distribution"), plotOutput("plot_states", height = "400px")))
    ),
    fluidRow(
      column(6, wellPanel(h5("Markov transitions"), plotOutput("plot_arbol_markov", height = "400px"))),
      column(6, wellPanel(h5("Game theory"), plotOutput("plot_juego_evolutivo", height = "400px")))
    ),
    fluidRow(
      column(6, wellPanel(h5("Markov trajectories"), plotOutput("plot_markov", height = "400px"))),
      column(6, wellPanel(h5("MPCS ranking"), plotOutput("plot_ranking", height = "400px")))
    ),
    fluidRow(
      column(12, wellPanel(
        h4("Sensitivity Analysis"),
        fluidRow(
          column(3, numericInput("sens_n_sim", "Simulations:", 10000, 100, 50000, 100), actionButton("run_sensitivity", "Run", class = "btn-primary w-100")),
          column(9, conditionalPanel("input.run_sensitivity > 0", plotOutput("sens_plot", height = "300px"), verbatimTextOutput("sens_summary")))
        )
      ))
    ),
    fluidRow(column(12, wellPanel(h4("Automatic interpretation"), uiOutput("interpretation_text"))))
  ),
  
  # ============================================================================
  # Pestaña 4: Reporte
  # ============================================================================
  nav_panel(
    "4. Report",
    fluidRow(
      column(6, wellPanel(h4("Download data"), downloadButton("download_csv", "Download CSV", class = "btn-primary w-100"))),
      column(6, wellPanel(h4("Generate report"), downloadButton("download_report", "Download HTML", class = "btn-danger w-100")))
    )
  )
)

# ==============================================================================
# SERVER — LOGICA DE LA APLICACION
# ==============================================================================

server <- function(input, output, session) {
  
  rv <- reactiveValues(data = NULL, results = NULL, results_df = NULL, plots = NULL, sens_resultado = NULL, bootstrap_resultado = NULL)
  
  observeEvent(input$file, {
    req(input$file); ext <- tools::file_ext(input$file$datapath)
    datos <- tryCatch({
      switch(ext, csv = read.csv(input$file$datapath), xlsx = readxl::read_excel(input$file$datapath), dta = haven::read_dta(input$file$datapath))
    }, error = function(e) { showNotification(paste("Error:", e$message), type = "error"); NULL })
    
    if (!is.null(datos)) {
      rv$data <- datos
      updateSelectInput(session, "markov_var", choices = names(datos), selected = if ("Estado_Markov" %in% names(datos)) "Estado_Markov" else names(datos)[1])
    }
  })
  
  observeEvent(input$load_demo, {
    if (file.exists("data/demo_data3.csv")) rv$data <- read.csv("data/demo_data3.csv") else rv$data <- generar_demo_data()
    updateSelectInput(session, "markov_var", choices = names(rv$data), selected = "Estado_Markov")
  })
  
  output$data_preview <- renderDT({ req(rv$data); datatable(head(rv$data, 5)) })
  output$data_stats <- renderPrint({ req(rv$data); cat("Rows:", nrow(rv$data), "Cols:", ncol(rv$data)) })
  
  output$graph_vars_ui <- renderUI({ req(rv$data); vars_num <- names(rv$data)[sapply(rv$data, is.numeric)]; checkboxGroupInput("graph_vars", "Variables:", choices = vars_num, selected = vars_num[1:min(5, length(vars_num))]) })
  output$group_var_ui <- renderUI({ req(rv$data); selectInput("group_var", "Group by:", choices = c("None", names(rv$data))) })
  output$markov_var_ui <- renderUI({ req(rv$data); selectInput("markov_var", "State column:", choices = names(rv$data)) })
  output$alpha_vars_ui <- renderUI({ req(rv$data); vars_num <- names(rv$data)[sapply(rv$data, is.numeric)]; selectInput("alpha_vars", "Vars for alpha:", choices = vars_num, multiple = TRUE, selected = vars_num[1]) })
  
  output$suma_ponderadores <- renderText({ total <- sum(input$w1, input$w2, input$w3); if (abs(total - 1) < 0.01) paste0(total, " OK") else paste0(total, " (must sum to 1)") })
  
  observeEvent(input$run_mpcs, {
    req(rv$data, input$graph_vars, input$markov_var)
    if (length(input$graph_vars) < 5) { output$validation_msg <- renderUI({ tags$div(class = "alert alert-danger", "Select at least 5 variables.") }); return() }
    output$validation_msg <- renderUI({ NULL })
    
    withProgress(message = 'Running...', value = 0.5, {
      data_analysis <- rv$data
      data_analysis$Group <- if (input$group_var == "None") "Global" else data_analysis[[input$group_var]]
      grupos <- unique(data_analysis$Group)
      results_list <- list()
      
      for (g in grupos) {
        sub <- data_analysis[data_analysis$Group == g, ]
        if (nrow(sub) < 30) next
        
        graph_res <- calcular_grafo(sub, input$graph_vars, input$threshold)
        orden_estados <- sort(unique(sub[[input$markov_var]]))
        
        # LÓGICA DE TASAS CORREGIDA
        tasas_avance <- NULL
        if (input$markov_rates_mode == "calibrated") {
          tasas_avance <- c(0.48, 0.89, 0.62)
        } else if (input$markov_rates_mode == "manual") {
          if (!is.null(input$tasas_avance_input) && input$tasas_avance_input != "") {
            t_num <- suppressWarnings(as.numeric(strsplit(trimws(input$tasas_avance_input), ",")[[1]]))
            if (length(t_num) == 4) {
              tasas_avance <- t_num
            } else {
              showNotification("Manual entry requires 4 rates. Using CALIBRATED values.", type = "warning")
              tasas_avance <- c(0.48, 0.89, 0.62) # Fallback seguro
            }
          } else {
            showNotification("Manual entry empty. Using CALIBRATED values.", type = "warning")
            tasas_avance <- c(0.48, 0.89, 0.62) # Fallback seguro
          }
        }
        
        markov_res <- calcular_markov(sub[[input$markov_var]], orden_estados, tasas_avance)
        alpha_val <- if(input$alpha_mode == "manual") input$alpha_manual else 0.60
        games_res <- calcular_juegos(markov_res$mat, input$R_factor, alpha_val)
        idx_res <- calcular_indice(graph_res$score, markov_res$score, games_res$score, input$w1, input$w2, input$w3, input$R_factor)
        
        results_list[[g]] <- list(grupo = g, n = nrow(sub), nodo = graph_res$optimal_node, I_M = markov_res$score, I_G = games_res$score, I_MPCS = idx_res$I_MPCS, k = idx_res$k, tipo = idx_res$nudge_type, graph = graph_res$graph, mat = markov_res$mat, sim = markov_res$sim_base)
      }
      
      rv$results <- results_list
      rv$results_df <- do.call(rbind, lapply(results_list, function(x) data.frame(Group = x$grupo, Node = x$nodo, I_Markov = round(x$I_M, 3), I_Games = round(x$I_G, 3), I_MPCS = round(x$I_MPCS, 3), Nudge = x$tipo)))
    })
  })
  
  output$results_table <- renderDT({ req(rv$results_df); datatable(rv$results_df) })
  output$plot_graph <- renderPlot({ req(rv$results); plot(rv$results[[1]]$graph) })
  output$plot_markov <- renderPlot({ req(rv$results); plot(rv$results[[1]]$sim[, 5], type = "l", ylim = c(0,1)) })
  output$download_csv <- downloadHandler(filename = function() "results.csv", content = function(file) write.csv(rv$results_df, file))
}

shinyApp(ui = ui, server = server)
