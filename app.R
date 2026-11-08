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
# GENERAR DATOS DE DEMOSTRACION INICIALES (SOLO ENDES)
# ==============================================================================
if (!file.exists("data/demo_data3.csv")) {
  dir.create("data", showWarnings = FALSE)
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
          
          # ==================================================================
          # CAMBIO 1: LÓGICA DE TASAS (Estimar vs Calibrado vs Manual)
          # ==================================================================
          hr(),
          h4("Markov Transition Rates"),
          radioButtons("markov_rates_mode", "Method to define transition probabilities:",
                       choices = c(
                         "Estimate from data (Equation X)" = "estimate",
                         "Use calibrated values (Manuscript)" = "calibrated",
                         "Enter manually (4 rates)" = "manual"
                       ),
                       selected = "estimate"),
          
          conditionalPanel(
            condition = "input.markov_rates_mode == 'estimate'",
            div(class = "alert alert-info", style = "font-size: 0.85em;",
                icon("info-circle"),
                HTML("<b>Estimation mode:</b> t1 (E1→E2) is fixed at 0.20 (diagnostic event). 
                     t2, t3, t4 are estimated from the observed state proportions using Eq. X.")
            )
          ),
          
          conditionalPanel(
            condition = "input.markov_rates_mode == 'calibrated'",
            div(class = "alert alert-success", style = "font-size: 0.85em;",
                icon("check-circle"),
                HTML("<b>Calibrated mode:</b> Uses t1=0.20, t2=0.48, t3=0.89, t4=0.62. 
                     This reproduces the exact values from the main manuscript.")
            )
          ),
          
          conditionalPanel(
            condition = "input.markov_rates_mode == 'manual'",
            textInput("tasas_avance_input", "Enter 4 rates (comma separated):",
                      placeholder = "e.g., 0.20, 0.48, 0.89, 0.62"),
            helpText("Provide rates for E1→E2, E2→E3, E3→E4, E4→E5.")
          ),
          
          # Modulo de Bootstrap para selección de umbral
          hr(),
          h4("Bootstrap threshold selection (optional)"),
          helpText("Evaluates the stability of the optimal node across different correlation thresholds."),
          fluidRow(
            column(6,
              numericInput("boot_n_sim", "Number of bootstrap replicates:",
                           value = 20, min = 10, max = 50, step = 5)
            ),
            column(6,
              actionButton("run_bootstrap", "Run Bootstrap",
                           class = "btn-info w-100", icon = icon("random"))
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
            column(6,
              radioButtons("alpha_mode", "Method for alpha:",
                           choices = c("Select variables" = "vars", "Enter manual value" = "manual", "Automatic" = "auto"),
                           selected = "auto")
            ),
            column(6,
              conditionalPanel(condition = "input.alpha_mode == 'vars'", uiOutput("alpha_vars_ui")),
              conditionalPanel(condition = "input.alpha_mode == 'manual'", numericInput("alpha_manual", "Alpha value (0-1):", value = 0.60, min = 0, max = 1, step = 0.01))
            )
          )
        )
      )
    ),
    fluidRow(
      column(
        width = 12,
        wellPanel(
          actionButton("run_mpcs", "Run MPCS", class = "btn-success btn-lg w-100", icon = icon("play")),
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
      column(width = 12, wellPanel(h4("Results table"), DTOutput("results_table")))
    ),
    fluidRow(
      column(width = 6, wellPanel(h5("Behavioral graph"), plotOutput("plot_graph", height = "500px"))),
      column(width = 6, wellPanel(h5("State distribution by group"), plotOutput("plot_states", height = "500px")))
    ),
    fluidRow(
      column(width = 6, wellPanel(h5("Markov transitions"), plotOutput("plot_arbol_markov", height = "450px"))),
      column(width = 6, wellPanel(h5("Game theory dynamics"), plotOutput("plot_juego_evolutivo", height = "450px")))
    ),
    fluidRow(
      column(width = 6, wellPanel(h5("Markov trajectories"), plotOutput("plot_markov", height = "400px"))),
      column(width = 6, wellPanel(h5("I_MPCS ranking"), plotOutput("plot_ranking", height = "400px")))
    ),
    fluidRow(
      column(
        width = 12,
        wellPanel(
          h4("Sensitivity Analysis"),
          fluidRow(
            column(width = 3,
              numericInput("sens_n_sim", "Number of simulations:", value = 10000, min = 100, max = 50000, step = 100),
              actionButton("run_sensitivity", "Run analysis", icon = icon("chart-line"), class = "btn-primary", style = "width: 100%;")
            ),
            column(width = 9,
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
    fluidRow(
      column(width = 12, wellPanel(h4("Automatic interpretation"), uiOutput("interpretation_text")))
    )
  ),
  
  # ============================================================================
  # Pestaña 4: Reporte y Descarga
  # ============================================================================
  nav_panel(
    "4. Report & Download",
    fluidRow(
      column(width = 6, wellPanel(h4("Download data"), downloadButton("download_csv", "Download CSV", class = "btn-primary w-100"))),
      column(width = 6, wellPanel(h4("Generate report"), downloadButton("download_report", "Download HTML Report", class = "btn-danger w-100")))
    )
  )
)

# ==============================================================================
# SERVER — LOGICA DE LA APLICACION
# ==============================================================================

server <- function(input, output, session) {
  
  rv <- reactiveValues(
    data = NULL, results = NULL, results_df = NULL, plots = NULL,
    sens_resultado = NULL, bootstrap_resultado = NULL, dominio = "salud"
  )
  
  observeEvent(input$file, {
    req(input$file)
    ext <- tools::file_ext(input$file$datapath)
    datos <- tryCatch({
      switch(ext, csv = read.csv(input$file$datapath), xlsx = readxl::read_excel(input$file$datapath),
             xls = readxl::read_excel(input$file$datapath), dta = haven::read_dta(input$file$datapath), NULL)
    }, error = function(e) { showNotification(paste("Error:", e$message), type = "error"); NULL })
    
    if (is.null(datos)) return()
    rv$data <- datos
    vars <- names(datos); vars_num <- names(datos)[sapply(datos, is.numeric)]
    updateSelectInput(session, "graph_vars", choices = vars_num, selected = vars_num[1:min(5, length(vars_num))])
    updateSelectInput(session, "markov_var", choices = vars, selected = if ("Estado_Markov" %in% vars) "Estado_Markov" else vars[length(vars)])
    updateSelectInput(session, "group_var", choices = c("None (Global)", vars), selected = if ("Region" %in% vars) "Region" else "None (Global)")
    updateSelectInput(session, "alpha_vars", choices = vars_num, selected = vars_num[1:min(2, length(vars_num))])
  })
  
  observeEvent(input$load_demo, {
    if (file.exists("data/demo_data3.csv")) {
      rv$data <- read.csv("data/demo_data3.csv")
    } else {
      rv$data <- generar_demo_data(1000, 123)
    }
    vars <- names(rv$data); vars_num <- names(rv$data)[sapply(rv$data, is.numeric)]
    updateSelectInput(session, "graph_vars", choices = vars_num, selected = vars_num[1:min(5, length(vars_num))])
    updateSelectInput(session, "markov_var", choices = vars, selected = "Estado_Markov")
    updateSelectInput(session, "group_var", choices = c("None (Global)", vars), selected = "Region")
    updateSelectInput(session, "alpha_vars", choices = vars_num, selected = c("Acceso_salud", "Tiene_seguro"))
  })
  
  output$data_preview <- renderDT({ req(rv$data); datatable(head(rv$data, 10), options = list(scrollX = TRUE, dom = 't')) })
  output$data_stats <- renderPrint({ req(rv$data); cat("Observations:", nrow(rv$data), "\nVariables:", ncol(rv$data), "\n") })
  
  output$graph_vars_ui <- renderUI({
    req(rv$data); vars_num <- names(rv$data)[sapply(rv$data, is.numeric)]
    checkboxGroupInput("graph_vars", "System variables:", choices = vars_num, selected = vars_num[1:min(5, length(vars_num))])
  })
  output$alpha_vars_ui <- renderUI({
    req(rv$data); vars_num <- names(rv$data)[sapply(rv$data, is.numeric)]
    selectInput("alpha_vars", "Variables for alpha:", choices = vars_num, multiple = TRUE, selected = vars_num[1:min(2, length(vars_num))])
  })
  output$markov_var_ui <- renderUI({
    req(rv$data); vars <- names(rv$data)
    default <- if ("Estado_Markov" %in% vars) "Estado_Markov" else vars[1]
    selectInput("markov_var", "Markov state variable:", choices = vars, selected = default)
  })
  output$group_var_ui <- renderUI({
    req(rv$data); vars <- c("None (Global)", names(rv$data))
    default <- if ("Region" %in% names(rv$data)) "Region" else "None (Global)"
    selectInput("group_var", "Grouping variable:", choices = vars, selected = default)
  })
  output$suma_ponderadores <- renderText({
    total <- sum(input$w1, input$w2, input$w3)
    if (abs(total - 1) < 0.01) paste0(total, " OK") else paste0(total, " (must sum to 1)")
  })
  
  detectar_orden_estados <- function(estados_unicos) {
    if (all(grepl("^E[0-9]+$", estados_unicos))) return(estados_unicos[order(as.numeric(gsub("E", "", estados_unicos)))])
    if (all(grepl("^[0-9]+$", estados_unicos))) return(as.character(sort(as.numeric(estados_unicos))))
    return(sort(estados_unicos))
  }
  
  calcular_alpha <- function(sub) {
    alpha <- 0.60
    if (input$alpha_mode == "manual") { alpha <- input$alpha_manual }
    else if (input$alpha_mode == "vars" && !is.null(input$alpha_vars)) {
      vars_alpha <- intersect(input$alpha_vars, names(sub))
      if (length(vars_alpha) > 0) {
        alpha_vals <- sapply(vars_alpha, function(v) {
          vals <- sub[[v]]; rng <- range(vals, na.rm = TRUE)
          if (rng[2] > rng[1]) (vals - rng[1]) / (rng[2] - rng[1]) else rep(0.5, length(vals))
        })
        alpha <- mean(alpha_vals, na.rm = TRUE)
      }
    } else {
      if ("Acceso_salud" %in% names(sub)) alpha <- mean(sub$Acceso_salud, na.rm = TRUE)
    }
    if (is.na(alpha) || is.nan(alpha)) alpha <- 0.60
    return(max(0, min(1, alpha)))
  }
  
  observeEvent(input$run_bootstrap, {
    req(rv$data)
    if (is.null(input$graph_vars) || length(input$graph_vars) < 5) {
      showNotification("Select at least 5 variables.", type = "error"); return()
    }
    withProgress(message = 'Running bootstrap...', value = 0.1, {
      resultado_boot <- seleccionar_umbral_bootstrap(
        datos = rv$data, variables = input$graph_vars,
        umbrales = c(0.05, 0.07, 0.10, 0.12, 0.15, 0.20),
        n_boot = input$boot_n_sim, seed = 123, criterio_jaccard = 0.70
      )
      rv$bootstrap_resultado <- resultado_boot
      incProgress(1.0, detail = "Completed!")
    })
  })
  
  output$bootstrap_results_ui <- renderUI({
    req(rv$bootstrap_resultado)
    res <- rv$bootstrap_resultado
    HTML(paste0("<div class='alert alert-info'><b>Optimal threshold:</b> ", round(res$umbral_optimo, 2), 
                "<br><b>Original node:</b> ", res$nodo_original, 
                "<br><b>Stability:</b> ", round(res$pct_nodo_mas_frecuente, 1), "%</div>"))
  })
  output$bootstrap_plot <- renderPlot({ req(rv$bootstrap_resultado); NULL })
  
  # ==========================================================================
  # EJECUTAR MPCS (CON NUEVA LÓGICA DE TASAS)
  # ==========================================================================
  observeEvent(input$run_mpcs, {
    req(rv$data)
    
    if (is.null(input$graph_vars) || length(input$graph_vars) < 5) {
      output$validation_msg <- renderUI({ tags$div(class = "alert alert-danger", "Select at least 5 variables.") }); return()
    }
    output$validation_msg <- renderUI({ NULL })
    
    withProgress(message = 'Running MPCS...', value = 0, {
      data_analysis <- rv$data
      if (input$group_var == "None (Global)") {
        data_analysis$Group <- "Global"; grupos <- "Global"
      } else {
        data_analysis$Group <- data_analysis[[input$group_var]]
        grupos <- unique(data_analysis$Group); grupos <- grupos[!is.na(grupos) & grupos != ""]
      }
      
      results_list <- list()
      for (i in seq_along(grupos)) {
        g <- grupos[i]; incProgress(1 / length(grupos), detail = paste("Processing:", g))
        sub <- data_analysis[data_analysis$Group == g, ]
        if (nrow(sub) < 30) next
        
        graph_data <- sub[, input$graph_vars, drop = FALSE]
        graph_res <- calcular_grafo(graph_data, input$graph_vars, input$threshold)
        if (is.null(graph_res$graph) || vcount(graph_res$graph) < 2) {
          graph_res$score <- 0.5; graph_res$optimal_node <- "No node"; graph_res$graph <- NULL
        }
        
        alpha_grupo <- calcular_alpha(sub)
        estados_unicos <- unique(sub[[input$markov_var]])
        estados_unicos <- estados_unicos[!is.na(estados_unicos) & estados_unicos != ""]
        orden_estados <- detectar_orden_estados(estados_unicos)
        
        # ======================================================================
        # CAMBIO 2: SELECCIÓN DE TASAS SEGÚN EL USUARIO
        # ======================================================================
        tasas_avance <- NULL
        if (input$markov_rates_mode == "calibrated") {
          tasas_avance <- c(0.48, 0.89, 0.62) # La función internamente fija t1=0.20
        } else if (input$markov_rates_mode == "manual") {
          if (!is.null(input$tasas_avance_input) && input$tasas_avance_input != "") {
            tasas_num <- suppressWarnings(as.numeric(strsplit(trimws(input$tasas_avance_input), ",")[[1]]))
            if (length(tasas_num) == 4) {
              tasas_avance <- tasas_num
            } else {
              showNotification("Manual entry requires 4 rates. Using estimation instead.", type = "warning")
            }
          }
        }
        # Si es "estimate", tasas_avance se queda NULL y calcular_markov lo estima.
        
        markov_res <- calcular_markov(
          sub[[input$markov_var]], 
          orden_estados = orden_estados,
          tasas_avance = tasas_avance,
          horizonte = 10
        )
        
        games_res <- calcular_juegos(markov_res$mat, input$R_factor, alpha = alpha_grupo)
        index_res <- calcular_indice(graph_res$score, markov_res$score, games_res$score, input$w1, input$w2, input$w3, input$R_factor)
        
        results_list[[as.character(g)]] <- list(
          grupo = g, n = nrow(sub), nodo_optimo = graph_res$optimal_node,
          I_grafo = graph_res$score, I_markov = markov_res$score, I_juegos = games_res$score,
          alpha = alpha_grupo, I_MPCS = index_res$I_MPCS, k = index_res$k, tipo = index_res$nudge_type,
          graph = graph_res$graph, graph_data = graph_data, markov_mat = markov_res$mat,
          sim_base = markov_res$sim_base, dist_actual = markov_res$dist_actual,
          T_base = markov_res$T_base, estados = colnames(markov_res$mat),
          tasas_usadas = markov_res$tasas_usadas # Guardamos las tasas usadas
        )
      }
      
      if (length(results_list) == 0) return()
      
      for (g in names(results_list)) {
        r <- results_list[[g]]
        if (!is.null(r$markov_mat) && !is.null(r$sim_base)) {
          P_n <- aplicar_nudge(r$markov_mat, r$k)
          sim_nudge <- r$sim_base
          for (j in 2:nrow(sim_nudge)) sim_nudge[j, ] <- sim_nudge[j-1, ] %*% P_n
          results_list[[g]]$sim_nudge <- sim_nudge
        }
      }
      
      rv$results <- results_list
      rv$results_df <- do.call(rbind, lapply(results_list, function(r) {
        data.frame(Grupo = r$grupo, n = r$n, Alpha = round(r$alpha, 3),
                   I_Markov_H10 = round(r$I_markov, 4), I_Games = round(r$I_juegos, 4),
                   I_MPCS = round(r$I_MPCS, 4), Nodo_Optimo = r$nodo_optimo,
                   k = round(r$k, 4), Tipo_Nudge = r$tipo, stringsAsFactors = FALSE)
      }))
      
      rv$plots <- generate_plots(data_analysis, input$graph_vars, input$markov_var, rv$results_df, results_list, input$threshold, rv$dominio)
    })
  })
  
  # ==========================================================================
  # FUNCIONES PARA GRAFICOS (Sin cambios)
  # ==========================================================================
  generate_plots <- function(data, graph_vars, markov_var, results, results_list, threshold, dominio = "salud") {
    p_graph <- p_states <- p_markov <- p_rank <- p_arbol <- p_juego <- NULL
    if (is.null(data) || nrow(data) == 0 || is.null(results) || nrow(results) == 0) return(list())
    
    first_group <- results$Grupo[1]
    r <- results_list[[as.character(first_group)]]
    
    if (!is.null(r$graph) && vcount(r$graph) > 0) {
      V(r$graph)$color <- ifelse(V(r$graph)$name == r$nodo_optimo, "#C0392B", "#F0DFC0")
      V(r$graph)$size <- ifelse(V(r$graph)$name == r$nodo_optimo, 25, 15)
      p_graph <- function() {
        plot(r$graph, layout = layout_with_fr(r$graph), vertex.label.cex = 0.8, vertex.label.color = "black", edge.color = "gray60", edge.width = 1.5, main = paste("Behavioral Graph -", first_group))
        legend("topright", legend = c("Optimal node", "Other nodes"), fill = c("#C0392B", "#F0DFC0"), cex = 0.8, bty = "n")
      }
    }
    
    if (!is.null(markov_var) && markov_var %in% names(data)) {
      p_states <- function() {
        if (!"Group" %in% names(data)) data$Group <- "Global"
        plot_data <- data[!is.na(data$Group) & !is.na(data[[markov_var]]), ]
        ggplot(plot_data, aes(x = .data[["Group"]], fill = .data[[markov_var]])) + geom_bar(position = "fill") + scale_fill_brewer(palette = "Set2") + theme_minimal() + labs(x = "Group", y = "Proportion", fill = "State") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
    }
    
    if (!is.null(results) && nrow(results) > 0) {
      p_rank <- function() {
        ggplot(results, aes(x = reorder(Grupo, I_MPCS), y = I_MPCS, fill = Tipo_Nudge)) + geom_col(width = 0.7) + coord_flip() + theme_minimal() + labs(x = "Group", y = "MPCS Index", fill = "Nudge Type") + scale_fill_manual(values = c("Informational" = "#74B3CE", "Structural" = "#2E86AB", "Normative" = "#E84855", "Systemic multi-nudge" = "#1A3A5C")) + geom_text(aes(label = round(I_MPCS, 3)), hjust = -0.2, size = 3.5) + theme(legend.position = "bottom")
      }
    }
    
    p_markov <- function() {
      if (!is.null(r$sim_base) && nrow(r$sim_base) > 0) {
        graficar_trayectorias_markov(sim_base = r$sim_base, sim_nudge = r$sim_nudge, estados = colnames(r$sim_base), titulo = paste("Markov Trajectories -", first_group))
      }
    }
    
    p_arbol <- function() {
      if (!is.null(r$markov_mat)) graficar_arbol_markov(P = r$markov_mat, estados = colnames(r$markov_mat), umbral_prob = 0.05, titulo = paste("Markov Transitions -", first_group))
    }
    
    p_juego <- function() {
      if (!is.null(r$alpha)) graficar_juego_evolutivo(alpha = r$alpha, p_star = r$I_juegos, titulo = paste("Replicator Dynamics -", first_group))
    }
    
    return(list(graph = p_graph, states = p_states, markov = p_markov, rank = p_rank, arbol = p_arbol, juego = p_juego))
  }
  
  output$results_table <- renderDT({ req(rv$results_df); datatable(rv$results_df, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE) %>% formatRound(columns = c("I_Markov_H10", "I_Games", "I_MPCS", "k"), digits = 4) })
  output$plot_graph <- renderPlot({ req(rv$plots); if (!is.null(rv$plots$graph)) rv$plots$graph() else plot(0, type = "n", axes = FALSE, xlab = "", ylab = "") })
  output$plot_states <- renderPlot({ req(rv$plots); if (!is.null(rv$plots$states)) rv$plots$states() })
  output$plot_markov <- renderPlot({ req(rv$plots); if (!is.null(rv$plots$markov)) rv$plots$markov() })
  output$plot_ranking <- renderPlot({ req(rv$plots); if (!is.null(rv$plots$rank)) rv$plots$rank() })
  output$plot_arbol_markov <- renderPlot({ req(rv$plots); if (!is.null(rv$plots$arbol)) rv$plots$arbol() })
  output$plot_juego_evolutivo <- renderPlot({ req(rv$plots); if (!is.null(rv$plots$juego)) rv$plots$juego() })
  
  # ==========================================================================
  # ANALISIS DE SENSIBILIDAD E INTERPRETACIÓN
  # ==========================================================================
  observeEvent(input$run_sensitivity, {
    req(rv$results_df)
    top_group <- rv$results_df$Grupo[which.max(rv$results_df$I_MPCS)]
    r <- rv$results[[as.character(top_group)]]
    if (is.null(r)) return()
    
    withProgress(message = 'Analyzing...', value = 0, {
      rv$sens_resultado <- analisis_sensibilidad(r$I_grafo, r$I_markov, r$I_juegos, n_sim = input$sens_n_sim, R = input$R_factor, seed = 123)
      incProgress(1, detail = "Completed")
    })
  })
  
  output$sens_plot <- renderPlot({
    req(rv$sens_resultado)
    df <- data.frame(I_MPCS = rv$sens_resultado$I_MPCS_sim, Tipo = rv$sens_resultado$tipo_sim)
    ggplot(df, aes(x = I_MPCS, fill = Tipo)) + geom_histogram(bins = 50, color = "white") + geom_vline(xintercept = rv$sens_resultado$I_MPCS_default, color = "red", linetype = "dashed") + theme_minimal()
  })
  
  output$sens_summary <- renderPrint({ req(rv$sens_resultado); cat("Mean:", rv$sens_resultado$stats$media, "\nAgreement:", rv$sens_resultado$pct_coincidencia, "%") })
  
  output$interpretation_text <- renderUI({
    req(rv$results_df)
    top_group <- rv$results_df[which.max(rv$results_df$I_MPCS), ]
    HTML(paste0("<p><b>Highest priority group:</b> ", top_group$Grupo, " (I_MPCS = ", round(top_group$I_MPCS, 4), ")</p>",
                "<p><b>Recommended nudge:</b> ", top_group$Tipo_Nudge, " with intensity k = ", round(top_group$k, 4), "</p>"))
  })
  
  output$download_csv <- downloadHandler(filename = function() paste0("MPCS_Results_", Sys.Date(), ".csv"), content = function(file) write.csv(rv$results_df, file, row.names = FALSE))
  output$download_report <- downloadHandler(filename = function() paste0("MPCS_Report_", Sys.Date(), ".html"), content = function(file) {
    # (Lógica de reporte omitida por brevedad, pero debe mantenerse igual a tu versión original)
    writeLines("<html><body><h1>MPCS Report</h1><p>See CSV for data.</p></body></html>", file)
  })
}

ui <- tagList(ui, tags$footer(class = "mpcs-footer", style = "text-align: center; padding: 20px; margin-top: 50px; background: #f8f9fa;", 
                              "MPCS Calculator v1.0 | ", tags$a(href = "https://github.com/Izela-meth/MPCS_APP", target = "_blank", "GitHub")))

shinyApp(ui = ui, server = server)
