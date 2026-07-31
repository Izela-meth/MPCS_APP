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
    "1. Carga de Datos",
    fluidRow(
      column(
        width = 4,
        wellPanel(
          h4("Cargar archivo"),
          fileInput("file", "Selecciona un archivo",
                    accept = c(".csv", ".xlsx", ".xls", ".dta"),
                    buttonLabel = "Examinar",
                    placeholder = "Ningun archivo seleccionado"),
          tags$small("Formatos soportados: CSV, Excel (.xlsx, .xls), Stata (.dta)"),
          hr(),
          h4("O usar datos de demostracion"),
          p("Carga datos simulados de salud (ENDES) para probar la aplicacion."),
          actionButton("load_demo", "Cargar Datos ENDES (Salud)", 
                       class = "btn-primary w-100",
                       icon = icon("heart"))
        )
      ),
      column(
        width = 8,
        wellPanel(
          h5("Vista previa (primeras 10 filas)"),
          DTOutput("data_preview"),
          hr(),
          h5("Estadisticas basicas"),
          verbatimTextOutput("data_stats")
        )
      )
    )
  ),
  
  # ============================================================================
  # Pestaña 2: Configuracion
  # ============================================================================
  nav_panel(
    "2. Configuracion",
    fluidRow(
      column(
        width = 6,
        wellPanel(
          h4("Variables del grafo"),
          helpText("Selecciona minimo 5 variables numericas para construir el grafo conductual."),
          uiOutput("graph_vars_ui"),
          hr(),
          h4("Variable de agrupacion (opcional)"),
          helpText("Si seleccionas una variable, el analisis se realizara por separado para cada grupo."),
          uiOutput("group_var_ui")
        )
      ),
      column(
        width = 6,
        wellPanel(
          h4("Configuracion de Markov"),
          helpText("Selecciona la columna que contiene los estados conductuales."),
          uiOutput("markov_var_ui"),
          hr(),
          h4("Parametros ajustables"),
          sliderInput("threshold", "Umbral de correlacion",
                      min = 0.05, max = 0.30, value = 0.10, step = 0.01,
                      post = tags$span("  (|r| > valor)")),
          hr(),
          h5("Ponderadores del Indice MPCS"),
          fluidRow(
            column(4, numericInput("w1", "Grafo (w1)", 
                                   value = 0.35, min = 0, max = 1, step = 0.05)),
            column(4, numericInput("w2", "Markov (w2)", 
                                   value = 0.40, min = 0, max = 1, step = 0.05)),
            column(4, numericInput("w3", "Juegos (w3)", 
                                   value = 0.25, min = 0, max = 1, step = 0.05))
          ),
          tags$small("Los ponderadores deben sumar 1. Actual: ", 
                     textOutput("suma_ponderadores", inline = TRUE)),
          hr(),
          sliderInput("R_factor", "Factor de recursos (R)",
                      min = 0, max = 1, value = 0.65, step = 0.05,
                      post = tags$span("  (mayor = mas recursos disponibles)")),
          hr(),
          h4("Indicador contextual (alpha)"),
          helpText("alpha representa la favorabilidad del contexto para la adopcion de la conducta."),
          fluidRow(
            column(6,
              radioButtons("alpha_mode", "Metodo para alpha:",
                           choices = c("Seleccionar variables" = "vars",
                                      "Ingresar valor manual" = "manual",
                                      "Automatico" = "auto"),
                           selected = "auto")
            ),
            column(6,
              conditionalPanel(
                condition = "input.alpha_mode == 'vars'",
                uiOutput("alpha_vars_ui")
              ),
              conditionalPanel(
                condition = "input.alpha_mode == 'manual'",
                numericInput("alpha_manual", "Valor de alpha (0-1):",
                             value = 0.60, min = 0, max = 1, step = 0.01)
              ),
              conditionalPanel(
                condition = "input.alpha_mode == 'auto'",
                helpText("La app usara el promedio de variables relacionadas con acceso/recurso.")
              )
            )
          ),
          div(
            class = "alert alert-info mt-2",
            style = "font-size: 0.85em;",
            icon("circle-info"),
            HTML(
              "<b>Nota metodologica.</b><br><br>
              <b>I_Markov:</b> Se calcula como la probabilidad acumulada en los estados 
              avanzados en un horizonte fijo de 10 periodos.<br><br>
              <b>alpha (contexto):</b> Si seleccionas variables, se normalizan a 0-1 y se promedian. 
              Si ingresas manual, usa ese valor. En modo automatico, busca variables comunes 
              como 'Acceso_salud' o usa la primera variable numerica."
            )
          )
        )
      )
    ),
    fluidRow(
      column(
        width = 12,
        wellPanel(
          actionButton("run_mpcs", "Ejecutar MPCS", 
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
    "3. Resultados",
    fluidRow(
      column(
        width = 12,
        wellPanel(
          h4("Tabla de resultados"),
          DTOutput("results_table")
        )
      )
    ),
    fluidRow(
      column(
        width = 6,
        wellPanel(
          h5("Grafo conductual"),
          plotOutput("plot_graph", height = "500px")
        )
      ),
      column(
        width = 6,
        wellPanel(
          h5("Distribucion de estados por grupo"),
          plotOutput("plot_states", height = "500px")
        )
      )
    ),
    fluidRow(
      column(
        width = 6,
        wellPanel(
          h5("Transiciones de Markov"),
          plotOutput("plot_arbol_markov", height = "450px")
        )
      ),
      column(
        width = 6,
        wellPanel(
          h5("Dinamica de Teoria de Juegos"),
          plotOutput("plot_juego_evolutivo", height = "450px")
        )
      )
    ),
    fluidRow(
      column(
        width = 6,
        wellPanel(
          h5("Trayectorias de Markov"),
          plotOutput("plot_markov", height = "400px")
        )
      ),
      column(
        width = 6,
        wellPanel(
          h5("Ranking de I_MPCS"),
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
          h4("Analisis de Sensibilidad"),
          p("Evalua la robustez de la recomendacion ante cambios en los pesos w1, w2, w3."),
          hr(),
          fluidRow(
            column(
              width = 3,
              numericInput("sens_n_sim", 
                           "Numero de simulaciones:", 
                           value = 10000, 
                           min = 100, 
                           max = 50000,
                           step = 100),
              br(),
              actionButton("run_sensitivity", 
                           "Ejecutar analisis",
                           icon = icon("chart-line"),
                           class = "btn-primary",
                           style = "width: 100%;")
            ),
            column(
              width = 9,
              conditionalPanel(
                condition = "input.run_sensitivity > 0",
                div(style = "padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
                    h5("Distribucion del Indice MPCS"),
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
          h4("Interpretacion automatica"),
          uiOutput("interpretation_text")
        )
      )
    )
  ),
  
  # ============================================================================
  # Pestaña 4: Reporte y Descarga
  # ============================================================================
  nav_panel(
    "4. Reporte y Descarga",
    fluidRow(
      column(
        width = 6,
        wellPanel(
          h4("Descargar datos"),
          p("Descarga la tabla de resultados en formato CSV."),
          downloadButton("download_csv", "Descargar tabla CSV", 
                         class = "btn-primary w-100",
                         icon = icon("file-csv"))
        )
      ),
      column(
        width = 6,
        wellPanel(
          h4("Generar reporte"),
          p("Descarga un reporte completo en formato HTML."),
          downloadButton("download_report", "Descargar reporte HTML", 
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
      showNotification("Formato no soportado.", type = "error")
      return()
    }
    
    rv$data <- datos
    rv$dominio <- "salud"
    showNotification(paste("Datos cargados:", nrow(datos), "filas,", ncol(datos), "columnas"), type = "message")
    
    vars <- names(datos)
    vars_num <- names(datos)[sapply(datos, is.numeric)]
    updateSelectInput(session, "graph_vars", choices = vars_num, selected = vars_num[1:min(5, length(vars_num))])
    updateSelectInput(session, "markov_var", choices = vars, selected = if ("Estado_Markov" %in% vars) "Estado_Markov" else vars[length(vars)])
    updateSelectInput(session, "group_var", choices = c("Ninguno (Global)", vars), selected = if ("Region" %in% vars) "Region" else "Ninguno (Global)")
    updateSelectInput(session, "alpha_vars", choices = vars_num, selected = vars_num[1:min(2, length(vars_num))])
  })
  
  # ==========================================================================
  # CARGA DE DATOS DE DEMOSTRACION (ENDES - Salud)
  # ==========================================================================
  observeEvent(input$load_demo, {
    showNotification("Cargando datos de demostracion (ENDES - Salud)...", type = "message")
    
    generar_y_cargar_demo <- function() {
      datos_demo <- generar_demo_data(1000, 123)
      rv$data <- datos_demo
      rv$dominio <- "salud"
      showNotification("Datos de demostracion generados exitosamente (n=1000).", type = "message")
    }
    
    if (file.exists("data/demo_data.csv")) {
      tryCatch({
        rv$data <- read.csv("data/demo_data.csv")
        rv$dominio <- "salud"
        showNotification("Datos de demostracion cargados desde archivo.", type = "message")
      }, error = function(e) {
        showNotification(paste("Error al leer archivo:", e$message), type = "error")
        generar_y_cargar_demo()
      })
    } else {
      showNotification("No se encontro archivo. Generando datos de demostracion...", type = "message")
      generar_y_cargar_demo()
    }
    
    vars <- names(rv$data)
    vars_num <- names(rv$data)[sapply(rv$data, is.numeric)]
    updateSelectInput(session, "graph_vars", choices = vars_num, selected = vars_num[1:min(5, length(vars_num))])
    updateSelectInput(session, "markov_var", choices = vars, selected = "Estado_Markov")
    updateSelectInput(session, "group_var", choices = c("Ninguno (Global)", vars), selected = "Region")
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
    cat("=== ESTADISTICAS DE LOS DATOS ===\n\n")
    cat("Observaciones (n):", nrow(rv$data), "\n")
    cat("Variables:", ncol(rv$data), "\n")
    cat("Valores faltantes:", sum(is.na(rv$data)), "\n\n")
    cat("--- Tipos de variables ---\n")
    print(table(sapply(rv$data, class)))
    cat("\n--- Primeras variables ---\n")
    print(names(rv$data)[1:min(10, ncol(rv$data))])
  })
  
  # ==========================================================================
  # UI DINAMICA
  # ==========================================================================
  output$graph_vars_ui <- renderUI({
    req(rv$data)
    vars_num <- names(rv$data)[sapply(rv$data, is.numeric)]
    if (length(vars_num) == 0) {
      return(div(class = "alert alert-warning", "No se encontraron variables numericas."))
    }
    checkboxGroupInput("graph_vars", "Variables del sistema:", choices = vars_num, 
                       selected = vars_num[1:min(5, length(vars_num))])
  })
  
  output$alpha_vars_ui <- renderUI({
    req(rv$data)
    vars_num <- names(rv$data)[sapply(rv$data, is.numeric)]
    if (length(vars_num) == 0) {
      return(helpText("No hay variables numericas disponibles para alpha."))
    }
    selectInput("alpha_vars", "Variables para alpha (selecciona 2-3):",
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
    selectInput("markov_var", "Variable de estados Markov:", choices = vars, selected = default)
  })
  
  output$group_var_ui <- renderUI({
    req(rv$data)
    vars <- c("Ninguno (Global)", names(rv$data))
    default <- if ("Region" %in% names(rv$data)) {
      "Region"
    } else if ("Grupo" %in% names(rv$data)) {
      "Grupo"
    } else {
      "Ninguno (Global)"
    }
    selectInput("group_var", "Variable de agrupacion:", choices = vars, selected = default)
  })
  
  output$suma_ponderadores <- renderText({
    total <- sum(input$w1, input$w2, input$w3)
    if (abs(total - 1) < 0.01) paste0(total, " OK") else paste0(total, " (debe sumar 1)")
  })
  
  # ==========================================================================
  # FUNCION PARA DETECTAR ORDEN LOGICO DE ESTADOS (MEJORADA)
  # ==========================================================================
  detectar_orden_estados <- function(estados_unicos) {
    
    # --- Caso 1: Estados de salud (E1, E2, E3, E4, E5) ---
    if (all(grepl("^E[0-9]+$", estados_unicos))) {
      return(estados_unicos[order(as.numeric(gsub("E", "", estados_unicos)))])
    }
    
    # --- Caso 2: Estados de participacion en aula ---
    orden_aula <- c("Nunca", "Rara_vez", "A_veces", "Casi_siempre")
    if (any(orden_aula %in% estados_unicos)) {
      # Devolver SOLO los que estan presentes, en el orden correcto
      return(orden_aula[orden_aula %in% estados_unicos])
    }
    
    # --- Caso 3: Estados con "Bajo", "Medio", "Alto" ---
    orden_bma <- c("Bajo", "Medio", "Alto")
    if (any(orden_bma %in% estados_unicos)) {
      return(orden_bma[orden_bma %in% estados_unicos])
    }
    
    # --- Caso 4: Estados numericos ---
    if (all(grepl("^[0-9]+$", estados_unicos))) {
      return(as.character(sort(as.numeric(estados_unicos))))
    }
    
    # --- Caso 5: Estados con "Nivel 1", "Nivel 2", "Nivel 3" ---
    if (all(grepl("^Nivel [0-9]+$", estados_unicos))) {
      return(estados_unicos[order(as.numeric(gsub("Nivel ", "", estados_unicos)))])
    }
    
    # --- Caso 6: Intentar detectar orden por frecuencia (fallback) ---
    # Si no se puede determinar, usar orden alfabetico CON ADVERTENCIA
    warning("No se pudo determinar el orden logico de estados. Usando orden alfabetico.")
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
  # EJECUTAR MPCS
  # ==========================================================================
  observeEvent(input$run_mpcs, {
    req(rv$data)
    
    if (is.null(input$graph_vars) || length(input$graph_vars) < 5) {
      showNotification("Seleccione al menos 5 variables.", type = "error")
      output$validation_msg <- renderUI({ 
        tags$div(class = "alert alert-danger mt-2", 
                 icon("exclamation-triangle"), 
                 " Seleccione al menos 5 variables.")
      })
      return()
    }
    output$validation_msg <- renderUI({ NULL })
    
    if (is.null(input$markov_var) || input$markov_var == "") {
      showNotification("Seleccione una variable de estados Markov.", type = "error")
      output$validation_msg <- renderUI({ 
        tags$div(class = "alert alert-danger mt-2", 
                 icon("exclamation-triangle"), 
                 " Seleccione una variable de estados Markov.")
      })
      return()
    }
    
    estados <- rv$data[[input$markov_var]]
    if (length(unique(na.omit(estados))) < 3) {
      showNotification("La variable de Markov debe tener al menos 3 estados.", type = "error")
      output$validation_msg <- renderUI({ 
        tags$div(class = "alert alert-danger mt-2", 
                 icon("exclamation-triangle"), 
                 " La variable de Markov debe tener al menos 3 estados.")
      })
      return()
    }
    output$validation_msg <- renderUI({ NULL })
    
    withProgress(message = 'Ejecutando MPCS...', value = 0, {
      data_analysis <- rv$data
      
      if (input$group_var == "Ninguno (Global)") {
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
        incProgress(1 / length(grupos), detail = paste("Procesando grupo:", g))
        
        sub <- data_analysis[data_analysis$Group == g, ]
        
        if (nrow(sub) < 30) {
          showNotification(paste("Grupo", g, "tiene menos de 30 observaciones. Saltando."), type = "message")
          next
        }
        
        # --- Grafo ---
        graph_data <- sub[, input$graph_vars, drop = FALSE]
        graph_res <- calcular_grafo(graph_data, input$graph_vars, input$threshold)
        
        if (is.null(graph_res$graph) || vcount(graph_res$graph) < 2) {
          showNotification(paste("Grupo", g, "no tiene suficientes conexiones en el grafo. Usando valor por defecto."), type = "message")
          graph_res$score <- 0.5
          graph_res$optimal_node <- "Sin nodo"
          graph_res$graph <- NULL
        }
        
        alpha_grupo <- calcular_alpha(sub)
        
        # ============================================================
        # [CORREGIDO] DETECTAR ORDEN LOGICO DE ESTADOS
        # ============================================================
        estados_unicos <- unique(sub[[input$markov_var]])
        estados_unicos <- estados_unicos[!is.na(estados_unicos) & estados_unicos != ""]
        orden_estados <- detectar_orden_estados(estados_unicos)
        
        # ============================================================
        # Markov CON ORDEN CORRECTO
        # ============================================================
        markov_res <- calcular_markov(
          sub[[input$markov_var]], 
          orden_estados = orden_estados,
          umbral_objetivo = 0.50, 
          horizonte = 10
        )
        
        # ============================================================
        # Juegos
        # ============================================================
        games_res <- calcular_juegos(markov_res$mat, input$R_factor, alpha = alpha_grupo)
        
        # ============================================================
        # Indice integrado
        # ============================================================
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
          nodo_optimo = if (is.null(graph_res$optimal_node)) "Sin nodo" else graph_res$optimal_node,
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
          # --- GUARDAR ESTADOS EN ORDEN CORRECTO ---
          estados = colnames(markov_res$mat)
        )
      }
      
      if (length(results_list) == 0) {
        showNotification("No se pudo procesar ningun grupo. Verifique que al menos un grupo tenga >30 observaciones.", type = "error")
        return()
      }
      
      # Calcular sim_nudge para cada grupo
      for (g in names(results_list)) {
        r <- results_list[[g]]
        if (!is.null(r$markov_mat) && !is.null(r$sim_base)) {
          P_n <- r$markov_mat
          for (i in 1:(nrow(P_n)-1)) {
            av <- P_n[i, i] * 0.4
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
          I_Juegos = round(r$I_juegos, 4),
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
      
      showNotification(paste("MPCS ejecutado correctamente para", nrow(results_df), "grupos."), type = "message")
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
             edge.width = 1.5, main = paste("Grafo del sistema -", first_group), cex.main = 0.9)
        legend("topright", legend = c("Nodo optimo", "Otros nodos"), 
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
                   annotate("text", x = 0.5, y = 0.5, label = "No hay datos"))
        }
        ggplot(plot_data, aes(x = .data[["Group"]], fill = .data[[markov_var]])) + 
          geom_bar(position = "fill") +
          scale_fill_brewer(palette = "Set2") + theme_minimal() +
          labs(x = "Grupo", y = "Proporcion", fill = "Estado") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
    }
    
    # --- Ranking ---
    if (!is.null(results) && nrow(results) > 0) {
      p_rank <- function() {
        ggplot(results, aes(x = reorder(Grupo, I_MPCS), y = I_MPCS, fill = Tipo_Nudge)) +
          geom_col(width = 0.7) + coord_flip() + theme_minimal(base_size = 12) +
          labs(x = "Grupo", y = "Indice MPCS", fill = "Tipo Nudge") +
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
        # Determinar etiqueta segun dominio
        y_label <- if (dominio == "educacion" || any(grepl("Participacion", colnames(r$sim_base)))) {
          "Probabilidad de participacion"
        } else {
          "Probabilidad de ocupacion"
        }
        return(graficar_trayectorias_markov(
          r$sim_base, 
          r$sim_nudge, 
          estados = colnames(r$sim_base),
          titulo = paste("Trayectorias de Markov -", first_group),
          y_label = y_label
        ))
      }
      return(ggplot() + theme_void() + 
               annotate("text", x = 0.5, y = 0.5, label = "No hay datos para trayectorias"))
    }
    
    # --- Arbol de Markov (CON NUEVA VERSION SIMPLE) ---
    p_arbol <- function() {
      if (!is.null(r$markov_mat)) {
        # Obtener estados en el orden correcto
        if (!is.null(r$estados) && length(r$estados) == ncol(r$markov_mat)) {
          estados_arbol <- r$estados
        } else if (!is.null(colnames(r$markov_mat))) {
          estados_arbol <- colnames(r$markov_mat)
        } else {
          estados_arbol <- paste0("E", 1:ncol(r$markov_mat))
        }
        
        return(graficar_arbol_markov(
          P = r$markov_mat,
          estados = estados_arbol,
          umbral_prob = 0.05,
          titulo = paste("Transiciones de Markov -", first_group)
        ))
      }
      return(ggplot() + theme_void() + 
               annotate("text", x = 0.5, y = 0.5, label = "No hay datos para arbol de Markov"))
    }
    
    # --- Juego evolutivo ---
    p_juego <- function() {
      if (!is.null(r$alpha)) {
        return(graficar_juego_evolutivo(
          alpha = r$alpha,
          p_star = r$I_juegos,
          titulo = paste("Dinamica Replicadora -", first_group)
        ))
      }
      return(ggplot() + theme_void() + 
               annotate("text", x = 0.5, y = 0.5, label = "No hay datos para teoria de juegos"))
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
      formatRound(columns = c("I_Markov_H10", "I_Juegos", "I_MPCS", "k"), digits = 4)
  })
  
  output$plot_graph <- renderPlot({
    req(rv$plots)
    if (!is.null(rv$plots$graph)) rv$plots$graph() else {
      plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(0, 0, "No se pudo generar el grafo.")
    }
  })
  
  output$plot_states <- renderPlot({
    req(rv$plots)
    if (!is.null(rv$plots$states)) rv$plots$states() else {
      ggplot() + theme_void() + 
        annotate("text", x = 0.5, y = 0.5, label = "No se pudo generar el grafico.")
    }
  })
  
  output$plot_markov <- renderPlot({
    req(rv$plots)
    if (!is.null(rv$plots$markov)) rv$plots$markov() else {
      ggplot() + theme_void() + 
        annotate("text", x = 0.5, y = 0.5, label = "No se pudo generar el grafico.")
    }
  })
  
  output$plot_ranking <- renderPlot({
    req(rv$plots)
    if (!is.null(rv$plots$rank)) rv$plots$rank() else {
      ggplot() + theme_void() + 
        annotate("text", x = 0.5, y = 0.5, label = "No se pudo generar el grafico.")
    }
  })
  
  output$plot_arbol_markov <- renderPlot({
    req(rv$plots)
    if (!is.null(rv$plots$arbol)) {
      rv$plots$arbol()
    } else {
      ggplot() + theme_void() + 
        annotate("text", x = 0.5, y = 0.5, label = "No se pudo generar el arbol de Markov.")
    }
  })
  
  output$plot_juego_evolutivo <- renderPlot({
    req(rv$plots)
    if (!is.null(rv$plots$juego)) rv$plots$juego() else {
      ggplot() + theme_void() + 
        annotate("text", x = 0.5, y = 0.5, label = "No se pudo generar el grafico de juegos.")
    }
  })
  
  # ==========================================================================
  # ANALISIS DE SENSIBILIDAD
  # ==========================================================================
  observeEvent(input$run_sensitivity, {
    
    if (is.null(rv$results_df) || nrow(rv$results_df) == 0) {
      showNotification("Primero ejecuta el MPCS para obtener resultados.", type = "error")
      return()
    }
    
    top_group_idx <- which.max(rv$results_df$I_MPCS)
    top_group <- rv$results_df$Grupo[top_group_idx]
    
    r <- rv$results[[as.character(top_group)]]
    
    if (is.null(r)) {
      showNotification("No se encontraron resultados para el grupo seleccionado.", type = "error")
      return()
    }
    
    I_Grafo <- r$I_grafo
    I_Markov <- r$I_markov
    I_Juegos <- r$I_juegos
    
    if (is.na(I_Grafo)) I_Grafo <- 0.5
    if (is.na(I_Markov)) I_Markov <- 0.5
    if (is.na(I_Juegos)) I_Juegos <- 0.5
    
    showNotification(
      paste("Ejecutando analisis de sensibilidad para el grupo:", top_group),
      type = "message",
      duration = 2
    )
    
    withProgress(message = 'Analizando sensibilidad...', value = 0, {
      
      resultado_sens <- analisis_sensibilidad(
        I_Grafo = I_Grafo,
        I_Markov = I_Markov,
        I_Juegos = I_Juegos,
        n_sim = input$sens_n_sim,
        R = input$R_factor,
        escala = 1.5,
        seed = 123
      )
      
      incProgress(0.8, detail = "Generando graficos...")
      
      rv$sens_resultado <- resultado_sens
      
      incProgress(1, detail = "Completado")
    })
    
    showNotification("Analisis de sensibilidad completado.", type = "message")
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
               label = paste0("Pesos por defecto\nI_MPCS = ", 
                             round(rv$sens_resultado$I_MPCS_default, 4)),
               hjust = 0, size = 3.5, color = "#C0392B") +
      scale_fill_manual(values = colores_filtrados) +
      labs(
        title = "Analisis de Sensibilidad del Indice MPCS",
        subtitle = paste0(rv$sens_resultado$n_sim, " combinaciones aleatorias de pesos w1, w2, w3"),
        x = "Indice MPCS",
        y = "Frecuencia",
        fill = "Tipo de nudge"
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
    
    cat("=== RESULTADOS DEL ANALISIS DE SENSIBILIDAD ===\n\n")
    cat(sprintf("Numero de simulaciones: %d\n", s$n_sim))
    cat(sprintf("Media de I_MPCS: %.4f\n", s$stats$media))
    cat(sprintf("Desviacion estandar: %.4f\n", s$stats$sd))
    cat(sprintf("Coeficiente de variacion: %.1f%%\n", s$stats$cv))
    cat(sprintf("Intervalo de confianza 95%%: [%.4f, %.4f]\n", 
                s$ic_inf, s$ic_sup))
    cat(sprintf("I_MPCS con pesos por defecto: %.4f\n", s$I_MPCS_default))
    cat(sprintf("Tipo de nudge por defecto: %s\n", s$tipo_default))
    cat(sprintf("\nTipo de nudge mas frecuente en simulaciones: %s (%.1f%%)\n",
                s$tipo_mas_frecuente, s$pct_mas_frecuente))
    cat(sprintf("Coincidencia con tipo por defecto: %.1f%%\n", 
                s$pct_coincidencia))
    
    cat("\n=== INTERPRETACION ===\n")
    if (s$pct_coincidencia >= 70) {
      cat("La recomendacion de nudge es ALTAMENTE ROBUSTA.\n")
      cat("  El tipo de nudge se mantiene en mas del 70% de las simulaciones.\n")
    } else if (s$pct_coincidencia >= 50) {
      cat("La recomendacion de nudge es ROBUSTA.\n")
      cat("  El tipo de nudge se mantiene en mas del 50% de las simulaciones.\n")
    } else {
      cat("La recomendacion de nudge es SENSIBLE a los pesos.\n")
      cat("  Se recomienda revisar la ponderacion o realizar un analisis adicional.\n")
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
      "<p><b>Resumen de resultados</b></p>",
      "<p>Se analizaron <b>", nrow(rv$results_df), " grupos</b> con un total de <b>", 
      sum(rv$results_df$n), " observaciones</b>.</p>",
      "<p><b>Dominio:</b> ",
      ifelse(any(grepl("Participacion", names(rv$data))), "Educacion - Participacion en el aula", "Salud - Adherencia a medicamentos"),
      "</p>",
      "<hr>",
      "<p><b>Grupo con mayor prioridad:</b> <span style='color:#C0392B;font-weight:bold;'>", 
      top_group$Grupo, "</span> (I_MPCS = ", round(top_group$I_MPCS, 4), ")</p>",
      "<p>El nodo optimo para la intervencion es <b>", top_group$Nodo_Optimo, "</b>.</p>",
      "<p>Se recomienda aplicar un <b>", top_group$Tipo_Nudge, 
      "</b> con intensidad k = ", round(top_group$k, 4), ".</p>",
      "<p><b>Masa critica estimada:</b> ", round(top_group$I_Juegos * 100, 1), 
      "% de adoptantes necesarios para que el cambio sea autosostenible.</p>",
      "<hr>",
      "<p><b>Grupo con menor prioridad:</b> <span style='color:#2ECC71;font-weight:bold;'>", 
      bottom_group$Grupo, "</span> (I_MPCS = ", round(bottom_group$I_MPCS, 4), ")</p>",
      "<p><b>Recomendacion:</b> Enfocar los esfuerzos en el grupo con mayor prioridad (", 
      top_group$Grupo, ") y aplicar un nudge de tipo <b>", top_group$Tipo_Nudge, 
      "</b> con intensidad ", round(top_group$k * 100, 1), "%.</p>",
      "</div>"
    ))
  })
  
  # ==========================================================================
  # DESCARGAS
  # ==========================================================================
  output$download_csv <- downloadHandler(
    filename = function() {
      paste0("MPCS_Resultados_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(rv$results_df)
      write.csv(rv$results_df, file, row.names = FALSE)
    }
  )
  
  output$download_report <- downloadHandler(
    filename = function() {
      paste0("MPCS_Reporte_", Sys.Date(), ".html")
    },
    content = function(file) {
      if (is.null(rv$results_df) || nrow(rv$results_df) == 0) {
        showNotification("No hay resultados para generar el reporte.", type = "error")
        return()
      }
      
      top_group <- rv$results_df[which.max(rv$results_df$I_MPCS), ]
      
      html_lines <- c(
        "<!DOCTYPE html>",
        "<html>",
        "<head>",
        "<meta charset='UTF-8'>",
        "<title>Reporte MPCS</title>",
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
        "<h1>Reporte del MPCS</h1>",
        "<p><b>Fecha de generacion:</b> ", format(Sys.Date(), "%d de %B de %Y"), 
        " a las ", format(Sys.time(), "%H:%M"), "</p>",
        "<p><b>Dominio:</b> ",
        ifelse(any(grepl("Participacion", names(rv$data))), "Educacion - Participacion en el aula", "Salud - Adherencia a medicamentos"),
        "</p>",
        "<hr>",
        "<h2>1. Resumen de resultados</h2>",
        "<p><b>Grupos analizados:</b> ", nrow(rv$results_df), "</p>",
        "<p><b>Total de observaciones:</b> ", sum(rv$results_df$n), "</p>",
        "<br>",
        "<h2>2. Tabla de resultados</h2>",
        "<table>",
        "<thead><tr><th>Grupo</th><th>n</th><th>Alpha</th><th>I_Markov_H10</th><th>I_Juegos</th><th>I_MPCS</th><th>Nodo optimo</th><th>k</th><th>Tipo Nudge</th></tr></thead>",
        "<tbody>"
      )
      
      for (i in 1:nrow(rv$results_df)) {
        row <- rv$results_df[i, ]
        is_top <- row$I_MPCS == max(rv$results_df$I_MPCS)
        bg_color <- if (is_top) " style='background-color: #FFEAA7;'" else ""
        icono <- if (is_top) " [PRIORITARIO]" else ""
        html_lines <- c(
          html_lines,
          paste0("<tr", bg_color, ">"),
          paste0("<td><b>", row$Grupo, icono, "</b></td>"),
          paste0("<td>", row$n, "</td>"),
          paste0("<td>", round(row$Alpha, 3), "</td>"),
          paste0("<td>", round(row$I_Markov_H10, 4), "</td>"),
          paste0("<td>", round(row$I_Juegos, 4), "</td>"),
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
        "<h2>3. Interpretacion de resultados</h2>",
        "<div class='highlight'>",
        "<p><b>Grupo con mayor prioridad de intervencion:</b></p>",
        "<ul>",
        "<li><b>Grupo:</b> ", top_group$Grupo, "</li>",
        "<li><b>Indice MPCS:</b> ", round(top_group$I_MPCS, 4), "</li>",
        "<li><b>Nodo optimo de intervencion:</b> ", top_group$Nodo_Optimo, "</li>",
        "<li><b>Intensidad recomendada (k):</b> ", round(top_group$k, 4), "</li>",
        "<li><b>Tipo de nudge recomendado:</b> <b>", top_group$Tipo_Nudge, "</b></li>",
        "<li><b>Masa critica:</b> ", round(top_group$I_Juegos * 100, 1), "% de adoptantes necesarios</li>",
        "</ul>",
        "</div>",
        "<h2>4. Citacion</h2>",
        "<p>MPCS: A method for the systematic design of behavioral nudges... (Zela Llanque, 2026).</p>",
        "<div class='footer'>",
        "<p>Reporte generado automaticamente por <b>MPCS Calculator</b></p>",
        "<p>Repositorio: <a href='https://github.com/Izela-meth/MPCS_APP' target='_blank'>",
        "https://github.com/Izela-meth/MPCS_APP</a></p>",
        "</div>",
        "</body>",
        "</html>"
      )
      
      writeLines(html_lines, file)
      showNotification("Reporte HTML generado correctamente.", type = "message")
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
    "MPCS Calculator v1.0 - Desarrollado por el equipo MPCS | ",
    tags$a(href = "https://github.com/Izela-meth/MPCS_APP", target = "_blank", "GitHub"),
    " | ",
    tags$a(href = "mailto:izela@unsa.edu.pe", "Contacto")
  )
)

# ==============================================================================
# EJECUTAR LA APLICACION
# ==============================================================================
shinyApp(
  ui = ui,
  server = server
)
