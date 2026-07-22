# ==============================================================================
# MPCS_APP — Aplicación Shiny del Modelo Predictivo de Cambio Conductual por Sistemas
# ==============================================================================
# Repositorio: https://github.com/Izela-meth/MPCS_APP
# ==============================================================================

# --- Cargar librerías ---
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

# --- Generar datos de demostración si no existen ---
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
                    placeholder = "Ningún archivo seleccionado"),
          tags$small("Formatos soportados: CSV, Excel (.xlsx, .xls), Stata (.dta)"),
          hr(),
          h4("O usar datos de demostración"),
          p("Carga los datos de la ENDES 2024 (simulados) para probar la aplicación."),
          actionButton("load_demo", "Cargar Datos Demo", 
                       class = "btn-primary w-100",
                       icon = icon("play"))
        )
      ),
      column(
        width = 8,
        wellPanel(
          h5("Vista previa (primeras 10 filas)"),
          DTOutput("data_preview"),
          hr(),
          h5("Estadísticas básicas"),
          verbatimTextOutput("data_stats")
        )
      )
    )
  ),
  
  # ============================================================================
  # Pestaña 2: Configuración
  # ============================================================================
  nav_panel(
    "2. Configuración",
    fluidRow(
      column(
        width = 6,
        wellPanel(
          h4("Variables del grafo"),
          helpText("Selecciona mínimo 5 variables numéricas para construir el grafo conductual."),
          uiOutput("graph_vars_ui"),
          hr(),
          h4("Variable de agrupación (opcional)"),
          helpText("Si seleccionas una variable, el análisis se realizará por separado para cada grupo."),
          uiOutput("group_var_ui")
        )
      ),
      column(
        width = 6,
        wellPanel(
          h4("Configuración de Markov"),
          helpText("Selecciona la columna que contiene los estados conductuales."),
          uiOutput("markov_var_ui"),
          hr(),
          h4("Parámetros ajustables"),
          sliderInput("threshold", "Umbral de correlación",
                      min = 0.05, max = 0.30, value = 0.10, step = 0.01,
                      post = tags$span("  (|r| > valor)")),
          hr(),
          h5("Ponderadores del Índice MPCS"),
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
                      post = tags$span("  (mayor = más recursos disponibles)"))
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
          h5("Grafo conductual (primer grupo)"),
          plotOutput("plot_graph", height = "500px")
        )
      ),
      column(
        width = 6,
        wellPanel(
          h5("Distribución de estados por grupo"),
          plotOutput("plot_states", height = "500px")
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
    fluidRow(
      column(
        width = 12,
        wellPanel(
          h4("Interpretación automática"),
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
# SERVER — LÓGICA DE LA APLICACIÓN
# ==============================================================================

server <- function(input, output, session) {
  
  rv <- reactiveValues(
    data = NULL,
    results = NULL,
    results_df = NULL,
    plots = NULL
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
    showNotification(paste("Datos cargados:", nrow(datos), "filas,", ncol(datos), "columnas"), type = "message")
    
    vars <- names(datos)
    vars_num <- names(datos)[sapply(datos, is.numeric)]
    updateSelectInput(session, "graph_vars", choices = vars_num, selected = vars_num[1:min(5, length(vars_num))])
    updateSelectInput(session, "markov_var", choices = vars, selected = if ("Estado_Markov" %in% vars) "Estado_Markov" else vars[length(vars)])
    updateSelectInput(session, "group_var", choices = c("Ninguno (Global)", vars), selected = if ("Region" %in% vars) "Region" else "Ninguno (Global)")
  })
  
  # ==========================================================================
  # CARGA DE DATOS DE DEMOSTRACIÓN
  # ==========================================================================
  observeEvent(input$load_demo, {
    showNotification("Cargando datos de demostración...", type = "message")
    
    generar_y_cargar_demo <- function() {
      datos_demo <- generar_demo_data(1000, 123)
      rv$data <- datos_demo
      showNotification("Datos de demostración generados exitosamente (n=1000).", type = "message")
    }
    
    if (file.exists("data/demo_data.csv")) {
      tryCatch({
        rv$data <- read.csv("data/demo_data.csv")
        showNotification("Datos de demostración cargados desde archivo.", type = "message")
      }, error = function(e) {
        showNotification(paste("Error al leer archivo:", e$message), type = "error")
        generar_y_cargar_demo()
      })
    } else {
      showNotification("No se encontró archivo. Generando datos de demostración...", type = "message")
      generar_y_cargar_demo()
    }
    
    vars <- names(rv$data)
    vars_num <- names(rv$data)[sapply(rv$data, is.numeric)]
    updateSelectInput(session, "graph_vars", choices = vars_num, selected = vars_num[1:min(5, length(vars_num))])
    updateSelectInput(session, "markov_var", choices = vars, selected = "Estado_Markov")
    updateSelectInput(session, "group_var", choices = c("Ninguno (Global)", vars), selected = "Region")
  })
  
  # ==========================================================================
  # SALIDAS DE VISTA PREVIA
  # ==========================================================================
  output$data_preview <- renderDT({
    req(rv$data)
    datatable(head(rv$data, 10), options = list(scrollX = TRUE, dom = 't', pageLength = 10), rownames = FALSE)
  })
  
  output$data_stats <- renderPrint({
    req(rv$data)
    cat("Observaciones (n):", nrow(rv$data), "\n")
    cat("Variables:", ncol(rv$data), "\n")
    cat("Valores faltantes:", sum(is.na(rv$data)), "\n\n")
    print(table(sapply(rv$data, class)))
  })
  
  # ==========================================================================
  # UI DINÁMICA
  # ==========================================================================
  output$graph_vars_ui <- renderUI({
    req(rv$data)
    vars_num <- names(rv$data)[sapply(rv$data, is.numeric)]
    if (length(vars_num) == 0) {
      return(div(class = "alert alert-warning", "⚠️ No se encontraron variables numéricas."))
    }
    checkboxGroupInput("graph_vars", "Variables del sistema:", choices = vars_num, selected = vars_num[1:min(5, length(vars_num))])
  })
  
  output$markov_var_ui <- renderUI({
    req(rv$data)
    vars <- names(rv$data)
    selectInput("markov_var", "Variable de estados Markov:", choices = vars, selected = if ("Estado_Markov" %in% vars) "Estado_Markov" else vars[1])
  })
  
  output$group_var_ui <- renderUI({
    req(rv$data)
    vars <- c("Ninguno (Global)", names(rv$data))
    selectInput("group_var", "Variable de agrupación:", choices = vars, selected = if ("Region" %in% vars) "Region" else "Ninguno (Global)")
  })
  
  output$suma_ponderadores <- renderText({
    total <- sum(input$w1, input$w2, input$w3)
    if (abs(total - 1) < 0.01) paste0(total, " ✅") else paste0(total, " ⚠️ (debe sumar 1)")
  })
  
  # ==========================================================================
  # EJECUTAR MPCS — VERSIÓN CORREGIDA
  # ==========================================================================
  observeEvent(input$run_mpcs, {
    req(rv$data)
    
    # --- Validaciones ---
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
      
      # --- Agrupación ---
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
        
        # --- 1. Grafo ---
        graph_data <- sub[, input$graph_vars, drop = FALSE]
        graph_res <- calcular_grafo(graph_data, input$graph_vars, input$threshold)
        
        if (is.null(graph_res$graph) || vcount(graph_res$graph) < 2) {
          showNotification(paste("Grupo", g, "no tiene suficientes conexiones en el grafo. Usando valor por defecto."), type = "message")
          graph_res$score <- 0.5
          graph_res$optimal_node <- "Sin nodo"
          graph_res$graph <- NULL
        }
        
        # --- 2. Markov ---
        markov_res <- calcular_markov(sub[[input$markov_var]], umbral_objetivo = 0.50)
        
        # --- 3. Juegos ---
        games_res <- calcular_juegos(markov_res$mat, input$R_factor)
        
        # --- 4. Índice MPCS ---
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
          I_MPCS = index_res$I_MPCS,
          k = index_res$k,
          tipo = index_res$nudge_type,
          graph = graph_res$graph,
          graph_data = graph_data,
          markov_mat = markov_res$mat,
          sim_base = markov_res$sim_base,
          dist_actual = markov_res$dist_actual,
          T_base = markov_res$T_base
        )
      }
      
      if (length(results_list) == 0) {
        showNotification("No se pudo procesar ningún grupo. Verifique que al menos un grupo tenga >30 observaciones.", type = "error")
        return()
      }
      
      rv$results <- results_list
      
      # --- Tabla de resultados ---
      results_df <- do.call(rbind, lapply(results_list, function(r) {
        data.frame(
          Grupo = r$grupo,
          n = r$n,
          I_MPCS = round(r$I_MPCS, 4),
          Nodo_Optimo = r$nodo_optimo,
          k = round(r$k, 4),
          Tipo_Nudge = r$tipo,
          stringsAsFactors = FALSE
        )
      }))
      
      rv$results_df <- results_df
      
      # --- Generar gráficos ---
      rv$plots <- generate_plots(
        data = data_analysis,
        graph_vars = input$graph_vars,
        markov_var = input$markov_var,
        results = results_df,
        results_list = results_list,
        threshold = input$threshold
      )
      
      showNotification(paste("MPCS ejecutado correctamente para", nrow(results_df), "grupos."), type = "message")
    })
  })
  
  # ==========================================================================
  # FUNCIONES PARA GRÁFICOS
  # ==========================================================================
  generate_plots <- function(data, graph_vars, markov_var, results, results_list, threshold) {
    p_graph <- NULL
    p_states <- NULL
    p_markov <- NULL
    p_rank <- NULL
    
    if (is.null(data) || nrow(data) == 0 || is.null(results) || nrow(results) == 0) {
      return(list(graph = NULL, states = NULL, markov = NULL, rank = NULL))
    }
    
    first_group <- results$Grupo[1]
    r <- results_list[[as.character(first_group)]]
    
    # --- Gráfico 1: Grafo ---
    if (!is.null(r$graph) && vcount(r$graph) > 0) {
      V(r$graph)$color <- ifelse(V(r$graph)$name == r$nodo_optimo, "#C0392B", "#F0DFC0")
      V(r$graph)$size <- ifelse(V(r$graph)$name == r$nodo_optimo, 25, 15)
      p_graph <- function() {
        plot(r$graph, layout = layout_with_fr(r$graph), vertex.label.cex = 0.8,
             vertex.label.color = "black", vertex.label.dist = 2, edge.color = "gray60",
             edge.width = 1.5, main = paste("Grafo del sistema —", first_group), cex.main = 0.9)
        legend("topright", legend = c("Nodo óptimo", "Otros nodos"), 
               fill = c("#C0392B", "#F0DFC0"), cex = 0.8, bty = "n")
      }
    }
    
    # --- Gráfico 2: Distribución de estados ---
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
          labs(x = "Grupo", y = "Proporción", fill = "Estado") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
    }
    
    # --- Gráfico 3: Ranking ---
    if (!is.null(results) && nrow(results) > 0) {
      p_rank <- function() {
        ggplot(results, aes(x = reorder(Grupo, I_MPCS), y = I_MPCS, fill = Tipo_Nudge)) +
          geom_col(width = 0.7) + coord_flip() + theme_minimal(base_size = 12) +
          labs(x = "Grupo", y = "Índice MPCS", fill = "Tipo Nudge") +
          scale_fill_manual(
            values = c("Informativo" = "#74B3CE", "Estructural" = "#2E86AB",
                       "Normativo" = "#E84855", "Sistémico multi-nudge" = "#1A3A5C")
          ) +
          geom_text(aes(label = round(I_MPCS, 3)), hjust = -0.2, size = 3.5) +
          theme(legend.position = "bottom")
      }
    }
    
    # --- Gráfico 4: Trayectorias de Markov ---
    p_markov <- function() {
      r <- results_list[[as.character(first_group)]]
      
      if (!is.null(r$sim_base) && nrow(r$sim_base) > 0) {
        m <- ncol(r$sim_base)
        if (m >= 2) {
          adh_base <- r$sim_base[, m] + r$sim_base[, max(1, m-1)]
          t <- 1:length(adh_base)
          
          P_n <- r$markov_mat
          if (!is.null(P_n)) {
            for (i in 1:(nrow(P_n)-1)) {
              av <- P_n[i, i] * r$k
              P_n[i, i] <- P_n[i, i] - av
              P_n[i, i+1] <- P_n[i, i+1] + av
              P_n[i, ] <- P_n[i, ] / sum(P_n[i, ])
            }
            sim_nudge <- r$sim_base
            for (j in 2:nrow(sim_nudge)) {
              sim_nudge[j, ] <- sim_nudge[j-1, ] %*% P_n
            }
            adh_nudge <- sim_nudge[, m] + sim_nudge[, max(1, m-1)]
          } else {
            adh_nudge <- adh_base * 1.3
          }
          
          df <- data.frame(
            Tiempo = c(t, t),
            Adherencia = c(adh_base, adh_nudge),
            Escenario = rep(c("Sin nudge", "Con nudge"), each = length(t))
          )
          
          return(ggplot(df, aes(x = Tiempo, y = Adherencia, color = Escenario)) +
                   geom_line(linewidth = 1.2) + 
                   geom_hline(yintercept = 0.50, linetype = "dashed", color = "gray50") +
                   scale_y_continuous(labels = scales::percent) + 
                   theme_minimal(base_size = 12) +
                   scale_color_manual(
                     values = c("Sin nudge" = "#E74C3C", "Con nudge" = "#2ECC71")
                   ) +
                   labs(x = "Período", y = "P(Adherencia ≥ 50%)") + 
                   theme(legend.position = "bottom"))
        }
      }
      
      # Fallback: simulación simple
      t <- 1:15
      df <- data.frame(
        Tiempo = t,
        Sin_Nudge = pnorm(t, mean = 8, sd = 3),
        Con_Nudge = pnorm(t, mean = 6, sd = 2.5)
      ) %>%
        pivot_longer(-Tiempo, names_to = "Escenario", values_to = "Prob_Adherencia")
      
      ggplot(df, aes(x = Tiempo, y = Prob_Adherencia, color = Escenario)) +
        geom_line(linewidth = 1.2) + 
        geom_hline(yintercept = 0.50, linetype = "dashed", color = "gray50") +
        scale_y_continuous(labels = scales::percent) + 
        theme_minimal(base_size = 12) +
        scale_color_manual(values = c("Sin_Nudge" = "#E74C3C", "Con_Nudge" = "#2ECC71")) +
        labs(x = "Período", y = "P(Adherencia ≥ 50%)") + 
        theme(legend.position = "bottom")
    }
    
    return(list(graph = p_graph, states = p_states, markov = p_markov, rank = p_rank))
  }
  
  # ==========================================================================
  # SALIDAS DE RESULTADOS
  # ==========================================================================
  output$results_table <- renderDT({
    req(rv$results_df)
    datatable(rv$results_df, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE) %>%
      formatRound(columns = c("I_MPCS", "k"), digits = 4)
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
        annotate("text", x = 0.5, y = 0.5, label = "No se pudo generar el gráfico.")
    }
  })
  
  output$plot_markov <- renderPlot({
    req(rv$plots)
    if (!is.null(rv$plots$markov)) rv$plots$markov() else {
      ggplot() + theme_void() + 
        annotate("text", x = 0.5, y = 0.5, label = "No se pudo generar el gráfico.")
    }
  })
  
  output$plot_ranking <- renderPlot({
    req(rv$plots)
    if (!is.null(rv$plots$rank)) rv$plots$rank() else {
      ggplot() + theme_void() + 
        annotate("text", x = 0.5, y = 0.5, label = "No se pudo generar el gráfico.")
    }
  })
  
  output$interpretation_text <- renderUI({
    req(rv$results_df)
    
    top_group <- rv$results_df[which.max(rv$results_df$I_MPCS), ]
    bottom_group <- rv$results_df[which.min(rv$results_df$I_MPCS), ]
    
    HTML(paste0(
      "<div class='well'>",
      "<p><b>📊 Resumen de resultados</b></p>",
      "<p>Se analizaron <b>", nrow(rv$results_df), " grupos</b> con un total de <b>", 
      sum(rv$results_df$n), " observaciones</b>.</p>",
      "<hr>",
      "<p><b>🔴 Grupo con mayor prioridad:</b> <span style='color:#C0392B;font-weight:bold;'>", 
      top_group$Grupo, "</span> (I_MPCS = ", round(top_group$I_MPCS, 4), ")</p>",
      "<p>El nodo óptimo para la intervención es <b>", top_group$Nodo_Optimo, "</b>.</p>",
      "<p>Se recomienda aplicar un <b>", top_group$Tipo_Nudge, 
      "</b> con intensidad k = ", round(top_group$k, 4), ".</p>",
      "<hr>",
      "<p><b>🟢 Grupo con menor prioridad:</b> <span style='color:#2ECC71;font-weight:bold;'>", 
      bottom_group$Grupo, "</span> (I_MPCS = ", round(bottom_group$I_MPCS, 4), ")</p>",
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
        "<h1>📊 Reporte del MPCS</h1>",
        "<p><b>Fecha de generación:</b> ", format(Sys.Date(), "%d de %B de %Y"), 
        " a las ", format(Sys.time(), "%H:%M"), "</p>",
        "<hr>",
        "<h2>📋 1. Resumen de resultados</h2>",
        "<p><b>Grupos analizados:</b> ", nrow(rv$results_df), "</p>",
        "<p><b>Total de observaciones:</b> ", sum(rv$results_df$n), "</p>",
        "<br>",
        "<h2>📊 2. Tabla de resultados</h2>",
        "<table>",
        "<thead><tr><th>Grupo</th><th>n</th><th>I_MPCS</th><th>Nodo óptimo</th><th>k</th><th>Tipo Nudge</th></tr></thead>",
        "<tbody>"
      )
      
      for (i in 1:nrow(rv$results_df)) {
        row <- rv$results_df[i, ]
        is_top <- row$I_MPCS == max(rv$results_df$I_MPCS)
        bg_color <- if (is_top) " style='background-color: #FFEAA7;'" else ""
        icono <- if (is_top) " 🔴" else ""
        html_lines <- c(
          html_lines,
          paste0("<tr", bg_color, ">"),
          paste0("<td><b>", row$Grupo, icono, "</b></td>"),
          paste0("<td>", row$n, "</td>"),
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
        "<h2>💡 3. Interpretación de resultados</h2>",
        "<div class='highlight'>",
        "<p><b>🔴 Grupo con mayor prioridad de intervención:</b></p>",
        "<ul>",
        "<li><b>Grupo:</b> ", top_group$Grupo, "</li>",
        "<li><b>Índice MPCS:</b> ", round(top_group$I_MPCS, 4), "</li>",
        "<li><b>Nodo óptimo de intervención:</b> ", top_group$Nodo_Optimo, "</li>",
        "<li><b>Intensidad recomendada (k):</b> ", round(top_group$k, 4), "</li>",
        "<li><b>Tipo de nudge recomendado:</b> <b>", top_group$Tipo_Nudge, "</b></li>",
        "</ul>",
        "</div>",
        "<h2>📚 4. Citación</h2>",
        "<p>MPCS: A Predictive Model of Systemic Behavioral Change... (Autor, año).</p>",
        "<div class='footer'>",
        "<p>Reporte generado automáticamente por <b>MPCS Calculator</b></p>",
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
# EJECUTAR LA APLICACIÓN CON FOOTER
# ==============================================================================

app <- shinyApp(ui = ui, server = server)

app <- tagList(
  app,
  tags$footer(
    class = "bg-light p-3 text-center small border-top mt-4",
    tags$b("Citación:"), 
    "MPCS: A Predictive Model of Systemic Behavioral Change... (Autor, año). ",
    tags$a("DOI del artículo", href = "#"), " | ",
    tags$a("Repositorio GitHub", href = "https://github.com/Izela-meth/MPCS_APP")
  )
)

app
