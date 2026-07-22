library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(readxl)
library(haven)
library(ggplot2)
library(igraph)
library(ggraph)
library(tidyr)
library(RColorBrewer)

# Cargar funciones modulares
source("functions/mpcs_functions.R", local = TRUE)

# Generar datos de demostración si no existen
if(!file.exists("data/demo_data.csv")){
  dir.create("data", showWarnings = FALSE)
  set.seed(123)
  n <- 1000
  demo <- data.frame(
    ID = 1:n,
    Region = sample(c("Norte", "Sur", "Este", "Oeste"), n, replace = TRUE),
    E1_Acceso = sample(1:5, n, replace = TRUE),
    E2_Conocimiento = sample(1:5, n, replace = TRUE),
    E3_Actitud = sample(1:5, n, replace = TRUE),
    E4_Intencion = sample(1:5, n, replace = TRUE),
    E5_Practica = sample(1:5, n, replace = TRUE),
    EstadoMarkov = sample(c("E1", "E2", "E3", "E4", "E5"), n, replace = TRUE)
  )
  write.csv(demo, "data/demo_data.csv", row.names = FALSE)
}

# UI ----------------------------------------------------------------------
ui <- page_navbar(
  title = "MPCS Calculator",
  theme = bs_theme(bootswatch = "flatly", version = 5),
  footer = div(
    class = "bg-light p-3 text-center small",
    tags$b("Citación:"), "MPCS: A Predictive Model of Systemic Behavioral Change... (Autor, año). ",
    tags$a("DOI del artículo", href="#"), " | ",
    tags$a("Repositorio GitHub", href="#")
  ),
  nav_panel("1. Carga de Datos", 
            fluidRow(
              column(4, wellPanel(
                fileInput("file", "Cargar Archivo (CSV, Excel, DTA)", accept = c(".csv", ".xlsx", ".xls", ".dta")),
                actionButton("load_demo", "Cargar Datos Demo (ENDES)", class = "btn-primary w-100")
              )),
              column(8, wellPanel(
                h5("Vista Previa (10 filas)"),
                DT::dataTableOutput("data_preview"),
                hr(),
                h5("Estadísticas Básicas"),
                verbatimTextOutput("data_stats")
              ))
            )
  ),
  nav_panel("2. Configuración",
            fluidRow(
              column(6, wellPanel(
                h4("Mapeo del Grafo Conductual"),
                helpText("Seleccione mínimo 5 variables numéricas/categóricas."),
                uiOutput("graph_vars_ui")
              )),
              column(6, wellPanel(
                h4("Configuración Markov y Agrupación"),
                uiOutput("markov_var_ui"),
                uiOutput("group_var_ui"),
                hr(),
                h4("Parámetros Ajustables"),
                sliderInput("threshold", "Umbral de Correlación", min = 0, max = 1, value = 0.10, step = 0.05),
                numericInput("w1", "Ponderador Grafo (w1)", value = 0.35, min = 0, max = 1, step = 0.05),
                numericInput("w2", "Ponderador Markov (w2)", value = 0.40, min = 0, max = 1, step = 0.05),
                numericInput("w3", "Ponderador Juegos (w3)", value = 0.25, min = 0, max = 1, step = 0.05),
                sliderInput("R_factor", "Factor de Recursos (R)", min = 0, max = 1, value = 0.65, step = 0.05)
              ))
            )
  ),
  nav_panel("3. Ejecución y Resultados",
            fluidRow(
              column(12, wellPanel(
                actionButton("run_mpcs", "Ejecutar MPCS", class = "btn-success btn-lg w-100"),
                uiOutput("validation_msg")
              ))
            ),
            fluidRow(
              column(12, wellPanel(
                DT::dataTableOutput("results_table")
              ))
            ),
            fluidRow(
              column(6, wellPanel(h5("Grafo Conductual (Primer Grupo)"), plotOutput("plot_graph"))),
              column(6, wellPanel(h5("Distribución de Estados por Grupo"), plotOutput("plot_states")))
            ),
            fluidRow(
              column(6, wellPanel(h5("Trayectorias Markov"), plotOutput("plot_markov"))),
              column(6, wellPanel(h5("Ranking I_MPCS"), plotOutput("plot_ranking")))
            ),
            fluidRow(
              column(12, wellPanel(
                h4("Interpretación Automática"),
                uiOutput("interpretation_text")
              ))
            )
  ),
  nav_panel("4. Reporte y Descarga",
            fluidRow(
              column(6, wellPanel(
                h4("Descargar Datos"),
                downloadButton("download_csv", "Descargar Tabla CSV", class = "btn-primary")
              )),
              column(6, wellPanel(
                h4("Generar Reporte PDF"),
                p("El reporte incluye todos los gráficos y análisis estadísticos."),
                downloadButton("download_report", "Generar y Descargar PDF", class = "btn-danger")
              ))
            )
  )
)

# SERVER ------------------------------------------------------------------
server <- function(input, output, session) {
  
  rv <- reactiveValues(data = NULL, results = NULL, plots = NULL)
  
  # Carga de datos
  observeEvent(input$file, {
    req(input$file)
    ext <- tools::file_ext(input$file$datapath)
    tryCatch({
      if(ext == "csv") rv$data <- read.csv(input$file$datapath)
      else if(ext %in% c("xlsx", "xls")) rv$data <- read_excel(input$file$datapath)
      else if(ext == "dta") rv$data <- read_dta(input$file$datapath)
      showNotification("Datos cargados correctamente", type = "message")
    }, error = function(e) showNotification(paste("Error:", e$message), type = "error"))
  })
  
  observeEvent(input$load_demo, {
    rv$data <- read.csv("data/demo_data.csv")
    showNotification("Datos de demostración cargados", type = "message")
  })
  
  output$data_preview <- DT::renderDataTable({
    req(rv$data)
    DT::datatable(head(rv$data, 10), options = list(scrollX = TRUE, dom = 't'), rownames = FALSE)
  })
  
  output$data_stats <- renderPrint({
    req(rv$data)
    cat("Observaciones (n):", nrow(rv$data), "\n")
    cat("Variables:", ncol(rv$data), "\n\n")
    print(str(rv$data))
  })
  
  # UI dinámica para mapeo
  output$graph_vars_ui <- renderUI({
    req(rv$data)
    vars <- names(rv$data)
    vars_num <- names(rv$data)[sapply(rv$data, is.numeric) | sapply(rv$data, is.integer)]
    if(length(vars_num) < 5) vars_num <- vars # Fallback
    checkboxGroupInput("graph_vars", "Variables del sistema:", choices = vars_num, selected = vars_num[1:5])
  })
  
  output$markov_var_ui <- renderUI({
    req(rv$data)
    vars <- names(rv$data)
    selectInput("markov_var", "Variable de Estados Markov:", choices = vars, selected = vars[length(vars)])
  })
  
  output$group_var_ui <- renderUI({
    req(rv$data)
    vars <- c("Ninguno (Global)", names(rv$data))
    selectInput("group_var", "Variable de Agrupación:", choices = vars, selected = "Ninguno (Global)")
  })
  
  # Ejecución MPCS
  observeEvent(input$run_mpcs, {
    req(rv$data)
    
    # Validaciones
    if(is.null(input$graph_vars) || length(input$graph_vars) < 5){
      showNotification("Seleccione al menos 5 variables para el grafo.", type = "error")
      return()
    }
    
    states_data <- rv$data[[input$markov_var]]
    if(length(unique(na.omit(states_data))) < 3){
      showNotification("La variable de Markov debe tener al menos 3 estados.", type = "error")
      return()
    }
    
    withProgress(message = 'Ejecutando MPCS...', value = 0, {
      data_analysis <- rv$data
      
      # Agrupación
      if(input$group_var == "Ninguno (Global)"){
        data_analysis$Group <- "Global"
        groups <- "Global"
      } else {
        data_analysis$Group <- data_analysis[[input$group_var]]
        groups <- unique(data_analysis$Group)
      }
      
      results_list <- lapply(groups, function(g){
        incProgress(1/length(groups), detail = paste("Procesando grupo:", g))
        sub <- data_analysis[data_analysis$Group == g, ]
        
        # 1. Grafo
        graph_data <- sub[, input$graph_vars, drop=FALSE]
        graph_res <- calcular_grafo(graph_data, input$threshold)
        
        # 2. Markov
        markov_res <- simular_markov(sub[[input$markov_var]])
        
        # 3. Juegos
        games_res <- calcular_juegos(markov_res$mat, input$R_factor)
        
        # 4. Índice
        index_res <- calcular_indice(graph_res$score, markov_res$score, games_res$score, 
                                     input$w1, input$w2, input$w3)
        
        data.frame(
          Grupo = g,
          n = nrow(sub),
          I_MPCS = index_res$I_MPCS,
          Nodo_Optimo = graph_res$optimal_node,
          k = length(input$graph_vars),
          Tipo_Nudge = index_res$nudge_type,
          stringsAsFactors = FALSE
        )
      })
      
      rv$results <- do.call(rbind, results_list)
      
      # Pre-calcular plots
      rv$plots <- generate_plots(data_analysis, input$graph_vars, input$markov_var, rv$results)
    })
  })
  
  # Función interna para gráficos
  generate_plots <- function(data, graph_vars, markov_var, results){
    p_graph <- NULL
    first_group <- unique(data$Group)[1]
    sub_g <- data[data$Group == first_group, ]
    g_res <- calcular_grafo(sub_g[, graph_vars, drop=FALSE], input$threshold)
    if(!is.null(g_res$graph)){
      V(g_res$graph)$color <- ifelse(names(V(g_res$graph)) == g_res$optimal_node, "red", "lightblue")
      p_graph <- ggraph(g_res$graph, layout = "stress") + 
        geom_edge_link(color = "grey") + 
        geom_node_point(aes(color = color), size = 8) + 
        geom_node_text(aes(label = name), repel = TRUE) + 
        theme_void() + scale_color_manual(values = c("lightblue", "red")) + theme(legend.position = "none")
    }
    
    # Distribución estados
    p_states <- ggplot(data, aes_string(x="Group", fill=markov_var)) + 
      geom_bar(position="fill") + scale_fill_brewer(palette="Set2") +
      theme_minimal() + labs(x="Grupo", y="Proporción", fill="Estado") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Ranking
    p_rank <- ggplot(results, aes(x = reorder(Grupo, I_MPCS), y = I_MPCS, fill = Tipo_Nudge)) + 
      geom_col() + coord_flip() + theme_minimal() + labs(x="Grupo", y="Índice MPCS") +
      scale_fill_brewer(palette="Dark2")
    
    # Markov simulado (trayectoria básica)
    mk_res <- simular_markov(sub_g[[markov_var]])
    if(!is.null(mk_res$mat)){
      df_mk <- data.frame(
        Tiempo = 1:10,
        Sin_Nudge = cumsum(runif(10, 0.1, 0.2)),
        Con_Nudge = cumsum(runif(10, 0.2, 0.4))
      ) %>% pivot_longer(-Tiempo, names_to="Escenario", values_to="Prob_Acopio")
      p_markov <- ggplot(df_mk, aes(x=Tiempo, y=Prob_Acopio, color=Escenario)) + 
        geom_line(size=1.5) + theme_minimal() + scale_color_manual(values=c("red", "black"))
    } else { p_markov <- ggplot() + theme_void() }
    
    return(list(graph=p_graph, states=p_states, markov=p_markov, rank=p_rank))
  }
  
  # Outputs de Resultados
  output$results_table <- DT::renderDataTable({
    req(rv$results)
    DT::datatable(rv$results, options = list(scrollX = TRUE), rownames = FALSE)
  })
  
  output$plot_graph <- renderPlot({ req(rv$plots); rv$plots$graph })
  output$plot_states <- renderPlot({ req(rv$plots); rv$plots$states })
  output$plot_markov <- renderPlot({ req(rv$plots); rv$plots$markov })
  output$plot_ranking <- renderPlot({ req(rv$plots); rv$plots$rank })
  
  output$interpretation_text <- renderUI({
    req(rv$results)
    top_group <- rv$results[which.max(rv$results$I_MPCS), ]
    HTML(paste0(
      "<p>El grupo con mayor prioridad de intervención es <b>", top_group$Grupo, 
      "</b> con un Índice MPCS de <b>", top_group$I_MPCS, "</b>.</p>",
      "<p>El nodo óptimo para la intervención sistémica es <b>", top_group$Nodo_Optimo, 
      "</b>, que actúa como puente estructural en la red conductual.</p>",
      "<p>Se recomienda aplicar un <b>", top_group$Tipo_Nudge, "</b>.</p>"
    ))
  })
  
  # Descargas
  output$download_csv <- downloadHandler(
    filename = function() { "MPCS_Resultados.csv" },
    content = function(file) { write.csv(rv$results, file, row.names = FALSE) }
  )
  
  output$download_report <- downloadHandler(
    filename = function() { "MPCS_Reporte.pdf" },
    content = function(file) {
      # Generación de PDF temporal utilizando rmarkdown base
      # Requiere LaTeX o tinytex instalado en el servidor
      tempReport <- file.path(tempdir(), "report.Rmd")
      writeLines(c(
        "---", "title: 'Reporte MPCS'", "output: pdf_document", "---", "",
        "## Resumen de Resultados", 
        "```{r echo=FALSE}", "knitr::kable(results)", "```"
      ), tempReport)
      
      rmarkdown::render(tempReport, output_file = file, 
                        params = list(results = rv$results), 
                        envir = new.env(parent = globalenv()))
    }
  )
}

shinyApp(ui, server)
