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

# ==============================================================================
# FUNCIONES PARA GENERAR DATOS DE DEMOSTRACIÓN (SOLO EJEMPLOS)
# ==============================================================================

#' Generar datos simulados de salud (ENDES demo) - SOLO EJEMPLO
generar_demo_data <- function(n = 1000, seed = 123) {
  set.seed(seed)
  data.frame(
    ID = 1:n,
    Region = sample(c("Norte", "Sur", "Este", "Oeste", "Centro"), n, replace = TRUE),
    Edad = sample(25:80, n, replace = TRUE),
    Sexo = sample(c("M", "F"), n, replace = TRUE),
    Acceso_salud = rbinom(n, 1, 0.35) * 0.5 + runif(n, 0, 0.5),
    Tiene_seguro = rbinom(n, 1, 0.60),
    Dx_HTA = rbinom(n, 1, 0.15),
    Dx_DM = rbinom(n, 1, 0.08),
    Adh_farma = runif(n, 0, 1),
    Control_PA = rbinom(n, 1, 0.40),
    IMC = rnorm(n, 26, 4),
    Estado_Markov = sample(c("E1", "E2", "E3", "E4", "E5"), n,
                           replace = TRUE,
                           prob = c(0.45, 0.18, 0.15, 0.12, 0.10))
  )
}

#' Generar datos simulados de aula (educación demo) - SOLO EJEMPLO
generar_datos_aula <- function(n = 200, seed = 456) {
  set.seed(seed)
  
  estados_base <- sample(1:4, n, replace = TRUE, prob = c(0.30, 0.35, 0.25, 0.10))
  
  nota_practica <- round(rnorm(n, mean = 12, sd = 3), 1)
  nota_practica <- pmax(0, pmin(20, nota_practica))
  
  nota_intervenciones <- round(rnorm(n, mean = 10 + 2 * (estados_base - 1), sd = 2.5), 1)
  nota_intervenciones <- pmax(0, pmin(20, nota_intervenciones))
  
  nota_exposicion <- round(rnorm(n, mean = 11 + 1.5 * (estados_base - 1), sd = 2.5), 1)
  nota_exposicion <- pmax(0, pmin(20, nota_exposicion))
  
  asistencia <- round(pmin(100, pmax(40, rnorm(n, mean = 75 + 5 * (estados_base - 1), sd = 10))), 1)
  asistencia_norm <- asistencia / 100 * 20
  eval_continua <- round(rowMeans(cbind(nota_practica, nota_intervenciones, nota_exposicion, asistencia_norm)), 1)
  eval_continua <- pmax(0, pmin(20, eval_continua))
  
  examen_parcial <- round(rnorm(n, mean = 10 + 1.5 * (estados_base - 1), sd = 3), 1)
  examen_parcial <- pmax(0, pmin(20, examen_parcial))
  
  nota_final <- round((eval_continua + examen_parcial) / 2, 1)
  nota_final <- pmax(0, pmin(20, nota_final))
  
  peer_influence <- round(pmin(4, pmax(1, estados_base + rnorm(n, mean = 0, sd = 0.6))), 0)
  peer_influence <- as.numeric(cut(peer_influence, breaks = c(0, 1.5, 2.5, 3.5, 5), labels = 1:4))
  
  teacher_encouragement <- round(pmin(10, pmax(0, rnorm(n, mean = 5 + 1.5 * (estados_base - 1), sd = 2))), 1)
  
  grupo <- sample(c("Sección A", "Sección B", "Sección C"), n, replace = TRUE, prob = c(0.40, 0.35, 0.25))
  
  data.frame(
    ID = 1:n,
    Grupo = grupo,
    Estado_Participacion = factor(estados_base, 
                                   levels = 1:4, 
                                   labels = c("Nunca", "Rara_vez", "A_veces", "Casi_siempre")),
    Nota_Practica = nota_practica,
    Nota_Intervenciones = nota_intervenciones,
    Nota_Exposicion = nota_exposicion,
    Asistencia = asistencia,
    Eval_Continua = eval_continua,
    Examen_Parcial = examen_parcial,
    Nota_Final = nota_final,
    Peer_Influence = peer_influence,
    Teacher_Encouragement = teacher_encouragement,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# FUNCIONES PARA GRÁFICOS ADICIONALES
# ==============================================================================

#' Graficar árbol de transición de Markov
graficar_arbol_markov <- function(P, estados, umbral_prob = 0.01, 
                                  titulo = "Árbol de Transición de Markov") {
  
  if (is.null(P) || is.null(estados) || nrow(P) < 2) {
    return(ggplot() + theme_void() + 
             annotate("text", x = 0.5, y = 0.5, label = "No hay datos para el árbol de Markov"))
  }
  
  m <- nrow(P)
  
  if (is.null(colnames(P))) {
    estados <- paste0("E", 1:m)
  } else {
    estados <- colnames(P)
  }
  
  aristas <- data.frame()
  for (i in 1:m) {
    for (j in 1:m) {
      if (P[i, j] >= umbral_prob && i != j) {
        nombre_from <- estados[i]
        nombre_to <- estados[j]
        if (nchar(nombre_from) > 15) nombre_from <- substr(nombre_from, 1, 12)
        if (nchar(nombre_to) > 15) nombre_to <- substr(nombre_to, 1, 12)
        
        aristas <- rbind(aristas, data.frame(
          from = nombre_from,
          to = nombre_to,
          prob = round(P[i, j] * 100, 1)
        ))
      }
    }
  }
  
  if (nrow(aristas) == 0) {
    return(ggplot() + theme_void() + 
             annotate("text", x = 0.5, y = 0.5, label = "No hay transiciones significativas"))
  }
  
  g <- igraph::graph_from_data_frame(aristas, directed = TRUE)
  
  tryCatch({
    layout <- igraph::layout_as_tree(g, root = 1, circular = FALSE)
  }, error = function(e) {
    layout <- igraph::layout_with_fr(g)
  })
  
  coords <- data.frame(
    x = layout[, 1],
    y = layout[, 2],
    name = igraph::V(g)$name
  )
  
  edge_data <- data.frame()
  for (i in 1:nrow(aristas)) {
    from_idx <- which(coords$name == aristas$from[i])
    to_idx <- which(coords$name == aristas$to[i])
    if (length(from_idx) > 0 && length(to_idx) > 0) {
      edge_data <- rbind(edge_data, data.frame(
        x = coords$x[from_idx],
        y = coords$y[from_idx],
        xend = coords$x[to_idx],
        yend = coords$y[to_idx],
        prob = aristas$prob[i]
      ))
    }
  }
  
  if (nrow(edge_data) == 0) {
    return(ggplot() + theme_void() + 
             annotate("text", x = 0.5, y = 0.5, label = "Error al generar el árbol"))
  }
  
  p <- ggplot() +
    geom_segment(
      data = edge_data,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
      color = "#2C3E50",
      linewidth = 0.8,
      alpha = 0.6
    ) +
    geom_text(
      data = edge_data,
      aes(x = (x + xend)/2, y = (y + yend)/2 + 0.05, 
          label = paste0(prob, "%")),
      size = 3.5,
      color = "#E74C3C",
      fontface = "bold"
    ) +
    geom_point(
      data = coords,
      aes(x = x, y = y),
      size = 20,
      color = "#3498DB",
      fill = "#D6EAF8",
      shape = 21,
      stroke = 1.5
    ) +
    geom_text(
      data = coords,
      aes(x = x, y = y, label = name),
      size = 3.5,
      fontface = "bold",
      color = "#1A3A5C"
    ) +
    labs(
      title = titulo,
      subtitle = paste0("Transiciones con probabilidad ≥ ", round(umbral_prob * 100, 0), "%"),
      x = "", y = ""
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#7F8C8D")
    )
  
  return(p)
}

#' Graficar dinámica replicadora (teoría de juegos)
graficar_juego_evolutivo <- function(alpha = 0.60, p_star = NULL, 
                                     titulo = "Dinámica Replicadora") {
  
  if (is.na(alpha) || is.null(alpha)) alpha <- 0.60
  
  a_AA <- 2.0
  a_AR <- -(0.5 + 0.5 * (1 - alpha))
  a_RA <- 0.5 + 0.5 * alpha
  a_RR <- 0.5 + 0.5 * (1 - alpha)
  
  if (is.null(p_star) || is.na(p_star)) {
    numerador <- a_RR - a_AR
    denominador <- a_AA - a_AR - a_RA + a_RR
    if (denominador != 0) {
      p_star <- numerador / denominador
      p_star <- max(0, min(1, p_star))
    } else {
      p_star <- 0.50
    }
  }
  
  p <- seq(0, 1, length.out = 100)
  f_A <- p * a_AA + (1 - p) * a_AR
  f_R <- p * a_RA + (1 - p) * a_RR
  dp_dt <- p * (1 - p) * (f_A - f_R)
  
  df <- data.frame(p = p, f_A = f_A, f_R = f_R, dp_dt = dp_dt)
  
  p1 <- ggplot(df, aes(x = p)) +
    geom_line(aes(y = f_A, color = "Adoptantes (A)"), linewidth = 1.2) +
    geom_line(aes(y = f_R, color = "Resistentes (R)"), linewidth = 1.2) +
    geom_vline(xintercept = p_star, linetype = "dashed", color = "#C0392B", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "#7F8C8D", linewidth = 0.5) +
    annotate("text", x = min(p_star + 0.05, 0.95), y = max(df$f_A) * 0.85,
             label = paste0("p* = ", round(p_star, 3)),
             color = "#C0392B", fontface = "bold", size = 4.5, hjust = 0) +
    annotate("rect", xmin = 0, xmax = p_star, ymin = -Inf, ymax = Inf,
             alpha = 0.08, fill = "#E74C3C") +
    annotate("rect", xmin = p_star, xmax = 1, ymin = -Inf, ymax = Inf,
             alpha = 0.08, fill = "#2ECC71") +
    annotate("text", x = p_star / 2, y = max(df$f_A) * 0.95,
             label = "Inestable", color = "#C0392B", size = 3.5, hjust = 0.5) +
    annotate("text", x = (1 + p_star) / 2, y = max(df$f_A) * 0.95,
             label = "Estable", color = "#27AE60", size = 3.5, hjust = 0.5) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_y_continuous(labels = scales::number) +
    scale_color_manual(
      name = "Estrategia",
      values = c("Adoptantes (A)" = "#2E86AB", "Resistentes (R)" = "#E84855")
    ) +
    labs(
      title = titulo,
      subtitle = paste0("α = ", round(alpha, 3), " | Masa crítica (p*) = ", round(p_star, 3)),
      x = "Proporción de adoptantes (p)",
      y = "Pago esperado",
      color = "Estrategia"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#7F8C8D")
    )
  
  p2 <- ggplot(df, aes(x = p, y = dp_dt)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "#7F8C8D", linewidth = 0.5) +
    geom_line(color = "#8E44AD", linewidth = 1.2) +
    geom_area(aes(fill = dp_dt > 0), alpha = 0.3) +
    geom_vline(xintercept = p_star, linetype = "dashed", color = "#C0392B", linewidth = 1) +
    annotate("text", x = min(p_star + 0.05, 0.95), y = max(abs(dp_dt), na.rm = TRUE) * 0.85,
             label = paste0("p* = ", round(p_star, 3)),
             color = "#C0392B", fontface = "bold", size = 4, hjust = 0) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_y_continuous(labels = scales::number) +
    scale_fill_manual(values = c("#E74C3C", "#2ECC71"), guide = "none") +
    labs(
      title = "Dinámica Replicadora (dp/dt)",
      subtitle = "Velocidad de cambio en la proporción de adoptantes",
      x = "Proporción de adoptantes (p)",
      y = "dp/dt (tasa de cambio)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#7F8C8D")
    )
  
  p_combinado <- p1 + p2 + 
    patchwork::plot_annotation(
      title = paste0("Análisis de Teoría de Juegos - MPCS"),
      subtitle = paste0("α = ", round(alpha, 3), " | p* = ", round(p_star, 3)),
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#7F8C8D")
      )
    )
  
  return(p_combinado)
}

#' Graficar trayectorias de Markov mejoradas
graficar_trayectorias_markov <- function(sim_base, sim_nudge = NULL, 
                                         estados = NULL, 
                                         titulo = "Trayectorias de Markov",
                                         y_label = "P(Adopción)") {
  
  if (is.null(sim_base)) {
    return(ggplot() + theme_void() + 
             annotate("text", x = 0.5, y = 0.5, label = "No hay datos para trayectorias"))
  }
  
  if (is.null(estados)) {
    estados <- colnames(sim_base)
    if (is.null(estados)) {
      estados <- paste0("E", 1:ncol(sim_base))
    }
  }
  
  df_base <- as.data.frame(sim_base)
  colnames(df_base) <- estados
  df_base$Periodo <- 0:(nrow(df_base) - 1)
  df_base$Escenario <- "Sin nudge"
  
  if (!is.null(sim_nudge) && nrow(sim_nudge) == nrow(sim_base)) {
    df_nudge <- as.data.frame(sim_nudge)
    if (ncol(df_nudge) == length(estados)) {
      colnames(df_nudge) <- estados
      df_nudge$Periodo <- 0:(nrow(df_nudge) - 1)
      df_nudge$Escenario <- "Con nudge"
      df_combined <- rbind(df_base, df_nudge)
    } else {
      df_combined <- df_base
    }
  } else {
    df_combined <- df_base
  }
  
  df_long <- tidyr::pivot_longer(
    df_combined,
    cols = all_of(estados),
    names_to = "Estado",
    values_to = "Probabilidad"
  )
  
  p <- ggplot(df_long, aes(x = Periodo, y = Probabilidad, 
                           color = Estado, linetype = Escenario)) +
    geom_line(linewidth = 1.1) +
    geom_point(data = df_long %>% filter(Periodo %% max(1, round(max(df_long$Periodo)/10)) == 0), 
               size = 1.5, alpha = 0.6) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, max(df_long$Periodo), by = max(1, round(max(df_long$Periodo)/5)))) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = titulo,
      x = "Período",
      y = y_label,
      color = "Estado",
      linetype = "Escenario"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#7F8C8D"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

# ==============================================================================
# GENERAR DATOS DE DEMOSTRACIÓN INICIALES (SOLO PARA DEMO)
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
                    placeholder = "Ningún archivo seleccionado"),
          tags$small("Formatos soportados: CSV, Excel (.xlsx, .xls), Stata (.dta)"),
          hr(),
          h4("O usar datos de demostración"),
          p("Carga datos simulados para probar la aplicación."),
          fluidRow(
            column(6,
              actionButton("load_demo", "Cargar Salud (ENDES)", 
                           class = "btn-primary w-100",
                           icon = icon("heart"))
            ),
            column(6,
              actionButton("load_aula", "Cargar Educación (Aula)", 
                           class = "btn-info w-100",
                           icon = icon("users"))
            )
          ),
          tags$small("Estos son solo ejemplos. La app funciona con cualquier dataset.")
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
                      post = tags$span("  (mayor = más recursos disponibles)")),
          hr(),
          h4("Indicador contextual (α)"),
          helpText("α representa la favorabilidad del contexto para la adopción de la conducta."),
          fluidRow(
            column(6,
              radioButtons("alpha_mode", "Método para α:",
                           choices = c("Seleccionar variables" = "vars",
                                      "Ingresar valor manual" = "manual",
                                      "Automático (primeras numéricas)" = "auto"),
                           selected = "auto")
            ),
            column(6,
              conditionalPanel(
                condition = "input.alpha_mode == 'vars'",
                uiOutput("alpha_vars_ui")
              ),
              conditionalPanel(
                condition = "input.alpha_mode == 'manual'",
                numericInput("alpha_manual", "Valor de α (0-1):",
                             value = 0.60, min = 0, max = 1, step = 0.01)
              ),
              conditionalPanel(
                condition = "input.alpha_mode == 'auto'",
                helpText("La app usará la primera variable numérica disponible.")
              )
            )
          ),
          div(
            class = "alert alert-info mt-2",
            style = "font-size: 0.85em;",
            icon("circle-info"),
            HTML(
              "<b>Nota metodológica.</b><br><br>
              <b>I_Markov:</b> Se calcula como la probabilidad acumulada en los estados 
              avanzados en un horizonte fijo de 10 períodos.<br><br>
              <b>α (contexto):</b> Si seleccionas variables, se normalizan a 0-1 y se promedian. 
              Si ingresas manual, usa ese valor. En modo automático, usa la primera variable numérica."
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
          h5("Distribución de estados por grupo"),
          plotOutput("plot_states", height = "500px")
        )
      )
    ),
    fluidRow(
      column(
        width = 6,
        wellPanel(
          h5("🌳 Árbol de Transición de Markov"),
          plotOutput("plot_arbol_markov", height = "400px")
        )
      ),
      column(
        width = 6,
        wellPanel(
          h5("🎮 Dinámica de Teoría de Juegos"),
          plotOutput("plot_juego_evolutivo", height = "400px")
        )
      )
    ),
    fluidRow(
      column(
        width = 6,
        wellPanel(
          h5("📈 Trayectorias de Markov"),
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
    # Análisis de Sensibilidad
    # ==========================================================================
    fluidRow(
      column(
        width = 12,
        wellPanel(
          h4("🔬 Análisis de Sensibilidad"),
          p("Evalúa la robustez de la recomendación ante cambios en los pesos w1, w2, w3."),
          hr(),
          fluidRow(
            column(
              width = 3,
              numericInput("sens_n_sim", 
                           "Número de simulaciones:", 
                           value = 10000, 
                           min = 100, 
                           max = 50000,
                           step = 100),
              br(),
              actionButton("run_sensitivity", 
                           "▶ Ejecutar análisis",
                           icon = icon("chart-line"),
                           class = "btn-primary",
                           style = "width: 100%;")
            ),
            column(
              width = 9,
              conditionalPanel(
                condition = "input.run_sensitivity > 0",
                div(style = "padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
                    h5("📈 Distribución del Índice MPCS"),
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
    # Interpretación
    # ==========================================================================
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
    plots = NULL,
    sens_resultado = NULL,
    dominio = "generico"
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
    rv$dominio <- "generico"
    showNotification(paste("Datos cargados:", nrow(datos), "filas,", ncol(datos), "columnas"), type = "message")
    
    vars <- names(datos)
    vars_num <- names(datos)[sapply(datos, is.numeric)]
    updateSelectInput(session, "graph_vars", choices = vars_num, selected = vars_num[1:min(5, length(vars_num))])
    updateSelectInput(session, "markov_var", choices = vars, selected = if ("Estado_Markov" %in% vars) "Estado_Markov" else vars[length(vars)])
    updateSelectInput(session, "group_var", choices = c("Ninguno (Global)", vars), selected = if ("Region" %in% vars) "Region" else "Ninguno (Global)")
    updateSelectInput(session, "alpha_vars", choices = vars_num, selected = vars_num[1:min(2, length(vars_num))])
  })
  
  # ==========================================================================
  # CARGA DE DATOS DE DEMOSTRACIÓN (SALUD) - SOLO EJEMPLO
  # ==========================================================================
  observeEvent(input$load_demo, {
    showNotification("Cargando datos de salud (ENDES demo)...", type = "message")
    
    datos_demo <- generar_demo_data(1000, 123)
    rv$data <- datos_demo
    rv$dominio <- "salud"
    
    showNotification(paste("Datos de salud generados:", nrow(datos_demo), "individuos"), type = "message")
    
    vars <- names(datos_demo)
    vars_num <- names(datos_demo)[sapply(datos_demo, is.numeric)]
    updateSelectInput(session, "graph_vars", choices = vars_num, selected = vars_num[1:min(5, length(vars_num))])
    updateSelectInput(session, "markov_var", choices = vars, selected = "Estado_Markov")
    updateSelectInput(session, "group_var", choices = c("Ninguno (Global)", vars), selected = "Region")
    updateSelectInput(session, "alpha_vars", choices = vars_num, selected = c("Acceso_salud", "Tiene_seguro"))
  })
  
  # ==========================================================================
  # CARGA DE DATOS DE DEMOSTRACIÓN (EDUCACIÓN) - SOLO EJEMPLO
  # ==========================================================================
  observeEvent(input$load_aula, {
    showNotification("Generando datos de educación (participación en aula)...", type = "message")
    
    datos_aula <- generar_datos_aula(n = 200, seed = 456)
    rv$data <- datos_aula
    rv$dominio <- "educacion"
    
    showNotification(paste("Datos de aula generados:", nrow(datos_aula), "estudiantes"), type = "message")
    
    vars <- names(datos_aula)
    vars_num <- names(datos_aula)[sapply(datos_aula, is.numeric)]
    updateSelectInput(session, "graph_vars", 
                      choices = vars_num, 
                      selected = c("Eval_Continua", "Examen_Parcial", "Nota_Final", 
                                   "Nota_Practica", "Nota_Intervenciones"))
    updateSelectInput(session, "markov_var", 
                      choices = vars, 
                      selected = "Estado_Participacion")
    updateSelectInput(session, "group_var", 
                      choices = c("Ninguno (Global)", vars), 
                      selected = "Grupo")
    updateSelectInput(session, "alpha_vars", 
                      choices = vars_num, 
                      selected = c("Peer_Influence", "Teacher_Encouragement"))
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
    checkboxGroupInput("graph_vars", "Variables del sistema:", choices = vars_num, 
                       selected = vars_num[1:min(5, length(vars_num))])
  })
  
  output$alpha_vars_ui <- renderUI({
    req(rv$data)
    vars_num <- names(rv$data)[sapply(rv$data, is.numeric)]
    if (length(vars_num) == 0) {
      return(helpText("No hay variables numéricas disponibles para α."))
    }
    selectInput("alpha_vars", "Variables para α (selecciona 2-3):",
                choices = vars_num, multiple = TRUE,
                selected = vars_num[1:min(2, length(vars_num))])
  })
  
  output$markov_var_ui <- renderUI({
    req(rv$data)
    vars <- names(rv$data)
    default <- if ("Estado_Participacion" %in% vars) {
      "Estado_Participacion"
    } else if ("Estado_Markov" %in% vars) {
      "Estado_Markov"
    } else {
      vars[1]
    }
    selectInput("markov_var", "Variable de estados Markov:", choices = vars, selected = default)
  })
  
  output$group_var_ui <- renderUI({
    req(rv$data)
    vars <- c("Ninguno (Global)", names(rv$data))
    default <- if ("Grupo" %in% names(rv$data)) {
      "Grupo"
    } else if ("Region" %in% names(rv$data)) {
      "Region"
    } else {
      "Ninguno (Global)"
    }
    selectInput("group_var", "Variable de agrupación:", choices = vars, selected = default)
  })
  
  output$suma_ponderadores <- renderText({
    total <- sum(input$w1, input$w2, input$w3)
    if (abs(total - 1) < 0.01) paste0(total, " ✅") else paste0(total, " ⚠️ (debe sumar 1)")
  })
  
  # ==========================================================================
  # FUNCIÓN PARA CALCULAR ALPHA (COMPLETAMENTE GENÉRICA)
  # ==========================================================================
  calcular_alpha <- function(sub) {
    alpha <- 0.60  # valor por defecto
    
    if (input$alpha_mode == "manual") {
      # Usar valor ingresado por el usuario
      alpha <- input$alpha_manual
      if (is.na(alpha) || alpha < 0 || alpha > 1) alpha <- 0.60
      
    } else if (input$alpha_mode == "vars" && !is.null(input$alpha_vars) && length(input$alpha_vars) > 0) {
      # Usar variables seleccionadas por el usuario
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
      # Modo automático: usar la primera variable numérica disponible
      vars_num <- names(sub)[sapply(sub, is.numeric)]
      if (length(vars_num) > 0) {
        vals <- sub[[vars_num[1]]]
        rng <- range(vals, na.rm = TRUE)
        if (rng[2] > rng[1]) {
          alpha <- mean((vals - rng[1]) / (rng[2] - rng[1]), na.rm = TRUE)
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
        
        # --- Alpha (contexto) - COMPLETAMENTE GENÉRICO ---
        alpha_grupo <- calcular_alpha(sub)
        
        # --- Markov ---
        markov_res <- calcular_markov(sub[[input$markov_var]], 
                                      umbral_objetivo = 0.50, 
                                      horizonte = 10)
        
        # --- Juegos ---
        games_res <- calcular_juegos(markov_res$mat, input$R_factor, alpha = alpha_grupo)
        
        # --- Índice integrado ---
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
          estados = colnames(markov_res$mat)
        )
      }
      
      if (length(results_list) == 0) {
        showNotification("No se pudo procesar ningún grupo. Verifique que al menos un grupo tenga >30 observaciones.", type = "error")
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
  # FUNCIONES PARA GRÁFICOS
  # ==========================================================================
  generate_plots <- function(data, graph_vars, markov_var, results, results_list, threshold, dominio = "generico") {
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
             edge.width = 1.5, main = paste("Grafo del sistema —", first_group), cex.main = 0.9)
        legend("topright", legend = c("Nodo óptimo", "Otros nodos"), 
               fill = c("#C0392B", "#F0DFC0"), cex = 0.8, bty = "n")
      }
    }
    
    # --- Distribución de estados ---
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
    
    # --- Ranking ---
    if (!is.null(results) && nrow(results) > 0) {
      p_rank <- function() {
        ggplot(results, aes(x = reorder(Grupo, I_MPCS), y = I_MPCS, fill = Tipo_Nudge)) +
          geom_col(width = 0.7) + coord_flip() + theme_minimal(base_size = 12) +
          labs(x = "Grupo", y = "Índice MPCS", fill = "Tipo Nudge") +
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
        # Determinar etiqueta del eje Y según dominio
        y_label <- if (dominio == "salud") {
          "P(Adherencia)"
        } else if (dominio == "educacion") {
          "P(Participación)"
        } else {
          "P(Adopción)"
        }
        return(graficar_trayectorias_markov(
          r$sim_base, 
          r$sim_nudge, 
          colnames(r$sim_base),
          titulo = paste("Trayectorias de Markov —", first_group),
          y_label = y_label
        ))
      }
      return(ggplot() + theme_void() + 
               annotate("text", x = 0.5, y = 0.5, label = "No hay datos para trayectorias"))
    }
    
    # --- Árbol de Markov ---
    p_arbol <- function() {
      if (!is.null(r$markov_mat)) {
        estados <- if (!is.null(r$estados)) r$estados else colnames(r$markov_mat)
        return(graficar_arbol_markov(
          r$markov_mat, 
          estados,
          umbral_prob = 0.05,
          titulo = paste("Árbol de Transición —", first_group)
        ))
      }
      return(ggplot() + theme_void() + 
               annotate("text", x = 0.5, y = 0.5, label = "No hay datos para árbol de Markov"))
    }
    
    # --- Juego evolutivo ---
    p_juego <- function() {
      if (!is.null(r$alpha)) {
        return(graficar_juego_evolutivo(
          alpha = r$alpha,
          p_star = r$I_juegos,
          titulo = paste("Dinámica Replicadora —", first_group)
        ))
      }
      return(ggplot() + theme_void() + 
               annotate("text", x = 0.5, y = 0.5, label = "No hay datos para teoría de juegos"))
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
  
  output$plot_arbol_markov <- renderPlot({
    req(rv$plots)
    if (!is.null(rv$plots$arbol)) rv$plots$arbol() else {
      ggplot() + theme_void() + 
        annotate("text", x = 0.5, y = 0.5, label = "No se pudo generar el árbol de Markov.")
    }
  })
  
  output$plot_juego_evolutivo <- renderPlot({
    req(rv$plots)
    if (!is.null(rv$plots$juego)) rv$plots$juego() else {
      ggplot() + theme_void() + 
        annotate("text", x = 0.5, y = 0.5, label = "No se pudo generar el gráfico de juegos.")
    }
  })
  
  # ==========================================================================
  # ANÁLISIS DE SENSIBILIDAD
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
      paste("Ejecutando análisis de sensibilidad para el grupo:", top_group),
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
      
      incProgress(0.8, detail = "Generando gráficos...")
      
      rv$sens_resultado <- resultado_sens
      
      incProgress(1, detail = "Completado")
    })
    
    showNotification("✅ Análisis de sensibilidad completado.", type = "message")
  })
  
  # Gráfico de sensibilidad
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
        title = "Análisis de Sensibilidad del Índice MPCS",
        subtitle = paste0(rv$sens_resultado$n_sim, " combinaciones aleatorias de pesos w1, w2, w3"),
        x = "Índice MPCS",
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
  
  # Resumen textual de sensibilidad
  output$sens_summary <- renderPrint({
    req(rv$sens_resultado)
    
    s <- rv$sens_resultado
    
    cat("=== RESULTADOS DEL ANÁLISIS DE SENSIBILIDAD ===\n\n")
    cat(sprintf("Número de simulaciones: %d\n", s$n_sim))
    cat(sprintf("Media de I_MPCS: %.4f\n", s$stats$media))
    cat(sprintf("Desviación estándar: %.4f\n", s$stats$sd))
    cat(sprintf("Coeficiente de variación: %.1f%%\n", s$stats$cv))
    cat(sprintf("Intervalo de confianza 95%%: [%.4f, %.4f]\n", 
                s$ic_inf, s$ic_sup))
    cat(sprintf("I_MPCS con pesos por defecto: %.4f\n", s$I_MPCS_default))
    cat(sprintf("Tipo de nudge por defecto: %s\n", s$tipo_default))
    cat(sprintf("\nTipo de nudge más frecuente en simulaciones: %s (%.1f%%)\n",
                s$tipo_mas_frecuente, s$pct_mas_frecuente))
    cat(sprintf("Coincidencia con tipo por defecto: %.1f%%\n", 
                s$pct_coincidencia))
    
    cat("\n=== INTERPRETACIÓN ===\n")
    if (s$pct_coincidencia >= 70) {
      cat("✅ La recomendación de nudge es ALTAMENTE ROBUSTA.\n")
      cat("   El tipo de nudge se mantiene en más del 70% de las simulaciones.\n")
    } else if (s$pct_coincidencia >= 50) {
      cat("✅ La recomendación de nudge es ROBUSTA.\n")
      cat("   El tipo de nudge se mantiene en más del 50% de las simulaciones.\n")
    } else {
      cat("⚠️ La recomendación de nudge es SENSIBLE a los pesos.\n")
      cat("   Se recomienda revisar la ponderación o realizar un análisis adicional.\n")
    }
  })
  
  # ==========================================================================
  # INTERPRETACIÓN
  # ==========================================================================
  output$interpretation_text <- renderUI({
    req(rv$results_df)
    
    top_group <- rv$results_df[which.max(rv$results_df$I_MPCS), ]
    bottom_group <- rv$results_df[which.min(rv$results_df$I_MPCS), ]
    
    # Determinar etiquetas según dominio
    if (rv$dominio == "educacion") {
      etiqueta_comportamiento <- "participación"
      etiqueta_nodo <- "variable educativa"
    } else if (rv$dominio == "salud") {
      etiqueta_comportamiento <- "adherencia"
      etiqueta_nodo <- "variable de salud"
    } else {
      etiqueta_comportamiento <- "comportamiento"
      etiqueta_nodo <- "variable del sistema"
    }
    
    HTML(paste0(
      "<div class='well'>",
      "<p><b>📊 Resumen de resultados</b></p>",
      "<p>Se analizaron <b>", nrow(rv$results_df), " grupos</b> con un total de <b>", 
      sum(rv$results_df$n), " observaciones</b>.</p>",
      "<p><b>Dominio:</b> ", ifelse(rv$dominio == "educacion", "Educación", 
                                   ifelse(rv$dominio == "salud", "Salud", "Genérico")), "</p>",
      "<hr>",
      "<p><b>🔴 Grupo con mayor prioridad:</b> <span style='color:#C0392B;font-weight:bold;'>", 
      top_group$Grupo, "</span> (I_MPCS = ", round(top_group$I_MPCS, 4), ")</p>",
      "<p>El nodo óptimo para la intervención es <b>", top_group$Nodo_Optimo, "</b> (", etiqueta_nodo, ").</p>",
      "<p>Se recomienda aplicar un <b>", top_group$Tipo_Nudge, 
      "</b> con intensidad k = ", round(top_group$k, 4), ".</p>",
      "<p><b>Masa crítica estimada:</b> ", round(top_group$I_Juegos * 100, 1), 
      "% de adoptantes necesarios para que el cambio sea autosostenible.</p>",
      "<hr>",
      "<p><b>🟢 Grupo con menor prioridad:</b> <span style='color:#2ECC71;font-weight:bold;'>", 
      bottom_group$Grupo, "</span> (I_MPCS = ", round(bottom_group$I_MPCS, 4), ")</p>",
      "<p><b>Recomendación:</b> Enfocar los esfuerzos en el grupo con mayor prioridad (", 
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
        "<h1>📊 Reporte del MPCS</h1>",
        "<p><b>Fecha de generación:</b> ", format(Sys.Date(), "%d de %B de %Y"), 
        " a las ", format(Sys.time(), "%H:%M"), "</p>",
        "<p><b>Dominio:</b> ", ifelse(rv$dominio == "educacion", "Educación", 
                                     ifelse(rv$dominio == "salud", "Salud", "Genérico")), "</p>",
        "<hr>",
        "<h2>📋 1. Resumen de resultados</h2>",
        "<p><b>Grupos analizados:</b> ", nrow(rv$results_df), "</p>",
        "<p><b>Total de observaciones:</b> ", sum(rv$results_df$n), "</p>",
        "<br>",
        "<h2>📊 2. Tabla de resultados</h2>",
        "<table>",
        "<thead><tr><th>Grupo</th><th>n</th><th>Alpha</th><th>I_Markov_H10</th><th>I_Juegos</th><th>I_MPCS</th><th>Nodo óptimo</th><th>k</th><th>Tipo Nudge</th></tr></thead>",
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
        "<h2>💡 3. Interpretación de resultados</h2>",
        "<div class='highlight'>",
        "<p><b>🔴 Grupo con mayor prioridad de intervención:</b></p>",
        "<ul>",
        "<li><b>Grupo:</b> ", top_group$Grupo, "</li>",
        "<li><b>Índice MPCS:</b> ", round(top_group$I_MPCS, 4), "</li>",
        "<li><b>Nodo óptimo de intervención:</b> ", top_group$Nodo_Optimo, "</li>",
        "<li><b>Intensidad recomendada (k):</b> ", round(top_group$k, 4), "</li>",
        "<li><b>Tipo de nudge recomendado:</b> <b>", top_group$Tipo_Nudge, "</b></li>",
        "<li><b>Masa crítica:</b> ", round(top_group$I_Juegos * 100, 1), "% de adoptantes necesarios</li>",
        "</ul>",
        "</div>",
        "<h2>📚 4. Citación</h2>",
        "<p>MPCS: A method for the systematic design of behavioral nudges... (Zela Llanque, 2026).</p>",
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
    "MPCS Calculator v1.0 — Desarrollado con ❤️ por el equipo MPCS | ",
    tags$a(href = "https://github.com/Izela-meth/MPCS_APP", target = "_blank", "GitHub"),
    " | ",
    tags$a(href = "mailto:izela@unsa.edu.pe", "Contacto")
  )
)

# ==============================================================================
# EJECUTAR LA APLICACIÓN
# ==============================================================================
shinyApp(
  ui = ui,
  server = server
)
