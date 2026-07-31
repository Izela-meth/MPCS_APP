# ============================================================================
# mpcs_functions.R — Funciones modulares del MPCS
# ============================================================================
# Este archivo contiene todas las funciones matemáticas del Modelo Predictivo
# de Cambio Conductual por Sistemas (MPCS)
# ============================================================================

# ============================================================================
# 1. calcular_grafo — Construye el grafo conductual y calcula centralidades
# ============================================================================

#' Calcular grafo conductual y centralidad
#'
#' @param datos data.frame con las variables del sistema
#' @param variables vector con nombres de columnas a incluir en el grafo
#' @param umbral valor mínimo de correlación para incluir aristas (defecto: 0.10)
#' @return lista con grafo, centralidad, nodo óptimo e índice de impacto
#' @export
calcular_grafo <- function(datos, variables, umbral = 0.10) {
  
  if (missing(datos) || missing(variables)) {
    warning("Se requieren datos y variables")
    return(list(
      graph = NULL,
      centralidad = NULL,
      optimal_node = NA,
      score = 0.5
    ))
  }
  
  if (length(variables) < 3) {
    warning("Se recomiendan al menos 5 variables para un grafo estable. Se usarán ", length(variables), " variables.")
  }
  
  df <- datos[, variables, drop = FALSE]
  df <- df[, colSums(is.na(df)) < nrow(df) * 0.5, drop = FALSE]
  
  if (ncol(df) < 3) {
    warning("Menos de 3 variables válidas después de limpiar datos faltantes")
    return(list(
      graph = NULL,
      centralidad = NULL,
      optimal_node = NA,
      score = 0.5
    ))
  }
  
  mat_cor <- tryCatch({
    cor(df, use = "pairwise.complete.obs", method = "spearman")
  }, error = function(e) {
    warning("Error al calcular la matriz de correlación: ", e$message)
    return(NULL)
  })
  
  if (is.null(mat_cor)) {
    return(list(
      graph = NULL,
      centralidad = NULL,
      optimal_node = NA,
      score = 0.5
    ))
  }
  
  aristas <- which(abs(mat_cor) > umbral & mat_cor != 1, arr.ind = TRUE)
  
  if (nrow(aristas) == 0) {
    warning("No se encontraron aristas con el umbral ", umbral)
    return(list(
      graph = NULL,
      centralidad = NULL,
      optimal_node = NA,
      score = 0.5
    ))
  }
  
  aristas_df <- as.data.frame(aristas) %>%
    dplyr::mutate(
      desde = rownames(mat_cor)[row],
      hasta = colnames(mat_cor)[col],
      peso = mat_cor[cbind(row, col)]
    ) %>%
    dplyr::filter(row < col) %>%
    dplyr::select(desde, hasta, peso)
  
  if (nrow(aristas_df) == 0) {
    warning("No se encontraron aristas no redundantes")
    return(list(
      graph = NULL,
      centralidad = NULL,
      optimal_node = NA,
      score = 0.5
    ))
  }
  
  g <- tryCatch({
    igraph::graph_from_data_frame(aristas_df, directed = FALSE)
  }, error = function(e) {
    warning("Error al crear el grafo: ", e$message)
    return(NULL)
  })
  
  if (is.null(g) || igraph::vcount(g) < 2) {
    return(list(
      graph = NULL,
      centralidad = NULL,
      optimal_node = NA,
      score = 0.5
    ))
  }
  
  grado_max <- max(igraph::degree(g))
  
  centr <- data.frame(
    Variable = igraph::V(g)$name,
    Grado = igraph::degree(g),
    Intermediacion = igraph::betweenness(g, normalized = TRUE)
  ) %>%
    dplyr::mutate(
      Impacto = round(0.60 * Intermediacion + 0.40 * (Grado / grado_max), 4)
    ) %>%
    dplyr::arrange(dplyr::desc(Impacto))
  
  return(list(
    graph = g,
    centralidad = centr,
    optimal_node = centr$Variable[1],
    score = centr$Impacto[1]
  ))
}

# ============================================================================
# 2. calcular_markov — Estima cadena de Markov (OPCIÓN B: HORIZONTE FIJO H=10)
# ============================================================================

#' Calcular cadena de Markov y convergencia (Opción B)
#'
#' @param estados vector con los estados de cada individuo
#' @param orden_estados vector con el orden progresivo de estados (opcional)
#' @param umbral_objetivo proporción para considerar convergencia (defecto: 0.50)
#' @param horizonte número de períodos para proyectar (defecto: 10)
#' @return lista con matriz P, simulación, índice de Markov y tiempo de referencia
#' @export
calcular_markov <- function(estados, orden_estados = NULL, 
                            umbral_objetivo = 0.50, horizonte = 10) {
  
  if (missing(estados)) {
    warning("Se requiere el vector de estados")
    return(list(
      mat = NULL,
      dist_actual = NULL,
      sim_base = NULL,
      T_base = horizonte,
      score = 0.50
    ))
  }
  
  estados_clean <- estados[!is.na(estados) & estados != ""]
  
  if (length(estados_clean) < 30) {
    warning("Pocos datos para estimar la cadena de Markov (n = ", length(estados_clean), ")")
    return(list(
      mat = NULL,
      dist_actual = NULL,
      sim_base = NULL,
      T_base = horizonte,
      score = 0.50
    ))
  }
  
  estados_unicos <- unique(estados_clean)
  
  if (is.null(orden_estados)) {
    if (all(grepl("^E[0-9]+$", estados_unicos))) {
      orden_estados <- estados_unicos[order(as.numeric(gsub("E", "", estados_unicos)))]
    } else {
      freq <- table(estados_clean)
      orden_estados <- names(sort(freq, decreasing = TRUE))
    }
  } else {
    faltantes <- setdiff(estados_unicos, orden_estados)
    if (length(faltantes) > 0) {
      orden_estados <- c(orden_estados, faltantes)
    }
  }
  
  m <- length(orden_estados)
  
  if (m < 3) {
    warning("Se necesitan al menos 3 estados para la cadena de Markov")
    return(list(
      mat = NULL,
      dist_actual = NULL,
      sim_base = NULL,
      T_base = horizonte,
      score = 0.50
    ))
  }
  
  freq <- table(estados_clean)
  dist_actual <- rep(0, m)
  names(dist_actual) <- orden_estados
  
  for (i in seq_along(orden_estados)) {
    if (orden_estados[i] %in% names(freq)) {
      dist_actual[i] <- freq[orden_estados[i]] / sum(freq)
    }
  }
  
  P <- matrix(0, nrow = m, ncol = m)
  colnames(P) <- orden_estados
  rownames(P) <- orden_estados
  
  for (i in 1:(m-1)) {
    prob_avance <- 0.30 + 0.20 * (1 - i/m)
    prob_quedarse <- 0.50 - 0.15 * (i/m)
    prob_retroceso <- 0.10 * (1 - i/m)
    
    P[i, i] <- prob_quedarse
    P[i, i+1] <- prob_avance
    
    if (i > 1) {
      P[i, i-1] <- prob_retroceso
    }
    
    resto <- 1 - sum(P[i, ])
    if (resto > 0) {
      otros <- setdiff(1:m, c(i, i+1, if (i > 1) i-1 else NULL))
      if (length(otros) > 0) {
        P[i, otros] <- resto / length(otros)
      }
    }
  }
  
  P[m, m] <- 0.85
  P[m, 1:(m-1)] <- (1 - 0.85) / (m - 1)
  
  P <- P / rowSums(P)
  
  simular <- function(P, v0, n = 30) {
    m <- nrow(P)
    dist <- matrix(0, n + 1, m)
    dist[1, ] <- v0
    for (t in 2:(n + 1)) {
      dist[t, ] <- dist[t-1, ] %*% P
    }
    colnames(dist) <- colnames(P)
    return(dist)
  }
  
  sim_base <- simular(P, dist_actual)
  
  idx_h <- min(horizonte + 1, nrow(sim_base))
  
  if (m >= 2) {
    score <- sim_base[idx_h, m] + sim_base[idx_h, max(1, m-1)]
  } else {
    score <- sim_base[idx_h, m]
  }
  
  score <- max(0, min(1, score))
  
  if (m >= 2) {
    objetivo <- sim_base[, m] + sim_base[, max(1, m-1)]
    T_base <- which(objetivo >= umbral_objetivo)[1]
    if (is.na(T_base) || is.infinite(T_base)) {
      T_base <- 20
    }
  } else {
    T_base <- 10
  }
  
  return(list(
    mat = P,
    dist_actual = dist_actual,
    sim_base = sim_base,
    T_base = T_base,
    score = score,
    horizonte = horizonte
  ))
}

# ============================================================================
# 3. calcular_juegos — Calcula masa crítica (TEORÍA DE JUEGOS DINÁMICA CON α)
# ============================================================================
 
#' Calcular teoría de juegos y masa crítica (dinámica con α)
#'
#' @param P matriz de transición de Markov (opcional, no se usa)
#' @param R_factor factor de recursos (defecto: 0.65) — OBSOLETO
#' @param alpha indicador de acceso a salud (0-1). Si no se proporciona, usa 0.60 por defecto.
#' @return lista con masa crítica e índice
#' @export
calcular_juegos <- function(P = NULL, R_factor = 0.65, alpha = NULL) {
  
  if (is.null(alpha) || is.na(alpha)) {
    alpha <- 0.60
  }
  
  a_AA <- 2.0
  a_AR <- -(0.5 + 0.5 * (1 - alpha))
  a_RA <- 0.5 + 0.5 * alpha
  a_RR <- 0.5 + 0.5 * (1 - alpha)
  
  numerador <- a_RR - a_AR
  denominador <- a_AA - a_AR - a_RA + a_RR
  p_star <- numerador / denominador
  p_star <- max(0, min(1, p_star))
  
  score <- p_star
  
  return(list(
    p_star = p_star,
    score = score,
    alpha = alpha,
    a_AA = a_AA,
    a_AR = a_AR,
    a_RA = a_RA,
    a_RR = a_RR
  ))
}

# ============================================================================
# 4. calcular_indice — Calcula Índice MPCS y tipo de nudge
# ============================================================================

#' Calcular Índice MPCS y tipo de nudge
#'
#' @param I_grafo índice del módulo de grafos
#' @param I_markov índice del módulo de Markov (Opción B: H=10)
#' @param I_juegos índice del módulo de juegos (dinámico con α)
#' @param w1 ponderador para grafos (defecto: 0.35)
#' @param w2 ponderador para Markov (defecto: 0.40)
#' @param w3 ponderador para juegos (defecto: 0.25)
#' @param R_factor factor de recursos (defecto: 0.65)
#' @return lista con I_MPCS, k, tipo de nudge
#' @export
calcular_indice <- function(I_grafo, I_markov, I_juegos, 
                            w1 = 0.35, w2 = 0.40, w3 = 0.25,
                            R_factor = 0.65) {
  
  total <- w1 + w2 + w3
  if (abs(total - 1) > 0.01) {
    warning("Los ponderadores no suman 1. Se normalizarán.")
    w1 <- w1 / total
    w2 <- w2 / total
    w3 <- w3 / total
  }
  
  I_MPCS <- w1 * I_grafo + w2 * I_markov + w3 * I_juegos
  I_MPCS <- max(0, min(1, I_MPCS))
  
  k <- min(1, I_MPCS * R_factor * 1.5)
  
  tipo <- dplyr::case_when(
    k < 0.25 ~ "Informational",
    k < 0.50 ~ "Structural",
    k < 0.75 ~ "Normative",
    TRUE ~ "Systemic multi-nudge"
  )
  
  return(list(
    I_MPCS = I_MPCS,
    k = k,
    nudge_type = tipo
  ))
}

# ============================================================================
# 5. aplicar_nudge — Aplica nudge a la matriz de transición
# ============================================================================

#' Aplicar nudge a la matriz de transición de Markov
#'
#' @param P matriz de transición original
#' @param k intensidad del nudge (0-1)
#' @return matriz de transición modificada
#' @export
aplicar_nudge <- function(P, k) {
  
  if (is.null(P) || !is.matrix(P)) {
    warning("Se requiere una matriz de transición válida")
    return(NULL)
  }
  
  if (k < 0 || k > 1) {
    warning("k debe estar entre 0 y 1. Se usará k = 0.4")
    k <- 0.4
  }
  
  P_n <- P
  m <- nrow(P)
  
  for (i in 1:(m-1)) {
    av <- P[i, i] * k
    P_n[i, i] <- P[i, i] - av
    P_n[i, i+1] <- P[i, i+1] + av
    P_n[i, ] <- P_n[i, ] / sum(P_n[i, ])
  }
  
  return(P_n)
}

# ============================================================================
# 6. generar_demo_data — Genera datos de salud (demostración ENDES)
# ============================================================================

#' Generar datos de demostración (salud - ENDES demo)
#'
#' @param n número de observaciones
#' @param seed semilla para reproducibilidad
#' @return data.frame con datos de demostración
#' @export
generar_demo_data <- function(n = 1000, seed = 123) {
  
  set.seed(seed)
  
  data.frame(
    ID = 1:n,
    Region = sample(c("Norte", "Sur", "Este", "Oeste", "Centro"), n, replace = TRUE),
    Edad = sample(25:80, n, replace = TRUE),
    Sexo = sample(c("M", "F"), n, replace = TRUE),
    Educacion = sample(0:4, n, replace = TRUE),
    Acceso_salud = rbinom(n, 1, 0.35) * 0.5 + runif(n, 0, 0.5),
    Tiene_seguro = rbinom(n, 1, 0.60),
    Dx_HTA = rbinom(n, 1, 0.15),
    Dx_DM = rbinom(n, 1, 0.08),
    Adh_farma = runif(n, 0, 1),
    HTA_Medida = rbinom(n, 1, 0.20),
    Control_PA = rbinom(n, 1, 0.40),
    IMC = rnorm(n, 26, 4),
    Obesidad_Abd = rbinom(n, 1, 0.30),
    Fuma = rbinom(n, 1, 0.15),
    Alcohol = rbinom(n, 1, 0.20),
    Dieta_Sana = runif(n, 0, 1),
    Actividad_Fisica = sample(0:2, n, replace = TRUE),
    Depresion_Bin = rbinom(n, 1, 0.10),
    Ansiedad_Bin = rbinom(n, 1, 0.15),
    Apoyo_Social = sample(1:5, n, replace = TRUE),
    Conocimiento_HTA = sample(0:10, n, replace = TRUE),
    Conocimiento_DM = sample(0:10, n, replace = TRUE),
    Estado_Markov = sample(c("E1", "E2", "E3", "E4", "E5"), n,
                           replace = TRUE,
                           prob = c(0.45, 0.18, 0.15, 0.12, 0.10))
  )
}

# ============================================================================
# 7. graficar_arbol_markov — Árbol de transición de Markov (CORREGIDO)
# ============================================================================

#' Graficar árbol de transición de Markov (VERSIÓN CORREGIDA)
#'
#' @param P matriz de transición (m x m)
#' @param estados vector con nombres de estados
#' @param umbral_prob probabilidad mínima para mostrar una flecha (default: 0.05)
#' @param titulo título del gráfico
#' @return objeto ggplot
#' @export
graficar_arbol_markov <- function(P, estados, umbral_prob = 0.05, 
                                  titulo = "Árbol de Transición de Markov") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Se requiere el paquete ggplot2")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Se requiere el paquete igraph")
  }
  
  if (is.null(P) || is.null(estados) || nrow(P) < 2) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + 
             ggplot2::annotate("text", x = 0.5, y = 0.5, 
                              label = "No hay datos para el árbol de Markov",
                              size = 5, color = "#7F8C8D"))
  }
  
  m <- nrow(P)
  
  if (is.null(colnames(P))) {
    estados <- paste0("E", 1:m)
  } else {
    estados <- colnames(P)
  }
  
  # --- Crear aristas con transiciones relevantes ---
  aristas <- data.frame()
  for (i in 1:m) {
    for (j in 1:m) {
      if (P[i, j] >= umbral_prob && i != j) {
        nombre_from <- estados[i]
        nombre_to <- estados[j]
        if (nchar(nombre_from) > 12) nombre_from <- substr(nombre_from, 1, 10)
        if (nchar(nombre_to) > 12) nombre_to <- substr(nombre_to, 1, 10)
        
        aristas <- rbind(aristas, data.frame(
          from = nombre_from,
          to = nombre_to,
          prob = round(P[i, j] * 100, 1),
          prob_raw = P[i, j]
        ))
      }
    }
  }
  
  if (nrow(aristas) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + 
             ggplot2::annotate("text", x = 0.5, y = 0.5, 
                              label = paste0("No hay transiciones ≥ ", round(umbral_prob*100, 0), "%"),
                              size = 5, color = "#7F8C8D"))
  }
  
  # --- Crear grafo ---
  g <- igraph::graph_from_data_frame(aristas, directed = TRUE)
  
  # --- Layout ---
  tryCatch({
    layout <- igraph::layout_as_tree(g, root = 1, circular = FALSE)
  }, error = function(e) {
    layout <- igraph::layout_with_fr(g, niter = 1000)
  })
  
  # --- Escalar layout ---
  layout[, 1] <- scale(layout[, 1]) * 2.5
  layout[, 2] <- scale(layout[, 2]) * 2.5
  
  coords <- data.frame(
    x = layout[, 1],
    y = layout[, 2],
    name = igraph::V(g)$name
  )
  
  # --- Crear edge_data ---
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
        prob = aristas$prob[i],
        prob_raw = aristas$prob_raw[i]
      ))
    }
  }
  
  if (nrow(edge_data) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + 
             ggplot2::annotate("text", x = 0.5, y = 0.5, 
                              label = "Error al generar el árbol",
                              size = 5, color = "#7F8C8D"))
  }
  
  # --- Tamaño de nodos según grado ---
  node_degrees <- table(c(as.character(aristas$from), as.character(aristas$to)))
  node_size <- ifelse(coords$name %in% names(node_degrees), 
                      18 + node_degrees[coords$name] * 2, 
                      20)
  node_size <- pmin(30, pmax(16, node_size))
  
  # --- Graficar ---
  p <- ggplot2::ggplot() +
    # Flechas de transición
    ggplot2::geom_curve(
      data = edge_data,
      aes(x = x, y = y, xend = xend, yend = yend),
      curvature = 0.15,
      arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
      color = "#2C3E50",
      linewidth = 0.7,
      alpha = 0.5
    ) +
    # Etiquetas de probabilidad (con fondo blanco para legibilidad)
    ggplot2::geom_label(
      data = edge_data,
      aes(x = (x + xend)/2, y = (y + yend)/2 + 0.08, 
          label = paste0(prob, "%")),
      size = 3.5,
      fill = "white",
      color = "#C0392B",
      fontface = "bold",
      label.size = 0.3,
      label.padding = unit(0.15, "lines"),
      alpha = 0.9
    ) +
    # Nodos (estados)
    ggplot2::geom_point(
      data = coords,
      aes(x = x, y = y),
      size = node_size,
      color = "#2C3E50",
      fill = "#D6EAF8",
      shape = 21,
      stroke = 1.5
    ) +
    # Etiquetas de nodos
    ggplot2::geom_text(
      data = coords,
      aes(x = x, y = y, label = name),
      size = 4,
      fontface = "bold",
      color = "#1A3A5C"
    ) +
    ggplot2::labs(
      title = titulo,
      subtitle = paste0("Transiciones con probabilidad ≥ ", round(umbral_prob * 100, 0), "%"),
      x = "", y = ""
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#7F8C8D"),
      plot.margin = margin(10, 20, 10, 20)
    )
  
  return(p)
}

# ============================================================================
# 8. graficar_juego_evolutivo — Dinámica replicadora (teoría de juegos)
# ============================================================================

#' Graficar dinámica replicadora (teoría de juegos)
#'
#' @param alpha indicador de contexto (0-1)
#' @param p_star masa crítica (calculada automáticamente si no se proporciona)
#' @param titulo título del gráfico
#' @return objeto ggplot
#' @export
graficar_juego_evolutivo <- function(alpha = 0.60, p_star = NULL, 
                                     titulo = "Dinámica Replicadora") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Se requiere el paquete ggplot2")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Se requiere el paquete patchwork")
  }
  
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

# ============================================================================
# 9. graficar_trayectorias_markov — Trayectorias de Markov (CORREGIDO)
# ============================================================================

#' Graficar trayectorias de Markov (VERSIÓN CORREGIDA)
#'
#' @param sim_base simulación base
#' @param sim_nudge simulación con nudge (opcional)
#' @param estados nombres de los estados
#' @param titulo título del gráfico
#' @param y_label etiqueta del eje Y
#' @return objeto ggplot
#' @export
graficar_trayectorias_markov <- function(sim_base, sim_nudge = NULL, 
                                         estados = NULL, 
                                         titulo = "Trayectorias de Markov",
                                         y_label = "Probabilidad de ocupación") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Se requiere el paquete ggplot2")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Se requiere el paquete tidyr")
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Se requiere el paquete RColorBrewer")
  }
  
  if (is.null(sim_base)) {
    return(ggplot() + theme_void() + 
             annotate("text", x = 0.5, y = 0.5, 
                     label = "No hay datos para trayectorias",
                     size = 5, color = "#7F8C8D"))
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
  
  n_estados <- length(estados)
  colores <- RColorBrewer::brewer.pal(min(n_estados, 8), "Set1")
  if (n_estados > 8) {
    colores <- colorRampPalette(colores)(n_estados)
  }
  
  p <- ggplot(df_long, aes(x = Periodo, y = Probabilidad, 
                           color = Estado, linetype = Escenario)) +
    geom_line(linewidth = 1.2) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1), 
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1)
    ) +
    scale_x_continuous(
      breaks = seq(0, max(df_long$Periodo), by = max(1, round(max(df_long$Periodo)/6)))
    ) +
    scale_color_manual(values = colores) +
    scale_linetype_manual(values = c("Sin nudge" = "solid", "Con nudge" = "dashed")) +
    labs(
      title = titulo,
      subtitle = "Muestra la proporción esperada de individuos en cada estado a lo largo del tiempo",
      x = "Período (unidades de tiempo)",
      y = y_label,
      color = "Estado conductual",
      linetype = "Escenario"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#7F8C8D"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#EAECEE", linewidth = 0.5)
    )
  
  return(p)
}

# ============================================================================
# 10. analisis_sensibilidad — Análisis de sensibilidad de pesos
# ============================================================================

#' Análisis de sensibilidad del Índice MPCS
#'
#' @param I_Grafo Índice del módulo de grafos (0-1)
#' @param I_Markov Índice del módulo de Markov (0-1)
#' @param I_Juegos Índice del módulo de juegos (0-1)
#' @param n_sim Número de simulaciones (defecto: 10000)
#' @param R Factor de recursos (defecto: 0.65)
#' @param escala Factor de escala (defecto: 1.5)
#' @param seed Semilla para reproducibilidad (defecto: 123)
#' @return Lista con resultados del análisis de sensibilidad
#' @export
analisis_sensibilidad <- function(I_Grafo, I_Markov, I_Juegos,
                                   n_sim = 10000, 
                                   R = 0.65, 
                                   escala = 1.5,
                                   seed = 123) {
  
  if (is.null(I_Grafo) || is.na(I_Grafo)) I_Grafo <- 0.5
  if (is.null(I_Markov) || is.na(I_Markov)) I_Markov <- 0.5
  if (is.null(I_Juegos) || is.na(I_Juegos)) I_Juegos <- 0.5
  
  set.seed(seed)
  pesos_sim <- matrix(NA, n_sim, 3)
  colnames(pesos_sim) <- c("w1", "w2", "w3")
  
  for (i in 1:n_sim) {
    w <- runif(3)
    w <- w / sum(w)
    pesos_sim[i, ] <- w
  }
  
  I_MPCS_sim <- apply(pesos_sim, 1, function(w) {
    w[1] * I_Grafo + w[2] * I_Markov + w[3] * I_Juegos
  })
  
  k_sim <- pmin(1, I_MPCS_sim * R * escala)
  
  tipo_sim <- dplyr::case_when(
    k_sim < 0.25 ~ "Informational",
    k_sim < 0.50 ~ "Structural",
    k_sim < 0.75 ~ "Normative",
    TRUE         ~ "Systemic multi-nudge"
  )
  
  w_default <- c(0.35, 0.40, 0.25)
  I_MPCS_default <- w_default[1] * I_Grafo +
                    w_default[2] * I_Markov +
                    w_default[3] * I_Juegos
  
  k_default <- min(1, I_MPCS_default * R * escala)
  
  tipo_default <- dplyr::case_when(
    k_default < 0.25 ~ "Informational",
    k_default < 0.50 ~ "Structural",
    k_default < 0.75 ~ "Normative",
    TRUE            ~ "Systemic multi-nudge"
  )
  
  freq_tipo <- table(tipo_sim)
  tipo_mas_frecuente <- names(freq_tipo)[which.max(freq_tipo)]
  pct_mas_frecuente <- max(freq_tipo) / n_sim * 100
  
  pct_coincidencia <- sum(tipo_sim == tipo_default) / n_sim * 100
  
  stats <- list(
    media = mean(I_MPCS_sim, na.rm = TRUE),
    mediana = median(I_MPCS_sim, na.rm = TRUE),
    sd = sd(I_MPCS_sim, na.rm = TRUE),
    cv = sd(I_MPCS_sim, na.rm = TRUE) / mean(I_MPCS_sim, na.rm = TRUE) * 100,
    min = min(I_MPCS_sim, na.rm = TRUE),
    max = max(I_MPCS_sim, na.rm = TRUE),
    q25 = quantile(I_MPCS_sim, 0.25, na.rm = TRUE),
    q75 = quantile(I_MPCS_sim, 0.75, na.rm = TRUE)
  )
  
  ic_inf <- quantile(I_MPCS_sim, 0.025, na.rm = TRUE)
  ic_sup <- quantile(I_MPCS_sim, 0.975, na.rm = TRUE)
  
  list(
    I_MPCS_sim = I_MPCS_sim,
    tipo_sim = tipo_sim,
    pesos_sim = pesos_sim,
    I_MPCS_default = I_MPCS_default,
    k_default = k_default,
    tipo_default = tipo_default,
    freq_tipo = freq_tipo,
    tipo_mas_frecuente = tipo_mas_frecuente,
    pct_mas_frecuente = pct_mas_frecuente,
    pct_coincidencia = pct_coincidencia,
    stats = stats,
    ic_inf = ic_inf,
    ic_sup = ic_sup,
    n_sim = n_sim,
    R = R,
    escala = escala
  )
}

# ============================================================================
# 11. grafico_sensibilidad — Gráfico de sensibilidad
# ============================================================================

#' Gráfico de sensibilidad
#'
#' @param resultado_sens Resultado de analisis_sensibilidad()
#' @param idioma "es" o "en"
#' @return Objeto ggplot
#' @export
grafico_sensibilidad <- function(resultado_sens, idioma = "es") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Se requiere el paquete ggplot2")
  }
  
  df <- data.frame(
    I_MPCS = resultado_sens$I_MPCS_sim,
    Tipo = resultado_sens$tipo_sim
  )
  
  if (idioma == "es") {
    titulo <- "Análisis de Sensibilidad del Índice MPCS"
    subtitulo <- paste("10,000 combinaciones aleatorias de pesos w1, w2, w3")
    eje_x <- "Índice MPCS"
    eje_y <- "Frecuencia"
    leyenda <- "Tipo de nudge"
    default_label <- "Pesos por defecto"
  } else {
    titulo <- "Sensitivity Analysis of MPCS Index"
    subtitulo <- paste("10,000 random combinations of weights w1, w2, w3")
    eje_x <- "MPCS Index"
    eje_y <- "Frequency"
    leyenda <- "Nudge type"
    default_label <- "Default weights"
  }
  
  colores <- c(
    "Informational" = "#74B3CE",
    "Structural" = "#2E86AB",
    "Normative" = "#E84855",
    "Systemic multi-nudge" = "#1A3A5C"
  )
  
  tipos_presentes <- unique(df$Tipo)
  colores_filtrados <- colores[names(colores) %in% tipos_presentes]
  
  bins <- 50
  hist_data <- hist(df$I_MPCS, breaks = bins, plot = FALSE)
  max_y <- max(hist_data$counts)
  
  ggplot2::ggplot(df, ggplot2::aes(x = I_MPCS, fill = Tipo)) +
    ggplot2::geom_histogram(bins = bins, color = "white", alpha = 0.85) +
    ggplot2::geom_vline(xintercept = resultado_sens$I_MPCS_default,
               color = "#C0392B", linetype = "dashed", linewidth = 1.2) +
    ggplot2::geom_vline(xintercept = resultado_sens$ic_inf,
               color = "#2C3E50", linetype = "dotted", linewidth = 0.8, alpha = 0.5) +
    ggplot2::geom_vline(xintercept = resultado_sens$ic_sup,
               color = "#2C3E50", linetype = "dotted", linewidth = 0.8, alpha = 0.5) +
    ggplot2::annotate("text", 
             x = resultado_sens$I_MPCS_default + 0.008, 
             y = max_y * 0.9,
             label = paste0(default_label, "\nI_MPCS = ", 
                           round(resultado_sens$I_MPCS_default, 4)),
             hjust = 0, size = 3.5, color = "#C0392B") +
    ggplot2::scale_fill_manual(values = colores_filtrados) +
    ggplot2::labs(
      title = titulo,
      subtitle = subtitulo,
      x = eje_x,
      y = eje_y,
      fill = leyenda
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# ============================================================================
# FIN DEL ARCHIVO
# ============================================================================
