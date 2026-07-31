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
#' @param umbral valor minimo de correlacion para incluir aristas (defecto: 0.10)
#' @return lista con grafo, centralidad, nodo optimo e indice de impacto
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
    warning("Se recomiendan al menos 5 variables para un grafo estable. Se usaran ", length(variables), " variables.")
  }
  
  df <- datos[, variables, drop = FALSE]
  df <- df[, colSums(is.na(df)) < nrow(df) * 0.5, drop = FALSE]
  
  if (ncol(df) < 3) {
    warning("Menos de 3 variables validas despues de limpiar datos faltantes")
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
    warning("Error al calcular la matriz de correlacion: ", e$message)
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
#' @param tasas_avance vector numérico de longitud m-1 con probabilidades de transición hacia adelante.
#'        Si es NULL, usa heurística genérica. Si se provee, construye la matriz con esas tasas.
#' @param umbral_objetivo proporción para considerar convergencia (defecto: 0.50)
#' @param horizonte número de períodos para proyectar (defecto: 10)
#' @return lista con matriz P, simulación, índice de Markov y tiempo de referencia
#' @export
calcular_markov <- function(estados, orden_estados = NULL, 
                            tasas_avance = NULL,
                            umbral_objetivo = 0.50, horizonte = 10) {
  
  # --- Validaciones ---
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
  
  # Limpiar datos
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
  
  # --- Determinar estados únicos ---
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
  
  # --- Calcular distribución actual ---
  freq <- table(estados_clean)
  dist_actual <- rep(0, m)
  names(dist_actual) <- orden_estados
  
  for (i in seq_along(orden_estados)) {
    if (orden_estados[i] %in% names(freq)) {
      dist_actual[i] <- freq[orden_estados[i]] / sum(freq)
    }
  }
  
  # --- Construir matriz de transición ---
  P <- matrix(0, nrow = m, ncol = m)
  colnames(P) <- orden_estados
  rownames(P) <- orden_estados
  
  # ============================================================
  # [NUEVO] Si se proveen tasas_avance, usarlas para construir la matriz
  # ============================================================
  if (!is.null(tasas_avance) && length(tasas_avance) == (m - 1)) {
    
    # --- Construcción con tasas_avance definidas por el usuario ---
    for (i in 1:(m-1)) {
      t <- tasas_avance[i]
      
      # Asegurar que t esté en [0,1]
      t <- max(0, min(1, t))
      
      # Inercia base: se reduce con el avance
      inercia <- 0.50 - 0.15 * (i/m)
      inercia <- max(0.10, min(0.80, inercia))
      
      # Retroceso base: se reduce con el avance
      retroceso <- 0.05 * (1 - i/m)
      retroceso <- max(0.01, min(0.20, retroceso))
      
      # ============================================================
      # [CORREGIDO] La tasa de avance se respeta EXACTAMENTE
      # Solo se ajustan inercia y retroceso para que quepan en 1 - t
      # ============================================================
      resto <- 1 - t
      
      # Escalar inercia para que quepa en el remanente
      if (inercia > resto) {
        inercia <- resto * 0.85  # deja 15% para retroceso
      }
      
      # Escalar retroceso para que quepa en el remanente
      if (retroceso > (resto - inercia)) {
        retroceso <- resto - inercia
      }
      retroceso <- max(0, retroceso)
      
      # Si sobra espacio, distribuir en inercia y retroceso proporcionalmente
      sobrante <- resto - inercia - retroceso
      if (sobrante > 0.001) {
        inercia <- inercia + sobrante * 0.7
        retroceso <- retroceso + sobrante * 0.3
      }
      
      # Asignar a la matriz
      P[i, i] <- inercia
      P[i, i+1] <- t
      
      if (i > 1) {
        P[i, i-1] <- retroceso
      }
      
      # Distribuir el resto en otros estados (si existe)
      resto_final <- 1 - sum(P[i, ])
      if (resto_final > 0.001) {
        otros <- setdiff(1:m, c(i, i+1, if (i > 1) i-1 else NULL))
        if (length(otros) > 0) {
          P[i, otros] <- resto_final / length(otros)
        }
      }
    }
    
    # Último estado: alta permanencia
    P[m, m] <- 0.85
    P[m, 1:(m-1)] <- (1 - 0.85) / (m - 1)
    
  } else {
    
    # ============================================================
    # [FALLBACK] Heurística genérica (cuando no se proveen tasas)
    # ============================================================
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
  }
  
  # --- Normalizar filas para asegurar que suman 1 ---
  P <- P / rowSums(P)
  
  # --- Simular cadena de Markov ---
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
  
  # ============================================================
  # [OPCIÓN B] I_Markov = adherencia en horizonte fijo H
  # ============================================================
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
# 3. calcular_juegos — Calcula masa critica (TEORIA DE JUEGOS DINAMICA CON α)
# ============================================================================
 
#' Calcular teoria de juegos y masa critica (dinamica con α)
#'
#' @param P matriz de transicion de Markov (opcional, no se usa)
#' @param R_factor factor de recursos (defecto: 0.65) — OBSOLETO
#' @param alpha indicador de acceso a salud (0-1). Si no se proporciona, usa 0.60 por defecto.
#' @return lista con masa critica e indice
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
# 4. calcular_indice — Calcula Indice MPCS y tipo de nudge
# ============================================================================

#' Calcular Indice MPCS y tipo de nudge
#'
#' @param I_grafo indice del modulo de grafos
#' @param I_markov indice del modulo de Markov (Opcion B: H=10)
#' @param I_juegos indice del modulo de juegos (dinamico con α)
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
    warning("Los ponderadores no suman 1. Se normalizaran.")
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
# 5. aplicar_nudge — Aplica nudge a la matriz de transicion
# ============================================================================

#' Aplicar nudge a la matriz de transicion de Markov
#'
#' @param P matriz de transicion original
#' @param k intensidad del nudge (0-1)
#' @return matriz de transicion modificada
#' @export
aplicar_nudge <- function(P, k) {
  
  if (is.null(P) || !is.matrix(P)) {
    warning("Se requiere una matriz de transicion valida")
    return(NULL)
  }
  
  if (k < 0 || k > 1) {
    warning("k debe estar entre 0 y 1. Se usara k = 0.4")
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
# 6. generar_demo_data — Genera datos de salud (demostracion)
# ============================================================================

#' Generar datos de demostracion (salud - ENDES demo)
#'
#' @param n numero de observaciones
#' @param seed semilla para reproducibilidad
#' @return data.frame con datos de demostracion
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
# 7. generar_datos_aula — Genera datos de educacion (demostracion)
# ============================================================================

#' Generar datos de demostracion (educacion - participacion en aula)
#'
#' @param n numero de estudiantes
#' @param seed semilla para reproducibilidad
#' @return data.frame con datos de demostracion
#' @export
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
  
  grupo <- sample(c("Seccion A", "Seccion B", "Seccion C"), n, replace = TRUE, prob = c(0.40, 0.35, 0.25))
  
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

# ============================================================================
# graficar_arbol_markov — Arbol de transicion de Markov (VERSION SIMPLE Y CLARA)
# ============================================================================

#' Graficar arbol de transicion de Markov (Version simple y clara)
#'
#' @param P matriz de transicion (m x m)
#' @param estados vector con nombres de estados
#' @param umbral_prob probabilidad minima para mostrar una flecha (default: 0.05)
#' @param titulo titulo del grafico
#' @return objeto ggplot
#' @export
graficar_arbol_markov <- function(P, estados = NULL, umbral_prob = 0.05, 
                                  titulo = "Markov Transitions") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Se requiere el paquete ggplot2")
  }
  
  if (is.null(P) || nrow(P) < 2) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + 
             ggplot2::annotate("text", x = 0.5, y = 0.5, 
                              label = "No data for Markov tree",
                              size = 6, color = "#7F8C8D", fontface = "bold"))
  }
  
  m <- nrow(P)
  
  if (!is.null(estados) && length(estados) == m) {
    estados_nombres <- estados
  } else if (!is.null(colnames(P)) && all(colnames(P) != "")) {
    estados_nombres <- colnames(P)
  } else {
    estados_nombres <- paste0("E", 1:m)
  }
  
  nodos <- data.frame(
    id = 1:m,
    nombre = estados_nombres,
    x = 1:m,
    y = 0
  )
  
  transiciones <- data.frame()
  for (i in 1:m) {
    for (j in 1:m) {
      if (P[i, j] >= umbral_prob && i != j) {
        offset_y <- if (j > i) 0.15 else -0.15
        transiciones <- rbind(transiciones, data.frame(
          from = i, to = j,
          from_x = i, to_x = j,
          from_y = offset_y, to_y = offset_y,
          prob = round(P[i, j] * 100, 1),
          tipo = if (j > i) "Advance" else "Regression",
          color = if (j > i) "#2E86AB" else "#E84855"
        ))
      }
    }
  }
  
  if (nrow(transiciones) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + 
             ggplot2::annotate("text", x = 0.5, y = 0.5, 
                              label = paste0("No transitions >= ", round(umbral_prob*100, 1), "%"),
                              size = 6, color = "#7F8C8D", fontface = "bold"))
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = transiciones,
      aes(x = from_x, xend = to_x, 
          y = from_y, yend = to_y,
          color = tipo,
          linewidth = prob),
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
      alpha = 0.8
    ) +
    ggplot2::scale_linewidth_continuous(range = c(0.5, 2.5), guide = "none") +
    ggplot2::scale_color_manual(
      values = c("Advance" = "#2E86AB", "Regression" = "#E84855"),
      name = "Type",
      labels = c("Advance" = "Advance \u2192", "Regression" = "Regression \u2190")
    ) +
    ggplot2::geom_label(
      data = transiciones,
      aes(x = (from_x + to_x) / 2, 
          y = (from_y + to_y) / 2 + 0.12,
          label = paste0(prob, "%")),
      size = 4.5, fontface = "bold", fill = "white",
      color = "#1A3A5C", label.size = 0.3,
      label.padding = unit(0.2, "lines"), alpha = 0.95,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = nodos, aes(x = x, y = y),
      size = 22, color = "#2C3E50", fill = "#F8F9F9",
      shape = 21, stroke = 2
    ) +
    ggplot2::geom_text(
      data = nodos, aes(x = x, y = y, label = nombre),
      size = 4, fontface = "bold", color = "#1A3A5C"
    ) +
    ggplot2::xlim(0.5, m + 0.5) + ggplot2::ylim(-0.6, 0.6) +
    ggplot2::labs(
      title = titulo,
      subtitle = paste0("Only transitions >= ", round(umbral_prob * 100, 1), "%"),
      x = "", y = "",
      color = ""
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "#1A3A5C"),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#5D6D7E"),
      plot.margin = margin(20, 30, 20, 30),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10)
    )
  
  return(p)
}

# ============================================================================
# 9. graficar_juego_evolutivo — Dinamica replicadora (teoria de juegos)
# ============================================================================

#' Graficar dinamica replicadora (teoria de juegos)
#'
#' @param alpha indicador de contexto (0-1)
#' @param p_star masa critica (calculada automaticamente si no se proporciona)
#' @param titulo titulo del grafico
#' @return objeto ggplot
#' @export
graficar_juego_evolutivo <- function(alpha = 0.60, p_star = NULL, 
                                     titulo = "Replicator Dynamics") {
  
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
    geom_line(aes(y = f_A, color = "Adopters (A)"), linewidth = 1.2) +
    geom_line(aes(y = f_R, color = "Resisters (R)"), linewidth = 1.2) +
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
             label = "Unstable", color = "#C0392B", size = 3.5, hjust = 0.5) +
    annotate("text", x = (1 + p_star) / 2, y = max(df$f_A) * 0.95,
             label = "Stable", color = "#27AE60", size = 3.5, hjust = 0.5) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_y_continuous(labels = scales::number) +
    scale_color_manual(
      name = "Strategy",
      values = c("Adopters (A)" = "#2E86AB", "Resisters (R)" = "#E84855")
    ) +
    labs(
      title = titulo,
      subtitle = paste0("α = ", round(alpha, 3), " | Critical mass (p*) = ", round(p_star, 3)),
      x = "Proportion of adopters (p)",
      y = "Expected payoff",
      color = "Strategy"
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
      title = "Replicator Dynamics (dp/dt)",
      subtitle = "Rate of change in the proportion of adopters",
      x = "Proportion of adopters (p)",
      y = "dp/dt (rate of change)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#7F8C8D")
    )
  
  p_combinado <- p1 + p2 + 
    patchwork::plot_annotation(
      title = paste0("Game Theory Analysis - MPCS"),
      subtitle = paste0("α = ", round(alpha, 3), " | p* = ", round(p_star, 3)),
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#7F8C8D")
      )
    )
  
  return(p_combinado)
}

# ============================================================================
# 10. graficar_trayectorias_markov — Trayectorias de Markov mejoradas
# ============================================================================

#' Graficar trayectorias de Markov mejoradas
#'
#' @param sim_base simulacion base
#' @param sim_nudge simulacion con nudge (opcional)
#' @param estados nombres de los estados
#' @param titulo titulo del grafico
#' @param y_label etiqueta del eje Y
#' @return objeto ggplot
#' @export
graficar_trayectorias_markov <- function(sim_base, sim_nudge = NULL, 
                                         estados = NULL, 
                                         titulo = "Markov Trajectories",
                                         y_label = "Occupation probability") {
  
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
                     label = "No data for trajectories",
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
  df_base$Escenario <- "Without nudge"
  
  if (!is.null(sim_nudge) && nrow(sim_nudge) == nrow(sim_base)) {
    df_nudge <- as.data.frame(sim_nudge)
    if (ncol(df_nudge) == length(estados)) {
      colnames(df_nudge) <- estados
      df_nudge$Periodo <- 0:(nrow(df_nudge) - 1)
      df_nudge$Escenario <- "With nudge"
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
    scale_linetype_manual(values = c("Without nudge" = "solid", "With nudge" = "dashed")) +
    labs(
      title = titulo,
      subtitle = "Expected proportion of individuals in each state over time",
      x = "Period (time units)",
      y = y_label,
      color = "Behavioral state",
      linetype = "Scenario"
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
# 11. format_report — Formatea resultados para reporte
# ============================================================================

#' Formatear resultados para reporte
#'
#' @param results data.frame con resultados del MPCS
#' @param plots lista de graficos generados (opcional)
#' @return lista formateada para reporte
#' @export
format_report <- function(results, plots = NULL) {
  
  if (missing(results) || is.null(results)) {
    warning("Se requieren resultados")
    return(NULL)
  }
  
  report <- list(
    resumen = results,
    timestamp = Sys.time(),
    n_grupos = nrow(results),
    n_total = sum(results$n),
    max_mpcs = max(results$I_MPCS),
    min_mpcs = min(results$I_MPCS),
    grupo_prioritario = results$Grupo[which.max(results$I_MPCS)],
    nodo_prioritario = results$Nodo_Optimo[which.max(results$I_MPCS)],
    tipo_prioritario = results$Tipo_Nudge[which.max(results$I_MPCS)]
  )
  
  if (!is.null(plots)) {
    report$plots <- plots
  }
  
  report$ranking <- results %>%
    dplyr::arrange(dplyr::desc(I_MPCS)) %>%
    dplyr::mutate(
      Prioridad = 1:n(),
      Nivel = dplyr::case_when(
        I_MPCS >= 0.75 ~ "High",
        I_MPCS >= 0.50 ~ "Medium",
        TRUE ~ "Low"
      )
    )
  
  return(report)
}

# ============================================================================
# 12. validar_datos — Valida la estructura de los datos de entrada
# ============================================================================

#' Validar la estructura de los datos de entrada
#'
#' @param datos data.frame a validar
#' @param min_filas numero minimo de filas (defecto: 30)
#' @param min_vars_num numero minimo de variables numericas (defecto: 5)
#' @return lista con resultado de validacion
#' @export
validar_datos <- function(datos, min_filas = 30, min_vars_num = 5) {
  
  if (missing(datos) || is.null(datos)) {
    return(list(
      valido = FALSE,
      errores = c("No se proporcionaron datos")
    ))
  }
  
  errores <- c()
  
  if (nrow(datos) < min_filas) {
    errores <- c(errores, paste("Se necesitan al menos", min_filas, "filas. Actual:", nrow(datos)))
  }
  
  vars_num <- names(datos)[sapply(datos, is.numeric)]
  if (length(vars_num) < min_vars_num) {
    errores <- c(errores, paste("Se necesitan al menos", min_vars_num, 
                                "variables numericas. Actual:", length(vars_num)))
  }
  
  total_na <- sum(is.na(datos))
  if (total_na > 0) {
    pct_na <- round(total_na / (nrow(datos) * ncol(datos)) * 100, 1)
    if (pct_na > 20) {
      errores <- c(errores, paste("Alto porcentaje de valores faltantes:", pct_na, "%"))
    }
  }
  
  return(list(
    valido = length(errores) == 0,
    errores = errores,
    n_filas = nrow(datos),
    n_cols = ncol(datos),
    vars_num = vars_num,
    pct_na = if (total_na > 0) round(total_na / (nrow(datos) * ncol(datos)) * 100, 1) else 0
  ))
}

# ============================================================================
# 13. analisis_sensibilidad — Analisis de sensibilidad de pesos
# ============================================================================

#' Analisis de sensibilidad del Indice MPCS
#'
#' @param I_Grafo Indice del modulo de grafos (0-1)
#' @param I_Markov Indice del modulo de Markov (0-1)
#' @param I_Juegos Indice del modulo de juegos (0-1)
#' @param n_sim Numero de simulaciones (defecto: 10000)
#' @param R Factor de recursos (defecto: 0.65)
#' @param escala Factor de escala (defecto: 1.5)
#' @param seed Semilla para reproducibilidad (defecto: 123)
#' @return Lista con resultados del analisis de sensibilidad
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
# 14. grafico_sensibilidad — Grafico de sensibilidad
# ============================================================================

#' Grafico de sensibilidad
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
    titulo <- "Analisis de Sensibilidad del Indice MPCS"
    subtitulo <- paste("10,000 combinaciones aleatorias de pesos w1, w2, w3")
    eje_x <- "Indice MPCS"
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
# 15. resumen_sensibilidad — Resumen textual de sensibilidad
# ============================================================================

#' Resumen textual de sensibilidad
#'
#' @param resultado_sens Resultado de analisis_sensibilidad()
#' @param idioma "es" o "en"
#' @return Texto formateado (impreso en consola)
#' @export
resumen_sensibilidad <- function(resultado_sens, idioma = "es") {
  
  if (idioma == "es") {
    cat("=== RESULTADOS DEL ANALISIS DE SENSIBILIDAD ===\n\n")
    cat(sprintf("Numero de simulaciones: %d\n", resultado_sens$n_sim))
    cat(sprintf("Media de I_MPCS: %.4f\n", resultado_sens$stats$media))
    cat(sprintf("Desviacion estandar: %.4f\n", resultado_sens$stats$sd))
    cat(sprintf("Coeficiente de variacion: %.1f%%\n", resultado_sens$stats$cv))
    cat(sprintf("Intervalo de confianza 95%%: [%.4f, %.4f]\n", 
                resultado_sens$ic_inf, resultado_sens$ic_sup))
    cat(sprintf("\nI_MPCS con pesos por defecto (w1=0.35, w2=0.40, w3=0.25): %.4f\n", 
                resultado_sens$I_MPCS_default))
    cat(sprintf("Tipo de nudge por defecto: %s\n", resultado_sens$tipo_default))
    cat(sprintf("\nTipo de nudge mas frecuente en simulaciones: %s (%.1f%%)\n",
                resultado_sens$tipo_mas_frecuente, resultado_sens$pct_mas_frecuente))
    cat(sprintf("Coincidencia con tipo por defecto: %.1f%%\n", 
                resultado_sens$pct_coincidencia))
    
    cat("\n=== INTERPRETACION ===\n")
    if (resultado_sens$pct_coincidencia >= 70) {
      cat("La recomendacion de nudge es ALTAMENTE ROBUSTA.\n")
      cat("  El tipo de nudge se mantiene en mas del 70% de las simulaciones.\n")
    } else if (resultado_sens$pct_coincidencia >= 50) {
      cat("La recomendacion de nudge es ROBUSTA.\n")
      cat("  El tipo de nudge se mantiene en mas del 50% de las simulaciones.\n")
    } else {
      cat("La recomendacion de nudge es SENSIBLE a los pesos.\n")
      cat("  Se recomienda revisar la ponderacion o realizar un analisis adicional.\n")
    }
    
  } else {
    cat("=== SENSITIVITY ANALYSIS RESULTS ===\n\n")
    cat(sprintf("Number of simulations: %d\n", resultado_sens$n_sim))
    cat(sprintf("Mean I_MPCS: %.4f\n", resultado_sens$stats$media))
    cat(sprintf("Standard deviation: %.4f\n", resultado_sens$stats$sd))
    cat(sprintf("Coefficient of variation: %.1f%%\n", resultado_sens$stats$cv))
    cat(sprintf("95%% Confidence interval: [%.4f, %.4f]\n", 
                resultado_sens$ic_inf, resultado_sens$ic_sup))
    cat(sprintf("\nI_MPCS with default weights (w1=0.35, w2=0.40, w3=0.25): %.4f\n", 
                resultado_sens$I_MPCS_default))
    cat(sprintf("Default nudge type: %s\n", resultado_sens$tipo_default))
    cat(sprintf("\nMost frequent nudge type in simulations: %s (%.1f%%)\n",
                resultado_sens$tipo_mas_frecuente, resultado_sens$pct_mas_frecuente))
    cat(sprintf("Agreement with default type: %.1f%%\n", 
                resultado_sens$pct_coincidencia))
    
    cat("\n=== INTERPRETATION ===\n")
    if (resultado_sens$pct_coincidencia >= 70) {
      cat("The nudge recommendation is HIGHLY ROBUST.\n")
      cat("  The nudge type remains in more than 70% of simulations.\n")
    } else if (resultado_sens$pct_coincidencia >= 50) {
      cat("The nudge recommendation is ROBUST.\n")
      cat("  The nudge type remains in more than 50% of simulations.\n")
    } else {
      cat("The nudge recommendation is SENSITIVE to weights.\n")
      cat("  Consider reviewing the weighting or conducting additional analysis.\n")
    }
  }
}

# ============================================================================
# 16. seleccionar_umbral_bootstrap — Selección de umbral por estabilidad
# ============================================================================
# ============================================================================
# 16. seleccionar_umbral_bootstrap — Selección de umbral por estabilidad
# ============================================================================

#' Seleccionar umbral de correlación por bootstrap stability (Jaccard similarity)
#'
#' @param datos data.frame con los datos
#' @param variables vector con nombres de columnas a incluir en el grafo
#' @param umbrales vector con umbrales a evaluar (defecto: c(0.05, 0.07, 0.10, 0.12, 0.15, 0.20))
#' @param n_boot número de remuestreos bootstrap (defecto: 20)
#' @param seed semilla para reproducibilidad (defecto: 123)
#' @param criterio_jaccard umbral mínimo de Jaccard para considerar estable (defecto: 0.70)
#' @return lista con umbral óptimo, estabilidad por umbral, y nodos óptimos
#' @export
seleccionar_umbral_bootstrap <- function(datos, variables, 
                                         umbrales = c(0.05, 0.07, 0.10, 0.12, 0.15, 0.20),
                                         n_boot = 20,
                                         seed = 123,
                                         criterio_jaccard = 0.70) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Se requiere el paquete igraph")
  }
  
  set.seed(seed)
  n <- nrow(datos)
  n_umbrales <- length(umbrales)
  
  resultados <- data.frame()
  
  # --- Bucle por cada umbral ---
  for (idx_u in seq_along(umbrales)) {
    umbral <- umbrales[idx_u]
    
    # --- 1. Grafo con datos originales ---
    res_original <- calcular_grafo(datos, variables, umbral = umbral)
    nodo_original <- res_original$optimal_node
    if (is.na(nodo_original) || is.null(nodo_original)) nodo_original <- "Sin_nodo"
    
    # --- 2. Bootstrap ---
    nodos_boot <- character()
    aristas_boot <- list()
    
    for (b in 1:n_boot) {
      # Remuestrear filas con reemplazo
      idx_boot <- sample(1:n, size = n, replace = TRUE)
      datos_boot <- datos[idx_boot, , drop = FALSE]
      
      # Construir grafo con el umbral actual
      res_boot <- calcular_grafo(datos_boot, variables, umbral = umbral)
      nodo_boot <- res_boot$optimal_node
      if (is.na(nodo_boot) || is.null(nodo_boot)) nodo_boot <- "Sin_nodo"
      nodos_boot <- c(nodos_boot, nodo_boot)
      
      # Guardar aristas si hay grafo
      if (!is.null(res_boot$graph)) {
        aristas_boot[[b]] <- igraph::as_edgelist(res_boot$graph)
      } else {
        aristas_boot[[b]] <- matrix(0, nrow = 0, ncol = 2)
      }
    }
    
    # --- 3. Calcular Jaccard para aristas ---
    jaccards <- numeric()
    for (b in 1:(n_boot - 1)) {
      for (c in (b + 1):n_boot) {
        E_b <- aristas_boot[[b]]
        E_c <- aristas_boot[[c]]
        
        if (nrow(E_b) == 0 && nrow(E_c) == 0) {
          jaccard <- 1
        } else if (nrow(E_b) == 0 || nrow(E_c) == 0) {
          jaccard <- 0
        } else {
          set_b <- apply(E_b, 1, function(x) paste(sort(x), collapse = "_"))
          set_c <- apply(E_c, 1, function(x) paste(sort(x), collapse = "_"))
          
          interseccion <- length(intersect(set_b, set_c))
          union <- length(union(set_b, set_c))
          jaccard <- interseccion / union
        }
        jaccards <- c(jaccards, jaccard)
      }
    }
    
    jaccard_promedio <- mean(jaccards, na.rm = TRUE)
    
    # --- 4. Estabilidad del nodo óptimo ---
    freq_nodos <- table(nodos_boot)
    nodo_mas_frecuente <- names(freq_nodos)[which.max(freq_nodos)]
    pct_nodo_mas_frecuente <- max(freq_nodos) / n_boot * 100
    
    # --- 5. Guardar resultados ---
    resultados <- rbind(resultados, data.frame(
      Umbral = umbral,
      Nodo_Original = nodo_original,
      Nodo_Mas_Frecuente = nodo_mas_frecuente,
      Pct_Nodo_Estable = round(pct_nodo_mas_frecuente, 1),
      Jaccard_Promedio = round(jaccard_promedio, 4),
      Cumple_Criterio = jaccard_promedio >= criterio_jaccard & pct_nodo_mas_frecuente >= 70
    ))
  }
  
  # --- 6. Seleccionar umbral óptimo ---
  resultados_filtrados <- resultados[resultados$Cumple_Criterio, ]
  
  if (nrow(resultados_filtrados) > 0) {
    idx_optimo <- which.max(resultados_filtrados$Jaccard_Promedio + resultados_filtrados$Pct_Nodo_Estable / 100)
    umbral_optimo <- resultados_filtrados$Umbral[idx_optimo]
  } else {
    idx_optimo <- which.max(resultados$Jaccard_Promedio)
    umbral_optimo <- resultados$Umbral[idx_optimo]
    warning("Ningún umbral cumple con el criterio de estabilidad. Seleccionado el de mayor Jaccard.")
  }
  
  return(list(
    umbral_optimo = umbral_optimo,
    resultados = resultados,
    criterio_jaccard = criterio_jaccard,
    n_boot = n_boot,
    nodo_original = resultados$Nodo_Original[which(resultados$Umbral == umbral_optimo)[1]],
    nodo_mas_frecuente = resultados$Nodo_Mas_Frecuente[which(resultados$Umbral == umbral_optimo)[1]]
  ))
}


# ============================================================================
# FIN DEL ARCHIVO
# ============================================================================
