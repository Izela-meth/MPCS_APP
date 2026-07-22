# ============================================================================
# mpcs_functions.R — Funciones modulares del MPCS
# ============================================================================
# Este archivo contiene todas las funciones matemáticas del Modelo Predictivo
# de Cambio Conductual por Sistemas (MPCS)
# ============================================================================
# Funciones incluidas:
#   1. calcular_grafo()     - Construye el grafo y calcula centralidades
#   2. calcular_markov()    - Estima cadena de Markov y convergencia
#   3. calcular_juegos()    - Calcula masa crítica (teoría de juegos)
#   4. calcular_indice()    - Calcula Índice MPCS y tipo de nudge
#   5. aplicar_nudge()      - Aplica nudge a matriz de transición
#   6. generar_demo_data()  - Genera datos de demostración
#   7. format_report()      - Formatea resultados para reporte
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
  
  # --- Validaciones ---
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
  
  # --- Seleccionar y limpiar datos ---
  df <- datos[, variables, drop = FALSE]
  
  # Eliminar columnas con > 50% de valores faltantes
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
  
  # --- Matriz de correlación ---
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
  
  # --- Crear aristas ---
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
  
  # --- Crear grafo ---
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
  
  # --- Calcular centralidades ---
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
  
  # --- Retornar resultados ---
  return(list(
    graph = g,
    centralidad = centr,
    optimal_node = centr$Variable[1],
    score = centr$Impacto[1]
  ))
}

# ============================================================================
# 2. calcular_markov — Estima cadena de Markov y tiempo de convergencia
# ============================================================================

#' Calcular cadena de Markov y tiempo de convergencia
#'
#' @param estados vector con los estados de cada individuo
#' @param orden_estados vector con el orden progresivo de estados (opcional)
#' @param umbral_objetivo proporción para considerar convergencia (defecto: 0.50)
#' @return lista con matriz P, simulación, tiempo de convergencia e índice
#' @export
calcular_markov <- function(estados, orden_estados = NULL, umbral_objetivo = 0.50) {
  
  # --- Validaciones ---
  if (missing(estados)) {
    warning("Se requiere el vector de estados")
    return(list(
      mat = NULL,
      dist_actual = NULL,
      sim_base = NULL,
      T_base = 10,
      score = 0.3
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
      T_base = 10,
      score = 0.3
    ))
  }
  
  # --- Determinar estados únicos ---
  estados_unicos <- unique(estados_clean)
  
  # Si no se especifica orden, intentar ordenar automáticamente
  if (is.null(orden_estados)) {
    # Intentar ordenar como E1, E2, E3, ...
    if (all(grepl("^E[0-9]+$", estados_unicos))) {
      orden_estados <- estados_unicos[order(as.numeric(gsub("E", "", estados_unicos)))]
    } else {
      # Ordenar por frecuencia (de más a menos común)
      freq <- table(estados_clean)
      orden_estados <- names(sort(freq, decreasing = TRUE))
    }
  } else {
    # Verificar que todos los estados estén en el orden
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
      T_base = 10,
      score = 0.3
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
  # En ausencia de datos longitudinales, se usan supuestos de progresión conservadores
  P <- matrix(0, nrow = m, ncol = m)
  colnames(P) <- orden_estados
  rownames(P) <- orden_estados
  
  for (i in 1:(m-1)) {
    # Probabilidad de avanzar al siguiente estado
    # Mayor probabilidad de avance desde estados tempranos
    prob_avance <- 0.30 + 0.20 * (1 - i/m)
    
    # Probabilidad de permanecer
    prob_quedarse <- 0.50 - 0.15 * (i/m)
    
    # Probabilidad de retroceder (pequeña)
    prob_retroceso <- 0.10 * (1 - i/m)
    
    P[i, i] <- prob_quedarse
    P[i, i+1] <- prob_avance
    
    if (i > 1) {
      P[i, i-1] <- prob_retroceso
    }
    
    # Distribuir el resto entre otros estados
    resto <- 1 - sum(P[i, ])
    if (resto > 0) {
      otros <- setdiff(1:m, c(i, i+1, if (i > 1) i-1 else NULL))
      if (length(otros) > 0) {
        P[i, otros] <- resto / length(otros)
      }
    }
  }
  
  # Último estado: alta permanencia (la conducta está consolidada)
  P[m, m] <- 0.85
  P[m, 1:(m-1)] <- (1 - 0.85) / (m - 1)
  
  # Normalizar filas para asegurar que suman 1
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
  
  # --- Calcular tiempo de convergencia ---
  # Se considera "convergencia" cuando los últimos dos estados alcanzan el umbral
  if (m >= 2) {
    objetivo <- sim_base[, m] + sim_base[, max(1, m-1)]
    T_base <- which(objetivo >= umbral_objetivo)[1]
    
    if (is.na(T_base) || is.infinite(T_base)) {
      T_base <- 20
    }
  } else {
    T_base <- 10
  }
  
  # --- Índice Markov (urgencia temporal) ---
  # Menor T = mayor urgencia = mayor I_Markov
  score <- 1 - exp(-T_base / 20)
  
  # --- Retornar resultados ---
  return(list(
    mat = P,
    dist_actual = dist_actual,
    sim_base = sim_base,
    T_base = T_base,
    score = score
  ))
}

# ============================================================================
# 3. calcular_juegos — Calcula masa crítica (teoría de juegos evolutiva)
# ============================================================================
 
#' Calcular teoría de juegos y masa crítica
#'
#' @param P matriz de transición de Markov (opcional)
#' @param R_factor factor de recursos (defecto: 0.65)
#' @return lista con masa crítica e índice
#' @export
calcular_juegos <- function(P = NULL, R_factor = 0.65) {
  
  # --- Matriz de pagos (fija por diseño del modelo) ---
  # Estos valores representan la interacción social entre adoptantes (A) y resistentes (R)
  a_AA <- 2.0   # Refuerzo mutuo entre adoptantes
  a_AR <- -0.5  # Presión social negativa sobre el adoptante aislado
  a_RA <- 1.0   # El resistente se beneficia del ejemplo ajeno
  a_RR <- 0.5   # Refuerzo mutuo de la inacción
  
  # --- Calcular masa crítica (p*) ---
  # p* es el umbral mínimo de adoptantes para que la conducta sea auto-sostenible.
  # Ecuación (Sección 3.5.3 / 2.4 del artículo): p* = (a_RR - a_AR) / (a_AA - a_AR - a_RA + a_RR)
  p_star <- (a_RR - a_AR) / (a_AA - a_AR - a_RA + a_RR)
  p_star <- max(0, min(1, p_star))
  
  # NOTA IMPORTANTE — R_factor NO se aplica aquí.
  # El artículo (Sección 3.5.4 / 2.5) especifica un único punto donde entra el
  # factor de disponibilidad de recursos: k = min(1, I_MPCS * R_factor * 1.5),
  # calculado en calcular_indice(). Aplicarlo también dentro de p* duplicaba
  # su efecto (doble conteo) e inflaba I_Juegos y, en cascada, I_MPCS.
  # R_factor se mantiene como argumento por compatibilidad de firma con las
  # llamadas existentes en server(), pero no se usa en este cálculo.
  
  # --- Índice de juegos ---
  # I_Juegos = p*, sin inversión ni transformación adicional (artículo, 3.5.3):
  # "Este valor de p* ingresa sin transformación adicional como I_Juegos".
  score <- p_star
  
  return(list(
    p_star = p_star,
    score = score
  ))
}

# ============================================================================
# 4. calcular_indice — Calcula Índice MPCS y tipo de nudge
# ============================================================================

#' Calcular Índice MPCS y tipo de nudge
#'
#' @param I_grafo índice del módulo de grafos
#' @param I_markov índice del módulo de Markov
#' @param I_juegos índice del módulo de juegos
#' @param w1 ponderador para grafos (defecto: 0.35)
#' @param w2 ponderador para Markov (defecto: 0.40)
#' @param w3 ponderador para juegos (defecto: 0.25)
#' @param R_factor factor de recursos (defecto: 0.65)
#' @return lista con I_MPCS, k, tipo de nudge
#' @export
calcular_indice <- function(I_grafo, I_markov, I_juegos, 
                            w1 = 0.35, w2 = 0.40, w3 = 0.25,
                            R_factor = 0.65) {
  
  # --- Verificar ponderadores ---
  # Asegurar que los ponderadores suman 1
  total <- w1 + w2 + w3
  if (abs(total - 1) > 0.01) {
    warning("Los ponderadores no suman 1. Se normalizarán.")
    w1 <- w1 / total
    w2 <- w2 / total
    w3 <- w3 / total
  }
  
  # --- Calcular I_MPCS ---
  I_MPCS <- w1 * I_grafo + w2 * I_markov + w3 * I_juegos
  
  # Asegurar que I_MPCS esté en [0, 1]
  I_MPCS <- max(0, min(1, I_MPCS))
  
  # --- Calcular k (intensidad del nudge) ---
  # k se calcula como una función del I_MPCS y el factor de recursos
  # Fórmula: k = min(1, I_MPCS * R_factor * 1.5)
  k <- min(1, I_MPCS * R_factor * 1.5)
  
  # --- Determinar tipo de nudge ---
  tipo <- dplyr::case_when(
    k < 0.25 ~ "Informativo",
    k < 0.50 ~ "Estructural",
    k < 0.75 ~ "Normativo",
    TRUE ~ "Sistémico multi-nudge"
  )
  
  # --- Retornar resultados ---
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
  
  # --- Validaciones ---
  if (is.null(P) || !is.matrix(P)) {
    warning("Se requiere una matriz de transición válida")
    return(NULL)
  }
  
  if (k < 0 || k > 1) {
    warning("k debe estar entre 0 y 1. Se usará k = 0.4")
    k <- 0.4
  }
  
  # --- Aplicar nudge ---
  P_n <- P
  m <- nrow(P)
  
  for (i in 1:(m-1)) {
    # Transferir probabilidad desde inercia al avance
    av <- P[i, i] * k
    P_n[i, i] <- P[i, i] - av
    P_n[i, i+1] <- P[i, i+1] + av
    P_n[i, ] <- P_n[i, ] / sum(P_n[i, ])
  }
  
  return(P_n)
}

# ============================================================================
# 6. generar_demo_data — Genera datos de demostración
# ============================================================================

#' Generar datos de demostración para probar la aplicación
#'
#' @param n número de observaciones
#' @param seed semilla para reproducibilidad
#' @return data.frame con datos de demostración
#' @export
generar_demo_data <- function(n = 1000, seed = 123) {
  
  # --- Fijar semilla ---
  set.seed(seed)
  
  # --- Generar datos ---
  data.frame(
    # Identificador
    ID = 1:n,
    
    # Variables demográficas
    Region = sample(c("Norte", "Sur", "Este", "Oeste", "Centro"), n, replace = TRUE),
    Edad = sample(25:80, n, replace = TRUE),
    Sexo = sample(c("M", "F"), n, replace = TRUE),
    Educacion = sample(0:4, n, replace = TRUE),
    
    # Variables de acceso y diagnóstico
    Acceso_Salud = rbinom(n, 1, 0.35),
    Tiene_Seguro = rbinom(n, 1, 0.60),
    Dx_HTA = rbinom(n, 1, 0.15),
    Dx_DM = rbinom(n, 1, 0.08),
    
    # Variables de adherencia y control
    Adh_Farma = runif(n, 0, 1),
    HTA_Medida = rbinom(n, 1, 0.20),
    Control_PA = rbinom(n, 1, 0.40),
    
    # Variables antropométricas
    IMC = rnorm(n, 26, 4),
    Obesidad_Abd = rbinom(n, 1, 0.30),
    
    # Variables de estilo de vida
    Fuma = rbinom(n, 1, 0.15),
    Alcohol = rbinom(n, 1, 0.20),
    Dieta_Sana = runif(n, 0, 1),
    Actividad_Fisica = sample(0:2, n, replace = TRUE),
    
    # Variables psicosociales
    Depresion_Bin = rbinom(n, 1, 0.10),
    Ansiedad_Bin = rbinom(n, 1, 0.15),
    Apoyo_Social = sample(1:5, n, replace = TRUE),
    
    # Variables de conocimiento
    Conocimiento_HTA = sample(0:10, n, replace = TRUE),
    Conocimiento_DM = sample(0:10, n, replace = TRUE),
    
    # Estado de Markov (conductual)
    Estado_Markov = sample(c("E1", "E2", "E3", "E4", "E5"), n,
                           replace = TRUE,
                           prob = c(0.45, 0.18, 0.15, 0.12, 0.10))
  )
}

# ============================================================================
# 7. format_report — Formatea resultados para reporte
# ============================================================================

#' Formatear resultados para reporte
#'
#' @param results data.frame con resultados del MPCS
#' @param plots lista de gráficos generados (opcional)
#' @return lista formateada para reporte
#' @export
format_report <- function(results, plots = NULL) {
  
  # --- Validaciones ---
  if (missing(results) || is.null(results)) {
    warning("Se requieren resultados")
    return(NULL)
  }
  
  # --- Formatear resumen ---
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
  
  # --- Añadir gráficos si están disponibles ---
  if (!is.null(plots)) {
    report$plots <- plots
  }
  
  # --- Añadir clasificación de grupos ---
  report$ranking <- results %>%
    dplyr::arrange(dplyr::desc(I_MPCS)) %>%
    dplyr::mutate(
      Prioridad = 1:n(),
      Nivel = dplyr::case_when(
        I_MPCS >= 0.75 ~ "Alta",
        I_MPCS >= 0.50 ~ "Media",
        TRUE ~ "Baja"
      )
    )
  
  return(report)
}

# ============================================================================
# 8. validar_datos — Valida la estructura de los datos de entrada
# ============================================================================

#' Validar la estructura de los datos de entrada
#'
#' @param datos data.frame a validar
#' @param min_filas número mínimo de filas (defecto: 30)
#' @param min_vars_num número mínimo de variables numéricas (defecto: 5)
#' @return lista con resultado de validación
#' @export
validar_datos <- function(datos, min_filas = 30, min_vars_num = 5) {
  
  # --- Validaciones básicas ---
  if (missing(datos) || is.null(datos)) {
    return(list(
      valido = FALSE,
      errores = c("No se proporcionaron datos")
    ))
  }
  
  errores <- c()
  
  # --- Verificar número de filas ---
  if (nrow(datos) < min_filas) {
    errores <- c(errores, paste("Se necesitan al menos", min_filas, "filas. Actual:", nrow(datos)))
  }
  
  # --- Verificar variables numéricas ---
  vars_num <- names(datos)[sapply(datos, is.numeric)]
  if (length(vars_num) < min_vars_num) {
    errores <- c(errores, paste("Se necesitan al menos", min_vars_num, 
                                "variables numéricas. Actual:", length(vars_num)))
  }
  
  # --- Verificar valores faltantes ---
  total_na <- sum(is.na(datos))
  if (total_na > 0) {
    pct_na <- round(total_na / (nrow(datos) * ncol(datos)) * 100, 1)
    if (pct_na > 20) {
      errores <- c(errores, paste("Alto porcentaje de valores faltantes:", pct_na, "%"))
    }
  }
  
  # --- Retornar resultado ---
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
# FIN DEL ARCHIVO
# ============================================================================
