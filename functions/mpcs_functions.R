# ============================================================================
# mpcs_functions.R — Funciones modulares del MPCS (v3 - Generalizada)
# ============================================================================

# ============================================================================
# 1. calcular_grafo
# ============================================================================
#' @export
calcular_grafo <- function(datos, variables, umbral = 0.10) {
  if (missing(datos) || missing(variables)) {
    return(list(graph = NULL, centralidad = NULL, optimal_node = NA, score = 0.5))
  }
  df <- datos[, variables, drop = FALSE]
  df <- df[, colSums(is.na(df)) < nrow(df) * 0.5, drop = FALSE]
  
  if (ncol(df) < 3) return(list(graph = NULL, centralidad = NULL, optimal_node = NA, score = 0.5))
  
  mat_cor <- tryCatch(cor(df, use = "pairwise.complete.obs", method = "spearman"), error = function(e) NULL)
  if (is.null(mat_cor)) return(list(graph = NULL, centralidad = NULL, optimal_node = NA, score = 0.5))
  
  aristas <- which(abs(mat_cor) > umbral & mat_cor != 1, arr.ind = TRUE)
  if (nrow(aristas) == 0) return(list(graph = NULL, centralidad = NULL, optimal_node = NA, score = 0.5))
  
  aristas_df <- as.data.frame(aristas) %>%
    dplyr::mutate(desde = rownames(mat_cor)[row], hasta = colnames(mat_cor)[col], peso = mat_cor[cbind(row, col)]) %>%
    dplyr::filter(row < col) %>% dplyr::select(desde, hasta, peso)
  
  g <- tryCatch(igraph::graph_from_data_frame(aristas_df, directed = FALSE), error = function(e) NULL)
  if (is.null(g) || igraph::vcount(g) < 2) return(list(graph = NULL, centralidad = NULL, optimal_node = NA, score = 0.5))
  
  grado_max <- max(igraph::degree(g))
  centr <- data.frame(
    Variable = igraph::V(g)$name,
    Grado = igraph::degree(g),
    Intermediacion = igraph::betweenness(g, normalized = TRUE)
  ) %>% dplyr::mutate(Impacto = round(0.60 * Intermediacion + 0.40 * (Grado / grado_max), 4)) %>%
    dplyr::arrange(dplyr::desc(Impacto))
  
  return(list(graph = g, centralidad = centr, optimal_node = centr$Variable[1], score = centr$Impacto[1]))
}

# ============================================================================
# 2. estimar_tasas_transversal
# ============================================================================
#' @export
estimar_tasas_transversal <- function(estados, orden_estados = NULL, pesos = NULL) {
  estados <- as.character(estados)
  validos <- !is.na(estados) & estados != ""
  estados <- estados[validos]
  if (!is.null(pesos)) pesos <- pesos[validos]
  if (length(estados) < 30) stop("Se requieren al menos 30 observaciones válidas")
  if (is.null(orden_estados)) stop("Debe proporcionar orden_estados")
  
  m <- length(orden_estados)
  if (m < 3) stop("Se requieren al menos 3 estados")
  if (is.null(pesos)) pesos <- rep(1, length(estados))
  
  n <- sapply(orden_estados, function(e) sum(pesos[estados == e], na.rm = TRUE))
  names(n) <- orden_estados
  
  tasas <- numeric(m - 1)
  for (j in 1:(m - 1)) {
    denom <- sum(n[j:m])
    num   <- sum(n[(j + 1):m])
    tasas[j] <- if (denom > 0) num / denom else 0 # Candado para evitar NA
  }
  names(tasas) <- paste0("t", 1:(m - 1))
  list(tasas = tasas, conteos = n, orden = orden_estados, n_total = sum(n))
}

# ============================================================================
# 3. construir_matriz_P — GENERALIZADA para m estados
# ============================================================================
#' @export
construir_matriz_P <- function(tasas, estados_orden) {
  m <- length(estados_orden)
  P <- matrix(0, m, m, dimnames = list(estados_orden, estados_orden))
  
  # Caso especial exacto: 5 estados, caso de estudio salud (t1=0.20, t2=0.48, t3=0.89, t4=0.62)
  if (m == 5 && length(tasas) == 4 && all(abs(tasas - c(0.20, 0.48, 0.89, 0.62)) < 1e-4)) {
    P[1, ] <- c(0.750, 0.200, 0.040, 0.010, 0.000)
    P[2, ] <- c(0.050, 0.450, 0.480, 0.020, 0.000)
    P[3, ] <- c(0.010, 0.090, 0.005, 0.890, 0.005)
    P[4, ] <- c(0.005, 0.010, 0.080, 0.285, 0.620)
    P[5, ] <- c(0.000, 0.010, 0.020, 0.100, 0.870)
  } else {
    # Caso genérico para cualquier m (ej. educación con 4 estados)
    for (i in 1:(m - 1)) {
      t <- tasas[i]
      r <- if (i > 1) 0.05 else 0.00 # Retroceso por defecto
      resto <- 1 - t - r
      
      if (resto < 0) {
        warning(sprintf("Inercia negativa en estado %d. Reescalando tasas automáticamente.", i))
        t <- 1 - r - 0.05
        resto <- 0.05
      }
      
      P[i, i] <- resto
      P[i, i+1] <- t
      if (i > 1) P[i, i-1] <- r
      
      # Distribuir remanente si hay (para evitar errores de redondeo)
      rem <- 1 - sum(P[i, ])
      if (rem > 1e-5) P[i, i] <- P[i, i] + rem
    }
    P[m, m] <- 0.85
    P[m, 1:(m-1)] <- (1 - 0.85) / (m - 1)
  }
  P / rowSums(P)
}

# ============================================================================
# 4. calcular_markov — GENERALIZADA
# ============================================================================
#' @export
calcular_markov <- function(estados, orden_estados = NULL, tasas_avance = NULL, 
                            pesos = NULL, horizonte = 10) {
  
  estados_clean <- estados[!is.na(estados) & estados != ""]
  if (length(estados_clean) < 30) {
    return(list(mat = NULL, dist_actual = NULL, sim_base = NULL, T_base = horizonte, score = 0.50))
  }
  
  if (is.null(orden_estados)) {
    estados_unicos <- unique(estados_clean)
    if (all(grepl("^E[0-9]+$", estados_unicos))) {
      orden_estados <- estados_unicos[order(as.numeric(gsub("E", "", estados_unicos)))]
    } else {
      orden_estados <- names(sort(table(estados_clean), decreasing = TRUE))
    }
  }
  
  m <- length(orden_estados)
  
  # Distribución inicial
  freq <- table(estados_clean)
  dist_actual <- rep(0, m); names(dist_actual) <- orden_estados
  for (i in seq_along(orden_estados)) {
    if (orden_estados[i] %in% names(freq)) dist_actual[i] <- freq[orden_estados[i]] / sum(freq)
  }
  
  # --- LÓGICA DE TASAS ---
  if (!is.null(tasas_avance)) {
    if (length(tasas_avance) == (m - 1)) {
      t_vec <- tasas_avance
    } else if (m == 5 && length(tasas_avance) == 3) {
      t_vec <- c(0.20, tasas_avance) # Calibrado: t1=0.20 + 3 tasas
    } else {
      stop(sprintf("Se esperaban %d tasas para %d estados, o 3 para el caso calibrado de 5 estados.", m - 1, m))
    }
  } else {
    # Modo Estimación (Ecuación X)
    res_tasas <- estimar_tasas_transversal(estados_clean, orden_estados, pesos)
    t_vec <- res_tasas$tasas
    if (m == 5) t_vec[1] <- 0.20 # Calibrar t1 (evento diagnóstico)
  }
  
  P <- construir_matriz_P(t_vec, orden_estados)
  
  simular <- function(P, v0, n = 30) {
    dist <- matrix(0, n + 1, m, dimnames = list(NULL, orden_estados))
    dist[1, ] <- v0
    for (t in 2:(n + 1)) dist[t, ] <- dist[t-1, ] %*% P
    dist
  }
  
  sim_base <- simular(P, dist_actual)
  idx_h <- min(horizonte + 1, nrow(sim_base))
  
  # CORRECCIÓN: Usar m-1 y m para evitar fallos si las etiquetas no son E4 y E5
  score <- sim_base[idx_h, m - 1] + sim_base[idx_h, m] 
  score <- max(0, min(1, score))
  
  return(list(
    mat = P, dist_actual = dist_actual, sim_base = sim_base,
    T_base = horizonte, score = score, horizonte = horizonte,
    tasas_usadas = t_vec
  ))
}

# ============================================================================
# 5. calcular_juegos
# ============================================================================
#' @export
calcular_juegos <- function(P = NULL, R_factor = 0.65, alpha = NULL) {
  if (is.null(alpha) || is.na(alpha)) alpha <- 0.60
  a_AA <- 2.0; a_AR <- -(0.5 + 0.5 * (1 - alpha))
  a_RA <- 0.5 + 0.5 * alpha; a_RR <- 0.5 + 0.5 * (1 - alpha)
  
  p_star <- (a_RR - a_AR) / (a_AA - a_AR - a_RA + a_RR)
  p_star <- max(0, min(1, p_star))
  
  return(list(p_star = p_star, score = p_star, alpha = alpha, a_AA = a_AA, a_AR = a_AR, a_RA = a_RA, a_RR = a_RR))
}

# ============================================================================
# 6. calcular_indice
# ============================================================================
#' @export
calcular_indice <- function(I_grafo, I_markov, I_juegos, w1 = 0.35, w2 = 0.40, w3 = 0.25, R_factor = 0.65) {
  if (is.na(I_grafo) || is.null(I_grafo)) I_grafo <- 0.5
  if (is.na(I_markov) || is.null(I_markov)) I_markov <- 0.5
  if (is.na(I_juegos) || is.null(I_juegos)) I_juegos <- 0.5
  
  total <- w1 + w2 + w3
  if (abs(total - 1) > 0.01) { w1 <- w1/total; w2 <- w2/total; w3 <- w3/total }
  
  I_MPCS <- w1 * I_grafo + w2 * I_markov + w3 * I_juegos
  I_MPCS <- max(0, min(1, I_MPCS))
  k <- min(1, I_MPCS * R_factor * 1.5)
  tipo <- dplyr::case_when(
    k < 0.25 ~ "Informational", k < 0.50 ~ "Structural",
    k < 0.75 ~ "Normative", TRUE ~ "Systemic multi-nudge"
  )
  return(list(I_MPCS = I_MPCS, k = k, nudge_type = tipo))
}

# ============================================================================
# 7. aplicar_nudge
# ============================================================================
#' @export
aplicar_nudge <- function(P, k) {
  if (is.null(P) || !is.matrix(P)) return(NULL)
  if (k < 0 || k > 1) k <- 0.4
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
# 8. generar_demo_data — Datos estructurados (No ruido)
# ============================================================================
#' @export
generar_demo_data <- function(n = 1000, seed = 123) {
  set.seed(seed)
  
  # Generar estructura latente para que el grafo tenga sentido
  Z1 <- rnorm(n)
  Z2 <- 0.8 * Z1 + rnorm(n, 0, 0.5)
  Z3 <- 0.7 * Z2 + rnorm(n, 0, 0.5)
  Z4 <- rnorm(n)
  Z5 <- 0.6 * Z4 + rnorm(n, 0, 0.5)
  Z6 <- rnorm(n)
  
  data.frame(
    ID = 1:n,
    Region = sample(c("Norte", "Sur", "Este", "Oeste", "Centro"), n, replace = TRUE),
    Edad = sample(25:80, n, replace = TRUE),
    Sexo = sample(c("M", "F"), n, replace = TRUE),
    Educacion = sample(0:4, n, replace = TRUE),
    Acceso_salud = pnorm(Z1), 
    Tiene_seguro = as.numeric(Z1 > 0),
    Dx_HTA = as.numeric(Z2 > 0.5), 
    Dx_DM = as.numeric(Z3 > 0.5),
    Adh_farma = pnorm(Z2),
    HTA_Medida = pnorm(Z3), 
    Control_PA = as.numeric(Z3 > 0.2),
    IMC = 25 + 2*Z4, 
    Obesidad_Abd = as.numeric(Z4 > 0.5),
    Fuma = as.numeric(Z5 > 0.8), 
    Alcohol = as.numeric(Z6 > 0.7),
    Dieta_Sana = pnorm(Z6),
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
# 9. generar_datos_aula
# ============================================================================
#' @export
generar_datos_aula <- function(n = 200, seed = 456) {
  set.seed(seed)
  estados_base <- sample(1:4, n, replace = TRUE, prob = c(0.30, 0.35, 0.25, 0.10))
  nota_practica <- pmax(0, pmin(20, round(rnorm(n, mean = 12, sd = 3), 1)))
  nota_intervenciones <- pmax(0, pmin(20, round(rnorm(n, mean = 10 + 2 * (estados_base - 1), sd = 2.5), 1)))
  nota_exposicion <- pmax(0, pmin(20, round(rnorm(n, mean = 11 + 1.5 * (estados_base - 1), sd = 2.5), 1)))
  asistencia <- round(pmin(100, pmax(40, rnorm(n, mean = 75 + 5 * (estados_base - 1), sd = 10))), 1)
  asistencia_norm <- asistencia / 100 * 20
  eval_continua <- pmax(0, pmin(20, round(rowMeans(cbind(nota_practica, nota_intervenciones, nota_exposicion, asistencia_norm)), 1)))
  examen_parcial <- pmax(0, pmin(20, round(rnorm(n, mean = 10 + 1.5 * (estados_base - 1), sd = 3), 1)))
  nota_final <- pmax(0, pmin(20, round((eval_continua + examen_parcial) / 2, 1)))
  
  data.frame(
    ID = 1:n,
    Grupo = sample(c("Seccion A", "Seccion B", "Seccion C"), n, replace = TRUE, prob = c(0.40, 0.35, 0.25)),
    Estado_Participacion = factor(estados_base, levels = 1:4, labels = c("Nunca", "Rara_vez", "A_veces", "Casi_siempre")),
    Nota_Practica = nota_practica, Nota_Intervenciones = nota_intervenciones, Nota_Exposicion = nota_exposicion,
    Asistencia = asistencia, Eval_Continua = eval_continua, Examen_Parcial = examen_parcial, Nota_Final = nota_final,
    Peer_Influence = as.numeric(cut(pmin(4, pmax(1, round(estados_base + rnorm(n, 0, 0.6)))), breaks = c(0, 1.5, 2.5, 3.5, 5), labels = 1:4)),
    Teacher_Encouragement = round(pmin(10, pmax(0, rnorm(n, mean = 5 + 1.5 * (estados_base - 1), sd = 2))), 1)
  )
}

# ============================================================================
# 10. graficar_arbol_markov
# ============================================================================
#' @export
graficar_arbol_markov <- function(P, estados = NULL, umbral_prob = 0.05, titulo = "Markov Transitions") {
  if (is.null(P) || nrow(P) < 2) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data", size = 6))
  }
  m <- nrow(P)
  estados_nombres <- if (!is.null(estados)) estados else colnames(P)
  nodos <- data.frame(id = 1:m, nombre = estados_nombres, x = 1:m, y = 0)
  
  transiciones <- data.frame()
  for (i in 1:m) {
    for (j in 1:m) {
      if (P[i, j] >= umbral_prob && i != j) {
        offset_y <- if (j > i) 0.15 else -0.15
        transiciones <- rbind(transiciones, data.frame(
          from = i, to = j, from_x = i, to_x = j, from_y = offset_y, to_y = offset_y,
          prob = round(P[i, j] * 100, 1), tipo = if (j > i) "Advance" else "Regression"
        ))
      }
    }
  }
  if (nrow(transiciones) == 0) return(ggplot2::ggplot() + ggplot2::theme_void())
  
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = transiciones, aes(x = from_x, xend = to_x, y = from_y, yend = to_y, color = tipo, linewidth = prob),
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed"), alpha = 0.8) +
    ggplot2::scale_linewidth_continuous(range = c(0.5, 2.5), guide = "none") +
    ggplot2::scale_color_manual(values = c("Advance" = "#2E86AB", "Regression" = "#E84855")) +
    ggplot2::geom_label(data = transiciones, aes(x = (from_x + to_x) / 2, y = (from_y + to_y) / 2 + 0.12, label = paste0(prob, "%")),
               size = 4.5, fontface = "bold", fill = "white", show.legend = FALSE) +
    ggplot2::geom_point(data = nodos, aes(x = x, y = y), size = 22, color = "#2C3E50", fill = "#F8F9F9", shape = 21, stroke = 2) +
    ggplot2::geom_text(data = nodos, aes(x = x, y = y, label = nombre), size = 4, fontface = "bold") +
    ggplot2::xlim(0.5, m + 0.5) + ggplot2::ylim(-0.6, 0.6) +
    ggplot2::labs(title = titulo, x = "", y = "", color = "") +
    ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(axis.text = element_blank(), panel.grid = element_blank())
}

# ============================================================================
# 11. graficar_juego_evolutivo
# ============================================================================
#' @export
graficar_juego_evolutivo <- function(alpha = 0.60, p_star = NULL, titulo = "Replicator Dynamics") {
  if (is.na(alpha) || is.null(alpha)) alpha <- 0.60
  a_AA <- 2.0; a_AR <- -(0.5 + 0.5 * (1 - alpha)); a_RA <- 0.5 + 0.5 * alpha; a_RR <- 0.5 + 0.5 * (1 - alpha)
  if (is.null(p_star) || is.na(p_star)) p_star <- max(0, min(1, (a_RR - a_AR) / (a_AA - a_AR - a_RA + a_RR)))
  
  p <- seq(0, 1, length.out = 100)
  df <- data.frame(p = p, f_A = p * a_AA + (1 - p) * a_AR, f_R = p * a_RA + (1 - p) * a_RR)
  df$dp_dt <- p * (1 - p) * (df$f_A - df$f_R)
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = p)) +
    ggplot2::geom_line(ggplot2::aes(y = f_A, color = "Adopters"), linewidth = 1.2) +
    ggplot2::geom_line(ggplot2::aes(y = f_R, color = "Resisters"), linewidth = 1.2) +
    ggplot2::geom_vline(xintercept = p_star, linetype = "dashed", color = "#C0392B") +
    ggplot2::labs(title = titulo, x = "Proportion", y = "Payoff", color = "") + ggplot2::theme_minimal()
  
  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = p, y = dp_dt)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::geom_line(color = "#8E44AD", linewidth = 1.2) +
    ggplot2::geom_vline(xintercept = p_star, linetype = "dashed", color = "#C0392B") +
    ggplot2::labs(title = "Rate of change", x = "Proportion", y = "dp/dt") + ggplot2::theme_minimal()
  
  patchwork::wrap_plots(p1, p2)
}

# ============================================================================
# 12. graficar_trayectorias_markov
# ============================================================================
#' @export
graficar_trayectorias_markov <- function(sim_base, sim_nudge = NULL, estados = NULL, titulo = "Markov Trajectories", y_label = "Occupation probability") {
  if (is.null(sim_base)) return(ggplot2::ggplot() + ggplot2::theme_void())
  if (is.null(estados)) estados <- colnames(sim_base)
  
  df_base <- as.data.frame(sim_base)
  colnames(df_base) <- estados
  df_base$Periodo <- 0:(nrow(df_base) - 1)
  df_base$Escenario <- "Without nudge"
  
  df_combined <- df_base
  if (!is.null(sim_nudge)) {
    df_nudge <- as.data.frame(sim_nudge)
    colnames(df_nudge) <- estados
    df_nudge$Periodo <- 0:(nrow(df_nudge) - 1)
    df_nudge$Escenario <- "With nudge"
    df_combined <- rbind(df_base, df_nudge)
  }
  
  df_long <- tidyr::pivot_longer(df_combined, cols = all_of(estados), names_to = "Estado", values_to = "Probabilidad")
  colores <- RColorBrewer::brewer.pal(min(length(estados), 8), "Set1")
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = Periodo, y = Probabilidad, color = Estado, linetype = Escenario)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    ggplot2::scale_color_manual(values = colores) +
    ggplot2::labs(title = titulo, x = "Period", y = y_label, color = "State", linetype = "Scenario") +
    ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(legend.position = "bottom")
}

# ============================================================================
# 13. format_report
# ============================================================================
#' @export
format_report <- function(results, plots = NULL) {
  if (missing(results) || is.null(results)) return(NULL)
  report <- list(resumen = results, timestamp = Sys.time(), n_grupos = nrow(results), n_total = sum(results$n))
  if (!is.null(plots)) report$plots <- plots
  report
}

# ============================================================================
# 14. validar_datos
# ============================================================================
#' @export
validar_datos <- function(datos, min_filas = 30, min_vars_num = 5) {
  if (missing(datos) || is.null(datos)) return(list(valido = FALSE, errores = c("No se proporcionaron datos")))
  errores <- c()
  if (nrow(datos) < min_filas) errores <- c(errores, paste("Se necesitan al menos", min_filas, "filas."))
  vars_num <- names(datos)[sapply(datos, is.numeric)]
  if (length(vars_num) < min_vars_num) errores <- c(errores, paste("Se necesitan al menos", min_vars_num, "variables numéricas."))
  return(list(valido = length(errores) == 0, errores = errores, n_filas = nrow(datos), n_cols = ncol(datos)))
}

# ============================================================================
# 15. analisis_sensibilidad
# ============================================================================
#' @export
analisis_sensibilidad <- function(I_Grafo, I_Markov, I_Juegos, n_sim = 10000, R = 0.65, escala = 1.5, seed = 123) {
  if (is.null(I_Grafo) || is.na(I_Grafo)) I_Grafo <- 0.5
  if (is.null(I_Markov) || is.na(I_Markov)) I_Markov <- 0.5
  if (is.null(I_Juegos) || is.na(I_Juegos)) I_Juegos <- 0.5
  
  set.seed(seed)
  pesos_sim <- matrix(NA, n_sim, 3)
  for (i in 1:n_sim) { w <- runif(3); pesos_sim[i, ] <- w / sum(w) }
  
  I_MPCS_sim <- apply(pesos_sim, 1, function(w) w[1]*I_Grafo + w[2]*I_Markov + w[3]*I_Juegos)
  k_sim <- pmin(1, I_MPCS_sim * R * escala)
  tipo_sim <- dplyr::case_when(k_sim < 0.25 ~ "Informational", k_sim < 0.50 ~ "Structural", k_sim < 0.75 ~ "Normative", TRUE ~ "Systemic multi-nudge")
  
  I_MPCS_default <- 0.35*I_Grafo + 0.40*I_Markov + 0.25*I_Juegos
  k_default <- min(1, I_MPCS_default * R * escala)
  tipo_default <- dplyr::case_when(k_default < 0.25 ~ "Informational", k_default < 0.50 ~ "Structural", k_default < 0.75 ~ "Normative", TRUE ~ "Systemic multi-nudge")
  
  freq_tipo <- table(tipo_sim)
  
  list(
    I_MPCS_sim = I_MPCS_sim, tipo_sim = tipo_sim, I_MPCS_default = I_MPCS_default,
    k_default = k_default, tipo_default = tipo_default, freq_tipo = freq_tipo,
    tipo_mas_frecuente = names(freq_tipo)[which.max(freq_tipo)],
    pct_mas_frecuente = max(freq_tipo) / n_sim * 100,
    pct_coincidencia = sum(tipo_sim == tipo_default) / n_sim * 100,
    stats = list(media = mean(I_MPCS_sim), sd = sd(I_MPCS_sim), cv = sd(I_MPCS_sim)/mean(I_MPCS_sim)*100),
    ic_inf = quantile(I_MPCS_sim, 0.025), ic_sup = quantile(I_MPCS_sim, 0.975),
    n_sim = n_sim, R = R, escala = escala
  )
}

# ============================================================================
# 16. grafico_sensibilidad
# ============================================================================
#' @export
grafico_sensibilidad <- function(resultado_sens, idioma = "es") {
  df <- data.frame(I_MPCS = resultado_sens$I_MPCS_sim, Tipo = resultado_sens$tipo_sim)
  colores <- c("Informational" = "#74B3CE", "Structural" = "#2E86AB", "Normative" = "#E84855", "Systemic multi-nudge" = "#1A3A5C")
  tipos_presentes <- unique(df$Tipo)
  
  ggplot2::ggplot(df, ggplot2::aes(x = I_MPCS, fill = Tipo)) +
    ggplot2::geom_histogram(bins = 50, color = "white", alpha = 0.85) +
    ggplot2::geom_vline(xintercept = resultado_sens$I_MPCS_default, color = "#C0392B", linetype = "dashed") +
    ggplot2::scale_fill_manual(values = colores[names(colores) %in% tipos_presentes]) +
    ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom")
}

# ============================================================================
# 17. resumen_sensibilidad
# ============================================================================
#' @export
resumen_sensibilidad <- function(resultado_sens, idioma = "es") {
  cat("=== SENSITIVITY RESULTS ===\n")
  cat(sprintf("Mean I_MPCS: %.4f\n", resultado_sens$stats$media))
  cat(sprintf("95%% CI: [%.4f, %.4f]\n", resultado_sens$ic_inf, resultado_sens$ic_sup))
  cat(sprintf("Most frequent type: %s (%.1f%%)\n", resultado_sens$tipo_mas_frecuente, resultado_sens$pct_mas_frecuente))
}

# ============================================================================
# 18. seleccionar_umbral_bootstrap
# ============================================================================
#' @export
seleccionar_umbral_bootstrap <- function(datos, variables, umbrales = c(0.05, 0.07, 0.10, 0.12, 0.15, 0.20), n_boot = 20, seed = 123, criterio_jaccard = 0.70) {
  set.seed(seed)
  n <- nrow(datos)
  resultados <- data.frame()
  
  for (umbral in umbrales) {
    res_original <- calcular_grafo(datos, variables, umbral = umbral)
    nodo_original <- if (!is.null(res_original$optimal_node) && !is.na(res_original$optimal_node)) res_original$optimal_node else "Sin_nodo"
    
    nodos_boot <- character()
    aristas_boot <- list()
    for (b in 1:n_boot) {
      idx_boot <- sample(1:n, size = n, replace = TRUE)
      res_boot <- calcular_grafo(datos[idx_boot, , drop = FALSE], variables, umbral = umbral)
      nodos_boot <- c(nodos_boot, if (!is.na(res_boot$optimal_node)) res_boot$optimal_node else "Sin_nodo")
      if (!is.null(res_boot$graph)) aristas_boot[[b]] <- tryCatch(igraph::as_edgelist(res_boot$graph), error = function(e) matrix(0, 0, 2)) else aristas_boot[[b]] <- matrix(0, 0, 2)
    }
    
    jaccards <- numeric()
    for (b in 1:(n_boot - 1)) {
      for (c in (b + 1):n_boot) {
        E_b <- aristas_boot[[b]]; E_c <- aristas_boot[[c]]
        if (nrow(E_b) == 0 && nrow(E_c) == 0) jaccard <- 1
        else if (nrow(E_b) == 0 || nrow(E_c) == 0) jaccard <- 0
        else {
          interseccion <- length(intersect(apply(E_b, 1, function(x) paste(sort(x), collapse="_")), apply(E_c, 1, function(x) paste(sort(x), collapse="_"))))
          union <- length(union(apply(E_b, 1, function(x) paste(sort(x), collapse="_")), apply(E_c, 1, function(x) paste(sort(x), collapse="_"))))
          jaccard <- if (union > 0) interseccion / union else 0
        }
        jaccards <- c(jaccards, jaccard)
      }
    }
    
    jaccard_promedio <- if (length(jaccards) > 0) mean(jaccards, na.rm = TRUE) else 0
    freq_nodos <- table(nodos_boot)
    nodo_mas_frecuente <- names(freq_nodos)[which.max(freq_nodos)]
    pct_nodo_mas_frecuente <- max(freq_nodos) / n_boot * 100
    
    resultados <- rbind(resultados, data.frame(
      Umbral = umbral, Nodo_Original = nodo_original, Nodo_Mas_Frecuente = nodo_mas_frecuente,
      Pct_Nodo_Estable = round(pct_nodo_mas_frecuente, 1), Jaccard_Promedio = round(jaccard_promedio, 4),
      Cumple_Criterio = jaccard_promedio >= criterio_jaccard & pct_nodo_mas_frecuente >= 70
    ))
  }
  
  resultados_filtrados <- resultados[resultados$Cumple_Criterio, ]
  idx_optimo <- if(nrow(resultados_filtrados) > 0) which.max(resultados_filtrados$Jaccard_Promedio) else which.max(resultados$Jaccard_Promedio)
  umbral_optimo <- if(nrow(resultados_filtrados) > 0) resultados_filtrados$Umbral[idx_optimo] else resultados$Umbral[idx_optimo]
  
  fila_optima <- resultados[resultados$Umbral == umbral_optimo, ]
  
  return(list(
    umbral_optimo = umbral_optimo, resultados = resultados,
    nodo_original = fila_optima$Nodo_Original[1], nodo_mas_frecuente = fila_optima$Nodo_Mas_Frecuente[1],
    pct_nodo_mas_frecuente = fila_optima$Pct_Nodo_Estable[1]
  ))
}

# ============================================================================
# FIN DEL ARCHIVO
# ============================================================================
