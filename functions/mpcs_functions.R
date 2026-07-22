# =======================================================
# Funciones Modulares del Modelo MPCS
# =======================================================

#' Calcular Grafo (Teoría de Grafos)
#' @param data_sub Dataframe con variables numéricas/categóricas codificadas
#' @param threshold Umbral de correlación
#' @return Lista con grafo, centralidades y nodo óptimo
calcular_grafo <- function(data_sub, threshold) {
  corr_matrix <- cor(data_sub, use = "pairwise.complete.obs")
  adj_matrix <- ifelse(abs(corr_matrix) > threshold & corr_matrix != 1, 1, 0)
  g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = NULL, diag = FALSE)
  
  if(igraph::vcount(g) == 0) return(list(score = 0, optimal_node = NA, graph = NULL))
  
  bet <- igraph::betweenness(g, normalized = TRUE)
  deg <- igraph::degree(g, normalized = TRUE)
  centralidad <- (bet + deg) / 2
  optimal_node <- names(which.max(centralidad))
  
  return(list(
    graph = g,
    centralities = centralidad,
    optimal_node = optimal_node,
    score = max(centralidad, na.rm = TRUE)
  ))
}

#' Simular Markov (Cadenas de Markov)
#' @param states Vector de estados
#' @return Lista con matriz de transición, convergencia y score
simular_markov <- function(states) {
  states <- factor(states)
  if(nlevels(states) < 3) return(list(score = 0, mat = NULL))
  
  trans_mat <- table(states[-length(states)], states[-1])
  trans_mat <- prop.table(trans_mat, 1)
  trans_mat[is.na(trans_mat)] <- 0
  
  # Tiempo de convergencia (aproximación simple)
  p <- trans_mat
  conv_time <- 0
  for(t in 1:100) {
    p <- p %*% trans_mat
    conv_time <- conv_time + 1
    if(max(abs(p - trans_mat)) < 0.001) break
  }
  
  # Score: mayor diversidad de transiciones (entropía)
  ent <- sum(trans_mat * log(trans_mat + 1e-10), na.rm = TRUE)
  score <- abs(ent) / log(nrow(trans_mat))
  
  return(list(mat = trans_mat, conv_time = conv_time, score = score))
}

#' Calcular Juegos (Teoría de Juegos Evolutiva)
#' @param markov_mat Matriz de transición de Markov
#' @param R Factor de recursos
#' @return Lista con masa crítica y score
calcular_juegos <- function(markov_mat, R) {
  if(is.null(markov_mat)) return(list(score = 0, critical_mass = NA))
  
  # Ecuación del replicador simplificada
  # Asumimos estrategia A (cambio) vs B (no cambio)
  payoffs <- matrix(c(1, 0.5, 0.5, 0), nrow=2) # Matriz de pagos genérica
  x <- markov_mat[1, ncol(markov_mat)] # Probabilidad de adoptar comportamiento
  dx_dt <- x * (1 - x) * ((payoffs[1,1] - payoffs[2,1]) * x + (payoffs[2,2] - payoffs[1,2]) * (1-x))
  
  critical_mass <- 1 - R # Masa crítica inversamente proporcional a recursos
  score <- abs(dx_dt)
  
  return(list(critical_mass = critical_mass, score = score))
}

#' Calcular Índice MPCS
#' @param graph_score Score del grafo
#' @param markov_score Score de Markov
#' @param game_score Score de Juegos
#' @param w1 Ponderador Grafo
#' @param w2 Ponderador Markov
#' @param w3 Ponderador Juegos
#' @return I_MPCS y tipo de nudge
calcular_indice <- function(graph_score, markov_score, game_score, w1, w2, w3) {
  I_MPCS <- (w1 * graph_score) + (w2 * markov_score) + (w3 * game_score)
  I_MPCS <- round(I_MPCS * 100, 2)
  
  nudge_type <- ifelse(I_MPCS > 75, "Arquitectura de Elección Fuerte (Structural)",
                       ifelse(I_MPCS > 50, "Nudge Social y de Default (Social)",
                              "Nudge Cognitivo Suave (Cognitive)"))
  
  return(list(I_MPCS = I_MPCS, nudge_type = nudge_type))
}
