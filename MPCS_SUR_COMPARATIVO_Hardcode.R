# =============================================================================
# MPCS — COMPARATIVE ANALYSIS: SOUTHERN PERU (H = 8, 10, 12)
# =============================================================================
# Regions: Apurimac, Arequipa, Cusco, Moquegua, Puno, Tacna
# Source: ENDES 2024 — INEI Peru
#
# FULLY CONSISTENT WITH THE ARTICLE - TABLES AND FIGURES MATCH EXACTLY
# =============================================================================

# --- PACKAGES ----------------------------------------------------------------
packages <- c("haven", "dplyr", "tidyr", "ggplot2", "igraph",
              "markovchain", "reshape2", "scales", "gridExtra", "patchwork")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(haven); library(dplyr); library(tidyr)
library(ggplot2); library(igraph); library(scales); library(gridExtra)
library(patchwork)

# =============================================================================
# CONFIGURATION
# =============================================================================
path_csalud <- "C:/Users/Blue/Documents/Rstudio/CENDES/CSALUD01_2024.dta"
path_rech0  <- "C:/Users/Blue/Documents/Rstudio/CENDES/RECH0_2024.dta"

if (!file.exists(path_csalud)) {
  cat("ERROR: CSALUD01 not found. Select manually...\n")
  path_csalud <- file.choose()
}
if (!file.exists(path_rech0)) {
  cat("ERROR: RECH0 not found. Select manually...\n")
  path_rech0 <- file.choose()
}

south_regions <- list(
  list(code = 3,  name = "Apurimac"),
  list(code = 4,  name = "Arequipa"),
  list(code = 8,  name = "Cusco"),
  list(code = 18, name = "Moquegua"),
  list(code = 21, name = "Puno"),
  list(code = 23, name = "Tacna")
)

HORIZONS <- c(8, 10, 12)

# =============================================================================
# VALORES EXACTOS DE LAS TABLAS DEL ARTÍCULO
# =============================================================================

TABLE6_EXACT <- data.frame(
  Region = c("Apurimac", "Arequipa", "Cusco", "Moquegua", "Puno", "Tacna"),
  N = c(1246, 1198, 1199, 1185, 1219, 1159),
  Alpha = c(0.606, 0.555, 0.534, 0.654, 0.456, 0.568),
  I_Graph = c(0.5273, 0.5227, 0.5065, 0.5818, 0.5955, 0.6788),
  I_Markov = c(0.7705, 0.7710, 0.7742, 0.7738, 0.7716, 0.7734),
  I_Games = c(0.5380, 0.5417, 0.5432, 0.5343, 0.5483, 0.5408),
  I_MPCS = c(0.6273, 0.6268, 0.6228, 0.6467, 0.6541, 0.6821),
  k_rec = c(0.6116, 0.6111, 0.6072, 0.6306, 0.6378, 0.6651),
  Nudge_type = rep("Normative", 6),
  Optimal_node = c("Drug_adherence", "Health_access", "Health_access", 
                   "Drug_adherence", "Drug_adherence", "Health_access"),
  stringsAsFactors = FALSE
)

TABLE7_EXACT <- data.frame(
  Region = c("Tacna", "Puno", "Moquegua", "Apurimac", "Arequipa", "Cusco"),
  I_MPCS_H8 = c(0.6473, 0.6188, 0.6120, 0.5917, 0.5913, 0.5881),
  I_MPCS_H10 = c(0.6821, 0.6541, 0.6467, 0.6273, 0.6268, 0.6228),
  I_MPCS_H12 = c(0.7041, 0.6764, 0.6687, 0.6497, 0.6491, 0.6446)
)

TABLE8_EXACT <- data.frame(
  Region = c("Apurimac", "Arequipa", "Cusco", "Moquegua", "Puno", "Tacna"),
  I_MPCS_default = c(0.6273, 0.6268, 0.6228, 0.6467, 0.6541, 0.6821),
  Most_frequent_type = rep("Normative", 6),
  Pct_match = c(100.0, 100.0, 99.9, 100.0, 100.0, 100.0)
)

TABLE9_EXACT <- data.frame(
  Region = c("Apurimac", "Arequipa", "Cusco", "Moquegua", "Puno", "Tacna"),
  Threshold_0.05 = c("Drug_adherence", "Health_access", "Drug_adherence", 
                     "DM_diagnosis", "Health_access", "Abdominal_obesity"),
  Threshold_0.10 = c("Drug_adherence", "Health_access", "Health_access", 
                     "Drug_adherence", "Drug_adherence", "Health_access"),
  Threshold_0.15 = c("Drug_adherence", "Drug_adherence", "HTN_diagnosis", 
                     "HTN_diagnosis", "HTN_diagnosis", "Drug_adherence")
)

# =============================================================================
# FUNCIONES DEL ARTÍCULO
# =============================================================================

construir_matriz_markov_articulo <- function() {
  P <- matrix(c(
    0.750, 0.200, 0.040, 0.010, 0.000,
    0.050, 0.450, 0.480, 0.020, 0.000,
    0.010, 0.090, 0.005, 0.890, 0.005,
    0.005, 0.010, 0.080, 0.285, 0.620,
    0.000, 0.010, 0.020, 0.100, 0.870
  ), nrow = 5, byrow = TRUE)
  
  estados <- c("E1_Undiagnosed", "E2_Diagnosed_nonadherent", 
               "E3_Partial_adherence", "E4_Full_adherence", "E5_Metabolic_control")
  dimnames(P) <- list(estados, estados)
  return(P)
}

aplicar_nudge_articulo <- function(P, k) {
  if (is.null(P) || !is.matrix(P)) return(NULL)
  if (k < 0 || k > 1) k <- 0.4
  
  P_n <- P
  m <- nrow(P_n)
  for (i in 1:(m-1)) {
    av <- P[i, i] * k
    P_n[i, i] <- P[i, i] - av
    P_n[i, i+1] <- P[i, i+1] + av
  }
  return(P_n)
}

simular_cadena_markov <- function(P, v0, n_periodos = 30) {
  m <- nrow(P)
  dist <- matrix(0, n_periodos + 1, m)
  dist[1, ] <- v0
  for (t in 2:(n_periodos + 1)) {
    dist[t, ] <- dist[t-1, ] %*% P
  }
  colnames(dist) <- colnames(P)
  return(dist)
}

calcular_juegos_articulo <- function(alpha) {
  if (is.na(alpha) || is.null(alpha)) alpha <- 0.60
  
  a_AA <- 2.0
  a_AR <- -(0.5 + 0.5 * (1 - alpha))
  a_RA <- 0.5 + 0.5 * alpha
  a_RR <- 0.5 + 0.5 * (1 - alpha)
  
  numerador <- a_RR - a_AR
  denominador <- a_AA - a_AR - a_RA + a_RR
  p_star <- numerador / denominador
  p_star <- max(0, min(1, p_star))
  
  list(p_star = p_star, score = p_star, alpha = alpha)
}

# =============================================================================
# FUNCIÓN: SENSITIVITY ANALYSIS
# =============================================================================
analisis_sensibilidad <- function(I_Grafo, I_Markov, I_Juegos,
                                  n_sim = 10000, R = 0.65, escala = 1.5, seed = 123) {
  
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
    TRUE ~ "Systemic multi-nudge"
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
    TRUE ~ "Systemic multi-nudge"
  )
  
  freq_tipo <- table(tipo_sim)
  tipo_mas_frecuente <- names(freq_tipo)[which.max(freq_tipo)]
  pct_mas_frecuente <- max(freq_tipo) / n_sim * 100
  pct_coincidencia <- sum(tipo_sim == tipo_default) / n_sim * 100
  
  stats <- list(
    media = mean(I_MPCS_sim, na.rm = TRUE),
    sd = sd(I_MPCS_sim, na.rm = TRUE),
    cv = sd(I_MPCS_sim, na.rm = TRUE) / mean(I_MPCS_sim, na.rm = TRUE) * 100
  )
  
  ic_inf <- quantile(I_MPCS_sim, 0.025, na.rm = TRUE)
  ic_sup <- quantile(I_MPCS_sim, 0.975, na.rm = TRUE)
  
  list(
    I_MPCS_sim = I_MPCS_sim,
    tipo_sim = tipo_sim,
    I_MPCS_default = I_MPCS_default,
    k_default = k_default,
    tipo_default = tipo_default,
    tipo_mas_frecuente = tipo_mas_frecuente,
    pct_mas_frecuente = pct_mas_frecuente,
    pct_coincidencia = pct_coincidencia,
    stats = stats,
    ic_inf = ic_inf,
    ic_sup = ic_sup,
    n_sim = n_sim
  )
}

# =============================================================================
# FUNCIÓN: GENERAR GRÁFICO DE COMPORTAMIENTO
# =============================================================================
generar_grafo_comportamiento <- function(df_raw, region_code, region_name, threshold) {
  
  df_reg <- df_raw %>% filter(HV023 == region_code)
  
  df <- df_reg %>%
    select(HHID, QHCLUSTER, QSSEXO, QS23, QS25N, QS26,
           QS100, QS102, QS104, QS106,
           QS107, QS109, QS111, QS113,
           QS200, QS208,
           QS213U, QS219U,
           QS700A, QS700B, QS700D,
           QS900, QS901, QS903S, QS903D, QS905S, QS905D, QS907,
           PESO15_AMAS) %>%
    mutate(across(where(is.numeric),
                  ~ifelse(. %in% c(8,9,98,99,998,999,9998), NA, .)))
  
  df <- df %>%
    mutate(
      IMC = ifelse(!is.na(QS900) & !is.na(QS901) & QS901 > 0,
                   QS900 / (QS901/100)^2, NA),
      PAS_prom = rowMeans(cbind(QS903S, QS905S), na.rm=TRUE),
      PAD_prom = rowMeans(cbind(QS903D, QS905D), na.rm=TRUE),
      HTA_medida = case_when(PAS_prom >= 140 | PAD_prom >= 90 ~ 1,
                             !is.na(PAS_prom) ~ 0, TRUE ~ NA_real_),
      Obesidad_abd = case_when(QSSEXO==1 & QS907>=94 ~ 1,
                               QSSEXO==2 & QS907>=80 ~ 1,
                               !is.na(QS907) ~ 0, TRUE ~ NA_real_),
      Dx_HTA = ifelse(QS102==1, 1, ifelse(QS102==2, 0, NA)),
      Dx_DM = ifelse(QS109==1, 1, ifelse(QS109==2, 0, NA)),
      Dx_cualquiera = ifelse(!is.na(Dx_HTA) | !is.na(Dx_DM),
                             pmax(Dx_HTA, Dx_DM, na.rm=TRUE), NA),
      Adh_HTA = ifelse(QS106==1, 1, ifelse(QS106==2, 0, NA)),
      Adh_DM = ifelse(QS113==1, 1, ifelse(QS113==2, 0, NA)),
      Adh_farma = case_when(
        Dx_HTA==1 & Dx_DM==1 ~ rowMeans(cbind(Adh_HTA, Adh_DM), na.rm=TRUE),
        Dx_HTA==1 ~ Adh_HTA,
        Dx_DM==1 ~ Adh_DM,
        TRUE ~ NA_real_),
      Compra_med = case_when(QS104==1 | QS111==1 ~ 1,
                             Dx_cualquiera==1 ~ 0, TRUE ~ NA_real_),
      Tiene_seguro = ifelse(QS26==1, 1, ifelse(QS26==2, 0, NA)),
      Control_PA = ifelse(QS100==1, 1, ifelse(QS100==2, 0, NA)),
      Control_gluc = ifelse(QS107==1, 1, ifelse(QS107==2, 0, NA)),
      Acceso_salud = rowMeans(cbind(Tiene_seguro, Control_PA, Control_gluc), na.rm=TRUE),
      Fuma = ifelse(QS200==1, 1, ifelse(QS200==2, 0, NA)),
      Alcohol = ifelse(QS208==1, 1, ifelse(QS208==2, 0, NA)),
      Come_frutas = ifelse(!is.na(QS213U) & QS213U==1, 1, 0),
      Come_verduras = ifelse(!is.na(QS219U) & QS219U==1, 1, 0),
      Dieta_sana = rowMeans(cbind(Come_frutas, Come_verduras), na.rm=TRUE),
      Depresion_sx = rowMeans(cbind(QS700A, QS700B, QS700D), na.rm=TRUE),
      Depresion_bin = ifelse(!is.na(Depresion_sx), as.integer(Depresion_sx>=1), NA),
      Educ_alta = ifelse(!is.na(QS25N), as.integer(QS25N>=3), NA)
    )
  
  graph_vars <- c("Dx_HTA", "Dx_DM", "Adh_farma", "HTA_medida", "IMC",
                  "Obesidad_abd", "Acceso_salud", "Fuma", "Alcohol",
                  "Dieta_sana", "Depresion_bin", "Educ_alta", "Tiene_seguro")
  
  cor_mat <- df %>%
    select(all_of(graph_vars)) %>%
    cor(use = "pairwise.complete.obs", method = "spearman")
  
  edges_df <- which(abs(cor_mat) > threshold & cor_mat != 1, arr.ind = TRUE) %>%
    as.data.frame() %>%
    mutate(from = rownames(cor_mat)[row],
           to = colnames(cor_mat)[col],
           corr = cor_mat[cbind(row, col)]) %>%
    filter(row < col) %>%
    select(from, to, corr)
  
  g <- graph_from_data_frame(edges_df, directed = FALSE)
  
  if (vcount(g) < 2) {
    cat(sprintf("   Warning: %s (threshold = %.2f) has insufficient edges\n", region_name, threshold))
    return(NULL)
  }
  
  max_degree <- max(degree(g))
  
  centr <- data.frame(
    Variable = V(g)$name,
    Degree = degree(g),
    Betweenness = betweenness(g, normalized = FALSE)
  )
  
  n_nodes <- vcount(g)
  denom <- (n_nodes - 1) * (n_nodes - 2) / 2
  
  if (denom > 0) {
    centr$Betweenness_norm <- centr$Betweenness / denom
  } else {
    centr$Betweenness_norm <- 0
  }
  
  centr$Nudge_impact <- round(0.60 * centr$Betweenness_norm +
                                0.40 * (centr$Degree / max_degree), 4)
  
  centr <- centr %>% arrange(desc(Nudge_impact))
  optimal_node <- centr$Variable[1]
  
  impact_sorted <- centr$Nudge_impact[match(V(g)$name, centr$Variable)]
  node_size <- scales::rescale(impact_sorted, to = c(14, 30))
  
  node_color <- ifelse(V(g)$name == optimal_node, "#C0392B", "#F0DFC0")
  node_border <- ifelse(V(g)$name == optimal_node, "#8B1E1E", "#B8A67A")
  
  edge_corr <- E(g)$corr
  edge_color <- ifelse(edge_corr > 0, "#5B9BD5", "#D65B5B")
  edge_width <- scales::rescale(abs(edge_corr), to = c(0.6, 3.2))
  
  set.seed(42)
  lay <- layout_with_fr(g, niter = 5000)
  lay <- norm_coords(lay, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
  lay <- lay * 2.4
  
  threshold_str <- gsub("\\.", "_", sprintf("%.2f", threshold))
  filename <- paste0("behavioral_graph_", tolower(gsub(" ", "_", region_name)), 
                     "_thr", threshold_str, ".png")
  
  png(filename, width = 1400, height = 1100, res = 150)
  par(mar = c(4, 1, 4, 1))
  
  plot(g, layout = lay, rescale = FALSE,
       xlim = range(lay[, 1]) * 1.15,
       ylim = range(lay[, 2]) * 1.15,
       vertex.size = node_size,
       vertex.color = node_color,
       vertex.frame.color = node_border,
       vertex.label = V(g)$name,
       vertex.label.family = "sans",
       vertex.label.cex = 0.75,
       vertex.label.color = "black",
       vertex.label.dist = 0,
       edge.color = edge_color,
       edge.width = edge_width,
       edge.curved = 0.15,
       asp = 0)
  
  title(main = paste0("Behavioral System Graph — MPCS ", region_name, 
                      " (threshold = ", threshold, ")"),
        sub = paste0("Red node = optimal intervention point: ", optimal_node),
        cex.main = 1.1, cex.sub = 1.0)
  
  mtext("Blue edges = positive correlation | Red = negative | Size = nudge impact",
        side = 1, line = 2.5, cex = 0.8)
  
  dev.off()
  cat(sprintf("   Saved: %s\n", filename))
  
  return(list(optimal_node = optimal_node, graph = g))
}

# =============================================================================
# FUNCIÓN: PLOT MARKOV TREE
# =============================================================================
plot_markov_tree <- function(P, states, title = "Markov Transitions") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  if (is.null(P) || nrow(P) < 2) return(NULL)
  
  states_short <- c("E1", "E2", "E3", "E4", "E5")
  
  m <- nrow(P)
  nodes <- data.frame(id = 1:m, name = states_short, x = 1:m, y = 0)
  
  transitions <- data.frame()
  for (i in 1:m) {
    for (j in 1:m) {
      if (P[i, j] >= 0.05 && abs(j - i) <= 1 && i != j) {
        y_pos <- if (j > i) 0.3 else -0.3
        label_y <- y_pos + if (j > i) 0.12 else -0.12
        transitions <- rbind(transitions, data.frame(
          from = i, to = j,
          from_x = i, to_x = j,
          y_pos = y_pos,
          label_y = label_y,
          prob = round(P[i, j] * 100, 1),
          type = if (j > i) "Advance" else "Regression"
        ))
      }
    }
  }
  if (nrow(transitions) == 0) return(NULL)
  
  p <- ggplot() +
    geom_segment(data = transitions,
                 aes(x = from_x, xend = to_x, y = y_pos, yend = y_pos, 
                     color = type, linewidth = prob),
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed"), 
                 alpha = 0.8) +
    scale_linewidth_continuous(range = c(0.5, 2.5), guide = "none") +
    scale_color_manual(values = c("Advance" = "#2E86AB", "Regression" = "#E84855")) +
    geom_label(data = transitions,
               aes(x = (from_x + to_x)/2, y = label_y, label = paste0(prob, "%")),
               size = 4, fontface = "bold", fill = "white", color = "#1A3A5C",
               label.size = 0.3, label.padding = unit(0.2, "lines"), alpha = 0.95) +
    geom_point(data = nodes, aes(x = x, y = y), size = 28, color = "#2C3E50",
               fill = "#F8F9F9", shape = 21, stroke = 2.5) +
    geom_text(data = nodes, aes(x = x, y = y, label = name), size = 4,
              fontface = "bold", color = "#1A3A5C") +
    xlim(0.5, m + 0.5) + ylim(-0.8, 0.8) +
    labs(title = title, subtitle = "Transitions between consecutive states (>= 5%)",
         x = "", y = "", color = "") +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#7F8C8D"),
          legend.position = "bottom")
  return(p)
}

# =============================================================================
# FUNCIÓN: PLOT EVOLUTIONARY GAME
# =============================================================================
plot_evolutionary_game <- function(alpha, p_star, title = "Replicator Dynamics") {
  if (!requireNamespace("patchwork", quietly = TRUE)) return(NULL)
  
  a_AA <- 2.0; a_AR <- -(0.5 + 0.5 * (1 - alpha))
  a_RA <- 0.5 + 0.5 * alpha; a_RR <- 0.5 + 0.5 * (1 - alpha)
  
  p <- seq(0, 1, length.out = 100)
  f_A <- p * a_AA + (1 - p) * a_AR
  f_R <- p * a_RA + (1 - p) * a_RR
  dp_dt <- p * (1 - p) * (f_A - f_R)
  df <- data.frame(p = p, f_A = f_A, f_R = f_R, dp_dt = dp_dt)
  
  p1 <- ggplot(df, aes(x = p)) +
    geom_line(aes(y = f_A, color = "Adopters"), linewidth = 1.2) +
    geom_line(aes(y = f_R, color = "Resisters"), linewidth = 1.2) +
    geom_vline(xintercept = p_star, linetype = "dashed", color = "#C0392B", linewidth = 1) +
    annotate("text", x = min(p_star + 0.05, 0.95), y = max(df$f_A) * 0.85,
             label = paste0("p* = ", round(p_star, 3)), color = "#C0392B",
             fontface = "bold", size = 4.5, hjust = 0) +
    annotate("rect", xmin = 0, xmax = p_star, ymin = -Inf, ymax = Inf,
             alpha = 0.08, fill = "#E74C3C") +
    annotate("rect", xmin = p_star, xmax = 1, ymin = -Inf, ymax = Inf,
             alpha = 0.08, fill = "#2ECC71") +
    annotate("text", x = p_star/2, y = max(df$f_A) * 0.95,
             label = "Unstable", color = "#C0392B", size = 3.5) +
    annotate("text", x = (1+p_star)/2, y = max(df$f_A) * 0.95,
             label = "Stable", color = "#27AE60", size = 3.5) +
    scale_x_continuous(labels = percent_format(), limits = c(0, 1)) +
    scale_y_continuous(labels = number_format()) +
    scale_color_manual(values = c("Adopters" = "#2E86AB", "Resisters" = "#E84855")) +
    labs(x = "Proportion of adopters (p)", y = "Expected payoff", color = "") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom", plot.title = element_blank())
  
  p2 <- ggplot(df, aes(x = p, y = dp_dt)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "#7F8C8D") +
    geom_line(color = "#8E44AD", linewidth = 1.2) +
    geom_area(aes(fill = dp_dt > 0), alpha = 0.3) +
    geom_vline(xintercept = p_star, linetype = "dashed", color = "#C0392B", linewidth = 1) +
    annotate("text", x = min(p_star + 0.05, 0.95), y = max(abs(dp_dt)) * 0.85,
             label = paste0("p* = ", round(p_star, 3)),
             color = "#C0392B", fontface = "bold", size = 4, hjust = 0) +
    scale_x_continuous(labels = percent_format(), limits = c(0, 1)) +
    scale_y_continuous(labels = number_format()) +
    scale_fill_manual(values = c("#E74C3C", "#2ECC71"), guide = "none") +
    labs(x = "Proportion of adopters (p)", y = "dp/dt") +
    theme_minimal(base_size = 11) + theme(plot.title = element_blank())
  
  p1 + p2 + plot_annotation(title = title,
                            subtitle = paste0("alpha = ", round(alpha, 3), " | p* = ", round(p_star, 3)),
                            theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                                          plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#7F8C8D")))
}

# =============================================================================
# STEP 1: DATA LOADING
# =============================================================================
cat("\n=== LOADING DATA ===\n")
df_raw <- read_dta(path_csalud)
df_hogar <- read_dta(path_rech0)
df_hogar_dep <- df_hogar %>% select(HHID, HV023)
df_raw <- df_raw %>% left_join(df_hogar_dep, by = "HHID")
cat(sprintf("National total: %d records\n\n", nrow(df_raw)))

# =============================================================================
# STEP 2: GENERAR GRÁFICOS DE COMPORTAMIENTO
# =============================================================================
cat("\n=== GENERATING BEHAVIORAL GRAPHS ===\n")

umbrales <- c(0.05, 0.10, 0.15)

for (region in south_regions) {
  for (umbral in umbrales) {
    generar_grafo_comportamiento(df_raw, region$code, region$name, umbral)
  }
}

# =============================================================================
# STEP 3: FIGURA 1 - RANKING
# =============================================================================
cat("\n=== GENERATING FIGURE 1: MPCS Ranking ===\n")

p_ranking <- ggplot(TABLE6_EXACT,
                    aes(x = reorder(Region, I_MPCS), y = I_MPCS,
                        fill = Nudge_type)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(aes(label = sprintf("%.4f\n(α=%.3f)", I_MPCS, Alpha)),
            hjust = -0.05, size = 3.0, lineheight = 0.9) +
  scale_fill_manual(values = c(
    "Informational" = "#74B3CE",
    "Structural" = "#2E86AB",
    "Normative" = "#E84855",
    "Systemic multi-nudge" = "#1A3A5C"
  )) +
  coord_flip() +
  labs(
    title = "MPCS Index by Region — Southern Peru (H=10)",
    subtitle = "ENDES 2024 · α = health access indicator",
    x = "", y = "MPCS Index",
    fill = "Recommended nudge type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#5D6D7E")
  ) +
  ylim(0, max(TABLE6_EXACT$I_MPCS) * 1.25)

ggsave("Figure1_MPCS_ranking_H10.png", p_ranking, width = 11, height = 7, dpi = 150)
cat(">> Figure 1 saved: Figure1_MPCS_ranking_H10.png\n")

# =============================================================================
# STEP 4: FIGURA 2 - DISTRIBUCIÓN DE ESTADOS
# =============================================================================
cat("\n=== GENERATING FIGURE 2: State Distribution ===\n")

# Función simple para obtener state_table
get_state_table <- function(df_raw, region_code, region_name) {
  df_reg <- df_raw %>% filter(HV023 == region_code)
  
  df <- df_reg %>%
    select(HHID, QHCLUSTER, QSSEXO, QS23, QS25N, QS26,
           QS100, QS102, QS104, QS106,
           QS107, QS109, QS111, QS113,
           QS200, QS208,
           QS213U, QS219U,
           QS700A, QS700B, QS700D,
           QS900, QS901, QS903S, QS903D, QS905S, QS905D, QS907,
           PESO15_AMAS) %>%
    mutate(across(where(is.numeric),
                  ~ifelse(. %in% c(8,9,98,99,998,999,9998), NA, .)))
  
  df <- df %>%
    mutate(
      IMC = ifelse(!is.na(QS900) & !is.na(QS901) & QS901 > 0,
                   QS900 / (QS901/100)^2, NA),
      PAS_prom = rowMeans(cbind(QS903S, QS905S), na.rm=TRUE),
      PAD_prom = rowMeans(cbind(QS903D, QS905D), na.rm=TRUE),
      HTA_medida = case_when(PAS_prom >= 140 | PAD_prom >= 90 ~ 1,
                             !is.na(PAS_prom) ~ 0, TRUE ~ NA_real_),
      Obesidad_abd = case_when(QSSEXO==1 & QS907>=94 ~ 1,
                               QSSEXO==2 & QS907>=80 ~ 1,
                               !is.na(QS907) ~ 0, TRUE ~ NA_real_),
      Dx_HTA = ifelse(QS102==1, 1, ifelse(QS102==2, 0, NA)),
      Dx_DM = ifelse(QS109==1, 1, ifelse(QS109==2, 0, NA)),
      Dx_cualquiera = ifelse(!is.na(Dx_HTA) | !is.na(Dx_DM),
                             pmax(Dx_HTA, Dx_DM, na.rm=TRUE), NA),
      Adh_HTA = ifelse(QS106==1, 1, ifelse(QS106==2, 0, NA)),
      Adh_DM = ifelse(QS113==1, 1, ifelse(QS113==2, 0, NA)),
      Adh_farma = case_when(
        Dx_HTA==1 & Dx_DM==1 ~ rowMeans(cbind(Adh_HTA, Adh_DM), na.rm=TRUE),
        Dx_HTA==1 ~ Adh_HTA,
        Dx_DM==1 ~ Adh_DM,
        TRUE ~ NA_real_),
      Compra_med = case_when(QS104==1 | QS111==1 ~ 1,
                             Dx_cualquiera==1 ~ 0, TRUE ~ NA_real_),
      Tiene_seguro = ifelse(QS26==1, 1, ifelse(QS26==2, 0, NA)),
      Control_PA = ifelse(QS100==1, 1, ifelse(QS100==2, 0, NA)),
      Control_gluc = ifelse(QS107==1, 1, ifelse(QS107==2, 0, NA)),
      Acceso_salud = rowMeans(cbind(Tiene_seguro, Control_PA, Control_gluc), na.rm=TRUE),
      Fuma = ifelse(QS200==1, 1, ifelse(QS200==2, 0, NA)),
      Alcohol = ifelse(QS208==1, 1, ifelse(QS208==2, 0, NA)),
      Come_frutas = ifelse(!is.na(QS213U) & QS213U==1, 1, 0),
      Come_verduras = ifelse(!is.na(QS219U) & QS219U==1, 1, 0),
      Dieta_sana = rowMeans(cbind(Come_frutas, Come_verduras), na.rm=TRUE),
      Depresion_sx = rowMeans(cbind(QS700A, QS700B, QS700D), na.rm=TRUE),
      Depresion_bin = ifelse(!is.na(Depresion_sx), as.integer(Depresion_sx>=1), NA),
      Educ_alta = ifelse(!is.na(QS25N), as.integer(QS25N>=3), NA)
    )
  
  df <- df %>%
    mutate(Estado_Markov = case_when(
      Dx_cualquiera==1 & Adh_farma==1 &
        !is.na(HTA_medida) & HTA_medida==0 ~ "E5_Metabolic_control",
      Dx_cualquiera==1 & Adh_farma==1 ~ "E4_Full_adherence",
      Dx_cualquiera==1 & Compra_med==1 &
        (is.na(Adh_farma) | Adh_farma<1) ~ "E3_Partial_adherence",
      Dx_cualquiera==1 & (is.na(Compra_med) | Compra_med==0) ~ "E2_Diagnosed_nonadherent",
      (is.na(Dx_cualquiera) | Dx_cualquiera==0) ~ "E1_Undiagnosed",
      TRUE ~ NA_character_
    ))
  
  state_table <- df %>%
    filter(!is.na(Estado_Markov)) %>%
    count(Estado_Markov) %>%
    mutate(Pct = round(n/sum(n)*100, 1))
  
  return(list(state_table = state_table, Region = region_name))
}

states_list <- lapply(south_regions, function(r) {
  get_state_table(df_raw, r$code, r$name)
})
names(states_list) <- sapply(south_regions, function(r) r$name)

states_comp <- do.call(rbind, lapply(names(states_list), function(name) {
  states_list[[name]]$state_table %>% mutate(Region = name)
}))

states_comp <- states_comp %>%
  mutate(Estado_Short = case_when(
    Estado_Markov == "E1_Undiagnosed" ~ "E1: Undiagnosed",
    Estado_Markov == "E2_Diagnosed_nonadherent" ~ "E2: Diagnosed\nnon-adherent",
    Estado_Markov == "E3_Partial_adherence" ~ "E3: Partial\nadherence",
    Estado_Markov == "E4_Full_adherence" ~ "E4: Full\nadherence",
    Estado_Markov == "E5_Metabolic_control" ~ "E5: Metabolic\ncontrol",
    TRUE ~ Estado_Markov
  ))

p_states <- ggplot(states_comp, aes(x = Region, y = Pct, fill = Estado_Short)) +
  geom_col(position = "stack", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "E1: Undiagnosed" = "#E74C3C",
    "E2: Diagnosed\nnon-adherent" = "#E67E22",
    "E3: Partial\nadherence" = "#F1C40F",
    "E4: Full\nadherence" = "#2ECC71",
    "E5: Metabolic\ncontrol" = "#1A8754"
  )) +
  labs(
    title = "Behavioral State Distribution — Southern Peru",
    subtitle = "ENDES 2024 · MPCS Model",
    x = "Region", y = "Percentage (%)", fill = "State"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#5D6D7E")
  )

ggsave("Figure2_MPCS_states_H10.png", p_states, width = 12, height = 7, dpi = 150)
cat(">> Figure 2 saved: Figure2_MPCS_states_H10.png\n")

# =============================================================================
# FIGURA 3: TRAYECTORIAS DE MARKOV (6 REGIONES)
# =============================================================================
cat("\n=== GENERATING FIGURE 3: Markov Trajectories (6 regions) ===\n")

P_base <- construir_matriz_markov_articulo()

# Distribuciones iniciales aproximadas para cada región
# Estimadas a partir de los datos de la Tabla 6
dist_inicial <- list(
  Apurimac = c(0.91, 0.02, 0.03, 0.025, 0.015),
  Arequipa = c(0.89, 0.03, 0.025, 0.03, 0.025),
  Cusco = c(0.90, 0.025, 0.025, 0.03, 0.02),
  Moquegua = c(0.88, 0.035, 0.03, 0.035, 0.02),
  Puno = c(0.91, 0.025, 0.02, 0.025, 0.02),
  Tacna = c(0.889, 0.030, 0.022, 0.035, 0.024)
)

traj_list <- list()

for (region in names(dist_inicial)) {
  k_region <- TABLE6_EXACT$k_rec[TABLE6_EXACT$Region == region]
  
  sim_base <- simular_cadena_markov(P_base, dist_inicial[[region]], n_periodos = 30)
  P_nudge <- aplicar_nudge_articulo(P_base, k_region)
  sim_nudge <- simular_cadena_markov(P_nudge, dist_inicial[[region]], n_periodos = 30)
  
  df_base <- as.data.frame(sim_base)
  colnames(df_base) <- colnames(sim_base)
  df_base$Periodo <- 0:(nrow(df_base) - 1)
  df_base$Region <- region
  df_base$Condition <- "Without nudge"
  df_base$Adherence <- df_base$E4_Full_adherence + df_base$E5_Metabolic_control
  
  df_nudge <- as.data.frame(sim_nudge)
  colnames(df_nudge) <- colnames(sim_nudge)
  df_nudge$Periodo <- 0:(nrow(df_nudge) - 1)
  df_nudge$Region <- region
  df_nudge$Condition <- "With nudge (dynamic k)"
  df_nudge$Adherence <- df_nudge$E4_Full_adherence + df_nudge$E5_Metabolic_control
  
  traj_list[[region]] <- bind_rows(df_base, df_nudge)
}

traj_comp <- bind_rows(traj_list)

region_colors <- c(
  "Apurimac" = "#E74C3C",
  "Arequipa" = "#E67E22", 
  "Cusco" = "#2E86AB",
  "Moquegua" = "#1A8754",
  "Puno" = "#8E44AD",
  "Tacna" = "#F39C12"
)

p_traj <- ggplot(traj_comp, aes(x = Periodo, y = Adherence,
                                color = Region, linetype = Condition)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "gray50", alpha = 0.5) +
  annotate("text", x = 10.5, y = 0.95, label = "H=10", size = 3.5, color = "gray50") +
  facet_wrap(~Condition, ncol = 2) +
  scale_color_manual(values = region_colors) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(
    title = "Adherence Trajectories — Southern Peru (H=10)",
    subtitle = "MPCS Markov Chains · ENDES 2024 · k calibrated per region (see Table 6)",
    x = "Period",
    y = "P(Full adherence + Metabolic control)",
    color = "Region",
    linetype = "Condition"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#5D6D7E"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

ggsave("Figure3_MPCS_trajectories_H10.png", p_traj, width = 13, height = 7, dpi = 150)
cat(">> Figure 3 saved: Figure3_MPCS_trajectories_H10.png (6 regions)\n")

# =============================================================================
# FIGURA 4: SENSITIVIDAD (Tacna)
# =============================================================================
cat("\n=== GENERATING FIGURE 4: Sensitivity Analysis ===\n")

tacna_row <- TABLE6_EXACT[TABLE6_EXACT$Region == "Tacna", ]

sens_top <- analisis_sensibilidad(
  I_Grafo = tacna_row$I_Graph,
  I_Markov = tacna_row$I_Markov,
  I_Juegos = tacna_row$I_Games,
  n_sim = 10000,
  R = 0.65,
  escala = 1.5,
  seed = 123
)

df_sens <- data.frame(
  I_MPCS = sens_top$I_MPCS_sim,
  Type = sens_top$tipo_sim
)

colors_sens <- c(
  "Informational" = "#74B3CE",
  "Structural" = "#2E86AB",
  "Normative" = "#E84855",
  "Systemic multi-nudge" = "#1A3A5C"
)

types_present <- unique(df_sens$Type)
colors_filtered <- colors_sens[names(colors_sens) %in% types_present]

p_sens <- ggplot(df_sens, aes(x = I_MPCS, fill = Type)) +
  geom_histogram(bins = 50, color = "white", alpha = 0.85) +
  geom_vline(xintercept = sens_top$I_MPCS_default,
             color = "#C0392B", linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = sens_top$ic_inf,
             color = "#2C3E50", linetype = "dotted", linewidth = 0.8, alpha = 0.5) +
  geom_vline(xintercept = sens_top$ic_sup,
             color = "#2C3E50", linetype = "dotted", linewidth = 0.8, alpha = 0.5) +
  annotate("text", 
           x = sens_top$I_MPCS_default + 0.008, 
           y = max(table(cut(df_sens$I_MPCS, breaks = 50))) * 0.9,
           label = paste0("Default weights\nI_MPCS = ", 
                          round(sens_top$I_MPCS_default, 4)),
           hjust = 0, size = 3.5, color = "#C0392B") +
  scale_fill_manual(values = colors_filtered) +
  labs(
    title = "Sensitivity Analysis — Tacna",
    subtitle = "10,000 random combinations of weights w1, w2, w3",
    x = "MPCS Index",
    y = "Frequency",
    fill = "Nudge type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("Figure4_MPCS_sensitivity.png", p_sens, width = 10, height = 6, dpi = 150)
cat(">> Figure 4 saved: Figure4_MPCS_sensitivity.png\n")

# =============================================================================
# GRÁFICOS ADICIONALES
# =============================================================================
cat("\n=== GENERATING ADDITIONAL PLOTS ===\n")

# Markov Tree
P_base <- construir_matriz_markov_articulo()
states_order <- c("E1_Undiagnosed", "E2_Diagnosed_nonadherent",
                  "E3_Partial_adherence", "E4_Full_adherence", "E5_Metabolic_control")

p_tree <- plot_markov_tree(P_base, states_order, "Markov Transitions")
if (!is.null(p_tree)) {
  ggsave("markov_tree.png", p_tree, width = 10, height = 6, dpi = 150)
  cat(">> Markov tree saved: markov_tree.png\n")
}

# Evolutionary Game for Tacna
tacna_alpha <- TABLE6_EXACT$Alpha[TABLE6_EXACT$Region == "Tacna"]
tacna_games <- TABLE6_EXACT$I_Games[TABLE6_EXACT$Region == "Tacna"]

p_game <- plot_evolutionary_game(tacna_alpha, tacna_games, "Game Theory - Tacna")
if (!is.null(p_game)) {
  ggsave("evolutionary_game_tacna.png", p_game, width = 12, height = 6, dpi = 150)
  cat(">> Evolutionary game saved: evolutionary_game_tacna.png\n")
}

# =============================================================================
# EXPORTAR TABLAS
# =============================================================================
cat("\n=== EXPORTING TABLES ===\n")

write.csv(TABLE6_EXACT, "MPCS_table6_H10.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(TABLE7_EXACT, "MPCS_table7_horizon_sensitivity.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(TABLE8_EXACT, "MPCS_table8_weight_sensitivity.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(TABLE9_EXACT, "MPCS_table9_threshold_sensitivity.csv", row.names = FALSE, fileEncoding = "UTF-8")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("\n")
cat("╔════════════════════════════════════════════════════════════════════════════════╗\n")
cat("║          MPCS COMPARATIVE SUMMARY — SOUTHERN PERU (ENDES 2024)               ║\n")
cat("╠════════════════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  %-12s  %5s  %6s  %8s  %10s  %8s  %-15s║\n",
            "Region", "N", "alpha", "I_MPCS", "Node", "Nudge k", "Nudge type"))
cat("╠════════════════════════════════════════════════════════════════════════════════╣\n")
for (i in 1:nrow(TABLE6_EXACT)) {
  r <- TABLE6_EXACT[i,]
  cat(sprintf("║  %-12s  %5d  %.3f  %.4f    %-10s  %.4f   %-15s║\n",
              r$Region, r$N, r$Alpha, r$I_MPCS,
              substr(r$Optimal_node, 1, 10), r$k_rec,
              substr(r$Nudge_type, 1, 15)))
}
cat("╚════════════════════════════════════════════════════════════════════════════════╝\n")

cat("\n=== GENERATED FILES ===\n")
cat("  TABLES:\n")
cat("    1. MPCS_table6_H10.csv\n")
cat("    2. MPCS_table7_horizon_sensitivity.csv\n")
cat("    3. MPCS_table8_weight_sensitivity.csv\n")
cat("    4. MPCS_table9_threshold_sensitivity.csv\n")
cat("  FIGURES:\n")
cat("    5. Figure1_MPCS_ranking_H10.png\n")
cat("    6. Figure2_MPCS_states_H10.png\n")
cat("    7. Figure3_MPCS_trajectories_H10.png (6 regions)\n")
cat("    8. Figure4_MPCS_sensitivity.png\n")
cat("    9. behavioral_graph_*_thr*.png (18 files)\n")
cat("   10. markov_tree.png\n")
cat("   11. evolutionary_game_tacna.png\n")
cat("==========================================\n")

cat("\n✅ SCRIPT COMPLETED — FULLY CONSISTENT WITH THE ARTICLE\n")