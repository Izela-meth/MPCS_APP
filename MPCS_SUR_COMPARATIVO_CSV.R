# =============================================================================
# MPCS — COMPARATIVE ANALYSIS: SOUTHERN PERU (H = 8, 10, 12)
# =============================================================================
# Regions: Apurimac, Arequipa, Cusco, Moquegua, Puno, Tacna
# Source: DEMO DATA (synthetic CSV mimicking ENDES 2024 structure)
#
# GENERALIZABLE VERSION — recalculates I_Markov from empirical distributions
# =============================================================================
# NOTE: This script recalculates all modular indices (I_Graph, I_Markov, 
# I_Games) from the synthetic data. Outputs may differ from published values
# by approximately 0.001. For bit-identical reproduction of Tables 6–10, 
# use the Hardcode version.
#
# INSTRUCTIONS FOR REVIEWERS:
# 1. Place this script and demo_data.csv in the same folder.
# 2. Run the entire script in R (RStudio or command line).
# 3. All tables and figures will be generated automatically.
# =============================================================================


# --- PACKAGES ----------------------------------------------------------------
packages <- c("dplyr", "tidyr", "ggplot2", "igraph",
              "markovchain", "reshape2", "scales", "gridExtra", "patchwork")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(dplyr); library(tidyr)
library(ggplot2); library(igraph); library(scales); library(gridExtra)
library(patchwork)

# =============================================================================
# CONFIGURATION
# =============================================================================
ruta_demo <- "demo_data.csv"

if (!file.exists(ruta_demo)) {
  cat("demo_data.csv not found in working directory.\n")
  cat("Please ensure demo_data.csv is in the same folder as this script.\n")
  stop("demo_data.csv missing")
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

# Default weights and parameters
W_DEFAULT <- c(w1 = 0.35, w2 = 0.40, w3 = 0.25)
R_DEFAULT <- 0.65
ESCALA_DEFAULT <- 1.5

# =============================================================================
# CORE FUNCTIONS
# =============================================================================

construir_matriz_markov <- function() {
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
  P <- P / rowSums(P)
  return(P)
}

aplicar_nudge <- function(P, k) {
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

calcular_I_Markov <- function(v0, P, H = 10) {
  sim <- simular_cadena_markov(P, v0, n_periodos = H)
  dist_h <- sim[H + 1, ]
  I_Markov <- dist_h["E4_Full_adherence"] + dist_h["E5_Metabolic_control"]
  return(max(0, min(1, I_Markov)))
}

calcular_juegos <- function(alpha) {
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
# DATA PREPARATION (from CSV)
# =============================================================================

cat("\n=== LOADING AND PREPARING DATA ===\n")
df_raw <- read.csv(ruta_demo, stringsAsFactors = FALSE, encoding = "UTF-8")
df_raw <- df_raw %>% mutate(across(-HHID, ~as.numeric(as.character(.))))
cat(sprintf("Total records: %d\n\n", nrow(df_raw)))

# --- Derive all behavioral variables ---
df_proc <- df_raw %>%
  mutate(
    IMC = ifelse(!is.na(QS900) & !is.na(QS901) & QS901 > 0,
                 QS900 / (QS901/100)^2, NA),
    PAS_prom = rowMeans(cbind(QS903S, QS905S), na.rm = TRUE),
    PAD_prom = rowMeans(cbind(QS903D, QS905D), na.rm = TRUE),
    HTA_medida = case_when(PAS_prom >= 140 | PAD_prom >= 90 ~ 1,
                           !is.na(PAS_prom) ~ 0, TRUE ~ NA_real_),
    Obesidad_abd = case_when(QSSEXO == 1 & QS907 >= 94 ~ 1,
                             QSSEXO == 2 & QS907 >= 80 ~ 1,
                             !is.na(QS907) ~ 0, TRUE ~ NA_real_),
    Dx_HTA = ifelse(QS102 == 1, 1, ifelse(QS102 == 2, 0, NA)),
    Dx_DM = ifelse(QS109 == 1, 1, ifelse(QS109 == 2, 0, NA)),
    Dx_cualquiera = ifelse(!is.na(Dx_HTA) | !is.na(Dx_DM),
                           pmax(Dx_HTA, Dx_DM, na.rm = TRUE), NA),
    Adh_HTA = ifelse(QS106 == 1, 1, ifelse(QS106 == 2, 0, NA)),
    Adh_DM = ifelse(QS113 == 1, 1, ifelse(QS113 == 2, 0, NA)),
    Adh_farma = case_when(
      Dx_HTA == 1 & Dx_DM == 1 ~ rowMeans(cbind(Adh_HTA, Adh_DM), na.rm = TRUE),
      Dx_HTA == 1 ~ Adh_HTA,
      Dx_DM == 1 ~ Adh_DM,
      TRUE ~ NA_real_
    ),
    Compra_med = case_when(QS104 == 1 | QS111 == 1 ~ 1,
                           Dx_cualquiera == 1 ~ 0, TRUE ~ NA_real_),
    Tiene_seguro = ifelse(QS26 == 1, 1, ifelse(QS26 == 2, 0, NA)),
    Control_PA = ifelse(QS100 == 1, 1, ifelse(QS100 == 2, 0, NA)),
    Control_gluc = ifelse(QS107 == 1, 1, ifelse(QS107 == 2, 0, NA)),
    Acceso_salud = rowMeans(cbind(Tiene_seguro, Control_PA, Control_gluc), na.rm = TRUE),
    Fuma = ifelse(QS200 == 1, 1, ifelse(QS200 == 2, 0, NA)),
    Alcohol = ifelse(QS208 == 1, 1, ifelse(QS208 == 2, 0, NA)),
    Come_frutas = ifelse(!is.na(QS213U) & QS213U == 1, 1, 0),
    Come_verduras = ifelse(!is.na(QS219U) & QS219U == 1, 1, 0),
    Dieta_sana = rowMeans(cbind(Come_frutas, Come_verduras), na.rm = TRUE),
    Depresion_sx = rowMeans(cbind(QS700A, QS700B, QS700D), na.rm = TRUE),
    Depresion_bin = ifelse(!is.na(Depresion_sx), as.integer(Depresion_sx >= 1), NA),
    Educ_alta = ifelse(!is.na(QS25N), as.integer(QS25N >= 3), NA)
  ) %>%
  mutate(Estado_Markov = case_when(
    Dx_cualquiera == 1 & Adh_farma == 1 & !is.na(HTA_medida) & HTA_medida == 0 ~ "E5_Metabolic_control",
    Dx_cualquiera == 1 & Adh_farma == 1 ~ "E4_Full_adherence",
    Dx_cualquiera == 1 & Compra_med == 1 & (is.na(Adh_farma) | Adh_farma < 1) ~ "E3_Partial_adherence",
    Dx_cualquiera == 1 & (is.na(Compra_med) | Compra_med == 0) ~ "E2_Diagnosed_nonadherent",
    (is.na(Dx_cualquiera) | Dx_cualquiera == 0) ~ "E1_Undiagnosed",
    TRUE ~ NA_character_
  ))

# =============================================================================
# MODULE 1: GRAPH THEORY (per region)
# =============================================================================

calcular_grafo_region <- function(df_reg, threshold = 0.10) {
  graph_vars <- c("Dx_HTA", "Dx_DM", "Adh_farma", "HTA_medida", "IMC",
                  "Obesidad_abd", "Acceso_salud", "Fuma", "Alcohol",
                  "Dieta_sana", "Depresion_bin", "Educ_alta", "Tiene_seguro")

  cor_mat <- df_reg %>%
    select(all_of(graph_vars)) %>%
    cor(use = "pairwise.complete.obs", method = "spearman")

  edges_df <- which(abs(cor_mat) > threshold & cor_mat != 1, arr.ind = TRUE) %>%
    as.data.frame() %>%
    mutate(from = rownames(cor_mat)[row],
           to = colnames(cor_mat)[col],
           corr = cor_mat[cbind(row, col)]) %>%
    filter(row < col) %>%
    select(from, to, corr)

  if (nrow(edges_df) < 1) {
    return(list(I_Graph = NA, optimal_node = NA, graph = NULL, centr = NULL))
  }

  g <- graph_from_data_frame(edges_df, directed = FALSE)
  n_nodes <- vcount(g)
  if (n_nodes < 2) {
    return(list(I_Graph = NA, optimal_node = NA, graph = g, centr = NULL))
  }

  max_degree <- max(degree(g))
  denom <- (n_nodes - 1) * (n_nodes - 2) / 2

  centr <- data.frame(
    Variable = V(g)$name,
    Degree = degree(g),
    Betweenness = betweenness(g, normalized = FALSE)
  )

  if (denom > 0) {
    centr$Betweenness_norm <- centr$Betweenness / denom
  } else {
    centr$Betweenness_norm <- 0
  }

  centr$Nudge_impact <- 0.60 * centr$Betweenness_norm +
    0.40 * (centr$Degree / max_degree)
  centr <- centr %>% arrange(desc(Nudge_impact))

  optimal_node <- centr$Variable[1]
  I_Graph <- centr$Nudge_impact[1]

  list(I_Graph = I_Graph, optimal_node = optimal_node, graph = g, centr = centr)
}

# =============================================================================
# MODULE 2: MARKOV CHAINS (per region)
# =============================================================================

calcular_markov_region <- function(df_reg, P, H = 10) {
  state_table <- df_reg %>%
    filter(!is.na(Estado_Markov)) %>%
    count(Estado_Markov) %>%
    mutate(prop = n / sum(n))

  estados <- c("E1_Undiagnosed", "E2_Diagnosed_nonadherent",
               "E3_Partial_adherence", "E4_Full_adherence", "E5_Metabolic_control")

  v0 <- sapply(estados, function(e) {
    p <- state_table %>% filter(Estado_Markov == e) %>% pull(prop)
    ifelse(length(p) == 0, 0, p)
  })
  v0 <- as.numeric(v0)
  names(v0) <- estados

  I_Markov <- calcular_I_Markov(v0, P, H)
  list(I_Markov = I_Markov, v0 = v0, state_table = state_table)
}

# =============================================================================
# MODULE 3: EVOLUTIONARY GAME THEORY (per region)
# =============================================================================

calcular_games_region <- function(df_reg) {
  alpha <- mean(df_reg$Acceso_salud, na.rm = TRUE)
  alpha <- max(0, min(1, alpha))
  res <- calcular_juegos(alpha)
  list(I_Games = res$p_star, alpha = alpha)
}

# =============================================================================
# MODULE 4: INTEGRATION
# =============================================================================

calcular_I_MPCS <- function(I_Graph, I_Markov, I_Games, 
                            w = W_DEFAULT, R = R_DEFAULT, escala = ESCALA_DEFAULT) {
  I_MPCS <- w[1] * I_Graph + w[2] * I_Markov + w[3] * I_Games
  k <- min(1, I_MPCS * R * escala)
  tipo <- case_when(
    k < 0.25 ~ "Informational",
    k < 0.50 ~ "Structural",
    k < 0.75 ~ "Normative",
    TRUE ~ "Systemic multi-nudge"
  )
  list(I_MPCS = I_MPCS, k = k, tipo = tipo)
}

# =============================================================================
# MAIN COMPUTATION: ALL REGIONS
# =============================================================================

cat("\n=== COMPUTING MPCS FOR ALL REGIONS ===\n")

P_base <- construir_matriz_markov()

results_list <- list()

for (region in south_regions) {
  reg_name <- region$name
  reg_code <- region$code

  cat(sprintf("Processing %s...\n", reg_name))

  df_reg <- df_proc %>% filter(HV023 == reg_code)
  n_reg <- nrow(df_reg)

  # Module 1: Graph
  grafo_res <- calcular_grafo_region(df_reg, threshold = 0.10)
  I_Graph <- grafo_res$I_Graph
  opt_node <- grafo_res$optimal_node

  # Module 2: Markov
  markov_res <- calcular_markov_region(df_reg, P_base, H = 10)
  I_Markov <- markov_res$I_Markov
  v0_reg <- markov_res$v0

  # Module 3: Games
  games_res <- calcular_games_region(df_reg)
  I_Games <- games_res$I_Games
  alpha_reg <- games_res$alpha

  # Module 4: Integration
  mpcs_res <- calcular_I_MPCS(I_Graph, I_Markov, I_Games)

  results_list[[reg_name]] <- list(
    Region = reg_name,
    Code = reg_code,
    N = n_reg,
    Alpha = round(alpha_reg, 3),
    I_Graph = round(I_Graph, 4),
    I_Markov = round(I_Markov, 4),
    I_Games = round(I_Games, 4),
    I_MPCS = round(mpcs_res$I_MPCS, 4),
    k = round(mpcs_res$k, 4),
    Optimal_node = opt_node,
    Nudge_type = mpcs_res$tipo,
    v0 = v0_reg,
    grafo = grafo_res,
    markov = markov_res,
    games = games_res
  )
}

# Assemble results table
TABLE6_COMP <- do.call(rbind, lapply(results_list, function(r) {
  data.frame(
    Region = r$Region,
    N = r$N,
    Alpha = r$Alpha,
    I_Graph = r$I_Graph,
    I_Markov = r$I_Markov,
    I_Games = r$I_Games,
    I_MPCS = r$I_MPCS,
    k_rec = r$k,
    Nudge_type = r$Nudge_type,
    Optimal_node = r$Optimal_node,
    stringsAsFactors = FALSE
  )
}))

# Sort by I_MPCS descending
TABLE6_COMP <- TABLE6_COMP %>% arrange(desc(I_MPCS))

cat("\n--- TABLE 6 (Computed) ---\n")
print(TABLE6_COMP)

# =============================================================================
# SENSITIVITY ANALYSIS FUNCTIONS
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

  tipo_sim <- case_when(
    k_sim < 0.25 ~ "Informational",
    k_sim < 0.50 ~ "Structural",
    k_sim < 0.75 ~ "Normative",
    TRUE ~ "Systemic multi-nudge"
  )

  w_default <- c(0.35, 0.40, 0.25)
  I_MPCS_default <- w_default[1] * I_Grafo + w_default[2] * I_Markov + w_default[3] * I_Juegos
  k_default <- min(1, I_MPCS_default * R * escala)

  tipo_default <- case_when(
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
# TABLE 7: Sensitivity on Horizon H
# =============================================================================

cat("\n=== GENERATING TABLE 7: Horizon Sensitivity ===\n")

TABLE7_COMP <- do.call(rbind, lapply(names(results_list), function(reg_name) {
  r <- results_list[[reg_name]]
  df_reg <- df_proc %>% filter(HV023 == r$Code)

  I_vals <- sapply(HORIZONS, function(H) {
    markov_h <- calcular_markov_region(df_reg, P_base, H)
    mpcs_h <- calcular_I_MPCS(r$I_Graph, markov_h$I_Markov, r$I_Games)
    round(mpcs_h$I_MPCS, 4)
  })

  data.frame(
    Region = reg_name,
    I_MPCS_H8 = I_vals[1],
    I_MPCS_H10 = I_vals[2],
    I_MPCS_H12 = I_vals[3],
    stringsAsFactors = FALSE
  )
})) %>% arrange(desc(I_MPCS_H10))

print(TABLE7_COMP)

# =============================================================================
# TABLE 8: Sensitivity on Weights
# =============================================================================

cat("\n=== GENERATING TABLE 8: Weight Sensitivity ===\n")

TABLE8_COMP <- do.call(rbind, lapply(names(results_list), function(reg_name) {
  r <- results_list[[reg_name]]
  sens <- analisis_sensibilidad(r$I_Graph, r$I_Markov, r$I_Games, n_sim = 10000)

  data.frame(
    Region = reg_name,
    I_MPCS_default = round(sens$I_MPCS_default, 4),
    Most_frequent_type = sens$tipo_mas_frecuente,
    Pct_match = round(sens$pct_coincidencia, 1),
    CV_pct = round(sens$stats$cv, 2),
    stringsAsFactors = FALSE
  )
}))

print(TABLE8_COMP)

# =============================================================================
# TABLE 9: Sensitivity on Correlation Threshold
# =============================================================================

cat("\n=== GENERATING TABLE 9: Threshold Sensitivity ===\n")

TABLE9_COMP <- do.call(rbind, lapply(names(results_list), function(reg_name) {
  r <- results_list[[reg_name]]
  df_reg <- df_proc %>% filter(HV023 == r$Code)

  thresholds <- c(0.05, 0.10, 0.15)
  opt_nodes <- sapply(thresholds, function(th) {
    g_res <- calcular_grafo_region(df_reg, threshold = th)
    ifelse(is.na(g_res$optimal_node), "NA", g_res$optimal_node)
  })

  data.frame(
    Region = reg_name,
    Threshold_0.05 = opt_nodes[1],
    Threshold_0.10 = opt_nodes[2],
    Threshold_0.15 = opt_nodes[3],
    stringsAsFactors = FALSE
  )
}))

print(TABLE9_COMP)

# =============================================================================
# TABLE 10: Complete Sensitivity (Markov + Game Theory parameters)
# =============================================================================

cat("\n=== GENERATING TABLE 10: Complete Sensitivity ===\n")

recalcular_mpcs_param <- function(fila, df_reg, t1 = 0.48, t2 = 0.89, t3 = 0.62,
                                  a_AA = 2.0, a_AR = NULL, a_RA = NULL, a_RR = NULL) {
  alpha <- mean(df_reg$Acceso_salud, na.rm = TRUE)
  alpha <- max(0, min(1, alpha))
  I_Grafo <- fila$I_Graph

  # Recalculate Markov
  estados <- c("E1_Undiagnosed", "E2_Diagnosed_nonadherent",
               "E3_Partial_adherence", "E4_Full_adherence", "E5_Metabolic_control")

  P <- matrix(c(
    0.750, 0.200, 0.040, 0.010, 0.000,
    0.050, 0.450, t1,    0.020, 0.000,
    0.010, 0.090, 0.005, t2,    0.005,
    0.005, 0.010, 0.080, 0.285, t3,
    0.000, 0.010, 0.020, 0.100, 0.870
  ), nrow = 5, byrow = TRUE, dimnames = list(estados, estados))
  P <- P / rowSums(P)

  markov_res <- calcular_markov_region(df_reg, P, H = 10)
  I_Markov <- markov_res$I_Markov

  # Recalculate Games
  if (is.null(a_AR)) a_AR <- -(0.5 + 0.5 * (1 - alpha))
  if (is.null(a_RA)) a_RA <- 0.5 + 0.5 * alpha
  if (is.null(a_RR)) a_RR <- 0.5 + 0.5 * (1 - alpha)

  p_star <- max(0, min(1, (a_RR - a_AR) / (a_AA - a_AR - a_RA + a_RR)))
  I_Games <- p_star

  # Combine
  w <- c(0.35, 0.40, 0.25)
  I_MPCS <- w[1] * I_Grafo + w[2] * I_Markov + w[3] * I_Games
  list(I_MPCS = I_MPCS, I_Markov = I_Markov, I_Games = I_Games)
}

TABLE10_COMP <- data.frame(
  Parameter = character(),
  Variation = character(),
  Mean_abs_delta = numeric(),
  Ranking_preserved = character(),
  Nudge_type_preserved = character(),
  stringsAsFactors = FALSE
)

# Markov rates
tasas_base <- c(t1 = 0.48, t2 = 0.89, t3 = 0.62)
variaciones <- c("minus20" = -0.20, "minus10" = -0.10, "plus10" = 0.10, "plus20" = 0.20)

for (param_name in names(tasas_base)) {
  for (v_name in names(variaciones)) {
    v_pct <- variaciones[v_name]
    t_nueva <- tasas_base[param_name] * (1 + v_pct)
    t_nueva <- max(0, min(1, t_nueva))

    I_MPCS_vals <- numeric(length(results_list))
    names(I_MPCS_vals) <- names(results_list)

    for (i in seq_along(results_list)) {
      reg_name <- names(results_list)[i]
      r <- results_list[[reg_name]]
      df_reg <- df_proc %>% filter(HV023 == r$Code)

      t1_i <- ifelse(param_name == "t1", t_nueva, tasas_base["t1"])
      t2_i <- ifelse(param_name == "t2", t_nueva, tasas_base["t2"])
      t3_i <- ifelse(param_name == "t3", t_nueva, tasas_base["t3"])

      res <- recalcular_mpcs_param(TABLE6_COMP[i, ], df_reg, t1 = t1_i, t2 = t2_i, t3 = t3_i)
      I_MPCS_vals[i] <- res$I_MPCS
    }

    delta <- abs(I_MPCS_vals - TABLE6_COMP$I_MPCS)
    mean_delta <- round(mean(delta), 3)

    rank_base <- order(TABLE6_COMP$I_MPCS, decreasing = TRUE)
    rank_new <- order(I_MPCS_vals, decreasing = TRUE)
    ranking_ok <- all(rank_base == rank_new)

    k_vals <- pmin(1, I_MPCS_vals * 0.65 * 1.5)
    tipos <- case_when(k_vals < 0.25 ~ "Informational", k_vals < 0.50 ~ "Structural",
                       k_vals < 0.75 ~ "Normative", TRUE ~ "Systemic multi-nudge")
    nudge_ok <- all(tipos == "Normative")

    TABLE10_COMP <- rbind(TABLE10_COMP, data.frame(
      Parameter = param_name,
      Variation = paste0(ifelse(v_pct > 0, "+", ""), v_pct * 100, "%"),
      Mean_abs_delta = mean_delta,
      Ranking_preserved = ifelse(ranking_ok, "Yes", "No"),
      Nudge_type_preserved = ifelse(nudge_ok, "Yes (100%)", "No"),
      stringsAsFactors = FALSE
    ))
  }
}

# Game theory payoffs
for (param_name in c("a_AA", "a_AR", "a_RA", "a_RR")) {
  for (v_pct in c(-0.20, 0.20)) {
    I_MPCS_vals <- numeric(length(results_list))
    names(I_MPCS_vals) <- names(results_list)

    for (i in seq_along(results_list)) {
      reg_name <- names(results_list)[i]
      r <- results_list[[reg_name]]
      df_reg <- df_proc %>% filter(HV023 == r$Code)
      alpha_reg <- mean(df_reg$Acceso_salud, na.rm = TRUE)

      a_AR_base <- -(0.5 + 0.5 * (1 - alpha_reg))
      a_RA_base <- 0.5 + 0.5 * alpha_reg
      a_RR_base <- 0.5 + 0.5 * (1 - alpha_reg)

      args <- list(fila = TABLE6_COMP[i, ], df_reg = df_reg)
      if (param_name == "a_AA") args$a_AA <- 2.0 * (1 + v_pct)
      if (param_name == "a_AR") args$a_AR <- a_AR_base * (1 + v_pct)
      if (param_name == "a_RA") args$a_RA <- a_RA_base * (1 + v_pct)
      if (param_name == "a_RR") args$a_RR <- a_RR_base * (1 + v_pct)

      res <- do.call(recalcular_mpcs_param, args)
      I_MPCS_vals[i] <- res$I_MPCS
    }

    delta <- abs(I_MPCS_vals - TABLE6_COMP$I_MPCS)
    mean_delta <- round(mean(delta), 3)

    rank_base <- order(TABLE6_COMP$I_MPCS, decreasing = TRUE)
    rank_new <- order(I_MPCS_vals, decreasing = TRUE)
    ranking_ok <- all(rank_base == rank_new)

    k_vals <- pmin(1, I_MPCS_vals * 0.65 * 1.5)
    tipos <- case_when(k_vals < 0.25 ~ "Informational", k_vals < 0.50 ~ "Structural",
                       k_vals < 0.75 ~ "Normative", TRUE ~ "Systemic multi-nudge")
    nudge_ok <- all(tipos == "Normative")

    TABLE10_COMP <- rbind(TABLE10_COMP, data.frame(
      Parameter = param_name,
      Variation = paste0(ifelse(v_pct > 0, "+", ""), v_pct * 100, "%"),
      Mean_abs_delta = mean_delta,
      Ranking_preserved = ifelse(ranking_ok, "Yes", "No"),
      Nudge_type_preserved = ifelse(nudge_ok, "Yes (100%)", "No"),
      stringsAsFactors = FALSE
    ))
  }
}

cat("\n--- TABLE 10 (Computed) ---\n")
print(TABLE10_COMP)

# =============================================================================
# PLOT FUNCTIONS
# =============================================================================

plot_markov_tree <- function(P, title = "Markov Transitions") {
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
          from = i, to = j, from_x = i, to_x = j,
          y_pos = y_pos, label_y = label_y,
          prob = round(P[i, j] * 100, 1),
          type = if (j > i) "Advance" else "Regression"
        ))
      }
    }
  }
  if (nrow(transitions) == 0) return(NULL)

  ggplot() +
    geom_segment(data = transitions,
                 aes(x = from_x, xend = to_x, y = y_pos, yend = y_pos,
                     color = type, linewidth = prob),
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed"), alpha = 0.8) +
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
}

plot_evolutionary_game <- function(alpha, p_star, title = "Replicator Dynamics") {
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
    annotate("rect", xmin = 0, xmax = p_star, ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "#E74C3C") +
    annotate("rect", xmin = p_star, xmax = 1, ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "#2ECC71") +
    annotate("text", x = p_star/2, y = max(df$f_A) * 0.95, label = "Unstable", color = "#C0392B", size = 3.5) +
    annotate("text", x = (1+p_star)/2, y = max(df$f_A) * 0.95, label = "Stable", color = "#27AE60", size = 3.5) +
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
             label = paste0("p* = ", round(p_star, 3)), color = "#C0392B", fontface = "bold", size = 4, hjust = 0) +
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
# FIGURE 1: MPCS Ranking
# =============================================================================

cat("\n=== GENERATING FIGURE 1: MPCS Ranking ===\n")

p_ranking <- ggplot(TABLE6_COMP,
                    aes(x = reorder(Region, I_MPCS), y = I_MPCS, fill = Nudge_type)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(aes(label = sprintf("%.4f\n(alpha=%.3f)", I_MPCS, Alpha)),
            hjust = -0.05, size = 3.0, lineheight = 0.9) +
  scale_fill_manual(values = c(
    "Informational" = "#74B3CE", "Structural" = "#2E86AB",
    "Normative" = "#E84855", "Systemic multi-nudge" = "#1A3A5C"
  )) +
  coord_flip() +
  labs(
    title = "MPCS Index by Region — Southern Peru (H=10)",
    subtitle = "Computed from synthetic demo data (ENDES 2024 structure)",
    x = "", y = "MPCS Index", fill = "Recommended nudge type"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "#5D6D7E")) +
  ylim(0, max(TABLE6_COMP$I_MPCS) * 1.25)

ggsave("Figure1_MPCS_ranking_H10_computed.png", p_ranking, width = 11, height = 7, dpi = 150)
cat(">> Figure 1 saved: Figure1_MPCS_ranking_H10_computed.png\n")

# =============================================================================
# FIGURE 2: State Distribution
# =============================================================================

cat("\n=== GENERATING FIGURE 2: State Distribution ===\n")

states_comp <- df_proc %>%
  filter(!is.na(Estado_Markov)) %>%
  group_by(Region = case_when(
    HV023 == 3 ~ "Apurimac", HV023 == 4 ~ "Arequipa", HV023 == 8 ~ "Cusco",
    HV023 == 18 ~ "Moquegua", HV023 == 21 ~ "Puno", HV023 == 23 ~ "Tacna"
  ), Estado_Markov) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Region) %>%
  mutate(Pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(Estado_Short = case_when(
    Estado_Markov == "E1_Undiagnosed" ~ "E1: Undiagnosed",
    Estado_Markov == "E2_Diagnosed_nonadherent" ~ "E2: Diagnosed\nnon-adherent",
    Estado_Markov == "E3_Partial_adherence" ~ "E3: Partial\nadherence",
    Estado_Markov == "E4_Full_adherence" ~ "E4: Full\nadherence",
    Estado_Markov == "E5_Metabolic_control" ~ "E5: Metabolic\ncontrol"
  ))

p_states <- ggplot(states_comp, aes(x = Region, y = Pct, fill = Estado_Short)) +
  geom_col(position = "stack", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "E1: Undiagnosed" = "#E74C3C", "E2: Diagnosed\nnon-adherent" = "#E67E22",
    "E3: Partial\nadherence" = "#F1C40F", "E4: Full\nadherence" = "#2ECC71",
    "E5: Metabolic\ncontrol" = "#1A8754"
  )) +
  labs(title = "Behavioral State Distribution — Southern Peru",
       subtitle = "Computed from synthetic demo data · MPCS Model",
       x = "Region", y = "Percentage (%)", fill = "State") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "#5D6D7E"))

ggsave("Figure2_MPCS_states_H10_computed.png", p_states, width = 12, height = 7, dpi = 150)
cat(">> Figure 2 saved: Figure2_MPCS_states_H10_computed.png\n")

# =============================================================================
# FIGURE 3: Markov Trajectories
# =============================================================================

cat("\n=== GENERATING FIGURE 3: Markov Trajectories ===\n")

traj_list <- list()

for (reg_name in names(results_list)) {
  r <- results_list[[reg_name]]
  k_region <- r$k

  sim_base <- simular_cadena_markov(P_base, r$v0, n_periodos = 30)
  P_nudge <- aplicar_nudge(P_base, k_region)
  sim_nudge <- simular_cadena_markov(P_nudge, r$v0, n_periodos = 30)

  df_base <- as.data.frame(sim_base)
  colnames(df_base) <- colnames(sim_base)
  df_base$Periodo <- 0:(nrow(df_base) - 1)
  df_base$Region <- reg_name
  df_base$Condition <- "Without nudge"
  df_base$Adherence <- df_base$E4_Full_adherence + df_base$E5_Metabolic_control

  df_nudge <- as.data.frame(sim_nudge)
  colnames(df_nudge) <- colnames(sim_nudge)
  df_nudge$Periodo <- 0:(nrow(df_nudge) - 1)
  df_nudge$Region <- reg_name
  df_nudge$Condition <- "With nudge (dynamic k)"
  df_nudge$Adherence <- df_nudge$E4_Full_adherence + df_nudge$E5_Metabolic_control

  traj_list[[reg_name]] <- bind_rows(df_base, df_nudge)
}

traj_comp <- bind_rows(traj_list)

region_colors <- c("Apurimac" = "#E74C3C", "Arequipa" = "#E67E22", "Cusco" = "#2E86AB",
                   "Moquegua" = "#1A8754", "Puno" = "#8E44AD", "Tacna" = "#F39C12")

p_traj <- ggplot(traj_comp, aes(x = Periodo, y = Adherence, color = Region, linetype = Condition)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "gray50", alpha = 0.5) +
  annotate("text", x = 10.5, y = 0.95, label = "H=10", size = 3.5, color = "gray50") +
  facet_wrap(~Condition, ncol = 2) +
  scale_color_manual(values = region_colors) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(title = "Adherence Trajectories — Southern Peru (H=10)",
       subtitle = "MPCS Markov Chains · Computed from demo data · k calibrated per region",
       x = "Period", y = "P(Full adherence + Metabolic control)", color = "Region", linetype = "Condition") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.box = "vertical",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "#5D6D7E"),
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank())

ggsave("Figure3_MPCS_trajectories_H10_computed.png", p_traj, width = 13, height = 7, dpi = 150)
cat(">> Figure 3 saved: Figure3_MPCS_trajectories_H10_computed.png\n")

# =============================================================================
# FIGURE 4: Sensitivity Analysis (Tacna)
# =============================================================================

cat("\n=== GENERATING FIGURE 4: Sensitivity Analysis ===\n")

tacna_res <- results_list[["Tacna"]]
sens_tacna <- analisis_sensibilidad(tacna_res$I_Graph, tacna_res$I_Markov, tacna_res$I_Games, n_sim = 10000)

df_sens <- data.frame(I_MPCS = sens_tacna$I_MPCS_sim, Type = sens_tacna$tipo_sim)

colors_sens <- c("Informational" = "#74B3CE", "Structural" = "#2E86AB",
                 "Normative" = "#E84855", "Systemic multi-nudge" = "#1A3A5C")
types_present <- unique(df_sens$Type)
colors_filtered <- colors_sens[names(colors_sens) %in% types_present]

p_sens <- ggplot(df_sens, aes(x = I_MPCS, fill = Type)) +
  geom_histogram(bins = 50, color = "white", alpha = 0.85) +
  geom_vline(xintercept = sens_tacna$I_MPCS_default, color = "#C0392B", linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = sens_tacna$ic_inf, color = "#2C3E50", linetype = "dotted", linewidth = 0.8, alpha = 0.5) +
  geom_vline(xintercept = sens_tacna$ic_sup, color = "#2C3E50", linetype = "dotted", linewidth = 0.8, alpha = 0.5) +
  annotate("text", x = sens_tacna$I_MPCS_default + 0.008,
           y = max(table(cut(df_sens$I_MPCS, breaks = 50))) * 0.9,
           label = paste0("Default weights\nI_MPCS = ", round(sens_tacna$I_MPCS_default, 4)),
           hjust = 0, size = 3.5, color = "#C0392B") +
  scale_fill_manual(values = colors_filtered) +
  labs(title = "Sensitivity Analysis — Tacna",
       subtitle = "10,000 random combinations of weights w1, w2, w3",
       x = "MPCS Index", y = "Frequency", fill = "Nudge type") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave("Figure4_MPCS_sensitivity_computed.png", p_sens, width = 10, height = 6, dpi = 150)
cat(">> Figure 4 saved: Figure4_MPCS_sensitivity_computed.png\n")

# =============================================================================
# BEHAVIORAL GRAPHS (all regions, all thresholds)
# =============================================================================

cat("\n=== GENERATING BEHAVIORAL GRAPHS ===\n")

for (region in south_regions) {
  for (umbral in c(0.05, 0.10, 0.15)) {
    reg_name <- region$name
    df_reg <- df_proc %>% filter(HV023 == region$code)

    g_res <- calcular_grafo_region(df_reg, threshold = umbral)
    if (is.null(g_res$graph)) next

    g <- g_res$graph
    centr <- g_res$centr
    optimal_node <- g_res$optimal_node

    impact_sorted <- centr$Nudge_impact[match(V(g)$name, centr$Variable)]
    node_size <- scales::rescale(impact_sorted, to = c(14, 30))
    node_color <- ifelse(V(g)$name == optimal_node, "#C0392B", "#F0DFC0")
    node_border <- ifelse(V(g)$name == optimal_node, "#8B1E1E", "#B8A67A")

    edges_df <- as.data.frame(as_edgelist(g))
    colnames(edges_df) <- c("from", "to")
    corr_vals <- sapply(1:nrow(edges_df), function(i) {
      cor_mat <- df_reg %>%
        select(Dx_HTA, Dx_DM, Adh_farma, HTA_medida, IMC, Obesidad_abd,
               Acceso_salud, Fuma, Alcohol, Dieta_sana, Depresion_bin, Educ_alta, Tiene_seguro) %>%
        cor(use = "pairwise.complete.obs", method = "spearman")
      from <- edges_df$from[i]
      to <- edges_df$to[i]
      cor_mat[from, to]
    })

    edge_color <- ifelse(corr_vals > 0, "#5B9BD5", "#D65B5B")
    edge_width <- scales::rescale(abs(corr_vals), to = c(0.6, 3.2))

    set.seed(42)
    lay <- layout_with_fr(g, niter = 5000)
    lay <- norm_coords(lay, ymin = -1, ymax = 1, xmin = -1, xmax = 1) * 2.4

    threshold_str <- gsub("\\.", "_", sprintf("%.2f", umbral))
    filename <- paste0("behavioral_graph_", tolower(gsub(" ", "_", reg_name)),
                       "_thr", threshold_str, "_computed.png")

    png(filename, width = 1400, height = 1100, res = 150)
    par(mar = c(4, 1, 4, 1))
    plot(g, layout = lay, rescale = FALSE,
         xlim = range(lay[, 1]) * 1.15, ylim = range(lay[, 2]) * 1.15,
         vertex.size = node_size, vertex.color = node_color,
         vertex.frame.color = node_border, vertex.label = V(g)$name,
         vertex.label.family = "sans", vertex.label.cex = 0.75,
         vertex.label.color = "black", vertex.label.dist = 0,
         edge.color = edge_color, edge.width = edge_width, edge.curved = 0.15, asp = 0)
    title(main = paste0("Behavioral System Graph — MPCS ", reg_name, " (threshold = ", umbral, ")"),
          sub = paste0("Red node = optimal intervention point: ", optimal_node),
          cex.main = 1.1, cex.sub = 1.0)
    mtext("Blue edges = positive correlation | Red = negative | Size = nudge impact",
          side = 1, line = 2.5, cex = 0.8)
    dev.off()
    cat(sprintf("   Saved: %s\n", filename))
  }
}

# =============================================================================
# ADDITIONAL PLOTS
# =============================================================================

cat("\n=== GENERATING ADDITIONAL PLOTS ===\n")

p_tree <- plot_markov_tree(P_base, "Markov Transitions")
if (!is.null(p_tree)) {
  ggsave("markov_tree_computed.png", p_tree, width = 10, height = 6, dpi = 150)
  cat(">> Markov tree saved: markov_tree_computed.png\n")
}

tacna_alpha <- results_list[["Tacna"]]$games$alpha
tacna_games <- results_list[["Tacna"]]$I_Games
p_game <- plot_evolutionary_game(tacna_alpha, tacna_games, "Game Theory - Tacna")
if (!is.null(p_game)) {
  ggsave("evolutionary_game_tacna_computed.png", p_game, width = 12, height = 6, dpi = 150)
  cat(">> Evolutionary game saved: evolutionary_game_tacna_computed.png\n")
}

# =============================================================================
# EXPORT TABLES
# =============================================================================

cat("\n=== EXPORTING TABLES ===\n")

write.csv(TABLE6_COMP, "MPCS_table6_H10_computed.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(TABLE7_COMP, "MPCS_table7_horizon_sensitivity_computed.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(TABLE8_COMP, "MPCS_table8_weight_sensitivity_computed.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(TABLE9_COMP, "MPCS_table9_threshold_sensitivity_computed.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(TABLE10_COMP, "MPCS_table10_complete_sensitivity_computed.csv", row.names = FALSE, fileEncoding = "UTF-8")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════════════════════╗\n")
cat("║     MPCS COMPUTED SUMMARY — SOUTHERN PERU (GENERALIZABLE VERSION)            ║\n")
cat("╠════════════════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  %-12s  %5s  %6s  %8s  %8s  %8s  %8s  %-15s║\n",
            "Region", "N", "alpha", "I_Graph", "I_Markov", "I_Games", "I_MPCS", "Nudge type"))
cat("╠════════════════════════════════════════════════════════════════════════════════╣\n")
for (i in 1:nrow(TABLE6_COMP)) {
  r <- TABLE6_COMP[i,]
  cat(sprintf("║  %-12s  %5d  %.3f  %.4f  %.4f  %.4f  %.4f  %-15s║\n",
              r$Region, r$N, r$Alpha, r$I_Graph, r$I_Markov, r$I_Games, r$I_MPCS,
              substr(r$Nudge_type, 1, 15)))
}
cat("╚════════════════════════════════════════════════════════════════════════════════╝\n")

cat("\n=== GENERATED FILES ===\n")
cat("  TABLES:\n")
cat("    1. MPCS_table6_H10_computed.csv\n")
cat("    2. MPCS_table7_horizon_sensitivity_computed.csv\n")
cat("    3. MPCS_table8_weight_sensitivity_computed.csv\n")
cat("    4. MPCS_table9_threshold_sensitivity_computed.csv\n")
cat("    5. MPCS_table10_complete_sensitivity_computed.csv\n")
cat("  FIGURES:\n")
cat("    6. Figure1_MPCS_ranking_H10_computed.png\n")
cat("    7. Figure2_MPCS_states_H10_computed.png\n")
cat("    8. Figure3_MPCS_trajectories_H10_computed.png\n")
cat("    9. Figure4_MPCS_sensitivity_computed.png\n")
cat("   10. behavioral_graph_*_thr*_computed.png (18 files)\n")
cat("   11. markov_tree_computed.png\n")
cat("   12. evolutionary_game_tacna_computed.png\n")
cat("==========================================\n")

cat("\n✅ SCRIPT COMPLETED — GENERALIZABLE VERSION\n")
cat("   All indices were computed from the synthetic demo_data.csv.\n")
cat("   Compare with Hardcode version to verify consistency (~0.001 diff expected).\n")
