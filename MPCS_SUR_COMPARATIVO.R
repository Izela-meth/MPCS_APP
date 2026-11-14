# =============================================================================
# MPCS — COMPARATIVE ANALYSIS SOUTHERN PERU (GENERALIZABLE)
# =============================================================================
# Regions: Apurímac, Arequipa, Cusco, Moquegua, Puno, Tacna
# Source: ENDES 2024 — INEI Peru
#
# NOTE: This is the generalizable implementation. The behavioral transition
# probabilities (E2→E3, E3→E4, E4→E5) are estimated from each region's data
# using estimar_tasas_transversal(). The E1→E2 transition is calibrated at
# 0.20 (diagnostic event, not behavioral transition; see manuscript Module 2).
# Outputs differ from the published values (Hardcode script, which uses
# calibrated rates t1=0.48, t2=0.89, t3=0.62) by 0.001–0.012, reflecting
# the difference between per-region estimates and calibrated rates. The
# nudge-type classification (Normative) is identical in both versions.

# --- PACKAGES ---
paquetes <- c("haven", "dplyr", "tidyr", "ggplot2", "igraph",
              "markovchain", "reshape2", "scales", "gridExtra", "patchwork")
for (p in paquetes) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(haven); library(dplyr); library(tidyr)
library(ggplot2); library(igraph); library(scales); library(gridExtra)
library(patchwork)

# =============================================================================
# CONFIGURATION
# =============================================================================
ruta_csalud <- "data/CSALUD01_2024.dta"
ruta_rech0  <- "data/RECH0_2024.dta"

if (!file.exists(ruta_csalud)) {
  cat("File not found. Please select CSALUD01_2024.dta manually...\n")
  ruta_csalud <- file.choose()
}
if (!file.exists(ruta_rech0)) {
  cat("File not found. Please select RECH0_2024.dta manually...\n")
  ruta_rech0 <- file.choose()
}

regiones_sur <- list(
  list(codigo = 3,  nombre = "Apurimac"),
  list(codigo = 4,  nombre = "Arequipa"),
  list(codigo = 8,  nombre = "Cusco"),
  list(codigo = 18, nombre = "Moquegua"),
  list(codigo = 21, nombre = "Puno"),
  list(codigo = 23, nombre = "Tacna")
)

# =============================================================================
# FUNCIÓN: estimar_tasas_transversal
# Devuelve m-1 tasas: t1 (E1→E2), t2 (E2→E3), t3 (E3→E4), t4 (E4→E5)
# =============================================================================
estimar_tasas_transversal <- function(estados, orden_estados = NULL,
                                      pesos = NULL) {
  estados <- as.character(estados)
  validos <- !is.na(estados) & estados != ""
  estados <- estados[validos]
  if (!is.null(pesos)) pesos <- pesos[validos]
  if (length(estados) < 30) stop("Se requieren al menos 30 observaciones válidas")
  if (is.null(orden_estados)) stop("Debe proporcionar orden_estados")
  m <- length(orden_estados)
  if (m < 3) stop("Se requieren al menos 3 estados")
  if (is.null(pesos)) pesos <- rep(1, length(estados))
  if (length(pesos) != length(estados)) stop("pesos y estados deben tener igual longitud")
  
  n <- sapply(orden_estados, function(e) sum(pesos[estados == e], na.rm = TRUE))
  names(n) <- orden_estados
  
  tasas <- numeric(m - 1)
  for (j in 1:(m - 1)) {
    denom <- sum(n[j:m])
    num   <- sum(n[(j + 1):m])
    tasas[j] <- if (denom > 0) num / denom else NA_real_
  }
  names(tasas) <- paste0("t", 1:(m - 1))
  
  list(tasas = tasas, conteos = n, orden = orden_estados, n_total = sum(n))
}

# =============================================================================
# FUNCIÓN: construir_matriz_P
# =============================================================================
construir_matriz_P <- function(t1, t2, t3, t4, estados_orden) {
  P <- matrix(0, 5, 5, dimnames = list(estados_orden, estados_orden))
  P[1, ] <- c(0,     t1,   0.04,  0.01,  0.00)
  P[2, ] <- c(0.05,  0,    t2,    0.02,  0.00)
  P[3, ] <- c(0.01,  0.09, 0,     t3,    0.005)
  P[4, ] <- c(0.005, 0.01, 0.08,  0,     t4)
  P[5, ] <- c(0.00,  0.01, 0.02,  0.10,  0.87)
  
  for (i in 1:5) {
    P[i, i] <- 1 - sum(P[i, -i])
    if (P[i, i] < 0) {
      P[i, -i] <- P[i, -i] / sum(P[i, -i]) * 0.99
      P[i, i] <- 0.01
    }
  }
  P
}

# =============================================================================
# FUNCIÓN: graficar_arbol_markov
# =============================================================================
graficar_arbol_markov <- function(P, estados, titulo = "Markov Transitions") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  if (is.null(P) || nrow(P) < 2) return(NULL)
  
  m <- nrow(P)
  nodos <- data.frame(id = 1:m, nombre = estados, x = 1:m, y = 0)
  
  transiciones <- data.frame()
  for (i in 1:m) {
    for (j in 1:m) {
      if (P[i, j] >= 0.05 && abs(j - i) <= 1 && i != j) {
        y_pos <- if (j > i) 0.3 else -0.3
        label_y <- y_pos + if (j > i) 0.12 else -0.12
        transiciones <- rbind(transiciones, data.frame(
          from = i, to = j,
          from_x = i, to_x = j,
          y_pos = y_pos,
          label_y = label_y,
          prob = round(P[i, j] * 100, 1),
          tipo = if (j > i) "Advance" else "Regression"
        ))
      }
    }
  }
  if (nrow(transiciones) == 0) return(NULL)
  
  p <- ggplot() +
    geom_segment(data = transiciones,
                 aes(x = from_x, xend = to_x, y = y_pos, yend = y_pos, color = tipo, linewidth = prob),
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed"), alpha = 0.8) +
    scale_linewidth_continuous(range = c(0.5, 2.5), guide = "none") +
    scale_color_manual(values = c("Advance" = "#2E86AB", "Regression" = "#E84855")) +
    geom_label(data = transiciones,
               aes(x = (from_x + to_x)/2, y = label_y, label = paste0(prob, "%")),
               size = 4, fontface = "bold", fill = "white", color = "#1A3A5C",
               label.size = 0.3, label.padding = unit(0.2, "lines"), alpha = 0.95) +
    geom_point(data = nodos, aes(x = x, y = y), size = 28, color = "#2C3E50",
               fill = "#F8F9F9", shape = 21, stroke = 2.5) +
    geom_text(data = nodos, aes(x = x, y = y, label = nombre), size = 4,
              fontface = "bold", color = "#1A3A5C") +
    xlim(0.5, m + 0.5) + ylim(-0.8, 0.8) +
    labs(title = titulo, subtitle = "Transitions between consecutive states (>= 5%)",
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
# FUNCIÓN: graficar_juego_evolutivo
# =============================================================================
graficar_juego_evolutivo <- function(alpha, p_star, titulo = "Replicator Dynamics") {
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
             label = paste0("p* = ", round(p_star, 3)), color = "#C0392B",
             fontface = "bold", size = 4, hjust = 0) +
    scale_x_continuous(labels = percent_format(), limits = c(0, 1)) +
    scale_y_continuous(labels = number_format()) +
    scale_fill_manual(values = c("#E74C3C", "#2ECC71"), guide = "none") +
    labs(x = "Proportion of adopters (p)", y = "dp/dt") +
    theme_minimal(base_size = 11) + theme(plot.title = element_blank())
  
  p1 + p2 + plot_annotation(title = titulo,
                            subtitle = paste0("α = ", round(alpha, 3), " | p* = ", round(p_star, 3)),
                            theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                                          plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#7F8C8D")))
}

# =============================================================================
# FUNCIÓN PRINCIPAL: procesar_region
# Acepta:
#   - tasas_forzadas = NULL  → estima 4 tasas
#   - tasas_forzadas de long 4 → usa esas 4 tasas
#   - tasas_forzadas de long 3 → usa 3 tasas + E1→E2=0.20 fija (estilo hardcode)
# =============================================================================
procesar_region <- function(df_raw, codigo_region, nombre_region,
                            horizonte = 10, umbral = 0.10,
                            tasas_forzadas = NULL,
                            a_AA = NULL, a_AR = NULL, a_RA = NULL, a_RR = NULL,
                            verbose = TRUE) {
  
  if (verbose) cat(sprintf(">> Processing: %s (HV023 = %d)...\n", nombre_region, codigo_region))
  
  df_reg <- df_raw %>% filter(HV023 == codigo_region)
  n_total <- nrow(df_reg)
  if (verbose) cat(sprintf("   n = %d individuals\n", n_total))
  
  df <- df_reg %>%
    select(HHID, QHCLUSTER, QSSEXO, QS23, QS25N, QS26,
           QS100, QS102, QS104, QS106, QS107, QS109, QS111, QS113,
           QS200, QS208, QS213U, QS219U,
           QS700A, QS700B, QS700D,
           QS900, QS901, QS903S, QS903D, QS905S, QS905D, QS907,
           PESO15_AMAS) %>%
    mutate(across(where(is.numeric), ~ifelse(. %in% c(8,9,98,99,998,999,9998), NA, .)))
  
  df <- df %>%
    mutate(
      IMC = ifelse(!is.na(QS900) & !is.na(QS901) & QS901 > 0, QS900 / (QS901/100)^2, NA),
      PAS_prom = rowMeans(cbind(QS903S, QS905S), na.rm=TRUE),
      PAD_prom = rowMeans(cbind(QS903D, QS905D), na.rm=TRUE),
      HTA_medida = case_when(PAS_prom >= 140 | PAD_prom >= 90 ~ 1, !is.na(PAS_prom) ~ 0, TRUE ~ NA_real_),
      Obesidad_abd = case_when(QSSEXO==1 & QS907>=94 ~ 1, QSSEXO==2 & QS907>=80 ~ 1, !is.na(QS907) ~ 0, TRUE ~ NA_real_),
      Dx_HTA = ifelse(QS102==1, 1, ifelse(QS102==2, 0, NA)),
      Dx_DM = ifelse(QS109==1, 1, ifelse(QS109==2, 0, NA)),
      Dx_cualquiera = ifelse(!is.na(Dx_HTA) | !is.na(Dx_DM), pmax(Dx_HTA, Dx_DM, na.rm=TRUE), NA),
      Adh_HTA = ifelse(QS106==1, 1, ifelse(QS106==2, 0, NA)),
      Adh_DM = ifelse(QS113==1, 1, ifelse(QS113==2, 0, NA)),
      Adh_farma = case_when(Dx_HTA==1 & Dx_DM==1 ~ rowMeans(cbind(Adh_HTA, Adh_DM), na.rm=TRUE),
                            Dx_HTA==1 ~ Adh_HTA, Dx_DM==1 ~ Adh_DM, TRUE ~ NA_real_),
      Compra_med = case_when(QS104==1 | QS111==1 ~ 1, Dx_cualquiera==1 ~ 0, TRUE ~ NA_real_),
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
      Dx_cualquiera==1 & Adh_farma==1 & !is.na(HTA_medida) & HTA_medida==0 ~ "E5_Control",
      Dx_cualquiera==1 & Adh_farma==1 ~ "E4_Adherencia_plena",
      Dx_cualquiera==1 & Compra_med==1 & (is.na(Adh_farma) | Adh_farma<1) ~ "E3_Adherencia_parcial",
      Dx_cualquiera==1 & (is.na(Compra_med) | Compra_med==0) ~ "E2_Sin_adherencia",
      (is.na(Dx_cualquiera) | Dx_cualquiera==0) ~ "E1_Sin_diagnostico",
      TRUE ~ NA_character_
    ))
  
  estados_orden <- c("E1_Sin_diagnostico","E2_Sin_adherencia","E3_Adherencia_parcial",
                     "E4_Adherencia_plena","E5_Control")
  
  tabla_est <- df %>%
    filter(!is.na(Estado_Markov)) %>%
    count(Estado_Markov) %>%
    mutate(Pct = round(n/sum(n)*100, 1)) %>%
    arrange(Estado_Markov)
  
  alpha <- mean(df$Acceso_salud, na.rm = TRUE)
  if (is.na(alpha) || is.nan(alpha)) alpha <- 0.60
  
  if (is.null(a_AA)) a_AA <- 2.0
  if (is.null(a_AR)) a_AR <- -(0.5 + 0.5 * (1 - alpha))
  if (is.null(a_RA)) a_RA <- 0.5 + 0.5 * alpha
  if (is.null(a_RR)) a_RR <- 0.5 + 0.5 * (1 - alpha)
  p_estrella <- max(0, min(1, (a_RR - a_AR) / (a_AA - a_AR - a_RA + a_RR)))
  I_juegos <- p_estrella
  
  vars_grafo <- c("Dx_HTA","Dx_DM","Adh_farma","HTA_medida","IMC",
                  "Obesidad_abd","Acceso_salud","Fuma","Alcohol",
                  "Dieta_sana","Depresion_bin","Educ_alta","Tiene_seguro")
  
  mat_cor <- df %>% select(all_of(vars_grafo)) %>% cor(use="pairwise.complete.obs", method="spearman")
  
  aristas_df <- which(abs(mat_cor) > umbral & mat_cor != 1, arr.ind=TRUE) %>%
    as.data.frame() %>%
    mutate(desde = rownames(mat_cor)[row], hasta = colnames(mat_cor)[col], peso = mat_cor[cbind(row,col)]) %>%
    filter(row < col) %>% select(desde, hasta, peso)
  
  g <- graph_from_data_frame(aristas_df, directed=FALSE)
  grado_max <- max(degree(g))
  
  centr <- data.frame(Variable = V(g)$name, Grado = degree(g), Intermediacion = betweenness(g, normalized = FALSE))
  n_nodos <- vcount(g)
  denom <- (n_nodos - 1) * (n_nodos - 2) / 2
  centr$Intermediacion_norm <- if (denom > 0) centr$Intermediacion / denom else 0
  centr$Impacto_nudge <- round(0.60 * centr$Intermediacion_norm + 0.40 * (centr$Grado / grado_max), 4)
  centr <- centr %>% arrange(desc(Impacto_nudge))
  
  nodo_optimo <- centr$Variable[1]
  indice_grafo <- centr$Impacto_nudge[1]
  
  # --- Tasas: forzadas (3 o 4) o estimadas ---
  if (!is.null(tasas_forzadas) && length(tasas_forzadas) == 3) {
    t1 <- 0.20
    t2 <- tasas_forzadas[1]
    t3 <- tasas_forzadas[2]
    t4 <- tasas_forzadas[3]
    if (verbose) cat(sprintf("   Tasas FORZADAS (3): E1→E2=0.20  E2→E3=%.4f  E3→E4=%.4f  E4→E5=%.4f\n",
                             t2, t3, t4))
  } else if (!is.null(tasas_forzadas) && length(tasas_forzadas) == 4) {
    t1 <- tasas_forzadas[1]
    t2 <- tasas_forzadas[2]
    t3 <- tasas_forzadas[3]
    t4 <- tasas_forzadas[4]
    if (verbose) cat(sprintf("   Tasas FORZADAS (4): E1→E2=%.4f  E2→E3=%.4f  E3→E4=%.4f  E4→E5=%.4f\n",
                             t1, t2, t3, t4))
  } else {
    res_tasas <- estimar_tasas_transversal(
      estados       = df$Estado_Markov,
      orden_estados = estados_orden,
      pesos         = if ("PESO15_AMAS" %in% names(df)) df$PESO15_AMAS else NULL
    )
    # t1 (E1→E2) es un evento de diagnóstico, no una transición de comportamiento.
    # La fórmula de proporciones secuenciales estima la probabilidad de estar
    # diagnosticado en la población total (~0.11-0.15), no la probabilidad
    # condicional a tener la condición. Se calibra a 0.20 (ENDES 2024).
    t1 <- 0.20
    t2 <- res_tasas$tasas[2]   # E2→E3
    t3 <- res_tasas$tasas[3]   # E3→E4
    t4 <- res_tasas$tasas[4]   # E4→E5
    if (verbose) cat(sprintf("   Tasas: t1(E1→E2)=0.2000 (calibrada)  t2=%.4f  t3=%.4f  t4=%.4f  (n=%.0f)\n",
                             t2, t3, t4, res_tasas$n_total))
  }
  
  dist_actual <- tabla_est %>% arrange(Estado_Markov) %>% pull(Pct) / 100
  if (length(dist_actual) < 5) {
    dist_completa <- setNames(rep(0,5), estados_orden)
    dist_completa[tabla_est$Estado_Markov] <- dist_actual
    dist_actual <- dist_completa
  }
  dist_actual <- dist_actual / sum(dist_actual)
  
  P_base <- construir_matriz_P(t1, t2, t3, t4, estados_orden)
  P_base <- P_base / rowSums(P_base)
  
  simular <- function(P, v0, n=25) {
    dist <- matrix(0, n+1, 5, dimnames=list(NULL, estados_orden))
    dist[1,] <- v0
    for (t in 2:(n+1)) dist[t,] <- dist[t-1,] %*% P
    cbind(as.data.frame(dist), Periodo=0:n)
  }
  
  sim_base <- simular(P_base, dist_actual)
  idx_h <- min(horizonte + 1, nrow(sim_base))
  I_markov <- max(0, min(1, sim_base[idx_h, "E4_Adherencia_plena"] + sim_base[idx_h, "E5_Control"]))
  
  w <- c(0.35, 0.40, 0.25)
  I_MPCS <- w[1]*indice_grafo + w[2]*I_markov + w[3]*I_juegos
  # k = intensidad del nudge (Eq. 13: k = min(1, I_MPCS * R * 1.5), R = 0.65)
  k_rec <- min(1, I_MPCS * 0.65 * 1.5)
  
  # La clasificación de la Tabla 4 se aplica sobre I_MPCS, no sobre k.
  tipo_nudge <- case_when(
    I_MPCS < 0.25 ~ "Informational",
    I_MPCS < 0.50 ~ "Structural",
    I_MPCS < 0.75 ~ "Normative",
    TRUE ~ "Systemic multi-nudge"
  )
  
  P_nudge <- P_base
  for (i in 1:(nrow(P_nudge)-1)) {
    av <- P_nudge[i, i] * k_rec
    P_nudge[i, i] <- P_nudge[i, i] - av
    P_nudge[i, i+1] <- P_nudge[i, i+1] + av
    P_nudge[i, ] <- P_nudge[i, ] / sum(P_nudge[i, ])
  }
  sim_nudge <- simular(P_nudge, dist_actual)
  
  if (verbose) {
    cat(sprintf("   α = %.3f | Optimal node: %s | I_MPCS = %.4f | Nudge: %s\n",
                alpha, nodo_optimo, I_MPCS, tipo_nudge))
  }
  
  return(list(
    Region = nombre_region,
    N = n_total,
    Alpha = round(alpha, 3),
    Pct_E1 = ifelse("E1_Sin_diagnostico" %in% tabla_est$Estado_Markov, tabla_est$Pct[tabla_est$Estado_Markov == "E1_Sin_diagnostico"], 0),
    Pct_E2 = ifelse("E2_Sin_adherencia" %in% tabla_est$Estado_Markov, tabla_est$Pct[tabla_est$Estado_Markov == "E2_Sin_adherencia"], 0),
    Nodo_optimo = nodo_optimo,
    I_Grafo = round(indice_grafo, 4),
    I_Markov = round(I_markov, 4),
    I_Juegos = round(I_juegos, 4),
    I_MPCS = round(I_MPCS, 4),
    k_rec = round(k_rec, 4),
    Tipo_nudge = tipo_nudge,
    sim_base = sim_base,
    sim_nudge = sim_nudge,
    dist_actual = dist_actual,
    tabla_estados = tabla_est,
    tasas = c(t1, t2, t3, t4),
    df = df
  ))
}

# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================
cat("\n=== LOADING DATA ===\n")
df_raw <- read_dta(ruta_csalud)
df_hogar <- read_dta(ruta_rech0)
df_hogar_dep <- df_hogar %>% select(HHID, HV023)
df_raw <- df_raw %>% left_join(df_hogar_dep, by = "HHID")
cat(sprintf("Total national: %d records\n\n", nrow(df_raw)))

# =============================================================================
# STEP 2: PROCESS REGIONS (estimación per-región de 4 tasas)
# =============================================================================
cat("\n=== PROCESSING 6 SOUTHERN REGIONS (H=10) ===\n\n")
resultados <- lapply(regiones_sur, function(r) {
  procesar_region(df_raw, r$codigo, r$nombre, horizonte = 10)
})
names(resultados) <- sapply(regiones_sur, function(r) r$nombre)

# =============================================================================
# STEP 3: COMPARATIVE TABLE
# =============================================================================
cat("\n=== COMPARATIVE TABLE — SOUTHERN PERU (H=10) ===\n\n")
tabla_comp <- do.call(rbind, lapply(resultados, function(r) {
  data.frame(
    Region = r$Region,
    N = r$N,
    Alpha = r$Alpha,
    Pct_E1_sin_dx = r$Pct_E1,
    Pct_E2_sin_adh = r$Pct_E2,
    t1_E1E2 = round(r$tasas[1], 4),
    t2_E2E3 = round(r$tasas[2], 4),
    t3_E3E4 = round(r$tasas[3], 4),
    t4_E4E5 = round(r$tasas[4], 4),
    Nodo_optimo = r$Nodo_optimo,
    I_Grafo = r$I_Grafo,
    I_Markov = r$I_Markov,
    I_Juegos = r$I_Juegos,
    I_MPCS = r$I_MPCS,
    k_rec = r$k_rec,
    Tipo_nudge = r$Tipo_nudge,
    stringsAsFactors = FALSE
  )
}))
print(tabla_comp)
write.csv(tabla_comp, "MPCS_comparative_table_south_H10.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("\n>> Table saved: MPCS_comparative_table_south_H10.csv\n")

# =============================================================================
# STEP 4: TRAJECTORIES PLOT
# =============================================================================
cat("\n=== GENERATING COMPARATIVE TRAJECTORIES ===\n")

tray_comp <- do.call(rbind, lapply(resultados, function(r) {
  df_base <- as.data.frame(r$sim_base)
  df_base$Periodo <- 0:(nrow(df_base) - 1)
  df_base$Region <- r$Region
  df_base$Condicion <- "Without nudge"
  df_base$Adherencia <- df_base$E4_Adherencia_plena + df_base$E5_Control
  
  df_nudge <- as.data.frame(r$sim_nudge)
  df_nudge$Periodo <- 0:(nrow(df_nudge) - 1)
  df_nudge$Region <- r$Region
  df_nudge$Condicion <- "With nudge"
  df_nudge$Adherencia <- df_nudge$E4_Adherencia_plena + df_nudge$E5_Control
  
  bind_rows(df_base, df_nudge)
}))

p_tray <- ggplot(tray_comp, aes(x = Periodo, y = Adherencia,
                                color = Region, linetype = Condicion)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "gray50", alpha = 0.5) +
  annotate("text", x = 10.5, y = 0.95, label = "H = 10", size = 3.5, color = "gray50") +
  facet_wrap(~Condicion) +
  scale_color_manual(values = c(
    "#E74C3C", "#E67E22", "#2E86AB", "#1A8754", "#8E44AD", "#F39C12"
  )) +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Adherence Trajectories — Southern Peru (H=10)",
       subtitle = "Markov Chains MPCS · ENDES 2024 · per-region estimated rates",
       x = "Period", y = "P(Full adherence + Control)",
       color = "Region", linetype = "Condition") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("MPCS_trajectories_comparative_H10.png", p_tray,
       width = 13, height = 7, dpi = 150)
cat(">> Plot saved: MPCS_trajectories_comparative_H10.png\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("\n")
cat("╔════════════════════════════════════════════════════════════════════════════════╗\n")
cat("║          COMPARATIVE MPCS SUMMARY — SOUTHERN PERU (ENDES 2024)                 ║\n")
cat("╠════════════════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  %-12s  %5s  %6s  %8s  %10s  %8s  %-20s║\n",
            "Region", "N", "α", "I_MPCS", "Node", "Nudge k", "Nudge type"))
cat("╠════════════════════════════════════════════════════════════════════════════════╣\n")
for (i in 1:nrow(tabla_comp)) {
  r <- tabla_comp[i,]
  cat(sprintf("║  %-12s  %5d  %.3f  %.4f    %-10s  %.4f   %-20s║\n",
              r$Region, r$N, r$Alpha, r$I_MPCS,
              substr(r$Nodo_optimo, 1, 10), r$k_rec,
              substr(r$Tipo_nudge, 1, 20)))
}
cat("╚════════════════════════════════════════════════════════════════════════════════╝\n")

cat("\n✅ SCRIPT COMPLETED\n")

# =============================================================================
# SENSITIVITY ANALYSIS
# =============================================================================
cat("\n=== GENERANDO FIGURAS Y TABLAS DE SENSIBILIDAD ===\n")

# ---------------------------------------------------------------------------
# Wrapper limpio
# ---------------------------------------------------------------------------
calcular_mpcs_numerico <- function(df_raw, codigo_region, nombre_region,
                                   horizonte = 10, umbral = 0.10,
                                   tasas_forzadas = NULL,
                                   a_AA = NULL, a_AR = NULL, a_RA = NULL, a_RR = NULL) {
  res <- procesar_region(df_raw, codigo_region, nombre_region,
                         horizonte = horizonte, umbral = umbral,
                         tasas_forzadas = tasas_forzadas,
                         a_AA = a_AA, a_AR = a_AR, a_RA = a_RA, a_RR = a_RR,
                         verbose = FALSE)
  list(
    Region = res$Region, N = res$N, Alpha = res$Alpha,
    Nodo_optimo = res$Nodo_optimo, I_Grafo = res$I_Grafo,
    I_Markov = res$I_Markov, I_Juegos = res$I_Juegos,
    I_MPCS = res$I_MPCS, k_rec = res$k_rec, Tipo_nudge = res$Tipo_nudge,
    tasas = res$tasas
  )
}

# =============================================================================
# FIGURA 1
# =============================================================================
cat("\n>> Generando Figura 1...\n")

p_fig1 <- ggplot(tabla_comp, aes(x = reorder(Region, I_MPCS), y = I_MPCS, fill = Tipo_nudge)) +
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
    title = "MPCS Index by region, southern Peru (H=10)",
    subtitle = "ENDES 2024 · α = health access indicator · per-region estimated rates",
    x = "", y = "MPCS Index",
    fill = "Recommended nudge type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#5D6D7E")
  ) +
  ylim(0, max(tabla_comp$I_MPCS) * 1.25)

ggsave("Figure1_MPCS_ranking_H10.png", p_fig1, width = 11, height = 7, dpi = 300)
cat("   Figura 1 guardada: Figure1_MPCS_ranking_H10.png\n")

# =============================================================================
# FIGURA 2
# =============================================================================
cat("\n>> Generando Figura 2...\n")

states_comp <- do.call(rbind, lapply(resultados, function(r) {
  r$tabla_estados %>% mutate(Region = r$Region)
}))

states_comp <- states_comp %>%
  mutate(Estado_Short = case_when(
    Estado_Markov == "E1_Sin_diagnostico" ~ "E1: Undiagnosed",
    Estado_Markov == "E2_Sin_adherencia" ~ "E2: Diagnosed\nnon-adherent",
    Estado_Markov == "E3_Adherencia_parcial" ~ "E3: Partial\nadherence",
    Estado_Markov == "E4_Adherencia_plena" ~ "E4: Full\nadherence",
    Estado_Markov == "E5_Control" ~ "E5: Metabolic\ncontrol",
    TRUE ~ Estado_Markov
  ))

p_fig2 <- ggplot(states_comp, aes(x = Region, y = Pct, fill = Estado_Short)) +
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

ggsave("Figure2_MPCS_states_H10.png", p_fig2, width = 12, height = 7, dpi = 300)
cat("   Figura 2 guardada: Figure2_MPCS_states_H10.png\n")

# =============================================================================
# VERIFICACIÓN DE TASAS
# =============================================================================
cat("\n=== VERIFICACIÓN DE TASAS ===\n")

estados_ord <- c("E1_Sin_diagnostico","E2_Sin_adherencia","E3_Adherencia_parcial",
                 "E4_Adherencia_plena","E5_Control")

n_pooled <- sapply(estados_ord, function(e) {
  sum(sapply(resultados, function(r) {
    idx <- which(r$tabla_estados$Estado_Markov == e)
    if (length(idx) == 0) 0 else r$tabla_estados$n[idx]
  }))
})
t1_u <- (n_pooled[3] + n_pooled[4] + n_pooled[5]) / (n_pooled[2] + n_pooled[3] + n_pooled[4] + n_pooled[5])
t2_u <- (n_pooled[4] + n_pooled[5]) / (n_pooled[3] + n_pooled[4] + n_pooled[5])
t3_u <- n_pooled[5] / (n_pooled[4] + n_pooled[5])
cat(sprintf("Test 1 (pooled, sin ponderar):   t1=%.4f  t2=%.4f  t3=%.4f\n", t1_u, t2_u, t3_u))

if (!is.null(resultados[[1]]$df)) {
  df_all <- do.call(rbind, lapply(resultados, function(r) r$df))
  res_w <- estimar_tasas_transversal(df_all$Estado_Markov,
                                     orden_estados = estados_ord,
                                     pesos = df_all$PESO15_AMAS)
  cat(sprintf("Test 2 (pooled, ponderado):      t1=%.4f  t2=%.4f  t3=%.4f  t4=%.4f\n",
              res_w$tasas[1], res_w$tasas[2], res_w$tasas[3], res_w$tasas[4]))
}

if (!is.null(resultados[[1]]$df)) {
  df_all <- do.call(rbind, lapply(resultados, function(r) r$df))
  df_dx  <- df_all %>% filter(Estado_Markov != "E1_Sin_diagnostico")
  res_dx <- estimar_tasas_transversal(df_dx$Estado_Markov,
                                      orden_estados = estados_ord[-1],
                                      pesos = df_dx$PESO15_AMAS)
  cat(sprintf("Test 3 (solo diagnosticados):    t1=%.4f  t2=%.4f  t3=%.4f\n",
              res_dx$tasas[1], res_dx$tasas[2], res_dx$tasas[3]))
}

cat(sprintf("Publicado (3 tasas):             t1=0.4800  t2=0.8900  t3=0.6200\n\n"))

# =============================================================================
# TABLA 7
# =============================================================================
cat("\n>> Calculando Tabla 7 (sensibilidad horizonte)...\n")

tabla_H <- data.frame()
for (H_val in c(8, 10, 12)) {
  for (r in regiones_sur) {
    res <- calcular_mpcs_numerico(df_raw, r$codigo, r$nombre, horizonte = H_val, umbral = 0.10)
    tabla_H <- rbind(tabla_H, data.frame(
      Region = res$Region,
      H = H_val,
      I_MPCS = res$I_MPCS,
      stringsAsFactors = FALSE
    ))
  }
}

tabla7 <- tabla_H %>%
  pivot_wider(names_from = H, values_from = I_MPCS, names_prefix = "H")

cat("\n--- Tabla 7 ---\n")
print(tabla7)
write.csv(tabla7, "MPCS_table7_horizon_sensitivity.csv", row.names = FALSE)

# =============================================================================
# TABLA 8 + FIGURA 4
# =============================================================================
cat("\n>> Calculando Tabla 8 / Figura 4...\n")

set.seed(123)
n_sim <- 10000
tabla8_rows <- list()
sim_por_region <- list()

for (nombre_region in names(resultados)) {
  res_reg <- resultados[[nombre_region]]
  
  I_G <- res_reg$I_Grafo
  I_M <- res_reg$I_Markov
  I_J <- res_reg$I_Juegos
  
  pesos_sim <- matrix(NA, n_sim, 3)
  for (i in 1:n_sim) {
    w <- runif(3)
    w <- w / sum(w)
    pesos_sim[i, ] <- w
  }
  
  I_MPCS_sim <- apply(pesos_sim, 1, function(w) {
    w[1]*I_G + w[2]*I_M + w[3]*I_J
  })
  # Clasificación sobre I_MPCS (Tabla 4), no sobre k.
  tipo_sim <- case_when(
    I_MPCS_sim < 0.25 ~ "Informational",
    I_MPCS_sim < 0.50 ~ "Structural",
    I_MPCS_sim < 0.75 ~ "Normative",
    TRUE ~ "Systemic multi-nudge"
  )
  
  I_MPCS_default <- 0.35*I_G + 0.40*I_M + 0.25*I_J
  tipo_default <- case_when(
    I_MPCS_default < 0.25 ~ "Informational",
    I_MPCS_default < 0.50 ~ "Structural",
    I_MPCS_default < 0.75 ~ "Normative",
    TRUE ~ "Systemic multi-nudge"
  )
  
  freq_tipo <- table(tipo_sim)
  tipo_mas_frec <- names(freq_tipo)[which.max(freq_tipo)]
  pct_match <- sum(tipo_sim == tipo_default) / n_sim * 100
  
  tabla8_rows[[nombre_region]] <- data.frame(
    Region = res_reg$Region,
    I_MPCS_default = round(I_MPCS_default, 4),
    Most_frequent_type = tipo_mas_frec,
    Match_pct = round(pct_match, 1),
    stringsAsFactors = FALSE
  )
  
  sim_por_region[[nombre_region]] <- I_MPCS_sim
}

tabla8 <- do.call(rbind, tabla8_rows)
rownames(tabla8) <- NULL

cat("\n--- Tabla 8 ---\n")
print(tabla8)
write.csv(tabla8, "MPCS_table8_weight_sensitivity.csv", row.names = FALSE)

df_sens <- data.frame(I_MPCS = I_MPCS_sim, Type = tipo_sim)
colores_sens <- c(
  "Informational" = "#74B3CE",
  "Structural" = "#2E86AB",
  "Normative" = "#E84855",
  "Systemic multi-nudge" = "#1A3A5C"
)
tipos_presentes <- unique(df_sens$Type)
colores_filtrados <- colores_sens[names(colores_sens) %in% tipos_presentes]

ic_inf <- quantile(I_MPCS_sim, 0.025)
ic_sup <- quantile(I_MPCS_sim, 0.975)

p_fig4 <- ggplot(df_sens, aes(x = I_MPCS, fill = Type)) +
  geom_histogram(bins = 50, color = "white", alpha = 0.85) +
  geom_vline(xintercept = I_MPCS_default, color = "#C0392B", linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = ic_inf, color = "#2C3E50", linetype = "dotted", linewidth = 0.8, alpha = 0.5) +
  geom_vline(xintercept = ic_sup, color = "#2C3E50", linetype = "dotted", linewidth = 0.8, alpha = 0.5) +
  annotate("text", x = I_MPCS_default + 0.008, 
           y = max(table(cut(df_sens$I_MPCS, breaks = 50))) * 0.9,
           label = paste0("Default weights\nI_MPCS = ", round(I_MPCS_default, 4)),
           hjust = 0, size = 3.5, color = "#C0392B") +
  scale_fill_manual(values = colores_filtrados) +
  labs(
    title = "Sensitivity Analysis",
    subtitle = "10,000 random combinations of weights w1, w2, w3",
    x = "MPCS Index", y = "Frequency",
    fill = "Nudge type"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave("Figure4_MPCS_sensitivity.png", p_fig4, width = 10, height = 6, dpi = 300)
cat("   Figura 4 guardada: Figure4_MPCS_sensitivity.png\n")

# =============================================================================
# TABLA 9
# =============================================================================
cat("\n>> Calculando Tabla 9...\n")

umbrales <- c(0.05, 0.10, 0.15)
tabla9 <- data.frame(Region = sapply(regiones_sur, function(r) r$nombre), stringsAsFactors = FALSE)

for (u in umbrales) {
  nodos <- character(length(regiones_sur))
  for (i in seq_along(regiones_sur)) {
    r <- regiones_sur[[i]]
    res <- calcular_mpcs_numerico(df_raw, r$codigo, r$nombre, horizonte = 10, umbral = u)
    nodos[i] <- ifelse(is.na(res$Nodo_optimo), "No node", as.character(res$Nodo_optimo))
  }
  tabla9[[paste0("Threshold_", gsub("\\.", "_", sprintf("%.2f", u)))]] <- nodos
}

cat("\n--- Tabla 9 ---\n")
print(tabla9)
write.csv(tabla9, "MPCS_table9_threshold_sensitivity.csv", row.names = FALSE)

# =============================================================================
# TABLA 10 (Markov): baseline per-región, 4 tasas
# =============================================================================
cat("\n>> Calculando Tabla 10 (sensibilidad tasas Markov)...\n")

variaciones <- c("minus20" = -0.20, "minus10" = -0.10, "plus10" = 0.10, "plus20" = 0.20)

tasas_base_reg <- do.call(rbind, lapply(resultados, function(r) r$tasas))
rownames(tasas_base_reg) <- names(resultados)
colnames(tasas_base_reg) <- c("t1", "t2", "t3", "t4")

tabla10 <- data.frame(
  Parameter = character(),
  Variation = character(),
  Mean_abs_delta = numeric(),
  Ranking_preserved = character(),
  Nudge_type_preserved = character(),
  stringsAsFactors = FALSE
)

for (param_name in c("t1", "t2", "t3", "t4")) {
  param_idx <- which(colnames(tasas_base_reg) == param_name)
  
  for (v_name in names(variaciones)) {
    v_pct <- variaciones[v_name]
    
    I_MPCS_vals <- numeric(length(regiones_sur))
    for (i in seq_along(regiones_sur)) {
      reg_name <- regiones_sur[[i]]$nombre
      tasas_reg <- tasas_base_reg[reg_name, ]
      tasas_reg[param_idx] <- tasas_reg[param_idx] * (1 + v_pct)
      tasas_reg <- pmin(pmax(tasas_reg, 0), 1)
      
      res <- calcular_mpcs_numerico(df_raw, regiones_sur[[i]]$codigo,
                                    regiones_sur[[i]]$nombre,
                                    horizonte = 10, umbral = 0.10,
                                    tasas_forzadas = as.numeric(tasas_reg))
      I_MPCS_vals[i] <- res$I_MPCS
    }
    
    delta <- abs(I_MPCS_vals - tabla_comp$I_MPCS)
    mean_delta <- round(mean(delta), 3)
    
    rank_base <- order(tabla_comp$I_MPCS, decreasing = TRUE)
    rank_new <- order(I_MPCS_vals, decreasing = TRUE)
    ranking_ok <- all(rank_base == rank_new)
    
    tipos <- case_when(I_MPCS_vals < 0.25 ~ "Informational", I_MPCS_vals < 0.50 ~ "Structural",
                       I_MPCS_vals < 0.75 ~ "Normative", TRUE ~ "Systemic multi-nudge")
    nudge_ok <- all(tipos == tabla_comp$Tipo_nudge)
    
    tabla10 <- rbind(tabla10, data.frame(
      Parameter = param_name,
      Variation = paste0(ifelse(v_pct > 0, "+", ""), v_pct * 100, "%"),
      Mean_abs_delta = mean_delta,
      Ranking_preserved = ifelse(ranking_ok, "Yes", "No"),
      Nudge_type_preserved = ifelse(nudge_ok, "Yes (100%)", "No"),
      stringsAsFactors = FALSE
    ))
  }
}

cat("\n--- Tabla 10 (Markov) ---\n")
print(tabla10)
write.csv(tabla10, "MPCS_table10_markov_sensitivity.csv", row.names = FALSE)
cat("\nNota: t1 (E1→E2) está calibrada a 0.20; t2, t3, t4 son estimadas per-región.\n")
cat("La sensibilidad de t1 mide el efecto de perturbaciones sobre la tasa calibrada.\n") 

# =============================================================================
# TABLA 10 (Game theory)
# =============================================================================
cat("\n>> Calculando Tabla 10 (sensibilidad payoffs game theory)...\n")

tabla10_gt <- data.frame(
  Parameter = character(),
  Variation = character(),
  Mean_abs_delta = numeric(),
  Ranking_preserved = character(),
  Nudge_type_preserved = character(),
  stringsAsFactors = FALSE
)

for (param_name in c("a_AA", "a_AR", "a_RA", "a_RR")) {
  for (v_pct in c(-0.20, 0.20)) {
    
    I_MPCS_vals <- numeric(length(regiones_sur))
    for (i in seq_along(regiones_sur)) {
      r <- regiones_sur[[i]]
      res_base <- calcular_mpcs_numerico(df_raw, r$codigo, r$nombre, 
                                         horizonte = 10, umbral = 0.10)
      alpha_reg <- res_base$Alpha
      
      a_AR_base <- -(0.5 + 0.5 * (1 - alpha_reg))
      a_RA_base <- 0.5 + 0.5 * alpha_reg
      a_RR_base <- 0.5 + 0.5 * (1 - alpha_reg)
      
      args_list <- list(
        df_raw = df_raw,
        codigo_region = r$codigo,
        nombre_region = r$nombre,
        horizonte = 10,
        umbral = 0.10
      )
      
      if (param_name == "a_AA") {
        args_list$a_AA <- 2.0 * (1 + v_pct)
      } else if (param_name == "a_AR") {
        args_list$a_AR <- a_AR_base * (1 + v_pct)
      } else if (param_name == "a_RA") {
        args_list$a_RA <- a_RA_base * (1 + v_pct)
      } else if (param_name == "a_RR") {
        args_list$a_RR <- a_RR_base * (1 + v_pct)
      }
      
      res <- do.call(calcular_mpcs_numerico, args_list)
      I_MPCS_vals[i] <- res$I_MPCS
    }
    
    delta <- abs(I_MPCS_vals - tabla_comp$I_MPCS)
    mean_delta <- round(mean(delta), 3)
    
    rank_base <- order(tabla_comp$I_MPCS, decreasing = TRUE)
    rank_new <- order(I_MPCS_vals, decreasing = TRUE)
    ranking_ok <- all(rank_base == rank_new)
    
    tipos <- case_when(I_MPCS_vals < 0.25 ~ "Informational", I_MPCS_vals < 0.50 ~ "Structural",
                       I_MPCS_vals < 0.75 ~ "Normative", TRUE ~ "Systemic multi-nudge")
    nudge_ok <- all(tipos == tabla_comp$Tipo_nudge)
    
    tabla10_gt <- rbind(tabla10_gt, data.frame(
      Parameter = param_name,
      Variation = paste0(ifelse(v_pct > 0, "+", ""), v_pct * 100, "%"),
      Mean_abs_delta = mean_delta,
      Ranking_preserved = ifelse(ranking_ok, "Yes", "No"),
      Nudge_type_preserved = ifelse(nudge_ok, "Yes (100%)", "No"),
      stringsAsFactors = FALSE
    ))
  }
}

cat("\n--- Tabla 10 (Game Theory) ---\n")
print(tabla10_gt)

tabla10_completa <- rbind(tabla10, tabla10_gt)
write.csv(tabla10_completa, "MPCS_table10_complete_sensitivity.csv", row.names = FALSE)
cat("\n>> Tabla 10 completa guardada: MPCS_table10_complete_sensitivity.csv\n")

# =============================================================================
# ROBUSTEZ: per-región (4 tasas) vs. calibrado (3 tasas + E1→E2=0.20)
# =============================================================================
cat("\n=== ROBUSTEZ: TASAS PER-REGIÓN vs. CALIBRADAS ===\n")

tasas_calibradas <- c(0.48, 0.89, 0.62)

tabla_robustez <- data.frame(
  Region = character(), N = integer(),
  t1_E1E2 = numeric(), t2_E2E3 = numeric(), t3_E3E4 = numeric(), t4_E4E5 = numeric(),
  I_MPCS_per_region = numeric(),
  I_MPCS_calibrado  = numeric(),
  delta = numeric(),
  stringsAsFactors = FALSE
)

for (reg in names(resultados)) {
  r <- resultados[[reg]]
  idx_reg <- which(sapply(regiones_sur, `[[`, "nombre") == reg)
  cod_reg <- regiones_sur[[idx_reg]]$codigo
  
  res_cal <- calcular_mpcs_numerico(df_raw, cod_reg, reg,
                                    horizonte = 10, umbral = 0.10,
                                    tasas_forzadas = tasas_calibradas)
  
  tabla_robustez <- rbind(tabla_robustez, data.frame(
    Region = reg, N = r$N,
    t1_E1E2 = round(r$tasas[1], 4),
    t2_E2E3 = round(r$tasas[2], 4),
    t3_E3E4 = round(r$tasas[3], 4),
    t4_E4E5 = round(r$tasas[4], 4),
    I_MPCS_per_region = round(r$I_MPCS, 4),
    I_MPCS_calibrado  = round(res_cal$I_MPCS, 4),
    delta = round(abs(r$I_MPCS - res_cal$I_MPCS), 4),
    stringsAsFactors = FALSE
  ))
}

print(tabla_robustez)
cat(sprintf("\nDelta máximo: %.4f   Delta medio: %.4f\n",
            max(tabla_robustez$delta), mean(tabla_robustez$delta)))
write.csv(tabla_robustez, "MPCS_supplementary_robustez_tasas.csv", row.names = FALSE)
cat(">> Guardada: MPCS_supplementary_robustez_tasas.csv\n")

cat("\n✅ TODAS LAS FIGURAS Y TABLAS GENERADAS CORRECTAMENTE\n")
