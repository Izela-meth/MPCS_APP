# =============================================================================
# MPCS — COMPARATIVE ANALYSIS SOUTHERN PERU (FINAL VERSION WITH OPTION B)
# =============================================================================
# Regions: Apurímac, Arequipa, Cusco, Moquegua, Puno, Tacna
# Source: ENDES 2024 — INEI Peru
# =============================================================================

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
# CONFIGURATION — CAMBIA LAS RUTAS AQUÍ
# =============================================================================

# --- Tus rutas reales ---
ruta_csalud <- "C:/Users/Blue/Documents/Rstudio/CENDES/CSALUD01_2024.dta"
ruta_rech0  <- "C:/Users/Blue/Documents/Rstudio/CENDES/RECH0_2024.dta"

# --- Verificar archivos ---
if (!file.exists(ruta_csalud)) {
  cat("ERROR: CSALUD01 no encontrado. Selecciona manualmente...\n")
  ruta_csalud <- file.choose()
}
if (!file.exists(ruta_rech0)) {
  cat("ERROR: RECH0 no encontrado. Selecciona manualmente...\n")
  ruta_rech0 <- file.choose()
}

# --- Regiones ---
regiones_sur <- list(
  list(codigo = 3,  nombre = "Apurimac"),
  list(codigo = 4,  nombre = "Arequipa"),
  list(codigo = 8,  nombre = "Cusco"),
  list(codigo = 18, nombre = "Moquegua"),
  list(codigo = 21, nombre = "Puno"),
  list(codigo = 23, nombre = "Tacna")
)

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
# =============================================================================
procesar_region <- function(df_raw, codigo_region, nombre_region, horizonte = 10, umbral = 0.10) {
  
  cat(sprintf(">> Processing: %s (HV023 = %d)...\n", nombre_region, codigo_region))
  
  df_reg <- df_raw %>% filter(HV023 == codigo_region)
  n_total <- nrow(df_reg)
  cat(sprintf("   n = %d individuals\n", n_total))
  
  # --- Variables ---
  df <- df_reg %>%
    select(HHID, QHCLUSTER, QSSEXO, QS23, QS25N, QS26,
           QS100, QS102, QS104, QS106, QS107, QS109, QS111, QS113,
           QS200, QS208, QS213U, QS219U,
           QS700A, QS700B, QS700D,
           QS900, QS901, QS903S, QS903D, QS905S, QS905D, QS907,
           PESO15_AMAS) %>%
    mutate(across(where(is.numeric), ~ifelse(. %in% c(8,9,98,99,998,999,9998), NA, .)))
  
  # --- Feature engineering ---
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
  
  # --- Markov states ---
  df <- df %>%
    mutate(Estado_Markov = case_when(
      Dx_cualquiera==1 & Adh_farma==1 & !is.na(HTA_medida) & HTA_medida==0 ~ "E5_Control",
      Dx_cualquiera==1 & Adh_farma==1 ~ "E4_Adherencia_plena",
      Dx_cualquiera==1 & Compra_med==1 & (is.na(Adh_farma) | Adh_farma<1) ~ "E3_Adherencia_parcial",
      Dx_cualquiera==1 & (is.na(Compra_med) | Compra_med==0) ~ "E2_Sin_adherencia",
      (is.na(Dx_cualquiera) | Dx_cualquiera==0) ~ "E1_Sin_diagnostico",
      TRUE ~ NA_character_
    ))
  
  estados_orden <- c("E1_Sin_diagnostico","E2_Sin_adherencia","E3_Adherencia_parcial","E4_Adherencia_plena","E5_Control")
  
  tabla_est <- df %>%
    filter(!is.na(Estado_Markov)) %>%
    count(Estado_Markov) %>%
    mutate(Pct = round(n/sum(n)*100, 1)) %>%
    arrange(Estado_Markov)
  
  # --- Alpha (context) ---
  alpha <- mean(df$Acceso_salud, na.rm = TRUE)
  if (is.na(alpha) || is.nan(alpha)) alpha <- 0.60
  
  a_AA <- 2.0; a_AR <- -(0.5 + 0.5 * (1 - alpha))
  a_RA <- 0.5 + 0.5 * alpha; a_RR <- 0.5 + 0.5 * (1 - alpha)
  p_estrella <- max(0, min(1, (a_RR - a_AR) / (a_AA - a_AR - a_RA + a_RR)))
  I_juegos <- p_estrella
  
  # --- Graph ---
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
  
  # --- Graph plot ---
  set.seed(42)
  lay <- layout_with_fr(g, niter = 5000)
  lay <- norm_coords(lay, ymin = -1, ymax = 1, xmin = -1, xmax = 1) * 2.4
  
  nombre_archivo <- paste0("behavioral_graph_", tolower(gsub(" ", "_", nombre_region)), ".png")
  png(nombre_archivo, width = 1400, height = 1100, res = 150)
  par(mar = c(4, 1, 4, 1))
  V(g)$color <- ifelse(V(g)$name == nodo_optimo, "#C0392B", "#F0DFC0")
  V(g)$size <- scales::rescale(centr$Impacto_nudge[match(V(g)$name, centr$Variable)], to = c(14, 30))
  E(g)$color <- ifelse(E(g)$peso > 0, "#5B9BD5", "#D65B5B")
  E(g)$width <- scales::rescale(abs(E(g)$peso), to = c(0.6, 3.2))
  
  plot(g, layout = lay, rescale = FALSE,
       xlim = range(lay[, 1]) * 1.15, ylim = range(lay[, 2]) * 1.15,
       vertex.frame.color = ifelse(V(g)$name == nodo_optimo, "#8B1E1E", "#B8A67A"),
       vertex.label = V(g)$name, vertex.label.cex = 0.75,
       vertex.label.color = "black", vertex.label.dist = 0,
       edge.curved = 0.15, asp = 0)
  title(main = paste0("Behavioral System Graph — MPCS ", nombre_region, " (ENDES 2024)"),
        sub = paste0("Red node = optimal intervention point: ", nodo_optimo), cex.main = 1.1, cex.sub = 1.0)
  mtext("Blue edges = positive correlation | Red = negative | Size = nudge impact",
        side = 1, line = 2.5, cex = 0.8)
  dev.off()
  cat(sprintf("   Graph saved: %s\n", nombre_archivo))
  
  # --- Markov rates ---
  t1 <- 0.48; t2 <- 0.89; t3 <- 0.62
  
  # --- Markov: distribution and matrix ---
  dist_actual <- tabla_est %>% arrange(Estado_Markov) %>% pull(Pct) / 100
  if (length(dist_actual) < 5) {
    dist_completa <- setNames(rep(0,5), estados_orden)
    dist_completa[tabla_est$Estado_Markov] <- dist_actual
    dist_actual <- dist_completa
  }
  dist_actual <- dist_actual / sum(dist_actual)
  
  P_base <- matrix(c(
    0.75, 0.20, 0.04, 0.01, 0.00,          # E1: 0.75+0.20+0.04+0.01+0.00 = 1.00
    0.05, 0.45, 0.48, 0.02, 0.00,          # E2: 0.05+0.45+0.48+0.02+0.00 = 1.00
    0.01, 0.09, 0.005, 0.89, 0.005,        # E3: 0.01+0.09+0.005+0.89+0.005 = 1.00 (t2=0.89 respetado)
    0.005, 0.01, 0.08, 0.285, 0.62,        # E4: 0.005+0.01+0.08+0.285+0.62 = 1.00
    0.00, 0.01, 0.02, 0.10, 0.87           # E5: 0.00+0.01+0.02+0.10+0.87 = 1.00
  ), nrow = 5, byrow = TRUE, dimnames = list(estados_orden, estados_orden))
  
  simular <- function(P, v0, n=25) {
    dist <- matrix(0, n+1, 5, dimnames=list(NULL, estados_orden))
    dist[1,] <- v0
    for (t in 2:(n+1)) dist[t,] <- dist[t-1,] %*% P
    cbind(as.data.frame(dist), Periodo=0:n)
  }
  
  sim_base <- simular(P_base, dist_actual)
  
  # --- I_Markov at horizon H ---
  idx_h <- min(horizonte + 1, nrow(sim_base))
  I_markov <- max(0, min(1, sim_base[idx_h, "E4_Adherencia_plena"] + sim_base[idx_h, "E5_Control"]))
  
  # --- MPCS Index ---
  w <- c(0.35, 0.40, 0.25)
  I_MPCS <- w[1]*indice_grafo + w[2]*I_markov + w[3]*I_juegos
  
  # --- k_rec (calculated BEFORE nudge simulation) ---
  k_rec <- min(1, I_MPCS * 0.65 * 1.5)
  
  # --- Nudge type ---
  tipo_nudge <- case_when(
    k_rec < 0.25 ~ "Informational",
    k_rec < 0.50 ~ "Structural",
    k_rec < 0.75 ~ "Normative",
    TRUE ~ "Systemic multi-nudge"
  )
  
  # --- Nudge simulation (USA k_rec, no 0.4 fijo) ---
  P_nudge <- P_base
  for (i in 1:(nrow(P_nudge)-1)) {
    av <- P_nudge[i, i] * k_rec
    P_nudge[i, i] <- P_nudge[i, i] - av
    P_nudge[i, i+1] <- P_nudge[i, i+1] + av
    P_nudge[i, ] <- P_nudge[i, ] / sum(P_nudge[i, ])
  }
  sim_nudge <- simular(P_nudge, dist_actual)
  
  cat(sprintf("   α = %.3f | Optimal node: %s | I_MPCS = %.4f | Nudge: %s\n",
              alpha, nodo_optimo, I_MPCS, tipo_nudge))
  
  # --- Additional plots ---
  p_arbol <- graficar_arbol_markov(P_base, estados_orden, paste("Markov Transitions -", nombre_region))
  if (!is.null(p_arbol)) {
    ggsave(paste0("markov_tree_", tolower(gsub(" ", "_", nombre_region)), ".png"),
           p_arbol, width = 10, height = 6, dpi = 150)
    cat(sprintf("   Markov tree saved: markov_tree_%s.png\n", tolower(gsub(" ", "_", nombre_region))))
  }
  
  p_juego <- graficar_juego_evolutivo(alpha, p_estrella, paste("Game Theory -", nombre_region))
  if (!is.null(p_juego)) {
    ggsave(paste0("game_theory_", tolower(gsub(" ", "_", nombre_region)), ".png"),
           p_juego, width = 12, height = 6, dpi = 150)
    cat(sprintf("   Game theory plot saved: game_theory_%s.png\n", tolower(gsub(" ", "_", nombre_region))))
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
    tabla_estados = tabla_est
  ))
}  # <--- ESTA LLAVE ES LA QUE FALTABA

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
# STEP 2: PROCESS REGIONS
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
