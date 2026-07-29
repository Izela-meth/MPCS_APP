# =============================================================================
# MPCS — ANÁLISIS COMPARATIVO SUR DEL PERÚ (VERSIÓN FINAL CON OPCIÓN B)
# =============================================================================
# Regiones: Apurímac, Arequipa, Cusco, Moquegua, Puno, Tacna
# Fuente: ENDES 2024 — INEI Perú
# Módulos: CSALUD01_2024.dta + RECH0_2024.dta
#
# MEJORAS INCORPORADAS:
#   1. Markov: I_Markov = adherencia en horizonte fijo (H=10) en lugar de T
#   2. Juegos: matriz de pagos DINÁMICA según α (acceso a salud) por región
#   3. Análisis de sensibilidad: por región (10,000 iteraciones)
#   4. Sensibilidad al horizonte H: H=8, 10, 12
#   5. Gráfico de sensibilidad para la región con mayor I_MPCS
# =============================================================================
# INSTRUCCIONES:
#   1. Verifica las rutas de los archivos en las líneas marcadas con <<<
#   2. Corre el script completo con: source("MPCS_SUR_COMPARATIVO.R")
#   3. Los resultados se guardan en tu carpeta de trabajo
# =============================================================================

# --- PAQUETES ----------------------------------------------------------------
paquetes <- c("haven", "dplyr", "tidyr", "ggplot2", "igraph",
              "markovchain", "reshape2", "scales", "gridExtra")
for (p in paquetes) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(haven); library(dplyr); library(tidyr)
library(ggplot2); library(igraph); library(scales); library(gridExtra)

# =============================================================================
# CONFIGURACIÓN — CAMBIA LAS RUTAS AQUÍ <<<
# =============================================================================
ruta_csalud <- "C:/Users/Blue/Documents/Rstudio/CENDES/CSALUD01_2024.dta"
ruta_rech0  <- "C:/Users/Blue/Documents/Rstudio/CENDES/RECH0_2024.dta"

# Regiones del sur con sus códigos HV023
regiones_sur <- list(
  list(codigo = 3,  nombre = "Apurimac"),
  list(codigo = 4,  nombre = "Arequipa"),
  list(codigo = 8,  nombre = "Cusco"),
  list(codigo = 18, nombre = "Moquegua"),
  list(codigo = 21, nombre = "Puno"),
  list(codigo = 23, nombre = "Tacna")
)

# =============================================================================
# PASO 1: CARGA DE DATOS
# =============================================================================
cat("\n=== CARGANDO DATOS ===\n")
df_raw   <- read_dta(ruta_csalud)
df_hogar <- read_dta(ruta_rech0)
df_hogar_dep <- df_hogar %>% select(HHID, HV023)
df_raw <- df_raw %>% left_join(df_hogar_dep, by = "HHID")
cat(sprintf("Total nacional: %d registros\n\n", nrow(df_raw)))

# =============================================================================
# FUNCIÓN PRINCIPAL — procesa una región y devuelve sus indicadores MPCS
# =============================================================================
procesar_region <- function(df_raw, codigo_region, nombre_region, horizonte = 10, umbral = 0.10) {
  
  cat(sprintf(">> Procesando: %s (HV023 = %d)...\n", nombre_region, codigo_region))
  
  # --- Filtrar región ---
  df_reg <- df_raw %>% filter(HV023 == codigo_region)
  n_total <- nrow(df_reg)
  cat(sprintf("   n = %d personas\n", n_total))
  
  # --- Seleccionar variables ---
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
  
  # --- Ingeniería de variables ---
  df <- df %>%
    mutate(
      IMC          = ifelse(!is.na(QS900) & !is.na(QS901) & QS901 > 0,
                            QS900 / (QS901/100)^2, NA),
      PAS_prom     = rowMeans(cbind(QS903S, QS905S), na.rm=TRUE),
      PAD_prom     = rowMeans(cbind(QS903D, QS905D), na.rm=TRUE),
      HTA_medida   = case_when(PAS_prom >= 140 | PAD_prom >= 90 ~ 1,
                               !is.na(PAS_prom) ~ 0, TRUE ~ NA_real_),
      Obesidad_abd = case_when(QSSEXO==1 & QS907>=94 ~ 1,
                               QSSEXO==2 & QS907>=80 ~ 1,
                               !is.na(QS907) ~ 0, TRUE ~ NA_real_),
      Dx_HTA       = ifelse(QS102==1, 1, ifelse(QS102==2, 0, NA)),
      Dx_DM        = ifelse(QS109==1, 1, ifelse(QS109==2, 0, NA)),
      Dx_cualquiera= ifelse(!is.na(Dx_HTA) | !is.na(Dx_DM),
                            pmax(Dx_HTA, Dx_DM, na.rm=TRUE), NA),
      Adh_HTA      = ifelse(QS106==1, 1, ifelse(QS106==2, 0, NA)),
      Adh_DM       = ifelse(QS113==1, 1, ifelse(QS113==2, 0, NA)),
      Adh_farma    = case_when(
        Dx_HTA==1 & Dx_DM==1 ~ rowMeans(cbind(Adh_HTA, Adh_DM), na.rm=TRUE),
        Dx_HTA==1 ~ Adh_HTA,
        Dx_DM==1  ~ Adh_DM,
        TRUE ~ NA_real_),
      Compra_med   = case_when(QS104==1 | QS111==1 ~ 1,
                               Dx_cualquiera==1 ~ 0, TRUE ~ NA_real_),
      Tiene_seguro = ifelse(QS26==1, 1, ifelse(QS26==2, 0, NA)),
      Control_PA   = ifelse(QS100==1, 1, ifelse(QS100==2, 0, NA)),
      Control_gluc = ifelse(QS107==1, 1, ifelse(QS107==2, 0, NA)),
      Acceso_salud = rowMeans(cbind(Tiene_seguro, Control_PA, Control_gluc), na.rm=TRUE),
      Fuma         = ifelse(QS200==1, 1, ifelse(QS200==2, 0, NA)),
      Alcohol      = ifelse(QS208==1, 1, ifelse(QS208==2, 0, NA)),
      Come_frutas  = ifelse(!is.na(QS213U) & QS213U==1, 1, 0),
      Come_verduras= ifelse(!is.na(QS219U) & QS219U==1, 1, 0),
      Dieta_sana   = rowMeans(cbind(Come_frutas, Come_verduras), na.rm=TRUE),
      Depresion_sx = rowMeans(cbind(QS700A, QS700B, QS700D), na.rm=TRUE),
      Depresion_bin= ifelse(!is.na(Depresion_sx), as.integer(Depresion_sx>=1), NA),
      Educ_alta    = ifelse(!is.na(QS25N), as.integer(QS25N>=3), NA)
    )
  
  # --- Estados de Markov ---
  df <- df %>%
    mutate(Estado_Markov = case_when(
      Dx_cualquiera==1 & Adh_farma==1 &
        !is.na(HTA_medida) & HTA_medida==0        ~ "E5_Control",
      Dx_cualquiera==1 & Adh_farma==1             ~ "E4_Adherencia_plena",
      Dx_cualquiera==1 & Compra_med==1 &
        (is.na(Adh_farma) | Adh_farma<1)          ~ "E3_Adherencia_parcial",
      Dx_cualquiera==1 & (is.na(Compra_med) |
                            Compra_med==0)         ~ "E2_Sin_adherencia",
      (is.na(Dx_cualquiera) | Dx_cualquiera==0)   ~ "E1_Sin_diagnostico",
      TRUE ~ NA_character_
    ))
  
  estados_orden <- c("E1_Sin_diagnostico","E2_Sin_adherencia",
                     "E3_Adherencia_parcial","E4_Adherencia_plena","E5_Control")
  
  tabla_est <- df %>%
    filter(!is.na(Estado_Markov)) %>%
    count(Estado_Markov) %>%
    mutate(Pct = round(n/sum(n)*100, 1)) %>%
    arrange(Estado_Markov)
  
  # ============================================================
  # [REVISADO] MARKOV — Tasas FIJAS de la literatura
  # ============================================================
  t1 <- 0.48  # Diagnosticado → Compra medicación
  t2 <- 0.89  # Compra → Toma correctamente
  t3 <- 0.62  # Adherente → Control metabólico
  
  # ============================================================
  # [REVISADO] TEORÍA DE JUEGOS — Dinámica con α
  # ============================================================
  alpha <- mean(df$Acceso_salud, na.rm = TRUE)
  if (is.na(alpha) || is.nan(alpha)) alpha <- 0.60
  
  a_AA <- 2.0
  a_AR <- -(0.5 + 0.5 * (1 - alpha))
  a_RA <- 0.5 + 0.5 * alpha
  a_RR <- 0.5 + 0.5 * (1 - alpha)
  
  numerador <- a_RR - a_AR
  denominador <- a_AA - a_AR - a_RA + a_RR
  p_estrella <- numerador / denominador
  p_estrella <- max(0, min(1, p_estrella))
  I_juegos <- p_estrella
  
  # --- Grafo conductual ---
  vars_grafo <- c("Dx_HTA","Dx_DM","Adh_farma","HTA_medida","IMC",
                  "Obesidad_abd","Acceso_salud","Fuma","Alcohol",
                  "Dieta_sana","Depresion_bin","Educ_alta","Tiene_seguro")
  
  mat_cor <- df %>%
    select(all_of(vars_grafo)) %>%
    cor(use="pairwise.complete.obs", method="spearman")
  
  aristas_df <- which(abs(mat_cor) > umbral & mat_cor != 1, arr.ind=TRUE) %>%
    as.data.frame() %>%
    mutate(desde = rownames(mat_cor)[row],
           hasta = colnames(mat_cor)[col],
           peso  = mat_cor[cbind(row,col)]) %>%
    filter(row < col) %>%
    select(desde, hasta, peso)
  
  g <- graph_from_data_frame(aristas_df, directed=FALSE)
  grado_max <- max(degree(g))
  
  centr <- data.frame(
    Variable       = V(g)$name,
    Grado          = degree(g),
    Intermediacion = betweenness(g, normalized = FALSE)
  )
  
  n_nodos <- vcount(g)
  denom <- (n_nodos - 1) * (n_nodos - 2) / 2
  
  if (denom > 0) {
    centr$Intermediacion_norm <- centr$Intermediacion / denom
  } else {
    centr$Intermediacion_norm <- 0
  }
  
  centr$Impacto_nudge <- round(0.60 * centr$Intermediacion_norm +
                                 0.40 * (centr$Grado / grado_max), 4)
  
  centr <- centr %>% arrange(desc(Impacto_nudge))
  
  nodo_optimo  <- centr$Variable[1]
  indice_grafo <- centr$Impacto_nudge[1]
  
  # --- Gráfico del grafo ---
  impacto_ordenado <- centr$Impacto_nudge[match(V(g)$name, centr$Variable)]
  tam_nodo <- scales::rescale(impacto_ordenado, to = c(14, 30))
  
  color_nodo <- ifelse(V(g)$name == nodo_optimo, "#C0392B", "#F0DFC0")
  color_borde <- ifelse(V(g)$name == nodo_optimo, "#8B1E1E", "#B8A67A")
  
  peso <- E(g)$peso
  color_arista <- ifelse(peso > 0, "#5B9BD5", "#D65B5B")
  ancho_arista <- scales::rescale(abs(peso), to = c(0.6, 3.2))
  
  set.seed(42)
  lay <- layout_with_fr(g, niter = 5000)
  lay <- norm_coords(lay, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
  lay <- lay * 2.4
  
  nombre_archivo <- paste0("grafo_conductual_", tolower(gsub(" ", "_", nombre_region)), ".png")
  png(nombre_archivo, width = 1400, height = 1100, res = 150)
  par(mar = c(4, 1, 4, 1))
  
  plot(
    g,
    layout = lay,
    rescale = FALSE,
    xlim = range(lay[, 1]) * 1.15,
    ylim = range(lay[, 2]) * 1.15,
    vertex.size = tam_nodo,
    vertex.color = color_nodo,
    vertex.frame.color = color_borde,
    vertex.label = V(g)$name,
    vertex.label.family = "sans",
    vertex.label.cex = 0.75,
    vertex.label.color = "black",
    vertex.label.dist = 0,
    edge.color = color_arista,
    edge.width = ancho_arista,
    edge.curved = 0.15,
    asp = 0
  )
  
  title(main = paste0("Grafo del Sistema Conductual — MPCS ", nombre_region, " (ENDES 2024)"),
        sub  = paste0("Nodo rojo = punto óptimo de intervención: ", nodo_optimo),
        cex.main = 1.1, cex.sub = 1.0)
  
  mtext("Aristas azules = correlación positiva | Rojas = negativa | Tamaño = impacto del nudge",
        side = 1, line = 2.5, cex = 0.8)
  
  dev.off()
  cat(sprintf("   Grafo guardado: %s\n", nombre_archivo))
  
  # --- Markov ---
  dist_actual <- tabla_est %>%
    arrange(Estado_Markov) %>%
    pull(Pct) / 100
  if (length(dist_actual) < 5) {
    estados_obs   <- tabla_est$Estado_Markov
    dist_completa <- setNames(rep(0,5), estados_orden)
    dist_completa[estados_obs] <- dist_actual
    dist_actual <- dist_completa
  }
  dist_actual <- dist_actual / sum(dist_actual)
  
  reg <- 0.05
  
  P_base <- matrix(c(
    0.75, 0.20, 0.04, 0.01, 0.00,
    reg, 1 - t1 - reg - 0.02, t1, 0.015, 0.005,
    0.01, 0.10, 1 - t2 - 0.10, t2, 0.01,
    0.005, 0.01, 0.08, 1 - t3 - 0.095, t3,
    0.00, 0.01, 0.02, 0.10, 0.87
  ), nrow = 5, byrow = TRUE, dimnames = list(estados_orden, estados_orden))
  
  P_base[P_base < 0] <- 0.001
  P_base <- P_base / rowSums(P_base)
  
  simular <- function(P, v0, n=25) {
    dist <- matrix(0, n+1, 5, dimnames=list(NULL, estados_orden))
    dist[1,] <- v0
    for (t in 2:(n+1)) dist[t,] <- dist[t-1,] %*% P
    cbind(as.data.frame(dist), Periodo=0:n)
  }
  
  aplicar_nudge <- function(P, k=0.3) {
    P_n <- P
    for (i in 1:(nrow(P)-1)) {
      av <- P[i,i]*k
      P_n[i,i]   <- P[i,i]   - av
      P_n[i,i+1] <- P[i,i+1] + av
      P_n[i,]    <- P_n[i,] / sum(P_n[i,])
    }
    P_n
  }
  
  sim_base   <- simular(P_base, dist_actual)
  sim_nudge  <- simular(aplicar_nudge(P_base, k=0.40), dist_actual)
  
  # ============================================================
  # [OPCIÓN B] I_Markov = adherencia en horizonte fijo H
  # ============================================================
  idx_h <- horizonte + 1
  
  if (nrow(sim_base) > idx_h) {
    adh_sin_nudge <- sim_base[idx_h, "E4_Adherencia_plena"] + 
      sim_base[idx_h, "E5_Control"]
    adh_con_nudge <- sim_nudge[idx_h, "E4_Adherencia_plena"] + 
      sim_nudge[idx_h, "E5_Control"]
  } else {
    ultimo <- nrow(sim_base)
    adh_sin_nudge <- sim_base[ultimo, "E4_Adherencia_plena"] + 
      sim_base[ultimo, "E5_Control"]
    adh_con_nudge <- sim_nudge[ultimo, "E4_Adherencia_plena"] + 
      sim_nudge[ultimo, "E5_Control"]
  }
  
  I_markov <- adh_sin_nudge
  I_markov <- max(0, min(1, I_markov))
  
  T_base <- horizonte
  T_nudge <- horizonte
  
  # --- Índice MPCS ---
  w <- c(0.35, 0.40, 0.25)
  I_MPCS <- w[1]*indice_grafo + w[2]*I_markov + w[3]*I_juegos
  k_rec  <- min(1, I_MPCS * 0.65 * 1.5)
  
  tipo_nudge <- case_when(
    k_rec < 0.25 ~ "Informativo",
    k_rec < 0.50 ~ "Estructural",
    k_rec < 0.75 ~ "Normativo",
    TRUE         ~ "Sistemico multi-nudge"
  )
  
  # Extraer % de cada estado
  pct_E1 <- tabla_est$Pct[tabla_est$Estado_Markov=="E1_Sin_diagnostico"]
  pct_E2 <- tabla_est$Pct[tabla_est$Estado_Markov=="E2_Sin_adherencia"]
  pct_E4 <- tabla_est$Pct[tabla_est$Estado_Markov=="E4_Adherencia_plena"]
  pct_E5 <- tabla_est$Pct[tabla_est$Estado_Markov=="E5_Control"]
  if (length(pct_E1)==0) pct_E1 <- 0
  if (length(pct_E2)==0) pct_E2 <- 0
  if (length(pct_E4)==0) pct_E4 <- 0
  if (length(pct_E5)==0) pct_E5 <- 0
  
  # --- Verificar cifras de adherencia (para el artículo) ---
  tasa_adh_entre_compradores <- mean(df$Adh_farma[df$Compra_med == 1], na.rm = TRUE)
  if (is.na(tasa_adh_entre_compradores) || is.nan(tasa_adh_entre_compradores)) {
    tasa_adh_entre_compradores <- 0.85
  }
  
  tasa_compra_entre_diagnosticados <- mean(df$Compra_med[df$Dx_cualquiera == 1], na.rm = TRUE)
  if (is.na(tasa_compra_entre_diagnosticados) || is.nan(tasa_compra_entre_diagnosticados)) {
    tasa_compra_entre_diagnosticados <- 0.50
  }
  
  cat(sprintf("   α = %.3f | Nodo óptimo: %s | I_MPCS = %.4f | Nudge: %s\n",
              alpha, nodo_optimo, I_MPCS, tipo_nudge))
  
  # Devolver todos los indicadores
  list(
    Region          = nombre_region,
    N               = n_total,
    Alpha           = round(alpha, 3),
    Pct_E1          = pct_E1,
    Pct_E2          = pct_E2,
    Pct_adherencia  = round(pct_E4 + pct_E5, 1),
    Tasa_compra_pct = round(tasa_compra_entre_diagnosticados * 100, 1),
    Tasa_adh_pct    = round(tasa_adh_entre_compradores * 100, 1),
    Nodo_optimo     = nodo_optimo,
    I_Grafo         = round(indice_grafo, 4),
    I_Markov        = round(I_markov, 4),
    I_Juegos        = round(I_juegos, 4),
    I_MPCS          = round(I_MPCS, 4),
    k_rec           = round(k_rec, 4),
    Tipo_nudge      = tipo_nudge,
    T_sin_nudge     = T_base,
    T_con_nudge     = T_nudge,
    Horizonte       = horizonte,
    sim_base        = sim_base,
    sim_nudge       = sim_nudge,
    dist_actual     = dist_actual,
    tabla_estados   = tabla_est
  )
}

# =============================================================================
# PASO 2: PROCESAR LAS 6 REGIONES (con horizonte H=10)
# =============================================================================
cat("\n=== PROCESANDO 6 REGIONES DEL SUR (H=10) ===\n\n")

resultados <- lapply(regiones_sur, function(r) {
  procesar_region(df_raw, r$codigo, r$nombre, horizonte = 10)
})
names(resultados) <- sapply(regiones_sur, function(r) r$nombre)

# =============================================================================
# PASO 3: TABLA COMPARATIVA PRINCIPAL
# =============================================================================
cat("\n=== TABLA COMPARATIVA — SUR DEL PERÚ (H=10) ===\n\n")

tabla_comp <- do.call(rbind, lapply(resultados, function(r) {
  data.frame(
    Region          = r$Region,
    N               = r$N,
    Alpha           = r$Alpha,
    Pct_E1_sin_dx   = r$Pct_E1,
    Pct_E2_sin_adh  = r$Pct_E2,
    Pct_adherencia  = r$Pct_adherencia,
    Tasa_compra_pct = r$Tasa_compra_pct,
    Tasa_adh_pct    = r$Tasa_adh_pct,
    Nodo_optimo     = r$Nodo_optimo,
    I_Grafo         = r$I_Grafo,
    I_Markov        = r$I_Markov,
    I_Juegos        = r$I_Juegos,
    I_MPCS          = r$I_MPCS,
    k_rec           = r$k_rec,
    Tipo_nudge      = r$Tipo_nudge,
    T_sin_nudge     = r$T_sin_nudge,
    T_con_nudge     = r$T_con_nudge,
    stringsAsFactors = FALSE
  )
}))

print(tabla_comp)

write.csv(tabla_comp, "MPCS_tabla_comparativa_sur_H10.csv",
          row.names=FALSE, fileEncoding="UTF-8")
cat("\n>> Tabla guardada: MPCS_tabla_comparativa_sur_H10.csv\n")

# =============================================================================
# PASO 4: ANÁLISIS DE SENSIBILIDAD AL HORIZONTE H (H=8, 10, 12)
# =============================================================================
cat("\n=== ANÁLISIS DE SENSIBILIDAD AL HORIZONTE H ===\n")

horizontes <- c(8, 10, 12)
ranking_por_h <- list()
resultados_por_h <- list()

for (h in horizontes) {
  cat(sprintf("\n>> Procesando con H = %d...\n", h))
  
  res_h <- lapply(regiones_sur, function(r) {
    procesar_region(df_raw, r$codigo, r$nombre, horizonte = h)
  })
  names(res_h) <- sapply(regiones_sur, function(r) r$nombre)
  resultados_por_h[[as.character(h)]] <- res_h
  
  ranking <- sapply(res_h, function(x) x$I_MPCS)
  ranking_por_h[[as.character(h)]] <- ranking
}

cat("\n=== COMPARACIÓN DE RANKINGS POR HORIZONTE ===\n")
ordenes <- list()
for (h in horizontes) {
  orden <- names(sort(ranking_por_h[[as.character(h)]], decreasing = TRUE))
  ordenes[[as.character(h)]] <- orden
  cat(sprintf("\nH = %d:\n", h))
  cat("  ", paste(orden, collapse = " > "), "\n")
  cat(sprintf("  I_MPCS: %s\n", 
              paste(round(sort(ranking_por_h[[as.character(h)]], decreasing = TRUE), 4), 
                    collapse = ", ")))
}

cat("\n=== ¿SE MANTIENE EL RANKING? ===\n")
if (identical(ordenes[["8"]], ordenes[["10"]]) && 
    identical(ordenes[["10"]], ordenes[["12"]])) {
  cat("✅ El ranking regional es IDÉNTICO para H=8, 10 y 12.\n")
  cat("   El orden es: ", paste(ordenes[["10"]], collapse = " > "), "\n")
} else {
  cat("⚠️ El ranking VARÍA ligeramente entre horizontes.\n")
  cat("   H=8:  ", paste(ordenes[["8"]], collapse = " > "), "\n")
  cat("   H=10: ", paste(ordenes[["10"]], collapse = " > "), "\n")
  cat("   H=12: ", paste(ordenes[["12"]], collapse = " > "), "\n")
}

tabla_H <- do.call(rbind, lapply(horizontes, function(h) {
  data.frame(
    Horizonte = h,
    Region = names(ranking_por_h[[as.character(h)]]),
    I_MPCS = round(ranking_por_h[[as.character(h)]], 4)
  )
}))
write.csv(tabla_H, "MPCS_sensibilidad_horizonte_H.csv", row.names = FALSE)
cat("\n>> Tabla guardada: MPCS_sensibilidad_horizonte_H.csv\n")

# =============================================================================
# PASO 5: ANÁLISIS DE SENSIBILIDAD POR REGIÓN (pesos aleatorios)
# =============================================================================
cat("\n=== ANÁLISIS DE SENSIBILIDAD POR REGIÓN (10,000 iteraciones) ===\n")

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
    k_sim < 0.25 ~ "Informativo",
    k_sim < 0.50 ~ "Estructural",
    k_sim < 0.75 ~ "Normativo",
    TRUE         ~ "Sistémico multi-nudge"
  )
  
  w_default <- c(0.35, 0.40, 0.25)
  I_MPCS_default <- w_default[1] * I_Grafo +
    w_default[2] * I_Markov +
    w_default[3] * I_Juegos
  
  k_default <- min(1, I_MPCS_default * R * escala)
  
  tipo_default <- dplyr::case_when(
    k_default < 0.25 ~ "Informativo",
    k_default < 0.50 ~ "Estructural",
    k_default < 0.75 ~ "Normativo",
    TRUE            ~ "Sistémico multi-nudge"
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

sensibilidad_por_region <- data.frame()

for (nombre_reg in names(resultados)) {
  r <- resultados[[nombre_reg]]
  
  I_Grafo <- r$I_Grafo
  I_Markov <- r$I_Markov
  I_Juegos <- r$I_Juegos
  
  if (is.na(I_Grafo)) I_Grafo <- 0.5
  if (is.na(I_Markov)) I_Markov <- 0.5
  if (is.na(I_Juegos)) I_Juegos <- 0.5
  
  res_sens <- analisis_sensibilidad(
    I_Grafo = I_Grafo,
    I_Markov = I_Markov,
    I_Juegos = I_Juegos,
    n_sim = 10000,
    R = 0.65,
    escala = 1.5,
    seed = 123
  )
  
  sensibilidad_por_region <- rbind(sensibilidad_por_region, data.frame(
    Region = nombre_reg,
    I_MPCS_default = round(res_sens$I_MPCS_default, 4),
    Tipo_default = res_sens$tipo_default,
    Tipo_mas_frecuente = res_sens$tipo_mas_frecuente,
    Pct_coincidencia = round(res_sens$pct_coincidencia, 1),
    CV = round(res_sens$stats$cv, 1)
  ))
}

cat("\n")
print(sensibilidad_por_region)

write.csv(sensibilidad_por_region, "MPCS_sensibilidad_por_region.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n>> Tabla guardada: MPCS_sensibilidad_por_region.csv\n")

# =============================================================================
# PASO 5b: GRÁFICO DE SENSIBILIDAD (para la región con mayor I_MPCS)
# =============================================================================
cat("\n=== GENERANDO GRÁFICO DE SENSIBILIDAD ===\n")

region_top <- tabla_comp$Region[which.max(tabla_comp$I_MPCS)]
cat(sprintf("   Región con mayor I_MPCS: %s\n", region_top))

r_top <- resultados[[region_top]]

I_Grafo_top <- r_top$I_Grafo
I_Markov_top <- r_top$I_Markov
I_Juegos_top <- r_top$I_Juegos

res_sens_top <- analisis_sensibilidad(
  I_Grafo = I_Grafo_top,
  I_Markov = I_Markov_top,
  I_Juegos = I_Juegos_top,
  n_sim = 10000,
  R = 0.65,
  escala = 1.5,
  seed = 123
)

df_sens <- data.frame(
  I_MPCS = res_sens_top$I_MPCS_sim,
  Tipo = res_sens_top$tipo_sim
)

colores_sens <- c(
  "Informativo" = "#74B3CE",
  "Estructural" = "#2E86AB",
  "Normativo" = "#E84855",
  "Sistémico multi-nudge" = "#1A3A5C"
)

tipos_presentes <- unique(df_sens$Tipo)
colores_filtrados <- colores_sens[names(colores_sens) %in% tipos_presentes]

p_sens <- ggplot(df_sens, aes(x = I_MPCS, fill = Tipo)) +
  geom_histogram(bins = 50, color = "white", alpha = 0.85) +
  geom_vline(xintercept = res_sens_top$I_MPCS_default,
             color = "#C0392B", linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = res_sens_top$ic_inf,
             color = "#2C3E50", linetype = "dotted", linewidth = 0.8, alpha = 0.5) +
  geom_vline(xintercept = res_sens_top$ic_sup,
             color = "#2C3E50", linetype = "dotted", linewidth = 0.8, alpha = 0.5) +
  annotate("text", 
           x = res_sens_top$I_MPCS_default + 0.008, 
           y = max(table(cut(df_sens$I_MPCS, breaks = 50))) * 0.9,
           label = paste0("Pesos por defecto\nI_MPCS = ", 
                          round(res_sens_top$I_MPCS_default, 4)),
           hjust = 0, size = 3.5, color = "#C0392B") +
  scale_fill_manual(values = colores_filtrados) +
  labs(
    title = paste("Análisis de Sensibilidad —", region_top),
    subtitle = "10,000 combinaciones aleatorias de pesos w1, w2, w3",
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

ggsave("MPCS_sensibilidad_histograma.png", p_sens,
       width = 10, height = 6, dpi = 150)
cat(">> Gráfico guardado: MPCS_sensibilidad_histograma.png\n")

cat("\n=== RESULTADOS DE SENSIBILIDAD PARA", region_top, "===\n")
cat(sprintf("   Media de I_MPCS: %.4f\n", res_sens_top$stats$media))
cat(sprintf("   Desviación estándar: %.4f\n", res_sens_top$stats$sd))
cat(sprintf("   Coeficiente de variación: %.1f%%\n", res_sens_top$stats$cv))
cat(sprintf("   I_MPCS con pesos por defecto: %.4f\n", res_sens_top$I_MPCS_default))
cat(sprintf("   Tipo de nudge más frecuente: %s (%.1f%%)\n",
            res_sens_top$tipo_mas_frecuente, res_sens_top$pct_mas_frecuente))

# =============================================================================
# PASO 6: GRÁFICOS
# =============================================================================

# --- Gráfico 1: Ranking de I_MPCS ---
p_ranking <- ggplot(tabla_comp,
                    aes(x = reorder(Region, I_MPCS), y = I_MPCS,
                        fill = Tipo_nudge)) +
  geom_col(width=0.7, color="white") +
  geom_text(aes(label=sprintf("%.4f\n(α=%.2f)", I_MPCS, Alpha)),
            hjust=-0.05, size=3.0) +
  scale_fill_manual(values=c(
    "Informativo"         = "#74B3CE",
    "Estructural"         = "#2E86AB",
    "Normativo"           = "#E84855",
    "Sistemico multi-nudge" = "#1A3A5C"
  )) +
  coord_flip() +
  labs(title    = "Índice MPCS por Región — Sur del Perú (H=10)",
       subtitle = "ENDES 2024 · α = indicador de acceso a salud",
       x = "", y = "Índice MPCS",
       fill = "Tipo de nudge recomendado") +
  theme_minimal(base_size=13) +
  theme(legend.position="bottom") +
  ylim(0, max(tabla_comp$I_MPCS) * 1.25)

ggsave("MPCS_ranking_regiones_H10.png", p_ranking,
       width=11, height=7, dpi=150)
cat(">> Gráfico guardado: MPCS_ranking_regiones_H10.png\n")

# --- Gráfico 2: Distribución de estados ---
estados_comp <- do.call(rbind, lapply(names(resultados), function(nombre) {
  r <- resultados[[nombre]]
  r$tabla_estados %>%
    mutate(Region = nombre)
}))

p_estados <- ggplot(estados_comp,
                    aes(x=Region, y=Pct, fill=Estado_Markov)) +
  geom_col(position="stack", color="white", linewidth=0.3) +
  scale_fill_manual(values=c(
    "E1_Sin_diagnostico"    = "#E74C3C",
    "E2_Sin_adherencia"     = "#E67E22",
    "E3_Adherencia_parcial" = "#F1C40F",
    "E4_Adherencia_plena"   = "#2ECC71",
    "E5_Control"            = "#1A8754"
  )) +
  labs(title    = "Distribución de Estados Conductuales — Sur del Perú",
       subtitle = "ENDES 2024 · Modelo MPCS",
       x="Región", y="Porcentaje (%)", fill="Estado") +
  theme_minimal(base_size=13) +
  theme(legend.position="bottom", axis.text.x=element_text(angle=20, hjust=1))

ggsave("MPCS_estados_por_region_H10.png", p_estados,
       width=12, height=7, dpi=150)
cat(">> Gráfico guardado: MPCS_estados_por_region_H10.png\n")

# --- Gráfico 3: Trayectorias Markov ---
tray_comp <- do.call(rbind, lapply(resultados, function(r) {
  bind_rows(
    r$sim_base  %>% mutate(Region=r$Region, Condicion="Sin nudge"),
    r$sim_nudge %>% mutate(Region=r$Region, Condicion="Con nudge (k=0.40)")
  ) %>%
    mutate(Adherencia = E4_Adherencia_plena + E5_Control)
}))

p_tray <- ggplot(tray_comp,
                 aes(x=Periodo, y=Adherencia,
                     color=Region, linetype=Condicion)) +
  geom_line(linewidth=1.1) +
  geom_vline(xintercept = 10, linetype="dashed", color="gray50", alpha=0.5) +
  annotate("text", x=10.5, y=0.95, label="H=10", size=3.5, color="gray50") +
  facet_wrap(~Condicion) +
  scale_color_manual(values=c(
    "#E74C3C","#E67E22","#2E86AB","#1A8754","#8E44AD","#F39C12")) +
  scale_y_continuous(labels=percent_format()) +
  labs(title    = "Trayectorias de Adherencia — Sur del Perú (H=10)",
       subtitle = "Cadenas de Markov MPCS · ENDES 2024",
       x="Período", y="P(Adherencia plena + Control)",
       color="Región", linetype="Condición") +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom")

ggsave("MPCS_trayectorias_comparativas_H10.png", p_tray,
       width=13, height=7, dpi=150)
cat(">> Gráfico guardado: MPCS_trayectorias_comparativas_H10.png\n")

# =============================================================================
# PASO 5c: ANÁLISIS DE SENSIBILIDAD AL UMBRAL (|r| > 0.05, 0.10, 0.15)
# =============================================================================
cat("\n=== ANÁLISIS DE SENSIBILIDAD AL UMBRAL ===\n")

umbrales <- c(0.05, 0.10, 0.15)
nodos_por_umbral <- list()
mpcs_por_umbral <- list()

for (u in umbrales) {
  cat(sprintf("\n>> Probando umbral = %.2f...\n", u))
  
  resultados_umbral <- lapply(regiones_sur, function(r) {
    procesar_region(df_raw, r$codigo, r$nombre, horizonte = 10, umbral = u)
  })
  names(resultados_umbral) <- sapply(regiones_sur, function(r) r$nombre)
  
  nodos <- sapply(resultados_umbral, function(x) x$Nodo_optimo)
  mpcs <- sapply(resultados_umbral, function(x) x$I_MPCS)
  
  nodos_por_umbral[[as.character(u)]] <- nodos
  mpcs_por_umbral[[as.character(u)]] <- mpcs
  
  cat(sprintf("   Nodos óptimos: %s\n", paste(nodos, collapse = ", ")))
  cat(sprintf("   I_MPCS (Tacna): %.4f\n", mpcs["Tacna"]))
}

# Comparar nodos
cat("\n=== COMPARACIÓN DE NODOS POR UMBRAL ===\n")
for (u in umbrales) {
  cat(sprintf("  Umbral %.2f: %s\n", u, paste(nodos_por_umbral[[as.character(u)]], collapse = ", ")))
}

# Verificar si el nodo se mantiene
if (identical(nodos_por_umbral[["0.05"]], nodos_por_umbral[["0.10"]]) &&
    identical(nodos_por_umbral[["0.10"]], nodos_por_umbral[["0.15"]])) {
  cat("\n✅ Los nodos óptimos son IDÉNTICOS para umbrales 0.05, 0.10 y 0.15.\n")
  cat("   La selección del nodo es ROBUSTA al umbral de correlación.\n")
} else {
  cat("\n⚠️ Los nodos óptimos VARÍAN entre umbrales.\n")
  cat("   Se recomienda discutir esta sensibilidad en Limitaciones.\n")
}

# Guardar resultados
tabla_umbral <- do.call(rbind, lapply(umbrales, function(u) {
  data.frame(
    Umbral = u,
    Region = names(nodos_por_umbral[[as.character(u)]]),
    Nodo_optimo = nodos_por_umbral[[as.character(u)]],
    I_MPCS = round(mpcs_por_umbral[[as.character(u)]], 4)
  )
}))
write.csv(tabla_umbral, "MPCS_sensibilidad_umbral.csv", row.names = FALSE)
cat("\n>> Tabla guardada: MPCS_sensibilidad_umbral.csv\n")

# =============================================================================
# RESUMEN FINAL EN CONSOLA
# =============================================================================
cat("\n")
cat("╔════════════════════════════════════════════════════════════════════════════════╗\n")
cat("║          RESUMEN COMPARATIVO MPCS — SUR DEL PERÚ (ENDES 2024)                  ║\n")
cat("╠════════════════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  %-12s  %5s  %6s  %8s  %10s  %8s  %-15s║\n",
            "Región", "N", "α", "I_MPCS", "Nodo", "Nudge k", "Tipo nudge"))
cat("╠════════════════════════════════════════════════════════════════════════════════╣\n")
for (i in 1:nrow(tabla_comp)) {
  r <- tabla_comp[i,]
  cat(sprintf("║  %-12s  %5d  %.3f  %.4f    %-10s  %.4f   %-15s║\n",
              r$Region, r$N, r$Alpha, r$I_MPCS,
              substr(r$Nodo_optimo,1,10), r$k_rec,
              substr(r$Tipo_nudge,1,15)))
}
cat("╚════════════════════════════════════════════════════════════════════════════════╝\n")

cat("\n=== ARCHIVOS GENERADOS ===\n")
cat("  1. MPCS_tabla_comparativa_sur_H10.csv\n")
cat("  2. MPCS_ranking_regiones_H10.png\n")
cat("  3. MPCS_estados_por_region_H10.png\n")
cat("  4. MPCS_trayectorias_comparativas_H10.png\n")
cat("  5. MPCS_sensibilidad_horizonte_H.csv\n")
cat("  6. MPCS_sensibilidad_por_region.csv\n")
cat("  7. MPCS_sensibilidad_histograma.png\n")
cat("  8. grafo_conductual_*.png (6 archivos, uno por región)\n")
cat("==========================================\n")

cat("\n✅ SCRIPT COMPLETADO\n")