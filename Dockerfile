# ============================================================================
# Dockerfile para MPCS Calculator
# Optimizado para Hugging Face Spaces y Docker Desktop
# ============================================================================
# Autor: [Tu nombre]
# Año: 2025
# Repositorio: https://github.com/Izela-meth/MPCS_APP
# ============================================================================
# INSTRUCCIONES:
#   Construir: docker build -t mpcs-calculator .
#   Ejecutar:  docker run -p 7860:7860 mpcs-calculator
#   Abrir:     http://localhost:7860
# ============================================================================

# --- Usar la imagen oficial de R con Shiny ---
FROM rocker/shiny:4.3.3

# --- Etiquetas para el contenedor ---
LABEL maintainer="[tu-email]" \
      version="1.0.0" \
      description="MPCS Calculator - Modelo Predictivo de Cambio Conductual por Sistemas"

# --- Instalar dependencias del sistema ---
RUN apt-get update && apt-get install -y \
    # Dependencias para paquetes R
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    # Dependencias para gráficos
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    # Dependencias para formatos de datos
    libudunits2-dev \
    libproj-dev \
    # Utilidades
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

# --- Instalar paquetes R necesarios ---
RUN R -e "install.packages(c( \
    'shiny', \
    'bslib', \
    'DT', \
    'dplyr', \
    'tidyr', \
    'readxl', \
    'haven', \
    'ggplot2', \
    'igraph', \
    'RColorBrewer', \
    'patchwork', \
    'rmarkdown', \
    'knitr', \
    'kableExtra' \
), repos='https://cloud.r-project.org/')"

# --- Crear directorio de la aplicación ---
WORKDIR /srv/shiny-server/

# --- Copiar archivos de la aplicación ---
COPY app.R .
COPY functions/ /srv/shiny-server/functions/
COPY data/ /srv/shiny-server/data/
COPY README.md /srv/shiny-server/README.md

# --- Asegurar permisos correctos ---
RUN chown -R shiny:shiny /srv/shiny-server

# --- Exponer el puerto (7860 para Hugging Face) ---
EXPOSE 7860

# --- Configurar variables de entorno ---
ENV PORT=7860 \
    SHINY_PORT=7860 \
    R_CONFIG_ACTIVE=production

# --- Comando para ejecutar la aplicación ---
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/', host='0.0.0.0', port=7860)"]