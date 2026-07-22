# Manual de Usuario — MPCS Calculator

## Índice

1. [Introducción](#introducción)
2. [Requisitos del sistema](#requisitos-del-sistema)
3. [Instalación](#instalación)
4. [Guía de uso paso a paso](#guía-de-uso-paso-a-paso)
5. [Interpretación de resultados](#interpretación-de-resultados)
6. [Preguntas frecuentes](#preguntas-frecuentes)

---

## Introducción

El **MPCS Calculator** es una aplicación web interactiva que implementa el **Modelo Predictivo de Cambio Conductual por Sistemas (MPCS)**. Esta herramienta está diseñada para investigadores, profesionales de la salud pública y tomadores de decisiones que necesitan diseñar intervenciones conductuales basadas en evidencia cuantitativa.

El MPCS integra tres herramientas matemáticas:
- **Teoría de Grafos:** Identifica el nodo óptimo de intervención en el sistema conductual.
- **Cadenas de Markov:** Estima el tiempo de convergencia y el efecto de un nudge.
- **Teoría de Juegos Evolutiva:** Calcula la masa crítica de adopción de la conducta.

---

## Requisitos del sistema

### Para usar la versión en línea (Hugging Face)
- Navegador web actualizado (Chrome, Firefox, Edge, Safari)
- Conexión a internet

### Para ejecutar localmente con Docker
- **Docker Desktop** (Windows/Mac) o **Docker Engine** (Linux)
- 2 GB de RAM disponible
- 1 GB de espacio en disco

### Para ejecutar directamente en R
- **R version 4.0 o superior**
- **RStudio** (recomendado)
- Paquetes: shiny, bslib, DT, dplyr, readxl, haven, ggplot2, igraph, tidyr, RColorBrewer, patchwork

---

## Instalación

### Opción A: Usar la versión en línea (recomendada)
1. Ve a: https://huggingface.co/spaces/tu-usuario/mpcs-calculator
2. La aplicación cargará automáticamente
3. No necesitas instalar nada

### Opción B: Instalar con Docker
```bash
# Clonar el repositorio
git clone https://github.com/Izela-meth/MPCS_APP.git
cd MPCS_APP

# Construir la imagen Docker
docker build -t mpcs-calculator .

# Ejecutar el contenedor
docker run -p 7860:7860 mpcs-calculator

# Abrir en el navegador
http://localhost:7860