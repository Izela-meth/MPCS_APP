# MPCS Calculator

## Modelo Predictivo de Cambio Conductual por Sistemas

[![Shiny](https://img.shields.io/badge/Shiny-1.7.1-blue)](https://shiny.rstudio.com/)
[![Docker](https://img.shields.io/badge/Docker-Ready-blue)](https://www.docker.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Hugging Face](https://img.shields.io/badge/Hugging%20Face-Spaces-orange)](https://huggingface.co/spaces)

---

## 📋 Descripción

**MPCS Calculator** es una aplicación web interactiva desarrollada en **R Shiny** que implementa el **Modelo Predictivo de Cambio Conductual por Sistemas (MPCS)**. Esta herramienta integra tres frameworks matemáticos:

| Módulo | Herramienta | Propósito |
|--------|-------------|-----------|
| **Módulo 1** | Teoría de Grafos | Identifica el nodo óptimo de intervención |
| **Módulo 2** | Cadenas de Markov | Estima tiempo de convergencia y efecto del nudge |
| **Módulo 3** | Teoría de Juegos Evolutiva | Calcula la masa crítica de adopción |

El MPCS genera un **Índice MPCS** ponderado que se traduce en un **tipo e intensidad de nudge** recomendada para intervenciones conductuales.

---

## 🚀 Probar la aplicación

### En línea (Hugging Face Spaces)

[https://huggingface.co/spaces/tu-usuario/mpcs-calculator](https://huggingface.co/spaces/tu-usuario/mpcs-calculator)

### Localmente con Docker

```bash
# 1. Clonar el repositorio
git clone https://github.com/Izela-meth/MPCS_APP.git
cd MPCS_APP

# 2. Construir la imagen
docker build -t mpcs-calculator .

# 3. Ejecutar el contenedor
docker run -p 7860:7860 mpcs-calculator

# 4. Abrir en el navegador
http://localhost:7860