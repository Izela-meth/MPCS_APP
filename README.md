# MPCS: Multi-Predictive Change System

[![Shiny App](https://img.shields.io/badge/Shiny-App-blue)](https://mpcs-calculator.onrender.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Reproducible workflow and interactive web application for the **Multi-Predictive 
Change System (MPCS)** model. MPCS integrates graph theory, Markov chains, and 
evolutionary game theory to design behavioral nudges for public health interventions.

## Overview

| Module | Framework | Purpose |
|--------|-----------|---------|
| **Module 1** | Graph Theory | Identifies the optimal intervention node |
| **Module 2** | Markov Chains | Estimates convergence time and nudge effect |
| **Module 3** | Evolutionary Game Theory | Calculates the critical adoption mass |

The MPCS generates a weighted **MPCS Index** that translates into a recommended 
**nudge type and intensity** for behavioral interventions.

## Repository Structure

MPCS/
├── MPCS_SUR_COMPARATIVO.R          # Main analysis script (ENDES 2024 data)
├── MPCS_SUR_COMPARATIVO_Hardcode.R # Reproduces exact manuscript figures/tables
├── functions/
│   └── mpcs_functions.R            # Core functions (graph, Markov, games)
├── app.R                           # Shiny application source
├── data/
│   ├── CSALUD01_2024.dta           # ENDES 2024 health module (user-provided)
│   ├── RECH0_2024.dta              # ENDES 2024 household module (user-provided)
│   └── demo_data.csv               # Synthetic demo dataset (no license needed)
├── outputs/                        # Generated figures and tables
├── README.md
└── LICENSE



## Interactive Web Application

A live version of the MPCS calculator is deployed on **Render**:

🔗 **[https://mpcs-calculator.onrender.com/](https://mpcs-calculator.onrender.com/)**

The Shiny app allows non-technical users to:
- Upload custom datasets
- Run the full MPCS workflow interactively
- Visualize behavioral graphs, Markov trajectories, and game dynamics
- Download results and sensitivity analyses

## Requirements

- **R** ≥ 4.2.0
- R packages: `haven`, `dplyr`, `tidyr`, `ggplot2`, `igraph`, `markovchain`, 
  `scales`, `gridExtra`, `patchwork`, `shiny`

Install dependencies:
```r
install.packages(c("haven", "dplyr", "tidyr", "ggplot2", "igraph", 
                   "markovchain", "scales", "gridExtra", "patchwork", "shiny"))

Usage
1. Full analysis from ENDES microdata
Place CSALUD01_2024.dta and RECH0_2024.dta in the data/ folder, then run: source("MPCS_SUR_COMPARATIVO.R")

This generates:
Figure1_MPCS_ranking_H10.png
Figure2_MPCS_states_H10.png
Figure3_MPCS_trajectories_H10.png
Figure4_MPCS_sensitivity.png
MPCS_table7_horizon_sensitivity.csv
MPCS_table8_weight_sensitivity.csv
MPCS_table9_threshold_sensitivity.csv
MPCS_table10_markov_sensitivity.csv

2. Exact manuscript reproduction
 source("MPCS_SUR_COMPARATIVO_Hardcode.R")

3. Local Shiny app
 shiny::runApp("app.R")

4. Demo without ENDES data
Use the synthetic demo_data.csv to test the workflow without licensed microdata

Numerical Reproducibility Note

The I_MPCS values reported in the manuscript (Table 6) were computed using the
ENDES 2024 microdata release. The provided analysis script (MPCS_SUR_COMPARATIVO.R)
reconstructs the full workflow from raw data and reproduces the same regional
ranking and nudge classifications.
Minor numerical differences (< 0.002 in I_MPCS, < 0.003 in I_Markov) may arise
due to iterative refinements in outlier handling and state-definition logic
during manuscript preparation. These differences do not affect the qualitative
conclusions (regional ranking, optimal nudge node, or nudge type classification)
and are consistent across all six regions.

License

This project is licensed under the MIT License — see the LICENSE file.


Citation


If you use this code or the MPCS model, please cite:
Zela Llanque, I. J. (2026). MPCS: A reproducible workflow integrating graph
theory, Markov chains, and evolutionary game theory for behavioral nudge design.
MethodsX, [Vol], [Pages]. https://doi.org/10.xxxx/xxxxx
