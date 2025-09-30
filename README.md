# scDEDS: Single-Cell Differential Equation Dynamical System for Gene Regulatory Network Inference

<img src="https://img.shields.io/badge/R%3E%3D-4.4.0-blue?style=flat&logo=R" /> <img src="https://img.shields.io/badge/Platform-Linux%20%7C%20Windows%20(WSL2)-lightgrey" /> <img src="https://img.shields.io/badge/License-MIT-yellow" />

## Overview

The `scDEDS` R package implements a novel framework for inferring Gene Regulatory Networks (GRNs) from single-cell multi-omics data (scRNA-seq and scATAC-seq). By modeling gene expression dynamics as a **Differential Equation Dynamical System (DEDS)**, `scDEDS` can predict context-specific regulatory interactions between Transcription Factors (TFs) and their Target Genes (TGs) across different cell states and pseudotemporal trajectories.

### Key Features

*   **Multi-omics Integration:** Seamlessly integrates scRNA-seq and scATAC-seq data to identify putative TGs based on chromatin accessibility.
*   **Pseudotemporal Ordering:** Utilizes pseudotime analysis to order cells and identify branching points (e.g., cell fate decisions).
*   **Dynamical System Modeling:** Employs a system of differential equations to model the regulatory strength (`theta_p`) of TF-TG pairs.
*   **Machine Learning Optimization:** Combines genetic algorithms (`GA` package) and gradient descent for robust parameter estimation and model training.
*   **Branch-Specific GRNs:** Constructs accurate, branch-specific GRNs, providing insights into dynamic regulatory changes during cellular processes.

## Installation

### 1. Install Dependencies
`scDEDS` requires several CRAN and Bioconductor packages. Please run the following code in your R session to install them.
