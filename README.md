# scDEDS: Single-Cell Differential Equation Dynamical System for Gene Regulatory Network Inference

<img src="https://img.shields.io/badge/R%3E%3D-4.4.0-blue?style=flat&logo=R" /> <img src="https://img.shields.io/badge/Platform-Linux%20%7C%20Windows%20(WSL2)-lightgrey" /> <img src="https://img.shields.io/badge/License-MIT-yellow" />

## Overview

The `scDEDS` R package implements a novel framework for inferring Gene Regulatory Networks (GRNs) from paired scRNA-seq and scATAC-seq data (the cells in the two datasets are identical; or, after using other methods to establish a one-to-one correspondence, the name of these matched cells is unified). By modeling gene expression dynamics as a **Differential Equation Dynamical System (DEDS)**, `scDEDS` can predict context-specific regulatory interactions between Transcription Factors (TFs) and their Target Genes (TGs) across different cell states and pseudotemporal trajectories. It is linked to paper *scDEDS: A Discrete Evolutionary Dynamical System for GRN Inference from Paired Single-Cell Multi-omics Data*.

### Key Features

*   **Multi-omics Integration:** Seamlessly integrates scRNA-seq and scATAC-seq data to identify putative TGs based on chromatin accessibility.
*   **Pseudotemporal Ordering:** Utilizes pseudotime analysis to order cells and identify branching points (e.g., cell fate decisions).
*   **Dynamical System Modeling:** Employs a system of differential equations to model the regulatory strength (`theta_p`) of TF-TG pairs.
*   **Machine Learning Optimization:** Combines genetic algorithms (`GA` package) and gradient descent for robust parameter estimation and model training.
*   **Branch-Specific GRNs:** Constructs accurate, branch-specific GRNs, providing insights into dynamic regulatory changes during cellular processes.

## Full Workflow

The complete analytical pipeline of `scDEDS` is structured into four main parts:

### Part 1: Data Preprocessing
*   **`identify_TGs()`:** Annotates chromatin fragments from scATAC-seq and identifies putative Target Genes based on promoter accessibility.
*   **`get_TGs_from_JASPAR2024()`:** Retrieves a list of Transcription Factors and their motifs from the JASPAR2024 database.
*   **`get_expression_and_activity_matrix()`:** Constructs count matrices for TGs/TFs expression and TG activity, creating a Seurat object for downstream analysis.

### Part 2: Pseudotime & Branching Analysis
*   **`get_interest_cell_type_data()`:** Subsets the data to focus on specific cell types of interest (e.g., CD14 Mono, CD4 Naive).
*   **`order_pseudotime_and_divide_branches()`:** Performs pseudotemporal ordering and partitions cells into distinct branches representing trajectories.
*   **`get_genes_pseudotime_info()`:** Maps gene expression and activity onto pseudotime.
*   **`cell_grouping()`:** Groups cells along pseudotime to reduce noise and facilitate the calculation of group-wise averages (`TFE_T`, `TGA_T`, `TGE_T`).

### Part 3: Standard GRN Construction
*   **`get_sGRN_by_TFBS_pwm_by_JASPAR2024()`:** Builds a prior, standard GRN (sGRN) (alse called as initial GRN, iGRN) by scanning TG promoters for TF binding sites (TFBS) using PWM models from JASPAR2024.
*   **`get_branch_sGRN()` & `get_sGRN()`:** Constructs branch-specific and overall cell-type-specific sGRNs, which serve as the ground truth for training the predictive model.

### Part 4: DEDS Predictive Model Training
*   **`spilt_dataset()`:** Splits the TF-TG pair data into training, validation, and test sets.
*   **Model Training:** The core of `scDEDS`. It fits a differential equation model to predict regulatory strength (`theta_p`).
    *   **`set_init_params()`,** **`set_init_params_lower()`,** **`set_init_params_upper()`:** Set initial parameters and constraints for optimization.
    *   The training loop uses a hybrid **Genetic Algorithm (GA)** and **Gradient Descent** approach to minimize the loss function and find optimal parameters.
    *   Model performance is evaluated using Cohen's **Kappa statistic** and other metrics (Precision, Recall, AUC) on the validation and test sets.

## Expected Output

The final output of `scDEDS` is a comprehensive list object (`interest_cell_type_branch_model_train`) containing:
*   **`best_pred_result`:** A dataframe for all TF-TG pairs with:
    *   `theta_s`: The standard regulatory strength from the sGRN.
    *   `theta_p`: The predicted regulatory strength from the DEDS model.
    *   `theta_s_bin`/`theta_p_bin`: Binarized version of the strengths (0/1).
    *   Dataset split information (training/validation/test).
*   **`evaluation`:** A dataframe summarizing model performance metrics (Accuracy, Kappa, F1, AUC, etc.) across the dataset splits.
*   **Model Parameters:** The optimized parameters for the differential equation model for each branch.

## Hardware Recommendations

The model training process is **extremely computationally intensive**.
*   **OS:** Linux is highly recommended. Windows users can use WSL2.
*   **CPU:** A high-core-count CPU (e.g., >= 40 cores) is strongly advised.
*   **RAM:** A big RAM is highly recommended (e.g. for data with 1,000 cells and 20,000 genes, >= 64 GB is recommended, and more is of course better).
*   **Time:** Training for one cell type may take **several days** (e.g., ~11 days for CD14 Mono in the paper on a powerful server with 200 GB RAM in WSL2).

## R Package

*   R package is in https://github.com/hth-buaa/scDEDS.

## Specific Guidelines (Code, Data, and Results Availability)

*   Specific Guidelines (Code, Data, and Results Availability) is in https://github.com/hth-buaa/scDEDS-code-data-and-result/tree/main.
*   See specific Guidelines in *article_experiment_R_code.txt* (it is also the code for experiment in the paper).
*   See the experiment results of the paper in folder *article experiment result* (including the data used in paper). 
*   See the code for figures and some mediate analysis results in folder *article figure and some analysis result file R code*.
