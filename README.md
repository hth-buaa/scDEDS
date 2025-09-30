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

r
Install from CRAN
install.packages("https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.0.3.tar.gz", repos = NULL, type = "source")
install.packages(c("dplyr", "ggplot2", "Signac", "SeuratObject", "Seurat", "VGAM", "ggrepel", "minpack.lm", "data.table", "tibble", "tidyr", "pROC", "DDRTree", "GA", "psych"))
Install from Bioconductor
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(c("GenomeInfoDb", "TFBSTools", "JASPAR2024", "IRanges", "ChIPseeker", "EnsDb.Hsapiens.v86", "Biobase", "BiocGenerics", "monocle", "Biostrings", "TxDb.Hsapiens.UCSC.hg38.knownGene", "BSgenome.Hsapiens.UCSC.hg38"))
### 2. Install `scDEDS`
Install the latest development version of `scDEDS` from GitHub:
r
devtools::install_github("hth-buaa/scDEDS")
library(scDEDS)
## Quick Start

This is a minimal example to demonstrate the core workflow using built-in PBMC multi-ome data.
r
library(scDEDS)
library(SeuratData) # For example data
Load example PBMC multi-ome data
SeuratData::InstallData("pbmcMultiome")
scRNAseq <- SeuratData::LoadData("pbmcMultiome", "pbmc.rna")
scATACseq <- SeuratData::LoadData("pbmcMultiome", "pbmc.atac")
1. Identify Target Genes (TGs) via chromatin accessibility
results_identify_TGs <- identify_TGs(
genome_for_anno = TxDb.Hsapiens.UCSC.hg38.knownGene,
promoter_range = 50,
scRNAseq = scRNAseq,
scATACseq = scATACseq,
annoDb = "org.Hs.eg.db"
)
2. Get TF information from JASPAR2024 database
TFs_in_JASPAR2024 <- get_TGs_from_JASPAR2024(
scRNAseq = scRNAseq,
species = "Homo sapiens"
)
3. Generate expression and activity matrices
Basic_Info <- get_expression_and_activity_matrix(
TGs = results_identify_TGs$TGs,
TFs = TFs_in_JASPAR2024,
scRNAseq = scRNAseq,
scATACseq = scATACseq,
promoter_range = 50,
path = "/path/to/your/fragments.tsv.gz" # Replace with your path
)
*Note: The full, complete workflow is extensive and computationally intensive. Please refer to the Full Workflow section and the provided test script for detailed instructions.*

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
*   **`get_sGRN_by_TFBS_pwm_by_JASPAR2024()`:** Builds a prior, standard GRN (sGRN) by scanning TG promoters for TF binding sites (TFBS) using PWM models from JASPAR2024.
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
*   **CPU:** A high-core-count CPU (>= 40 cores) is strongly advised.
*   **RAM:** At least 128 GB of RAM is required. 256 GB or more is preferable for larger datasets.
*   **Time:** Training for one cell type can take **several days** (e.g., ~11 days for CD14 Mono on a powerful server).

## Citation

If you use `scDEDS` in your research, please cite:

> *[Publication/Preprint Title]*
> *[Authors]*
> *[Journal/Repository], [Year].*
> *DOI/BioRxiv link to be provided.*

## License

This project is licensed under the MIT License.

## Contact

For questions, bug reports, or contributions, please open an issue on the [GitHub repository](https://github.com/hth-buaa/scDEDS) or contact the maintainer.
