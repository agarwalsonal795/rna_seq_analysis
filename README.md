**Exploring the Transcriptomic Landscape of Mouse Gastrulation : A Single Cell RNA-Seq Analysis Challenge**

## Overview

This repository contains a Jupyter notebook (`rna_eda.ipynb`) for analyzing single-cell RNA sequencing (scRNA-seq) datasets from mouse gastrulation studies, as described in Argelaguet et al. (2019) and Pijuan-Sala et al. (2019). The analysis focuses on preprocessing, exploratory data analysis (EDA), dataset integration using scVI, and differential expression analysis to study gene expression changes across developmental stages (E4.5, E5.5, E6.5, E7.5). The workflow identifies key differentially expressed genes and provides biological insights into the gastrulation process.

The datasets include:
- **Argelaguet Batch 1**: Covers stages E4.5 to E7.5, capturing the exit from pluripotency to lineage specification.
- **Pijuan-Sala Batch 1 and Batch 2**: Focus on early gastrulation (E6.5 and E7.5).

## Dependencies

To run the notebook, ensure you have the following Python packages installed in your environment:

- `scanpy` (>=1.9.0): For single-cell data analysis.
- `scvi-tools` (>=1.0.0): For dataset integration using scVI.
- `pandas` (>=1.5.0): For data manipulation.
- `numpy` (>=1.23.0): For numerical operations.
- `matplotlib` (>=3.7.0): For plotting.
- `seaborn` (>=0.12.0): For enhanced visualizations.
- `scipy` (>=1.10.0): For sparse matrix operations and statistical tests.
- `torch` (>=2.0.0): For GPU support in scVI.
- `lightning` (>=2.0.0): Underlying dependency for scVI.
- `adjustText` (>=0.8): For adjusting text labels in plots.
- `sklearn` (>=1.2.0): For silhouette scores and data scaling.

You can install these dependencies using `pip`:

```bash
pip install scanpy scvi-tools pandas numpy matplotlib seaborn scipy torch lightning adjustText scikit-learn
```

Additionally, ensure you have Jupyter installed to run the notebook:

```bash
pip install jupyter
```

**Hardware Requirements**:
- At least 16 GB of RAM is recommended for handling large scRNA-seq datasets, alternatively can use AWS Cloud.

## Directory Structure

The notebook expects the following directory structure:

```
project_root/
│
├── data/
│   ├── Arg_batch1.parquet
│   ├── Arg_batch1_meta.parquet
│   ├── PJ_batch1.parquet
│   ├── PJ_batch1_meta.parquet
│   ├── PJ_batch2.parquet
│   ├── PJ_batch2_meta.parquet
│   └── mart_export.csv
│
└── rna_eda.ipynb
```

Ensure the `data/` directory contains the required .csv or .parquet files before running the notebook.

## Workflow

The notebook (`rna_eda.ipynb`) is organized into four main parts, each addressing a specific step in the scRNA-seq analysis pipeline:

### Part 1: Data Download and Preprocessing

**Objective**: Load, preprocess, and quality-control the scRNA-seq datasets to prepare them for analysis.

1. **Load Datasets & Create AnnData Objects**:
   - Loads expression data and metadata for Argelaguet Batch 1, Pijuan-Sala Batch 1, and Batch 2.
   - Converts CSV files to Parquet format for efficiency.
   - Creates `AnnData` objects (`adata_arg`, `adata_pj1`, `adata_pj2`) with cells as rows and genes as columns.
   - Filters Argelaguet Batch 1 by `pass_rnaQC` and removes cells with missing `celltype` annotations.

2. **Map Gene Names**:
   - Maps Ensembl IDs to gene symbols using the `mart_export.csv` file for better interpretability.
   - Verifies mapping quality (>96% genes mapped) and handles duplicates.

3. **Quality Control (QC)**:
   - Computes QC metrics: total counts (`n_counts`), number of expressed genes (`n_genes`), and mitochondrial fraction (`mt_frac%`).
   - Visualizes distributions using violin and scatter plots to identify batch differences.
   - Applies dynamic thresholds (5th/99th percentiles) to filter cells based on `n_counts`, `n_genes`, and `mt_frac` (<20%).
   - Filters genes expressed in fewer than 20 cells to reduce noise.

4. **Normalization**:
   - Applies counts per million (CPM) normalization and log-transformation to make expression values comparable across cells.
   - Computes size factors and visualizes their distribution.

5. **Highly Variable Genes (HVG)**:
   - Identifies the top 4,000 highly variable genes (HVGs) using the `cell_ranger` method to reduce dimensionality.
   - Visualizes HVG selection with dispersion plots.

6. **Analyze Differences Between Datasets**:
   - Comments on differences in developmental stages (Argelaguet: E4.5–E7.5; Pijuan-Sala: E6.5–E7.5) and total mRNA counts, highlighting the need for batch correction.

**Output**:
- Preprocessed `AnnData` objects (`adata_arg`, `adata_pj1`, `adata_pj2`) with filtered cells and genes, normalized counts, and HVGs.
- QC plots (violin, scatter, histograms) saved in the `figures/` directory.

### Part 2: Exploratory Data Analysis (EDA)

**Objective**: Perform dimensionality reduction and visualize the data to explore cell clusters and batch effects before integration.

1. **Parameter Optimization**:
   - Optimizes the number of principal components (`n_pcs`) using a PCA knee plot (targeting 70% cumulative variance).
   - Optimizes `n_neighbors` for UMAP using silhouette scores based on `celltype` labels.

2. **Dimensionality Reduction**:
   - Applies PCA, t-SNE, UMAP, Diffusion Maps, and Draw Graph to visualize the data for each dataset separately.
   - Performs Louvain clustering at resolutions 0.5 and 1.0 to identify cell clusters.

3. **Visualization**:
   - Generates plots for each dimensionality reduction method, colored by `celltype`, `stage`, and Louvain clusters.
   - Creates stacked bar plots to show cluster composition by cell type.
   - Discusses Louvain clusters, highlighting their alignment with developmental stages and cell types (e.g., Epiblast at E4.5, Mesoderm at E7.5).

**Output**:
- UMAP, t-SNE, PCA, Diffusion Map, and Draw Graph plots saved in the `figures/` directory.
- Stacked bar plots showing cluster composition.
- A markdown section discussing Louvain clusters and their biological relevance.

### Part 3: Dataset Integration with scVI

**Objective**: Integrate the datasets to correct for batch effects using scVI, a deep learning-based method.

1. **Combine Datasets**:
   - Concatenates the raw `AnnData` objects into a single `adata_combined` object, treating each dataset as a batch.

2. **Setup scVI**:
   - Adds mitochondrial fraction (`mt_frac`) as a continuous covariate.
   - Sets up scVI with `batch_key='batch'` and `continuous_covariate_keys=['mt_frac']`.

3. **Hyperparameter Tuning**:
   - Tests various scVI hyperparameters (`n_latent`, `n_layers`, `dispersion`, `gene_likelihood`) to find the best configuration.
   - Evaluates models using silhouette scores and training loss.
   - Visualizes results with bar plots to compare configurations.

4. **Train Final Model**:
   - Trains the final scVI model with the optimal hyperparameters for 200 epochs.
   - Saves the trained model for reproducibility.

5. **Interpret Latent Space**:
   - Extracts the latent representation and computes UMAP.
   - Visualizes the integrated space colored by `batch`, `celltype`, and `stage`.
   - Compares with pre-integration UMAPs to assess batch correction.

**Output**:
- Integrated `AnnData` object (`adata_combined`) with a latent representation (`X_scVI`).
- Hyperparameter tuning plots and integrated UMAP plots saved in the `figures/` directory.
- Trained scVI model saved in `scvi_model_final`.

### Part 4: Differential Expression Analysis

**Objective**: Identify differentially expressed genes for each developmental stage and interpret their biological significance.

1. **scVI DE Method**:
   - Performs DE analysis using scVI’s `differential_expression` method, leveraging the integrated latent space.
   - Generates volcano plots (using log fold change and Bayes factor) and a heatmap of normalized expression for the top 5 genes per stage.

2. **Scanpy Wilcoxon Method**:
   - Performs DE analysis using Scanpy’s `rank_genes_groups` with the Wilcoxon rank-sum test on raw counts.
   - Generates volcano plots (using log fold change and adjusted p-values) and a heatmap for the top 5 genes per stage.

3. **Biological Interpretation**:
   - Lists the top 5 differentially expressed genes per stage for both methods.
   - Provides biological insights into the top genes, linking them to gastrulation processes (e.g., pluripotency, germ layer formation).

**Output**:
- Volcano plots and heatmaps for both scVI and Wilcoxon methods, saved in the `figures/` directory.
- A markdown section with biological insights into the top genes.

## Running the Code

1. **Setup Environment**:
   - Create a Python environment (e.g., using `conda` or `venv`).
   - Install the dependencies listed above.

2. **Prepare Data**:
   - Ensure the `data/` directory contains the required files as described in the "Directory Structure" section.

3. **Run the Notebook**:
   - Launch Jupyter Notebook:
     ```bash
     jupyter notebook
     ```
   - Open `rna_eda.ipynb` and run the cells sequentially.

4. **Outputs**:
   - Plots (e.g., QC, UMAP, volcano, heatmaps) will be saved in the `figures/` directory.
   - The trained scVI model will be saved in `scvi_model_final`.

## Notes
- **File Paths**: Ensure the `figures/` directory exists or is writable, as plots are saved there.
- **Font Issues**: If the DejaVu Sans font is unavailable, the code will use the default font, which may affect plot appearance.
- **Training Time**: scVI training (Part 3) can be time-consuming with `max_epochs=200`. Reduce epochs (e.g., to 50) for faster testing.

## References

- Argelaguet et al. (2019). *Nature*, 576. DOI: [10.1038/s41586-019-1825-8](https://doi.org/10.1038/s41586-019-1825-8).
- Pijuan-Sala et al. (2019). *Nature*, 566. DOI: [10.1038/s41586-019-0933-9](https://doi.org/10.1038/s41586-019-0933-9).
- Lopez et al. (2018). *Nature Methods*, 15. DOI: [10.1038/s41592-018-0229-2](https://doi.org/10.1038/s41592-018-0229-2).
- Svensson (2020). *Nature Biotechnology*, 38. DOI: [10.1038/s41587-019-0379-5](https://doi.org/10.1038/s41587-019-0379-5).
- Theis Lab Single-Cell Tutorial: [Mouse Intestinal Epithelium](https://github.com/theislab/single-cell-tutorial).
