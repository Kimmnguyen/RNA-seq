# Drought Tolerance Prediction Using Motif and Expression Analysis

This repository contains all the code, data, and analysis scripts used in the study of drought tolerance prediction across Arabidopsis accessions. The project integrates gene expression data, motif occurrence analysis, and machine learning models to identify drought-tolerant plant accessions based on sequence and transcriptomic features.

## 📁 Repository Contents

### 🔬 Data Files
- `TPM10.txt`: Original gene expression dataset, filtered to include genes with TPM > 10.
- `gene_motifs dataset.csv`: Merged dataset linking genes to significant motifs.
- `17_accession_drought_category and.csv`: Classification of 17 accessions based on observed drought responses.
- `xgboost_predictions_drought_tolerance.csv`: Predicted drought tolerance scores from XGBoost model.
- `Accession Distribution.html`: Visualization of drought response distributions across accessions.

### 📜 Scripts
- `Accession_XGboosting.py`: Python script for building and evaluating an XGBoost model to predict drought tolerance.
- `DSeq2.R`: Differential expression analysis script using the TPM-filtered gene set.
- `Motifs matching by Fimo and MCAST.R`: Detects regulatory motifs using FIMO and MCAST tools.
- `Visualization for clustering.R`: Generates clustering plots based on motif or expression profiles.
- `motif occur.R`: Analyzes motif occurrences across different accessions.
- `geographical & rainfall.R`: Integrates environmental data like rainfall and geography into accession profiles.

### 📄 Documentation
- `README.md`: This file. Provides overview, usage instructions, and descriptions.

## 🚀 Getting Started

### Prerequisites
- R (with packages: `ggplot2`, `dplyr`, `DESeq2`, etc.)
- Python 3 (with packages: `xgboost`, `pandas`, `sklearn`, `matplotlib`, etc.)
- MEME Suite (for FIMO and MCAST)

### Running the Analysis

1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/your-repo-name.git
   cd your-repo-name
