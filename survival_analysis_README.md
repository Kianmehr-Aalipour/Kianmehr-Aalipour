# Univariate Survival Analysis with Multiple Testing Correction

This repository contains R code for performing univariate survival analysis on TCGA (The Cancer Genome Atlas) data with cell type signature scores, including proper multiple testing correction.

## Overview

The analysis performs the following steps:

1. **Data Preparation**: Extracts clinical survival data from TCGA SummarizedExperiment object
2. **Score Integration**: Merges survival data with cell type signature scores
3. **Univariate Analysis**: Performs Kaplan-Meier survival analysis for each cell type
4. **Multiple Testing Correction**: Applies False Discovery Rate (FDR) correction using Benjamini-Hochberg method
5. **Visualization**: Generates survival plots for each cell type
6. **Results Export**: Saves summary statistics and corrected p-values

## Requirements

### R Packages

#### CRAN Packages:
- `survival` - Core survival analysis functions
- `survminer` - Enhanced survival plots and statistics  
- `ggplot2` - Plotting system
- `dplyr` - Data manipulation

#### Bioconductor Packages:
- `SummarizedExperiment` - For handling TCGA data structures

### Installation

Run the requirements script to install all necessary packages:

```r
source("requirements.R")
```

Or install manually:

```r
# CRAN packages
install.packages(c("survival", "survminer", "ggplot2", "dplyr"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")
```

## Input Data Requirements

The script expects the following objects to be available in the R environment:

1. **`tcga_data`**: A SummarizedExperiment object containing TCGA data with clinical information in `colData()` including:
   - `patient`: Patient identifiers
   - `vital_status`: Vital status ("Alive" or "Dead")
   - `days_to_last_follow_up`: Days to last follow-up for living patients
   - `days_to_death`: Days to death for deceased patients

2. **`estimated_scores`**: A matrix or data frame with cell type signature scores where:
   - Rows are samples (with TCGA barcode format)
   - Columns are cell types
   - Values are numeric scores

3. **`base_path`**: Character string specifying the base directory for saving outputs

## Usage

### Method 1: As part of a larger pipeline
If you have the required objects already loaded:

```r
source("step4_univariate_survival_analysis.R")
```

### Method 2: Standalone execution
If running independently, load your data first:

```r
# Load your data
load("tcga_data.RData")
load("estimated_scores.RData")
base_path <- "/path/to/output/directory"

# Run the analysis
source("step4_univariate_survival_analysis.R")
```

## Cell Types Analyzed

The script analyzes the following cell types (if available in the scores data):

- T-Cells
- B-Cells  
- CD8+ T-Cells
- NK Cells
- Macrophages
- Monocytes
- Neutrophils
- Mast Cells
- Dendritic Cells
- Fibroblasts
- Endothelial Cells
- AT1 Cells
- AT2 Cells
- Plasma Cells
- Regulatory T-Cells (Tregs)

## Output Files

The script generates the following outputs:

### 1. Survival Plots
- **Location**: `{base_path}/Survival_Plots_Univariate/`
- **Format**: PNG files (300 DPI, 8×7 inches)
- **Content**: Kaplan-Meier curves with risk tables for each cell type
- **Naming**: `Survival_Plot_{CellType}.png`

### 2. Summary Results Table
- **File**: `{base_path}/Univariate_Survival_Summary_with_FDR.csv`
- **Columns**:
  - `CellType`: Name of the cell type
  - `P_Value`: Raw p-value from log-rank test
  - `Adjusted_P_Value_FDR`: FDR-corrected p-value (Benjamini-Hochberg)
  - `Significant_Raw`: Boolean indicating significance at p < 0.05
  - `Significant_FDR`: Boolean indicating significance after FDR correction

## Statistical Methods

### Survival Analysis
- **Method**: Kaplan-Meier estimator with log-rank test
- **Stratification**: Median split (High Score vs. Low Score)
- **Test**: Log-rank test for comparing survival curves

### Multiple Testing Correction
- **Method**: False Discovery Rate (FDR) correction
- **Implementation**: Benjamini-Hochberg procedure via `p.adjust(method = "BH")`
- **Significance Level**: α = 0.05

## Interpretation

### P-values
- **Raw P-value**: Significance of survival difference for individual cell type
- **FDR-corrected P-value**: Significance after controlling for multiple comparisons

### Survival Curves
- **High Score**: Samples with cell type scores ≥ median
- **Low Score**: Samples with cell type scores < median
- **Risk Tables**: Number of patients at risk at different time points

## Example Output

```
Cell types available for analysis: 12
Cell types to be tested:
  - T-Cells
  - B-Cells
  - CD8+ T-Cells
  - NK Cells
  - Macrophages
  ...

--- Summary of Univariate Survival Analysis Results (with FDR) ---
     CellType P_Value Adjusted_P_Value_FDR Significant_Raw Significant_FDR
1    T-Cells   0.0234               0.1404           TRUE           FALSE
2    B-Cells   0.0456               0.1825           TRUE           FALSE
...

Summary Statistics:
  Total cell types tested: 12
  Significant at p < 0.05: 3
  Significant after FDR correction: 1
```

## Notes

- The script includes error handling for individual cell type analyses
- Missing or invalid data points are automatically excluded
- The analysis assumes TCGA barcode format for sample identifiers
- Plots use a consistent color scheme (gold for high scores, blue for low scores)

## Troubleshooting

### Common Issues

1. **Missing columns in clinical data**: Ensure your TCGA data has the required clinical variables
2. **Sample ID mismatch**: Check that sample IDs in scores match TCGA barcodes (first 12 characters)
3. **No significant results**: This may indicate weak associations or need for different stratification methods

### Contact

For questions or issues, please contact: aalipourez@gmail.com