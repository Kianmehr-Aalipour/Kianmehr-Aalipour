# Univariate Survival Analysis with Multiple Testing Correction

This repository contains R scripts for performing univariate survival analysis on cell type signature scores with FDR (False Discovery Rate) correction for multiple testing.

## Overview

The analysis performs Kaplan-Meier survival analysis for multiple cell types, comparing high vs. low signature scores and correcting for multiple testing using the Benjamini-Hochberg method.

## Files

### Core Analysis Scripts

1. **`step4_univariate_survival_analysis.R`** - Main analysis script
   - Performs univariate survival analysis for multiple cell types
   - Creates Kaplan-Meier plots
   - Applies FDR correction for multiple testing
   - Generates summary tables

2. **`install_packages.R`** - Package installation helper
   - Installs all required R packages
   - Handles both CRAN and Bioconductor packages
   - Verifies successful installation

3. **`generate_sample_data.R`** - Sample data generator
   - Creates realistic sample TCGA-like data for testing
   - Generates clinical data and cell type scores
   - Useful for testing the analysis pipeline

## Required R Packages

### CRAN Packages
- `survival` - Survival analysis functions
- `survminer` - Enhanced survival plots
- `ggplot2` - Data visualization
- `dplyr` - Data manipulation

### Bioconductor Packages
- `SummarizedExperiment` - Bioconductor data structure

## Usage

### Option 1: Using Your Own Data

If you have your own TCGA data and cell type scores:

```r
# 1. Ensure you have the required objects in your environment:
#    - tcga_data: SummarizedExperiment object with clinical data
#    - estimated_scores: Matrix of cell type scores
#    - base_path: Directory path for saving results (optional)

# 2. Install required packages
source("install_packages.R")

# 3. Run the survival analysis
source("step4_univariate_survival_analysis.R")
```

### Option 2: Using Sample Data

For testing or demonstration purposes:

```r
# 1. Install required packages
source("install_packages.R")

# 2. Generate sample data
source("generate_sample_data.R")

# 3. Run the survival analysis
source("step4_univariate_survival_analysis.R")
```

## Expected Data Structure

### tcga_data (SummarizedExperiment)
The clinical data should contain these columns in `colData(tcga_data)`:
- `patient` - Patient identifiers
- `vital_status` - "Alive" or "Dead"
- `days_to_last_follow_up` - Days to last follow-up (for alive patients)
- `days_to_death` - Days to death (for deceased patients)

### estimated_scores (Matrix)
- Rows: Sample identifiers (TCGA barcode format)
- Columns: Cell type names
- Values: Numeric scores (typically 0-1 range)

## Output Files

The analysis generates several output files:

### Plots
- **`Survival_Plots_Univariate/`** directory containing:
  - Individual Kaplan-Meier plots for each cell type (PNG format)
  - High-resolution plots (300 DPI) suitable for publication

### Data Files
- **`Univariate_Survival_Summary_with_FDR.csv`** - Main results table containing:
  - Cell type names
  - Raw p-values
  - FDR-adjusted p-values
  - Significance indicators

## Cell Types Analyzed

The analysis focuses on these key cell types:
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

*Note: Only cell types present in your `estimated_scores` data will be analyzed.*

## Statistical Methods

### Survival Analysis
- **Method**: Kaplan-Meier survival estimation
- **Stratification**: Median split (high vs. low scores)
- **Statistical test**: Log-rank test
- **Plots**: Include confidence intervals and risk tables

### Multiple Testing Correction
- **Method**: Benjamini-Hochberg (FDR) correction
- **Significance threshold**: p < 0.05 (both raw and adjusted)

## Interpretation

### Results Table Columns
- **CellType**: Name of the cell type analyzed
- **P_Value**: Raw p-value from log-rank test
- **Adjusted_P_Value_FDR**: FDR-corrected p-value
- **Significant_Raw**: Significance at p < 0.05 (before correction)
- **Significant_FDR**: Significance at p < 0.05 (after FDR correction)

### Plot Interpretation
- **High Score**: Samples with cell type scores ≥ median
- **Low Score**: Samples with cell type scores < median
- **P-value**: Displayed on plot (raw p-value)
- **Risk Table**: Shows number of samples at risk over time

## Troubleshooting

### Common Issues

1. **Missing packages**: Run `install_packages.R` first
2. **Data format errors**: Ensure your data matches the expected structure
3. **No significant results**: This is normal - not all cell types may be prognostic
4. **Plot saving errors**: Check write permissions in the output directory

### Data Requirements
- Minimum sample size: ~50 samples recommended
- Clinical data completeness: Patients need either death date or follow-up date
- Score data: Should be numeric and not all identical values

## Example Output

```
SUMMARY OF UNIVARIATE SURVIVAL ANALYSIS RESULTS (WITH FDR CORRECTION)
======================================================================
     CellType P_Value Adjusted_P_Value_FDR Significant_Raw Significant_FDR
      T-Cells   0.023                0.087              Yes              No
    B-Cells   0.156                0.234              No               No
CD8+ T-Cells   0.045                0.112              Yes              No
...

ADDITIONAL SUMMARY STATISTICS:
--------------------------------------------------
Total cell types tested: 15
Significant at p < 0.05 (raw): 3
Significant at p < 0.05 (FDR-corrected): 1
```

## Contact

For questions or issues, please refer to the original analysis pipeline or contact the development team.

---

*This analysis is designed for research purposes. Please validate results with independent datasets and consider biological relevance alongside statistical significance.*