# Enhanced Forest Plot for Survival Analysis

## Overview
This repository contains improved forest plot visualizations for multivariate Cox regression analysis. The enhancements provide better visual design, additional information display, and multiple output formats.

## Key Improvements

### 1. Visual Enhancements
- **Color-coded significance**: Variables are colored based on statistical significance and effect direction
  - Red: Significant risk factors (HR > 1, p < 0.05)
  - Green: Significant protective factors (HR < 1, p < 0.05)  
  - Gray: Non-significant associations
- **Better typography**: Improved font sizes, bold text for important elements
- **Professional styling**: Clean, publication-ready appearance
- **Informative subtitles**: Include sample size and model performance metrics

### 2. Additional Information Display
- **Hazard ratios and confidence intervals**: Displayed as text alongside the plot
- **P-values**: Clearly shown with appropriate formatting (< 0.001 for very small values)
- **Model performance metrics**: C-index, R-squared, AIC/BIC
- **Variable importance ranking**: Based on effect size (absolute log HR)

### 3. Multiple Output Formats
- **Enhanced ggforest**: Improved version of the original survminer plot
- **Custom ggplot**: Fully customizable forest plot with text annotations
- **Professional forestplot**: Using the forestplot package for publication-quality figures
- **Summary table**: CSV export of all results

### 4. Additional Diagnostics
- **Proportional hazards testing**: Automatic testing of Cox model assumptions
- **Schoenfeld residual plots**: Generated if assumptions are violated
- **Comprehensive model summary**: Saved as RDS file for future analysis

## Quick Implementation

### Option 1: Simple Replacement
Replace your existing forest plot code:

```r
# OLD CODE:
forest_plot_adjusted <- ggforest(multi_cox_model_adjusted, data = survival_analysis_data)
ggsave(file.path(base_path, "Multivariate_Forest_Plot_ADJUSTED.png"), plot = forest_plot_adjusted, width = 12, height = 10)

# NEW CODE:
source("forest_plot_replacement.R")
```

### Option 2: Full Enhanced Workflow
For comprehensive analysis with all diagnostics:

```r
source("improved_step5_forest_plot.R")
```

### Option 3: Custom Function Approach
For maximum flexibility:

```r
source("improved_forest_plot.R")
enhanced_results <- enhanced_forest_plot_workflow(
  multi_cox_model_adjusted = your_cox_model,
  survival_analysis_data = your_data,
  base_path = your_output_path
)
```

## Files Generated

### Plots
- `*_enhanced_ggforest.png`: Improved version of original ggforest
- `*_custom_forest.png`: Custom designed forest plot with annotations
- `*_professional_forest.png`: Publication-quality forestplot
- `Schoenfeld_Residuals_Plot.png`: Diagnostic plot (if needed)

### Data Files
- `Forest_Plot_Summary_Table.csv`: Tabular summary of results
- `Variable_Importance_Ranking.csv`: Variables ranked by effect size
- `Model_Performance_Metrics.txt`: Model performance summary
- `Proportional_Hazards_Test.txt`: Assumption testing results
- `Comprehensive_Analysis_Summary.rds`: Complete analysis object

## Required Packages

```r
# Core packages
library(ggplot2)
library(survminer) 
library(survival)
library(dplyr)

# Additional packages for enhancements
library(broom)        # For model tidying
library(forestplot)   # For professional forest plots
```

## Features Comparison

| Feature | Original ggforest | Enhanced Version |
|---------|------------------|------------------|
| Basic forest plot | ✓ | ✓ |
| Color coding | ✗ | ✓ |
| HR/CI text display | ✗ | ✓ |
| P-value display | Limited | ✓ |
| Model metrics | ✗ | ✓ |
| Multiple formats | ✗ | ✓ |
| Summary table | ✗ | ✓ |
| Variable ranking | ✗ | ✓ |
| Assumption testing | ✗ | ✓ |
| Publication ready | Limited | ✓ |

## Customization Options

### Colors
Modify the color scheme by changing:
```r
scale_color_manual(
  values = c(
    "Significant Risk" = "#your_color",
    "Significant Protective" = "#your_color", 
    "Non-significant" = "#your_color"
  )
)
```

### Variable Labels
Customize variable names in the `variable_clean` case_when statement:
```r
variable_clean = case_when(
  term == "your_variable" ~ "Your Custom Label",
  # ... other mappings
)
```

### Plot Dimensions
Adjust plot size:
```r
ggsave(..., width = your_width, height = your_height, dpi = 300)
```

## Tips for Best Results

1. **Variable Names**: Ensure your variable names are descriptive and will display well
2. **Sample Size**: Include sample size information in plot titles/subtitles
3. **Model Performance**: Always report C-index and other performance metrics
4. **Assumptions**: Check proportional hazards assumptions and report violations
5. **Effect Sizes**: Consider reporting both statistical significance and clinical significance
6. **Multiple Comparisons**: Apply appropriate corrections if testing many variables

## Troubleshooting

### Common Issues
1. **Package installation**: Some packages may need to be installed
2. **Variable names**: Clean variable names for better display
3. **Large datasets**: May need to adjust plot dimensions for many variables
4. **Memory issues**: Use data.table or similar for very large datasets

### Performance Tips
- Use parallel processing for large models
- Save intermediate results as RDS files
- Consider variable selection to reduce model complexity

## Citation
If you use these enhancements in your research, please consider citing the original packages:
- survminer: Kassambara et al.
- forestplot: Gordon & Lumley
- ggplot2: Wickham et al.
- survival: Therneau & Grambsch