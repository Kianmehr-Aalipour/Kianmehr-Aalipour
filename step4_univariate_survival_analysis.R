# ---
# STEP 4: UNIVARIATE SURVIVAL ANALYSIS (WITH MULTIPLE TESTING CORRECTION)
# ---

# Load required libraries
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
})

cat("Performing univariate survival analysis for key cell types...\n")

# Note: This script assumes that the following variables are available from previous steps:
# - tcga_data: SummarizedExperiment object with TCGA data
# - estimated_scores: Matrix/data frame with cell type scores
# - base_path: Base directory path for saving outputs

# If running this script independently, you'll need to load these objects first
# Example:
# load("tcga_data.RData")
# load("estimated_scores.RData")
# base_path <- getwd()

# 4.1: Prepare data for survival analysis
clinical_data_raw <- as.data.frame(SummarizedExperiment::colData(tcga_data))

# Select only the necessary columns for survival
survival_info <- clinical_data_raw[, c("patient", "vital_status", "days_to_last_follow_up", "days_to_death")]
colnames(survival_info) <- c("submitter_id", "status", "time_to_last_contact", "time_to_death")

# Create unified time variable
survival_info$time <- ifelse(survival_info$status == "Alive", 
                           survival_info$time_to_last_contact, 
                           survival_info$time_to_death)

# Remove cases with missing time data
survival_info <- survival_info[!is.na(survival_info$time), ]

# Create numeric status variable (0 = censored/alive, 1 = event/dead)
survival_info$status_numeric <- ifelse(survival_info$status == "Alive", 0, 1)

# Prepare scores data frame
scores_df <- as.data.frame(estimated_scores)
scores_df$submitter_id <- substr(rownames(scores_df), 1, 12)

# Merge survival data with scores
survival_analysis_data <- merge(survival_info, scores_df, by = "submitter_id")

cat("Merged dataset contains", nrow(survival_analysis_data), "samples\n")

# 4.2: Define list of cell types to test
key_cell_types_of_interest <- c(
  "T-Cells", "B-Cells", "CD8+ T-Cells", "NK Cells", "Macrophages", "Monocytes", 
  "Neutrophils", "Mast Cells", "Dendritic Cells", "Fibroblasts", 
  "Endothelial Cells", "AT1 Cells", "AT2 Cells", "Plasma Cells", "Regulatory T-Cells (Tregs)"
)

available_scores <- colnames(estimated_scores)
valid_cell_types_to_test <- intersect(key_cell_types_of_interest, available_scores)

cat("Cell types available for analysis:", length(valid_cell_types_to_test), "\n")
cat("Cell types to be tested:\n")
for (ct in valid_cell_types_to_test) {
  cat("  -", ct, "\n")
}

# 4.3: Create output directory and perform survival analysis
plot_dir <- file.path(base_path, "Survival_Plots_Univariate")
if (!dir.exists(plot_dir)) { 
  dir.create(plot_dir, recursive = TRUE)
  cat("Created directory:", plot_dir, "\n")
}

# Initialize results storage
univariate_results <- list()

# Loop through each cell type and perform survival analysis
for (cell_type in valid_cell_types_to_test) {
  cat("  Analyzing survival for:", cell_type, "\n")
  
  # Dichotomize based on median score
  median_score <- median(survival_analysis_data[[cell_type]], na.rm = TRUE)
  survival_analysis_data$strata <- ifelse(survival_analysis_data[[cell_type]] >= median_score, 
                                        "High Score", "Low Score")
  
  # Perform survival analysis
  tryCatch({
    # Fit Kaplan-Meier survival curves
    fit <- survfit(Surv(time, status_numeric) ~ strata, data = survival_analysis_data)
    
    # Calculate p-value using log-rank test
    p_val <- surv_pvalue(fit)$pval
    univariate_results[[cell_type]] <- p_val
    
    # Create Kaplan-Meier plot
    km_plot <- ggsurvplot(
      fit, 
      data = survival_analysis_data, 
      pval = TRUE, 
      risk.table = TRUE,
      title = paste("Survival Analysis based on", cell_type, "Signature Score"),
      legend.title = paste(cell_type, "Score"), 
      xlab = "Time (days)",
      ylab = "Survival Probability",
      palette = c("#E7B800", "#2E9FDF"),  # Custom colors for high/low groups
      ggtheme = theme_minimal(),
      risk.table.height = 0.3
    )
    
    # Save plot
    plot_filename <- file.path(plot_dir, paste0("Survival_Plot_", gsub("[^A-Za-z0-9]", "_", cell_type), ".png"))
    png(filename = plot_filename, width = 8, height = 7, units = "in", res = 300)
    print(km_plot)
    dev.off()
    
    cat("    P-value:", sprintf("%.4f", p_val), "- Plot saved\n")
    
  }, error = function(e) {
    cat("    Error analyzing", cell_type, ":", e$message, "\n")
    univariate_results[[cell_type]] <- NA
  })
}

# 4.4: Create and save the results table with multiple testing correction
cat("\nApplying multiple testing correction...\n")

# Extract p-values and remove any NA values
original_p_values <- unlist(univariate_results)
valid_p_values <- original_p_values[!is.na(original_p_values)]

if (length(valid_p_values) > 0) {
  # Apply FDR correction (Benjamini-Hochberg)
  adjusted_p_values_fdr <- p.adjust(original_p_values, method = "BH")
  
  # Create results summary table
  survival_summary_df <- data.frame(
    CellType = names(original_p_values),
    P_Value = original_p_values,
    Adjusted_P_Value_FDR = adjusted_p_values_fdr,
    Significant_Raw = ifelse(is.na(original_p_values), NA, original_p_values < 0.05),
    Significant_FDR = ifelse(is.na(adjusted_p_values_fdr), NA, adjusted_p_values_fdr < 0.05),
    stringsAsFactors = FALSE
  )
  
  # Sort by raw p-value
  survival_summary_df <- survival_summary_df[order(survival_summary_df$P_Value, na.last = TRUE), ]
  
  # Display results
  cat("\n--- Summary of Univariate Survival Analysis Results (with FDR) ---\n")
  print(survival_summary_df)
  
  # Save results to CSV
  output_file <- file.path(base_path, "Univariate_Survival_Summary_with_FDR.csv")
  write.csv(survival_summary_df, file = output_file, row.names = FALSE)
  cat("\nResults saved to:", output_file, "\n")
  
  # Print summary statistics
  n_significant_raw <- sum(survival_summary_df$Significant_Raw, na.rm = TRUE)
  n_significant_fdr <- sum(survival_summary_df$Significant_FDR, na.rm = TRUE)
  
  cat("\nSummary Statistics:\n")
  cat("  Total cell types tested:", nrow(survival_summary_df), "\n")
  cat("  Significant at p < 0.05:", n_significant_raw, "\n")
  cat("  Significant after FDR correction:", n_significant_fdr, "\n")
  
  if (n_significant_fdr > 0) {
    cat("\nCell types significant after FDR correction:\n")
    significant_types <- survival_summary_df[survival_summary_df$Significant_FDR == TRUE & !is.na(survival_summary_df$Significant_FDR), ]
    for (i in 1:nrow(significant_types)) {
      cat(sprintf("  - %s (p = %.4f, FDR = %.4f)\n", 
                  significant_types$CellType[i], 
                  significant_types$P_Value[i], 
                  significant_types$Adjusted_P_Value_FDR[i]))
    }
  }
  
} else {
  cat("No valid p-values obtained from survival analysis.\n")
}

cat("\nUnivariate survival analysis completed!\n")