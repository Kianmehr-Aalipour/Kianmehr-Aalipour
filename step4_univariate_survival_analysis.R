# =============================================================================
# STEP 4: UNIVARIATE SURVIVAL ANALYSIS (WITH MULTIPLE TESTING CORRECTION)
# =============================================================================

# Load required libraries
suppressMessages({
  library(SummarizedExperiment)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
})

cat("Performing univariate survival analysis for key cell types...\n")

# =============================================================================
# 4.1: PREPARE DATA FOR SURVIVAL ANALYSIS
# =============================================================================

# Assuming tcga_data and estimated_scores are loaded from previous steps
# If not available, you may need to load them from saved files:
# load("tcga_data.RData")
# load("estimated_scores.RData")

# Extract clinical data from SummarizedExperiment object
clinical_data_raw <- as.data.frame(SummarizedExperiment::colData(tcga_data))

# Select only the necessary columns for survival
survival_info <- clinical_data_raw[, c("patient", "vital_status", "days_to_last_follow_up", "days_to_death")]
colnames(survival_info) <- c("submitter_id", "status", "time_to_last_contact", "time_to_death")

# Create unified time variable
survival_info$time <- ifelse(survival_info$status == "Alive", 
                           survival_info$time_to_last_contact, 
                           survival_info$time_to_death)

# Remove rows with missing time data
survival_info <- survival_info[!is.na(survival_info$time), ]

# Create numeric status variable (0 = alive, 1 = dead)
survival_info$status_numeric <- ifelse(survival_info$status == "Alive", 0, 1)

# Prepare scores dataframe
scores_df <- as.data.frame(estimated_scores)
scores_df$submitter_id <- substr(rownames(scores_df), 1, 12)

# Merge survival data with scores
survival_analysis_data <- merge(survival_info, scores_df, by = "submitter_id")

cat("Survival analysis data prepared with", nrow(survival_analysis_data), "samples\n")

# =============================================================================
# 4.2: DEFINE LIST OF CELL TYPES TO TEST
# =============================================================================

# Define key cell types of interest
key_cell_types_of_interest <- c(
  "T-Cells", "B-Cells", "CD8+ T-Cells", "NK Cells", "Macrophages", "Monocytes", 
  "Neutrophils", "Mast Cells", "Dendritic Cells", "Fibroblasts", 
  "Endothelial Cells", "AT1 Cells", "AT2 Cells", "Plasma Cells", "Regulatory T-Cells (Tregs)"
)

# Get available scores from the estimated_scores object
available_scores <- colnames(estimated_scores)

# Find intersection of desired cell types and available scores
valid_cell_types_to_test <- intersect(key_cell_types_of_interest, available_scores)

cat("Found", length(valid_cell_types_to_test), "valid cell types to test:\n")
cat(paste(valid_cell_types_to_test, collapse = ", "), "\n\n")

# =============================================================================
# 4.3: PERFORM SURVIVAL ANALYSIS AND CREATE PLOTS
# =============================================================================

# Create directory for survival plots
# Assuming base_path is defined, otherwise use current directory
if (!exists("base_path")) {
  base_path <- getwd()
}

plot_dir <- file.path(base_path, "Survival_Plots_Univariate")
if (!dir.exists(plot_dir)) { 
  dir.create(plot_dir, recursive = TRUE)
  cat("Created directory:", plot_dir, "\n")
}

# Initialize results storage
univariate_results <- list()

# Loop through each valid cell type
for (cell_type in valid_cell_types_to_test) {
  cat("  Analyzing survival for:", cell_type, "\n")
  
  # Calculate median score for stratification
  median_score <- median(survival_analysis_data[[cell_type]], na.rm = TRUE)
  
  # Create high/low score strata based on median
  survival_analysis_data$strata <- ifelse(survival_analysis_data[[cell_type]] >= median_score, 
                                        "High Score", "Low Score")
  
  # Perform survival analysis
  fit <- survfit(Surv(time, status_numeric) ~ strata, data = survival_analysis_data)
  
  # Calculate p-value
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
    conf.int = TRUE,
    palette = c("#E7B800", "#2E9FDF"),
    ggtheme = theme_minimal()
  )
  
  # Save plot
  plot_filename <- file.path(plot_dir, paste0("Survival_Plot_", gsub("[^A-Za-z0-9]", "_", cell_type), ".png"))
  
  # Save plot with error handling
  tryCatch({
    png(filename = plot_filename, width = 8, height = 7, units = "in", res = 300)
    print(km_plot)
    dev.off()
    cat("    Plot saved:", basename(plot_filename), "\n")
  }, error = function(e) {
    cat("    Error saving plot for", cell_type, ":", e$message, "\n")
    if (dev.cur() != 1) dev.off()  # Close any open graphics device
  })
}

# =============================================================================
# 4.4: MULTIPLE TESTING CORRECTION AND RESULTS SUMMARY
# =============================================================================

# Extract p-values
original_p_values <- unlist(univariate_results)

# Apply FDR correction (Benjamini-Hochberg)
adjusted_p_values_fdr <- p.adjust(original_p_values, method = "BH")

# Create results summary dataframe
survival_summary_df <- data.frame(
  CellType = names(original_p_values),
  P_Value = original_p_values,
  Adjusted_P_Value_FDR = adjusted_p_values_fdr,
  Significant_Raw = ifelse(original_p_values < 0.05, "Yes", "No"),
  Significant_FDR = ifelse(adjusted_p_values_fdr < 0.05, "Yes", "No"),
  stringsAsFactors = FALSE
)

# Sort by adjusted p-value
survival_summary_df <- survival_summary_df[order(survival_summary_df$Adjusted_P_Value_FDR), ]

# Display results
cat("\n" + paste(rep("=", 70), collapse = "") + "\n")
cat("SUMMARY OF UNIVARIATE SURVIVAL ANALYSIS RESULTS (WITH FDR CORRECTION)\n")
cat(paste(rep("=", 70), collapse = "") + "\n")
print(survival_summary_df, row.names = FALSE)

# Save results to CSV
results_filename <- file.path(base_path, "Univariate_Survival_Summary_with_FDR.csv")
write.csv(survival_summary_df, file = results_filename, row.names = FALSE)
cat("\nResults saved to:", results_filename, "\n")

# Additional summary statistics
cat("\n" + paste(rep("-", 50), collapse = "") + "\n")
cat("ADDITIONAL SUMMARY STATISTICS:\n")
cat(paste(rep("-", 50), collapse = "") + "\n")
cat("Total cell types tested:", nrow(survival_summary_df), "\n")
cat("Significant at p < 0.05 (raw):", sum(survival_summary_df$Significant_Raw == "Yes"), "\n")
cat("Significant at p < 0.05 (FDR-corrected):", sum(survival_summary_df$Significant_FDR == "Yes"), "\n")

if (sum(survival_summary_df$Significant_FDR == "Yes") > 0) {
  cat("\nCell types significant after FDR correction:\n")
  significant_types <- survival_summary_df[survival_summary_df$Significant_FDR == "Yes", "CellType"]
  cat(paste(significant_types, collapse = ", "), "\n")
}

cat("\nUnivariate survival analysis completed successfully!\n")