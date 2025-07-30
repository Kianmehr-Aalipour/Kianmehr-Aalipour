# ---
# STEP 5: MULTIVARIATE SURVIVAL ANALYSIS (ADJUSTED FOR CLINICAL COVARIATES) - IMPROVED
# ---
cat("\nPerforming multivariate survival analysis adjusted for clinical covariates...\n")

# Load required libraries for enhanced forest plots
library(ggplot2)
library(survminer)
library(survival)
library(dplyr)
library(forestplot)
library(broom)

# Source the enhanced forest plot functions
source("improved_forest_plot.R")

# 5.1: Get additional clinical data and merge it
cat("  Fetching and merging clinical covariates (age, gender, stage)...\n")
clinical_covariates <- clinical_data_raw[, c("patient", "age_at_index", "gender", "ajcc_pathologic_stage")]
colnames(clinical_covariates)[1] <- "submitter_id"

# Merge with the main survival analysis data frame
survival_analysis_data <- merge(survival_analysis_data, clinical_covariates, by = "submitter_id", all.x = TRUE)

# 5.2: Clean and simplify tumor stage
# Grouping stages is crucial for model stability.
survival_analysis_data$tumor_stage_simple <- case_when(
  grepl("Stage IV", survival_analysis_data$ajcc_pathologic_stage, ignore.case = TRUE) ~ "Stage IV",
  grepl("Stage III", survival_analysis_data$ajcc_pathologic_stage, ignore.case = TRUE) ~ "Stage III",
  grepl("Stage II", survival_analysis_data$ajcc_pathologic_stage, ignore.case = TRUE) ~ "Stage II",
  grepl("Stage I", survival_analysis_data$ajcc_pathologic_stage, ignore.case = TRUE) ~ "Stage I",
  TRUE ~ NA_character_
)
survival_analysis_data <- survival_analysis_data[!is.na(survival_analysis_data$tumor_stage_simple), ]
survival_analysis_data$tumor_stage_simple <- factor(survival_analysis_data$tumor_stage_simple, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
cat("  Simplified tumor stage. Patients remaining after removing NA stages:", nrow(survival_analysis_data), "\n")

# 5.3: Build the comprehensive multivariate Cox model formula
cols_for_model <- make.names(valid_cell_types_to_test)
colnames(survival_analysis_data) <- make.names(colnames(survival_analysis_data))

# Formula now includes cell scores AND clinical covariates
cox_formula_adjusted <- as.formula(paste(
  "Surv(time, status_numeric) ~", 
  paste(cols_for_model, collapse = " + "),
  "+ age_at_index + gender + tumor_stage_simple"
))

# 5.4: Fit the adjusted model
multi_cox_model_adjusted <- coxph(cox_formula_adjusted, data = survival_analysis_data)

cat("\n--- Summary of ADJUSTED Multivariate Cox Model ---\n")
print(summary(multi_cox_model_adjusted))

# 5.5: Save the adjusted model summary to a new file
model_summary_text_adj <- capture.output(summary(multi_cox_model_adjusted))
writeLines(model_summary_text_adj, file.path(base_path, "Multivariate_Cox_Model_ADJUSTED_Summary.txt"))

# 5.6: ENHANCED FOREST PLOT VISUALIZATION
cat("Generating ENHANCED forest plots for the ADJUSTED multivariate analysis...\n")

# Replace the basic ggforest with our enhanced version
enhanced_results <- enhanced_forest_plot_workflow(
  multi_cox_model_adjusted = multi_cox_model_adjusted,
  survival_analysis_data = survival_analysis_data,
  base_path = base_path
)

# Additional model diagnostics and visualizations
cat("\n--- Enhanced Model Diagnostics ---\n")

# 1. Check proportional hazards assumption
cat("Testing proportional hazards assumption...\n")
ph_test <- cox.zph(multi_cox_model_adjusted)
print(ph_test)

# Save PH test results
ph_test_output <- capture.output(print(ph_test))
writeLines(ph_test_output, file.path(base_path, "Proportional_Hazards_Test.txt"))

# Plot Schoenfeld residuals if any violations
if(any(ph_test$table[,"p"] < 0.05)) {
  cat("  Warning: Proportional hazards assumption may be violated for some variables.\n")
  cat("  Generating Schoenfeld residual plots...\n")
  
  png(file.path(base_path, "Schoenfeld_Residuals_Plot.png"), width = 12, height = 8, units = "in", res = 300)
  plot(ph_test, main = "Schoenfeld Residuals Test for Proportional Hazards")
  dev.off()
}

# 2. Model performance summary
model_performance <- list(
  concordance = enhanced_results$model_metrics$concordance,
  rsquared = enhanced_results$model_metrics$rsquared,
  lr_pvalue = enhanced_results$model_metrics$lr_pvalue,
  aic = AIC(multi_cox_model_adjusted),
  bic = BIC(multi_cox_model_adjusted)
)

cat("\n--- Model Performance Metrics ---\n")
cat("Concordance Index (C-index):", round(model_performance$concordance, 3), "\n")
cat("R-squared:", round(model_performance$rsquared, 3), "\n")
cat("AIC:", round(model_performance$aic, 1), "\n")
cat("BIC:", round(model_performance$bic, 1), "\n")
cat("Likelihood Ratio Test p-value:", format.pval(model_performance$lr_pvalue, digits = 3), "\n")

# Save model performance metrics
performance_text <- c(
  "=== Model Performance Metrics ===",
  paste("Concordance Index (C-index):", round(model_performance$concordance, 3)),
  paste("R-squared:", round(model_performance$rsquared, 3)),
  paste("AIC:", round(model_performance$aic, 1)),
  paste("BIC:", round(model_performance$bic, 1)),
  paste("Likelihood Ratio Test p-value:", format.pval(model_performance$lr_pvalue, digits = 3)),
  "",
  "=== Interpretation ===",
  "C-index: Measure of discriminative ability (0.5 = random, 1.0 = perfect)",
  "R-squared: Proportion of variation explained by the model",
  "AIC/BIC: Model selection criteria (lower is better)",
  "LR test: Overall model significance"
)
writeLines(performance_text, file.path(base_path, "Model_Performance_Metrics.txt"))

# 3. Variable importance ranking
cat("\n--- Variable Importance Ranking ---\n")
variable_importance <- enhanced_results$plots$forest_data %>%
  mutate(
    importance_score = abs(log(HR)),
    risk_direction = ifelse(HR > 1, "Risk Factor", "Protective Factor")
  ) %>%
  arrange(desc(importance_score)) %>%
  select(variable_clean, HR, p_value, importance_score, risk_direction)

print(variable_importance)

# Save variable importance
write.csv(variable_importance, file.path(base_path, "Variable_Importance_Ranking.csv"), row.names = FALSE)

# 4. Create a comprehensive results summary
comprehensive_summary <- list(
  model_formula = deparse(cox_formula_adjusted),
  sample_size = nrow(survival_analysis_data),
  number_of_events = sum(survival_analysis_data$status_numeric),
  follow_up_time = list(
    median = median(survival_analysis_data$time),
    range = range(survival_analysis_data$time)
  ),
  model_performance = model_performance,
  significant_variables = enhanced_results$summary_table[enhanced_results$summary_table$`P-value` != "" & 
                                                        enhanced_results$summary_table$Significance %in% c("*", "**", "***"), ],
  ph_assumption_test = ph_test$table
)

# Save comprehensive summary as RDS for future use
saveRDS(comprehensive_summary, file.path(base_path, "Comprehensive_Analysis_Summary.rds"))

cat("\n--- ENHANCED Prognostic analysis workflow complete. ---\n")
cat("Enhanced forest plots (3 versions) are saved with improved formatting.\n")
cat("Summary table and model performance metrics are saved.\n")
cat("Proportional hazards assumption has been tested.\n")
cat("Variable importance ranking has been generated.\n")
cat("Comprehensive analysis summary saved as RDS file.\n")

# Return the enhanced results for further use if needed
return(enhanced_results)