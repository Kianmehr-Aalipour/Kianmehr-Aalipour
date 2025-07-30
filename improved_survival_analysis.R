# ---
# STEP 5: MULTIVARIATE SURVIVAL ANALYSIS (ADJUSTED FOR CLINICAL COVARIATES)
# ---
cat("\nPerforming multivariate survival analysis adjusted for clinical covariates...\n")

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

# 5.6: Create improved forest plots with multiple options
cat("Generating improved forest plots for the ADJUSTED multivariate analysis...\n")

# Load required libraries for improved plotting
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

# Function to create an improved forest plot
create_improved_forest_plot <- function(cox_model, data, 
                                       title = "Multivariate Cox Regression Analysis",
                                       subtitle = "Adjusted for Clinical Covariates",
                                       ref_line = 1,
                                       conf_level = 0.95) {
  
  # Extract model results
  model_summary <- summary(cox_model)
  coef_table <- model_summary$conf.int
  
  # Create data frame for plotting
  forest_data <- data.frame(
    variable = rownames(coef_table),
    hr = coef_table[, 1],
    lower_ci = coef_table[, 3],
    upper_ci = coef_table[, 4],
    p_value = model_summary$coefficients[, 5]
  )
  
  # Add significance levels
  forest_data$significance <- case_when(
    forest_data$p_value < 0.001 ~ "***",
    forest_data$p_value < 0.01 ~ "**",
    forest_data$p_value < 0.05 ~ "*",
    TRUE ~ ""
  )
  
  # Create variable labels with significance
  forest_data$variable_label <- paste0(forest_data$variable, " ", forest_data$significance)
  
  # Add confidence interval text
  forest_data$ci_text <- sprintf("%.2f (%.2f-%.2f)", 
                                forest_data$hr, 
                                forest_data$lower_ci, 
                                forest_data$upper_ci)
  
  # Add p-value text
  forest_data$p_text <- sprintf("p = %.3f", forest_data$p_value)
  
  # Color coding based on significance
  forest_data$color_group <- case_when(
    forest_data$p_value < 0.001 ~ "Highly Significant",
    forest_data$p_value < 0.01 ~ "Significant", 
    forest_data$p_value < 0.05 ~ "Marginally Significant",
    TRUE ~ "Not Significant"
  )
  
  # Create the forest plot
  p <- ggplot(forest_data, aes(x = hr, y = reorder(variable_label, hr))) +
    # Add reference line
    geom_vline(xintercept = ref_line, linetype = "dashed", color = "gray50", alpha = 0.7) +
    
    # Add confidence intervals
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, color = color_group), 
                   height = 0.3, size = 1.2) +
    
    # Add hazard ratio points
    geom_point(aes(color = color_group), size = 3, shape = 16) +
    
    # Color scale
    scale_color_manual(values = c("Highly Significant" = "#D32F2F",
                                 "Significant" = "#FF5722", 
                                 "Marginally Significant" = "#FF9800",
                                 "Not Significant" = "#757575")) +
    
    # Scale and labels
    scale_x_continuous(breaks = seq(0, max(forest_data$upper_ci, na.rm = TRUE) + 0.5, 0.5),
                      limits = c(0, max(forest_data$upper_ci, na.rm = TRUE) + 0.5)) +
    
    # Labels
    labs(title = title,
         subtitle = subtitle,
         x = "Hazard Ratio (95% CI)",
         y = NULL,
         color = "Significance Level") +
    
    # Theme
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50"),
      axis.text.y = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, size = 0.5)
    )
  
  return(p)
}

# Function to create a comprehensive forest plot with model diagnostics
create_comprehensive_forest_plot <- function(cox_model, data,
                                           title = "Comprehensive Multivariate Cox Analysis",
                                           subtitle = "Adjusted for Clinical Covariates") {
  
  # Extract model results
  model_summary <- summary(cox_model)
  coef_table <- model_summary$conf.int
  
  # Create data frame for plotting
  forest_data <- data.frame(
    variable = rownames(coef_table),
    hr = coef_table[, 1],
    lower_ci = coef_table[, 3],
    upper_ci = coef_table[, 4],
    p_value = model_summary$coefficients[, 5],
    se = model_summary$coefficients[, 3]
  )
  
  # Add significance levels
  forest_data$significance <- case_when(
    forest_data$p_value < 0.001 ~ "***",
    forest_data$p_value < 0.01 ~ "**",
    forest_data$p_value < 0.05 ~ "*",
    TRUE ~ ""
  )
  
  # Create variable labels with significance
  forest_data$variable_label <- paste0(forest_data$variable, " ", forest_data$significance)
  
  # Add confidence interval text
  forest_data$ci_text <- sprintf("%.2f (%.2f-%.2f)", 
                                forest_data$hr, 
                                forest_data$lower_ci, 
                                forest_data$upper_ci)
  
  # Add p-value text
  forest_data$p_text <- sprintf("p = %.3f", forest_data$p_value)
  
  # Color coding based on significance
  forest_data$color_group <- case_when(
    forest_data$p_value < 0.001 ~ "Highly Significant",
    forest_data$p_value < 0.01 ~ "Significant", 
    forest_data$p_value < 0.05 ~ "Marginally Significant",
    TRUE ~ "Not Significant"
  )
  
  # Create the main forest plot with enhanced styling
  p1 <- ggplot(forest_data, aes(x = hr, y = reorder(variable_label, hr))) +
    # Add reference line
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
    
    # Add confidence intervals with gradient
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, color = color_group), 
                   height = 0.3, size = 1.2) +
    
    # Add hazard ratio points
    geom_point(aes(color = color_group), size = 3, shape = 16) +
    
    # Add text annotations for HR values
    geom_text(aes(label = sprintf("%.2f", hr), x = upper_ci + 0.1), 
              size = 3, hjust = 0) +
    
    # Color scale
    scale_color_manual(values = c("Highly Significant" = "#D32F2F",
                                 "Significant" = "#FF5722", 
                                 "Marginally Significant" = "#FF9800",
                                 "Not Significant" = "#757575")) +
    
    # Scale and labels
    scale_x_continuous(breaks = seq(0, max(forest_data$upper_ci, na.rm = TRUE) + 0.5, 0.5),
                      limits = c(0, max(forest_data$upper_ci, na.rm = TRUE) + 0.5)) +
    
    # Labels
    labs(title = title,
         subtitle = subtitle,
         x = "Hazard Ratio (95% CI)",
         y = NULL,
         color = "Significance Level",
         caption = paste("Model fit: Concordance =", 
                        round(concordance(cox_model)$concordance, 3),
                        "| Likelihood ratio test p =",
                        sprintf("%.3f", anova(cox_model)[2,4]))) +
    
    # Theme
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50"),
      plot.caption = element_text(size = 10, color = "gray60", hjust = 0.5),
      axis.text.y = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, size = 0.5)
    )
  
  return(p1)
}

# Generate multiple forest plot versions
cat("  Creating basic improved forest plot...\n")
forest_plot_basic <- create_improved_forest_plot(multi_cox_model_adjusted, survival_analysis_data)

cat("  Creating comprehensive forest plot with diagnostics...\n")
forest_plot_comprehensive <- create_comprehensive_forest_plot(multi_cox_model_adjusted, survival_analysis_data)

# Save the improved plots
cat("  Saving improved forest plots...\n")
ggsave(file.path(base_path, "Multivariate_Forest_Plot_ADJUSTED_BASIC.png"), 
       plot = forest_plot_basic, width = 12, height = 10, dpi = 300)

ggsave(file.path(base_path, "Multivariate_Forest_Plot_ADJUSTED_COMPREHENSIVE.png"), 
       plot = forest_plot_comprehensive, width = 14, height = 12, dpi = 300)

# Also save the original ggforest for comparison
cat("  Creating original ggforest plot for comparison...\n")
forest_plot_original <- ggforest(multi_cox_model_adjusted, data = survival_analysis_data)
ggsave(file.path(base_path, "Multivariate_Forest_Plot_ADJUSTED_ORIGINAL.png"), 
       plot = forest_plot_original, width = 12, height = 10, dpi = 300)

cat("\n--- Prognostic analysis workflow complete. ---\n")
cat("Univariate plots and summary (with FDR) are saved.\n")
cat("ADJUSTED multivariate results are saved with multiple forest plot versions:\n")
cat("  - Original ggforest plot\n")
cat("  - Basic improved forest plot\n")
cat("  - Comprehensive forest plot with diagnostics\n")