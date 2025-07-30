# Improved Forest Plot for Multivariate Survival Analysis
# This script enhances the basic ggforest() with better visualization and customization

library(ggplot2)
library(survminer)
library(survival)
library(dplyr)
library(forestplot)
library(broom)

# Function to create an enhanced forest plot
create_enhanced_forest_plot <- function(cox_model, data, title = "Multivariate Cox Regression Analysis", 
                                      filename = NULL, width = 14, height = 10) {
  
  cat("Creating enhanced forest plot...\n")
  
  # Extract model results
  model_summary <- summary(cox_model)
  model_tidy <- broom::tidy(cox_model, conf.int = TRUE, exponentiate = TRUE)
  
  # Prepare data for forest plot
  forest_data <- model_tidy %>%
    mutate(
      variable = term,
      HR = estimate,
      CI_lower = conf.low,
      CI_upper = conf.high,
      p_value = p.value,
      significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**", 
        p_value < 0.05 ~ "*",
        p_value < 0.1 ~ ".",
        TRUE ~ ""
      ),
      # Format HR and CI for display
      HR_text = sprintf("%.2f", HR),
      CI_text = sprintf("(%.2f-%.2f)", CI_lower, CI_upper),
      p_text = case_when(
        p_value < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", p_value)
      ),
      # Create variable labels (clean up variable names)
      variable_clean = case_when(
        # Cell type variables (assuming they start with specific patterns)
        grepl("^[A-Z].*cell|^[A-Z].*Cell", variable) ~ gsub("\\.", " ", variable),
        variable == "age_at_index" ~ "Age at diagnosis",
        variable == "genderMale" ~ "Gender (Male vs Female)",
        variable == "tumor_stage_simpleStage.II" ~ "Stage II (vs Stage I)",
        variable == "tumor_stage_simpleStage.III" ~ "Stage III (vs Stage I)", 
        variable == "tumor_stage_simpleStage.IV" ~ "Stage IV (vs Stage I)",
        TRUE ~ gsub("\\.", " ", variable)
      ),
      # Color coding based on significance and effect direction
      color_group = case_when(
        p_value < 0.05 & HR > 1 ~ "Significant Risk",
        p_value < 0.05 & HR < 1 ~ "Significant Protective", 
        p_value >= 0.05 ~ "Non-significant",
        TRUE ~ "Non-significant"
      )
    ) %>%
    arrange(desc(HR))
  
  # Method 1: Enhanced ggforest plot
  enhanced_ggforest <- ggforest(
    cox_model, 
    data = data,
    main = title,
    cpositions = c(0.02, 0.22, 0.4),
    fontsize = 1.0,
    refLabel = "Reference",
    noDigits = 3
  ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "gray90", size = 0.5),
      panel.grid.minor = element_blank(),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    labs(
      subtitle = paste("n =", nrow(data), "patients"),
      caption = paste("Model: Cox Proportional Hazards\n",
                     "Concordance index:", round(model_summary$concordance[1], 3))
    )
  
  # Method 2: Custom ggplot forest plot with more control
  custom_forest <- ggplot(forest_data, aes(y = reorder(variable_clean, HR))) +
    # Confidence intervals
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = color_group), 
                   height = 0.3, size = 1.2, alpha = 0.8) +
    # Point estimates
    geom_point(aes(x = HR, color = color_group, size = abs(log(HR))), 
               alpha = 0.9) +
    # Reference line at HR = 1
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 0.8) +
    # Custom color scheme
    scale_color_manual(
      name = "Association",
      values = c(
        "Significant Risk" = "#d73027",
        "Significant Protective" = "#1a9850", 
        "Non-significant" = "#737373"
      )
    ) +
    scale_size_continuous(range = c(2, 5), guide = "none") +
    scale_x_log10(
      breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8),
      labels = c("0.1", "0.25", "0.5", "1", "2", "4", "8")
    ) +
    labs(
      title = title,
      subtitle = paste("n =", nrow(data), "patients |", 
                      "C-index:", round(model_summary$concordance[1], 3)),
      x = "Hazard Ratio (95% CI)",
      y = "",
      caption = "Dashed line represents HR = 1 (no effect)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text.y = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "gray90", size = 0.3),
      panel.grid.major.x = element_line(color = "gray95", size = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(20, 20, 20, 20),
      strip.text = element_text(face = "bold")
    )
  
  # Add text annotations for HR values and p-values
  custom_forest_with_text <- custom_forest +
    # HR values
    geom_text(aes(x = max(forest_data$CI_upper) * 1.1, 
                  label = paste0(HR_text, " ", CI_text)), 
              size = 3, hjust = 0, fontface = "bold") +
    # P-values
    geom_text(aes(x = max(forest_data$CI_upper) * 1.4, 
                  label = paste0("p=", p_text, " ", significance)), 
              size = 3, hjust = 0) +
    # Extend x-axis to accommodate text
    coord_cartesian(xlim = c(min(forest_data$CI_lower) * 0.8, 
                            max(forest_data$CI_upper) * 1.6))
  
  # Method 3: Professional forestplot using forestplot package
  # Prepare data for forestplot package
  forestplot_data <- forest_data %>%
    select(variable_clean, HR_text, CI_text, p_text) %>%
    add_row(variable_clean = "Variable", HR_text = "HR", CI_text = "95% CI", p_text = "P-value", .before = 1)
  
  # Create matrix for forestplot
  tabletext <- as.matrix(forestplot_data)
  hr_values <- c(NA, forest_data$HR)
  ci_lower <- c(NA, forest_data$CI_lower)
  ci_upper <- c(NA, forest_data$CI_upper)
  
  # Professional forest plot
  professional_forest <- function() {
    forestplot(
      tabletext,
      mean = hr_values,
      lower = ci_lower,
      upper = ci_upper,
      title = title,
      xlab = "Hazard Ratio (95% CI)",
      hrzl_lines = list("2" = gpar(lty = 1, lwd = 2)),
      txt_gp = fpTxtGp(
        label = gpar(cex = 1.1),
        ticks = gpar(cex = 1.0),
        xlab = gpar(cex = 1.2, fontface = "bold"),
        title = gpar(cex = 1.4, fontface = "bold")
      ),
      col = fpColors(
        box = c("red", "darkgreen")[ifelse(forest_data$p_value < 0.05, 
                                          ifelse(forest_data$HR > 1, 1, 2), 1)],
        line = "darkblue",
        summary = "royalblue"
      ),
      boxsize = 0.3,
      colgap = unit(6, "mm"),
      lineheight = unit(8, "mm"),
      cex = 1.1,
      graphwidth = unit(4, "inches")
    )
  }
  
  # Save plots
  if (!is.null(filename)) {
    # Save enhanced ggforest
    ggsave(paste0(gsub("\\.png$", "", filename), "_enhanced_ggforest.png"), 
           plot = enhanced_ggforest, width = width, height = height, dpi = 300)
    
    # Save custom forest plot
    ggsave(paste0(gsub("\\.png$", "", filename), "_custom_forest.png"), 
           plot = custom_forest_with_text, width = width, height = height, dpi = 300)
    
    # Save professional forest plot
    png(paste0(gsub("\\.png$", "", filename), "_professional_forest.png"), 
        width = width, height = height, units = "in", res = 300)
    professional_forest()
    dev.off()
    
    cat("Forest plots saved with enhanced formatting!\n")
  }
  
  # Return list of plots
  return(list(
    enhanced_ggforest = enhanced_ggforest,
    custom_forest = custom_forest_with_text,
    forest_data = forest_data,
    professional_forest_func = professional_forest
  ))
}

# Function to create a summary table alongside the forest plot
create_forest_summary_table <- function(cox_model, filename = NULL) {
  
  model_tidy <- broom::tidy(cox_model, conf.int = TRUE, exponentiate = TRUE)
  
  summary_table <- model_tidy %>%
    mutate(
      Variable = case_when(
        grepl("^[A-Z].*cell|^[A-Z].*Cell", term) ~ gsub("\\.", " ", term),
        term == "age_at_index" ~ "Age at diagnosis",
        term == "genderMale" ~ "Gender (Male vs Female)",
        term == "tumor_stage_simpleStage.II" ~ "Stage II (vs Stage I)",
        term == "tumor_stage_simpleStage.III" ~ "Stage III (vs Stage I)", 
        term == "tumor_stage_simpleStage.IV" ~ "Stage IV (vs Stage I)",
        TRUE ~ gsub("\\.", " ", term)
      ),
      `Hazard Ratio` = sprintf("%.3f", estimate),
      `95% CI` = sprintf("(%.3f - %.3f)", conf.low, conf.high),
      `P-value` = case_when(
        p.value < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", p.value)
      ),
      Significance = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**", 
        p.value < 0.05 ~ "*",
        p.value < 0.1 ~ ".",
        TRUE ~ ""
      )
    ) %>%
    select(Variable, `Hazard Ratio`, `95% CI`, `P-value`, Significance)
  
  if (!is.null(filename)) {
    write.csv(summary_table, paste0(gsub("\\.png$", "", filename), "_summary_table.csv"), 
              row.names = FALSE)
    cat("Summary table saved!\n")
  }
  
  return(summary_table)
}

# Enhanced version of your original code
enhanced_forest_plot_workflow <- function(multi_cox_model_adjusted, survival_analysis_data, base_path) {
  
  cat("Generating ENHANCED forest plot for the ADJUSTED multivariate analysis...\n")
  
  # Create enhanced forest plots
  forest_plots <- create_enhanced_forest_plot(
    cox_model = multi_cox_model_adjusted,
    data = survival_analysis_data,
    title = "Multivariate Cox Regression Analysis\n(Adjusted for Clinical Covariates)",
    filename = file.path(base_path, "Multivariate_Forest_Plot_ADJUSTED.png"),
    width = 14,
    height = 10
  )
  
  # Create summary table
  summary_table <- create_forest_summary_table(
    cox_model = multi_cox_model_adjusted,
    filename = file.path(base_path, "Multivariate_Forest_Plot_ADJUSTED.png")
  )
  
  # Print summary table
  cat("\n--- Forest Plot Summary Table ---\n")
  print(summary_table)
  
  # Display model performance metrics
  model_summary <- summary(multi_cox_model_adjusted)
  cat("\n--- Model Performance ---\n")
  cat("Concordance Index:", round(model_summary$concordance[1], 3), "\n")
  cat("R-squared:", round(model_summary$rsq["rsq"], 3), "\n")
  cat("Likelihood ratio test p-value:", 
      format.pval(model_summary$logtest["pvalue"], digits = 3), "\n")
  
  return(list(
    plots = forest_plots,
    summary_table = summary_table,
    model_metrics = list(
      concordance = model_summary$concordance[1],
      rsquared = model_summary$rsq["rsq"],
      lr_pvalue = model_summary$logtest["pvalue"]
    )
  ))
}

cat("Enhanced forest plot functions loaded successfully!\n")
cat("Use enhanced_forest_plot_workflow() to replace your existing forest plot code.\n")