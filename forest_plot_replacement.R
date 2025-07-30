# 5.6: Visualize the results with an improved forest plot
cat("Generating improved forest plot for the ADJUSTED multivariate analysis...\n")

# Load required libraries
library(ggplot2)
library(dplyr)

# Create improved forest plot function
create_improved_forest_plot <- function(cox_model, data, 
                                       title = "Multivariate Cox Regression Analysis",
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
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
    
    # Add confidence intervals
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
  
  return(p)
}

# Generate the improved forest plot
forest_plot_adjusted <- create_improved_forest_plot(multi_cox_model_adjusted, survival_analysis_data)

# Save the improved plot
ggsave(file.path(base_path, "Multivariate_Forest_Plot_ADJUSTED_IMPROVED.png"), 
       plot = forest_plot_adjusted, width = 14, height = 12, dpi = 300)