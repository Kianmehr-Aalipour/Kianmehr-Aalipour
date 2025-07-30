# REPLACE THIS SECTION:
# forest_plot_adjusted <- ggforest(multi_cox_model_adjusted, data = survival_analysis_data)
# ggsave(file.path(base_path, "Multivariate_Forest_Plot_ADJUSTED.png"), plot = forest_plot_adjusted, width = 12, height = 10)

# WITH THIS ENHANCED VERSION:

# Load additional required libraries
if (!require(broom)) install.packages("broom")
if (!require(forestplot)) install.packages("forestplot")
library(broom)
library(forestplot)

# Create enhanced forest plot data
model_tidy <- broom::tidy(multi_cox_model_adjusted, conf.int = TRUE, exponentiate = TRUE)

forest_data <- model_tidy %>%
  mutate(
    variable_clean = case_when(
      grepl("^[A-Z].*cell|^[A-Z].*Cell", term) ~ gsub("\\.", " ", term),
      term == "age_at_index" ~ "Age at diagnosis",
      term == "genderMale" ~ "Gender (Male vs Female)",
      term == "tumor_stage_simpleStage.II" ~ "Stage II (vs Stage I)",
      term == "tumor_stage_simpleStage.III" ~ "Stage III (vs Stage I)", 
      term == "tumor_stage_simpleStage.IV" ~ "Stage IV (vs Stage I)",
      TRUE ~ gsub("\\.", " ", term)
    ),
    color_group = case_when(
      p.value < 0.05 & estimate > 1 ~ "Significant Risk",
      p.value < 0.05 & estimate < 1 ~ "Significant Protective", 
      TRUE ~ "Non-significant"
    ),
    HR_text = sprintf("%.2f", estimate),
    CI_text = sprintf("(%.2f-%.2f)", conf.low, conf.high),
    p_text = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  arrange(desc(estimate))

# Method 1: Enhanced ggforest (improved version of original)
model_summary <- summary(multi_cox_model_adjusted)
enhanced_ggforest <- ggforest(
  multi_cox_model_adjusted, 
  data = survival_analysis_data,
  main = "Multivariate Cox Regression Analysis\n(Adjusted for Clinical Covariates)",
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
    panel.grid.major.y = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(
    subtitle = paste("n =", nrow(survival_analysis_data), "patients | C-index =", 
                    round(model_summary$concordance[1], 3)),
    caption = "Model: Cox Proportional Hazards Regression"
  )

# Method 2: Custom forest plot with better visual design
custom_forest <- ggplot(forest_data, aes(y = reorder(variable_clean, estimate))) +
  # Confidence intervals
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = color_group), 
                 height = 0.3, size = 1.2, alpha = 0.8) +
  # Point estimates
  geom_point(aes(x = estimate, color = color_group, size = abs(log(estimate))), 
             alpha = 0.9) +
  # Reference line at HR = 1
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 0.8) +
  # Custom colors
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
  # Add HR and p-value text
  geom_text(aes(x = max(forest_data$conf.high) * 1.1, 
                label = paste0(HR_text, " ", CI_text)), 
            size = 3, hjust = 0, fontface = "bold") +
  geom_text(aes(x = max(forest_data$conf.high) * 1.4, 
                label = paste0("p=", p_text)), 
            size = 3, hjust = 0) +
  coord_cartesian(xlim = c(min(forest_data$conf.low) * 0.8, 
                          max(forest_data$conf.high) * 1.6)) +
  labs(
    title = "Multivariate Cox Regression Analysis",
    subtitle = paste("Adjusted for Clinical Covariates | n =", nrow(survival_analysis_data), 
                    "patients | C-index =", round(model_summary$concordance[1], 3)),
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
    plot.margin = margin(20, 20, 20, 20)
  )

# Save both versions
ggsave(file.path(base_path, "Multivariate_Forest_Plot_ADJUSTED_enhanced.png"), 
       plot = enhanced_ggforest, width = 14, height = 10, dpi = 300)

ggsave(file.path(base_path, "Multivariate_Forest_Plot_ADJUSTED_custom.png"), 
       plot = custom_forest, width = 14, height = 10, dpi = 300)

# Create and save summary table
summary_table <- forest_data %>%
  select(variable_clean, HR_text, CI_text, p_text) %>%
  rename(
    Variable = variable_clean,
    `Hazard Ratio` = HR_text,
    `95% CI` = CI_text,
    `P-value` = p_text
  )

write.csv(summary_table, file.path(base_path, "Forest_Plot_Summary_Table.csv"), row.names = FALSE)

# Print summary
cat("Enhanced forest plots saved:\n")
cat("  - Enhanced ggforest version\n")
cat("  - Custom designed version with HR/CI values\n")
cat("  - Summary table (CSV)\n")
print(summary_table)