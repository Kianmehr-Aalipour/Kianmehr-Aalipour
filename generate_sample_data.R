# =============================================================================
# GENERATE SAMPLE DATA FOR TESTING SURVIVAL ANALYSIS
# =============================================================================

cat("Generating sample data for survival analysis testing...\n")

# Load required libraries
suppressMessages({
  library(SummarizedExperiment)
  library(dplyr)
})

# Set seed for reproducibility
set.seed(42)

# =============================================================================
# GENERATE SAMPLE CLINICAL DATA
# =============================================================================

# Number of samples
n_samples <- 500

# Generate sample patient IDs (TCGA format)
patient_ids <- paste0("TCGA-", 
                     sprintf("%02d", sample(1:99, n_samples, replace = TRUE)), "-",
                     sprintf("%04d", sample(1000:9999, n_samples, replace = TRUE)))

# Generate clinical data
clinical_data <- data.frame(
  patient = patient_ids,
  vital_status = sample(c("Alive", "Dead"), n_samples, replace = TRUE, prob = c(0.6, 0.4)),
  days_to_last_follow_up = sample(100:3000, n_samples, replace = TRUE),
  days_to_death = sample(50:2500, n_samples, replace = TRUE),
  age_at_diagnosis = sample(30:80, n_samples, replace = TRUE),
  gender = sample(c("Male", "Female"), n_samples, replace = TRUE),
  tumor_stage = sample(c("Stage I", "Stage II", "Stage III", "Stage IV"), 
                      n_samples, replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2)),
  stringsAsFactors = FALSE
)

# Ensure consistency: if alive, no days_to_death; if dead, no days_to_last_follow_up
clinical_data$days_to_death[clinical_data$vital_status == "Alive"] <- NA
clinical_data$days_to_last_follow_up[clinical_data$vital_status == "Dead"] <- NA

# =============================================================================
# GENERATE SAMPLE CELL TYPE SCORES
# =============================================================================

# Define cell types (matching the ones in the analysis)
cell_types <- c(
  "T-Cells", "B-Cells", "CD8+ T-Cells", "NK Cells", "Macrophages", "Monocytes", 
  "Neutrophils", "Mast Cells", "Dendritic Cells", "Fibroblasts", 
  "Endothelial Cells", "AT1 Cells", "AT2 Cells", "Plasma Cells", "Regulatory T-Cells (Tregs)"
)

# Generate sample IDs (TCGA barcode format)
sample_ids <- paste0(patient_ids, "-01A-11R-1234-07")

# Generate cell type scores (beta distribution to simulate realistic scores)
estimated_scores <- matrix(
  rbeta(n_samples * length(cell_types), shape1 = 2, shape2 = 5),
  nrow = n_samples,
  ncol = length(cell_types),
  dimnames = list(sample_ids, cell_types)
)

# Add some correlation structure to make it more realistic
# Some cell types should be correlated
estimated_scores[, "CD8+ T-Cells"] <- estimated_scores[, "T-Cells"] * 0.7 + 
                                      rbeta(n_samples, 1, 3) * 0.3

estimated_scores[, "Regulatory T-Cells (Tregs)"] <- estimated_scores[, "T-Cells"] * 0.5 + 
                                                    rbeta(n_samples, 1, 4) * 0.5

# Introduce some survival-associated signal for demonstration
# Make high T-cell scores associated with better survival
survival_effect <- estimated_scores[, "T-Cells"]
clinical_data$survival_modifier <- survival_effect

# Adjust survival times based on T-cell scores (higher scores = longer survival)
clinical_data$days_to_death[!is.na(clinical_data$days_to_death)] <- 
  clinical_data$days_to_death[!is.na(clinical_data$days_to_death)] + 
  (survival_effect[!is.na(clinical_data$days_to_death)] * 500)

clinical_data$days_to_last_follow_up[!is.na(clinical_data$days_to_last_follow_up)] <- 
  clinical_data$days_to_last_follow_up[!is.na(clinical_data$days_to_last_follow_up)] + 
  (survival_effect[!is.na(clinical_data$days_to_last_follow_up)] * 300)

# Remove the modifier column
clinical_data$survival_modifier <- NULL

# =============================================================================
# CREATE SUMMARIZEDEXPERIMENT OBJECT
# =============================================================================

# Create a simple gene expression matrix (not used in survival analysis but needed for SE object)
gene_expression <- matrix(
  rnorm(n_samples * 100, mean = 10, sd = 2),
  nrow = 100,
  ncol = n_samples,
  dimnames = list(paste0("Gene_", 1:100), sample_ids)
)

# Create SummarizedExperiment object
tcga_data <- SummarizedExperiment(
  assays = list(counts = gene_expression),
  colData = clinical_data
)

# =============================================================================
# SAVE DATA
# =============================================================================

# Save objects
save(tcga_data, file = "tcga_data.RData")
save(estimated_scores, file = "estimated_scores.RData")

# Also save as separate files for inspection
write.csv(clinical_data, "sample_clinical_data.csv", row.names = FALSE)
write.csv(estimated_scores, "sample_cell_type_scores.csv")

cat("Sample data generated successfully!\n")
cat("Files created:\n")
cat("  - tcga_data.RData (SummarizedExperiment object)\n")
cat("  - estimated_scores.RData (cell type scores matrix)\n")
cat("  - sample_clinical_data.csv (clinical data for inspection)\n")
cat("  - sample_cell_type_scores.csv (cell type scores for inspection)\n")

cat("\nData summary:\n")
cat("  - Number of samples:", n_samples, "\n")
cat("  - Number of cell types:", length(cell_types), "\n")
cat("  - Vital status distribution:\n")
print(table(clinical_data$vital_status))

cat("\nYou can now run the survival analysis with:\n")
cat("  source('step4_univariate_survival_analysis.R')\n")