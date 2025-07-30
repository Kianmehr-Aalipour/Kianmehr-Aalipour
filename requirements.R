# R Package Requirements for Univariate Survival Analysis
# Install these packages before running the analysis

# CRAN packages
required_packages <- c(
  "survival",      # Core survival analysis functions
  "survminer",     # Enhanced survival plots and statistics
  "ggplot2",       # Plotting system
  "dplyr"          # Data manipulation
)

# Bioconductor packages
bioc_packages <- c(
  "SummarizedExperiment"  # For handling TCGA data structures
)

# Function to install missing packages
install_if_missing <- function(packages, bioc = FALSE) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    }
  }
}

# Install packages
cat("Installing CRAN packages...\n")
install_if_missing(required_packages, bioc = FALSE)

cat("Installing Bioconductor packages...\n")
install_if_missing(bioc_packages, bioc = TRUE)

cat("All required packages installed successfully!\n")