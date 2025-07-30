# =============================================================================
# INSTALL REQUIRED PACKAGES FOR SURVIVAL ANALYSIS
# =============================================================================

cat("Installing required R packages for survival analysis...\n")

# List of required packages
required_packages <- c(
  "survival",
  "survminer", 
  "ggplot2",
  "dplyr"
)

# Bioconductor packages
bioc_packages <- c(
  "SummarizedExperiment"
)

# Function to install CRAN packages
install_cran_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing", pkg, "from CRAN...\n")
      install.packages(pkg, dependencies = TRUE, repos = "https://cran.r-project.org/")
    } else {
      cat(pkg, "is already installed.\n")
    }
  }
}

# Function to install Bioconductor packages
install_bioc_packages <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager", repos = "https://cran.r-project.org/")
  }
  
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing", pkg, "from Bioconductor...\n")
      BiocManager::install(pkg, update = FALSE)
    } else {
      cat(pkg, "is already installed.\n")
    }
  }
}

# Install packages
cat("\n--- Installing CRAN packages ---\n")
install_cran_packages(required_packages)

cat("\n--- Installing Bioconductor packages ---\n")
install_bioc_packages(bioc_packages)

cat("\n--- Verifying installations ---\n")
all_packages <- c(required_packages, bioc_packages)
failed_packages <- c()

for (pkg in all_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "loaded successfully\n")
  } else {
    cat("✗", pkg, "failed to load\n")
    failed_packages <- c(failed_packages, pkg)
  }
}

if (length(failed_packages) == 0) {
  cat("\n🎉 All packages installed and loaded successfully!\n")
} else {
  cat("\n⚠️  The following packages failed to install:\n")
  cat(paste(failed_packages, collapse = ", "), "\n")
  cat("Please install them manually.\n")
}

cat("\nPackage installation complete.\n")