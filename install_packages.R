# Install packages required by the downstream TORCH analysis scripts.
# Run this script explicitly once; analysis scripts never install packages.

cran_packages <- c(
  "DescTools",
  "dplyr",
  "ggplot2",
  "pROC",
  "scales",
  "survival",
  "survminer",
  "tidyr"
)

missing_cran <- cran_packages[!vapply(
  cran_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]

if (length(missing_cran) > 0L) {
  install.packages(missing_cran)
}

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE)
}

message("Required packages are installed.")

