# .Rprofile

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load required packages
pkgs <- c("Seurat", "tidyverse", "Matrix", "hdf5r")

for (pkg in pkgs) {
  if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}
