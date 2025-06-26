# .Rprofile

# Set static variables
SCRIPT_DIR <- "sc_crohns/clustering"
SEURAT_OBJECT_PATH <- "project-area/data/crohns_scrnaseq/clustering_output"
MATRIX_DIR <- "project-area/data/crohns_scrnaseq/crohns_samples"

# Set clustering parameters
disps <- c(0.25, 0.5, 0.75, 1)
n_features <- c(500, 1000, 2000, 3000)
res <- c(0.2, 0.5, 0.8, 1.0, 1.2)
neighbors <- c(10, 20, 30, 50)

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load required packages
pkgs <- c("Seurat", "tidyverse", "Matrix", "hdf5r", "future", "future.apply", 
          "parallel", "cluster", "parallelDist", "factoextra", "cowplot", 
          "patchwork", "DoubletFinder")

for (pkg in pkgs) {
  if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Set parallelisation
set.seed(1234)
options(future.globals.maxSize = 60 * 1024^3, future.seed = TRUE)
max_workers <- parallelly::availableCores()
cores <- min(max_workers, detectCores() - 1)
plan("multisession", workers = cores)
