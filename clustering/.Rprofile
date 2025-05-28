# .Rprofile

# Set static variables
SCRIPT_DIR <- "sc_crohns/clustering"
SEURAT_OBJECT_LOC <- "project-area/data/crohns_scrnaseq/clustering_output"
MATRIX_DIR <- "project-area/data/crohns_scrnaseq/crohns_samples"
OUTDIR <- SEURAT_OBJECT_LOC

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load required packages
pkgs <- c("Seurat", "tidyverse", "Matrix", "hdf5r", "future", "future.apply", 
          "parallel", "cluster", "parallelDist", "factoextra", "cowplot")

for (pkg in pkgs) {
  if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Set parallelisation
options(future.globals.maxSize = 60 * 1024^3)
max_workers <- parallelly::availableCores()
cores <- min(max_workers, detectCores() - 1)
plan("multisession", workers = cores)

source(file.path(SCRIPT_DIR, "1_clustering_functions.R"))
