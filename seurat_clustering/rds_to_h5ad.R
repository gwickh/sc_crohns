# Load required libraries
library(Seurat)
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Missing package: Seurat", call. = FALSE)
}

library(SeuratDisk)
if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  stop("Missing package: SeuratDisk", call. = FALSE)
}

library(hdf5r)
if (!requireNamespace("hdf5r", quietly = TRUE)) {
  stop("Missing package: hdf5r", call. = FALSE)
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: rds_to_h5ad.R <input.rds> [output.h5ad]\n", call. = FALSE)
}

infile  <- args[1]
outfile <- if (length(args) >= 2) args[2] else sub("\\.rds$", ".h5ad", infile, ignore.case = TRUE)

# Read Seurat object from RDS file
cat("Reading Seurat object from:", infile, "\n")
obj <- readRDS(infile)
if (!inherits(obj, "Seurat")) {
  stop(sprintf(
    "RDS is not a Seurat object",
    paste(class(obj), collapse = ", ")
  ), call. = FALSE)
}

obj0 <- obj  # For summary

# Prepare Seurat object for conversion
obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])    # Join layers into assay to avoid issues with conversion

counts <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts")
obj[["RNA"]] <- Seurat::CreateAssayObject(counts = counts)  # Recreate RNA assay with counts only

# Convert to h5seurat then h5ad
tmp_h5s <- sub("\\.h5ad$", ".h5Seurat", outfile, ignore.case = TRUE)
cat("Writing temporary h5Seurat file:", tmp_h5s, "\n")
SeuratDisk::SaveH5Seurat(obj, filename = tmp_h5s, overwrite = TRUE)

cat("Converting to h5ad file \n")
SeuratDisk::Convert(tmp_h5s, dest = "h5ad", overwrite = TRUE)

cat("Wrote:", outfile, "\n")

# Print summary
cat("\n=== Seurat ===\n")
for (a in Assays(obj0)) {
  ass <- obj0[[a]]
  cat("assay:", a, "\n")
  lay <- tryCatch(SeuratObject::Layers(ass), error = function(e) NULL)
  if (!is.null(lay)) {
    cat("  layers:\n")
    cat(paste0("    - ", lay), sep = "\n")
    cat("\n")
  } else if (inherits(ass, "Assay")) {
    cat("  slots:\n")
    cat("    - counts\n    - data\n    - scale.data\n")
  }
}

cat("\n=== H5AD ===\n")
f <- H5File$new(outfile, "r")
on.exit(try(f$close_all(), silent=TRUE), add=TRUE)

cat("top:\n")
cat(paste0("  - ", names(f)), sep="\n")
cat("\n")

for (gname in intersect(c("obs","var","obsm","varm","layers","uns","raw"), names(f))) {
  g <- f[[gname]]
  cat(gname, ":\n", sep="")
  if (inherits(g, "H5Group")) {
    cat(paste0("  - ", names(g)), sep="\n")
  } else {
    cat("  - (dataset)\n")
  }
  cat("\n")
}