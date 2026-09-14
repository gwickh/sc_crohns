library(Seurat)
library(SeuratObject)
library(Matrix)
library(CellChat)
library(patchwork)
library(igraph)
library(future)
library(org.Hs.eg.db)
library(AnnotationDbi)

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 16 * 1024^3)

future::plan("sequential")


# Paths and settings
output_dir <- path.expand(
    "~/project-area/data/crohns_scrnaseq/10c_14n_analysis/cellchat_output"
)

input_file <- file.path(
    output_dir,
    "adata_low_conf_hc_celltype_anchor_profiles.Rds"
)

diagnosis_column <- "diagnosis"
min_cells <- 10

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
}

if (file.access(output_dir, mode = 2) != 0) {
    stop("Output directory is not writable: ", output_dir)
}


# Load Seurat object
seurat_file <- readRDS(input_file)

if (!inherits(seurat_file, "Seurat")) {
    stop(
        "Expected a Seurat object, but found: ",
        paste(class(seurat_file), collapse = ", ")
    )
}

message("Available assays: ", paste(Assays(seurat_file), collapse = ", "))


# Rename the expression assay to RNA
if (
    "originalexp" %in% Assays(seurat_file) &&
    !"RNA" %in% Assays(seurat_file)
) {
    seurat_file <- RenameAssays(
        object = seurat_file,
        originalexp = "RNA"
    )
}

if (!"RNA" %in% Assays(seurat_file)) {
    stop(
        "RNA assay not found. Available assays: ",
        paste(Assays(seurat_file), collapse = ", ")
    )
}

DefaultAssay(seurat_file) <- "RNA"

message("Default assay: ", DefaultAssay(seurat_file))
message(
    "RNA layers: ",
    paste(Layers(seurat_file[["RNA"]]), collapse = ", ")
)

seurat_file <- subset(seurat_file, subset = final_curated_labels != "low_conf")

# Extract metadata
metadata <- seurat_file[[]]

if (
    "sample_id" %in% colnames(metadata) &&
    !"samples" %in% colnames(metadata)
) {
    metadata[["samples"]] <- metadata[["sample_id"]]
}

if (!diagnosis_column %in% colnames(metadata)) {
    possible_columns <- grep(
        "diagn",
        colnames(metadata),
        value = TRUE,
        ignore.case = TRUE
    )

    stop(
        "Metadata column '", diagnosis_column, "' was not found. ",
        "Possible diagnosis columns: ",
        paste(possible_columns, collapse = ", ")
    )
}

# Prefer final_curated_labels, but fall back to curated_labels
if ("final_curated_labels" %in% colnames(metadata)) {
    curated_column <- "final_curated_labels"
} else if ("curated_labels" %in% colnames(metadata)) {
    curated_column <- "curated_labels"

    warning(
        "'final_curated_labels' was not found; using 'curated_labels'."
    )
} else {
    stop(
        "Neither 'final_curated_labels' nor 'curated_labels' ",
        "was found in the metadata."
    )
}

grouping_columns <- c(curated_column, "category")

missing_columns <- setdiff(
    c(diagnosis_column, grouping_columns),
    colnames(metadata)
)

if (length(missing_columns) > 0) {
    stop(
        "Missing metadata columns: ",
        paste(missing_columns, collapse = ", ")
    )
}

message(
    "CellChat grouping columns: ",
    paste(grouping_columns, collapse = ", ")
)


# Extract normalized expression data
if (!"data" %in% Layers(seurat_file[["RNA"]])) {
    stop(
        "The RNA assay does not contain a 'data' layer. ",
        "CellChat requires normalized expression data. Available layers: ",
        paste(Layers(seurat_file[["RNA"]]), collapse = ", ")
    )
}

data_input <- GetAssayData(
    object = seurat_file,
    assay = "RNA",
    layer = "data"
)

message(
    "Expression matrix before gene conversion: ",
    nrow(data_input), " genes x ",
    ncol(data_input), " cells"
)

if (nrow(data_input) <= 5000) {
    warning(
        "The expression matrix contains only ",
        nrow(data_input),
        " genes. Check that the full-gene H5AD was converted to this RDS."
    )
}

if (is.null(rownames(data_input))) {
    stop("The expression matrix does not have gene names.")
}

if (is.null(colnames(data_input))) {
    stop("The expression matrix does not have cell names.")
}


# Convert Ensembl identifiers to gene symbols
original_gene_names <- rownames(data_input)

# Remove Ensembl version suffixes, for example ENSG000001.12
gene_ids_without_version <- sub(
    "\\.[0-9]+$",
    "",
    original_gene_names
)

is_ensembl <- grepl(
    "^ENSG[0-9]+$",
    gene_ids_without_version
)

gene_symbols <- original_gene_names

if (any(is_ensembl)) {
    ensembl_to_map <- unique(
        gene_ids_without_version[is_ensembl]
    )

    mapped_symbols <- AnnotationDbi::mapIds(
        x = org.Hs.eg.db,
        keys = ensembl_to_map,
        keytype = "ENSEMBL",
        column = "SYMBOL",
        multiVals = "first"
    )

    gene_symbols[is_ensembl] <- unname(
        mapped_symbols[
            gene_ids_without_version[is_ensembl]
        ]
    )
}

valid_gene <- (
    !is.na(gene_symbols) &
    nzchar(gene_symbols)
)

message(
    "Dropping ",
    sum(!valid_gene),
    " features without valid gene symbols"
)

data_input <- data_input[
    valid_gene,
    ,
    drop = FALSE
]

gene_symbols <- gene_symbols[valid_gene]


# For duplicated symbols, retain the row with the highest mean expression
gene_means <- Matrix::rowMeans(data_input)

gene_order <- order(
    gene_means,
    decreasing = TRUE,
    na.last = TRUE
)

data_input <- data_input[
    gene_order,
    ,
    drop = FALSE
]

gene_symbols <- gene_symbols[gene_order]

unique_gene <- !duplicated(gene_symbols)

message(
    "Dropping ",
    sum(!unique_gene),
    " duplicated gene-symbol rows"
)

data_input <- data_input[
    unique_gene,
    ,
    drop = FALSE
]

rownames(data_input) <- gene_symbols[unique_gene]

message(
    "Expression matrix after gene conversion: ",
    nrow(data_input), " genes x ",
    ncol(data_input), " cells"
)


# Match expression columns to metadata rows
common_cells <- colnames(data_input)[
    colnames(data_input) %in% rownames(metadata)
]

if (length(common_cells) == 0) {
    stop(
        "No matching cell names were found between the expression ",
        "matrix and Seurat metadata."
    )
}

if (length(common_cells) < ncol(data_input)) {
    warning(
        ncol(data_input) - length(common_cells),
        " expression-matrix cells have no matching metadata row."
    )
}

data_input <- data_input[
    ,
    common_cells,
    drop = FALSE
]

metadata <- metadata[
    common_cells,
    ,
    drop = FALSE
]

if (!identical(colnames(data_input), rownames(metadata))) {
    stop("Expression matrix and metadata cell order do not match.")
}


# Prepare diagnosis values
metadata[[diagnosis_column]] <- trimws(
    as.character(metadata[[diagnosis_column]])
)

valid_diagnosis <- (
    !is.na(metadata[[diagnosis_column]]) &
    nzchar(metadata[[diagnosis_column]])
)

diagnosis_type <- unique(
    metadata[[diagnosis_column]][valid_diagnosis]
)

if (length(diagnosis_type) == 0) {
    stop(
        "No valid diagnosis values were found in metadata column '",
        diagnosis_column,
        "'."
    )
}

message(
    "Diagnoses: ",
    paste(diagnosis_type, collapse = ", ")
)

print(
    table(
        metadata[[diagnosis_column]],
        useNA = "ifany"
    )
)

# Prepare CellChat database
CellChatDB.use <- subsetDB(CellChatDB.human)

# Safe filename helper
safe_filename <- function(x) {
    x <- gsub(
        "[^A-Za-z0-9_.-]+",
        "_",
        as.character(x)
    )

    x <- gsub("^_+|_+$", "", x)

    if (!nzchar(x)) {
        x <- "unnamed"
    }

    x
}


# Run CellChat function
run_cellchat <- function(celltypes) {
    if (!celltypes %in% colnames(metadata)) {
        stop("Metadata column not found: ", celltypes)
    }

    results <- list()

    for (samp in diagnosis_type) {
        message("")
        message("Diagnosis: ", samp)
        message("Grouping column: ", celltypes)

        diagnosis_values <- metadata[[diagnosis_column]]
        grouping_values <- trimws(
            as.character(metadata[[celltypes]])
        )

        selected <- (
            !is.na(diagnosis_values) &
            diagnosis_values == samp &
            !is.na(grouping_values) &
            nzchar(grouping_values)
        )

        cells <- rownames(metadata)[selected]

        message("Selected cells: ", length(cells))

        if (length(cells) == 0) {
            warning(
                "No cells selected for ",
                celltypes,
                " in ",
                samp
            )

            next
        }

        data_subset <- data_input[
            ,
            cells,
            drop = FALSE
        ]

        metadata_subset <- metadata[
            cells,
            ,
            drop = FALSE
        ]

        metadata_subset[[celltypes]] <- droplevels(
            factor(
                as.character(metadata_subset[[celltypes]])
            )
        )

        group_counts <- table(
            metadata_subset[[celltypes]],
            useNA = "ifany"
        )

        print(group_counts)

        if (sum(group_counts >= min_cells) == 0) {
            warning(
                "No cell groups contain at least ",
                min_cells,
                " cells for ",
                celltypes,
                " in ",
                samp,
                ". Skipping."
            )

            next
        }

        safe_group <- safe_filename(celltypes)
        safe_diagnosis <- safe_filename(samp)

        final_file <- file.path(
            output_dir,
            paste0(
                "cellchat_object_",
                safe_group,
                "_",
                safe_diagnosis,
                ".Rds"
            )
        )

        checkpoint_file <- file.path(
            output_dir,
            paste0(
                "cellchat_checkpoint_",
                safe_group,
                "_",
                safe_diagnosis,
                ".Rds"
            )
        )

        error_file <- file.path(
            output_dir,
            paste0(
                "cellchat_error_",
                safe_group,
                "_",
                safe_diagnosis,
                ".txt"
            )
        )

        result <- tryCatch(
            {
                message("Creating CellChat object")

                cellChat <- createCellChat(
                    object = data_subset,
                    meta = metadata_subset,
                    group.by = celltypes
                )

                cellChat@DB <- CellChatDB.use

                message("Subsetting signaling data")
                cellChat <- subsetData(cellChat)

                message("Identifying overexpressed genes")
                cellChat <- identifyOverExpressedGenes(cellChat)

                message("Identifying overexpressed interactions")
                cellChat <- identifyOverExpressedInteractions(cellChat)

                message("Computing communication probabilities")
                cellChat <- computeCommunProb(
                    cellChat,
                    type = "triMean"
                )

                message(
                    "Filtering interactions with fewer than ",
                    min_cells,
                    " cells"
                )

                cellChat <- filterCommunication(
                    cellChat,
                    min.cells = min_cells
                )

                message("Computing pathway probabilities")
                cellChat <- computeCommunProbPathway(cellChat)

                message("Aggregating CellChat network")
                cellChat <- aggregateNet(cellChat)

                # Save before centrality calculations
                saveRDS(
                    cellChat,
                    checkpoint_file
                )

                message("Checkpoint saved: ", checkpoint_file)

                message("Computing pathway centrality")
                cellChat <- netAnalysis_computeCentrality(
                    cellChat,
                    slot.name = "netP"
                )

                message("Computing interaction centrality")
                cellChat <- netAnalysis_computeCentrality(
                    cellChat,
                    slot.name = "net"
                )

                message("Saving final object: ", final_file)

                saveRDS(
                    cellChat,
                    final_file
                )

                if (!file.exists(final_file)) {
                    stop("Final output file was not created.")
                }

                message(
                    "Saved successfully: ",
                    round(file.info(final_file)$size / 1024^2, 2),
                    " MB"
                )

                cellChat
            },
            error = function(e) {
                error_message <- paste0(
                    "CellChat failed for grouping column '",
                    celltypes,
                    "' and diagnosis '",
                    samp,
                    "': ",
                    conditionMessage(e)
                )

                message(error_message)

                writeLines(
                    c(
                        error_message,
                        "",
                        paste(capture.output(traceback()), collapse = "\n")
                    ),
                    con = error_file
                )

                NULL
            }
        )

        if (!is.null(result)) {
            results[[as.character(samp)]] <- result
        }
    }

    results
}


# Run CellChat
cellchat_results <- setNames(
    lapply(
        grouping_columns,
        run_cellchat
    ),
    grouping_columns
)

# Report generated files
output_files <- list.files(
    output_dir,
    pattern = "^cellchat_(object|checkpoint|error)_",
    full.names = TRUE
)

message("Generated files:")
print(output_files)