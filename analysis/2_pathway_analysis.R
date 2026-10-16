suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(fgsea)
})

# Print warnings where they occur in batch-job logs.
options(warn = 1)


# Configuration
project_dir <- path.expand(
  "~/project-area/data/crohns_scrnaseq/10c_14n_analysis/pseudobulk_dea"
)

input_dir <- file.path(project_dir, "degs")
output_dir <- file.path(project_dir, "pathway_enrichment")

min_gene_set_size <- 15
max_gene_set_size <- 500
n_simple_permutations <- 10000
max_pathways_in_plot <- 30

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir)
}


# Read pseudobulk DEA files
files <- list.files(
  input_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No CSV files found in: ", input_dir)
}

required_columns <- c(
  "gene_name",
  "cell_type",
  "case",
  "reference",
  "wald_statistic",
  "log2FC",
  "q_value"
)

de <- map_dfr(files, function(file) {
  x <- read_csv(
    file,
    show_col_types = FALSE
  )

  missing_columns <- setdiff(
    required_columns,
    names(x)
  )

  if (length(missing_columns) > 0) {
    stop(
      basename(file),
      " is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  mutate(
    x,
    source_file = basename(file)
  )
})

if (!is.numeric(de$wald_statistic)) {
  stop("wald_statistic must be numeric")
}

comparisons <- de |>
  filter(
    !is.na(cell_type),
    !is.na(case),
    !is.na(reference)
  ) |>
  distinct(
    cell_type,
    case,
    reference
  )

if (nrow(comparisons) == 0) {
  stop("No valid cell-type comparisons were found")
}

# Load offline gene-set files
gene_set_files <- c(
  Hallmark = file.path(
    project_dir,
    "msigdb_2026.1_human_hallmark.csv"
  ),
  Reactome = file.path(
    project_dir,
    "msigdb_2026.1_human_reactome.csv"
  ),
  GO_BP = file.path(
    project_dir,
    "msigdb_2026.1_human_go_bp.csv"
  )
)

missing_gene_set_files <- gene_set_files[
  !file.exists(gene_set_files)
]

if (length(missing_gene_set_files) > 0) {
  stop(
    "Missing gene-set files:\n",
    paste(
      missing_gene_set_files,
      collapse = "\n"
    )
  )
}

pathway_collections <- map(
  gene_set_files,
  function(file) {
    gene_set_table <- read_csv(
      normalizePath(file, mustWork = TRUE),
      show_col_types = FALSE
    )

    missing_columns <- setdiff(
      c("gs_name", "gene_symbol"),
      names(gene_set_table)
    )

    if (length(missing_columns) > 0) {
      stop(
        basename(file),
        " is missing: ",
        paste(missing_columns, collapse = ", ")
      )
    }

    gene_set_table <- gene_set_table |>
      filter(
        !is.na(gs_name),
        gs_name != "",
        !is.na(gene_symbol),
        gene_symbol != ""
      ) |>
      distinct(
        gs_name,
        gene_symbol
      )

    pathways <- split(
      gene_set_table$gene_symbol,
      gene_set_table$gs_name
    )

    pathways[lengths(pathways) > 0]
  }
)

message(
  "Loaded gene-set collections: ",
  paste(
    paste0(
      names(pathway_collections),
      " (",
      lengths(pathway_collections),
      " pathways)"
    ),
    collapse = ", "
  )
)

# Run GSEA
run_one_gsea <- function(
  cell_type,
  case,
  reference,
  pathways,
  collection_name
) {
  x <- de |>
    filter(
      .data$cell_type == .env$cell_type,
      .data$case == .env$case,
      .data$reference == .env$reference,
      !is.na(gene_name),
      gene_name != "",
      is.finite(wald_statistic)
    )

  if (nrow(x) == 0) {
    warning(
      "No valid genes for ",
      cell_type,
      ": ",
      case,
      " vs ",
      reference
    )

    return(tibble::tibble())
  }

  # fgsea requires unique gene identifiers. For duplicate symbols,
  # retain the row with the strongest absolute Wald statistic.
  ranked_table <- x |>
    group_by(gene_name) |>
    slice_max(
      order_by = abs(wald_statistic),
      n = 1,
      with_ties = FALSE
    ) |>
    ungroup()

  duplicate_symbols_removed <- nrow(x) - nrow(ranked_table)

  ranks <- setNames(
    ranked_table$wald_statistic,
    ranked_table$gene_name
  )

  ranks <- sort(
    ranks,
    decreasing = TRUE
  )

  if (length(ranks) < min_gene_set_size) {
    warning(
      "Too few ranked genes for ",
      cell_type,
      ": ",
      case,
      " vs ",
      reference
    )

    return(tibble::tibble())
  }

  pathway_genes <- unique(
    unlist(
      pathways,
      use.names = FALSE
    )
  )

  mapped_genes <- sum(
    names(ranks) %in% pathway_genes
  )

  if (mapped_genes < min_gene_set_size) {
    warning(
      "Low gene-set mapping for ",
      collection_name,
      " in ",
      cell_type,
      ": only ",
      mapped_genes,
      " genes mapped"
    )
  }

  gsea_result <- fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = min_gene_set_size,
    maxSize = max_gene_set_size,
    eps = 0,
    nPermSimple = n_simple_permutations
  )

  if (nrow(gsea_result) == 0) {
    warning(
      "No pathways passed size filters for ",
      collection_name,
      " in ",
      cell_type
    )

    return(tibble::tibble())
  }

  tibble::as_tibble(gsea_result) |>
    mutate(
      collection = collection_name,
      cell_type = cell_type,
      case = case,
      reference = reference,
      comparison = paste(case, "vs", reference),
      n_ranked_genes = length(ranks),
      n_mapped_genes = mapped_genes,
      n_duplicate_symbols_removed = duplicate_symbols_removed,
      leading_edge = map_chr(
        leadingEdge,
        ~ paste(.x, collapse = ";")
      )
    ) |>
    select(-leadingEdge)
}

results <- imap_dfr(
  pathway_collections,
  function(pathways, collection_name) {
    message("Running ", collection_name, "...")

    pmap_dfr(
      comparisons,
      function(cell_type, case, reference) {
        run_one_gsea(
          cell_type = cell_type,
          case = case,
          reference = reference,
          pathways = pathways,
          collection_name = collection_name
        )
      }
    )
  }
)

if (nrow(results) == 0) {
  stop("GSEA returned no results")
}

# fgsea's padj is calculated within each collection × cell-type run.
results <- results |>
  rename(
    q_value_within_collection_cell_type = padj
  ) |>
  group_by(
    cell_type,
    case,
    reference
  ) |>
  mutate(
    q_value_across_collections_cell_type = p.adjust(
      pval,
      method = "BH"
    )
  ) |>
  ungroup() |>
  mutate(
    q_value_global = p.adjust(
      pval,
      method = "BH"
    )
  ) |>
  arrange(
    collection,
    cell_type,
    q_value_within_collection_cell_type
  )

# Write combined outputs
write_csv(
  results,
  file.path(
    output_dir,
    "all_collections_gsea_all_results.csv"
  )
)

run_summary <- results |>
  group_by(
    collection,
    cell_type,
    case,
    reference
  ) |>
  summarise(
    pathways_tested = n(),
    significant_pathways = sum(
      q_value_within_collection_cell_type < 0.05,
      na.rm = TRUE
    ),
    missing_NES_results = sum(
      !is.finite(NES) |
        !is.finite(q_value_within_collection_cell_type)
    ),
    n_ranked_genes = first(n_ranked_genes),
    n_mapped_genes = first(n_mapped_genes),
    .groups = "drop"
  )

write_csv(
  run_summary,
  file.path(
    output_dir,
    "gsea_run_summary.csv"
  )
)


# Write one set of outputs per collection
results_by_collection <- split(
  results,
  results$collection
)

iwalk(
  results_by_collection,
  function(collection_results, collection_name) {
    file_prefix <- tolower(
      gsub(
        "[^A-Za-z0-9]+",
        "_",
        collection_name
      )
    )

    # All results
    write_csv(
      collection_results,
      file.path(
        output_dir,
        paste0(
          file_prefix,
          "_gsea_all_results.csv"
        )
      )
    )

    # Significant results using the collection-specific FDR
    significant_results <- collection_results |>
      filter(
        is.finite(NES),
        is.finite(
          q_value_within_collection_cell_type
        ),
        q_value_within_collection_cell_type < 0.05
      )

    write_csv(
      significant_results,
      file.path(
        output_dir,
        paste0(
          file_prefix,
          "_gsea_significant.csv"
        )
      )
    )

    # NES matrix
    nes_matrix <- collection_results |>
      mutate(
        population_contrast = paste0(
          cell_type,
          " (",
          case,
          " vs ",
          reference,
          ")"
        )
      ) |>
      select(
        pathway,
        population_contrast,
        NES
      ) |>
      distinct() |>
      pivot_wider(
        names_from = population_contrast,
        values_from = NES
      )

    write_csv(
      nes_matrix,
      file.path(
        output_dir,
        paste0(
          file_prefix,
          "_NES_matrix.csv"
        )
      )
    )

    # Select the 30 strongest pathways that are significant
    # in at least one cell population.
    pathway_order <- collection_results |>
      filter(
        is.finite(NES),
        is.finite(
          q_value_within_collection_cell_type
        )
      ) |>
      group_by(pathway) |>
      summarise(
        minimum_q_value = min(
          q_value_within_collection_cell_type
        ),
        significant_somewhere = any(
          q_value_within_collection_cell_type < 0.05
        ),
        .groups = "drop"
      ) |>
      filter(significant_somewhere) |>
      arrange(minimum_q_value) |>
      slice_head(n = max_pathways_in_plot)

    if (nrow(pathway_order) > 0) {
      plot_data <- collection_results |>
        semi_join(
          pathway_order,
          by = "pathway"
        ) |>
        filter(
          is.finite(NES),
          is.finite(
            q_value_within_collection_cell_type
          )
        ) |>
        mutate(
          pathway = factor(
            pathway,
            levels = rev(
              pathway_order$pathway
            )
          )
        )

      pathway_plot <- ggplot(
        plot_data,
        aes(
          x = cell_type,
          y = pathway,
          colour = NES,
          size = -log10(
            pmax(
              q_value_within_collection_cell_type,
              1e-300
            )
          )
        )
      ) +
        geom_point() +
        scale_colour_gradient2(
          low = "#2166AC",
          mid = "white",
          high = "#B2182B",
          midpoint = 0
        ) +
        scale_size_continuous(
          range = c(2, 8)
        ) +
        labs(
          title = paste(
            collection_name,
            "pathway enrichment"
          ),
          subtitle = paste(
            "Top",
            max_pathways_in_plot,
            "pathways significant in at least one population"
          ),
          x = "Cell population",
          y = NULL,
          colour = "NES",
          size = "-log10(FDR)"
        ) +
        theme_bw(base_size = 10) +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
          )
        )

      ggsave(
        file.path(
          output_dir,
          paste0(
            file_prefix,
            "_gsea_dotplot.png"
          )
        ),
        pathway_plot,
        width = 11,
        height = max(
          6,
          0.25 * nrow(pathway_order) + 2
        ),
        dpi = 300
      )
    } else {
      message(
        "No significant ",
        collection_name,
        " pathways to plot"
      )
    }
  }
)