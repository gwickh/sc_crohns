#!/usr/bin/env Rscript

# Visualise adaptively summarised Reactome pathway GSEA results.
#
# Normal batch use:
#   Rscript 7_Reactome_visualisation.R \
#     <rollup_batch_root> [visualisation_root] \
#     [max_overview_terms] [max_graph_terms] [graph_parent_hops]
#
# The roll-up root should contain scope_run_manifest.csv and one result
# directory per scope. The script renders the complete existing figure suite
# for every completed scope, uses one symmetric NES colour range throughout,
# and adds cross-scope comparison tables and figures. Arguments 6-9 are
# reserved for internal worker calls made by the batch driver.

options(stringsAsFactors = FALSE, warn = 1)
set.seed(14)

# Offline Reactome hierarchy support, embedded in scripts 6/7/8.
# Official format: ReactomePathways.txt = stable ID, name, species;
# ReactomePathwaysRelation.txt = parent stable ID, child stable ID.
reactome_normalise_name <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT", sub = "")
  x <- gsub("^REACTOME_", "", toupper(x))
  paste0("REACTOME_", gsub("^_+|_+$", "", gsub("[^A-Z0-9]+", "_", x)))
}

reactome_extract_id <- function(x) {
  x <- as.character(x)
  result <- rep(NA_character_, length(x))
  hit <- !is.na(x) & grepl("R-HSA-[0-9]+", x)
  result[hit] <- sub(".*?(R-HSA-[0-9]+).*", "\\1", x[hit], perl = TRUE)
  result
}

load_reactome_hierarchy <- function(resources_dir = NA_character_, snapshot_file = NA_character_) {
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Install the R package igraph.")
  snapshot <- !is.na(snapshot_file) && nzchar(snapshot_file)
  if (snapshot) {
    if (!file.exists(snapshot_file)) stop("Hierarchy snapshot not found: ", snapshot_file)
    h <- readRDS(snapshot_file)
    if (!is.list(h) || !all(c("terms", "edges", "source") %in% names(h))) {
      stop("Invalid Reactome snapshot: expected terms, edges and source.")
    }
    if (!all(c("reactome_id", "reactome_term") %in% names(h$terms))) {
      stop("Snapshot is not a Reactome hierarchy exported by these scripts.")
    }
    source_files <- if (!is.null(h$source_files)) h$source_files else character()
    signature <- if (!is.null(h$signature)) h$signature else unname(tools::md5sum(snapshot_file))
  } else {
    if (is.na(resources_dir)) stop("An offline Reactome resources directory is required.")
    paths <- file.path(path.expand(resources_dir),
      c("ReactomePathways.txt", "ReactomePathwaysRelation.txt"))
    if (any(!file.exists(paths))) stop(
      "Missing offline Reactome hierarchy file(s): ", paste(paths[!file.exists(paths)], collapse = ", "),
      "\nDownload the two files on a networked machine and transfer them to this directory."
    )
    table <- utils::read.delim(paths[[1]], header = FALSE, sep = "\t", quote = "",
      comment.char = "", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE)
    if (ncol(table) < 3L) stop("ReactomePathways.txt must have ID, name and species columns.")
    names(table)[1:3] <- c("reactome_id", "reactome_term", "species")
    terms <- unique(table[grepl("^R-HSA-[0-9]+$", table$reactome_id) &
      table$species == "Homo sapiens", c("reactome_id", "reactome_term", "species")])
    if (!nrow(terms)) stop("No human stable pathway IDs found in ReactomePathways.txt.")
    if (anyDuplicated(terms$reactome_id)) stop("Human Reactome IDs have conflicting names.")
    relation <- utils::read.delim(paths[[2]], header = FALSE, sep = "\t", quote = "",
      comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
    if (ncol(relation) != 2L) stop("ReactomePathwaysRelation.txt must have two columns: parent, child.")
    names(relation) <- c("parent", "child")
    edges <- unique(relation[relation$parent %in% terms$reactome_id &
      relation$child %in% terms$reactome_id, ])
    edges$relationship <- "has_subpathway"
    version_file <- file.path(resources_dir, "reactome_release.txt")
    release <- if (file.exists(version_file)) paste(readLines(version_file, warn = FALSE), collapse = " ") else "unspecified release"
    source_files <- stats::setNames(paths, basename(paths))
    signature <- paste(unname(tools::md5sum(paths)), collapse = ":")
    h <- list(terms = terms, edges = edges,
      source = paste("Offline Reactome human pathway hierarchy;", release))
  }
  root <- "REACTOME_ROOT"
  terms <- as.data.frame(h$terms, stringsAsFactors = FALSE)
  terms <- terms[terms$reactome_id != root, , drop = FALSE]
  if (!"species" %in% names(terms)) terms$species <- "Homo sapiens"
  terms <- terms[grepl("^R-HSA-[0-9]+$", terms$reactome_id),
    c("reactome_id", "reactome_term", "species"), drop = FALSE]
  edges <- as.data.frame(h$edges, stringsAsFactors = FALSE)
  if (!all(c("parent", "child") %in% names(edges))) stop("Hierarchy edge columns missing.")
  edges <- unique(edges[edges$parent %in% terms$reactome_id &
    edges$child %in% terms$reactome_id, c("parent", "child"), drop = FALSE])
  edges$relationship <- rep("has_subpathway", nrow(edges))
  top_ids <- setdiff(terms$reactome_id, edges$child)
  terms <- rbind(data.frame(reactome_id = root, reactome_term = "Reactome human pathways",
    species = "Homo sapiens"), terms)
  edges <- rbind(data.frame(parent = rep(root, length(top_ids)), child = top_ids,
    relationship = rep("structural_root", length(top_ids))), edges)
  edges <- unique(edges[order(edges$child, edges$parent), ])
  g <- igraph::graph_from_data_frame(edges[, c("parent", "child")], directed = TRUE,
    vertices = data.frame(name = terms$reactome_id))
  if (!igraph::is_dag(g)) stop("Reactome pathway/subpathway hierarchy contains a cycle.")
  parent_map <- split(edges$parent, edges$child)
  min_depth <- max_depth <- stats::setNames(rep(NA_integer_, nrow(terms)), terms$reactome_id)
  min_depth[[root]] <- max_depth[[root]] <- 0L
  ancestors <- stats::setNames(vector("list", nrow(terms)), terms$reactome_id)
  for (id in igraph::as_ids(igraph::topo_sort(g, mode = "out"))) {
    parents <- unique(parent_map[[id]])
    ancestors[[id]] <- unique(c(parents, unlist(ancestors[parents], use.names = FALSE)))
    if (id == root) next
    min_depth[[id]] <- min(min_depth[parents]) + 1L
    max_depth[[id]] <- max(max_depth[parents]) + 1L
  }
  real_ids <- setdiff(terms$reactome_id, root)
  counts <- table(unlist(lapply(real_ids, function(id) c(id, ancestors[[id]])), use.names = FALSE))
  descendant_count <- as.numeric(counts[terms$reactome_id])
  descendant_count[!is.finite(descendant_count)] <- 1
  # Structural IC depends on the complete pathway hierarchy, not GO annotations
  # or this study's significance counts. Root IC is exactly zero.
  ic <- stats::setNames(-log(pmin(descendant_count / length(real_ids), 1)), terms$reactome_id)
  ic[[root]] <- 0
  terms$ontology <- "Reactome"
  terms$pathway_key <- reactome_normalise_name(terms$reactome_term)
  list(schema_version = "reactome_hierarchy_v1", terms = terms, edges = edges,
    root_id = root, top_level_ids = top_ids, parent_map = parent_map, ancestors = ancestors,
    depths = data.frame(reactome_id = terms$reactome_id,
      min_depth = unname(min_depth[terms$reactome_id]),
      max_depth = unname(max_depth[terms$reactome_id])),
    ic = ic, source = h$source, source_files = source_files, signature = signature,
    ic_definition = "-ln(number of real descendant pathways including self / all real human pathways)")
}

reactome_structural_similarity <- function(hierarchy, ids) {
  ids <- as.character(ids)
  if (!all(ids %in% hierarchy$terms$reactome_id)) stop("Unknown IDs requested for structural similarity.")
  similarity <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  lineage <- lapply(ids, function(id) unique(c(id, hierarchy$ancestors[[id]])))
  ic <- hierarchy$ic
  for (i in seq_along(ids)) {
    similarity[i, i] <- 1
    if (i == length(ids)) next
    for (j in seq.int(i + 1L, length(ids))) {
      shared <- intersect(lineage[[i]], lineage[[j]])
      denominator <- unname(ic[ids[[i]]] + ic[ids[[j]]])
      value <- if (!length(shared) || !is.finite(denominator) || denominator <= 0) 0 else {
        2 * max(ic[shared], na.rm = TRUE) / denominator
      }
      similarity[i, j] <- similarity[j, i] <- max(0, min(1, value))
    }
  }
  similarity
}

args <- commandArgs(trailingOnly = TRUE)

default_input_dir <- path.expand(paste0(
  "~/project-area/data/crohns_scrnaseq/10c_14n_analysis/",
  "pseudobulk_dea/pathway_enrichment/reactome_adaptive_summary_by_compartment_uniform"
))

input_dir <- if (length(args) >= 1 && nzchar(args[[1]])) {
  path.expand(args[[1]])
} else {
  default_input_dir
}

output_dir <- if (length(args) >= 2 && nzchar(args[[2]])) {
  path.expand(args[[2]])
} else {
  file.path(input_dir, "visualisations")
}

parse_integer_arg <- function(position, default, minimum = 1L, maximum = Inf) {
  if (length(args) < position || !nzchar(args[[position]])) return(default)
  value <- suppressWarnings(as.integer(args[[position]]))
  if (is.na(value) || value < minimum || value > maximum) {
    stop(
      "Argument ", position, " must be an integer between ", minimum,
      " and ", maximum, "."
    )
  }
  value
}

max_overview_terms <- parse_integer_arg(3L, 30L, 6L, 50L)
max_graph_terms <- parse_integer_arg(4L, 18L, 3L, 30L)
graph_parent_hops <- parse_integer_arg(5L, 2L, 1L, 4L)
run_mode <- if (length(args) >= 6L && nzchar(args[[6L]])) {
  args[[6L]]
} else {
  "batch"
}
scope_id <- if (length(args) >= 7L && nzchar(args[[7L]])) {
  args[[7L]]
} else {
  "whole_atlas"
}
scope_label <- if (length(args) >= 8L && nzchar(args[[8L]])) {
  args[[8L]]
} else {
  "Whole atlas"
}
shared_nes_colour_limit <- if (length(args) >= 9L && nzchar(args[[9L]])) {
  suppressWarnings(as.numeric(args[[9L]]))
} else {
  NA_real_
}

if (!run_mode %in% c("batch", "worker")) {
  stop("Argument 6 must be either 'batch' or 'worker'.")
}
max_graph_nodes <- 120L
max_dechaining_parents <- 8L
max_descendants_per_parent <- 6L
top_terms_per_cell <- 12L
supplement_terms_per_page <- 40L

required_packages <- c(
  "readr", "dplyr", "tidyr", "purrr", "ggplot2", "stringr",
  "scales", "igraph", "ggrepel"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(stringr)
  library(scales)
  library(igraph)
  library(ggrepel)
})

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument) == 1L) {
  normalizePath(sub("^--file=", "", script_argument), mustWork = TRUE)
} else {
  stop("Could not determine this script's path for compartment batch execution.")
}

# Backward-compatible single-directory use: if the caller supplies one completed
# roll-up directory rather than a batch root, render it as a whole-atlas worker.
if (
  run_mode == "batch" &&
  file.exists(file.path(input_dir, "adaptive_reactome_final_by_cell.csv")) &&
  !any(file.exists(file.path(
    input_dir, c("scope_run_manifest.csv", "scope_manifest.csv")
  )))
) {
  run_mode <- "worker"
  scope_id <- "whole_atlas"
  scope_label <- "Whole atlas"
}

save_batch_plot_pair <- function(plot, directory, stem, width, height) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    file.path(directory, paste0(stem, ".png")), plot,
    width = width, height = height, dpi = 300, bg = "white",
    limitsize = FALSE
  )
  ggsave(
    file.path(directory, paste0(stem, ".pdf")), plot,
    width = width, height = height, bg = "white", limitsize = FALSE
  )
}

# Ignore legacy custom biological annotations when older result tables are read.
drop_custom_annotations <- function(data) {
  select(data, -matches("priority|^theme_type$|^max_granular$", ignore.case = TRUE))
}

run_compartment_visualisation_batch <- function() {
  manifest_candidates <- file.path(
    input_dir, c("scope_run_manifest.csv", "scope_manifest.csv")
  )
  manifest_path <- manifest_candidates[file.exists(manifest_candidates)][1]
  if (is.na(manifest_path)) {
    stop(
      "No scope manifest was found in ", input_dir,
      ". Expected scope_run_manifest.csv (preferred) or scope_manifest.csv."
    )
  }

  scope_manifest <- read_csv(
    manifest_path, show_col_types = FALSE, progress = FALSE
  )

  first_existing_column <- function(data, candidates, required = TRUE) {
    found <- candidates[candidates %in% names(data)]
    if (length(found) > 0L) return(found[[1L]])
    if (required) {
      stop(
        "Scope manifest is missing one of these columns: ",
        paste(candidates, collapse = ", ")
      )
    }
    NA_character_
  }

  id_column <- first_existing_column(
    scope_manifest, c("scope_id", "analysis_scope", "scope")
  )
  label_column <- first_existing_column(
    scope_manifest, c("scope_label", "analysis_scope_label", "label"),
    required = FALSE
  )
  status_column <- first_existing_column(
    scope_manifest, c("status", "run_status"), required = FALSE
  )
  output_column <- first_existing_column(
    scope_manifest,
    c("output_dir", "scope_output_dir", "results_dir"),
    required = FALSE
  )

  scopes <- scope_manifest |>
    mutate(
      scope_id = as.character(.data[[id_column]]),
      scope_label = if (!is.na(label_column)) {
        as.character(.data[[label_column]])
      } else {
        stringr::str_to_sentence(stringr::str_replace_all(scope_id, "_", " "))
      },
      rollup_status = if (!is.na(status_column)) {
        tolower(as.character(.data[[status_column]]))
      } else {
        "completed"
      },
      scope_type = if ("scope_type" %in% names(scope_manifest)) {
        as.character(.data$scope_type)
      } else {
        if_else(scope_id == "whole_atlas", "whole_atlas", "compartment")
      },
      n_requested_cell_types = if (
        "n_requested_cell_types" %in% names(scope_manifest)
      ) {
        as.integer(.data$n_requested_cell_types)
      } else {
        NA_integer_
      },
      n_matched_cell_types = if (
        "n_matched_cell_types" %in% names(scope_manifest)
      ) {
        as.integer(.data$n_matched_cell_types)
      } else {
        NA_integer_
      },
      n_missing_cell_types = if (
        "n_missing_cell_types" %in% names(scope_manifest)
      ) {
        as.integer(.data$n_missing_cell_types)
      } else {
        NA_integer_
      },
      matched_cell_types = if (
        "matched_cell_types" %in% names(scope_manifest)
      ) {
        as.character(.data$matched_cell_types)
      } else {
        NA_character_
      },
      missing_cell_types = if (
        "missing_cell_types" %in% names(scope_manifest)
      ) {
        as.character(.data$missing_cell_types)
      } else {
        NA_character_
      },
      scope_input_dir = if (!is.na(output_column)) {
        path.expand(as.character(.data[[output_column]]))
      } else {
        file.path(.env$input_dir, "scopes", scope_id)
      }
    ) |>
    mutate(
      scope_input_dir = if_else(
        dir.exists(scope_input_dir),
        scope_input_dir,
        if_else(
          dir.exists(file.path(.env$input_dir, "scopes", scope_id)),
          file.path(.env$input_dir, "scopes", scope_id),
          file.path(.env$input_dir, scope_id)
        )
      ),
      scope_output_dir = file.path(.env$output_dir, scope_id),
      has_final_by_cell = file.exists(file.path(
        scope_input_dir, "adaptive_reactome_final_by_cell.csv"
      )),
      has_final_terms = file.exists(file.path(
        scope_input_dir, "adaptive_reactome_final_significant_terms.csv"
      )),
      eligible = (
        grepl("^completed", rollup_status) |
          rollup_status %in% c("success", "succeeded")
      ) &
        has_final_by_cell & has_final_terms
    )

  if (!any(scopes$eligible)) {
    stop("No completed scope contains the two final Reactome result tables.")
  }

  eligible_scopes <- scopes |> filter(eligible)

  scoped_cell_results <- purrr::pmap_dfr(
    eligible_scopes |> select(scope_id, scope_label, scope_input_dir),
    function(scope_id, scope_label, scope_input_dir) {
      read_csv(
        file.path(scope_input_dir, "adaptive_reactome_final_by_cell.csv"),
        show_col_types = FALSE, progress = FALSE
      ) |>
        drop_custom_annotations() |>
        mutate(
          scope_id = .env$scope_id,
          scope_label = .env$scope_label
        )
    }
  )

  global_nes_limit <- scoped_cell_results |>
    filter(
      !is.na(evidence_NES), is.finite(evidence_NES),
      is.na(direction) | direction != "mixed"
    ) |>
    summarise(value = max(abs(evidence_NES), na.rm = TRUE)) |>
    pull(value)

  if (
    length(global_nes_limit) != 1L || !is.finite(global_nes_limit) ||
    global_nes_limit <= 0
  ) {
    global_nes_limit <- 1
  }

  run_records <- vector("list", nrow(scopes))
  for (index in seq_len(nrow(scopes))) {
    current <- scopes[index, ]
    if (!current$eligible) {
      run_records[[index]] <- current |>
        transmute(
          scope_id, scope_label, scope_type,
          n_requested_cell_types, n_matched_cell_types, n_missing_cell_types,
          matched_cell_types, missing_cell_types,
          scope_input_dir, scope_output_dir, rollup_status,
          visualisation_status = "skipped",
          exit_status = NA_integer_,
          note = "Missing completed final-term inputs"
        )
      next
    }

    dir.create(current$scope_output_dir, recursive = TRUE, showWarnings = FALSE)
    worker_args <- c(
      shQuote(script_path),
      shQuote(current$scope_input_dir),
      shQuote(current$scope_output_dir),
      as.character(max_overview_terms),
      as.character(max_graph_terms),
      as.character(graph_parent_hops),
      "worker",
      shQuote(current$scope_id),
      shQuote(current$scope_label),
      format(global_nes_limit, digits = 16, scientific = FALSE)
    )
    exit_status <- system2(
      file.path(R.home("bin"), "Rscript"), args = worker_args
    )

    run_records[[index]] <- current |>
      transmute(
        scope_id, scope_label, scope_type,
        n_requested_cell_types, n_matched_cell_types, n_missing_cell_types,
        matched_cell_types, missing_cell_types,
        scope_input_dir, scope_output_dir, rollup_status,
        visualisation_status = if_else(
          exit_status == 0L, "completed", "failed"
        ),
        exit_status = as.integer(exit_status),
        note = if_else(
          exit_status == 0L, "", "Worker visualization returned a non-zero status"
        )
      )
  }

  visualisation_manifest <- bind_rows(run_records)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(
    visualisation_manifest,
    file.path(output_dir, "visualisation_scope_manifest.csv")
  )

  successful_scopes <- visualisation_manifest |>
    filter(visualisation_status == "completed")
  if (nrow(successful_scopes) == 0L) {
    stop("Every eligible scope visualization worker failed.")
  }

  cross_scope_dir <- file.path(output_dir, "cross_scope")
  dir.create(cross_scope_dir, recursive = TRUE, showWarnings = FALSE)

  successful_ids <- successful_scopes$scope_id
  cross_cell <- scoped_cell_results |> filter(scope_id %in% successful_ids)
  cross_terms <- purrr::pmap_dfr(
    successful_scopes |>
      select(
        scope_id, scope_label, scope_input_dir, n_matched_cell_types
      ),
    function(
      scope_id, scope_label, scope_input_dir, n_matched_cell_types
    ) {
      data <- read_csv(
        file.path(scope_input_dir, "adaptive_reactome_final_significant_terms.csv"),
        show_col_types = FALSE, progress = FALSE
      ) |> drop_custom_annotations()
      n_scope_cells <- suppressWarnings(as.integer(n_matched_cell_types))
      if (
        length(n_scope_cells) != 1L || !is.finite(n_scope_cells) ||
        n_scope_cells < 1L
      ) {
        n_scope_cells <- cross_cell |>
          filter(.data$scope_id == .env$scope_id) |>
          summarise(value = n_distinct(cell_type)) |>
          pull(value)
      }
      data |>
        mutate(
          scope_id = .env$scope_id,
          scope_label = .env$scope_label,
          n_scope_cell_types = n_scope_cells,
          significant_cell_fraction = pmin(
            1, n_significant_cell_types / pmax(1, n_scope_cells)
          )
        )
    }
  )
  write_csv(
    cross_terms,
    file.path(cross_scope_dir, "cross_scope_final_terms_long.csv")
  )

  cross_counts <- cross_terms |>
    count(
      scope_id, scope_label, overall_direction, resolution,
      name = "n_final_terms"
    )
  write_csv(
    cross_counts,
    file.path(cross_scope_dir, "cross_scope_final_term_counts.csv")
  )

  scope_order <- successful_scopes$scope_label
  counts_plot <- cross_counts |>
    mutate(scope_label = factor(scope_label, levels = rev(scope_order))) |>
    ggplot(aes(x = scope_label, y = n_final_terms, fill = resolution)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~overall_direction, scales = "free_x") +
    scale_fill_brewer(palette = "Set2", name = "Resolution") +
    labs(
      title = "Final non-redundant Reactome terms by analysis scope",
      subtitle = "Counts are separated by enrichment direction and retained resolution",
      x = NULL, y = "Number of final Reactome terms",
      caption = paste(
        "Counts reflect each independently recomputed roll-up; compartments with",
        "fewer represented cell types should not be compared as if sampling were equal."
      )
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  save_batch_plot_pair(
    counts_plot, cross_scope_dir, "01_final_term_counts_by_scope", 12, 7
  )

  term_sets <- split(cross_terms$display_reactome_id, cross_terms$scope_id)
  scope_keys <- successful_scopes$scope_id
  jaccard_long <- tidyr::expand_grid(
    scope_id_1 = scope_keys, scope_id_2 = scope_keys
  ) |>
    rowwise() |>
    mutate(
      intersection_n = length(intersect(
        unique(term_sets[[scope_id_1]]), unique(term_sets[[scope_id_2]])
      )),
      union_n = length(union(
        unique(term_sets[[scope_id_1]]), unique(term_sets[[scope_id_2]])
      )),
      jaccard = if_else(union_n == 0L, NA_real_, intersection_n / union_n)
    ) |>
    ungroup() |>
    left_join(
      successful_scopes |> select(scope_id_1 = scope_id, label_1 = scope_label),
      by = "scope_id_1"
    ) |>
    left_join(
      successful_scopes |> select(scope_id_2 = scope_id, label_2 = scope_label),
      by = "scope_id_2"
    )
  write_csv(
    jaccard_long,
    file.path(cross_scope_dir, "cross_scope_Reactome_ID_jaccard_long.csv")
  )
  jaccard_matrix <- jaccard_long |>
    select(scope_id_1, scope_id_2, jaccard) |>
    pivot_wider(names_from = scope_id_2, values_from = jaccard)
  write_csv(
    jaccard_matrix,
    file.path(cross_scope_dir, "cross_scope_Reactome_ID_jaccard_matrix.csv")
  )

  jaccard_plot <- jaccard_long |>
    mutate(
      label_1 = factor(label_1, levels = rev(scope_order)),
      label_2 = factor(label_2, levels = scope_order)
    ) |>
    ggplot(aes(x = label_2, y = label_1, fill = jaccard)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = if_else(is.na(jaccard), "NA", sprintf("%.2f", jaccard))), size = 3) +
    scale_fill_gradient(
      low = "white", high = "#3F6C8E", limits = c(0, 1),
      na.value = "#EFEFEF", name = "Jaccard"
    ) +
    coord_equal() +
    labs(
      title = "Overlap of final Reactome identifiers across scopes",
      subtitle = "Jaccard = shared Reactome IDs divided by the union of Reactome IDs",
      x = NULL, y = NULL,
      caption = paste(
        "This measures overlap between independently selected summaries, not",
        "similarity of the underlying transcriptomes."
      )
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  save_batch_plot_pair(
    jaccard_plot, cross_scope_dir, "02_final_Reactome_ID_jaccard", 10, 8
  )

  shared_term_rank <- cross_terms |>
    group_by(display_reactome_id) |>
    summarise(
      display_term = dplyr::first(display_term),
      n_scopes = n_distinct(scope_id),
      max_fraction = max(significant_cell_fraction, na.rm = TRUE),
      best_q_value = min(best_member_q_value, na.rm = TRUE),
      max_abs_nes = max(largest_member_abs_NES, na.rm = TRUE),
      .groups = "drop"
    ) |>
    filter(n_scopes >= 2L) |>
    arrange(desc(n_scopes), desc(max_fraction), best_q_value, desc(max_abs_nes)) |>
    slice_head(n = 40L)

  shared_terms <- cross_terms |>
    semi_join(shared_term_rank, by = "display_reactome_id") |>
    left_join(
      shared_term_rank |> select(display_reactome_id, n_scopes),
      by = "display_reactome_id"
    ) |>
    mutate(
      signed_effect = case_when(
        overall_direction == "positive" ~ largest_member_abs_NES,
        overall_direction == "negative" ~ -largest_member_abs_NES,
        TRUE ~ 0
      ),
      term_label = stringr::str_wrap(
        paste0(display_term, " (", display_reactome_id, ")"), width = 45
      )
    )
  write_csv(
    shared_terms,
    file.path(cross_scope_dir, "cross_scope_shared_terms_selected.csv")
  )

  if (nrow(shared_terms) > 0L) {
    shared_term_order <- shared_term_rank |>
      mutate(
        term_label = stringr::str_wrap(
          paste0(display_term, " (", display_reactome_id, ")"), width = 45
        )
      ) |>
      pull(term_label)
    shared_plot <- shared_terms |>
      mutate(
        scope_label = factor(scope_label, levels = scope_order),
        term_label = factor(term_label, levels = rev(shared_term_order))
      ) |>
      ggplot(aes(x = scope_label, y = term_label)) +
      geom_point(
        aes(fill = signed_effect, size = significant_cell_fraction),
        shape = 21, colour = "#252525", stroke = 0.35
      ) +
      scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0, limits = c(-global_nes_limit, global_nes_limit),
        oob = scales::squish, name = "Signed max |NES|"
      ) +
      scale_size_continuous(
        range = c(2, 8), limits = c(0, 1),
        name = "Fraction of scope\ncell types significant"
      ) +
      labs(
        title = "Final Reactome terms retained in multiple analysis scopes",
        subtitle = "Up to 40 terms, prioritised by scope recurrence and within-scope breadth",
        x = NULL, y = NULL,
        caption = paste(
          "Blank combinations mean the term was not retained in that scope's final",
          "summary; they do not prove absence of biological activity."
        )
      ) +
      theme_bw(base_size = 9) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        legend.position = "right"
      )
    save_batch_plot_pair(
      shared_plot, cross_scope_dir, "03_shared_final_terms_dotplot",
      max(12, 1.1 * length(scope_order) + 7),
      max(8, 0.25 * n_distinct(shared_terms$display_reactome_id) + 3)
    )
  }

  cross_manifest <- tibble(
    artifact = c(
      "cross_scope_final_terms_long.csv",
      "cross_scope_final_term_counts.csv",
      "01_final_term_counts_by_scope.png; 01_final_term_counts_by_scope.pdf",
      "cross_scope_Reactome_ID_jaccard_long.csv; cross_scope_Reactome_ID_jaccard_matrix.csv",
      "02_final_Reactome_ID_jaccard.png; 02_final_Reactome_ID_jaccard.pdf",
      "cross_scope_shared_terms_selected.csv",
      if (nrow(shared_terms) > 0L) {
        "03_shared_final_terms_dotplot.png; 03_shared_final_terms_dotplot.pdf"
      } else {
        "03_shared_final_terms_dotplot (not created: no Reactome ID occurred in >=2 scopes)"
      }
    ),
    interpretation = c(
      "Combined final, non-redundant term list with scope labels and within-scope breadth.",
      "Final-term counts by scope, direction and retained resolution.",
      "Visual comparison of the sizes and composition of the final scope summaries.",
      "Exact final Reactome-ID intersections, unions and Jaccard similarities.",
      "Overlap of selected summaries; this is not transcriptomic similarity.",
      "The recurrent terms used in the cross-scope dot plot.",
      "Recurrence, direction, effect strength and cell-type breadth for shared final terms."
    )
  )
  write_csv(
    cross_manifest,
    file.path(cross_scope_dir, "cross_scope_figure_manifest.csv")
  )

  batch_parameters <- tibble(
    parameter = c(
      "input_root", "output_root", "scope_manifest", "n_scopes_eligible",
      "n_scopes_visualised", "shared_symmetric_NES_colour_limit",
      "max_overview_terms", "max_graph_terms", "graph_parent_hops"
    ),
    value = as.character(c(
      input_dir, output_dir, manifest_path, nrow(eligible_scopes),
      nrow(successful_scopes), global_nes_limit, max_overview_terms,
      max_graph_terms, graph_parent_hops
    ))
  )
  write_csv(
    batch_parameters,
    file.path(output_dir, "visualisation_batch_parameters.csv")
  )

  failed_n <- sum(visualisation_manifest$visualisation_status == "failed")
  message("Created scope visualisations in: ", output_dir)
  message("Successfully visualised scopes: ", nrow(successful_scopes))
  message("Failed scope workers: ", failed_n)
  if (failed_n > 0L) {
    stop(
      failed_n, " scope visualization worker(s) failed. See ",
      file.path(output_dir, "visualisation_scope_manifest.csv")
    )
  }
}

if (run_mode == "batch") {
  run_compartment_visualisation_batch()
  quit(save = "no", status = 0L)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

colour_positive <- "#B2182B"
colour_negative <- "#2166AC"
colour_mixed <- "#7B5EA7"
colour_missing <- "#EFEFEF"
colour_connector <- "#D9D9D9"
colour_text <- "#222222"

direction_colours <- c(
  positive = colour_positive,
  negative = colour_negative,
  mixed = colour_mixed
)

resolution_shapes <- c(
  granular = 21,
  direct = 22,
  ontology_rollup = 24,
  semantic_group = 23
)

theme_reactome <- theme_bw(base_size = 10, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#E7E7E7", linewidth = 0.25),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", colour = colour_text),
    plot.subtitle = element_text(colour = "#4A4A4A"),
    plot.caption = element_text(colour = "#5A5A5A", hjust = 0),
    axis.text = element_text(colour = colour_text),
    legend.position = "right"
  )

data_files <- c(
  display_mapping = "adaptive_display_mapping.csv",
  final_by_cell = "adaptive_reactome_final_by_cell.csv",
  final_terms = "adaptive_reactome_final_significant_terms.csv",
  global_selection = "adaptive_reactome_global_granular_selection.csv",
  pruning_audit = "adaptive_reactome_hierarchy_pruning_audit.csv",
  nes_matrix = "adaptive_Reactome_NES_matrix.csv",
  fdr_matrix = "adaptive_Reactome_source_FDR_matrix.csv",
  summary_long = "adaptive_reactome_summary_long.csv",
  term_metrics = "reactome_term_metrics.csv",
  reduction_counts = "reduction_counts.csv"
)

data_paths <- file.path(input_dir, data_files)
names(data_paths) <- names(data_files)

missing_files <- data_paths[!file.exists(data_paths)]
if (length(missing_files) > 0) {
  stop(
    "The following required output files were not found:\n",
    paste(" -", missing_files, collapse = "\n")
  )
}

read_output <- function(key) {
  read_csv(data_paths[[key]], show_col_types = FALSE, progress = FALSE) |>
    drop_custom_annotations()
}

display_mapping <- read_output("display_mapping")
final_by_cell <- read_output("final_by_cell")
final_terms <- read_output("final_terms")
global_selection <- read_output("global_selection")
pruning_audit <- read_output("pruning_audit")
term_metrics <- read_output("term_metrics")
reduction_counts <- read_output("reduction_counts")

validate_columns <- function(data, columns, table_name) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    stop(
      table_name, " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

validate_columns(
  final_terms,
  c(
    "display_reactome_id", "display_term", "resolution", "display_max_depth",
    "overall_direction", "n_significant_cell_types", "supporting_reactome_ids",
    "best_member_q_value", "largest_member_abs_NES"
  ),
  "adaptive_reactome_final_significant_terms.csv"
)

validate_columns(
  final_by_cell,
  c(
    "cell_type", "display_reactome_id", "display_term", "direction", "resolution",
    "evidence_source", "evidence_reactome_id", "evidence_NES", "evidence_q_value",
    "evidence_leading_edge"
  ),
  "adaptive_reactome_final_by_cell.csv"
)

validate_columns(
  display_mapping,
  c("reactome_id"),
  "adaptive_display_mapping.csv"
)

validate_columns(
  global_selection,
  c(
    "member_direction",
    "global_granular_selected"
  ),
  "adaptive_reactome_global_granular_selection.csv"
)

validate_columns(
  pruning_audit,
  c(
    "pruning_scope", "removed_reactome_id", "removed_term",
    "retained_descendant_reactome_ids", "direction", "n_significant_cell_types"
  ),
  "adaptive_reactome_hierarchy_pruning_audit.csv"
)

validate_columns(
  reduction_counts,
  c("metric", "value"),
  "reduction_counts.csv"
)

validate_columns(
  term_metrics,
  c(
    "reactome_id", "reactome_term", "max_abs_nes", "evidence_pass", "max_depth",
    "information_content", "best_q", "n_significant_cell_types",
    "keep_granular"
  ),
  "reactome_term_metrics.csv"
)

if (anyDuplicated(final_terms$display_reactome_id)) {
  stop(
    "The final significant-term table contains duplicated display_reactome_id values. ",
    "Resolve these before plotting the final list."
  )
}

if (anyDuplicated(final_by_cell[c("cell_type", "display_reactome_id")])) {
  stop(
    "The cell-specific table contains duplicated cell_type/display_reactome_id pairs."
  )
}

as_logical_flag <- function(x) {
  if (is.logical(x)) return(replace_na(x, FALSE))
  toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES", "Y")
}

collapse_non_missing <- function(x) {
  x <- unique(str_trim(as.character(x)))
  x <- x[!is.na(x) & nzchar(x) & x != "NA"]
  if (length(x) == 0) NA_character_ else paste(sort(x), collapse = "; ")
}

scale_01 <- function(x) {
  x <- as.numeric(x)
  finite <- is.finite(x)
  result <- rep(0, length(x))
  if (!any(finite)) return(result)
  limits <- range(x[finite], na.rm = TRUE)
  if (diff(limits) == 0) {
    result[finite] <- 0.5
  } else {
    result[finite] <- (x[finite] - limits[[1]]) / diff(limits)
  }
  result
}

safe_neg_log10 <- function(x, cap = 50) {
  x <- suppressWarnings(as.numeric(x))
  result <- -log10(pmax(x, 1e-300))
  pmin(result, cap)
}

wrap_term <- function(term, reactome_id, width = 42) {
  paste0(str_wrap(term, width = width), "\n", reactome_id)
}

wrap_caption <- function(...) {
  str_wrap(paste(...), width = 145)
}

final_terms_plot <- final_terms |>
  mutate(
    q_strength = safe_neg_log10(best_member_q_value),
    breadth_component = scale_01(n_significant_cell_types),
    effect_component = scale_01(largest_member_abs_NES),
    q_component = scale_01(q_strength),
    depth_component = scale_01(display_max_depth),
    presentation_score =
      0.30 * breadth_component +
      0.30 * effect_component +
      0.25 * q_component +
      0.15 * depth_component,
    display_label = wrap_term(display_term, display_reactome_id)
  )

select_overview_terms <- function(data, n_terms) {
  n_broad <- min(ceiling(n_terms / 2), nrow(data))

  selected_broad <- data |>
    arrange(
      desc(n_significant_cell_types), desc(q_strength),
      desc(largest_member_abs_NES)
    ) |>
    slice_head(n = n_broad) |>
    transmute(display_reactome_id, selection_reason = "Broad across cell populations")

  selected_so_far <- selected_broad
  n_remaining <- max(0L, n_terms - nrow(selected_so_far))

  selected_strong <- data |>
    anti_join(selected_so_far, by = "display_reactome_id") |>
    arrange(
      desc(largest_member_abs_NES), desc(q_strength),
      desc(display_max_depth)
    ) |>
    slice_head(n = n_remaining) |>
    transmute(display_reactome_id, selection_reason = "Strong or specific remaining signal")

  bind_rows(selected_so_far, selected_strong) |>
    distinct(display_reactome_id, .keep_all = TRUE) |>
    left_join(data, by = "display_reactome_id")
}

overview_terms <- select_overview_terms(final_terms_plot, max_overview_terms)

write_csv(
  overview_terms |>
    select(
      display_reactome_id, display_term, selection_reason,
      overall_direction, resolution, n_significant_cell_types,
      best_member_q_value, largest_member_abs_NES, display_max_depth,
      presentation_score
    ),
  file.path(output_dir, "visualised_terms.csv")
)

final_cell_plot <- final_by_cell |>
  semi_join(final_terms_plot, by = "display_reactome_id") |>
  left_join(
    final_terms_plot |>
      select(
        display_reactome_id, global_direction = overall_direction,
        global_resolution = resolution,
        global_n_significant_cell_types = n_significant_cell_types,
        global_presentation_score = presentation_score
      ),
    by = "display_reactome_id"
  ) |>
  mutate(
    direction = if_else(
      direction %in% c("positive", "negative", "mixed"),
      direction,
      "mixed"
    ),
    evidence_mode = if_else(
      evidence_source == "display_term" & evidence_reactome_id == display_reactome_id,
      "Direct displayed term",
      "Supporting member"
    ),
    plot_log_fdr = safe_neg_log10(evidence_q_value, cap = 20),
    signed_nes = if_else(direction == "mixed", NA_real_, evidence_NES),
    display_label = wrap_term(display_term, display_reactome_id)
  )

# Reuse one symmetric NES colour range in every dot plot and heatmap so the
# same red/blue shade always represents the same effect size across figures.
nes_colour_limit <- final_cell_plot |>
  filter(direction != "mixed", is.finite(evidence_NES)) |>
  summarise(limit = max(abs(evidence_NES), na.rm = TRUE)) |>
  pull(limit)

if (length(nes_colour_limit) != 1 || !is.finite(nes_colour_limit) ||
    nes_colour_limit <= 0) {
  nes_colour_limit <- 1
}

if (
  length(shared_nes_colour_limit) == 1L &&
  is.finite(shared_nes_colour_limit) && shared_nes_colour_limit > 0
) {
  nes_colour_limit <- shared_nes_colour_limit
}

scoped_title <- function(title) {
  paste0(scope_label, ": ", title)
}

cell_order_from_data <- function(data, term_ids) {
  cells <- sort(unique(final_cell_plot$cell_type))
  matrix_data <- data |>
    filter(display_reactome_id %in% term_ids, direction != "mixed") |>
    select(display_reactome_id, cell_type, evidence_NES) |>
    pivot_wider(names_from = cell_type, values_from = evidence_NES)

  if (nrow(matrix_data) < 2 || ncol(matrix_data) < 3) return(cells)

  ids <- matrix_data$display_reactome_id
  matrix_data <- as.data.frame(matrix_data[, setdiff(names(matrix_data), "display_reactome_id")])
  matrix_values <- as.matrix(matrix_data)
  rownames(matrix_values) <- ids
  matrix_values[!is.finite(matrix_values)] <- 0

  if (ncol(matrix_values) < 2) return(cells)

  clustered <- tryCatch(
    colnames(matrix_values)[hclust(dist(t(matrix_values)))$order],
    error = function(e) cells
  )

  c(clustered, setdiff(cells, clustered))
}

order_terms_for_plot <- function(term_ids) {
  final_terms_plot |>
    filter(display_reactome_id %in% term_ids) |>
    arrange(
      factor(overall_direction, levels = c("positive", "mixed", "negative")),
      desc(n_significant_cell_types), desc(presentation_score)
    ) |>
    pull(display_label)
}

figure_manifest <- tibble(
  scope_id = character(),
  scope_label = character(),
  artifact = character(),
  title = character(),
  type = character(),
  source_tables = character(),
  interpretation = character(),
  status = character()
)

record_artifact <- function(
  artifact, title, type, source_tables, interpretation,
  status = "created"
) {
  figure_manifest <<- bind_rows(
    figure_manifest,
    tibble(
      scope_id = scope_id,
      scope_label = scope_label,
      artifact = artifact,
      title = title,
      type = type,
      source_tables = source_tables,
      interpretation = interpretation,
      status = status
    )
  )
}

save_plot_pair <- function(
  plot, stem, title, width, height, source_tables, interpretation
) {
  png_name <- paste0(stem, ".png")
  pdf_name <- paste0(stem, ".pdf")
  full_title <- scoped_title(title)
  plot <- plot + labs(title = full_title)

  ggsave(
    file.path(output_dir, png_name), plot,
    width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(output_dir, pdf_name), plot,
    width = width, height = height, bg = "white", limitsize = FALSE
  )

  record_artifact(
    paste(png_name, pdf_name, sep = "; "), full_title, "PNG and PDF",
    source_tables, interpretation
  )
}

make_dotplot <- function(term_ids, title, subtitle, stem) {
  plot_data <- final_cell_plot |>
    filter(display_reactome_id %in% term_ids)

  term_order <- order_terms_for_plot(term_ids)
  cell_order <- cell_order_from_data(plot_data, term_ids)

  plot_data <- plot_data |>
    mutate(
      display_label = factor(display_label, levels = rev(term_order)),
      cell_type = factor(cell_type, levels = cell_order)
    )

  p <- ggplot(plot_data, aes(x = cell_type, y = display_label)) +
    geom_point(
      data = filter(plot_data, direction != "mixed"),
      aes(
        fill = evidence_NES,
        size = plot_log_fdr,
        shape = evidence_mode
      ),
      colour = "#252525", stroke = 0.35, alpha = 0.95
    ) +
    geom_point(
      data = filter(plot_data, direction == "mixed"),
      aes(size = plot_log_fdr, shape = evidence_mode),
      fill = colour_mixed, colour = "#252525", stroke = 0.35
    ) +
    scale_fill_gradient2(
      low = colour_negative, mid = "white", high = colour_positive,
      midpoint = 0,
      limits = c(-nes_colour_limit, nes_colour_limit),
      oob = squish,
      name = "Evidence NES"
    ) +
    scale_size_continuous(
      range = c(1.8, 7),
      name = expression(paste("Capped ", -log[10], "(FDR)"))
    ) +
    scale_shape_manual(
      values = c("Direct displayed term" = 21, "Supporting member" = 24),
      name = "Evidence source"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Cell population",
      y = NULL,
      caption = wrap_caption(
        "Red: positive NES (enriched toward Crohn's disease); blue: negative NES",
        "(enriched toward Normal). Purple marks mixed within-cell evidence. Blank",
        "combinations were not retained as significant and are not zero effects."
      )
    ) +
    theme_reactome +
    theme(
      axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 7.8),
      panel.grid.major = element_blank()
    )

  save_plot_pair(
    p, stem, title,
    width = max(12, 0.36 * length(cell_order) + 5),
    height = max(7, 0.31 * length(term_ids) + 3),
    source_tables = paste(
      "adaptive_reactome_final_significant_terms.csv",
      "adaptive_reactome_final_by_cell.csv",
      "adaptive_display_mapping.csv",
      sep = "; "
    ),
    interpretation = paste(
      "Each point is a retained significant cell-population/Reactome-term result.",
      "Colour is signed NES, size is FDR strength, and shape indicates whether",
      "the statistic came directly from the displayed term or a supporting member."
    )
  )
}

make_heatmap <- function(term_ids, title, subtitle, stem) {
  plot_data <- final_cell_plot |>
    filter(display_reactome_id %in% term_ids)

  term_order <- order_terms_for_plot(term_ids)
  cell_order <- cell_order_from_data(plot_data, term_ids)

  full_grid <- expand_grid(
    display_label = factor(rev(term_order), levels = rev(term_order)),
    cell_type = factor(cell_order, levels = cell_order)
  )

  plot_data <- plot_data |>
    mutate(
      display_label = factor(display_label, levels = rev(term_order)),
      cell_type = factor(cell_type, levels = cell_order)
    )

  p <- ggplot() +
    geom_tile(
      data = full_grid,
      aes(x = cell_type, y = display_label),
      fill = colour_missing, colour = "white", linewidth = 0.18
    ) +
    geom_tile(
      data = filter(plot_data, direction != "mixed"),
      aes(x = cell_type, y = display_label, fill = evidence_NES),
      colour = "white", linewidth = 0.18
    ) +
    geom_tile(
      data = filter(plot_data, direction == "mixed"),
      aes(x = cell_type, y = display_label),
      fill = colour_mixed, colour = "white", linewidth = 0.18
    ) +
    geom_point(
      data = filter(plot_data, evidence_q_value <= 0.01),
      aes(x = cell_type, y = display_label),
      shape = 4, size = 0.8, stroke = 0.3, colour = "#111111"
    ) +
    scale_fill_gradient2(
      low = colour_negative, mid = "white", high = colour_positive,
      midpoint = 0,
      limits = c(-nes_colour_limit, nes_colour_limit),
      oob = squish,
      name = "Evidence NES"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Cell population",
      y = NULL,
      caption = wrap_caption(
        "Grey = no retained significant result; purple = mixed within-cell evidence;",
        "x = FDR <= 0.01. Missing cells are not interpreted as zero NES or biological absence."
      )
    ) +
    theme_reactome +
    theme(
      axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 7.8),
      panel.grid = element_blank()
    )

  save_plot_pair(
    p, stem, title,
    width = max(12, 0.36 * length(cell_order) + 5),
    height = max(7, 0.30 * length(term_ids) + 3),
    source_tables = paste(
      "adaptive_reactome_final_significant_terms.csv",
      "adaptive_reactome_final_by_cell.csv",
      "adaptive_display_mapping.csv",
      sep = "; "
    ),
    interpretation = paste(
      "The heatmap shows the direction and strength of retained significant evidence",
      "for each selected term and cell population. It does not impute unreported cells."
    )
  )
}

overview_ids <- overview_terms$display_reactome_id

make_dotplot(
  overview_ids,
  "Adaptive Reactome enrichment overview",
  "Balanced selection of broadly recurrent and strong specific terms",
  "01_overview_dotplot"
)

make_heatmap(
  overview_ids,
  "Adaptive Reactome NES heatmap",
  "The same auditable term selection used in the overview dot plot",
  "02_overview_heatmap"
)

label_terms <- final_terms_plot |>
  arrange(desc(presentation_score)) |>
  slice_head(n = 16L) |>
  pull(display_reactome_id)

p_breadth <- ggplot(
  final_terms_plot,
  aes(
    x = n_significant_cell_types,
    y = largest_member_abs_NES,
    fill = overall_direction,
    shape = resolution,
    size = q_strength
  )
) +
  geom_point(colour = "#242424", stroke = 0.35, alpha = 0.82) +
  geom_text_repel(
    data = filter(final_terms_plot, display_reactome_id %in% label_terms),
    aes(label = str_wrap(display_term, 28)),
    size = 2.7, max.overlaps = Inf, box.padding = 0.35,
    point.padding = 0.15, min.segment.length = 0,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = direction_colours, name = "Overall direction") +
  scale_shape_manual(values = resolution_shapes, name = "Resolution") +
  scale_size_continuous(range = c(2, 7), name = "Capped -log10(FDR)") +
  labs(
    title = "Breadth versus enrichment strength of final Reactome terms",
    subtitle = "Terms are assessed using the same evidence and specificity criteria",
    x = "Number of cell populations with significant evidence",
    y = "Largest absolute supporting NES",
    caption = wrap_caption(
      "Upper right: strong, broadly recurrent processes. Upper left: strong but",
      "cell-restricted processes. NES describes coordinated enrichment, not proven activation."
    )
  ) +
  theme_reactome

save_plot_pair(
  p_breadth,
  "05_breadth_vs_strength",
  "Breadth versus enrichment strength of final Reactome terms",
  width = 11, height = 8,
  source_tables = paste(
    "adaptive_reactome_final_significant_terms.csv",
    "adaptive_display_mapping.csv",
    sep = "; "
  ),
  interpretation = paste(
    "Separates widespread signals from strong cell-restricted signals using",
    "recurrence and supporting enrichment strength."
  )
)

# Similarity is computed from sets of retained Reactome-ID + direction labels. It is
# intentionally not computed by treating absent NES values as zero.
cell_names <- sort(unique(final_cell_plot$cell_type))
cell_sets <- setNames(
  lapply(cell_names, function(cell) {
    final_cell_plot |>
      filter(cell_type == cell) |>
      transmute(key = paste(display_reactome_id, direction, sep = "|")) |>
      pull(key) |>
      unique()
  }),
  cell_names
)

jaccard_matrix <- matrix(
  NA_real_, nrow = length(cell_names), ncol = length(cell_names),
  dimnames = list(cell_names, cell_names)
)

for (i in seq_along(cell_names)) {
  for (j in seq_along(cell_names)) {
    union_size <- length(union(cell_sets[[i]], cell_sets[[j]]))
    intersection_size <- length(intersect(cell_sets[[i]], cell_sets[[j]]))
    jaccard_matrix[i, j] <- if (union_size == 0) NA_real_ else {
      intersection_size / union_size
    }
  }
}

similarity_order <- if (length(cell_names) > 2) {
  cell_names[hclust(as.dist(1 - replace(jaccard_matrix, is.na(jaccard_matrix), 0)))$order]
} else {
  cell_names
}

jaccard_long <- as.data.frame(as.table(jaccard_matrix)) |>
  as_tibble() |>
  rename(cell_type_1 = Var1, cell_type_2 = Var2, jaccard = Freq) |>
  mutate(
    cell_type_1 = factor(cell_type_1, levels = rev(similarity_order)),
    cell_type_2 = factor(cell_type_2, levels = similarity_order)
  )

write_csv(
  jaccard_long |>
    mutate(across(starts_with("cell_type_"), as.character)),
  file.path(output_dir, "cell_programme_jaccard_matrix_long.csv")
)

p_jaccard <- ggplot(
  jaccard_long,
  aes(x = cell_type_2, y = cell_type_1, fill = jaccard)
) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradient(
    low = "white", high = "#006D77", limits = c(0, 1),
    na.value = colour_missing, name = "Jaccard"
  ) +
  coord_equal() +
  labs(
    title = "Shared significant Reactome programmes between cell populations",
    subtitle = "Jaccard similarity of retained Reactome-term + direction sets",
    x = NULL, y = NULL,
    caption = wrap_caption(
      "This measures overlap in significant enriched programmes, not overall",
      "transcriptomic similarity. Power and annotation redundancy also influence overlap."
    )
  ) +
  theme_reactome +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 7.5),
    axis.text.y = element_text(size = 7.5),
    panel.grid = element_blank()
  )

save_plot_pair(
  p_jaccard,
  "06_cell_programme_jaccard_heatmap",
  "Shared significant Reactome programmes between cell populations",
  width = 11, height = 10,
  source_tables = "adaptive_reactome_final_by_cell.csv; adaptive_reactome_final_significant_terms.csv",
  interpretation = paste(
    "Cells close in this heatmap share more retained significant Reactome terms in the",
    "same direction; this is not a whole-transcriptome similarity measure."
  )
)

cell_burden <- final_cell_plot |>
  count(cell_type, direction, name = "n_terms") |>
  group_by(cell_type) |>
  mutate(total = sum(n_terms)) |>
  ungroup() |>
  mutate(cell_type = reorder(cell_type, total))

p_burden <- ggplot(
  cell_burden,
  aes(x = n_terms, y = cell_type, fill = direction)
) +
  geom_col(width = 0.76, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = direction_colours, name = "Direction") +
  labs(
    title = "Retained significant Reactome-term burden by cell population",
    x = "Number of final retained terms",
    y = NULL,
    caption = wrap_caption(
      "Counts reflect statistical power, sample variability and Reactome hierarchy structure",
      "as well as biology; they should not be read as a direct measure of disease severity."
    )
  ) +
  theme_reactome

save_plot_pair(
  p_burden,
  "07_cell_term_burden",
  "Retained significant Reactome-term burden by cell population",
  width = 10, height = 8.5,
  source_tables = "adaptive_reactome_final_by_cell.csv; adaptive_reactome_final_significant_terms.csv",
  interpretation = paste(
    "Shows how many final programmes are represented in each cell population and",
    "their enrichment direction; counts are sensitive to power and annotation."
  )
)

resolution_summary <- final_terms_plot |>
  count(resolution, overall_direction, name = "n_terms") |>
  complete(
    resolution,
    overall_direction = c("positive", "negative", "mixed"),
    fill = list(n_terms = 0L)
  )

p_resolution <- ggplot(
  resolution_summary,
  aes(x = overall_direction, y = resolution, fill = n_terms)
) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = comma(n_terms)), fontface = "bold", size = 4) +
  scale_fill_gradient(low = "white", high = "#3F6C8E", name = "Final terms") +
  labs(
    title = "Composition of the final non-redundant Reactome list",
    x = "Overall direction across cell populations",
    y = "Presentation resolution",
    caption = wrap_caption(
      "Granular terms retain specific evidence. Ontology roll-ups and semantic groups",
      "summarise less distinctive or redundant signals at a broader level."
    )
  ) +
  theme_reactome +
  theme(panel.grid = element_blank())

save_plot_pair(
  p_resolution,
  "08_resolution_direction_summary",
  "Composition of the final non-redundant Reactome list",
  width = 8.5, height = 5.5,
  source_tables = "adaptive_reactome_final_significant_terms.csv",
  interpretation = paste(
    "Counts final terms by presentation resolution and overall enrichment direction,",
    "making the balance between specificity and roll-up explicit."
  )
)

count_value <- function(metric_name) {
  value <- reduction_counts |>
    filter(metric == metric_name) |>
    pull(value)
  if (length(value) != 1) {
    stop("Expected exactly one reduction_counts value for: ", metric_name)
  }
  as.numeric(value)
}

reduction_plot_data <- tribble(
  ~unit, ~stage, ~stage_order, ~value,
  "Cell x term evidence rows", "Mapped significant input", 1,
  count_value("mapped_significant_input_rows"),
  "Cell x term evidence rows", "After semantic-cluster selection", 2,
  count_value("compact_member_rows_after_cluster_selection"),
  "Cell x term evidence rows", "Before within-cell parent pruning", 3,
  count_value("cell_term_rows_before_parent_pruning"),
  "Cell x term evidence rows", "After within-cell parent pruning", 4,
  count_value("cell_term_rows_after_parent_pruning"),
  "Unique Reactome terms", "Mapped candidate Reactome IDs", 1,
  count_value("mapped_candidate_reactome_ids"),
  "Unique Reactome terms", "Before global parent pruning", 2,
  count_value("unique_terms_before_global_parent_pruning"),
  "Unique Reactome terms", "Final non-redundant list", 3,
  count_value("final_unique_significant_terms")
) |>
  group_by(unit) |>
  mutate(
    baseline = first(value),
    percent_of_baseline = 100 * value / baseline,
    label = paste0(comma(value), "  (", number(percent_of_baseline, accuracy = 0.1), "%)"),
    stage_panel = paste(unit, stage, sep = "___")
  ) |>
  ungroup() |>
  mutate(stage_panel = factor(stage_panel, levels = rev(stage_panel)))

p_reduction <- ggplot(
  reduction_plot_data,
  aes(x = value, y = stage_panel)
) +
  geom_col(fill = "#4C78A8", width = 0.7) +
  geom_text(
    aes(label = label), hjust = -0.08, size = 3.1,
    colour = colour_text
  ) +
  facet_grid(unit ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub(".*___", "", x)) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.28))) +
  labs(
    title = "How the adaptive summary reduces Reactome redundancy",
    subtitle = "Separate panels preserve the distinction between row counts and unique terms",
    x = "Count",
    y = NULL,
    caption = wrap_caption(
      "Percentages are relative to the first stage within each panel. The two panels",
      "use different units and should not be interpreted as one continuous funnel."
    )
  ) +
  theme_reactome +
  theme(
    strip.text.y = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

save_plot_pair(
  p_reduction,
  "09_reduction_summary",
  "How the adaptive summary reduces Reactome redundancy",
  width = 10.5, height = 7.5,
  source_tables = "reduction_counts.csv",
  interpretation = paste(
    "Tracks evidence rows and unique Reactome terms separately through selection and",
    "parent-pruning stages; the final list contains no parent-child duplicates."
  )
)

selection_plot_data <- global_selection |>
  mutate(
    selected = if_else(
      as_logical_flag(global_granular_selected), "Selected", "Not selected"
    )
  ) |>
  count(member_direction, selected, name = "n_terms")

p_selection <- ggplot(
  selection_plot_data,
  aes(x = member_direction, y = n_terms, fill = selected)
) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(
    aes(label = comma(n_terms)),
    position = position_dodge(width = 0.7), vjust = -0.25, size = 3
  ) +
  scale_fill_manual(
    values = c("Selected" = "#2A9D8F", "Not selected" = "#BDBDBD"),
    name = "Global granular status"
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Global granular-term selection audit",
    x = "Direction stratum",
    y = "Candidate Reactome terms",
    caption = wrap_caption(
      "Selection preserves the strongest informative granular representatives within",
      "each semantic/direction group; unselected terms may still contribute to a roll-up."
    )
  ) +
  theme_reactome +
  theme(panel.grid.major.x = element_blank())

save_plot_pair(
  p_selection,
  "10_global_granular_selection_audit",
  "Global granular-term selection audit",
  width = 9, height = 5.8,
  source_tables = "adaptive_reactome_global_granular_selection.csv",
  interpretation = paste(
    "Shows the number of candidate terms selected or not selected for granular",
    "presentation, stratified by enrichment direction."
  )
)

term_metric_plot_data <- term_metrics |>
  mutate(
    evidence_pass = as_logical_flag(evidence_pass),
    keep_granular = as_logical_flag(keep_granular),
    presentation_score =
      0.40 * scale_01(pmin(safe_neg_log10(best_q), 20)) +
      0.25 * scale_01(pmin(max_abs_nes, 3)) +
      0.15 * scale_01(log1p(n_significant_cell_types)) +
      0.20 * scale_01(information_content),
    presentation_outcome = case_when(
      keep_granular ~ "Kept granular",
      evidence_pass ~ "Evidence passed; summarised",
      TRUE ~ "Lower granular evidence"
    )
  ) |>
  filter(
    is.finite(max_depth),
    is.finite(information_content),
    is.finite(presentation_score)
  )

metric_labels <- term_metric_plot_data |>
  arrange(desc(presentation_score), desc(max_abs_nes)) |>
  slice_head(n = 10L)

p_metric_landscape <- ggplot(
  term_metric_plot_data,
  aes(
    x = max_depth,
    y = information_content,
    colour = presentation_score,
    shape = presentation_outcome
  )
) +
  geom_point(alpha = 0.48, size = 1.8) +
  geom_text_repel(
    data = metric_labels,
    aes(label = str_wrap(reactome_term, 28)),
    size = 2.6, max.overlaps = Inf, min.segment.length = 0,
    box.padding = 0.35, point.padding = 0.15,
    show.legend = FALSE
  ) +
  scale_colour_gradient(
    low = "#BDBDBD", high = "#007C83",
    name = "Presentation score"
  ) +
  scale_shape_manual(
    values = c(
      "Kept granular" = 16,
      "Evidence passed; summarised" = 15,
      "Lower granular evidence" = 1
    ),
    name = "Adaptive outcome"
  ) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(
    title = "Reactome specificity and adaptive presentation landscape",
    subtitle = "All Reactome candidates use the same evidence and specificity criteria",
    x = "Maximum Reactome depth",
    y = "Information content",
    caption = str_wrap(
      paste(
        "Deeper and higher-information terms are more specific. Presentation score also",
        "uses statistical strength, effect size, recurrence and information content, so depth alone",
        "does not determine whether a term is retained granularly or summarised."
      ),
      width = 135
    )
  ) +
  theme_reactome

save_plot_pair(
  p_metric_landscape,
  "10b_reactome_term_metric_landscape",
  "Reactome specificity and adaptive presentation landscape",
  width = 11, height = 8,
  source_tables = "reactome_term_metrics.csv",
  interpretation = paste(
    "Shows how ontology depth and information content relate to the adaptive",
    "presentation outcome; it also demonstrates that hierarchy depth was not used",
    "as the sole selection criterion."
  )
)

layout_dag <- function(nodes, edges) {
  graph <- graph_from_data_frame(
    edges |>
      select(parent, child),
    directed = TRUE,
    vertices = nodes |>
      transmute(name = reactome_id)
  )

  layout_result <- layout_with_sugiyama(graph)$layout
  coordinates <- tibble(
    reactome_id = V(graph)$name,
    x = layout_result[, 1],
    y = -layout_result[, 2]
  )

  edge_coordinates <- edges |>
    left_join(
      coordinates |>
        rename(parent = reactome_id, x_parent = x, y_parent = y),
      by = "parent"
    ) |>
    left_join(
      coordinates |>
        rename(child = reactome_id, x_child = x, y_child = y),
      by = "child"
    )

  list(
    nodes = nodes |> left_join(coordinates, by = "reactome_id"),
    edges = edge_coordinates,
    graph = graph
  )
}

plot_dag <- function(
  nodes, edges, title, subtitle, caption, stem,
  source_tables, edge_type_values
) {
  if (nrow(nodes) == 0 || nrow(edges) == 0) {
    record_artifact(
      paste0(stem, ".png; ", stem, ".pdf"), scoped_title(title), "PNG and PDF",
      source_tables, "No connected nodes were available.", "skipped"
    )
    return(invisible(NULL))
  }

  laid_out <- layout_dag(nodes, edges)

  p <- ggplot() +
    geom_segment(
      data = laid_out$edges,
      aes(
        x = x_parent, y = y_parent,
        xend = x_child, yend = y_child,
        linetype = relationship
      ),
      colour = "#737373", linewidth = 0.45, alpha = 0.75,
      arrow = grid::arrow(length = grid::unit(0.09, "inches"), type = "closed")
    ) +
    geom_point(
      data = laid_out$nodes,
      aes(
        x = x, y = y, fill = fill_category,
        shape = node_role, size = node_size
      ),
      colour = "#303030", stroke = 0.5
    ) +
    geom_text_repel(
      data = filter(laid_out$nodes, show_label),
      aes(x = x, y = y, label = label),
      size = 2.6, max.overlaps = Inf, min.segment.length = 0,
      box.padding = 0.3, point.padding = 0.2, seed = 14,
      colour = colour_text
    ) +
    scale_fill_manual(
      values = c(
        "Retained: Crohn enriched" = colour_positive,
        "Retained: Normal enriched" = colour_negative,
        "Retained: mixed" = colour_mixed,
        "Significant parent removed" = "#F2C66D",
        "Ontology connector" = colour_connector
      ),
      name = "Node meaning",
      drop = FALSE
    ) +
    scale_shape_manual(
      values = c(
        "Retained final term" = 21,
        "Significant pruned parent" = 22,
        "Ontology connector" = 24
      ),
      name = "Node role",
      drop = FALSE
    ) +
    scale_size_continuous(range = c(2.5, 8), name = "Significant cells") +
    scale_linetype_manual(values = edge_type_values, name = "Relationship") +
    coord_equal(clip = "off") +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL, y = NULL,
      caption = caption
    ) +
    theme_void(base_size = 10, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "#4A4A4A"),
      plot.caption = element_text(colour = "#5A5A5A", hjust = 0),
      legend.position = "right",
      plot.margin = margin(12, 30, 12, 12)
    )

  save_plot_pair(
    p, stem, title,
    width = 14, height = 10,
    source_tables = source_tables,
    interpretation = caption
  )

  invisible(laid_out)
}

# A focused graph of the global pruning decisions. These edges are transitive
# "strict ancestor of" audit relationships, not necessarily immediate Reactome edges.
global_pruning <- pruning_audit |>
  filter(pruning_scope == "global_unique_list") |>
  mutate(
    n_descendants = if_else(
      is.na(retained_descendant_reactome_ids) | retained_descendant_reactome_ids == "",
      0L,
      str_count(retained_descendant_reactome_ids, fixed(";")) + 1L
    )
  ) |>
  arrange(desc(n_descendants), best_member_q_value)

chosen_pruning <- global_pruning |>
  slice_head(n = max_dechaining_parents)

dechaining_edges <- chosen_pruning |>
  select(
    parent = removed_reactome_id,
    parent_term = removed_term,
    parent_direction = direction,
    parent_n_significant_cells = n_significant_cell_types,
    total_retained_descendants = n_descendants,
    retained_descendant_reactome_ids
  ) |>
  separate_rows(retained_descendant_reactome_ids, sep = ";") |>
  transmute(
    parent,
    parent_term,
    parent_direction,
    parent_n_significant_cells,
    total_retained_descendants,
    child = str_trim(retained_descendant_reactome_ids),
    relationship = "strict ancestor of (possibly indirect)"
  ) |>
  filter(!is.na(child), child != "") |>
  distinct() |>
  left_join(
    final_terms_plot |>
      transmute(
        child = display_reactome_id,
        child_term = display_term,
        child_direction = overall_direction,
        child_n_significant_cells = as.numeric(n_significant_cell_types),
        child_presentation_score = presentation_score
      ),
    by = "child"
  ) |>
  mutate(
    child_term = coalesce(child_term, child),
    child_direction = coalesce(child_direction, "mixed"),
    child_n_significant_cells = coalesce(child_n_significant_cells, 1),
    child_presentation_score = coalesce(child_presentation_score, 0)
  ) |>
  group_by(parent) |>
  arrange(
    desc(child_n_significant_cells),
    desc(child_presentation_score),
    .by_group = TRUE
  ) |>
  slice_head(n = max_descendants_per_parent) |>
  mutate(
    child_rank = row_number(),
    parent_y = (n() + 1) / 2,
    parent_label = paste0(
      str_wrap(first(parent_term), 42), "\n", first(parent),
      "  (", first(total_retained_descendants), " retained descendants total)"
    )
  ) |>
  ungroup()

if (nrow(dechaining_edges) > 0) {
  parent_nodes <- dechaining_edges |>
    distinct(
      parent, parent_term, parent_direction,
      parent_n_significant_cells, parent_y, parent_label
    ) |>
    transmute(
      reactome_id = parent,
      term = parent_term,
      direction = parent_direction,
      n_significant_cells = as.numeric(parent_n_significant_cells),
      node_role = "Significant pruned parent"
    )

  child_nodes <- dechaining_edges |>
    transmute(
      reactome_id = child,
      term = child_term,
      direction = child_direction,
      n_significant_cells = child_n_significant_cells,
      node_role = "Retained final term"
    )

  dechaining_nodes <- bind_rows(parent_nodes, child_nodes) |>
    distinct(reactome_id, node_role, .keep_all = TRUE)

  write_csv(dechaining_nodes, file.path(output_dir, "dechaining_graph_nodes.csv"))
  write_csv(dechaining_edges, file.path(output_dir, "dechaining_graph_edges.csv"))

  p_dechaining <- ggplot() +
    geom_segment(
      data = dechaining_edges,
      aes(x = 0, y = parent_y, xend = 1, yend = child_rank),
      colour = "#777777", linewidth = 0.45, linetype = "dotted",
      arrow = grid::arrow(length = grid::unit(0.07, "inches"), type = "closed")
    ) +
    geom_point(
      data = dechaining_edges |>
        distinct(parent, parent_label, parent_y),
      aes(x = 0, y = parent_y),
      shape = 22, size = 5, fill = "#F2C66D", colour = "#303030",
      stroke = 0.6
    ) +
    geom_point(
      data = dechaining_edges,
      aes(
        x = 1, y = child_rank, fill = child_direction,
        size = child_n_significant_cells
      ),
      shape = 21, colour = "#303030", stroke = 0.45
    ) +
    geom_text(
      data = dechaining_edges,
      aes(
        x = 1.08, y = child_rank,
        label = paste0(str_wrap(child_term, 27), "\n", child)
      ),
      hjust = 0, size = 2.45, lineheight = 0.9, colour = colour_text
    ) +
    facet_wrap(~parent_label, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = direction_colours, name = "Retained direction") +
    scale_size_continuous(range = c(2.5, 6), name = "Significant cells") +
    scale_x_continuous(
      limits = c(-0.15, 2.35),
      breaks = c(0, 1),
      labels = c("Removed parent", "Retained descendant")
    ) +
    labs(
      title = "How redundant parent terms were removed",
      subtitle = paste(
        max_dechaining_parents, " largest global pruning decisions; up to ",
        max_descendants_per_parent, " broad retained",
        "descendants shown per parent"
      ),
      x = NULL, y = NULL,
      caption = wrap_caption(
        "Gold squares are significant parent terms removed because more specific",
        "descendants were retained. Dotted arrows mean strict ancestor-of and may span",
        "several direct Reactome edges; this is a decision audit, not the literal Reactome DAG."
      )
    ) +
    theme_void(base_size = 9, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(colour = "#4A4A4A"),
      plot.caption = element_text(colour = "#5A5A5A", hjust = 0),
      strip.text = element_text(face = "bold", size = 8, lineheight = 0.95),
      strip.background = element_rect(fill = "#F3F3F3", colour = "#D0D0D0"),
      axis.text.x = element_text(size = 7.5, colour = colour_text),
      legend.position = "right",
      panel.spacing = grid::unit(1.1, "lines"),
      plot.margin = margin(10, 35, 10, 10)
    )

  save_plot_pair(
    p_dechaining,
    "11_dechaining_audit_graph",
    "How redundant parent terms were removed",
    width = 14, height = 13,
    source_tables = paste(
      "adaptive_reactome_hierarchy_pruning_audit.csv",
      "adaptive_reactome_final_significant_terms.csv",
      sep = "; "
    ),
    interpretation = paste(
      "Each facet links one removed significant parent to a readable subset of",
      "retained, more specific descendants. Edges are transitive ancestry summaries,",
      "not necessarily direct Reactome edges."
    )
  )
}

collect_parent_neighbourhood <- function(seed_ids, edges, n_hops = 2L) {
  nodes <- tibble(reactome_id = unique(seed_ids), ancestor_hops = 0L)
  frontier <- nodes$reactome_id

  for (hop in seq_len(n_hops)) {
    new_parents <- edges |>
      filter(child %in% frontier) |>
      pull(parent) |>
      unique() |>
      setdiff(nodes$reactome_id)

    if (length(new_parents) == 0) break

    nodes <- bind_rows(
      nodes,
      tibble(reactome_id = new_parents, ancestor_hops = hop)
    )
    frontier <- new_parents
  }

  graph_edges <- edges |>
    filter(parent %in% nodes$reactome_id, child %in% nodes$reactome_id)

  list(nodes = nodes, edges = graph_edges)
}

snapshot_candidates <- c(
  file.path(input_dir, "reactome_hierarchy_snapshot.rds"),
  file.path(dirname(dirname(input_dir)), "scopes", "whole_atlas", "reactome_hierarchy_snapshot.rds")
)
snapshot_file <- snapshot_candidates[file.exists(snapshot_candidates)][1]
hierarchy <- tryCatch({
  if (!is.na(snapshot_file)) {
    load_reactome_hierarchy(snapshot_file = snapshot_file)
  } else {
    resources <- Sys.getenv("REACTOME_RESOURCES_DIR",
      unset = path.expand(paste0("~/project-area/data/crohns_scrnaseq/10c_14n_analysis/",
        "pseudobulk_dea/pathway_enrichment/resources/reactome")))
    load_reactome_hierarchy(resources_dir = resources)
  }
}, error = function(e) {
  warning("Reactome hierarchy context unavailable: ", conditionMessage(e))
  NULL
})
reactome_hierarchy_available <- !is.null(hierarchy)
reactome_hierarchy_source <- if (reactome_hierarchy_available) hierarchy$source else NA_character_
reactome_dag_created <- FALSE

if (reactome_hierarchy_available) {
  reactome_edges <- as_tibble(hierarchy$edges)

  seed_candidates <- final_terms_plot |>
    mutate(in_overview = display_reactome_id %in% overview_ids) |>
    arrange(
      desc(in_overview),
      desc(presentation_score)
    ) |>
    pull(display_reactome_id) |>
    unique()

  n_seed <- min(max_graph_terms, length(seed_candidates))
  dag_context <- NULL

  repeat {
    dag_seed_ids <- head(seed_candidates, n_seed)
    dag_context <- collect_parent_neighbourhood(
      dag_seed_ids, reactome_edges, graph_parent_hops
    )

    if (nrow(dag_context$nodes) <= max_graph_nodes || n_seed <= 3L) break
    n_seed <- n_seed - 1L
  }

  reactome_annotation <- as_tibble(hierarchy$terms) |>
    filter(reactome_id %in% dag_context$nodes$reactome_id) |>
    select(reactome_id, reactome_term)

  global_pruned_nodes <- global_pruning |>
    transmute(
      reactome_id = removed_reactome_id,
      pruned_term = removed_term,
      pruned_direction = direction,
      pruned_n_cells = as.numeric(n_significant_cell_types)
    ) |>
    distinct(reactome_id, .keep_all = TRUE)

  dag_nodes <- dag_context$nodes |>
    left_join(reactome_annotation, by = "reactome_id") |>
    left_join(
      final_terms_plot |>
        transmute(
          reactome_id = display_reactome_id,
          retained_term = display_term,
          retained_direction = overall_direction,
          retained_n_cells = as.numeric(n_significant_cell_types)
        ),
      by = "reactome_id"
    ) |>
    left_join(global_pruned_nodes, by = "reactome_id") |>
    mutate(
      node_role = case_when(
        reactome_id %in% dag_seed_ids ~ "Retained final term",
        !is.na(pruned_term) ~ "Significant pruned parent",
        TRUE ~ "Ontology connector"
      ),
      term = coalesce(retained_term, pruned_term, reactome_term, reactome_id),
      direction = coalesce(retained_direction, pruned_direction),
      node_size = case_when(
        node_role == "Retained final term" ~ pmax(retained_n_cells, 1),
        node_role == "Significant pruned parent" ~ pmax(pruned_n_cells, 1),
        TRUE ~ 1
      ),
      fill_category = case_when(
        node_role == "Retained final term" & direction == "positive" ~
          "Retained: Crohn enriched",
        node_role == "Retained final term" & direction == "negative" ~
          "Retained: Normal enriched",
        node_role == "Retained final term" ~ "Retained: mixed",
        node_role == "Significant pruned parent" ~ "Significant parent removed",
        TRUE ~ "Ontology connector"
      )
    )

  dag_graph <- graph_from_data_frame(
    dag_context$edges |> select(parent, child),
    directed = TRUE,
    vertices = dag_nodes |> transmute(name = reactome_id)
  )

  if (!is_dag(dag_graph)) {
    stop(
      "The selected Reactome pathway/subpathway subgraph contains a cycle; cannot use a DAG layout."
    )
  }

  selected_descendant_count <- vapply(
    dag_nodes$reactome_id,
    function(node_id) {
      reachable <- as_ids(subcomponent(dag_graph, node_id, mode = "out"))
      as.integer(sum(dag_seed_ids %in% reachable))
    },
    integer(1)
  )

  dag_nodes <- dag_nodes |>
    mutate(
      selected_descendant_count = selected_descendant_count,
      show_label = node_role == "Retained final term" |
        selected_descendant_count >= 2L,
      label = if_else(show_label, str_wrap(term, 25), NA_character_)
    )

  selected_missing_from_reactome_hierarchy <- setdiff(dag_seed_ids, reactome_annotation$reactome_id)
  if (length(selected_missing_from_reactome_hierarchy) > 0) {
    warning(
      "Selected Reactome IDs absent from the offline Reactome hierarchy: ",
      paste(selected_missing_from_reactome_hierarchy, collapse = ", ")
    )
    write_lines(
      selected_missing_from_reactome_hierarchy,
      file.path(output_dir, "reactome_dag_selected_ids_missing_from_hierarchy.txt")
    )
  }

  write_csv(dag_nodes, file.path(output_dir, "reactome_dag_nodes.csv"))
  write_csv(dag_context$edges, file.path(output_dir, "reactome_dag_edges.csv"))

  dag_plot_result <- plot_dag(
    dag_nodes,
    dag_context$edges,
    "Reactome hierarchy context for selected final terms",
    paste0(
      length(dag_seed_ids), " retained terms plus up to ", graph_parent_hops,
      " direct ancestor levels; ", reactome_hierarchy_source
    ),
    paste(
      "Coloured circles are retained significant processes. Gold squares are significant",
      "ancestors removed as redundant; grey triangles supply ontology context only.",
      "Solid arrows are direct parent/subpathway relationships; dashed arrows link the",
      "synthetic structural root. Arrows run from broader to narrower pathways. Multiple parents",
      "are retained because Reactome is a directed acyclic graph, not a tree."
    ),
    "12_reactome_dag_context",
    paste(
      "adaptive_reactome_final_significant_terms.csv",
      "adaptive_reactome_hierarchy_pruning_audit.csv",
      reactome_hierarchy_source,
      sep = "; "
    ),
    c("has_subpathway" = "solid", "structural_root" = "dashed")
  )

  reactome_dag_created <- !is.null(dag_plot_result)
} else {
  warning(
    "Offline Reactome hierarchy files/snapshot are unavailable. The hierarchy DAG was skipped; ",
    "all other figures were created."
  )
  record_artifact(
    "12_reactome_dag_context.png; 12_reactome_dag_context.pdf",
    scoped_title("Reactome hierarchy context for selected final terms"),
    "PNG and PDF",
    "ReactomePathways.txt; ReactomePathwaysRelation.txt",
    paste(
      "Supply the offline Reactome hierarchy resources or rerun script6 to create its snapshot.",
      "Then rerun this script for the direct pathway/subpathway hierarchy DAG."
    ),
    "skipped: missing Reactome hierarchy resources"
  )
}

# Representative leading-edge recurrence for a small subset of overview terms.
leading_edge_term_ids <- overview_terms |>
  arrange(desc(presentation_score)) |>
  slice_head(n = 12L) |>
  pull(display_reactome_id)

leading_edge_counts <- final_cell_plot |>
  filter(display_reactome_id %in% leading_edge_term_ids) |>
  select(display_reactome_id, display_term, cell_type, evidence_leading_edge) |>
  filter(!is.na(evidence_leading_edge), evidence_leading_edge != "") |>
  separate_rows(evidence_leading_edge, sep = ";") |>
  transmute(
    display_reactome_id,
    display_term,
    cell_type,
    gene = str_trim(evidence_leading_edge)
  ) |>
  filter(gene != "") |>
  distinct(display_reactome_id, display_term, cell_type, gene) |>
  count(display_reactome_id, display_term, gene, name = "n_cell_populations")

if (nrow(leading_edge_counts) > 0) {
  top_leading_genes <- leading_edge_counts |>
    group_by(gene) |>
    summarise(
      n_terms = n_distinct(display_reactome_id),
      total_cell_occurrences = sum(n_cell_populations),
      .groups = "drop"
    ) |>
    arrange(desc(n_terms), desc(total_cell_occurrences), gene) |>
    slice_head(n = 30L) |>
    pull(gene)

  leading_edge_plot_data <- leading_edge_counts |>
    filter(gene %in% top_leading_genes) |>
    mutate(
      term_label = wrap_term(display_term, display_reactome_id, width = 42),
      term_label = factor(
        term_label,
        levels = rev(order_terms_for_plot(leading_edge_term_ids))
      ),
      gene = factor(gene, levels = rev(top_leading_genes))
    )

  write_csv(
    leading_edge_plot_data |>
      mutate(across(c(term_label, gene), as.character)),
    file.path(output_dir, "leading_edge_recurrence.csv")
  )

  p_leading <- ggplot(
    leading_edge_plot_data,
    aes(x = gene, y = term_label, fill = n_cell_populations)
  ) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_gradient(
      low = "#E8F1F2", high = "#005F73",
      name = "Cell populations"
    ) +
    labs(
      title = "Recurrent genes in representative leading edges",
      subtitle = "Top genes across a focused set of final terms",
      x = "Leading-edge gene",
      y = NULL,
      caption = str_wrap(
        paste(
          "Counts show the number of cell populations in which a gene occurred in the",
          "representative evidence leading edge. Roll-ups may have additional supporting",
          "members, so this is evidence provenance rather than a complete pathway gene map."
        ),
        width = 155
      )
    ) +
    theme_reactome +
    theme(
      axis.text.x = element_text(angle = 55, hjust = 1, size = 7.5),
      axis.text.y = element_text(size = 7.4),
      panel.grid = element_blank()
    )

  save_plot_pair(
    p_leading,
    "13_leading_edge_recurrence",
    "Recurrent genes in representative leading edges",
    width = 13, height = 8,
    source_tables = "adaptive_reactome_final_by_cell.csv; adaptive_reactome_final_significant_terms.csv",
    interpretation = paste(
      "Highlights recurrent representative leading-edge genes while preserving the",
      "distinction between evidence provenance and complete pathway membership."
    )
  )
}

# Multi-page cell-specific supplement: strongest final terms in each cell.
per_cell_pdf <- file.path(output_dir, "14_top_terms_by_cell.pdf")
pdf(per_cell_pdf, width = 10.5, height = 8.2, onefile = TRUE)

for (cell in sort(unique(final_cell_plot$cell_type))) {
  cell_data <- final_cell_plot |>
    filter(cell_type == cell) |>
    arrange(evidence_q_value, desc(abs(evidence_NES))) |>
    slice_head(n = top_terms_per_cell) |>
    mutate(
      plot_nes = if_else(direction == "mixed", 0, evidence_NES),
      display_label = factor(
        wrap_term(display_term, display_reactome_id, width = 42),
        levels = rev(wrap_term(display_term, display_reactome_id, width = 42))
      )
    )

  p_cell <- ggplot(
    cell_data,
    aes(x = plot_nes, y = display_label)
  ) +
    geom_vline(xintercept = 0, colour = "#9E9E9E", linewidth = 0.4) +
    geom_segment(
      aes(x = 0, xend = plot_nes, yend = display_label, colour = direction),
      linewidth = 0.8, alpha = 0.7
    ) +
    geom_point(
      aes(
        fill = direction,
        size = plot_log_fdr,
        shape = evidence_mode
      ),
      colour = "#222222", stroke = 0.35
    ) +
    scale_colour_manual(values = direction_colours, guide = "none") +
    scale_fill_manual(values = direction_colours, name = "Direction") +
    scale_shape_manual(
      values = c("Direct displayed term" = 21, "Supporting member" = 24),
      name = "Evidence source"
    ) +
    scale_size_continuous(range = c(2.5, 7), name = "Capped -log10(FDR)") +
    labs(
      title = paste0(scope_label, ": Top retained Reactome terms: ", cell),
      subtitle = "Ranked first by evidence FDR, then by absolute NES",
      x = "Evidence NES",
      y = NULL,
      caption = wrap_caption(
        "Positive values are enriched toward Crohn's disease; negative values toward",
        "Normal. Mixed terms are shown at zero because one representative NES cannot",
        "summarise conflicting within-cell supporting evidence."
      )
    ) +
    theme_reactome +
    theme(axis.text.y = element_text(size = 8))

  print(p_cell)
}

dev.off()

record_artifact(
  basename(per_cell_pdf),
  scoped_title("Top retained Reactome terms for every cell population"),
  "Multi-page PDF",
  "adaptive_reactome_final_by_cell.csv; adaptive_reactome_final_significant_terms.csv",
  paste(
    "One page per cell population preserves strong cell-specific results that may",
    "not enter the cross-cell overview."
  )
)

# Multi-page all-term heatmap supplement. Each page uses a different subset of
# rows but the same cell ordering within that page. No absent value is imputed.
all_term_order <- final_terms_plot |>
  arrange(
    desc(n_significant_cell_types), desc(presentation_score)
  ) |>
  pull(display_reactome_id)

all_term_pages <- split(
  all_term_order,
  ceiling(seq_along(all_term_order) / supplement_terms_per_page)
)

all_heatmap_pdf <- file.path(output_dir, "15_all_final_terms_heatmap.pdf")
pdf(all_heatmap_pdf, width = 14, height = 12, onefile = TRUE)

for (page_index in seq_along(all_term_pages)) {
  ids <- all_term_pages[[page_index]]
  page_data <- final_cell_plot |>
    filter(display_reactome_id %in% ids)
  term_order <- order_terms_for_plot(ids)
  cell_order <- cell_order_from_data(page_data, ids)

  page_grid <- expand_grid(
    display_label = factor(rev(term_order), levels = rev(term_order)),
    cell_type = factor(cell_order, levels = cell_order)
  )

  page_data <- page_data |>
    mutate(
      display_label = factor(display_label, levels = rev(term_order)),
      cell_type = factor(cell_type, levels = cell_order)
    )

  p_page <- ggplot() +
    geom_tile(
      data = page_grid,
      aes(x = cell_type, y = display_label),
      fill = colour_missing, colour = "white", linewidth = 0.12
    ) +
    geom_tile(
      data = filter(page_data, direction != "mixed"),
      aes(x = cell_type, y = display_label, fill = evidence_NES),
      colour = "white", linewidth = 0.12
    ) +
    geom_tile(
      data = filter(page_data, direction == "mixed"),
      aes(x = cell_type, y = display_label),
      fill = colour_mixed, colour = "white", linewidth = 0.12
    ) +
    scale_fill_gradient2(
      low = colour_negative, mid = "white", high = colour_positive,
      midpoint = 0,
      limits = c(-nes_colour_limit, nes_colour_limit),
      oob = squish,
      name = "Evidence NES"
    ) +
    labs(
      title = paste0(
        scope_label, ": All final non-redundant Reactome terms - page ", page_index,
        " of ", length(all_term_pages)
      ),
      x = "Cell population", y = NULL,
      caption = "Grey = no retained significant result; purple = mixed within-cell evidence."
    ) +
    theme_reactome +
    theme(
      axis.text.x = element_text(angle = 50, hjust = 1, size = 7),
      axis.text.y = element_text(size = 6.3),
      panel.grid = element_blank()
    )

  print(p_page)
}

dev.off()

record_artifact(
  basename(all_heatmap_pdf),
  scoped_title("All final non-redundant Reactome terms"),
  "Multi-page PDF",
  "adaptive_reactome_final_by_cell.csv; adaptive_reactome_final_significant_terms.csv",
  paste(
    "Supplementary heatmap containing every final term without crowding the main",
    "overview. Grey cells are unreported retained evidence, not zero effects."
  )
)

# Inventory all upstream outputs. The three pre-final matrix/long files are
# deliberately retained as audit inputs rather than used for the main figures;
# their row counts are scope-specific and they may repeat display Reactome IDs.
source_inventory <- tibble(
  scope_id = scope_id,
  scope_label = scope_label,
  key = names(data_paths),
  file = unname(data_files),
  path = unname(data_paths),
  size_mb = as.numeric(file.info(unname(data_paths))$size) / 1024^2,
  use_in_figures = c(
    "Member-to-display mapping provenance",
    "Primary cell-specific evidence",
    "Primary final term list",
    "Granular-selection audit",
    "Parent-pruning audit and graph",
    "Pre-final audit only; not used as primary heatmap",
    "Pre-final audit only; not used as primary heatmap",
    "Pre-final audit only; not used as primary heatmap",
    "Term-level metric provenance",
    "Reduction summary"
  )
)

write_csv(source_inventory, file.path(output_dir, "source_file_inventory.csv"))

parameters <- tibble(
  parameter = c(
    "scope_id", "scope_label", "run_mode",
    "input_dir", "output_dir", "max_overview_terms", "max_graph_terms",
    "graph_parent_hops", "max_graph_nodes", "max_dechaining_parents",
    "max_descendants_per_parent",
    "top_terms_per_cell", "supplement_terms_per_page", "Reactome_hierarchy_source",
    "Reactome_DAG_created", "symmetric_NES_colour_limit"
  ),
  value = c(
    scope_id, scope_label, run_mode,
    input_dir, output_dir, max_overview_terms, max_graph_terms,
    graph_parent_hops, max_graph_nodes, max_dechaining_parents,
    max_descendants_per_parent,
    top_terms_per_cell, supplement_terms_per_page,
    if_else(is.na(reactome_hierarchy_source), "not installed", reactome_hierarchy_source),
    as.character(reactome_dag_created), nes_colour_limit
  )
)

write_csv(parameters, file.path(output_dir, "visualisation_parameters.csv"))
write_csv(figure_manifest, file.path(output_dir, "figure_manifest.csv"))

message("Created visualisations in: ", output_dir)
message("Main overview terms: ", nrow(overview_terms))
message("Final cell-term evidence rows plotted: ", nrow(final_cell_plot))
message("Reactome DAG created: ", reactome_dag_created)
