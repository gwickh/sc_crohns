#!/usr/bin/env Rscript

# Adaptive, hierarchy-aware summarisation of Reactome preranked GSEA results.
#
# This script performs no enrichment testing. It uses already-computed GSEA
# results, calculates term-level evidence and redundancy metrics, and chooses
# an adaptive display level:
#   * evidence-supported, distinct terms retain selected granular detail;
#   * redundant or less informative terms are mapped to a verified Reactome ancestor
#     where possible, using the same criteria for every Reactome branch;
#   * a structural-similarity representative is explicitly labelled as a group
#     representative rather than a parent when no useful common ancestor exists.
# The primary output in each analysis scope is a unique Reactome-ID list. Granular
# representatives are capped once per semantic cluster/direction, and strict
# Reactome pathway ancestors are removed when a retained descendant is available.
# By default, this script runs the computation once for the whole atlas and
# once for each predefined biological compartment. Each run is independent:
# candidate terms, evidence metrics, semantic groups, roll-ups and hierarchy
# pruning are all recalculated after restricting the input cell populations.
#
# Usage (all arguments are optional):
#   Rscript 6_Reactome_hierarchy_traversal.R \
#     significant_results.csv \
#     all_results.csv \
#     output_root_directory \
#     reactome_id_mapping.csv \
#     reactome_resources_directory
#
# Use the literal string "NA" to omit an optional file argument.
# Argument 5 contains offline ReactomePathways.txt and ReactomePathwaysRelation.txt.
# Arguments 6 and 7 are internal worker arguments used by the batch driver.

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

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

project_dir <- path.expand(
  "~/project-area/data/crohns_scrnaseq/10c_14n_analysis/pseudobulk_dea"
)

config <- list(
  significant_file = file.path(
    project_dir, "pathway_enrichment", "reactome_gsea_significant.csv"
  ),
  all_results_file = file.path(
    project_dir, "pathway_enrichment", "reactome_gsea_all_results.csv"
  ),
  output_dir = file.path(
    project_dir, "pathway_enrichment", "reactome_adaptive_summary_by_compartment_uniform"
  ),

  # Preferred optional mapping columns:
  #   pathway (or gs_name), reactome_id (or gs_exact_source), reactome_term (optional)
  # Without an explicit map, stable IDs and Reactome pathway names are matched
  # against the supplied offline hierarchy, with an ambiguity/exclusion audit.
  reactome_id_mapping_file = NA_character_,
  resources_dir = file.path(project_dir, "pathway_enrichment", "resources", "reactome"),

  fdr_column = "q_value_within_collection_cell_type",
  fdr_threshold = 0.05,
  min_abs_nes = 1.0,
  min_leading_edge_genes = 3L,
  min_leading_edge_fraction = 0.05,

  # This is hierarchy-derived structural Lin similarity. It is NOT GO/GOSemSim
  # semantic similarity. Threshold sensitivity outputs are retained for review.
  semantic_method = "StructuralLin",
  semantic_threshold = 0.70,
  sensitivity_thresholds = c(0.60, 0.70, 0.80),
  max_terms_for_similarity = 4500L,

  min_shared_cell_types_for_correlation = 4L,
  min_shared_significant_cell_types_for_leading_edge = 2L,
  min_profile_correlation = 0.60,
  min_direction_concordance = 0.75,
  min_leading_edge_jaccard = 0.30,
  min_leading_edge_overlap = 0.50,

  # The final reporting list uses a stricter gate than statistical
  # significance. Terms that fail this gate are retained at a coarser level;
  # they are not discarded. These display thresholds apply to every Reactome branch.
  final_granular_fdr = 0.01,
  final_granular_min_abs_nes = 1.50,
  final_granular_min_depth = 4,
  final_max_granular_per_cluster = 1L,

  # Clusters are rolled up only to a genuine common ancestor.
  # Depth is descriptive because Reactome is a DAG; information content is also used.
  min_rollup_depth = 2,
  target_rollup_depth = 3.5,
  target_rollup_information_content = 4.0,

  minimum_mapping_fraction = 0.90,
  max_plot_terms = 50L,
  cache_similarity_matrix = TRUE
)

args <- commandArgs(trailingOnly = TRUE)
arg_or_default <- function(position, default) {
  if (length(args) < position || identical(args[[position]], "")) {
    return(default)
  }
  if (toupper(args[[position]]) == "NA") {
    return(NA_character_)
  }
  path.expand(args[[position]])
}

config$significant_file <- arg_or_default(1, config$significant_file)
config$all_results_file <- arg_or_default(2, config$all_results_file)
config$output_dir <- arg_or_default(3, config$output_dir)
config$reactome_id_mapping_file <- arg_or_default(4, config$reactome_id_mapping_file)
config$resources_dir <- arg_or_default(5, config$resources_dir)
if (is.na(config$resources_dir)) config$resources_dir <- file.path(
  project_dir, "pathway_enrichment", "resources", "reactome"
)

# -----------------------------------------------------------------------------
# Analysis scopes and batch driver
# -----------------------------------------------------------------------------

# These labels are matched after Unicode transliteration, whitespace cleanup
# and case folding. The spelling found in the input is retained in every
# scientific output. In particular, CD8 T mem is not treated as an alias for a
# naive CD8 T-cell population.
scope_definitions <- list(
  whole_atlas = list(
    label = "Whole atlas",
    type = "atlas",
    cell_types = character()
  ),
  t_lymphocytes = list(
    label = "T-lymphocytes",
    type = "compartment",
    cell_types = c(
      "CD4 T cell", "CD8 T cell", "Naïve CD8 T cell", "Tfh",
      "NK T cell", "NK cell"
    )
  ),
  b_lymphocytes = list(
    label = "B-lymphocytes",
    type = "compartment",
    cell_types = c("B cells_cl_8", "B cells_cl_20", "B cells_cl_21")
  ),
  plasma_cells = list(
    label = "Plasma cells",
    type = "compartment",
    cell_types = c("Plasma cells_cl_4", "Plasma cells_cl_16")
  ),
  endothelial_pericytes = list(
    label = "Endothelial cells & related pericytes",
    type = "compartment",
    cell_types = c("Endothelial", "Pericyte")
  ),
  myofibroblasts_stromal = list(
    label = "Myofibroblasts & stromal cells",
    type = "compartment",
    cell_types = c(
      "Myofibroblast", "Mesenchymal", "ADAMDEC1+ stromal_cl_17",
      "ADAMDEC1+ stromal_cl_19", "NPY+ Stromal",
      "Transitional Stromal C3+"
    )
  ),
  macrophages_monocytes = list(
    label = "Macrophages & monocytes",
    type = "compartment",
    cell_types = c(
      "Macrophages_cl_34", "Macrophages_cl_14", "Macrophages_cl_2",
      "Macrophages_cl_3", "Mast cell", "Monocytes"
    )
  )
)

normalise_cell_type_key <- function(x) {
  x <- enc2utf8(as.character(x))
  x <- gsub("\u00a0", " ", x, fixed = TRUE)
  x <- gsub(
    "[\u2018\u2019\u201a\u201b\u201c\u201d\u201e\u201f]", "", x,
    perl = TRUE
  )
  transliterated <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "")
  transliterated[is.na(transliterated)] <- x[is.na(transliterated)]
  transliterated <- gsub("[[:space:]]+", " ", transliterated)
  tolower(trimws(transliterated))
}

collapse_manifest_values <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  if (length(x) == 0L) "" else paste(x, collapse = ";")
}

read_cell_type_column <- function(path, label) {
  if (is.na(path) || !file.exists(path)) {
    if (identical(label, "complete")) return(character())
    stop(label, " input does not exist: ", path)
  }
  x <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8"
  )
  if (!"cell_type" %in% names(x)) {
    stop(label, " input lacks a cell_type column: ", path)
  }
  as.character(x$cell_type)
}

resolve_script_path <- function() {
  file_argument <- grep(
    "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
  )
  if (length(file_argument) == 0L) {
    stop("Batch mode must be run with Rscript, not source()")
  }
  normalizePath(sub("^--file=", "", file_argument[[1]]), mustWork = TRUE)
}

analysis_scope_id <- if (length(args) >= 6L && nzchar(args[[6]])) {
  args[[6]]
} else {
  "all_scopes"
}

shared_similarity_cache_file <- if (
  length(args) >= 7L && nzchar(args[[7]]) && toupper(args[[7]]) != "NA"
) {
  path.expand(args[[7]])
} else {
  NA_character_
}

if (identical(analysis_scope_id, "all_scopes")) {
  output_root <- config$output_dir
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  scopes_root <- file.path(output_root, "scopes")
  shared_cache_dir <- file.path(output_root, "_shared_cache")
  dir.create(scopes_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(shared_cache_dir, recursive = TRUE, showWarnings = FALSE)

  significant_cell_types <- read_cell_type_column(
    config$significant_file, "significant"
  )
  complete_cell_types <- read_cell_type_column(
    config$all_results_file, "complete"
  )
  available_cell_types <- sort(unique(c(
    significant_cell_types, complete_cell_types
  )))
  available_cell_types <- available_cell_types[
    !is.na(available_cell_types) & nzchar(trimws(available_cell_types))
  ]

  available_keys <- normalise_cell_type_key(available_cell_types)
  collision_table <- split(available_cell_types, available_keys)
  collision_table <- collision_table[
    vapply(collision_table, function(z) length(unique(z)) > 1L, logical(1))
  ]
  if (length(collision_table) > 0L) {
    stop(
      "Cell-type normalisation is ambiguous for: ",
      paste(names(collision_table), collapse = ", ")
    )
  }
  available_lookup <- stats::setNames(available_cell_types, available_keys)

  membership_rows <- list()
  run_rows <- list()
  for (scope_id in names(scope_definitions)) {
    spec <- scope_definitions[[scope_id]]
    requested <- if (identical(scope_id, "whole_atlas")) {
      available_cell_types
    } else {
      spec$cell_types
    }
    requested_keys <- normalise_cell_type_key(requested)
    matched <- unname(available_lookup[requested_keys])
    present <- !is.na(matched)

    membership_rows[[scope_id]] <- data.frame(
      scope_id = scope_id,
      scope_label = spec$label,
      scope_type = spec$type,
      requested_cell_type = requested,
      requested_match_key = requested_keys,
      matched_input_cell_type = ifelse(present, matched, NA_character_),
      present_in_significant = ifelse(
        present, matched %in% significant_cell_types, FALSE
      ),
      present_in_complete_results = ifelse(
        present, matched %in% complete_cell_types, FALSE
      ),
      match_status = ifelse(present, "matched", "missing"),
      stringsAsFactors = FALSE
    )

    matched_values <- unique(matched[present])
    missing_values <- requested[!present]
    n_significant_scope_rows <- sum(significant_cell_types %in% matched_values)
    can_run <- length(matched_values) > 0L && n_significant_scope_rows > 0L
    run_rows[[scope_id]] <- data.frame(
      scope_id = scope_id,
      scope_label = spec$label,
      scope_type = spec$type,
      status = if (can_run) "pending" else "skipped",
      scope_completeness = if (
        length(missing_values) == 0L
      ) "complete" else "partial",
      status_reason = if (can_run) {
        "awaiting worker run"
      } else if (length(matched_values) == 0L) {
        "no requested cell populations were found in the inputs"
      } else {
        "no significant input rows were available for the matched populations"
      },
      output_dir = file.path(scopes_root, scope_id),
      n_requested_cell_types = length(requested),
      n_matched_cell_types = length(matched_values),
      n_missing_cell_types = length(missing_values),
      n_significant_input_rows = n_significant_scope_rows,
      matched_cell_types = collapse_manifest_values(matched_values),
      missing_cell_types = collapse_manifest_values(missing_values),
      exit_status = NA_integer_,
      final_unique_terms = NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  scope_membership <- do.call(rbind, membership_rows)
  rownames(scope_membership) <- NULL
  scope_manifest <- do.call(rbind, run_rows)
  rownames(scope_manifest) <- NULL
  utils::write.csv(
    scope_membership,
    file.path(output_root, "scope_cell_type_manifest.csv"),
    row.names = FALSE,
    na = ""
  )
  utils::write.csv(
    scope_manifest,
    file.path(output_root, "scope_run_manifest.csv"),
    row.names = FALSE,
    na = ""
  )

  assigned_compartment_cells <- unique(unlist(lapply(
    scope_definitions[names(scope_definitions) != "whole_atlas"],
    function(z) z$cell_types
  )))
  assigned_keys <- normalise_cell_type_key(assigned_compartment_cells)
  unassigned <- available_cell_types[
    !available_keys %in% assigned_keys
  ]
  utils::write.csv(
    data.frame(cell_type = unassigned, stringsAsFactors = FALSE),
    file.path(output_root, "unassigned_atlas_cell_types.csv"),
    row.names = FALSE,
    na = ""
  )

  if (identical(Sys.getenv("REACTOME_ROLLUP_PLAN_ONLY"), "1")) {
    message(
      "Scope plan written without running workers because ",
      "REACTOME_ROLLUP_PLAN_ONLY=1: ", output_root
    )
    quit(save = "no", status = 0L)
  }

  script_path <- resolve_script_path()
  rscript_path <- file.path(R.home("bin"), "Rscript")
  shared_cache <- file.path(
    shared_cache_dir,
    paste0("semantic_similarity_", config$semantic_method, ".rds")
  )

  for (scope_id in names(scope_definitions)) {
    manifest_index <- match(scope_id, scope_manifest$scope_id)
    if (scope_manifest$status[[manifest_index]] == "skipped") next

    scope_output <- scope_manifest$output_dir[[manifest_index]]
    dir.create(scope_output, recursive = TRUE, showWarnings = FALSE)
    worker_arguments <- c(
      config$significant_file,
      ifelse(is.na(config$all_results_file), "NA", config$all_results_file),
      scope_output,
      ifelse(is.na(config$reactome_id_mapping_file), "NA", config$reactome_id_mapping_file),
      config$resources_dir,
      scope_id,
      shared_cache
    )
    log_file <- file.path(scope_output, "rollup_console.log")
    message("Running Reactome roll-up for ", scope_definitions[[scope_id]]$label, "...")
    exit_status <- system2(
      rscript_path,
      args = c(shQuote(script_path), shQuote(worker_arguments)),
      stdout = log_file,
      stderr = log_file
    )

    scope_manifest$exit_status[[manifest_index]] <- exit_status
    scope_manifest$status[[manifest_index]] <- if (
      identical(as.integer(exit_status), 0L)
    ) {
      if (scope_manifest$n_missing_cell_types[[manifest_index]] > 0L) {
        "completed_partial"
      } else {
        "completed"
      }
    } else {
      "failed"
    }
    scope_manifest$status_reason[[manifest_index]] <- if (
      identical(as.integer(exit_status), 0L)
    ) {
      if (scope_manifest$n_missing_cell_types[[manifest_index]] > 0L) {
        "completed; one or more requested populations were absent"
      } else {
        "completed successfully"
      }
    } else {
      "worker returned a non-zero exit status"
    }

    final_file <- file.path(
      scope_output, "adaptive_reactome_final_significant_terms.csv"
    )
    if (identical(as.integer(exit_status), 0L) && file.exists(final_file)) {
      scope_manifest$final_unique_terms[[manifest_index]] <- nrow(
        utils::read.csv(final_file, stringsAsFactors = FALSE)
      )
    }
    utils::write.csv(
      scope_manifest,
      file.path(output_root, "scope_run_manifest.csv"),
      row.names = FALSE,
      na = ""
    )
  }

  failed_scopes <- scope_manifest$scope_id[scope_manifest$status == "failed"]
  if (length(failed_scopes) > 0L) {
    stop(
      "One or more scope runs failed: ", paste(failed_scopes, collapse = ", "),
      ". See each scope's rollup_console.log."
    )
  }
  message(
    "All available Reactome roll-up scopes complete: ", output_root,
    "\nRun manifest: scope_run_manifest.csv",
    "\nCell-type audit: scope_cell_type_manifest.csv"
  )
  quit(save = "no", status = 0L)
}

if (!analysis_scope_id %in% names(scope_definitions)) {
  stop(
    "Unknown worker scope '", analysis_scope_id, "'. Expected one of: ",
    paste(names(scope_definitions), collapse = ", ")
  )
}

config$analysis_scope <- analysis_scope_id
config$analysis_scope_label <- scope_definitions[[analysis_scope_id]]$label
config$analysis_scope_type <- scope_definitions[[analysis_scope_id]]$type
config$shared_similarity_cache_file <- shared_similarity_cache_file

# The synthetic super-root is structural context, never a reporting pathway.
uninformative_parent_ids <- "REACTOME_ROOT"

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------

required_packages <- c(
  "readr", "dplyr", "tidyr", "purrr", "ggplot2", "igraph"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them in the offline environment before running this script."
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else min(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else median(x)
}

collapse_unique_text <- function(x, separator = ";") {
  x <- as.character(x)
  x <- trimws(x[!is.na(x) & nzchar(trimws(x))])
  x <- sort(unique(x))
  if (length(x) == 0) NA_character_ else paste(x, collapse = separator)
}

split_delimited_values <- function(x, separator = ";") {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (length(x) == 0) return(character())

  values <- unlist(
    strsplit(x, separator, fixed = TRUE),
    use.names = FALSE
  )
  values <- trimws(values)
  sort(unique(values[nzchar(values)]))
}

collapse_delimited_values <- function(x, separator = ";") {
  values <- split_delimited_values(x, separator = separator)
  if (length(values) == 0) NA_character_ else paste(values, collapse = separator)
}

rescale_01 <- function(x) {
  out <- rep(0, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  limits <- range(x[ok])
  if (diff(limits) == 0) {
    out[ok] <- 1
  } else {
    out[ok] <- (x[ok] - limits[[1]]) / diff(limits)
  }
  out
}

split_gene_string <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return(character())
  unique(trimws(Filter(nzchar, strsplit(x, ";", fixed = TRUE)[[1]])))
}

gene_jaccard <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  union_n <- length(union(a, b))
  if (union_n == 0) return(NA_real_)
  length(intersect(a, b)) / union_n
}

gene_overlap_coefficient <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  denominator <- min(length(a), length(b))
  if (denominator == 0) return(NA_real_)
  length(intersect(a, b)) / denominator
}

normalise_msigdb_name <- function(x) {
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "")
  x <- toupper(x)
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  paste0("REACTOME_", x)
}

extract_reactome_id <- function(x) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  hit <- !is.na(x) & grepl("R-HSA-[0-9]+", x)
  out[hit] <- sub(".*(R-HSA-[0-9]+).*", "\\1", x[hit])
  out
}

first_existing_column <- function(columns, candidates) {
  hit <- candidates[candidates %in% columns]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

one_or_na <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x) == 1) x else NA_character_
}

read_gsea_results <- function(path, label) {
  if (is.na(path) || !file.exists(path)) {
    stop(label, " file does not exist: ", path)
  }

  x <- readr::read_csv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    na = c("", "NA", "NaN")
  )

  required <- c(
    "pathway", "collection", "pval", "ES", "NES", "size", "cell_type",
    "case", "reference", "comparison", "leading_edge", config$fdr_column
  )
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }

  numeric_columns <- intersect(
    unique(c(
      "pval", "q_value_within_collection_cell_type", "log2err", "ES", "NES",
      "size", "n_ranked_genes", "n_mapped_genes",
      "n_duplicate_symbols_removed", "q_value_across_collections_cell_type",
      "q_value_global", config$fdr_column
    )),
    names(x)
  )

  x <- x |>
    mutate(
      across(all_of(numeric_columns), ~ suppressWarnings(as.numeric(.x))),
      leading_edge = ifelse(is.na(leading_edge), "", leading_edge),
      leading_edge_n = vapply(
        leading_edge,
        function(z) length(split_gene_string(z)),
        integer(1)
      ),
      leading_edge_fraction = ifelse(
        is.finite(size) & size > 0,
        leading_edge_n / size,
        NA_real_
      ),
      evidence_q = .data[[config$fdr_column]],
      evidence_pass_row = is.finite(evidence_q) &
        evidence_q <= config$fdr_threshold &
        is.finite(NES),
      granular_evidence_pass_row = evidence_pass_row &
        abs(NES) >= config$min_abs_nes &
        leading_edge_n >= config$min_leading_edge_genes &
        is.finite(leading_edge_fraction) &
        leading_edge_fraction >= config$min_leading_edge_fraction
    )

  identity_columns <- c(
    "pathway", "collection", "cell_type", "case", "reference", "comparison"
  )
  if (any(vapply(
    x[identity_columns],
    function(z) any(is.na(z) | trimws(as.character(z)) == ""),
    logical(1)
  ))) {
    stop(label, " contains missing values in required identity columns")
  }

  duplicate_keys <- x |>
    count(pathway, cell_type, case, reference, comparison, name = "n") |>
    filter(n > 1)

  if (nrow(duplicate_keys) > 0) {
    readr::write_csv(
      duplicate_keys,
      file.path(config$output_dir, paste0(label, "_duplicate_keys.csv"))
    )
    stop(label, " contains duplicate pathway/cell-type/comparison keys")
  }

  x
}

# -----------------------------------------------------------------------------
# Read inputs
# -----------------------------------------------------------------------------

message("Reading significant Reactome results...")
significant_input <- read_gsea_results(config$significant_file, "significant_results")

have_full_results <- !is.na(config$all_results_file) &&
  file.exists(config$all_results_file)

if (have_full_results) {
  message("Reading complete Reactome results for full cell-type profiles...")
  all_results <- read_gsea_results(config$all_results_file, "all_results")
  profile_source <- "complete_results"
} else {
  warning(
    "Complete results were not found. Cell-type profiles will contain only ",
    "significant rows; parent evidence and correlations may be incomplete."
  )
  all_results <- significant_input
  profile_source <- "significant_results_only"
}

if (any(is.na(significant_input$collection) |
        toupper(significant_input$collection) != "REACTOME")) {
  stop("The significant input must contain REACTOME results only")
}
if (any(is.na(all_results$collection) | toupper(all_results$collection) != "REACTOME")) {
  stop("The complete-results input must contain REACTOME results only")
}

# Resolve and apply the scope before selecting candidate terms. This is what
# makes recurrence, profile similarity, semantic representatives, roll-ups and
# the final hierarchy antichain specific to the current compartment.
available_cell_types <- sort(unique(c(
  as.character(significant_input$cell_type),
  as.character(all_results$cell_type)
)))
available_keys <- normalise_cell_type_key(available_cell_types)
collision_table <- split(available_cell_types, available_keys)
collision_table <- collision_table[
  vapply(collision_table, function(z) length(unique(z)) > 1L, logical(1))
]
if (length(collision_table) > 0L) {
  stop(
    "Cell-type normalisation is ambiguous for: ",
    paste(names(collision_table), collapse = ", ")
  )
}
available_lookup <- stats::setNames(available_cell_types, available_keys)

requested_cell_types <- if (identical(config$analysis_scope, "whole_atlas")) {
  available_cell_types
} else {
  scope_definitions[[config$analysis_scope]]$cell_types
}
requested_keys <- normalise_cell_type_key(requested_cell_types)
matched_cell_types <- unname(available_lookup[requested_keys])
matched_flags <- !is.na(matched_cell_types)
included_cell_types <- sort(unique(matched_cell_types[matched_flags]))
missing_cell_types <- requested_cell_types[!matched_flags]

if (length(included_cell_types) == 0L) {
  stop(
    "No requested cell populations were found for scope ",
    config$analysis_scope, " (", config$analysis_scope_label, ")"
  )
}

scope_membership_audit <- tibble(
  scope_id = config$analysis_scope,
  scope_label = config$analysis_scope_label,
  scope_type = config$analysis_scope_type,
  requested_cell_type = requested_cell_types,
  requested_match_key = requested_keys,
  matched_input_cell_type = ifelse(
    matched_flags, matched_cell_types, NA_character_
  ),
  present_in_significant = ifelse(
    matched_flags,
    matched_cell_types %in% unique(significant_input$cell_type),
    FALSE
  ),
  present_in_complete_results = ifelse(
    matched_flags,
    matched_cell_types %in% unique(all_results$cell_type),
    FALSE
  ),
  match_status = ifelse(matched_flags, "matched", "missing")
)
readr::write_csv(
  scope_membership_audit,
  file.path(config$output_dir, "scope_cell_type_audit.csv")
)

significant_input <- significant_input |>
  filter(normalise_cell_type_key(cell_type) %in%
           normalise_cell_type_key(included_cell_types))
all_results <- all_results |>
  filter(normalise_cell_type_key(cell_type) %in%
           normalise_cell_type_key(included_cell_types))

if (nrow(significant_input) == 0L) {
  stop(
    "The scope contains no significant input rows: ",
    config$analysis_scope, " (", config$analysis_scope_label, ")"
  )
}
if (nrow(all_results) == 0L) {
  stop(
    "The scope contains no complete-result rows: ",
    config$analysis_scope, " (", config$analysis_scope_label, ")"
  )
}

if (length(missing_cell_types) > 0L) {
  warning(
    "Scope ", config$analysis_scope_label,
    " is partial. Requested populations not present in the input: ",
    paste(missing_cell_types, collapse = ", ")
  )
}

# Retain the atlas-wide four-cell minimum in every scope. Correlation based on
# only two or three cell populations is unstable, so small compartments leave
# profile correlation unavailable and rely on semantic and leading-edge
# evidence instead of silently relaxing the original criterion.

config$requested_cell_types <- requested_cell_types
config$included_cell_types <- included_cell_types
config$missing_cell_types <- missing_cell_types
config$scope_completeness <- if (
  length(missing_cell_types) == 0L
) "complete" else "partial"

write_scope_metadata <- function(
  status,
  n_candidate_pathways = NA_integer_,
  n_final_unique_terms = NA_integer_,
  n_final_by_cell_rows = NA_integer_
) {
  metadata <- tibble(
    scope_id = config$analysis_scope,
    scope_label = config$analysis_scope_label,
    scope_type = config$analysis_scope_type,
    scope_completeness = config$scope_completeness,
    status = status,
    requested_cell_types = collapse_manifest_values(
      config$requested_cell_types
    ),
    included_cell_types = collapse_manifest_values(
      config$included_cell_types
    ),
    missing_cell_types = collapse_manifest_values(
      config$missing_cell_types
    ),
    n_requested_cell_types = length(config$requested_cell_types),
    n_included_cell_types = length(config$included_cell_types),
    n_missing_cell_types = length(config$missing_cell_types),
    n_significant_input_rows = nrow(significant_input),
    n_complete_result_rows = nrow(all_results),
    n_candidate_pathways = n_candidate_pathways,
    n_final_unique_terms = n_final_unique_terms,
    n_final_by_cell_rows = n_final_by_cell_rows,
    min_shared_cell_types_for_correlation =
      config$min_shared_cell_types_for_correlation
  )
  readr::write_csv(
    metadata,
    file.path(config$output_dir, "run_scope_metadata.csv")
  )
  invisible(metadata)
}

write_scope_metadata("running")

message(
  "Analysis scope: ", config$analysis_scope_label,
  " [", config$analysis_scope, "]",
  "\nIncluded cell populations (", length(included_cell_types), "): ",
  paste(included_cell_types, collapse = ", ")
)

contrast_table <- bind_rows(
  significant_input |>
    distinct(case, reference, comparison) |>
    mutate(input = "significant"),
  all_results |>
    distinct(case, reference, comparison) |>
    mutate(input = "complete")
) |>
  distinct(case, reference, comparison)

if (nrow(contrast_table) != 1) {
  stop(
    "This script expects one case/reference comparison per run, but found ",
    nrow(contrast_table),
    ". Split the input by comparison and run each comparison separately."
  )
}

# The chosen FDR criterion, rather than the filename, determines candidate terms.
candidate_pathways <- significant_input |>
  filter(
    is.finite(evidence_q),
    evidence_q <= config$fdr_threshold,
    is.finite(NES)
  ) |>
  distinct(pathway) |>
  pull(pathway)

if (length(candidate_pathways) == 0) {
  stop("No candidate pathways pass the configured evidence threshold")
}

missing_from_complete <- setdiff(candidate_pathways, unique(all_results$pathway))
if (length(missing_from_complete) > 0) {
  stop(
    length(missing_from_complete),
    " candidate pathways are absent from the complete-results file. ",
    "Use matching result files from the same run."
  )
}

# -----------------------------------------------------------------------------
# Map MSigDB pathway labels to Reactome IDs
# -----------------------------------------------------------------------------

message("Loading the offline human Reactome hierarchy...")
hierarchy <- load_reactome_hierarchy(config$resources_dir)
reactome_meta <- as_tibble(hierarchy$terms) |> filter(reactome_id != hierarchy$root_id)
saveRDS(hierarchy, file.path(config$output_dir, "reactome_hierarchy_snapshot.rds"))
readr::write_csv(as_tibble(hierarchy$edges), file.path(config$output_dir, "reactome_hierarchy_edges.csv"))

all_pathways <- tibble(pathway = sort(unique(all_results$pathway)))
embedded_mapping <- if ("reactome_id" %in% names(all_results)) {
  all_results |> transmute(pathway,
    input_stable_id = reactome_extract_id(reactome_id)) |> distinct()
} else {
  all_pathways |> mutate(input_stable_id = NA_character_)
}

if (!is.na(config$reactome_id_mapping_file) && file.exists(config$reactome_id_mapping_file)) {
  message("Using supplied pathway-to-Reactome-ID mapping...")
  supplied_map <- readr::read_csv(
    config$reactome_id_mapping_file,
    show_col_types = FALSE,
    progress = FALSE
  )

  pathway_column <- first_existing_column(
    names(supplied_map), c("pathway", "gs_name")
  )
  id_column <- first_existing_column(
    names(supplied_map), c("reactome_id", "gs_exact_source")
  )

  if (is.na(pathway_column) || is.na(id_column)) {
    stop(
      "Reactome mapping file must contain pathway/gs_name and reactome_id/gs_exact_source"
    )
  }

  raw_mapping <- supplied_map |>
    transmute(
      pathway = .data[[pathway_column]],
      reactome_id = extract_reactome_id(as.character(.data[[id_column]]))
    ) |>
    filter(!is.na(pathway)) |>
    distinct()

  mapping_method <- "supplied_reactome_id"
} else {
  message(
    "No explicit mapping supplied; matching stable IDs and offline Reactome pathway names..."
  )
  raw_mapping <- all_pathways |>
    mutate(pathway_key = reactome_normalise_name(pathway),
      embedded_id = reactome_extract_id(pathway)) |>
    left_join(reactome_meta |> select(pathway_key, reactome_id), by = "pathway_key") |>
    left_join(embedded_mapping, by = "pathway") |>
    mutate(reactome_id = coalesce(input_stable_id, embedded_id, reactome_id)) |>
    select(pathway, reactome_id)
  mapping_method <- "normalised_term_name"
}

pathway_mapping <- all_pathways |>
  left_join(raw_mapping, by = "pathway") |>
  group_by(pathway) |>
  summarise(
    n_reactome_matches = n_distinct(reactome_id[!is.na(reactome_id)]),
    matched_reactome_ids = paste(sort(unique(reactome_id[!is.na(reactome_id)])), collapse = ";"),
    reactome_id = one_or_na(reactome_id),
    .groups = "drop"
  ) |>
  mutate(
    provisional_mapping_status = case_when(
      n_reactome_matches == 0 ~ "unmapped",
      n_reactome_matches == 1 ~ "mapped",
      TRUE ~ "ambiguous"
    ),
    mapping_method = mapping_method
  ) |>
  left_join(reactome_meta |> select(reactome_id, reactome_term, ontology), by = "reactome_id") |>
  mutate(
    mapping_status = case_when(
      provisional_mapping_status == "unmapped" ~ "unmapped",
      provisional_mapping_status == "ambiguous" ~ "ambiguous",
      is.na(reactome_term) | is.na(ontology) | ontology != "Reactome" ~
        "invalid_or_non_human_pathway",
      TRUE ~ "mapped"
    )
  ) |>
  select(-provisional_mapping_status)

if (
  nrow(pathway_mapping) != nrow(all_pathways) ||
    anyDuplicated(pathway_mapping$pathway)
) {
  stop("Reactome-ID mapping must contain exactly one audit row per pathway")
}

readr::write_csv(
  pathway_mapping,
  file.path(config$output_dir, "reactome_id_mapping_qc.csv")
)
readr::write_csv(
  pathway_mapping |> filter(mapping_status != "mapped"),
  file.path(config$output_dir, "unmapped_or_ambiguous_terms.csv")
)

candidate_mapping <- pathway_mapping |>
  filter(pathway %in% candidate_pathways)
mapping_fraction <- mean(candidate_mapping$mapping_status == "mapped")
omitted_candidate_pathways <- sum(candidate_mapping$mapping_status != "mapped")

readr::write_csv(
  candidate_mapping |> filter(mapping_status != "mapped"),
  file.path(config$output_dir, "candidate_mapping_failures.csv")
)

message(
  sprintf(
    "Mapped %.1f%% of candidate pathways (%d/%d)",
    100 * mapping_fraction,
    sum(candidate_mapping$mapping_status == "mapped"),
    nrow(candidate_mapping)
  )
)

if (!is.finite(mapping_fraction) ||
    mapping_fraction < config$minimum_mapping_fraction) {
  stop(
    "Reactome-ID mapping coverage is below ",
    100 * config$minimum_mapping_fraction,
    "%. Inspect reactome_id_mapping_qc.csv and supply an exact mapping containing ",
    "gs_exact_source/Reactome IDs from the matching MSigDB release."
  )
}

if (omitted_candidate_pathways > 0) {
  warning(
    omitted_candidate_pathways,
    " significant candidate pathways were excluded because their Reactome-ID ",
    "mapping was not unique and current. See candidate_mapping_failures.csv."
  )
}

all_results <- all_results |>
  select(-any_of(c("reactome_id", "reactome_term", "mapping_status", "mapping_method")))
significant_input <- significant_input |>
  select(-any_of(c("reactome_id", "reactome_term", "mapping_status", "mapping_method")))

all_results_mapped <- all_results |>
  left_join(
    pathway_mapping |>
      select(pathway, reactome_id, reactome_term, mapping_status, mapping_method),
    by = "pathway"
  )

significant_mapped <- significant_input |>
  left_join(
    pathway_mapping |>
      select(pathway, reactome_id, reactome_term, mapping_status, mapping_method),
    by = "pathway"
  )

candidate_ids <- candidate_mapping |>
  filter(mapping_status == "mapped") |>
  pull(reactome_id) |>
  unique() |>
  sort()

if (length(candidate_ids) > config$max_terms_for_similarity) {
  stop(
    "There are ", length(candidate_ids),
    " candidate Reactome IDs, exceeding max_terms_for_similarity=",
    config$max_terms_for_similarity,
    ". Tighten the evidence filter or increase the limit after checking memory."
  )
}

# -----------------------------------------------------------------------------
# Build the Reactome pathway/subpathway DAG
# -----------------------------------------------------------------------------

message("Building the Reactome parent graph...")
parent_map <- hierarchy$parent_map
root_reactome_id <- hierarchy$root_id

ancestor_cache <- new.env(parent = emptyenv())
depth_cache <- new.env(parent = emptyenv())

get_ancestors <- function(reactome_id, trail = character()) {
  if (exists(reactome_id, envir = ancestor_cache, inherits = FALSE)) {
    return(get(reactome_id, envir = ancestor_cache, inherits = FALSE))
  }
  if (reactome_id %in% trail) {
    stop("Cycle detected in filtered Reactome graph at ", reactome_id)
  }

  parents <- parent_map[[reactome_id]]
  parents <- unique(parents[!is.na(parents)])
  if (length(parents) == 0) {
    answer <- character()
  } else {
    answer <- unique(c(
      parents,
      unlist(
        lapply(parents, get_ancestors, trail = c(trail, reactome_id)),
        use.names = FALSE
      )
    ))
  }
  assign(reactome_id, answer, envir = ancestor_cache)
  answer
}

# Return a strict Reactome antichain flag for a set of display terms. A term is
# marked for removal only when it is a verified ancestor of another term in
# the same biological context. Reactome depth alone is not used because Reactome is a DAG.
get_hierarchy_pruning_flags <- function(reactome_ids) {
  ids <- sort(unique(reactome_ids[!is.na(reactome_ids) & nzchar(reactome_ids)]))

  if (length(ids) == 0) {
    return(tibble(
      display_reactome_id = character(),
      hierarchy_parent_removed = logical(),
      retained_descendant_reactome_ids = character()
    ))
  }

  strict_descendants <- setNames(
    lapply(ids, function(parent_id) {
      setdiff(ids, parent_id)[
        vapply(
          setdiff(ids, parent_id),
          function(child_id) parent_id %in% get_ancestors(child_id),
          logical(1)
        )
      ]
    }),
    ids
  )

  hierarchy_parent_removed <- lengths(strict_descendants) > 0
  retained_ids <- ids[!hierarchy_parent_removed]

  retained_descendants <- lapply(ids, function(parent_id) {
    retained_ids[
      vapply(
        retained_ids,
        function(child_id) {
          !identical(parent_id, child_id) &
            parent_id %in% get_ancestors(child_id)
        },
        logical(1)
      )
    ]
  })

  tibble(
    display_reactome_id = ids,
    hierarchy_parent_removed = unname(hierarchy_parent_removed),
    retained_descendant_reactome_ids = vapply(
      retained_descendants,
      collapse_unique_text,
      character(1)
    )
  )
}

get_depth_range <- function(reactome_id, trail = character()) {
  if (exists(reactome_id, envir = depth_cache, inherits = FALSE)) {
    return(get(reactome_id, envir = depth_cache, inherits = FALSE))
  }
  if (identical(reactome_id, root_reactome_id)) {
    answer <- c(min_depth = 0, max_depth = 0)
    assign(reactome_id, answer, envir = depth_cache)
    return(answer)
  }
  if (reactome_id %in% trail) {
    stop("Cycle detected while calculating Reactome depth at ", reactome_id)
  }

  parents <- parent_map[[reactome_id]]
  parents <- unique(parents[!is.na(parents)])
  if (length(parents) == 0) {
    answer <- c(min_depth = NA_real_, max_depth = NA_real_)
  } else {
    parent_depths <- do.call(
      rbind,
      lapply(parents, get_depth_range, trail = c(trail, reactome_id))
    )
    valid <- is.finite(parent_depths[, "min_depth"]) &
      is.finite(parent_depths[, "max_depth"])
    if (!any(valid)) {
      answer <- c(min_depth = NA_real_, max_depth = NA_real_)
    } else {
      answer <- c(
        min_depth = 1 + min(parent_depths[valid, "min_depth"]),
        max_depth = 1 + max(parent_depths[valid, "max_depth"])
      )
    }
  }
  assign(reactome_id, answer, envir = depth_cache)
  answer
}

candidate_ancestors <- setNames(
  lapply(candidate_ids, get_ancestors),
  candidate_ids
)

# -----------------------------------------------------------------------------
# Information content and term-level metrics
# -----------------------------------------------------------------------------

message("Using hierarchy-derived structural information content...")
ic_vector <- hierarchy$ic

candidate_depths <- do.call(rbind, lapply(candidate_ids, get_depth_range))
candidate_depths <- tibble(
  reactome_id = candidate_ids,
  min_depth = candidate_depths[, "min_depth"],
  max_depth = candidate_depths[, "max_depth"]
)

evidence_rows <- all_results_mapped |>
  filter(
    reactome_id %in% candidate_ids,
    evidence_pass_row,
    evidence_q <= config$fdr_threshold
  )

# Use the complete table for profiles but the configured evidence threshold for
# counts and leading-edge summaries.
term_metrics <- all_results_mapped |>
  filter(reactome_id %in% candidate_ids) |>
  group_by(reactome_id) |>
  summarise(
    pathway = first(pathway),
    reactome_term = first(reactome_term),
    best_q = safe_min(evidence_q[evidence_pass_row]),
    max_abs_nes = safe_max(abs(NES[evidence_pass_row])),
    median_abs_nes = safe_median(abs(NES[evidence_pass_row])),
    n_tested_cell_types = n_distinct(cell_type[is.finite(NES)]),
    n_significant_cell_types = n_distinct(
      cell_type[evidence_pass_row & evidence_q <= config$fdr_threshold]
    ),
    n_positive_significant = n_distinct(
      cell_type[evidence_pass_row & evidence_q <= config$fdr_threshold & NES > 0]
    ),
    n_negative_significant = n_distinct(
      cell_type[evidence_pass_row & evidence_q <= config$fdr_threshold & NES < 0]
    ),
    n_granular_evidence_cell_types = n_distinct(
      cell_type[granular_evidence_pass_row]
    ),
    max_leading_edge_n = safe_max(
      leading_edge_n[evidence_pass_row & evidence_q <= config$fdr_threshold]
    ),
    max_leading_edge_fraction = safe_max(
      leading_edge_fraction[
        evidence_pass_row & evidence_q <= config$fdr_threshold
      ]
    ),
    .groups = "drop"
  ) |>
  mutate(
    dominant_direction = case_when(
      n_positive_significant > 0 & n_negative_significant > 0 ~ "mixed",
      n_positive_significant > 0 ~ "positive",
      n_negative_significant > 0 ~ "negative",
      TRUE ~ "none"
    ),
    direction_consistency = ifelse(
      n_positive_significant + n_negative_significant > 0,
      pmax(n_positive_significant, n_negative_significant) /
        (n_positive_significant + n_negative_significant),
      NA_real_
    ),
    evidence_pass = n_granular_evidence_cell_types > 0
  ) |>
  left_join(candidate_depths, by = "reactome_id") |>
  mutate(
    information_content = unname(ic_vector[reactome_id]),
    information_content = ifelse(
      is.finite(information_content), information_content, NA_real_
    )
  )

if (!setequal(candidate_ids, term_metrics$reactome_id)) {
  stop(
    "Candidate Reactome IDs and term-metric Reactome IDs differ; check that the two input ",
    "files and Reactome mapping come from matching releases."
  )
}

# Leading-edge union and recurrent genes across significant cell populations.
leading_edge_long <- evidence_rows |>
  select(reactome_id, cell_type, leading_edge) |>
  tidyr::separate_rows(leading_edge, sep = ";") |>
  mutate(leading_edge = trimws(leading_edge)) |>
  filter(!is.na(leading_edge), leading_edge != "") |>
  distinct(reactome_id, cell_type, leading_edge)

leading_edge_metrics <- leading_edge_long |>
  count(reactome_id, leading_edge, name = "n_cell_types") |>
  group_by(reactome_id) |>
  summarise(
    leading_edge_union = paste(sort(unique(leading_edge)), collapse = ";"),
    n_leading_edge_union = n_distinct(leading_edge),
    recurrent_leading_edge = paste(
      sort(unique(leading_edge[n_cell_types >= 2])), collapse = ";"
    ),
    n_recurrent_leading_edge = n_distinct(leading_edge[n_cell_types >= 2]),
    .groups = "drop"
  )

term_metrics <- term_metrics |>
  left_join(leading_edge_metrics, by = "reactome_id") |>
  mutate(
    leading_edge_union = ifelse(
      is.na(leading_edge_union), "", leading_edge_union
    ),
    recurrent_leading_edge = ifelse(
      is.na(recurrent_leading_edge), "", recurrent_leading_edge
    ),
    n_leading_edge_union = ifelse(
      is.na(n_leading_edge_union), 0L, n_leading_edge_union
    ),
    n_recurrent_leading_edge = ifelse(
      is.na(n_recurrent_leading_edge), 0L, n_recurrent_leading_edge
    )
  )

# -----------------------------------------------------------------------------
# Evidence-based representative scoring
# -----------------------------------------------------------------------------

# Presentation score chooses representatives. It is not a statistical test.
# Every Reactome branch uses the same evidence, recurrence and information-content
# weights; no named branch or biological topic receives a score adjustment.
term_metrics <- term_metrics |>
  mutate(
    capped_log_q = pmin(-log10(pmax(best_q, 1e-300)), 20),
    presentation_score =
      0.40 * rescale_01(capped_log_q) +
      0.25 * rescale_01(pmin(max_abs_nes, 3)) +
      0.15 * rescale_01(log1p(n_significant_cell_types)) +
      0.20 * rescale_01(information_content)
  )

# -----------------------------------------------------------------------------
# Semantic similarity and clusters
# -----------------------------------------------------------------------------

message(
  "Calculating structural similarity for ", length(candidate_ids),
  " candidate Reactome terms..."
)

local_similarity_cache_file <- file.path(
  config$output_dir,
  paste0("semantic_similarity_", config$semantic_method, ".rds")
)

cache_candidates <- unique(c(
  if (!is.na(config$shared_similarity_cache_file)) {
    config$shared_similarity_cache_file
  } else {
    character()
  },
  local_similarity_cache_file
))

cache_is_valid <- FALSE
cache_source_file <- NA_character_
if (config$cache_similarity_matrix) {
  for (candidate_cache_file in cache_candidates) {
    if (!file.exists(candidate_cache_file)) next
    cached <- tryCatch(
      readRDS(candidate_cache_file),
      error = function(e) NULL
    )
    cached_ids <- if (is.list(cached)) as.character(cached$ids) else character()
    cached_matrix <- if (is.list(cached)) cached$similarity else NULL
    cache_is_valid <- is.list(cached) &&
      identical(cached$method, config$semantic_method) &&
      identical(cached$hierarchy_signature, hierarchy$signature) &&
      is.matrix(cached_matrix) &&
      all(candidate_ids %in% cached_ids) &&
      all(candidate_ids %in% rownames(cached_matrix)) &&
      all(candidate_ids %in% colnames(cached_matrix))
    if (cache_is_valid) {
      semantic_similarity <- cached_matrix[
        candidate_ids, candidate_ids, drop = FALSE
      ]
      cache_source_file <- candidate_cache_file
      message(
        "Using cached structural-similarity matrix (",
        length(cached_ids), " cached Reactome IDs; ",
        length(candidate_ids), " used in this scope): ",
        candidate_cache_file
      )
      break
    }
  }
}

if (!cache_is_valid) {
  if (length(candidate_ids) == 1L) {
    semantic_similarity <- matrix(
      1,
      nrow = 1,
      ncol = 1,
      dimnames = list(candidate_ids, candidate_ids)
    )
  } else {
    semantic_similarity <- reactome_structural_similarity(hierarchy, candidate_ids)
  }
  semantic_similarity[!is.finite(semantic_similarity)] <- 0
  semantic_similarity <- (semantic_similarity + t(semantic_similarity)) / 2
  semantic_similarity <- pmin(pmax(semantic_similarity, 0), 1)
  diag(semantic_similarity) <- 1

  if (config$cache_similarity_matrix) {
    # The whole-atlas matrix is a reusable superset. A compartment that has to
    # recompute writes only its local cache and never replaces that superset.
    similarity_cache_file <- if (
      identical(config$analysis_scope, "whole_atlas") &&
        !is.na(config$shared_similarity_cache_file)
    ) {
      config$shared_similarity_cache_file
    } else {
      local_similarity_cache_file
    }
    dir.create(
      dirname(similarity_cache_file), recursive = TRUE, showWarnings = FALSE
    )
    saveRDS(
      list(
        ids = candidate_ids,
        method = config$semantic_method,
        hierarchy_signature = hierarchy$signature,
        hierarchy_source = hierarchy$source,
        similarity = semantic_similarity
      ),
      similarity_cache_file,
      compress = FALSE
    )
    cache_source_file <- similarity_cache_file
  }
}

if (is.null(dim(semantic_similarity)) ||
    nrow(semantic_similarity) != length(candidate_ids) ||
    ncol(semantic_similarity) != length(candidate_ids)) {
  stop("Structural-similarity matrix dimensions do not match candidate Reactome IDs")
}

if (is.null(rownames(semantic_similarity)) ||
    is.null(colnames(semantic_similarity))) {
  rownames(semantic_similarity) <- candidate_ids
  colnames(semantic_similarity) <- candidate_ids
} else {
  if (!setequal(rownames(semantic_similarity), candidate_ids) ||
      !setequal(colnames(semantic_similarity), candidate_ids)) {
    stop("Structural-similarity matrix names do not match candidate Reactome IDs")
  }
  semantic_similarity <- semantic_similarity[candidate_ids, candidate_ids, drop = FALSE]
}

semantic_similarity[!is.finite(semantic_similarity)] <- 0
semantic_similarity <- (semantic_similarity + t(semantic_similarity)) / 2
semantic_similarity <- pmin(pmax(semantic_similarity, 0), 1)
diag(semantic_similarity) <- 1

representative_scores <- setNames(
  term_metrics$presentation_score,
  term_metrics$reactome_id
)
representative_scores[!is.finite(representative_scores)] <- 0

reduce_at_threshold <- function(threshold) {
  if (length(candidate_ids) == 1L) {
    return(tibble(
      reactome_id = candidate_ids,
      semantic_representative = candidate_ids
    ))
  }

  clustering <- stats::hclust(stats::as.dist(1 - semantic_similarity), method = "complete")
  clusters <- stats::cutree(clustering, h = 1 - threshold)
  members <- split(candidate_ids, clusters[candidate_ids])
  representatives <- lapply(members, function(ids) ids[order(-representative_scores[ids], ids)][[1]])
  lookup <- character()
  for (group in names(members)) {
    lookup <- c(lookup, setNames(rep(representatives[[group]], length(members[[group]])), members[[group]]))
  }
  tibble(reactome_id = candidate_ids,
    semantic_representative = unname(lookup[candidate_ids]))
}

threshold_labels <- gsub("\\.", "_", as.character(config$sensitivity_thresholds))
cluster_tables <- map(
  config$sensitivity_thresholds,
  reduce_at_threshold
)

semantic_clusters <- term_metrics |> select(reactome_id)
for (i in seq_along(cluster_tables)) {
  column_name <- paste0("representative_t", threshold_labels[[i]])
  semantic_clusters <- semantic_clusters |>
    left_join(
      cluster_tables[[i]] |>
        rename(!!column_name := semantic_representative),
      by = "reactome_id"
    )
}

main_threshold_label <- paste0(
  "representative_t", gsub("\\.", "_", as.character(config$semantic_threshold))
)
if (!main_threshold_label %in% names(semantic_clusters)) {
  stop("semantic_threshold must be included in sensitivity_thresholds")
}

representative_columns <- grep(
  "^representative_t", names(semantic_clusters), value = TRUE
)
semantic_clusters$representative_stable <- apply(
  semantic_clusters[, representative_columns, drop = FALSE],
  1,
  function(z) {
    z <- unique(z[!is.na(z)])
    length(z) == 1
  }
)
semantic_clusters$semantic_representative <-
  semantic_clusters[[main_threshold_label]]

if (any(is.na(semantic_clusters$semantic_representative))) {
  stop("Structural clustering returned missing representatives for one or more candidate terms")
}

# Complete-linkage structural similarity defines the clusters. Within each cluster, use a fixed
# direction-specific member as the comparison/display representative so that a
# positive-only term is never assessed against a negative-only representative.
semantic_clusters <- semantic_clusters |>
  rename(structural_cluster_representative = semantic_representative) |>
  left_join(
    term_metrics |>
      select(
        reactome_id, cluster_direction = dominant_direction, presentation_score
      ),
    by = "reactome_id"
  ) |>
  group_by(structural_cluster_representative, cluster_direction) |>
  arrange(desc(presentation_score), reactome_id, .by_group = TRUE) |>
  mutate(semantic_representative = first(reactome_id)) |>
  ungroup() |>
  select(-cluster_direction, -presentation_score)

# -----------------------------------------------------------------------------
# NES-profile and leading-edge similarity to the semantic representative
# -----------------------------------------------------------------------------

cell_types <- sort(unique(all_results_mapped$cell_type))

profile_rows <- all_results_mapped |>
  filter(reactome_id %in% candidate_ids) |>
  arrange(reactome_id, cell_type, is.na(evidence_q), evidence_q) |>
  group_by(reactome_id, cell_type) |>
  slice(1) |>
  ungroup() |>
  select(reactome_id, cell_type, NES)

nes_matrix <- matrix(
  NA_real_,
  nrow = length(candidate_ids),
  ncol = length(cell_types),
  dimnames = list(candidate_ids, cell_types)
)

row_index <- match(profile_rows$reactome_id, rownames(nes_matrix))
column_index <- match(profile_rows$cell_type, colnames(nes_matrix))
nes_matrix[cbind(row_index, column_index)] <- profile_rows$NES

evidence_rows_sorted <- evidence_rows |> arrange(reactome_id, cell_type)
evidence_by_term <- split(evidence_rows_sorted, evidence_rows_sorted$reactome_id)

leading_edge_by_term_cell <- lapply(evidence_by_term, function(df) {
  by_cell <- split(df$leading_edge, df$cell_type)
  lapply(by_cell, function(strings) {
    unique(unlist(lapply(strings, split_gene_string), use.names = FALSE))
  })
})

pair_profile_metrics <- function(reactome_id, representative) {
  if (is.na(representative) ||
      !reactome_id %in% rownames(nes_matrix) ||
      !representative %in% rownames(nes_matrix)) {
    return(c(
      shared_cell_types = 0,
      profile_correlation = NA_real_,
      direction_concordance = NA_real_
    ))
  }

  a <- nes_matrix[reactome_id, ]
  b <- nes_matrix[representative, ]
  ok <- is.finite(a) & is.finite(b)
  n_shared <- sum(ok)
  correlation <- NA_real_
  if (n_shared >= config$min_shared_cell_types_for_correlation &&
      stats::sd(a[ok]) > 0 && stats::sd(b[ok]) > 0) {
    correlation <- suppressWarnings(
      stats::cor(a[ok], b[ok], method = "spearman")
    )
  }

  direction_ok <- ok & a != 0 & b != 0
  direction_concordance <- if (any(direction_ok)) {
    mean(sign(a[direction_ok]) == sign(b[direction_ok]))
  } else {
    NA_real_
  }

  c(
    shared_cell_types = n_shared,
    profile_correlation = correlation,
    direction_concordance = direction_concordance
  )
}

pair_leading_edge_metrics <- function(reactome_id, representative) {
  a <- leading_edge_by_term_cell[[reactome_id]]
  b <- leading_edge_by_term_cell[[representative]]
  if (is.null(a) || is.null(b)) {
    return(c(
      shared_significant_cell_types = 0,
      median_leading_edge_jaccard = NA_real_,
      median_leading_edge_overlap = NA_real_,
      union_leading_edge_jaccard = NA_real_
    ))
  }

  shared <- intersect(names(a), names(b))
  per_cell_jaccard <- vapply(
    shared,
    function(cell) gene_jaccard(a[[cell]], b[[cell]]),
    numeric(1)
  )
  per_cell_overlap <- vapply(
    shared,
    function(cell) gene_overlap_coefficient(a[[cell]], b[[cell]]),
    numeric(1)
  )

  c(
    shared_significant_cell_types = length(shared),
    median_leading_edge_jaccard = safe_median(per_cell_jaccard),
    median_leading_edge_overlap = safe_median(per_cell_overlap),
    union_leading_edge_jaccard = gene_jaccard(
      unique(unlist(a, use.names = FALSE)),
      unique(unlist(b, use.names = FALSE))
    )
  )
}

semantic_clusters <- semantic_clusters |>
  mutate(
    semantic_similarity_to_representative = map2_dbl(
      reactome_id,
      semantic_representative,
      function(id, parent) {
        if (is.na(parent) ||
            !id %in% rownames(semantic_similarity) ||
            !parent %in% colnames(semantic_similarity)) {
          return(NA_real_)
        }
        semantic_similarity[id, parent]
      }
    )
  )

profile_pair_table <- map2_dfr(
  semantic_clusters$reactome_id,
  semantic_clusters$semantic_representative,
  function(id, representative) {
    values <- pair_profile_metrics(id, representative)
    tibble(
      reactome_id = id,
      shared_cell_types = unname(values[["shared_cell_types"]]),
      profile_correlation = unname(values[["profile_correlation"]]),
      direction_concordance_to_representative = unname(
        values[["direction_concordance"]]
      )
    )
  }
)

leading_edge_pair_table <- map2_dfr(
  semantic_clusters$reactome_id,
  semantic_clusters$semantic_representative,
  function(id, representative) {
    values <- pair_leading_edge_metrics(id, representative)
    tibble(
      reactome_id = id,
      shared_significant_cell_types = unname(
        values[["shared_significant_cell_types"]]
      ),
      median_leading_edge_jaccard = unname(
        values[["median_leading_edge_jaccard"]]
      ),
      median_leading_edge_overlap = unname(
        values[["median_leading_edge_overlap"]]
      ),
      union_leading_edge_jaccard = unname(
        values[["union_leading_edge_jaccard"]]
      )
    )
  }
)

term_metrics <- term_metrics |>
  left_join(semantic_clusters, by = "reactome_id") |>
  left_join(profile_pair_table, by = "reactome_id") |>
  left_join(leading_edge_pair_table, by = "reactome_id") |>
  mutate(
    profile_similarity_available =
      shared_cell_types >= config$min_shared_cell_types_for_correlation &
      is.finite(profile_correlation) &
      is.finite(direction_concordance_to_representative),
    leading_edge_similarity_available =
      shared_significant_cell_types >=
        config$min_shared_significant_cell_types_for_leading_edge &
      (
        is.finite(median_leading_edge_jaccard) |
        is.finite(median_leading_edge_overlap)
      ),
    profile_concordant = profile_similarity_available &
      profile_correlation >= config$min_profile_correlation &
      direction_concordance_to_representative >=
        config$min_direction_concordance,
    leading_edge_concordant = leading_edge_similarity_available &
      (
        coalesce(
          median_leading_edge_jaccard >= config$min_leading_edge_jaccard,
          FALSE
        ) |
        coalesce(
          median_leading_edge_overlap >= config$min_leading_edge_overlap,
          FALSE
        )
      ),
    redundancy_supported =
      coalesce(
        semantic_similarity_to_representative >= config$semantic_threshold,
        FALSE
      ) &
      (coalesce(profile_concordant, FALSE) |
         coalesce(leading_edge_concordant, FALSE)),
    similarity_uncertain =
      !profile_similarity_available & !leading_edge_similarity_available,
    signal_distinct_from_representative =
      (profile_similarity_available & !profile_concordant) |
      (leading_edge_similarity_available & !leading_edge_concordant),
    direction_stratum = dominant_direction
  ) |>
  mutate(
    keep_granular = evidence_pass &
      (
        reactome_id == semantic_representative |
        signal_distinct_from_representative |
        similarity_uncertain |
        !redundancy_supported
      )
  ) |>
  ungroup()

# -----------------------------------------------------------------------------
# Choose true common ancestors for semantic clusters
# -----------------------------------------------------------------------------

choose_common_ancestor <- function(member_ids) {
  member_ids <- unique(member_ids)
  if (length(member_ids) == 0L) return(NA_character_)
  lineages <- lapply(member_ids, get_ancestors)
  common <- Reduce(intersect, lineages)
  common <- setdiff(common, uninformative_parent_ids)
  if (length(common) == 0) return(NA_character_)

  candidate_table <- map_dfr(common, function(id) {
    depth <- get_depth_range(id)
    tibble(
      reactome_id = id,
      min_depth = depth[["min_depth"]],
      max_depth = depth[["max_depth"]],
      depth_midpoint = mean(depth, na.rm = TRUE),
      information_content = unname(ic_vector[id])
    )
  }) |>
    filter(
      is.finite(depth_midpoint),
      min_depth >= config$min_rollup_depth
    ) |>
    mutate(
      information_content = ifelse(
        is.finite(information_content), information_content, NA_real_
      ),
      depth_distance = abs(depth_midpoint - config$target_rollup_depth),
      ic_distance = abs(
        information_content - config$target_rollup_information_content
      ),
      ancestor_selection_score = depth_distance + 0.25 * ifelse(
        is.finite(ic_distance), ic_distance, 4
      )
    ) |>
    arrange(
      ancestor_selection_score, depth_distance, ic_distance,
      desc(information_content), reactome_id
    )

  if (nrow(candidate_table) == 0) NA_character_ else candidate_table$reactome_id[[1]]
}

cluster_ancestors <- term_metrics |>
  filter(!keep_granular) |>
  group_by(semantic_representative, direction_stratum) |>
  summarise(
    cluster_common_ancestor = choose_common_ancestor(reactome_id),
    n_cluster_terms = n_distinct(reactome_id),
    .groups = "drop"
  )

# A separate fallback is calculated across every semantic cluster.
# It does not alter the audit/display mapping above; it is used only when a
# term lacks sufficient evidence to be shown granularly in a particular cell.
fallback_cluster_ancestors <- term_metrics |>
  group_by(semantic_representative, direction_stratum) |>
  summarise(
    fallback_common_ancestor = choose_common_ancestor(reactome_id),
    .groups = "drop"
  )

# If a cluster has no useful strict common ancestor, use an informative strict
# ancestor of the individual term before falling back to a semantic label.
individual_rollup_ancestors <- tibble(
  reactome_id = term_metrics$reactome_id,
  individual_rollup_ancestor = vapply(
    term_metrics$reactome_id,
    function(id) choose_common_ancestor(id),
    character(1)
  )
)

term_metrics <- term_metrics |>
  left_join(
    cluster_ancestors,
    by = c("semantic_representative", "direction_stratum")
  ) |>
  left_join(
    fallback_cluster_ancestors,
    by = c("semantic_representative", "direction_stratum")
  ) |>
  left_join(
    individual_rollup_ancestors,
    by = "reactome_id"
  )

display_mapping <- term_metrics |>
  mutate(
    display_reactome_id = case_when(
      keep_granular ~ reactome_id,
      !is.na(cluster_common_ancestor) ~ cluster_common_ancestor,
      TRUE ~ semantic_representative
    ),
    mapping_type = case_when(
      keep_granular ~ "granular_self",
      !is.na(cluster_common_ancestor) ~ "ontology_ancestor",
      TRUE ~ "semantic_representative"
    ),
    granularity = ifelse(
      keep_granular, "granular", "high_level"
    ),
    coarse_display_reactome_id = case_when(
      !is.na(fallback_common_ancestor) ~ fallback_common_ancestor,
      !is.na(individual_rollup_ancestor) ~ individual_rollup_ancestor,
      TRUE ~ semantic_representative
    ),
    coarse_mapping_type = case_when(
      !is.na(fallback_common_ancestor) ~ "ontology_ancestor",
      !is.na(individual_rollup_ancestor) ~ "ontology_ancestor",
      TRUE ~ "semantic_representative"
    ),
    coarse_rollup_source = case_when(
      !is.na(fallback_common_ancestor) ~ "cluster_common_ancestor",
      !is.na(individual_rollup_ancestor) ~ "individual_ancestor",
      TRUE ~ "semantic_representative"
    ),
    coarse_is_strict_ancestor = map2_lgl(
      reactome_id,
      coarse_display_reactome_id,
      function(member_id, coarse_id) {
        !is.na(coarse_id) & coarse_id %in% get_ancestors(member_id)
      }
    ),
    rollup_basis = case_when(
      keep_granular ~ "distinct_or_uncertain_retention",
      !evidence_pass ~ "low_granular_evidence",
      redundancy_supported ~ "redundancy_supported",
      TRUE ~ "review_required"
    ),
    decision_reason = case_when(
      keep_granular & reactome_id == semantic_representative ~
        "evidence-supported representative retained as a granular candidate",
      keep_granular & similarity_uncertain ~
        "term retained because redundancy evidence was insufficient",
      keep_granular ~
        "term retained because its NES or leading-edge pattern differed",
      !evidence_pass & !is.na(cluster_common_ancestor) ~
        paste(
          "term with low granular evidence moved to a strict",
          "Reactome ancestor; this is a presentation roll-up, not redundancy evidence"
        ),
      !is.na(cluster_common_ancestor) ~
        "semantic cluster collapsed to a verified common Reactome ancestor",
      !evidence_pass ~
        paste(
          "low-granular-evidence term grouped by semantic label;",
          "review because no informative strict common ancestor was available"
        ),
      TRUE ~
        "no informative common Reactome ancestor; semantic representative used"
    )
  ) |>
  left_join(
    reactome_meta |> select(display_reactome_id = reactome_id, display_term = reactome_term),
    by = "display_reactome_id"
  ) |>
  left_join(
    reactome_meta |>
      select(coarse_display_reactome_id = reactome_id, coarse_display_term = reactome_term),
    by = "coarse_display_reactome_id"
  ) |>
  mutate(
    display_term = ifelse(
      is.na(display_term),
      paste0("Semantic group represented by ", semantic_representative),
      display_term
    ),
    coarse_display_term = ifelse(
      is.na(coarse_display_term),
      paste0("Semantic group represented by ", coarse_display_reactome_id),
      coarse_display_term
    ),
    display_key = paste(
      mapping_type,
      direction_stratum,
      semantic_representative,
      display_reactome_id,
      sep = "|"
    )
  )

# Assert that every claimed ontology roll-up is genuinely ancestral.
ancestry_check <- pmap_lgl(
  display_mapping |>
    select(reactome_id, display_reactome_id, mapping_type),
  function(reactome_id, display_reactome_id, mapping_type) {
    if (mapping_type != "ontology_ancestor") {
      return(TRUE)
    }
    display_reactome_id %in% get_ancestors(reactome_id)
  }
)

if (!all(ancestry_check)) {
  readr::write_csv(
    display_mapping[!ancestry_check, ],
    file.path(config$output_dir, "failed_ancestry_assertions.csv")
  )
  stop("At least one claimed ontology roll-up failed the ancestry assertion")
}

if (
  nrow(display_mapping) != length(candidate_ids) ||
    anyDuplicated(display_mapping$reactome_id)
) {
  stop("Display mapping must contain exactly one row per candidate Reactome term")
}

# -----------------------------------------------------------------------------
# Select one fixed tested source term for every display row
# -----------------------------------------------------------------------------

tested_reactome_ids <- all_results_mapped |>
  filter(!is.na(reactome_id), is.finite(NES)) |>
  distinct(reactome_id) |>
  pull(reactome_id)

display_sources <- display_mapping |>
  group_by(
    display_key, display_reactome_id, display_term, mapping_type,
    direction_stratum
  ) |>
  arrange(
    desc(n_tested_cell_types), desc(presentation_score), reactome_id,
    .by_group = TRUE
  ) |>
  summarise(
    source_reactome_id = if (first(display_reactome_id) %in% reactome_id &&
                       first(display_reactome_id) %in% tested_reactome_ids) {
      first(display_reactome_id)
    } else {
      first(reactome_id)
    },
    source_type = ifelse(
      first(display_reactome_id) %in% reactome_id &
        first(display_reactome_id) %in% tested_reactome_ids,
      "display_term",
      "representative_member"
    ),
    n_collapsed_terms = n_distinct(reactome_id),
    .groups = "drop"
  ) |>
  left_join(
    reactome_meta |> select(source_reactome_id = reactome_id, source_term = reactome_term),
    by = "source_reactome_id"
  ) |>
  left_join(
    term_metrics |>
      select(
        source_reactome_id = reactome_id,
        source_n_tested_cell_types = n_tested_cell_types
      ),
    by = "source_reactome_id"
  )

if (anyDuplicated(display_sources$display_key)) {
  stop("Each display key must have exactly one fixed source term")
}

display_mapping <- display_mapping |>
  left_join(
    display_sources |>
      select(
        display_key, source_reactome_id, source_term, source_type,
        source_n_tested_cell_types, n_collapsed_terms
      ),
    by = "display_key"
  )

# -----------------------------------------------------------------------------
# Output annotated data, mapping, summaries and plot
# -----------------------------------------------------------------------------

readr::write_csv(
  term_metrics |> arrange(best_q, desc(presentation_score), reactome_id),
  file.path(config$output_dir, "reactome_term_metrics.csv")
)

readr::write_csv(
  semantic_clusters |> arrange(semantic_representative, reactome_id),
  file.path(config$output_dir, "semantic_cluster_membership.csv")
)

mapping_output <- display_mapping |>
  select(
    pathway, reactome_id, reactome_term, min_depth, max_depth, information_content,
    dominant_direction, direction_consistency, direction_stratum,
    n_granular_evidence_cell_types,
    structural_cluster_representative, semantic_representative,
    semantic_similarity_to_representative,
    representative_stable, profile_correlation,
    direction_concordance_to_representative,
    median_leading_edge_jaccard, median_leading_edge_overlap,
    profile_similarity_available, leading_edge_similarity_available,
    profile_concordant, leading_edge_concordant, redundancy_supported,
    similarity_uncertain, signal_distinct_from_representative,
    evidence_pass, presentation_score, keep_granular,
    display_key, display_reactome_id, display_term,
    mapping_type, granularity, rollup_basis, decision_reason,
    cluster_common_ancestor, fallback_common_ancestor,
    individual_rollup_ancestor,
    coarse_display_reactome_id, coarse_display_term, coarse_mapping_type,
    coarse_rollup_source, coarse_is_strict_ancestor,
    source_reactome_id, source_term,
    source_type, source_n_tested_cell_types, n_collapsed_terms
  )

if (anyDuplicated(mapping_output$reactome_id)) {
  stop("The exported member-to-display mapping has duplicate Reactome IDs")
}

readr::write_csv(
  mapping_output |> arrange(display_term, reactome_term),
  file.path(config$output_dir, "adaptive_display_mapping.csv")
)

annotated_results <- all_results_mapped |>
  left_join(
    mapping_output |>
      select(
        reactome_id, direction_stratum, display_key,
        display_reactome_id, display_term, mapping_type, granularity,
        rollup_basis, decision_reason, source_reactome_id, source_term, source_type,
        n_collapsed_terms, keep_granular,
        min_depth, max_depth, information_content, semantic_representative,
        presentation_score,
        signal_distinct_from_representative,
        coarse_display_reactome_id, coarse_display_term, coarse_mapping_type,
        coarse_rollup_source, coarse_is_strict_ancestor
      ),
    by = "reactome_id"
  )

if (nrow(annotated_results) != nrow(all_results_mapped)) {
  stop("Annotating results changed the number of GSEA rows")
}

readr::write_csv(
  annotated_results,
  file.path(config$output_dir, "reactome_results_annotated.csv")
)

# Descriptive member evidence is labelled explicitly and is never copied to a
# display ancestor or semantic representative as though it had been tested.
member_by_cell <- annotated_results |>
  filter(!is.na(display_key), evidence_pass_row) |>
  group_by(display_key, cell_type) |>
  summarise(
    best_member_q = safe_min(evidence_q),
    largest_member_abs_NES = safe_max(abs(NES)),
    n_significant_members = n_distinct(
      reactome_id[evidence_q <= config$fdr_threshold]
    ),
    n_positive_members = n_distinct(
      reactome_id[evidence_q <= config$fdr_threshold & NES > 0]
    ),
    n_negative_members = n_distinct(
      reactome_id[evidence_q <= config$fdr_threshold & NES < 0]
    ),
    .groups = "drop"
  )

if (anyDuplicated(member_by_cell[c("display_key", "cell_type")])) {
  stop("Grouped member evidence is not unique by display key and cell type")
}

source_results <- all_results_mapped |>
  filter(!is.na(reactome_id)) |>
  inner_join(
    display_sources |>
      select(display_key, source_reactome_id, source_type),
    by = c("reactome_id" = "source_reactome_id")
  ) |>
  arrange(display_key, cell_type, is.na(evidence_q), evidence_q) |>
  group_by(display_key, cell_type) |>
  slice(1) |>
  ungroup() |>
  transmute(
    display_key,
    cell_type,
    source_pathway = pathway,
    source_NES = NES,
    source_ES = ES,
    source_p_value = pval,
    source_q_value = evidence_q,
    source_leading_edge = leading_edge
  )

if (anyDuplicated(source_results[c("display_key", "cell_type")])) {
  stop("Source results are not unique by display key and cell type")
}

summary_grid <- tidyr::crossing(
  display_key = display_sources$display_key,
  cell_type = cell_types
)

adaptive_summary <- summary_grid |>
  left_join(display_sources, by = "display_key") |>
  left_join(source_results, by = c("display_key", "cell_type")) |>
  left_join(member_by_cell, by = c("display_key", "cell_type")) |>
  mutate(
    case = contrast_table$case[[1]],
    reference = contrast_table$reference[[1]],
    comparison = contrast_table$comparison[[1]],
    source_direction = case_when(
      source_NES > 0 ~ "positive",
      source_NES < 0 ~ "negative",
      TRUE ~ NA_character_
    ),
    source_significant = is.finite(source_q_value) &
      source_q_value <= config$fdr_threshold,
    member_evidence_present = is.finite(best_member_q) &
      best_member_q <= config$fdr_threshold,
    member_evidence_without_source_significance =
      member_evidence_present & !source_significant
  )

expected_summary_rows <- nrow(display_sources) * length(cell_types)
if (
  nrow(adaptive_summary) != expected_summary_rows ||
    anyDuplicated(adaptive_summary[c("display_key", "cell_type")])
) {
  stop("Adaptive summary is not one row per display key and cell type")
}

readr::write_csv(
  adaptive_summary,
  file.path(config$output_dir, "adaptive_reactome_summary_long.csv")
)

# -----------------------------------------------------------------------------
# Primary compact significant-term output
# -----------------------------------------------------------------------------

# The complete grid above is useful for matrices and audit, but it is not the
# final biological list. Build the final list from significant tested rows so
# every displayed term has explicit evidence provenance. A member is shown at
# its own Reactome ID only when it passes the row-level granular evidence rule and
# was selected for granular retention. Otherwise it is moved to the best
# available verified ontology ancestor; semantic grouping is
# used only when no informative strict ancestor exists.
final_member_rows_raw <- annotated_results |>
  filter(
    !is.na(reactome_id),
    !is.na(display_key),
    evidence_pass_row,
    is.finite(evidence_q),
    evidence_q <= config$fdr_threshold,
    is.finite(NES)
  )

# Different MSigDB pathway labels can occasionally resolve to the same Reactome ID.
# They are aliases of one ontology concept and must not consume multiple
# granular-selection slots. Keep the strongest tested row for each
# Reactome-ID x cell-type x comparison combination.
final_member_evidence <- final_member_rows_raw |>
  arrange(
    reactome_id,
    cell_type,
    comparison,
    evidence_q,
    desc(abs(NES)),
    desc(leading_edge_fraction),
    pathway
  ) |>
  distinct(reactome_id, cell_type, comparison, .keep_all = TRUE) |>
  mutate(
    member_direction = case_when(
      NES > 0 ~ "positive",
      NES < 0 ~ "negative",
      TRUE ~ "zero"
    ),
    granular_candidate_here = granular_evidence_pass_row &
      is.finite(min_depth) &
      min_depth >= config$final_granular_min_depth &
      coalesce(keep_granular, FALSE) &
      evidence_q <= config$final_granular_fdr &
      abs(NES) >= config$final_granular_min_abs_nes &
      (
        reactome_id == semantic_representative |
          coalesce(signal_distinct_from_representative, FALSE)
      )
  )

# Select granular representatives once across the complete experiment rather
# than independently in every cell type. This prevents a semantic cluster from
# contributing up to n_cell_types near-duplicate sibling terms to the final
# list. Cell-specific evidence is retained when the selected term is displayed.
global_granular_selection <- final_member_evidence |>
  filter(granular_candidate_here) |>
  group_by(
    member_direction,
    semantic_representative,
    reactome_id,
    reactome_term
  ) |>
  summarise(
    selection_max_depth = safe_max(max_depth),
    selection_information_content = safe_max(information_content),
    selection_presentation_score = safe_max(presentation_score),
    n_candidate_cell_types = n_distinct(cell_type),
    best_candidate_q_value = safe_min(evidence_q),
    largest_candidate_abs_NES = safe_max(abs(NES)),
    .groups = "drop"
  ) |>
  group_by(member_direction, semantic_representative) |>
  arrange(
    desc(selection_max_depth),
    desc(selection_information_content),
    desc(n_candidate_cell_types),
    desc(selection_presentation_score),
    best_candidate_q_value,
    desc(largest_candidate_abs_NES),
    reactome_id,
    .by_group = TRUE
  ) |>
  mutate(
    global_granular_rank = row_number(),
    global_granular_limit = config$final_max_granular_per_cluster,
    global_granular_selected =
      global_granular_rank <= global_granular_limit
  ) |>
  ungroup()

readr::write_csv(
  global_granular_selection |>
    arrange(
      member_direction,
      semantic_representative,
      global_granular_rank,
      reactome_id
    ),
  file.path(
    config$output_dir,
    "adaptive_reactome_global_granular_selection.csv"
  )
)

final_member_rows <- final_member_evidence |>
  left_join(
    global_granular_selection |>
      select(
        member_direction,
        semantic_representative,
        reactome_id,
        global_granular_rank,
        global_granular_limit,
        global_granular_selected,
        n_candidate_cell_types,
        best_candidate_q_value,
        largest_candidate_abs_NES
      ),
    by = c(
      "member_direction",
      "semantic_representative",
      "reactome_id"
    )
  ) |>
  mutate(
    global_granular_selected = coalesce(
      global_granular_selected,
      FALSE
    ),
    granular_display_here =
      granular_candidate_here & global_granular_selected
  ) |>
  group_by(cell_type, member_direction, semantic_representative) |>
  mutate(
    cluster_has_granular_display = any(granular_display_here),
    cluster_reactome_ids = collapse_unique_text(reactome_id),
    cluster_terms = collapse_unique_text(reactome_term)
  ) |>
  ungroup() |>
  # Once a cluster has an evidence-supported granular winner, do not also emit
  # a coarser or semantic label for the same cluster. Non-winning members stay
  # available in reactome_results_annotated.csv and are counted below for audit.
  filter(granular_display_here | !cluster_has_granular_display) |>
  mutate(
    final_display_reactome_id = case_when(
      granular_display_here ~ reactome_id,
      !is.na(coarse_display_reactome_id) ~ coarse_display_reactome_id,
      TRUE ~ display_reactome_id
    ),
    direct_display_evidence = reactome_id == final_display_reactome_id,
    final_mapping_type = case_when(
      granular_display_here ~ "granular_self",
      direct_display_evidence ~ "direct_self",
      !is.na(coarse_mapping_type) ~ coarse_mapping_type,
      TRUE ~ mapping_type
    ),
    final_rollup_source = case_when(
      granular_display_here ~ "direct_granular_evidence",
      direct_display_evidence ~ "direct_tested_term",
      !is.na(coarse_rollup_source) ~ coarse_rollup_source,
      TRUE ~ "existing_display_mapping"
    ),
    final_is_strict_ontology_rollup =
      !direct_display_evidence & coalesce(coarse_is_strict_ancestor, FALSE)
  ) |>
  filter(!final_display_reactome_id %in% uninformative_parent_ids) |>
  left_join(
    reactome_meta |>
      select(final_display_reactome_id = reactome_id, final_display_term = reactome_term),
    by = "final_display_reactome_id"
  ) |>
  mutate(
    final_display_term = ifelse(
      is.na(final_display_term),
      paste0("Semantic group represented by ", final_display_reactome_id),
      final_display_term
    )
  )

if (nrow(final_member_rows) == 0) {
  stop("No mapped significant member rows were available for the final list")
}

# Select one real tested result as the evidence example for every displayed
# term and cell. Prefer a directly tested display term; otherwise use the
# strongest supporting member. Statistics always retain the evidence_* prefix
# so an ancestor never appears to have inherited a child's NES or FDR.
best_cell_evidence <- final_member_rows |>
  arrange(
    cell_type,
    final_display_reactome_id,
    desc(direct_display_evidence),
    evidence_q,
    desc(abs(NES)),
    desc(leading_edge_fraction),
    reactome_id
  ) |>
  group_by(cell_type, final_display_reactome_id) |>
  slice_head(n = 1) |>
  ungroup() |>
  transmute(
    cell_type,
    display_reactome_id = final_display_reactome_id,
    evidence_source = ifelse(
      direct_display_evidence,
      "display_term",
      "supporting_member"
    ),
    evidence_reactome_id = reactome_id,
    evidence_term = reactome_term,
    evidence_pathway = pathway,
    evidence_NES = NES,
    evidence_ES = ES,
    evidence_p_value = pval,
    evidence_q_value = evidence_q,
    evidence_leading_edge_n = leading_edge_n,
    evidence_leading_edge_fraction = leading_edge_fraction,
    evidence_leading_edge = leading_edge
  )

cell_term_summary <- final_member_rows |>
  group_by(cell_type, final_display_reactome_id, final_display_term) |>
  summarise(
    case = one_or_na(case),
    reference = one_or_na(reference),
    comparison = one_or_na(comparison),
    n_significant_supporting_terms = n_distinct(reactome_id),
    n_granular_supporting_terms = n_distinct(
      reactome_id[granular_display_here]
    ),
    n_positive_supporting_terms = n_distinct(
      reactome_id[member_direction == "positive"]
    ),
    n_negative_supporting_terms = n_distinct(
      reactome_id[member_direction == "negative"]
    ),
    best_member_q_value = safe_min(evidence_q),
    largest_member_abs_NES = safe_max(abs(NES)),
    supporting_reactome_ids = collapse_unique_text(reactome_id),
    supporting_terms = collapse_unique_text(reactome_term),
    cluster_reactome_ids = collapse_delimited_values(cluster_reactome_ids),
    cluster_terms = collapse_delimited_values(cluster_terms),
    semantic_cluster_ids = collapse_unique_text(semantic_representative),
    selected_granular_reactome_ids = collapse_unique_text(
      reactome_id[granular_display_here]
    ),
    mapping_types = collapse_unique_text(final_mapping_type),
    rollup_sources = collapse_unique_text(final_rollup_source),
    has_direct_display_evidence = any(direct_display_evidence),
    has_granular_display_evidence = any(
      granular_display_here & direct_display_evidence
    ),
    all_rollups_are_strict_ancestors = all(
      granular_display_here |
        direct_display_evidence |
        final_is_strict_ontology_rollup |
        final_mapping_type == "semantic_representative"
    ),
    uses_semantic_fallback = any(
      final_mapping_type == "semantic_representative"
    ),
    .groups = "drop"
  ) |>
  rename(
    display_reactome_id = final_display_reactome_id,
    display_term = final_display_term
  ) |>
  mutate(
    n_cluster_significant_terms = vapply(
      cluster_reactome_ids,
      function(x) length(split_delimited_values(x)),
      integer(1)
    ),
    direction = case_when(
      n_positive_supporting_terms > 0 & n_negative_supporting_terms > 0 ~
        "mixed",
      n_positive_supporting_terms > 0 ~ "positive",
      n_negative_supporting_terms > 0 ~ "negative",
      TRUE ~ "zero"
    ),
    resolution = case_when(
      has_granular_display_evidence ~ "granular",
      has_direct_display_evidence ~ "direct",
      uses_semantic_fallback ~ "semantic_group",
      TRUE ~ "ontology_rollup"
    )
  ) |>
  left_join(
    best_cell_evidence,
    by = c("cell_type", "display_reactome_id")
  )

# Count cluster members not represented elsewhere in the same cell. Comparing
# sets, rather than subtracting row counts, avoids labelling a second retained
# granular winner as "not displayed".
cell_displayed_id_sets <- cell_term_summary |>
  group_by(cell_type) |>
  summarise(
    all_displayed_reactome_ids_in_cell = collapse_unique_text(display_reactome_id),
    .groups = "drop"
  )

cell_term_summary <- cell_term_summary |>
  left_join(cell_displayed_id_sets, by = "cell_type") |>
  mutate(
    n_cluster_terms_not_displayed = map2_int(
      cluster_reactome_ids,
      all_displayed_reactome_ids_in_cell,
      function(cluster_ids, displayed_ids) {
        length(setdiff(
          split_delimited_values(cluster_ids),
          split_delimited_values(displayed_ids)
        ))
      }
    )
  ) |>
  select(-all_displayed_reactome_ids_in_cell)

display_depth_table <- map_dfr(
  sort(unique(cell_term_summary$display_reactome_id)),
  function(id) {
    depth <- get_depth_range(id)
    tibble(
      display_reactome_id = id,
      display_min_depth = unname(depth[["min_depth"]]),
      display_max_depth = unname(depth[["max_depth"]]),
      display_information_content = unname(ic_vector[id])
    )
  }
) |>
  mutate(
    display_information_content = ifelse(
      is.finite(display_information_content),
      display_information_content,
      NA_real_
    )
  )

cell_term_summary <- cell_term_summary |>
  left_join(display_depth_table, by = "display_reactome_id")

if (
  anyDuplicated(cell_term_summary[c("cell_type", "display_reactome_id")]) ||
    any(!cell_term_summary$all_rollups_are_strict_ancestors) ||
    any(cell_term_summary$display_reactome_id %in% uninformative_parent_ids)
) {
  stop(
    "Compact Reactome rows failed uniqueness, ancestry or informative-label checks"
  )
}

# Remove a displayed parent when a more specific displayed descendant exists
# in the same cell population. Positive and negative support are combined into
# direction='mixed' where necessary, so this cannot create duplicate Reactome terms.
cell_hierarchy_flags <- cell_term_summary |>
  group_split(cell_type, .keep = TRUE) |>
  map_dfr(function(x) {
    get_hierarchy_pruning_flags(x$display_reactome_id) |>
      mutate(cell_type = x$cell_type[[1]], .before = 1)
  })

cell_term_summary_all <- cell_term_summary |>
  left_join(
    cell_hierarchy_flags,
    by = c("cell_type", "display_reactome_id")
  )

final_by_cell <- cell_term_summary_all |>
  filter(!hierarchy_parent_removed) |>
  mutate(
    scope_id = config$analysis_scope,
    scope_label = config$analysis_scope_label,
    scope_completeness = config$scope_completeness,
    .before = 1
  ) |>
  arrange(
    cell_type,
    desc(resolution == "granular"),
    evidence_q_value,
    desc(largest_member_abs_NES),
    display_reactome_id
  )

if (
  anyDuplicated(final_by_cell[c("cell_type", "display_reactome_id")]) ||
    any(!is.finite(final_by_cell$best_member_q_value)) ||
    any(final_by_cell$best_member_q_value > config$fdr_threshold)
) {
  stop("Final by-cell Reactome list failed uniqueness or significance checks")
}

remaining_cell_parent_flags <- final_by_cell |>
  group_split(cell_type, .keep = TRUE) |>
  map_lgl(function(x) {
    any(get_hierarchy_pruning_flags(
      x$display_reactome_id
    )$hierarchy_parent_removed)
  })

if (any(remaining_cell_parent_flags)) {
  stop("A parent-child pair remains in the final by-cell Reactome list")
}

readr::write_csv(
  final_by_cell |>
    select(-hierarchy_parent_removed, -retained_descendant_reactome_ids),
  file.path(config$output_dir, "adaptive_reactome_final_by_cell.csv")
)

# Aggregate the cell-specific frontier to one row per Reactome term. This is the
# designated final result for reporting: exact Reactome IDs are unique, directions
# and cell populations are retained in explicit list/count fields, and a final
# scope-wide antichain removes any parent whose retained child appears elsewhere.
final_by_cell <- final_by_cell |>
  mutate(
    cell_evidence = sprintf(
      "%s[%s,q=%.4g,NES=%.4g]",
      cell_type,
      direction,
      evidence_q_value,
      evidence_NES
    )
  )

representative_global_evidence <- final_by_cell |>
  arrange(
    display_reactome_id,
    desc(evidence_source == "display_term"),
    evidence_q_value,
    desc(abs(evidence_NES)),
    cell_type
  ) |>
  group_by(display_reactome_id) |>
  slice_head(n = 1) |>
  ungroup() |>
  transmute(
    display_reactome_id,
    representative_evidence_cell_type = cell_type,
    representative_evidence_source = evidence_source,
    representative_evidence_reactome_id = evidence_reactome_id,
    representative_evidence_term = evidence_term,
    representative_evidence_pathway = evidence_pathway,
    representative_evidence_NES = evidence_NES,
    representative_evidence_ES = evidence_ES,
    representative_evidence_p_value = evidence_p_value,
    representative_evidence_q_value = evidence_q_value,
    representative_evidence_leading_edge_n = evidence_leading_edge_n,
    representative_evidence_leading_edge_fraction =
      evidence_leading_edge_fraction,
    representative_evidence_leading_edge = evidence_leading_edge
  )

global_term_summary <- final_by_cell |>
  group_by(display_reactome_id, display_term) |>
  summarise(
    case = one_or_na(case),
    reference = one_or_na(reference),
    comparison = one_or_na(comparison),
    resolution = case_when(
      any(resolution == "granular") ~ "granular",
      any(resolution == "direct") ~ "direct",
      any(resolution == "semantic_group") ~ "semantic_group",
      TRUE ~ "ontology_rollup"
    ),
    display_min_depth = first(display_min_depth),
    display_max_depth = first(display_max_depth),
    display_information_content = first(display_information_content),
    overall_direction = case_when(
      any(n_positive_supporting_terms > 0) &
        any(n_negative_supporting_terms > 0) ~ "mixed",
      any(n_positive_supporting_terms > 0) ~ "positive",
      any(n_negative_supporting_terms > 0) ~ "negative",
      TRUE ~ "zero"
    ),
    n_significant_cell_types = n_distinct(cell_type),
    significant_cell_types = collapse_unique_text(cell_type),
    positive_cell_types = collapse_unique_text(
      cell_type[n_positive_supporting_terms > 0]
    ),
    negative_cell_types = collapse_unique_text(
      cell_type[n_negative_supporting_terms > 0]
    ),
    mixed_direction_cell_types = collapse_unique_text(
      cell_type[direction == "mixed"]
    ),
    n_significant_member_cell_results = sum(n_significant_supporting_terms),
    supporting_reactome_ids = collapse_delimited_values(supporting_reactome_ids),
    supporting_terms = collapse_delimited_values(supporting_terms),
    cluster_reactome_ids = collapse_delimited_values(cluster_reactome_ids),
    cluster_terms = collapse_delimited_values(cluster_terms),
    semantic_cluster_ids = collapse_delimited_values(semantic_cluster_ids),
    selected_granular_reactome_ids = collapse_delimited_values(
      selected_granular_reactome_ids
    ),
    mapping_types = collapse_delimited_values(mapping_types),
    rollup_sources = collapse_delimited_values(rollup_sources),
    best_member_q_value = safe_min(best_member_q_value),
    largest_member_abs_NES = safe_max(largest_member_abs_NES),
    cell_type_evidence = paste(sort(cell_evidence), collapse = ";"),
    .groups = "drop"
  ) |>
  mutate(
    n_unique_supporting_reactome_terms = vapply(
      supporting_reactome_ids,
      function(x) length(split_delimited_values(x)),
      integer(1)
    ),
    n_unique_cluster_reactome_terms = vapply(
      cluster_reactome_ids,
      function(x) length(split_delimited_values(x)),
      integer(1)
    )
  ) |>
  left_join(representative_global_evidence, by = "display_reactome_id")

global_hierarchy_flags <- get_hierarchy_pruning_flags(
  global_term_summary$display_reactome_id
)

global_term_summary_all <- global_term_summary |>
  left_join(global_hierarchy_flags, by = "display_reactome_id")

final_significant_terms <- global_term_summary_all |>
  filter(!hierarchy_parent_removed) |>
  select(-hierarchy_parent_removed, -retained_descendant_reactome_ids) |>
  mutate(
    scope_id = config$analysis_scope,
    scope_label = config$analysis_scope_label,
    scope_completeness = config$scope_completeness,
    .before = 1
  ) |>
  arrange(
    desc(resolution == "granular"),
    best_member_q_value,
    desc(largest_member_abs_NES),
    display_reactome_id
  )

if (
  anyDuplicated(final_significant_terms$display_reactome_id) ||
    any(get_hierarchy_pruning_flags(
      final_significant_terms$display_reactome_id
    )$hierarchy_parent_removed) ||
    any(final_significant_terms$n_significant_cell_types < 1L)
) {
  stop("Final unique Reactome list failed uniqueness, hierarchy or evidence checks")
}

readr::write_csv(
  final_significant_terms,
  file.path(config$output_dir, "adaptive_reactome_final_significant_terms.csv")
)

cell_pruning_audit <- cell_term_summary_all |>
  filter(hierarchy_parent_removed) |>
  transmute(
    pruning_scope = "within_cell_type",
    context = cell_type,
    removed_reactome_id = display_reactome_id,
    removed_term = display_term,
    retained_descendant_reactome_ids,
    removal_reason = paste(
      "A more specific significant display term connected by strict pathway",
      "ancestry was retained in the same cell population"
    ),
    n_significant_cell_types = 1L,
    significant_cell_types = cell_type,
    direction,
    supporting_reactome_ids,
    cell_type_evidence = NA_character_,
    best_member_q_value,
    largest_member_abs_NES
  )

global_pruning_audit <- global_term_summary_all |>
  filter(hierarchy_parent_removed) |>
  transmute(
    pruning_scope = "global_unique_list",
    context = config$analysis_scope,
    removed_reactome_id = display_reactome_id,
    removed_term = display_term,
    retained_descendant_reactome_ids,
    removal_reason = paste(
      "A more specific significant display term connected by strict pathway",
      "ancestry was retained in the final scope-wide list"
    ),
    n_significant_cell_types,
    significant_cell_types,
    direction = overall_direction,
    supporting_reactome_ids,
    cell_type_evidence,
    best_member_q_value,
    largest_member_abs_NES
  )

hierarchy_pruning_audit <- bind_rows(
  cell_pruning_audit,
  global_pruning_audit
) |>
  arrange(pruning_scope, context, removed_term, removed_reactome_id)

readr::write_csv(
  hierarchy_pruning_audit,
  file.path(config$output_dir, "adaptive_reactome_hierarchy_pruning_audit.csv")
)

nes_matrix_output <- adaptive_summary |>
  mutate(
    display_label = paste0(
      display_term, " (", display_reactome_id, ")",
      " {", direction_stratum, "}"
    )
  ) |>
  select(
    display_key, display_reactome_id, display_term, mapping_type,
    direction_stratum, source_reactome_id, source_type, display_label, cell_type,
    source_NES
  ) |>
  tidyr::pivot_wider(names_from = cell_type, values_from = source_NES)

readr::write_csv(
  nes_matrix_output,
  file.path(config$output_dir, "adaptive_Reactome_NES_matrix.csv")
)

source_q_matrix_output <- adaptive_summary |>
  select(
    display_key, display_reactome_id, display_term, mapping_type,
    direction_stratum, source_reactome_id, source_type, cell_type, source_q_value
  ) |>
  tidyr::pivot_wider(names_from = cell_type, values_from = source_q_value)

readr::write_csv(
  source_q_matrix_output,
  file.path(config$output_dir, "adaptive_Reactome_source_FDR_matrix.csv")
)

collapse_qc <- mapping_output |>
  group_by(
    display_key, display_reactome_id, display_term, mapping_type,
    direction_stratum, source_reactome_id, source_term, source_type
  ) |>
  summarise(
    n_terms = n_distinct(reactome_id),
    median_semantic_similarity = safe_median(
      semantic_similarity_to_representative[
        reactome_id != semantic_representative
      ]
    ),
    median_profile_correlation = safe_median(
      profile_correlation[reactome_id != semantic_representative]
    ),
    median_direction_concordance = safe_median(
      direction_concordance_to_representative[
        reactome_id != semantic_representative
      ]
    ),
    group_median_leading_edge_jaccard = safe_median(
      median_leading_edge_jaccard[reactome_id != semantic_representative]
    ),
    group_median_leading_edge_overlap = safe_median(
      median_leading_edge_overlap[reactome_id != semantic_representative]
    ),
    n_similarity_uncertain = sum(similarity_uncertain, na.rm = TRUE),
    n_mixed_direction_terms = sum(dominant_direction == "mixed", na.rm = TRUE),
    all_representatives_stable = all(representative_stable %in% TRUE),
    .groups = "drop"
  ) |>
  mutate(
    review_flag =
      mapping_type == "semantic_representative" |
      (!is.na(median_semantic_similarity) &
         median_semantic_similarity < config$semantic_threshold) |
      (!is.na(median_profile_correlation) &
         median_profile_correlation < config$min_profile_correlation) |
      (!is.na(median_direction_concordance) &
         median_direction_concordance < config$min_direction_concordance) |
      (!is.na(group_median_leading_edge_jaccard) &
         group_median_leading_edge_jaccard <
           config$min_leading_edge_jaccard) |
      (!is.na(group_median_leading_edge_overlap) &
         group_median_leading_edge_overlap <
           config$min_leading_edge_overlap) |
      n_similarity_uncertain > 0 |
      n_mixed_direction_terms > 0 |
      !all_representatives_stable
  )

readr::write_csv(
  collapse_qc,
  file.path(config$output_dir, "collapse_qc.csv")
)

reduction_counts <- tibble(
  metric = c(
    "significant_input_rows",
    "candidate_pathways",
    "mapped_candidate_pathways",
    "omitted_candidate_pathways_mapping",
    "mapped_candidate_reactome_ids",
    "granular_terms_retained_in_display_mapping",
    "terms_rolled_up",
    "display_groups",
    "adaptive_summary_grid_rows",
    "display_groups_flagged_for_review",
    "mapped_significant_input_rows",
    "unique_significant_reactome_cell_rows",
    "significant_reactome_alias_rows_removed",
    "global_granular_reactome_direction_selections",
    "compact_member_rows_after_cluster_selection",
    "cell_term_rows_before_parent_pruning",
    "cell_parent_rows_removed",
    "cell_term_rows_after_parent_pruning",
    "unique_terms_before_global_parent_pruning",
    "global_parent_terms_removed",
    "final_unique_significant_terms",
    "final_parent_child_pairs"
  ),
  value = c(
    nrow(significant_input),
    length(candidate_pathways),
    sum(candidate_mapping$mapping_status == "mapped"),
    omitted_candidate_pathways,
    length(candidate_ids),
    sum(display_mapping$keep_granular, na.rm = TRUE),
    sum(!display_mapping$keep_granular, na.rm = TRUE),
    nrow(display_sources),
    nrow(adaptive_summary),
    sum(collapse_qc$review_flag, na.rm = TRUE),
    nrow(final_member_rows_raw),
    nrow(final_member_evidence),
    nrow(final_member_rows_raw) - nrow(final_member_evidence),
    sum(global_granular_selection$global_granular_selected, na.rm = TRUE),
    nrow(final_member_rows),
    nrow(cell_term_summary_all),
    sum(cell_term_summary_all$hierarchy_parent_removed, na.rm = TRUE),
    nrow(final_by_cell),
    nrow(global_term_summary_all),
    sum(global_term_summary_all$hierarchy_parent_removed, na.rm = TRUE),
    nrow(final_significant_terms),
    sum(
      get_hierarchy_pruning_flags(
        final_significant_terms$display_reactome_id
      )$hierarchy_parent_removed,
      na.rm = TRUE
    )
  )
)

readr::write_csv(
  reduction_counts,
  file.path(config$output_dir, "reduction_counts.csv")
)

output_dictionary <- tibble(
  file = c(
    "adaptive_reactome_final_significant_terms.csv",
    "adaptive_reactome_final_by_cell.csv",
    "adaptive_reactome_global_granular_selection.csv",
    "adaptive_reactome_hierarchy_pruning_audit.csv",
    "reactome_id_mapping_qc.csv",
    "unmapped_or_ambiguous_terms.csv",
    "candidate_mapping_failures.csv",
    "reactome_term_metrics.csv",
    "semantic_cluster_membership.csv",
    "adaptive_display_mapping.csv",
    "reactome_results_annotated.csv",
    "adaptive_reactome_summary_long.csv",
    "adaptive_Reactome_NES_matrix.csv",
    "adaptive_Reactome_source_FDR_matrix.csv",
    "collapse_qc.csv",
    "reduction_counts.csv",
    "adaptive_Reactome_dotplot.png"
  ),
  meaning = c(
    paste(
      "Primary reporting output: one row per unique significant displayed Reactome",
      "term across the analysis. Exact duplicates and strict parent/subpathway",
      "ancestor-descendant pairs are removed; granular semantic siblings are",
      "capped scope-wide. Cell populations and directions are aggregated."
    ),
    paste(
      "Cell-specific significant Reactome frontier: one row per cell type and Reactome",
      "ID, with mixed directions combined and no within-cell strict pathway/",
      "pathway ancestor pairs. evidence_* columns identify the tested source."
    ),
    paste(
      "Scope-wide ranking within each semantic cluster and direction, showing",
      "which evidence-supported granular Reactome terms passed the cluster cap."
    ),
    paste(
      "Terms removed for strict parent/subpathway hierarchy redundancy, with the",
      "retained descendant IDs and the removed term's evidence provenance."
    ),
    "Pathway-to-Reactome-ID mapping and mapping status for every tested pathway.",
    "Terms that need mapping review; no such term enters the summary.",
    "Significant candidate pathways omitted because mapping was not usable.",
    "Depth range, IC, evidence, direction and redundancy metrics by Reactome term.",
    "Hierarchy-similarity sensitivity clusters and direction-specific representatives.",
    "Auditable member-to-display mapping and the reason for each decision.",
    "Complete original GSEA table with mapping and display annotations.",
    paste(
      "Complete display-group x cell-type audit grid, including",
      "non-significant combinations. source_* values belong to one fixed",
      "tested term; this is not the primary final result."
    ),
    "Cell-type matrix of the fixed source term's NES; no synthetic parent NES.",
    "Cell-type matrix of the same fixed source term's within-cell-type FDR.",
    "Group-level similarity checks and review flags.",
    "Counts before and after adaptive summarisation.",
    paste(
      "Dot plot in which colour, size and alpha use one fixed tested source",
      "term; nonsignificant source results are faint. FDR is within Reactome",
      "and cell type."
    )
  )
)

readr::write_csv(
  output_dictionary,
  file.path(config$output_dir, "output_dictionary.csv")
)

# Plot the strongest evidence-supported display rows while retaining the complete
# table in CSV outputs.
plot_selection <- adaptive_summary |>
  group_by(
    display_key, display_reactome_id, display_term,
    direction_stratum, mapping_type
  ) |>
  summarise(
    best_source_q = safe_min(source_q_value),
    n_source_significant_cell_types = sum(source_significant, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(is.finite(best_source_q)) |>
  arrange(
    best_source_q,
    desc(n_source_significant_cell_types)
  ) |>
  slice_head(n = config$max_plot_terms)

plot_label_table <- plot_selection |>
  mutate(
    display_label = paste0(
      display_term, " (", display_reactome_id, ")",
      " {", direction_stratum, "}",
      " <", mapping_type, ">"
    ),
    display_label = vapply(
      display_label,
      function(x) paste(strwrap(x, width = 55), collapse = "\n"),
      character(1)
    )
  )

plot_data <- adaptive_summary |>
  filter(
    display_key %in% plot_selection$display_key,
    is.finite(source_NES),
    is.finite(source_q_value)
  ) |>
  left_join(
    plot_label_table |> select(display_key, display_label),
    by = "display_key"
  ) |>
  mutate(
    significance_size = pmin(-log10(pmax(source_q_value, 1e-300)), 20)
  )

if (nrow(plot_data) > 0) {
  plot_data$display_key <- factor(
    plot_data$display_key,
    levels = rev(plot_label_table$display_key)
  )

  display_labels <- setNames(
    plot_label_table$display_label,
    plot_label_table$display_key
  )

  adaptive_plot <- ggplot(
    plot_data,
    aes(
      x = cell_type,
      y = display_key,
      colour = source_NES,
      size = significance_size,
      shape = source_type,
      alpha = source_significant
    )
  ) +
    geom_point() +
    scale_colour_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0
    ) +
    scale_y_discrete(labels = display_labels) +
    scale_shape_manual(
      values = c(display_term = 16, representative_member = 17)
    ) +
    scale_alpha_manual(
      values = c(`TRUE` = 0.95, `FALSE` = 0.20),
      name = "Within-cell-type\nReactome FDR < 0.05"
    ) +
    scale_size_continuous(range = c(1, 6)) +
    labs(
      x = "Cell population",
      y = NULL,
      colour = "Source NES",
      size = "-log10(within-cell-type\nReactome FDR)",
      shape = "Statistic source",
      caption = paste(
        "Colour, size and alpha use the same fixed tested source term.",
        "FDR is within Reactome x cell type; grouped-member evidence is",
        "descriptive and reported separately."
      )
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )

  ggsave(
    file.path(config$output_dir, "adaptive_Reactome_dotplot.png"),
    adaptive_plot,
    width = 15,
    height = max(8, 0.28 * n_distinct(plot_data$display_key) + 2),
    dpi = 300,
    limitsize = FALSE
  )
}

# -----------------------------------------------------------------------------
# Reproducibility metadata
# -----------------------------------------------------------------------------

parameter_table <- tibble(
  parameter = c(
    "significant_file", "all_results_file", "profile_source",
    "reactome_id_mapping_file", "fdr_column", "fdr_threshold",
    "min_abs_nes", "min_leading_edge_genes",
    "min_leading_edge_fraction", "semantic_method", "semantic_threshold",
    "sensitivity_thresholds", "min_shared_cell_types_for_correlation",
    "min_shared_significant_cell_types_for_leading_edge",
    "min_profile_correlation",
    "min_direction_concordance", "min_leading_edge_jaccard",
    "min_leading_edge_overlap", "min_rollup_depth",
    "target_rollup_depth", "target_rollup_information_content",
    "final_granular_fdr", "final_granular_min_abs_nes",
    "final_granular_min_depth", "final_max_granular_per_cluster",
    "minimum_mapping_fraction",
    "mapping_fraction"
  ),
  value = as.character(c(
    config$significant_file,
    ifelse(is.na(config$all_results_file), "NA", config$all_results_file),
    profile_source,
    ifelse(is.na(config$reactome_id_mapping_file), "NA", config$reactome_id_mapping_file),
    config$fdr_column,
    config$fdr_threshold,
    config$min_abs_nes,
    config$min_leading_edge_genes,
    config$min_leading_edge_fraction,
    config$semantic_method,
    config$semantic_threshold,
    paste(config$sensitivity_thresholds, collapse = ";"),
    config$min_shared_cell_types_for_correlation,
    config$min_shared_significant_cell_types_for_leading_edge,
    config$min_profile_correlation,
    config$min_direction_concordance,
    config$min_leading_edge_jaccard,
    config$min_leading_edge_overlap,
    config$min_rollup_depth,
    config$target_rollup_depth,
    config$target_rollup_information_content,
    config$final_granular_fdr,
    config$final_granular_min_abs_nes,
    config$final_granular_min_depth,
    config$final_max_granular_per_cluster,
    config$minimum_mapping_fraction,
    mapping_fraction
  ))
)

parameter_table <- bind_rows(
  parameter_table,
  tibble(
    parameter = c(
      "analysis_scope", "analysis_scope_label", "analysis_scope_type",
      "scope_completeness", "requested_cell_types", "included_cell_types",
      "missing_cell_types", "n_requested_cell_types",
      "n_included_cell_types", "n_missing_cell_types",
      "n_scope_significant_rows", "n_scope_complete_result_rows",
      "shared_similarity_cache_file", "semantic_similarity_cache_source"
    ),
    value = as.character(c(
      config$analysis_scope,
      config$analysis_scope_label,
      config$analysis_scope_type,
      config$scope_completeness,
      collapse_manifest_values(config$requested_cell_types),
      collapse_manifest_values(config$included_cell_types),
      collapse_manifest_values(config$missing_cell_types),
      length(config$requested_cell_types),
      length(config$included_cell_types),
      length(config$missing_cell_types),
      nrow(significant_input),
      nrow(all_results),
      ifelse(
        is.na(config$shared_similarity_cache_file),
        "NA",
        config$shared_similarity_cache_file
      ),
      ifelse(is.na(cache_source_file), "not_cached", cache_source_file)
    ))
  )
)

readr::write_csv(
  parameter_table,
  file.path(config$output_dir, "run_parameters.csv")
)

input_files <- unique(c(
  unname(hierarchy$source_files),
  config$significant_file,
  if (have_full_results) config$all_results_file else character(),
  if (!is.na(config$reactome_id_mapping_file) && file.exists(config$reactome_id_mapping_file)) {
    config$reactome_id_mapping_file
  } else character()
))

checksum_table <- tibble(
  file = input_files,
  md5 = unname(tools::md5sum(input_files))
)
readr::write_csv(
  checksum_table,
  file.path(config$output_dir, "input_checksums.csv")
)

reactome_metadata <- tibble(
  name = c("source", "hierarchy_signature", "information_content_definition",
    "similarity_method", "clustering_method", "depth_convention"),
  value = c(hierarchy$source, hierarchy$signature, hierarchy$ic_definition,
    "Structural Lin; not GO semantic similarity", "complete linkage",
    "synthetic all-pathways root=0; biological top levels=1")
)
readr::write_csv(reactome_metadata, file.path(config$output_dir, "Reactome_hierarchy_metadata.csv"))

writeLines(
  capture.output(sessionInfo()),
  file.path(config$output_dir, "sessionInfo.txt")
)

write_scope_metadata(
  status = if (
    identical(config$scope_completeness, "complete")
  ) "completed" else "completed_partial",
  n_candidate_pathways = length(candidate_pathways),
  n_final_unique_terms = nrow(final_significant_terms),
  n_final_by_cell_rows = nrow(final_by_cell)
)

message(
  "Adaptive Reactome summary complete for ", config$analysis_scope_label,
  ": ", config$output_dir,
  "\nPrimary compact output: adaptive_reactome_final_significant_terms.csv",
  " (", nrow(final_significant_terms), " unique Reactome terms)",
  "\nCell-specific output: adaptive_reactome_final_by_cell.csv",
  if (length(config$missing_cell_types) > 0L) {
    paste0(
      "\nPartial scope; missing requested populations: ",
      paste(config$missing_cell_types, collapse = ", ")
    )
  } else {
    ""
  }
)
