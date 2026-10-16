#!/usr/bin/env Rscript

# Adaptive, hierarchy-aware summarisation of GO-BP preranked GSEA results.
#
# This script performs no enrichment testing. It uses already-computed GSEA
# results, calculates term-level evidence and redundancy metrics, and chooses
# an adaptive display level:
#   * biologically prioritised fibrosis branches retain selected granular terms;
#   * peripheral/redundant terms are mapped to a verified GO ancestor where
#     possible;
#   * an rrvgo representative is explicitly labelled as semantic rather than
#     as an ontology parent when no useful common ancestor exists.
# The primary output is a globally unique GO-ID list. Granular representatives
# are capped once per semantic cluster/direction, and strict is_a/part_of
# ancestors are removed when a retained descendant is available.
#
# Usage (all arguments are optional):
#   Rscript summarize_go_bp_adaptive.R \
#     significant_results.csv \
#     all_results.csv \
#     output_directory \
#     go_id_mapping.csv \
#     priority_branches.csv
#
# Use the literal string "NA" to omit an optional file argument.

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

project_dir <- path.expand(
  "~/project-area/data/crohns_scrnaseq/10c_14n_analysis/pseudobulk_dea"
)

config <- list(
  significant_file = file.path(
    project_dir, "pathway_enrichment", "go_bp_gsea_significant.csv"
  ),
  all_results_file = file.path(
    project_dir, "pathway_enrichment", "go_bp_gsea_all_results.csv"
  ),
  output_dir = file.path(
    project_dir, "pathway_enrichment", "go_bp_adaptive_summary"
  ),

  # Preferred optional mapping columns:
  #   pathway (or gs_name), go_id (or gs_exact_source), go_term (optional)
  # The supplied MSigDB CSV lacks GO IDs, so the default fallback reconstructs
  # IDs from GO.db term names and writes a complete mapping audit.
  go_id_mapping_file = NA_character_,

  # Optional CSV columns: theme, root_go_id, max_granular
  priority_file = NA_character_,

  fdr_column = "q_value_within_collection_cell_type",
  fdr_threshold = 0.05,
  min_abs_nes = 1.0,
  min_leading_edge_genes = 3L,
  min_leading_edge_fraction = 0.05,

  semantic_method = "Rel",
  semantic_threshold = 0.70,
  sensitivity_thresholds = c(0.60, 0.70, 0.80),
  max_terms_for_similarity = 4500L,

  min_shared_cell_types_for_correlation = 4L,
  min_shared_significant_cell_types_for_leading_edge = 2L,
  min_profile_correlation = 0.60,
  min_direction_concordance = 0.75,
  min_leading_edge_jaccard = 0.30,
  min_leading_edge_overlap = 0.50,

  default_max_granular_per_cluster = 3L,
  hard_max_granular_per_cluster = 6L,

  # The final reporting list uses a stricter gate than statistical
  # significance. Terms that fail this gate are retained at a coarser level;
  # they are not discarded. Priority branches retain their existing evidence
  # gate, while peripheral branches use these stricter display thresholds.
  final_granular_fdr = 0.01,
  final_granular_min_abs_nes = 1.50,
  final_granular_min_depth = 4,
  final_max_peripheral_per_cluster = 1L,

  # Peripheral clusters are rolled up only to a genuine common ancestor.
  # Depth is descriptive because GO is a DAG; information content is also used.
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
config$go_id_mapping_file <- arg_or_default(4, config$go_id_mapping_file)
config$priority_file <- arg_or_default(5, config$priority_file)

# These are deliberately excluded as uninformative roll-up labels.
uninformative_parent_ids <- c(
  "GO:0008150", # biological_process
  "GO:0009987", # cellular process
  "GO:0065007"  # biological regulation
)

# Built-in defaults can be replaced with priority_file. Multiple roots can share
# a custom display theme without being claimed to share a GO parent.
default_priority_branches <- data.frame(
  theme = c(
    "Extracellular matrix organization",
    "Extracellular matrix organization",
    "Collagen turnover",
    "Collagen turnover",
    "Fibroblast activation and expansion",
    "Fibroblast activation and expansion",
    "Fibroblast activation and expansion",
    "Fibroblast activation and expansion",
    "Fibroblast activation and expansion",
    "TGF-beta response and signalling",
    "TGF-beta response and signalling",
    "TGF-beta response and signalling",
    "Wound healing and tissue repair",
    "Wound healing and tissue repair",
    "Cell-matrix adhesion and mechanotransduction",
    "Cell-matrix adhesion and mechanotransduction",
    "Cell-matrix adhesion and mechanotransduction"
  ),
  root_go_id = c(
    "GO:0030198",
    "GO:1903053",
    "GO:0032963",
    "GO:0010712",
    "GO:0072537",
    "GO:0048144",
    "GO:0048145",
    "GO:0036446",
    "GO:1904760",
    "GO:0071559",
    "GO:0007179",
    "GO:0017015",
    "GO:0042060",
    "GO:0061041",
    "GO:0031589",
    "GO:0010810",
    "GO:0009612"
  ),
  max_granular = rep(3L, 17),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------

required_packages <- c(
  "readr", "dplyr", "tidyr", "purrr", "ggplot2",
  "AnnotationDbi", "GO.db", "org.Hs.eg.db", "GOSemSim", "rrvgo"
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
  paste0("GOBP_", x)
}

extract_go_id <- function(x) {
  out <- rep(NA_character_, length(x))
  hit <- !is.na(x) & grepl("GO:[0-9]{7}", x)
  out[hit] <- sub(".*(GO:[0-9]{7}).*", "\\1", x[hit])
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

message("Reading significant GO-BP results...")
significant_input <- read_gsea_results(config$significant_file, "significant_results")

have_full_results <- !is.na(config$all_results_file) &&
  file.exists(config$all_results_file)

if (have_full_results) {
  message("Reading complete GO-BP results for full cell-type profiles...")
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
        significant_input$collection != "GO_BP")) {
  stop("The significant input must contain GO_BP results only")
}
if (any(is.na(all_results$collection) | all_results$collection != "GO_BP")) {
  stop("The complete-results input must contain GO_BP results only")
}

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
# Map MSigDB pathway labels to GO IDs
# -----------------------------------------------------------------------------

message("Loading GO term metadata...")
go_meta <- AnnotationDbi::select(
  GO.db::GO.db,
  keys = AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID"),
  keytype = "GOID",
  columns = c("TERM", "ONTOLOGY")
) |>
  as_tibble() |>
  transmute(go_id = GOID, go_term = TERM, ontology = ONTOLOGY) |>
  filter(ontology == "BP", !is.na(go_term)) |>
  distinct(go_id, .keep_all = TRUE) |>
  mutate(pathway_key = normalise_msigdb_name(go_term))

all_pathways <- tibble(pathway = sort(unique(all_results$pathway)))

if (!is.na(config$go_id_mapping_file) && file.exists(config$go_id_mapping_file)) {
  message("Using supplied pathway-to-GO-ID mapping...")
  supplied_map <- readr::read_csv(
    config$go_id_mapping_file,
    show_col_types = FALSE,
    progress = FALSE
  )

  pathway_column <- first_existing_column(
    names(supplied_map), c("pathway", "gs_name")
  )
  id_column <- first_existing_column(
    names(supplied_map), c("go_id", "gs_exact_source")
  )

  if (is.na(pathway_column) || is.na(id_column)) {
    stop(
      "GO mapping file must contain pathway/gs_name and go_id/gs_exact_source"
    )
  }

  raw_mapping <- supplied_map |>
    transmute(
      pathway = .data[[pathway_column]],
      go_id = extract_go_id(as.character(.data[[id_column]]))
    ) |>
    filter(!is.na(pathway)) |>
    distinct()

  mapping_method <- "supplied_go_id"
} else {
  message(
    "No exact GO-ID mapping supplied; reconstructing IDs from GO.db term names..."
  )
  raw_mapping <- all_pathways |>
    left_join(
      go_meta |> select(pathway = pathway_key, go_id),
      by = "pathway"
    )
  mapping_method <- "normalised_term_name"
}

pathway_mapping <- all_pathways |>
  left_join(raw_mapping, by = "pathway") |>
  group_by(pathway) |>
  summarise(
    n_go_matches = n_distinct(go_id[!is.na(go_id)]),
    matched_go_ids = paste(sort(unique(go_id[!is.na(go_id)])), collapse = ";"),
    go_id = one_or_na(go_id),
    .groups = "drop"
  ) |>
  mutate(
    provisional_mapping_status = case_when(
      n_go_matches == 0 ~ "unmapped",
      n_go_matches == 1 ~ "mapped",
      TRUE ~ "ambiguous"
    ),
    mapping_method = mapping_method
  ) |>
  left_join(go_meta |> select(go_id, go_term, ontology), by = "go_id") |>
  mutate(
    mapping_status = case_when(
      provisional_mapping_status == "unmapped" ~ "unmapped",
      provisional_mapping_status == "ambiguous" ~ "ambiguous",
      is.na(go_term) | is.na(ontology) | ontology != "BP" ~
        "invalid_or_non_bp",
      TRUE ~ "mapped"
    )
  ) |>
  select(-provisional_mapping_status)

if (
  nrow(pathway_mapping) != nrow(all_pathways) ||
    anyDuplicated(pathway_mapping$pathway)
) {
  stop("GO-ID mapping must contain exactly one audit row per pathway")
}

readr::write_csv(
  pathway_mapping,
  file.path(config$output_dir, "go_id_mapping_qc.csv")
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
    "GO-ID mapping coverage is below ",
    100 * config$minimum_mapping_fraction,
    "%. Inspect go_id_mapping_qc.csv and supply an exact mapping containing ",
    "gs_exact_source/GO IDs from the matching MSigDB release."
  )
}

if (omitted_candidate_pathways > 0) {
  warning(
    omitted_candidate_pathways,
    " significant candidate pathways were excluded because their GO-ID ",
    "mapping was not unique and current. See candidate_mapping_failures.csv."
  )
}

all_results_mapped <- all_results |>
  left_join(
    pathway_mapping |>
      select(pathway, go_id, go_term, mapping_status, mapping_method),
    by = "pathway"
  )

significant_mapped <- significant_input |>
  left_join(
    pathway_mapping |>
      select(pathway, go_id, go_term, mapping_status, mapping_method),
    by = "pathway"
  )

candidate_ids <- candidate_mapping |>
  filter(mapping_status == "mapped") |>
  pull(go_id) |>
  unique() |>
  sort()

if (length(candidate_ids) > config$max_terms_for_similarity) {
  stop(
    "There are ", length(candidate_ids),
    " candidate GO IDs, exceeding max_terms_for_similarity=",
    config$max_terms_for_similarity,
    ". Tighten the evidence filter or increase the limit after checking memory."
  )
}

# -----------------------------------------------------------------------------
# Build a conservative GO-BP DAG using is_a and part_of only
# -----------------------------------------------------------------------------

message("Building the GO-BP parent graph...")
raw_parent_map <- as.list(GO.db::GOBPPARENTS)

clean_parent_vector <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(character())
  # GO.db returns parent GO IDs as vector values and relationship labels as
  # names (for example, name "is_a" with value "GO:0048308").
  relation <- names(x)
  if (is.null(relation)) return(character())
  relation <- gsub("[^a-z]", "", tolower(relation))
  unique(unname(x[relation %in% c("isa", "partof") & !is.na(x)]))
}

parent_map <- lapply(raw_parent_map, clean_parent_vector)
root_go_id <- "GO:0008150"

ancestor_cache <- new.env(parent = emptyenv())
depth_cache <- new.env(parent = emptyenv())

get_ancestors <- function(go_id, trail = character()) {
  if (exists(go_id, envir = ancestor_cache, inherits = FALSE)) {
    return(get(go_id, envir = ancestor_cache, inherits = FALSE))
  }
  if (go_id %in% trail) {
    stop("Cycle detected in filtered GO graph at ", go_id)
  }

  parents <- parent_map[[go_id]]
  parents <- unique(parents[!is.na(parents)])
  if (length(parents) == 0) {
    answer <- character()
  } else {
    answer <- unique(c(
      parents,
      unlist(
        lapply(parents, get_ancestors, trail = c(trail, go_id)),
        use.names = FALSE
      )
    ))
  }
  assign(go_id, answer, envir = ancestor_cache)
  answer
}

# Return a strict GO antichain flag for a set of display terms. A term is
# marked for removal only when it is a verified ancestor of another term in
# the same biological context. GO depth alone is not used because GO is a DAG.
get_hierarchy_pruning_flags <- function(go_ids) {
  ids <- sort(unique(go_ids[!is.na(go_ids) & nzchar(go_ids)]))

  if (length(ids) == 0) {
    return(tibble(
      display_go_id = character(),
      hierarchy_parent_removed = logical(),
      retained_descendant_go_ids = character()
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
    display_go_id = ids,
    hierarchy_parent_removed = unname(hierarchy_parent_removed),
    retained_descendant_go_ids = vapply(
      retained_descendants,
      collapse_unique_text,
      character(1)
    )
  )
}

get_depth_range <- function(go_id, trail = character()) {
  if (exists(go_id, envir = depth_cache, inherits = FALSE)) {
    return(get(go_id, envir = depth_cache, inherits = FALSE))
  }
  if (identical(go_id, root_go_id)) {
    answer <- c(min_depth = 0, max_depth = 0)
    assign(go_id, answer, envir = depth_cache)
    return(answer)
  }
  if (go_id %in% trail) {
    stop("Cycle detected while calculating GO depth at ", go_id)
  }

  parents <- parent_map[[go_id]]
  parents <- unique(parents[!is.na(parents)])
  if (length(parents) == 0) {
    answer <- c(min_depth = NA_real_, max_depth = NA_real_)
  } else {
    parent_depths <- do.call(
      rbind,
      lapply(parents, get_depth_range, trail = c(trail, go_id))
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
  assign(go_id, answer, envir = depth_cache)
  answer
}

candidate_ancestors <- setNames(
  lapply(candidate_ids, get_ancestors),
  candidate_ids
)

# -----------------------------------------------------------------------------
# Information content and term-level metrics
# -----------------------------------------------------------------------------

message("Calculating GO information content...")
semantic_data <- GOSemSim::godata(
  "org.Hs.eg.db",
  ont = "BP",
  computeIC = TRUE,
  keytype = "ENTREZID"
)

ic_vector <- tryCatch(
  methods::slot(semantic_data, "IC"),
  error = function(e) numeric()
)

candidate_depths <- do.call(rbind, lapply(candidate_ids, get_depth_range))
candidate_depths <- tibble(
  go_id = candidate_ids,
  min_depth = candidate_depths[, "min_depth"],
  max_depth = candidate_depths[, "max_depth"]
)

evidence_rows <- all_results_mapped |>
  filter(
    go_id %in% candidate_ids,
    evidence_pass_row,
    evidence_q <= config$fdr_threshold
  )

# Use the complete table for profiles but the configured evidence threshold for
# counts and leading-edge summaries.
term_metrics <- all_results_mapped |>
  filter(go_id %in% candidate_ids) |>
  group_by(go_id) |>
  summarise(
    pathway = first(pathway),
    go_term = first(go_term),
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
  left_join(candidate_depths, by = "go_id") |>
  mutate(
    information_content = unname(ic_vector[go_id]),
    information_content = ifelse(
      is.finite(information_content), information_content, NA_real_
    )
  )

if (!setequal(candidate_ids, term_metrics$go_id)) {
  stop(
    "Candidate GO IDs and term-metric GO IDs differ; check that the two input ",
    "files and GO mapping come from matching releases."
  )
}

# Leading-edge union and recurrent genes across significant cell populations.
leading_edge_long <- evidence_rows |>
  select(go_id, cell_type, leading_edge) |>
  tidyr::separate_rows(leading_edge, sep = ";") |>
  mutate(leading_edge = trimws(leading_edge)) |>
  filter(!is.na(leading_edge), leading_edge != "") |>
  distinct(go_id, cell_type, leading_edge)

leading_edge_metrics <- leading_edge_long |>
  count(go_id, leading_edge, name = "n_cell_types") |>
  group_by(go_id) |>
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
  left_join(leading_edge_metrics, by = "go_id") |>
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
# Priority-branch membership
# -----------------------------------------------------------------------------

if (!is.na(config$priority_file) && file.exists(config$priority_file)) {
  priority_branches <- readr::read_csv(
    config$priority_file,
    show_col_types = FALSE,
    progress = FALSE
  )
  missing_priority_columns <- setdiff(
    c("theme", "root_go_id", "max_granular"), names(priority_branches)
  )
  if (length(missing_priority_columns) > 0) {
    stop(
      "Priority file is missing: ",
      paste(missing_priority_columns, collapse = ", ")
    )
  }
} else {
  priority_branches <- as_tibble(default_priority_branches)
}

priority_branches <- priority_branches |>
  mutate(max_granular = as.integer(max_granular)) |>
  left_join(
    go_meta |> select(root_go_id = go_id, root_go_term = go_term),
    by = "root_go_id"
  )

if (any(is.na(priority_branches$theme) | priority_branches$theme == "") ||
    any(is.na(priority_branches$root_go_id)) ||
    any(is.na(priority_branches$max_granular) |
        priority_branches$max_granular < 1L |
        priority_branches$max_granular >
          config$hard_max_granular_per_cluster)) {
  stop(
    "Each priority row needs a theme, GO root, and max_granular between 1 and ",
    config$hard_max_granular_per_cluster
  )
}

if (any(is.na(priority_branches$root_go_term))) {
  warning(
    "Some priority roots are absent from this GO.db version; see ",
    "priority_branches_used.csv"
  )
}

readr::write_csv(
  priority_branches,
  file.path(config$output_dir, "priority_branches_used.csv")
)

assign_priority <- function(go_id) {
  lineage <- unique(c(go_id, candidate_ancestors[[go_id]]))
  hits <- priority_branches |>
    filter(root_go_id %in% lineage, !is.na(root_go_term))

  if (nrow(hits) == 0) {
    return(tibble(
      go_id = go_id,
      priority_branch = FALSE,
      priority_root_go_id = NA_character_,
      priority_theme = NA_character_,
      all_priority_themes = NA_character_,
      max_granular = config$default_max_granular_per_cluster
    ))
  }

  hit_depths <- map_dfr(hits$root_go_id, function(id) {
    d <- get_depth_range(id)
    tibble(root_go_id = id, root_depth = mean(d, na.rm = TRUE))
  })
  hits <- hits |>
    left_join(hit_depths, by = "root_go_id") |>
    arrange(desc(root_depth), root_go_id)
  chosen <- hits[1, ]

  tibble(
    go_id = go_id,
    priority_branch = TRUE,
    priority_root_go_id = chosen$root_go_id,
    priority_theme = chosen$theme,
    all_priority_themes = paste(sort(unique(hits$theme)), collapse = ";"),
    max_granular = chosen$max_granular
  )
}

priority_membership <- map_dfr(candidate_ids, assign_priority)
term_metrics <- term_metrics |>
  left_join(priority_membership, by = "go_id")

# Presentation score chooses representatives. It is not a statistical test.
term_metrics <- term_metrics |>
  mutate(
    capped_log_q = pmin(-log10(pmax(best_q, 1e-300)), 20),
    presentation_score =
      0.40 * rescale_01(capped_log_q) +
      0.25 * rescale_01(pmin(max_abs_nes, 3)) +
      0.15 * rescale_01(log1p(n_significant_cell_types)) +
      0.20 * rescale_01(information_content) +
      ifelse(priority_branch, 2, 0)
  )

# -----------------------------------------------------------------------------
# Semantic similarity and clusters
# -----------------------------------------------------------------------------

message(
  "Calculating semantic similarity for ", length(candidate_ids),
  " candidate GO terms..."
)

similarity_cache_file <- file.path(
  config$output_dir,
  paste0("semantic_similarity_", config$semantic_method, ".rds")
)

cache_is_valid <- FALSE
if (config$cache_similarity_matrix && file.exists(similarity_cache_file)) {
  cached <- readRDS(similarity_cache_file)
  cache_is_valid <- is.list(cached) &&
    identical(cached$ids, candidate_ids) &&
    identical(cached$method, config$semantic_method) &&
    identical(cached$go_db_version, as.character(packageVersion("GO.db"))) &&
    identical(
      cached$gosemsim_version, as.character(packageVersion("GOSemSim"))
    ) &&
    identical(
      cached$org_db_version, as.character(packageVersion("org.Hs.eg.db"))
    )
  if (cache_is_valid) {
    semantic_similarity <- cached$similarity
    message("Using cached semantic-similarity matrix")
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
    semantic_similarity <- GOSemSim::mgoSim(
      candidate_ids,
      candidate_ids,
      semData = semantic_data,
      measure = config$semantic_method,
      combine = NULL
    )
    semantic_similarity <- as.matrix(semantic_similarity)
  }
  semantic_similarity[!is.finite(semantic_similarity)] <- 0
  semantic_similarity <- (semantic_similarity + t(semantic_similarity)) / 2
  semantic_similarity <- pmin(pmax(semantic_similarity, 0), 1)
  diag(semantic_similarity) <- 1

  if (config$cache_similarity_matrix) {
    saveRDS(
      list(
        ids = candidate_ids,
        method = config$semantic_method,
        go_db_version = as.character(packageVersion("GO.db")),
        gosemsim_version = as.character(packageVersion("GOSemSim")),
        org_db_version = as.character(packageVersion("org.Hs.eg.db")),
        similarity = semantic_similarity
      ),
      similarity_cache_file,
      compress = FALSE
    )
  }
}

if (is.null(dim(semantic_similarity)) ||
    nrow(semantic_similarity) != length(candidate_ids) ||
    ncol(semantic_similarity) != length(candidate_ids)) {
  stop("Semantic-similarity matrix dimensions do not match candidate GO IDs")
}

if (is.null(rownames(semantic_similarity)) ||
    is.null(colnames(semantic_similarity))) {
  rownames(semantic_similarity) <- candidate_ids
  colnames(semantic_similarity) <- candidate_ids
} else {
  if (!setequal(rownames(semantic_similarity), candidate_ids) ||
      !setequal(colnames(semantic_similarity), candidate_ids)) {
    stop("Semantic-similarity matrix names do not match candidate GO IDs")
  }
  semantic_similarity <- semantic_similarity[candidate_ids, candidate_ids]
}

semantic_similarity[!is.finite(semantic_similarity)] <- 0
semantic_similarity <- (semantic_similarity + t(semantic_similarity)) / 2
semantic_similarity <- pmin(pmax(semantic_similarity, 0), 1)
diag(semantic_similarity) <- 1

representative_scores <- setNames(
  term_metrics$presentation_score,
  term_metrics$go_id
)
representative_scores[!is.finite(representative_scores)] <- 0

reduce_at_threshold <- function(threshold) {
  if (length(candidate_ids) == 1L) {
    return(tibble(
      go_id = candidate_ids,
      semantic_representative = candidate_ids
    ))
  }

  reduced <- rrvgo::reduceSimMatrix(
    semantic_similarity,
    scores = representative_scores,
    threshold = threshold,
    orgdb = "org.Hs.eg.db"
  ) |>
    as_tibble()

  if (!all(c("go", "parent") %in% names(reduced))) {
    stop("Unexpected rrvgo output: expected columns 'go' and 'parent'")
  }
  if (anyDuplicated(reduced$go) ||
      !setequal(as.character(reduced$go), candidate_ids)) {
    stop("rrvgo output does not contain exactly one row per candidate GO ID")
  }

  reduced |>
    transmute(
      go_id = as.character(go),
      semantic_representative = as.character(parent)
    )
}

threshold_labels <- gsub("\\.", "_", as.character(config$sensitivity_thresholds))
cluster_tables <- map(
  config$sensitivity_thresholds,
  reduce_at_threshold
)

semantic_clusters <- term_metrics |> select(go_id)
for (i in seq_along(cluster_tables)) {
  column_name <- paste0("representative_t", threshold_labels[[i]])
  semantic_clusters <- semantic_clusters |>
    left_join(
      cluster_tables[[i]] |>
        rename(!!column_name := semantic_representative),
      by = "go_id"
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
  stop("rrvgo returned missing representatives for one or more candidate terms")
}

# rrvgo defines the semantic clusters. Within each cluster, use a fixed
# direction-specific member as the comparison/display representative so that a
# positive-only term is never assessed against a negative-only representative.
semantic_clusters <- semantic_clusters |>
  rename(rrvgo_representative = semantic_representative) |>
  left_join(
    term_metrics |>
      select(
        go_id, cluster_direction = dominant_direction, presentation_score
      ),
    by = "go_id"
  ) |>
  group_by(rrvgo_representative, cluster_direction) |>
  arrange(desc(presentation_score), go_id, .by_group = TRUE) |>
  mutate(semantic_representative = first(go_id)) |>
  ungroup() |>
  select(-cluster_direction, -presentation_score)

# -----------------------------------------------------------------------------
# NES-profile and leading-edge similarity to the semantic representative
# -----------------------------------------------------------------------------

cell_types <- sort(unique(all_results_mapped$cell_type))

profile_rows <- all_results_mapped |>
  filter(go_id %in% candidate_ids) |>
  arrange(go_id, cell_type, is.na(evidence_q), evidence_q) |>
  group_by(go_id, cell_type) |>
  slice(1) |>
  ungroup() |>
  select(go_id, cell_type, NES)

nes_matrix <- matrix(
  NA_real_,
  nrow = length(candidate_ids),
  ncol = length(cell_types),
  dimnames = list(candidate_ids, cell_types)
)

row_index <- match(profile_rows$go_id, rownames(nes_matrix))
column_index <- match(profile_rows$cell_type, colnames(nes_matrix))
nes_matrix[cbind(row_index, column_index)] <- profile_rows$NES

evidence_rows_sorted <- evidence_rows |> arrange(go_id, cell_type)
evidence_by_term <- split(evidence_rows_sorted, evidence_rows_sorted$go_id)

leading_edge_by_term_cell <- lapply(evidence_by_term, function(df) {
  by_cell <- split(df$leading_edge, df$cell_type)
  lapply(by_cell, function(strings) {
    unique(unlist(lapply(strings, split_gene_string), use.names = FALSE))
  })
})

pair_profile_metrics <- function(go_id, representative) {
  if (is.na(representative) ||
      !go_id %in% rownames(nes_matrix) ||
      !representative %in% rownames(nes_matrix)) {
    return(c(
      shared_cell_types = 0,
      profile_correlation = NA_real_,
      direction_concordance = NA_real_
    ))
  }

  a <- nes_matrix[go_id, ]
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

pair_leading_edge_metrics <- function(go_id, representative) {
  a <- leading_edge_by_term_cell[[go_id]]
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
      go_id,
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
  semantic_clusters$go_id,
  semantic_clusters$semantic_representative,
  function(id, representative) {
    values <- pair_profile_metrics(id, representative)
    tibble(
      go_id = id,
      shared_cell_types = unname(values[["shared_cell_types"]]),
      profile_correlation = unname(values[["profile_correlation"]]),
      direction_concordance_to_representative = unname(
        values[["direction_concordance"]]
      )
    )
  }
)

leading_edge_pair_table <- map2_dfr(
  semantic_clusters$go_id,
  semantic_clusters$semantic_representative,
  function(id, representative) {
    values <- pair_leading_edge_metrics(id, representative)
    tibble(
      go_id = id,
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
  left_join(semantic_clusters, by = "go_id") |>
  left_join(profile_pair_table, by = "go_id") |>
  left_join(leading_edge_pair_table, by = "go_id") |>
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
  group_by(semantic_representative, priority_theme, direction_stratum) |>
  arrange(desc(priority_branch), desc(presentation_score), .by_group = TRUE) |>
  mutate(
    priority_rank_in_cluster = ifelse(
      priority_branch, cumsum(priority_branch), NA_integer_
    ),
    keep_granular = priority_branch & evidence_pass &
      (
        priority_rank_in_cluster <= max_granular |
        (
          (signal_distinct_from_representative | similarity_uncertain) &
            priority_rank_in_cluster <= config$hard_max_granular_per_cluster
        )
      ),
    keep_distinct_peripheral = !priority_branch & evidence_pass &
      (
        signal_distinct_from_representative |
        similarity_uncertain |
        !redundancy_supported
      )
  ) |>
  ungroup()

# -----------------------------------------------------------------------------
# Choose true common ancestors for peripheral semantic clusters
# -----------------------------------------------------------------------------

choose_common_ancestor <- function(member_ids) {
  member_ids <- unique(member_ids)
  lineages <- lapply(member_ids, get_ancestors)
  common <- Reduce(intersect, lineages)
  common <- setdiff(common, uninformative_parent_ids)
  if (length(common) == 0) return(NA_character_)

  candidate_table <- map_dfr(common, function(id) {
    depth <- get_depth_range(id)
    tibble(
      go_id = id,
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
      desc(information_content), go_id
    )

  if (nrow(candidate_table) == 0) NA_character_ else candidate_table$go_id[[1]]
}

cluster_ancestors <- term_metrics |>
  filter(!priority_branch, !keep_distinct_peripheral) |>
  group_by(semantic_representative, direction_stratum) |>
  summarise(
    cluster_common_ancestor = choose_common_ancestor(go_id),
    n_cluster_terms = n_distinct(go_id),
    .groups = "drop"
  )

# A separate fallback is calculated across every peripheral semantic cluster.
# It does not alter the audit/display mapping above; it is used only when a
# term lacks sufficient evidence to be shown granularly in a particular cell.
fallback_cluster_ancestors <- term_metrics |>
  filter(!priority_branch) |>
  group_by(semantic_representative, direction_stratum) |>
  summarise(
    fallback_common_ancestor = choose_common_ancestor(go_id),
    .groups = "drop"
  )

# If a cluster has no useful strict common ancestor, use an informative strict
# ancestor of the individual term before falling back to a semantic label.
individual_rollup_ancestors <- tibble(
  go_id = term_metrics$go_id,
  individual_rollup_ancestor = vapply(
    term_metrics$go_id,
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
    by = "go_id"
  )

display_mapping <- term_metrics |>
  mutate(
    display_go_id = case_when(
      keep_granular | keep_distinct_peripheral ~ go_id,
      priority_branch ~ priority_root_go_id,
      !is.na(cluster_common_ancestor) ~ cluster_common_ancestor,
      TRUE ~ semantic_representative
    ),
    mapping_type = case_when(
      keep_granular ~ "granular_self",
      keep_distinct_peripheral ~ "distinct_peripheral_self",
      priority_branch & go_id == priority_root_go_id ~ "priority_root_self",
      priority_branch ~ "priority_branch_rollup",
      !is.na(cluster_common_ancestor) ~ "ontology_ancestor",
      TRUE ~ "semantic_representative"
    ),
    granularity = ifelse(
      keep_granular | keep_distinct_peripheral, "granular", "high_level"
    ),
    coarse_display_go_id = case_when(
      priority_branch & !is.na(priority_root_go_id) ~ priority_root_go_id,
      !priority_branch & !is.na(fallback_common_ancestor) ~
        fallback_common_ancestor,
      !is.na(individual_rollup_ancestor) ~ individual_rollup_ancestor,
      TRUE ~ semantic_representative
    ),
    coarse_mapping_type = case_when(
      priority_branch & !is.na(priority_root_go_id) ~
        "priority_branch_rollup",
      !priority_branch & !is.na(fallback_common_ancestor) ~
        "ontology_ancestor",
      !is.na(individual_rollup_ancestor) ~ "ontology_ancestor",
      TRUE ~ "semantic_representative"
    ),
    coarse_rollup_source = case_when(
      priority_branch & !is.na(priority_root_go_id) ~ "priority_root",
      !priority_branch & !is.na(fallback_common_ancestor) ~
        "cluster_common_ancestor",
      !is.na(individual_rollup_ancestor) ~ "individual_ancestor",
      TRUE ~ "semantic_representative"
    ),
    coarse_is_strict_ancestor = map2_lgl(
      go_id,
      coarse_display_go_id,
      function(member_id, coarse_id) {
        !is.na(coarse_id) & coarse_id %in% get_ancestors(member_id)
      }
    ),
    theme_type = ifelse(
      priority_branch, "custom_fibrosis_priority", NA_character_
    ),
    rollup_basis = case_when(
      keep_granular ~ "priority_granular_retention",
      keep_distinct_peripheral ~ "distinct_or_uncertain_retention",
      !evidence_pass ~ "low_granular_evidence",
      redundancy_supported ~ "redundancy_supported",
      priority_branch ~ "priority_rank_cap",
      TRUE ~ "review_required"
    ),
    decision_reason = case_when(
      keep_granular & signal_distinct_from_representative ~
        "priority branch; evidence passed; distinct NES/leading-edge pattern",
      keep_granular ~
        "priority branch; evidence passed; top-ranked within semantic cluster",
      keep_distinct_peripheral & similarity_uncertain ~
        "peripheral term retained because redundancy evidence was insufficient",
      keep_distinct_peripheral ~
        "peripheral term retained because its NES or leading-edge pattern differed",
      priority_branch & go_id == priority_root_go_id ~
        "term is itself the configured high-level priority root",
      priority_branch & !evidence_pass ~
        "priority term rolled to its root because granular evidence did not pass",
      priority_branch & redundancy_supported ~
        "priority term was redundant and below the granular retention rank",
      priority_branch &
        (signal_distinct_from_representative | similarity_uncertain) ~
        "priority term was distinct or uncertain but exceeded the hard rank cap",
      priority_branch ~
        "priority term exceeded the granular retention rank",
      !evidence_pass & !is.na(cluster_common_ancestor) ~
        paste(
          "peripheral term with low granular evidence moved to a strict",
          "GO ancestor; this is a presentation roll-up, not redundancy evidence"
        ),
      !is.na(cluster_common_ancestor) ~
        "peripheral semantic cluster collapsed to a verified common GO ancestor",
      !evidence_pass ~
        paste(
          "low-granular-evidence peripheral term grouped by semantic label;",
          "review because no informative strict common ancestor was available"
        ),
      TRUE ~
        "no informative common GO ancestor; semantic representative used"
    )
  ) |>
  left_join(
    go_meta |> select(display_go_id = go_id, display_term = go_term),
    by = "display_go_id"
  ) |>
  left_join(
    go_meta |>
      select(coarse_display_go_id = go_id, coarse_display_term = go_term),
    by = "coarse_display_go_id"
  ) |>
  mutate(
    display_term = ifelse(
      is.na(display_term),
      paste0("Semantic group represented by ", semantic_representative),
      display_term
    ),
    coarse_display_term = ifelse(
      is.na(coarse_display_term),
      paste0("Semantic group represented by ", coarse_display_go_id),
      coarse_display_term
    ),
    display_key = paste(
      mapping_type,
      ifelse(is.na(priority_theme), "", priority_theme),
      direction_stratum,
      semantic_representative,
      display_go_id,
      sep = "|"
    )
  )

# Assert that every claimed ontology/priority roll-up is genuinely ancestral.
ancestry_check <- pmap_lgl(
  display_mapping |>
    select(go_id, display_go_id, mapping_type),
  function(go_id, display_go_id, mapping_type) {
    if (!mapping_type %in% c("ontology_ancestor", "priority_branch_rollup")) {
      return(TRUE)
    }
    display_go_id %in% get_ancestors(go_id)
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
    anyDuplicated(display_mapping$go_id)
) {
  stop("Display mapping must contain exactly one row per candidate GO term")
}

# -----------------------------------------------------------------------------
# Select one fixed tested source term for every display row
# -----------------------------------------------------------------------------

tested_go_ids <- all_results_mapped |>
  filter(!is.na(go_id), is.finite(NES)) |>
  distinct(go_id) |>
  pull(go_id)

display_sources <- display_mapping |>
  group_by(
    display_key, display_go_id, display_term, mapping_type, priority_theme,
    direction_stratum
  ) |>
  arrange(
    desc(n_tested_cell_types), desc(presentation_score), go_id,
    .by_group = TRUE
  ) |>
  summarise(
    source_go_id = if (first(display_go_id) %in% go_id &&
                       first(display_go_id) %in% tested_go_ids) {
      first(display_go_id)
    } else {
      first(go_id)
    },
    source_type = ifelse(
      first(display_go_id) %in% go_id &
        first(display_go_id) %in% tested_go_ids,
      "display_term",
      "representative_member"
    ),
    n_collapsed_terms = n_distinct(go_id),
    .groups = "drop"
  ) |>
  left_join(
    go_meta |> select(source_go_id = go_id, source_term = go_term),
    by = "source_go_id"
  ) |>
  left_join(
    term_metrics |>
      select(
        source_go_id = go_id,
        source_n_tested_cell_types = n_tested_cell_types
      ),
    by = "source_go_id"
  )

if (anyDuplicated(display_sources$display_key)) {
  stop("Each display key must have exactly one fixed source term")
}

display_mapping <- display_mapping |>
  left_join(
    display_sources |>
      select(
        display_key, source_go_id, source_term, source_type,
        source_n_tested_cell_types, n_collapsed_terms
      ),
    by = "display_key"
  )

# -----------------------------------------------------------------------------
# Output annotated data, mapping, summaries and plot
# -----------------------------------------------------------------------------

readr::write_csv(
  term_metrics |> arrange(desc(priority_branch), best_q),
  file.path(config$output_dir, "go_term_metrics.csv")
)

readr::write_csv(
  semantic_clusters |> arrange(semantic_representative, go_id),
  file.path(config$output_dir, "semantic_cluster_membership.csv")
)

mapping_output <- display_mapping |>
  select(
    pathway, go_id, go_term, min_depth, max_depth, information_content,
    dominant_direction, direction_consistency, direction_stratum,
    n_granular_evidence_cell_types,
    priority_branch, priority_theme, all_priority_themes, theme_type,
    priority_root_go_id, max_granular,
    rrvgo_representative, semantic_representative,
    semantic_similarity_to_representative,
    representative_stable, profile_correlation,
    direction_concordance_to_representative,
    median_leading_edge_jaccard, median_leading_edge_overlap,
    profile_similarity_available, leading_edge_similarity_available,
    profile_concordant, leading_edge_concordant, redundancy_supported,
    similarity_uncertain, signal_distinct_from_representative,
    evidence_pass, presentation_score, keep_granular,
    keep_distinct_peripheral,
    display_key, display_go_id, display_term,
    mapping_type, granularity, rollup_basis, decision_reason,
    cluster_common_ancestor, fallback_common_ancestor,
    individual_rollup_ancestor,
    coarse_display_go_id, coarse_display_term, coarse_mapping_type,
    coarse_rollup_source, coarse_is_strict_ancestor,
    source_go_id, source_term,
    source_type, source_n_tested_cell_types, n_collapsed_terms
  )

if (anyDuplicated(mapping_output$go_id)) {
  stop("The exported member-to-display mapping has duplicate GO IDs")
}

readr::write_csv(
  mapping_output |> arrange(display_term, go_term),
  file.path(config$output_dir, "adaptive_display_mapping.csv")
)

annotated_results <- all_results_mapped |>
  left_join(
    mapping_output |>
      select(
        go_id, priority_theme, theme_type, direction_stratum, display_key,
        display_go_id, display_term, mapping_type, granularity,
        rollup_basis, decision_reason, source_go_id, source_term, source_type,
        n_collapsed_terms, keep_granular, keep_distinct_peripheral,
        min_depth, max_depth, information_content, semantic_representative,
        priority_branch, max_granular, presentation_score,
        signal_distinct_from_representative,
        coarse_display_go_id, coarse_display_term, coarse_mapping_type,
        coarse_rollup_source, coarse_is_strict_ancestor
      ),
    by = "go_id"
  )

if (nrow(annotated_results) != nrow(all_results_mapped)) {
  stop("Annotating results changed the number of GSEA rows")
}

readr::write_csv(
  annotated_results,
  file.path(config$output_dir, "go_results_annotated.csv")
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
      go_id[evidence_q <= config$fdr_threshold]
    ),
    n_positive_members = n_distinct(
      go_id[evidence_q <= config$fdr_threshold & NES > 0]
    ),
    n_negative_members = n_distinct(
      go_id[evidence_q <= config$fdr_threshold & NES < 0]
    ),
    .groups = "drop"
  )

if (anyDuplicated(member_by_cell[c("display_key", "cell_type")])) {
  stop("Grouped member evidence is not unique by display key and cell type")
}

source_results <- all_results_mapped |>
  filter(!is.na(go_id)) |>
  inner_join(
    display_sources |>
      select(display_key, source_go_id, source_type),
    by = c("go_id" = "source_go_id")
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
  file.path(config$output_dir, "adaptive_go_summary_long.csv")
)

# -----------------------------------------------------------------------------
# Primary compact significant-term output
# -----------------------------------------------------------------------------

# The complete grid above is useful for matrices and audit, but it is not the
# final biological list. Build the final list from significant tested rows so
# every displayed term has explicit evidence provenance. A member is shown at
# its own GO ID only when it passes the row-level granular evidence rule and
# was selected for granular retention. Otherwise it is moved to the best
# available priority root or verified ontology ancestor; semantic grouping is
# used only when no informative strict ancestor exists.
final_member_rows_raw <- annotated_results |>
  filter(
    !is.na(go_id),
    !is.na(display_key),
    evidence_pass_row,
    is.finite(evidence_q),
    evidence_q <= config$fdr_threshold,
    is.finite(NES)
  )

# Different MSigDB pathway labels can occasionally resolve to the same GO ID.
# They are aliases of one ontology concept and must not consume multiple
# granular-selection slots. Keep the strongest tested row for each
# GO-ID x cell-type x comparison combination.
final_member_evidence <- final_member_rows_raw |>
  arrange(
    go_id,
    cell_type,
    comparison,
    evidence_q,
    desc(abs(NES)),
    desc(leading_edge_fraction),
    pathway
  ) |>
  distinct(go_id, cell_type, comparison, .keep_all = TRUE) |>
  mutate(
    member_direction = case_when(
      NES > 0 ~ "positive",
      NES < 0 ~ "negative",
      TRUE ~ "zero"
    ),
    granular_candidate_here = granular_evidence_pass_row &
      is.finite(min_depth) &
      min_depth >= config$final_granular_min_depth &
      (
        (
          coalesce(priority_branch, FALSE) &
            coalesce(keep_granular, FALSE)
        ) |
          (
            !coalesce(priority_branch, FALSE) &
              coalesce(keep_distinct_peripheral, FALSE) &
              evidence_q <= config$final_granular_fdr &
              abs(NES) >= config$final_granular_min_abs_nes &
              (
                go_id == semantic_representative |
                  coalesce(signal_distinct_from_representative, FALSE)
              )
          )
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
    priority_theme,
    go_id,
    go_term
  ) |>
  summarise(
    priority_branch = any(coalesce(priority_branch, FALSE)),
    max_granular = safe_max(max_granular),
    selection_max_depth = safe_max(max_depth),
    selection_information_content = safe_max(information_content),
    selection_presentation_score = safe_max(presentation_score),
    n_candidate_cell_types = n_distinct(cell_type),
    best_candidate_q_value = safe_min(evidence_q),
    largest_candidate_abs_NES = safe_max(abs(NES)),
    .groups = "drop"
  ) |>
  group_by(
    member_direction,
    semantic_representative,
    priority_theme
  ) |>
  arrange(
    desc(coalesce(priority_branch, FALSE)),
    desc(selection_max_depth),
    desc(selection_information_content),
    desc(n_candidate_cell_types),
    desc(selection_presentation_score),
    best_candidate_q_value,
    desc(largest_candidate_abs_NES),
    go_id,
    .by_group = TRUE
  ) |>
  mutate(
    global_granular_rank = row_number(),
    global_granular_limit = if (any(coalesce(priority_branch, FALSE))) {
      valid_limits <- max_granular[
        coalesce(priority_branch, FALSE) & is.finite(max_granular)
      ]
      if (length(valid_limits) == 0) {
        config$default_max_granular_per_cluster
      } else {
        max(valid_limits)
      }
    } else {
      config$final_max_peripheral_per_cluster
    },
    global_granular_selected =
      global_granular_rank <= global_granular_limit
  ) |>
  ungroup()

readr::write_csv(
  global_granular_selection |>
    arrange(
      member_direction,
      semantic_representative,
      priority_theme,
      global_granular_rank,
      go_id
    ),
  file.path(
    config$output_dir,
    "adaptive_go_global_granular_selection.csv"
  )
)

final_member_rows <- final_member_evidence |>
  left_join(
    global_granular_selection |>
      select(
        member_direction,
        semantic_representative,
        priority_theme,
        go_id,
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
      "priority_theme",
      "go_id"
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
  group_by(
    cell_type,
    member_direction,
    semantic_representative,
    priority_theme
  ) |>
  mutate(
    cluster_has_granular_display = any(granular_display_here),
    cluster_go_ids = collapse_unique_text(go_id),
    cluster_terms = collapse_unique_text(go_term)
  ) |>
  ungroup() |>
  # Once a cluster has an evidence-supported granular winner, do not also emit
  # a coarser or semantic label for the same cluster. Non-winning members stay
  # available in go_results_annotated.csv and are counted below for audit.
  filter(granular_display_here | !cluster_has_granular_display) |>
  mutate(
    final_display_go_id = case_when(
      granular_display_here ~ go_id,
      !is.na(coarse_display_go_id) ~ coarse_display_go_id,
      TRUE ~ display_go_id
    ),
    direct_display_evidence = go_id == final_display_go_id,
    final_mapping_type = case_when(
      granular_display_here & coalesce(keep_granular, FALSE) ~
        "granular_self",
      granular_display_here & coalesce(keep_distinct_peripheral, FALSE) ~
        "distinct_peripheral_self",
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
  filter(!final_display_go_id %in% uninformative_parent_ids) |>
  left_join(
    go_meta |>
      select(final_display_go_id = go_id, final_display_term = go_term),
    by = "final_display_go_id"
  ) |>
  mutate(
    final_display_term = ifelse(
      is.na(final_display_term),
      paste0("Semantic group represented by ", final_display_go_id),
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
    final_display_go_id,
    desc(direct_display_evidence),
    evidence_q,
    desc(abs(NES)),
    desc(leading_edge_fraction),
    go_id
  ) |>
  group_by(cell_type, final_display_go_id) |>
  slice_head(n = 1) |>
  ungroup() |>
  transmute(
    cell_type,
    display_go_id = final_display_go_id,
    evidence_source = ifelse(
      direct_display_evidence,
      "display_term",
      "supporting_member"
    ),
    evidence_go_id = go_id,
    evidence_term = go_term,
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
  group_by(cell_type, final_display_go_id, final_display_term) |>
  summarise(
    case = one_or_na(case),
    reference = one_or_na(reference),
    comparison = one_or_na(comparison),
    n_significant_supporting_terms = n_distinct(go_id),
    n_granular_supporting_terms = n_distinct(
      go_id[granular_display_here]
    ),
    n_positive_supporting_terms = n_distinct(
      go_id[member_direction == "positive"]
    ),
    n_negative_supporting_terms = n_distinct(
      go_id[member_direction == "negative"]
    ),
    best_member_q_value = safe_min(evidence_q),
    largest_member_abs_NES = safe_max(abs(NES)),
    supporting_go_ids = collapse_unique_text(go_id),
    supporting_terms = collapse_unique_text(go_term),
    cluster_go_ids = collapse_delimited_values(cluster_go_ids),
    cluster_terms = collapse_delimited_values(cluster_terms),
    semantic_cluster_ids = collapse_unique_text(semantic_representative),
    selected_granular_go_ids = collapse_unique_text(
      go_id[granular_display_here]
    ),
    mapping_types = collapse_unique_text(final_mapping_type),
    rollup_sources = collapse_unique_text(final_rollup_source),
    priority_themes = collapse_unique_text(priority_theme),
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
    display_go_id = final_display_go_id,
    display_term = final_display_term
  ) |>
  mutate(
    n_cluster_significant_terms = vapply(
      cluster_go_ids,
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
    by = c("cell_type", "display_go_id")
  )

# Count cluster members not represented elsewhere in the same cell. Comparing
# sets, rather than subtracting row counts, avoids labelling a second retained
# granular winner as "not displayed".
cell_displayed_id_sets <- cell_term_summary |>
  group_by(cell_type) |>
  summarise(
    all_displayed_go_ids_in_cell = collapse_unique_text(display_go_id),
    .groups = "drop"
  )

cell_term_summary <- cell_term_summary |>
  left_join(cell_displayed_id_sets, by = "cell_type") |>
  mutate(
    n_cluster_terms_not_displayed = map2_int(
      cluster_go_ids,
      all_displayed_go_ids_in_cell,
      function(cluster_ids, displayed_ids) {
        length(setdiff(
          split_delimited_values(cluster_ids),
          split_delimited_values(displayed_ids)
        ))
      }
    )
  ) |>
  select(-all_displayed_go_ids_in_cell)

display_depth_table <- map_dfr(
  sort(unique(cell_term_summary$display_go_id)),
  function(id) {
    depth <- get_depth_range(id)
    tibble(
      display_go_id = id,
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
  left_join(display_depth_table, by = "display_go_id")

if (
  anyDuplicated(cell_term_summary[c("cell_type", "display_go_id")]) ||
    any(!cell_term_summary$all_rollups_are_strict_ancestors) ||
    any(cell_term_summary$display_go_id %in% uninformative_parent_ids)
) {
  stop(
    "Compact GO rows failed uniqueness, ancestry or informative-label checks"
  )
}

# Remove a displayed parent when a more specific displayed descendant exists
# in the same cell population. Positive and negative support are combined into
# direction='mixed' where necessary, so this cannot create duplicate GO terms.
cell_hierarchy_flags <- cell_term_summary |>
  group_split(cell_type, .keep = TRUE) |>
  map_dfr(function(x) {
    get_hierarchy_pruning_flags(x$display_go_id) |>
      mutate(cell_type = x$cell_type[[1]], .before = 1)
  })

cell_term_summary_all <- cell_term_summary |>
  left_join(
    cell_hierarchy_flags,
    by = c("cell_type", "display_go_id")
  )

final_by_cell <- cell_term_summary_all |>
  filter(!hierarchy_parent_removed) |>
  arrange(
    cell_type,
    desc(resolution == "granular"),
    evidence_q_value,
    desc(largest_member_abs_NES),
    display_go_id
  )

if (
  anyDuplicated(final_by_cell[c("cell_type", "display_go_id")]) ||
    any(!is.finite(final_by_cell$best_member_q_value)) ||
    any(final_by_cell$best_member_q_value > config$fdr_threshold)
) {
  stop("Final by-cell GO list failed uniqueness or significance checks")
}

remaining_cell_parent_flags <- final_by_cell |>
  group_split(cell_type, .keep = TRUE) |>
  map_lgl(function(x) {
    any(get_hierarchy_pruning_flags(
      x$display_go_id
    )$hierarchy_parent_removed)
  })

if (any(remaining_cell_parent_flags)) {
  stop("A parent-child pair remains in the final by-cell GO list")
}

readr::write_csv(
  final_by_cell |>
    select(-hierarchy_parent_removed, -retained_descendant_go_ids),
  file.path(config$output_dir, "adaptive_go_final_by_cell.csv")
)

# Aggregate the cell-specific frontier to one row per GO term. This is the
# designated final result for reporting: exact GO IDs are unique, directions
# and cell populations are retained in explicit list/count fields, and a final
# global antichain removes any parent whose retained child appears elsewhere.
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
    display_go_id,
    desc(evidence_source == "display_term"),
    evidence_q_value,
    desc(abs(evidence_NES)),
    cell_type
  ) |>
  group_by(display_go_id) |>
  slice_head(n = 1) |>
  ungroup() |>
  transmute(
    display_go_id,
    representative_evidence_cell_type = cell_type,
    representative_evidence_source = evidence_source,
    representative_evidence_go_id = evidence_go_id,
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
  group_by(display_go_id, display_term) |>
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
    supporting_go_ids = collapse_delimited_values(supporting_go_ids),
    supporting_terms = collapse_delimited_values(supporting_terms),
    cluster_go_ids = collapse_delimited_values(cluster_go_ids),
    cluster_terms = collapse_delimited_values(cluster_terms),
    semantic_cluster_ids = collapse_delimited_values(semantic_cluster_ids),
    selected_granular_go_ids = collapse_delimited_values(
      selected_granular_go_ids
    ),
    priority_themes = collapse_delimited_values(priority_themes),
    mapping_types = collapse_delimited_values(mapping_types),
    rollup_sources = collapse_delimited_values(rollup_sources),
    best_member_q_value = safe_min(best_member_q_value),
    largest_member_abs_NES = safe_max(largest_member_abs_NES),
    cell_type_evidence = paste(sort(cell_evidence), collapse = ";"),
    .groups = "drop"
  ) |>
  mutate(
    n_unique_supporting_go_terms = vapply(
      supporting_go_ids,
      function(x) length(split_delimited_values(x)),
      integer(1)
    ),
    n_unique_cluster_go_terms = vapply(
      cluster_go_ids,
      function(x) length(split_delimited_values(x)),
      integer(1)
    )
  ) |>
  left_join(representative_global_evidence, by = "display_go_id")

global_hierarchy_flags <- get_hierarchy_pruning_flags(
  global_term_summary$display_go_id
)

global_term_summary_all <- global_term_summary |>
  left_join(global_hierarchy_flags, by = "display_go_id")

final_significant_terms <- global_term_summary_all |>
  filter(!hierarchy_parent_removed) |>
  select(-hierarchy_parent_removed, -retained_descendant_go_ids) |>
  arrange(
    desc(!is.na(priority_themes)),
    desc(resolution == "granular"),
    best_member_q_value,
    desc(largest_member_abs_NES),
    display_go_id
  )

if (
  anyDuplicated(final_significant_terms$display_go_id) ||
    any(get_hierarchy_pruning_flags(
      final_significant_terms$display_go_id
    )$hierarchy_parent_removed) ||
    any(final_significant_terms$n_significant_cell_types < 1L)
) {
  stop("Final unique GO list failed uniqueness, hierarchy or evidence checks")
}

readr::write_csv(
  final_significant_terms,
  file.path(config$output_dir, "adaptive_go_final_significant_terms.csv")
)

cell_pruning_audit <- cell_term_summary_all |>
  filter(hierarchy_parent_removed) |>
  transmute(
    pruning_scope = "within_cell_type",
    context = cell_type,
    removed_go_id = display_go_id,
    removed_term = display_term,
    retained_descendant_go_ids,
    removal_reason = paste(
      "A more specific significant display term connected by strict is_a",
      "or part_of ancestry was retained in the same cell population"
    ),
    n_significant_cell_types = 1L,
    significant_cell_types = cell_type,
    direction,
    supporting_go_ids,
    cell_type_evidence = NA_character_,
    best_member_q_value,
    largest_member_abs_NES
  )

global_pruning_audit <- global_term_summary_all |>
  filter(hierarchy_parent_removed) |>
  transmute(
    pruning_scope = "global_unique_list",
    context = "all_cell_types",
    removed_go_id = display_go_id,
    removed_term = display_term,
    retained_descendant_go_ids,
    removal_reason = paste(
      "A more specific significant display term connected by strict is_a",
      "or part_of ancestry was retained in the final cross-cell-type list"
    ),
    n_significant_cell_types,
    significant_cell_types,
    direction = overall_direction,
    supporting_go_ids,
    cell_type_evidence,
    best_member_q_value,
    largest_member_abs_NES
  )

hierarchy_pruning_audit <- bind_rows(
  cell_pruning_audit,
  global_pruning_audit
) |>
  arrange(pruning_scope, context, removed_term, removed_go_id)

readr::write_csv(
  hierarchy_pruning_audit,
  file.path(config$output_dir, "adaptive_go_hierarchy_pruning_audit.csv")
)

nes_matrix_output <- adaptive_summary |>
  mutate(
    display_label = paste0(
      display_term, " (", display_go_id, ")",
      ifelse(
        is.na(priority_theme) | priority_theme == "",
        "",
        paste0(" [", priority_theme, "]")
      ),
      " {", direction_stratum, "}"
    )
  ) |>
  select(
    display_key, display_go_id, display_term, mapping_type, priority_theme,
    direction_stratum, source_go_id, source_type, display_label, cell_type,
    source_NES
  ) |>
  tidyr::pivot_wider(names_from = cell_type, values_from = source_NES)

readr::write_csv(
  nes_matrix_output,
  file.path(config$output_dir, "adaptive_GO_NES_matrix.csv")
)

source_q_matrix_output <- adaptive_summary |>
  select(
    display_key, display_go_id, display_term, mapping_type, priority_theme,
    direction_stratum, source_go_id, source_type, cell_type, source_q_value
  ) |>
  tidyr::pivot_wider(names_from = cell_type, values_from = source_q_value)

readr::write_csv(
  source_q_matrix_output,
  file.path(config$output_dir, "adaptive_GO_source_FDR_matrix.csv")
)

collapse_qc <- mapping_output |>
  group_by(
    display_key, display_go_id, display_term, mapping_type, priority_theme,
    direction_stratum, source_go_id, source_term, source_type
  ) |>
  summarise(
    n_terms = n_distinct(go_id),
    median_semantic_similarity = safe_median(
      semantic_similarity_to_representative[
        go_id != semantic_representative
      ]
    ),
    median_profile_correlation = safe_median(
      profile_correlation[go_id != semantic_representative]
    ),
    median_direction_concordance = safe_median(
      direction_concordance_to_representative[
        go_id != semantic_representative
      ]
    ),
    group_median_leading_edge_jaccard = safe_median(
      median_leading_edge_jaccard[go_id != semantic_representative]
    ),
    group_median_leading_edge_overlap = safe_median(
      median_leading_edge_overlap[go_id != semantic_representative]
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
    "mapped_candidate_go_ids",
    "priority_candidate_terms",
    "granular_priority_terms_retained",
    "distinct_peripheral_terms_retained",
    "terms_rolled_up",
    "display_groups",
    "adaptive_summary_grid_rows",
    "display_groups_flagged_for_review",
    "mapped_significant_input_rows",
    "unique_significant_go_cell_rows",
    "significant_go_alias_rows_removed",
    "global_granular_go_direction_selections",
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
    sum(term_metrics$priority_branch, na.rm = TRUE),
    sum(display_mapping$keep_granular, na.rm = TRUE),
    sum(display_mapping$keep_distinct_peripheral, na.rm = TRUE),
    sum(
      !display_mapping$keep_granular &
        !display_mapping$keep_distinct_peripheral,
      na.rm = TRUE
    ),
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
        final_significant_terms$display_go_id
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
    "adaptive_go_final_significant_terms.csv",
    "adaptive_go_final_by_cell.csv",
    "adaptive_go_global_granular_selection.csv",
    "adaptive_go_hierarchy_pruning_audit.csv",
    "go_id_mapping_qc.csv",
    "unmapped_or_ambiguous_terms.csv",
    "candidate_mapping_failures.csv",
    "priority_branches_used.csv",
    "go_term_metrics.csv",
    "semantic_cluster_membership.csv",
    "adaptive_display_mapping.csv",
    "go_results_annotated.csv",
    "adaptive_go_summary_long.csv",
    "adaptive_GO_NES_matrix.csv",
    "adaptive_GO_source_FDR_matrix.csv",
    "collapse_qc.csv",
    "reduction_counts.csv",
    "adaptive_GO_dotplot.png"
  ),
  meaning = c(
    paste(
      "Primary reporting output: one row per unique significant displayed GO",
      "term across the analysis. Exact duplicates and strict is_a/part_of",
      "ancestor-descendant pairs are removed; granular semantic siblings are",
      "capped globally. Cell populations and directions are aggregated."
    ),
    paste(
      "Cell-specific significant GO frontier: one row per cell type and GO",
      "ID, with mixed directions combined and no within-cell strict is_a/",
      "part_of ancestor pairs. evidence_* columns identify the tested source."
    ),
    paste(
      "Global ranking within each semantic cluster and direction, showing",
      "which evidence-supported granular GO terms passed the cluster cap."
    ),
    paste(
      "Terms removed for strict is_a/part_of hierarchy redundancy, with the",
      "retained descendant IDs and the removed term's evidence provenance."
    ),
    "Pathway-to-GO-ID mapping and mapping status for every tested pathway.",
    "Terms that need mapping review; no such term enters the summary.",
    "Significant candidate pathways omitted because mapping was not usable.",
    "Custom fibrosis-priority roots used for this run.",
    "Depth range, IC, evidence, direction and redundancy metrics by GO term.",
    "rrvgo sensitivity clusters and direction-specific representatives.",
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
      "term; nonsignificant source results are faint. FDR is within GO-BP",
      "and cell type."
    )
  )
)

readr::write_csv(
  output_dictionary,
  file.path(config$output_dir, "output_dictionary.csv")
)

# Plot the highest-priority/evidence display rows while retaining the complete
# table in CSV outputs.
plot_selection <- adaptive_summary |>
  group_by(
    display_key, display_go_id, display_term, priority_theme,
    direction_stratum, mapping_type
  ) |>
  summarise(
    is_priority = any(!is.na(priority_theme) & priority_theme != ""),
    best_source_q = safe_min(source_q_value),
    n_source_significant_cell_types = sum(source_significant, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(is.finite(best_source_q)) |>
  arrange(
    desc(is_priority), best_source_q,
    desc(n_source_significant_cell_types)
  ) |>
  slice_head(n = config$max_plot_terms)

plot_label_table <- plot_selection |>
  mutate(
    display_label = paste0(
      display_term, " (", display_go_id, ")",
      ifelse(
        is.na(priority_theme) | priority_theme == "",
        "",
        paste0(" [", priority_theme, "]")
      ),
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
      name = "Within-cell-type\nGO FDR < 0.05"
    ) +
    scale_size_continuous(range = c(1, 6)) +
    labs(
      x = "Cell population",
      y = NULL,
      colour = "Source NES",
      size = "-log10(within-cell-type\nGO FDR)",
      shape = "Statistic source",
      caption = paste(
        "Colour, size and alpha use the same fixed tested source term.",
        "FDR is within GO-BP x cell type; grouped-member evidence is",
        "descriptive and reported separately."
      )
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )

  ggsave(
    file.path(config$output_dir, "adaptive_GO_dotplot.png"),
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
    "go_id_mapping_file", "priority_file", "fdr_column", "fdr_threshold",
    "min_abs_nes", "min_leading_edge_genes",
    "min_leading_edge_fraction", "semantic_method", "semantic_threshold",
    "sensitivity_thresholds", "min_shared_cell_types_for_correlation",
    "min_shared_significant_cell_types_for_leading_edge",
    "min_profile_correlation",
    "min_direction_concordance", "min_leading_edge_jaccard",
    "min_leading_edge_overlap", "min_rollup_depth",
    "target_rollup_depth", "target_rollup_information_content",
    "default_max_granular_per_cluster",
    "hard_max_granular_per_cluster",
    "final_granular_fdr", "final_granular_min_abs_nes",
    "final_granular_min_depth", "final_max_peripheral_per_cluster",
    "minimum_mapping_fraction",
    "mapping_fraction"
  ),
  value = as.character(c(
    config$significant_file,
    ifelse(is.na(config$all_results_file), "NA", config$all_results_file),
    profile_source,
    ifelse(is.na(config$go_id_mapping_file), "NA", config$go_id_mapping_file),
    ifelse(is.na(config$priority_file), "built_in", config$priority_file),
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
    config$default_max_granular_per_cluster,
    config$hard_max_granular_per_cluster,
    config$final_granular_fdr,
    config$final_granular_min_abs_nes,
    config$final_granular_min_depth,
    config$final_max_peripheral_per_cluster,
    config$minimum_mapping_fraction,
    mapping_fraction
  ))
)

readr::write_csv(
  parameter_table,
  file.path(config$output_dir, "run_parameters.csv")
)

input_files <- unique(c(
  config$significant_file,
  if (have_full_results) config$all_results_file else character(),
  if (!is.na(config$go_id_mapping_file) && file.exists(config$go_id_mapping_file)) {
    config$go_id_mapping_file
  } else character(),
  if (!is.na(config$priority_file) && file.exists(config$priority_file)) {
    config$priority_file
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

go_metadata <- tryCatch(
  AnnotationDbi::metadata(GO.db::GO.db),
  error = function(e) data.frame(name = "metadata_error", value = conditionMessage(e))
)
readr::write_csv(
  as_tibble(go_metadata),
  file.path(config$output_dir, "GO_db_metadata.csv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(config$output_dir, "sessionInfo.txt")
)

message(
  "Adaptive GO-BP summary complete: ", config$output_dir,
  "\nPrimary compact output: adaptive_go_final_significant_terms.csv",
  " (", nrow(final_significant_terms), " unique GO terms)",
  "\nCell-specific output: adaptive_go_final_by_cell.csv"
)
