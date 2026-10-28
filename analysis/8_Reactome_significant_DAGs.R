#!/usr/bin/env Rscript
# All significant Reactome terms, without adaptive roll-up or term selection.
# Offline: uses ReactomePathways.txt and ReactomePathwaysRelation.txt, or a
# snapshot exported by script6 or script8. No network calls or Bioconductor.
#
# Rscript 8_Reactome_significant_DAGs.R --help
# Rscript 8_Reactome_significant_DAGs.R --resources /path/to/resources/reactome
#
# R >= 4.1. Install once before transferring the R environment to the HPC:
# install.packages(c("readr", "dplyr", "igraph", "ggplot2", "jsonlite", "svglite"))
#
# Nodes: one significant Reactome ID per scope. Colour: arithmetic mean NES across
# significant cell populations (one value per Reactome ID/cell/contrast). Node AREA:
# number of significant cell populations. A small grey square is the structural
# Reactome root if it is not significant. All multiple-parent branches are retained.
# Edges point broader -> narrower. Solid: a direct parent-subpathway edge. Dashed:
# a real ontology path whose intermediate terms are absent from the significant
# set in THIS scope. One shortest witness path is exported per endpoint pair.
# Depth: longest path from Reactome root in the FULL parent-subpathway ontology, root=0.
# This is a DAG depth convention, not a measure of biological importance.

options(stringsAsFactors = FALSE, warn = 1)

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

default_project <- path.expand(paste0(
  "~/project-area/data/crohns_scrnaseq/10c_14n_analysis/",
  "pseudobulk_dea/pathway_enrichment"
))

parse_options <- function(args) {
  opt <- list(
    input = file.path(default_project, "reactome_gsea_significant.csv"),
    output = file.path(default_project, "reactome_significant_DAGs"),
    mapping = NA_character_, hierarchy_rds = NA_character_,
    resources = Sys.getenv("REACTOME_RESOURCES_DIR", unset = file.path(default_project, "resources", "reactome")),
    fdr_column = "q_value_within_collection_cell_type", fdr = 0.05,
    label_top = 60L, include_root = TRUE, self_test = FALSE, help = FALSE
  )
  index <- 1L
  while (index <= length(args)) {
    flag <- args[[index]]
    if (flag %in% c("--help", "-h", "--self-test", "--no-root-anchor")) {
      if (flag %in% c("--help", "-h")) opt$help <- TRUE
      if (flag == "--self-test") opt$self_test <- TRUE
      if (flag == "--no-root-anchor") opt$include_root <- FALSE
      index <- index + 1L
      next
    }
    key <- switch(flag,
      "--input" = "input", "--output" = "output", "--mapping" = "mapping",
      "--hierarchy-rds" = "hierarchy_rds", "--resources" = "resources", "--fdr-column" = "fdr_column",
      "--fdr" = "fdr", "--label-top" = "label_top", NA_character_
    )
    if (is.na(key) || index == length(args)) {
      stop("Unknown option or missing value: ", flag, ". Use --help.")
    }
    opt[[key]] <- args[[index + 1L]]
    index <- index + 2L
  }
  for (key in c("input", "output", "mapping", "hierarchy_rds", "resources")) {
    if (!is.na(opt[[key]])) opt[[key]] <- path.expand(opt[[key]])
  }
  opt$fdr <- suppressWarnings(as.numeric(opt$fdr))
  if (length(opt$fdr) != 1L || !is.finite(opt$fdr) ||
      opt$fdr <= 0 || opt$fdr > 1) stop("--fdr must be in (0, 1].")
  opt$label_top <- if (tolower(as.character(opt$label_top)) == "all") {
    Inf
  } else suppressWarnings(as.integer(opt$label_top))
  if (is.na(opt$label_top) || opt$label_top < 0) {
    stop("--label-top must be a non-negative integer or 'all'.")
  }
  opt
}

print_help <- function() {
  cat(paste(
    "Usage: Rscript 8_Reactome_significant_DAGs.R [options]",
    "--input FILE           Original significant Reactome GSEA CSV",
    "--output DIRECTORY     Output root (default: reactome_significant_DAGs)",
    "--mapping FILE         Optional pathway/gs_name + reactome_id/gs_exact_source map",
    "--hierarchy-rds FILE    Reuse a reactome_hierarchy_snapshot.rds from script6 or script8",
    "--resources DIRECTORY  Offline ReactomePathways.txt and ReactomePathwaysRelation.txt",
    "--fdr-column NAME      Default: q_value_within_collection_cell_type",
    "--fdr NUMBER           Default: 0.05 (existing FDR; no enrichment is rerun)",
    "--label-top N|all      Static labels only; default 60. All nodes are plotted.",
    "--no-root-anchor       Omit the grey structural Reactome root (may disconnect branches)",
    "--self-test            Run synthetic DAG/path tests and exit",
    "",
    "Outputs: one comparison folder, then whole_atlas/, compartments/<id>/,",
    "cell_types/<id>/. Each contains DAG.html, DAG.svg, DAG.pdf, nodes.csv,",
    "edges.csv and significant_cell_evidence.csv. HTML is offline and self-contained.",
    "No Reactome terms are pruned or rolled up. Missing ontology ancestors are compressed",
    "into dashed links, with one shortest witness path in edges.csv.",
    sep = "\n"
  ), "\n")
}

normalise_cell <- function(x) {
  x <- enc2utf8(as.character(x))
  x <- gsub("\u00a0", " ", x, fixed = TRUE)
  x <- gsub("[\u2018\u2019\u201a\u201b\u201c\u201d\u201e\u201f]", "", x)
  converted <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "")
  converted[is.na(converted)] <- x[is.na(converted)]
  tolower(trimws(gsub("[[:space:]]+", " ", converted)))
}

slug <- function(x) {
  x <- normalise_cell(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unnamed")
}

collapse_values <- function(x) paste(sort(unique(x[!is.na(x)])), collapse = ";")

normalise_reactome_name <- reactome_normalise_name
extract_reactome_id <- reactome_extract_id

compartments <- list(
  t_lymphocytes = list(label = "T-lymphocytes", cells = c(
    "CD4 T cell", "CD8 T cell", "Na\u00efve CD8 T cell", "Tfh", "NK T cell", "NK cell"
  )),
  b_lymphocytes = list(label = "B-lymphocytes", cells = c(
    "B cells_cl_8", "B cells_cl_20", "B cells_cl_21"
  )),
  plasma_cells = list(label = "Plasma cells", cells = c(
    "Plasma cells_cl_4", "Plasma cells_cl_16"
  )),
  endothelial_pericytes = list(label = "Endothelial cells & related pericytes",
    cells = c("Endothelial", "Pericyte")),
  myofibroblasts_stromal = list(label = "Myofibroblasts & stromal cells", cells = c(
    "Myofibroblast", "Mesenchymal", "ADAMDEC1+ stromal_cl_17",
    "ADAMDEC1+ stromal_cl_19", "NPY+ Stromal", "Transitional Stromal C3+"
  )),
  macrophages_monocytes = list(label = "Macrophages & monocytes", cells = c(
    "Macrophages_cl_34", "Macrophages_cl_14", "Macrophages_cl_2",
    "Macrophages_cl_3", "Mast cell", "Monocytes"
  ))
)

# ---------------------------------------------------------------------------
# Graph algorithms: these also run in --self-test without external hierarchy files.
# ---------------------------------------------------------------------------

calculate_depths <- function(terms, edges, root = "REACTOME_ROOT") {
  graph <- igraph::graph_from_data_frame(
    edges[, c("parent", "child")], directed = TRUE,
    vertices = data.frame(name = terms$reactome_id)
  )
  if (!igraph::is_dag(graph)) stop("The parent-subpathway ontology contains a cycle.")
  if (!root %in% terms$reactome_id) stop("Reactome root is missing from ontology.")
  parent_map <- split(edges$parent, edges$child)
  min_depth <- max_depth <- stats::setNames(rep(NA_integer_, nrow(terms)), terms$reactome_id)
  min_depth[[root]] <- max_depth[[root]] <- 0L
  for (id in igraph::as_ids(igraph::topo_sort(graph, mode = "out"))) {
    if (id == root) next
    parents <- unique(parent_map[[id]])
    valid <- parents[is.finite(min_depth[parents]) & is.finite(max_depth[parents])]
    if (length(valid)) {
      min_depth[[id]] <- min(min_depth[valid]) + 1L
      max_depth[[id]] <- max(max_depth[valid]) + 1L
    }
  }
  data.frame(reactome_id = terms$reactome_id,
    depth_min = unname(min_depth[terms$reactome_id]),
    depth_max = unname(max_depth[terms$reactome_id]))
}

# Every upward branch stops at its FIRST visible ancestor. Memoisation is
# scope-specific: a node that is hidden in one scope may be significant in another.
compress_reactome_edges <- function(visible_ids, ontology_edges, significant_ids) {
  visible_ids <- unique(visible_ids)
  parent_rows <- split(ontology_edges, ontology_edges$child)
  memo <- new.env(parent = emptyenv())
  visiting <- new.env(parent = emptyenv())
  choose_route <- function(routes, key, proposed) {
    existing <- routes[[key]]
    if (is.null(existing) || length(proposed$ids) < length(existing$ids) ||
        (length(proposed$ids) == length(existing$ids) &&
         paste(proposed$ids, collapse = ";") < paste(existing$ids, collapse = ";"))) {
      routes[[key]] <- proposed
    }
    routes
  }
  frontier <- function(id) {
    if (id %in% visible_ids) {
      return(stats::setNames(list(list(ids = id, relations = character())), id))
    }
    if (exists(id, memo, inherits = FALSE)) return(get(id, memo))
    if (exists(id, visiting, inherits = FALSE)) stop("Cycle during ancestor traversal.")
    assign(id, TRUE, visiting)
    answer <- list()
    rows <- parent_rows[[id]]
    if (!is.null(rows)) for (j in seq_len(nrow(rows))) {
      routes <- frontier(rows$parent[[j]])
      for (ancestor in names(routes)) {
        route <- routes[[ancestor]]
        answer <- choose_route(answer, ancestor, list(
          ids = c(route$ids, id),
          relations = c(route$relations, rows$relationship[[j]])
        ))
      }
    }
    rm(list = id, envir = visiting)
    assign(id, answer, memo)
    answer
  }
  records <- list()
  index <- 0L
  for (child in visible_ids) {
    routes <- list()
    rows <- parent_rows[[child]]
    if (!is.null(rows)) for (j in seq_len(nrow(rows))) {
      candidates <- frontier(rows$parent[[j]])
      for (ancestor in names(candidates)) {
        route <- candidates[[ancestor]]
        routes <- choose_route(routes, ancestor, list(
          ids = c(route$ids, child),
          relations = c(route$relations, rows$relationship[[j]])
        ))
      }
    }
    for (parent in names(routes)) {
      route <- routes[[parent]]
      hops <- length(route$ids) - 1L
      hidden <- if (hops > 1L) route$ids[2L:hops] else character()
      if (any(hidden %in% visible_ids)) stop("A compressed edge crosses a visible term.")
      index <- index + 1L
      records[[index]] <- data.frame(
        parent_reactome_id = parent, child_reactome_id = child,
        edge_type = if (hops == 1L) "direct" else "via_nonsignificant_ancestors",
        line_style = if (hops == 1L && parent %in% significant_ids &&
          child %in% significant_ids) "solid" else "dashed",
        ontology_hops = hops, n_hidden_intermediates = length(hidden),
        hidden_reactome_ids = paste(hidden, collapse = ";"),
        one_shortest_path_reactome_ids = paste(route$ids, collapse = ";"),
        one_shortest_path_relationships = paste(route$relations, collapse = ";"),
        structural_root_link = !parent %in% significant_ids,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(records)) return(data.frame(
    parent_reactome_id = character(), child_reactome_id = character(), edge_type = character(),
    line_style = character(), ontology_hops = integer(), n_hidden_intermediates = integer(),
    hidden_reactome_ids = character(), one_shortest_path_reactome_ids = character(),
    one_shortest_path_relationships = character(), structural_root_link = logical()
  ))
  do.call(rbind, records)
}

run_self_tests <- function() {
  edge <- function(p, c) data.frame(parent = p, child = c, relationship = "has_subpathway")
  check <- function(edges, visible, expected) {
    result <- compress_reactome_edges(visible, edges, visible)
    actual <- paste(result$parent_reactome_id, result$child_reactome_id, result$line_style, sep = ":")
    stopifnot(setequal(actual, expected))
    invisible(result)
  }
  chain <- edge(c("A", "B"), c("B", "C"))
  check(chain, c("A", "B", "C"), c("A:B:solid", "B:C:solid"))
  bridge <- check(chain, c("A", "C"), "A:C:dashed")
  stopifnot(bridge$ontology_hops == 2L, bridge$hidden_reactome_ids == "B")
  check(edge(c("A", "X", "B", "Y"), c("X", "C", "Y", "C")),
    c("A", "B", "C"), c("A:C:dashed", "B:C:dashed"))
  check(edge(c("A", "A", "X"), c("C", "X", "C")),
    c("A", "C"), "A:C:solid")
  check(edge(c("A", "B", "A", "X"), c("B", "C", "X", "C")),
    c("A", "B", "C"), c("A:B:solid", "B:C:solid", "A:C:dashed"))
  check(edge(c("R", "R"), c("A", "B")), c("A", "B"), character())
  anchored <- compress_reactome_edges(c("R", "A", "B"),
    edge(c("R", "R"), c("A", "B")), c("A", "B"))
  stopifnot(all(anchored$line_style == "dashed"))
  depth_edges <- edge(c("R", "R", "A", "B", "C"), c("A", "C", "B", "C", "D"))
  depths <- calculate_depths(data.frame(reactome_id = c("R", "A", "B", "C", "D")),
    depth_edges, root = "R")
  stopifnot(depths$depth_min[depths$reactome_id == "C"] == 1L,
    depths$depth_max[depths$reactome_id == "C"] == 3L,
    depths$depth_max[depths$reactome_id == "D"] == 4L)
  cat("All synthetic multiple-parent, hidden-path, root and full-depth tests passed.\n")
}

# HTML_HELPER_BEGIN
write_dag_html <- function(
  nodes, edges, file, title, colour_limit, size_limit, scope_note
) {
  payload <- list(
    title = title, note = scope_note, colourLimit = colour_limit,
    sizeLimit = size_limit, nodes = nodes, edges = edges
  )
  payload_json <- jsonlite::toJSON(
    payload, dataframe = "rows", auto_unbox = TRUE, na = "null",
    null = "null", digits = 12
  )
  # Prevent a source label from closing the script element containing JSON.
  payload_json <- gsub("<", "\\u003c", payload_json, fixed = TRUE)
  template <- r"---(<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Reactome significant-term DAG</title>
<style>
*{box-sizing:border-box}body{margin:0;background:#f6f8fb;color:#172536;font:14px system-ui,-apple-system,sans-serif}header{padding:17px 22px 12px;background:white;border-bottom:1px solid #dce2eb}h1{font-size:21px;margin:0 0 5px;line-height:1.25}#note{color:#526176;font-size:12px;max-width:1350px;white-space:pre-wrap}.toolbar{display:flex;flex-wrap:wrap;gap:7px;align-items:center;padding:11px 18px;background:#fff;border-bottom:1px solid #dce2eb}input{width:min(370px,65vw);padding:8px;border:1px solid #bbc7d4;border-radius:5px;font:inherit}button{border:1px solid #bbc7d4;border-radius:5px;background:#fff;padding:8px 11px;color:#21364e;cursor:pointer}button:hover{background:#edf3fa}button:focus-visible,input:focus-visible{outline:2px solid #3b79b4;outline-offset:2px}.layout{display:grid;grid-template-columns:minmax(0,1fr) 320px;height:calc(100vh - 182px);min-height:470px}#view{position:relative;overflow:hidden;background:white;touch-action:none}canvas{display:block;width:100%;height:100%;cursor:grab}canvas.dragging{cursor:grabbing}aside{padding:15px 16px;border-left:1px solid #dce2eb;overflow:auto;background:#fcfdff}h2{font-size:15px;margin:0 0 10px}#details{font-size:12px;white-space:pre-wrap;line-height:1.5;overflow-wrap:anywhere;margin-bottom:20px}.legend{font-size:12px;line-height:1.5;color:#445268}.gradient{height:13px;border:1px solid #ccd4df;background:linear-gradient(to right,#2166ac,#fff,#b2182b);margin-top:7px}.range{display:flex;justify-content:space-between;font-size:11px}.linekey{display:inline-block;width:31px;border-top:2px solid #68798c;margin-right:8px;vertical-align:middle}.dashed{border-top-style:dashed}.rootkey{display:inline-block;background:#b8c1cb;width:9px;height:9px;margin-right:9px}.hint{position:absolute;bottom:10px;left:12px;padding:5px 8px;background:#ffffffdc;color:#4c5d72;border-radius:4px;font-size:11px;pointer-events:none}#tooltip{position:absolute;display:none;max-width:330px;pointer-events:none;white-space:pre-wrap;background:#12283eef;color:white;border-radius:5px;padding:9px 11px;font-size:12px;line-height:1.4;box-shadow:0 3px 15px #0003}#status{font-size:12px;color:#526176;margin-left:6px}#sizes{display:flex;align-items:end;gap:17px;min-height:42px;margin:8px 0}.sizeitem{display:flex;flex-direction:column;align-items:center;gap:2px;font-size:10px}.circlekey{display:block;border-radius:100%;border:1px solid #4d637a;background:#d6e3f1}p{margin:8px 0}.small{font-size:11px;color:#607086}@media(max-width:850px){.layout{grid-template-columns:minmax(0,1fr);height:auto}#view{height:65vh;min-height:430px}aside{border-left:0;border-top:1px solid #dce2eb}.toolbar{padding:9px}.hint{max-width:70%}}
</style></head><body>
<header><h1 id="title"></h1><div id="note"></div></header>
<div class="toolbar"><input id="search" type="search" placeholder="Find Reactome ID or term name" aria-label="Find Reactome ID or term name"><button id="find">Find</button><button id="fit">Fit all</button><button id="plus" aria-label="Zoom in">+</button><button id="minus" aria-label="Zoom out">−</button><button id="clear">Clear selection</button><button id="labels" aria-pressed="false">Labels: auto</button><button id="download">Download PNG view</button><span id="status" role="status"></span></div>
<div class="layout"><div id="view"><canvas id="graph" tabindex="0" aria-label="Reactome hierarchy graph; use search to select terms, arrow keys to pan and plus or minus to zoom"></canvas><div id="tooltip"></div><div class="hint">Scroll to zoom · drag to pan · click a term or link · double-click to fit · arrows show broader → narrower</div></div>
<aside><h2>Selection</h2><div id="details">Choose a node or link to inspect its evidence or ontology path.</div><div class="legend"><h2>How to read this graph</h2><div>Node colour: mean NES</div><div class="gradient"></div><div class="range"><span id="low"></span><span>0</span><span id="high"></span></div><p>Mean NES includes significant cell populations only. Positive and negative values can cancel; inspect their counts in the selection panel.</p><p>Node area: number of significant cell populations</p><div id="sizes"></div><p>All significant nodes have count 1 in a single-cell-type graph.</p><p><span class="linekey"></span>Direct parent/subpathway connection</p><p><span class="linekey dashed"></span>Connection through hidden non-significant ancestors, or a structural root connection</p><p><span class="rootkey"></span>Grey square: synthetic all-Reactome root used only as an ontology anchor if it is not significant.</p><p>Layers use the longest root-to-term path in the full ontology. Root depth = 0. Terms can have several parents.</p><p class="small">Every significant mapped term is present. Search and zoom expose dense regions. This file is self-contained and works offline.</p></div></aside></div>
<script id="dag-data" type="application/json">__DAG_PAYLOAD__</script>
<script>
(() => {
'use strict';
const D=JSON.parse(document.getElementById('dag-data').textContent), nodes=D.nodes||[], rawEdges=D.edges||[];
document.getElementById('title').textContent=D.title;document.title=D.title;
document.getElementById('note').textContent=D.note;
const cv=document.getElementById('graph'),view=document.getElementById('view'),ctx=cv.getContext('2d'),tooltip=document.getElementById('tooltip'),detail=document.getElementById('details'),status=document.getElementById('status');
const lim=Math.max(1e-8,+D.colourLimit||1),sizeLim=Math.max(1,+D.sizeLimit||1);
document.getElementById('low').textContent=String(-lim.toFixed(2));document.getElementById('high').textContent=lim.toFixed(2);
const legendCounts=[...new Set([1,Math.max(1,Math.round(sizeLim/2)),sizeLim])];
legendCounts.forEach(n=>{const item=document.createElement('div');item.className='sizeitem';const dot=document.createElement('span');dot.className='circlekey';const diam=24*Math.sqrt(n/sizeLim);dot.style.width=diam+'px';dot.style.height=diam+'px';const lab=document.createElement('span');lab.textContent=String(n);item.append(dot,lab);document.getElementById('sizes').append(item);});
const byId=new Map(nodes.map(n=>[n.reactome_id,n]));
const depthCounts=new Map();nodes.forEach(n=>depthCounts.set(+n.depth_max,(depthCounts.get(+n.depth_max)||0)+1));
const xs=nodes.map(n=>+n.x||0),ys=nodes.map(n=>+n.y||0),minX=Math.min(...xs,0),maxX=Math.max(...xs,1),maxY=Math.max(...ys,0),minY=Math.min(...ys,0);
const worldWidth=Math.max(1700,Math.max(1,...depthCounts.values())*24),worldHeight=Math.max(550,(maxY-minY)*105);
nodes.forEach(n=>{n.wx=80+((+n.x||0)-minX)/Math.max(1,maxX-minX)*worldWidth;n.wy=65+(maxY-(+n.y||0))/Math.max(1,maxY-minY)*worldHeight;});
const edges=rawEdges.map(e=>({...e,a:byId.get(e.parent_reactome_id),b:byId.get(e.child_reactome_id)})).filter(e=>e.a&&e.b);
const adjacent=new Map(nodes.map(n=>[n.reactome_id,new Set()]));edges.forEach(e=>{adjacent.get(e.a.reactome_id).add(e.b.reactome_id);adjacent.get(e.b.reactome_id).add(e.a.reactome_id);});
const sigCount=nodes.filter(n=>n.is_significant).length;
let width=1,height=1,k=1,tx=0,ty=0,selected=null,selectedEdge=null,hover=null,searchMatches=new Set(),allLabels=false,drag=null,needsDraw=false;
function num(v,d=3){return v==null||!Number.isFinite(+v)?'NA':(+v).toPrecision(d);}
function nesColour(v){if(v==null||!Number.isFinite(+v))return '#b8c1cb';const z=Math.max(-1,Math.min(1,+v/lim)),t=Math.abs(z),a=[255,255,255],b=z<0?[33,102,172]:[178,24,43];return 'rgb('+a.map((x,i)=>Math.round(x+(b[i]-x)*t)).join(',')+')';}
function radius(n){return n.is_significant?12*Math.sqrt(Math.max(0,+n.n_significant_cell_types)/sizeLim)*Math.max(.35,Math.min(2.5,k)):3*Math.max(.65,Math.min(2,k));}
function screen(n){return [n.wx*k+tx,n.wy*k+ty];}
function schedule(){if(!needsDraw){needsDraw=true;requestAnimationFrame(()=>{needsDraw=false;draw();});}}
function fit(targetNodes=nodes){if(!targetNodes.length)return;const x=targetNodes.map(n=>n.wx),y=targetNodes.map(n=>n.wy),loX=Math.min(...x),hiX=Math.max(...x),loY=Math.min(...y),hiY=Math.max(...y);k=Math.min((width-95)/Math.max(160,hiX-loX),(height-75)/Math.max(160,hiY-loY),3);k=Math.max(.015,k);tx=width/2-(hiX+loX)/2*k;ty=height/2-(hiY+loY)/2*k;schedule();}
function zoom(f,cx=width/2,cy=height/2){const next=Math.max(.012,Math.min(9,k*f));tx=cx-(cx-tx)*next/k;ty=cy-(cy-ty)*next/k;k=next;schedule();}
function neighbours(){const s=new Set();if(selected){s.add(selected.reactome_id);adjacent.get(selected.reactome_id).forEach(id=>s.add(id));}return s;}
function edgePath(e){const [ax,ay]=screen(e.a),[bx,by]=screen(e.b);return {ax,ay,bx,by};}
function draw(){
 ctx.clearRect(0,0,width,height);ctx.fillStyle='#fff';ctx.fillRect(0,0,width,height);
 const near=neighbours(),active=!!selected||!!selectedEdge;
 ctx.font='10px system-ui';ctx.textAlign='left';
 const layers=new Map();nodes.forEach(n=>{if(!layers.has(+n.depth_max))layers.set(+n.depth_max,n);});
 [...layers].sort((a,b)=>a[0]-b[0]).forEach(([depth,n])=>{const y=screen(n)[1];if(y<0||y>height)return;ctx.strokeStyle='#edf0f5';ctx.lineWidth=1;ctx.setLineDash([]);ctx.beginPath();ctx.moveTo(0,y);ctx.lineTo(width,y);ctx.stroke();ctx.fillStyle='#6c7d92';ctx.fillText('Depth '+depth,7,Math.max(11,y-5));});
 for(const e of edges){const p=edgePath(e);if(Math.max(p.ax,p.bx)<0||Math.min(p.ax,p.bx)>width||Math.max(p.ay,p.by)<0||Math.min(p.ay,p.by)>height)continue;const isActive=e===selectedEdge||(selected&&(e.a===selected||e.b===selected));ctx.globalAlpha=active?(isActive?.85:.045):.22;ctx.strokeStyle=isActive?'#263e59':'#6f8094';ctx.lineWidth=isActive?1.8:.65;ctx.setLineDash(e.line_style==='dashed'?[5,4]:[]);ctx.beginPath();ctx.moveTo(p.ax,p.ay);ctx.lineTo(p.bx,p.by);ctx.stroke();if(isActive||k>.65){const dx=p.bx-p.ax,dy=p.by-p.ay,len=Math.hypot(dx,dy);if(len>13){const r=radius(e.b)+2,bx=p.bx-dx/len*r,by=p.by-dy/len*r;ctx.setLineDash([]);ctx.fillStyle=ctx.strokeStyle;ctx.beginPath();ctx.moveTo(bx,by);ctx.lineTo(bx-dx/len*6-dy/len*3,by-dy/len*6+dx/len*3);ctx.lineTo(bx-dx/len*6+dy/len*3,by-dy/len*6-dx/len*3);ctx.closePath();ctx.fill();}}}
 ctx.setLineDash([]);
 for(const n of nodes){const [x,y]=screen(n),r=radius(n);if(x<-30||x>width+30||y<-30||y>height+30)continue;ctx.globalAlpha=active?(selectedEdge?(n===selectedEdge.a||n===selectedEdge.b?1:.19):(near.has(n.reactome_id)?1:.19)):1;ctx.fillStyle=nesColour(n.mean_NES);ctx.strokeStyle=searchMatches.has(n.reactome_id)?'#00a67d':n===selected?'#071c35':'#364e67';ctx.lineWidth=searchMatches.has(n.reactome_id)||n===selected?2.3:.65;ctx.beginPath();if(n.is_significant)ctx.arc(x,y,Math.max(.15,r),0,Math.PI*2);else ctx.rect(x-r,y-r,r*2,r*2);ctx.fill();ctx.stroke();if(allLabels||n===selected||searchMatches.has(n.reactome_id)||(k>1.3&&nodes.length<450)){ctx.font=(n===selected?'bold ':'')+'11px system-ui';ctx.fillStyle='#1b3048';ctx.textAlign='left';const label=n.reactome_id+' '+n.reactome_term;ctx.fillText(label.length>95?label.slice(0,92)+'…':label,x+r+4,y+3);}}
 ctx.globalAlpha=1;status.textContent=sigCount.toLocaleString()+' significant nodes · '+edges.length.toLocaleString()+' links · '+Math.round(k*100)+'% zoom'+(searchMatches.size?' · '+searchMatches.size+' matches':'');
}
function nodeAt(x,y){let best=null,dist=Infinity;for(const n of nodes){const [sx,sy]=screen(n),d=Math.hypot(x-sx,y-sy);if(d<=Math.max(5,radius(n)+2)&&d<dist){best=n;dist=d;}}return best;}
function edgeAt(x,y){let best=null,dist=6;for(const e of edges){const p=edgePath(e),dx=p.bx-p.ax,dy=p.by-p.ay,ll=dx*dx+dy*dy;if(ll===0)continue;const t=Math.max(0,Math.min(1,((x-p.ax)*dx+(y-p.ay)*dy)/ll));if(t<.025||t>.975)continue;const d=Math.hypot(x-p.ax-t*dx,y-p.ay-t*dy);if(d<dist){dist=d;best=e;}}return best;}
function nodeText(n){if(!n.is_significant)return n.reactome_id+' — '+n.reactome_term+'\nStructural root anchor; no significance or NES assigned.\nFull-ontology depth: '+n.depth_max;return n.reactome_id+' — '+n.reactome_term+'\n\nMean NES: '+num(n.mean_NES)+'\nSignificant cell populations: '+n.n_significant_cell_types+'\nPositive / negative NES cells: '+n.n_positive_cells+' / '+n.n_negative_cells+'\nMean |NES|: '+num(n.mean_abs_NES)+'\nNES range: '+num(n.min_NES)+' to '+num(n.max_NES)+'\nBest source FDR: '+num(n.best_q)+'\nOntology depth min / max: '+n.depth_min+' / '+n.depth_max+'\n\nSignificant populations:\n'+String(n.significant_cell_types||'').split(';').join('\n');}
function edgeText(e){return (e.line_style==='solid'?'Solid direct ontology link':'Dashed ontology connection')+'\n'+e.a.reactome_id+' → '+e.b.reactome_id+'\n\nOntology edges in witness path: '+e.ontology_hops+'\nPath (broader → narrower):\n'+String(e.one_shortest_path_reactome_ids).split(';').join(' → ')+'\n\nRelationships in path:\n'+String(e.one_shortest_path_relationships).split(';').join(' → ')+(e.hidden_reactome_ids?'\n\nHidden, non-significant intermediate terms:\n'+String(e.hidden_reactome_ids).split(';').join('\n'):'')+(e.structural_root_link?'\n\nThe parent is a structural root anchor.':'')+'\n\nThe exported edge table retains one shortest valid witness path per endpoint pair.';}
function choose(n,e){selected=n||null;selectedEdge=e||null;detail.textContent=n?nodeText(n):e?edgeText(e):'Choose a node or link to inspect its evidence or ontology path.';schedule();}
cv.addEventListener('pointerdown',ev=>{const r=cv.getBoundingClientRect();drag={x:ev.clientX-r.left,y:ev.clientY-r.top,tx,ty,moved:false};cv.setPointerCapture(ev.pointerId);cv.classList.add('dragging');tooltip.style.display='none';});
cv.addEventListener('pointermove',ev=>{const r=cv.getBoundingClientRect(),x=ev.clientX-r.left,y=ev.clientY-r.top;if(drag){if(Math.hypot(x-drag.x,y-drag.y)>3)drag.moved=true;tx=drag.tx+x-drag.x;ty=drag.ty+y-drag.y;schedule();return;}const n=nodeAt(x,y),e=n?null:edgeAt(x,y);hover=n||e;if(!hover){tooltip.style.display='none';return;}tooltip.textContent=n?(n.reactome_id+' — '+n.reactome_term+'\n'+(n.is_significant?'Mean NES '+num(n.mean_NES)+' · '+n.n_significant_cell_types+' significant cell populations':'Structural Reactome root anchor')+'\nClick for details'):(e.line_style==='solid'?'Direct Reactome edge':'Dashed connection')+' · '+e.ontology_hops+' ontology edges\n'+e.a.reactome_id+' → '+e.b.reactome_id+'\nClick for hidden path';tooltip.style.display='block';tooltip.style.left=Math.min(x+12,Math.max(5,width-345))+'px';tooltip.style.top=Math.min(y+12,height-tooltip.offsetHeight-8)+'px';});
cv.addEventListener('pointerup',ev=>{if(drag&&!drag.moved){const r=cv.getBoundingClientRect(),x=ev.clientX-r.left,y=ev.clientY-r.top,n=nodeAt(x,y);choose(n,n?null:edgeAt(x,y));}drag=null;cv.classList.remove('dragging');});
cv.addEventListener('pointercancel',()=>{drag=null;cv.classList.remove('dragging');});
cv.addEventListener('pointerleave',()=>{tooltip.style.display='none';});
cv.addEventListener('wheel',ev=>{ev.preventDefault();const r=cv.getBoundingClientRect();zoom(Math.exp(-ev.deltaY*.0015),ev.clientX-r.left,ev.clientY-r.top);},{passive:false});
cv.addEventListener('dblclick',()=>fit());
cv.addEventListener('keydown',ev=>{if(['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','+','-','=','Escape'].includes(ev.key))ev.preventDefault();if(ev.key==='ArrowLeft')tx+=45;if(ev.key==='ArrowRight')tx-=45;if(ev.key==='ArrowUp')ty+=45;if(ev.key==='ArrowDown')ty-=45;if(ev.key==='+'||ev.key==='=')zoom(1.35);if(ev.key==='-')zoom(1/1.35);if(ev.key==='Escape')choose(null,null);schedule();});
function find(){const q=document.getElementById('search').value.trim().toLowerCase();searchMatches=new Set();if(!q){choose(null,null);schedule();return;}const matches=nodes.filter(n=>(n.reactome_id+' '+n.reactome_term).toLowerCase().includes(q));matches.forEach(n=>searchMatches.add(n.reactome_id));if(matches.length===1){choose(matches[0],null);fit([matches[0],...[...adjacent.get(matches[0].reactome_id)].map(id=>byId.get(id))]);}else if(matches.length){choose(null,null);fit(matches);detail.textContent=matches.length+' matching terms.\n\n'+matches.slice(0,80).map(n=>n.reactome_id+' — '+n.reactome_term).join('\n')+(matches.length>80?'\n…':'');}else{detail.textContent='No matching Reactome ID or term name.';}schedule();}
document.getElementById('find').onclick=find;document.getElementById('search').addEventListener('keydown',ev=>{if(ev.key==='Enter')find();});
document.getElementById('fit').onclick=()=>fit();document.getElementById('plus').onclick=()=>zoom(1.4);document.getElementById('minus').onclick=()=>zoom(1/1.4);
document.getElementById('clear').onclick=()=>{searchMatches.clear();document.getElementById('search').value='';choose(null,null);};
document.getElementById('labels').onclick=ev=>{allLabels=!allLabels;ev.currentTarget.textContent=allLabels?'Labels: all':'Labels: auto';ev.currentTarget.setAttribute('aria-pressed',String(allLabels));schedule();};
document.getElementById('download').onclick=()=>{draw();const a=document.createElement('a');a.download='Reactome_DAG_view.png';a.href=cv.toDataURL('image/png');a.click();};
function resize(){const bounds=view.getBoundingClientRect();width=Math.max(1,bounds.width);height=Math.max(1,bounds.height);const ratio=window.devicePixelRatio||1;cv.width=Math.round(width*ratio);cv.height=Math.round(height*ratio);ctx.setTransform(ratio,0,0,ratio,0,0);schedule();}
new ResizeObserver(resize).observe(view);resize();fit();
})();
</script></body></html>)---"
  marker <- "__DAG_PAYLOAD__"
  marker_position <- regexpr(marker, template, fixed = TRUE)[[1]]
  before <- substr(template, 1L, marker_position - 1L)
  after <- substr(template, marker_position + nchar(marker), nchar(template))
  writeLines(paste0(before, payload_json, after), con = file, useBytes = TRUE)
  invisible(file)
}
# HTML_HELPER_END

# MAIN_PIPELINE_BEGIN

load_ontology <- function(snapshot_file, resources_dir) {
  h <- load_reactome_hierarchy(resources_dir = resources_dir, snapshot_file = snapshot_file)
  h$depths <- dplyr::rename(h$depths, depth_min = min_depth, depth_max = max_depth)
  h
}

map_pathways <- function(input, ontology, mapping_file) {
  pathways <- sort(unique(input$pathway))
  exact <- data.frame(pathway = character(), reactome_id = character())
  input_id_columns <- intersect(c("reactome_id", "REACTOME_ID", "gs_exact_source"), names(input))
  if (length(input_id_columns)) for (column in input_id_columns) {
    exact <- dplyr::bind_rows(exact, data.frame(
      pathway = input$pathway, reactome_id = extract_reactome_id(input[[column]])
    ))
  }
  embedded <- split(exact$reactome_id, exact$pathway)
  supplied <- list()
  if (!is.na(mapping_file)) {
    if (!file.exists(mapping_file)) stop("Mapping file does not exist: ", mapping_file)
    table <- readr::read_csv(mapping_file, show_col_types = FALSE, progress = FALSE)
    path_col <- intersect(c("pathway", "gs_name"), names(table))[1L]
    id_col <- intersect(c("reactome_id", "gs_exact_source", "REACTOME_ID"), names(table))[1L]
    if (is.na(path_col) || is.na(id_col)) {
      stop("Mapping file needs pathway/gs_name and reactome_id/gs_exact_source columns.")
    }
    map_ids <- extract_reactome_id(table[[id_col]])
    # Preserve ambiguity recorded by a previous mapping audit.
    if ("matched_reactome_ids" %in% names(table)) {
      problematic <- if ("n_reactome_matches" %in% names(table)) {
        !is.na(table$n_reactome_matches) & table$n_reactome_matches > 1L
      } else rep(FALSE, nrow(table))
      for (i in which(problematic)) {
        ids <- unlist(strsplit(as.character(table$matched_reactome_ids[[i]]), ";", fixed = TRUE))
        exact <- dplyr::bind_rows(exact, data.frame(
          pathway = as.character(table[[path_col]][[i]]), reactome_id = extract_reactome_id(ids)
        ))
      }
    }
    supplied <- split(map_ids, as.character(table[[path_col]]))
  }
  name_lookup <- split(ontology$terms$reactome_id, normalise_reactome_name(ontology$terms$reactome_term))
  records <- lapply(pathways, function(pathway) {
    embedded_ids <- unique(stats::na.omit(c(embedded[[pathway]],
      exact$reactome_id[exact$pathway == pathway])))
    supplied_ids <- unique(stats::na.omit(supplied[[pathway]]))
    literal_id <- extract_reactome_id(pathway)
    if (!is.na(literal_id)) embedded_ids <- unique(c(embedded_ids, literal_id))
    if (length(embedded_ids)) {
      ids <- embedded_ids
      method <- "embedded_or_literal_Reactome_ID"
    } else if (length(supplied_ids)) {
      ids <- supplied_ids
      method <- "supplied_mapping"
    } else {
      ids <- unique(name_lookup[[pathway]])
      if (!length(ids)) ids <- unique(name_lookup[[normalise_reactome_name(pathway)]])
      method <- "normalised_Reactome_db_term_name"
    }
    ids <- sort(unique(ids[!is.na(ids)]))
    status <- if (!length(ids)) "unmapped" else if (length(ids) > 1L) {
      "ambiguous"
    } else if (!ids %in% ontology$terms$reactome_id) "invalid_or_non_human_pathway" else "mapped"
    data.frame(pathway = pathway, reactome_id = if (length(ids) == 1L) ids else NA_character_,
      n_reactome_matches = length(ids), matched_reactome_ids = paste(ids, collapse = ";"),
      mapping_method = method, mapping_status = status)
  })
  dplyr::bind_rows(records) |>
    dplyr::left_join(ontology$terms, by = "reactome_id") |>
    dplyr::left_join(ontology$depths, by = "reactome_id") |>
    dplyr::mutate(mapping_status = dplyr::if_else(
      mapping_status == "mapped" & !is.finite(depth_max),
      "unreachable_from_Reactome_root", mapping_status
    ))
}

make_scopes <- function(cells) {
  cells <- sort(unique(cells))
  keys <- normalise_cell(cells)
  if (anyDuplicated(keys)) stop("Distinct input cell labels have ambiguous normalised names.")
  scopes <- list(list(id = "whole_atlas", label = "Whole atlas", type = "atlas",
    cells = cells, requested = cells, missing = character(), directory = "whole_atlas"))
  for (id in names(compartments)) {
    spec <- compartments[[id]]
    matched <- match(normalise_cell(spec$cells), keys)
    scopes[[length(scopes) + 1L]] <- list(id = id, label = spec$label,
      type = "compartment", cells = cells[matched[!is.na(matched)]],
      requested = spec$cells, missing = spec$cells[is.na(matched)],
      directory = file.path("compartments", id))
  }
  cell_slugs <- make.unique(slug(cells), sep = "_")
  for (j in seq_along(cells)) scopes[[length(scopes) + 1L]] <- list(
    id = paste0("cell_", cell_slugs[[j]]), label = cells[[j]], type = "cell_type",
    cells = cells[[j]], requested = cells[[j]], missing = character(),
    directory = file.path("cell_types", cell_slugs[[j]])
  )
  scopes
}

summarise_nodes <- function(evidence, ontology, scope) {
  evidence |>
    dplyr::group_by(reactome_id) |>
    dplyr::summarise(
      mean_NES = mean(NES), mean_abs_NES = mean(abs(NES)),
      min_NES = min(NES), max_NES = max(NES),
      sd_NES = if (dplyr::n() > 1L) stats::sd(NES) else NA_real_,
      n_significant_cell_types = dplyr::n_distinct(cell_type),
      n_positive_cells = sum(NES > 0), n_negative_cells = sum(NES < 0),
      best_q = min(evidence_q),
      significant_cell_types = collapse_values(cell_type),
      significant_pathway_labels = collapse_values(pathway),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      positive_fraction = n_positive_cells / n_significant_cell_types,
      mixed_direction = n_positive_cells > 0L & n_negative_cells > 0L,
      is_significant = TRUE,
      scope_id = scope$id, scope_label = scope$label,
      scope_type = scope$type,
      n_scope_cell_types_in_significant_input = length(scope$cells)
    ) |>
    dplyr::left_join(ontology$terms, by = "reactome_id") |>
    dplyr::left_join(ontology$depths, by = "reactome_id") |>
    dplyr::arrange(reactome_id)
}

layout_nodes <- function(nodes, edges) {
  if (!nrow(nodes)) return(nodes)
  graph <- igraph::graph_from_data_frame(
    edges[, c("parent_reactome_id", "child_reactome_id")], directed = TRUE,
    vertices = data.frame(name = nodes$reactome_id)
  )
  if (!igraph::is_dag(graph)) stop("Compressed display graph is not a DAG.")
  if (nrow(edges)) {
    parent_depth <- nodes$depth_max[match(edges$parent_reactome_id, nodes$reactome_id)]
    child_depth <- nodes$depth_max[match(edges$child_reactome_id, nodes$reactome_id)]
    if (any(parent_depth >= child_depth)) stop("An ontology edge violates maximum-depth ordering.")
  }
  if (nrow(nodes) < 2L || nrow(edges) == 0L) {
    nodes$x <- seq_len(nrow(nodes))
  } else {
    positions <- tryCatch(
      igraph::layout_with_sugiyama(graph, layers = nodes$depth_max, maxiter = 60L)$layout,
      error = function(e) {
        warning("Sugiyama ordering failed; using deterministic per-depth ordering: ",
          conditionMessage(e))
        NULL
      }
    )
    nodes$x <- if (is.null(positions)) {
      ave(seq_len(nrow(nodes)), nodes$depth_max, FUN = seq_along)
    } else positions[, 1L]
  }
  # Depth is fixed from the full ontology; layout only controls horizontal order.
  nodes$y <- -nodes$depth_max
  nodes
}

write_static_dag <- function(nodes, edges, directory, title, colour_limit,
                             size_limit, label_top, scope_note) {
  plot_edges <- edges |>
    dplyr::left_join(nodes |> dplyr::select(parent_reactome_id = reactome_id, x, y),
      by = "parent_reactome_id") |>
    dplyr::left_join(nodes |> dplyr::select(child_reactome_id = reactome_id, xend = x, yend = y),
      by = "child_reactome_id")
  labelled <- nodes |>
    dplyr::filter(is_significant) |>
    dplyr::arrange(best_q, dplyr::desc(n_significant_cell_types), dplyr::desc(abs(mean_NES)))
  labelled <- if (is.infinite(label_top)) labelled else utils::head(labelled, label_top)
  labelled$plot_label <- paste0(labelled$reactome_id, " ",
    vapply(labelled$reactome_term, function(x) paste(strwrap(x, 32), collapse = "\n"), character(1)))
  depths <- seq.int(min(nodes$depth_max), max(nodes$depth_max))
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = -depths, colour = "#e7e9ec", linewidth = 0.25) +
    ggplot2::geom_segment(data = plot_edges,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linetype = line_style),
      colour = "#667585", alpha = 0.24, linewidth = 0.23,
      arrow = grid::arrow(length = grid::unit(0.035, "inches"), type = "closed")) +
    ggplot2::geom_point(data = nodes |> dplyr::filter(is_significant),
      ggplot2::aes(x = x, y = y, fill = mean_NES, size = n_significant_cell_types),
      shape = 21, colour = "#303945", stroke = 0.16, alpha = 0.94) +
    ggplot2::geom_point(data = nodes |> dplyr::filter(!is_significant),
      ggplot2::aes(x = x, y = y), shape = 22, size = 2.6,
      fill = "#c6ccd3", colour = "#536171") +
    ggplot2::geom_text(data = labelled,
      ggplot2::aes(x = x, y = y, label = plot_label),
      size = 1.6, vjust = -0.7, lineheight = 0.9, check_overlap = TRUE) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
      midpoint = 0, limits = c(-colour_limit, colour_limit), name = "Mean NES\n(significant cells)") +
    ggplot2::scale_size_area(max_size = 6.5, limits = c(0, size_limit),
      breaks = unique(c(1, pretty(c(1, size_limit), n = 4))), name = "Significant\ncell types") +
    ggplot2::scale_linetype_manual(values = c(solid = "solid", dashed = "dashed"),
      breaks = c("solid", "dashed"), labels = c("Direct Reactome edge", "Hidden ancestors / root link"),
      name = "Connection") +
    ggplot2::scale_y_continuous(breaks = -depths, labels = depths,
      expand = ggplot2::expansion(mult = c(0.05, 0.08))) +
    ggplot2::scale_x_continuous(breaks = NULL, expand = ggplot2::expansion(mult = 0.025)) +
    ggplot2::labs(title = title,
      subtitle = sprintf("%s significant Reactome terms; %s displayed connections. %s",
        sum(nodes$is_significant), nrow(edges), scope_note),
      x = NULL, y = "Full Reactome depth (longest root path)",
      caption = paste(
        "All mapped significant terms are plotted; static labels are limited for legibility.",
        "Use DAG.html to search, zoom and inspect every term and connection.",
        "Arrows run broader to narrower. Grey square: structural Reactome root; no NES assigned.",
        "White can reflect opposing NES values across cells; inspect positive/negative counts."
      )) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(size = 7, hjust = 0),
      axis.text.y = ggplot2::element_text(size = 8), legend.position = "bottom")
  width <- 25
  height <- max(12, min(23, 5 + length(depths) * 0.58))
  ggplot2::ggsave(file.path(directory, "DAG.pdf"), plot = p,
    device = grDevices::pdf, width = width, height = height, limitsize = FALSE,
    bg = "white", useDingbats = FALSE)
  ggplot2::ggsave(file.path(directory, "DAG.svg"), plot = p,
    device = svglite::svglite, width = width, height = height, limitsize = FALSE,
    bg = "white")
  invisible(p)
}

write_readme <- function(file) {
  writeLines(c(
    "SIGNIFICANT Reactome DAGS",
    "",
    "INPUT AND SCOPES",
    "Uses the original significant GSEA results, with the chosen existing FDR <= threshold.",
    "It does not use adaptive roll-up outputs or select a subset of significant Reactome terms.",
    "Separate graphs are written for every case/reference comparison, the whole atlas,",
    "each requested compartment and each cell type occurring in the significant input.",
    "Compartment members match by case, whitespace and accent normalisation only.",
    "Missing requested cell labels are reported. CD8 T mem is never substituted for naive CD8.",
    "With a significant-only input, absence cannot distinguish non-significance from not tested.",
    "",
    "NODE COLOUR AND AREA",
    "mean_NES is the arithmetic mean over significant cell types in the current scope.",
    "There is one contribution per Reactome ID/cell type/contrast, so cells receive equal weight.",
    "Conflicting NES values for aliases of the same Reactome pathway/cell/contrast cause an explicit error.",
    "Node AREA represents n_significant_cell_types. All cell-type graph counts equal 1.",
    "One symmetric colour scale and one count-to-area scale are used across the entire run.",
    "Positive and negative NES can cancel: a white node may have strong opposing effects.",
    "Inspect n_positive_cells, n_negative_cells, mean_abs_NES, min_NES and max_NES.",
    "The displayed mean is conditional on significance; it is not a pooled GSEA or meta-analysis.",
    "",
    "ONTOLOGY AND CONNECTIONS",
    "Only Reactome pathway/subpathway relationships are used; multiple parents are retained.",
    "Arrows are drawn from broader parent to narrower child for reading down the page.",
    "Depth is the longest path from All Reactome pathways (REACTOME_ROOT, depth 0) in the full ontology.",
    "The ontology is a DAG, so min and max root depth can differ; both are exported.",
    "Layer number is descriptive and is not a biological importance or information-content score.",
    "For each significant child, walk upward along every branch and stop at its first visible ancestor.",
    "Solid line: one direct ontology edge between two significant terms in this scope.",
    "Dashed line: at least one intermediate ancestor is absent from the significant set in this scope.",
    "The grey structural-root square is an optional structural anchor and has no significance or NES.",
    "Links to that non-significant root are also dashed, even when directly connected.",
    "A path is stopped at any visible intermediate term. Alternate parent branches are preserved.",
    "No generic transitive reduction is applied; a direct edge takes precedence for identical endpoints.",
    "edges.csv gives one shortest witness path, its relation sequence, hops and hidden Reactome IDs.",
    "It does not enumerate every possible ontology path, which can be very numerous.",
    "With --no-root-anchor, valid disconnected components may remain; no sibling edges are invented.",
    "",
    "FILES",
    "index.html: links to every completed atlas, compartment and cell-type DAG.",
    "DAG.html: offline canvas viewer with search, hover, selection, pan/zoom and depth guides.",
    "DAG.pdf / DAG.svg: all-node static overview; only top labels are shown by default.",
    "nodes.csv: every plotted node, mean/count metrics, source cell labels, depth and coordinates.",
    "edges.csv: every plotted connection with auditable original ontology paths.",
    "significant_cell_evidence.csv: the per-cell significant input values behind the means/counts.",
    "scope_metadata.csv: requested/matched/missing populations and plotting scales.",
    "scope_run_manifest.csv: completion/skip/failure and counts for every graph.",
    "reactome_id_mapping_audit.csv / omitted_input_rows.csv: all mapping and filtering exclusions.",
    "reactome_hierarchy_snapshot.rds: exact ontology used, reusable with --hierarchy-rds on an offline machine.",
    "ontology_depths.csv / ontology_edges.csv: the full Reactome depth/reference-edge tables.",
    "input_checksums.csv / run_parameters.csv / sessionInfo.txt: reproducibility records.",
    "",
    "ONTOLOGY SOURCE",
    "Use a offline Reactome hierarchy or snapshot compatible with the MSigDB release when possible; renamed/obsolete",
    "labels are audited rather than assigned speculative Reactome IDs.",
    "Hierarchy documentation: https://reactome.org/2-uncategorised/61-pathway-information"
  ), file, useBytes = TRUE)
}

write_graph_index <- function(manifest, output_root) {
  escape_html <- function(x) {
    x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    gsub('"', "&quot;", x, fixed = TRUE)
  }
  rows <- vapply(seq_len(nrow(manifest)), function(i) {
    row <- manifest[i, ]
    relative <- substring(row$html, nchar(output_root) + 2L)
    label <- escape_html(row$scope_label)
    if (grepl("^completed", row$status)) {
      label <- paste0('<a href="', escape_html(relative), '">', label, '</a>')
    }
    paste0("<tr><td>", escape_html(row$comparison), "</td><td>",
      escape_html(row$scope_type), "</td><td>", label, "</td><td>",
      row$n_significant_Reactome_terms, "</td><td>", escape_html(row$status),
      "</td><td>", escape_html(row$missing_cell_types_in_significant_input), "</td></tr>")
  }, character(1))
  writeLines(c(
    '<!doctype html><html lang="en"><meta charset="utf-8"><title>Significant Reactome DAGs</title>',
    '<style>body{font:14px system-ui;color:#243447;max-width:1500px;margin:35px auto;padding:0 22px}table{border-collapse:collapse;width:100%}td,th{text-align:left;padding:9px;border-bottom:1px solid #dce3e9}th{background:#eef3f7}a{color:#17639b}p{line-height:1.6}</style>',
    '<h1>Significant Reactome Biological Process DAGs</h1>',
    '<p>Open a scope below. Every mapped significant term is included. Colour is mean NES across significant cell populations; node area is their count. Layers show maximum full-ontology depth. Dashed links pass through hidden ancestors or join the structural Reactome root. The viewers work offline.</p>',
    '<table><thead><tr><th>Comparison</th><th>Scope</th><th>Graph</th><th>Significant Reactome terms</th><th>Status</th><th>Requested populations absent from significant input</th></tr></thead><tbody>',
    rows, '</tbody></table><p>See README.txt and scope_run_manifest.csv for interpretation and audit details.</p></html>'
  ), file.path(output_root, "index.html"), useBytes = TRUE)
}

main <- function(opt) {
  if (!file.exists(opt$input)) stop("Input CSV does not exist: ", opt$input)
  dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
  opt$output <- normalizePath(opt$output, mustWork = TRUE)
  write_readme(file.path(opt$output, "README.txt"))
  input <- readr::read_csv(opt$input, show_col_types = FALSE, progress = FALSE)
  required <- c("pathway", "cell_type", "NES", "case", "reference", opt$fdr_column)
  missing <- setdiff(required, names(input))
  if (length(missing)) stop("Input is missing: ", paste(missing, collapse = ", "))
  if ("collection" %in% names(input)) {
    wrong <- !is.na(input$collection) & input$collection != "Reactome"
    if (any(wrong)) stop("Input contains non-Reactome collections. Supply Reactome results only.")
  }
  for (column in c("pathway", "cell_type", "case", "reference")) {
    if (any(is.na(input[[column]]) | !nzchar(trimws(as.character(input[[column]]))))) {
      stop("Missing identity values in column: ", column)
    }
  }
  input$NES <- suppressWarnings(as.numeric(input$NES))
  input$evidence_q <- suppressWarnings(as.numeric(input[[opt$fdr_column]]))
  if (any(is.finite(input$evidence_q) & (input$evidence_q < 0 | input$evidence_q > 1))) {
    stop("FDR values must be between zero and one.")
  }
  if (!"comparison" %in% names(input)) {
    input$comparison <- paste(input$case, "vs", input$reference)
  }
  if (any(is.na(input$comparison) | !nzchar(trimws(input$comparison)))) {
    stop("The comparison column contains missing labels.")
  }
  contrasts <- input |>
    dplyr::distinct(case, reference, comparison) |>
    dplyr::arrange(case, reference, comparison) |>
    dplyr::mutate(contrast_id = sprintf("c%02d_%s", dplyr::row_number(),
      slug(paste(case, "vs", reference))))
  input <- input |>
    dplyr::left_join(contrasts, by = c("case", "reference", "comparison")) |>
    dplyr::mutate(input_row = dplyr::row_number(),
      passes_input_filter = is.finite(NES) & is.finite(evidence_q) & evidence_q <= opt$fdr)
  filtered <- input |> dplyr::filter(passes_input_filter)
  if (!nrow(filtered)) stop("No finite significant rows pass the chosen FDR threshold.")
  readr::write_csv(contrasts |>
    dplyr::filter(!contrast_id %in% filtered$contrast_id),
    file.path(opt$output, "comparisons_without_significant_rows.csv"))
  contrasts <- contrasts |> dplyr::filter(contrast_id %in% filtered$contrast_id)

  message("Loading the offline Reactome ontology...")
  ontology <- load_ontology(opt$hierarchy_rds, opt$resources)
  saveRDS(ontology, file.path(opt$output, "reactome_hierarchy_snapshot.rds"), compress = TRUE)
  readr::write_csv(ontology$edges, file.path(opt$output, "ontology_edges.csv"))
  readr::write_csv(ontology$terms |> dplyr::left_join(ontology$depths, by = "reactome_id"),
    file.path(opt$output, "ontology_depths.csv"))

  mapping <- map_pathways(filtered, ontology, opt$mapping)
  readr::write_csv(mapping, file.path(opt$output, "reactome_id_mapping_audit.csv"))
  readr::write_csv(mapping |> dplyr::filter(mapping_status != "mapped"),
    file.path(opt$output, "unmapped_or_ambiguous_terms.csv"))
  map_columns <- mapping |> dplyr::select(pathway, mapped_reactome_id = reactome_id,
    mapping_status, mapping_method)
  annotated <- input |> dplyr::left_join(map_columns, by = "pathway")
  omitted <- annotated |>
    dplyr::filter(!passes_input_filter | is.na(mapping_status) | mapping_status != "mapped") |>
    dplyr::mutate(omission_reason = dplyr::case_when(
      !passes_input_filter ~ "nonfinite_NES_or_FDR_or_did_not_pass_threshold",
      TRUE ~ mapping_status
    ))
  readr::write_csv(omitted, file.path(opt$output, "omitted_input_rows.csv"))
  evidence <- annotated |>
    dplyr::filter(passes_input_filter, mapping_status == "mapped") |>
    dplyr::mutate(reactome_id = mapped_reactome_id) |>
    dplyr::select(-mapped_reactome_id)
  if (!nrow(evidence)) stop("No significant rows have usable current human Reactome IDs.")
  n_omitted_terms <- sum(mapping$mapping_status != "mapped")
  if (n_omitted_terms) warning(n_omitted_terms,
    " pathway labels have no unique, reachable human pathway mapping; see mapping audit.")

  aliases <- evidence |>
    dplyr::group_by(contrast_id, reactome_id, cell_type) |>
    dplyr::summarise(n_rows = dplyr::n(), NES_range = max(NES) - min(NES),
      pathway_labels = collapse_values(pathway), .groups = "drop") |>
    dplyr::filter(n_rows > 1L)
  readr::write_csv(aliases, file.path(opt$output, "duplicate_Reactome_cell_audit.csv"))
  if (any(aliases$NES_range > 1e-8)) stop(
    "Different NES values map to the same Reactome ID/cell/contrast. Resolve conflicting ",
    "aliases before averaging; see duplicate_Reactome_cell_audit.csv."
  )
  # Identical-NES duplicates contribute only once; choose smallest q deterministically.
  evidence <- evidence |>
    dplyr::arrange(evidence_q, pathway, input_row) |>
    dplyr::distinct(contrast_id, reactome_id, cell_type, .keep_all = TRUE)
  colour_limit <- max(abs(evidence$NES))
  if (!is.finite(colour_limit) || colour_limit <= 0) colour_limit <- 1
  size_limit <- evidence |>
    dplyr::count(contrast_id, reactome_id) |>
    dplyr::summarise(n = max(n)) |>
    dplyr::pull(n)

  manifest_rows <- membership_rows <- list()
  run_index <- membership_index <- 0L
  for (contrast_index in seq_len(nrow(contrasts))) {
    contrast <- contrasts[contrast_index, ]
    contrast_input <- filtered |> dplyr::filter(contrast_id == contrast$contrast_id)
    contrast_evidence <- evidence |> dplyr::filter(contrast_id == contrast$contrast_id)
    scopes <- make_scopes(contrast_input$cell_type)
    for (scope in scopes) {
      directory <- file.path(opt$output, contrast$contrast_id, scope$directory)
      dir.create(directory, recursive = TRUE, showWarnings = FALSE)
      matched <- scope$cells[match(normalise_cell(scope$requested), normalise_cell(scope$cells))]
      membership_index <- membership_index + 1L
      membership_rows[[membership_index]] <- data.frame(
        contrast_id = contrast$contrast_id, scope_id = scope$id, scope_label = scope$label,
        scope_type = scope$type, requested_cell_type = scope$requested,
        matched_cell_type = matched, present = !is.na(matched)
      )
      scoped <- contrast_evidence |> dplyr::filter(cell_type %in% scope$cells)
      nodes <- if (nrow(scoped)) {
        summarise_nodes(scoped, ontology, scope)
      } else data.frame(reactome_id = character())
      n_significant_nodes <- nrow(nodes)
      scope_note <- if (length(scope$missing)) {
        paste0("Partial scope: ", length(scope$cells), "/", length(scope$requested),
          " requested populations present in the significant input.")
      } else paste0(length(scope$cells), " cell populations in this scope.")
      scope_metadata <- data.frame(
        contrast_id = contrast$contrast_id, comparison = contrast$comparison,
        case = contrast$case, reference = contrast$reference,
        scope_id = scope$id, scope_label = scope$label, scope_type = scope$type,
        n_requested_cell_types = length(scope$requested),
        n_matched_cell_types_in_significant_input = length(scope$cells),
        requested_cell_types = collapse_values(scope$requested),
        matched_cell_types = collapse_values(scope$cells),
        missing_cell_types_in_significant_input = collapse_values(scope$missing),
        scope_note = scope_note, symmetric_NES_colour_limit = colour_limit,
        node_count_size_limit = size_limit, mean_NES_population = "significant cells only",
        depth_definition = "maximum root path in full parent/subpathway hierarchy",
        ontology_source = as.character(ontology$source)
      )
      readr::write_csv(scope_metadata, file.path(directory, "scope_metadata.csv"))
      status <- "skipped_no_mapped_significant_terms"
      error_message <- ""
      n_edges <- n_dashed <- 0L
      if (n_significant_nodes) {
        message("Building DAG: ", contrast$comparison, " / ", scope$label,
          " (", n_significant_nodes, " significant terms)")
        if (opt$include_root && !"REACTOME_ROOT" %in% nodes$reactome_id) {
          root <- data.frame(reactome_id = "REACTOME_ROOT", reactome_term = "All Reactome pathways",
            mean_NES = NA_real_, mean_abs_NES = NA_real_, min_NES = NA_real_,
            max_NES = NA_real_, sd_NES = NA_real_, n_significant_cell_types = 0L,
            n_positive_cells = 0L, n_negative_cells = 0L, best_q = NA_real_,
            significant_cell_types = "", significant_pathway_labels = "",
            positive_fraction = NA_real_, mixed_direction = FALSE,
            is_significant = FALSE, scope_id = scope$id, scope_label = scope$label,
            scope_type = scope$type, n_scope_cell_types_in_significant_input = length(scope$cells),
            depth_min = 0L, depth_max = 0L)
          nodes <- dplyr::bind_rows(root, nodes)
        }
        tryCatch({
          edges <- compress_reactome_edges(nodes$reactome_id, ontology$edges,
            nodes$reactome_id[nodes$is_significant])
          nodes <- layout_nodes(nodes, edges)
          edges <- edges |>
            dplyr::mutate(scope_id = scope$id, scope_label = scope$label,
              contrast_id = contrast$contrast_id) |>
            dplyr::left_join(nodes |> dplyr::select(parent_reactome_id = reactome_id,
              parent_term = reactome_term), by = "parent_reactome_id") |>
            dplyr::left_join(nodes |> dplyr::select(child_reactome_id = reactome_id,
              child_term = reactome_term), by = "child_reactome_id")
          readr::write_csv(nodes, file.path(directory, "nodes.csv"))
          readr::write_csv(edges, file.path(directory, "edges.csv"))
          readr::write_csv(scoped, file.path(directory, "significant_cell_evidence.csv"))
          title <- paste(scope$label, contrast$comparison, sep = " | ")
          write_dag_html(nodes, edges, file.path(directory, "DAG.html"), title,
            colour_limit, size_limit, scope_note)
          write_static_dag(nodes, edges, directory, title, colour_limit, size_limit,
            opt$label_top, scope_note)
          expected_files <- file.path(directory,
            c("DAG.html", "DAG.pdf", "DAG.svg", "nodes.csv", "edges.csv"))
          if (any(!file.exists(expected_files)) ||
              any(file.info(expected_files)$size == 0L)) {
            stop("One or more graph exports are missing or empty.")
          }
          n_edges <- nrow(edges)
          n_dashed <- sum(edges$line_style == "dashed")
          status <- if (length(scope$missing)) "completed_partial" else "completed"
        }, error = function(e) {
          status <<- "failed"
          error_message <<- conditionMessage(e)
          warning("Scope failed: ", scope$label, ": ", error_message)
        })
      }
      run_index <- run_index + 1L
      manifest_rows[[run_index]] <- dplyr::bind_cols(scope_metadata, data.frame(
        status = status, error_message = error_message,
        n_significant_Reactome_terms = n_significant_nodes,
        n_displayed_nodes_including_root = if (n_significant_nodes) nrow(nodes) else 0L,
        n_edges = n_edges, n_dashed_edges = n_dashed,
        output_directory = directory,
        html = if (grepl("^completed", status)) file.path(directory, "DAG.html") else ""
      ))
      readr::write_csv(dplyr::bind_rows(manifest_rows),
        file.path(opt$output, "scope_run_manifest.csv"))
    }
  }
  readr::write_csv(dplyr::bind_rows(membership_rows),
    file.path(opt$output, "scope_cell_type_manifest.csv"))
  write_graph_index(dplyr::bind_rows(manifest_rows), opt$output)
  readr::write_csv(contrasts, file.path(opt$output, "comparisons.csv"))
  parameters <- data.frame(parameter = names(opt),
    value = vapply(opt, as.character, character(1)))
  parameters <- dplyr::bind_rows(parameters, data.frame(
    parameter = c("ontology_source", "colour_limit", "size_limit",
      "original_input_rows", "significant_input_rows", "unique_mapped_Reactome_cell_rows",
      "unmapped_or_ambiguous_pathway_labels"),
    value = as.character(c(ontology$source, colour_limit, size_limit, nrow(input),
      nrow(filtered), nrow(evidence), n_omitted_terms))))
  readr::write_csv(parameters, file.path(opt$output, "run_parameters.csv"))
  input_files <- unique(c(opt$input, stats::na.omit(c(opt$mapping, opt$hierarchy_rds))))
  readr::write_csv(data.frame(file = input_files, md5 = unname(tools::md5sum(input_files))),
    file.path(opt$output, "input_checksums.csv"))
  writeLines(capture.output(utils::sessionInfo()), file.path(opt$output, "sessionInfo.txt"))
  failed <- dplyr::bind_rows(manifest_rows) |> dplyr::filter(status == "failed")
  if (nrow(failed)) stop(nrow(failed), " DAG scope(s) failed; inspect scope_run_manifest.csv.")
  message("Significant Reactome DAGs completed: ", opt$output,
    "\nOpen index.html in a browser. Full results: scope_run_manifest.csv")
}

opt <- parse_options(commandArgs(trailingOnly = TRUE))
if (opt$help) {
  print_help()
} else {
  required_packages <- if (opt$self_test) "igraph" else c(
    "readr", "dplyr", "igraph", "ggplot2", "jsonlite", "svglite"
  )
  missing <- required_packages[!vapply(required_packages, requireNamespace,
    logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  if (opt$self_test) run_self_tests() else main(opt)
}
