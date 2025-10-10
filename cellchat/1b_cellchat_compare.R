library(CellChat)
library(Seurat)
library(Matrix)
library(patchwork)
library(igraph)
library(future)
library(NMF)
library(ggalluvial)
library(ComplexHeatmap)
library(scales)

options(stringsAsFactors = FALSE)

# Vars
path <- "dev/project-area/data/crohns_scrnaseq/cellchat_output"

# load files
cellchat.crohns <- readRDS(
  file.path(path, "cellchat_object_curated_Crohn's Disease.Rds")
)
cellchat.normal <- readRDS(
  file.path(path, "cellchat_object_curated_Normal.Rds")
)

crohns_celltypes = levels(cellchat.crohns@idents)
normal_celltypes = levels(cellchat.normal@idents)

object.list <- list(Crohns = cellchat.crohns, Normal = cellchat.normal)

cellchat_merged <- mergeCellChat(
  object.list, add.names = names(object.list), cell.prefix = TRUE
)

# total interactions barplot ----------------------------------------------
pdf(file.path(path, paste0("total_interactions.pdf")), width = 2, height = 4)
gg1 <- compareInteractions(
  cellchat_merged, 
  show.legend = F, 
  group = c(2, 1))
gg1
dev.off()


# heatmap of interactions -------------------------------------------------
plot_interaction_heatmap <- function(filename, cellchat_obj) {
  pdf(file.path(path, paste0(filename, ".pdf")), width = 10, height = 6)
  
  gg1 <- netVisual_heatmap(
    cellchat_obj,
    comparison = c("Normal", "Crohns")
  )
  gg2 <- netVisual_heatmap(
    cellchat_obj, 
    measure = "weight",
    comparison = c("Normal", "Crohns")
  )
  draw(gg1+gg2)
  decorate_heatmap_body(
    gg2@name, {grid.text(
        "Destination (Receiver)", 
        x = unit(-1, "lines"), 
        y = unit(-10, "lines"), 
        gp = gpar(fontsize = 10))}
  )
  dev.off()
}

plot_interaction_heatmap("crohns_interactions_heatmap", cellchat.crohns)
plot_interaction_heatmap("normal_interactions_heatmap", cellchat.normal)
plot_interaction_heatmap("relative_interactions_heatmap", cellchat_merged)

# interaction strength plots ----------------------------------------------
plot_pathways_heatmap <- function(filename, cellchat_obj, slot_name) {
  cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj,  slot.name = slot_name)
  
  pdf(file.path(path, paste0(filename, ".pdf")), width = 10, height = 12)
  
  gg1 <- netAnalysis_signalingRole_heatmap(
    cellchat_obj, 
    pattern = "outgoing",
    width = 8,
    height = 18,
    slot.name = slot_name
  )
  gg2 <- netAnalysis_signalingRole_heatmap(
    cellchat_obj, 
    pattern = "incoming",
    width = 8,
    height = 18,
    slot.name = slot_name
  )
  
  draw(gg1+gg2, gap = unit(0.5, "cm"))
  decorate_heatmap_body(
    gg2@name, {grid.text(
      "Cell type", 
      x = unit(-3, "lines"), 
      y = unit(-6, "lines"), 
      gp = gpar(fontsize = 10))}
  )
  decorate_heatmap_body(
    gg2@name, {grid.text(
      "Signaling pathway", 
      x = unit(-35, "lines"), 
      y = unit(22.5, "lines"),
      rot = 90,   
      gp = gpar(fontsize = 10))}
  )
  dev.off()
}

plot_pathways_heatmap("crohns_pathways_heatmap", cellchat.crohns, "netP")
plot_pathways_heatmap("normal_pathways_heatmap", cellchat.normal, "netP")

plot_pathways_heatmap("crohns_LR_pairs_heatmap", cellchat.crohns, "net")
plot_pathways_heatmap("normal_LR_pairs_heatmap", cellchat.normal, "net")


# dotplots -----------------------------------------------------------------
selectK(cellchat.normal, pattern = "outgoing") # 5 patterns
selectK(cellchat.normal, pattern = "incoming") # 6 patterns
selectK(cellchat.crohns, pattern = "outgoing") # 8 patterns
selectK(cellchat.crohns, pattern = "incoming") # 6 patterns

plot_pathways_dotplot <- function(filename, cellchat_obj, k_in, k_out) {
  pdf(file.path(path, paste0(filename, "_outgoing_higher_order.pdf")), width = 10, height = 10)
  cellchat.out <- identifyCommunicationPatterns(
    cellchat_obj, 
    pattern = "outgoing", 
    k = k_out,
    width = 4,
    height = 15
  )
  dev.off()
  
  pdf(file.path(path, paste0(filename, "_incoming_higher_order.pdf")), width = 10, height = 10)
  cellchat.in <- identifyCommunicationPatterns(
    cellchat_obj, 
    pattern = "incoming", 
    k = k_in,
    width = 4,
    height = 15
  )
  dev.off()
  
  gg1 <- netAnalysis_dot(
    cellchat.out, 
    pattern = "outgoing"
  )
  gg2 <- netAnalysis_dot(
    cellchat.in, 
    pattern = "incoming"
  )
  pdf(file.path(path, paste0(filename, "_pathways_dotplot.pdf")), width = 24, height = 6)
  print((gg1 + gg2) + plot_layout(widths = c(1,1)))
  dev.off()
}

plot_pathways_dotplot("normal_interactions", cellchat.normal, 6, 5)
plot_pathways_dotplot("crohns_interactions", cellchat.crohns, 6, 8)

# manifold learning
pdf(file.path(path, "UMAP.pdf"), width = 10, height = 6)
cellchat.normal <- computeNetSimilarity(cellchat.normal, type = "functional", thresh = 0.25)
cellchat.normal <- netEmbedding(cellchat.normal, type = "functional")
cellchat.normal <- netClustering(cellchat.normal, type = "functional")
umap1 <- netVisual_embedding(
  cellchat.normal, 
  type = "functional", 
  pathway.remove.show = FALSE,
  label.size = 3.5,
  title = "Normal"
  )

cellchat.crohns <- computeNetSimilarity(cellchat.crohns, type = "functional", thresh = 0.25)
cellchat.crohns <- netEmbedding(cellchat.crohns, type = "functional")
cellchat.crohns <- netClustering(cellchat.crohns, type = "functional")
umap2 <- netVisual_embedding(
  cellchat.crohns, 
  type = "functional", 
  pathway.remove.show = FALSE,
  label.size = 3.5,
  title = "Crohn's Disease")
umap1 + umap2
dev.off()

cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "functional", thresh = 0.25)
cellchat_merged <- netEmbedding(cellchat_merged, type = "functional")
cellchat_merged <- netClustering(cellchat_merged, type = "functional")
umap <- netVisual_embeddingPairwise(
  cellchat_merged, 
  slot.name = "netP",
  type = "functional",
  label.size = 3.5,
  comparison = c(1,2)
)

umap_crohns <- subset(umap@data, group == "Crohns")
umap_normal <- subset(umap@data, group == "Normal")

pdf(file.path(path, "UMAP.pdf"), width = 10, height = 6)
gg1 <- ggplot(umap_crohns, aes(x = x, y = y, color = clusters)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = "Clusters",
    title = "Crohn's Disease"
  )

gg2 <- ggplot(umap_normal, aes(x = x, y = y, color = clusters)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = "Clusters",
    title = "Normal"
  )
gg1 + gg2

dev.off()

# identify significant interactions
get_enriched_interactions <- function(filename, cellchat_obj) {
  cellchat_obj <- rankNetPairwise(cellchat_obj)
  enriched_interactions <- list()
  idx <- 1L
  for (k in seq_along(crohns_celltypes)) {
    for (j in seq_along(crohns_celltypes)) {
      res <- identifyEnrichedInteractions(
        from        = crohns_celltypes[k],
        to          = crohns_celltypes[j],
        object      = cellchat_obj,
        bidirection = FALSE,
        pair.only   = TRUE,
        pairLR.use0 = NULL,
        thresh      = 0.05
      )
      if (!is.null(res) && nrow(res) > 0) {
        if ("interaction_name" %in% names(res)) {
          res <- res %>%
            separate(interaction_name, into = c("ligand","receptor"), sep = "\\|", remove = FALSE)
        }
        
        keep_cols <- intersect(c("ligand","receptor"), names(res))
        tmp <- res[, keep_cols, drop = FALSE]
        tmp$from <- crohns_celltypes[k]
        tmp$to   <- crohns_celltypes[j]
        
        tmp <- tmp[, c("from","to","ligand","receptor")]
        
        enriched_interactions[[idx]] <- tmp
        idx <- idx + 1L
      }
    }
  }
  enriched_interactions <- bind_rows(enriched_interactions)
  write.csv(
    enriched_interactions, 
    file.path(path, paste0(filename, "_enriched_interactions.csv"))
  )
}

get_enriched_interactions("crohns", cellchat.crohns)
get_enriched_interactions("normal", cellchat.normal)


# networks ----------------------------------------------------------------
draw_color_legend <- function(
    x = 0, y = 0.5, height = 0.5, width = 0.1,
    col = colorRampPalette(c("blue", "red"))(100),
    labels = c("Minimum", "Maximum"),
    title = "Comm. Prob.") {
  
  n <- length(col)
  y_seq <- seq(y, y + height, length.out = n + 1)
  
  for (i in 1:n) {
    rect(x, y_seq[i], x + width, y_seq[i + 1], col = col[i], border = NA)
  }
  
  # Labels
  text(x + width + 0.05, y, labels[1], adj = 0, cex = 0.8)
  text(x + width + 0.05, y + height, labels[2], adj = 0, cex = 0.8)
  text(x + width + 0.05, y + height / 2, title, srt = 90, adj = 0.5, font = 2, cex = 0.9)
}

# Network plot generator
generate_network <- function() {
  pathways_crohns <- cellchat_merged@netP$Crohns[1]$pathways
  pathways_normal <- cellchat_merged@netP$Normal[1]$pathways
  pathways <- unique(c(pathways_crohns, pathways_normal))
  
  for (pathway in pathways) {
    comms_crohns <- try(subsetCommunication(cellchat.crohns, slot.name = "netP", signaling = pathway), silent = TRUE)
    comms_normal <- try(subsetCommunication(cellchat.normal, slot.name = "netP", signaling = pathway), silent = TRUE)
    
    # Skip pathway if no interactions in either
    if ((inherits(comms_crohns, "try-error") || nrow(comms_crohns) < 1) &&
        (inherits(comms_normal, "try-error") || nrow(comms_normal) < 1)) {
      next
    }
    
    # Construct edge lists
    edges1 <- if (!inherits(comms_crohns, "try-error") && nrow(comms_crohns) > 0) {
      df <- comms_crohns[, c("source", "target", "prob")]
      colnames(df) <- c("from", "to", "weight")
      df
    } else {
      data.frame(from = character(0), to = character(0), weight = numeric(0))
    }
    
    edges2 <- if (!inherits(comms_normal, "try-error") && nrow(comms_normal) > 0) {
      df <- comms_normal[, c("source", "target", "prob")]
      colnames(df) <- c("from", "to", "weight")
      df
    } else {
      data.frame(from = character(0), to = character(0), weight = numeric(0))
    }
    
    # Common node set
    all_nodes <- union(c(edges1$from, edges1$to), c(edges2$from, edges2$to))
    gr1 <- graph_from_data_frame(edges1, directed = TRUE, vertices = data.frame(name = all_nodes))
    gr2 <- graph_from_data_frame(edges2, directed = TRUE, vertices = data.frame(name = all_nodes))
    
    layout_fixed <- layout_with_fr(gr1)
    
    # Edge widths (scaled)
    E(gr1)$width <- if (ecount(gr1) > 0) rescale(E(gr1)$weight, to = c(1, 6)) else numeric(0)
    E(gr2)$width <- if (ecount(gr2) > 0) rescale(E(gr2)$weight, to = c(1, 6)) else numeric(0)
    
    # Edge colors
    edge_colors <- colorRampPalette(c("blue", "red"))(100)
    if (ecount(gr1) > 0) {
      probs1 <- rescale(E(gr1)$weight, to = c(1, 100))
      E(gr1)$color <- edge_colors[round(probs1)]
    } else {
      E(gr1)$color <- NA
    }
    
    if (ecount(gr2) > 0) {
      probs2 <- rescale(E(gr2)$weight, to = c(1, 100))
      E(gr2)$color <- edge_colors[round(probs2)]
    } else {
      E(gr2)$color <- NA
    }
    
    pdf(file.path(path, paste0("network_", pathway, ".pdf")), width = 14, height = 5)
    layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 0.5))
    par(mar = c(0, 0, 4, 0), oma = c(0, 0, 3, 0))
    par(xpd = TRUE)  # allow labels and edges outside bounds
    
    # Plot Normal
    plot(
      gr2,
      layout = layout_fixed,
      vertex.label.cex = 0.8,
      vertex.label.family = "Helvetica",
      vertex.size = 20,
      edge.arrow.size = 1,
      main = "Normal"
    )
    
    # Plot Crohn's
    plot(
      gr1,
      layout = layout_fixed,
      vertex.label.cex = 0.8,
      vertex.label.family = "Helvetica",
      vertex.size = 20,
      edge.arrow.size = 1,
      main = "Crohn's Disease"
    )
    
    # Legend panel
    plot.new()
    draw_color_legend()
    
    # Main title
    mtext(
      paste0(pathway, " Communication Network"),
      outer = TRUE,
      cex = 1.2,
      line = 1.5,
      adj = 0.35
    )
    
    dev.off()
  }
}

generate_network()