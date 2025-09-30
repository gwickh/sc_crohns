library(CellChat)
library(Seurat)
library(Matrix)
library(patchwork)
library(igraph)
library(future)
library(NMF)
library(ggalluvial)
library(ComplexHeatmap)

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


# dotplot -----------------------------------------------------------------
selectK(cellchat.normal, pattern = "outgoing") # 6 patterns
selectK(cellchat.normal, pattern = "incoming") # 7 patterns
selectK(cellchat.crohns, pattern = "outgoing") # 5 patterns
selectK(cellchat.crohns, pattern = "incoming") # 7 patterns

plot_pathways_dotplot <- function(filename, cellchat_obj, k_in, k_out) {
  pdf(file.path(path, paste0(filename, "_higher_order.pdf")), width = 10, height = 12)
  
  cellchat.out <- identifyCommunicationPatterns(
    cellchat_obj, 
    pattern = "outgoing", 
    k = k_out,
    width = 4,
    height = 15
  )
  gg1 <- netAnalysis_dot(
    cellchat.out, 
    pattern = "outgoing"
  )
  
  cellchat.in <- identifyCommunicationPatterns(
    cellchat_obj, 
    pattern = "incoming", 
    k = k_in,
    width = 4,
    height = 15
  )
  gg2 <- netAnalysis_dot(
    cellchat.in, 
    pattern = "incoming"
  )
  (gg1 | gg2) + plot_layout(widths = c(1,1))
  dev.off()
}

plot_pathways_dotplot("normal_interactions", cellchat.normal, 6, 5)
plot_pathways_dotplot("crohns_interactions", cellchat.crohns, 6, 5)

# dotplots
cellchat.normal_out <- identifyCommunicationPatterns(
  cellchat.normal, 
  pattern = "outgoing", 
  k = 5
)
cellchat.normal_in <- identifyCommunicationPatterns(
  cellchat.normal, 
  pattern = "incoming", 
  k = 6
)

pdf(file.path(path, "normal_pathways_dotplot.pdf"), width = 24, height = 6)
gg1 <- netAnalysis_dot(cellchat.normal_out, pattern = "outgoing")
gg2 <- netAnalysis_dot(cellchat.normal_in, pattern = "incoming")
(gg1 | gg2) + plot_layout(widths = c(1,1))
dev.off()

cellchat.crohns_out <- identifyCommunicationPatterns(
  cellchat.crohns, 
  pattern = "outgoing", 
  k = 5
)
cellchat.crohns_in <- identifyCommunicationPatterns(
  cellchat.crohns, 
  pattern = "incoming", 
  k = 6
)

pdf(file.path(path, "crohns_pathways_dotplot.pdf"), width = 24, height = 6)
gg1 <- netAnalysis_dot(cellchat.crohns_out, pattern = "outgoing")
gg2 <- netAnalysis_dot(cellchat.crohns_in, pattern = "incoming")
(gg1 | gg2) + plot_layout(widths = c(1,1))
dev.off()

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
