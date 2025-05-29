ClusterGeneMarkersByPairs <- function(object, out.dir, id = "", clusters = NULL, test.use = "wilcox", min.pct = 0.1, logFC = 0.25, latent.vars = NULL) {
  
  if (id != "") Idents(object) <- id
  if (length(clusters) == 0) clusters <- levels(object@active.ident)
  
  for (i in 1:(length(clusters)-1)) {    
    for (j in (i+1):length(clusters)) {
      deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], ".tsv", sep="")
      GeneMarkersTable(object, out.name = deg.table, ident.1 = clusters[i], ident.2 = clusters[j], test.use = test.use, min.pct = min.pct, logFC = logFC, latent.vars = latent.vars)            
    }
  }
  
  return()
}

ClusterGeneMarkersVsAll <- function(object, out.dir, id = "", clusters = NULL, test.use = "wilcox", min.pct = 0.1, logFC = 0.25, latent.vars = NULL) {
  
  if (id != "") Idents(object) <- id
  if (length(clusters) == 0) clusters <- levels(object@active.ident)
  
  for (i in 1:length(clusters)) {
    deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-all.tsv", sep="")
    GeneMarkersTable(object, out.name = deg.table, ident.1 = clusters[i], ident.2 = setdiff(clusters, c(clusters[i])), test.use = test.use, min.pct = min.pct, logFC = logFC, latent.vars = latent.vars)    
  }
  
  return()
}

FilterClusterGeneMarkersPairs <- function(object, out.dir, id = "", clusters = NULL, test.use = "wilcox", logFC.filt = 1, adjpval.filt = 0.1, min.pct = 0.1, num = 100, genes.use = NULL, heatmap = FALSE, volcano = FALSE) {
  
  if (id != "") Idents(object) <- id
  if (length(clusters) == 0) clusters <- levels(Idents(object))
  
  dir.create(paste(out.dir, "heatmap", sep="/"), showWarnings = FALSE)
  dir.create(paste(out.dir, "volcano", sep="/"), showWarnings = FALSE)
  
  for (i in 1:(length(clusters)-1)) {
    for (j in (i+1):length(clusters)) {
      deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], ".tsv", sep="")
      filt.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], "-filtered.tsv", sep="")
      if (file.exists(deg.table)) {
        FilterClusterGeneMarkers(deg.table, filt.table, logFC.filt = logFC.filt, adjpval.filt = adjpval.filt, min.pct = min.pct, num = num, genes.use = genes.use)
        if (heatmap) {
          cells <- colnames(GetAssayData(object))[Idents(object) %in% sort(c(clusters[i], clusters[j]))]
          plot.name <- paste(out.dir, "/heatmap/heatmap_DEG_", test.use, "_cl", clusters[i], "-", clusters[j], ".pdf", sep="")
          HeatmapPlot(object, filt.table, plot.name, cells = cells)
        }
        if (volcano) {
          volcano.plot <- paste(out.dir, "/volcano/volcano_DEG_", test.use, "_cl", clusters[i], "-", clusters[j], ".pdf", sep="")
          title <- paste("cl", clusters[i], "-", clusters[j], sep="")
          VolcanoPlotFilter(deg.table, filt.table, plot.name = volcano.plot, title = title, genes.use = genes.use)
          
          volcano_rev.plot <- paste(out.dir, "/volcano/volcano_DEG_", test.use, "_cl", clusters[j], "-", clusters[i], ".pdf", sep="")
          title_rev <- paste("cl", clusters[j], "-", clusters[i], sep="")
          VolcanoPlotFilter(deg.table, filt.table, plot.name = volcano_rev.plot, title = title_rev, revert = TRUE, genes.use = genes.use)
        }
      }
    }
  }
  
  FilterClusterGeneMarkersAll <- function(object, out.dir, id = "", clusters = NULL, test.use = "wilcox", logFC.filt = 1, adjpval.filt = 0.1, min.pct = 0.1, num = 100, genes.use = NULL, heatmap = FALSE, volcano = FALSE) {
    
    dir.create(paste(out.dir, "heatmap", sep="/"))
    dir.create(paste(out.dir, "volcano", sep="/"))
    
    if (id != "") Idents(object) <- id
    if (length(clusters) == 0) clusters <- levels(Idents(object))
    cells <- colnames(GetAssayData(object))[Idents(object) %in% clusters]
    
    for (i in 1:length(clusters)) {
      deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-all.tsv", sep="")
      filt.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-all-filtered.tsv", sep="")
      if (file.exists(deg.table)) {
        FilterClusterGeneMarkers(deg.table, filt.table, logFC.filt = logFC.filt, adjpval.filt = adjpval.filt, min.pct = min.pct, num = num, genes.use = genes.use)
        if (heatmap) {
          plot.name <- paste(out.dir, "/heatmap/heatmap_DEG_", test.use, "_cl", clusters[i], "-all.pdf", sep="")
          HeatmapPlot(object, filt.table, plot.name, cells = cells, width = 9, height = 7)
        }
        if (volcano) {
          volcano.plot <- paste(out.dir, "/volcano/volcano_DEG_", test.use, "_cl", clusters[i], "-all.pdf", sep="")
          title <- paste("cl", clusters[i], "-all", sep="")
          VolcanoPlotFilter(deg.table, filt.table, plot.name = volcano.plot, title = title, genes.use = genes.use)
        }
      }
    }
    
    return()
  }