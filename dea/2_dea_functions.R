
# --------------------- FUNCTION: Create Gene Markers table ------------------
# Identify marker genes (DEGs) between  groups of cells (ident.1 vs. ident.2)
# and summarize their expression levels
GeneMarkersTable <- function(
    object, 
    out.name, 
    ident.1, 
    ident.2, 
    test.use = "wilcox", 
    min.pct = 0.1, 
    logFC = 0.25, 
    latent.vars = NULL
  ) {
  
  # Get mean expression values for cells with identity ident.1 and ident.2
  cells.1 <- WhichCells(object = object, idents = ident.1)
  exp1 <- apply(
    GetAssayData(object = object, assay = "RNA")[, cells.1, drop = F], 
    1, 
    function(x) mean(x = expm1(x = x)) # reverse log-normalization 
  )
  cells.2 <- WhichCells(object = object, idents = ident.2)
  exp2 <- apply(
    GetAssayData(object = object, assay = "RNA")[, cells.2, drop = F],
    1, 
    function(x) mean(x = expm1(x = x))
  )
  
  # Determine DE between the two groups
  cluster_markers <- FindMarkers(
    object = object, 
    assay = "RNA", 
    ident.1 = ident.1, 
    ident.2 = ident.2, 
    min.pct = min.pct, 
    logfc.threshold = logFC, 
    test.use = test.use, 
    latent.vars = latent.vars
  )
  
  # Combine into a result table and write out
  df <- data.frame(
    geneID = rownames(cluster_markers),
    pct.1 = cluster_markers$pct.1, pct.2 = cluster_markers$pct.2,
    avg_exp.1 = exp1[rownames(cluster_markers)], 
    avg_exp.2 = exp2[rownames(cluster_markers)],
    avg_log2FC = cluster_markers$avg_log2FC,
    p_val = cluster_markers$p_val,
    p_val_adj = cluster_markers$p_val_adj
  )
  write.table(df, file=out.name, sep="\t", quote=FALSE, row.names = FALSE)
  
  return()
}

# --------------------- FUNCTION: Perform Pairwise DEA ------------------------
# perform pairwise DEA between all unique pairs of clusters 
ClusterGeneMarkersByPairs <- function(
    object, 
    out.dir, 
    id = "", 
    clusters = NULL, 
    test.use = "wilcox", 
    min.pct = 0.1, 
    logFC = 0.25, 
    latent.vars = NULL
) {
  # Set cluster identity class if id is specified
  if (id != "") { Idents(object) <- id }
  # If clusters is not provided, use all cluster levels in the active identity
  if (length(clusters) == 0) { clusters <- levels(object@active.ident) }
  
  # iterate over all unique combinations of clusters with no repeats or 
  # self-comparisons 
  for (i in 1:(length(clusters)-1)) {    
    for (j in (i+1):length(clusters)) {
      deg.table <- paste(
        out.dir, "/DEG_", 
        test.use, "_cl", 
        clusters[i], "-", 
        clusters[j], 
        ".tsv", 
        sep=""
      )
      GeneMarkersTable(
        object, 
        out.name = deg.table, 
        ident.1 = clusters[i], 
        ident.2 = clusters[j], 
        test.use = test.use, 
        min.pct = min.pct, 
        logFC = logFC, 
        latent.vars = latent.vars
      )            
    }
  }
  return()
}

# --------------------- FUNCTION: Perform One vs All DEA ----------------------
ClusterGeneMarkersVsAll <- function(
  object, 
  out.dir, 
  id = "", 
  clusters = NULL, 
  test.use = "wilcox", 
  min.pct = 0.1, 
  logFC = 0.25, 
  latent.vars = NULL
) {
  
  if (id != "") Idents(object) <- id
  if (length(clusters) == 0) clusters <- levels(object@active.ident)
  
  for (i in 1:length(clusters)) {
    deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-all.tsv", sep="")
    GeneMarkersTable(
      object, 
      out.name = deg.table, 
      ident.1 = clusters[i], 
      ident.2 = setdiff(clusters, c(clusters[i])), 
      test.use = test.use, 
      min.pct = min.pct, 
      logFC = logFC, 
      latent.vars = latent.vars
    )    
  }
  return()
}

# ----- FUNCTION: Filter most significant and biologically relevant DEGs ------
FilterClusterGeneMarkers <- function(
    input.table, 
    output.table, 
    logFC.filt = 1, 
    adjpval.filt = 0.1, 
    min.pct = 0.1, 
    num = 100, 
    sort = FALSE, 
    genes.use = NULL, 
    genes.filter = NULL
) {
  m <- read.table(input.table, sep="\t", header=TRUE)
  
  #  keeps genes that meet all criteria:
  #     Absolute log2 fold-change > 1 (default)
  #     Adjusted p-value < 0.1 (default)
  #     Gene expressed in > 10% of cells in either cluster 1 or cluster 2
  fold_change_filter <- abs(m$avg_log2FC) > logFC.filt
  pvalue_filter <- m$p_val_adj < adjpval.filt
  expression_filter <- (m$pct.1 > min.pct | m$pct.2 > min.pct)
  
  v <- m[which(fold_change_filter & pvalue_filter & expression_filter),]
  
  # exclude undefined p-values and HIV genes
  v <- v[!is.na(v[,1]),]
  
  # exclude gene.filter genes / keep only specified genes
  if (!is.null(genes.filter)) v <- v[!(v$geneID %in% genes.filter),]
  if (!is.null(genes.use)) v <- v[v$geneID %in% genes.use,]
  
  # print the top num DEG, according to the value of ad. p-value
  num <- min(num, nrow(v))
  v <- v[1:num,]
  if (sort) {
    v <- v[order(v$avg_log2FC, decreasing=TRUE),]
  }
  write.table(v, file=output.table,  sep="\t", quote=FALSE, row.names = FALSE)
  
  return()
}

# --------------------- FUNCTION: Generate heatmap ----------------------
HeatmapPlot <- function(
  object, 
  filt.table, 
  plot.name, 
  cells = NULL, 
  width = 5, 
  height = 7
) {
  table <- read.table(filt.table, sep="\t", header = TRUE)
  genes <- as.character(table$geneID)
  if (length(genes) == 0) next
  
  pdf(plot.name, width = width, height = height)    
  print(
    DoHeatmap(
      object = object, 
      cells = cells, 
      features = genes, 
      label = TRUE, 
      size = 2
    )
  )
  dev.off()
  
  return()
}


# --------------------- FUNCTION: Generate volcano plot ----------------------
VolcanoPlotFilter <- function(
  all.table, 
  filt.table, 
  plot.name, 
  title = "", 
  subtitle = "", 
  genes.use = NULL, 
  genes.filter = NULL, 
  revert = FALSE
) {
  filt <- read.table(filt.table, sep="\t", header=TRUE)
  
  label.genes <- color.genes <- as.character(filt$geneID) # otherwise they are seen as factors!!
  
  VolcanoPlot(
    all.table = all.table, 
    color.genes = color.genes, 
    label.genes = label.genes, 
    plot.name = plot.name, 
    title = title, 
    subtitle = subtitle, 
    genes.use = genes.use, 
    genes.filter = genes.filter, 
    revert = revert
  )
  
  return()
}

# --------------------- FUNCTION: Plot pairwise DEGs  ----------------------
FilterClusterGeneMarkersPairs <- function(
    object, 
    out.dir, 
    id = "", 
    clusters = NULL, 
    test.use = "wilcox", 
    logFC.filt = 1, 
    adjpval.filt = 0.1, 
    min.pct = 0.1, 
    num = 100, 
    genes.use = NULL, 
    heatmap = FALSE, 
    volcano = FALSE
) {
  
  if (id != "") Idents(object) <- id
  if (length(clusters) == 0) clusters <- levels(Idents(object))
  
  dir.create(paste(out.dir, "heatmap", sep="/"), showWarnings = FALSE)
  dir.create(paste(out.dir, "volcano", sep="/"), showWarnings = FALSE)
  
  for (i in 1:(length(clusters)-1)) {
    for (j in (i+1):length(clusters)) {
      deg.table <- paste0(
        out.dir, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], ".tsv"
      )
      filt.table <- paste0(
        out.dir, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], "-filtered.tsv"
      )
      if (file.exists(deg.table)) {
        FilterClusterGeneMarkers(
          deg.table, 
          filt.table, 
          logFC.filt = logFC.filt, 
          adjpval.filt = adjpval.filt, 
          min.pct = min.pct, 
          num = num, 
          genes.use = genes.use
        )
        if (heatmap) {
          cells <- colnames(
            GetAssayData(object))[Idents(object) %in% sort(
              c(clusters[i], clusters[j])
            )]
          plot.name <- paste0(
            out.dir, "/heatmap/heatmap_DEG_", test.use, "_cl", clusters[i], "-", clusters[j], ".png"
          )
          HeatmapPlot(object, filt.table, plot.name, cells = cells)
        }
        if (volcano) {
          volcano.plot <- paste0(
            out.dir, "/volcano/volcano_DEG_", test.use, "_cl", clusters[i], "-", clusters[j], ".png"
          )
          title <- paste0(
            "cl", clusters[i], "-", clusters[j]
          )
          VolcanoPlotFilter(
            deg.table, 
            filt.table, 
            plot.name = volcano.plot, 
            title = title, 
            genes.use = genes.use
          )
          
          volcano_rev.plot <- paste0(
            out.dir, "/volcano/volcano_DEG_", test.use, "_cl", clusters[j], "-", clusters[i], ".png"
          )
          title_rev <- paste0("cl", clusters[j], "-", clusters[i])
          VolcanoPlotFilter(
            deg.table, 
            filt.table, 
            plot.name = volcano_rev.plot, 
            title = title_rev, 
            revert = TRUE, 
            genes.use = genes.use
          )
        }
      }
    }
  }
}

# --------------------- FUNCTION: Plot DEGs One vs All  ----------------------
FilterClusterGeneMarkersAll <- function(
  object, 
  out.dir, 
  id = "", 
  clusters = NULL, 
  test.use = "wilcox", 
  logFC.filt = 1, 
  adjpval.filt = 0.1, 
  min.pct = 0.1, 
  num = 100, 
  genes.use = NULL, 
  heatmap = FALSE, 
  volcano = FALSE
) {
  
  dir.create(paste(out.dir, "heatmap", sep="/"))
  dir.create(paste(out.dir, "volcano", sep="/"))
  
  if (id != "") Idents(object) <- id
  if (length(clusters) == 0) clusters <- levels(Idents(object))
  
  cells <- colnames(GetAssayData(object))[Idents(object) %in% clusters]
  
  for (i in 1:length(clusters)) {
    deg.table <- paste0(
      out.dir, "/DEG_", test.use, "_cl", clusters[i], "-all.tsv"
    )
    filt.table <- paste0(
      out.dir, "/DEG_", test.use, "_cl", clusters[i], "-all-filtered.tsv"
    )
    if (file.exists(deg.table)) {
      FilterClusterGeneMarkers(
        deg.table, 
        filt.table, 
        logFC.filt = logFC.filt, 
        adjpval.filt = adjpval.filt, 
        min.pct = min.pct, 
        num = num, 
        genes.use = genes.use
      )
      if (heatmap) {
        plot.name <- paste0(
          out.dir, "/heatmap/heatmap_DEG_", test.use, "_cl", clusters[i], "-all.pdf"
        )
        HeatmapPlot(
          object, 
          filt.table, 
          plot.name, 
          cells = cells, 
          width = 9, 
          height = 7
        )
      }
      if (volcano) {
        volcano.plot <- paste0(
          out.dir, "/volcano/volcano_DEG_", test.use, "_cl", clusters[i], "-all.pdf"
        )
        title <- paste0(
          "cl", clusters[i], "-all"
        )
        VolcanoPlotFilter(
          deg.table, 
          filt.table, 
          plot.name = volcano.plot, 
          title = title, 
          genes.use = genes.use
        )
      }
    }
  }
  
  return()
}