#' Read single-cell data into a Seurat object
#'
#' Supports 10X Matrix, TXT, or RDS formats for raw count import.
#'
#' @param sc_path Character. Path to the data directory or file.
#' @param project_name Character. Project name for Seurat object. Defaults to basename(sc_path) if NULL.
#' @param method Character. Data format to import: "10X_matrix", "txt", or "rds". Default: "10X_matrix".
#' @param gene.column Integer. Column index for gene names in 10X data. Default: 1.
#' @param min.cells Integer. Minimum cells per feature for CreateSeuratObject. Default: 3.
#' @param min.features Integer. Minimum features per cell for CreateSeuratObject. Default: 300.
#' @param verbose Logical. Whether to print `leo.log` messages. Default: TRUE.
#' @param ... Additional arguments passed to CreateSeuratObject().
#'
#' @return A Seurat object of raw counts.
#'
#' @importFrom Seurat Read10X CreateSeuratObject
#' @importFrom data.table fread
#' @importFrom leo.basic leo_log
#' @export
read_sc_data <- function(sc_path,
                         project_name = NULL,
                         method = c("10x", "txt", "rds"),
                         gene.column = 1,
                         min.cells = 3,
                         min.features = 300,
                         verbose = TRUE,
                         ...) {
  project_name <- if (is.null(project_name)) basename(sc_path) else project_name
  method <- match.arg(method)

  if (method == "10x") {
    counts <- Read10X(sc_path, gene.column = gene.column)
  } else if (method == "txt") {
    dt <- data.table::fread(sc_path, data.table = FALSE)
    rownames(dt) <- dt[[1]]; counts <- as.matrix(dt[, -1])
  } else if (method == "rds") {
    seurat_obj <- readRDS(sc_path)
  } else {
    leo.basic::leo_log("Unknown method for reading data", level = "danger", verbose = verbose); stop()
  }

  seurat_obj <- CreateSeuratObject(counts, project = project_name, min.cells = min.cells, min.features = min.features, ...)
  leo.basic::leo_log(paste0("Imported ", ncol(seurat_obj), " cells into Seurat object."), level = "success", verbose = verbose)
  return(seurat_obj)
}


#' Standard seurat normalize_and_scale
#'
#' Performs data normalization, finds variable features, scales data, runs PCA,
#' constructs graph-based neighbors and clusters, and generates UMAP embeddings.
#'
#' @param seurat_obj A Seurat object.
#' @param normalize_method Character. Normalization method:
#'  - "classic" (LogNormalize + FindVariableFeatures + ScaleData)
#'  - "sctransform" (SCTransform).
#'  - "sctransform&regress" (SCTransform + regress out unwanted sources [currently only percent.mt]).
#' @param conserve.memory Logical. Whether to conserve memory in SCTransform. Default: FALSE.
#' @param cluster_resolution Numeric. Resolution for FindClusters. Default: 1.0.
#' @param n_hv_gene Integer. Number of top variable features to label. Default: 10.
#' @param plot_dir Character. Directory to save variable feature plot. Default: NULL (i.e., not plot).
#' @param verbose Logical. Whether to print `leo.log` messages. Default: TRUE.
#'
#' @return A Seurat object after normalization, scaling, clustering, and UMAP embedding.
#'
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData SCTransform DefaultAssay RunPCA FindNeighbors FindClusters RunUMAP VariableFeaturePlot LabelPoints
#' @importFrom ggplot2 ggsave
#' @importFrom leo.basic leo_log
#' @export
seurat_standard_normalize_and_scale <- function(seurat_obj, normalize_method = c("classic","sctransform", "sctransform&regress"),
                                                conserve.memory = F, cluster_resolution = 1.0, plot_dir = NULL, n_hv_gene = 10, verbose = TRUE){
  # normalize and make UMAP
  normalize_method <- match.arg(normalize_method)
  if (normalize_method == "classic"){
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000) # GetAssayData(seurat_obj,slot="data",assay="RNA")
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
    seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj)) # GetAssayData(seurat_obj,slot="scale.data",assay="RNA")
  } else if (normalize_method == "sctransform"){
    seurat_obj <- SCTransform(seurat_obj, variable.features.n = 2000, conserve.memory = conserve.memory)
    DefaultAssay(seurat_obj) <- "SCT" # default for seurat v5 already
  } else if (normalize_method == "sctransform&regress"){
    seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")
    seurat_obj <- SCTransform(seurat_obj, variable.features.n = 2000, conserve.memory = conserve.memory,
                              vars.to.regress = "percent.mt")
    DefaultAssay(seurat_obj) <- "SCT" # default for seurat v5 already
  } else {
    stop("Unknown normalization method")
  }

  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj)) # ElbowPlot(seurat_obj, ndims=20, reduction="pca")
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = cluster_resolution) # BuildClusterTree(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

  if (!is.null(plot_dir)){
    top_genes <- head(VariableFeatures(seurat_obj), n_hv_gene)
    p <- VariableFeaturePlot(seurat_obj) %>%
      LabelPoints(points = top_genes, repel = TRUE, xnudge = 0, ynudge = 0)
    save_path <- file.path(plot_dir, "variable_features.pdf")
    ggsave(save_path, plot = p, width = 8, height = 3, create.dir = T)
    leo.basic::leo_log("vst_plot save to >>> {.path {save_path}}", level = "success", verbose = verbose)
  }

  return(seurat_obj)
}

#' Seurat basic QC
#'
#' Calculates mitochondrial and hemoglobin percentages, filters cells by QC thresholds,
#' and optionally plots QC metrics before and after filtering.
#'
#' @param seurat_obj A Seurat object.
#' @param out_path Character. Directory to save QC plots. Required.
#' @param plot Logical. Whether to generate QC plots. Default: TRUE.
#' @param verbose Logical. Whether to print `leo.log` messages. Default: TRUE.
#'
#' @return A filtered Seurat object.
#'
#' @importFrom Seurat PercentageFeatureSet
#' @importFrom leo.basic leo_log
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom cowplot plot_grid
#' @export
seurat_basic_qc <- function(seurat_obj, out_path = NULL, save_plot = T, verbose = TRUE){
  if (save_plot & is.null(out_path)) {stop("output path for qc plot is required")}
  leo.basic::leo_log("Performing basic QC on the Seurat object", verbose = verbose)

  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  hb_genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
  hb_matches <- match(hb_genes, rownames(seurat_obj@assays$RNA))
  hb_genes_present <- hb_genes[!is.na(hb_matches)]
  seurat_obj[["percent.HB"]] <- PercentageFeatureSet(seurat_obj, features = hb_genes_present)
  if (save_plot) {
    p_before <- plot_qc(seurat_obj, out_path, prefix = "before_", save_plot = F) +
      plot_annotation(
        title = "Before QC",
        theme = theme(
          plot.title = element_text(
            color = "blue", face  = "bold",
            size  = 18, hjust = 0.5))
        )
  }

  n1 <- ncol(seurat_obj)
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 &  nFeature_RNA < 7500 &
                                            nCount_RNA < 10000 & percent.mt < 10 & percent.HB < 1)
  n2 <- ncol(seurat_obj); filtered_cells <- n1 - n2
  leo.basic::leo_log("{filtered_cells} cells were filtered", level = "success", verbose = verbose)

  if (save_plot) {
    p_after <- plot_qc(seurat_obj, out_path, prefix = "after_", save_plot = F) +
      plot_annotation(
        title = "After QC",
        theme = theme(
          plot.title = element_text(
            color = "blue", face  = "bold",
            size  = 18, hjust = 0.5))
      )
    p <- cowplot::plot_grid(p_before, p_after, ncol = 1,
                            align       = "v",
                            rel_heights = c(1, 1)) +
      plot_annotation(
        title = paste0("QC before and after filtering: ", filtered_cells, " cells removed"),
        theme = theme(
          plot.title = element_text(
            color = "red", face  = "bold.italic",
            size  = 18, hjust = 0.5))
      )

    save_path <- file.path(out_path, "qc_plot.pdf")
    ggsave(save_path, plot = p, width = 18, height = 15, create.dir = T)
    leo.basic::leo_log("QC plots save to >>> {.path {save_path}}", level = "success", verbose = verbose)
  }
  return(seurat_obj)
}

#' Plot QC metrics
#'
#' Generates violin and scatter plots for QC features and saves them as PDFs.
#'
#' @param seurat_obj A Seurat object with QC metrics.
#' @param out_path Character. Directory to save plots.
#' @param prefix Character. Filename prefix. Default: "before_".
#' @param verbose Logical. Whether to print `leo.log` messages. Default: TRUE.
#'
#' @return Invisible NULL.
#'
#' @importFrom Seurat VlnPlot FeatureScatter
#' @importFrom ggplot2 ggsave theme
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom leo.basic leo_log
#' @export
plot_qc <- function(seurat_obj, out_path, prefix = "before_", save_plot = T, verbose = TRUE, ...) {
  p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"),
                pt.size = 0, ncol = 4) &
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.HB")
  p2 <- (plot1 + plot2 + plot3) + plot_layout(guides = "collect", nrow = 1) &
    theme(legend.position = "right") &
    guides(color = guide_legend(ncol = 2, byrow = TRUE))
  # p2 <- p2 + plot_annotation(title = "Scatter Plots of QC Metrics")

  p3 <- p1/p2 + plot_layout(ncol = 1, heights = c(1, 1))

  if (save_plot) {
    save <- file.path(out_path, paste0("QC_Plot_", prefix, "combined.pdf"))
    ggsave(filename = save, plot = p3, width = 20, height = 10, create.dir=T)
    leo.basic::leo_log("QC plots for {prefix} save to >>> {.path {save}}", level = "success", verbose = verbose)
  }
  return(p3)
}


#' Remove doublets using DoubletFinder
#'
#' Estimates doublet rate (if not provided), finds optimal pK, and removes
#' heterotypic doublets from a Seurat object.
#'
#' @param seurat_obj Seurat object (after normalization & PCA).
#' @param doublet_rate Numeric. Expected doublet rate (0–1). If NULL, looked up via \code{doublet_rate_dictionary}.
#' @param PCs Integer vector. PCs to use. Default: 1:20.
#' @param nPCs_for_pK Integer vector. PCs to use for pK estimation. Default: 1:10.
#' @param pN Numeric. Proportion of artificial doublets. Default: 0.25.
#' @param sct Logical. Whether data is SCTransformed. Default: FALSE.
#' @param out_path Character. Directory to save UMAP plots. Default: "."
#' @param slim Logical. Whether to slim the Seurat object to only leave RNA count. Default: TRUE.
#' @param verbose Logical. Whether to print `leo.log` messages. Default: TRUE.
#'
#' @return Seurat object with doublets removed.
#'
#' @note
#' - PCs: a vector of statistically significant PCs to use
#' - pN: the number of artificially generated doublets (default = 0.25); robust, normally do not need fine-tune
#' - pK: PC neighborhood size used to compute network; need fine-tune
#' - nExp: threshold used to make doublet/singlet call; based on empirical ratios
#' @importFrom DoubletFinder paramSweep summarizeSweep find.pK doubletFinder
#' @importFrom Seurat DimPlot DietSeurat
#' @importFrom ggplot2 ggsave
#' @export
doublet_removal <- function(seurat_obj=sc_matrix, out_path,
                            pN = 0.25,
                            PCs = 1:20,
                            nPCs_for_pK = 1:10,
                            doublet_rate = NULL,
                            sct = FALSE,
                            slim = TRUE,
                            verbose = TRUE){
  # function for doubletfinder based doublet removal
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    leo.basic::leo_log("Please install DoubletFinder", level = "danger")
    leo.basic::leo_log("remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')")
    stop("DoubletFinder not installed")
  }

  # look up for estimated doublet rate if not provided
  n_cells <- ncol(seurat_obj)
  if (is.null(doublet_rate)) {
    doublet_rate <- doublet_rate_dictionary(n_cells)
    leo.basic::leo_log(sprintf("Doublet rate not provided, using dictionary: %0.3f", doublet_rate), level = "warning", verbose = verbose)
  }

  # estimate doublets and homotypic proportion
  nExp    <- round(doublet_rate * n_cells) # default ~7.5% doublet rate
  hom.prop <- modelHomotypic(seurat_obj@meta.data$seurat_clusters)
  nExp.adj <- round(nExp * (1 - hom.prop))
  leo.basic::leo_log(sprintf("Expected doublets: %d (adjusted: %d)", nExp, nExp.adj), level = "info", verbose = verbose)

  # determine optimal pK
  sweep     <- paramSweep(seurat_obj, PCs = nPCs_for_pK, sct = sct)
  stats     <- summarizeSweep(sweep, GT = FALSE)
  bcmvn     <- find.pK(stats)
  mpK       <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  leo.basic::leo_log(sprintf("Optimal pK: %s", mpK), level = "info", verbose = verbose)

  # run DoubletFinder
  seurat_obj <- DoubletFinder::doubletFinder(
    seurat_obj,
    PCs   = PCs,
    pN    = pN,
    pK    = mpK,
    nExp  = nExp.adj,
    reuse.pANN = NULL,
    sct   = sct
  )

  # rename & clean metadata
  df.col <- grep("^DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
  seurat_obj$doublet.class <- seurat_obj[[df.col]]; seurat_obj[[df.col]] <- NULL
  pann    <- grep("^pANN", names(seurat_obj@meta.data), value = TRUE)
  seurat_obj$pANN <- seurat_obj[[pann]]; seurat_obj[[pann]] <- NULL

  # visualize & subset
  seurat_obj@meta.data$doublet.class <- factor(seurat_obj$doublet.class, levels = c("Singlet", "Doublet"))
  dc <- DimPlot(seurat_obj, reduction = "umap", group.by = "doublet.class", cols = c("gold", "black")) +
    ggtitle("Doublet distribution") & theme(plot.title = element_text(hjust = 0, face = "bold"))

  seurat_obj <- subset(seurat_obj, subset = doublet.class != "Doublet")
  post_cluster <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
    ggtitle("Cluster after doublet removal") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  p <- cowplot::plot_grid(dc, post_cluster, nrow = 1)
  ggsave(plot = p, filename = file.path(out_path, paste0("doublets_filter.pdf")), width = 10, height = 5)
  leo.basic::leo_log("DoubletFinder Done!", level = "success", verbose = verbose)

  # slim if needed
  if (slim) {
    leo.basic::leo_log("Sliming the seruat object!", level = "info", verbose = verbose)
    seurat_obj <-  DietSeurat(seurat_obj, assays = "RNA", layers = "counts",
                              dimreducs = NULL, graphs = NULL, misc = FALSE)
    gc()
  }
  return(seurat_obj)
}

#' Lookup expected doublet rate based on cells loaded
#'
#' Estimates the 10x Genomics multiplet rate (as a fraction) by finding the closest
#' “# of Cells Loaded” bracket in the standard guideline table.
#' The ratios can refer to: https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/76
#' We currently did not find any corresponding ratio for other single cell sequancing platforms.
#'
#' @param n_cells Integer. Number of cells loaded onto the 10x Chromium.
#' @return Numeric. Expected multiplet (doublet) rate.
#' @export
doublet_rate_dictionary <- function(n_cells) {
  # Named vector: cells loaded → multiplet rate
  dict <- c(
    `800`  = 0.004,
    `1600` = 0.008,
    `3200` = 0.016,
    `4800` = 0.023,
    `6400` = 0.031,
    `8000` = 0.039,
    `9600` = 0.046,
    `11200`= 0.054,
    `12800`= 0.061,
    `14400`= 0.069,
    `16000`= 0.076
  )
  keys    <- as.numeric(names(dict))
  nearest <- keys[which.min(abs(keys - n_cells))]
  return(unname(dict[as.character(nearest)]))
}

#' get basic summary of a seurat object
#'
#' @param seurat_obj Seurat object
#' @param assay character. assay name, defaults to active assay
#' @param slot character. data slot ("counts", "data" or "scale.data")
#' @return named list with total_genes, total_cells, avg_genes_per_cell, genes_per_cell
#' @importFrom Seurat DefaultAssay GetAssayData
#' @importFrom leo.basic leo_log
#' @export
seurat_basic_info <- function(seurat_obj, assay = NULL, slot = "counts") {
  if (is.null(assay)) assay <- DefaultAssay(seurat_obj)
  mat <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  total_genes       <- nrow(mat)
  total_cells       <- ncol(mat)
  genes_per_cell    <- colSums(mat > 0)
  avg_genes_per_cell <- mean(genes_per_cell)

  leo.basic::leo_log(sprintf("total genes: %d", total_genes))
  leo.basic::leo_log(sprintf("total cells: %d", total_cells))
  leo.basic::leo_log(sprintf("avg genes per cell: %.2f", avg_genes_per_cell))

  list(
    total_genes        = total_genes,
    total_cells        = total_cells,
    avg_genes_per_cell = avg_genes_per_cell
  )
}
