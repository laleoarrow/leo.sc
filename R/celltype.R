#' Calculate ROIE
#'
#' Calculate the Ro/e value from the given crosstab
#'
#' @param crosstab the contingency table of given distribution
#'
#' @return matrix of Ro/e values
#' @export
#' @note This function is from https://github.com/yuyang3/pan-B/blob/main/Figure4.R (line 545)
ROIE <- function(crosstab){
  # helper function
  divMatrix <- function(m1, m2){
    ## Divide each element in turn in two same dimension matrixes
    ## Returns:
    ## a matrix with the same dimension, row names and column names as m1.
    ## result[i,j] = m1[i,j] / m2[i,j]
    dim_m1 <- dim(m1) # m1: the first matrix
    dim_m2 <- dim(m2) # m2: the second matrix
    if( sum(dim_m1 == dim_m2) == 2 ){
      div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
      row.names(div.result) <- row.names(m1)
      colnames(div.result) <- colnames(m1)
      for(i in 1:dim_m1[1]){
        for(j in 1:dim_m1[2]){
          div.result[i,j] <- m1[i,j] / m2[i,j]
        }
      }
      return(div.result)
    }
    else{
      warning("The dimensions of m1 and m2 are different")
    }
  }
  # main function
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}

#' Compute Ro/e and draw heatmap
#'
#' Build a Ro/e (observed / expected) matrix from single-cell metadata and,
#' optionally, save a heatmap.
#' @param srt              seurat object
#' @param filter_col       column used to subset (optional)
#' @param filter_criteria  values kept in `filter_col`
#' @param anno_col         row group column; default `"cell_anno"`
#' @param group_col        column group column; default `"Stage1"`
#' @param group_col_order  column order in heatmap; auto if `NULL` (based on the level of `group_col` if it is a factor)
#' @param plot             save heatmap if `TRUE`
#' @param col_fun          colour scale: Defaut orange style or `redblue`
#' @param plot_path        PDF path; default `"./ROIE.pdf"`
#' @param width,height     device size; auto when `NULL`
#' @param fontsize         text size inside heatmap
#' @param heatmap_anno     `"num"`, `"+++"`, or `"none"`
#' @param sym_break        numeric breakpoints; **only two presets supported** --> `c(-Inf,.1,1,2,3,Inf)` or `c(-Inf,0,.2,.8,1,Inf)`
#' @param ...              passed to `ComplexHeatmap::Heatmap()`
#'
#' @return Ro/e matrix
#'
#' @examples
#' \dontrun{
#' roe <- leo.ROIE(all, filter_col = 'lineage', filter_criteria = c('DC'), fontsize = 4,
#'                 anno_col = 'cell_anno', group_col = 'Stage1', col_fun = "redblue",
#'                 plot = T, plot_path = "./figure/sc/celltype_analysis/dc.roe.pdf",
#'                 heatmap_anno = "+++", border = NA, rect_gp = gpar(col = NA))
#' head(roe)
#' }
#' @export
leo.ROIE <- function(srt, filter_col = NULL, filter_criteria = NULL,
                     anno_col = "cell_anno", group_col = "Stage1", group_col_order = NULL,
                     plot = F, plot_path = NULL, col_fun = NULL,
                     width = NULL, height = NULL, fontsize = 6,
                     heatmap_anno = c("num", "+++", "none"),
                     sym_break = c(-Inf, .1, 1, 2, 3, Inf), ...) {
  meta <- srt@meta.data
  if (!is.null(filter_col) && !is.null(filter_criteria)) {
    if (class(filter_criteria) != "vector") filter_criteria <- as.vector(filter_criteria)
    meta <- meta[meta[[filter_col]] %in% filter_criteria, ]
    leo.basic::leo_log("Filter {.emph {filter_col}} for {.emph {filter_criteria}} --> {nrow(meta)} cell{?s} kept")
  }
  if (nrow(meta) == 0) return("Stopped for no cells left after filtering. Please check your filter criteria.")

  if (is.null(group_col_order)) {
    if (is.factor(meta[[group_col]])) {
      group_col_order <- levels(meta[[group_col]])
    } else {
      group_col_order <- sort(unique(meta[[group_col]]))
    }
  }
  meta[[anno_col]] = as.character(meta[[anno_col]])
  meta[[group_col]] = as.character(meta[[group_col]])

  summary <- table(meta[[anno_col]], meta[[group_col]])
  roe <- as.data.frame(ROIE(summary))
  roe <- roe[, group_col_order]

  if (plot) {
    require(ComplexHeatmap);require(circlize)
    if (is.null(col_fun)) col_fun <- circlize::colorRamp2(c(0, 1, 1.5, 2), c("#fffde7", "#ffe0b2","#ff9800", "#e65100"))
    if (identical(col_fun, "redblue")) col_fun <- circlize::colorRamp2(c(0, 1, 1.5, 2), c("#83CEF3", "white","#B45C5E", "#671A19")) # https://doi.org/10.1038/s41591-023-02371-y

    plot_path <- ifelse(is.null(plot_path), "./ROIE.pdf", plot_path)
    leo.basic::leo_log("Save heatmap ➜ {.path {plot_path}}")
    if (!dir.exists(dirname(plot_path))) dir.create(dirname(plot_path), recursive = TRUE)

    if (is.null(width)) width <- 1 + ncol(roe) * .2
    if (is.null(height)) height <- .75 + nrow(roe) * .15

    plus_sym <- function(v) {
      brks <- sym_break
      lbls <- c("-", "+/-", "+", "++", "+++")
      as.character(cut(v, breaks = brks, labels = lbls, right = TRUE, include.lowest = TRUE))
    }
    cell_fun <- switch(heatmap_anno,
                       num  = function(j, i, x, y, w, h, fill)
                         grid::grid.text(sprintf("%.1f", roe[i, j]), x, y,
                                         gp = grid::gpar(fontsize = fontsize)),
                       "+++" = function(j, i, x, y, w, h, fill)
                         grid::grid.text(plus_sym(roe[i, j]), x, y,
                                         gp = grid::gpar(fontsize = fontsize)),
                       none  = NULL)
    if (heatmap_anno == "+++") {
      if (identical(sym_break, c(-Inf, .1, 1, 2, 3, Inf))) {
        leo.basic::leo_log("Annotate the roe with: +++ (>3); ++ (2–3); + (1–2); +/- (0.1–1); - (<0.1)")
        leo.basic::leo_log("Refer to: https://www.sciencedirect.com/science/article/pii/S2352396424004341?#fig4")
      }
      if (identical(sym_break, c(-Inf, 0, 0.2, 0.8, 1, Inf))) {
        leo.basic::leo_log("Annotate the roe with: +++, Ro/e > 1; ++, 0.8 < Ro/e ≤ 1; +, 0.2 ≤ Ro/e ≤ 0.8; ±, 0 < Ro/e < 0.2; −, Ro/e = 0")
        leo.basic::leo_log("Refer to: ")
      }
    }

    pdf(plot_path, width = width, height = height)
    ht <- ComplexHeatmap::Heatmap(
      roe, col = col_fun, cluster_rows = T, cluster_columns = F,
      clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D", # https://github.com/yuyang3/pan-B/blob/main/Figure4.R
      column_names_gp = grid::gpar(fontsize = 6), row_names_gp = grid::gpar(fontsize = 6),
      cell_fun = cell_fun, name = "Ro/e", ...
    )
    draw(ht, heatmap_legend_side = "right")
    dev.off()
    return(list(heatmap = ht, roe = roe))
  } else {
    return(roe)
  }
}

#' Run Augur analysis
#'
#' We added functionality to subset srt obj so that you can compare differet group's results
#'
#' @param srt Seurat object
#' @param subset_col Character. Column in meta.data to subset by (e.g. "Stage1")
#' @param subset_value Character vector. Values to retain in `subset_col`
#' @param label_level Character vector. Factor level order for `subset_col` (default: same as `subset_value`)
#' @param label_col Character. Column in meta.data used as label for Augur
#' @param cell_type_col Character. Column in meta.data representing cell types
#' @param n_threads Integer. Number of threads for parallel Augur computation
#' @param return Character. If `"plot"` (default), returns a list of plot and data; otherwise returns raw Augur object
#'
#' @import Seurat
#' @importFrom Augur calculate_auc
#' @importFrom ggplot2 aes geom_point geom_segment scale_color_manual theme element_text
#'
#' @return If `return = "plot"`, a list with `plot` (ggplot object) and `dat` (Augur result); else, Augur result only
#' @export
leo.augur <- function(srt, subset_col = NULL, subset_value = c("Control", "Inactive"),
                      label_level = NULL, label_col = "Stage1", cell_type_col = "cell_anno",
                      n_threads = 8, return = "plot"){
  if (class(subset_value) != "vector") subset_value <- as.vector(subset_value)
  if (is.null(label_level)) label_level <- unique(srt@meta.data[[label_col]])
  if (!is.null(subset_col)) {
    srt <- subset_srt(srt, subset_col, subset_value)
    leo_log("After subset {subset_col} for **{subset_value}**, {ncol(srt)} cell{?s} remains.")
  }
  if (!is.null(levels)) {
    if (length(unique(srt[[label_col]][[1]])) != length(label_level)) stop("The number of levels does not match.")
    srt@meta.data[[label_col]] <- factor(srt@meta.data[[label_col]], levels = label_level)
  }
  leo_log("Calculating AUC using {label_col} as label_col and {cell_type_col} as cell_type_col with {n_threads} thread{?s}")
  augur <- calculate_auc(srt,
                         label_col = label_col,
                         cell_type_col = cell_type_col,
                         n_threads = n_threads)
  if(return == "plot"){
    p <- plot_lollipop(augur) +
      # geom_segment(aes(xend = cell_type, yend = 0.5), size = 1) +
      geom_point(size = 3, aes(color = cell_type)) +
      scale_color_manual(values = cell_anno_color) +
      theme(
        plot.title      = element_text(face = 'bold', hjust = 0, size = 14, color = 'black'),
        legend.title    = element_text(size = 10, face = 'bold', color = 'black'),
        legend.text     = element_text(size = 9, color = 'black'),
        axis.title      = element_text(size = 10, face = 'bold', color = 'black'),
        axis.text       = element_text(size = 8, color = 'black'),
        axis.text.x     = element_text(size = 8, color = 'black'),
        axis.text.y     = element_text(size = 8, color = 'black'),
        panel.grid.minor= element_blank(),
        legend.position ='right'
      )
    return(list(plot = p, dat = augur))
  } else {
    return(augur)
  }
}


#' Run Milo differential abundance workflow on a Seurat object
#'
#' Lightweight wrapper to build a Milo object, define neighbourhoods, count cells,
#' and test differential abundance with optional batch-aware design.
#'
#' @param all Seurat object
#' @param sample Character. Meta column used as sample identifier (default: "orig.ident")
#' @param group Character. Meta column encoding the biological group/condition (default: "Stage1")
#' @param group_level Character vector. Desired factor levels for `group` (controls ordering)
#' @param cell_type Ignored in current version (placeholder for future stratification)
#' @param milo_mode Character, either `"fast"` or other. Controls neighbourhood refinement and testing behavior
#' @param batch NULL or character. If non-NULL, meta column name used as batch covariate
#'
#' @return A list with components `da_results` (Milo differential abundance result) and `milo_obj` (constructed Milo object)
#' @examples
#' \dontrun{
#' res <- leo.milo(all)
#' milo_obj   <- res$milo_obj
#' da_results <- res$da_results
#' }
#' @importFrom miloR Milo buildGraph makeNhoods plotNhoodSizeHist countCells calcNhoodDistance testNhoods buildNhoodGraph
#' @importFrom dplyr distinct
#' @importFrom SingleCellExperiment reducedDimNames
#' @import Seurat
#' @export
leo.milo <- function(all, sample = "orig.ident", group = "case_ctrl",
                     group_level = c("ctrl", "vkh"), # maybe it is better to use only two group?
                     cell_type = "cell_anno", milo_mode = "fast", batch = NULL) { # set to NULL if no batch correction required
  ec <- leo.basic::leo_log; ec("Tutorial: https://www.bioconductor.org/packages/devel/bioc/vignettes/miloR/inst/doc/milo_gastrulation.html?")
  require(miloR); require(SingleCellExperiment)
  if (class(all$orig.ident) != "character") all$orig.ident <- as.character(all$orig.ident)
  if (!identical(levels(all@meta.data[[group]]), group_level)) all@meta.data[[group]] <- factor(all@meta.data[[group]], levels = group_level)

  ec("Converting to a SingleCellExperiment obj...")
  all2 <- as.SingleCellExperiment(all); ec("Created a SingleCellExperiment object")
  n_d <- min(50, ncol(reducedDim(all2, "PCA")))

  all2 <- miloR::Milo(all2); ec("Created a Milo object")
  all2 <- miloR::buildGraph(all2, k = 30, d = n_d); ec("Constructed KNN graph")

  if (milo_mode == "fast") {
    all2 <- miloR::makeNhoods(all2, prop = 0.05, k = 30, d = n_d, refined = TRUE, refinement_scheme="graph"); ec("(fast mode --> refinement_scheme='graph') Defined representative neighbourhoods on the KNN graph")
  } else {
    all2 <- miloR::makeNhoods(all2, prop = 0.1, k = 30, d = n_d, refined = TRUE); ec("Defined representative neighbourhoods on the KNN graph")
  }
  ec("Ploting neighbourhood size histogram --> As a rule of thumb we want to have an average neighbourhood size over 5 x N_samples = {5*length(unique(all@meta.data[[sample]]))}")
  miloR::plotNhoodSizeHist(all2)

  all2 <- miloR::countCells(all2, sample = sample, meta.data = data.frame(colData(all2))); ec("Cells in neighbourhoods Counted")

  # Defining experimental design
  if (is.null(batch)) {
    exp_design <- data.frame(colData(all2))[,c(sample, group)]; ec("The design is NOT batch-aware")
    colnames(exp_design) <- c("sample", "group")
  } else {
    exp_design <- data.frame(colData(all2))[,c(sample, group, batch)]; ec("The design is batch-aware")
    exp_design[[batch]] <- as.factor(exp_design[[batch]])
    colnames(exp_design) <- c("sample", "group", "batch")
  }
  exp_design$sample <- as.factor(exp_design$sample)
  exp_design <- distinct(exp_design)
  rownames(exp_design) <- exp_design$sample

  # Computing neighbourhood connectivity
  if (mode != "fast") {
    ec("Next, we will be computing neighbourhood connectivity")
    ec("It took ~10 hours to calculate for a ~40w cells dataset. I almost thought it has been stucked.")
    ec("This is a computationally intensive step, it may take a while to run. Take a walk or coffee.")
    all2 <- miloR::calcNhoodDistance(all2, d = n_d); ec("Neighbourhood connectivity computed")
  }

  # Testing
  ec("Now we can do the DA test, explicitly defining our experimental design")
  if (mode == "fast") {
    if (is.null(batch)) {
      da_results <- miloR::testNhoods(all2, design = ~ group,
                                      design.df = exp_design,
                                      fdr.weighting="graph-overlap")
    } else {
      da_results <- miloR::testNhoods(all2, design = ~ batch + group,
                                      design.df = exp_design,
                                      fdr.weighting="graph-overlap")
    }
  } else {
    if (is.null(batch)) {
      da_results <- miloR::testNhoods(all2,
                                      design = ~ group,
                                      design.df = exp_design)
    } else {
      da_results <- miloR::testNhoods(all2,
                                      design = ~ batch + group,
                                      design.df = exp_design)
    }
  }
  # ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
  # ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  #   geom_point() +
  #   geom_hline(yintercept = 1)
  all2 <- miloR::buildNhoodGraph(all2); ec("Now you can visualize the milo results.")
  return(list(da_results = da_results, milo_obj = all2))
}

