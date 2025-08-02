#' Calculate ROIE
#'
#' Calculate the Ro/e value from the given crosstab
#'
#' @param crosstab the contingency table of given distribution
#'
#' @return A matrix of Ro/e values
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
#' @param srt              Seurat object
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
#' @param sym_break        numeric breakpoints; **only two presets supported**
#'  supported:** `c(-Inf,.1,1,2,3,Inf)` or `c(-Inf,0,.2,.8,1,Inf)`.
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

    if (is.null(width)) width <- 2 + ncol(roe) * .2
    if (is.null(height)) height <- 1.5 + nrow(roe) * .15

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
