#' Draw alluvial bars with optional custom palette
#'
#' Quickly visualise stacked proportions (e.g. cell-type composition over
#' conditions) as an alluvial plot.
#'
#' @param df           A long-format data frame.
#' @param x_col        Column mapped to the x-axis, default `"Group"`.
#' @param weight_col   Numeric column summed within each x–stratum, default `"Percentage"`.
#' @param stratum_col  Column defining each stacked segment, default `"Cluster"`.
#' @param palette      Optional colour vector; unnamed = applied by order, named = matched by `stratum_col`.
#' @param width        Numeric; width of each stratum (default 0.3).
#' @param x_angle      Rotation angle for x-axis labels. Default `0`. Accepts `0`, `45`, or `90`.
#' @param border_size  Line width for stratum borders. Default `0.5`.
#'
#' @return A `ggplot` object.
#' @importFrom dplyr mutate group_by summarise filter arrange ungroup
#' @importFrom ggplot2 ggplot aes scale_x_discrete scale_y_continuous labs
#'   theme theme_classic element_text guides guide_legend scale_fill_manual
#'   scale_fill_discrete expansion geom_text
#' @importFrom ggalluvial geom_flow geom_stratum
#' @importFrom scales percent_format
#' @importFrom stats setNames
#' @examples
#' library(dplyr)
#' example_data <- tibble(
#'   Group      = rep(c("Ctrl","10","20","30"), each = 7),
#'   Cluster    = rep(c("Il1b+","Cxcl9+","Spp1+","Folr2+",
#'                      "Clps+","Mki67+","Marco+"), times = 4),
#'   Percentage = c(
#'      5,10,15,20,20,20,10,
#'     25,20,15,15,10,10, 5,
#'     30,20,15,10,10,10, 5,
#'     35,25,15,10, 5, 5, 5)
#' )
#'
#' # default palette
#' plot_alluvial(example_data)
#'
#' # custom palette (unnamed vector)
#' my_cols <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C",
#'              "#FB9A99","#E31A1C","#FDBF6F")
#' plot_alluvial(example_data, palette = my_cols)
#'
#' # custom palette (named vector)
#' named_cols <- c(
#'   "Il1b+"  = "#BDD7EE",
#'   "Cxcl9+" = "#6FA8DC",
#'   "Spp1+"  = "#C6E0B4",
#'   "Folr2+" = "#93C47D",
#'   "Clps+"  = "#F4B6C2",
#'   "Mki67+" = "#E06666",
#'   "Marco+" = "#F9CB9C")
#' plot_alluvial(example_data, palette = named_cols)
#' @export
plot_alluvial <- function(df, x_col = "Group", weight_col = "Percentage",
                          stratum_col = "Cluster", width = 0.3, border_size = 0.5,
                          x_angle = 0, palette = NULL) {
  df <- dplyr::mutate(df, !!x_col := factor(.data[[x_col]], levels = unique(.data[[x_col]])))

  # adjust parameters
  fill_scale <- if (!is.null(palette)) {
    if (is.null(names(palette)))
      palette <- stats::setNames(palette, sort(unique(df[[stratum_col]])))
    ggplot2::scale_fill_manual(values = palette, name = NULL)
  } else { ggplot2::scale_fill_discrete(name = NULL) }
  hv <- switch(as.character(x_angle),
               "0"  = list(h = .5, v = 2),
               "90" = list(h = 1,  v = .5),
               "45" = list(h = 1,  v = 1),
               list(h = .5, v = 1))

  p <- ggplot2::ggplot(df,
    ggplot2::aes(x = .data[[x_col]], y = .data[[weight_col]],
                 stratum  = .data[[stratum_col]],
                 alluvium = .data[[stratum_col]],
                 fill     = .data[[stratum_col]])) +
    ggalluvial::geom_flow(stat = "alluvium", lode.guidance = "frontback",
                          color = "grey70", width = width) +
    ggalluvial::geom_stratum(width = width, color = "black", size = border_size) +
    fill_scale +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(.1, .1))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.0, .0)),
                                labels = scales::percent_format(scale = 1)) +
    ggplot2::labs(x = "Group", y = "Percentage (%)") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 14, color = 'black'),
      legend.title    = ggplot2::element_text(face = "bold", size = 12, color = 'black'),
      legend.text     = ggplot2::element_text(size = 10),
      axis.title      = ggplot2::element_text(face = "bold", size = 12, color = 'black'),
      axis.text       = ggplot2::element_text(size = 10, color = 'black'),
      axis.title.x    = ggplot2::element_blank(),
      axis.text.x     = ggplot2::element_text(hjust = hv$h, vjust = hv$v, size = 10, angle = x_angle),
      axis.ticks.x    = ggplot2::element_blank(),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "right") +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(color = NA, size = 0)))

  return(p)
}

#' Alluvial plot from a Seurat object
#'
#' Aggregate cell counts by any two metadata fields and draw an alluvial plot
#' with \code{\link{plot_alluvial}}. Extra arguments are passed straight to
#' \code{plot_alluvial()}.
#'
#' @param obj         A \pkg{Seurat} object.
#' @param group_col   Metadata field mapped to the x-axis (e.g. sample, time).
#' @param cluster_col Metadata field defining strata (e.g. cell type / ident).
#' @param ...         Additional arguments forwarded to \code{plot_alluvial()}.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' library(Seurat)
#' data("pbmc_small")
#' pbmc_small$Group   <- pbmc_small$orig.ident   # mock group
#' pbmc_small$Cluster <- Idents(pbmc_small)      # use idents
#'
#' ## default colours
#' plot_alluvial_sc(pbmc_small)
#'
#' ## custom palette (unnamed)
#' plot_alluvial_sc(
#'   pbmc_small,
#'   palette = c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072",
#'               "#80B1D3","#FDB462","#B3DE69"))
#' @export
plot_alluvial_sc <- function(obj,
                             group_col   = "Group",
                             cluster_col = "Cluster",
                             ...) {

  stopifnot(inherits(obj, "Seurat"))

  meta <- obj[[]] |>
    dplyr::count(.data[[group_col]], .data[[cluster_col]], name = "n") |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::mutate(Percentage = 100 * n / sum(n), .after = n) |>
    dplyr::ungroup()

  plot_alluvial(
    meta,
    x_col       = group_col,
    weight_col  = "Percentage",
    stratum_col = cluster_col,
    ...)
}

#' Plot gene-weighted density
#'
#' Quickly visualize expression gene-weighted density of one or more features
#' based on \pkg{Nebulosa}.
#'
#' @param data     A Seurat object with a "harmony" reduction and UMAP computed on it.
#' @param features Character vector of feature names to plot.
#' @param ncol     Number of columns in the output layout (default 2).
#'
#' @return A patchwork object arranging density plots in a grid.
#'
#' @importFrom Nebulosa plot_density
#' @importFrom ggplot2 coord_fixed theme_void
#' @importFrom patchwork wrap_plots
#' @examples
#' # Assuming 'all' is a Seurat object with harmony UMAP already run
#' library(Seurat)
#' library(Nebulosa)
#' library(patchwork)
#' # Plot density for CD8A and CD8B in two columns
#' p <- plot_gw_density(all, features = c("CD8A", "CD8B"), ncol = 2)
#' print(p)
#' @export
plot_gw_density <- function(data, features, ncol = 2) {
  missing <- setdiff(features, rownames(data))
  if (length(missing)) leo_log("Features not found: ", paste(missing, collapse = ", "), level = "danger")
  plots <- lapply(features, function(gene) {
    Nebulosa::plot_density(
      object    = data,
      reduction = "umap.harmony",
      features  = gene,
      method    = "wkde",
      pal       = "magma",
      size      = 0.2,
      raster    = TRUE
    ) +
      ggplot2::coord_fixed() +
      ggplot2::theme_void()
  })
  patchwork::wrap_plots(plots, ncol = ncol)
}

#' Plot Highlight a cluster on a dimensional reduction
#'
#' @param obj A Seurat object
#' @param cluster_id Cluster identifier to highlight
#' @param reduction Dimensional reduction name (defaults to active reduction)
#' @param group.by Metadata column or NULL to use identities
#' @param highlight.col Color for highlighted cluster
#' @param other.col Color for other cells
#' @param pt.size Point size
#' @param dpi DPI for rasterization if >1e5 cells
#' @return ggplot2 object
#' @importFrom Seurat Embeddings Idents
#' @importFrom SeuratObject DefaultDimReduc
#' @importFrom ggplot2 ggplot aes labs theme_void geom_point
#' @importFrom ggrastr geom_point_rast
#' @importFrom leo.basic leo_log
#' @export
plot_highlight_cluster <- function(obj, cluster_id, reduction = NULL, group.by = NULL,
                                   highlight.col = "red", other.col = "grey80",
                                   pt.size = 1, dpi = 300) {
  if (is.null(reduction)) reduction <- SeuratObject::DefaultDimReduc(obj)
  coords <- Seurat::Embeddings(obj, reduction)[, 1:2]
  grp <- if (is.null(group.by)) Seurat::Idents(obj) else obj@meta.data[[group.by]]

  df <- data.frame(Dim1 = coords[,1], Dim2 = coords[,2], grp = as.character(grp), stringsAsFactors = FALSE)
  df_other <- df[df$grp != cluster_id, ]
  df_high  <- df[df$grp == cluster_id, ]
  use_raster <- nrow(df) > 1e5
  if (use_raster) {
    leo.basic::leo_log("More than 100,000 cells detected, using rasterized points for performance.")
    p <- ggplot2::ggplot() +
      ggrastr::geom_point_rast(data = df_other, aes(Dim1,Dim2),
                               color = other.col, size = pt.size, raster.dpi = dpi) +
      ggrastr::geom_point_rast(data = df_high, aes(Dim1,Dim2),
                               color = highlight.col, size = pt.size, raster.dpi = dpi)
  } else {
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(data = df_other, aes(Dim1,Dim2),
                          color = other.col, size = pt.size) +
      ggplot2::geom_point(data = df_high, aes(Dim1,Dim2),
                          color = highlight.col, size = pt.size)
  }
  p + ggplot2::labs(title = paste0("Cluster ", cluster_id, " highlighted")) +
    ggplot2::theme_void()
}
