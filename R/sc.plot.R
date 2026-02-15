#' Draw alluvial bars with optional custom palette
#'
#' Quickly visualise stacked proportions (e.g. cell-type composition over
#' conditions) as an alluvial plot.
#'
#' @param df           A long-format data frame.
#' @param x_col        Column mapped to the x-axis, default `"Group"`.
#' @param weight_col   Numeric column summed within each x-stratum, default `"Percentage"`.
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
#' @importFrom rlang sym
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
  df[[x_col]] <- factor(df[[x_col]], levels = unique(df[[x_col]]))

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
    ggplot2::aes(x = !!rlang::sym(x_col), y = !!rlang::sym(weight_col),
                 stratum  = !!rlang::sym(stratum_col),
                 alluvium = !!rlang::sym(stratum_col),
                 fill     = !!rlang::sym(stratum_col))) +
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
#' @param return      return a plot (set "plot") or a list with plot and data (set "both").
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
plot_alluvial_sc <- function(obj, group_col = "Group", cluster_col = "Cluster",
                             return = "plot", ...) {

  stopifnot(inherits(obj, "Seurat"))

  dat <- meta <- obj[[]]
  dat <- dplyr::count(dat, .data[[group_col]], .data[[cluster_col]], name = "n")
  dat <- dplyr::group_by(dat, .data[[group_col]])
  dat <- dplyr::mutate(dat, Percentage = 100 * n / sum(n), .after = n)
  dat <- dplyr::ungroup(dat)

  p <- plot_alluvial(dat, x_col = group_col,
                     weight_col  = "Percentage",
                     stratum_col = cluster_col, ...)

  if (return == "plot") return(p)
  return(list(plot = p, data = dat))
}

#' Plot gene-weighted density
#'
#' Quickly visualize expression gene-weighted density of one or more features
#' based on \pkg{Nebulosa}.
#'
#' @param data A Seurat object with a "harmony" reduction and UMAP computed on it.
#' @param features Character vector of feature names to plot.
#' @param ncol Number of columns in the output layout (default 2).
#' @param reduction Name of the reduction to use for plotting (default "umap.harmony").
#' @param joint Return joint density plot? By default FALSE
#' @param size Size of the geom to be plotted (e.g. point size)
#' @param pal Choose from Nebulosa's palettes, e.g. "magma", "inferno", "plasma", "viridis", "cividis".
#' @param combine Passed to \code{Nebulosa::plot_density()}.
#' @param adjustment Numeric value to adjust the bandwidth of the kernel density estimation (default 1).
#' @param ... Additional arguments passed to \code{Nebulosa::plot_density()}.
#'
#' @return A patchwork object arranging density plots in a grid.
#'
#' @importFrom Nebulosa plot_density
#' @importFrom ggplot2 coord_fixed theme_void
#' @importFrom patchwork wrap_plots
#' @importFrom leo.basic leo_log
#' @examples
#' # Assuming 'all' is a Seurat object with harmony UMAP already run
#' library(Seurat)
#' library(Nebulosa)
#' library(patchwork)
#' data <- SeuratObject::pbmc_small
#' # Plot density for CD8A and CD8B in two columns
#' plot_gw_density(data, features = c("CD3D", "CD3E"), reduction = "tsne", ncol = 2)
#' # plot joint density for CD3D and CD3E
#' plot_gw_density(data, features = c("CD3D", "CD3E"), reduction = "tsne",
#'                 joint = TRUE, combine = FALSE)
#' @export
plot_gw_density <- function(data, features, reduction = "umap.harmony",
                            size = 0.2, pal = "magma", ncol = 2,
                            joint = FALSE, combine = TRUE, adjustment = 1, ...) {
  # Check features
  available_features <- rownames(data)
  missing <- setdiff(features, available_features)
  if (length(missing) > 0) {
     features <- intersect(features, available_features)
     leo.basic::leo_log("Features not found: {paste(missing, collapse = ', ')}", level = "warning")
     if (length(features) == 0) stop("No valid features found.")
  }

  if (is.null(reduction)) reduction <- SeuratObject::DefaultDimReduc(data)

  # Check Nebulosa availability and compatibility
  use_fallback <- FALSE
  if (!requireNamespace("Nebulosa", quietly = TRUE)) {
    leo.basic::leo_log("Nebulosa not installed. Using basic ggplot2 fallback.", level = "warning")
    use_fallback <- TRUE
  }
  
  # Try running Nebulosa::plot_density
  if (!use_fallback) {
    tryCatch({
      # Attempt to run Nebulosa - test run on first feature if not joint, or all if joint
      # We just return the result directly if it works
      if (!joint) {
        plots <- lapply(features, function(gene) {
           p <- Nebulosa::plot_density(data, features = gene, reduction = reduction,
                                  method = "wkde", pal = pal, size = size,
                                  raster = TRUE) +
             ggplot2::coord_fixed() + ggplot2::theme_void()
           return(p)
        })
        if (combine) return(patchwork::wrap_plots(plots, ncol = ncol))
        return(plots)
      } else {
        leo.basic::leo_log("Plot joint density for {length(features)} features.")
        x <- Nebulosa::plot_density(data, features = features, reduction = reduction,
                                    joint = TRUE, method = "wkde", pal = pal,
                                    size = size, raster = TRUE, combine = TRUE)
        return(x)
      }
    }, error = function(e) {
      if (grepl("slot", e$message, ignore.case = TRUE) || grepl("FetchData", e$message, ignore.case = TRUE)) {
         leo.basic::leo_log("Nebulosa failed (likely Seurat v5 incompatibility). Using robust fallback.", level = "warning")
         # Logic continues below to fallback
      } else {
         stop(e) # Re-throw unexpected errors
      }
    })
    # If we reached here without returning, an error occurred in tryCatch that we caught
    use_fallback <- TRUE
  }

  # Fallback Implementation (Seurat v5 compatible)
  if (use_fallback) {
      coords <- Seurat::Embeddings(data, reduction = reduction)[, 1:2]
      colnames(coords) <- c("Dim1", "Dim2")

      make_plot <- function(feat_name, expr_vals) {
        df <- data.frame(coords, Expression = expr_vals)
        p <- ggplot2::ggplot(df, ggplot2::aes(x = Dim1, y = Dim2)) +
          ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(density), weight = .data$Expression),
                                   geom = "raster", contour = FALSE, adjust = adjustment, n = 200) +
          ggplot2::scale_fill_viridis_c(option = pal) +
          ggplot2::coord_fixed() +
          ggplot2::theme_void() +
          ggplot2::labs(title = feat_name) +
          ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5))
        return(p)
      }

      if (isTRUE(joint)) {
        expr_data <- Seurat::FetchData(data, vars = features, layer = "data")
        joint_expr <- apply(expr_data, 1, prod)
        title <- paste(features, collapse = " & ")
        return(make_plot(title, joint_expr))
      } else {
        plots <- lapply(features, function(f) {
          expr <- Seurat::FetchData(data, vars = f, layer = "data")[, 1]
          make_plot(f, expr)
        })
        if (isTRUE(combine)) return(patchwork::wrap_plots(plots, ncol = ncol))
        return(plots)
      }
  }
}

#' Highlight a cluster
#'
#' This function highlight a cluster on a dimensional reduction in NPG palette with highlight on top.
#'
#' @param obj Seurat object.
#' @param cluster_id Value to highlight (in Idents(obj) or `obj[[group.by]]`).
#' @param reduction Dimred name; default active reduction.
#' @param group.by Metadata column; NULL uses Idents.
#' @param highlight.col Highlight color; NULL -> ggsci NPG red.
#' @param other.col Background color (default "grey80").
#' @param pt.size Point size.
#' @param pt.shape Point shape (ggplot2 pch), default 16.
#' @param raster Logical or NULL; NULL auto-enable when >1e5 cells.
#' @param dpi Single numeric DPI for raster geoms (default 300).
#' @param legend Show legend (default TRUE).
#' @param legend_labels Named vector to rename legend entries; must include
#'   c("highlight","other"). Defaults to c(highlight=cluster_id, other="other").
#' @param legend_breaks Character vector to set legend order/visibility; accepts keys
#'   in c("highlight","other") or their display labels in `legend_labels`.
#'   Default c("highlight","other"). Use "highlight" only to hide "other".
#' @param legend_pos Legend position (inside); NULL uses ggplot2 default (outside, right).
#' @return ggplot object.
#'
#' @importFrom Seurat Embeddings Idents NoAxes
#' @importFrom SeuratObject DefaultDimReduc
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_void theme
#' @importFrom ggrastr geom_point_rast
#' @importFrom ggsci pal_npg
#' @importFrom leo.basic leo_log
#' @examples
#' # load demo data
#' data("pbmc_small", package = "SeuratObject")
#'
#' # 1) Basic: highlight first cluster on PCA (NPG red, other gray)
#' plot_highlight_cluster(pbmc_small, cluster_id = levels(Seurat::Idents(pbmc_small))[1],
#'                        reduction = "pca", pt.size = 0.5, pt.shape = 16, raster = FALSE, dpi = 300)
#'
#' # 2) Use UMAP (active reduction will auto-detect; here we set explicitly if available)
#' \donttest{
#' pbmc_small <- Seurat::RunUMAP(pbmc_small, reduction = "pca", dims = 1:10)
#' plot_highlight_cluster(pbmc_small, cluster_id = "0", reduction = "umap",
#'                        pt.size = 0.6, pt.shape = 16, raster = FALSE, dpi = 300)
#' }
#'
#' # 3) Custom colors (manual override): highlight blue, background light gray
#' plot_highlight_cluster(pbmc_small, cluster_id = "0", reduction = "pca",
#'                        highlight.col = "#377EB8", other.col = "grey85",
#'                        pt.size = 0.7, pt.shape = 16, raster = FALSE, dpi = 300)
#'
#' # 4) Hide legend entirely
#' plot_highlight_cluster(pbmc_small, cluster_id = "0", reduction = "pca",
#'                        legend = FALSE, pt.size = 0.6, pt.shape = 16, raster = FALSE, dpi = 300)
#'
#' # 5) Rename legend entries & order them (show highlight then other)
#' plot_highlight_cluster(pbmc_small, cluster_id = "0", reduction = "pca",
#'                        legend_labels = c(highlight = "Yes", other = "No"),
#'                        legend_breaks = c("highlight","other"),
#'                        pt.size = 0.6, pt.shape = 16, raster = FALSE, dpi = 300)
#'
#' # 6) Show only the highlight entry in legend (hide "other")
#' plot_highlight_cluster(pbmc_small, cluster_id = "0", reduction = "pca",
#'                        legend_labels = c(highlight = "Yes", other = "No"),
#'                        legend_breaks = "highlight",
#'                        pt.size = 0.6, pt.shape = 16, raster = FALSE, dpi = 300)
#'
#' # 7) Big data style: force rasterization and set DPI
#' plot_highlight_cluster(pbmc_small, cluster_id = "0", reduction = "pca",
#'                        raster = TRUE, dpi = 500, pt.size = 0.6, pt.shape = 16)
#'
#' # 8) Use a metadata column for grouping (e.g., Stage1), highlight "Control"
#' #    (replace "Stage1" with your metadata column that contains "Control")
#' pbmc_small$Stage1 <- ifelse(as.character(Seurat::Idents(pbmc_small)) == "0", "Control", "Other")
#' plot_highlight_cluster(pbmc_small, cluster_id = "Control", group.by = "Stage1",
#'                        reduction = "pca", pt.size = 0.6, pt.shape = 16,
#'                        raster = FALSE, dpi = 300)
#' # 9) Same as above but customize legend position (inside plot area, e.g. top-right)
#' plot_highlight_cluster(pbmc_small, cluster_id = "Control", group.by = "Stage1",
#'                        reduction = "pca", pt.size = 0.6, pt.shape = 16,legend_pos = c(0.8,0.8),
#'                        raster = FALSE, dpi = 300)
#' @export
plot_highlight_cluster <- function(obj, cluster_id, reduction = NULL, group.by = NULL,
                                   highlight.col = NULL, other.col = "grey80",
                                   pt.size = 0.6, pt.shape = 16,
                                   raster = NULL, dpi = 300,
                                   legend = TRUE, legend_labels = NULL, legend_breaks = NULL, legend_pos = NULL) {
  if (!is.null(group.by) && !group.by %in% colnames(obj@meta.data)) stop("`group.by` not found in meta.data.")
  if (is.null(reduction)) reduction <- SeuratObject::DefaultDimReduc(obj)

  groups <- if (is.null(group.by)) Seurat::Idents(obj) else obj[[group.by, drop = TRUE]]
  cluster_label <- as.character(cluster_id)
  cells_highlight <- names(groups)[as.character(groups) == cluster_label]
  if (length(cells_highlight) == 0) stop("No cells match `cluster_id` under current grouping.")

  coords <- Seurat::Embeddings(obj, reduction)[, 1:2]
  df <- data.frame(cell = rownames(coords), Dim1 = coords[, 1], Dim2 = coords[, 2], stringsAsFactors = FALSE)
  df_other <- df[!df$cell %in% cells_highlight, , drop = FALSE]
  df_high  <- df[ df$cell %in% cells_highlight, , drop = FALSE]

  raster_use <- if (is.null(raster)) nrow(df) > 1e5 else raster
  dpi_num <- as.numeric(dpi[1]); if (length(dpi_num) != 1 || is.na(dpi_num)) stop("`dpi` must be a single numeric.")

  if (is.null(highlight.col)) {
    highlight.col <- if (requireNamespace("ggsci", quietly = TRUE)) ggsci::pal_npg("nrc")(10)[1] else "#DE2D26"
  }
  color_map <- c(other = other.col, highlight = highlight.col)

  if (is.null(legend_labels)) legend_labels <- c(highlight = cluster_label, other = "other")
  if (!all(c("highlight","other") %in% names(legend_labels))) stop("`legend_labels` must be named with c('highlight','other').")

  to_keys <- function(x) {
    keys <- ifelse(x %in% names(legend_labels), x, names(legend_labels)[match(x, unname(legend_labels))])
    keys[!is.na(keys) & keys %in% c("highlight","other")]
  }
  if (is.null(legend_breaks)) legend_breaks <- c("highlight","other")
  legend_breaks <- to_keys(legend_breaks)

  leo.basic::leo_log("Highlight --> <{cluster_id}>: {nrow(df_high)}/{nrow(df)} cells on {reduction}; raster={raster_use}@{dpi_num}dpi.")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Dim1, y = Dim2))
  if (isTRUE(raster_use) && requireNamespace("ggrastr", quietly = TRUE)) {
    p <- p +
      ggrastr::geom_point_rast(data = df_other, ggplot2::aes(color = "other"),
                               size = pt.size, shape = pt.shape, raster.dpi = dpi_num) +
      ggrastr::geom_point_rast(data = df_high,  ggplot2::aes(color = "highlight"),
                               size = pt.size, shape = pt.shape, raster.dpi = dpi_num)
  } else {
    if (isTRUE(raster_use)) leo.basic::leo_log("ggrastr not installed; fallback to vector points.", level = "warning")
    p <- p +
      ggplot2::geom_point(data = df_other, ggplot2::aes(color = "other"),
                          size = pt.size, shape = pt.shape) +
      ggplot2::geom_point(data = df_high,  ggplot2::aes(color = "highlight"),
                          size = pt.size, shape = pt.shape)
  }

  p <- p + ggplot2::scale_color_manual(values = color_map, breaks = legend_breaks,
                                       labels = unname(legend_labels[legend_breaks]), name = NULL) +
    ggplot2::theme_void() + Seurat::NoAxes()

  if (legend) {
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, shape = pt.shape)))
  } else {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  if (!is.null(legend_pos) && legend) {
    p <- p +
      ggplot2::theme(legend.position = "inside",
                     legend.position.inside = legend_pos,
                     legend.background = ggplot2::element_rect(fill = "transparent"),
                     legend.box.background = ggplot2::element_rect(color = "transparent")
                     )
  }

  leo.basic::leo_log("Done.", level = "success")
  return(p)
}

#' Different-effect-variable Beeswarm Plot
#'
#' Visualize continuous effect values (e.g., log2FC/beta/FC) across groups with a reference dashed line.
#' Non-significant points (by p-value and/or deadband) are drawn in gray; significant points use a diverging palette.
#'
#' @param df data.frame/tibble containing grouping and effect columns
#' @param group.by character, column name for grouping
#' @param effect_col character, column name for effect (e.g., "logFC", "beta")
#' @param p_col character or NULL, column name for p-values; if provided, p >= p_thresh is treated as non-significant
#' @param p_thresh numeric, p-value threshold for significance (default 0.05)
#' @param effect_thresh numeric, reference threshold for dashed line and color midpoint (default 0)
#' @param pal_color named vector c(low, mid, high) for diverging palette
#'   (default c(low="#5062A7", mid="white", high="#BC4B59"))
#' @param log2fc_limits NULL or numeric length-2 c(L, R); if set, color scale limits
#'   are c(effect_thresh-L, effect_thresh+R)
#' @param insignificant_color character, color for non-significant/gray-zone points (default "gray80")
#' @param deadband NULL or non-negative numeric; if set, |effect - effect_thresh| <= deadband will be gray
#' @param flip_coord logical, flip coordinates to show groups vertically (default TRUE)
#' @param point_size numeric, point size (default 1)
#' @param seed NULL or integer, for reproducible quasirandom placement
#' @param ... extra args passed to ggbeeswarm::geom_quasirandom()
#'
#' @return ggplot object
#' @export
#'
#' @importFrom dplyr mutate %>% filter
#' @importFrom ggplot2 ggplot aes geom_hline scale_color_gradient2 coord_flip labs theme_classic
#' @importFrom ggplot2 guides xlab ylab
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom stats na.omit
#' @importFrom utils head
#' @importFrom tibble tibble
#'
#' @examples
#' # ---- Example 1: MiloR-like DA results ----
#'  set.seed(1)
#'  milo_df <- tibble::tibble(
#'    Nhood = paste0("n", seq_len(1200)),
#'    `Cell Type` = sample(paste0("CT", 1:6), 1200, replace = TRUE),
#'    logFC = rnorm(1200, sd = 1.2),
#'    SpatialFDR = runif(1200)
#'  )
#'  # Visualize logFC by cell types; non-sig: SpatialFDR >= 0.1; dashed line at 0
#'  p1 <- plot_dbee(milo_df, group.by = "Cell Type", effect_col = "logFC",
#'                  p_col = "SpatialFDR", p_thresh = 0.1, effect_thresh = 0,
#'                  log2fc_limits = c(-.1, .1), deadband = 0.1, point_size = 2, seed = 42)
#'  print(p1)
#'
#'  # ---- Example 2: scRNA-seq DEG-like results ----
#'  set.seed(123)
#'  deg_df <- tibble::tibble(
#'    gene = paste0("G", 1:900),
#'    cluster = sample(paste0("C", 1:5), 900, replace = TRUE),
#'    log2FC = rnorm(900, mean = rep(seq(-0.4, 0.4, length.out = 5), each = 180), sd = 1),
#'    p_val_adj = pmin(runif(900)^2, 1)
#'  )
#'  # Visualize log2FC by cluster; non-sig: p_val_adj >= 0.05; dashed line at 0
#'  p2 <- plot_dbee(deg_df, group.by = "cluster", effect_col = "log2FC",
#'                  p_col = "p_val_adj", p_thresh = 0.05, effect_thresh = 0,
#'                  pal_color = c(low = "#2C7BB6", mid = "#FFFFBF", high = "#D7191C"), flip_coord = FALSE,
#'                  log2fc_limits = NULL, deadband = 0.05, point_size = 2, seed = 7)
#'  print(p2)
plot_dbee <- function(df, group.by, effect_col, p_col = NULL, p_thresh = 0.05, effect_thresh = 0,
                      pal_color = c(low = "#5062A7", mid = "white", high = "#BC4B59"),
                      log2fc_limits = NULL, insignificant_color = "gray80", deadband = NULL,
                      flip_coord = TRUE, point_size = 1, seed = NULL, ...) {
  # basic checks
  if (!is.data.frame(df)) stop("`df` must be a data.frame / tibble.")
  if (!group.by %in% names(df)) stop(glue::glue("{group.by} not found in `df`."))
  if (!effect_col %in% names(df)) stop(glue::glue("{effect_col} not found in `df`."))
  if (!all(c("low","mid","high") %in% names(pal_color))) stop("`pal_color` must have names: low, mid, high.")
  if (!is.null(log2fc_limits)) if (length(log2fc_limits) != 2 || any(!is.finite(log2fc_limits))) stop("`log2fc_limits` must be NULL or numeric length-2.")
  if (!is.null(p_col) && !p_col %in% names(df)) stop(glue::glue("{p_col} not found in `df`."))
  if (!is.null(seed)) set.seed(seed)

  # logging
  n_rows <- nrow(df); n_groups <- length(unique(df[[group.by]]))
  leo.basic::leo_log("plot_dbee(): start; rows={n_rows}, groups={n_groups}")

  # data prep
  df2 <- dplyr::mutate(df, group_by = .data[[group.by]], effect = as.numeric(.data[[effect_col]]))
  if (!is.factor(df2$group_by)) df2$group_by <- factor(df2$group_by, levels = unique(df2$group_by))

  # non-significant / gray-zone flags
  in_gray <- rep(FALSE, nrow(df2))
  if (!is.null(p_col)) {
    pvals <- suppressWarnings(as.numeric(df2[[p_col]]))
    in_gray <- in_gray | (!is.na(pvals) & pvals >= p_thresh)
  }
  if (!is.null(deadband) && is.finite(deadband) && deadband >= 0) {
    in_gray <- in_gray | (abs(df2$effect - effect_thresh) <= deadband)
  }
  df2$in_gray <- in_gray

  # color scale limits
  limits <- NULL
  if (!is.null(log2fc_limits)) limits <- c(effect_thresh - log2fc_limits[1], effect_thresh + log2fc_limits[2])

  # plot
  p <- ggplot2::ggplot(df2, ggplot2::aes(x = group_by, y = effect))
  if (any(df2$in_gray)) p <- p + ggbeeswarm::geom_quasirandom(data = df2[df2$in_gray, , drop = FALSE],
                                                              color = insignificant_color, size = point_size, ...)
  p <- p + ggbeeswarm::geom_quasirandom(data = df2[!df2$in_gray, , drop = FALSE],
                                        ggplot2::aes(color = .data$effect), size = point_size, ...)
  p <- p + ggplot2::geom_hline(yintercept = effect_thresh, linetype = "dashed")
  p <- p + ggplot2::scale_color_gradient2(low = pal_color["low"], mid = pal_color["mid"], high = pal_color["high"],
                                          midpoint = effect_thresh, limits = limits)
  if (isTRUE(flip_coord)) p <- p + ggplot2::coord_flip()
  p <- p + ggplot2::labs(x = group.by, y = effect_col, color = effect_col) + ggplot2::theme_classic()

  leo.basic::leo_log("plot_dbee(): done", level = "success")
  return(p)
}
