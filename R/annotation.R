#' Format marker-gene lists for (bulk) upload
#'
#' Convert a Seurat/Scanpy **marker table** into plain-text lines of the form
#' `cluster_<id>:geneA,geneB,…`, one line per cluster.
#' Such text can be pasted into web tools or saved to a file.
#'
#' @param markers_tbl  A tibble / data.frame **containing at least** the
#'   columns `cluster`, `gene`, and the column named in `order_by`.
#' @param top_n        Integer. Number of genes to keep for each cluster
#'   (default **10**).
#' @param order_by     Column used to rank genes *within* each cluster.
#'   Accepts either a **character string** (e.g. `"avg_log2FC"`) or a
#'   *bare* column name (unquoted). Default `"avg_log2FC"`.
#' @param verbose      Logical. Print `leo_log()` messages? Default **TRUE**.
#'
#' @return A length-one character vector: newline-separated text ready to
#'   copy or write.
#'
#' @examples
#' txt <- format_markers_for_upload(top10, top_n = 15, order_by = "p_val_adj")
#' cat(txt)
#'
#' ## tidy-eval style (bare name):
#' txt2 <- format_markers_for_upload(top10, order_by = avg_log2FC)
#'
#' @importFrom dplyr arrange desc group_by slice_head summarise pull first
#' @importFrom rlang enquo quo_is_symbol get_expr sym
#' @importFrom leo.basic leo_log
#' @export
format_markers_for_upload <- function(markers_tbl,
                                      top_n    = 10,
                                      order_by = "avg_log2FC",
                                      verbose  = TRUE) {

  # ---- Tidy-eval: allow bare names OR character strings -----------------
  order_quo <- rlang::enquo(order_by)
  order_sym <- if (rlang::quo_is_symbol(order_quo)) {
    rlang::get_expr(order_quo)          # bare column
  } else {
    if (!order_by %in% colnames(markers_tbl)) {
      leo.basic::leo_log(
        sprintf("Column '%s' not found in `markers_tbl`.", order_by),
        level   = "danger", verbose = verbose
      )
      stop("`order_by` column is missing.", call. = FALSE)
    }
    rlang::sym(order_by)                # character → symbol
  }

  # ---- Main pipeline ----------------------------------------------------
  txt <- markers_tbl |>
    dplyr::arrange(cluster, dplyr::desc(!!order_sym)) |>
    dplyr::group_by(cluster) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::summarise(
      line = paste0(
        "cluster_", dplyr::first(cluster), ":",
        paste(gene, collapse = ",")
      ),
      .groups = "drop"
    ) |>
    dplyr::pull(line) |>
    paste(collapse = "\n")

  leo.basic::leo_log("Website for uploading: http://xteam.xbio.top/ACT/index.jsp",
                     level = "info", verbose = verbose)
  leo.basic::leo_log("Marker list formatted ✔", level = "success", verbose = verbose)
  return(txt)
}

#' Sort string-based cluster labels numerically and refactor
#'
#' @param seurat_obj A Seurat object.
#' @param cluster_col A character string. The name of the metadata column containing cluster labels as strings.
#' @return The Seurat object with the specified cluster column re-factored so that its levels are sorted numerically.
#' @examples
#' # seurat_obj <- sort_string_numeric_clusters(seurat_obj, "cluster_labels")
#' @importFrom leo.basic leo_log
#' @export
sort_string_numeric_clusters <- function(seurat_obj, cluster_col) {
  clust_vec <- as.character(seurat_obj[[cluster_col]])
  new_levels <- as.character(sort(as.numeric(unique(clust_vec))))
  leo.basic::leo_log("ℹ️ New ordered levels for '", cluster_col, "': ", paste(new_levels, collapse = ", "))
  seurat_obj[[cluster_col]] <- factor(clust_vec, levels = new_levels)
  return(seurat_obj)
}

#' Get cluster counts sorted descending
#'
#' @param seurat_obj A Seurat object.
#' @param cluster_col String. Name of the metadata column containing cluster labels.
#' @return A data.frame with columns `cluster` and `count`, sorted by `count` descending.
#' @importFrom dplyr count arrange desc
#' @importFrom rlang sym
#' @importFrom leo.basic leo_log
#' @examples
#' # Suppose your Seurat object is `all` and cluster column is "harmony_clusters"
#' # get_cluster_counts(all, "harmony_clusters")
#' @export
get_cluster_counts <- function(seurat_obj, cluster_col) {
  meta_df <- seurat_obj@meta.data
  counts_df <- meta_df %>%
    dplyr::count(!!rlang::sym(cluster_col), name = "count") %>%
    dplyr::arrange(dplyr::desc(count)) %>%
    stats::setNames(c("cluster", "count"))
  return(counts_df)
}


#' Filter clusters by size (proportion or absolute count)
#'
#' Reports cluster counts (descending) then returns the labels of clusters meeting
#' either a minimum fraction of total cells or a user‐provided absolute cell count.
#'
#' @param seurat_obj A Seurat object.
#' @param cluster_col Character. Metadata column with cluster labels.
#' @param prop_threshold Numeric. Minimum fraction of total cells (e.g., 0.001 = 0.1%).
#'   Ignored if \code{abs_threshold} is supplied. Default: 0.001.
#' @param abs_threshold Numeric. Minimum absolute cell count. Overrides
#'   \code{prop_threshold} if not \code{NULL}. Default: \code{NULL}.
#' @param verbose Logical. Whether to print \code{leo.log} messages. Default: TRUE.
#' @return A character vector of cluster labels to keep.
#' @seealso \code{\link{get_cluster_counts}}
#' @importFrom dplyr count arrange desc filter pull
#' @importFrom rlang sym
#' @examples
#' # Keep clusters comprising ≥0.5% of cells
#' keep_pct <- filter_clusters_by_size(all,
#'                                     "harmony_clusters",
#'                                     prop_threshold = 0.005)
#' # Keep clusters with at least 100 cells
#' keep_abs <- filter_clusters_by_size(all,
#'                                     "harmony_clusters",
#'                                     abs_threshold = 100)
#' @export
filter_clusters_by_size <- function(seurat_obj,
                                    cluster_col,
                                    prop_threshold = 0.001,
                                    abs_threshold = NULL,
                                    verbose = TRUE) {
  # 1. report counts
  counts_df <- get_cluster_counts(seurat_obj, cluster_col)
  leo.basic::leo_log("Cluster counts (descending):", level = "info", verbose = verbose)
  print(counts_df)

  # 2. determine cutoff
  total <- ncol(seurat_obj)
  cutoff <- if (!is.null(abs_threshold)) abs_threshold else prop_threshold * total
  leo.basic::leo_log(
    paste0("Cutoff = ", cutoff,
           if (is.null(abs_threshold)) paste0(" (", prop_threshold*100, "%)") else " cells"),
    level = "info", verbose = verbose
  )

  # 3. filter clusters
  keep <- counts_df %>%
    dplyr::filter(count >= cutoff) %>%
    dplyr::pull(!!rlang::sym(cluster_col))
  leo.basic::leo_log(paste0("Kept clusters: ", paste(keep, collapse = ", ")),
                     level = "success", verbose = verbose)

  return(keep)
}

#' Highlight one cluster on UMAP
#'
#' @param obj           A Seurat object
#' @param cluster_id    The cluster to highlight (matching levels of `group.by`)
#' @param reduction     Reduction slot (default "umap.harmony")
#' @param group.by      Metadata column with cluster IDs (default "harmony_clusters")
#' @param highlight.col Color for the highlighted cluster (default "red")
#' @param other.col     Color for all other clusters (default "grey80")
#' @param pt.size       Point size (default 1)
#' @param raster        Logical, rasterize points for speed (default TRUE)
#' @param dpi           Rasterization DPI (default 72)
#' @return ggplot2 object
#' @importFrom Seurat Embeddings
#' @importFrom ggplot2 ggplot aes labs theme_void
#' @importFrom ggrastr geom_point_rast
#' @examples
#' # Highlight cluster 34 on the harmony UMAP
#' highlightCluster(all, cluster_id = "34")
#' @export
highlightCluster <- function(obj,
                             cluster_id,
                             reduction     = "umap.harmony",
                             group.by      = "harmony_clusters",
                             highlight.col = "red",
                             other.col     = "grey80",
                             pt.size       = 1,
                             raster        = TRUE,
                             dpi           = 72) {
  coords <- Seurat::Embeddings(obj, reduction)[, 1:2]
  df <- data.frame(
    Dim1 = coords[,1],
    Dim2 = coords[,2],
    grp  = obj@meta.data[[group.by]]
  )
  df_other <- df[df$grp != cluster_id, ]
  df_high  <- df[df$grp == cluster_id, ]

  geom_fn <- if (raster) {
    function(...) ggrastr::geom_point_rast(..., raster.dpi = dpi)
  } else {
    ggplot2::geom_point
  }

  ggplot2::ggplot() +
    geom_fn(data = df_other, aes(Dim1, Dim2),
            color = other.col, size = pt.size) +
    geom_fn(data = df_high,  aes(Dim1, Dim2),
            color = highlight.col, size = pt.size) +
    labs(title = paste0("Cluster ", cluster_id, " highlighted")) +
    theme_void()
}
