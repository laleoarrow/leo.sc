#' Format marker-gene lists for (bulk) upload
#'
#' Convert a Seurat **marker table** into plain-text lines of the form `cluster_<id>:geneA,geneB,…`, one line per cluster.
#'
#' @param markers_tbl  A tibble / data.frame **containing at least** the columns `cluster`, `gene`, and the column named in `order_by`.
#' @param top_n        Integer. Number of genes to keep for each cluster (default **10**).
#' @param order_by     Column used to rank genes *within* each cluster. Accepts either a **character string** (e.g. `"avg_log2FC"`) or a *bare* column name (unquoted). Default `"avg_log2FC"`.
#' @param verbose      Logical. Print `leo_log()` messages? Default **TRUE**.
#'
#' @return A length-one character vector
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
  order_quo <- rlang::enquo(order_by)
  order_sym <- if (rlang::quo_is_symbol(order_quo)) {
    rlang::get_expr(order_quo) # bare column
  } else {
    if (!order_by %in% colnames(markers_tbl)) {
      leo.basic::leo_log(
        sprintf("Column '%s' not found in `markers_tbl`.", order_by),
        level   = "danger", verbose = verbose
      )
      stop("`order_by` column is missing.", call. = FALSE)
    }
    rlang::sym(order_by) # character → symbol
  }

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


#' Filter Seurat clusters by percentage or absolute cell count
#'
#' Filters clusters in a Seurat object based on either a percentage of total cells or an absolute count threshold.
#'
#' @param seurat_obj A Seurat object.
#' @param cluster_col Character. Name of the metadata column for clusters.
#' @param pct_threshold Numeric. Minimum fraction of total cells to keep a cluster (default 0.001).
#' @param abs_threshold Numeric or NULL. If provided, minimum absolute cell count to keep a cluster (overrides pct_threshold).
#' @return Character vector of cluster IDs to keep.
#'
#' @importFrom dplyr rename arrange filter pull desc
#' @import leo.basic
#' @export
filter_clusters_by_percent_or_cell_count <- function(seurat_obj, cluster_col,
                                                     pct_threshold = 0.001,
                                                     abs_threshold = NULL) {
  total_cells <- ncol(seurat_obj)
  cell_num_threshold <- if (!is.null(abs_threshold)) abs_threshold else pct_threshold*total_cells
  leo.basic::leo_log("Total cells: ", total_cells)
  leo.basic::leo_log("Cell threshold: ", cell_num_threshold)

  cell_counts <- as.data.frame(table(seurat_obj[[cluster_col, drop = TRUE]])) %>%
    dplyr::rename(cluster = Var1, cell_count = Freq) %>%
    dplyr::arrange(desc(cell_count))
  keep_clusters <- cell_counts %>%
    dplyr::filter(cell_count >= cell_num_threshold) %>%
    dplyr::pull(cluster)
  drop_clusters <- cell_counts %>%
    dplyr::filter(cell_count < cell_num_threshold) %>%
    dplyr::pull(cluster)
  leo.basic::leo_log("✅ Kept clusters (", length(keep_clusters), "): ", paste(keep_clusters, collapse = ", "))
  leo.basic::leo_log("❌ Dropped clusters (", length(drop_clusters), "): ", paste(drop_clusters, collapse = ", "))
  return(keep_clusters)
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

#' Calculate and plot ROGUE index for a Seurat object, with inline filtering
#'
#' Filters raw counts by minimum cells and genes, computes entropy and ROGUE metrics.
#'
#' @param obj A Seurat object with `celltype_broad` and `orig.ident` metadata.
#' @param assay Assay name containing raw counts (default "RNA").
#' @param layer Layer name for raw counts (default "counts").
#' @param min.cells Minimum cells per gene for filtering (default 10).
#' @param min.genes Minimum genes per cell for filtering (default 10).
#' @param anno_label Cell type annotation label (default "celltype_broad").
#' @param sample_label Sample label for grouping (default "orig.ident").
#' @param ... Additional arguments passed to `rogue()`.or `mc_rogue()`.
#'
#' @return An object returned by `rogue()`, invisibly.
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom Matrix colSums rowSums
#' @importFrom ROGUE SE_fun SEplot CalculateRogue rogue rogue.boxplot
#' @importFrom leo.basic leo_log
#' @export
calcROGUE <- function(obj, assay = "RNA", layer = "counts", min.cells = 10, min.genes = 10,
                      anno_label = "celltype_broad", sample_label = "orig.ident", ...) {
  matr.filter2 <- function(expr, min.cells, min.genes) {
    # This function enhances `matr.filter` in ROGUE
    if (!inherits(expr, "matrix") && !inherits(expr, "dgCMatrix")) {
      stop("`expr` must be a matrix or dgCMatrix, but got: ", class(expr), call. = FALSE)
    }
    if (length(dim(expr)) != 2) {
      stop("`expr` must be 2-dimensional, but dim(expr) = ", paste(dim(expr), collapse = "x"), call. = FALSE)
    }
    requireNamespace("Matrix", quietly = TRUE)
    gene_count <- Matrix::colSums(expr > 0, na.rm = T)
    cell_count <- Matrix::rowSums(expr > 0, na.rm = T)
    keep_cells <- cell_count >= min.cells
    keep_genes <- gene_count >= min.genes
    expr[keep_cells, keep_genes, drop = FALSE]
  }

  leo.basic::leo_log("Filtering expression matrix with min.cells = {min.cells} and min.genes = {min.genes}")
  expr <- GetAssayData(obj, assay = "RNA", slot="counts") %>%
    matr.filter2(.,  min.cells, min.genes) %>%
    as.matrix(); gc()

  leo.basic::leo_log("Filtered expression matrix and begin calculating the expression entropy for each gene")
  ent.res <- SE_fun(expr); SEplot(ent.res)

  leo.basic::leo_log("Calculating overall ROGUE value")
  val <- CalculateRogue(ent.res, platform = "UMI"); gc()
  leo.basic::leo_log(paste0("Overall ROGUE value: ", round(val, 4)), level = "success")

  leo.basic::leo_log("Calculating per-cluster/sample ROGUE")
  res <- mc_rogue(expr,
                  labels = obj@meta.data[[anno_label]],
                  samples = obj@meta.data[[sample_label]],
                  ...)
  leo.basic::leo_log("Done! `rogue.boxplot(res)` to visualize.", level = "success")
  return(res)
}

# multi-core version of rogue calculation
mc_rogue <- function(expr, labels, samples, platform = "UMI", k = NULL,
                     min.cell.n = 10, remove.outlier.n = 2, span = 0.5, r = 1,
                     filter = F, min.cells = 10, min.genes = 10, mt.method = "fdr",
                     mc.cores = 1){
  requireNamespace("pbmcapply", quietly = TRUE)
  clusters <- unique(labels)
  meta <- tibble(CellID = 1:ncol(expr), ct = labels, sample = samples)
  sample.rogue <- function(meta, cluster) {
    tmp <- meta %>% dplyr::filter(ct == cluster)
    s <- unique(samples)
    rogue <- c()
    for (i in 1:length(s)) {
      index1 <- tmp %>% dplyr::filter(sample == s[i]) %>%
        dplyr::pull(CellID)
      if (length(index1) >= min.cell.n) {
        tmp.matr <- expr[, index1]
        if (isTRUE(filter)) {
          print("Filtering out low-abundance genes and low-quality cells")
          tmp.matr <- matr.filter(tmp.matr, min.cells = min.cells,
                                  min.genes = min.genes)
        }
        else {
          tmp.matr <- tmp.matr
        }
        tmp.res <- SE_fun(tmp.matr, span = span, r = r)
        tmp.res <- ent.toli(tmp.res, tmp.matr, span = span,
                            r = r, n = remove.outlier.n)
        rogue[i] <- CalculateRogue(tmp.res, platform = platform,
                                   k = k)
      }
      else {
        rogue[i] <- NA
      }
    }
    return(rogue)
  }
  # res <- list()
  # for (i in 1:length(clusters)) {
    # res[[i]] <- sample.rogue(meta, clusters[i])
  # }
  res <- pbmcapply::pbmclapply(clusters, sample.rogue,
                               meta = meta,
                               mc.cores = mc.cores)
  res.tibble <- Reduce(rbind, res) %>% as.matrix() %>% t() %>%
    as.data.frame()
  colnames(res.tibble) <- clusters
  rownames(res.tibble) <- unique(samples)
  return(res.tibble)
}
