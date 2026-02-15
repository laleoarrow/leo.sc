#' Format marker-gene lists for (bulk) upload
#'
#' Convert a Seurat **marker table** into plain-text lines of the form
#' `cluster_<id>:geneA,geneB,…`, one line per cluster.
#'
#' @param markers_tbl  A tibble / data.frame **containing at least** the columns
#'   `cluster`, `gene`, and the column named in `order_by`.
#' @param top_n        Integer. Number of genes to keep for each cluster (default **10**).
#' @param order_by     Column used to rank genes *within* each cluster. Accepts either a
#'   **character string** (e.g. `"avg_log2FC"`) or a *bare* column name (unquoted).
#'   Default `"avg_log2FC"`.
#' @param verbose      Logical. Print `leo_log()` messages? Default **TRUE**.
#'
#' @return A length-one character vector
#'
#' @examples
#' top10 <- data.frame(
#'   cluster = rep(1:2, each = 2),
#'   gene = c("GeneA", "GeneB", "GeneC", "GeneD"),
#'   p_val_adj = c(0.01, 0.05, 0.001, 0.02),
#'   avg_log2FC = c(1.5, 1.2, 2.0, 1.8),
#'   avg_logFC = c(1.5, 1.2, 2.0, 1.8)
#' )
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
    dplyr::pull("line") |>
    paste(collapse = "\n")

  leo.basic::leo_log("Website for uploading: http://xteam.xbio.top/ACT/index.jsp",
                     level = "info", verbose = verbose)
  leo.basic::leo_log("Backup website: http://biocc.hrbmu.edu.cn/ACT/index.jsp",
                     level = "info", verbose = verbose)
  leo.basic::leo_log("Marker list formatted", level = "success", verbose = verbose)
  return(txt)
}

#' Sort string-based cluster labels numerically and refactor
#'
#' @param seurat_obj A Seurat object.
#' @param cluster_col A character string. The name of the metadata column 
#'                    containing cluster labels as strings.
#' @return The Seurat object with the specified cluster column re-factored so 
#'         that its levels are sorted numerically.
#' @examples
#' # seurat_obj <- sort_string_numeric_clusters(seurat_obj, "cluster_labels")
#' @importFrom leo.basic leo_log
#' @export
sort_string_numeric_clusters <- function(seurat_obj, cluster_col) {
  clust_vec   <- as.character(seurat_obj@meta.data[[cluster_col]])
  new_levels  <- as.character(sort(as.numeric(unique(clust_vec))))
  leo.basic::leo_log("[i] New ordered levels for '", cluster_col, "': ",
                     paste(new_levels, collapse = ", "))
  seurat_obj@meta.data[[cluster_col]] <- factor(clust_vec, levels = new_levels)
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
#' Filters clusters in a Seurat object based on either a percentage of total cells
#' or an absolute count threshold.
#'
#' @param seurat_obj A Seurat object.
#' @param cluster_col Character. Name of the metadata column for clusters.
#' @param pct_threshold Numeric. Minimum fraction of total cells to keep a cluster
#'                      (default 0.001).
#' @param abs_threshold Numeric or NULL. If provided, minimum absolute cell count
#'                      to keep a cluster (overrides pct_threshold).
#' @return Character vector of cluster IDs to keep.
#'
#' @importFrom dplyr rename arrange filter pull desc
#' @importFrom leo.basic leo_log
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
  leo.basic::leo_log("⬇️ The number of cell in each clusters️"); print(sort(table(seurat_obj[[cluster_col]]), decreasing = T))
  return(keep_clusters)
}

#' Calculate and plot ROGUE index for a Seurat object, with inline filtering
#'
#' Filters raw counts by minimum cells and genes, computes entropy and ROGUE metrics.
#'
#' @param obj A Seurat object with `celltype_broad` and `orig.ident` metadata.
#' @param assay Assay name containing raw counts (default "RNA").
#' @param layer Layer name for raw counts (default "counts").
#' @param downsample Number of cells to downsample per group (default 3000).
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
#' @export
calcROGUE <- function(obj, assay = "RNA", layer = "counts", downsample = 3000,
                      min.cells = 10, min.genes = 10,
                      anno_label = "celltype_broad", sample_label = "orig.ident", ...) {
  t0 <- Sys.time()
  requireNamespace("Matrix", quietly = TRUE)
  requireNamespace("ROGUE", quietly = TRUE)
  if (!is.null(downsample)) {
    leo.basic::leo_log("Downsampling to", downsample, "cells per", anno_label, "group")
    Idents(obj) <- obj[[anno_label]][,1]
    obj <- subset(obj, downsample = downsample)
  }

  matr.filter2 <- function(expr, min.cells, min.genes) {
    # This function enhances `matr.filter` in ROGUE
    if (!inherits(expr, "matrix") && !inherits(expr, "dgCMatrix")) {
      stop("`expr` must be a matrix or dgCMatrix, but got: ", class(expr), call. = FALSE)
    }
    if (length(dim(expr)) != 2) {
      stop("`expr` must be 2-dimensional, but dim(expr) = ", paste(dim(expr), collapse = "x"), call. = FALSE)
    }
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
  leo.basic::leo_log("Done! `rogue.boxplot(rogue.res)` to visualize.", level = "success")
  leo.basic::leo_time_elapsed(t0)
  return(res)
}

# multi-core version of rogue calculation
mc_rogue <- function(expr, labels, samples, platform = "UMI", k = NULL,
                     min.cell.n = 10, remove.outlier.n = 2, span = 0.5, r = 1,
                     filter = FALSE, min.cells = 10, min.genes = 10, mt.method = "fdr",
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
          tmp.matr <- ROGUE:::matr.filter(tmp.matr, min.cells = min.cells,
                                  min.genes = min.genes)
        }
        else {
          tmp.matr <- tmp.matr
        }
        tmp.res <- SE_fun(tmp.matr, span = span, r = r)
        tmp.res <- ROGUE:::ent.toli(tmp.res, tmp.matr, span = span,
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

  # 20250927: update the logics here ---
  # res.tibble <- Reduce(rbind, res) %>% as.matrix() %>% t() %>%
  #   as.data.frame()
  # # print(length(clusters))
  # # print(length(unique(samples)))
  # colnames(res.tibble) <- clusters
  # rownames(res.tibble) <- unique(samples)
  # 20250927: replace the above to the below---
  mat <- do.call(cbind, res)
  colnames(mat) <- clusters
  rownames(mat) <- s
  res.tibble <- as.data.frame(mat)
  # ---
  return(res.tibble)
}

#' Score sc_obj with signature list
#'
#' This function calculates module scores for a Seurat object based on a list of gene sets.
#'
#' @param sc_obj A Seurat object.
#' @param signature_list list. A list of gene sets, where each element is a character vector of gene names.
#' @param seed An integer seed for reproducibility. Default is 1.
#'
#' @returns A Seurat object with module scores added to the metadata.
#' @export
#' @importFrom Seurat AddModuleScore
#' @importFrom leo.basic leo_log
#' @examples
#' \dontrun{
#' sc_obj <- score_signature(sc_obj, signature_list))
#' }
score_signature <- function(sc_obj, signature_list, seed = 1) {
  for (nm in names(signature_list)) {
    leo.basic::leo_log("Calculating {nm}")
    old_names <- colnames(sc_obj@meta.data)
    sc_obj <- AddModuleScore(sc_obj, features = list(signature_list[[nm]]),
                             name = nm, ctrl = 100, seed = seed, search = FALSE)
    # AddModuleScore gives "nm1" not what we expected---nm. So change it back to nm.
    new_names <- setdiff(colnames(sc_obj@meta.data), old_names)
    if (length(new_names) != 1) {
      leo.basic::leo_log("Clear old calculated signature scores for {nm} before re-calculating.", level = "warning")
      leo.basic::leo_log("Use: sc_obj@meta.data = sc_obj@meta.data[,-c(xxx)]", level = "warning")
      stop()
    }

    colnames(sc_obj@meta.data)[match(new_names, colnames(sc_obj@meta.data))] <- nm
    colnames(sc_obj@meta.data)[match(old_names, colnames(sc_obj@meta.data))] <- old_names
  }
  return(sc_obj)
}

#' Plot heatmap of signature scores
#'
#' This function generates a heatmap of signature scores. You should do this after score_signature().
#'
#' @param sc_obj A Seurat object.
#' @param signature_list list. A list of gene sets, where each element is a character vector of gene names.
#' @param group A character string specifying the metadata column to group by. (e.g., group = "RNA_snn_res.0.6")
#' @param group_prefix A character string to prefix the column names in the heatmap. Default is NULL.
#' @param scale Whether to scale the matrix. Options are "none", "row", or "column". Default is "none".
#' @param signature_cat Named vector. c("rowname1" = "category1", "rowname2" = "category2"). Provide to each row a category.
#' @param signature_cat_col Named vector. c("category1" = "color1", "category2" = "color2"). Provide to each category a color.
#' @param save_path Path to save the heatmap PDF (Only support pdf). Default is "./signature.pdf".
#' @param width Width of the saved heatmap PDF. Default is 6.
#' @param height Height of the saved heatmap PDF. Default is 6.
#' @param heatmap_title Character. Title for the heatmap. Default: NULL.
#'
#' @returns heatmap obj.
#' @importFrom dplyr mutate group_by summarise across
#' @importFrom tibble column_to_rownames
#' @importFrom matrixStats rowSds colSds
#' @importFrom grDevices adjustcolor hcl.colors
#' @importFrom grid unit gpar
#' @importFrom ComplexHeatmap Heatmap rowAnnotation Legend draw anno_block
#' @importFrom circlize colorRamp2
#' @export
#' @examples
#' \dontrun{
#' signature_category <- c("Naïve"                        = "Differentiation",
#' "Activation/Effector function" = "Differentiation",
#' "Exhaustion"                   = "Differentiation",
#' "TCR Signaling"                = "Function",
#' "Cytotoxicity"                 = "Function",
#' "Cytokine/Cytokine receptor"   = "Function",
#' "Chemokine/Chemokine receptor" = "Function",
#' "Senescence"                   = "Function",
#' "Anergy"                       = "Function",
#' "NFKB Signaling"               = "Function",
#' "Stress response"              = "Function",
#' "MAPK Signaling"               = "Function",
#' "Adhesion"                     = "Function",
#' "IFN Response"                 = "Function",
#' "Oxidative phosphorylation"    = "Metabolism",
#' "Glycolysis"                   = "Metabolism",
#' "Fatty acid metabolism"        = "Metabolism",
#' "Pro-apoptosis"                = "Apoptosis",
#' "Anti-apoptosis"               = "Apoptosis")
#' signature_category_color <- c("Differentiation" = "#E49446",
#'                               "Function"        = "#4BA5DB",
#'                               "Metabolism"      = "#815F43",
#'                               "Apoptosis"       = "#B12F3A")
#' plot_score_signature_heatmap(cd8, cd8_signature, group = "RNA_snn_res.0.6",
#'                              group_prefix = "CD8_c", scale = "row",
#'                              signature_cat = signature_category,
#'                              signature_cat_col = signature_category_color,
#'                              save_path = "./figure/sc/annotation/tnk/cd8/demo.signature.pdf",
#'                              width = 6, height = 6)
#'}
plot_score_signature_heatmap <- function(sc_obj, signature_list, group, group_prefix = NULL,
                                         scale = "none", signature_cat, signature_cat_col, heatmap_title = NULL,
                                         save_path = "./signature.pdf", width = 6, height = 6){
  require(ComplexHeatmap); require(circlize)

  mat <- sc_obj@meta.data[, names(signature_list)] %>%
    mutate(cluster = sc_obj@meta.data[[group]]) %>%
    group_by(cluster) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    tibble::column_to_rownames("cluster") %>%
    t()

  if (!is.null(group_prefix)) colnames(mat) <- paste0(group_prefix, colnames(mat))
  if ("row" %in% scale) { # adpated from pheatmap's scale method
    if (any(is.na(mat))) {
      mat = (mat - rowMeans(mat, na.rm = TRUE))/rowSds(mat, na.rm = TRUE)
    }
    else {
      mat = t(scale(t(mat)))
    }
  }
  else if ("column" %in% scale) {
    if (any(is.na(mat))) {
      mat = t((t(mat) - colMeans(mat, na.rm = TRUE))/colSds(mat, na.rm = TRUE))
    }
    else {
      mat = scale(mat)
    }
  }

  rowname_col <- signature_cat[rownames(mat)]
  rowname_colors <- grDevices::adjustcolor(signature_cat_col[rowname_col], alpha.f = .5)
  names(rowname_colors) <- rownames(mat)

  left_ha <- rowAnnotation(Category = anno_block(gp = gpar(fill = signature_cat_col)),
                           show_annotation_name = FALSE,
                           width = unit(4, "mm"))

  ht <- Heatmap(mat, name = "Signature\nscore",
                col = circlize::colorRamp2(c(-1, 0, 1), c("#5DB3E6","white","#E7717D")),
                border = "black",
                rect_gp = gpar(col = "black", lwd = 1),
                cluster_rows = FALSE, cluster_columns = TRUE,
                left_annotation = left_ha,
                row_split = factor(signature_cat[rownames(mat)],
                                   levels = names(signature_cat_col)),
                row_gap = unit(3, "pt"),
                row_title_gp = gpar(fontsize = 0),
                row_names_gp = gpar(col = "black", fontsize = 10, lwd = 0, lineheight = 0,
                                    fill = rowname_colors),
                row_names_side = "right", row_names_rot = 0,
                column_names_rot= 90,
                heatmap_legend_param = list(title = "Signature\nscore")
  )
  lgd_cat <- Legend(title = "Category", labels = names(signature_cat_col),
                    legend_gp = gpar(fill = signature_cat_col, col = NA), ncol = 1)
  pdf(save_path, width = width, height = height)
  draw(ht, padding = unit(c(6, 4, 4, 4), "mm"), merge_legend = TRUE, annotation_legend_list = list(lgd_cat))
  if (!is.null(heatmap_title)) { grid.text(heatmap_title,
                                           x = unit(0, "npc") + unit(3, "mm"),
                                           y = unit(1, "npc") - unit(2, "mm"),
                                           just = c("left", "top"),
                                           gp = gpar(fontsize = 14, fontface = "bold"))
    }
  dev.off()
}

#' Locate most distinguishing markers between two clusters
#'
#' @param sc_obj       Seurat object
#' @param ident1       First cluster id (string)
#' @param ident2       Second cluster id (string)
#' @param assay        Assay name, default "RNA"
#' @param test.use     DE test, one of "wilcox","roc","MAST"; default "wilcox"
#' @param pval.adj     Adjusted p‐value cutoff, default 0.05
#' @param logfc        |log2FC| cutoff, default 0.5
#' @param min.pct1     Minimum detection pct in marker group, default 0.5
#' @param max.pct2     Maximum detection pct in other group, default 0.2
#' @param return_top_n Number of top markers to return for each group, default 1.
#'
#' @return named vector c(marker1,marker2)
#' @export
#' @importFrom Seurat FindMarkers Idents DefaultAssay
#' @importFrom dplyr arrange slice_head
#' @examples
#' \dontrun{
#' tops <- locate_most_different_g_in_2_group(pbmc, "seurat_clusters", "4","5")
#' }
locate_most_different_g_in_2_group <- function(sc_obj, ident.1, ident.2, assay = "RNA",
                                               test.use = "wilcox", pval.adj = 0.05,
                                               logfc = 0.5, min.pct1 = 0.5, max.pct2 = 0.2,
                                               return_top_n = 1) {
  DefaultAssay(sc_obj) <- assay
  deg <- FindMarkers(sc_obj, ident.1 = ident.1, ident.2 = ident.2,
                     assay = assay, test.use = test.use,
                     logfc.threshold = logfc, min.pct = min.pct1)
  # marker for ident.1
  m1 <- deg[ deg$p_val_adj < pval.adj & abs(deg$avg_log2FC) > logfc &
               deg$pct.1 >= min.pct1 & deg$pct.2 <= max.pct2, ]
  leo.basic::leo_log("Found ", nrow(m1), " markers for ident.1: ", ident.1)
  # marker for ident.2
  m2 <- deg[ deg$p_val_adj < pval.adj & abs(deg$avg_log2FC) > logfc &
               deg$pct.2 >= min.pct1 & deg$pct.1 <= max.pct2, ]
  leo.basic::leo_log("Found ", nrow(m2), " markers for ident.2: ", ident.2)

  if (nrow(m1) == 0 && nrow(m2) == 0) {
    leo.basic::leo_log("No markers found for either ident.1 or ident.2.", level = "warning")
    leo.basic::leo_log("Returning the DEG for you to mannually check", level = "warning")
    return(deg)
  }

  if (return_top_n == 1) {
    top1 <- if(nrow(m1)) rownames(m1)[which.max(m1$avg_log2FC)] else NA
    top2 <- if(nrow(m2)) rownames(m2)[which.min(m2$avg_log2FC)] else NA
    return(list(ident.1 = top1, ident.2 = top2))
  }

  m1 <- m1 %>% arrange(desc(avg_log2FC)) %>% slice_head(n = return_top_n)
  m2 <- m2 %>% arrange(avg_log2FC) %>% slice_head(n = return_top_n)
  return(list(ident.1 = m1, ident.2 = m2))
}

#' Give me marker!
#'
#' This function gives the most distinguishing marker for a cluster
#'
#' @param markers A data frame or tibble containing marker genes with columns
#'   `cluster`, `pct.1`, and `pct.2`.
#' @param cluster_of_interest The cluster number for which to find the most
#'   distinguishing marker. Default is 1.
#'
#' @returns The most distinguishing marker for the specified cluster.
#' @export
gimme_marker <- function(markers, cluster_of_interest = 1) {
  tmp <- markers %>% filter(cluster == cluster_of_interest & pct.1>0.5 & pct.2<0.5) %>% arrange(pct.2 - pct.1) %>% head()
  return(tmp)
}
