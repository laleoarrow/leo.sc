#' Get specified metadata columns from a Seurat object
#'
#' @param srt Seurat object.
#' @param pattern Pattern to match in column names, default is "RNA_snn_res".
#'
#' @return A character vector of column names that match the pattern.
#' @export
#'
#' @examples
#' library(Seurat)
#' # Create a dummy Seurat object with clustering columns
#' srt <- SeuratObject::pbmc_small
#' srt$RNA_snn_res.0.8 <- Idents(srt)
#' metadata_get_colnames(srt)
metadata_get_colnames <- function(srt, pattern = "RNA_snn_res") {
  res_cols <- base::grep(base::paste0("^", pattern, "\\."), base::colnames(srt@meta.data), value = TRUE)
  if (pattern == "RNA_snn_res") {
    # Sort resolutions numerically
    res_nums <- base::as.numeric(base::sub("^RNA_snn_res\\.", "", res_cols))
    res_cols <- res_cols[base::order(res_nums)]
  }
  return(res_cols)
}


#' Keep only selected metadata columns in a Seurat object
#'
#' @param seu Seurat object.
#' @param ... Bare metadata names or strings (e.g., cell_anno, "cell_anno").
#' @param cols Optional character vector of metadata column names to keep.
#' @param select A tidyselect expression (e.g., starts_with("predicted")).
#' @param drop_na Logical; if TRUE, drop cells with NA in the kept columns.
#' @param ignore_missing Logical; if TRUE, ignore missing columns and notify.
#'
#' @return Seurat object with trimmed meta.data.
#'
#' @importFrom rlang dots_list is_string is_symbol as_string enexpr is_null
#' @importFrom tidyselect eval_select
#' @importFrom stats complete.cases
#' @importFrom cli cli_alert_info
#' @export
metadata_keep <- function(seu, ..., cols = NULL, select = NULL, drop_na = FALSE, ignore_missing = FALSE) {
  stopifnot(inherits(seu, "Seurat"))
  md <- seu@meta.data
  md_names <- colnames(md)

  dots <- rlang::dots_list(..., .ignore_empty = "all")
  explicit_cols <- character(0)
  if (length(dots) > 0) {
    for (i in seq_along(dots)) {
      if (rlang::is_string(dots[[i]])) explicit_cols <- c(explicit_cols, dots[[i]])
      if (rlang::is_symbol(dots[[i]])) explicit_cols <- c(explicit_cols, rlang::as_string(dots[[i]]))
    }
  }

  if (!is.null(cols)) {
    if (!is.character(cols)) stop("`cols` must be a character vector.")
    explicit_cols <- unique(c(explicit_cols, cols))
  }

  missing_cols <- setdiff(explicit_cols, md_names)
  if (length(missing_cols) > 0) {
    msg <- sprintf("metadata_keep(): not in meta.data (ignored): %s", paste(missing_cols, collapse = ", "))
    if (!ignore_missing) stop(msg)
    cli::cli_alert_info("{msg}")
  }
  keep_explicit <- intersect(explicit_cols, md_names)

  select_expr <- rlang::enexpr(select)
  keep_select <- character(0)
  if (!rlang::is_null(select_expr)) {
    keep_select <- names(tidyselect::eval_select(select_expr, md))
  }

  cols_keep <- unique(c(keep_explicit, keep_select))
  if (length(cols_keep) == 0) stop("metadata_keep(): no metadata columns selected to keep.")

  if (drop_na) {
    keep_cells <- rownames(md)[stats::complete.cases(md[, cols_keep, drop = FALSE])]
    n_drop <- nrow(md) - length(keep_cells)
    if (n_drop > 0) cli::cli_alert_info("metadata_keep(): dropping {n_drop} cells with NA in kept metadata.")
    seu <- Seurat::subset(seu, cells = keep_cells)
    md <- seu@meta.data
  }

  seu@meta.data <- md[, cols_keep, drop = FALSE]
  leo_log("metadata_keep(): kept {length(cols_keep)} metadata cols; drop_na={drop_na}")
  seu
}

#' Drop specified metadata columns from a Seurat object
#'
#' @param seu Seurat object.
#' @param ... Bare metadata names or strings (e.g., seurat_clusters, "seurat_clusters").
#' @param cols Optional character vector of metadata column names to drop.
#' @param select A tidyselect expression (e.g., starts_with("predicted")).
#' @param ignore_missing Logical; if TRUE, ignore missing columns and notify.
#'
#' @return Seurat object with columns removed from meta.data.
#'
#' @importFrom rlang dots_list is_string is_symbol as_string enexpr is_null
#' @importFrom tidyselect eval_select
#' @importFrom cli cli_alert_info
#' @export
metadata_drop <- function(seu, ..., cols = NULL, select = NULL, ignore_missing = TRUE) {
  stopifnot(inherits(seu, "Seurat"))
  md <- seu@meta.data
  md_names <- colnames(md)

  dots <- rlang::dots_list(..., .ignore_empty = "all")
  explicit_cols <- character(0)
  if (length(dots) > 0) {
    for (i in seq_along(dots)) {
      if (rlang::is_string(dots[[i]])) explicit_cols <- c(explicit_cols, dots[[i]])
      if (rlang::is_symbol(dots[[i]])) explicit_cols <- c(explicit_cols, rlang::as_string(dots[[i]]))
    }
  }

  if (!is.null(cols)) {
    if (!is.character(cols)) stop("`cols` must be a character vector.")
    explicit_cols <- unique(c(explicit_cols, cols))
  }

  missing_cols <- setdiff(explicit_cols, md_names)
  if (length(missing_cols) > 0) {
    msg <- sprintf("metadata_drop(): not in meta.data (ignored): %s", paste(missing_cols, collapse = ", "))
    if (!ignore_missing) stop(msg)
    cli::cli_alert_info("{msg}")
  }
  drop_explicit <- intersect(explicit_cols, md_names)

  select_expr <- rlang::enexpr(select)
  drop_select <- character(0)
  if (!rlang::is_null(select_expr)) {
    drop_select <- tryCatch(
      names(tidyselect::eval_select(select_expr, md)),
      error = function(e) {
        if (!ignore_missing) stop(conditionMessage(e))
        cli::cli_alert_info("metadata_drop(): selector error ignored: {conditionMessage(e)}")
        character(0)
      }
    )
  }

  cols_drop <- unique(c(drop_explicit, drop_select))
  if (length(cols_drop) == 0) {
    cli::cli_alert_info("metadata_drop(): no metadata columns matched; unchanged.")
    return(seu)
  }

  keep_cols <- setdiff(md_names, cols_drop)
  seu@meta.data <- md[, keep_cols, drop = FALSE]
  leo_log("metadata_drop(): dropped {length(cols_drop)} metadata cols; remaining={length(keep_cols)}")
  seu
}
