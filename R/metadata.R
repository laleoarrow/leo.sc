#' get colnames of srt obj meta data name
#'
#' @param srt Seurat object
#' @param pattern Pattern to match in column names, default is "RNA_snn_res"
#'
#' @returns A vector of column names that match the pattern
#' @export
#'
#' @examples
#' library(Seurat)
#' # Create a dummy Seurat object with clustering columns
#' srt <- SeuratObject::pbmc_small
#' srt$RNA_snn_res.0.8 <- Idents(srt)
#' get_meta_colnames(srt)
get_meta_colnames <- function(srt, pattern = "RNA_snn_res") {
  res_cols <- grep("^RNA_snn_res\\.", colnames(srt@meta.data), value = TRUE)
  if (pattern == "RNA_snn_res") {
    # Sort resolutions numerically
    res_nums <- as.numeric(sub("^RNA_snn_res\\.", "", res_cols))
    res_cols <- res_cols[order(res_nums)]
  }
  return(res_cols)
}
