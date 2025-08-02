#' get colnames of srt obj meta data name
#'
#' @param srt Seurat object
#' @param pattern Pattern to match in column names, default is "RNA_snn_res"
#'
#' @returns A vector of column names that match the pattern
#' @export
#'
#' @examples
#' srt %>% get_meta_colnames()
get_meta_colnames <- function(srt, pattern = "RNA_snn_res") {
  res_cols <- grep("^RNA_snn_res\\.", colnames(srt@meta.data), value = TRUE)
  if (pattern == "RNA_snn_res") res_cols <- res_cols |> (\(x) x[order(as.numeric(sub("^RNA_snn_res\\.", "", x)))])()
  return(res_cols)
}
