#' Inspection of gene in multi-batch
#'
#' Batch-wise gene statistics checking, average expression, and common-gene subsetting
#'
#' @param srt        A Seurat object.
#' @param batch_col  Metadata column indicating batch/sample. Default \code{"Batch"}.
#' @param slot       Expression slot to inspect: \code{"counts"} for presence; \code{"data"} for averaged expression.
#' @param genes      Character vector of genes (only for \code{check_gene_avg_in_multi_batch}).
#' @param common_genes Character vector of shared genes (only for \code{subset_common_gene_in_multi_batch}).
#' @param assay      Assay to subset. Default \code{"RNA"}.
#'
#' @name   check_gene_stats_in_multi_batch
#' @rdname check_gene_stats_in_multi_batch
#'
#' @return
#' * \code{check_gene_stats_in_multi_batch}: list with a tibble \code{stats} and two
#'   character vectors \code{common}, \code{drop}.
#' * \code{check_gene_avg_in_multi_batch}: tibble of mean expression (rows = genes,
#'   cols = batches).
#' * \code{subset_common_gene_in_multi_batch}: a Seurat object containing only the
#'   genes shared by all batches.
#'
#' @examples
#' \dontrun{
#' ## --------------------------------------------------------------------
#' ## Tutorial: reconcile gene sets across multiple batches
#' ## --------------------------------------------------------------------
#'
#' # Assume `cd4` is a Seurat object with metadata column "Batch"
#'
#' # Step 1 - inspect gene overlap/uniqueness across batches
#' res <- check_gene_stats_in_multi_batch(cd4)  # returns stats, common, drop
#' res$stats                                    # view the tibble summary
#'
#' # Step 2 - examine average expression of genes missing from >=1 batch
#' avg <- check_gene_avg_in_multi_batch(cd4, res$drop); head(avg); colSums(avg[-1])
#'
#' # Step 3 - subset the Seurat object to the intersection gene set
#' cd4 <- subset_common_gene_in_multi_batch(cd4, res$common)
#' }
#' @export
check_gene_stats_in_multi_batch <- function(srt, batch_col = "Batch", slot = "counts") {
  gs   <- SplitObject(srt, split.by = batch_col) |>
    lapply(\(x) rownames(x)[Matrix::rowSums(GetAssayData(x, slot = slot) > 0) > 0])
  cmn  <- Reduce(intersect, gs)
  uni  <- Reduce(union, gs)
  n_unq <- lengths(lapply(seq_along(gs), \(i) setdiff(gs[[i]], Reduce(union, gs[-i]))))
  stats <- tibble::tibble(!!batch_col := names(gs),
                          n_genes = lengths(gs),
                          n_unique = n_unq,
                          n_common = length(cmn))
  list(stats = stats, common = cmn, drop = setdiff(uni, cmn))
}

#' @rdname check_gene_stats_in_multi_batch
#' @export
check_gene_avg_in_multi_batch <- function(srt, genes, batch_col = "Batch", slot = "counts") {
  m <- AverageExpression(srt, features = genes, group.by = batch_col, slot = slot, verbose = FALSE)[[DefaultAssay(srt)]]
  miss <- setdiff(genes, rownames(m))
  if (length(miss)) m <- rbind(m, matrix(0, length(miss), ncol(m), dimnames = list(miss, colnames(m))))
  tibble::as_tibble(m, rownames = "gene")
}

#' @rdname check_gene_stats_in_multi_batch
#' @export
subset_common_gene_in_multi_batch <- function(srt, common_genes, assay = "RNA") {
  DefaultAssay(srt) <- assay
  subset(srt, features = intersect(common_genes, rownames(srt)))
}
