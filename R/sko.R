#' In-silico Knock-Out / Knock-In Analysis (single gene)
#'
#' Performs a virtual knock-out or knock-in simulation for a target gene
#' in a single-cell dataset. The function identifies cells with high or
#' low expression of the gene, matches group sizes and cell-type composition,
#' and performs differential expression and enrichment analysis to infer
#' potential functional effects.
#'
#' @param all A Seurat object containing expression data and cell metadata.
#' @param gene Character(1). The target gene symbol (must exist in the expression matrix).
#' @param sko_mode One of \code{c("ko", "ki")}. Specifies virtual knock-out ("ko") or knock-in ("ki") mode.
#' @param cell_col Character(1). The metadata column indicating cell type.
#' @param filter_cell_threshold Integer. Minimum number of cells per cell type to retain.
#' @param pct_threshold Numeric. Fraction of expressing cells to extract per group
#'   (ignored if \code{abs_threshold} is set).
#' @param abs_threshold Integer. Absolute number of cells per group to extract;
#'   overrides \code{pct_threshold} if provided.
#' @param deg_method Character. Differential expression method: \code{"default"} or \code{"MAST"}.
#' @param enrichment_method Character vector. Enrichment methods to run; default is \code{c("ORA", "GSEA")}.
#' @param enrichment_bg Character vector. Background databases for enrichment;
#'   default is \code{c("GO", "KEGG", "MKEGG", "Reactome")}.
#'
#' @details
#' The pipeline includes:
#' \enumerate{
#'   \item Filtering cells expressing the target gene.
#'   \item Retaining cell types with sufficient cell counts.
#'   \item Selecting top-expressing cells (high group) and matched low-expressing cells (low group).
#'   \item Running differential expression between the two groups using \code{Seurat::FindMarkers}.
#'   \item Performing enrichment analysis via \code{leo.basic::leo_enrich}.
#' }
#'
#' @return Invisibly returns a list containing:
#' \itemize{
#'   \item \code{high_cells}: Vector of selected high-expression cells.
#'   \item \code{low_cells}: Vector of selected low-expression cells.
#'   \item \code{cell_rato}: Data frame summarizing cell-type composition.
#'   \item \code{deg_results}: Differential expression results (data frame).
#'   \item \code{enrichments}: Enrichment analysis results.
#' }
#'
#' @seealso \code{\link[Seurat]{FindMarkers}}, \code{\link[leo.basic]{leo_enrich}}
#'
#' @examples
#' \dontrun{
#' res <- silico_ko(
#'   all = srt, gene = "LRRK2",
#'   sko_mode = "ko", cell_col = "cell_anno",
#'   pct_threshold = 0.1, deg_method = "default"
#' )
#' }
#'
#' @export
#' @importFrom Seurat GetAssayData FindMarkers Idents
#' @importFrom tibble tibble rownames_to_column
#' @importFrom dplyr count filter pull mutate arrange desc select slice_head everything case_when n
#' @importFrom rlang sym
#' @importFrom stats setNames
#' @importFrom leo.basic leo_log leo_enrich format_geneList
silico_ko <- function(all, gene, sko_mode = c("ko", "ki"),
                      cell_col = "cell_anno", filter_cell_threshold = 10,
                      pct_threshold = 0.1,
                      abs_threshold = NULL, # if set, will override pct_threshold
                      deg_method = "default",
                      enrichment_method = c("ORA", "GSEA"),
                      enrichment_bg = c("GO", "KEGG", "MKEGG", "Reactome")){
  # Here we develop a new pipeline to estimate the effect of the virtual knock-down for the genes of interests.
  if (!gene %in% rownames(all)) return(leo_log("Gene {gene} not found in the dataset.", level = "danger"))
  if (length(gene) > 1) return(leo_log("Only one gene can be processed at a time.", level = "danger"))

  leo_log <- leo.basic::leo_log; start_time <- Sys.time()
  sko_mode <- match.arg(sko_mode, c("ko", "ki"))
  if (sko_mode == "ko") leo_log("Performing in-silico knock-down for: {gene}")
  if (sko_mode == "ki") leo_log("Performing in-silico knock-in for: {gene}")

  # Step 1: Pre-filter -------
  # Filter out cells with zero expression for the target gene
  leo_log("Step 1: Filter out cells with zero expression for the target gene")
  expr_vals <- Seurat::GetAssayData(all, slot = "data")[gene, ]
  cells_express <- names(expr_vals)[expr_vals > 0]
  n_total <- length(cells_express)

  if (n_total == 0)  return(leo_log("No cells express {gene}, aborting.", level = "danger"))
  if (n_total < 100) leo_log("Only {n_total} cell{?s} express {gene}, results may be unstable.", level = "warning")
  all2 <- subset(all, cells = cells_express); leo_log("Found {n_total} cell{?s} expressing {gene}.", level = "info")

  # Filter out also the cell type with only too few number of cells
  cell_counts_in_type <- all2@meta.data %>%
    dplyr::select(!!rlang::sym(cell_col)) %>%
    dplyr::count(!!rlang::sym(cell_col)) %>%
    dplyr::filter(n > filter_cell_threshold) %>%
    dplyr::pull(!!rlang::sym(cell_col)) %>% as.vector()
  leo_log("Found {length(cell_counts_in_type)} cell type{?s} with more than {filter_cell_threshold} cells --> applying filtering")
  all2 <- subset_srt(all2, subset_col = cell_col, subset_value = cell_counts_in_type)

  # Step 2: Determine number of high-expression cells to extract ------
  n_extract <- if (!is.null(abs_threshold)) abs_threshold else floor(n_total * pct_threshold)
  if (n_extract < 100) {n_extract <- 100; leo_log("n_extract < 100, set to 100 to ensure minimal power", level = "warning")}  # ensure at least 100 cell
  leo_log("Step 2: extract {.emph {n_extract}} cell{?s} with high {.emph {gene}} expression for analysis")

  # Step 3: Get top expressing cells ------
  gene_vals <- Seurat::GetAssayData(all2, slot = "data")[gene, ]
  df_gene   <- tibble::tibble(cell = names(gene_vals),
                              gene = gene,
                              expr = gene_vals,
                              cell_type = all2@meta.data[, cell_col]
  ) %>% dplyr::arrange(dplyr::desc(expr))

  # We set a quota for each cell type in case we can not locate bottom cell
  leo_log("Calculating quota for each cell type based on 50% of expressing cells")
  quota_tbl <- df_gene %>%
    dplyr::count(cell_type, name = "total") %>%
    dplyr::mutate(max_keep = floor(total / 2))
  quota_max <- sum(quota_tbl$max_keep)
  if (n_extract > quota_max) {  # check if set n_extract exceeds the quota
    leo_log("Requested {n_extract} exceeds per-type 50% cap (max {quota_max}).",
            "Automatically reset n_extract = {quota_max} (~50% of expressing cells).", level = "warning")
    n_extract <- quota_max
    if (n_extract < 100)
      leo_log("After adjustment, n_extract = {n_extract} < 100. Statistical power limited.", level = "warning")
  } else { leo_log("Quota max ({quota_max}) for each cell type satisfied n_extract ({n_extract})") }

  # extraction
  greedy_quota_select <- function(df, quota_vec, n_target) {
    # Agrs:
    #   @df: data frame with at least columns 'cell', 'expr', 'cell_type', ordered by 'expr' (ps: this func go down from line 1)
    #   @quota_vec: named vector: name is cell type, value is the max number of cells to keep for that cell type
    #   @n_target: total number of cells to select
    # Returns:
    #   A character vector of selected cell names
    select <- character(0); q_left <- quota_vec
    for (i in seq_len(nrow(df))) {
      ct <- as.character(df$cell_type[i])
      if (is.na(q_left[ct])) next
      if (q_left[ct] > 0 && length(select) < n_target) {
        select <- c(select, df$cell[i])
        q_left[ct] <- q_left[ct] - 1
      }
      if (length(select) == n_target) break
    }
    select
  }
  high_cells <- greedy_quota_select(
    df = df_gene,
    quota_vec = stats::setNames(quota_tbl$max_keep, quota_tbl$cell_type),
    n_target = n_extract
  )

  # Echo final high group composition by cell type
  high_cells_composition <- tibble::tibble(cell_type = all2@meta.data[high_cells, cell_col]) %>%
    dplyr::count(cell_type, name = "n") %>%
    dplyr::mutate(pct = round(100 * n / sum(n), 1)) %>%
    dplyr::arrange(dplyr::desc(n))
  leo_log("High group composition (N = {sum(high_cells_composition$n)}): {paste0(high_cells_composition$cell_type, ': ', high_cells_composition$n, ' (', high_cells_composition$pct, '%)', collapse = ', ')}")


  # step 4: Get bottom expressing cells  ------
  # We extract same number of cells while keep the same cell type ratio as high expression group
  leo_log("Locating bottom {.emph {n_extract}} cell{?s} with low expression of {.emph {gene}}")
  df_gene_low <- df_gene %>% dplyr::filter(!cell %in% high_cells) %>% dplyr::arrange(expr)
  quota_low <- stats::setNames(high_cells_composition$n, high_cells_composition$cell_type)

  low_cells <- greedy_quota_select(
    df = df_gene_low,
    quota_vec = quota_low,
    n_target  = n_extract
  )

  low_cells_composition <- df_gene_low %>%
    dplyr::filter(cell %in% low_cells) %>%
    dplyr::count(cell_type, name = "n") %>%
    dplyr::mutate(pct = round(100 * n / sum(n), 1)) %>%
    dplyr::arrange(dplyr::desc(n))

  # defensive coding (in case)
  if (!identical(low_cells_composition, high_cells_composition)) {
    leo_log("The low_cells_composition and high_cells_composition is not the same; check if it is expected!", level = "warning")
  }
  if (length(intersect(high_cells, low_cells)) > 0) {
    leo_log("High and low cells overlap, this is not expected! Check", level = "danger")
  }
  # 若低表达数量不足，同步缩减两组 in case
  # if (length(low_cells) < n_extract) {
  #   leo_log("Low group shortage: only {length(low_cells)} – down-scaling both groups.", level = "warning")
  #   n_extract <- length(low_cells); high_cells <- head(high_cells, n_extract)
  # }

  # step 5: DEG for high/ low in-silico knock-out group -----
  all2@meta.data <- all2@meta.data %>%
    mutate(sko_group = dplyr::case_when(
      rownames(.) %in% high_cells ~ "high",
      rownames(.) %in% low_cells  ~ "low",
      T ~ "other"))
  all2 <- subset_srt(all2, "sko_group", c("high", "low"))
  all2@meta.data$sko_group <- factor(all2@meta.data$sko_group, level = c("high", "low"))
  Seurat::Idents(all2) <- all2@meta.data$sko_group

  # Ensure ident.1 (log2FC>1) stands for the group of interests.
  ident.1 <- case_when(sko_mode == "ko"  ~ "low",
                       sko_mode == "ki"  ~ "high")
  ident.2 <- case_when(sko_mode == "ko"  ~ "high",
                       sko_mode == "ki"  ~ "low")

  deg_results <- switch(
    deg_method,
    "default" = Seurat::FindMarkers(all2, ident.1 = ident.1, ident.2 = ident.2),
    "MAST"    = Seurat::FindMarkers(all2, ident.1 = ident.1, ident.2 = ident.2, test.use = "MAST")
  )
  deg_results <- deg_results %>%
    dplyr::mutate(ident.1 = ident.1, ident.2 = ident.2, sko_gene = gene) %>%
    dplyr::select(sko_gene, ident.1, ident.2, dplyr::everything())

  up  <- sum(deg_results$avg_log2FC > 0 & deg_results$p_val_adj < 0.05)
  down<- sum(deg_results$avg_log2FC < 0 & deg_results$p_val_adj < 0.05)
  leo_log("DEG done –-> more expressed in [{ident.1}]: {up}; more expressed in [{ident.2}]: {down}")

  # step 6: Enrichment analysis ------
  leo_log("Step 6: Enrichment analysis for {.emph {gene}} in {sko_mode} mode")
  ORA_gene <-  deg_results %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05) %>%
    dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
    dplyr::slice_head(n=100) %>%
    dplyr::pull(gene)
  GSEA_geneList <- deg_results %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::select(gene, avg_log2FC) %>%
    dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
    leo.basic::format_geneList()

  enrichment_res <- leo.basic::leo_enrich(ORA_gene, GSEA_geneList,
                                          simplify =T, input = "SYMBOL",
                                          method = enrichment_method,
                                          background = enrichment_bg)

  # return
  leo_time_elapsed(start_time)

  invisible(list(high_cells  = high_cells,
                 low_cells   = low_cells,
                 cell_rato   = high_cells_composition,
                 deg_results = deg_results,
                 enrichments = enrichment_res))
}
