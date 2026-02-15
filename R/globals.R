#' @importFrom utils globalVariables
utils::globalVariables(c(
  # sko.R
  "cell_type", "sko_gene", "avg_log2FC", "cell", "cluster", "expr", "total", "max_keep", "n", "pct",
  "ident.1", "ident.2", "p_val_adj",

  # basic.R / plot_qc
  "nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB", "doublet.class", "seurat_clusters",

  # sc.plot.R
  "group_by", ".data", "Percentage",

  # celltype.R
  "cell_anno_color",

  # other possible reporting
  "Var1", "Freq",
  
  # added to fix check
  "Dim1", "Dim2",
  ":=", ".",
  "Cell Type", "cell_count", "CellID", "colData", "ct",
  "effect", "gene", "nhoods",
  "pct.1", "pct.2", "s",
  "sc_matrix", "n_rm", "n_tot"
))
