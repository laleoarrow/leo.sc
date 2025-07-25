% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic.R
\name{seurat_basic_info}
\alias{seurat_basic_info}
\title{get basic summary of a seurat object}
\usage{
seurat_basic_info(seurat_obj, assay = NULL, slot = "counts")
}
\arguments{
\item{seurat_obj}{Seurat object}

\item{assay}{character. assay name, defaults to active assay}

\item{slot}{character. data slot ("counts", "data" or "scale.data")}
}
\value{
named list with total_genes, total_cells, avg_genes_per_cell, genes_per_cell
}
\description{
get basic summary of a seurat object
}
