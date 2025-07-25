% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotation.R
\name{calcROGUE}
\alias{calcROGUE}
\title{Calculate and plot ROGUE index for a Seurat object, with inline filtering}
\usage{
calcROGUE(
  obj,
  assay = "RNA",
  layer = "counts",
  downsample = 3000,
  min.cells = 10,
  min.genes = 10,
  anno_label = "celltype_broad",
  sample_label = "orig.ident",
  ...
)
}
\arguments{
\item{obj}{A Seurat object with \code{celltype_broad} and \code{orig.ident} metadata.}

\item{assay}{Assay name containing raw counts (default "RNA").}

\item{layer}{Layer name for raw counts (default "counts").}

\item{downsample}{Number of cells to downsample per group (default 3000).}

\item{min.cells}{Minimum cells per gene for filtering (default 10).}

\item{min.genes}{Minimum genes per cell for filtering (default 10).}

\item{anno_label}{Cell type annotation label (default "celltype_broad").}

\item{sample_label}{Sample label for grouping (default "orig.ident").}

\item{...}{Additional arguments passed to \code{rogue()}.or \code{mc_rogue()}.}
}
\value{
An object returned by \code{rogue()}, invisibly.
}
\description{
Filters raw counts by minimum cells and genes, computes entropy and ROGUE metrics.
}
