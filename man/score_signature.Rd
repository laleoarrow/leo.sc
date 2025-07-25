% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotation.R
\name{score_signature}
\alias{score_signature}
\title{Score sc_obj with signature list}
\usage{
score_signature(sc_obj, signature_list, seed = 1)
}
\arguments{
\item{sc_obj}{A Seurat object.}

\item{signature_list}{list. A list of gene sets, where each element is a character vector of gene names.}

\item{seed}{An integer seed for reproducibility. Default is 1.}
}
\value{
A Seurat object with module scores added to the metadata.
}
\description{
This function calculates module scores for a Seurat object based on a list of gene sets.
}
\examples{
\dontrun{
sc_obj <- score_signature(sc_obj, signature_list))
}
}
