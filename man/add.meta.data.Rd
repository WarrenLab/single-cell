% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seurat_helpers.R
\name{add.meta.data}
\alias{add.meta.data}
\title{Add metadata to a Seurat object from a data frame}
\usage{
add.meta.data(seurat, meta.data)
}
\arguments{
\item{seurat}{A seurat object.}

\item{meta.data}{A table containing gem-group-level metadata to be added to
\code{seurat@meta.data}. Rows must be in order of gem group, as output by
\code{cellranger aggregate}.}
}
\value{
The seurat object with metadata added.
}
\description{
When creating a Seurat object with, for example, \code{Read10X}, no metadata
is loaded automatically, even though \code{cellranger aggregate} gives you
a nice aggregation csv. The cell barcodes just contain a numerical suffix to
indicate which library they're from. It's useful to have all the metadata
associated with each library as part of the Seurat object, so this function
parses each cell barcode to figure out what gem group the cell came from,
and then applies gem-group-level metadata contained in a table to the cell.
}
\examples{
meta.data <- read.csv('aggregation.csv')
seurat <- add.meta.data(seurat, meta.data)
}
