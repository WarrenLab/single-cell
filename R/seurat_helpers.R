#' Add metadata to a Seurat object from a data frame
#'
#' When creating a Seurat object with, for example, \code{Read10X}, no metadata
#' is loaded automatically, even though \code{cellranger aggregate} gives you
#' a nice aggregation csv. The cell barcodes just contain a numerical suffix to
#' indicate which library they're from. It's useful to have all the metadata
#' associated with each library as part of the Seurat object, so this function
#' parses each cell barcode to figure out what gem group the cell came from,
#' and then applies gem-group-level metadata contained in a table to the cell.
#'
#' @param seurat A seurat object.
#' @param meta.data A table containing gem-group-level metadata to be added to
#'   \code{seurat@meta.data}. Rows must be in order of gem group, as output by
#'   \code{cellranger aggregate}.
#' @return The seurat object with metadata added.
#' @examples
#' meta.data <- read.csv('aggregation.csv')
#' seurat <- add.meta.data(seurat, meta.data)
add.meta.data <- function(seurat, meta.data) {
    gem.group <- as.numeric(sapply(strsplit(rownames(
        seurat@meta.data), split="-"), "[[", 2))
    df <- data.frame(
        lapply(
            colnames(meta.data),
            function(c) sapply(gem.group, function(g) meta.data[g, c])
        ),
        row.names = rownames(seurat@meta.data)
    )
    names(df) <- colnames(meta.data)
    return(AddMetaData(seurat, df))
}

#' Find conserved markers for each cluster
#'
#' Find conserved markers for each cluster. Basically, this does what
#' \link{\code{Seurat::FindAllMarkers}} does, except it uses
#' \link{\code{Seurat::FindConservedMarkers}} inside instead of plain old
#' \link{\code{Seurat::FindMarkers}}.
#'
#' @param seurat A seurat object.
#' @param grouping.var The name of the metadata variable to use to separate the
#'   cells into two conditions.
#' @param test.use String specifying the test to pass to \code{FindMarkers}.
#' @return A \code{data.frame} of the outputs of all calls to
#'   \code{FindConservedMarkers}, row-bound together with the added columns
#'   \code{cluster} and \code{feature}.
FindAllConservedMarkers <- function(
    seurat,
    grouping.var,
    test.use = 'wilcox'
    ) {
    Reduce(function(df, cluster) {
        df2 <- FindConservedMarkers(
            seurat,
            ident.1 = cluster,
            grouping.var = grouping.var,
            only.pos = TRUE,
            test.use = test.use
        )
        df2$cluster <- cluster
        df2$feature <- rownames(df2)
        rbind(df, df2, make.row.names = FALSE)
    }, levels(seurat$seurat_clusters), data.frame())
}
