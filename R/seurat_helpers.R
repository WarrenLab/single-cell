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

# this does what FindAllMarkers does, except it uses FindConservedMarkers
# inside instead of FindMarkers.
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
