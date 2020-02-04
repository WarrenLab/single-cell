#!/usr/bin/env Rscript

library(optparse)
library(Seurat)
library(dplyr)

add.meta.data <- function(seurat, meta.data) {
    gem.group <- as.numeric(sapply(strsplit(rownames(
        spleen@meta.data), split="-"), "[[", 2))
    df <- data.frame(
        lapply(
            function(c) sapply(gem.group, function(g) meta.data[g, c]),
            colnames(meta.data)
        ),
        row.names = rownames(seurat@meta.data)
    )
    names(df) <- colnames(meta.data)
    return AddMetaData(seurat, df)
}

# this does what FindAllMarkers does, except it uses FindConservedMarkers
# inside instead of FindMarkers.
FindAllConservedMarkers <- function(
    seurat,
    grouping.var,
    test.use = 'wilcox',
    ) {
    Reduce(function(df, cluster) {
        df2 <- FindConservedMarkers(
            seurat,
            ident.1 = cluster,
            grouping.var = grouping.var,
            test.use = test.use,
        )
        df2$cluster <- cluster
        rbind(df, df2)
    }, levels(seurat@active.ident), data.frame())
}

options <- commandArgs(trailingOnly = TRUE)
feature.bc.matrix <- options[1]
aggregation <- options[2]
# TODO add these as optparse options
min.nfeature <- 200
max.nfeature <- 5000
max.percent.mt <- 5
nfeatures <- 2000
integrate <- TRUE
group.var <- "challenge_status"

# load Cell Ranger output
seurat.data <- Read10X(data.dir = feature.bc.matrix)
seurat <- CreateSeuratObject(seurat.data, min.cells = 3, min.features = 200)

# add metadata to Seurat object
meta.data <- read.csv(aggregation)
seurat <- add.meta.data(seurat, meta.data)

# TODO make some diagnostic plots to help with filtering

# do some basic filtering
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat <- subset(
    seurat,
    subset = nFeature_RNA > min.nfeature & nFeature_RNA < max.nFeature &
        percent.mt < max.percent.mt
)

# to integrate or not to integrate?
if (integrate) {
    seurat.list <- SplitObject(seurat, split.by = group.var)
    seurat.list <- lapply(seurat.list, function (x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method="vst",
                                  nfeatures=nfeatures)
    })
    seurat.anchors <- FindIntegrationAnchors(seurat.list)
    seurat <- IntegrateData(seurat.anchors)
} else {
    seurat <- NormalizeData(seurat)
    seurat <- VindVariableFeatures(
        seurat, selection.method = "vst", nfeatures = nfeatures)
}

# do several super-standard analysis steps
seurat <- FindClusters(FindNeighbors(RunUMAP(RunPCA(
    ScaleData(seurat)), dims=1:30), dims=1:30))

# TODO make some UMAP plots
# TODO find conserved markers
# TODO make a heatmap
