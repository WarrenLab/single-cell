#!/usr/bin/env Rscript

library(optparse)
library(Seurat)
library(dplyr)
library(ggplot2)

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
max.nfeature <- 2000
max.percent.mt <- 5
nfeatures <- 2000
integrate <- TRUE
group.var <- "challenge_status"
output.dir <- "."
num.pcs <- 20

# load Cell Ranger output
seurat.data <- Read10X(data.dir = feature.bc.matrix)
seurat <- CreateSeuratObject(seurat.data, min.cells = 3, min.features = 200)

# add metadata to Seurat object
meta.data <- read.csv(aggregation)
seurat <- add.meta.data(seurat, meta.data)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

# make a diagnostic plot to help with filtering
p <- ggplot(seurat@meta.data, aes(nCount_RNA, nFeature_RNA, color=percent.mt))
p <- p + geom_point(size=0.2)
p <- p + geom_hline(yintercept=min.nfeature, color='red')
p <- p + geom_hline(yintercept=max.nfeature, color='red')
ggsave(paste(output.dir, 'feature_plot.pdf', sep='/'), plot=p)

# do some basic filtering
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
        x <- FindVariableFeatures(x, nfeatures=nfeatures)
    })
    seurat.anchors <- FindIntegrationAnchors(seurat.list)
    seurat <- IntegrateData(seurat.anchors)
} else {
    seurat <- NormalizeData(seurat)
    seurat <- FindVariableFeatures(seurat, nfeatures = nfeatures)
}

# do several super-standard analysis steps
seurat <- FindClusters(FindNeighbors(RunUMAP(RunPCA(
    ScaleData(seurat)), dims=1:num.pcs)))

# make an elbow plot to help choose number of PCs
p <- ElbowPlot(seurat) + geom_vline(xintercept = num.pcs)
ggsave(paste(output.dir, 'elbow.pdf', sep='/'), plot=p)

# make some UMAP plots
p <- DimPlot(seurat, group.by="library_id")
ggsave(paste(output.dir, 'umap.batches.pdf', sep='/'), plot=p)
if (group.var == '') {
    p <- DimPlot(seurat, label=TRUE) + NoLegend()
} else {
    p <- DimPlot(seurat, group.by=group.var)
    ggsave(paste(output.dir, 'umap.groups.pdf', sep='/'), plot=p)
    p <- DimPlot(seurat, split.by=group.var, label=TRUE) + NoLegend()
}
ggsave(paste(output.dir, 'umap.clusters.pdf', sep='/'), plot=p)

# TODO find conserved markers
# TODO make a heatmap

