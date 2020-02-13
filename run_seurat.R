#!/usr/bin/env Rscript

library(argparser)
library(Seurat)
library(dplyr)
library(ggplot2)

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
            test.use = test.use
        )
        df2$cluster <- cluster
        rbind(df, df2)
    }, levels(seurat@active.ident), data.frame())
}

ParseArguments <- function() {
    p <- arg_parser('Run Seurat')
    p <- add_argument(p, '--min-nfeature', default=200,
                      help='minimum number of features to use a cell')
    p <- add_argument(p, '--max-nfeature', default=5000,
                      help='maximum number of features to use a cell')
    p <- add_argument(p, '--max-percent-mt', default=6.0,
                      help='max % of counts from MT genes to use a cell')
    p <- add_argument(p, '--nfeatures', default=2000,
                      help='use top N features [2000]')
    p <- add_argument(p, '--integrate', flag=TRUE,
                      help='integrate dataset divided by --group-var')
    p <- add_argument(p, '--group-var', help='metadata variable to group by')
    p <- add_argument(p, '--nonlinear', flag=TRUE, help='use SCTransform')
    p <- add_argument(p, '--output-dir', default='.',
                      help='output directory for plots and tables')
    p <- add_argument(p, '--num-pcs', default=20,
                      help="number of principal components to use")
    p <- add_argument(p, 'feature-matrix',
                      help='folder containing feature matrix')
    p <- add_argument(p, 'aggregation',
                      help='csv containing metadata, output by Cell Ranger')
    return(parse_args(p))
}

argv <- ParseArguments()

# load Cell Ranger output
seurat.data <- Read10X(data.dir = argv$feature_matrix)
seurat <- CreateSeuratObject(seurat.data, min.cells = 3, min.features = 100)

# add metadata to Seurat object
meta.data <- read.csv(argv$aggregation)
seurat <- add.meta.data(seurat, meta.data)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^RP[LS]")

# make a diagnostic plot to help with filtering
p <- ggplot(seurat@meta.data, aes(nCount_RNA, nFeature_RNA, color=percent.mt))
p <- p + geom_point(size=0.2)
p <- p + geom_hline(yintercept=argv$min_nfeature, color='red')
p <- p + geom_hline(yintercept=argv$max_nfeature, color='red')
ggsave(paste(argv$output_dir, 'feature_plot.pdf', sep='/'), plot=p)

# do some basic filtering
seurat <- subset(
    seurat,
    subset = (nFeature_RNA > argv$min_nfeature
              & nFeature_RNA < argv$max_nfeature
              & percent.mt < argv$max_percent_mt)
)

# normalize, find variable features, and scale. Integrate datasets if
# requested, and use non-linear normalization if requested.
if (argv$integrate) {
    seurat.list <- SplitObject(seurat, split.by = argv$group_var)
    if (argv$nonlinear) {
        seurat.list <- lapply(seurat.list, function (x) {
            SCTransform(
                x, vars.to.regress = c("percent.mt", "percent.ribo"),
                do.scale = TRUE,
                variable.features.n = nfeatures)
        })
        seurat.features <- SelectIntegrationFeatures(seurat.list,
                                                     nfeatures=argv$nfeatures)
        options(future.globals.maxSize = 891289600)
        seurat.list <- PrepSCTIntegration(seurat.list,
                                          anchor.features = seurat.features)
        seurat.anchors <- FindIntegrationAnchors(
            seurat.list, normalization.method = "SCT",
            anchor.features = seurat.features)
        seurat <- IntegrateData(seurat.anchors, normalization.method = "SCT")
    } else {
        seurat.list <- lapply(seurat.list, function (x) {
            FindVariableFeatures(NormalizatData(x), nfeatures=argv$nfeatures)
        })
        seurat.anchors <- FindIntegrationAnchors(seurat.list)
        seurat <- IntegrateData(seurat.anchors)
    }
} else {
    if (argv$nonlinear) {
        seurat <- SCTransform(
            seurat,
            vars.to.regress = c("percent.mt", "percent.ribo"),
            do.scale = TRUE)
    } else {
        seurat <- NormalizeData(seurat)
        seurat <- FindVariableFeatures(seurat, nfeatures = argv$nfeatures)
        seurat <- ScaleData(seurat)
    }
}

# do several super-standard analysis steps
seurat <- FindClusters(FindNeighbors(RunUMAP(RunPCA(seurat),
                                             dims=1:argv$num_pcs)))

# make an elbow plot to help choose number of PCs
p <- ElbowPlot(seurat, ndims = argv$num_pcs + 10)
p <- p + geom_vline(xintercept = argv$num_pcs)
ggsave(paste(argv$output_dir, 'elbow.pdf', sep='/'), plot=p)

# make some UMAP plots
p <- DimPlot(seurat, group.by="library_id")
ggsave(paste(argv$output_dir, 'umap.batches.pdf', sep='/'), plot=p)
if (is.na(argv$group_var)) {
    p <- DimPlot(seurat, label=TRUE) + NoLegend()
} else {
    # TODO increase width of this plot
    p <- DimPlot(seurat, group.by=argv$group_var)
    ggsave(paste(argv$output_dir, 'umap.groups.pdf', sep='/'), plot=p)
    p <- DimPlot(seurat, split.by=argv$group_var, label=TRUE) + NoLegend()
}
ggsave(paste(argv$output_dir, 'umap.clusters.pdf', sep='/'), plot=p)

# find biomarkers for each cluster
DefaultAssay(seurat) <- 'RNA' # always do DE analysis on raw counts
if (argv$integrated)
    all.markers <- FindAllConservedMarkers(seurat, grouping.var=argv_group_var)
else
    all.markers <- FindAllMarkers(seurat)

all.markers$feature <- rownames(all.markers)
top5 <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = -max_pval)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = -max_pval)

# make a heatmap # TODO change dimensions
if (argv$integrated) DefaultAssay(seurat) <- 'integrated'
p <- DoHeatmap(spleen, features=top5$feature) + NoLegend()
ggsave(paste(argv$output_dir, 'markers_heatmap.pdf', sep=','), plot=p)

# TODO differential expression per cell type between groups

