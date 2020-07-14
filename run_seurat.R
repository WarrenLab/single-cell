#!/usr/bin/env Rscript

library(argparser)

ParseArguments <- function() {
    p <- arg_parser('Run Seurat')
    p <- add_argument(p, '--min-nfeature', default=200,
                      help='minimum number of features to use a cell')
    p <- add_argument(p, '--max-nfeature', default=5000,
                      help='maximum number of features to use a cell')
    p <- add_argument(p, '--max-percent-mt', default=20.0,
                      help='max % of counts from MT genes to use a cell')
    p <- add_argument(p, '--max-percent-ribo', default=100.0,
                      help='max % of counts from ribosomal genes to use a cell')
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
                      help='folder or h5 containing feature matrix')
    p <- add_argument(p, '--aggregation',
                      help='csv containing metadata, output by Cell Ranger')
    return(parse_args(p))
}

argv <- ParseArguments()

library(Seurat)
library(dplyr)
library(ggplot2)
library(warrenlabSC)

# create output directory
dir.create(argv$output_dir)

# load Cell Ranger output
if (endsWith(argv$feature_matrix, '.h5')) {
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
        stop("Please install package \"hdf5r\" to use 10x h5 as input.",
             call. = FALSE)
    }
    seurat.data <- Read10X_h5(argv$feature_matrix)
} else {
    seurat.data <- Read10X(data.dir = argv$feature_matrix)
}
seurat <- CreateSeuratObject(seurat.data, min.cells = 3, min.features = 100)

# add metadata to Seurat object
if (!is.na(argv$aggregation)) {
    meta.data <- read.csv(argv$aggregation)
    seurat <- add.meta.data(seurat, meta.data)
}
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^(MT|Mt|mt)-")
seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^(RP[LS]|rp[ls])")

# make a diagnostic plot to help with filtering
p <- ggplot(seurat@meta.data, aes(nCount_RNA, nFeature_RNA, color=percent.mt))
p <- p + geom_point(size=0.2)
p <- p + geom_hline(yintercept=argv$min_nfeature, color='red')
p <- p + geom_hline(yintercept=argv$max_nfeature, color='red')
ggsave(file.path(argv$output_dir, 'feature_plot.pdf'), plot=p)

if (!is.na(argv$aggregation)) {
    p <- VlnPlot(
        seurat,
        features = c('percent.mt', 'percent.ribo', 'nFeature_RNA', 'nCount_RNA'),
        group.by = 'library_id', ncol = 1, pt.size = 0
    )
    ggsave(file.path(argv$output_dir, 'violin_plot.pdf'), plot = p,
           height = 21, width = 7)
} else {
    p <- VlnPlot(
        seurat,
        features = c('percent.mt', 'percent.ribo', 'nFeature_RNA', 'nCount_RNA'),
        ncol = 2, pt.size = 0
    )
    ggsave(file.path(argv$output_dir, 'violin_plot.pdf'), plot = p)
}

# do some basic filtering
seurat <- subset(
    seurat,
    subset = (nFeature_RNA > argv$min_nfeature
              & nFeature_RNA < argv$max_nfeature
              & percent.mt < argv$max_percent_mt
              & percent.ribo < argv$max_percent_ribo)
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
                variable.features.n = argv$nfeatures)
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
            FindVariableFeatures(NormalizeData(x), nfeatures=argv$nfeatures)
        })
        seurat.anchors <- FindIntegrationAnchors(seurat.list)
        seurat <- IntegrateData(seurat.anchors)
        DefaultAssay(seurat) <- "integrated"
        seurat <- ScaleData(seurat)
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

# save seurat object for easy loading later
saveRDS(seurat, file = file.path(argv$output_dir, 'seurat.rds'))

# make an elbow plot to help choose number of PCs
p <- ElbowPlot(seurat, ndims = argv$num_pcs + 10)
p <- p + geom_vline(xintercept = argv$num_pcs)
ggsave(file.path(argv$output_dir, 'elbow.pdf'), plot=p)

# make some UMAP plots
if (!is.na(argv$aggregation)) {
    p <- DimPlot(seurat, group.by="library_id")
    ggsave(file.path(argv$output_dir, 'umap.batches.pdf'), plot=p)
}
if (is.na(argv$group_var)) {
    p <- DimPlot(seurat, label=TRUE) + NoLegend()
    ggsave(file.path(argv$output_dir, 'umap.clusters.pdf'), plot=p)
} else {
    p <- DimPlot(seurat, group.by=argv$group_var)
    ggsave(file.path(argv$output_dir, 'umap.groups.pdf'), plot=p)
    p <- DimPlot(seurat, split.by=argv$group_var, label=TRUE) + NoLegend()
    ggsave(file.path(argv$output_dir, 'umap.clusters.pdf'), width = 14, plot=p)
}

# find biomarkers for each cluster
DefaultAssay(seurat) <- 'RNA' # always do DE analysis on raw counts
if (argv$integrate) {
    all.markers <- FindAllConservedMarkers(seurat, grouping.var=argv$group_var)
    top5 <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = -max_pval)
} else {
    all.markers <- FindAllMarkers(seurat)
    top5 <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
}

write.csv(all.markers, file = file.path(argv$output_dir, "all_markers.csv"),
          quote = FALSE)
write.csv(top5, file = file.path(argv$output_dir, "top5.csv"), quote = FALSE)

# make a heatmap
if (argv$integrate) {
    # if we did integrated normalization, there might be missing genes in the
    # combined set, so we need to have a plain normalized slot to make a heatmap
    seurat <- ScaleData(NormalizeData(seurat))
    p <- DoHeatmap(seurat, features = top5$feature) + NoLegend()
} else {
    p <- DoHeatmap(seurat, features = top5$gene) + NoLegend()
}
ggsave(file.path(argv$output_dir, 'markers_heatmap.pdf'),
       plot=p, width=10, height=20)

if (!is.na(argv$group_var)) {
    # differential expression per cell type between groups
    per.cluster.DE.edgeR <- FindAllClusterDE(
        seurat,
        argv$group_var,
        method=RunEdgeR
    )
    write.csv(per.cluster.DE.edgeR,
              file.path(argv$output_dir, 'per_cluster_DE.edgeR.csv'),
              quote = FALSE)
}

