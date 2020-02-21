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

library(Seurat)
library(dplyr)
library(ggplot2)
library(edgeR)
library(MAST)

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
    }, levels(seurat$seurat_clusters), data.frame())
}

calc.cdr <- function(counts) scale(colMeans(as.matrix(counts) > 0))

# this is based on https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLFDetRate.R
RunEdgeR <- function(
    seurat,
    cluster,
    grouping.var
) {
    this.cluster.only <- subset(seurat, seurat_clusters == cluster)
    counts <- this.cluster.only@assays$RNA@counts
    group <- as.matrix(this.cluster.only[[grouping.var]])

    dge <- DGEList(counts, group = t(group))
    # perform TMM normalization
    dge <- edgeR::calcNormFactors(dge)
    # calculate cellular detection rate to use as covariate
    cdr <- calc.cdr(counts)
    design <- model.matrix(~ cdr + group)
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    qlf <- glmQLFTest(fit)
    tt <- topTags(qlf, n = Inf)

    return(data.frame(p = tt$table$PValue,
                      fold.change = tt$table$logFC,
                      row.names = rownames(tt$table)))
}

# based on https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_MASTcpmDetRate.R
RunMAST <- function(
    seurat,
    cluster,
    grouping.var
) {
    this.cluster.only <- subset(seurat, seurat_clusters == cluster)
    counts <- this.cluster.only@assays$RNA@counts
    group <- as.matrix(this.cluster.only[[grouping.var]])
    group.names <- levels(factor(group[,1]))

    cdr <- calc.cdr(counts)
    dge <- DGEList(counts)
    dge <- edgeR::calcNormFactors(dge)
    cpms <- edgeR::cpm(dge)
    sca <- FromMatrix(exprsArray = log2(cpms + 1),
                      cData = data.frame(group = group, cdr = cdr))
    zlmdata <- zlm(~cdr + group, sca)
    mast <- lrTest(zlmdata, "group")

    idx <- 1:nrow(mast[,'hurdle',])
    names(idx) <- rownames(mast[,'hurdle',])
    fold.change <- sapply(
        idx, function(i) log2(mean(cpms[i, group == group.names[1]]) + 1)
                       - log2(mean(cpms[i, group == group.names[2]]) + 1)
        )


    return(data.frame(p = mast[, 'hurdle', "Pr(>Chisq)"],
                      fold.change = fold.change,
                      row.names = rownames(mast[, 'hurdle',])))
}

# based on https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_Wilcoxon.R
RunWilcoxon <- function(
    seurat,
    cluster,
    grouping.var
) {
    this.cluster.only <- subset(seurat, seurat_clusters == cluster)
    counts <- this.cluster.only@assays$RNA@counts
    group <- this.cluster.only[[grouping.var]]
    group.names <- levels(factor(group[,1]))

    dge <- DGEList(counts)
    dge <- edgeR::calcNormFactors(dge)
    cpms <- edgeR::cpm(dge)
    idx <- 1:nrow(cpms)
    names(idx) <- rownames(cpms)
    wilcox.p <- sapply(
        idx, function(i) wilcox.test(cpms[i, ] ~ group)$p.value)
    fold.change <- sapply(
        idx, function(i) log2(mean(cpms[i, group == group.names[1]]) + 1)
                       - log2(mean(cpms[i, group == group.names[2]]) + 1)
        )
    return(data.frame(p = wilcox.p,
                      fold.change = fold.change,
                      row.names = names(wilcox.p)))
}

FindAllClusterDE <- function(
    seurat,
    grouping.var,
    method = RunEdgeR
) {
    Reduce(function(df, cluster) {
        df2 <- method(
            seurat,
            cluster = cluster,
            grouping.var = grouping.var
        )
        df2$cluster <- cluster
        rbind(df, df2)
    }, levels(seurat$seurat_clusters), data.frame())
}

# create output directory
dir.create(argv$output_dir)

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
ggsave(file.path(argv$output_dir, 'feature_plot.pdf'), plot=p)

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
p <- DimPlot(seurat, group.by="library_id")
ggsave(file.path(argv$output_dir, 'umap.batches.pdf'), plot=p)
if (is.na(argv$group_var)) {
    p <- DimPlot(seurat, label=TRUE) + NoLegend()
} else {
    p <- DimPlot(seurat, group.by=argv$group_var)
    ggsave(file.path(argv$output_dir, 'umap.groups.pdf'),
           width=14, plot=p)
    p <- DimPlot(seurat, split.by=argv$group_var, label=TRUE) + NoLegend()
}
ggsave(file.path(argv$output_dir, 'umap.clusters.pdf'), plot=p)

# find biomarkers for each cluster
DefaultAssay(seurat) <- 'RNA' # always do DE analysis on raw counts
if (argv$integrate) {
    all.markers <- FindAllConservedMarkers(seurat, grouping.var=argv$group_var)
} else {
    all.markers <- FindAllMarkers(seurat)
}

all.markers$feature <- rownames(all.markers)
top5 <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = -max_pval)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = -max_pval)

# make a heatmap
if (argv$integrate) DefaultAssay(seurat) <- 'integrated'
p <- DoHeatmap(seurat, features=top5$feature) + NoLegend()
ggsave(file.path(argv$output_dir, 'markers_heatmap.pdf'),
       plot=p, width=10, height=20)

# differential expression per cell type between groups
per.cluster.DE.edgeR <- FindAllClusterDE(seurat,
                                         argv$group_var, method=RunEdgeR)
write.csv(per.cluster.DE.edgeR,
          file.path(argv$output_dir, 'per_cluster_DE.edgeR.csv'))
per.cluster.DE.MAST <- FindAllClusterDE(seurat, argv$group_var, method=RunMAST)
write.csv(per.cluster.DE.MAST,
          file.path(argv$output_dir, 'per_cluster_DE.MAST.csv'))
per.cluster.DE.wilcoxon <- FindAllClusterDE(seurat, argv$group_var,
                                            method=RunWilcoxon)
write.csv(per.cluster.DE.wilcoxon,
          file.path(argv$output_dir, 'per_cluster_DE.wilcox.csv'))

