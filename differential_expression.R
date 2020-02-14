# this file contains functions for performing differential expression on
# Seurat objects

library(edgeR)

# this is based on https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLFDetRate.R
RunEdgeR <- function(
    seurat,
    cluster,
    grouping.var
) {
    this.cluster.only <- subset(seurat, seurat_clusters == cluster)
    counts <- this.cluster.only@assays$RNA@counts
    group <- this.cluster.only[[grouping.var]]

    dge <- DGEList(counts, group = t(group))
    # perform TMM normalization
    dge <- calcNormFactors(dge)
    # calculate cellular detection rate to use as covariate
    cdr <- scale(colMeans(as.matrix(counts) > 0))
    design <- model.matrix(~ cdr + group)
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    qlf <- glmQLFTest(fit)
    tt <- topTags(qlf, n = Inf)

    return(tt$table)
}
