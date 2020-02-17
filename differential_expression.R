# this file contains functions for performing differential expression on
# Seurat objects

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

    library(edgeR)
    dge <- DGEList(counts, group = t(group))
    # perform TMM normalization
    dge <- calcNormFactors(dge)
    # calculate cellular detection rate to use as covariate
    cdr <- calc.cdr(counts)
    design <- model.matrix(~ cdr + group)
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    qlf <- glmQLFTest(fit)
    tt <- topTags(qlf, n = Inf)

    return(tt$table)
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

    library(MAST)
    library(edgeR)
    cdr <- calc.cdr(counts)
    dge <- DGEList(counts)
    dge <- edgeR::calcNormFactors(dge)
    cpms <- edgeR::cpm(dge)
    sca <- FromMatrix(exprsArray = log2(cpms + 1),
                      cData = data.frame(group = group, cdr = cdr))
    zlmdata <- zlm(~cdr + group, sca)
    mast <- lrTest(zlmdata, "group")

    return(mast[, 'hurdle',])
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

    # TODO calculate tpm! this can be done using Seurat's NormalizeData, but
    # need to make sure we get the linearly scaled data out of the object
    tmm <- edgeR::calcNormFactors(tpm)
    tpm.tmm <- edgeR::cpm(tpm, lib.size = tmm * colSums(tpm))
    idx <- 1:nrow(tpm.tmm)
    names(idx) <- rownames(tpm.tmm)
    wilcox_p <- sapply(idx, function(i) {
        wilcox.test(tpm.tmm[i, ] ~ group)$p.value
    })
    # TODO calculate fold changes too, unless wilcox.test actually returns this
}
