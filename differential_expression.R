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
    group.names <- levels(factor(group[,1]))

    library(edgeR)
    dge <- DGEList(counts)
    dge <- calcNormFactors(dge)
    cpms <- cpm(dge)
    idx <- 1:nrow(cpms)
    names(idx) <- rownames(cpms)
    wilcox.p <- sapply(
        idx, function(i) wilcox.test(cpms[i, ] ~ group)$p.value)
    fold.change <- sapply(
        idx, function(i) log2(mean(cpms[i, group == group.names[1]]) + 1)
                       - log2(mean(cpms[i, group == group.names[2]]) + 1)
        )
}

