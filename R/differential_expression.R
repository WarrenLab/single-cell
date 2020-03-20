calc.cdr <- function(counts) scale(colMeans(as.matrix(counts) > 0))

# this is based on https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLFDetRate.R
RunEdgeR <- function(
    seurat,
    cluster,
    grouping.var
) {
    this.cluster.only <- seurat[, which(seurat[['seurat_clusters']] == cluster)]
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
    this.cluster.only <- seurat[, which(seurat[['seurat_clusters']] == cluster)]
    counts <- this.cluster.only@assays$RNA@counts
    group <- as.matrix(this.cluster.only[[grouping.var]])
    group.names <- levels(factor(group[,1]))

    cdr <- calc.cdr(counts)
    dge <- DGEList(counts)
    dge <- edgeR::calcNormFactors(dge)
    cpms <- edgeR::cpm(dge)
    sca <- FromMatrix(
        exprsArray = log2(cpms + 1),
        cData = data.frame(
            wellKey = rownames(group),
            group = group,
            cdr = cdr
        )
    )
    zlmdata <- zlm(as.formula(paste("~cdr +", grouping.var)), sca)
    mast <- lrTest(zlmdata, grouping.var)

    idx <- 1:nrow(mast[,'hurdle',])
    names(idx) <- rownames(mast[,'hurdle',])
    fold.change <- sapply(
        idx, function(i) log2(mean(cpms[i, group == group.names[2]]) + 1)
                       - log2(mean(cpms[i, group == group.names[1]]) + 1)
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
    this.cluster.only <- seurat[, which(seurat[['seurat_clusters']] == cluster)]
    counts <- this.cluster.only@assays$RNA@counts
    group <- as.matrix(this.cluster.only[[grouping.var]])
    group.names <- levels(factor(group[,1]))

    dge <- DGEList(counts)
    dge <- edgeR::calcNormFactors(dge)
    cpms <- edgeR::cpm(dge)
    idx <- 1:nrow(cpms)
    names(idx) <- rownames(cpms)
    wilcox.p <- sapply(
        idx, function(i) wilcox.test(cpms[i, ] ~ group)$p.value)
    fold.change <- sapply(
        idx, function(i) log2(mean(cpms[i, group == group.names[2]]) + 1)
                       - log2(mean(cpms[i, group == group.names[1]]) + 1)
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
        df2$feature <- rownames(df2)
        rbind(df, df2, make.row.names = FALSE)
    }, levels(seurat$seurat_clusters), data.frame())
}
