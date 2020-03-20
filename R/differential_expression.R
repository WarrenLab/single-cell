#' Calculate the cellular detection rate
#'
#' Given a feature X cell table of raw counts, calculate the cellular detection
#' rate for each cell. The cellular detection rate of a cell is a number between
#' 0 and 1 representing the fraction of the total number of features measured
#' that have a count > 0 in that cell. It is useful as a covariate in
#' differential expression analyses.
#'
#' @param counts A feature X cell table of raw counts.
#' @return A vector of cellular detection rate per cell
#' @examples
#' cdr <- calc.cdr(seurat@assays$RNA@counts)
calc.cdr <- function(counts) scale(colMeans(as.matrix(counts) > 0))

#' Run edgeR to perform DE analysis
#'
#' Run edgeR to assign a p-value to each feature for the null hypothesis that
#' this feature is not significantly differently expressed between the cells
#' of a given cluster between the two levels of a condition. This gives edgeR
#' raw counts, and uses cellular detection rate as a covariate. Based on
#' \link{https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLFDetRate.R}
#'
#' @param seurat A seurat object.
#' @param cluster The number of a cluster to perform DE analysis on
#' @param grouping.var The name of the metadata variable used to group the
#'   samples into two groups to compare.
#' @return A \code{data.frame} with row names as feature names, and the columns
#'   \code{p} and \code{fold.change}.
#' @examples
#' de.results <- RunEdgeR(seurat, 0, 'challenge_status')
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

#' Run MAST to perform DE analysis.
#'
#' Run MAST to assign a p-value to each feature for the null hypothesis that
#' this feature is not significantly differently expressed in the cells of a
#' given cluster between two levels of a condition. This function gives MAST
#' the CPM (counts per million) normalized expression level, and includes
#' cellular detection rate as a covariate. Based on
#' \link{https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_MASTcpmDetRate.R}
#'
#' @inheritParams RunEdgeR
#' @return A \code{data.frame} with row names as feature names, and the columns
#'   \code{p} and \code{fold.change}.
#' @examples
#' de.results <- RunMAST(seurat, 0, 'challenge_status')
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

#' Run Wilcoxon tests to perform DE analysis.
#'
#' Use Wilcoxon tests to assign a p-value to each feature for the null
#' hypothesis that this feature is not significantly differently expressed in
#' the cells of a given cluster between two levels of a condition. This function
#' uses the CPM (counts per million) normalized expression levels as input to
#' the Wilcoxon test. Based on
#' \link{https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_Wilcoxon.R}
#'
#' @inheritParams RunEdgeR
#' @return A \code{data.frame} with row names as feature names, and the columns
#'   \code{p} and \code{fold.change}.
#' @examples
#' de.results <- RunMAST(seurat, 0, 'challenge_status')
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

#' Run differential expression analysis on every cluster
#'
#' Use a given function to perform differential expression analysis between two
#' conditions on every cluster in a Seurat object.
#'
#' @param method A function that takes a seurat object, a cluster number, and a
#'   grouping variable, and returns a \code{data.frame} with row names as
#'   feature names, and the columns \code{p} and \code{fold.change}. See, for
#'   example, \code{\link{RunEdgeR}}
#' @inheritParams RunEdgeR
#' @return A \code{data.frame} containing one row for each cluster/feature
#'   combination, with columns \code{p}, \code{fold.change}, \code{cluster}, and
#'   \code{feature}.
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
