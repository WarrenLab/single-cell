#' Calculate the cellular detection rate
#'
#' Given a feature X cell table of raw counts, calculate the scaled cellular detection
#' rate for each cell. The cellular detection rate of a cell is a number between
#' 0 and 1 representing the fraction of the total number of features measured
#' that have a count > 0 in that cell. It is useful as a covariate in
#' differential expression analyses. This function returns a vector of scaled CDRs.
#'
#' @param counts A feature X cell table of raw counts.
#' @return A vector of cellular detection rate per cell
#' @examples
#' cdr <- calc.cdr(seurat@assays$RNA@counts)
#' @export
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
#' @param cluster.var The name of the metadata variable containing cluster
#'   assignments ('seurat_clusters' by default)
#' @return A \code{data.frame} with row names as feature names, and the columns
#'   \code{p} and \code{fold.change}.
#' @examples
#' de.results <- RunEdgeR(seurat, 0, 'challenge_status')
#' @export
RunEdgeR <- function(
    seurat,
    cluster,
    grouping.var,
    cluster.var = 'seurat_clusters'
) {
    this.cluster.only <- seurat[, which(seurat[[cluster.var]] == cluster)]
    counts <- this.cluster.only@assays$RNA@counts
    group <- as.matrix(this.cluster.only[[grouping.var]])

    dge <- edgeR::DGEList(counts, group = t(group))
    # perform TMM normalization
    dge <- edgeR::calcNormFactors(dge)
    # calculate cellular detection rate to use as covariate
    cdr <- calc.cdr(counts)
    design <- model.matrix(~ cdr + group)
    dge <- edgeR::estimateDisp(dge, design = design)
    fit <- edgeR::glmQLFit(dge, design = design)
    qlf <- edgeR::glmQLFTest(fit)
    tt <- edgeR::topTags(qlf, n = Inf)

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
#' @export
RunMAST <- function(
    seurat,
    cluster,
    grouping.var,
    cluster.var = 'seurat_clusters'
) {
    this.cluster.only <- seurat[, which(seurat[[cluster.var]] == cluster)]
    counts <- this.cluster.only@assays$RNA@counts
    group <- as.matrix(this.cluster.only[[grouping.var]])
    group.names <- levels(factor(group[,1]))

    cdr <- calc.cdr(counts)
    dge <- edgeR::DGEList(counts)
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
    zlmdata <- MAST::zlm(as.formula(paste("~cdr +", grouping.var)), sca)
    mast <- MAST::lrTest(zlmdata, grouping.var)

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
#' @export
RunWilcoxon <- function(
    seurat,
    cluster,
    grouping.var,
    cluster.var = 'seurat_clusters'
) {
    this.cluster.only <- seurat[, which(seurat[[cluster.var]] == cluster)]
    counts <- this.cluster.only@assays$RNA@counts
    group <- as.matrix(this.cluster.only[[grouping.var]])
    group.names <- levels(factor(group[,1]))

    dge <- edgeR::DGEList(counts)
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
#' @export
FindAllClusterDE <- function(
    seurat,
    grouping.var,
    method = RunEdgeR,
    cluster.var = 'seurat_clusters'
) {
    Reduce(function(df, cluster) {
        df2 <- method(
            seurat,
            cluster = cluster,
            grouping.var = grouping.var,
            cluster.var = cluster.var
        )
        df2$cluster <- cluster
        df2$feature <- rownames(df2)
        rbind(df, df2, make.row.names = FALSE)
    }, Filter(
        function(c) is_contrastable(seurat, grouping.var, c, cluster.var),
        levels(seurat@meta.data[,cluster.var])
    ), data.frame())
}

#' Make a heatmap of the most differentially expressed genes in a cluster
#'
#' Given a seurat object, a table of differential expression fold changes and
#' p-values, a cluster number, and a grouping variable name, make a heatmap
#' that shows per-cell expression of the genes with the smallest differential
#' expression p-values, ordered by fold change.
#'
#' @param seurat A seurat object with defined clusters
#' @param diff.exp A table of differential expression p-values and fold changes
#'   for all clusters, with columns \code{feature}, \code{cluster},
#'   \code{fold.change}, and \code{p}.
#' @param cluster The number of the cluster to plot expression for
#' @param grouping.var The name of the metadata variable in \code{seurat} that
#'   cells were grouped by in differential expression analysis
#' @param min.fc Minimum absolute value of fold change to include a feature in
#'   heatmap
#' @param max.p Maximum differential expression p-value to include a feature in
#'   heatmap
#' @param max.features Maximum number of features to include in the heatmap
#' @return A ggplot object containing a heatmap with expression of at most
#'   \code{max.features} genes with the smallest differential expression
#'   p-values, ordered by fold change, for all cells in \code{cluster}, grouped
#'   by \code{grouping.var}.
#' @export
make.diffexp.heatmap <- function(
    seurat,
    diff.exp,
    cluster,
    grouping.var,
    min.fc = 1,
    max.p = 1e-3,
    max.features = 30,
    cluster.var = 'seurat_clusters'
) {
    selected.rows <- diff.exp[diff.exp$cluster == cluster
                              & abs(diff.exp$fold.change) > min.fc
                              & diff.exp$p < max.p, ]
    if (dim(selected.rows)[2] == 0) stop("empty heatmap!")
    selected.rows <- selected.rows[order(selected.rows$fold.change), ]
    features <- dplyr::top_n(selected.rows, n = max.features, wt = -p)$feature
    this.cluster.only <- seurat[, which(seurat[[cluster.var]] == cluster)]
    Seurat::DoHeatmap(this.cluster.only,
                      features = features,
                      group.by = grouping.var) + NoLegend()
}

is_contrastable <- function(seurat, grouping.var, cluster, cluster.var) {
    this.cluster.only <- seurat[, which(seurat[[cluster.var]] == cluster)]
    group <- as.matrix(this.cluster.only[[grouping.var]])
    return(length(levels(factor(group[,1]))) != 1)
}
