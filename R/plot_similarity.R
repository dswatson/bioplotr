#' Similarity Matrix Heatmap
#'
#' This function displays pairwise distances between samples as a heatmap.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples. It is strongly recommended that data be
#'   filtered and normalized prior to plotting. Raw counts stored in \code{
#'   \link[edgeR]{DGEList}} or \code{\link[DESeq2]{DESeqDataSet}} objects are
#'   automatically extracted and transformed to the log2-CPM scale, with a
#'   warning.
#' @param group Optional character or factor vector of length equal to sample
#'   size. Alternatively, a data frame or list of such vectors, optionally
#'   named.
#' @param covar Optional continuous covariate of length equal to sample size.
#'   Alternatively, a data frame or list of such vectors, optionally named.
#' @param top Optional number (if > 1) or proportion (if < 1) of top probes to
#'   be used for distance calculations.
#' @param filter_method String specifying whether to apply a \code{"pairwise"}
#'   or \code{"common"} filter if \code{top} is non-\code{NULL}. See Details.
#' @param dist Distance measure to be used. Currently supports \code{
#'   'euclidean'}, \code{'pearson'}, \code{'MI'}, or \code{'KLD'}. See Details.
#' @param hclustfun The agglomeration method to be used for hierarchical
#'   clustering. Supports any method available in \code{\link[stats]{hclust}}.
#' @param pal_group String specifying the color palette to use if \code{group}
#'   is non-\code{NULL}, or a vector of such strings with length equal to the
#'   number of vectors passed to \code{group}. Options include \code{'ggplot'},
#'   as well as the complete collection of \code{
#'   \href{https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html}{
#'   ggsci}} palettes, which can be identified by name (e.g., \code{'npg'},
#'   \code{'aaas'}, etc.). Alternatively, any character vector of colors with
#'   length equal to the cumulative number of levels in \code{group}.
#' @param pal_covar String specifying the color palette to use if \code{covar}
#'   is non-\code{NULL}, or a vector of such strings with length equal to the
#'   number of vectors passed to \code{covar}. Options include \code{'blues'},
#'   \code{'greens'}, \code{'purples'}, \code{'greys'}, \code{'oranges'}, and
#'   \code{'reds'}. Alternatively, any character vector of colors representing
#'   a smooth gradient, or a list of such vectors with length equal to the
#'   number of continuous variables to visualize.
#' @param pal_tiles String specifying the color palette to use for heatmap
#'   tiles. Options include \code{'RdBu'} for red to blue gradient, \code{
#'   'GrRd'} for green to red gradient, and \code{'BuYl'} for blue to yellow
#'   gradient. Alternatively, any user-supplied color palette is acceptable.
#' @param title Optional plot title.
#'
#' @details
#' Similarity matrices are a valuable tool for exploratory data analysis. A
#' hierarchical clustering dendrogram atop the figure helps identify potential
#' clusters and/or outliers in the data. Annotation tracks can help investigate
#' associations with phenotypic features.
#'
#' Different distance metrics and agglomeration methods can lead to different
#' results. The default settings, which perform average linkage hierarchical
#' clustering on a Euclidean distance matrix, are mathematically straightforward
#' and commonly used for omic EDA. Complete linkage is also fairly common.
#'
#' Pearson distance, defined as 1 - the Pearson correlation, is another popular
#' method for evaluating sample similarity. Mutual information and
#' Kullback-Leibler divergence are more complicated distance metrics that
#' require some simplifying assumptions to be efficiently applied to continuous
#' data distributions. See \code{\link[bioDist]{MIdist}} and \code{
#' \link[bioDist]{KLdist.matrix}} for more info.
#'
#' The \code{top} argument optionally filters data using either probewise
#' variance (if \code{filter_method = "common"}) or the leading fold change
#' method of Smyth et al. (if \code{filter_method = "pairwise"}). See \code{
#' \link{plot_mds}} for more details.
#'
#' @examples
#' mat <- matrix(rnorm(5000), nrow = 1000, ncol = 5)
#' plot_similarity(mat, title = 'Nothin' Doin'')
#'
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet()
#' dds <- rlog(dds)
#' plot_similarity(dds, group = colData(dds)$condition,
#'                 title = 'Somethin' Cookin'')
#'
#' @export
#' @importFrom matrixStats rowMedians
#' @importFrom RColorBrewer brewer.pal
#'

plot_similarity <- function(dat,
                            group = NULL,
                            covar = NULL,
                              top = NULL,
                    filter_method = 'pairwise',
                             dist = 'euclidean',
                        hclustfun = 'average',
                        pal_group = 'npg',
                        pal_covar = 'blues',
                        pal_tiles = 'RdBu',
                            title = NULL) {

  # Preliminaries
  if (!is.null(group)) {
    group <- format_features(dat, group, var_type = 'Categorical')
    grp_cols <- track_cols(group, pal_group, var_type = 'Categorical')
  } else {
    grp_cols <- NULL
  }
  if (!is.null(covar)) {
    covar <- format_features(dat, covar, var_type = 'Continuous')
    cov_cols <- track_cols(covar, pal_covar, var_type = 'Continuous')
  } else {
    cov_cols <- NULL
  }
  if (!is.null(c(group, covar))) {
    anno <- c(group, covar)
    ann_cols <- c(grp_cols, cov_cols)
  } else {
    anno <- NULL
  }
  if (!dist %in% c('euclidean', 'pearson', 'MI', 'KLD')) {
    stop('dist must be one of "euclidean", "pearson", "MI", or "KLD".')
  }
  if (!hclustfun %in% c('ward.D', 'ward.D2', 'single', 'complete', 'average',
                        'mcquitty', 'median', 'centroid')) {
    stop('hclustfun not recognized. See ?hclust for available options.')
  }
  if (pal_tiles == 'RdBu') {
    pal_tiles <- colorRampPalette(brewer.pal(10L, 'RdBu'))(n = 256L)
  } else if (pal_tiles == 'GrRd') {
    pal_tiles <- colorRampPalette(c('green', 'black', 'red'))(n = 256L)
  } else if (pal_tiles == 'BuYl') {
    pal_tiles <- colorRampPalette(c('blue', 'grey', 'yellow'))(n = 256L)
  }
  if (is.null(title)) {
    title <- 'Sample Similarity Matrix'
  }

  # Tidy data
  dat <- matrixize(dat)
  dat <- dat - rowMedians(dat)                   # Median center data
  dm <- dist_mat(dat, top, filter_method, dist)

  # Build plot
  require(NMF)
  if (is.null(anno)) {
    aheatmap(dm, col = pal_tiles, Rowv = FALSE, revC = TRUE, main = title,
             distfun = function(x) as.dist(x), hclustfun = hclustfun,
             border_color = 'grey60')
  } else {
    aheatmap(dm, col = pal_tiles, Rowv = FALSE, revC = TRUE, main = title,
             distfun = function(x) as.dist(x), hclustfun = hclustfun,
             annCol = anno, annColors = ann_cols, border_color = 'grey60')
  }

}


