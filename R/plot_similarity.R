#' Similarity Matrix Heatmap
#'
#' This function displays pairwise distances between samples as a heatmap.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples. It is strongly recommended that data be
#'   normalized and filtered prior to plotting the sample similarity matrix.
#'   For count data, this means undergoing some sort of variance stabilizing
#'   transformation, such as \code{\link[edgeR]{cpm} (with \code{log = TRUE}),
#'   \link[DESeq2]{vst}, \link[DESeq2]{rlog}}, etc.
#' @param feat Optional character, factor, numeric, or logical vector of length
#'   equal to sample size. Alternatively, a data frame or list of such vectors,
#'   optionally named. Values are used to color one or several annotation tracks
#'   atop the heatmap.
#' @param dist Distance measure to be used. Currently supports \code{"euclidean",
#'   "pearson", "MI",} or \code{"KLD"}. See Details.
#' @param hclustfun The agglomeration method to be used for hierarchical clustering.
#'   See \code{\link[stats]{hclust}} for available options.
#' @param title Optional plot title.
#'
#' @details
#' Similarity matrices are a valuable tool for exploratory data analysis. A
#' hierarchical clustering dendrogram atop the figure helps identify potential
#' clusters and/or outliers in the data. Annotation tracks can help investigate
#' associations with phenotypic features.
#'
#' Different distance measures and agglomeration methods can lead to different
#' results. The default settings, which perform average linkage hierarchical
#' clustering on a Euclidean distance matrix, are mathematically straightforward
#' and commonly used for omic EDA. Complete linkage is also fairly common for
#' hierarchical clustering, while other options (single linkage, Ward, etc.) are less
#' informative.
#'
#' Pearson distance, defined as 1 - the Pearson correlation, is another popular
#' method for evaluating sample similarity. Mutual information and Kullback-Leibler
#' divergence are more complicated distance metrics that require some simplifying
#' assumptions to be efficiently applied to continuous data distributions. See
#' \code{\link[bioDist]{MIdist}} and \code{\link[bioDist]{KLdist.matrix}}
#' for more info.
#'
#' @examples
#' mat <- matrix(rnorm(5000), nrow = 1000, ncol = 5)
#' plot_similarity(mat, title = "Nothin' Doin'")
#'
#' library(DESeq2)
#' mat <- cbind(matrix(rnbinom(5000, mu = 4, size = 1), nrow = 1000, ncol = 5),
#'              matrix(rnbinom(5000, mu = 4, size = 10), nrow = 1000, ncol = 5))
#' mat <- rlog(mat)
#' grp <- rep(c("A", "B"), each = 5)
#' plot_similarity(mat, feat = grp, title = "Somethin' Cookin'")
#'
#' @export
#' @importFrom purrr map_lgl
#' @importFrom limma getEAWP
#' @importFrom wordspace dist.matrix
#' @importFrom bioDist MIdist KLdist.matrix
#' @importFrom NMF aheatmap
#' @import RColorBrewer
#'

plot_similarity <- function(dat,
                            feat = NULL,
                            dist = 'euclidean',
                       hclustfun = 'average',
                           title = NULL) {

  # Preliminaries
  if (is.data.frame(feat)) {
    feat <- as.list(feat)
  } else if (!is.list(feat)) {
    feat <- list('Variable' = feat)
  } else {
    if (is.null(names(feat))) {
      if (length(feat) == 1L) {
        names(feat) <- 'Variable'
      } else {
        names(feat) <- paste('Variable', seq_along(feat))
      }
    }
  }
  if (any(map_lgl(seq_along(feat), function(j) {
    length(feat[[j]]) != ncol(dat)
  }))) {
    stop('feat length must match number of samples in dat.')
  }
  if (any(map_lgl(seq_along(feat), function(j) {
    if (is.numeric(feat[[j]])) var(feat[[j]]) == 0L
    else length(unique(feat[[j]])) == 1L
  }))) {
    stop('feat is invariant.')
  }
  if (!dist %in% c('euclidean', 'pearson', 'MI', 'KLD')) {
    stop('dist must be one of "euclidean", "pearson", "MI", or "KLD".')
  }
  if (!hclustfun %in% c('ward.D', 'ward.D2', 'single', 'complete', 'average',
                        'mcquitty', 'median', 'centroid')) {
    stop('hclustfun must be one of "ward.D", "ward.D2", "single", "complete", ',
         '"average", "mcquitty", "median", or "centroid". See ?hclust.')
  }
  if (is.null(title)) {
    title <- 'Sample Similarity Matrix'
  }

  # Tidy data
  dat <- getEAWP(dat)$expr
  keep <- rowSums(is.finite(dat)) == ncol(dat)
  dat <- dat[keep, , drop = FALSE]
  if (dist == 'euclidean') {
    dm <- dist.matrix(t(dat), method = 'euclidean')
  } else if (dist == 'pearson') {
    dm <- 1 - cor(mat)
  } else if (dist == 'MI') {
    dm <- as.matrix(MIdist(t(dat)))
  } else {
    dm <- as.matrix(KLdist.matrix(t(dat)))
  }

  # Plot
  rb <- colorRampPalette(brewer.pal(10, 'RdBu'))(n = 256)
  if (is.null(feat)) {
    aheatmap(dm, col = rb, Rowv = FALSE, main = title,
             distfun = function(x) as.dist(x), hclustfun = hclustfun)
  } else {
    aheatmap(dm, col = rb, Rowv = FALSE, main = title,
             distfun = function(x) as.dist(x), hclustfun = hclustfun,
             annCol = feat)
  }

}


# Use shiny to toggle between distance measures and agglomeration methods
