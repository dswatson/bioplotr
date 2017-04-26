#' Similarity Matrix Heatmap
#'
#' This function displays pairwise distances between samples as a heatmap.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples. It is strongly recommended that data be
#'   normalized and filtered prior to plotting the sample similarity matrix.
#'   For count data, this means undergoing some sort of variance stabilizing
#'   transformation, such as \code{\link[edgeR]{cpm}} (with \code{log = TRUE}),
#'   \code{\link[DESeq2]{vst}}, \code{\link[DESeq2]{rlog}}, etc.
#' @param anno Optional character, factor, numeric, or logical vector of length
#'   equal to sample size. Alternatively, a data frame or list of such vectors,
#'   optionally named. Values are used to color one or several annotation tracks
#'   atop the heatmap.
#' @param top Optional number (if > 1) or proportion (if < 1) of top probes to be used
#'   for distance calculations. See Details.
#' @param dist Distance measure to be used. Currently supports \code{"euclidean",
#'   "pearson", "MI",} or \code{"KLD"}. See Details.
#' @param hclustfun The agglomeration method to be used for hierarchical clustering.
#'   See \code{\link[stats]{hclust}} for available options.
#' @param col Color palette to use for heatmap tiles. Preset options include \code{
#'   "RdBu"} for red to blue gradient, \code{"GrRd"} for green to red gradient, and
#'   \code{"BuYl"} for blue to yellow gradient. Alternatively, any user-supplied
#'   color palette is acceptable.
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
#' The \code{top} argument optionally filters probes using the leading fold change
#' method of Smyth et al. See \code{\link{plot_mds}} for more details.
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
#' plot_similarity(mat, anno = grp, title = "Somethin' Cookin'")
#'
#' @export
#' @importFrom purrr map_lgl
#' @importFrom edgeR calcNormFactors cpm
#' @importFrom DESeq2 sizeFactors normalizationFactors estimateSizeFactors
#' @importFrom SummarizedExperiment assay
#' @importFrom limma getEAWP
#' @importFrom wordspace dist.matrix
#' @importFrom bioDist MIdist KLdist.matrix
#' @importFrom NMF aheatmap
#' @import RColorBrewer
#'

plot_similarity <- function(dat,
                            anno = NULL,
                             top = NULL,
                            dist = 'euclidean',
                       hclustfun = 'average',
                             col = 'RdBu',
                           title = NULL) {

  # Preliminaries
  if (!is.null(anno)) {
    if (is.data.frame(anno)) {
      anno <- as.list(anno)
    } else if (!is.list(anno)) {
      anno <- list('Variable' = anno)
    } else {
      if (is.null(names(anno))) {
        if (length(anno) == 1L) {
          names(anno) <- 'Variable'
        } else {
          names(anno) <- paste('Variable', seq_along(anno))
        }
      }
    }
    if (!all(map_lgl(seq_along(anno), function(j) {
      length(anno[[j]]) == ncol(dat)
    }))) {
      stop('anno length must match number of samples in dat.')
    }
    if (any(map_lgl(seq_along(anno), function(j) {
      if (is.numeric(anno[[j]])) var(anno[[j]]) == 0L
      else length(unique(anno[[j]])) == 1L
    }))) {
      stop('anno is invariant.')
    }
  }
  if (!dist %in% c('euclidean', 'pearson', 'MI', 'KLD')) {
    stop('dist must be one of "euclidean", "pearson", "MI", or "KLD".')
  }
  if (!hclustfun %in% c('ward.D', 'ward.D2', 'single', 'complete', 'average',
                        'mcquitty', 'median', 'centroid')) {
    stop('hclustfun must be one of "ward.D", "ward.D2", "single", "complete", ',
         '"average", "mcquitty", "median", or "centroid". See ?hclust.')
  }
  if (col == 'RdBu') {
    col <- colorRampPalette(brewer.pal(10L, 'RdBu'))(n = 256L)
  } else if (col == 'GrRd') {
    col <- colorRampPalette(c('green', 'black', 'red'))(n = 256L)
  } else if (col == 'BuYl') {
    col <- colorRampPalette(c('blue', 'grey', 'yellow'))(n = 256L)
  }
  if (is.null(title)) {
    title <- 'Sample Similarity Matrix'
  }

  # Tidy data
  if (is(dat, 'DGEList')) {
    keep <- rowSums(dat$counts) > 1L             # Minimal count filter
    dat <- dat[keep, ]
    if (is.null(dat$samples$norm.factors) |      # Calculate size factors
        all(dat$samples$norm.factors == 1L)) {
      dat <- calcNormFactors(dat)
    }
    dat <- cpm(dat, log = TRUE, prior.count = 1L)
    warning('Transforming raw counts to log2-CPM scale.')
  } else if (is(dat, 'DESeqDataSet')) {
    if (is.null(sizeFactors(dat)) & is.null(normalizationFactors(dat))) {
      dat <- estimateSizeFactors(dat)            # Normalize counts
    }
    dat <- counts(dat, normalized = TRUE)
    keep <- rowMeans(dat) > 0L                   # Minimal count filter
    dat <- dat[keep, , drop = FALSE]
    dat <- cpm(dat, log = TRUE, prior.count = 1L)
    warning('Transforming raw counts to log2-CPM scale.')
  } else if (is(dat, 'DESeqTransform')) {
    dat <- assay(dat)
  } else {
    dat <- getEAWP(dat)$expr
    keep <- rowSums(is.finite(dat)) == ncol(dat)
    dat <- dat[keep, , drop = FALSE]
  }
  dat <- sweep(dat, 1L, apply(dat, 1L, median))  # Median center data
  if (!is.null(top)) {
    if (top > 1L) {
      if (top > nrow(dat)) {
        warning(paste('top exceeds nrow(dat), at least after removing probes with
                      missing values and/or applying a minimal expression filter.
                      Proceeding with the complete', nrow(dat), 'x', ncol(dat), 'matrix.'))
        top <- NULL
      }
    } else {
      top <- round(top * nrow(dat))
    }
  }
  if (dist == 'euclidean') {
    if (is.null(top)) {
      dm <- dist.matrix(t(dat), method = 'euclidean')
    } else {
      dm <- matrix(0L, nrow = ncol(dat), ncol = ncol(dat))
      for (i in 2L:ncol(dat)) {
        for (j in 1L:(i - 1L)) {
          tops <- order((dat[, i] - dat[, j])^2, decreasing = TRUE)[seq_len(top)]
          dm[i, j] <- sqrt(sum((dat[tops, i] - dat[tops, j])^2))
        }
      }
      dm <- pmax(dm, t(dm))
    }
  } else if (dist == 'pearson') {
    if (is.null(top)) {
      dm <- 1 - cor(dat)
    } else {
      dm <- matrix(0L, nrow = ncol(dat), ncol = ncol(dat))
      for (i in 2L:ncol(dat)) {
        for (j in 1L:(i - 1L)) {
          tops <- order((dat[, i] - dat[, j])^2, decreasing = TRUE)[seq_len(top)]
          dm[i, j] <- 1 - cor(dat[tops, i], dat[tops, j])
        }
      }
      dm <- pmax(dm, t(dm))
    }
  } else if (dist == 'MI') {
    if (is.null(top)) {
      dm <- as.matrix(MIdist(t(dat)))
    } else {
      dm <- matrix(0L, nrow = ncol(dat), ncol = ncol(dat))
      for (i in 2L:ncol(dat)) {
        for (j in 1L:(i - 1L)) {
          tops <- order((dat[, i] - dat[, j])^2, decreasing = TRUE)[seq_len(top)]
          dm[i, j] <- max(as.matrix(MIdist(t(dat[tops, c(i, j)]))))
        }
      }
      dm <- pmax(dm, t(dm))
    }
  } else if (dist == 'KLD') {
    if (is.null(top)) {
      dm <- as.matrix(KLdist.matrix(t(dat)))
    } else {
      dm <- matrix(0L, nrow = ncol(dat), ncol = ncol(dat))
      for (i in 2L:ncol(dat)) {
        for (j in 1L:(i - 1L)) {
          tops <- order((dat[, i] - dat[, j])^2, decreasing = TRUE)[seq_len(top)]
          dm[i, j] <- max(as.matrix(KLdist.matrix(t(dat[tops, c(i, j)]))))
        }
      }
      dm <- pmax(dm, t(dm))
    }
  }

  # Plot
  if (is.null(anno)) {
    aheatmap(dm, col = col, Rowv = FALSE, revC = TRUE, main = title,
             distfun = function(x) as.dist(x), hclustfun = hclustfun,
             border_color = 'grey60')
  } else {
    aheatmap(dm, col = col, Rowv = FALSE, revC = TRUE, main = title,
             distfun = function(x) as.dist(x), hclustfun = hclustfun,
             annCol = anno, border_color = 'grey60')
  }

}


