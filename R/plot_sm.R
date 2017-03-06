#' Similarity Matrix Heatmap
#'
#' This function displays the pairwise scaled Euclidean distance between samples
#' as a heatmap.
#'
#' @param dat Omic data matrix with rows corresponding to probes and columns
#'   to samples. \code{NA} values are silently removed.
#' @param feat Optional character, factor, numeric, or logical vector of length
#'   equal to sample size. Alternatively, a data frame or list of such vectors,
#'   optionally named. Values are used to color one or several annotation tracks
#'   atop the heatmap.
#' @param main Optional plot title.
#'
#' @details
#' Similarity matrices are a valuable tool for exploratory data analysis. A
#' hierarchical clustering dendrogram atop the figure helps identify potential
#' outliers and/or clusters in the data.
#'
#' @examples
#' mat <- matrix(rnorm(5000), nrow = 1000, ncol = 5)
#' plot_sm(mat, main = "Nothin' Doin'")
#'
#' library(DESeq2)
#' mat <- cbind(matrix(rnbinom(5000, mu = 4, size = 1), nrow = 1000, ncol = 5),
#'              matrix(rnbinom(5000, mu = 4, size = 10), nrow = 1000, ncol = 5))
#' mat <- rlog(mat)
#' grp <- rep(c("A", "B"), each = 5)
#' plot_sm(mat, feat = grp, main = "Somethin' Cookin'")
#'
#' @export
#' @importFrom purrr map_lgl
#' @importFrom wordspace dist.matrix
#' @importFrom NMF aheatmap
#' @import RColorBrewer
#'

plot_sm <- function(dat,
                    feat = NULL,
                    main = NULL) {

  # Preliminaries
  if (is.data.frame(feat)) feat <- as.list(feat)
  else if (!is.list(feat)) feat <- list('Variable' = feat)
  else {
    if (is.null(names(feat))) {
      if (length(feat) == 1L) names(feat) <- 'Variable'
      else names(feat) <- paste('Variable', seq_along(feat))
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
  if (is.null(main)) main <- 'Sample Similarity Matrix'

  # Tidy data
  dat <- getEAWP(dat)$expr
  keep <- rowSums(is.finite(dat)) == ncol(dat)
  dat <- dat[keep, , drop = FALSE]
  dm <- dist.matrix(dat, method = 'euclidean', byrow = FALSE)
  rb <- colorRampPalette(brewer.pal(10, 'RdBu'))(n = 256)

  # Plot
  if (is.null(feat)) {
    aheatmap(dm, col = rb, Rowv = FALSE, main = main,
             distfun = function(x) as.dist(x), hclustfun = 'average')
  } else {
    aheatmap(dm, col = rb, Rowv = FALSE, main = main,
             distfun = function(x) as.dist(x), hclustfun = 'average',
             annCol = feat)
  }

}


