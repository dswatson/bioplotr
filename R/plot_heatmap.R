#' Omic Heatmap
#'
#' This function visualizes a probe by sample matrix as a heatmap.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples.
#' @param feat Optional character, factor, numeric, or logical vector of length
#'   equal to sample size. Alternatively, a data frame or list of such vectors,
#'   optionally named. Values are used to color one or several annotation tracks
#'   atop the heatmap.
#' @param col Color scheme to use for heatmap tiles. Options are \code{"rb"} (for
#'   red to blue gradient) or \code{"gr"} (for green to red gradient).
#' @param main Optional plot title.
#'
#' @details
#' Heatmaps are a common and intuitive way to display the values of an omic data
#' matrix, especially after top probes have been selected for closer investigation.
#' Hierarchical clustering dendrograms cluster both the rows and the columns,
#' revealing latent structure in the data. Annotation tracks atop the figure may be
#' used to investigate associations with phenotypic features.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 10), nrow = 1000, ncol = 10)
#' grp <- rep(c("A", "B"), each = 5)
#' plot_heatmap(mat, feat = grp)
#'
#' @export
#' @importFrom purrr map_lgl
#' @importFrom limma getEAWP
#' @importFrom NMF aheatmap
#' @import RColorBrewer
#'

plot_heatmap <- function(dat,
                         feat = NULL,
                         col = 'rb',
                         main = NULL) {

  # Preliminaries
  if (col == 'rb') {
    cols <- colorRampPalette(brewer.pal(10L, 'RdBu'))(n = 256L)
  } else if (col == 'gr') {
    cols <- colorRampPalette(c('green', 'black', 'red'))(n = 256L)
  } else {
    stop('col must be one of "rb" or "gr".')
  }
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
  if (is.null(main)) {
    main <- 'Omic Heatmap'
  }

  # Tidy data
  dat <- getEAWP(dat)$expr
  keep <- rowSums(is.finite(dat)) == ncol(dat)
  dat <- dat[keep, , drop = FALSE]

  # Plot
  if (is.null(feat)) {
    aheatmap(dat, distfun = 'pearson', scale = 'row', col = cols,
             hclustfun = 'average', main = main)
  } else {
    aheatmap(dat, distfun = 'pearson', scale = 'row', col = cols,
             hclustfun = 'average', main = main, annCol = feat)
  }

}


