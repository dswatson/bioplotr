#' Omic Heatmap
#'
#' This function visualizes a probe by sample matrix as a heatmap.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples.
#' @param group Optional character or factor vector of length equal to sample
#'   size. Alternatively, a data frame or list of such vectors, optionally
#'   named.
#' @param covar Optional continuous covariate of length equal to sample size.
#'   Alternatively, a data frame or list of such vectors, optionally named.
#' @param dist Distance measure to be used. Currently supports any method
#'   available in \code{\link[stats]{dist}} or \code{\link[stats]{cor}}.
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
#' Heatmaps are a common and intuitive way to display the values of an omic data
#' matrix, especially after top probes have been selected for closer
#' investigation. Hierarchical clustering dendrograms cluster both the rows and
#' the columns, revealing latent structure in the data. Annotation tracks atop
#' the figure may be used to investigate associations with phenotypic features.
#'
#' @examples
#' mat <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
#' grp <- rep(c("A", "B"), each = 5)
#' plot_heatmap(mat, group = grp)
#'
#' @export
#' @importFrom purrr map_lgl
#' @importFrom RColorBrewer brewer.pal
#'

plot_heatmap <- function(dat,
                         group = NULL,
                         covar = NULL,
                          dist = 'pearson',
                     hclustfun = 'average',
                     pal_group = 'd3',
                     pal_covar = 'blues',
                     pal_tiles = 'RdBu',
                         title = NULL) {

  # Preliminaries
  if (!is.null(group)) {
    group <- format_features(dat, group, var_type = 'Categorical')
    grp_cols <- track_cols(group, var_type = 'Categorical')
  } else {
    grp_cols <- NULL
  }
  if (!is.null(covar)) {
    covar <- format_features(dat, covar, var_type = 'Continuous')
    cov_cols <- track_cols(covar, var_type = 'Continuous')
  } else {
    cov_cols <- NULL
  }
  if (!is.null(c(group, covar))) {
    anno <- c(group, covar)
    ann_cols <- c(grp_cols, cov_cols)
  }
  if (!dist %in% c('euclidean', 'maximum', 'manhattan', 'canberra', 'binary',
                   'minkowski', 'pearson', 'kendall', 'spearman')) {
    stop('dist measure not recognized. See ?dist and ?cor for available ',
         'options.')
  }
  if (!hclustfun %in% c('ward.D', 'ward.D2', 'single', 'complete', 'average',
                        'mcquitty', 'median', 'centroid')) {
    stop('hclustfun not recognized. See ?hclust for available options.')
  }
  if (pal_tiles == 'RdBu') {
    pal_tiles <- rev(colorRampPalette(brewer.pal(10L, 'RdBu'))(n = 256L))
  } else if (pal_tiles == 'GrRd') {
    pal_tiles <- colorRampPalette(c('green', 'black', 'red'))(n = 256L)
  } else if (pal_tiles == 'BuYl') {
    pal_tiles <- colorRampPalette(c('blue', 'grey', 'yellow'))(n = 256L)
  }
  if (is.null(title)) {
    title <- 'Omic Heatmap'
  }

  # Tidy data
  dat <- matrixize(dat)

  # Plot
  require(NMF)
  if (is.null(anno)) {
    aheatmap(dat, distfun = dist, scale = 'row', col = pal_tiles,
             hclustfun = hclustfun, main = title, border_color = 'grey60')
  } else {
    aheatmap(dat, distfun = dist, scale = 'row', col = pal_tiles,
             hclustfun = hclustfun, main = title, annCol = anno,
             annColors = ann_cols, border_color = 'grey60')
  }

}


