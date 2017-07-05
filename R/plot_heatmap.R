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
#' @param dist Distance measure to be used. Supports all methods available in
#'   \code{\link[stats]{dist}} and \code{\link[stats]{cor}}.
#' @param hclustfun The agglomeration method to be used for hierarchical
#'   clustering. Supports all methods available in \code{\link[stats]{hclust}}.
#' @param pal_group String specifying the color palette to use if \code{group}
#'   is non-\code{NULL}, or a vector of such strings with length equal to the
#'   number of vectors passed to \code{group}. Options include \code{"ggplot"},
#'   all qualitative color schemes available in \code{RColorBrewer}, and the
#'   complete collection of \code{\href{http://bit.ly/2bxnuGB}{ggsci}} palettes.
#'   Alternatively, a character vector of colors with length equal to the
#'   cumulative number of levels in \code{group}.
#' @param pal_covar String specifying the color palette to use if \code{covar}
#'   is non-\code{NULL}, or a vector of such strings with length equal to the
#'   number of vectors passed to \code{covar}. Options include all sequential
#'   color schemes available in \code{RColorBrewer}. Alternatively, a
#'   character vector of colors representing a smooth gradient, or a list of
#'   such vectors with length equal to the number of continuous variables to
#'   visualize.
#' @param pal_tiles String specifying the color palette to use for heatmap
#'   tiles. Options include all diverging color schemes available in \code{
#'   RColorBrewer}. Alternatively, a character vector of at least two colors.
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
#' @importFrom RColorBrewer brewer.pal
#'

plot_heatmap <- function(dat,
                         group = NULL,
                         covar = NULL,
                          dist = 'pearson',
                     hclustfun = 'average',
                     pal_group = 'npg',
                     pal_covar = 'Blues',
                     pal_tiles = 'RdBu',
                         title = NULL) {

  # Preliminaries
  if (!(group %>% is.null)) {
    group <-  dat %>% format_features(group, var_type = 'Categorical')
    grp_cols <- group %>% track_cols(pal_group, var_type = 'Categorical')
  } else {
    grp_cols <- NULL
  }
  if (!(covar %>% is.null)) {
    covar <- format_features(dat, covar, var_type = 'Continuous')
    cov_cols <- track_cols(covar, pal_covar, var_type = 'Continuous')
  } else {
    cov_cols <- NULL
  }
  if (!c(group, covar) %>% is.null) {
    anno <- c(group, covar)
    ann_cols <- c(grp_cols, cov_cols)
  }
  d <- c('euclidean', 'maximum', 'manhattan', 'canberra', 'minkowski',
         'pearson', 'kendall', 'spearman')
  if (!dist %in% d) {
    stop('dist must be one of ', stringify(d, 'or'), '.')
  }
  hclusts <- c('ward.D', 'ward.D2', 'single', 'complete', 'average',
               'mcquitty', 'median', 'centroid')
  if (!hclustfun %in% hclusts) {
    stop('hclustfun must be one of ', stringify(hclusts, 'or'), '. ',
         'See ?hclust.')
  }
  pal_tiles <- colorize(pal_tiles, var_type = 'Continuous')
  if (title %>% is.null) {
    title <- 'Omic Heatmap'
  }

  # Tidy data
  dat <- matrixize(dat)

  # Plot
  suppressPackageStartupMessages(require(NMF))
  if (anno %>% is.null) {
    aheatmap(dat, distfun = dist, hclustfun = hclustfun, scale = 'row',
             col = pal_tiles, main = title, border_color = 'grey60')
  } else {
    aheatmap(dat, distfun = dist, hclustfun = hclustfun, scale = 'row',
             annCol = anno, annColors = ann_cols, col = pal_tiles,
             main = title, border_color = 'grey60')
  }

}


### REPLACE AHEATMAP WITH PHEATMAP ###

