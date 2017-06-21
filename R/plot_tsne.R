#' t-SNE Plot
#'
#' This function plots a low-dimensional projection of an omic data matrix using
#' t-distributed stochastic neighbor embedding.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples. It is strongly recommended that data be
#'   filtered and normalized prior to plotting. Raw counts stored in \code{
#'   \link[edgeR]{DGEList}} or \code{\link[DESeq2]{DESeqDataSet}} objects are
#'   automatically extracted and transformed to the log2-CPM scale, with a
#'   warning.
#' @param group Optional character or factor vector of length equal to sample
#'   size, or up to two such vectors organized into a list or data frame. Supply
#'   legend title(s) by passing a named list or data frame.
#' @param covar Optional continuous covariate. If non-\code{NULL}, then plot can
#'   render at most one \code{group} variable. Supply legend title by passing
#'   a named list or data frame.
#' @param top Optional number (if > 1) or proportion (if < 1) of top probes to
#'   be used for t-SNE.
#' @param filter_method String specifying whether to apply a \code{"pairwise"}
#'   or \code{"common"} filter if \code{top} is non-\code{NULL}. See Details.
#' @param dims Vector specifying which dimensions to plot. Must be of length
#'   two unless \code{D3 = TRUE}.
#' @param perplexity How many nearest neighbors should the algorithm consider
#'   when building projections?
#' @param theta Speed/accuracy tradeoff of the Barnes-Hut algorithm. See Details.
#' @param max_iter Maximum number of iterations over which to minimize the loss
#'   function. See Details.
#' @param label Label data points by sample name? Defaults to \code{FALSE}
#'   unless \code{group} and \code{covar} are both \code{NULL}. If \code{TRUE},
#'   then plot can render at most one phenotypic feature.
#' @param pal_group String specifying the color palette to use if \code{group}
#'   is non-\code{NULL}. Options include \code{"ggplot"}, as well as the
#'   complete collection of \code{
#'   \href{https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html}{
#'   ggsci}} palettes, which can be identified by name (e.g., \code{"npg"},
#'   \code{"aaas"}, etc.). Alternatively, any character vector of colors with
#'   length equal to the number of levels in \code{group}.
#' @param pal_covar String specifying the color palette to use if \code{covar}
#'   is non-\code{NULL}. Options include \code{"blues"}, \code{"greens"}, \code{
#'   "purples"}, \code{"greys"}, \code{"oranges"}, and \code{"reds"}.
#'   Alternatively, any character vector of colors representing a smooth
#'   gradient.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"right"}, \code{
#'   "left"}, \code{"top"}, \code{"bottom"}, \code{"topright"}, \code{
#'   "topleft"}, \code{"bottomright"}, or \code{"bottomleft"}.
#' @param hover Show sample name by hovering mouse over data point? If \code{
#'   TRUE}, the plot is rendered in HTML and will either open in your browser's
#'   graphic display or appear in the RStudio viewer.
#' @param D3 Render plot in three dimensions?
#' @param ... Additional arguments to be passed to \code{\link[Rtsne]{Rtsne}}.
#'
#' @details
#' This function plots the samples of an omic data matrix in a two- or
#' three-dimensional stochastic neighbor subspace. t-SNE is a popular machine
#' learning method for unsupervised cluster detection. It can also aid in
#' spotting potential outliers, and generally helps to visualize the latent
#' structure of a data set.
#'
#' \code{plot_tsne} relies on a C++ implementation of the Barnes-Hut algorithm,
#' which vastly accelerates the original t-SNE projection method. An exact t-SNE
#' plot may be rendered by setting \code{theta = 0}. Briefly, the algorithm
#' computes samplewise similarities based on distances in the original \emph{
#' p}-dimensional space (where \emph{p} = the number of probes); generates a
#' low-dimensional embedding of the samples based on the user-defined \code{
#' perplexity} parameter; and iteratively minimizes the Kullback-Leibler
#' divergence between these two distributions using an efficient tree search.
#' See \code{\link[Rtsne]{Rtsne}} for more details. A thorough introduction to
#' and explication of the original t-SNE method and the Barnes-Hut approximation
#' may be found in the references below.
#'
#' The \code{top} argument optionally filters data using either probewise
#' variance (if \code{filter_method = "common"}) or the leading fold change
#' method of Smyth et al. (if \code{filter_method = "pairwise"}). See \code{
#' \link{plot_mds}} for more details.
#'
#' @references
#' van der Maaten, L.J.P. (2014).
#' \href{http://www.jmlr.org/papers/volume15/vandermaaten14a/source/vandermaaten14a.pdf}{
#' Accelerating t-SNE using Tree-Based Algorithms}. \emph{Journal of Machine
#' Learning Research}, \emph{15}: 3221-3245.
#'
#' van der Maaten, L.J.P. & Hinton, G.E. (2008).
#' \href{http://www.jmlr.org/papers/volume9/vandermaaten08a/vandermaaten08a.pdf}{
#' Visualizing High-Dimensional Data Using t-SNE}. \emph{Journal of Machine
#' Learning Research}, \emph{9}: 2579-2605.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 5), nrow = 1000, ncol = 5)
#' plot_tsne(mat)
#'
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet()
#' dds <- rlog(dds)
#' plot_tsne(dds, group = colData(dds)$condition)
#'
#' @seealso
#' \code{\link[Rtsne]{Rtsne}}, \code{\link{plot_pca}}, \code{\link{plot_mds}}
#'
#' @export
#' @importFrom Rtsne Rtsne
#' @import dplyr
#' @import ggplot2
#'

plot_tsne <- function(dat,
                      group = NULL,
                      covar = NULL,
                        top = NULL,
              filter_method = 'pairwise',
                       dims = c(1, 2),
                 perplexity = ncol(dat) / 4,
                      theta = 0.1,
                   max_iter = 1000,
                      label = FALSE,
                  pal_group = 'npg',
                  pal_covar = 'blues',
                      title = NULL,
                     legend = 'right',
                      hover = FALSE,
                         D3 = FALSE, ...) {

  # Preliminaries
  if (ncol(dat) < 3L) {
    stop(paste('dat includes only', ncol(dat), 'samples;',
               'need at least 3 for t-SNE.'))
  }
  if (!(group %>% is.null)) {
    group <- dat %>% format_features(group, var_type = 'Categorical')
    if (length(group) > 2L) {
      stop('Plot can render at most two categorical features.')
    }
    if (length(group) == 2L && !(covar %>% is.null)) {
      stop('Plot can render at most one categorical feature when a continuous ',
           'covariate is also supplied.')
    }
    group_cols <- colorize(pal = pal_group, var_type = 'Categorical',
                           n = length(levels(group[[1L]])))
  }
  if (!(covar %>% is.null)) {
    covar <- dat %>% format_features(covar, var_type = 'Continuous')
    if (length(covar) != 1L) {
      stop('Plot can render at most one continuous feature.')
    }
    covar_cols <- colorize(pal = pal_covar, var_type = 'Continuous')
  }
  if (!(c(group, covar) %>% is.null)) {
    features <- c(covar, group)
    feature_names <- names(features)
    names(features) <- paste0('Feature', seq_along(features))
  } else {
    features <- NULL
  }
  if (length(dims) > 2L && !D3) {
    stop('dims must be of length 2 when D3 = FALSE.')
  } else if (length(dims) > 3L) {
    stop('dims must be a vector of length <= 3.')
  }
  if (label && length(features) == 2L) {
    stop('If label is TRUE, then plot can render at most one phenotypic ',
         'feature.')
  }
  if (title %>% is.null) {
    title <- 't-SNE'
  }
  if (!legend %in% c('right', 'left', 'top', 'bottom',
                     'topright', 'topleft', 'bottomright', 'bottomleft')) {
    stop('legend must be one of "right", "left", "top", "bottom", ',
         '"topright", "topleft", "bottomright", or "bottomleft".')
  }

  # Tidy data
  dat <- matrixize(dat)
  if (rownames(dat) %>% is.null) {
    rownames(dat) <- seq_len(nrow(dat))
  }
  if (colnames(dat) %>% is.null) {
    colnames(dat) <- paste0('Sample', seq_len(ncol(dat)))
  }
  dm <- dist_mat(dat, top, filter_method, dist = 'euclidean') %>% as.dist(.)
  tsne <- Rtsne(dm, perplexity = perplexity, dims = max(dims),
                theta = theta, max_iter = max_iter, check_duplicates = FALSE,
                is_distance = TRUE, ...)
  tsne <- tsne$Y                                           # t-SNE
  df <- data_frame(Sample = colnames(dat))                 # Melt
  if (length(dims) == 2L) {
    df <- df %>% mutate(PC1 = tsne[, min(dims)],
                        PC2 = tsne[, max(dims)])
  } else {
    other <- setdiff(dims, c(min(dims), max(dims)))
    df <- df %>% mutate(PC1 = tsne[, min(dims)],
                        PC2 = tsne[, other],
                        PC3 = tsne[, max(dims)])
  }
  if (!(features %>% is.null)) {
    df <- df %>% cbind(tbl_df(features))
  }

  # Build plot
  if (top %>% is.null) {
    xlab <- paste('t-SNE Dim', min(dims))
    ylab <- paste('t-SNE Dim', max(dims))
  } else {
    xlab <- paste('Leading logFC, t-SNE Dim', min(dims))
    ylab <- paste('Leading logFC, t-SNE Dim', max(dims))
  }
  embed(df, group, covar, group_cols, covar_cols, feature_names,
        label, title, xlab, ylab, legend, hover, D3)

}

# Set seed?

