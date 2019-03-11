#' UMAP Plot
#'
#' This function plots a low-dimensional projection of an omic data matrix using
#' the uniform manifold approximation and projection algorithm. 
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
#' @param dist Distance measure to be used. Supports all methods available in
#'   \code{\link[stats]{dist}}, \code{Rfast::\link[Rfast]{Dist}}, and \code{
#'   \link[vegan]{vegdist}}, as well as those implemented in the \code{bioDist} 
#'   package. See Details.
#' @param p Power of the Minkowski distance.
#' @param top Optional number (if > 1) or proportion (if < 1) of top probes to
#'   be used for UMAP.
#' @param filter_method String specifying whether to apply a \code{"pairwise"}
#'   or \code{"common"} filter if \code{top} is non-\code{NULL}.
#' @param dims Vector specifying which dimensions to plot. Must be of length
#'   two unless \code{D3 = TRUE}.
#' @param n_neighbors Number of nearest neighbors for the UMAP algorithm. See
#'   \code{\link[umap]{umap.defaults}} for details.
#' @param label Label data points by sample name? Defaults to \code{FALSE}
#'   unless \code{group} and \code{covar} are both \code{NULL}. If \code{TRUE},
#'   then plot can render at most one phenotypic feature.
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
#' @param size Point size. 
#' @param alpha Point transparency.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
#' @param hover Show sample name by hovering mouse over data point? If \code{
#'   TRUE}, the plot is rendered in HTML and will either open in your browser's
#'   graphic display or appear in the RStudio viewer.
#' @param D3 Render plot in three dimensions?
#' @param ... Additional arguments to be passed to \code{\link[umap]{umap}}.
#'
#' @details
#' UMAP is a novel machine learning method for visualizing high-dimensional
#' datasets. It is designed to preserve local structure and aids in revealing
#' unsupervised clusters. UMAP is constructed from a theoretical framework based 
#' on Riemannian geometry and algebraic topology. See McInnes et al., 2018 for
#' details.
#'
#' The \code{umap} function can operate directly on a distance matrix.
#' Available distance measures include: \code{"euclidean"}, \code{"maximum"},
#' \code{"manhattan"}, \code{"canberra"}, \code{"minkowski"}, \code{"cosine"},
#' \code{"pearson"}, \code{"kendall"}, \code{"spearman"}, \code{"bray"}, \code{
#' "kulczynski"}, \code{"jaccard"}, \code{"gower"}, \code{"altGower"}, \code{
#' "morisita"}, \code{"horn"}, \code{"mountford"}, \code{"raup"}, \code{
#' "binomial"}, \code{"chao"}, \code{"cao"}, \code{"mahalanobis"}, \code{"MI"},
#' or \code{"KLD"}. Some distance measures are unsuitable for certain types of
#' data. See \code{\link{dist_mat}} for more details on these methods and links
#' to documentation on each.
#'
#' The \code{top} argument optionally filters data using either probewise
#' variance (if \code{filter_method = "common"}) or the leading fold change
#' method of Smyth et al. (if \code{filter_method = "pairwise"}). See \code{
#' \link{plot_mds}} for more details.
#'
#' @references
#' McInnes, L., Healy, J. & Melville, J. (2018).
#' \href{https://arxiv.org/abs/1802.03426}{UMAP: Uniform Manifold Approximation 
#' and Projection for Dimension Reduction}.
#' \emph{arXiv preprint: 1802.03426v2}.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 20), nrow = 1000, ncol = 20)
#' plot_umap(mat)
#'
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet(m = 20)
#' dds <- rlog(dds)
#' plot_umap(dds, group = colData(dds)$condition)
#'
#' @seealso
#' \code{\link[umap]{umap}}, \code{\link{plot_pca}}, \code{\link{plot_mds}}
#'
#' @export
#' @importFrom umap umap
#' @importFrom umap umap.defaults
#' @import dplyr
#' @import ggplot2
#'

plot_umap <- function(dat,
                      group = NULL,
                      covar = NULL,
                       dist = 'euclidean',
                          p = 2L,
                        top = NULL,
              filter_method = 'pairwise',
                       dims = c(1L, 2L),
                n_neighbors = 15L,
                      label = FALSE,
                  pal_group = 'npg',
                  pal_covar = 'Blues',
                       size = NULL,
                      alpha = NULL,
                      title = 'UMAP',
                     legend = 'right',
                      hover = FALSE,
                         D3 = FALSE, ...) {

  # Preliminaries
  if (!dat %>% is('dist')) {
    n <- ncol(dat)
    if (n < 3L) {
      stop('dat includes only ', ncol(dat), ' samples; ',
           'need at least 3 for UMAP.')
    }
    d <- c('euclidean', 'maximum', 'manhattan', 'canberra', 'minkowski',
           'bhattacharyya', 'hellinger', 'kullback_leibler', 'cosine', 
           'bray', 'kulczynski', 'jaccard', 'gower', 'altGower', 'morisita', 
           'horn', 'mountford', 'raup' , 'binomial', 'chao', 'cao',
           'mahalanobis', 'pearson', 'kendall', 'spearman', 'MI')
    if (!dist %in% d) {
      stop('dist must be one of ', stringify(d, 'or'), '.')
    }
    if (!filter_method %in% c('pairwise', 'common')) {
      stop('filter_method must be either "pairwise" or "common".')
    }
  }
  if (!group %>% is.null) {
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
  } else {
    group_cols <- NULL
  }
  if (!(covar %>% is.null)) {
    covar <- dat %>% format_features(covar, var_type = 'Continuous')
    if (length(covar) != 1L) {
      stop('Plot can render at most one continuous feature.')
    }
    covar_cols <- colorize(pal = pal_covar, var_type = 'Continuous')
  } else {
    covar_cols <- NULL
  }
  if (!c(group, covar) %>% is.null) {
    features <- c(covar, group)
    feature_names <- names(features)
    names(features) <- paste0('Feature', seq_along(features))
  } else {
    features <- feature_names <- NULL
  }
  if (length(dims) > 2L && !D3) {
    stop('dims must be of length 2 when D3 = FALSE.')
  } else if (length(dims) > 3L) {
    stop('dims must be a vector of length <= 3.')
  }
  if (n_neighbors > n) {
    n_neighbors <- n - 1L
    warning('Number of neighbors must be smaller than sample size n. ',
            'Resetting n_neighbors to n - 1.')
  }
  if (label && length(features) == 2L) {
    stop('If label is TRUE, then plot can render at most one phenotypic ',
         'feature.')
  }
  loc <- c('bottom', 'left', 'top', 'right',
           'bottomright', 'bottomleft', 'topleft', 'topright')
  if (!legend %in% loc) {
    stop('legend must be one of ', stringify(loc, 'or'), '.')
  }

  # Tidy data
  dat <- matrixize(dat)
  config <- list(...)
  custom.config <- umap.defaults
  for (x in seq_along(config)) {
    name_x <- names(config)[x]
    custom.config[[name_x]] <- config[[x]]
  }
  custom.config$n_neighbors <- n_neighbors
  proj <- umap(t(dat), custom.config)                      # UMAP
  df <- tibble(Sample = colnames(dat))                     # Melt
  if (length(dims) == 2L) {
    df <- df %>% mutate(PC1 = proj$layout[, min(dims)],
                        PC2 = proj$layout[, max(dims)])
  } else {
    other <- setdiff(dims, c(min(dims), max(dims)))
    df <- df %>% mutate(PC1 = proj$layout[, min(dims)],
                        PC2 = proj$layout[, other],
                        PC3 = proj$layout[, max(dims)])
  }
  if (!features %>% is.null) {
    df <- df %>% bind_cols(as_tibble(features))
  }

  # Build plot
  xlab <- paste('UMAP Dim', min(dims))
  ylab <- paste('UMAP Dim', max(dims))
  embed(df, group, covar, group_cols, covar_cols, feature_names,
        label, size, alpha, title, xlab, ylab, legend, hover, D3)

}

# Extend to other distance functions
