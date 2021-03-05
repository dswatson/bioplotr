#' LLE Plot
#'
#' This function plots a low-dimensional projection of an omic data matrix using
#' locally linear embedding.
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
#' @param k How many nearest neighbors should the algorithm consider when
#'   building projections? See Details.
#' @param top Optional number (if > 1) or proportion (if < 1) of most variable
#'   probes to be used for LLE.
#' @param dims Vector specifying which principal components to plot. Must be of
#'   length two unless \code{D3 = TRUE}.
#' @param label Label data points by sample name? Defaults to \code{FALSE}
#'   unless \code{group} and \code{covar} are both \code{NULL}. If \code{TRUE},
#'   then plot can render at most one phenotypic feature.
#' @param pal_group String specifying the color palette to use if \code{group}
#'   is non-\code{NULL}, or a vector of such strings with length equal to the
#'   number of vectors passed to \code{group}. Options include \code{"ggplot"}, 
#'   all qualitative color schemes available in \code{\href{
#'   https://bit.ly/2ipuEjn}{RColorBrewer}}, and the complete collection of 
#'   \code{\href{http://bit.ly/2bxnuGB}{ggsci}} palettes. Alternatively, a 
#'   character vector of colors with length equal to the cumulative number of 
#'   levels in \code{group}.
#' @param pal_covar String specifying the color palette to use if \code{covar}
#'   is non-\code{NULL}, or a vector of such strings with length equal to the
#'   number of vectors passed to \code{covar}. Options include the complete 
#'   collection of \code{\href{https://bit.ly/2n7D6tF}{viridis}} palettes, as 
#'   well as all sequential color schemes available in \code{\href{
#'   https://bit.ly/2ipuEjn}{RColorBrewer}}. Alternatively, a character vector 
#'   of colors representing a smooth gradient, or a list of such vectors with 
#'   length equal to the number of continuous variables to visualize.
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
#'
#' @details
#' Locally linear embedding is a feature extraction method that projects a
#' high-dimensional dataset in two or three dimensions. It prioritizes local
#' over global structure, and is designed to capture nonlinearities that would
#' be impossible to detect using PCA or classical MDS.
#'
#' By default, the number of nearest neighbors \emph{k} used in the LLE
#' projection is fixed at one quarter the sample size. This is an arbitrary
#' value that may be unsuitable for some datasets. The parameter may be
#' optimized using Kayo's method (2006). See See \code{lle::\link[lle]{calc_k}}.
#'
#' @references
#' Roweis, S.T. & Saul, L.K. (2000).
#' \href{http://science.sciencemag.org/content/290/5500/2323}{Nonlinear
#' Dimensionality Reduction by Locally Linear Embedding}. \emph{Science,
#' 290}(5500): 2323-2326.
#'
#' Kayo, O. (2006).
#' \href{http://jultika.oulu.fi/files/isbn9514280415.pdf}{Locally Linear
#' Embedding Algorithm: Extensions and Applications} (Unpublished doctoral
#' dissertation). University of Oulu, Finland.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 5), nrow = 1000, ncol = 5)
#' plot_lle(mat)
#'
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet()
#' dds <- rlog(dds)
#' plot_lle(dds, group = colData(dds)$condition)
#'
#' @seealso
#' \code{\link{plot_pca}}, \code{\link{plot_tsne}}
#'
#' @export
#' @importFrom lle lle
#' @import dplyr
#' @import ggplot2
#'

plot_lle <- function(
  dat,
      group = NULL,
      covar = NULL,
          k = ncol(dat) / 4L,
        top = NULL,
       dims = c(1L, 2L),
      label = FALSE,
  pal_group = 'npg',
  pal_covar = 'Blues',
       size = NULL,
      alpha = NULL,
      title = 'LLE',
     legend = 'right',
      hover = FALSE,
         D3 = FALSE
) {

  # Preliminaries
  if (ncol(dat) < 3L) {
    stop('dat includes only ', ncol(dat), ' samples; need at least 3 for LLE.')
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
  if (!covar %>% is.null) {
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
  if (!k %>% is.numeric && k != 'test') {
    stop('k must be numeric or "test".')
  }
  if (k >= ncol(dat)) {
    stop('k must be less than the sample size of the dataset.')
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
  locations <- c('bottom', 'left', 'top', 'right',
                 'bottomright', 'bottomleft', 'topleft', 'topright')
  legend <- match.arg(legend, locations)

  # Tidy data
  dat <- matrixize(dat)
  if (!top %>% is.null) {                        # Filter by variance?
    dat <- var_filt(dat, top, robust = FALSE)
  }
  capture.output(y <- lle(t(dat), m = max(dims), k = k, p = 1L)$Y)        # LLE
  df <- tibble(Sample = colnames(dat))                                    # Melt
  if (length(dims) == 2L) {
    df <- df %>% mutate(PC1 = y[, min(dims)],
                        PC2 = y[, max(dims)])
  } else {
    other <- setdiff(dims, c(min(dims), max(dims)))
    df <- df %>% mutate(PC1 = y[, min(dims)],
                        PC2 = y[, other],
                        PC3 = y[, max(dims)])
  }
  if (!(features %>% is.null)) {
    df <- df %>% bind_cols(as_tibble(features))
  }

  # Build plot
  xlab <- paste('LLE Dim', min(dims))
  ylab <- paste('LLE Dim', max(dims))
  embed(df, group, covar, group_cols, covar_cols, feature_names,
        label, size, alpha, title, xlab, ylab, legend, hover, D3)

}




