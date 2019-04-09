#' ICA Plot
#'
#' This function plots a low-dimensional projection of an omic data matrix using
#' independent component analysis.
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
#' @param top Optional number (if > 1) or proportion (if < 1) of most variable
#'   probes to be used for ICA.
#' @param dims Vector specifying which independent components to plot. Must be 
#'   of length two unless \code{D3 = TRUE}.
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
#' This function plots the samples of an omic data matrix in a two- or
#' three-dimensional independent component subspace. ICA is an easy and
#' popular projection method that can aid in identifying clusters, detecting
#' outliers, and visualizing the latent structure of a dataset.
#' 
#' ICA is a general method for separating a multivariate signal into additive
#' subcomponents. ICA algorithms differ in their objective functions and 
#' optimization methods. \code{plot_ica} relies on Cardoso's JADE ICA algorithm,
#' as implemented in the \code{\link[JADE]{JADE}} package. 
#' 
#' By default, \code{plot_ica} decomposes the complete \code{dat} matrix. 
#' Limit the ICA to only the most variable probes by using the \code{top} 
#' argument.
#'
#' @references
#' Miettinen, J., Nordhausen, K. & Taskinen, S. (2017). 
#' \href{https://www.jstatsoft.org/article/view/v076i02}{Blind Source Separation 
#' Based on Joint Diagonalization in R: The Packages JADE and BSSasymp}. \emph{
#' Journal of Statistical Software}, \emph{76}: 1–31.
#' 
#' @examples
#' n <- 10L
#' p <- 1000L
#' x <- matrix(rnorm(n * p), ncol = n)
#' plot_ica(x)
#'
#' @seealso
#' \code{\link[JADE]{JADE}}, \code{\link{plot_pca}}, \code{\link{plot_tsne}}
#'
#' @export
#' @importFrom JADE JADE
#' @importFrom Rfast rowmeans
#' @import dplyr
#' @import ggplot2
#'

plot_ica <- function(dat,
                     group = NULL,
                     covar = NULL,
                       top = NULL,
                      dims = c(1L, 2L),
                     label = FALSE,
                 pal_group = 'npg',
                 pal_covar = 'Blues',
                      size = NULL,
                     alpha = NULL,
                     title = 'ICA',
                    legend = 'right',
                     hover = FALSE,
                        D3 = FALSE) {
  
  # Preliminaries
  if (ncol(dat) < 3L) {
    stop('dat includes only ', ncol(dat), ' samples; need at least 3 for ICA.')
  }
  if (!group %>% is.null) {
    group <- dat %>% format_features(group, var_type = 'Categorical')
    if (length(group) > 2L) {
      stop('Plot can render at most two categorical features.')
    }
    if (length(group) == 2L && !covar %>% is.null) {
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
  if (length(dims) > 2L && !D3) {
    stop('dims must be of length 2 when D3 = FALSE.')
  } else if (length(dims) > 3L) {
    stop('dims must be a vector of length <= 3.')
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
  if (!top %>% is.null) {                        # Filter by variance?
    dat <- var_filt(dat, top, robust = FALSE)
  }                                              # Mean center the probes
  dat <- dat - rowmeans(dat)
  ica <- JADE(dat, n.comp = max(dims))           # ICA
  df <- tibble(Sample = colnames(dat))           # Melt
  if (length(dims) == 2L) {
    df <- df %>% mutate(PC1 = ica$A[, min(dims)],
                        PC2 = ica$A[, max(dims)])
  } else {
    other <- setdiff(dims, c(min(dims), max(dims)))
    df <- df %>% mutate(PC1 = ica$A[, min(dims)],
                        PC2 = ica$A[, other],
                        PC3 = ica$A[, max(dims)])
  }
  if (!features %>% is.null) {
    df <- df %>% bind_cols(as_tibble(features))
  }
  
  # Build plot
  xlab <- paste0('IC', min(dims))
  ylab <- paste0('IC', max(dims))
  embed(df, group, covar, group_cols, covar_cols, feature_names,
        label, size, alpha, title, xlab, ylab, legend, hover, D3)
  
}

# Options for other decomposed matrices?


