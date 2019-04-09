#' Kernel PCA Plot
#'
#' This function plots a low-dimensional projection of an omic data matrix using
#' kernel principal component analysis.
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
#' @param kernel The kernel generating function, if using KPCA. Options include
#'   \code{"rbfdot"}, \code{"polydot"}, \code{"tanhdot"}, \code{"vanilladot"}, 
#'   \code{"laplacedot"}, \code{"besseldot"}, \code{"anovadot"}, and 
#'   \code{"splinedot"}. To run normal PCA, set to \code{NULL}. See Details.
#' @param kpar A named list of arguments setting parameters for the kernel
#'   function. If \code{NULL}, \code{kernlab} defaults are used. See Details.
#' @param top Optional number (if > 1) or proportion (if < 1) of most variable
#'   probes to be used for KPCA.
#' @param dims Vector specifying which dimensions to plot. Must be of length
#'   two unless \code{D3 = TRUE}.
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
#' This function maps an omic data matrix into a Hilbert space with a
#' user-selected kernel function. The resulting kernel matrix is then projected
#' into a principal component subspace, which is visualized.
#'
#' Kernel PCA is a nonlinear variant of traditional PCA. Kernel methods are
#' designed to uncover subtle structures in complex datasets, and are the basis
#' of many supervised and unsupervised learning techniques, including support
#' vector machines and spectral clustering. For more details on kernel functions
#' and their input parameters, see \code{kernlab::\link[kernlab]{dots}}.
#'
#' By default, \code{plot_kpca} maps the entire \code{dat} matrix into a
#' corresponding Hilbert space. Limit the analysis to only the most variable
#' probes by using the \code{top} argument.
#'
#' @references
#' Schölkopf, B., Smola, A. & Müller, K. (1998).
#' \href{https://www.mitpressjournals.org/doi/abs/10.1162/089976698300017467}{
#' Nonlinear Component Analysis as a Kernel Eigenvalue Problem}. \emph{Neural
#' Computation}, \emph{10}(5), 1299-1319.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 5), nrow = 1000, ncol = 5)
#' plot_kpca(mat)
#'
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet()
#' dds <- rlog(dds)
#' plot_kpca(dds, group = colData(dds)$condition)
#'
#' @seealso
#' \code{\link{plot_pca}}, \code{kernlab::\link[kernlab]{dots}}
#'
#' @export
#' @importFrom purrr map_chr
#' @importFrom kernlab rbfdot
#' @importFrom kernlab polydot
#' @importFrom kernlab tanhdot
#' @importFrom kernlab vanilladot
#' @importFrom kernlab laplacedot
#' @importFrom kernlab besseldot
#' @importFrom kernlab anovadot
#' @importFrom kernlab splinedot
#' @importFrom kernlab kernelMatrix
#' @importFrom kernlab kpca
#' @importFrom kernlab eig
#' @importFrom kernlab rotated
#' @import dplyr
#' @import ggplot2
#'

plot_kpca <- function(dat,
                      group = NULL,
                      covar = NULL,
                     kernel = 'rbfdot',
                       kpar = list('sigma' = 1e-4),
                        top = NULL,
                       dims = c(1L, 2L),
                      label = FALSE,
                  pal_group = 'npg',
                  pal_covar = 'Blues',
                       size = NULL,
                      alpha = NULL,
                      title = 'Kernel PCA',
                     legend = 'right',
                      hover = FALSE,
                         D3 = FALSE) {

  # Preliminaries
  if (ncol(dat) < 3L) {
    stop('dat includes only ', ncol(dat), ' samples; need at least 3 for KPCA.')
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
  if (!(c(group, covar) %>% is.null)) {
    features <- c(covar, group)
    feature_names <- names(features)
    names(features) <- paste0('Feature', seq_along(features))
  } else {
    features <- feature_names <- NULL
  }
  kernels <- c('rbfdot', 'polydot', 'tanhdot', 'vanilladot', 'laplacedot', 
               'besseldot', 'anovadot', 'splinedot')
  if (!kernel %in% kernels) {
    stop('kernel must be one of ', stringify(kernels, 'or'), '. ', 
         'For more info, see ?plot_kpca or ?kernlab::dots.')
  }
  if (length(dims) > 2L & !D3) {
    stop('dims must be of length 2 when D3 = FALSE.')
  } else if (length(dims) > 3L) {
    stop('dims must be a vector of length <= 3.')
  }
  if (label & length(features) == 2L) {
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
  }
  if (kernel == 'rbfdot') {                      # Initialize kernel function
    if (kpar %>% is.null) {
      kpar <- list(sigma = 1L)
    }
    kf <- rbfdot(unlist(kpar))
  } else if (kernel == 'polydot') {
    if (kpar %>% is.null) {
      kpar <- list(degree = 1L, scale = 1L, offset = 1L)
    }
    kf <- polydot(unlist(kpar))
  } else if (kernel == 'tanhdot') {
    if (kpar %>% is.null) {
      kpar <- list(scale = 1L, offset = 1L)
    }
    kf <- tanhdot(unlist(kpar))
  } else if (kernel == 'vanilladot') {
    kf <- vanilladot()
  } else if (kernel == 'laplacedot') {
    if (kpar %>% is.null) {
      kpar <- list(sigma = 1L)
    }
    kf <- laplacedot(unlist(kpar))
  } else if (kernel == 'besseldot') {
    if (kpar %>% is.null) {
      kpar <- list(sigma = 1L, order = 1L, degree = 1L)
    }
    kf <- besseldot(unlist(kpar))
  } else if (kernel == 'anovadot') {
    if (kpar %>% is.null) {
      kpar <- list(sigma = 1L, degree = 1L)
    }
    kf <- anovadot(unlist(kpar))
  } else if (kernel == 'splinedot') {
    kf <- splinedot()
  }
  k_mat <- kernelMatrix(kernel = kf, x = t(dat))
  pca <- kpca(k_mat)                             # PCA, % variance explained
  pve <- seq_len(max(dims)) %>% map_chr(function(pc) {
    p <- as.numeric(eig(pca)[pc] / sum(eig(pca)) * 100L)
    paste0('KPC', pc, ' (', round(p, 2L), '%)')
  })
  df <- tibble(Sample = colnames(dat))           # Melt
  if (length(dims) == 2L) {
    df <- df %>% mutate(PC1 = rotated(pca)[, min(dims)],
                        PC2 = rotated(pca)[, max(dims)])
  } else {
    other <- setdiff(dims, c(min(dims), max(dims)))
    df <- df %>% mutate(PC1 = rotated(pca)[, min(dims)],
                        PC2 = rotated(pca)[, other],
                        PC3 = rotated(pca)[, max(dims)])
  }
  if (!(features %>% is.null)) {
    df <- df %>% bind_cols(as_tibble(features))
  }

  # Build plot
  xlab <- pve[min(dims)]
  ylab <- pve[max(dims)]
  embed(df, group, covar, group_cols, covar_cols, feature_names,
        label, size, alpha, title, xlab, ylab, legend, hover, D3)

}

# Extend filter_method argument?

