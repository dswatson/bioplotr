#' PCA Plot
#'
#' This function plots a low-dimensional projection of an omic data matrix using
#' principal component analysis.
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
#'   probes to be used for PCA.
#' @param pcs Vector specifying which principal components to plot. Must be of
#'   length two unless \code{D3 = TRUE}.
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
#' three-dimensional principal component subspace. Axis labels include the
#' percentage of variance explained by each component. PCA is an easy and
#' popular projection method that can aid in identifying clusters, detecting
#' outliers, and visualizing the latent structure of a dataset.
#'
#' By default, \code{plot_pca} performs singular value decomposition on the
#' complete \code{dat} matrix. Limit the PCA to only the most variable probes by
#' using the \code{top} argument.
#'
#' @references
#' Hotelling, H. (1933).
#' \href{http://psycnet.apa.org/journals/edu/24/6/417/}{Analysis of a complex of
#' variables into principal components}. \emph{Journal of Educational
#' Psychology}, \emph{24}(6): 414:441.
#'
#' Pearson, K. (1901).
#' \href{http://www.tandfonline.com/doi/abs/10.1080/14786440109462720}{On lines
#' and planes of closest fit to systems of points in space}. \emph{Philosophical
#' Magazine}, \emph{2}(11): 559–572.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 5), nrow = 1000, ncol = 5)
#' plot_pca(mat)
#'
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet()
#' dds <- rlog(dds)
#' plot_pca(dds, group = colData(dds)$condition)
#'
#' @seealso
#' \code{\link[DESeq2]{plotPCA}}, \code{\link{plot_mds}},
#' \code{\link{plot_tsne}}
#'
#' @export
#' @importFrom purrr map_chr
#' @importFrom matrixStats rowVars
#' @import dplyr
#' @import ggplot2
#'

plot_pca <- function(dat,
                     group = NULL,
                     covar = NULL,
                       top = NULL,
                       pcs = c(1, 2),
                     label = FALSE,
                 pal_group = 'npg',
                 pal_covar = 'Blues',
                     title = NULL,
                    legend = 'right',
                     hover = FALSE,
                        D3 = FALSE) {

  # Preliminaries
  if (ncol(dat) < 3L) {
    stop('dat includes only ', ncol(dat), ' samples; need at least 3 for PCA.')
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
  if (length(pcs) > 2L && !D3) {
    stop('pcs must be of length 2 when D3 = FALSE.')
  } else if (length(pcs) > 3L) {
    stop('pcs must be a vector of length <= 3.')
  }
  if (label && length(features) == 2L) {
    stop('If label is TRUE, then plot can render at most one phenotypic ',
         'feature.')
  }
  if (title %>% is.null) {
    title <- 'PCA'
  }
  loc <- c('bottom', 'left', 'top', 'right',
           'bottomright', 'bottomleft', 'topleft', 'topright')
  if (!legend %in% loc) {
    stop('legend must be one of ', stringify(loc, 'or'), '.')
  }

  # Tidy data
  dat <- matrixize(dat)
  if (!(top %>% is.null)) {                      # Filter by variance?
    dat <- var_filt(dat, top, robust = FALSE)
  }
  pca <- prcomp(t(dat))                          # PCA, % variance explained
  pve <- seq_len(max(pcs)) %>% map_chr(function(pc) {
    p <- pca$sdev[pc]^2L / sum(pca$sdev^2L) * 100L
    paste0('PC', pc, ' (', round(p, 2L), '%)')
  })
  df <- data_frame(Sample = colnames(dat))       # Melt
  if (length(pcs) == 2L) {
    df <- df %>% mutate(PC1 = pca$x[, min(pcs)],
                        PC2 = pca$x[, max(pcs)])
  } else {
    other <- setdiff(pcs, c(min(pcs), max(pcs)))
    df <- df %>% mutate(PC1 = pca$x[, min(pcs)],
                        PC2 = pca$x[, other],
                        PC3 = pca$x[, max(pcs)])
  }
  if (!(features %>% is.null)) {
    df <- df %>% bind_cols(as_tibble(features))
  }

  # Build plot
  xlab <- pve[min(pcs)]
  ylab <- pve[max(pcs)]
  embed(df, group, covar, group_cols, covar_cols, feature_names,
        label, title, xlab, ylab, legend, hover, D3)

}


