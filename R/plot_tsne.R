#' t-SNE Plot
#'
#' This function plots a low-dimensional projection of an omic data matrix using
#' t-distributed stochastic neighbor embedding.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples. It is strongly recommended that data be
#'   normalized and filtered prior to running t-SNE. For count data, this means
#'   undergoing some sort of variance stabilizing transformation, such as
#'   \code{\link{lcpm}, \link[DESeq2]{vst}, \link[DESeq2]{rlog}}, etc.
#' @param group Optional character or factor vector of length equal to sample size,
#'   or up to two such vectors organized into a list or data frame. Supply legend
#'   title(s) by passing a named list or data frame.
#' @param covar Optional continuous covariate. If non-\code{NULL}, then function can
#'   take at most one \code{group} variable. Supply legend title by passing a named
#'   list or data frame.
#' @param top Optional number (if > 1) or proportion (if < 1) of top probes to be used
#'   for t-SNE. See Details.
#' @param dims Vector specifying which dimensions to plot. Must be of length
#'   two unless \code{D3 = TRUE}.
#' @param perplexity How many nearest neighbors should the algorithm consider when
#'   building projections?
#' @param theta Speed/accuracy tradeoff of the Barnes-Hut algorithm. See Details.
#' @param max_iter Maximum number of iterations over which to minimize the loss
#'   function. See Details.
#' @param label Label data points by sample name? Defaults to \code{FALSE} unless
#'   \code{covar = NULL}. If \code{TRUE}, then plot can render at most one covariate.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside",
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
#' @param hover Show sample name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#' @param D3 Render plot in three dimensions?
#' @param ... Additional arguments to be passed to \code{\link[Rtsne]{Rtsne}}.
#'
#' @details
#' This function plots the samples of an omic data matrix in a two- or
#' three-dimensional stochastic neighbor subspace. t-SNE is a popular machine learning
#' method for unsupervised cluster detection. It can also aid in spotting potential
#' outliers, and generally helps to visualize the latent structure of a data set.
#'
#' The \code{top} argument optionally filters probes using the leading fold change
#' method of Smyth et al. See \code{\link{plot_mds}} for more details.
#'
#' \code{plot_tsne} relies on a C++ implementation of the Barnes-Hut algorithm, which
#' vastly accelerates the original t-SNE projection method. An exact t-SNE plot may
#' be rendered by setting \code{theta = 0}. Briefly, the algorithm computes
#' samplewise similarities based on distances in the original \emph{p}-dimensional
#' space (where \emph{p} = the number of probes); generates a low-dimensional
#' embedding of the samples based on the user-defined \code{perplexity} parameter;
#' and iteratively minimizes the Kullback-Leibler divergence between these two
#' distributions using an efficient tree search. See \code{\link[Rtsne]{Rtsne}} for
#' more details. A thorough introduction to and explication of the original t-SNE
#' method and the Barnes-Hut approximation may be found in the references below.
#'
#' @references
#' van der Maaten, L.J.P. (2014).
#' \href{http://www.jmlr.org/papers/volume15/vandermaaten14a/source/vandermaaten14a.pdf}{Accelerating
#' t-SNE using Tree-Based Algorithms}. \emph{Journal of Machine Learning Research},
#' \emph{15}: 3221-3245.
#'
#' van der Maaten, L.J.P. & Hinton, G.E. (2008).
#' \href{http://www.jmlr.org/papers/volume9/vandermaaten08a/vandermaaten08a.pdf}{Visualizing
#' High-Dimensional Data Using t-SNE}. \emph{Journal of Machine Learning Research},
#' \emph{9}: 2579-2605.
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
#' \code{\link[Rtsne]{Rtsne}, \link{plot_pca}, \link{plot_mds}}
#'
#' @export
#' @import dplyr
#' @importFrom purrr map
#' @importFrom limma getEAWP
#' @importFrom wordspace dist.matrix
#' @importFrom Rtsne Rtsne
#' @importFrom ggsci scale_color_d3 pal_d3
#' @import ggplot2
#' @import plotly
#'

plot_tsne <- function(dat,
                      group = NULL,
                      covar = NULL,
                        top = NULL,
                       dims = c(1, 2),
                 perplexity = ncol(dat) / 4,
                      theta = 0.1,
                   max_iter = 1000,
                      label = FALSE,
                      title = NULL,
                     legend = 'outside',
                      hover = FALSE,
                         D3 = FALSE,
                        ...) {

  # Preliminaries
  if (ncol(dat) < 3L) {
    stop(paste('dat includes only', ncol(dat), 'samples; need at least 3 for t-SNE.'))
  }
  dat <- getEAWP(dat)$expr
  keep <- rowSums(is.finite(dat)) == ncol(dat)
  dat <- dat[keep, , drop = FALSE]
  if (!is.null(group)) {
    if (is.data.frame(group)) {
      group <- as.list(group)
    } else if (!is.list(group)) {
      group <- list(group)
    }
    if (length(group) > 2L) {
      stop('Plot can render at most two categorical variables.')
    }
    if (length(group) == 2L) {
      if (!is.null(covar)) {
        stop('Plot can render at most one categorical variable when a continuous ',
             'covariate is also supplied.')
      }
      if (is.null(names(group))) {
        names(group) <- paste('Factor', seq_len(2))
      }
    } else {
      if (is.null(names(group))) {
        names(group) <- 'Group'
      }
    }
    for (i in seq_along(group)) {
      if (length(group[[i]]) != ncol(dat)) {
        stop(paste(names(group)[i], 'not equal to sample size.'))
      }
      if (length(unique(group[[i]])) == 1L) {
        warning(paste(names(group)[i], 'is invariant.'))
      }
      if (!(is.character(group[[i]]) || is.factor(group[[i]]))) {
        warning(paste0(names(group)[i], ' is of class ', class(group[[i]]), '; ',
                       'coercing to factor.'))
        group[[i]] <- as.factor(group[[i]])
      }
    }
  } else if (!is.null(covar)) {
    if (is.data.frame(covar)) {
      covar <- as.list(covar)
    } else if (!is.list(covar)) {
      covar <- list(covar)
    }
    if (length(covar) != 1L) {
      stop('Plot can render at most one continuous covariate.')
    }
    if (length(covar[[1]] != ncol(dat))) {
      stop('covar must be of length equal to sample size.')
    }
    if (!is.numeric(covar[[1]])) {
      stop('covar must be of class numeric.')
    }
    if (var(covar[[1]]) == 0L) {
      warning('covar is invariant.')
    }
    if (is.null(names(covar))) {
      names(covar) <- 'Feature'
    }
  } else {
    features <- NULL
  }
  if (!is.null(c(group, covar))) {
    features <- c(group, covar)
    feature_names <- names(features)
    names(features) <- paste0('Feature', seq_along(features))
  }
  if (!is.null(top)) {
    if (top > 1L) {
      if (top > nrow(dat)) {
        warning(paste('top exceeds nrow(dat), at least after removing probes with
                      infinite or missing values. Proceeding with the complete',
                      nrow(dat), 'x', ncol(dat), 'matrix.'))
        top <- NULL
      }
    } else {
      top <- round(top * nrow(dat))
    }
  }
  if (length(dims) > 2L && !D3) {
    stop('dims must be of length 2 when D3 = FALSE.')
  } else if (length(dims) > 3L) {
    stop('dims must be a vector of length <= 3.')
  }
  if (label && !is.null(features) && length(features) == 2L) {
    stop('If label is TRUE, then plot can render at most one phenotypic feature.')
  }
  if (is.null(title)) {
    title <- 't-SNE'
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright"')
  }

  # Tidy data
  set.seed(123)
  if (is.null(rownames(dat))) {
    rownames(dat) <- seq_len(nrow(dat))
  }
  if (is.null(colnames(dat))) {
    colnames(dat) <- paste0('Sample', seq_len(ncol(dat)))
  }
  if (is.null(top)) {                                      # Distance matrix
    dm <- dist.matrix(t(dat), method = 'euclidean')
  } else {
    dm <- matrix(nrow = ncol(dat), ncol = ncol(dat))
    top_idx <- nrow(dat) - top + 1L
    for (i in 2L:ncol(dat)) {
      for (j in 1L:(i - 1L)) {
        dm[i, j] <- sqrt(sum(sort.int((dat[, i] - dat[, j])^2L,
                                      partial = top_idx)[top_idx:nrow(dat)]))
      }
    }
  }
  tsne <- Rtsne(as.dist(dm), perplexity = perplexity, dims = max(dims), theta = theta,
                max_iter = max_iter, check_duplicates = FALSE, is_distance = TRUE, ...)
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
  if (!is.null(features)) {
    df <- df %>% cbind(tbl_df(features))
  }

  # Build plot
  sample <- sample_ptsize(df)
  alpha <- sample_alpha(df)
  if (!D3) {
    p <- ggplot(df, aes(PC1, PC2)) +
      geom_hline(yintercept = 0L, color = 'grey') +
      geom_vline(xintercept = 0L, color = 'grey') +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    if (is.null(top)) {
      p <- p + labs(title = title,
                        x = paste('Dim', min(dims)),
                        y = paste('Dim', max(dims)))
    } else {
      p <- p + labs(title = title,
                        x = paste('Leading logFC Dim', min(dims)),
                        y = paste('Leading logFC Dim', max(dims)))
    }
    if (is.null(features)) {
      if (label) {
        p <- p + geom_text(aes(label = Sample), alpha = alpha)
      } else {
        suppressWarnings(
          p <- p + geom_point(aes(text = Sample), size = size, alpha = alpha)
        )
      }
    } else if (length(features) == 1L) {
      if (label) {
        p <- p + geom_text(aes(label = Sample, color = Feature1),
                           alpha = alpha)
      } else {
        if (!is.null(covar)) {
          p <- p + geom_point(aes(text = Sample, color = Feature1),
                              size = size, alpha = alpha)
        } else {
          p <- p + geom_point(aes(text = Sample, color = Feature1, shape = Feature1),
                              size = size, alpha = alpha)
        }
      }
      p <- p + guides(color = guide_legend(title = feature_names[1]),
                      shape = guide_legend(title = feature_names[1]))
    } else {
      p <- p + geom_point(aes(text = Sample, color = Feature1, shape = Feature2),
                          size = size, alpha = alpha) +
        guides(color = guide_legend(title = feature_names[1]),
               shape = guide_legend(title = feature_names[2]))
    }
    p <- p + scale_color_d3()
    if (legend == 'bottomleft') {                          # Locate legend
      p <- p + theme(legend.justification = c(0.01, 0.01),
                     legend.position = c(0.01, 0.01))
    } else if (legend == 'bottomright') {
      p <- p + theme(legend.justification = c(0.99, 0.01),
                     legend.position = c(0.99, 0.01))
    } else if (legend == 'topleft') {
      p <- p + theme(legend.justification = c(0.01, 0.99),
                     legend.position = c(0.01, 0.99))
    } else if (legend == 'topright') {
      p <- p + theme(legend.justification = c(0.99, 0.99),
                     legend.position = c(0.99, 0.99))
    }
    if (!hover) {
      print(p)
    } else {
      if (legend == 'outside') {
        p <- ggplotly(p, tooltip = 'text', height = 525, width = 600)
      } else {
        p <- ggplotly(p, tooltip = 'text', height = 600, width = 600)
      }
      print(p)
    }
  } else {
    # ???
  }

}


# Use gganimate, tweenr, and shiny to:
# 1) filter probes
# 2) filter samples
# 3) change PCs
# 4) tweak perplexity, theta

