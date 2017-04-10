#' PCA Plot
#'
#' This function plots a low-dimensional projection of an omic data matrix using
#' principal component analysis.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to probes
#'   and columns to samples. It is strongly recommended that data be filtered and
#'   normalized prior to running PCA. For count data, this means undergoing some sort
#'   of variance stabilizing transformation, such as\code{\link[edgeR]{cpm}} (with
#'   \code{log = TRUE}), \code{\link[DESeq2]{vst}, \link[DESeq2]{rlog}}, etc. Count
#'   matrices stored in \code{\link[edgeR]{DGEList}} or \code{\link[DESeq2]{
#'   DESeqDataSet}} objects are automatically extracted and transformed to the
#'   log2-CPM scale, with a warning.
#' @param group Optional character or factor vector of length equal to sample size,
#'   or up to two such vectors organized into a list or data frame. Supply legend
#'   title(s) by passing a named list or data frame.
#' @param covar Optional continuous covariate. If non-\code{NULL}, then function can
#'   take at most one \code{group} variable. Supply legend title by passing a named
#'   list or data frame.
#' @param top Optional number (if > 1) or proportion (if < 1) of most variable probes
#'   to be used for PCA.
#' @param pcs Vector specifying which principal components to plot. Must be of length
#'   two unless \code{D3 = TRUE}.
#' @param label Label data points by sample name? Defaults to \code{FALSE} unless
#'   \code{covar = NULL}. If \code{TRUE}, then plot can render at most one covariate.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside", "bottomleft",
#'   "bottomright", "topleft",} or \code{"topright"}.
#' @param hover Show sample name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#' @param D3 Render plot in three dimensions?
#'
#' @details
#' This function plots the samples of an omic data matrix in a two- or
#' three-dimensional principal component subspace. Axis labels include the percentage
#' of variance explained by each component. PCA is an easy and popular method for
#' unsupervised cluster detection. It can also aid in spotting potential outliers,
#' and generally helps to visualize the latent structure of a data set.
#'
#' By default, \code{plot_pca} performs singular value decomposition on the complete
#' \code{dat} matrix. Limit the PCA to only the most variable probes by using the
#' \code{top} argument.
#'
#' @references
#' Hotelling, H. (1933). \href{http://psycnet.apa.org/journals/edu/24/6/417/}{Analysis
#' of a complex of variables into principal components}. \emph{Journal of Educational
#' Psychology}, \emph{24}(6): 414:441.
#'
#' Pearson, K. (1901). \href{http://www.tandfonline.com/doi/abs/10.1080/14786440109462720}{On
#' lines and planes of closest fit to systems of points in space}. \emph{Philosophical
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
#' \code{\link[DESeq2]{plotPCA}, \link{plot_mds}, \link{plot_tsne}}
#'
#' @export
#' @importFrom edgeR calcNormFactors cpm
#' @importFrom DESeq2 estimateSizeFactors counts
#' @importFrom SummarizedExperiment assay
#' @importFrom limma getEAWP
#' @importFrom purrr map_chr
#' @importFrom matrixStats rowVars
#' @importFrom ggsci scale_color_d3 pal_d3
#' @import dplyr
#' @import ggplot2
#' @import plotly
#'

plot_pca <- function(dat,
                     group = NULL,
                     covar = NULL,
                       top = NULL,
                       pcs = c(1, 2),
                     label = FALSE,
                     title = NULL,
                    legend = 'outside',
                     hover = FALSE,
                        D3 = FALSE) {

  # Preliminaries
  if (ncol(dat) < 3L) {
    stop(paste('dat includes only', ncol(dat), 'samples; need at least 3 for PCA.'))
  }
  if (is(dat, 'DGEList')) {
    keep <- rowSums(dat$counts) > 1L             # Minimal count filter
    dat <- dat[keep, ]
    if (is.null(dat$samples$norm.factors) |      # Calculate size factors
        all(dat$samples$norm.factors == 1L)) {
      dat <- calcNormFactors(dat)
    }
    dat <- cpm(dat, log = TRUE, prior.count = 1L)
    warning('Transforming raw counts to log2-CPM scale.')
  } else if (is(dat, 'DESeqDataSet')) {
    if (is.null(sizeFactors(dat)) & is.null(normalizationFactors(dat))) {
      dat <- estimateSizeFactors(dat)            # Normalize counts
    }
    dat <- counts(dat, normalized = TRUE)
    keep <- rowMeans(dat) > 0L                   # Minimal count filter
    dat <- dat[keep, , drop = FALSE]
    dat <- cpm(dat, log = TRUE, prior.count = 1L)
    warning('Transforming raw counts to log2-CPM scale.')
  } else if (is(dat, 'DESeqTransform')) {
    dat <- assay(dat)
  } else {
    dat <- getEAWP(dat)$expr
    keep <- rowSums(is.finite(dat)) == ncol(dat)
    dat <- dat[keep, , drop = FALSE]
  }
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
                      missing values and/or applying a minimal expression filter.
                      Proceeding with the complete', nrow(dat), 'x', ncol(dat), 'matrix.'))
      }
    } else {
      top <- round(top * nrow(dat))
    }
    vars <- rowVars(dat)
    keep <- order(vars, decreasing = TRUE)[seq_len(min(top, nrow(dat)))]
    dat <- dat[keep, , drop = FALSE]
  }
  if (length(pcs) > 2L & !D3) {
    stop('pcs must be of length 2 when D3 = FALSE.')
  } else if (length(pcs) > 3L) {
    stop('pcs must be a vector of length <= 3.')
  }
  if (label & !is.null(features) & length(features) == 2L) {
    stop('If label is TRUE, then plot can render at most one phenotypic feature.')
  }
  if (is.null(title)) {
    title <- 'PCA'
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright"')
  }

  # Tidy data
  if (is.null(rownames(dat))) {
    rownames(dat) <- seq_len(nrow(dat))
  }
  if (is.null(colnames(dat))) {
    colnames(dat) <- paste0('Sample', seq_len(ncol(dat)))
  }
  pca <- prcomp(t(dat))                          # PCA, % variance explained
  pve <- map_chr(seq_len(max(pcs)), function(pc) {
    p <- round(pca$sdev[pc]^2L / sum(pca$sdev^2L) * 100L, 2L)
    paste0('PC', pc, ' (', p, '%)')
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
  if (!is.null(features)) {
    df <- df %>% cbind(tbl_df(features))
  }

  # Build plot
  size <- sample_ptsize(df)
  alpha <- sample_alpha(df)
  if (!D3) {
    p <- ggplot(df, aes(PC1, PC2)) +
      geom_hline(yintercept = 0L, color = 'grey') +
      geom_vline(xintercept = 0L, color = 'grey') +
      labs(title = title, x = pve[min(pcs)], y = pve[max(pcs)]) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    if (is.null(features)) {
      if (label) {
        suppressWarnings(
          p <- p + geom_text(aes(label = Sample), alpha = alpha)
        )
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
          suppressWarnings(
            p <- p + geom_point(aes(text = Sample, color = Feature1),
                                size = size, alpha = alpha)
          )
        } else {
          suppressWarnings(
            p <- p + geom_point(aes(text = Sample, color = Feature1, shape = Feature1),
                                size = size, alpha = alpha)
          )
        }
      }
      p <- p + guides(color = guide_legend(title = feature_names[1]),
                      shape = guide_legend(title = feature_names[1]))
    } else {
      suppressWarnings(
        p <- p + geom_point(aes(text = Sample, color = Feature1, shape = Feature2),
                            size = size, alpha = alpha) +
          guides(color = guide_legend(title = feature_names[1]),
                 shape = guide_legend(title = feature_names[2]))
      )
    }
    p <- p + scale_color_d3()
    p <- locate_legend(p, legend)
    gg_out(p, hover, legend)
  } else {
    ### REWRITE ###
    # symbls <- c(16, 17, 15, 3, 7, 8)           # This would be right if plotly worked
    symbls <- c(16, 18, 15, 3, 7, 8)
    p <- plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3,
                 text = ~Sample, color = ~Group, symbol = ~Group,
                 colors = pal_d3()(length(unique(df$Group))),
                 symbols = symbls[1:length(unique(df$Group))],
                 type = 'scatter3d', mode = 'markers',
                 alpha = 0.85, hoverinfo = 'text', marker = list(size = 5)) %>%
      layout(hovermode = 'closest', title = title, scene = list(
        xaxis = list(title = pve[min(pcs)]),
        yaxis = list(title = pve[other]),
        zaxis = list(title = pve[max(pcs)])))
    print(p)
  }

}


# Use gganimate, tweenr, and shiny to:
# 1) filter probes
# 2) filter samples
# 3) change PCs

