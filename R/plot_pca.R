#' PCA Plot
#'
#' This function plots a low-dimensional projection of an omic data matrix using
#' principal component analysis.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples. It is strongly recommended that data be
#'   normalized and filtered prior to running PCA. For count data, this means
#'   undergoing some sort of variance stabilizing transformation, such as
#'   \code{\link[edgeR]{cpm} (with \code{log = TRUE}), \link[DESeq2]{vst},
#'   \link[DESeq2]{rlog}}, etc.
#' @param covar Optional vector of length equal to sample size, or up to two such
#'   vectors organized into a list or data frame. If passing two covariates, only one
#'   may be continuous. Supply legend title(s) by passing a named list or data frame.
#' @param top Optional number (if > 1) or proportion (if < 1) of most variable probes
#'   to be used for PCA.
#' @param pcs Vector specifying which principal components to plot. Must be of length
#'   two unless \code{D3 = TRUE}.
#' @param label Label data points by sample name? Defaults to \code{FALSE} unless
#'   \code{covar = NULL}. If \code{TRUE}, then plot can render at most one covariate.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside",
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
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
#' Hotelling, H. (1933). "Analysis of a complex of variables into principal components."
#' \emph{Journal of Educational Psychology}, \emph{24}(6): 414:441.
#' \url{http://psycnet.apa.org/journals/edu/24/6/417/}
#'
#' Pearson, K. (1901). "On lines and planes of closest fit to systems of points in
#' space." \emph{Philosophical Magazine}, \emph{2}(11): 559–572.
#' \url{http://www.tandfonline.com/doi/abs/10.1080/14786440109462720}
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 5), nrow = 1000, ncol = 5)
#' plot_pca(mat)
#'
#' library(DESeq2)
#' mat <- cbind(matrix(rnbinom(5000, mu = 4, size = 1), nrow = 1000, ncol = 5),
#'              matrix(rnbinom(5000, mu = 4, size = 10), nrow = 1000, ncol = 5))
#' mat <- rlog(mat)
#' grp <- gl(n = 2, k = 5, labels = c("A", "B"))
#' plot_pca(mat, covar = grp)
#'
#' @seealso
#' \code{\link{plot_mds}}, \code{\link[DESeq2]{plotPCA}}
#'
#' @export
#' @importFrom limma getEAWP
#' @import dplyr
#' @importFrom purrr map map_chr
#' @importFrom matrixStats rowVars
#' @import ggplot2
#' @import plotly
#' @importFrom scales hue_pal
#'

plot_pca <- function(dat,
                     covar = NULL,
                       top = NULL,
                       pcs = c(1, 2),
                     label = FALSE,
                      main = NULL,
                    legend = 'outside',
                     hover = FALSE,
                        D3 = FALSE) {

  # Preliminaries
  if (ncol(dat) < 3L) {
    stop(paste('dat includes only', ncol(dat), 'samples; need at least 3 for PCA.'))
  }
  dat <- getEAWP(dat)$expr
  keep <- rowSums(is.finite(dat)) == ncol(dat)
  dat <- dat[keep, , drop = FALSE]
  if (!is.null(covar)) {
    if (is.data.frame(covar)) covar <- as.list(covar)
    else if (!is.list(covar)) covar <- list(covar)
    if (length(covar) > 2L) {
      stop('covar cannot contain more than two covariates.')
    }
    for (i in seq_along(covar)) {
      if (length(covar[[i]]) != ncol(dat)) {
        stop('Covariate(s) must be of length equal to sample size.')
      }
      if (is.numeric(covar[[i]]) && var(covar[[i]]) == 0L) {
        warning('Continuous feature is invariant.')
      } else if (!is.numeric(covar[[i]]) && length(unique(covar[[i]])) == 1L) {
        warning('Grouping factor is invariant.')
      }
    }
    nums <- as.logical(map(covar, is.numeric))
    if (sum(nums) == 2L) {
      stop('Only one continuous covariate can be plotted at a time.')
    }
    if (any(nums)) {
      cont_cov <- TRUE
      if (which(nums) == 2L) covar <- covar[c(2, 1)]
    } else cont_cov <- FALSE
    if (!is.null(names(covar))) covars <- names(covar)
    else {
      if (cont_cov) covars <- c('Feature', 'Group')
      else {
        if (length(covar) == 1L) covars <- 'Group'
        else covars <- c('Factor 1', 'Factor 2')
      }
    }
    names(covar) <- paste0('Feature', seq_along(covar))
  }
  if (!is.null(top)) {
    if (top > 1L) {
      if (top > nrow(dat)) {
        warning(paste('top exceeds nrow(dat), at least after removing probes with infinite',
                      'or missing values. Proceeding with the complete', nrow(dat), 'x',
                      ncol(dat), 'matrix.'))
      }
    } else top <- round(top * nrow(dat))
    vars <- rowVars(dat)
    dat <- dat[rev(order(vars))[seq_len(top)], , drop = FALSE]
  }
  if (length(pcs) > 2L && !D3) {
    stop('pcs must be of length 2 when D3 = FALSE.')
  } else if (length(pcs) > 3L) {
    stop('pcs must be a vector of length <= 3.')
  }
  if (label && !is.null(covar) && length(covar) == 2L) {
    stop('If label is TRUE, then plot can render at most one covariate.')
  }
  if (is.null(main)) main <- 'PCA'
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
  pca <- prcomp(t(dat))                          # PCA
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
  if (!is.null(covar)) {
    covar <- tbl_df(covar) %>% mutate(Sample = colnames(dat))
    df <- df %>% inner_join(covar, by = 'Sample')
  }

  # Build plot
  if (!D3) {
    p <- ggplot(df, aes(PC1, PC2)) +
      geom_hline(yintercept = 0L, color = 'grey') +
      geom_vline(xintercept = 0L, color = 'grey') +
      labs(title = main, x = pve[min(pcs)], y = pve[max(pcs)]) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    if (is.null(covar)) {
      if (label) p <- p + geom_text(aes(label = Sample), alpha = 0.85)
      else suppressWarnings(p <- p + geom_point(aes(text = Sample)))
    } else if (ncol(covar) == 2L) {
      if (label) {
        p <- p + geom_text(aes(label = Sample, color = Feature1),
                           alpha = 0.85)
      } else {
        if (cont_cov) {
          suppressWarnings(
            p <- p + geom_point(aes(text = Sample, color = Feature1),
                                alpha = 0.85)
          )
        } else {
          suppressWarnings(
            p <- p + geom_point(aes(text = Sample, color = Feature1, shape = Feature1),
                                alpha = 0.85)
          )
        }
      }
      p <- p + guides(color = guide_legend(title = covars[1]),
                      shape = guide_legend(title = covars[1]))
    } else {
      suppressWarnings(
        p <- p + geom_point(aes(text = Sample, color = Feature1, shape = Feature2),
                            alpha = 0.85) +
          guides(color = guide_legend(title = covars[1]),
                 shape = guide_legend(title = covars[2]))
      )
    }
    if (legend == 'bottomleft') {             # Locate legend
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
    if (!hover) print(p)
    else {
      if (legend == 'outside') {
        p <- ggplotly(p, tooltip = 'text', height = 525, width = 600)
      } else {
        p <- ggplotly(p, tooltip = 'text', height = 600, width = 600)
      }
      print(p)
    }
  } else {
    ### REWRITE ###
    # symbls <- c(16, 17, 15, 3, 7, 8)      # This would be right if plotly worked
    symbls <- c(16, 18, 15, 3, 7, 8)
    p <- plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3,
                 text = ~Sample, color = ~Group, symbol = ~Group,
                 colors = hue_pal()(length(unique(df$Group))),
                 symbols = symbls[1:length(unique(df$Group))],
                 type = 'scatter3d', mode = 'markers',
                 alpha = 0.85, hoverinfo = 'text', marker = list(size = 5)) %>%
      layout(hovermode = 'closest', title = main, scene = list(
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
