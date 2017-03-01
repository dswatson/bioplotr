#' MDS Plot
#'
#' This function plots the principal coordinates of an omic data matrix.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples. It is strongly recommended that data be
#'   normalized and filtered prior to running MDS.
#' @param covar Optional vector of length equal to sample size, or up to two such
#'   vectors organized into a list or data frame. If passing two covariates, only one
#'   may be continuous. Supply legend title(s) by passing a named list or data frame.
#' @param top Optional number (if > 1) or proportion (if < 1) of top probes to be used
#'   for MDS. See Details.
#' @param pcs Vector specifying which principal coordinates to plot. Must be of length
#'   2 unless \code{D3 = TRUE}.
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
#' This function plots the samples of an omic data matrix in a two- or three-
#' dimensional principal coordinate subspace. MDS is an easy and popular method for
#' unsupervised cluster detection. It can also aid in spotting potential outliers
#' and generally helps to visualize the latent structure of a data set.
#'
#' The \code{top} argument filters probes using the leading fold change method of
#' Smyth et al. (See \code{\link[limma]{plotMDS}}. Pairwise Euclidean distances are
#' calculated using the most differentially expressed probes between the two samples.
#' This method is appropriate when different molecular pathways are relevant for
#' distinguishing different pairs of samples. To run MDS on the complete data, set
#' \code{top = NULL}. This is functionally equivalent to running PCA on the full
#' matrix. See \code{\link{plot_pca}}.
#'
#' @examples
#' mat <- matrix(rnorm(5000), nrow = 1000, ncol = 5)
#' plot_mds(mat)
#'
#' library(DESeq2)
#' mat <- cbind(matrix(rnbinom(5000, mu = 4, size = 1), nrow = 1000, ncol = 5),
#'              matrix(rnbinom(5000, mu = 4, size = 10), nrow = 1000, ncol = 5))
#' mat <- rlog(mat)
#' grp <- gl(n = 2, k = 5, labels = c("A", "B"))
#' plot_pca(mat, covar = grp)
#'
#' @seealso
#' \code{\link{plot_pca}} \code{\link[limma]{plotMDS}}
#'
#' @export
#' @import dplyr
#' @importFrom purrr map
#' @importFrom limma getEAWP
#' @importFrom wordspace dist.matrix
#' @import ggplot2
#' @import plotly
#' @importFrom scales hue_pal
#'

plot_mds <- function(dat,
                     covar = NULL,
                       top = 500,
                       pcs = c(1, 2),
                     label = FALSE,
                      main = NULL,
                    legend = 'outside',
                     hover = FALSE,
                        D3 = FALSE) {

  # Preliminaries
  if (ncol(dat) < 3L) {
    stop(paste('dat includes only', ncol(dat), 'samples; need at least 3 for MDS.'))
  }
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
      if (which(nums) == 2L) covar <- covar[c(2L, 1L)]
      else cont_cov <- FALSE
    }
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
        warning('top exceeds nrow(dat), at least after removing probes with infinite ',
                'or missing values. Proceeding with the complete matrix.')
      }
    } else top <- round(top * nrow(dat))
  }
  if (length(pcs) > 2L && !D3) {
    stop('pcs must be of length 2 when D3 = FALSE.')
  } else if (length(pcs) > 3L) {
    stop('pcs must be a vector of length <= 3.')
  }
  if (label && length(covars) == 2L) {
    stop('If label is TRUE, then plot can render at most one covariate.')
  }
  if (is.null(main)) main <- 'MDS'
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright"')
  }

  # Tidy data
  dat <- getEAWP(dat)$expr
  keep <- rowSums(is.finite(dat)) == ncol(dat)
  dat <- dat[keep, , drop = FALSE]
  if (is.null(rownames(dat))) {
    rownames(dat) <- seq_len(nrow(dat))
  }
  if (is.null(colnames(dat))) {
    colnames(dat) <- paste0('Sample', seq_len(ncol(dat)))
  }
  if (is.null(top)) {                                            # Distance matrix
    dm <- dist.matrix(t(dat), method = 'euclidean')
    dimnames(dm) <- list(colnames(dat), colnames(dat))
  } else {
    dm <- matrix(0, nrow = ncol(dat), ncol = ncol(dat),
                 dimnames = list(colnames(dat), colnames(dat)))
    top_idx <- nrow(dat) - top + 1L
    for (i in 2L:ncol(dat)) {
      for (j in 1L:(i - 1L)) {
        dm[i, j] <- sqrt(sum(sort.int((dat[, i] - dat[, j])^2L,
                                      partial = top_idx)[top_idx:nrow(dat)]))
      }
    }
  }
  mds <- suppressWarnings(cmdscale(as.dist(dm), k = max(pcs)))    # MDS
  if (length(pcs) == 2L) {                                        # Melt
    df <- df %>% mutate(PC1 = mds[, min(pcs)],
                        PC2 = mds[, max(pcs)])
  } else {
    other <- setdiff(pcs, c(min(pcs), max(pcs)))
    df <- df %>% mutate(PC1 = mds[, min(pcs)],
                        PC2 = mds[, other],
                        PC3 = mds[, max(pcs)])
  }
  if (!is.null(covar)) {
    covar <- tbl_df(covar) %>% mutate(Sample = colnames(dat))
    df <- inner_join(df, covar, by = 'Sample')
  }

  # Build plot
  if (!D3) {
    p <- ggplot(df, aes(PC1, PC2)) +
      geom_hline(yintercept = 0L, size = 0.2) +
      geom_vline(xintercept = 0L, size = 0.2) +
      labs(title = main, x = paste('PC', min(pcs)), y = paste0('PC', max(pcs))) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    if (ncol(covar) == 2) {
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
      p <- p + guides(color = guide_legend(title = covars[1L]),
                      shape = guide_legend(title = covars[1L]))
    } else {
      suppressWarnings(
        p <- p + geom_point(aes(text = Sample, color = Feature1, shape = Feature2),
                            alpha = 0.85) +
          guides(color = guide_legend(title = covars[1L]),
                 shape = guide_legend(title = covars[2L]))
      )
    }
    if (legend == 'bottomleft') {                            # Locate legend
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

# Some way to check that covar vector is ordered the same as cols of dat?
