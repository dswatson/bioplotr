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
#' calculated using the \code{top} most differentially expressed probes between the
#' two samples. This method is appropriate when different molecular pathways are
#' relevant for distinguishing different pairs of samples. To run MDS on the complete
#' data, set \code{top = NULL}. This is functionally equivalent to running PCA
#' on the full matrix. See \code{\link{plot_pca}}.
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
#' @importFrom limma getEAWP
#' @import dplyr
#' @importFrom wordspace dist.matrix
#' @import ggplot2
#' @import plotly
#' @importFrom scales hue_pal
#'

plot_mds <- function(dat,
                     covar = NULL,
                       top = 500,
                     label = FALSE,
                      main = NULL,
                    legend = 'outside',
                     hover = FALSE,
                        D3 = FALSE) {

  # Preliminaries
  dat <- getEAWP(dat)
  dat <- dat$expr
  bad <- rowSums(is.finite(dat)) < ncol(dat)
  if (any(bad)) {
    dat <- dat[!bad, , drop = FALSE]
  }
  if (ncol(dat) < 3) {
    stop(paste('dat includes only', ncol(dat), 'samples; need at least 3 for MDS.'))
  }
  if (is.null(rownames(dat))) {
    rownames(dat) <- seq_len(nrow(dat))
  }
  if (is.null(colnames(dat))) {
    colnames(dat) <- paste0('Sample', seq_len(ncol(dat)))
  }
  if (!is.null(covar)) {
    if (is.data.frame(covar)) {
      covar <- as.list(covar)
    } else if (!is.list(covar)) {
      covar <- list(covar)
    }
    if (length(covar) > 2) {
      stop('covar cannot contain more than two covariates.')
    }
    for (i in seq_along(covar)) {
      if (length(covar[[i]]) != ncol(dat)) {
        stop('Covariate(s) must be of length equal to sample size.')
      }
      if (is.numeric(covar[[i]]) && var(covar[[i]]) == 0) {
        warning('Continuous feature is invariant.')
      } else if (!is.numeric(covar[[i]]) && length(unique(covar[[i]])) == 1) {
        warning('Grouping factor is invariant.')
      }
    }
    nums <- as.logical(map(covar, is.numeric))
    if (sum(nums) == 2) {
      stop('Only one continuous covariate can be plotted at a time.')
    }
    if (any(nums)) {
      cont_cov <- TRUE
      if (which(nums) == 2) {
        covar <- covar[c(2, 1)]
      }
    } else {
      cont_cov <- FALSE
    }
    if (!is.null(names(covar))) {
      covars <- names(covar)
    } else {
      if (cont_cov) {
        covars <- c('Feature', 'Group')
      } else {
        if (length(covar) == 1) {
          covars <- 'Group'
        } else {
          covars <- c('Factor 1', 'Factor 2')
        }
      }
    }
    names(covar) <- paste0('Feature', seq_along(covar))
    covar <- tbl_df(covar) %>% mutate(Sample = colnames(dat))
  }
  if (!is.null(top)) {
    if (top > 1) {
      if (top > nrow(dat)) {
        warning('top exceeds nrow(dat), at least after removing probes with infinite ',
                'or missing values. Proceeding with the complete matrix.')
      }
    } else {
      top <- round(top * nrow(dat))
    }
  } else {
    top <- nrow(dat)
  }
  if (label && length(covars) == 2) {
    stop('If label is TRUE, then plot can render at most one covariate.')
  }
  if (is.null(main)) {
    main <- 'MDS'
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright"')
  }

  # Tidy data
  if (is.null(top)) {                                      # Distance matrix
    dm <- dist.matrix(t(dat), method = 'euclidean')
    dimnames(dm) <- list(colnames(dat), colnames(dat))
  } else {
    dm <- matrix(0, nrow = ncol(dat), ncol = ncol(dat),
                 dimnames = list(colnames(dat), colnames(dat)))
    top_idx <- nrow(dat) - top + 1
    for (i in 2:ncol(dat)) {
      for (j in 1:(i - 1)) {
        if (is.null(top)) {
          dm[i, j] <- sqrt(sum(sort.int((dat[, i] - dat[, j])^2,
                                        partial = top_idx)[top_idx:nrow(dat)]))
        }
      }
    }
  }
  mds <- suppressWarnings(cmdscale(as.dist(dm), k = 3))    # MDS
  df <- data_frame(Sample = colnames(dat),                 # Melt
                   PC1 = mds[, 1],
                   PC2 = mds[, 2],
                   PC3 = mds[, 3])
  if (!is.null(covar)) {
    df <- inner_join(df, covar, by = 'Sample')
  }

  # Build plot
  p <- ggplot(df, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, size = 0.2) +
    geom_vline(xintercept = 0, size = 0.2) +
    labs(title = main, x = 'PC1', y = 'PC2') +
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

  # Output
  if (!D3) {
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
        xaxis = list(title = pve[1]),
        yaxis = list(title = pve[2]),
        zaxis = list(title = pve[3])))
    print(p)
  }

}

# Argument to pick which PCs to plot
