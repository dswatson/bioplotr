#' MDS Plot
#'
#' This function plots a low-dimensional projection of an omic data matrix using
#' multi-dimensional scaling.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples. It is strongly recommended that data be
#'   normalized and filtered prior to running MDS. For count data, this means
#'   undergoing some sort of variance stabilizing transformation, such as
#'   \code{\link{lcpm}, \link[DESeq2]{vst}, \link[DESeq2]{rlog}}, etc.
#' @param group Optional character or factor vector of length equal to sample size,
#'   or up to two such vectors organized into a list or data frame. Supply legend
#'   title(s) by passing a named list or data frame.
#' @param covar Optional continuous covariate. If non-\code{NULL}, then function can
#'   take at most one \code{group} variable. Supply legend title by passing a named
#'   list or data frame.
#' @param top Optional number (if > 1) or proportion (if < 1) of top probes to be used
#'   for MDS. See Details.
#' @param pcs Vector specifying which principal coordinates to plot. Must be of length
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
#' three-dimensional principal coordinate subspace. MDS is an easy and popular method
#' for unsupervised cluster detection. It can also aid in spotting potential outliers,
#' and generally helps to visualize the latent structure of a data set.
#'
#' The \code{top} argument filters probes using the leading fold change method of
#' Smyth et al. (See \code{limma::\link[limma]{plotMDS}}). Pairwise Euclidean
#' distances are calculated using the most differentially expressed probes between the
#' two samples. This method is appropriate when different molecular pathways are
#' relevant for distinguishing different pairs of samples. To run MDS on the complete
#' data, set \code{top = NULL}. This is functionally equivalent to running PCA on the
#' full matrix. See \code{\link{plot_pca}}.
#'
#' @references
#' Cox, T.F. & Cox, M.A.A. (2001). \emph{Multidimensional Scaling.} Second edition.
#' Chapman and Hall.
#'
#' Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., & Smyth, G.K. (2015).
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/25605792}{limma powers differential
#' expression analyses for RNA-sequencing and microarray studies}. \emph{Nucleic
#' Acids Res.}, emph{43}(7): e47.
#'
#' Torgerson, W.S. (1958). \emph{Theory and Methods of Scaling.} New York: Wiley.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 5), nrow = 1000, ncol = 5)
#' plot_mds(mat)
#'
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet()
#' dds <- rlog(dds)
#' plot_mds(mat, group = colData(dds)$condition)
#'
#' @seealso
#' \code{\link[limma]{plotMDS}, \link{plot_pca}, \link{plot_tsne}}
#'
#' @export
#' @import dplyr
#' @importFrom purrr map
#' @importFrom limma getEAWP
#' @importFrom wordspace dist.matrix
#' @importFrom ggsci scale_color_d3 pal_d3
#' @import ggplot2
#' @import plotly
#'

plot_mds <- function(dat,
                     group = NULL,
                     covar = NULL,
                       top = 500,
                       pcs = c(1, 2),
                     label = FALSE,
                     title = NULL,
                    legend = 'outside',
                     hover = FALSE,
                        D3 = FALSE) {

  # Preliminaries
  if (ncol(dat) < 3L) {
    stop(paste('dat includes only', ncol(dat), 'samples; need at least 3 for MDS.'))
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
  if (length(pcs) > 2L && !D3) {
    stop('pcs must be of length 2 when D3 = FALSE.')
  } else if (length(pcs) > 3L) {
    stop('pcs must be a vector of length <= 3.')
  }
  if (label && !is.null(features) && length(features) == 2L) {
    stop('If label is TRUE, then plot can render at most one phenotypic feature.')
  }
  if (is.null(title)) {
    title <- 'MDS'
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
  if (is.null(top)) {                                             # Distance matrix
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
  mds <- suppressWarnings(cmdscale(as.dist(dm), k = max(pcs)))    # MDS
  df <- data_frame(Sample = colnames(dat))                        # Melt
  if (length(pcs) == 2L) {
    df <- df %>% mutate(PC1 = mds[, min(pcs)],
                        PC2 = mds[, max(pcs)])
  } else {
    other <- setdiff(pcs, c(min(pcs), max(pcs)))
    df <- df %>% mutate(PC1 = mds[, min(pcs)],
                        PC2 = mds[, other],
                        PC3 = mds[, max(pcs)])
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
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    if (is.null(top)) {
      p <- p + labs(title = title,
                        x = paste0('PC', min(pcs)),
                        y = paste0('PC', max(pcs)))
    } else {
      p <- p + labs(title = title,
                        x = paste('Leading logFC Dim', min(pcs)),
                        y = paste('Leading logFC Dim', max(pcs)))
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

