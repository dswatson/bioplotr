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
#'   is non-\code{NULL}. Options include \code{"ggplot"}, as well as the
#'   complete collection of \code{
#'   \href{https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html}{
#'   ggsci}} palettes, which can be identified by name (e.g., \code{"npg"},
#'   \code{"aaas"}, etc.). Alternatively, any character vector of colors with
#'   length equal to the cumulative number of levels in \code{group}.
#' @param pal_covar String specifying the color palette to use if \code{covar}
#'   is non-\code{NULL}. Options include \code{"blues"}, \code{"greens"}, \code{
#'   "purples"}, \code{"greys"}, \code{"oranges"}, and \code{"reds"}.
#'   Alternatively, any character vector of colors representing a smooth
#'   gradient.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"right"}, \code{
#'   "left"}, \code{"top"}, \code{"bottom"}, \code{"bottomright"},
#'   \code{"bottomleft"}, \code{"topright"}, or \code{"topleft"}.
#' @param hover Show sample name by hovering mouse over data point? If \code{
#'   TRUE}, the plot is rendered in HTML and will either open in your browser's
#'   graphic display or appear in the RStudio viewer.
#' @param D3 Render plot in three dimensions?
#'
#' @details
#' This function plots the samples of an omic data matrix in a two- or
#' three-dimensional principal component subspace. Axis labels include the
#' percentage of variance explained by each component. PCA is an easy and
#' popular method for unsupervised cluster detection. It can also aid in
#' spotting potential outliers, and generally helps to visualize the latent
#' structure of a data set.
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
                 pal_group = 'd3',
                 pal_covar = 'blues',
                     title = NULL,
                    legend = 'outside',
                     hover = FALSE,
                        D3 = FALSE) {

  # Preliminaries
  if (ncol(dat) < 3L) {
    stop(paste('dat includes only', ncol(dat), 'samples;',
               'need at least 3 for PCA.'))
  }
  if (!is.null(group)) {
    group <- format_features(dat, group, var_type = 'Categorical')
    if (length(group) > 2L) {
      stop('Plot can render at most two categorical features.')
    }
    if (length(group) == 2L & !is.null(covar)) {
      stop('Plot can render at most one categorical feature when a continuous ',
           'covariate is also supplied.')
    }
  }
  if (!is.null(covar)) {
    covar <- format_features(dat, covar, var_type = 'Continuous')
    if (length(covar) != 1L) {
      stop('Plot can render at most one continuous feature.')
    }
  }
  if (!is.null(c(group, covar))) {
    features <- c(covar, group)
    feature_names <- names(features)
    names(features) <- paste0('Feature', seq_along(features))
  } else {
    features <- NULL
  }
  if (length(pcs) > 2L & !D3) {
    stop('pcs must be of length 2 when D3 = FALSE.')
  } else if (length(pcs) > 3L) {
    stop('pcs must be a vector of length <= 3.')
  }
  if (label & length(features) == 2L) {
    stop('If label is TRUE, then plot can render at most one phenotypic ',
         'feature.')
  }
  if (length(pal_group) == 1L & !is.color(pal_group)) {
    if (!pal_group %in% c('ggplot', 'npg', 'aaas', 'nejm', 'lancet', 'jco',
                          'ucscgb', 'd3', 'locuszoom', 'igv', 'uchicago',
                          'startrek', 'futurama', 'rickandmorty', 'simpsons',
                          'gsea')) {
      stop('pal_group not recognized.')
    }
  } else {
    if (!all(is.color(pal_group))) {
      stop('When passing multiple strings to pal_group, each must denote a ',
           'valid color in R.')
    }
    if (length(levels(group[[1]])) != length(pal_group)) {
      stop('When passing individual colors to pal_group, length(pal_group) ',
           'must equal the number of unique groups.')
    }
  }
  if (length(pal_covar) == 1L & !is.color(pal_covar)) {
    if (!pal_covar %in% c('blues', 'greens', 'purples',
                          'greys', 'oranges', 'reds')) {
      stop('pal_covar not recognized.')
    }
  } else {
    if (!all(is.color(pal_covar))) {
      stop('When passing multiple strings to pal_covar, each must denote a ',
           'valid color in R.')
    }
  }
  if (is.null(title)) {
    title <- 'PCA'
  }
  if (!legend %in% c('right', 'left', 'top', 'bottom', 'bottomright',
                     'bottomleft', 'topright', 'topleft')) {
    stop('legend must be one of "right", "left", "top", "bottom", ',
         '"bottomright", "bottomleft", "topright", or "topleft".')
  }

  # Tidy data
  dat <- matrixize(dat)
  if (is.null(rownames(dat))) {
    rownames(dat) <- seq_len(nrow(dat))
  }
  if (is.null(colnames(dat))) {
    colnames(dat) <- paste0('Sample', seq_len(ncol(dat)))
  }
  if (!is.null(top)) {                           # Filter by variance?
    if (top > 1L) {
      if (top > nrow(dat)) {
        warning(paste('top exceeds nrow(dat), at least after removing probes',
                      'with missing values and/or applying a minimal expression',
                      'filter. Proceeding with the complete', nrow(dat), 'x',
                      ncol(dat), 'matrix.'))
      }
    } else {
      top <- round(top * nrow(dat))
    }
    vars <- rowVars(dat)
    keep <- order(vars, decreasing = TRUE)[seq_len(min(top, nrow(dat)))]
    dat <- dat[keep, , drop = FALSE]
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
  size <- pt_size(df)
  alpha <- pt_alpha(df)
  if (!D3) {
    p <- ggplot(df, aes(PC1, PC2)) +
      geom_hline(yintercept = 0L, color = 'grey') +
      geom_vline(xintercept = 0L, color = 'grey') +
      labs(title = title, x = pve[min(pcs)], y = pve[max(pcs)]) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    if (is.null(features)) {
      p <- p + geom_text(aes(label = Sample), alpha = alpha,
                         hjust = 'inward', vjust = 'inward')
    } else if (length(features) == 1L) {
      if (label) {
        p <- p + geom_text(aes(label = Sample, color = Feature1), alpha = alpha,
                           hjust = 'inward', vjust = 'inward')
      } else {
        if (!is.null(covar)) {
          suppressWarnings(
            p <- p + geom_point(aes(text = Sample, color = Feature1),
                                size = size, alpha = alpha) +
              labs(color = feature_names[1L])
          )
        } else {
          suppressWarnings(
            p <- p + geom_point(aes(text = Sample, color = Feature1, shape = Feature1),
                                size = size, alpha = alpha) +
              labs(color = feature_names[1L], shape = feature_names[1L])
          )
        }
      }
    } else {
      suppressWarnings(
        p <- p + geom_point(aes(text = Sample, color = Feature1, shape = Feature2),
                            size = size, alpha = alpha) +
          labs(color = feature_names[1L], shape = feature_names[2L])
      )
      if (is.null(covar)) {
        p <- p + guides(color = guide_legend(order = 1L),
                        shape = guide_legend(order = 2L))
      } else {
        p <- p + guides(color = guide_colorbar(order = 1L),
                        shape = guide_legend(order = 2L))
      }
    }
    if (is.null(covar)) {
      p <- p + scale_color_manual(values = colorize(pal_group,
                                                    length(levels(group[[1L]])),
                                                    'Categorical'))
    } else {
      p <- p + scale_color_gradientn(colors = colorize(pal_covar,
                                                       var_type = 'Continuous'))
    }
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


