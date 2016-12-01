#' Plot the principal components of an omic data matrix
#'
#' @param dat Omic data matrix with rows corresponding to probes and columns
#'   to samples. It is strongly recommended that data be filtered and
#'   normalized prior to running PCA.
#' @param group Optional factor or character vector of length equal to sample size.
#'   Levels are used to color and shape points. If supplied, legend title defaults to
#'   "Group". Override this feature by passing a named list instead.
#' @param label Label data points by sample name? Defaults to \code{FALSE} unless
#'   \code{group = NULL}.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside",
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
#' @param hover Show sample name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#' @param D3 Render the plot in three dimensions?
#'
#' @details
#' This function plots the samples of an omic data matrix in a two- or three-
#' dimensional principal component subspace. The numbers printed along the axis
#' labels indicate the percentage of variance explained by each component. PCA is an
#' easy and popular method for unsupervised cluster detection. It can also aid in
#' spotting potential outliers and generally helps to visualize the latent
#' structure of a data set.
#'
#' @examples
#' mat <- matrix(rnorm(5000), nrow = 1000, ncol = 5)
#' plot_pca(mat)
#'
#' library(DESeq2)
#' mat <- cbind(matrix(rnbinom(5000, mu = 4, size = 1), nrow = 1000, ncol = 5),
#'              matrix(rnbinom(5000, mu = 4, size = 10), nrow = 1000, ncol = 5))
#' mat <- rlog(mat)
#' grp <- rep(c("A", "B"), each = 5)
#' plot_pca(mat, group = grp)
#'
#' @export
#' @import dplyr
#' @importFrom purrr map_dbl
#' @import ggplot2
#' @import plotly
#' @importFrom scales hue_pal
#'

plot_pca <- function(dat,
                     group  = NULL,
                     label  = FALSE,
                     main   = NULL,
                     legend = 'outside',
                     hover  = FALSE,
                     D3     = FALSE) {

  # Preliminaries
  if (is.null(group)) {
    group <- rep(1, ncol(dat))
  } else {
    if (!is.list(group)) {
      if (!is.character(group) & !is.factor(group)) {
        stop('group must be a character or factor variable.')
      }
    } else {
      if (!is.character(group[[1]]) & !is.factor(group[[1]])) {
        stop('group must be a character or factor variable.')
      }
    }
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright",
  "topleft", or "topright".')
  }
  if (is.null(colnames(dat))) {
    sample <- paste0('Sample', 1:ncol(dat))
  } else {
    sample <- colnames(dat)
  }
  if (is.null(main)) {
    main <- 'PCA'
  }

  # PCA
  pca <- prcomp(t(dat), center = TRUE, scale. = TRUE)
  ve <- function(pc) round(pca$sdev[pc]^2 / sum(pca$sdev^2) * 100, 2)
  vars <- map_dbl(1:3, ve)

  # Tidy
  df <- data_frame(PC1    = pca$x[, 1],
                   PC2    = pca$x[, 2],
                   PC3    = pca$x[, 3],
                   Sample = sample)
  if (!is.list(group)) {
    df <- df %>% mutate(Group = group)
  } else {
    df <- df %>% mutate(Group = group[[1]])
  }

  # Basic plot
  p <- ggplot(df, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, size = 0.2) +
    geom_vline(xintercept = 0, size = 0.2) +
    labs(title = main,
         x = paste0('PC1 (', vars[1], '%)'),
         y = paste0('PC2 (', vars[2], '%)')) +
    theme_bw() +
    theme(plot.title = element_text(hjust = .5))

  # Sample labels
  if (label == TRUE) {
    if (!is.numeric(df$Group)) {
      p <- p + geom_text(aes(label = Sample, color = Group), alpha = 0.85)
    } else {
      p <- p + geom_text(aes(label = Sample), alpha = 0.85)
    }
  } else {
    if (!is.numeric(df$Group)) {
      suppressWarnings(p <- p + geom_point(aes(text  = Sample,
                                               color = Group,
                                               shape = Group),
                                           alpha = 0.85))
    } else {
      p <- p + geom_text(aes(label = Sample))
    }
  }

  # Named list?
  if (!is.null(names(group)) & !is.numeric(df$Group)) {
    p <- p + guides(color = guide_legend(title = names(group)),
                    shape = guide_legend(title = names(group)))
  }

  # Locate legend
  if (legend == 'bottomleft') {
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
  if (D3 == FALSE) {
    if (hover == FALSE) {
      print(p)
    } else {
      p <- ggplotly(p, tooltip = 'text', height = 600, width = 600)
      print(p)
    }
  } else {
    # symbls <- c(16, 17, 15, 3, 7, 8)  # This would be right if plotly worked
    symbls <- c(16, 18, 15, 3, 7, 8)
    p <- plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3,
                 text = ~Sample, color = ~Group, symbol = ~Group,
                 colors = hue_pal()(length(unique(df$Group))),
                 symbols = symbls[1:length(unique(df$Group))],
                 type = 'scatter3d', mode = 'markers',
                 alpha = 0.85, hoverinfo = 'text', marker = list(size = 5)) %>%
      layout(hovermode = 'closest', title = main,
             xaxis = list(title = paste0('PC1 (', vars[1], '%)')),
             yaxis = list(title = paste0('PC2 (', vars[2], '%)')),
             zaxis = list(title = paste0('PC3 (', vars[3], '%)')))
    print(p)
  }

}


