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
#' @param hover Show sample name by hovering mouse over data point? Renders
#'   the plot in HTML, which opens the graphic display in your browser. The plot can
#'   also be embedded in an HTML doc using Rmarkdown so long as \code{knitr = TRUE}
#'   and you set the code chunk option \code{plotly = TRUE}.
#' @param D3 Render the plot in three dimensions? Like \code{hover}, this creates a
#'   plotly object that opens in your browser or can be embedded in an HTML doc.
#' @param knitr Set this to \code{TRUE} if you want to embed a plotly object (i.e.,
#'   the \code{plot_pca} output when \code{hover = TRUE} or \code{D3 = TRUE}) in
#'   an HTML doc. Make sure to set \code{plotly = TRUE} in the corresponding code
#'   chunk options.
#'
#' @details
#' This function plots the samples of an omic data matrix in a two- or three-
#' dimensional principal component subspace. The numbers printed along the axis
#' labels indicate the percentage of variance explained by each component. PCA is an
#' easy and popular method for unsupervised cluster detection. It can also aid in
#' spotting potential outliers and generally  help to visualize the latent
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
#' @import ggplot2
#' @import plotly
#'

plot_pca <- function(dat,
                     group  = NULL,
                     label  = FALSE,
                     main   = NULL,
                     legend = 'outside',
                     hover  = FALSE,
                     D3     = FALSE,
                     knitr  = FALSE) {

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
  if (is.null(colnames(dat))) {
    sample <- paste0('Sample', seq_along(1:ncol(dat)))
  } else {
    sample <- colnames(dat)
  }
  if (is.null(main)) {
    main <- 'PCA'
  }

  # PCA
  pca <- prcomp(t(dat), center = TRUE, scale. = TRUE)
  var1 <- round(pca$sdev[1]^2 / sum(pca$sdev^2) * 100, 2)
  var2 <- round(pca$sdev[2]^2 / sum(pca$sdev^2) * 100, 2)

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
         x = paste0('PC1 (', var1, '%)'),
         y = paste0('PC2 (', var2, '%)')) +
    theme_bw() +
    theme(plot.title = element_text(hjust = .5))

  # Sample labels
  if (label == TRUE) {
    if (!is.numeric(df$Group)) {
      p <- p + geom_text(aes(label = Sample, color = Group, alpha = 0.85))
    } else {
      p <- p + geom_text(aes(label = Sample, alpha = 0.85))
    }
  } else {
    if (!is.numeric(df$Group)) {
      suppressWarnings(p <- p + geom_point(aes(text  = paste('Sample:', Sample),
                                               color = Group, shape = Group),
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
      if (knitr == FALSE) {
        p <- ggplotly(p, tooltip = 'text', width = 600, height = 500)
        print(p)
      } else {
        p <- ggplotly(p, tooltip = 'text', width = 600, height = 500,
                      session = 'knitr')
        print(p)
      }
    }
  } else {
      p <- plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3,
                   text = ~Sample, color = ~Group, symbol = ~Group,
                   colors = c('#E41A1C', '#4DAF4A', '#377EB8'),
                   symbols = c('circle', 'triangle-up', 'square'),
                   type = 'scatter3d', mode = 'markers',
                   alpha = 0.85, hoverinfo = 'text', marker = list(size = 5)) %>%
        layout(title = main, hovermode = 'closest', scene = list(
          xaxis = list(title = paste0('PC1 (', var1, '%)')),
          yaxis = list(title = paste0('PC2 (', var2, '%)')),
          zaxis = list(title = paste0('PC3 (', var3, '%)'))
        ))
      print(p)
  }

}


