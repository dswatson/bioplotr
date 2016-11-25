#' Create box plots by sample
#'
#' @param dat Omic data matrix with rows corresponding to probes and columns
#'   to samples.
#' @param group Optional factor or character vector of length equal to sample size.
#'   Levels are used to color box plots. If supplied, legend title defaults to
#'   "Group". Override this feature by passing a named list instead.
#' @param type Optional string specifying omic data type. Currently supports
#'   \code{"microarray", "RNA-seq",} and \code{"methylation"}.
#' @param ylab Label for y-axis. If left \code{NULL}, this defaults to
#'   log expression for \code{type = "microarray"}, log CPM for
#'   \code{type = "RNA-seq"}, and Beta for \code{type = "methylation"}. At least
#'   one of \code{type} or \code{xlab} must be specified.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside",
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
#' @param hover Show sample name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#'
#' @details
#' This function displays each sample's omic data distribution as a box and whisker
#' plot. This is especially helpful when contrasting pre- and post-normalization
#' matrices. It may additionally be used to inspect for batch effects or other
#' associations with phenotypic features by using the \code{group} argument.
#'
#' @examples
#' mat <- matrix(rnorm(5000), nrow = 1000, ncol = 5)
#' plot_box(mat, type = "microarray")
#'
#' library(edgeR)
#' mat <- cbind(matrix(rnbinom(5000, mu = 4, size = 1), nrow = 1000, ncol = 5),
#'              matrix(rnbinom(5000, mu = 4, size = 10), nrow = 1000, ncol = 5))
#' mat <- calcNormFactors(DGEList(mat))
#' mat <- cpm(mat, log = TRUE, prior.count = 0.5)
#' batch <- rep(c("A", "B"), each = 5)
#' plot_box(mat, group = batch, type = "RNA-seq")
#'
#' @export
#' @importFrom tidyr gather
#' @import dplyr
#' @import ggplot2
#' @importFrom plotly ggplotly
#'

plot_box <- function(dat,
                     group  = NULL,
                     type   = NULL,
                     ylab   = NULL,
                     main   = NULL,
                     legend = 'outside',
                     hover  = FALSE) {

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
  if (is.null(ylab)) {
    if (is.null(type)) {
      stop('Either data type or ylab must be provided.')
    } else if (type == 'microarray') {
      ylab <- expression('log'[2]*' Expression')
    } else if (type == 'RNA-seq') {
      ylab <- expression('log'[2]*' Counts Per Million')
    } else if (type == 'methylation') {
      ylab <- 'Beta'
    } else if (!type %in% c('microarray', 'RNA-seq', 'methylation')) {
      stop('type must be one of "microarray", "RNA-seq",  "methylation", or NULL.')
    }
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright",
  "topleft", or "topright".')
  }
  if (is.null(main)) {
    if (is.null(group)) {
      main <- 'Expression By Sample'
    } else {
      if (is.null(names(group))) {
        main <- 'Expression By Group'
      } else {
        main <- paste('Expression By', names(group))
      }
    }
  }

  # Tidy
  df <- gather(as_data_frame(dat), Sample, Expression)
  if (!is.null(group)) {
    if (!is.list(group)) {
      df <- df %>% mutate(Group = rep(group, each = nrow(dat)))
    } else {
      df <- df %>% mutate(Group = rep(group[[1]], each = nrow(dat)))
    }
  } else {
    df <- df %>% mutate(Group = rep(1, nrow(df)))
  }
  df <- df %>%
    arrange(Group) %>%
    mutate(Sample = factor(Sample, levels = unique(Sample)))

  # Basic plot
  p <- ggplot(df, aes(Sample, Expression)) +
    labs(title = main, x = 'Sample', y = ylab) +
    theme_bw() +
    theme(plot.title = element_text(hjust = .5),
          axis.text.x = element_text(angle = 45, hjust = 1))
  if (!is.null(group)) {
    p <- p + suppressWarnings(geom_boxplot(aes(text = Sample,
                                               fill = Group)))
  } else {
    p <- p + suppressWarnings(geom_boxplot(text = Sample))
  }

  # Named list?
  if (!is.null(names(group))) {
    p <- p + guides(fill = guide_legend(title = names(group)))
  }

  # Legend location
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
  if (hover == FALSE) {
    print(p)
  } else {
    p <- ggplotly(p, tooltip = 'text', height = 600, width = 600)
    print(p)
  }

}

# Fun fact: plotly won't display text for boxplots:
# https://community.plot.ly/t/boxplot-hoverinfo-text-not-display/1959

