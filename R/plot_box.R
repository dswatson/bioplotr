#' Box Plots by Sample
#'
#' This function displays each sample's omic data distribution as a box and whisker
#' plot.
#'
#' @param dat Omic data matrix or matrix-like object with rows corresponding to
#'   probes and columns to samples.
#' @param group Optional character or factor vector of length equal to sample size.
#'   Numeric or logical vectors will be silently coerced to factor. Levels are used
#'   to color density curves. If supplied, legend title defaults to "Group". Override
#'   this feature by passing a named list instead.
#' @param ylab Optional label for y-axis.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside", "bottomleft",
#'   "bottomright", "topleft",} or \code{"topright"}.
#' @param hover Show sample name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#'
#' @details
#' Box plots are an intuitive way to visualize an omic data distribution. They are
#' especially helpful when contrasting pre- and post-normalization matrices.
#' \code{plot_box} may additionally be used to inspect for batch effects
#' or associations with phenotypic factors by using the \code{group} argument.
#'
#' @examples
#' mat <- matrix(rnorm(1000 * 5), nrow = 1000, ncol = 5)
#' plot_box(mat)
#'
#' mat <- cbind(matrix(rnbinom(1000 * 5, size = 1, mu = 4), nrow = 1000, ncol = 5),
#'              matrix(rnbinom(1000 * 5, size = 3, mu = 4), nrow = 1000, ncol = 5))
#' mat <- lcpm(mat)
#' batch <- rep(c("A", "B"), each = 5)
#' plot_box(mat, group = batch, ylab = "Normalized Counts")
#'
#' @export
#' @importFrom limma getEAWP
#' @importFrom tidyr gather
#' @importFrom ggsci scale_fill_d3
#' @import dplyr
#' @import ggplot2
#' @importFrom plotly ggplotly
#'

plot_box <- function(dat,
                     group = NULL,
                      type = NULL,
                      ylab = NULL,
                     title = NULL,
                    legend = 'outside',
                     hover = FALSE) {

  # Preliminaries
  if (!is.null(group)) {
    if (!is.list(group)) {
      group <- list(group)
    }
    if (length(group) > 1L) {
      stop('group cannot be a list of length > 1.')
    }
    if (length(group[[1]]) != ncol(dat)) {
      stop('group length must match number of samples in dat.')
    }
    if (length(unique(group[[1]])) == 1L) {
      warning('group is invariant.')
    }
    group[[1]] <- as.factor(group[[1]])
    if (is.null(names(group))) {
      names(group) <- 'Group'
    }
  }
  if (is.null(xlab)) {
    xlab <- 'Value'
  }
  if (is.null(main)) {
    if (is.null(group)) {
      main <- 'Density By Sample'
    } else {
      main <- paste('Density By', names(group))
    }
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright" ',
         '"topleft", or "topright".')
  }

  # Tidy data
  dat <- getEAWP(dat)$expr
  keep <- rowSums(is.finite(dat)) == ncol(dat)
  dat <- dat[keep, , drop = FALSE]
  df <- gather(tbl_df(dat), Sample, Expression) %>%
    mutate(Sample = factor(Sample, levels = unique(Sample)))
  if (!is.null(group)) {
    df <- df %>%
      mutate(Group = rep(group[[1]], each = nrow(dat))) %>%
      arrange(Group)
  }

  # Build plot
  p <- ggplot(df, aes(Sample, Expression, text = Sample)) +
    labs(title = title, x = 'Sample', y = ylab) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45L, hjust = 1L))
  if (!is.null(group)) {                         # Fill by group?
    p <- p + geom_boxplot(aes(fill = Group)) +
      guides(fill = guide_legend(title = names(group))) +
      scale_fill_d3()
  } else {
    p <- p + geom_boxplot()
  }
  if (legend == 'bottomleft') {                  # Locate legend
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

}

# Fun fact: plotly won't display text for boxplots:
# https://community.plot.ly/t/boxplot-hoverinfo-text-not-display/1959

# Use gganimate, tweenr, and shiny to toggle btw matrices

