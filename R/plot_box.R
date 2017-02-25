#' Box Plots by Sample
#'
#' This function displays each sample's omic data distribution as a box and
#' whicker plot.
#'
#' @param dat Omic data matrix with rows corresponding to probes and columns
#'   to samples.
#' @param group Optional character or factor vector of length equal to sample size.
#'   Levels are used to color box plots. If supplied, legend title defaults to
#'   "Group". Override this feature by passing a named list instead.
#' @param ylab Optional label for y-axis.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside",
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
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
#' plot_box(mat, ylab = "Normalized Expression")
#'
#' library(edgeR)
#' mat <- cbind(matrix(rnbinom(1000 * 5, mu = 4, size = 1), nrow = 1000, ncol = 5),
#'              matrix(rnbinom(1000 * 5, mu = 4, size = 5), nrow = 1000, ncol = 5))
#' mat <- calcNormFactors(DGEList(mat))
#' mat <- cpm(mat, log = TRUE)
#' batch <- rep(c("A", "B"), each = 5)
#' plot_box(mat, group = batch, ylab = "Normalized Counts")
#'
#' @export
#' @importFrom tidyr gather
#' @import dplyr
#' @import ggplot2
#' @importFrom plotly ggplotly
#'

plot_box <- function(dat,
                     group = NULL,
                      type = NULL,
                      ylab = NULL,
                      main = NULL,
                    legend = 'outside',
                     hover = FALSE) {

  # Preliminaries
  dat <- getEAWP(dat)
  dat <- dat$expr
  bad <- rowSums(is.finite(dat)) < ncol(dat)
  if (any(bad)) {
    dat <- dat[!bad, , drop = FALSE]
  }
  if (is.null(group)) {
    group <- list(rep(1, times = ncol(dat)))
  } else {
    if (!is.list(group)) {
      group <- list(group)
    }
    if (!is.character(group[[1]]) & !is.factor(group[[1]])) {
      stop('group must be a character or factor variable.')
    }
    if (length(group) > 1) {
      stop('group cannot be a list of length > 1.')
    }
    if (length(group[[1]]) != ncol(dat)) {
      stop('group length must match number of samples in dat.')
    }
    if (length(unique(group[[1]])) == 1) {
      warning('group is invariant.')
    }
  }
  if (is.null(ylab)) {
    ylab <- 'Value'
  }
  if (is.null(main)) {
    if (is.numeric(group[[1]])) {
      main <- 'Expression By Sample'
    } else {
      if (is.null(names(group))) {
        main <- 'Expression By Group'
      } else {
        main <- paste('Expression By', names(group))
      }
    }
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright" ',
         '"topleft", or "topright".')
  }

  # Tidy data
  df <- gather(tbl_df(dat), Sample, Expression) %>%
    mutate(Group = rep(group[[1]], each = nrow(dat))) %>%
    arrange(Group) %>%
    mutate(Sample = factor(Sample, levels = unique(Sample)))

  # Build plot
  suppressWarnings(
    p <- ggplot(df, aes(Sample, Expression, text = Sample)) +
      labs(title = main, x = 'Sample', y = ylab) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1))
  )
  if (!is.numeric(group[[1]])) {  # Fill by group?
    p <- p + geom_boxplot(aes(fill = Group))
  } else {
    p <- p + geom_boxplot()
  }
  if (!is.null(names(group))) {   # Named list?
    p <- p + guides(fill = guide_legend(title = names(group)))
  }
  if (legend == 'bottomleft') {   # Locate legend
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
    p <- ggplotly(p, tooltip = 'text', height = 600, width = 600)
    print(p)
  }

}

# Fun fact: plotly won't display text for boxplots:
# https://community.plot.ly/t/boxplot-hoverinfo-text-not-display/1959

