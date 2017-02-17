#' Lorenz Curve(s)
#'
#' This functions plots empirical Lorenz curves for one or several vectors.
#'
#' @param dat Numeric vector, or several such vectors organized into a list
#'   (optionally named). May also be a data frame of such vectors, however in that
#'   case each must be of equal length. Lorenz curves may be plotted and Gini
#'   coefficients calculated for vectors with negative values, but results are
#'   easiest to interpret when data are exclusively positive.
#' @param xlab Optional label for x-axis.
#' @param ylab Optional label for y-axis.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside", "bottomleft",
#'   "bottomright", "topleft",} or \code{"topright"}.
#' @param hover Show vector name by hovering mouse over Lorenz curve? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#'
#' @details
#' The Lorenz curve of a distribution plots its cumulative proportion of observations
#' against its cumulative proportion of values. The extent to which the curve sags
#' below the straight diagonal line indicates variable's degree of inequality. This
#' is measured by the Gini coefficient, which represents the ratio of the area between
#' the line of perfect equality and the distribution's Lorenz curve to the total area
#' under the diagonal line. The statistic has range [0, 1] for non-negative data, with
#' higher coefficients corresponding to more unequal distributions.
#'
#' Lorenz curves are common in economics, where Gini coefficients are used to measure
#' the inequality of national income distributions. They also have sevel omic
#' applications, e.g. to visualize the degree of biodiversity in a microbiome, or the
#' spread of RNA-sequencing counts across libraries.
#'
#' @examples
#' x <- runif(100)
#' plot_lorenz(x)
#'
#' X <- list("x1" = runif(100), "x2" = rpois(200, lambda = 5))
#' plot_lorenz(X)
#'
#' @export
#' @import dplyr
#' @importFrom purrr map map_chr
#' @import ggplot2
#' @importFrom scales hue_pal
#' @importFrom plotly ggplotly
#'

plot_lorenz <- function(dat,
                        main = NULL,
                        xlab = NULL,
                        ylab = NULL,
                      legend = 'topleft',
                       hover = FALSE) {

  # Preliminaries
  if (is.data.frame(dat)) {
    dat <- as.list(dat)
  } else if (!is.list(dat)) {
    dat <- list(dat)
  }
  if (is.null(names(dat))) {
    names(dat) <- paste0('x', seq_along(dat))
  }
  for (i in seq_along(dat)) {
    if (!is.numeric(dat[[i]])) {
      stop('dat must be a numeric vector, or several such vectors organized into ',
           'a list or data frame.')
    }
    if (min(dat[[i]] < 0)) {
      warning('Lorenz curves and Gini coefficients for data with negative values ',
              'should be interpreted with caution.')
    }
    if (length(na.omit(dat[[i]])) <= 2) {
      warning('Lorenz curves and Gini coefficients for vectors of length <= 2 ',
              'should be interpreted with caution.')
    }
  }
  if (is.null(xlab)) {
    xlab <- 'Cumulative Proportion of Observations'
  }
  if (is.null(ylab)) {
    ylab <- 'Cumulative Proportion of Values'
  }
  if (is.null(main)) {
    if (length(dat) == 1) {
      main <- 'Lorenz Curve'
    } else {
      main <- 'Lorenz Curves'
    }
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright".')
  }

  # Tidy data
  dfs <- map(seq_along(dat), function(i) {
    x <- sort(na.omit(dat[[i]]))
    n <- rep(1, length(x))
    p <- cumsum(n) / sum(n)
    L <- cumsum(x) / sum(x)
    p <- c(0, p)
    L <- c(0, L)
    df <- data_frame(Title = names(dat)[i],
                Proportion = p,
                    Lorenz = L)
    return(df)
  })

  # Plot
  gini <- function(x) {          # Calculate Gini coefficient
    x <- sort(na.omit(x))
    n <- length(x)
    g <- (2 * sum(x * seq_len(n)) / sum(x) - (n + 1)) / n
    return(g)
  }
  leg <- function(i) {           # Print Gini coefficient
    paste0(names(dat)[i], ', Gini = ', round(gini(dat[[i]]), 2))
  }
  p <- ggplot() +
    geom_abline(intercept = 0, slope = 1) +
    labs(title = main, x = xlab, y = ylab) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (length(dat) > 1) {         # Multiple curves?
    for (i in seq_along(dat)) {
      suppressWarnings(
        p <- p + geom_path(data = dfs[[i]], aes(Proportion, Lorenz,
                                                text = Title, color = Title)) +
          scale_colour_manual(name = 'Data',
                            labels = map_chr(seq_along(dat), leg),
                            values = hue_pal()(length(dat)))
      )
    }
  } else {
    p <- p + geom_path(data = dfs[[1]], aes(Proportion, Lorenz)) +
      scale_colour_manual(name = 'Data',
                        labels = map_chr(seq_along(dat), leg),
                        values = hue_pal()(length(dat)))
  }
  if (legend == 'bottomleft') {  # Locate legend
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


