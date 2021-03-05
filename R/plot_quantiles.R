#' Plot Quantiles
#'
#' This function plots the quantiles of two vectors against each other as either
#' a QQ or an MD plot.
#'
#' @param x Vector of numeric values.
#' @param y Second vector of numeric values for comparison.
#' @param method Plot quantiles against quantiles (\code{method = "QQ"}) or mean
#'   quantiles against difference in quantiles (\code{method = "MD"})?
#' @param pts Number of points to plot.
#' @param main Optional plot title.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#'
#' @details
#' Quantile-quantile (QQ) and mean-difference (MD) plots visualize the
#' relationship between two numeric vectors, or between expected and observed 
#' values for given variable. They are a quick and easy alternative to 
#' scatterplots when variables are of unequal length. They can slso help 
#' visually assess parametric goodness of fit.
#'
#' @examples
#' x1 <- rnorm(100)
#' x2 <- runif(500)
#' plot_quantiles(x1, x2)
#' plot_quantiles(x1, x2, method = "MD")
#'
#' @export
#' @importFrom purrr keep
#' @import dplyr
#' @import ggplot2
#'

plot_quantiles <- function(
  x,
  y,
  method = c('QQ', 'MD'),
     pts = 1000L,
   title = ifelse(method == 'QQ', 'QQ Plot', 'MD Plot'),
    xlab = NULL,
    ylab = NULL
) {

  # Preliminaries
  x <- x %>% keep(is.finite)
  if (length(x) < 1L) {
    stop('x must have at least one finite, non-missing value.')
  }
  y <- y %>% keep(is.finite)
  if (length(y) < 1L) {
    stop('y must have at least one finite, non-missing value.')
  }
  method <- match.arg(method, c('QQ', 'MD'))
  if (xlab %>% is.null) {
    xlab <- 'X'
  }
  if (ylab %>% is.null) {
    ylab <- 'Y'
  }

  # Tidy data
  x <- quantile(x, probs = seq(0L, 1L, length.out = pts))
  y <- quantile(y, probs = seq(0L, 1L, length.out = pts))
  if (method == 'QQ') {
    df <- tibble(X = x, Y = y)
  } else {
    df <- tibble(X = (x + y) / 2L, Y = x - y)
  }

  # Build plot
  size <- pt_size(df)
  p <- ggplot(df, aes(X, Y)) +
    geom_point(size = size) +
    labs(title = title, x = xlab, y = ylab) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (method == 'QQ') {
    p <- p + geom_abline(intercept = 0L, slope = 1L, color = 'red', size = 0.2)
  } else {
    p <- p + geom_hline(yintercept = 0L, color = 'red', size = 0.2)
  }

  # Output
  print(p)

}


# Extend to theoretical distros?

