#' Lorenz Curve(s)
#'
#' This functions plots empirical Lorenz curves for one or several vectors.
#'
#' @param dat Numeric vector, or several such vectors organized into a list,
#'   optionally named. May also be a data frame of such vectors, however in that
#'   case each must be of equal length. Data may include negative values, but if
#'   so a warning will be issued to proceed with caution.
#' @param pal_curves String specifying the color palette to use when plotting
#'   multiple vectors. Options include \code{"ggplot"}, as well as the complete
#'   collection of \code{
#'   \href{https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html}{
#'   ggsci}} palettes, which can be identified by name (e.g., \code{"npg"},
#'   \code{"aaas"}, etc.). Alternatively, a character vector of colors with
#'   length equal to the number of vectors in \code{dat}.
#' @param title Optional plot title.
#' @param leg.txt Optional legend title.
#' @param legend Legend position. Must be one of \code{"right"}, \code{
#'   "left"}, \code{"top"}, \code{"bottom"}, \code{"topright"}, \code{
#'   "topleft"}, \code{"bottomright"}, or \code{"bottomleft"}.
#' @param hover Show vector name by hovering mouse over Lorenz curve? If \code{
#'   TRUE}, the plot is rendered in HTML and will either open in your browser's
#'   graphic display or appear in the RStudio viewer.
#'
#' @details
#' The Lorenz curve of a distribution plots its cumulative proportion of
#' observations against its cumulative proportion of values. The extent to which
#' the curve sags below the straight diagonal line indicates the variable's
#' degree of inequality. This is measured by the Gini coefficient, which
#' represents the ratio of the area between the line of perfect equality and the
#' distribution's Lorenz curve to the total area under the diagonal line. The
#' statistic has range [0, 1] for non-negative data, with higher coefficients
#' corresponding to more unequal distributions.
#'
#' Lorenz curves are common in economics, where Gini coefficients are often used
#' to measure the inequality of national income distributions. They also have
#' several omic applications, e.g. to visualize the degree of biodiversity in a
#' microbiome, or the spread of RNA-sequencing counts across libraries.
#'
#' @examples
#' x <- runif(100)
#' plot_lorenz(x)
#'
#' X <- list("Uniform" = runif(100),
#'           "Poisson" = rpois(200, lambda = 5))
#' plot_lorenz(X)
#'
#' @export
#' @importFrom dplyr data_frame
#' @importFrom purrr map map_chr
#' @import ggplot2
#'

plot_lorenz <- function(dat,
                        pal_curves = 'npg',
                             title = NULL,
                           leg.txt = NULL,
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
      stop('dat must be a numeric vector, or several such vectors organized ',
           'into a list or data frame.')
    }
    if (min(dat[[i]] < 0L)) {
      warning('Lorenz curves and Gini coefficients for data with negative ',
              'values should be interpreted with caution.')
    }
    if (length(na.omit(dat[[i]])) <= 2L) {
      warning('Lorenz curves and Gini coefficients for vectors of length <= 2 ',
              'should be interpreted with caution.')
    }
  }
  if (length(dat) > 1L) {
    cols <- colorize(pal_curves, var_type = 'Categorical', n = length(dat))
  }
  if (is.null(title)) {
    if (length(dat) == 1L) {
      title <- 'Lorenz Curve'
    } else {
      title <- 'Lorenz Curves'
    }
  }
  if (is.null(leg.txt)) {
    leg.txt <- 'Data'
  }
  if (!legend %in% c('right', 'left', 'top', 'bottom',
                     'topright', 'topleft', 'bottomright', 'bottomleft')) {
    stop('legend must be one of "right", "left", "top", "bottom", ',
         '"topright", "topleft", "bottomright", or "bottomleft".')
  }

  # Tidy data
  dfs <- map(seq_along(dat), function(i) {
    x <- sort(dat[[i]][is.finite(dat[[i]])])
    n <- rep(1L, length(x))
    p <- c(0L, cumsum(n) / sum(n))
    L <- c(0L, cumsum(x) / sum(x))
    df <- data_frame(Title = names(dat)[i],
                Proportion = p,
                    Lorenz = L)
    return(df)
  })

  # Build plot
  gini <- function(x) {                     # Calculate Gini coefficient
    x <- as.numeric(sort(x[is.finite(x)]))
    n <- length(x)
    g <- (2L * sum(x * seq_len(n)) / sum(x) - (n + 1L)) / n
    return(round(g, 2L))
  }
  p_gin <- function(i) {                    # Print Gini coefficient
    paste0(names(dat)[i], ', Gini = ', gini(dat[[i]]))
  }
  p <- ggplot() +
    geom_abline(intercept = 0L, slope = 1L, color = 'grey') +
    labs(title = title,
             x = 'Cumulative Proportion of Observations',
             y = 'Cumulative Proportion of Values') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (length(dat) > 1L) {                   # Multiple curves?
    for (i in seq_along(dat)) {
      suppressWarnings(
        p <- p + geom_path(data = dfs[[i]], aes(Proportion, Lorenz,
                                                text = Title, color = Title))
      )
    }
    p <- p + scale_color_manual(name = leg.txt,
                              labels = map_chr(seq_along(dat), p_gin),
                              values = cols)
  } else {
    p <- p + geom_path(data = dfs[[1]], aes(Proportion, Lorenz, color = Title)) +
      scale_color_manual(name = leg.txt,
                       labels = map_chr(seq_along(dat), p_gin),
                       values = 'black')
  }

  # Output
  gg_out(p, hover, legend)

}

# Use gganimate, tweenr, and shiny to toggle between vectors
