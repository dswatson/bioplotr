#' Lorenz Curve(s)
#'
#' This functions plots empirical Lorenz curves for one or several vectors.
#'
#' @param dat Numeric vector, or several such vectors organized into a list,
#'   optionally named. May also be a data frame of such vectors, however in that
#'   case each must be of equal length. Data may include negative values, but if
#'   so a warning will be issued to proceed with caution.
#' @param pal_curves String specifying the color palette to use when plotting
#'   multiple vectors. Options include \code{"ggplot"}, all qualitative color
#'   schemes available in \code{RColorBrewer}, and the complete collection of
#'   \code{\href{http://bit.ly/2bxnuGB}{ggsci}} palettes. Alternatively, a
#'   character vector of colors with length equal to the number of vectors
#'   in \code{dat}.
#' @param title Optional plot title.
#' @param leg.txt Optional legend title.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
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
#' @importFrom purrr map keep map_chr
#' @import dplyr
#' @import ggplot2
#'

plot_lorenz <- function(dat,
                        pal_curves = 'npg',
                             title = NULL,
                           leg.txt = NULL,
                            legend = 'topleft',
                             hover = FALSE) {

  # Preliminaries
  if (dat %>% is.data.frame) {
    dat <- dat %>% as.list(.)
  } else if (!(dat %>% is.list)) {
    dat <- list(dat)
  }
  if (names(dat) %>% is.null) {
    names(dat) <- paste0('x', seq_along(dat))
  }
  for (i in seq_along(dat)) {
    if (!(dat[[i]] %>% is.numeric)) {
      stop('dat must be a numeric vector, or several such vectors organized ',
           'into a list or data frame.')
    }
    if (min(dat[[i]]) < 0L) {
      warning('Lorenz curves and Gini coefficients for data with negative ',
              'values should be interpreted with caution.')
    }
    if (length(dat[[i]] %>% na.omit(.)) <= 2L) {
      warning('Lorenz curves and Gini coefficients for vectors of length <= 2 ',
              'should be interpreted with caution.')
    }
  }
  if (length(dat) > 1L) {
    cols <- colorize(pal_curves, var_type = 'Categorical', n = length(dat))
  }
  if (title %>% is.null) {
    if (length(dat) == 1L) {
      title <- 'Lorenz Curve'
    } else {
      title <- 'Lorenz Curves'
    }
  }
  if (leg.txt %>% is.null) {
    leg.txt <- 'Data'
  }
  loc <- c('bottom', 'left', 'top', 'right',
           'bottomright', 'bottomleft', 'topleft', 'topright')
  if (!legend %in% loc) {
    stop('legend must be one of ', stringify(loc, 'or'), '.')
  }

  # Tidy data
  dfs <- seq_along(dat) %>% map(function(j) {
    x <- dat[[j]] %>%
      keep(is.finite) %>%
      sort(.)
    n <- rep(1L, length(x))
    p <- c(0L, cumsum(n) / sum(n))
    L <- c(0L, cumsum(x) / sum(x))
    data_frame(Title = names(dat)[j],
          Proportion = p,
              Lorenz = L) %>%
      return(.)
  })

  # Build plot
  gini <- function(x) {                     # Calculate Gini coefficient
    x <- x %>%
      keep(is.finite) %>%
      sort(.) %>%
      as.numeric(.)
    n <- length(x)
    g <- (2L * sum(x * seq_len(n)) / sum(x) - (n + 1L)) / n
    return(g)
  }
  p_gin <- function(i) {                    # Print Gini coefficient
    paste0(names(dat)[i], ', Gini = ',
           gini(dat[[i]]) %>% round(2L))
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
                              labels = seq_along(dat) %>% map_chr(p_gin),
                              values = cols)
  } else {
    p <- p + geom_path(data = dfs[[1L]], aes(Proportion, Lorenz, color = Title)) +
      scale_color_manual(name = leg.txt,
                       labels = seq_along(dat) %>% map_chr(p_gin),
                       values = 'black')
  }

  # Output
  gg_out(p, hover, legend)

}

# Use gganimate, tweenr, and shiny to toggle between vectors
