#' Probability Calibration Curve(s)
#'
#' This functions plots probability calibration curves for one or several
#' classifiers.
#'
#' @param obs Vector of observed outcomes. Must be dichotomous. Can be numeric,
#'   logical, character, or factor. If numeric, \code{obs} must be coded \code{
#'   1} or \code{0}. If character or factor, a warning will be issued clarifying
#'   that the first level is assumed to be the reference.
#' @param pred Vector of predicted probabilities, or several such vectors
#'   organized into a data frame or list, optionally named. Must be numeric on
#'   \code{[0, 1]}.
#' @param pal_curves String specifying the color palette to use when plotting
#'   multiple vectors. Options include \code{"ggplot"}, all qualitative color
#'   schemes available in \code{RColorBrewer}, and the complete collection of
#'   \code{\href{http://bit.ly/2bxnuGB}{ggsci}} palettes. Alternatively, a
#'   character vector of colors with length equal to the number of vectors
#'   in \code{dat}.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
#' @param hover Show predictor name by hovering mouse over ROC curve? If \code{
#'   TRUE}, the plot is rendered in HTML and will either open in your browser's
#'   graphic display or appear in the RStudio viewer.
#'
#' @details
#' Calibration curves are a quick and easy way to evaluate a classifier's fit to
#' the data. This function allows one or several models to be plotted in the
#' same figure, with points sized by the number of observations that fall within
#' the corresponding bin.
#'
#' @examples
#' x1 <- runif(1000)
#' y <- rbinom(1000, size = 1, prob = x1)
#' plot_calibration(obs = y, pred = x1)
#'
#' x2 <- rbeta(1000, shape1 = 5/2, shape2 = 3/2)
#' plot_calibration(obs = y, pred = list("Good" = x1, "Bad" = x2))
#'
#' @export
#' @importFrom purrr some map_lgl map map_dbl map_df
#' @importFrom dplyr data_frame
#' @import ggplot2
#'

plot_calibration <- function(obs,
                             pred,
                             pal_curves = 'npg',
                                  title = NULL,
                                 legend = 'right',
                                  hover = FALSE) {

  # Preliminaries
  obs <- format_binom(obs, vec_type = 'obs')
  pred <- format_binom(pred, vec_type = 'pred', n = length(obs))
  if (pred %>% some(function(m) max(m > 1L) || min(m < 0L))) {
    stop('pred values must be on [0, 1].')
  }
  if (length(pred) > 1L) {
    cols <- colorize(pal_curves, var_type = 'Categorical', n = length(pred))
  }
  if (title %>% is.null) {
    if (length(pred) == 1L) {
      title <- 'Calibration Curve'
    } else {
      title <- 'Calibration Curves'
    }
  }
  loc <- c('bottom', 'left', 'top', 'right',
           'bottomright', 'bottomleft', 'topleft', 'topright')
  if (!legend %in% loc) {
    stop('legend must be one of ', stringify(loc, 'or'), '.')
  }

  # Tidy data
  brks <- seq(from = 0.05, to = 1L, by = 0.05)
  bin <- pred %>% map(function(m) {
    seq_along(m) %>% map_dbl(function(i) which.max(m[i] <= brks))
  })
  exp_grps <- seq_along(pred) %>% map(function(m) split(pred[[m]], bin[[m]]))
  obs_grps <- seq_along(bin) %>% map(function(x) split(obs, bin[[x]]))
  df <- seq_along(pred) %>% map_df(function(m) {
    data_frame(Y = obs_grps[[m]] %>% map_dbl(mean),
               X = exp_grps[[m]] %>% map_dbl(mean),
       Frequency = exp_grps[[m]] %>% map_dbl(length),
      Classifier = names(pred)[m])
  })

  # Build plot
  p <- ggplot(df, aes(X, Y)) +
    geom_abline(intercept = 0L, slope = 1L, color = 'grey') +
    scale_size(range = c(1L, 5L)) +
    labs(title = title, x = 'Expected Probability', y = 'Observed Probability') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (length(pred) > 1L) {                       # Multiple curves?
    suppressWarnings(
      p <- p + geom_point(aes(size = Frequency, color = Classifier,
                              text = Classifier)) +
        geom_path(aes(text = Classifier, color = Classifier)) +
        scale_color_manual(values = cols)
    )
  } else {
    p <- p + geom_point(aes(size = Frequency)) +
      geom_path()
  }

  # Output
  gg_out(p, hover, legend)

}

# Use gganimate, tweenr, and shiny to toggle btw classifiers
