#' Probability Calibration Curve(s)
#'
#' This functions plots probability calibration curves for one or several classifiers.
#'
#' @param obs Vector of observed outcomes. Must be dichotomous. Can be numeric,
#'   logical, character, or factor. If numeric, \code{obs} must be coded \code{1}
#'   or \code{0}. If character or factor, a warning will be issued clarifying that
#'   the first level is assumed to be the reference.
#' @param pred Vector of predicted probabilities, or several such vectors organized
#'   into a data frame or list, optionally named. Must be numeric on \code{(0, 1)}.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside", "bottomleft",
#'   "bottomright", "topleft",} or \code{"topright"}.
#' @param hover Show predictor name by hovering mouse over ROC curve? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#'
#' @details
#' Calibration curves are a quick and easy way to evaluate a classifier's fit to the
#' data. This function allows one or several models to be plotted in the same figure,
#' with points sized by the number of observations that fall within the corresponding
#' bin.
#'
#' @examples
#' x <- runif(1000)
#' y <- rbinom(1000, size = 1, prob = x)
#' plot_calibration(obs = y, pred = x)
#'
#' x2 <- rbeta(1000, shape1 = 5/2, shape2 = 3/2)
#' plot_calibration(obs = y, pred = list("Good" = x, "Bad" = x2))
#'
#' @export
#' @importFrom purrr map map_dbl map_df
#' @importFrom dplyr data_frame
#' @importFrom ggsci scale_color_d3
#' @import ggplot2
#' @importFrom plotly ggplotly
#'

plot_calibration <- function(obs,
                             pred,
                             main = NULL,
                           legend = 'outside',
                            hover = FALSE) {

  # Preliminaries
  if (is.character(obs)) {
    obs <- as.factor(obs)
  }
  if (is.factor(obs)) {
    if (length(levels(obs)) != 2L) {
      stop('Response must be dichotomous.')
    } else {
      warning('A positive outcome is hereby defined as obs == "', levels(obs)[1], '". ',
              'To change this to obs == "', levels(obs)[2], '", either relevel the ',
              'factor or recode response as numeric (1/0).')
      obs <- ifelse(obs == levels(obs)[1], 1L, 0L)
    }
  }
  if (is.logical(obs)) {
    obs <- ifelse(obs, 1L, 0L)
  }
  if (!all(obs %in% c(0L, 1L))) {
    stop('A numeric response can only take on values of 0 or 1.')
  }
  if (var(obs) == 0L) {
    stop('Response is invariant.')
  }
  if (is.data.frame(pred)) {
    pred <- as.list(pred)
  } else if (!is.list(pred)) {
    pred <- list(pred)
  }
  pred <- map(pred, function(x) x <- x[is.finite(x)])
  if (is.null(names(pred))) {
    names(pred) <- paste0('M', seq_along(pred))
  }
  for (x in seq_along(pred)) {
    if (!is.numeric(pred[[x]])) {
      stop('pred must be a numeric vector, or several such vectors organized into ',
           'a list or data frame.')
    }
    if (max(pred[[x]] > 1L || min(pred[[x]] < 0L))) {
      stop('pred values must be on (0, 1).')
    }
    if (length(obs) != length(pred[[x]])) {
      stop('obs and pred vectors must be of equal length.')
    }
  }
  if (is.null(main)) {
    if (length(pred) == 1L) {
      main <- 'Calibration Curve'
    } else {
      main <- 'Calibration Curves'
    }
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright".')
  }

  # Tidy data
  breaks <- seq(from = 0.05, to = 1L, by = 0.05)
  bin <- map(pred, function(x) {
    map_dbl(seq_along(x), function(i) which.max(x[i] <= breaks))
  })
  exp_grps <- map(seq_along(pred), function(x) {
    split(pred[[x]], bin[[x]])
  })
  obs_grps <- map(seq_along(bin), function(x) split(obs, bin[[x]]))
  df <- map_df(seq_along(pred), function(x) {
    data_frame(Y = map_dbl(obs_grps[[x]], mean),
               X = map_dbl(exp_grps[[x]], mean),
            Freq = map_dbl(exp_grps[[x]], length),
      Classifier = names(pred)[x])
  })

  # Build plot
  p <- ggplot(df, aes(X, Y)) +
    geom_abline(intercept = 0L, slope = 1L, color = 'grey') +
    scale_size(range = c(1L, 5L)) +
    labs(title = main,
         x = 'Expected Probability',
         y = 'Observed Probability') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (length(pred) > 1L) {       # Multiple curves?
    suppressWarnings(
      p <- p + geom_point(aes(size = Freq, color = Classifier, text = Classifier)) +
        geom_path(aes(color = Classifier, text = Classifier)) +
        scale_color_d3()
    )
  } else {
    p <- p + geom_point(aes(size = Freq)) +
      geom_path()
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
    if (legend == 'outside') {
      p <- ggplotly(p, tooltip = 'text', height = 525, width = 600)
    } else {
      p <- ggplotly(p, tooltip = 'text', height = 600, width = 600)
    }
    print(p)
  }

}

# Use gganimate, tweenr, and shiny to toggle btw classifiers
