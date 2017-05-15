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
#' @param pal String specifying the color palette to use when plotting multiple
#'   curves. Options include \code{"ggplot"}, as well as the complete collection
#'   of \code{
#'   \href{https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html}{
#'   ggsci}} palettes, which can be identified by name (e.g., \code{"npg"},
#'   \code{"aaas"}, etc.). Alternatively, a character vector of colors with
#'   length equal to the number of vectors in \code{pred}.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"right"}, \code{
#'   "left"}, \code{"top"}, \code{"bottom"}, \code{"bottomright"},
#'   \code{"bottomleft"}, \code{"topright"}, or \code{"topleft"}.
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
#' @importFrom purrr map_lgl map map_dbl map_df
#' @importFrom dplyr data_frame
#' @import ggplot2
#'

plot_calibration <- function(obs,
                             pred,
                             pal = 'd3',
                           title = NULL,
                          legend = 'right',
                           hover = FALSE) {

  # Preliminaries
  obs <- format_binom(obs, vec_type = 'obs')
  pred <- format_binom(pred, vec_type = 'pred', n = length(obs))
  if (any(map_lgl(seq_along(pred), function(m) {
    max(pred[[m]] > 1L | min(pred[[m]] < 0L))
  }))) {
    stop('pred values must be on [0, 1].')
  }
  if (length(pred) > 1L) {
    cols <- colorize(pal, var_type = 'Categorical', n = length(pred))
  }
  if (is.null(title)) {
    if (length(pred) == 1L) {
      title <- 'Calibration Curve'
    } else {
      title <- 'Calibration Curves'
    }
  }
  if (!legend %in% c('right', 'left', 'top', 'bottom',
                     'bottomright', 'bottomleft', 'topright', 'topleft')) {
    stop('legend must be one of "right", "left", "top", "bottom", ',
         '"bottomright", "bottomleft", "topright", or "topleft".')
  }

  # Tidy data
  brks <- seq(from = 0.05, to = 1L, by = 0.05)
  bin <- map(pred, function(m) {
    map_dbl(seq_along(m), function(i) which.max(m[i] <= brks))
  })
  exp_grps <- map(seq_along(pred), function(m) {
    split(pred[[m]], bin[[m]])
  })
  obs_grps <- map(seq_along(bin), function(x) split(obs, bin[[x]]))
  df <- map_df(seq_along(pred), function(m) {
    data_frame(Y = map_dbl(obs_grps[[m]], mean),
               X = map_dbl(exp_grps[[m]], mean),
       Frequency = map_dbl(exp_grps[[m]], length),
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
