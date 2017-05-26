#' Reicever Operating Characteristic Curve(s)
#'
#' This functions plots ROC curves for one or several classifiers.
#'
#' @param obs Vector of observed outcomes. Must be dichotomous. Can be logical,
#'   numeric, character, or factor. If numeric, \code{obs} must be coded \code{
#'   1} or \code{0}. If character or factor, a warning will be issued clarifying
#'   that the first level is assumed to be the reference.
#' @param pred Vector of predicted values, or several such vectors organized
#'   into a data frame or list, optionally named. Must be numeric. Common
#'   examples include the probabilities output by a logistic model, or the
#'   expression levels of a particular biomarker.
#' @param pal_curves String specifying the color palette to use when plotting
#'   multiple curves. Options include \code{"ggplot"}, as well as the complete
#'   collection of \code{
#'   \href{https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html}{
#'   ggsci}} palettes, which can be identified by name (e.g., \code{"npg"},
#'   \code{"aaas"}, etc.). Alternatively, a character vector of colors with
#'   length equal to the number of vectors in \code{pred}.
#' @param title Optional plot title.
#' @param leg.txt Optional legend title.
#' @param legend Legend position. Must be one of \code{"right"}, \code{
#'   "left"}, \code{"top"}, \code{"bottom"}, \code{"topright"}, \code{
#'   "topleft"}, \code{"bottomright"}, or \code{"bottomleft"}.
#' @param hover Show predictor name by hovering mouse over ROC curve? If \code{
#'   TRUE}, the plot is rendered in HTML and will either open in your browser's
#'   graphic display or appear in the RStudio viewer.
#'
#' @details
#' ROC curves plot the false positive rate (i.e., 1 - specificity) against the
#' true positive rate (i.e., sensitivity) for a given classifier and vector of
#' observations. The area under the ROC curve (AUC) is a common performance
#' metric for binary classifiers. The grey diagonal line across the plot
#' represents the performance of a theoretical random classifier.
#'
#' @examples
#' y <- rbinom(1000, size = 1, prob = 0.5)
#' x1 <- rnorm(1000, mean = y)
#' plot_roc(obs = y, pred = x1)
#'
#' x2 <- rnorm(1000, mean = y, sd = 2)
#' plot_roc(obs = y, pred = list("Better" = x1, "Worse" = x2))
#'
#' @export
#' @importFrom precrec evalmod
#' @importFrom purrr map_df map_chr
#' @import dplyr
#' @import ggplot2
#'

plot_roc <- function(obs,
                     pred,
                     pal_curves = 'npg',
                          title = NULL,
                        leg.txt = NULL,
                         legend = 'bottomright',
                          hover = FALSE) {

  # Preliminaries
  obs <- format_binom(obs, vec_type = 'obs')
  pred <- format_binom(pred, vec_type = 'pred', n = length(obs))
  if (length(pred) > 1L) {
    cols <- colorize(pal_curves, var_type = 'Categorical', n = length(pred))
  }
  if (is.null(title)) {
    if (length(pred) == 1L) {
      title <- 'ROC Curve'
    } else {
      title <- 'ROC Curves'
    }
  }
  if (is.null(leg.txt)) {
    leg.txt <- 'Classifier'
  }
  if (!legend %in% c('right', 'left', 'top', 'bottom',
                     'topright', 'topleft', 'bottomright', 'bottomleft')) {
    stop('legend must be one of "right", "left", "top", "bottom", ',
         '"topright", "topleft", "bottomright", or "bottomleft".')
  }

  # Tidy data
  rocs <- evalmod(scores = pred, labels = obs)$rocs
  df <- map_df(seq_along(pred), function(m) {
    data_frame(FPR = rocs[[m]]$x,
               TPR = rocs[[m]]$y,
        Classifier = names(pred)[m]) %>%
      return()
  })

  # Build plot
  p_auc <- function(m) {                         # Print AUC
    paste0(names(pred)[m], ', AUC = ', round(attr(rocs[[m]], 'auc'), 2L))
  }
  p <- ggplot(df, aes(FPR, TPR)) +
    lims(x = c(0L, 1L), y = c(0L, 1L)) +
    geom_abline(intercept = 0L, slope = 1L, linetype = 'dashed', color = 'grey') +
    labs(title = title, x = 'False Positive Rate', y = 'True Positive Rate') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (length(pred) > 1L) {                       # Multiple curves?
    suppressWarnings(
      p <- p + geom_line(aes(text = Classifier,
                            color = Classifier)) +
        scale_color_manual(name = leg.txt,
                         breaks = names(pred),
                         labels = map_chr(seq_along(pred), p_auc),
                         values = cols)

    )
  } else {
    p <- p + geom_line(aes(color = Classifier)) +
      scale_color_manual(name = leg.txt,
                       labels = map_chr(seq_along(pred), p_auc),
                       values = 'black')
  }

  # Output
  gg_out(p, hover, legend)

}

# Confidence intervals?
