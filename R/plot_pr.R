#' Precision-Recall Curve(s)
#'
#' This function plots PR curves for one or several classifiers.
#'
#' @param obs Vector of observed outcomes. Must be dichotomous. Can be numeric,
#'   character, factor, or logical. If numeric, \code{obs} must be coded \code{1
#'   } or \code{0}. If character or factor, a warning will be issued clarifying
#'   that the first level is assumed to be the reference.
#' @param pred Vector of predicted values, or several such vectors organized
#'   into a data frame or list, optionally named. Must be numeric. Common
#'   examples include the probabilities output by a logistic model, or the
#'   expression levels of a particular biomarker.
#' @param pal_curves String specifying the color palette to use when plotting
#'   multiple vectors. Options include \code{"ggplot"}, all qualitative color
#'   schemes available in \code{RColorBrewer}, and the complete collection of
#'   \code{\href{http://bit.ly/2bxnuGB}{ggsci}} palettes. Alternatively, a
#'   character vector of colors with length equal to the number of vectors
#'   in \code{pred}.
#' @param title Optional plot title.
#' @param leg.txt Optional legend title.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
#' @param hover Show predictor name by hovering mouse over PR curve? If \code{
#'   TRUE}, the plot is rendered in HTML and will either open in your browser's
#'   graphic display or appear in the RStudio viewer.
#'
#' @details
#' PR curves plot the precision (i.e., positive predictive value) against the
#' recall (i.e., true positive rate/sensitivity) for a given classifier and
#' vector of observations. The area under the PR curve (AUC) is a useful
#' performance metric for binary classifiers, especially in cases of extreme
#' class imbalance, which is typical in omic contexts (Saito & Rehmsmeier,
#' 2015). The grey horizontal line represents the performance of a theoretical
#' random classifier. Interpolations for tied \code{pred} values are computed
#' using the nonlinear method of Davis & Goadrich (2006).
#'
#' @references
#' Davis, J. & Goadrich, M. (2006).
#' \href{http://pages.cs.wisc.edu/~jdavis/davisgoadrichcamera2.pdf}{The
#' Relationship Between Precision-Recall and ROC Curves}. In \emph{Proceedings
#' of the 23rd International Conference on Machine Learning}, pp. 223-240. New
#' York: ACM.
#'
#' Saito, T. & Rehmsmeier, M. (2015). \href{http://bit.ly/2vrw32I}{
#' The Precision-Recall Plot Is More Informative than the ROC Plot When
#' Evaluating Binary Classifiers on Imbalanced Datasets}. \emph{PLoS ONE,
#' 10}(3): e0118432.
#'
#' @examples
#' y <- rbinom(1000, size = 1, prob = 0.1)
#' x1 <- rnorm(1000, mean = y)
#' plot_pr(obs = y, pred = x1)
#'
#' x2 <- rnorm(1000, mean = y, sd = 2)
#' plot_pr(obs = y, pred = list("Better" = x1, "Worse" = x2))
#'
#' @export
#' @importFrom precrec evalmod
#' @importFrom purrr map_df map_chr
#' @import dplyr
#' @import ggplot2
#'

plot_pr <- function(obs,
                    pred,
                    pal_curves = 'npg',
                         title = NULL,
                       leg.txt = NULL,
                        legend = 'topright',
                         hover = FALSE) {

  # Preliminaries
  obs <- format_binom(obs, vec_type = 'obs')
  pred <- format_binom(pred, vec_type = 'pred', n = length(obs))
  if (length(pred) > 1L) {
    cols <- colorize(pal_curves, var_type = 'Categorical', n = length(pred))
  }
  if (title %>% is.null) {
    if (length(pred) == 1L) {
      title <- 'Precision-Recall Curve'
    } else {
      title <- 'Precision-Recall Curves'
    }
  }
  if (leg.txt %>% is.null) {
    leg.txt <- 'Classifier'
  }
  loc <- c('bottom', 'left', 'top', 'right',
           'bottomright', 'bottomleft', 'topleft', 'topright')
  if (!legend %in% loc) {
    stop('legend must be one of ', stringify(loc, 'or'), '.')
  }

  # Tidy data
  prcs <- evalmod(scores = pred, labels = obs)$prcs
  df <- seq_along(pred) %>%
    map_df(~ data_frame(Recall = prcs[[.x]]$x,
                     Precision = prcs[[.x]]$y,
                    Classifier = names(pred)[.x]))

  # Build plot
  p_auc <- function(m) {                         # Print AUC
    auc <- prcs[[m]] %>%
      attr('auc') %>%
      round(2L)
    paste0(names(pred)[m], ', AUPR = ', auc)
  }
  p <- ggplot(df, aes(Recall, Precision)) +
    lims(x = c(0L, 1L), y = c(0L, 1L)) +
    geom_hline(yintercept = sum(obs) / length(obs),
               linetype = 'dashed', color = 'grey') +
    labs(title = title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (length(pred) > 1L) {                       # Multiple curves?
    suppressWarnings(
      p <- p + geom_line(aes(text = Classifier,
                            color = Classifier)) +
        scale_color_manual(name = leg.txt,
                         breaks = names(pred),
                         labels = seq_along(pred) %>% map_chr(p_auc),
                         values = cols)
    )
  } else {
    p <- p + geom_line(aes(color = Classifier)) +
      scale_color_manual(name = leg.txt,
                       labels = seq_along(pred) %>% map_chr(p_auc),
                       values = 'black')
  }

  # Output
  gg_out(p, hover, legend)

}


