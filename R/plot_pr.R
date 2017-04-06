#' Precision-Recall Curve(s)
#'
#' This function plots PR curves for one or several classifiers.
#'
#' @param obs Vector of observed outcomes. Must be dichotomous. Can be numeric,
#'   character, factor, or logical. If numeric, \code{obs} must be coded \code{1}
#'   or \code{0}. If character or factor, a warning will be issued clarifying that
#'   the first level is assumed to be the reference.
#' @param pred Vector of predicted values, or several such vectors organized into a
#'   data frame or list, optionally named. Must be numeric. Common examples include
#'   the probabilities output by a logistic model, or the expression levels of a
#'   particular biomarker.
#' @param title Optional plot title.
#' @param leg.txt Optional legend title.
#' @param legend Legend position. Must be one of \code{"outside", "bottomleft",
#'   "bottomright", "topleft",} or \code{"topright"}.
#' @param hover Show predictor name by hovering mouse over PR curve? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#'
#' @details
#' PR curves plot the precision (i.e., positive predictive value) against the recall
#' (i.e., true positive rate/sensitivity) for a given classifier and vector of
#' observations. The area under the PR curve (AUC) is a useful performance metric for
#' binary classifiers, especially in cases of extreme class imbalance, which is
#' typical in omic contexts (Saito & Rehmsmeier, 2015). The grey horizontal line
#' represents the performance of a theoretical random classifier. Interpolations
#' for tied \code{pred} values are computed using the nonlinear method of Davis &
#' Goadrich (2006).
#'
#' @references
#' Davis, J. & Goadrich, M. (2006).
#' \href{http://pages.cs.wisc.edu/~jdavis/davisgoadrichcamera2.pdf}{The Relationship
#' Between Precision-Recall and ROC Curves}. In \emph{Proceedings of the 23rd
#' International Conference on Machine Learning}, pp. 223-240. New York: ACM.
#'
#' Saito, T. & Rehmsmeier, M. (2015).
#' \href{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0118432}{The
#' Precision-Recall Plot Is More Informative than the ROC Plot When Evaluating Binary
#' Classifiers on Imbalanced Datasets}. \emph{PLoS ONE, 10}(3): e0118432.
#'
#' @examples
#' y <- rbinom(1000, size = 1, prob = 0.1)
#' x1 <- rnorm(1000, mean = y)
#' plot_pr(obs = y, pred = x1)
#'
#' x2 <- rnorm(1000, mean = y, sd = 2)
#' plot_pr(obs = y, pred = list("x1" = x1, "x2" = x2))
#'
#' @export
#' @import dplyr
#' @importFrom precrec evalmod
#' @importFrom purrr map_df map_chr
#' @import ggplot2
#' @importFrom ggsci pal_d3
#' @importFrom plotly ggplotly
#'

plot_pr <- function(obs,
                    pred,
                    title = NULL,
                  leg.txt = NULL,
                   legend = 'topright',
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
              'factor or recode response as logical or numeric (1/0).')
      obs <- ifelse(obs == levels(obs)[1], TRUE, FALSE)
    }
  }
  if (is.numeric(obs)) {
    if (!all(obs %in% c(0L, 1L))) {
      stop('A numeric response can only take on values of 0 or 1.')
    } else {
      obs <- ifelse(obs == 1L, TRUE, FALSE)
    }
  }
  if (any(is.na(obs))) {
    stop('obs cannot contain missing values.')
  }
  if (all(obs) | all(!obs)) {
    stop('obs is invariant.')
  }
  if (is.data.frame(pred)) {
    pred <- as.list(pred)
  } else if (!is.list(pred)) {
    pred <- list(pred)
  }
  for (x in seq_along(pred)) {
    if (!is.numeric(pred[[x]])) {
      stop('pred must be a numeric vector, or several such vectors organized into ',
           'a list or data frame.')
    }
    if (length(obs) != length(pred[[x]])) {
      stop('obs and pred vectors must be of equal length.')
    }
    if (any(!is.finite(pred[[x]]))) {
      stop('pred cannot contain missing or infinite values.')
    }
  }
  if (is.null(names(pred))) {
    names(pred) <- paste0('M', seq_along(pred))
  }
  if (is.null(title)) {
    if (length(pred) == 1L) {
      title <- 'Precision-Recall Curve'
    } else {
      title <- 'Precision-Recall Curves'
    }
  }
  if (is.null(leg.txt)) {
    leg.txt <- 'Classifier'
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright".')
  }

  # Tidy data
  prcs <- evalmod(scores = pred, labels = obs)$prcs
  df <- map_df(seq_along(pred), function(m) {
    data_frame(Recall = prcs[[m]]$x,
            Precision = prcs[[m]]$y,
           Classifier = names(pred)[m]) %>%
      return()
  })

  # Build plot
  p_auc <- function(m) {                         # Print AUC
    paste0(names(pred)[m], ', AUC = ', round(attr(prcs[[m]], 'auc'), 2L))
  }
  p <- ggplot(df, aes(Recall, Precision)) +
    lims(x = c(0L, 1L), y = c(0L, 1L)) +
    geom_hline(yintercept = sum(obs) / length(obs), linetype = 2L, color = 'grey') +
    labs(title = title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (length(pred) > 1L) {                       # Multiple curves?
    suppressWarnings(
      p <- p + geom_line(aes(text = Classifier,
                            group = Classifier,
                            color = Classifier)) +
        scale_color_manual(name = leg.txt,
                         labels = map_chr(seq_along(pred), p_auc),
                         values = pal_d3()(length(pred)))
    )
  } else {
    p <- p + geom_line(aes(color = Classifier)) +
      scale_color_manual(name = leg.txt,
                       labels = map_chr(seq_along(pred), p_auc),
                       values = 'black')
  }
  p <- locate_legend(p, legend)

  # Output
  gg_out(p, hover, legend)

}


