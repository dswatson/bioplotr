#' Reicever Operating Characteristic Curve(s)
#'
#' This functions plots ROC curves for one or several classifiers.
#'
#' @param obs Vector of observed outcomes. Must be dichotomous. Can be logical,
#'   numeric, character, or factor. If numeric, \code{obs} must be coded \code{1}
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
#' @param hover Show predictor name by hovering mouse over ROC curve? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#'
#' @details
#' ROC curves plot the false positive rate (i.e., 1 - specificity) against the true
#' positive rate (i.e., sensitivity) for a given classifier and vector of observations.
#' The area under the ROC curve (AUC) is a common performance metric for binary
#' classifiers. The grey diagonal line across the plot represents the performance of
#' a theoretical random classifier.
#'
#' @examples
#' y <- rbinom(100, size = 1, prob = 0.5)
#' x1 <- rnorm(100, mean = y)
#' plot_roc(obs = y, pred = x1)
#'
#' x2 <- rnorm(100, mean = y, sd = 2)
#' plot_roc(obs = y, pred = list("x1" = x1, "x2" = x2))
#'
#' @export
#' @import dplyr
#' @importFrom precrec evalmod
#' @importFrom purrr map_df map_chr
#' @import ggplot2
#' @importFrom ggsci pal_d3
#' @importFrom plotly ggplotly
#'

plot_roc <- function(obs,
                     pred,
                     title = NULL,
                   leg.txt = NULL,
                    legend = 'bottomright',
                     hover = FALSE) {

  # Preliminaries
  if (is.character(obs)) {
    obs <- as.factor(obs)
  }
  if (is.factor(obs)) {
    if (length(levels(obs)) != 2L) {
      stop('obs must be dichotomous.')
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
  for (m in seq_along(pred)) {
    if (!is.numeric(pred[[m]])) {
      stop('pred must be a numeric vector, or several such vectors organized into ',
           'a list or data frame.')
    }
    if (length(obs) != length(pred[[m]])) {
      stop('obs and pred vectors must be of equal length.')
    }
    if (any(!is.finite(pred[[m]]))) {
      stop('pred cannot contain missing or infinite values.')
    }
  }
  if (is.null(names(pred))) {
    names(pred) <- paste0('M', seq_along(pred))
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
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright".')
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
    geom_abline(intercept = 0L, slope = 1L, linetype = 2L, color = 'grey') +
    labs(title = title, x = 'False Positive Rate', y = 'True Positive Rate') +
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

# Confidence intervals?
