#' Reicever Operating Characteristic Curve(s)
#'
#' This functions plots ROC curves for one or several classifiers.
#'
#' @param obs Vector of observed outcomes. Must be dichotomous. Can be numeric,
#'   logical, character, or factor. If numeric, \code{obs} must be coded \code{1}
#'   or \code{0}. If character or factor, a warning will be issued clarifying that
#'   the first level is assumed to be the reference.
#' @param pred Vector of predicted values, or several such vectors organized into a
#'   data frame or list, optionally named. Must be numeric. Common examples include
#'   the probabilities output by a logistic model, or the expression levels of a
#'   particular biomarker.
#' @param main Optional plot title.
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
#' classifiers.
#'
#' @examples
#' y <- rbinom(100, size = 1, prob = 0.5)
#' x <- rnorm(100, mean = y, sd = 0.5)
#' plot_roc(obs = y, pred = x)
#'
#' x2 <- rnorm(100, mean = y, sd = 2)
#' plot_roc(obs = y, pred = list("x1" = x, "x2" = x2))
#'
#' @export
#' @import dplyr
#' @importFrom purrr map map_df map_chr
#' @import ggplot2
#' @importFrom limma auROC
#' @importFrom scales hue_pal
#' @importFrom plotly ggplotly
#'

plot_roc <- function(obs,
                     pred,
                     main = NULL,
                   legend = 'bottomright',
                    hover = FALSE) {

  # Preliminaries
  if (is.character(obs)) obs <- as.factor(obs)
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
  if (is.logical(obs)) obs <- ifelse(obs, 1L, 0L)
  if (!all(obs %in% c(0L, 1L))) {
    stop('A numeric response can only take on values of 0 or 1.')
  }
  if (var(obs) == 0L) {
    stop('Response is invariant.')
  }
  if (is.data.frame(pred)) pred <- as.list(pred)
  else if (!is.list(pred)) pred <- list(pred)
  pred <- map(pred, function(x) x <- x[is.finite(x)])
  if (is.null(names(pred))) {
    names(pred) <- paste0('M', seq_along(pred))
  }
  for (x in seq_along(pred)) {
    if (!is.numeric(pred[[x]])) {
      stop('pred must be a numeric vector, or several such vectors organized into ',
           'a list or data frame.')
    }
    if (length(obs) != length(pred[[x]])) {
      stop('obs and pred vectors must be of equal length.')
    }
  }
  if (is.null(main)) {
    if (length(pred) == 1L) main <- 'ROC Curve'
    else main <- 'ROC Curves'
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright".')
  }

  # Tidy data
  originate <- function(tbl) {
    data_frame(TPR = 0L,
               FPR = 0L,
        Classifier = tbl$Classifier[1]) %>%
      rbind(tbl) %>%
      return()
  }
  df <- map_df(seq_along(pred), function(x) {
    data_frame(Y = obs,
               X = pred[[x]],
      Classifier = names(pred)[x]) %>%
      arrange(desc(X)) %>%
      mutate(TPR = cumsum(Y == 1L) / sum(Y == 1L),
             FPR = cumsum(Y == 0L) / sum(Y == 0L)) %>%
      select(TPR, FPR, Classifier) %>%
      originate() %>%
      return()
  })

  # Build plot
  p_auc <- function(x) {           # Print AUC
    paste0(names(pred)[x], ', AUC = ', round(auROC(obs, pred[[x]]), 2L))
  }
  p <- ggplot(df, aes(FPR, TPR)) +
    geom_abline(intercept = 0L, slope = 1L, color = 'grey') +
    labs(title = main,
             x = 'False Positive Rate',
             y = 'True Positive Rate') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (length(pred) > 1L) {        # Multiple curves?
    suppressWarnings(
      p <- p + geom_step(aes(text = Classifier,
                            group = Classifier,
                            color = Classifier)) +
        scale_colour_manual(name = 'Classifier',
                          labels = map_chr(seq_along(pred), p_auc),
                          values = hue_pal()(length(pred)))
    )
  } else {
    p <- p + geom_step(aes(color = Classifier)) +
      scale_colour_manual(name = 'Classifier',
                        labels = map_chr(seq_along(pred), p_auc),
                        values = 'black')
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


