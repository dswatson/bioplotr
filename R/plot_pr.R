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
#' typical in omic contexts (Saito & Rehmsmeier, 2015).
#'
#' @references
#' Saito, T. & Rehmsmeier, M. (2015).
#' \href{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0118432}{The
#' Precision-Recall Plot Is More Informative than the ROC Plot When Evaluating Binary
#' Classifiers on Imbalanced Datasets}. \emph{PLoS ONE, 10}(3): e0118432.
#'
#' @examples
#' y <- rbinom(1000, size = 1, prob = 0.1)
#' x <- rnorm(1000, mean = y, sd = 0.5)
#' plot_pr(obs = y, pred = x)
#'
#' x2 <- rnorm(1000, mean = y, sd = 2)
#' plot_pr(obs = y, pred = list("x1" = x, "x2" = x2))
#'
#' @export
#' @import dplyr
#' @importFrom purrr map map_df map_chr
#' @import ggplot2
#' @importFrom PRROC pr.curve
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
    if (length(obs) != length(pred[[x]])) {
      stop('obs and pred vectors must be of equal length.')
    }
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
  df <- map_df(seq_along(pred), function(x) {
    data_frame(Y = obs,
               X = pred[[x]],
      Classifier = names(pred)[x]) %>%
      arrange(desc(X)) %>%
      mutate(TPR = cumsum(Y == 1L) / sum(Y == 1L),
             PPV = cumsum(Y == 1L) / (cumsum(Y == 1L) + cumsum(Y == 0L))) %>%
      return()
  })

  # Build plot
  p_auc <- function(x) {                         # Print AUC
    pos <- df %>% filter(Classifier == names(pred)[x], Y == 1L)
    neg <- df %>% filter(Classifier == names(pred)[x], Y == 0L)
    txt <- paste0(names(pred)[x], ', AUC = ',
                  round(pr.curve(pos$X, neg$X)$auc.integral, 2L))
    return(txt)
  }
  p <- ggplot(df, aes(TPR, PPV)) +
    labs(title = title, x = 'Recall', y = 'Precision') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (length(pred) > 1L) {                       # Multiple curves?
    p <- p + geom_line(aes(text = Classifier,
                          group = Classifier,
                          color = Classifier)) +
      scale_color_manual(name = leg.txt,
                       labels = map_chr(seq_along(pred), p_auc),
                       values = pal_d3()(length(pred)))
  } else {
    p <- p + geom_point(size = 0.1) +
      geom_line(aes(color = Classifier)) +
      scale_color_manual(name = leg.txt,
                       labels = map_chr(seq_along(pred), p_auc),
                       values = 'black')
  }
  if (legend == 'bottomleft') {                  # Locate legend
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
# Slash maybe do cumulative from left to right?
