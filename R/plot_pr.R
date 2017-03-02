#' Precision-Recall Curve(s)
#'
#' This function plots PR curves for one or several classifiers.
#'
#' @param obs Vector of observed outcomes. Must be dichotomous. Can be numeric,
#'   character, factor, or logical. If numeric, \code{obs} must be coded \code{1}
#'   or \code{0}. If character or factor, a warning will be issued clarifying that
#'   the first level is assumed to be the reference.
#' @param pred Vector of predicted values, or several such vectors organized into
#'   a list or data frame. Must be numeric. Common examples include the probabilities
#'   output by a logistic model, or the expression levels of a particular biomarker.
#' @param main Optional plot title.
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
#' binary classifiers, especially when the prevalence of the outcome in question is
#' relatively low.
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
#' @importFrom scales hue_pal
#' @importFrom plotly ggplotly
#'

plot_pr <- function(obs,
                    pred,
                    main = NULL,
                  legend = 'bottomright',
                   hover = FALSE) {

  # Preliminaries
  if (is.character(obs)) {
    obs <- as.factor(obs)
  }
  if (is.factor(obs)) {
    if (length(levels(obs)) != 2L) {
      stop('Response must be dichotomous.')
    } else {
      warning('A positive outcome is hereby defined as obs == "', levels(obs)[1L], '". ',
              'To change this to obs == "', levels(obs)[2L], '", either relevel the ',
              'factor or recode response as numeric (1/0).')
      obs <- ifelse(obs == levels(obs)[1], 1, 0)
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
  if (is.null(names(pred))) {
      names(pred) <- paste0('M', seq_along(pred))
  }
  for (i in seq_along(pred)) {
    if (!is.numeric(pred[[i]])) {
      stop('pred must be a numeric vector, or several such vectors organized into ',
           'a list or data frame.')
    }
    if (length(obs) != length(pred[[i]])) {
      stop('obs and pred vectors must be of equal length.')
    }
  }
  if (is.null(main)) {
    if (length(pred) == 1L) main <- 'Precision-Recall Curve'
    else main <- 'Precision-Recall Curves'
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright".')
  }

  # Tidy data
  pred <- map(pred, function(x) {
    x <- x[is.finite(x)]
  })
  df <- map_df(seq_along(pred), function(i) {
    data_frame(Y = obs,
               X = pred[[i]],
      Classifier = names(pred)[i]) %>%
      arrange(desc(X)) %>%
      mutate(TPR = cumsum(Y == 1L) / sum(Y == 1L),
             PPV = cumsum(Y == 1L) / (cumsum(Y == 1L) + cumsum(Y == 0L))) %>%
      return()
  })

  # Build plot
  leg <- function(i) {                      # Print AUC
    pos <- df %>% filter(Classifier == names(pred)[i], Y == 1L)
    neg <- df %>% filter(Classifier == names(pred)[i], Y == 0L)
    txt <- paste0(names(pred)[i], ', AUC = ',
                  round(pr.curve(pos$X, neg$X)$auc.integral, 2L))
    return(txt)
  }
  p <- ggplot(df, aes(TPR, PPV)) +
    labs(title = main, x = 'Recall', y = 'Precision') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (length(pred) > 1L) {                  # Multiple curves?
    suppressWarnings(
      p <- p + geom_line(aes(text = Classifier,
                            group = Classifier,
                            color = Classifier)) +
        scale_colour_manual(name = 'Classifier',
                          labels = map_chr(seq_along(pred), leg),
                          values = hue_pal()(length(pred)))
    )
  } else {
    p <- p + geom_point(size = 0.1) +
      geom_line(aes(color = Classifier)) +
      scale_colour_manual(name = 'Classifier',
                        labels = map_chr(seq_along(pred), leg),
                        values = 'black')
  }
  if (legend == 'bottomleft') {             # Locate legend
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


