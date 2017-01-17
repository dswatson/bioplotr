#' Plot PR curves for one or several classifiers
#'
#' @param obs Vector of observed outcomes. Must be dichotomous. Can be numeric,
#'   logical, character, or factor. If numeric, \code{obs} must be coded \code{1}
#'   or \code{0}. If character or factor, a warning will be issued clarifying that
#'   the first level is assumed to be the reference.
#' @param pred Vector of predicted values or a list of such vectors, optionally named.
#'   Must be numeric. Common examples include the probabilities output by a logistic
#'   model, or the expression levels of a particular biomarker.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside", "bottomleft",
#'   "bottomright", "topleft",} or \code{"topright"}.
#' @param hover Show predictor name by hovering mouse over PR curve? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#'
#' @details
#' This function plots one or several precision-recall (PR) curves. PR curves plot the
#' precision (i.e., positive predictive value) against the recall (i.e., true positive
#' rate/sensitivity) for a given classifier and vector of observations. The area under
#' the PR curve (AUC) is a useful performance metric for binary classifiers, especially
#' when the prevalence of the outcome in question is relatively low.
#'
#' @examples
#' y <- rbinom(300, size = 1, prob = 0.1)
#' x <- rnorm(300, mean = y, sd = 0.5)
#' plot_pr(obs = y, pred = x)
#'
#' x2 <- rnorm(300, mean = y, sd = 2)
#' plot_pr(obs = y, pred = list("x1" = x, "x2" = x2))
#'
#' @export
#' @import dplyr
#' @importFrom purrr map_df map_chr
#' @import ggplot2
#' @importFrom PRROC pr.curve
#' @importFrom scales hue_pal
#' @importFrom plotly ggplotly
#'

plot_pr <- function(obs,
                    pred,
                    main   = NULL,
                    legend = 'outside',
                    hover  = FALSE) {

  # Preliminaries
  if (is.list(pred)) {
    for (i in seq_along(pred)) {
      if (length(obs) != length(pred[[i]])) {
        stop('obs and pred vectors must be of equal length.')
      }
    }
    if (is.null(names(pred))) {
      names(pred) <- paste0('M', seq_along(pred))
    }
  } else {
    if (length(obs) != length(pred)) {
      stop('obs and pred vectors must be of equal length.')
    }
    pred <- list('M1' = pred)
  }
  if (is.character(obs)) {
    obs <- as.factor(obs)
  }
  if (is.factor(obs)) {
    if (length(levels(obs)) > 2) {
      stop('Response must be dichotomous.')
    } else {
      warning('Response vector is character or factor. A positive outcome is hereby ',
              'defined as obs == "', levels(obs)[1], '". To change this to obs == "',
              levels(obs)[2], '", either relevel the factor or recode response as ',
              'numeric (1/0).')
      obs <- ifelse(obs == levels(obs)[1], 1, 0)
    }
  }
  if (is.logical(obs)) {
    obs <- ifelse(obs, 1, 0)
  }
  if (!all(obs %in% c(0, 1))) {
    stop('A numeric response can only take on values of 1 or 0.')
  }
  if (var(obs) == 0) {
    stop('Response is invariant.')
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright", ',
         '"topleft", or "topright".')
  }
  if (is.null(main)) {
    if (length(pred) == 1) {
      main <- 'Precision-Recall Curve'
    } else {
      main <- 'Precision-Recall Curves'
    }
  }

  # Tidy
  rates <- function(i) {
    data_frame(Y = obs,
               X = pred[[i]],
               Classifier = names(pred)[i]) %>%
      arrange(desc(X)) %>%
      mutate(TPR = cumsum(Y == 1) / sum(Y == 1),
             PPV = cumsum(Y == 1) / (cumsum(Y == 1) + cumsum(Y == 0))) %>%
      return()
  }
  df <- map_df(seq_along(pred), rates)

  # Plot
  leg <- function(i) {
    pos <- df %>%
      filter(Classifier == names(pred)[i],
             Y == 1)
    neg <- df %>%
      filter(Classifier == names(pred)[i],
             Y == 0)
    txt <- paste0(names(pred)[i], ', AUC = ',
                  round(pr.curve(pos$X, neg$X)$auc.integral, 2))
    return(txt)
  }
  p <- ggplot(df, aes(TPR, PPV)) +
    labs(title = main,
         x = 'Recall',
         y = 'Precision') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (length(pred) > 1) {
    p <- p + geom_point(aes(color = Classifier), size = 0.1) +
      suppressWarnings(geom_line(aes(text  = Classifier,
                                     group = Classifier,
                                     color = Classifier))) +
      scale_colour_manual(name   = 'Classifier',
                          labels = map_chr(seq_along(pred), leg),
                          values = hue_pal()(length(pred)))
  } else {
    p <- p + geom_point(size = 0.1) +
      geom_line(aes(color = Classifier)) +
      scale_colour_manual(name   = 'Classifier',
                          labels = map_chr(seq_along(pred), leg),
                          values = 'black')
  }

  # Legend location
  if (legend == 'bottomleft') {
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
  if (hover == FALSE) {
    print(p)
  } else {
    p <- ggplotly(p, tooltip = 'text', height = 600, width = 600)
    print(p)
  }

}


