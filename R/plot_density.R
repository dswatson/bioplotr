#' Create density plots by sample
#'
#' @param dat Omic data matrix with rows corresponding to probes and columns
#'   to samples.
#' @param group Optional character or factor vector of length equal to sample size.
#'   Levels are used to color density curves. If supplied, legend title defaults to
#'   "Group". Override this feature by passing a named list instead.
#' @param xlab Optional label for x-axis.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside",
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
#' @param hover Show sample name by hovering mouse over data point? If \code{TRUE},
#'   the plot is rendered in HTML and will either open in your browser's graphic
#'   display or appear in the RStudio viewer.
#'
#' @details
#' This function displays each sample's omic data distribution as a density curve.
#' It is especially helpful when contrasting pre- and post-normalization matrices.
#' It may additionally be used to inspect for batch effects or other associations
#' with phenotypic features by using the \code{group} argument.
#'
#' @examples
#' mat <- matrix(rnorm(5000), nrow = 1000, ncol = 5)
#' plot_density(mat, xlab = "Normalized Expression")
#'
#' library(DESeq2)
#' mat <- cbind(matrix(rnbinom(5000, mu = 4, size = 1), nrow = 1000, ncol = 5),
#'              matrix(rnbinom(5000, mu = 4, size = 10), nrow = 1000, ncol = 5))
#' mat <- rlog(mat)
#' batch <- gl(n = 2, k = 5, labels = c("A", "B"))
#' plot_density(mat, group = batch, xlab = "Normalized Counts")
#'
#' @export
#' @importFrom tidyr gather
#' @import dplyr
#' @import ggplot2
#' @importFrom plotly ggplotly
#'

plot_density <- function(dat,
                         group  = NULL,
                         xlab   = NULL,
                         main   = NULL,
                         legend = 'outside',
                         hover  = FALSE) {

  # Preliminaries
  if (is.null(group)) {
    group <- list(rep(1, times = ncol(dat)))
  } else {
    if (!is.list(group)) {
      group <- list(group)
    }
    if (!is.character(group[[1]]) & !is.factor(group[[1]])) {
      stop('group must be a character or factor variable')
    }
    if (length(group) > 1) {
      stop('group cannot be a list of length > 1')
    }
    if (length(group[[1]]) != ncol(dat)) {
      stop('group length must match number of samples in dat')
    }
    if (length(unique(group[[1]])) == 1) {
      warning('group is invariant')
    }
  }
  if (is.null(xlab)) {
    xlab <- 'Value'
  } else {
    if (!is.character(xlab) | length(xlab) != 1) {
      stop('xlab must be NULL or a single string')
    }
  }
  if (is.null(main)) {
    if (is.numeric(group[[1]])) {
      main <- 'Density By Sample'
    } else {
      if (is.null(names(group))) {
        main <- 'Density By Group'
      } else {
        main <- paste('Density By', names(group))
      }
    }
  }
  if (!legend %in% c('outside', 'bottomleft', 'bottomright', 'topleft', 'topright')) {
    stop('legend must be one of "outside", "bottomleft", "bottomright" ',
         '"topleft", or "topright"')
  }

  # Tidy
  df <- gather(as_data_frame(dat), Sample, Value) %>%
    mutate(Group = rep(group[[1]], each = nrow(dat)))

  # Basic plot
  p <- ggplot(df, aes(Value, group = Sample)) +
    labs(title = main,
         x = xlab,
         y = 'Density') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (!is.numeric(group)) {
    p <- p + suppressWarnings(geom_path(stat = 'density',
                                        aes(text = Sample, color = Group)))
  } else {
    p <- p + geom_path(stat = 'density', aes(color = Sample))
  }

  # Named list?
  if (!is.null(names(group))) {
    p <- p + guides(color = guide_legend(title = names(group)))
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
    p <- ggplotly(p, tooltip = 'text', height = 600, width = 650)
    print(p)
  }

}


