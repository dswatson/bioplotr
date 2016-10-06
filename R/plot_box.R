#' Create box plots by sample
#'
#' @param dat Omic data matrix with rows corresponding to probes and columns
#'   to samples.
#' @param group Optional factor or character vector of length equal to sample size.
#'   Levels are used to color box plots.
#' @param type Optional string specifying omic data type. Currently supports
#'   \code{"microarray", "RNA-seq",} and \code{"methylation"}.
#' @param ylab Label for y-axis. If left \code{NULL}, this defaults to
#'   log expression for \code{type = "microarray"}, log CPM for
#'   \code{type = "RNA-seq"}, and Beta for \code{type = "methylation"}. At least
#'   one of \code{type} or \code{xlab} must be specified.
#' @param main Optional plot title.
#' @param legend Legend position. Must be one of \code{"outside",
#'   "bottomleft", "bottomright", "topleft",} or \code{"topright"}.
#'
#' @details
#' This function displays each sample's distribution of expression or methylation
#' values in a single plot. It is especially helpful when contrasting pre- and
#' post-normalization matrices. It may additionally be used to inspect for batch effects
#' or other associations with phenotypic categories using the \code{group} argument.
#'
#' @examples
#' mat <- matrix(rnorm(5000), nrow = 1000, ncol = 5)
#' plot_box(mat, type = "microarray")
#'
#' library(edgeR)
#' mat <- cbind(matrix(rnbinom(5000, mu = 4, size = 1), nrow = 1000, ncol = 5),
#'              matrix(rnbinom(5000, mu = 4, size = 10), nrow = 1000, ncol = 5))
#' mat <- calcNormFactors(DGEList(mat))
#' mat <- cpm(mat, log = TRUE, prior.count = 0.5)
#' batch <- rep(c("A", "B"), each = 5)
#' plot_box(mat, group = batch, type = "RNA-seq")
#'
#' @export
#' @importFrom tidyr gather
#' @import dplyr
#' @import ggplot2
#'

plot_box <- function(dat,
                     group  = NULL,
                     type   = NULL,
                     ylab   = NULL,
                     main   = NULL,
                     legend = 'outside') {

  df <- gather(as_tibble(dat), Sample, Expression) %>%
    mutate(Sample = factor(Sample, levels = unique(Sample)))

  if (!is.null(group)) {
    df <- mutate(df, Group = rep(group, each = nrow(dat)))
  }
  if (is.null(main) & is.null(group)) {
    main <- 'Expression By Sample'
  }
  if (is.null(main) & !is.null(group)) {
    main <- 'Expresion By Group'
  }
  if (is.null(type) & is.null(ylab)) {
    stop('Either data type or ylab must be provided')
  }
  if (type == 'microarray') {
    ylab <- expression('log'[2]*' Expression')
  }
  if (type == 'RNA-seq') {
    ylab <- expression('log'[2]*' Counts Per Million')
  }
  if (type == 'methylation') {
    ylab <- 'Beta'
  }

  p <- ggplot(df, aes(Sample, Expression)) +
    labs(title = main, y = ylab) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (!is.null(group)) {
    p <- p + geom_boxplot(aes(fill = Group))
  } else {
    p <- p + geom_boxplot()
  }

  if (legend == 'bottomleft') {
    p <- p + theme(legend.justification = c(0, 0), legend.position = c(0, 0))
  } else if (legend == 'bottomright') {
    p <- p + theme(legend.justification = c(1, 0), legend.position = c(1, 0))
  } else if (legend == 'topleft') {
    p <- p + theme(legend.justification = c(0, 1), legend.position = c(0, 1))
  } else if (legend == 'topright') {
    p <- p + theme(legend.justification = c(1, 1), legend.position = c(1, 1))
  }

  print(p)

}


